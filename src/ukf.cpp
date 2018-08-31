#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include <fstream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;
  n_x_ = 5;
  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = M_PI/12.0;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  // augmented state dimension
  n_aug_ = 7;
  aug_size_ = 2 * n_aug_ + 1;
  // Weights of sigma points
  weights_ = VectorXd(aug_size_);
  // Sigma point spreading parameter
  lambda_ = 3 - n_aug_;
  // Predicted sigma points
  Xsig_pred_ = MatrixXd(n_x_, aug_size_);

  // State covariance matrix (start with identity, then will see)
  P_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 1, 0,
        0, 0, 0, 0, 1;

  H_laser_ = MatrixXd(2, n_x_);
  H_laser_ << 1, 0, 0, 0, 0,
              0, 1, 0, 0, 0;
  // Measurement covariance matrix R - laser
  R_laser_ = MatrixXd(2,2);
  R_laser_ << std_laspx_ * std_laspx_, 0,
              0, std_laspy_ * std_laspy_;

  //laser_file_ = "NIS_laser.txt";

  NIS_laser_file_.open(laser_file_, std::ofstream::out | std::ofstream::trunc);
  NIS_radar_file_.open(radar_file_, std::ofstream::out | std::ofstream::trunc);
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  // Initialize the state x_ with first measurement.
  if(!is_initialized_)
  {
    cout << "UKF: " << endl;

    x_ << 0, 0, 0, 0, 0;
    if(meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
      cout << "Initialization with radar measurement.." << endl;
      float ro = meas_package.raw_measurements_(0);
      float phi = meas_package.raw_measurements_(1);
      float ro_dot = meas_package.raw_measurements_(2);
      if(fabs(phi) > M_PI)
      {
        phi -= round(phi / (2 * M_PI)) * (2 * M_PI);
      }
      float p_x = ro * cos(phi);
      float p_y = ro * sin(phi);
      float v_x = ro_dot * cos(phi);
      float v_y = ro_dot * sin(phi);
      float v = sqrt(v_x * v_x + v_y * v_y);

      x_ << p_x, p_y, v, phi, 0;
      cout << "Initialized!" << endl;
    }
    else // LIDAR measurement
    {
      cout << "Initialization with laser measurement.." << endl;
      x_ << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1), 0, 0, 0;
      cout << "Initialized!" << endl;
    }

    is_initialized_ = true;
    time_us_ = meas_package.timestamp_;
    //return;
  }
  // Calculate delta t
  double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
  // Predict state
  Prediction(dt);
  // Measurement update
  if(meas_package.sensor_type_ == MeasurementPackage::RADAR)
  {
    UpdateRadar(meas_package);
  }
  else
  {
    UpdateLidar(meas_package);
  }
  time_us_ = meas_package.timestamp_;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  // Augmented mean state
  //cout << "Prediction started.." << endl;
  VectorXd x_aug = VectorXd(n_aug_);

  // Augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  // Sigma points
  MatrixXd Xsig_aug = MatrixXd(n_aug_, aug_size_);

  // initialize augmented mean state
  x_aug.head(n_x_) = x_;
  x_aug(n_x_) = 0;
  x_aug(n_x_ + 1) = 0;

  // initialize augmented state covariance
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(n_x_,n_x_) = std_a_ * std_a_;
  P_aug(n_x_ + 1, n_x_ + 1) = std_yawdd_ * std_yawdd_;
  // create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  // create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i(0); i< n_aug_; ++i)
  {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }
  // predict sigma points
  for (int i(0); i < aug_size_; ++i)
  {
    double p_x = Xsig_aug(0, i);
    double p_y = Xsig_aug(1, i);
    double v = Xsig_aug(2, i);
    double yaw = Xsig_aug(3, i);
    double yawd = Xsig_aug(4, i);
    double nu_a = Xsig_aug(5, i);
    double nu_yawdd = Xsig_aug(6, i);

    double px_pred;
    double py_pred;

    //avoid division by zero
    if (fabs(yawd) > 0.0001) {
      px_pred = p_x + v/yawd * ( sin (yaw + yawd * delta_t) - sin(yaw));
      py_pred = p_y + v/yawd * ( cos(yaw) - cos(yaw + yawd * delta_t) );
    }
    else {
      px_pred = p_x + v * delta_t * cos(yaw);
      py_pred = p_y + v * delta_t * sin(yaw);
    }
    double v_p = v;
    double yaw_p = yaw + yawd * delta_t;
    double yawd_p = yawd;

    //add noise
    px_pred = px_pred + 0.5 * nu_a * delta_t * delta_t * cos(yaw);
    py_pred = py_pred + 0.5 * nu_a * delta_t * delta_t * sin(yaw);
    v_p = v_p + nu_a * delta_t;

    yaw_p = yaw_p + 0.5 * nu_yawdd * delta_t * delta_t;
    yawd_p = yawd_p + nu_yawdd * delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_pred;
    Xsig_pred_(1,i) = py_pred;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }
  // set weights
  weights_.fill(0.5/(n_aug_+lambda_));
  weights_(0) = lambda_/(lambda_+n_aug_);

  //predicted state mean
  x_.fill(0.0);
  for (int i(0); i < aug_size_; ++i) {
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i(0); i < aug_size_; ++i) {  //iterate over sigma points
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    if(fabs(x_diff(3)) > M_PI)
    {
      x_diff(3) -= round(x_diff(3)/(2.0 * M_PI)) * (2.0 * M_PI);
    }
    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }
  //cout << "Prediction finished!" << endl;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  // measurement dimension
  //cout << "Laser update started.." << endl;
  VectorXd z_pred = H_laser_ * x_;
  VectorXd y = meas_package.raw_measurements_ - z_pred;
  MatrixXd Ht = H_laser_.transpose();
  MatrixXd S = H_laser_ * P_ * Ht + R_laser_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;
  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_laser_) * P_;
  // calculate NIS and append to file
  double NIS_laser = y.transpose()*Si*y;

  // write NIS value to file, will analyse it later with python
  ofstream NIS_laser_file;
  NIS_laser_file.open (laser_file_, ios::app);
  //auto str = std::to_string(NIS_laser);
  NIS_laser_file << NIS_laser << ";";
  NIS_laser_file.close();
  //cout << "Laser update finished!" << endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  // measurement dimension
  //cout << "Radar update started.." << endl;
  int n_z_= 3;
  // measurements
  VectorXd z = VectorXd(n_z_);
  z << meas_package.raw_measurements_(0),
       meas_package.raw_measurements_(1),
       meas_package.raw_measurements_(2);

  // Sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z_, aug_size_);

  // transform sigma points into measurement space
  for (int i(0); i < aug_size_; ++i) {  //2n+1 simga points

    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw) * v;
    double v2 = sin(yaw) * v;

    // measurement model
    Zsig(0,i) = sqrt(p_x * p_x + p_y * p_y);
    double phi = atan2(p_y, p_x);
    /*
    if(fabs(phi) > M_PI)
    {
      phi -= round(phi/(2.0 * M_PI)) * (2.0 * M_PI);
    }*/
    Zsig(1,i) = phi;
    Zsig(2,i) = (p_x * v1 + p_y * v2 ) / sqrt(p_x * p_x + p_y * p_y);
  }

  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_);
  z_pred.fill(0.0);
  for (int i(0); i < aug_size_; ++i) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z_, n_z_);
  S.fill(0.0);
  for (int i(0); i < aug_size_; ++i) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    if(fabs(z_diff(1)) > M_PI)
    {
      z_diff(1) -= round(z_diff(1)/(2.0 * M_PI)) * (2.0 * M_PI);
    }

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z_,n_z_);
  R <<    std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0,std_radrd_*std_radrd_;

  S = S + R;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i(0); i < aug_size_; ++i) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    if(fabs(z_diff(1)) > M_PI)
    {
      z_diff(1) -= round(z_diff(1)/(2.0 * M_PI)) * (2.0 * M_PI);
    }

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    if(fabs(x_diff(3)) > M_PI)
    {
      x_diff(3) -= round(x_diff(3)/(2.0 * M_PI)) * (2.0 * M_PI);
    }

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  if(fabs(z_diff(1)) > M_PI)
  {
    z_diff(1) -= round(z_diff(1)/(2.0 * M_PI)) * (2.0 * M_PI);
  }

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

  //NIS calculation for Radar
  double NIS_Radar = z_diff.transpose()*S.inverse()*z_diff;
  ofstream NIS_radar_file;
  NIS_radar_file.open (radar_file_, ios::app);
  NIS_radar_file << NIS_Radar << ";";
  NIS_radar_file.close();
  //cout << "Radar update finished!" << endl;
}
