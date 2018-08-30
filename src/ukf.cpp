#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

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
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;
  
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
    }
    else // LIDAR measurement
    {
      x_ << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1), 0, 0, 0;
    }

    is_initialized_ = true;
    time_us_ = meas_package.timestamp_;
    return;
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
      px_pred = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
      py_pred = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
      px_pred = p_x + v*delta_t*cos(yaw);
      py_pred = p_y + v*delta_t*sin(yaw);
    }
    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_pred = px_pred + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_pred = py_pred + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

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
  for (int i(0); i < aug_size_ + 1; i++) {
    x_ = x_+ weights_(i) * Xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i(0); i < aug_size_; i++) {  //iterate over sigma points
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    if(fabs(x_diff(3)) > M_PI)
    {
      x_diff(3) -= round(x_diff(3)/(2.0 * M_PI)) * (2.0 * M_PI);
    }
    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }
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
}
