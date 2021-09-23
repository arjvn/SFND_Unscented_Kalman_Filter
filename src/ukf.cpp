#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include <fstream>


using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // State dimension
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = 7;

  // Sigma point spreading parameter
  lambda_ = 3.0 - n_aug_;

  // Measurement space lidar and radar state size
  n_z_lidar_ = 2; // px and py
  n_z_radar_ = 3; // rho, phi and rhodot

  // initial state vector
  x_ = VectorXd(n_x_);
  x_.fill(0.0);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);
  P_.fill(0.0);
  P_ << 4, 0, 0,  0,  0,
      0, 4, 0,  0,  0,
      0, 0, 25, 0,  0,
      0, 0, 0,  10, 0,
      0, 0, 0,  0,  0.04;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3.0;
  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.4;

  // Augmented Sigma matrix of size 7 x 14
  Xsig_aug_ = Eigen::MatrixXd(n_aug_, 2*n_aug_ + 1);
  Xsig_aug_.fill(0.0);

  // Sigma matrix of size 5 x 14 which stores all of the sigma points in the predicted state
  Xsig_pred_ = Eigen::MatrixXd(n_x_, 2*n_aug_ + 1);
  Xsig_pred_.fill(0.0);

  // weights_ of sigma points
  weights_ = VectorXd(2*n_aug_+1);

  time_us_ = 0;

  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

  // LiDAR measurement noise
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  //Radar measurement noise omega k (wk)
  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
   * End DO NOT MODIFY section for measurement noise values
   */

  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */

  Lidar_sensor_noise = MatrixXd(n_z_lidar_, n_z_lidar_);
  Lidar_sensor_noise << (std_laspx_*std_laspx_), 0,
              0, (std_laspy_*std_laspy_);
  H_lidar_ = MatrixXd(n_z_lidar_, n_x_);
  H_lidar_ << 1, 0, 0, 0, 0,
              0, 1, 0, 0, 0;

  Radar_SensorNoise_ = MatrixXd(n_z_radar_, n_z_radar_);
  Radar_SensorNoise_ << std_radr_*std_radr_, 0, 0,
                        0, std_radphi_*std_radphi_, 0,
                        0, 0,std_radrd_*std_radrd_;

  std::stringstream ss;
  ss << "UKF_log" << ".txt";
  log_file_ = std::make_shared<std::ofstream>(ss.str(), std::ios::out);
  *log_file_ << "--------------- UKF Log ---------------"<< '\n'<< '\n' << '\n';

  std::stringstream nis_l;
  nis_l << "nis_l" << ".txt";
  nis_l_ = std::make_shared<std::ofstream>(nis_l.str(), std::ios::out);

  std::stringstream nis_r;
  nis_r << "nis_r" << ".txt";
  nis_r_ = std::make_shared<std::ofstream>(nis_r.str(), std::ios::out);
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

  // --------------- UKF Log --------------- //
  // std::ofstream log;
  // log.open ("UKF_log.txt");
  // --------------------------------------- //

  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.

   * Tools.cpp calls processmeasurement - hence used to call prediction and update steps
   */

   // * Prep data for initlisation
   // * Create switch case for Radar and LiDAR
  if (!is_initialized_) {
    double x;
    double y;
    switch (meas_package.sensor_type_) {
     case MeasurementPackage::LASER:
      *log_file_ << "--------------- Initilisation with Lidar Data ---------------" << '\n';
      *log_file_ << "LiDAR data input: " << meas_package.timestamp_ << '\n' << '\n';
      // initilise lidar data for the first time
      x = meas_package.raw_measurements_(0);
      y = meas_package.raw_measurements_(1);
      x_ << x, y, 0.0, 0.0, 0.0;
      *log_file_ << "Data: " << x <<", " << y << '\n';
      break;
     case MeasurementPackage::RADAR:
      *log_file_ << "--------------- Initilisation with Radar Data ---------------" << '\n';
      *log_file_ << "Radar data input: " << meas_package.timestamp_ << '\n' << '\n';
      double range = meas_package.raw_measurements_(0);
      double azimuth = meas_package.raw_measurements_(1);
      // do we need these
      double speed = meas_package.raw_measurements_(2);
      double x_radar = range * std::cos(azimuth);
      double y_radar = range * std::sin(azimuth);
      double vx = speed * std::cos(azimuth);
      double vy = speed * std::sin(azimuth);
      double v = std::sqrt(vx * vx + vy * vy);
      x_ << x_radar, y_radar, v, range, speed;
      *log_file_ << "Data: " << range << ", " << azimuth << '\n';
      *log_file_ << "In cartisian: " << x_radar << ", " << y_radar << '\n';
      break;
    }
  is_initialized_ = true;
  return;
  }
  else {
    double delta_t = (meas_package.timestamp_ - time_us_) * 1.0e-6; // in seconds for prediction

    Prediction(delta_t);
    time_us_ = meas_package.timestamp_;
    switch(meas_package.sensor_type_) {
      case MeasurementPackage::LASER:
        UpdateLidar(meas_package);
        break;
      case MeasurementPackage::RADAR:
        UpdateRadar(meas_package);
        break;
    }
  }
}

void UKF::Prediction(double delta_t) {
  /**
  * TODO: Complete this function! Estimate the object's location.
  * Modify the state vector, x_. Predict sigma points, the state,
  * and the state covariance matrix.
  */

  // initialise augmented sigma points (7 x 15)
  calcAugmentedSigmaPoints();

  // get sigma points of predicted space (5 x 15)
  // What is delta t here??
  calcSigmaPointsPredictedSpace(delta_t);

  // predict mean and covariance from sigma points
  calcMeanCovarSigmaPredictedSpace();

}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
/**
 * TODO: Complete this function! Use lidar data to update the belief
 * about the object's position. Modify the state vector, x_, and
 * covariance, P_.
 * You can also calculate the lidar NIS, if desired.
 */

  // Predicted (aka expected) measurement.
  VectorXd z_pred = H_lidar_ * x_;

  // Actual measurement.
  VectorXd z = VectorXd(n_z_lidar_);
  z(0) = meas_package.raw_measurements_(0);
  z(1) = meas_package.raw_measurements_(1);

  // Innovation (no need to normalize any angle, as we only have position data here).
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_lidar_.transpose();
  MatrixXd S = H_lidar_ * P_ * Ht + Lidar_sensor_noise;
  MatrixXd S_inverse = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * S_inverse;

  //new estimate
  x_ = x_ + (K * y);
  x_(3) = normalizeAngle(x_(3));
  MatrixXd I = MatrixXd::Identity(n_x_, n_x_);
  P_ = (I - K * H_lidar_) * P_;

  // Compute and export the NIS (normalized innovation squared) value.
  double nis = y.transpose() * S_inverse * y;

  // *log_file_ << "--------------- Update with LiDAR Log ---------------" << '\n';
  // *log_file_ << "LiDAR data input: " << meas_package.timestamp_ << '\n';
  // *log_file_ << "NIS LiDAR: " << nis << '\n';
  // *log_file_ << nis << '\n';

}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief
   * about the object's position. Modify the state vector, x_, and
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */

   // linear error hence no need to generate augmented sigma points again.
   // transform the existing augmented sigma points into the measurement space
   // the new matrix that holds the sigma points in the measurement space is calligrafic Z
   // the measurement space of the radar is 3 dimensional (rho, phi and rhodot)
   // Hence calligrafic Z is a 3 x 15 vector

  MatrixXd Zsig = MatrixXd(n_z_radar_, 2 * n_aug_ + 1);
  Zsig.fill(0.0);
  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_radar_);
  z_pred.fill(0.0);
  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z_radar_,n_z_radar_);
  S.fill(0.0);

  // transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
     Zsig(0, i) = std::pow(std::pow(Xsig_pred_(0, i), 2) + std::pow(Xsig_pred_(1, i), 2), 0.5);
     Zsig(1, i) = atan2(Xsig_pred_(1, i), Xsig_pred_(0, i));
     Zsig(2, i) = (Xsig_pred_(0, i)*std::cos(Xsig_pred_(3, i))*Xsig_pred_(2, i) + Xsig_pred_(1, i)*std::sin(Xsig_pred_(3, i))*Xsig_pred_(2, i)) / std::pow((std::pow(Xsig_pred_(0, i), 2) + std::pow(Xsig_pred_(1, i), 2)), 0.5);
  }

  // Next we need to calculate the mean and covariance of the sigma points in the measurement space
  // Remember to account for the measurement covariance noise R | we do this instead of augmentation
  // because the noise does not have a non-linear effect on the measurement, its purely additive

  // calculate mean predicted measurement
  for (int i=0; i < 2*n_aug_+1; ++i) {
   z_pred += weights_(i) * Zsig.col(i);
  }

  // calculate innovation covariance matrix S
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // 2n+1 simga points
   // residual
   VectorXd z_diff = Zsig.col(i) - z_pred;

   // angle normalization
   while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
   while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

   S += weights_(i) * z_diff * z_diff.transpose();
  }

  S = S + Radar_SensorNoise_;

  // std::cout << "S UpdateRadar" << '\n';
  // std::cout << S << '\n';


  // Calculate the cross correlation matrix.
  MatrixXd Tc = MatrixXd(n_x_, n_z_radar_);
  Tc.fill(0.0);
  for (int ii=0; ii<2*n_aug_ + 1; ii++) {
    VectorXd x_dev = Xsig_pred_.col(ii) - x_;
    x_dev(3) = normalizeAngle(x_dev(3));
    VectorXd z_dev = Zsig.col(ii) - z_pred;
    z_dev(1) = normalizeAngle(z_dev(1));
    Tc += weights_(ii) * (x_dev * z_dev.transpose());
  }

  // Calculate the Kalman gain K.
  MatrixXd K = Tc * S.inverse();

  // Create a vector to hold the actual measurement.
  VectorXd z = VectorXd(n_z_radar_);
  z(0) = meas_package.raw_measurements_(0);
  z(1) = meas_package.raw_measurements_(1);
  z(2) = meas_package.raw_measurements_(2);

  // Update the state mean x and the covariance matrix P.
  VectorXd z_diff = z - z_pred;
  z_diff(1) = normalizeAngle(z_diff(1));
  x_ += K*z_diff;
  x_(3) = normalizeAngle(x_(3));
  P_ -= K*S*K.transpose();

  // Compute and export the NIS (normalized innovation squared) value.
  double nis = z_diff.transpose() * S.inverse() * z_diff;

  *log_file_ << "--------------- Update with Radar Log ---------------" << '\n';
  *log_file_ << "Radar data input: " << meas_package.timestamp_ << '\n';
  // *log_file_ << "NIS Radar: " << nis << '\n';
  *log_file_ << nis << '\n';

  // std::cout << "x_ UpdateRadar" << '\n';
  // std::cout << x_ << '\n';
  // std::cout << "P_ UpdateRadar" << '\n';
  // std::cout << P_ << '\n';
}

void UKF::calcAugmentedSigmaPoints() {
  lambda_ = 3 - n_aug_;

  // create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_; // sort initial P out pls

  // Process noise covariance matrix Q
  P_aug(5, 5) = std::pow(std_a_, 2);
  P_aug(6,6) = (std_yawdd_);

  // create augmented state vector
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.head(5) = x_;
  x_aug(5) = 0.0;
  x_aug(6) = 0.0;

  // create sigma point matrix
  Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  MatrixXd A_ = P_aug.llt().matrixL();
  // std::cout << "x_aug" << "\n" << x_aug << '\n';
  // std::cout << "P_aug" << "\n" << P_aug << '\n';
  // int pause;
  // std::cin >> pause; //used to pause code
  // std::cout << "A_aug" << "\n" << A_ << '\n';
  A_ = A_ * std::sqrt(lambda_ + n_aug_);

  // Create augmented sigma points for the state.
  Xsig_aug_.fill(0.0);
  for (int cc=0; cc<n_aug_; cc++) {
      Xsig_aug_.col(cc+1) = A_.col(cc);
      Xsig_aug_.col(cc+1+n_aug_) = A_.col(cc) * (-1.0);
  }
  for (int cc=0; cc<2*n_aug_ + 1; cc++) {
      Xsig_aug_.col(cc) += x_aug;
  }
}

void UKF::calcSigmaPointsPredictedSpace(double delta_t) {

  for (int i = 0; i < Xsig_aug_.cols(); i++) {
    VectorXd x_k = VectorXd(5);
    VectorXd term1 = VectorXd(5);
    VectorXd term2 = VectorXd(5);

    x_k(0) = Xsig_aug_(0, i);
    x_k(1) = Xsig_aug_(1, i);
    x_k(2) = Xsig_aug_(2, i);
    x_k(3) = Xsig_aug_(3, i);
    x_k(4) = Xsig_aug_(4, i);

    // Term 1: deterministic portion of model
    // Term 2: stocastic portion of model
    if (Xsig_aug_(4, i) < 0.001) {
      term1(0) = Xsig_aug_(2, i)*(std::cos(Xsig_aug_(3, i)))*delta_t;
      term1(1) = Xsig_aug_(2, i)*(std::sin(Xsig_aug_(3, i)))*delta_t;

    }
    else {
      term1(0) = Xsig_aug_(2, i)/Xsig_aug_(4, i)*(std::sin(Xsig_aug_(3, i) + Xsig_aug_(4, i)*delta_t) - std::sin(Xsig_aug_(3, i)));
      term1(1) = Xsig_aug_(2, i)/Xsig_aug_(4, i)*(-std::cos(Xsig_aug_(3, i) + Xsig_aug_(4, i)*delta_t) + std::cos(Xsig_aug_(3, i)));
      term1(3) = Xsig_aug_(4, i) * delta_t;
    }
    term1(2) = 0;
    term1(3) = Xsig_aug_(4, i) * delta_t;
    term1(4) = 0;

    term2(0) = 0.5*std::pow(delta_t, 2)*std::cos(Xsig_aug_(3, i))*Xsig_aug_(5, i);
    term2(1) = 0.5*std::pow(delta_t, 2)*std::sin(Xsig_aug_(3, i))*Xsig_aug_(5, i);
    term2(2) = delta_t*Xsig_aug_(5, i);
    term2(3) = 0.5*std::pow(delta_t, 2)*Xsig_aug_(6, i);
    term2(4) = delta_t*Xsig_aug_(6, i);
    Xsig_pred_.col(i) = x_k + term1 + term2;
    // std::cout << "x_k: " << '\n' << x_k << '\n';
    // std::cout << "term1: " << '\n' << term1 << '\n';
    // std::cout << "term2: " << '\n' << term2 << '\n';
    // std::cout << "Xsig_pred_: " << '\n' << Xsig_pred_ << '\n';
  }
  // std::cout << "Xsig_pred_ = " << std::endl << Xsig_pred_ << std::endl;
  // non zero rhodot

}

void UKF::calcMeanCovarSigmaPredictedSpace() {
  // predict mean
  x_.fill(0.0);
  weights_(0) = lambda_/(lambda_+n_aug_);
  for (int i = 1; i < 2*n_aug_+1; i++) {
    weights_(i) = 0.5/(lambda_+n_aug_);
  }
  for (int i = 0; i < 2*n_aug_+1; i++) {
    x_ = x_ + weights_(i)*Xsig_pred_.col(i);
  }
  // predict covar
  P_.fill(0.0); // important to clear previous values
  for (int i = 0; i < 2*n_aug_+1; i++) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    x_diff(3) = normalizeAngle(x_diff(3));
    P_ += weights_(i) * x_diff * x_diff.transpose();
  }

  // std::cout << "x_ calcMeanCovarSigmaPredictedSpace" << '\n';
  // std::cout << x_ << '\n';
  // std::cout << "P_ calcMeanCovarSigmaPredictedSpace" << '\n';
  // std::cout << P_ << '\n';
}

double UKF::normalizeAngle(double angle) const {
  constexpr double two_pi = 2.0 * M_PI;
  while (angle > M_PI) {
    angle -= two_pi;
  }
  while (angle < -M_PI) {
    angle += two_pi;
  }

  return angle;
}
