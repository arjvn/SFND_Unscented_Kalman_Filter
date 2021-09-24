#ifndef UKF_H
#define UKF_H

#include <fstream>
#include <string>
#include <memory>

#include "Eigen/Dense"
#include "measurement_package.h"

class UKF {
 public:
  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(MeasurementPackage meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(MeasurementPackage meas_package);

  /**
     * \brief Compute augmented sigma points based on the mean state x and the state error covariance P.
     * \details The resulting sigma points are stored in Xsig_aug_ a 7 by 15 matrix.
     */
  void calcAugmentedSigmaPoints();

  /**
    * \brief Predict the sigma points using the state transition model.
    * \details moves the sigma points into the predicted state space which are now stored in Xsig_pred_ a 5 by 15 matrix.
    *           handles 2 cases: zero rhodot and non zero rho_dot
    *           model is made up of two portions stocastic and determinsitic.
    *           only 3 terms of the stocastic model differ between zero and nonzero rhodot
    *           refer to lesson 4 chapter 20 for more detail and equations
    */
  void calcSigmaPointsPredictedSpace(double delta_t);

  void calcMeanCovarSigmaPredictedSpace();

  double normalizeAngle(double angle) const;

  // initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  // if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  // if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  // state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  Eigen::VectorXd x_;

  // state covariance matrix
  Eigen::MatrixXd P_;

  // Augmented sigma points matrix
  Eigen::MatrixXd Xsig_aug_;

  // predicted sigma points matrix
  Eigen::MatrixXd Xsig_pred_;

  // Radar measuement noise matrix R_radar
  Eigen::MatrixXd R_Lidar_SensorNoise_;

  // Lidar measurement matrix H_lidar
  Eigen::MatrixXd H_lidar_;

  // Radar measuement noise matrix R_radar
  Eigen::MatrixXd R_Radar_SensorNoise_;

  // time when the state is true, in us
  long long time_us_;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  // Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  // Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  // Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  // Radar measurement noise standard deviation radius in m
  double std_radr_;

  // Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  // Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  // Weights of sigma points
  Eigen::VectorXd weights_;

  // State dimension
  int n_x_;

  // Augmented state dimension
  int n_aug_;

  // Lidar and radar state dimension measurement state
  int n_z_lidar_ = 2; // px and py
  int n_z_radar_ = 3; // rho, phi and rhodot

  // Sigma point spreading parameter
  double lambda_;

  // Log File
  std::shared_ptr<std::ofstream> log_file_{nullptr};
  std::shared_ptr<std::ofstream> nis_l_{nullptr};
  std::shared_ptr<std::ofstream> nis_r_{nullptr};
};

#endif  // UKF_H
