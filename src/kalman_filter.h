#ifndef KALMAN_FILTER_H_
#define KALMAN_FILTER_H_
#include "Eigen/Dense"

using namespace Eigen;

class KalmanFilter {
public:

  VectorXd x_; // state vector
  MatrixXd F_; // state transistion matrix
  MatrixXd P_; // state covariance matrix
  MatrixXd Q_; // process covariance matrix
  MatrixXd H_; // measurement matrix
  MatrixXd R_; // measurement covariance matrix

  /**
   * Constructor
   */
  KalmanFilter();

  /**
   * Destructor
   */
  virtual ~KalmanFilter();

  void Initialize(VectorXd &x, MatrixXd &F, MatrixXd &P, MatrixXd &Q);

  /**
   * Prediction Predicts the state and the state covariance
   * using the process model
   * @param delta_T Time between k and k+1 in s
   */
  void Predict();

  /**
   * Updates the state by using standard Kalman Filter equations
   * @param z The measurement at k+1
   */
  void Update(const VectorXd &y);

};

#endif /* KALMAN_FILTER_H_ */
