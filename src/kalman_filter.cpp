#include "kalman_filter.h"
#include <iostream>

using namespace std;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Initialize(VectorXd &x, MatrixXd &F, MatrixXd &P,
                              MatrixXd &Q) {
  x_ = x;
  F_ = F;
  P_ = P;
  Q_ = Q;
}

void KalmanFilter::Predict() {
  // State estimate prediction
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &y) {
  // State estimate update using Kalman Filter equations
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd K = P_ * Ht * S.inverse();

  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  x_ = x_ + K * y;
  P_ = (I - K * H_) * P_;
}
