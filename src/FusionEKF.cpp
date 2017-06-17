#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <math.h>
#include <iostream>

using namespace std;
using namespace Eigen;

// Constructor
AjitFusion::AjitFusion() {
  is_initialized_ = false;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  R_radar_ = MatrixXd(3, 3);
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  H_laser_ = MatrixXd(2, 4);
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  // process noise characteristics
  S_ax = 9;
  S_ay = 9;

}

// Destructor
AjitFusion::~AjitFusion() {}

VectorXd AjitFusion::NormalizeAngle(const Eigen::VectorXd& y) {
	VectorXd y_norm(3);

	double angle = y[1];
	
	while (angle > M_PI){
		angle -= 2*M_PI;
	}
	
	while (angle < -M_PI) {
		angle += 2*M_PI;
	}
	
	y_norm << y[0], angle, y[2];
	return y_norm;
}

void AjitFusion::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {

    double rho, phi, rho_dot, px, py, vx, vy;
    VectorXd x(4);
    MatrixXd F(4, 4);
    MatrixXd P(4, 4);
    MatrixXd Q(4, 4);

    // Decide how to initialize depending on what type of measurement we have
    if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      x << measurement_pack.raw_measurements_, 0, 0; // Assuming starting at 0 velocity
      P << 0.5, 0, 0, 0,
           0, 0.5, 0, 0,
           0, 0, 100, 0,
           0, 0, 0, 100;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      P << 10, 0, 0, 0,
           0, 10, 0, 0,
           0, 0, 0.5, 0,
           0, 0, 0, 0.5;
      rho = measurement_pack.raw_measurements_(0);
      phi = measurement_pack.raw_measurements_(1);
      rho_dot = measurement_pack.raw_measurements_(2);
      px = rho * cos(phi);
      py = rho * sin(phi);
      vx = 0;
      vy = 0;
      x << px, py, vx, vy;
    }

    ekf_.Initialize(x, F, P, Q); // Perform initiliaziation of kf object
    previous_timestamp_ = measurement_pack.timestamp_;
    is_initialized_ = true;

    return;

  }

  /*****************************************************************************
  *  Prediction
  ****************************************************************************/
  // Calculate elapsed time
  dt_ = (measurement_pack.timestamp_ - previous_timestamp_)/1000000.0; // Assumed seconds

  previous_timestamp_ = measurement_pack.timestamp_;

  // Update the F, Q matrices given elapsed time
  ekf_.F_ << 1, 0, dt_, 0,
            0, 1, 0, dt_,
            0, 0, 1, 0,
            0, 0, 0, 1;

  double dt2 = dt_ * dt_;
  double dt3 = dt2 * dt_;
  double dt4 = dt3 * dt_;
  double dt4_4 = dt4/4;
  double dt3_2 = dt3/2; 


  ekf_.Q_ << dt4_4 * S_ax, 0, dt3_2 * S_ax, 0,
            0, dt4_4 * S_ay, 0, dt3_2 * S_ay,
            dt3_2 * S_ax, 0, dt2 * S_ax, 0,
            0, dt3_2 * S_ay, 0, dt2 * S_ay;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/
  VectorXd y;
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar H, R matrices
    double rho_pred, phi_pred, rhodot_pred;

    // sqrt(px^2 + py^2)
    rho_pred = sqrt(pow(ekf_.x_[0], 2) + pow(ekf_.x_[1], 2));

    // arctan(py/px)
    phi_pred = 0.0;
    if (fabs(ekf_.x_[0]) > 0.001) {
      phi_pred  = atan2(ekf_.x_[1], ekf_.x_[0]);
    }

    // (px*vx + py*vy)/rho_pred
    rhodot_pred = 0.0;
    if (fabs(rho_pred) > 0.001) {
      rhodot_pred = (ekf_.x_[0] * ekf_.x_[2] + ekf_.x_[1] * ekf_.x_[3]) / rho_pred;
    }

    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);;
    ekf_.R_ = R_radar_;
    VectorXd z_pred(3);
    z_pred << rho_pred, phi_pred, rhodot_pred;
    y = measurement_pack.raw_measurements_ - z_pred;
    y = NormalizeAngle(y);   // advised by reviewer
    ekf_.Update(y);

  }
  else {

    // Laser H, R matrices
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    y = measurement_pack.raw_measurements_ - ekf_.H_ * ekf_.x_;
    ekf_.Update(y);

  }

}
