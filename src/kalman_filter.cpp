#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
VectorXd convert_cartesian_to_polar(const VectorXd& v);
// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  x_ = F_ * x_;
	P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  VectorXd z_pred = H_ * x_;
	VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();  
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd K = P_ * Ht * S.inverse(); 
  MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
  x_ = x_ + (K * y);
	P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  VectorXd z_pred = convert_cartesian_to_polar(x_);
	VectorXd y = z - z_pred;
	y(1) = atan2(sin(y(1)), cos(y(1)));
  MatrixXd Ht = H_.transpose();  
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd K = P_ * Ht * S.inverse(); 
  MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
  x_ = x_ + (K * y);
	P_ = (I - K * H_) * P_;
}

VectorXd convert_cartesian_to_polar(const VectorXd& v){

  float THRESH = 0.0001;
  VectorXd polar_vector(3);

  float px = v(0);
  float py = v(1);
  float vx = v(2);
  float vy = v(3);

  float rho = sqrt( px * px + py * py);

	float phi = 0.0;
	phi = atan2(py, px);
	float drho = ( px * vx + py * vy )/rho;
  polar_vector << rho, phi, drho;
  return polar_vector;
}