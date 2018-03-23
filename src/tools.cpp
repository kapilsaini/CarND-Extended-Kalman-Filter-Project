#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  
  VectorXd rmse(4);
	rmse << 0.0, 0.0, 0.0, 0.0;
  for (int i = 0; i < estimations.size() ; i++) {
    VectorXd diff = ground_truth[i] - estimations[i];
    diff = diff.array() * diff.array();
    rmse += diff;
  }
  
  rmse = rmse / estimations.size();
  rmse = rmse.array().sqrt();
  
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  MatrixXd Hj = MatrixXd::Zero(3, 4);
  float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);
  
  float pxy = px*px + py*py;
  float pxy_sqrt = sqrt(pxy);
  
	float THRESH = 0.0001;
	
	float Hj01 = px/pxy_sqrt;
	float Hj02 = py/pxy_sqrt;
	Hj << Hj01, Hj02, 0.0, 0.0,
				-py/pxy, px/pxy, 0.0, 0.0,
				 py*(vx*py - vy*px)/(pxy*pxy_sqrt), px*(vy*px - vx*py)/(pxy*pxy_sqrt), Hj01, Hj02;
  return Hj;
}