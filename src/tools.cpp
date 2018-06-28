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
  rmse << 0,0,0,0;
  
  if ((estimations.size() == 0) ||
      (ground_truth.size() == 0) ||
      (estimations.size() != ground_truth.size()))
  {
    std::cout << "CalculateRMSE - Error - Input vector dimensions are not as expected";
    std::cout << "(e:"<< estimations.size() << ", g:" << ground_truth.size() << ")" << std::endl;
	    
    return rmse;
  }

  VectorXd res;
  for(int i=0; i < estimations.size(); ++i){
    res = estimations[i] - ground_truth[i];
    res = res.array() * res.array();
    rmse += res;
		
   }

   rmse = rmse/estimations.size();
   rmse = rmse.array().sqrt();

   return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {

  MatrixXd Hj(3,4);
  
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  //check division by zero
  if ((px == 0) && (py == 0))
  {
     Hj << 0, 0, 0, 0,
	   0, 0, 0, 0,
	   0, 0, 0, 0;
     std::cout << "CalculateJacobian() - Error - Division by Zero" << std::endl;
	    
     return Hj;
  }

  //compute the Jacobian matrix
  float pxpysq = px*px + py*py;
  float sqrtpxpysq = sqrt(pxpysq);
  float sqrt3pxpysq = pxpysq * sqrt(pxpysq);
	
  Hj <<  px / sqrtpxpysq, py / sqrtpxpysq, 0, 0,
	 -py / pxpysq,     px / pxpysq,     0, 0,
	 py * (vx*py - vy*px) / sqrt3pxpysq, px * (vy*px - vx*py) / sqrt3pxpysq, px / sqrtpxpysq, py / sqrtpxpysq;

  return Hj;
}
