#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using namespace std;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
    
  */
  	if (estimations.size() == 0 || ground_truth.size()!=estimations.size()) 
	{
	   
	   cout<<"Invalid ground trusth or estimation data input for RMSE"<<"\n";

     return rmse;
	    
	}
  	for(unsigned int i=0; i < estimations.size(); ++i){

		VectorXd residual = estimations[i] - ground_truth[i];

		//coefficient-wise multiplication
		residual = residual.array()*residual.array();
		rmse += residual;
	}

	//calculate the mean
	rmse = rmse/estimations.size();

	//calculate the squared root
	rmse = rmse.array().sqrt();

	//return the result
	return rmse;


}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
  MatrixXd Hj(3,4);
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);
  
  if (((px*px)+(py*py)) == 0)
    {
        cout<<"divide by zero"<<"\n";
        
    }
    
  else{
    
	//check division by zero
      float var1 = px/(sqrt((px*px)+(py*py)));
      float var2 = -(py/((px*px)+(py*py)));
      float var3 = py*(vx*py - vy*px)/pow((px*px)+(py*py),3/2);
      float var4 = py/(sqrt((px*px)+(py*py)));
      float var5 = px/((px*px)+(py*py));
      float var6 = px*(vy*px - vx*py)/pow((px*px)+(py*py),3/2);
      float var7 = 0;
      float var8 = 0;
      float var9 = px/(sqrt((px*px)+(py*py)));
      float var10 = 0;
      float var11 = 0;
      float var12 = py/(sqrt((px*px)+(py*py)));
      
      //compute the Jacobian matrix
      
      Hj << var1,var4,var7,var10,
            var2,var5,var8,var11,
            var3,var6,var9,var12;
    }
	return Hj;
}
