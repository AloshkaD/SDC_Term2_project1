#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <math.h>
using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  H_laser_ << 1, 0, 0, 0,
        0, 1, 0, 0;
  /*
  Hj_ << 1, 1, 0, 0,
        1, 1, 0, 0,
        1, 1, 0, 0,
        1, 1, 1, 1;
  */
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0,
          0, 1, 0, 1,
          0, 0, 1, 0,
          0, 0, 0, 1;
  ekf_.P_ = MatrixXd(4, 4);
  //acceleration noise, from the course it 
  //was advised to set it to 9
  int noise_ax = 9;
  int noise_ay = 9;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/

  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    //cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    //previous_timestamp_ = measurement_pack.timestamp_;
    // Zero initialization for the first measurment 
    ekf_.x_ << 1, 1, 1, 1;
    //ekf_.x_ << 0, 0, 0, 0;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float theta = measurement_pack.raw_measurements_[1];
      float px = measurement_pack.raw_measurements_[0]*cos(theta);
      float py = measurement_pack.raw_measurements_[0]*sin(theta);
      if(px == 0 or py == 0){
        return;
      }
      ekf_.x_ << px, py, 0, 0;
          //ekf_.x_(0) = measurement_pack.raw_measurements_(0)*cos(measurement_pack.raw_measurements_(1));
	  	    //ekf_.x_(1) = measurement_pack.raw_measurements_(0)*sin(measurement_pack.raw_measurements_(1));
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state with 0,0,0,0 value
      */
      //float px_laser = measurement_pack.raw_measurements_[0];
      //float py_laser = measurement_pack.raw_measurements_[1];

      //if(px_laser == 0 or py_laser == 0){
        //return;
            ekf_.x_(0) = measurement_pack.raw_measurements_(0);
	          ekf_.x_(1) = measurement_pack.raw_measurements_(1);
      }
      //ekf_.x_ << px_laser, py_laser, 0, 0;
    

    
    ekf_.F_ = MatrixXd(4, 4);
    ekf_.F_ << 1, 0, 0, 0,
               0, 1, 0, 0,
               0, 0, 1, 0,
               0, 0, 0, 1;
    // done initializing, no need to predict or update
    previous_timestamp_ = measurement_pack.timestamp_;
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;


  //if(dt > 0.001){


    float dt_2 = dt * dt;
    float dt_3 = dt_2 * dt;
    float dt_4 = dt_3 * dt;

    ekf_.F_(0, 2) = dt;
    ekf_.F_(1, 3) = dt;

    ekf_.Q_ = MatrixXd(4, 4);
    ekf_.Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
        0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
        dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
        0, dt_3/2*noise_ay, 0, dt_2*noise_ay;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */
  //}
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
        //MatrixXd z_predicted (3,1);
        //Hj_ << tools.CalculateJacobian(ekf_.x_);
        //ekf_.H_ = Hj_;
        //ekf_.R_ = R_radar_;

        //Convert from cartesian to polar notation prior calling the updating function
        float x = ekf_.x_(0);

        float y = ekf_.x_(1);
        float vx = ekf_.x_(2);
        float vy = ekf_.x_(3);
        
        ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
        ekf_.R_ = R_radar_;
        //h(x')
        
        float rho = sqrt(x*x + y*y);
        float theta = atan2(y,x);
        float radial_velocity = (x*vx + y*vy)/rho;
        VectorXd z_pred = VectorXd(3);
        z_pred << rho, theta, radial_velocity;


        ekf_.UpdateRadar(measurement_pack.raw_measurements_, z_pred);
  } else {
    // Laser updates
    	  ekf_.H_ = H_laser_;
	      ekf_.R_ = R_laser_;

	      ekf_.UpdateLaser(measurement_pack.raw_measurements_);
  }

  // print the output
  //cout << "x_ = " << ekf_.x_ << endl;
  //cout << "P_ = " << ekf_.P_ << endl;
}

