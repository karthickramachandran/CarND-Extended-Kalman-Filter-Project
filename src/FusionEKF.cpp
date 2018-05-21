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
                0,      0.0225;

    //measurement covariance matrix - radar
    R_radar_ << 0.09,   0,      0,
                0,      0.0009, 0,
                0,      0,      0.09;

    noise_ax = 9;
    noise_ay = 9;

    ekf_.P_ = MatrixXd(4, 4);
    ekf_.P_ << 1, 0, 0, 0,
               0, 1, 0, 0,
               0, 0, 1000, 0,
               0, 0, 0, 1000;

    H_laser_ << 1, 0, 0, 0,
                0, 1, 0, 0;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack)
{
   /*****************************************************************************
   *  Initialization
   ****************************************************************************/
    if (!is_initialized_) {
        // first measurement
        ekf_.x_ = VectorXd(4);

        ekf_.x_ << 1, 1, 1, 1;

        ekf_.F_ = MatrixXd(4, 4);

        ekf_.F_ <<  1,0,1,0,
                0,1,0,0,
                0,0,1,0,
                0,0,0,1;

        ekf_.P_ = MatrixXd(4, 4);

        ekf_.P_ <<  1,0,0,0,
                0,1,0,0,
                0,0,1000,0,
                0,0,0,1000;



        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
        {
            //std::cout<<"RADAR - measurement_pack.raw_measurements_ : "<<measurement_pack.raw_measurements_<<std::endl;

            float range         = measurement_pack.raw_measurements_[0];
            float theta         = measurement_pack.raw_measurements_[1];
            float range_data    = measurement_pack.raw_measurements_[2];

            float px= range*cos(theta);
            float py= range*sin(theta);
            float vx=range_data*cos(theta);
            float vy=range_data*sin(theta);

            if (px<0.0001)
                px=0.0001;
            if (py<0.0001)
                py=0.0001;

            ekf_.x_ << px, py, vx, vy;

        }
        else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER)
        {
        //    std::cout<<"LIDAR - measurement_pack.raw_measurements_ : "<<measurement_pack.raw_measurements_<<std::endl;

            float px= measurement_pack.raw_measurements_[0];
            float py= measurement_pack.raw_measurements_[1];
            float vx=0;
            float vy=0;

            ekf_.x_ << px, py, vx, vy;
            std::cout<<"Updated the Laser Init"<<std::endl;
        }

        is_initialized_ = true;
        previous_timestamp_ =measurement_pack.timestamp_;
        std::cout<<"previous_timestamp_ is initialized"<<std::endl;
        return;
    }

   /*****************************************************************************
   *  Prediction
   ****************************************************************************/

    double dt_1 = (measurement_pack.timestamp_ -previous_timestamp_ )/1000000.0;
    double dt_2=dt_1*dt_1;
    double dt_3=dt_2*dt_1;
    double dt_4=dt_2*dt_2;

    ekf_.F_  =Eigen::MatrixXd(4,4);

    ekf_.F_<<   1, 0, dt_1, 0,
                0, 1, 0,    dt_1,
                0, 0, 1,    0,
                0, 0, 0,    1;

    ekf_.Q_  =Eigen::MatrixXd(4,4);  

    ekf_.Q_ <<      0.25*dt_4*noise_ax,     0,                      0.5*dt_3*noise_ax,  0,
                    0,                      0.25*dt_4*noise_ay,     0,                  0.5*dt_3*noise_ay,
                    0.5*dt_3*noise_ax,      0,                      dt_2*noise_ax,      0,
                    0,                      0.5*dt_3*noise_ay,      0,                  dt_2*noise_ay          ;

    ekf_.Predict();

    std::cout<<"Prediction step is scuccessfully completed"<<std::endl;
    /*****************************************************************************
   *  Update
   ****************************************************************************/

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
    {
        // Radar updates
        ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
        ekf_.R_ = R_radar_;;
        ekf_.UpdateEKF(measurement_pack.raw_measurements_);
    }
    else
    {
        // Laser updates
        ekf_.H_ = H_laser_;
        ekf_.R_ = R_laser_;;
        ekf_.Update(measurement_pack.raw_measurements_);
    }

    // print the output
    cout << "x_ = " << ekf_.x_ << endl;
    cout << "P_ = " << ekf_.P_ << endl;

    previous_timestamp_ =measurement_pack.timestamp_;
}
