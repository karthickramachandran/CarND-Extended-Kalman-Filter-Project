#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth)
{
    Eigen::VectorXd rmse(4);
    rmse << 0,0,0,0;
    if(estimations.size() != ground_truth.size() || estimations.size() == 0)
    {
            std::cerr << "Wrong size: Size is not equal to the ground truth or equal to zero" << endl;
            return rmse;
    }

    for(uint i=0; i < estimations.size(); ++i)
    {
        Eigen::VectorXd tmp = estimations[i] - ground_truth[i];

        tmp = tmp.array()*tmp.array();
        rmse += tmp;
    }

    rmse = rmse/estimations.size();

    rmse = rmse.array().sqrt();

    return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state)
{
    Eigen::MatrixXd H_(3,4);

    if (x_state.size()!=4)
    {
        std::cerr<<"Must have size of 4"<<std::endl;
        return H_;
    }
    float px=x_state(0);
    float py=x_state(1);
    float vx=x_state(2);
    float vy=x_state(3);

    float denom = px*px+py*py;

    float r_30 = (py*((vx*py)-(vy*px)))/pow(denom, 1.5);
    float r_31 = (px*((vy*px)-(vx*py)))/pow(denom, 1.5);

    if (denom< 0.0001)
    {
        std::cout << "Error: Jacobian Function -> Division by Zero" << std::endl;
        return H_;
    }

    std::cout<<"Calulcation jacobian now ...."<<std::endl;
    H_ <<   px/sqrt(denom), py/sqrt(denom), 0,      0,
            -py/denom,    px/denom,   0,      0,
            r_30,       r_31,   px/sqrt(denom),   py/sqrt(denom);

    return H_;

}
double Tools::NormalizeAngle(double angle) {
    /*
    This method normalizes angles to be between -pi and pi
    */
    double a = fmod(angle + M_PI, 2 * M_PI);
    return a >= 0 ? (a - M_PI) : (a + M_PI);
}
