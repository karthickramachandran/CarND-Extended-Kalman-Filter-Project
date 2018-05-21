#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in)
{
    x_ = x_in;
    P_ = P_in;
    F_ = F_in;
    H_ = H_in;
    R_ = R_in;
    Q_ = Q_in;
    std::cout<<"Kalman Init is performed "<<std::endl;
}

void KalmanFilter::Predict()
{
    /// p'=F*P*Ft+Q
    x_          = F_ * x_;
    MatrixXd Ft = F_.transpose();
    P_          = F_ * P_ * Ft + Q_;

}

void KalmanFilter::Update(const VectorXd &z)
{
    /// P = (I - K*H)*P'
    /// Find I, K
    /// Laser
    VectorXd z_pred = H_ * x_;
    VectorXd y      = z - z_pred; /// y = z - H*x', where x'=F*x+u
    MatrixXd Ht     = H_.transpose();
    MatrixXd S      = H_ * P_ * Ht + R_;
    MatrixXd Si     = S.inverse();
    MatrixXd PHt    = P_ * Ht;
    MatrixXd K      = PHt * Si;/// K=P'*Ht*Si

    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z)
{
    /// P = (I - K*H)*P'
    /// Find I, K


    /// Radar
    Tools tools;

    /// Convert Cartesian coordinates to Polar
    float px, py, vx, vy;
    px = x_(0);
    py = x_(1);
    vx = x_(2);
    vy = x_(3);

    float pos=sqrt(px*px+py*py);
    float theta, pos_rate;

    if (fabs(pos) <0.0001)
    {
        theta = 0;
        pos_rate= 0;
    }
    else
    {
        theta = std::atan2(py, px);
        pos_rate = (px*vx+py*vy)/pos;
    }

    VectorXd z_pred = VectorXd(3);
    z_pred << pos, theta, pos_rate;


    VectorXd y      = z - z_pred; /// y= z-H*x'

    /// Normalize
    y(1)=tools.NormalizeAngle(y(1));

    MatrixXd Ht     = H_.transpose();
    MatrixXd S      = H_ * P_ * Ht + R_;
    MatrixXd Si     = S.inverse();
    MatrixXd PHt    = P_ * Ht;
    MatrixXd K      = PHt * Si;/// K= P'*Ht*Si

    //new estimate
    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_; /// P =(I-K*H)*P'
}
