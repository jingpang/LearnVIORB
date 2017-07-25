#ifndef IMUDATA_H
#define IMUDATA_H

#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <vector>

namespace ORB_SLAM2
{

using namespace Eigen;

class IMUData
{

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    typedef std::vector<IMUData, Eigen::aligned_allocator<IMUData> > vector_t;

    // covariance of measurement
    static Matrix3d _gyrMeasCov;
    static Matrix3d _accMeasCov;
    static Matrix3d getGyrMeasCov(void) {return _gyrMeasCov;}
    static Matrix3d getAccMeasCov(void) {return _accMeasCov;}

    // covariance of bias random walk
    static Matrix3d _gyrBiasRWCov;
    static Matrix3d _accBiasRWCov;
    static Matrix3d getGyrBiasRWCov(void) {return _gyrBiasRWCov;}
    static Matrix3d getAccBiasRWCov(void) {return _accBiasRWCov;}

    static double _gyrBiasRw2;
    static double _accBiasRw2;
    static double getGyrBiasRW2(void) {return _gyrBiasRw2;}
    static double getAccBiasRW2(void) {return _accBiasRw2;}


    IMUData(const double& gx, const double& gy, const double& gz,
            const double& ax, const double& ay, const double& az,
            const double& t);
    //IMUData(const IMUData& imu);

    // Raw data of imu's
    Vector3d _g;    //gyr data
    Vector3d _a;    //acc data
    double _t;      //timestamp
};

}

#endif // IMUDATA_H
