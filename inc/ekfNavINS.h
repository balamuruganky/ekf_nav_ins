/*
*) Refactor the code to remove reduentent part and improve the readabilty. 
*) Compiled for Linux with C++14 standard
Copyright (c) 2021 Balamurugan Kandan.
MIT License; See LICENSE.md for complete details
Author: 2021 Balamurugan Kandan
*/

/*
Updated to be a class, use Eigen, and compile as an Arduino library.
Added methods to get gyro and accel bias. Added initialization to
estimated angles rather than assuming IMU is level.

Copyright (c) 2016 - 2019 Regents of the University of Minnesota and Bolder Flight Systems Inc.
MIT License; See LICENSE.md for complete details
Author: Brian Taylor
*/

/*
Addapted from earlier version
Copyright 2011 Regents of the University of Minnesota. All rights reserved.
Original Author: Adhika Lie
*/

#pragma once

#include <stdint.h>
#include <math.h>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <tuple>
#include <mutex>
#include <shared_mutex>
// ++++++++++++++++++++++++++++++++++++++++++++++++
// error characteristics of navigation parameters
// ++++++++++++++++++++++++++++++++++++++++++++++++
// Correlation time or time constant
constexpr float TAU_A = 100.0f;
// Correlati1on time or time constant
constexpr float TAU_G = 50.0f;
// Initial set of covariance
constexpr float P_P_INIT = 10.0f;
constexpr float P_V_INIT = 1.0f;
constexpr float P_A_INIT = 0.34906f;
constexpr float P_HDG_INIT = 3.14159f;
constexpr float P_AB_INIT = 0.9810f;
constexpr float P_GB_INIT = 0.01745f;
// acceleration due to gravity
constexpr float G = 9.807f;
// major eccentricity squared
constexpr double ECC2 = 0.0066943799901;
// earth semi-major axis radius (m)
constexpr double EARTH_RADIUS = 6378137.0;

struct gpsPosData {
    uint64_t pos_time;
    double lat;
    double lon;
    double alt;
};
using gpsPosDataPtr = std::shared_ptr<gpsPosData>;

struct gpsVelData {
    uint64_t vel_time;
    double vN;
    double vE;
    double vD;
};
using gpsVelDataPtr = std::shared_ptr<gpsVelData>;

struct imuData {
    uint64_t imu_time;
    float gyroX;
    float gyroY;
    float gyroZ;
    float acclX;
    float acclY;
    float acclZ;
};
using imuDataPtr = std::shared_ptr<imuData>;

struct magData {
    uint64_t mag_time;
    float hX;
    float hY;
    float hZ;
};
using magDataPtr = std::shared_ptr<magData>;

struct ekfState {
    uint64_t timestamp;
    Eigen::Vector3d lla;       // in radians
    // Velocities
    Eigen::Vector3d velNED;
    // The imu data.
    Eigen::Matrix<float,3,1> linear;
    Eigen::Matrix<float,3,1> angular;
    // Eular angles
    Eigen::Matrix<float,4,1> quat;       // Quaternion
    Eigen::Vector3d accl_bias;  // The bias of the acceleration sensor.
    Eigen::Vector3d gyro_bias; // The bias of the gyroscope sensor.
    // Covariance.
    Eigen::Matrix<float, 15, 15> cov;
};

/*
struct ImuData {
    double timestamp;      // In second.
    Eigen::Vector3d acc;   // Acceleration in m/s^2
    Eigen::Vector3d gyro;  // Angular velocity in radian/s.
};
using ImuDataPtr = std::shared_ptr<ImuData>;

struct MagData {
    double timestamp;      // In second.
    Eigen::Vector3d mag;   // All same unit. For eg. mT
};
using MagDataPtr = std::shared_ptr<MagData>;

struct GpsPosData {
    double timestamp;     // In second.
    Eigen::Vector3d lla;  // Latitude in degree, longitude in degree, and altitude in meter.
    Eigen::Matrix3d cov;  // Covariance in m^2.
};
using GpsPosDataPtr = std::shared_ptr<GpsPosData>;

struct State {
    double timestamp;
    Eigen::Vector3d lla;       // WGS84 position.
    Eigen::Vector3d acc_bias;  // The bias of the acceleration sensor.
    Eigen::Vector3d gyro_bias; // The bias of the gyroscope sensor.
    // Covariance.
    Eigen::Matrix<double, 15, 15> cov;
    // The imu data.
    ImuDataPtr imu_data_ptr;
    // The mag data
    MagDataPtr mag_data_ptr;
};
*/

class ekfNavINS {
  public:
    // ekf_update
    void ekf_update(uint64_t time/*, unsigned long TOW*/,   /* Time, Time of the week from GPS */
                    double vn, double ve, double vd,    /* Velocity North, Velocity East, Velocity Down */
                    double lat, double lon, double alt, /* GPS latitude, GPS longitude, GPS/Barometer altitude */
                    float p, float q, float r,          /* Gyro P, Q and R  */
                    float ax, float ay, float az,       /* Accelarometer X, Y and Z */
                    float hx, float hy, float hz        /* Magnetometer HX, HY, HZ */ );
    // returns whether the INS has been initialized
    bool initialized()          { return initialized_; }
    // returns the pitch angle, rad
    float getPitch_rad()        { return theta; }
    // returns the roll angle, rad
    float getRoll_rad()         { return phi; }
    // returns the heading angle, rad
    float getHeadingConstrainAngle180_rad()      { return constrainAngle180(psi); }
    float getHeading_rad()      { return psi; }
    // returns the INS latitude, rad
    double getLatitude_rad()    { return lat_ins; }
    // returns the INS longitude, rad
    double getLongitude_rad()   { return lon_ins; }
    // returns the INS altitude, m
    double getAltitude_m()      { return alt_ins; }
    // returns the INS north velocity, m/s
    double getVelNorth_ms()     { return vn_ins; }
    // returns the INS east velocity, m/s
    double getVelEast_ms()      { return ve_ins; }
    // returns the INS down velocity, m/s
    double getVelDown_ms()      { return vd_ins; }
    // returns the INS ground track, rad
    float getGroundTrack_rad()  { return atan2f((float)ve_ins,(float)vn_ins); }
    // returns the gyro bias estimate in the x direction, rad/s
    float getGyroBiasX_rads()   { return gbx; }
    // returns the gyro bias estimate in the y direction, rad/s
    float getGyroBiasY_rads()   { return gby; }
    // returns the gyro bias estimate in the z direction, rad/s
    float getGyroBiasZ_rads()   { return gbz; }
    // returns the accel bias estimate in the x direction, m/s/s
    float getAccelBiasX_mss()   { return abx; }
    // returns the accel bias estimate in the y direction, m/s/s
    float getAccelBiasY_mss()   { return aby; }
    // returns the accel bias estimate in the z direction, m/s/s
    float getAccelBiasZ_mss()   { return abz; }
    // return pitch, roll and yaw
    std::tuple<float,float,float> getPitchRollYaw(float ax, float ay, float az, float hx, float hy, float hz);
    //
    // ROS
    //
    bool imuDataUpdateEKF(const imuDataPtr imu, ekfState* ekfOut);
    void magDataUpdateEKF(const magDataPtr mag);
    void gpsPosDataUpdateEKF(const gpsPosDataPtr pos);
    void gpsVelDataUpdateEKF(const gpsVelDataPtr vel);
    void setAcclNoise (float acclNoise) { SIG_W_A = acclNoise; }
    void setAcclBias  (float acclBias)  { SIG_A_D = acclBias;  }
    void setGyroNoise (float gyroNoise) { SIG_W_G = gyroNoise; }
    void setGyroBias  (float gyroBias)  { SIG_G_D = gyroBias;  }
    void setStdDevGpsPosNE(float stdDevGpsPosNE)    { SIG_GPS_P_NE = stdDevGpsPosNE; }
    void setStdDevGpsPosD(float stdDevGpsPosD)      { SIG_GPS_P_D  = stdDevGpsPosD;  }
    void setStdDevGpsVelNE(float stdDevGpsVelNE)    { SIG_GPS_V_NE = stdDevGpsVelNE; }
    void setStdDevGpsVelD(float stdDevGpsVelD)      { SIG_GPS_V_D  = stdDevGpsVelD;  }

  private:
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // error characteristics of navigation parameters
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Std dev of Accelerometer Wide Band Noise (m/s^2)
    float SIG_W_A = 0.05f;
    // Std dev of Accelerometer Markov Bias
    float SIG_A_D = 0.01f;
    // Std dev of gyro output noise (rad/s)
    float SIG_W_G = 0.00175f;
    // Std dev of correlated gyro bias
    float SIG_G_D = 0.00025;
    // GPS measurement noise std dev (m)
    float SIG_GPS_P_NE = 3.0f;
    float SIG_GPS_P_D = 6.0f;
    // GPS measurement noise std dev (m/s)
    float SIG_GPS_V_NE = 0.5f;
    float SIG_GPS_V_D = 1.0f;
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////// member variables /////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    gpsPosDataPtr pGpsPosDat;
    gpsVelDataPtr pGpsVelDat;
    imuDataPtr pImuDat;
    magDataPtr pMagDat;
    mutable std::shared_mutex shMutex;
    // initialized
    bool initialized_ = false;
    bool is_gps_pos_initialized = false;
    bool is_gps_vel_initialized = false;
    bool is_imu_initialized = false;
    bool is_mag_initialized = false;
    // timing
    uint64_t _tprev;
    //float _dt;
    unsigned long previousTOW;
    // estimated attitude
    float phi, theta, psi;
    // estimated NED velocity
    double vn_ins, ve_ins, vd_ins;
    // estimated location
    double lat_ins, lon_ins, alt_ins;
    // magnetic heading corrected for roll and pitch angle
    float Bxc, Byc;
    // accelerometer bias
    float abx = 0.0, aby = 0.0, abz = 0.0;
    // gyro bias
    float gbx = 0.0, gby = 0.0, gbz = 0.0;
    // earth radius at location
    double Re, Rn, denom;
    // State matrix
    Eigen::Matrix<float,15,15> Fs = Eigen::Matrix<float,15,15>::Identity();
    // State transition matrix
    Eigen::Matrix<float,15,15> PHI = Eigen::Matrix<float,15,15>::Zero();
    // Covariance matrix
    Eigen::Matrix<float,15,15> P = Eigen::Matrix<float,15,15>::Zero();
    // For process noise transformation
    Eigen::Matrix<float,15,12> Gs = Eigen::Matrix<float,15,12>::Zero();
    Eigen::Matrix<float,12,12> Rw = Eigen::Matrix<float,12,12>::Zero();
    // Process noise matrix
    Eigen::Matrix<float,15,15> Q = Eigen::Matrix<float,15,15>::Zero();
    // Gravity model
    Eigen::Matrix<float,3,1> grav = Eigen::Matrix<float,3,1>::Zero();
    // Rotation rate
    Eigen::Matrix<float,3,1> om_ib = Eigen::Matrix<float,3,1>::Zero();
    // Specific force
    Eigen::Matrix<float,3,1> f_b = Eigen::Matrix<float,3,1>::Zero();
    // DCM
    Eigen::Matrix<float,3,3> C_N2B = Eigen::Matrix<float,3,3>::Zero();
    // DCM transpose
    Eigen::Matrix<float,3,3> C_B2N = Eigen::Matrix<float,3,3>::Zero();
    // Temporary to get dxdt
    Eigen::Matrix<float,3,1> dx = Eigen::Matrix<float,3,1>::Zero();
    Eigen::Matrix<double,3,1> dxd = Eigen::Matrix<double,3,1>::Zero();
    // Estimated INS
    Eigen::Matrix<double,3,1> estmimated_ins = Eigen::Matrix<double,3,1>::Zero();
    // NED velocity INS
    Eigen::Matrix<double,3,1> V_ins = Eigen::Matrix<double,3,1>::Zero();
    // LLA INS
    Eigen::Matrix<double,3,1> lla_ins = Eigen::Matrix<double,3,1>::Zero();
    // NED velocity GPS
    Eigen::Matrix<double,3,1> V_gps = Eigen::Matrix<double,3,1>::Zero();
    // LLA GPS
    Eigen::Matrix<double,3,1> lla_gps = Eigen::Matrix<double,3,1>::Zero();
    // Position ECEF INS
    Eigen::Matrix<double,3,1> pos_ecef_ins = Eigen::Matrix<double,3,1>::Zero();
    // Position NED INS
    Eigen::Matrix<double,3,1> pos_ned_ins = Eigen::Matrix<double,3,1>::Zero();
    // Position ECEF GPS
    Eigen::Matrix<double,3,1> pos_ecef_gps = Eigen::Matrix<double,3,1>::Zero();
    // Position NED GPS
    Eigen::Matrix<double,3,1> pos_ned_gps = Eigen::Matrix<double,3,1>::Zero();
    // Quat
    Eigen::Matrix<float,4,1> quat = Eigen::Matrix<float,4,1>::Zero();
    // dquat
    Eigen::Matrix<float,4,1> dq = Eigen::Matrix<float,4,1>::Zero();
    // difference between GPS and INS
    Eigen::Matrix<float,6,1> y = Eigen::Matrix<float,6,1>::Zero();
    // GPS measurement noise
    Eigen::Matrix<float,6,6> R = Eigen::Matrix<float,6,6>::Zero();
    Eigen::Matrix<float,15,1> x = Eigen::Matrix<float,15,1>::Zero();
    // Kalman Gain
    Eigen::Matrix<float,15,6> K = Eigen::Matrix<float,15,6>::Zero();
    Eigen::Matrix<float,6,15> H = Eigen::Matrix<float,6,15>::Zero();
    // skew symmetric
    Eigen::Matrix<float,3,3> sk(Eigen::Matrix<float,3,1> w);

    //////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////// member functions /////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    // ekf_init
    void ekf_init(uint64_t time, 
                 double vn,double ve,double vd, 
                 double lat,double lon,double alt,
                 float p,float q,float r,
                 float ax,float ay,float az,
                 float hx,float hy, float hz);
    // lla rate
    Eigen::Matrix<double,3,1> llarate(Eigen::Matrix<double,3,1> V, Eigen::Matrix<double,3,1> lla);
    Eigen::Matrix<double,3,1> llarate(Eigen::Matrix<double,3,1> V, double lat, double alt);
    // lla to ecef
    Eigen::Matrix<double,3,1> lla2ecef(Eigen::Matrix<double,3,1> lla);
    // ecef to ned
    Eigen::Matrix<double,3,1> ecef2ned(Eigen::Matrix<double,3,1> ecef, Eigen::Matrix<double,3,1> pos_ref);
    // quaternion to dcm
    Eigen::Matrix<float,3,3> quat2dcm(Eigen::Matrix<float,4,1> q);
    // quaternion multiplication
    Eigen::Matrix<float,4,1> qmult(Eigen::Matrix<float,4,1> p, Eigen::Matrix<float,4,1> q);
    // maps angle to +/- 180
    float constrainAngle180(float dta);
    // maps angle to 0-360
    float constrainAngle360(float dta);
    // Returns Radius - East West and Radius - North South
    constexpr std::pair<double, double> earthradius(double lat);
    // Yaw, Pitch, Roll to Quarternion
    Eigen::Matrix<float,4,1> toQuaternion(float yaw, float pitch, float roll);
    // Quarternion to Yaw, Pitch, Roll
    std::tuple<float, float, float> toEulerAngles(Eigen::Matrix<float,4,1> quat);
    // Update Jacobian matrix
    void updateJacobianMatrix();
    // Update Process Noise and Covariance Time
    void updateProcessNoiseCovarianceTime(float _dt);
    // Update Gyro and Accelerometer Bias
    void updateBias(float ax,float ay,float az,float p,float q, float r);
    // Update 15 states after KF state update
    void update15statesAfterKF();
    // Update differece between predicted and calculated GPS and IMU values
    void updateCalculatedVsPredicted();
    void ekf_update(uint64_t time);
    void updateINS();
};
