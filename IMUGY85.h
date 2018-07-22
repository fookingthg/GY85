/*
	IMUGY85.h - IMU library for GY85 9DOF IMU (ADXL345 accelerometer, ITG3200 gyroscope and HMC5883L magnetometer) 
	using open source Madgwick quaternion filter
	Requires I2Cdev, ADXL345, HMC5883L, ITG3200 libraries (https://github.com/jrowberg/i2cdevlib/tree/master/Arduino)
	By thg(https://github.com/fookingthg). Jul 2018
	Based on multiple libraries of Jeff Rowberg(https://github.com/jrowberg) and others
*/
#ifndef IMUGY85_h
#define IMUGY85_h

#include <math.h>
#include <Wire.h>
#include "I2Cdev.h"
#include "ADXL345.h"
#include "HMC5883L.h"
#include "ITG3200.h"
#include <string.h>

class IMUGY85 {
	enum Ascale {
	  AFS_2G = 0,
	  AFS_4G,
	  AFS_8G,
	  AFS_16G
	};

	enum Gscale {
	  GFS_250DPS = 0,
	  GFS_500DPS,
	  GFS_1000DPS,
	  GFS_2000DPS,
	  GFS_CUSTOM
	};

	enum Mscale {
	  MFS_14BITS = 0, // 0.6 mG per LSB
	  MFS_16BITS      // 0.15 mG per LSB
	};

	uint8_t Gscale = GFS_CUSTOM;
	uint8_t Ascale = AFS_16G;
	uint8_t Mscale = MFS_16BITS; // Choose either 14-bit or 16-bit magnetometer resolution
	float aRes, gRes, mRes;      // scale resolutions per LSB for the sensors

	int16_t accelCount[3];  // Stores the 16-bit signed accelerometer sensor output
	int16_t gyroCount[3];   // Stores the 16-bit signed gyro sensor output
	int16_t magCount[3];    // Stores the 16-bit signed magnetometer sensor output
	float magCalibration[3] = {1, 1, 1}, magbias[3] = {0, 0, 0}; // Factory mag calibration and mag bias
	float magbias_max[3] = {0, 0, 0}; 
	float magbias_min[3] = {0, 0, 0};
	float gyroBias[3] = {0, 0, 0}, accelBias[3] = {0, 0, 0};      // Bias corrections for gyro and accelerometer


	float GyroMeasError = PI * (10.0f / 180.0f);   // gyroscope measurement error in rads/s (start at 40 deg/s)
	float GyroMeasDrift = PI * (0.0f  / 180.0f);   // gyroscope measurement drift in rad/s/s (start at 0.0 deg/s/s)
	float beta = sqrt(3.0f / 4.0f) * GyroMeasError;   // compute beta
	float zeta = sqrt(3.0f / 4.0f) * GyroMeasDrift;   // compute zeta, the other free parameter in the Madgwick scheme usually set to a small or zero value
	#define Kp 2.0f * 5.0f // these are the free parameters in the Mahony filter and fusion scheme, Kp for proportional feedback, Ki for integral
	#define Ki 0.1f


	double deltat = 0.0f, sum = 0.0f;        // integration interval for both filter schemes
	uint32_t lastUpdate = 0, firstUpdate = 0; // used to calculate integration interval
	uint32_t Now = 0;        // used to calculate integration interval

	float ax, ay, az, gx, gy, gz, mx, my, mz; // variables to hold latest sensor data values 
	float q[4] = {1.0f, 0.0f, 0.0f, 0.0f};    // vector to hold quaternion
	float eInt[3] = {0.0f, 0.0f, 0.0f};       // vector to hold integral error for Mahony method

	ADXL345 accel;
	ITG3200 gyro;
	HMC5883L mag;
	double pitch=0, roll=0, yaw=0;

    public:
    	IMUGY85();
    	void init();
    	void update();
    	double getRoll();
		double getPitch();
		double getYaw();
		double getRawYaw();
		double getAcceleration(double *a1, double *a2, double *a3);
		double getGyro(double *m1, double *m2, double *m3);
	private:
		void computeEuler();
		void getAres();
		void getGres();
		void getMres();
		void MadgwickQuaternionUpdate(float ax, float ay, float az, float gx, float gy, float gz, float mx, float my, float mz);
};

#endif