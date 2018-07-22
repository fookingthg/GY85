/*
	IMUGY85.cpp - IMU library for GY85 9DOF IMU (ADXL345 accelerometer, ITG3200 gyroscope and HMC5883L magnetometer) 
	using open source Madgwick quaternion filter
	Requires I2Cdev, ADXL345, HMC5883L, ITG3200 libraries (https://github.com/jrowberg/i2cdevlib/tree/master/Arduino)
	By thg(https://github.com/fookingthg). Jul 2018
	Based on multiple libraries of Jeff Rowberg(https://github.com/jrowberg)
*/
#include "IMUGY85.h"

IMUGY85::IMUGY85()
{

}

void IMUGY85::init()
{
	Wire.begin();
	accel.initialize();
	gyro.initialize();
	mag.initialize();
}

void IMUGY85::update()
{  

	accel.getAcceleration(&accelCount[0], &accelCount[1], &accelCount[2]);
	gyro.getRotation(&gyroCount[0], &gyroCount[1], &gyroCount[2]);
	mag.getHeading(&magCount[0], &magCount[1], &magCount[2]);

	getAres();

	ax = (float)accelCount[0]*aRes; // - accelBias[0];  // get actual g value, this depends on scale being set
	ay = (float)accelCount[1]*aRes; // - accelBias[1];   
	az = (float)accelCount[2]*aRes; // - accelBias[2];  

	getGres();

	gx = (float)gyroCount[0]*gRes;  // get actual gyro value, this depends on scale being set
	gy = (float)gyroCount[1]*gRes;  
	gz = (float)gyroCount[2]*gRes;   

	getMres();  

	// check the magnetic max. and min. Values
	if (magCount[0] > magbias_max[0]) magbias_max[0] = magCount[0];
	if (magCount[1] > magbias_max[1]) magbias_max[1] = magCount[1];
	if (magCount[2] > magbias_max[2]) magbias_max[2] = magCount[2];

	if (magCount[0] < magbias_min[0]) magbias_min[0] = magCount[0];
	if (magCount[1] < magbias_min[1]) magbias_min[1] = magCount[1];
	if (magCount[2] < magbias_min[2]) magbias_min[2] = magCount[2];

	magbias[0] = ((magbias_max[0] + magbias_min[0])/2);  // User environmental x-axis correction in milliGauss, should be automatically calculated
	magbias[1] = ((magbias_max[1] + magbias_min[1])/2);  // User environmental x-axis correction in milliGauss
	magbias[2] = ((magbias_max[2] + magbias_min[2])/2);  // User environmental x-axis correction in milliGauss 

	mx = (float)(magCount[0]) ;  // get actual magnetometer value, this depends on scale being set
	my = (float)(magCount[1]) ;  
	mz = (float)(magCount[2]);   

	//    mx = (float)(magCount[0]- magbias[0])*mRes*magCalibration[0] ;  // get actual magnetometer value, this depends on scale being set
	//    my = (float)(magCount[1]- magbias[1])*mRes*magCalibration[1] ;  
	//    mz = (float)(magCount[2]- magbias[2])*mRes*magCalibration[2] ;   

	Now = micros();
	deltat = ((Now - lastUpdate)/1000000.0f); // set integration time by time elapsed since last filter update
	lastUpdate = Now;
	MadgwickQuaternionUpdate(ax, ay, az, (gx)*PI/180.0f, (gy)*PI/180.0f, (gz)*PI/180.0f,  mx,  my, mz); // change mx & my, and negate az to have all axes in one direction!!!
	computeEuler();
}

double IMUGY85::getRoll()
{
    return degrees(roll);
}

double IMUGY85::getPitch()
{
    return degrees(pitch);
}

double IMUGY85::getYaw()
{
    double indegrees = degrees(yaw);
    return indegrees<0?indegrees+360:indegrees;
}

double IMUGY85::getRawYaw()
{
    return degrees(yaw);
}

double IMUGY85::getAcceleration(double *a1, double *a2, double *a3)
{
    *a1 = ax;
    *a2 = ay;
    *a3 = az;
}

double IMUGY85::getGyro(double *m1, double *m2, double *m3)
{
    *m1 = gx;
    *m2 = gy;
    *m3 = gz;
}

void IMUGY85::computeEuler()
{
    double q2sqr = q[2] * q[2];
    double t0 = -2.0 * (q2sqr + q[3] * q[3]) + 1.0;
    double t1 = +2.0 * (q[1] * q[2] + q[0] * q[3]);
    double t2 = -2.0 * (q[1] * q[3] - q[0] * q[2]);
    double t3 = +2.0 * (q[2] * q[3] + q[0] * q[1]);
    double t4 = -2.0 * (q[1] * q[1] + q2sqr) + 1.0;

    t2 = t2 > 1.0 ? 1.0 : t2;
    t2 = t2 < -1.0 ? -1.0 : t2;

    pitch = asin(t2);
    roll = atan2(t3, t4);
    yaw = atan2(t1, t0);
}

// Implementation of Sebastian Madgwick's "...efficient orientation filter for... inertial/magnetic sensor arrays"
// (see http://www.x-io.co.uk/category/open-source/ for examples and more details)
// which fuses acceleration, rotation rate, and magnetic moments to produce a quaternion-based estimate of absolute
// device orientation -- which can be converted to yaw, pitch, and roll. Useful for stabilizing quadcopters, etc.
// The performance of the orientation filter is at least as good as conventional Kalman-based filtering algorithms
// but is much less computationally intensive---it can be performed on a 3.3 V Pro Mini operating at 8 MHz!
void IMUGY85::MadgwickQuaternionUpdate(float ax, float ay, float az, float gx, float gy, float gz, float mx, float my, float mz)
{

	float q1 = q[0], q2 = q[1], q3 = q[2], q4 = q[3];   // short name local variable for readability
	float norm;
	float hx, hy, _2bx, _2bz;
	float s1, s2, s3, s4;
	float qDot1, qDot2, qDot3, qDot4;

	// Auxiliary variables to avoid repeated arithmetic
	float _2q1mx;
	float _2q1my;
	float _2q1mz;
	float _2q2mx;
	float _4bx;
	float _4bz;
	float _2q1 = 2.0f * q1;
	float _2q2 = 2.0f * q2;
	float _2q3 = 2.0f * q3;
	float _2q4 = 2.0f * q4;
	float _2q1q3 = 2.0f * q1 * q3;
	float _2q3q4 = 2.0f * q3 * q4;
	float q1q1 = q1 * q1;
	float q1q2 = q1 * q2;
	float q1q3 = q1 * q3;
	float q1q4 = q1 * q4;
	float q2q2 = q2 * q2;
	float q2q3 = q2 * q3;
	float q2q4 = q2 * q4;
	float q3q3 = q3 * q3;
	float q3q4 = q3 * q4;
	float q4q4 = q4 * q4;

	// Normalise accelerometer measurement
	norm = sqrt(ax * ax + ay * ay + az * az);
	if (norm == 0.0f) return; // handle NaN
	norm = 1.0f/norm;
	ax *= norm;
	ay *= norm;
	az *= norm;

	// Normalise magnetometer measurement
	norm = sqrt(mx * mx + my * my + mz * mz);
	if (norm == 0.0f) return; // handle NaN
	norm = 1.0f/norm;
	mx *= norm;
	my *= norm;
	mz *= norm;

	// Reference direction of Earth's magnetic field
	_2q1mx = 2.0f * q1 * mx;
	_2q1my = 2.0f * q1 * my;
	_2q1mz = 2.0f * q1 * mz;
	_2q2mx = 2.0f * q2 * mx;
	hx = mx * q1q1 - _2q1my * q4 + _2q1mz * q3 + mx * q2q2 + _2q2 * my * q3 + _2q2 * mz * q4 - mx * q3q3 - mx * q4q4;
	hy = _2q1mx * q4 + my * q1q1 - _2q1mz * q2 + _2q2mx * q3 - my * q2q2 + my * q3q3 + _2q3 * mz * q4 - my * q4q4;
	_2bx = sqrt(hx * hx + hy * hy);
	_2bz = -_2q1mx * q3 + _2q1my * q2 + mz * q1q1 + _2q2mx * q4 - mz * q2q2 + _2q3 * my * q4 - mz * q3q3 + mz * q4q4;
	_4bx = 2.0f * _2bx;
	_4bz = 2.0f * _2bz;

	// Gradient decent algorithm corrective step
	s1 = -_2q3 * (2.0f * q2q4 - _2q1q3 - ax) + _2q2 * (2.0f * q1q2 + _2q3q4 - ay) - _2bz * q3 * (_2bx * (0.5f - q3q3 - q4q4) + _2bz * (q2q4 - q1q3) - mx) + (-_2bx * q4 + _2bz * q2) * (_2bx * (q2q3 - q1q4) + _2bz * (q1q2 + q3q4) - my) + _2bx * q3 * (_2bx * (q1q3 + q2q4) + _2bz * (0.5f - q2q2 - q3q3) - mz);
	s2 = _2q4 * (2.0f * q2q4 - _2q1q3 - ax) + _2q1 * (2.0f * q1q2 + _2q3q4 - ay) - 4.0f * q2 * (1.0f - 2.0f * q2q2 - 2.0f * q3q3 - az) + _2bz * q4 * (_2bx * (0.5f - q3q3 - q4q4) + _2bz * (q2q4 - q1q3) - mx) + (_2bx * q3 + _2bz * q1) * (_2bx * (q2q3 - q1q4) + _2bz * (q1q2 + q3q4) - my) + (_2bx * q4 - _4bz * q2) * (_2bx * (q1q3 + q2q4) + _2bz * (0.5f - q2q2 - q3q3) - mz);
	s3 = -_2q1 * (2.0f * q2q4 - _2q1q3 - ax) + _2q4 * (2.0f * q1q2 + _2q3q4 - ay) - 4.0f * q3 * (1.0f - 2.0f * q2q2 - 2.0f * q3q3 - az) + (-_4bx * q3 - _2bz * q1) * (_2bx * (0.5f - q3q3 - q4q4) + _2bz * (q2q4 - q1q3) - mx) + (_2bx * q2 + _2bz * q4) * (_2bx * (q2q3 - q1q4) + _2bz * (q1q2 + q3q4) - my) + (_2bx * q1 - _4bz * q3) * (_2bx * (q1q3 + q2q4) + _2bz * (0.5f - q2q2 - q3q3) - mz);
	s4 = _2q2 * (2.0f * q2q4 - _2q1q3 - ax) + _2q3 * (2.0f * q1q2 + _2q3q4 - ay) + (-_4bx * q4 + _2bz * q2) * (_2bx * (0.5f - q3q3 - q4q4) + _2bz * (q2q4 - q1q3) - mx) + (-_2bx * q1 + _2bz * q3) * (_2bx * (q2q3 - q1q4) + _2bz * (q1q2 + q3q4) - my) + _2bx * q2 * (_2bx * (q1q3 + q2q4) + _2bz * (0.5f - q2q2 - q3q3) - mz);
	norm = sqrt(s1 * s1 + s2 * s2 + s3 * s3 + s4 * s4);    // normalise step magnitude
	norm = 1.0f/norm;
	s1 *= norm;
	s2 *= norm;
	s3 *= norm;
	s4 *= norm;

	qDot1 = 0.5f * (-q2 * gx - q3 * gy - q4 * gz) - beta * s1;
	qDot2 = 0.5f * (q1 * gx + q3 * gz - q4 * gy) - beta * s2;
	qDot3 = 0.5f * (q1 * gy - q2 * gz + q4 * gx) - beta * s3;
	qDot4 = 0.5f * (q1 * gz + q2 * gy - q3 * gx) - beta * s4;

	q1 += qDot1 * deltat;
	q2 += qDot2 * deltat;
	q3 += qDot3 * deltat;
	q4 += qDot4 * deltat;
	norm = sqrt(q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4);
	norm = 1.0f/norm;
	q[0] = q1 * norm;
	q[1] = q2 * norm;
	q[2] = q3 * norm;
	q[3] = q4 * norm;
}

void IMUGY85::getMres() {
	switch (Mscale)
	{
		case MFS_14BITS:
			mRes = 10.*4219./8190.; // Proper scale to return milliGauss
			break;
		case MFS_16BITS:
			mRes = 10.*4219./32760.0; // Proper scale to return milliGauss
			break;
	}
}

void IMUGY85::getGres() {
	switch (Gscale)
	{
		case GFS_250DPS:
			gRes = 250.0/32768.0;
			break;
		case GFS_500DPS:
			gRes = 500.0/32768.0;
			break;
		case GFS_1000DPS:
			gRes = 1000.0/32768.0;
			break;
		case GFS_2000DPS:
			gRes = 2000.0/32768.0;
			break;
		case GFS_CUSTOM:
			gRes = 2400.0/32768.0;
			break;
	}
}

void IMUGY85::getAres() {
	switch (Ascale)
	{
		case AFS_2G:
			aRes = 2.0/32768.0;
			break;
		case AFS_4G:
			aRes = 4.0/32768.0;
			break;
		case AFS_8G:
			aRes = 8.0/32768.0;
			break;
		case AFS_16G:
			aRes = 16.0/32768.0;
			break;
	}
}