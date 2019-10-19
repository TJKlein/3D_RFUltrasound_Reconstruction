#ifndef DATATOOLS_H
#define DATATOOLS_H

#include <Eigen/StdVector>
#include <iostream>
#include <fstream>
#include <time.h>
#include <string>
//#include "Common/DataTypes.h"
#include <Eigen/Core>
#include <Eigen/Geometry>
//#include <vector>

struct TrackingData
{
	float quaternion[4];
	float translation[3];
	float scale[3];
	float error;
	float timestamp;
	int isValid;
	std::string targetID;
};

class DataTools
{
	public:
		DataTools();
		~DataTools();
		std::string computeElapsedTime(const clock_t &startTimer, const clock_t &stopTimer);
		bool getTargetNamesFromFile(std::string filename, std::vector<std::string> &targetNameVector);
		bool writeTrackingDataToFile(std::vector<Eigen::Transform3f, Eigen::aligned_allocator<Eigen::Transform3f>> &_data, std::string targetName, std::string filename);
		bool getTrackingDataFromFileRelative(std::string filename, std::string target, std::string referenceTarget, std::vector<TrackingData> &trackingVector);
		bool getTrackingDataFromFile(std::string filename, std::string target, std::vector<TrackingData> &trackingVector);
		bool getVideoTimestampsFromFile(std::string filename, std::vector<float> &timestampVector);
		//TJK_removed_CAMP: void eigen2CAMP(const Eigen::Matrix4f &eigMat, CAMP::Matrix4<float> &campMat);
		//TJK_removed_CAMP: void eigen2CAMP(const Eigen::Matrix3f &eigMat, CAMP::Matrix3<float> &ampMat);
		/*void eigen2CAMP(const Eigen::Matrix4f &eigMat, CAMP::Matrix4<float> &campMat);
		void eigen2CAMP(const Eigen::Matrix3f &eigMat, CAMP::Matrix3<float> &ampMat);*/
}; 

#endif