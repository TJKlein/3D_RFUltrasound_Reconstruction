#ifndef RECONSTRUCTION_H
#define RECONSTRUCTION_H
#include <Eigen/StdVector>
#include <hash_map>

#include <Eigen\Core>
#include <Eigen/Geometry>
#include <Eigen/Array>
#include <Eigen/QR>
#include <Eigen/LU>
#include <Eigen/SVD>
#include <windows.h>
//#include <boost/serialization/base_object.hpp>
//#include <boost/archive/text_oarchive.hpp> 

#define MATLAB_SUPPORT 0

class ColorLUT
{
private:
	std::vector<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f>>  m_LUT;
	std::vector<float> m_velocity;
	stdext::hash_map<DWORD, float> m_hashMapVelocity;
	float m_maxVelocity;
public:
	bool colorToVelocity(const int &R, const int &G, const int &B, float &velocity);
	ColorLUT(std::vector<int> interpolation, std::vector<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f>> color, const float &maxVelocity);
	~ColorLUT();
};


//typedef Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, int> DynamicMatrixf;

class VelocityData
{
public:
	//friend class boost::serialization::access;
	Eigen::MatrixXf normalMatrix;
	Eigen::VectorXf velocityVector;
	Eigen::VectorXi frameVector;
	float approximationError;

	/*template<class Data>
	friend void serialize(Data &ar, const unsigned int version) 
	{ 
		ar & age_; 
	} */

};

class Reconstruction
{
public:
	std::vector<VelocityData*> m_measurementField;
std::vector<Eigen::Vector3f> m_normalPlane;
private:
	
	std::vector<Eigen::Vector3f> m_spacePlane;
	Eigen::Vector3f *m_velocityField;
	Eigen::Vector3f *m_velocityFieldUpdated;
	bool *m_velocityFieldValid;

	// backup data structures
	
	Eigen::Vector3f *m_velocityFieldBackup;
	bool *m_velocityFieldValidBackup;


	Eigen::Transform3f m_velocityField2World, m_world2velocityField;
	int m_dimX, m_dimY, m_dimZ;
	float m_scaleX, m_scaleY, m_scaleZ;
	float m_maxError;
	float m_usScaleX, m_usScaleY;

	int m_roiHeight, m_roiWidth;

	bool m_calibrationMatrixSet;
	Eigen::Transform3f m_calibrationMatrix;
	Eigen::Vector4f m_centroid;
	Eigen::Vector3f m_minVector;
	Eigen::Vector3f m_maxVector;

	
public:

	EIGEN_MAKE_ALIGNED_OPERATOR_NEW 

		void setMaxError(const float &error);
		float getMaxError();
		VelocityData* getMeasurementsByCoordinates(const int &x, const int &y, const int &z);
		VelocityData* getMeasurementsByIndex(int index);
		Eigen::Vector3f getVelocityVectorByCoordinates(const int &x, const int &y, const int &z);
		Eigen::Vector3f getVelocityVectorByIndex(const int &index);
		void setVelocityVectorByIndex(const int &index, Eigen::Vector3f &_vec);

		float partialDerivativeComponent(int x, int y, int z, const int &velocityComponent, const int &x1, const int &x2, const float &hx, const float &hy);
		bool updateDeltaVelocity(const int &x, const int &y, const int &z, const Eigen::Vector3f &velocity);
		bool updateDeltaVelocity(const int &index, const Eigen::Vector3f &velocity);
		bool getApproximationTemporaryError(const int &x, const int &y, const int &z, Eigen::VectorXf &error);
		bool getProjection(const int &x, const int &y, const int &z, Eigen::VectorXf &solutionVector);
		bool getApproximationTemporaryError(const int &index, Eigen::VectorXf &error);
		bool getProjection(const int &index, Eigen::VectorXf &solutionVector);
		bool getApproximationError(const int &x, const int &y, const int &z, Eigen::VectorXf &error);
		bool getApproximationError(const int &index, Eigen::VectorXf &error);
		void performVectorFieldUpdate();
		bool updateVelocity(const int &x, const int &y, const int &z, const Eigen::Vector3f &velocity);
		bool getVelocityFieldValidity(const int &x, const int &y, const int &z);
		bool getVelocityFieldValidity(const int &index);
		bool setVelocityFieldValidity(const int &x, const int &y, const int &z, const bool &validity);
		bool getIndex(const int &x, const int &y, const int &z, int &index);
		bool getIndex(const Eigen::Vector3f& worldPosition, int &index);
		void getCoordinatesFromIndex(const int &index, int &x, int &y, int &z);
		Eigen::Vector3f* getVelocityField() { return m_velocityField; }
		Eigen::Vector3f* getTemporaryVelocityField() { return m_velocityFieldUpdated; }
	std::vector<Eigen::Vector3f>* getNormalPlane() { return &m_normalPlane; }
	std::vector<Eigen::Vector3f>* getSpacePlane() { return &m_spacePlane; }

	Eigen::Vector3f computeDivergence(const int &x, const int &y, const int &z, const float &physicalX, const float &physicalY, const float &physicalZ);
	Eigen::Vector3f computeLaplacian(const int &x, const int &y, const int &z, const float &physicalX, const float &physicalY, const float &physicalZ);
	
#if MATLAB_SUPPORT
	void saveNormalPlaneAsMAT();
	void saveMatlabMatrix1Dim(const Eigen::MatrixXf &data, const int &width, const int &height, char* name, char* filename);
	void saveMatlabMatrix3Dim(std::vector<Eigen::Vector3f> &data, const int &width, const int &height, char* name, char* filename);
#endif

	Reconstruction();
	~Reconstruction();

	/**
	\brief Creates a backup of the vector field
	**/
	void createFieldBackup();

	/**
	\bried Restores a previously created backup of the vector field
	**/
	void restoreFieldBackup();

	/**
	\brief Returns the precomputed normal of the corresponding pixel in the region of interest of the ultrasound image
	\param x X coordinate within the region of interest of the ultrasound image
	\param y Y coordinate within the region of interest of the ultrasound image
	\param normal Precomputed normal value will be written here if successful
	**/
	bool getNormal(int x, int y, Eigen::Vector3f &normal);

	/**
	\brief Returns the precomputed space position of the corresponding pixel in the region of interest of the ultrasound image
	\param x X coordinate within the region of interest of the ultrasound image
	\param y Y coordinate within the region of interest of the ultrasound image
	\param space Precomputed value will be written here if successful
	**/
	bool getSpace(int x, int y, Eigen::Vector3f &space);

	/**
	\brief Computes the normal vectors and space positions for each pixel in the region of interest of the ultrasound image
	\param apexX X coordinate position of the ultrasound apex in mm
	\param apexY Y coordinate position of the ultrasound apex in mm
	\param minRoiX Minimum X coordinate of the region of interest in the ultrasound image
	\param minRoiY Minimum Y coordinate of the region of interest in the ultrasound image
	\param maxRoiX Maximum X coordinate of the region of interest in the ultrasound image
	\param maxRoiY Maximum Y coordinate of the region of interest in the ultrasound image
	**/
	bool preComputePlanes(const float &apexX, const float &apexY, const int &minRoiX, const int &minRoiY, const int &maxRoiX, const int &maxRoiY);
	

	bool approximateVelocityVector(const Eigen::Vector3f &worldPosition, Eigen::VectorXf &velocity);
	bool approximateVelocityVector(const int &x, const int &y, const int &z, Eigen::VectorXf &velocity, float &error, float &conditionNumber, int &measurements);
	bool getMeasurementData(const Eigen::Vector3f &position, Eigen::MatrixXf &normalMatrix, Eigen::VectorXf &velocityVector);

	/**
	\brief Sets the calibration data of the ultrasound transducer relative to the tracking target.
	\param calibrationMatrix Transformation from the image plane (mm) to the tracking target. Origin in the image is assumed to be the left upper corner of the image.
	\param scaleX Scaling parameter in X direction defining the number of pixels per mm
	\param scaleY Scaling parameter in Y direction defining the number of pixels per mm
	**/
	void setCalibration(Eigen::Matrix4f &calibrationMatrix, const float &scaleX, const float &scaleY);

	void generateVelocityField(const int &dimX, const int &dimY, const int &dimZ, const float &scaleX, const float &scaleY, const float &scaleZ, const Eigen::Matrix4f &origin);

	bool getClosestVelocity(const Eigen::Vector3f &worldPosition, Eigen::Vector3f &velocity);

	bool getClosestMeasurements(const Eigen::Vector3f &worldPosition, VelocityData * velocityData);
	
	bool pushMeasurement(const Eigen::Vector3f &normal, const float &velocity, const int &index, const int &frame);

	bool setVelocity(const int &dimX, const int &dimY, const int &dimZ, const Eigen::Vector3f &velocity, bool validity = true);

};


#endif