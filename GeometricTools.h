#ifndef GEOMETRICTOOLS_H
#define GEOMETRICTOOLS_H
#include <Eigen/StdVector>
#include <Eigen\Core>
#include <Eigen/Geometry>
#include <Eigen/Array>
#include <Eigen/QR>
#include <Eigen/LU>
//TJK_removed_CAMP: #include "Common/DataTypes.h"

#define PI 3.1415926536


class CostFunction {
public:
	/// Virtual destructor definition
	virtual ~CostFunction() {}

	/// Method called by optimization algorithms
	virtual double evaluate(double *params) = 0;

	/// Method called by optimization algorithms for vector-valued results
	virtual void evaluate(double *params, double *result) { result[0] = evaluate(params); }

	/// Method called for parallel optimization
	virtual void evaluateParallel(int n, double *params, double *results) { }

	/// Indicates if parallel optimization is supported
	virtual bool parallel() { return false; }

};


class Optimizer {
public:
	enum Success {
		CONVERGED_FUNC = 0,
		CONVERGED_PARS = 1,
		FAILED_ITER    = 2,
		FAILED_EVAL    = 3,
		UNKNOWN        = 4,
		WRONG_PARAMS   = 5,
		STEP_INF	   = 6,
		STEP_NAN	   = 7,
		STEP_NOT_SOLVED = 8
	};

	enum OptimizerType {
		BESTNEIGHBOR = 0,
		POWELLBRENT  = 1,
		EXHAUSTIVE   = 2,
		GRADIENT	 = 3,
		GAUSSNEWTON  = 4,
		LEVMAR       = 5,
		GENETIC      = 6,
		SIMPLEX		 = 7,
		BESTNEIGHBORFEATURED = 8,
		LEVMARWEIGHTED = 9
	};

	/// Constructor with number of dimensions of the cost function
	Optimizer(int dimensions, CostFunction *cf = 0);
	virtual ~Optimizer();

	/// Performs the optimization
	virtual int run() = 0;

	/// Factory method for convenience
	//static Optimizer *createOptimizer(int type, int dimensions);

	/// Factory method for convenience
	static Optimizer *createOptimizer(OptimizerType type, int paramDimensions, int resultDimensions=1);

	/// @name Access Methods
	//@{
	virtual void setCostFunction(CostFunction *cf) { m_cf = cf; }
	virtual void setParTol(double x) { m_ptol = x; }
	virtual void setFunTol(double x) { m_ftol = x; }
	virtual void setParameters(double *params) { m_params = params; }
	virtual void setStepSizes(double *stepsizes);
	virtual void setMaxIter(int n) { m_maxiter = n; }
	virtual void setMaxEval(int n) { m_maxeval = n; }
	virtual void setMinimize(bool flag) { m_minimize = flag; }
	//virtual void setParamSequence(double **ps) { m_paramSequence = ps; }
	virtual void setCbIter(void (*cbIter)(void*, void*), void *param) { m_cbIter = cbIter; m_paramCbIter = param; }
	void setParamsOrder(int *order);


	virtual double getParTol() { return m_ptol; }
	virtual double getFunTol() { return m_ftol; }
	virtual double getCurrentValue() {return m_currval; }
	virtual double getInitialValue() {return m_initval; }
	virtual double* getParameters() { return m_params; }
	virtual int getIter() { return m_iter; }
	virtual int getEval() { return m_eval; }
	virtual int getMaxIter() { return m_maxiter; }
	virtual int getMaxEval() { return m_maxeval; }
	virtual bool getRunning() { return m_running; }
	virtual int getDimensions() { return m_dimensions; }
//	virtual void* getCbIter() { return m_cbIter; }
	virtual void abort() { m_abort = true; }
	virtual int getType(){ return m_type;}
	Success getSuccess() { return m_success; }

	//@}

protected:
	/// Evaluate the cost function while scaling the parameters with the step sizes
	double evaluateScaled(double *params);

	double  m_initval;		///< initial cost function value
	double  m_currval;		///< actual cost function value
	double  m_ptol;			///< parameter tolerance
	double  m_ftol;			///< cost function tolerance
	double* m_params;		///< pointer to the parameters
	double* m_stepsizes;	///< step size for each parameter
	int     m_dimensions;	///< number of dimensions/parameters
	int     m_maxiter;		///< maximum number of iterations
	int     m_maxeval;		///< maximum number of evaluations
	int     m_iter;			///< current number of iterations
	int     m_eval;			///< current number of evaluations
	Success m_success;		///< result state of the last optimization
	bool    m_minimize;		///< if true, we are seeking a minimum
	bool    m_running;		///< is the optimization running
	bool    m_abort;		///< are we supposed to abort
	CostFunction* m_cf;		///< the cost function to minimize/maximize

	/**
	* m_cbIter, callback function
	* - if set by the user should be called by every
	*   implementation after each iteration
	* - the first parameter should be a this pointer (such that the callback
	*   has access to the optimizer)
	* - the second parameter is the pointer m_ParamCbIter that will be set 
	*   by the user and has to be passed to the function always!
	*   (this is due to the fact that the callback is a static function and
	*   in most cases wants to operate on an object or data structure)
	*
	* possible use: each parameter update can be tracked and for example
	* be used to visualize the optimization process
	*/
    void (*m_cbIter)(void*, void*);
	void *m_paramCbIter;

	int  m_type;			///< type of optimizer
	int m_paramsOrder[6];		///< order of parameters to be evaluated

};


class OptimizerBestNeighbor : public Optimizer {
public:
	OptimizerBestNeighbor(int dimensions);
	int run();
};

struct BBExtent {
	float dimX;
	float dimY;
	float dimZ;
};

class Line
{
private:
	Eigen::Vector3f m_startPoint;
	Eigen::Vector3f m_endPoint;
	Eigen::Vector3f m_direction;
	Eigen::Vector3f m_directionNormalized;
	float m_directionDotProduct;

public:
	Line(Eigen::Vector3f &startPoint, Eigen::Vector3f &endPoint);
	Eigen::Vector3f& getStartPoint() { return m_startPoint; }
	Eigen::Vector3f& getEndPoint() { return m_endPoint; }

	void distanceToPoint(const Eigen::Vector3f &point, float &distance, float &param, Eigen::Vector3f &intersectionPoint); 
	void distanceToPointMAX(const Eigen::Vector3f &point, float &distance, float &param, Eigen::Vector3f &intersectionPoint); 

};

class Plane
{
private:
	// m_plane = [n_vec d] = [n_x n_y n_z d]
	// n_vec is the normal vector on the plane, d is the distance to the origin
	Eigen::Vector4f m_planeEigen;
	Eigen::Vector3f m_normalEigen;
	//TJK_removed_CAMP CAMP::Vector4<float> m_planeCAMP;
	//TJK_removed_CAMP CAMP::Vector3<float> m_normalCAMP;
	float m_distance;
	float m_denominator;
	float m_denominatorINV;
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	// define the plane by a point on the plane and a normal vector
	Plane(Eigen::Vector3f point, Eigen::Vector3f normal);
	//TJK_removed_CAMP Plane(CAMP::Vector3<float> point, CAMP::Vector3<float> normal);

	/**
	\brief Defines a plane given a transformation matrix, describing the position and orienation of the plane origin

	\param image2World The mapping which transforms from image to world coordinate space
	*/
	Plane(const Eigen::Matrix4f& image2World);
	//TJK_removed_CAMP Plane(const CAMP::Matrix4<float>& image2World);
	~Plane();

	/**

	\brief Transforms a point from world to local coordinate system using the transformation H, which describes the mapping from image to world space.

	H^(-1)*worldPoint = point in pixel coordinate space

	\param H_worlde2image Homography that transforms 3f world coordinates to 3f image coordinates (z=0); the origin is assumed in the top left corner of the image
	\param worldPoint The point given in world coordinates that we wish to transform into local coordinate space
	\param imagePoint Reference to the world point transformed in local coordinate space

	*/
	void worldToImage(Eigen::Matrix4f &H_world2image, Eigen::Vector3f &worldPoint, Eigen::Vector3f &imagePoint);

	/**
	
	\brief Computes the intersection of the plane with a ray defined by point and direction vector

	\param point Defines a point on the ray
	\param direction Defines the direction vector of the ray
	\param result Defines the intersection point
	\result True if the line intersects with the ray, False otherwise
	*/
	bool intersectWithRay(const Eigen::Vector3f &point, const Eigen::Vector3f &direction, Eigen::Vector3f &intersection);
	bool intersectWithRay(const Eigen::Vector3f &point, const Eigen::Vector3f &direction, Eigen::Vector3f &intersection, float &distance);

	float getOrthogonalDistance(const Eigen::Vector3f &point, Eigen::Vector3f &intersection);


	//TJK_removed_CAMP: bool intersectWithRay(const CAMP::Vector3<float> &point, const CAMP::Vector3<float> &direction, CAMP::Vector3<float> &intersection);
	//TJK_removed_CAMP: bool intersectWithRay(const CAMP::Vector3<float> &point, const CAMP::Vector3<float> &direction, CAMP::Vector3<float> &intersection, float &distance);

	//TJK_removed_CAMP: float getOrthogonalDistance(const CAMP::Vector3<float> &point, CAMP::Vector3<float> &intersection);

};




class GeometricTools
{
private:
public:
	GeometricTools();
	~GeometricTools();
	void determineExtent(std::vector<Eigen::Vector3f> data, Eigen::Matrix4f &transformation, BBExtent &extent);
	void computeCentroid(std::vector<Eigen::Vector3f> data, Eigen::Matrix3f rotation, Eigen::Vector3f &centroid, Eigen::Vector3f &originInWorld, BBExtent& extent);

	void computeOrientedBoundingBox(std::vector<Eigen::Vector3f> data, BBExtent &extent, Eigen::Transform3f &trans);
	void computeManualBoundingBox(Eigen::Matrix3f &rotation, std::vector<Eigen::Vector3f> data, BBExtent &extent, Eigen::Transform3f &trans);

	void getEulerAngles(const Eigen::Quaternionf &quat, float &angleX, float &angleY, float &angleZ, bool radians = true);
};


class BoundingBox
{

protected:
	Eigen::Vector3f coordinateWorld2Local(Eigen::Vector3f point);
};


class OrientedBoundingBox : public BoundingBox
{
private:
	Eigen::Vector3f m_centroid;
	Eigen::Transform3f m_world2local;
	Eigen::Transform3f m_local2world;
	Eigen::Vector3f m_halfLengths;
	BBExtent m_extent;
public:
	/*OrientedBoundingBox(Eigen::Vector3f centroid, Eigen::Matrix3f orientation, BBExtent extent) 
	{
		m_centroid = centroid;
		m_orientation = orientation;
		m_extent = extent;
	}*/

	/**
	\brief Oriented bounding box
	\param local2world Transformation from local to world coordinate system
	\param extent Lengths of the box from the origin
	**/
	OrientedBoundingBox(Eigen::Transform3f &local2world, BBExtent extent)
	{
		m_local2world = local2world;
		m_world2local = local2world.inverse();
		m_extent = extent;
		m_halfLengths << extent.dimX/2.0f,extent.dimY/2.0f, extent.dimZ/2.0f;
		m_centroid = m_halfLengths;
	}
	~OrientedBoundingBox() {}
	bool lineSegmentIntersection(Line *line); //(const Eigen::Vector3f &point1, const Eigen::Vector3f &point2);
};

#endif