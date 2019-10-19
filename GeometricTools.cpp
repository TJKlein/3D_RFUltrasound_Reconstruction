#include "GeometricTools.h"
#include <stdio.h>
#include <iostream>
/*
#include "Numerics/OptimizerBestNeighbor.h"
#include "Numerics/OptimizerPowellBrent.h"
#include "Numerics/OptimizerSimplex.h"
#include "Numerics/CostFunction.h"
#include "Common/DataTypes.h"
*/


Optimizer::Optimizer(int dimensions, CostFunction *cf) {
	m_dimensions = dimensions;
	m_params = NULL;
	
	// init stepsizes
	m_stepsizes = new double[dimensions];
	for (int i=0; i<dimensions; i++) {
		m_stepsizes[i] = 1;
	}

	m_cf = cf;
	m_minimize = false;
	m_running = false;
	m_abort = false;
	m_iter = 0;
	m_maxiter = 100000;
	m_maxeval = 1000000;
	m_success = UNKNOWN;
	m_ftol = 0.00001;
	m_ptol = 0.00001;
	m_cbIter = NULL;
	m_currval = -1;
}


Optimizer::~Optimizer() {
	delete [] m_stepsizes;
}


void Optimizer::setStepSizes(double *stepsizes) {
	for (int i = 0; i < m_dimensions; i++) m_stepsizes[i] = stepsizes[i];
}

 //This is giving problems... Sorry, but I needed to comment it, Thomas
void Optimizer::setParamsOrder(int *paramsOrder) {
	for (int i = 0; i < 6; i++) m_paramsOrder[i] = paramsOrder[i];
}


Optimizer *Optimizer::createOptimizer(OptimizerType type, int paramDimensions, int resultDimensions) {
	
	/*
	For optimizers 0, 1, 2, 3, 6, 7, 8 it is only required the dimensions of parameters
	For optimizers 4, 5, 9 it is required the dimensions of parameters AND the dimensions of result
	*/

	switch (type) {
		case BESTNEIGHBOR: return new OptimizerBestNeighbor(paramDimensions); break;
		/*case EXHAUSTIVE: return new OptimizerScan(paramDimensions,resultDimensions); break;
		case SIMPLEX: return new OptimizerSimplex(paramDimensions); break;
		case BESTNEIGHBORFEATURED: return new OptimizerBestNeighborFeatured(paramDimensions); break;
		case LEVMAR: return new OptimizerLevenbergMarquardt(resultDimensions,paramDimensions); break;
		case GAUSSNEWTON: return new OptimizerGaussNewton(resultDimensions,paramDimensions); break;
		case LEVMARWEIGHTED: return new OptimizerLMWeighted(resultDimensions,paramDimensions); break;
		case POWELLBRENT: return new OptimizerPowellBrent(paramDimensions); break;
		case GRADIENT: return new OptimizerGradient(paramDimensions); break;
		case GENETIC: return new OptimizerGenetic(paramDimensions); break;*/
		default: return NULL; break;
	}
}

double Optimizer::evaluateScaled(double *params) {
	for (int i = 0; i < m_dimensions; i++)
		m_params[i] = params[i] * m_stepsizes[i];
	return (m_minimize)? m_cf->evaluate(m_params): -m_cf->evaluate(m_params);
}


OptimizerBestNeighbor::OptimizerBestNeighbor(int dimensions)
	: Optimizer(dimensions) {}

int OptimizerBestNeighbor::run() {
	m_running = true;
	m_abort = false;
	m_iter = 0;
	m_eval = 0;
	int i;
	bool downscaled = false, better = false;
	double lastvalue = 0, currvalue = 0, scale = 1.0;
	double *workparams = new double[(m_cf->parallel()? 2*m_dimensions: 1) * m_dimensions];
	double *workresults = new double[2*m_dimensions];
	for (i = 0; i < m_dimensions; i++) workparams[i] = m_params[i];

	// determine smallest step
/*	double smallest_step = m_stepsizes[0];
	for (i = 1; i < m_dimensions; i++)
		if ((m_stepsizes[i]) && (m_stepsizes[i] < smallest_step))
			smallest_step = m_stepsizes[i]; */

	// evaluate the current position first
	double c = m_cf->evaluate(workparams);
	m_currval = c;
	m_initval = c;
	// negate result if we are minimizing
	currvalue = (m_minimize)? -c: c;

	// main optimization loop
	while ((m_iter < m_maxiter) && (m_eval < m_maxeval)
	&& ((m_iter == 0) || (downscaled) || (currvalue - lastvalue > m_ftol))) {

		// save last result
		lastvalue = currvalue;
		int bestpos = 0;

		// evaluate neighborhood

		// For Thorben: Parallel Evaluation
		if (m_cf->parallel()) {
			int i, nsteps = 0;
			// copy current parameters into all 2n pars for parallelization
			for (i = 0; i < 2*m_dimensions*m_dimensions; i++)
				workparams[i] = m_params[i%m_dimensions];

			for (i = 0; i < 2*m_dimensions; i++) {
				// neglect degrees of freedom, if the stepsize if zero
				if (m_stepsizes[i%m_dimensions]) {
					// move into specific direction
					workparams[m_dimensions*nsteps+i%m_dimensions]
					+= m_stepsizes[i%m_dimensions] * scale * (i < m_dimensions? 1.0: -1.0);
					nsteps++;
				}
			}
			m_cf->evaluateParallel(nsteps, workparams, workresults);
			m_eval += nsteps;
			// figure out the best result
			for (i = 0; i < nsteps; i++) {
				if (m_minimize) workresults[i] = -workresults[i];
				if (i) {
					if (currvalue < workresults[i]) {
						currvalue = workresults[i]; bestpos = i;
					}
				} else currvalue = workresults[i];
			}
		} else {
		// sequential evaluation
			for (int i = 0; i < 2*m_dimensions; i++) {

				// neglect degrees of freedom, if the stepsize if zero
				if (m_stepsizes[i%m_dimensions]) {
					// move into specific direction
					workparams[i%m_dimensions] += m_stepsizes[i%m_dimensions]
					* scale * (i < m_dimensions? 1.0: -1.0);

					// calculate cost function
					double c = m_cf->evaluate(workparams);
					// negate result if we are minimizing
					workresults[i] = (m_minimize)? -c: c;
					m_eval++;

					// reset work parameter
					workparams[i%m_dimensions] = m_params[i%m_dimensions];
				}
				else workresults[i] = lastvalue;

				// check if result is better than the last one
				if (i) {
					if (currvalue < workresults[i]) {
						currvalue = workresults[i]; bestpos = i;
					}
				} else currvalue = workresults[i];
			}
		}

		// If best result is better, proceed
		if (currvalue > lastvalue) {

			// Try combination of different DOFs
			// Only for those who have different signs in each directions
/*			int combined = 0;
			for (int j = 0; j < m_dimensions; j++) {
				if ((m_stepsizes[j])
				&& (((workresults[j] > lastvalue) && (workresults[j+m_dimensions] < lastvalue))
				||  ((workresults[j] < lastvalue) && (workresults[j+m_dimensions] > lastvalue)))) {
					workparams[j] += m_stepsizes[j] * m_scale * ((workresults[j] > lastvalue)? 1.0: -1.0);
					combined++;
				}
			}
			// Evaluate combination of the steps!
			double c;
			if (combined > 1) {
				c = m_cf->evaluate(workparams);
				if (m_minimize) c = -c;
			}
				// is this even better?
			if ((combined > 1) && (c > currvalue)) {
				for (int j = 0; j < m_dimensions; j++)
					m_params[j] = workparams[j];
				m_costFunctionValue = (m_minimize)? -c: c;
				better = true;
				downscaled = false;
				sprintf(last_direction, "combined");
			} else { */
				// set this direction permanently
				m_params[bestpos%m_dimensions] += m_stepsizes[bestpos%m_dimensions]
				 * scale * (bestpos < m_dimensions? 1.0: -1.0);
				for (int j = 0; j < m_dimensions; j++) workparams[j] = m_params[j];
				m_currval = (m_minimize)? -currvalue: currvalue;
				better = true;
				downscaled = false;
//			}
		}
		// otherwise scale the step size down
		else {
			currvalue = lastvalue;
			scale *= 0.5;
			downscaled = true;
			if (/*smallest_step * */ scale < m_ptol) break;
		}

		if (m_cbIter != NULL) {
			m_cbIter(this, m_paramCbIter);
		}

		m_iter++;
		if (m_abort) break;
	}

	if (m_iter >= m_maxiter) m_success = FAILED_ITER;
	else if (m_eval >= m_maxeval) m_success = FAILED_EVAL;
	else if (currvalue - lastvalue <= m_ftol) m_success = CONVERGED_FUNC;
	else if (/* smallest_step * */ scale < m_ptol) m_success = CONVERGED_PARS;
	else m_success = UNKNOWN;

	delete [] workparams;
	delete [] workresults;
	m_running = false;
	return better;
}




/*float Plane::getOrthogonalDistance(const CAMP::Vector3<float> &point, CAMP::Vector3<float> &intersection)
{
	float _distance = 0.0;
	intersectWithRay(point, m_normalCAMP,intersection, _distance);

	return _distance;
}*/



float Plane::getOrthogonalDistance(const Eigen::Vector3f &point, Eigen::Vector3f &intersection)
{
	float _distance = 0.0;
	intersectWithRay(point, m_normalEigen,intersection, _distance);

	return _distance;
}

/*
bool Plane::intersectWithRay(const CAMP::Vector3<float> &point, const CAMP::Vector3<float> &direction, CAMP::Vector3<float> &intersection, float &distance)
{
	if ( intersectWithRay(point, direction, intersection) )
	{
		distance = (intersection-point).norm();
		return true;
	}
	else
		return false;
}*/


bool Plane::intersectWithRay(const Eigen::Vector3f &point, const Eigen::Vector3f &direction, Eigen::Vector3f &intersection, float &distance)
{
	if ( intersectWithRay(point, direction, intersection) )
	{
		distance = (intersection-point).norm();
		return true;
	}
	else
		return false;
}

bool Plane::intersectWithRay(const Eigen::Vector3f &point, const Eigen::Vector3f &direction, Eigen::Vector3f &intersection)
{
  const float _nominator = -(m_normalEigen.dot(point) + m_distance);
  intersection.noalias() = (_nominator * m_denominatorINV * direction);
  intersection += point;

 return true;
}

/*
bool Plane::intersectWithRay(const CAMP::Vector3<float> &point, const CAMP::Vector3<float> &direction, CAMP::Vector3<float> &intersection)
{
	float _nominator = -(CAMP::dot(m_normalCAMP,point) + m_distance);

	// check if ray and plane are orthogonal
	//if ( m_denominator <= 0.0000001 )
	///	return false;

	intersection = _nominator * m_denominatorINV * direction + point;

	return true;
}*/

void Plane::worldToImage(Eigen::Matrix4f &H_world2image, Eigen::Vector3f &worldPoint, Eigen::Vector3f &imagePoint)
{
	Eigen::Vector4f _tempHomogenous;
	_tempHomogenous << worldPoint(0),worldPoint(1),worldPoint(2),1.0;

	_tempHomogenous=H_world2image*_tempHomogenous;

	imagePoint << _tempHomogenous(0),_tempHomogenous(1),_tempHomogenous(2);
}


Plane::Plane(const Eigen::Matrix4f& image2world)
{

	Eigen::Vector3f _origin;
	Eigen::Vector4f _tempNormal;
	Eigen::Vector4f _tempOrigin;

	_tempOrigin << 0,0,0,1;
	_tempNormal << 0,0,1,1;

	_tempOrigin = (image2world * _tempOrigin);
	_tempNormal	= (image2world * _tempNormal)-_tempOrigin;

	m_normalEigen << _tempNormal(0),_tempNormal(1),_tempNormal(2);
	m_normalEigen = m_normalEigen.normalized();

	_origin << _tempOrigin(0),_tempOrigin(1),_tempOrigin(2);


	// d = m_plane(3), distance to the origin
	m_planeEigen(3) = -m_normalEigen.dot( _origin );
	m_distance = m_planeEigen(3);

	// plug in the normal vector components
	m_planeEigen(0) = m_normalEigen(0);
	m_planeEigen(1) = m_normalEigen(1);
	m_planeEigen(2) = m_normalEigen(2);

	m_denominator = m_normalEigen.dot( m_normalEigen);
	m_denominatorINV = 1.0/m_denominator;

	//TJK_removed_CAMP: m_planeCAMP = CAMP::Vector4<float>(m_planeEigen(0), m_planeEigen(1), m_planeEigen(2), m_planeEigen(3));
	//TJK_removed_CAMP: m_normalCAMP = CAMP::Vector3<float>(m_normalEigen(0), m_normalEigen(1), m_normalEigen(2));
}


Plane::Plane(Eigen::Vector3f point, Eigen::Vector3f normal)
{
	// d = m_plane(3), distance to the origin
	m_planeEigen(3) = -normal.dot( point );
	m_distance = m_planeEigen(3);

	// plug in the normal vector components
	m_planeEigen(0) = normal(0);
	m_planeEigen(1) = normal(1);
	m_planeEigen(2) = normal(2);

	
	m_normalEigen = normal;

	//TJK_removed_CAMP: m_planeCAMP = CAMP::Vector4<float>(m_planeEigen(0), m_planeEigen(1), m_planeEigen(2), m_planeEigen(3));
	//TJK_removed_CAMP: m_normalCAMP = CAMP::Vector3<float>(normal(0), normal(1), normal(2));

	m_denominator = m_normalEigen.dot( m_normalEigen);
	m_denominatorINV = 1.0/m_denominator;
}

/*
Plane::Plane(const CAMP::Matrix4<float>& image2world)
{

	CAMP::Vector3<float> _origin;
	CAMP::Vector4<float> _tempNormal(0,0,1,1);
	CAMP::Vector4<float> _tempOrigin(0,0,0,1);


	_tempOrigin = (image2world * _tempOrigin);
	_tempNormal	= (image2world * _tempNormal)-_tempOrigin;

	m_normalCAMP = CAMP::Vector3<float>( _tempNormal.c_array()[0],_tempNormal.c_array()[1],_tempNormal.c_array()[2]);
	m_normalCAMP.normalize();

	_origin = CAMP::Vector3<float>( _tempOrigin.c_array()[0],_tempOrigin.c_array()[1],_tempOrigin.c_array()[2]);


	// d = m_plane(3), distance to the origin
	m_planeCAMP.c_array()[3] = CAMP::dot(-m_normalCAMP, _origin );
	m_distance = m_planeCAMP.c_array()[3];

	// plug in the normal vector components
	m_planeCAMP.c_array()[0] = m_normalCAMP.c_array()[0];
	m_planeCAMP.c_array()[1] = m_normalCAMP.c_array()[1];
	m_planeCAMP.c_array()[2] = m_normalCAMP.c_array()[2];

	m_denominator = CAMP::dot(m_normalCAMP, m_normalCAMP);
	m_denominatorINV = 1.0/m_denominator;

	m_planeEigen << m_planeCAMP.c_array()[0],m_planeCAMP.c_array()[1], m_planeCAMP.c_array()[2], m_planeCAMP.c_array()[3];
	m_normalEigen << m_normalCAMP.c_array()[0],m_normalCAMP.c_array()[1],m_normalCAMP.c_array()[2];
}


Plane::Plane(CAMP::Vector3<float> point, CAMP::Vector3<float> normal)
{
	// d = m_plane(3), distance to the origin
	m_planeCAMP.c_array()[3] =CAMP::dot(-normal, point );
	m_distance = m_planeCAMP.c_array()[3];

	// plug in the normal vector components
	m_planeCAMP.c_array()[0] = normal.c_array()[0];
	m_planeCAMP.c_array()[1] = normal.c_array()[1];
	m_planeCAMP.c_array()[2] = normal.c_array()[2];

	
	m_normalCAMP = normal;

	m_denominator = CAMP::dot(m_normalCAMP, m_normalCAMP);
	m_denominatorINV = 1.0/m_denominator;

	m_planeEigen << m_planeCAMP.c_array()[0],m_planeCAMP.c_array()[1], m_planeCAMP.c_array()[2], m_planeCAMP.c_array()[3];
	m_normalEigen << m_normalCAMP.c_array()[0],m_normalCAMP.c_array()[1],m_normalCAMP.c_array()[2];
}*/


Plane::~Plane()
{
}



class OBB_Computation : CostFunction
{
private:
	std::vector<Eigen::Vector3f> m_data;
public:

	OBB_Computation(std::vector<Eigen::Vector3f> &data) { m_data = data;}
	double evaluate(double* params)
	{
		//Eigen::Vector3f _trans;
		//Eigen::Quaternionf _rot = Eigen::Quaternionf( (float)params[0], (float)params[1], (float)params[2], (float)params[3]);

		Eigen::Matrix3f _temp;
		_temp = Eigen::AngleAxisf((float)params[2], Eigen::Vector3f::UnitZ()) * Eigen::AngleAxisf((float)params[1], Eigen::Vector3f::UnitY()) * Eigen::AngleAxisf((float)params[0], Eigen::Vector3f::UnitX());

		float minX,minY,minZ;
		float maxX,maxY,maxZ;

		minX = std::numeric_limits<float>::max();
		minY = std::numeric_limits<float>::max();
		minZ = std::numeric_limits<float>::max();
		maxX = -std::numeric_limits<float>::max();
		maxY = -std::numeric_limits<float>::max();
		maxZ = -std::numeric_limits<float>::max();


		for (int i=0; i<m_data.size();i++)
		{
			Eigen::Vector3f _transformed = _temp*m_data[i];


			if ( _transformed(0) > maxX )
				maxX = _transformed(0);
			if ( _transformed(1) > maxY )
				maxY = _transformed(1);
			if ( _transformed(2) > maxZ )
				maxZ = _transformed(2);
			if ( _transformed(0) < minX )
				minX = _transformed(0);
			if ( _transformed(1) < minY )
				minY = _transformed(1);
			if ( _transformed(2) < minZ )
				minZ = _transformed(2);
		}

		// return the volume of the bounding box
		double volume =  (maxX-minX)*(maxY-minY)*(maxZ-minZ);

		if ( volume < 0 )
			volume =  std::numeric_limits<double>::max();

		return volume;
	}
};


GeometricTools::GeometricTools()
{
}

GeometricTools::~GeometricTools()
{
}

void GeometricTools::determineExtent(std::vector<Eigen::Vector3f> data, Eigen::Matrix4f &transformation, BBExtent &extent)
{
	
	

	float minX,minY,minZ;
	float maxX,maxY,maxZ;

	minX = std::numeric_limits<float>::max();
	minY = std::numeric_limits<float>::max();
	minZ = std::numeric_limits<float>::max();
	maxX = -std::numeric_limits<float>::max();
	maxY = -std::numeric_limits<float>::max();
	maxZ = -std::numeric_limits<float>::max();

	int index[6];
	for (int i=0; i<data.size();i++)
		{
			Eigen::Vector4f _tmpVec;
			_tmpVec << data[i].x(), data[i].y(),data[i].z(), 1.0;
			_tmpVec = transformation*_tmpVec;
			Eigen::Vector3f _transformed;
			_transformed << _tmpVec.x(),_tmpVec.y(),_tmpVec.z();

			if ( _transformed(0) > maxX )
			{
				maxX = _transformed(0);
				index[0] = i;
			}
			if ( _transformed(1) > maxY )
			{
				maxY = _transformed(1);
				index[1] = i;
			}
			if ( _transformed(2) > maxZ )
			{			
				maxZ = _transformed(2);
				index[2] = i;
			}
			if ( _transformed(0) < minX )
			{
				minX = _transformed(0);
				index[3] = i;
			}
			if ( _transformed(1) < minY )
			{
				minY = _transformed(1);
				index[4] = i;
			}
			if ( _transformed(2) < minZ )
			{
				minZ = _transformed(2);
				index[5] = i;
			}
	}

	Eigen::Vector3f _bbCorner[8];

	_bbCorner[0] << minX,minY,minZ;
	_bbCorner[1] << maxX,minY,minZ;
	_bbCorner[2] << maxX,maxY,minZ;
	_bbCorner[3] << minX,maxY,minZ;
	
	_bbCorner[4] << minX,minY,maxZ;
	_bbCorner[5] << maxX,minY,maxZ;
	_bbCorner[6] << maxX,maxY,maxZ;
	_bbCorner[7] << minX,maxY,maxZ;

	extent.dimX = maxX - minX;
	extent.dimY = maxY - minY;
	extent.dimZ = maxZ - minZ;
}

void GeometricTools::computeCentroid(std::vector<Eigen::Vector3f> data, Eigen::Matrix3f rotation, Eigen::Vector3f &centroid, Eigen::Vector3f &originInWorld, BBExtent &extent)
{
	
	Eigen::Matrix3f _rotation;
	_rotation = rotation.inverse();

	float minX,minY,minZ;
	float maxX,maxY,maxZ;

	minX = std::numeric_limits<float>::max();
	minY = std::numeric_limits<float>::max();
	minZ = std::numeric_limits<float>::max();
	maxX = -std::numeric_limits<float>::max();
	maxY = -std::numeric_limits<float>::max();
	maxZ = -std::numeric_limits<float>::max();

	int index[6];
	for (int i=0; i<data.size();i++)
		{
			Eigen::Vector3f _transformed = rotation*data[i];


			if ( _transformed(0) > maxX )
			{
				maxX = _transformed(0);
				index[0] = i;
			}
			if ( _transformed(1) > maxY )
			{
				maxY = _transformed(1);
				index[1] = i;
			}
			if ( _transformed(2) > maxZ )
			{			
				maxZ = _transformed(2);
				index[2] = i;
			}
			if ( _transformed(0) < minX )
			{
				minX = _transformed(0);
				index[3] = i;
			}
			if ( _transformed(1) < minY )
			{
				minY = _transformed(1);
				index[4] = i;
			}
			if ( _transformed(2) < minZ )
			{
				minZ = _transformed(2);
				index[5] = i;
			}
	}

	Eigen::Vector3f _bbCorner[8];

	_bbCorner[0] << minX,minY,minZ;
	_bbCorner[1] << maxX,minY,minZ;
	_bbCorner[2] << maxX,maxY,minZ;
	_bbCorner[3] << minX,maxY,minZ;
	
	_bbCorner[4] << minX,minY,maxZ;
	_bbCorner[5] << maxX,minY,maxZ;
	_bbCorner[6] << maxX,maxY,maxZ;
	_bbCorner[7] << minX,maxY,maxZ;


	Eigen::Vector3f _centroid;
	_centroid << 0,0,0;

	for (int i=0;i<8;i++)
	{

		_centroid += _bbCorner[i];
	}

	_centroid = _centroid / 8.0;

	// centroid in BB coordinate system
	centroid = /*_rotation **/ _centroid;

	originInWorld << minX,minY,minZ;

	extent.dimX = maxX - minX;
	extent.dimY = maxY - minY;
	extent.dimZ = maxZ - minZ;


}

void GeometricTools::computeOrientedBoundingBox(std::vector<Eigen::Vector3f> data, BBExtent &extent, Eigen::Transform3f &trans)
{
	OBB_Computation* _costfunction = new OBB_Computation(data);

	double params[]={0.0,0.0,0.0};

	Optimizer* optimizer = new OptimizerBestNeighbor(3); 

	optimizer->setCostFunction((CostFunction*)_costfunction);
	optimizer->setMinimize(true);
	optimizer->setMaxIter(1000000);
	optimizer->setParameters(params);

	optimizer->run();

	// the optimizer returns the rotation parameters that transform from WORLD to LOCAL
	// Local means the optimized bounding box (axis aligned)
	// _rot orientation is world to local
	
	Eigen::Matrix3f _rot;
	_rot = Eigen::AngleAxisf((float)params[2], Eigen::Vector3f::UnitZ()) * Eigen::AngleAxisf((float)params[1], Eigen::Vector3f::UnitY()) * Eigen::AngleAxisf((float)params[0], Eigen::Vector3f::UnitX());


	//std::cout << "Params: " << params[2] << " / " << params[1] << " / " << params[0] << std::endl;

	Eigen::Vector3f _originInWorld;
	Eigen::Vector3f _centroid;
	computeCentroid(data, _rot, _centroid, _originInWorld, extent);

	// compute the transformation from local to world coordinate system
	// required for voxel traversal
	// W_T_L (Local => World)

	trans.setIdentity();
	trans.rotate(Eigen::Quaternionf( _rot ).inverse());
	trans.translate(1.0*_originInWorld);

	// origin of trans is the lowest point
	// now we have to transform by transforming to another corner by reversing the Y-direction

	Eigen::Matrix4f _translate;
	Eigen::Matrix4f _mirror;

	_translate << 1.0, 0.0, 0.0, 0.0,		0.0, 1.0, 0.0, extent.dimY,	0.0, 0.0, 1.0, 0.0,		0.0, 0.0, 0.0, 1.0;
	_mirror << 1.0, 0.0, 0.0, 0.0,		0.0, -1.0, 0.0, 0.0,	0.0, 0.0, 1.0, 0.0,		0.0, 0.0, 0.0, 1.0;
	

	// bug inside :(
	trans =   trans.matrix() * _translate * _mirror;
	delete _costfunction;
}

void GeometricTools::computeManualBoundingBox(Eigen::Matrix3f &_rot, std::vector<Eigen::Vector3f> data, BBExtent &extent, Eigen::Transform3f &trans)
{
	OBB_Computation* _costfunction = new OBB_Computation(data);

	

	//std::cout << "Params: " << params[2] << " / " << params[1] << " / " << params[0] << std::endl;

	Eigen::Vector3f _originInWorld;
	Eigen::Vector3f _centroid;
	computeCentroid(data, _rot.inverse(), _centroid, _originInWorld, extent);

	// compute the transformation from local to world coordinate system
	// required for voxel traversal
	// W_T_L (Local => World)

	trans.setIdentity();
	trans.rotate(Eigen::Quaternionf( _rot ));
	trans.translate(1.0*_originInWorld);

	// origin of trans is the lowest point
	// now we have to transform by transforming to another corner by reversing the Y-direction

	Eigen::Matrix4f _translate;
	Eigen::Matrix4f _mirror;

	_translate << 1.0, 0.0, 0.0, 0.0,		0.0, 1.0, 0.0, extent.dimY,	0.0, 0.0, 1.0, 0.0,		0.0, 0.0, 0.0, 1.0;
	_mirror << 1.0, 0.0, 0.0, 0.0,		0.0, -1.0, 0.0, 0.0,	0.0, 0.0, 1.0, 0.0,		0.0, 0.0, 0.0, 1.0;
	

	// bug inside :(
	trans =   trans.matrix() * _translate * _mirror;
	delete _costfunction;
}


void GeometricTools::getEulerAngles(const Eigen::Quaternionf &quat, float &angleX, float &angleY, float &angleZ, bool radians)
{

 // conversion to Euler angles

 float _q0 = quat.w();
 float _q1 = quat.x();
 float _q2 = quat.y();
 float _q3 = quat.z();

 angleX = atan2(2.0*(_q0*_q1+_q2*_q3),(1.0-2.0*(_q1*_q1+_q2*_q2)));
 angleY = asin(2.0*(_q0*_q2-_q3*_q1));
 angleZ = atan2(2.0*(_q0*_q3+_q1*_q2),(1.0-2.0*(_q2*_q2+_q3*_q3)));

 if ( !radians )
 {
	 angleX*=180.0/PI;
	 angleY*=180.0/PI;
	 angleZ*=180.0/PI;
 }
}


Line::Line(Eigen::Vector3f &startPoint, Eigen::Vector3f &endPoint)
{
	m_startPoint = startPoint;
	m_direction = endPoint - startPoint;
	m_endPoint = endPoint;
	// precompute the dot product, which will be needed all the time later on

	m_directionDotProduct = m_direction.dot(m_direction);

	m_directionNormalized = m_direction.normalized();
}


void Line::distanceToPoint(const Eigen::Vector3f &point, float &distance, float &param, Eigen::Vector3f &intersectionPoint)
{

	Eigen::Vector3f _w = point - m_startPoint;

	double _c1 = _w.dot(m_direction);
	double _c2 = m_directionDotProduct;
	if ( _c1 <= 0 )
	{
		distance=  _w.norm();
		param = 0.0;
		intersectionPoint = m_startPoint;
	}
	else if ( _c2 <= _c1 )
	{
		distance = (point - m_endPoint).norm();
		param = 1.0;
		intersectionPoint = m_endPoint;
	}
	else
	{
		param = _c1 / _c2;
		intersectionPoint = m_startPoint + param * m_direction;
		distance = (intersectionPoint - point).norm();
	}
}


void Line::distanceToPointMAX(const Eigen::Vector3f &point, float &distance, float &param, Eigen::Vector3f &intersectionPoint)
{

	Eigen::Vector3f _w = point - m_startPoint;

	double _c1 = _w.dot(m_direction);
	double _c2 = m_directionDotProduct;
	if ( _c1 <= 0 )
	{
		
		distance=  _w.cwiseAbs().maxCoeff();
		param = 0.0;
		intersectionPoint = m_startPoint;
	}
	else if ( _c2 <= _c1 )
	{
		distance  = (point - m_endPoint).cwiseAbs().maxCoeff();
		param = 1.0;
		intersectionPoint = m_endPoint;
	}
	else
	{
		param = _c1 / _c2;
		intersectionPoint = m_startPoint + param * m_direction;
		distance = (intersectionPoint - point).cwiseAbs().maxCoeff();
	}
}


bool OrientedBoundingBox::lineSegmentIntersection(Line *line) //const Eigen::Vector3f &point1, const Eigen::Vector3f &point2)
{
	float EPSILON = 1.0E-1;
	Eigen::Vector3f _localPoint1 = m_world2local * line->getStartPoint(); //m_world2local * point1;
	Eigen::Vector3f _localPoint2 = m_world2local * line->getEndPoint(); //m_world2local * point2;


	//Eigen::Vector3f c = m_centroid;//(b.min + b.max) * 0.5f; // Box center-point
	Eigen::Vector3f e = m_halfLengths;
	Eigen::Vector3f m = (_localPoint1 + _localPoint2) * 0.5f; // Segment midpoint
	Eigen::Vector3f d = _localPoint2 - m; // Segment halflength vector
	m = m - m_centroid; // Translate box and segment to origin
	// Try world coordinate axes as separating axes
	float adx = abs(d.x());
	if (abs(m.x()) > e.x() + adx) return false;
	float ady = abs(d.y());
	if (abs(m.y()) > e.y() + ady) return false;
	float adz = abs(d.z());
	if (abs(m.z()) > e.z() + adz) return false;
	// Add in an epsilon term to counteract arithmetic errors when segment is
	// (near) parallel to a coordinate axis 
	adx += EPSILON; ady += EPSILON; adz += EPSILON;

	// Try cross products of segment direction vector with coordinate axes
	if (abs(m.y() * d.z() - m.z() * d.y()) > e.y() * adz + e.z() * ady) return false;
	if (abs(m.z() * d.x() - m.x() * d.z()) > e.x() * adz + e.z() * adx) return false;
	if (abs(m.x() * d.y() - m.y() * d.x()) > e.x() * ady + e.y() * adx) return false;
	// No separating axis found; segment must be overlapping AABB
	return true;
}