#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>

#include "compounding.h"



#include <windows.h>
#include <process.h>

#include <queue>
#include <iostream>
#include <sstream>
#include <string>
#include <hash_map>
#include <Winbase.h>
#include <time.h>
#include <algorithm>
#include <stdio.h>
#include <math.h>
#include <omp.h>

#include <io.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <share.h>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include "GraphicTools.h"
#include "Reconstruction.h"
#include "RandomNumbers.h"
#include <omp.h>
#include "RFCompounding.h"
#include "XMLTools.h"
#include "IntensityCompounding.h"

 namespace matlab {
#include <matrix.h>
	#include "mex.h"
	 #include "mat.h"
}

#define DEBUG_TXT_OUTPUT 0
#define DEBUG_MHD_OUTPUT 0



// CAMP Image depenencies
//#include "Common/CAMPImageIO.h"
//#include "Common/ImageBase.h"
//#include "Common/DataTypes.h"



struct WeightedMedianData
{
	float weight;
	float intensity;
};

bool WeightedMedianDataComparator (WeightedMedianData i,WeightedMedianData j) { return (i.intensity<j.intensity); }


double round(double r) {
    return (r > 0.0) ? floor(r + 0.5) : ceil(r - 0.5);
	
}

/*
float round(float r) {
    return (r > 0.0f) ? floor(r + 0.5f) : ceil(r - 0.5f);
}

*/


float round(float number, int digits)
{
	float v[] = { 1, 10, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8 };  // mgl. verlängern
	return floor(number * v[digits] + 0.5) / v[digits];
}


/*#define _CRTDBG_MAP_ALLOC

#include <stdlib.h>

#include <crtdbg.h>
*/

// Note: Within the CAMP_COMMON <vector> has to be replaced with <Eigen/StdVector> solve the unresolved external issue

USING_PART_OF_NAMESPACE_EIGEN

	/**
	\brief A class for precomputing or bining the distances according the radius. Allows to find the appropriate radius for backward warping given a distance.
	It precomputes the maximum distances for each radius (up to a certain maximum) and then determines the appropriate radius by maximum comparision (to avoid operations which are potentially more computationally expensive)
	**/
class RadiusSelect
{
private:
	float *m_distances;
	int m_maxIndex;
public:
	/**
	\brief Constructor for fast radius determination class
	\param threshold Maximum distance for a point to be in range for backward warping.
	\param scale Contains the pixel to space (mm) conversion scale.
	\param maxIndex The maximum number of radii that is possible for backward warping.
	**/
	RadiusSelect(float  threshold, float scale, int maxIndex);
	~RadiusSelect();
	/**
	\brief Returns the backward warping radius, which contains pixel values which have a distance below the specified threshold value.
	\param distance Distance of a point
	\return Radius of the circle containing points of interest.
	**/
	inline int findRadius(float distance);
};

RadiusSelect::~RadiusSelect()
{
	delete[] m_distances;
}

RadiusSelect::RadiusSelect(float threshold, float scale, int maxIndex) :
m_distances(0)
{
	// m_distances is going to contain the maximum orthogonal distance for its indices
	// e.g. m_distance[1] contains the maximum distance for pixel radius 1
	if ( m_distances) delete[] m_distances;
	m_distances = new float[maxIndex];

	m_maxIndex = maxIndex;

	// now compute the maximum distances
	for (int i=1;i<=maxIndex;i++)
	{
		m_distances[i-1]=sqrt(pow(static_cast<float>(i)*scale,2.0f) - pow(threshold,2.0f));
	}
}


inline int RadiusSelect::findRadius(float distance)
{
	// look into every array element for the distance,
	// if it exceeds the specified one, then we have found the radius
	for(int i=0;i<m_maxIndex;i++)
	{
		if ( m_distances[i] > distance)
			return i;
	}
	return m_maxIndex;
}

class USFrame
{
public:
	// video image
	int videoFrame;
	int vectorIndex;
	//CAMP::Matrix4<float>/*<float,4,4,Eigen::DontAlign>*/ image2World;
	//CAMP::Matrix4<float>/*<float,4,4,Eigen::DontAlign>*/ world2Image;
	Eigen::Transform3f image2World;
	Eigen::Transform3f world2Image;
	Eigen::Transform3f local2World; // already in mm

	//Eigen::Transform3f image2WorldEigen;
	//Eigen::Transform3f world2ImageEigen;
	Plane *plane;
	float distance;
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};
/*
class USFrameComparator
{
bool reverse;

public:
USFrameComparator(const bool& revparam=false)	{ reverse=revparam; }
bool operator() (USFrame* lhs, USFrame* rhs) const
{
if (reverse) 
return (lhs->distance < rhs->distance);
else 
return (lhs->distance > rhs->distance);
}
};*/


/**
\brief This class describes how a voxel volume is to be traversed. By consecquitively calling nextPosition, new position data is generated.
*/
class TraversalScheme
{
private:
	int m_columnCounter, m_rowCounter;
	bool m_forward, m_downward;
	int m_X, m_Y, m_Z;
	int m_voxelsX, m_voxelsY, m_voxelsZ;
	std::vector<Eigen::Vector2i> m_traversalList[6];
	int m_offsetX;
	int m_offsetY;
	int m_offsetZ;
public:
	/**
	\brief Constructor for the traversal scheme
	\param voxelsX Number of voxels in dimension X
	\param voxelsY Number of voxels in dimension Y
	\param voxelsZ Number of voxels in dimenion Z
	\param offsetX Voxel X position in 3D space where to start position (default = 0)
	\param offsetY Voxel Y position in 3D space where to start position (default = 0)
	\param offsetZ Voxel Z poistion in 3D space where to start position (default = 0)
	**/
	TraversalScheme(const int &voxelsX, const int &voxelsY, const int &voxelsZ, int offsetX = 0, int offsetY = 0, int offsetZ = 0);

	/**
	\brief Prior to use of traversal scheme, the start position has to be reset.
	**/
	void init();

	/** 
	\brief Generates the next 3f position in the voxel volume according to the traversal rules.
	\param posX X coordinate of the next position
	\param posY Y coordinate of the next position
	\param posZ Z coordinate of the next position
	\return True if the returned position is valid, false otherwise (e.g. traversal is over)
	**/
	inline bool nextPosition(int &posX, int &posY, int &posZ);
};




TraversalScheme::TraversalScheme(const int &voxelsX, const int &voxelsY, const int &voxelsZ, int offsetX, int offsetY, int offsetZ):
m_forward(true),
	m_downward(true),
	m_X(0),
	m_Y(0),
	m_Z(0),
	m_columnCounter(0),
	m_rowCounter(0)
{
	m_voxelsX = voxelsX;
	m_voxelsY = voxelsY;
	m_voxelsZ = voxelsZ;

	m_offsetX = offsetX;
	m_offsetY = offsetY;
	m_offsetZ = offsetZ;
}

void TraversalScheme::init() 
{
	m_rowCounter = 0;
	m_columnCounter = 0;
	m_downward = true;
	m_forward = true;

	m_X = m_offsetX;
	m_Y = m_offsetY;
	m_Z = m_offsetZ;
}

inline bool TraversalScheme::nextPosition(int &posX, int &posY, int &posZ)
{
	// check if traversal is still in range of the volume, otherwise return false
	if ( m_X >= m_voxelsX || m_Y >= m_voxelsY || m_Z >= m_voxelsZ )
		return false;
	posX = m_X;
	posY = m_Y;
	posZ = m_Z;


	// check if you have to reverse the X direction (after following voxelsX steps along X direction)
	if ( ++m_columnCounter == m_voxelsX )
	{


		m_columnCounter = 0;

		// change direction
		m_forward = !m_forward;

		// check if you have to reverse the Y direction (after following voxelsY steps along Y direction)
		if ( ++m_rowCounter == m_voxelsY )
		{
			m_Z++;
			m_rowCounter = 0;

			// change direction
			m_downward = !m_downward;
		}
		else
		{
			// increase or decrease Y coordinate, depending on the mode: forward or backward
			if ( m_downward)
				m_Y++;
			else
				m_Y--;
		}
	}
	else
	{
		// inrease or decrease X coordinate, depending on the mode:  forward or backward
		if ( m_forward )
			m_X++;
		else
			m_X--;
	}


	// traversal has not exceeded the volume limits
	return true;
}

/**
\brief Element containing information about the ultrasound slices, which are inserted into the rotation queue.
**/
template <typename T>
struct QueueElement
{
	// each vector represents a bin, containing slices of the same distance in the rotation queue
	std::vector<T*>* m_sliceVector;
	// pointer to the next element in the rotation queue
	QueueElement* next;
};

/**
\brief The RotationQueue class defines a queue, in which ultrasound slices are sorted according to their distance to a voxel.
**/
template <typename T>
class RotationQueue
{
private:
	QueueElement<T> *m_firstElement,*m_lastElement;
	int m_elements;
	float m_step;
	int m_startIndex;
	int m_averageElements;
	QueueElement<T> **m_queue;
	std::vector<std::vector<T*>*> m_recyclingContainer; 
public:
	/**
	\brief Constructor for the rotation queue.
	\param step The step size of the rotation queue. Extent divided by step equals the number of bins of the rotation queue.
	\param extent The spatial extent of the volume. Extent divided by step equals the number of bins of the rotation queue.
	\param slices The number of slices to be inserted into the rotation queue. Might be used for statistical reasons (can be removed).
	**/
	RotationQueue(float step, float extent, int slices);
	~RotationQueue();

	/**
	\brief Pushes a pointer to an ultrasound slice to the rotation queue into the appropriate bin according to the distance.
	\param frame Pointer to the ultrasound frame.
	\param distance Distance of the ultrasound frame to the current voxel.
	**/
	inline void push(T* frame, float distance);

	/**
	\brief Pops the minimum bin vector (containing the ultrasound slices)
	\return Vector of ultrasound frames contained in the minimum distance bin.
	**/
	inline std::vector<T*>* pop();

	void recycle(std::vector<T*> *obj);
};


template  <typename T> void RotationQueue<T>::recycle(std::vector<T*> *obj)
{
	m_recyclingContainer.push_back(obj);
}

template  <typename T> inline std::vector<T*>* RotationQueue<T>::pop()
{
	// grab the pointer to the first bin in the rotation queue
	std::vector<T*> *_frameVec = m_firstElement->m_sliceVector;

	QueueElement<T> *_temp = m_firstElement;
	m_firstElement = m_firstElement->next;

	// now we have to generate a new bin, as we return one and the hole needs to be filled
	std::vector<T*> *_tempVec;
	if ( m_recyclingContainer.size() > 0 )
	{
		_tempVec = m_recyclingContainer.back();
		m_recyclingContainer.pop_back();
	}
	else
		_tempVec = new std::vector<T*>;
	// reserve memory for the slices to be re-inserted
	// reserve at least the average slice number, or the last number of slices
	//_tempVec->reserve(std::max((int)_frameVec->size(), m_averageElements));
	_temp->m_sliceVector = _tempVec;

	// the new bin will be the last one and doesn't have a successor yet, so its next is NULL
	_temp->next = NULL;

	// the former last element will now the last but one element, so its next has to be adapted
	m_lastElement->next = _temp;
	m_lastElement = _temp;

	// as we are working on a fixed size array, the index of the first element shifts
	m_startIndex++;

	return _frameVec;
}

template  <typename T> inline void RotationQueue<T>::push(T* frame, float distance)
{
	//compute the appropriate index of the ultrasound frame according to its distance
	int index = std::max((int)floor((distance/m_step)),0);

	// m_queue contains the pointers to the rotation queue elements
	// every pop operations shifts the startIndex one position further, which needs to be
	// accounted for when adding new elements, to find the right bin
	// e.g. when only pushing without poppiong m_startIndex equals 0
	// however, when poping once, the first element in the rotation queue will not be on position 0
	// in m_queue anymore, but one 1.
	m_queue[(index+m_startIndex) % m_elements]->m_sliceVector->push_back(frame);
}

template  <typename T> RotationQueue<T>::~RotationQueue()
{
	QueueElement<T> *_current = m_firstElement;
	while( _current->next != NULL)
	{
		QueueElement<T> *_temp = _current;
		_current = _current->next;

		if ( _temp->m_sliceVector )
			delete _temp->m_sliceVector;
		if ( _temp )
			delete _temp;
	}
	if (m_queue)
		delete[] m_queue;

	while(m_recyclingContainer.size() > 0 )
	{

		std::vector<T*> *_tempVec = m_recyclingContainer.back();
		m_recyclingContainer.pop_back();
		delete _tempVec;
	}

	delete _current;
}

template  <typename T> RotationQueue<T>::RotationQueue(float step, float extent, int slices)
{ 
	m_startIndex = 0;
	m_queue = 0;

	// determine the number of bin the rotation queue is going to have
	m_elements = (int)ceil(extent/step);
	// average elements per slot with a safety margin (25% added)
	m_averageElements = (int)ceil(((float)m_elements / (float)slices)*1.25);

	// if constructor is called several times, remove the old queue
	if ( m_queue) delete[] m_queue;

	// this array has to point to vectors represening the bins
	// however, the index of the first element shifts with each pop operation
	m_queue = new QueueElement<T>*[m_elements];

	m_step = step;

	QueueElement<T> *last = NULL;
	// generate emptry vectors for each bin
	for (int i=0;i<m_elements;i++)
	{
		QueueElement<T> *_elem = new QueueElement<T>();
		std::vector<T*> *_vec = new std::vector<T*>;
		_elem->next = last;
		_elem->m_sliceVector = _vec;

		last = _elem;

		// set the pointers for the first and last element
		if ( i == 0 )
			m_lastElement = _elem;
		else if ( i == m_elements-1 )
			m_firstElement = _elem;
	}
	// now set the next pointers to neighboring elements properly
	QueueElement<T> *_temp = m_firstElement;
	for (int i=0;i<m_elements;i++)
	{
		m_queue[i] = _temp;
		_temp = _temp->next;
	}
}



void CompoundVolumeThread::run()
{

	// determine the maximum dimension of the voxel spacing
	float _maxDim = m_voxelSize;
	float _voxelDiameter = sqrt(m_physicalX*m_physicalX+m_physicalY*m_physicalY+m_physicalZ*m_physicalZ)/2.0f;

	// initialize the rotation queue
	RotationQueue<USFrame> _rotQueue(_maxDim,m_extent.dimX+m_extent.dimY+m_extent.dimZ, m_imageCounter);

	Eigen::Matrix4f _scalingMatrix;
	_scalingMatrix <<	m_scaleX,0,0,0, 
		0,m_scaleY,0,0, 
		0,0,1,0, 
		0,0,0,1;

	// now we set up a data structure containing information about the position of the ultrasound planes in space
	for (int i=0; i<m_matrixDataTracking->size();i++)
	{
		USFrame* _tempFrame = new USFrame();
		// the tracking data specifies where the plane lies in space, here the origin is assumed in the upper left corner of the image
		Plane *_plane = new Plane((*m_matrixDataTracking)[i]);
		_tempFrame->distance = (float)0.0;
		// _matrixDataTracking.at(x) transforms from unit [mm] on the image to world [mm]
		// but we typically need it for pixels, thus we multiply it with the scaling factors

		//TJK_27_07_09: m_dataTools.eigen2CAMP(m_matrixDataTracking[i] * _scalingMatrix, _tempFrame->image2World);
		_tempFrame->image2World = ((*m_matrixDataTracking)[i] * _scalingMatrix);

		// compute all the relevant information describing the ultrasound plane
		//TJK_27_07_09: _tempFrame->world2Image = _tempFrame->image2World.getInverse();
		_tempFrame->world2Image = _tempFrame->image2World.inverse();
		//_tempFrame->image2WorldEigen = _matrixDataTracking.at(i) * _scalingMatrix;
		//_tempFrame->world2ImageEigen = _tempFrame->image2WorldEigen.inverse();
		_tempFrame->vectorIndex = (int)floor((float)i/NUMBER_OF_SLICES); 
		_tempFrame->videoFrame = i % (int)(NUMBER_OF_SLICES);
		_tempFrame->plane = _plane;

		// set the ultrasound frame at initial distance 0.0
		_rotQueue.push(_tempFrame,0.0);
	}

	int test = m_matrixDataTracking->size();

	RadiusSelect _radiusSelect(m_voxelSize, std::max(m_scaleX, m_scaleY), 50);


	// choose the distance threshold, that is the maximum distance a US slice can have
	// in order to be considered for contributing intensity to a voxel

	float _distanceScalar = 0.0f;
	switch(m_maxDistanceScalar)
	{
	case 0: { _distanceScalar = DISTANCE_SCALAR; break;}
	case 1: { _distanceScalar = 1.0; break;}
	case 2: { _distanceScalar = 1.1; break;}
	case 3: { _distanceScalar = 1.2; break;}
	case 4: { _distanceScalar = 1.3; break;}
	case 5: { _distanceScalar = 1.4; break;}
	case 6: { _distanceScalar = 1.5; break;}
	};

	float _threshold = 0.0f;
	switch(m_maxDistanceUnit)
	{
	case 0: { _threshold = _maxDim; break;}
	case 1: { _threshold = _voxelDiameter; break;}
	};

	// now we compute the actual threshold value from the UI selection
	_threshold = _threshold * _distanceScalar;

	// precompute the squared distance to speed-up the computation
	//float _thresholdSQR = pow(_threshold,2.0f);
	stdext::hash_map<int, std::vector<Segment>*> _hashMap;
	stdext::hash_map<int, std::vector<Segment>*>::iterator _hashMapIterator;



	TraversalScheme _traversalScheme(m_voxelsX, m_voxelsY, m_voxelsZ, m_offsetX, m_offsetY, m_offsetZ);
	_traversalScheme.init();

	int _x,_y,_z;
	int _traversalCount = 0;

	// Gaussian backward
	// sigma is only needed as squared, so precalculate to save computation time
	float _sigmaSQR = pow(m_gaussianSigma,2.0f);

	// Inverse distance weighted smoothness factor
	float _inverseWeightedSmoothness = (m_inverseDistanceFactor);

	Eigen::Vector3f _worldIntersection;


	GraphicTools _graphicTools;

	// traverse through the volume

	//ui.progressBar->setRange(0,100);
	// as to avoid too frequen updates of the GUI we do only updates when values change
	//int _maxVoxels =  m_voxelsX*m_voxelsY*m_voxelsZ;
	//int _1percent = (int)_maxVoxels/100.0;
	int _percentCounter = 0;
	//int _percent = 0;
	//ui.progressBar->setValue(_percent);
	//qApp->processEvents();
	while(_traversalScheme.nextPosition(_x,_y,_z))
	{
		_percentCounter++;
		_traversalCount++;

		//  get the memory addess where the intensity to be set in the MHD data is located
		COMPOUNDING_DATA_TYPE *_outputImageIntensityPtr = static_cast<COMPOUNDING_DATA_TYPE*>(m_targetImage->GetScalarPointer(_x, _y, _z));

		if ( _percentCounter == m_1percent)
		{
			emit percentEvent();
			_percentCounter = 0;
		}
		//	_percent++;
		//	_percentCounter = 0;
		//	ui.progressBar->setValue(_percent);
		//	char buffer[50];
		//	sprintf(buffer,"TASScomp - %d%%",_percent);
		//	this->setWindowTitle(QString(buffer));
		//	qApp->processEvents();
		//}


		// compute the world coordinates from the voxel coordinates
		const Eigen::Vector3f _local(static_cast<float>(_x)*m_physicalX+m_physicalX/2.0,static_cast<float>(_y)*m_physicalY+m_physicalY/2.0,static_cast<float>(_z)*m_physicalZ+m_physicalZ/2.0);

		const Eigen::Vector3f _world = (m_boundingBoxTransformation*_local);



		// initalize the values which determine the color of the voxel (gaussian backward warping, inverse distance weighting)
		float _intensityNominator = 0.0;
		float _intensityDenominator = 0.0;

		// needed for sigmoid
		std::vector<float> _confidence;
		std::vector<float> _confidenceIntensity;


ADAPTIVE_REPEAT:

		int _numberOfVoxelContributors = 0;

		// pop the minimum distance bin vector from the rotation queue
		std::vector<USFrame*>* _usFrameQueue = _rotQueue.pop();
		int _oldFrameSize = _usFrameQueue->size();

		// this vector is needed for determining the mean color
		std::vector<float> _intensityVector;

		std::vector<WeightedMedianData> _weightedMedianVector;
		float _weightedIntensitySum = 0.0f;

		// this variables are needed for nearest neighbor
		float _minDistance = std::numeric_limits<float>::max();
		float _minDistanceIntensity;

		// for each slice within the vector returned form the queue do a backward warping
		while( _usFrameQueue->size() > 0 )
		{

			USFrame *_frame = _usFrameQueue->back();

			_usFrameQueue->pop_back();

			Plane *_plane = _frame->plane;


			// determine the orthogonal distance from the voxel to the ultrasound plane
			const float _distance = _plane->getOrthogonalDistance (_world, _worldIntersection );

			// determine the point of minimum distance to the voxel on the plane (in pixel coordinates)
			const Eigen::Vector3f _localIntersection = (_frame->world2Image*_worldIntersection);

			// point is relevant for the voxel, its intensity has to be considered / it has to be within the image
			// NOTE: the _roiMin/Max is because of the use of the cropped image area within the video frame
			if ( _distance < _threshold && _localIntersection[0] >= m_roiMinX && _localIntersection[0] <= m_roiMaxX && _localIntersection[1] >= m_roiMinY && _localIntersection[1] <= m_roiMaxY)
			{

				// update the distance
				_frame->distance = _distance;

				// and push it into the temporary queue, as the point is not to be considered again for the same voxel
				_rotQueue.push(_frame, _distance);




				// grab the video image

				// very slow code
				//HBITMAP bmp;
				//m_videoSegmentationTools.getFrame(_frame.videoFrame, bmp);
				//float _tempIntensity = m_videoSegmentationTools.getIntensity(bmp, (int)_localIntersection(0), (int)_localIntersection(1));




				// now compute the radius in the plane that is within distance specified by the threshold
				// that's simply pythagoras 
				//int _radius = static_cast<int>(sqrt ( _distance * _distance + _thresholdSQR ) / _scaleX);
				// instead of pythagoras a faster version is to take the precomputed distance values

				int _radius = 0; 
				switch(m_compoundingMode )
				{
				case GAUSSIAN: 
					{
						// determine the pixel area contributing to backward-warping
						_radius = _radiusSelect.findRadius(_distance);
						break;
					}
				case INVERSE_DISTANCE:
					{
						// determine the pixel area contributing to backward-warping
						_radius = _radiusSelect.findRadius(_distance);
						break;
					}
				case NEAREST_NEIGHBOR:
					{
						// assume _radius = 0, so on the slice, only the pixel with min.
						// orthogonal distance is selected
						break;
					}
				case MEDIAN:
					{
						// assume _radius = 0, so on the slice, only the pixel with min.
						// orthogonal distance is selected
						break;
					}
				case WEIGHTED_MEDIAN:
					{
						// assume _radius = 0, so on the slice, only the pixel with min.
						// orthogonal distance is selected
						break;
					}
				};

				std::vector<Segment> *_SegmentVector;

				// check the hash map if for the radius there exist already Segment information
				// the Segment contains the discretization of the circle for a given radius
				// e.g. K pixels above the center, move from -M to +M (relative position)
				// it is like the active edge table in the rasterizer
				_hashMapIterator = _hashMap.find(_radius);


				// in case there is no Segment information available, we have to generate it
				if ( _hashMapIterator == _hashMap.end() )
				{
					_SegmentVector = new std::vector<Segment>;

					// compute the relative position Segments of the circle
					// rasterization table for a circle of given radius (relative coordinates => 0,0 as first parameters)
					_graphicTools.bresenhamCircle(0,0, _radius, _SegmentVector);

					// place the information in the hashmap
					_hashMap[_radius]=_SegmentVector;
				} // Segment information is available, grab it
				else
				{
					// discretization was computed previously, so we simply have to acquire the pointer
					_SegmentVector = _hashMapIterator->second;
				}

				// now we look up follow the information contained in the vector (visit all the pixels, which are in range)
				// cache the size of the vector to save computations
				std::vector<Segment>::iterator _itBegin = _SegmentVector->begin();
				std::vector<Segment>::iterator _itEnd = _SegmentVector->end();

				for(std::vector<Segment>::iterator it = _itBegin; it != _itEnd;++it )
				{
					Segment _Segment = *it;

					// as we have only relative positional information from the mid-point/bresenham
					// we have to transform them in absolute positions on the image
					// we have to adapt the coordinates, as we have cropped the image in memory (that's why using _roiMin*)
					_Segment.Y			+= static_cast<int>(_localIntersection[1])- m_roiMinY;
					_Segment.beginX	+= static_cast<int>(_localIntersection[0])- m_roiMinX;
					_Segment.endX		+= static_cast<int>(_localIntersection[0])- m_roiMinX;
					// only consider points which are within the image
					if ( _Segment.Y < m_roiHeight && _Segment.Y >= 0 )
					{

						int _beginX, _endX;
						bool _draw = true;

						// if the Segment is outside the window, forget about it
						if ( _Segment.beginX >= m_roiWidth )
							_draw = false;
						if ( _Segment.endX < 0 )
							_draw = false;

						// if Segment-range is within the area of interest backward-warp it
						if ( _draw && ((_Segment.beginX	 >= 0 && _Segment.beginX <= m_roiWidth-1) ||( _Segment.endX >= 0 && _Segment.endX <= m_roiWidth-1)) )
						{
							// clip to the boundary, define the rasterization limits (the bounding box of the clipping mask)
							_beginX = std::max( std::min( _Segment.beginX, m_roiWidth-1), 0);
							_endX = std::min(std::max( _Segment.endX, 0 ), m_roiWidth-1 );


							IplImage *_img = (*m_imageVector)[_frame->vectorIndex][(_frame->videoFrame)];
#ifdef IMAGE_FILTER
							IplImage *_filter = (*m_filterVector)[_frame->vectorIndex][(_frame->videoFrame)];
#endif
							// get the pointers to the ultrasound image and the cropping mask
							unsigned char *_imgPointer= &(((unsigned char*)(_img->imageData + _Segment.Y*_img->widthStep))[_beginX]);
							unsigned char *_croppingMaskPointer = &(((unsigned char*)(m_croppingMask->imageData + _Segment.Y*m_croppingMask->widthStep))[_beginX]);
							// now pixel by pixel backward warping (given the radius)
							#ifdef IMAGE_FILTER
								unsigned char *_filterPointer= &(((unsigned char*)(_filter->imageData + _Segment.Y*_filter->widthStep))[_beginX]);
							#endif
							for (int x_pos=_beginX;x_pos <= _endX; x_pos++)
							{

								// now apply the cropping
								// get the intensity of the cropping mask
								const float _maskIntensity = (float)*_croppingMaskPointer++;

								// get the intensity of the ultrasound image
								const float _tempIntensity = (float)*_imgPointer++;	

								#ifdef IMAGE_FILTER
								const float _filterIntensity = (float)*_filterPointer++;	
								#endif
								
								// if the intensity equals white, then we are in non-cropped area
								if ( _maskIntensity > 0.0)
								{

									//if ( _filterIntensity > 0.0 && _filterIntensity < 255.0 )
								//		int debug = 5;
									// here we go: we have a pixel, whose itensity is contributing to the voxel intensity!
									_numberOfVoxelContributors ++;
									// get the intensity of the ultrasound image
									//const float _tempIntensity = (float)*_imgPointer++;			

									//cvSaveImage("c:\\Resultant.jpg",m_imageVector.at(_frame->videoFrame));
									// compute the distance for the point => _distance is just the orthogonal distance to the plane


									switch(m_compoundingMode )
									{
									case SIGMOID:
										{
											_confidenceIntensity.push_back(_tempIntensity);
											#ifdef IMAGE_FILTER
											_confidence.push_back(_filterIntensity);
											#endif
											break;
										}
									case GAUSSIAN: 
										{
											// Gaussian backward warping

											// determine the position of the pixel in world coordinates (mm)
											const Eigen::Vector3f _temp(x_pos+m_roiMinX, _Segment.Y+m_roiMinY, 0.0);
											const Eigen::Vector3f _worldPoint = (_frame->image2World * _temp);

											// determine squared euclidean distance of the point to the pixel in world coordinates
											const float _euclideanDistanceSQR =  pow(static_cast<float>((_world - _worldPoint).norm()),2.0f);

											// determine the interpolation ratio
											const float _tempRatio = exp(-(_euclideanDistanceSQR)/(_sigmaSQR));

											// add the intensity contribution by weighting distance to the nominator and denominator
											// which after considering all relevant pixels will by division supply the gaussian interpolated intensity
											_intensityDenominator += _tempRatio;
											#ifdef IMAGE_FILTER
											
											_intensityNominator += _tempIntensity * _tempRatio * _filterIntensity/255.f;

#else
											_intensityNominator += _tempIntensity * _tempRatio;// * _filterIntensity/255.f;

#endif								
											break;
										}
									case INVERSE_DISTANCE: 
										{
											// Inverse distance backward warping

											// determine the position of the pixel in world coordinates (mm)
											const Eigen::Vector3f _temp(x_pos+m_roiMinX, _Segment.Y+m_roiMinY, 0.0);
											const Eigen::Vector3f _worldPoint = (_frame->image2World * _temp);

											// determine euclidean distance by power of smoothness factor of the point to the pixel in world coordinates
											const float _euclideanDistanceSmoothed =  pow(static_cast<float>((_world - _worldPoint).norm()),2.0f);


											// add the intensity contribution by weighting distance to the nominator and denominator
											// which after considering all relevant pixels will by division (and denominator * N) give the inverse distance measure
											_intensityDenominator += _euclideanDistanceSmoothed;
											_intensityNominator += _tempIntensity * _euclideanDistanceSmoothed; 
											break;
										}
									case NEAREST_NEIGHBOR:
										{
											if ( _distance < _minDistance )
											{
												_minDistance = _distance;
												_minDistanceIntensity = _tempIntensity;
											}
											break;
										}
									case MEDIAN:
										{
											_intensityVector.push_back(_tempIntensity);
											break;
										}
									case WEIGHTED_MEDIAN:
										{
											WeightedMedianData _wmd;

											_wmd.weight = 1.0 - _distance / _threshold;
											_wmd.intensity = _tempIntensity;
											_weightedMedianVector.push_back(_wmd);
											_weightedIntensitySum += _wmd.weight;
											break;
										}
									}
								}
							}
						}
					}

				}

			}
			else // point too far away or not within the region of interest, but distance was updated so we have to update it in the rotation queue
			{
				// update the distance
				_frame->distance = _distance;

				// and push it into the temporary queue, as the point is not to be considered again for the same voxel
				_rotQueue.push(_frame, _distance);
			}

		}
		// instead of delete put empty vectors back to rotation queue, so we save new operations
		_usFrameQueue->reserve(_oldFrameSize);
		_rotQueue.recycle(_usFrameQueue);

		// now set the intensity of the voxel
		// small box for showing the origin and the orientation
		//if ( _x <=25 && _y <=12 && _z <=5 )
		//{
		//	_intensityDenominator = 1.0;
		//	_intensityNominator = 255.0;
		//}



		// calculate the intensity of the voxel given the specified method

		switch(m_compoundingMode )
		{
		case SIGMOID:
			{
				NumericTools _nr;
				float _confidenceVar, _confidenceMean, _confidenceMax, _confidenceMin;
				if (_confidence.size() > 0 ) {
					_nr.variance(_confidence,& _confidenceVar,& _confidenceMean, &_confidenceMax, &_confidenceMin);

				for(int k=0;k<_confidence.size();k++)
					 {
						 // Tassilo's approach
						 float _tmp = exp(-pow(_confidence[k]-_confidenceMax,2.0f)/pow(_confidenceMax,2.0f)); 
						 // Athansios' approach
						/* float _sigmaSQR = 1.0;
						 float _tmp = exp(-pow(_confidence[k]-_confidenceMax,2.0f)/(2.0*_sigmaSQR)); */

						 _intensityNominator += _tmp*_confidenceIntensity[k];
						_intensityDenominator += _tmp;

						 // multiplicative
						// _intensityNominator +=  _confidenceIntensity[k] * (_confidence[k]/255.0f);
					 }
				_intensityDenominator = static_cast<float>(_confidence.size());
				*_outputImageIntensityPtr = static_cast<COMPOUNDING_DATA_TYPE>( _intensityNominator / _intensityDenominator );
				}

			break;
			}
		case GAUSSIAN: 
			{
				// avoid division by zero
				if ( _intensityDenominator > 0.0 && _intensityNominator > 0.0 ) // avoid division by ZERO
				{
					// CAMP Image replaced
					//compounding->setVoxelFloat( _intensityNominator / _intensityDenominator,_x+_y*m_voxelsX+_z*m_voxelsX*m_voxelsY);
					*_outputImageIntensityPtr = static_cast<COMPOUNDING_DATA_TYPE>( _intensityNominator / _intensityDenominator );
				}
				break;
			}

		case INVERSE_DISTANCE: 
			{
				// avoid division by zero
				if ( _intensityDenominator > 0.0 && _intensityNominator > 0.0 ) // avoid division by ZERO
				{
					// CAMP Image replaced
					//compounding->setVoxelFloat( _intensityNominator / ((float)(_numberOfVoxelContributors)*_intensityDenominator),_x+_y*m_voxelsX+_z*m_voxelsX*m_voxelsY);
					*_outputImageIntensityPtr = static_cast<COMPOUNDING_DATA_TYPE>( _intensityNominator / ((float)(_numberOfVoxelContributors)*_intensityDenominator) );
				}
				break;
			}

		case NEAREST_NEIGHBOR:
			{
				// check if there was a slice within range
				if ( _minDistance != std::numeric_limits<float>::max() )
				{
					// CAMP Image replaced
					//compounding->setVoxelFloat( _minDistanceIntensity,_x+_y*m_voxelsX+_z*m_voxelsX*m_voxelsY);
					*_outputImageIntensityPtr = static_cast<COMPOUNDING_DATA_TYPE>( _minDistanceIntensity );
					// reset the minimum distance for the next run
					_minDistance = std::numeric_limits<float>::max();
				}
				break;
			}
		case MEDIAN:
			{
				int _vecSize = _intensityVector.size();
				if ( _vecSize> 1 )
				{
					// sort the intensities and take the mean
					sort (_intensityVector.begin(), _intensityVector.end());
					// CAMP Image replaced
					//compounding->setVoxelFloat( _intensityVector[round((float)_intensityVector.size()/2.0,2)],_x+_y*m_voxelsX+_z*m_voxelsX*m_voxelsY);
					*_outputImageIntensityPtr = static_cast<COMPOUNDING_DATA_TYPE>( _intensityVector[round((float)_intensityVector.size()/2.0,2)] );
				}
				else if ( _vecSize == 1) 
				{
					// we only have one intensity...no need to sort, just take it
					// CAMP Image replaced
					//compounding->setVoxelFloat( _intensityVector[0],_x+_y*m_voxelsX+_z*m_voxelsX*m_voxelsY);
					*_outputImageIntensityPtr = static_cast<COMPOUNDING_DATA_TYPE>( _intensityVector[0] );
				}
				else
				{
					// do nothing, leave the voxel black
				}
				// reset the intensity vector for the next run
				_intensityVector.clear();
				break;
			}
		case WEIGHTED_MEDIAN:
			{
				int _vecSize = _weightedMedianVector.size();
				if ( _vecSize > 1 )
				{
					float _tempLast = 0.0f;
					float _half = _weightedIntensitySum / 2.0f;

					for(int i=1;i<_vecSize;i++)
					{
						float _w1 = _weightedMedianVector[i-1].weight+_tempLast;
						float _w2 = _weightedMedianVector[i].weight+_w1;
						_tempLast = _w1;

						if ( _w1 <=  _half && _half <= _w2)
						{
							// CAMP Image replaced
							//compounding->setVoxelFloat( _weightedMedianVector[i].intensity,_x+_y*m_voxelsX+_z*m_voxelsX*m_voxelsY);
							*_outputImageIntensityPtr = static_cast<COMPOUNDING_DATA_TYPE>(  _weightedMedianVector[i].intensity );
							break;
						}
					}
				}
				else if ( _vecSize == 1)
				{
					// CAMP Image replaced
					//compounding->setVoxelFloat( _weightedMedianVector[0].intensity,_x+_y*m_voxelsX+_z*m_voxelsX*m_voxelsY);
					*_outputImageIntensityPtr = static_cast<COMPOUNDING_DATA_TYPE>( _weightedMedianVector[0].intensity );
					break;
				}
				else
				{
					// leave it black
				}
				// reset the intensity vector for the next run
				_weightedMedianVector.clear();
				break;
			}
		}



	} // end of while

	// free the memory

	stdext::hash_map<int, std::vector<Segment>*>::iterator _hashBegin = _hashMap.begin();
	stdext::hash_map<int, std::vector<Segment>*>::iterator _hashEnd = _hashMap.end();

	for(_hashMapIterator = _hashBegin; _hashMapIterator!= _hashEnd; ++_hashMapIterator)
	{
		if ( _hashMapIterator->second )
			delete _hashMapIterator->second;
	}


}

void CompoundVolumeThread::setTransformation(Eigen::Transform3f& boundingBoxTransformation)
{
	m_boundingBoxTransformation = boundingBoxTransformation;
}

void CompoundVolumeThread::setUserConfiguration(int maxDistanceScalar, int maxDistanceUnit, float gaussianSigma, float inverseDistanceFactor, int _1percent)
{
	m_maxDistanceScalar = maxDistanceScalar;
	m_maxDistanceUnit = maxDistanceUnit;
	m_gaussianSigma = gaussianSigma;
	m_inverseDistanceFactor = inverseDistanceFactor;
	m_1percent = _1percent;

}

void CompoundVolumeThread::setROI(int roiMaxX, int roiMaxY, int roiMinX, int roiMinY, int roiHeight, int roiWidth)
{
	m_roiMaxX = roiMaxX;
	m_roiMaxY = roiMaxY;
	m_roiMinX = roiMinX;
	m_roiMinY = roiMinY;
	m_roiHeight = roiHeight;
	m_roiWidth = roiWidth;
}

void CompoundVolumeThread::setDataPointers(std::vector<std::vector<IplImage*>> *imageVector, std::vector<std::vector<IplImage*>> *filterVector, std::vector< Eigen::Matrix4f, Eigen::aligned_allocator<Eigen::Matrix4f> > *matrixDataTracking, std::vector<Eigen::Transform3f, Eigen::aligned_allocator<Eigen::Transform3f> > *trackingData)
{
	m_imageVector = imageVector;
	m_matrixDataTracking = matrixDataTracking;
	m_trackingData = trackingData;
	#ifdef IMAGE_FILTER
	m_filterVector = filterVector;
#endif
}


void CompoundVolumeThread::setAuxiliaryPointers(VideoTools *videoTools, IplImage* croppingMask)
{
	m_videoTools = videoTools;
	m_croppingMask = croppingMask;
}

void CompoundVolumeThread::setCompoundingMode(COMPOUNDING_MODE mode)
{
	m_compoundingMode = mode;
}

void CompoundVolumeThread::setTargetImage(vtkImageData* image)
{
	m_targetImage = image;
}

void CompoundVolumeThread::setDataSize(const int dimX, const int dimY, const int dimZ, const float physicalX, const float physicalY, const float physicalZ, BBExtent extent, int offsetX, int offsetY, int offsetZ, int imageCounter)
{

	m_imageCounter = imageCounter;

	m_voxelsX = dimX;
	m_voxelsY = dimY;
	m_voxelsZ = dimZ;

	m_physicalX = physicalX;
	m_physicalY = physicalY;
	m_physicalZ = physicalZ;

	m_voxelSize = std::min(std::min(physicalX,physicalY),physicalZ);
	m_extent = extent;

	m_offsetX = offsetX;
	m_offsetY = offsetY;
	m_offsetZ = offsetZ;
}

void CompoundVolumeThread::setUSResolution(float scaleX, float scaleY)
{
	m_scaleX = scaleX;
	m_scaleY = scaleY;
}


Compounding::Compounding(QWidget *parent, Qt::WindowFlags flags)
	: QMainWindow(parent, flags),
	m_guiVideoReady(false),
	m_guiTrackingReady(false),
	m_guiDataLoaded(false),
	m_guiRFCalibrationReady(false),
	m_guiRFDataReady(false),
	m_guiRFTrackingReady(false),
	m_logCompress(true), // initial state is true
	m_croppingMask(0),
	m_voxelSize(0.5),
	m_imageCounter(0),
	m_recon(0),
	m_colLUT(0),
	m_percentage(0),
	m_imageForVisualization(0),
	m_masterImage(0),
	m_trackingIndex(0),
	m_subvolumeWriting(NULL),
	m_rfScanlineIndex(0),
	m_settingsIndex(0),
	m_fileIndex(0),
	m_rfImage(0),
	m_rfPreviewImage(0),
	m_xmlTools(0),
	m_threadingFactor(4.0/3.0), // how many more threads than cores should be created (for optimized CPU usage)
	m_rfFrameCounter(0)
{

	m_mhdFileName = "testCompounding.mhd";
	m_rfRFDfilename = "testRFCompounding.xml";
	//hMutex = CreateMutex(NULL, FALSE, NULL);
	/*
	BBExtent _extent;

	_extent.dimX = 100;
	_extent.dimY = 100;
	_extent.dimZ = 100;

	float _rot[]={0.5,0.2,0.7};
	Eigen::Vector3f _trans;
	_trans << 0.0,0.0,0.0;
	Eigen::Matrix3f _temp;
	_temp = Eigen::AngleAxisf((float)_rot[2], Eigen::Vector3f::UnitZ()) * Eigen::AngleAxisf((float)_rot[1], Eigen::Vector3f::UnitY()) * Eigen::AngleAxisf((float)_rot[0], Eigen::Vector3f::UnitX());


	Eigen::Transform3f _test;
	_test.setIdentity();
	_test.rotate(_temp);
	_test.translation() = _trans;

	OrientedBoundingBox *_obb = new OrientedBoundingBox(_test, _extent);

	Eigen::Vector3f _point1;
	_point1 << 5.0, 0.0, 0.0;
	Eigen::Vector3f _point2;
	_point2 << 150.0, 0.0, 0.0;

	bool _intersecting = _obb->lineSegmentIntersection(_point1, _point2);

	delete _obb;*/
	/*

	RFFileStream *_rfs = new RFFileStream();

	//_rfs->openFile("testFloat.xml");
	_rfs->openFile("exampleRF.xml");
	//int16 *data = reinterpret_cast<int16*>(_rfs->getFrameByIndex(1));

	float *data = reinterpret_cast<float*>(_rfs->getFrameByIndex(0));

	RFData *_rfData = new RFData();

	_rfData->setData(data,128,3072);

	float *_val = _rfData->getScanline(0);



	std::cout << "WTF: " << *(_val) << std::endl;
	{
	float *_ptr = _val;
	float _last = -1;
	for (int i=0; i<3072;i++)
	{
	if ( _last != *_ptr)
	{
	std::cout << i << " => " << *_ptr << std::endl;
	_last = *_ptr;
	}
	_ptr++;
	}
	}

	//data = reinterpret_cast<int16*>(_rfs->getNextFrame(0));
	_rfs->closeFile();

	delete _rfs;*/
	/*RFProcessing *_rp = new RFProcessing(2,50, 6000, 19000, 40000,20);

	float *_test;//{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
	_test = new float[100];
	for(int i=0;i<100;i++)
	_test[i]=(float)i;

	float *_ptr = 0;
	RFData *_rf = new RFData();
	_rf->setData(_test,2,50);
	RFData *_rf2 = _rp->bandpassFilter(_rf);
	_ptr = _rp->envelopeDetection(_test);

	float *_p = &_rf2->getData()[0];
	for (int i=0;i<50;i++)
	{
	std::cout << i << ". " << *_p << std::endl;
	_p++;
	}
	delete[] _ptr;
	delete[] _test;
	delete _rp;
	delete _rf2;

	*/
	/*
	int16 _test[10];

	_test[0].value = 100;
	_test[1].value = 200;
	_test[2].value = 300;
	_test[3].value = 400;
	_test[4].value = 500;
	_test[5].value = 600;
	_test[6].value = 700;
	_test[7].value = 800;
	_test[8].value = 900;
	_test[9].value = 1000;

	char _data[20] = {_test[0].byte.low, _test[0].byte.high, _test[1].byte.low, _test[1].byte.high, _test[2].byte.low, _test[2].byte.high,
	_test[3].byte.low, _test[3].byte.high,_test[4].byte.low, _test[4].byte.high,_test[5].byte.low, _test[5].byte.high,_test[6].byte.low, _test[6].byte.high,
	_test[7].byte.low, _test[7].byte.high,_test[8].byte.low, _test[8].byte.high,_test[9].byte.low, _test[9].byte.high};

	_rp->convert(_data, sizeof(char)*2*10);

	delete _rp;*/
	/*int16 _test;
	//_test.value = -32767;
	//_test.byte.low = 255;
	//_test.byte.hi = 128;
	_test.value = 32767;
	std::cout << "low: " << (int)_test.byte.low << std::endl;
	std::cout << "high: " << (int)_test.byte.hi << std::endl;
	std::cout << "value: " << (int)_test.value << std::endl;

	_test.value = -32767;
	std::cout << "low: " << (int)_test.byte.low << std::endl;
	std::cout << "high: " << (int)_test.byte.hi << std::endl;
	std::cout << "value: " << (int)_test.value << std::endl;*/
	// TJK test code: can be removed! 
	/*Eigen::Vector3f _apex(1, 0, 0);
	Eigen::Vector3f _target(9.1, 0, 0);
	Eigen::Vector3f _velocity(5, 2, 3);
	std::cout << "Jacobian: " << jacobianNormalUncertainty(_target, _apex, _velocity) << std::endl;
	*/
	m_xmlConfiguration.loadFile(std::string(COMPOUNDING_CONFIGURATION));

	// load settings from previous run

	// check if we previously loaded some data from a specific directory, if not use the one from the XML file
	if ( m_lastPath.length() == 0 )
	{
		std::string _path;
		m_xmlConfiguration.handleData(false,_path, std::string(TARGET_DIRECTORY_TAG));
		m_lastTargetPath = QString(_path.c_str());
	}

	// check if we previously loaded some data from a specific directory, if not use the one from the XML file
	if ( m_lastCalibrationPath.length() == 0 )
	{
		std::string _path;
		m_xmlConfiguration.handleData(false,_path, std::string(CALIBRATION_DIRECTORY_TAG));
		m_lastCalibrationPath = QString(_path.c_str());
	}

	// check if we previously loaded some data from a specific directory, if not use the one from the XML file
	if ( m_lastPath.length() == 0 )
	{
		std::string _path;
		m_xmlConfiguration.handleData(false,_path, std::string(DATA_DIRECTORY_TAG));
		m_lastPath = QString(_path.c_str());
	}


	ui.setupUi(this);
	ui.voxelSizeSlider->setRange(0,150);

	// VIDEO IMAGE VISUALIZATION
	// working stuff
	glWidget = new GLWidget;
	//ui.gridLayout->addWidget(glWidget);

	// OFFLINE VISUALIZATION (needed for rendering the cropping mask)
	m_offlineVisualization = new SlicesVisualization;


	// SLICES VISUALIZATION
	m_slicesVisualization = new SlicesVisualization;
	ui.gridLayout->addWidget(m_slicesVisualization);
	m_slicesVisualization->setFocusPolicy(Qt::FocusPolicy::StrongFocus);

	// VELOCITY VISUALIZATION
	m_velocityVisualization = new VelocityVisualization;
	ui.gridLayout2->addWidget(m_velocityVisualization);
	m_velocityVisualization->setFocusPolicy(Qt::FocusPolicy::StrongFocus);


	m_videoTools = new VideoTools();
	m_numericTools = new NumericTools();


	// initialize the vectors to store the images
	for (int i=0;i<NUMBER_OF_BLOCKS;i++)
	{
		std::vector<IplImage*> _temp1;
		m_imageVector.push_back(_temp1);

		std::vector<IplImage*> _temp2;
		m_auxImageVector.push_back(_temp2);

#ifdef IMAGE_FILTER
		std::vector<IplImage*> _temp3;
		m_filterVector.push_back(_temp3);
#endif
	}

	// set the default MHD file name
	m_mhdFileName = ui.outputFileName->text();

	// Doppler specific initialization
	// set the color LUT
	// set up the look-up table for the doppler velocity compounding with exemplary color values
	std::vector<int> _width;
	/* _width.push_back(20);
	_width.push_back(39);
	_width.push_back(5);
	_width.push_back(5);
	_width.push_back(39);
	_width.push_back(20);*/
	_width.push_back(55);
	_width.push_back(255);
	_width.push_back(255);
	_width.push_back(255);
	_width.push_back(255);
	_width.push_back(55);

	std::vector<Eigen::Vector3f,Eigen::aligned_allocator<Eigen::Vector3f>> _color;
	_color.push_back(Eigen::Vector3f(0.0,0.0,200.0));
	_color.push_back(Eigen::Vector3f(0.0,0.0,255.0));
	_color.push_back(Eigen::Vector3f(0.0,255.0,255.0));
	_color.push_back(Eigen::Vector3f(0.0,0.0,0.0));
	_color.push_back(Eigen::Vector3f(255.0,255.0,0.0));
	_color.push_back(Eigen::Vector3f(255.0,0.0,0.0));
	_color.push_back(Eigen::Vector3f(200.0,0.0,0.0));
	/*_color.push_back(Eigen::Vector3f(0.0,0.0,255.0));
	_color.push_back(Eigen::Vector3f(0.0,255.0,255.0));
	_color.push_back(Eigen::Vector3f(0.0,255.0,0.0));
	_color.push_back(Eigen::Vector3f(0.0,0.0,0.0));
	_color.push_back(Eigen::Vector3f(255.0,255.0,0.0));
	_color.push_back(Eigen::Vector3f(255.0,0.0,255.0));
	_color.push_back(Eigen::Vector3f(255.0,0.0,0.0));
	*/

	// initialize the look-up table
	// assume maximum velocity 9.88
	m_colLUT = new ColorLUT(_width, _color, 9.882117688026186);


	// in order for this to work (that is showing the documentation), the workspace has to be set to where the source code is
	QString _address;
	_address = "file:///"+QDir::currentPath()+"/Documentation/Compounding.html";
	ui.webView->load(_address);

	ui.textEdit->append(QString(_address));
	std::cout << _address.toLatin1().data() << std::endl;

}


void Compounding::updatePercentage()
{
	m_percentage++;
	std::stringstream _tmp;
	_tmp << "TASScomp - " << m_percentage << "%";
	this->setWindowTitle(QString(_tmp.str().c_str()));
	ui.progressBar->setValue(m_percentage);
	//qApp->processEvents();
	//QApplication::processEvents();

}

void Compounding::compound()
{
	bool _doppler = false;

	if ( ui.MasterFrameBModeRadio->isChecked())
		determineBoundingBoxFromMasterFrame();
	else if ( ui.ManualBModeRadio->isChecked() )
		determineMinimumBoundingBoxManualFromBMode();
	else
		determineMinimumBoundingBoxFromBMode();

	// identify the sender
	QPushButton *_tempButton = (QPushButton *)sender();
	if ( _tempButton != NULL)
	if ( strcmp(_tempButton->objectName(). toLatin1().data(),"DopplerCompoundButton") == 0 )
	{
		// we want to reconstruct Doppler data, this need to be handled as a special case
		_doppler = true;
	}

	m_slicesVisualization->defineOrientedBoundingBox(m_extent, m_boundingBoxTransformation);

	determineVoxelSize();

	// refresh the widget once every minute (to speed up computation)
	glWidget->stopTimer();

	// special handling of Doppler compounding
	if ( !_doppler )
	{
		// compound the volume with specified compounding mode such as gaussian distance weighting, nearest neighbor, mean
		compoundVolume((COMPOUNDING_MODE)ui.reconstructionModeComboBox->currentIndex ());
	}
	else // Dopple reconstruction
	{
		// grab parameters from the 
		float _tau = atof(ui.TauEdit->text(). toLatin1().data());
		float _lambda1 = atof(ui.Lambda1Edit->text(). toLatin1().data());
		float _lambda2 = atof(ui.Lambda2Edit->text(). toLatin1().data());
		int _iterations = atoi(ui.IterationsEdit->text(). toLatin1().data());
		compoundVolumeVelocity(_tau, _lambda1, _lambda2, _iterations);
	}
	// refresh the widget every couple of ms, as we don't need CPU performance so much anymore...let's waste it
	glWidget->restartTimer();

	// free the XML data handler
	if ( m_xmlTools ) 
	{
		delete m_xmlTools;
		m_xmlTools = 0;
	}
}

void Compounding::updateVoxelDimension(int position)
{
	m_voxelSize = (float)position * 0.05;
	determineVoxelSize();
}

Compounding::~Compounding()
{
	if ( m_xmlTools ) delete m_xmlTools;
	if ( m_rfImage != 0)
		DeleteObject(m_rfImage);
	if ( m_rfPreviewImage != 0 )
		DeleteObject(m_rfPreviewImage);
	{
		int _elements = m_subvolumeWriting.size();
		for(int i=0;i<_elements;i++)
		{
			CompoundRFSubvolumeWritingThread* _thread = m_subvolumeWriting.back();

			if (_thread) delete _thread;

			m_subvolumeWriting.pop_back();
		}
	}
	for (int j=0; j<m_rfFileStreamVector.size();j++)
	{
		RFFileStream *_rfFileStream = m_rfFileStreamVector.back();
		m_rfFileStreamVector.pop_back();

		if ( _rfFileStream) delete _rfFileStream;
	}
	if ( m_imageForVisualization ) DeleteObject(m_imageForVisualization);
	if ( m_masterImage ) DeleteObject(m_masterImage);

	if ( m_videoTools) delete m_videoTools;
	if ( m_numericTools) delete m_numericTools;
	// save the configuration back to the XML file
	m_xmlConfiguration.saveToFile(std::string(COMPOUNDING_CONFIGURATION));

	delete m_slicesVisualization, glWidget, m_offlineVisualization;

	//clear the scanline vector
	int _elements = m_rfDataVector.size(); 
	for(int j=0;j < _elements;j++)
	{
		// data conversion
#ifdef USE_USHORT16
		StandardRF<unsigned short> *_rf = reinterpret_cast<StandardRF<unsigned short>*>(m_rfDataVector.back());
		m_rfDataVector.pop_back();
		// data conversion
		unsigned short *_data = _rf->getData();

		//if  (_data) delete[] _data; // free the image data
		delete _rf;
#endif

#ifdef USE_FLOAT32 || USE_FLOAT5X32
		StandardRF<float> *_rf = reinterpret_cast<StandardRF<float>*>(m_rfDataVector.back());
		m_rfDataVector.pop_back();
		// data conversion
		float *_data = _rf->getData();

		//if  (_data) delete[] _data; // free the image data
		delete _rf;
#endif
	}
	// all pointers have been cleared now empty the vector
	m_rfDataVector.clear();

	//clear the US image vector if it contains images
	for (int j=0;j < m_imageVector.size();j++)
	{
		for (int i=0;i < m_imageVector[j].size();i++)
		{
			IplImage* _img = m_imageVector[j][i];

			if ( _img )
				cvReleaseImage(&_img);

		}
	}
#ifdef IMAGE_FILTER
	// clear the auxilary image vector if it contains stuff
	for (int j=0;j < m_filterVector.size();j++)
	{
		for (int i=0;i < m_filterVector[j].size();i++)
		{
			IplImage* _img = m_filterVector[j][i];

			if ( _img )
				cvReleaseImage(&_img);
		}
	}
#endif
	// clear the auxilary image vector if it contains stuff
	for (int j=0;j < m_auxImageVector.size();j++)
	{
		for (int i=0;i < m_auxImageVector[j].size();i++)
		{
			IplImage* _img = m_auxImageVector[j][i];

			if ( _img )
				cvReleaseImage(&_img);
		}
	}

	if ( m_croppingMask )
		cvReleaseImage(&m_croppingMask);
	//_CrtDumpMemoryLeaks();


	// Doppler related variables
	if ( m_recon ) 
		delete m_recon;

	// free the Doppler LUT
	if ( m_colLUT )
		delete m_colLUT;
}

void Compounding::guiUpdated()
{
	if ( m_guiVideoReady && m_guiTrackingReady && m_guiConfigurationReady)
	{
		ui.loadDataButton->setDisabled(false);
		ui.LoadDopplerButton->setDisabled(false);
		//ui.dataScrollBar->setDisabled(false);
	} 
	else
	{
		ui.loadDataButton->setDisabled(true);
		ui.LoadDopplerButton->setDisabled(true);
	}
	// RF specific GUI
	if ( m_guiRFCalibrationReady && m_guiRFDataReady && m_guiRFTrackingReady)
	{
		ui.loadDataRFButton->setDisabled(false);
	}
	else
	{
		ui.loadDataRFButton->setDisabled(true);
	}
}

void Compounding::selectConfigurationFile()
{
	// initialisation of the output values
	//QString configurationFilter  = tr("Configuration files (*.xml)");		
	QStringList filterList;// = ((QStringList)configurationFilter);
	filterList << "*.xml";
	// creation of a new file dialog
	QFileDialog* fileDialog = new QFileDialog(this);
	fileDialog->setWindowTitle("Loading Data...");
	fileDialog->setNameFilters(filterList);

	// check if we previously loaded some data from a specific directory, if not use the one from the XML file
	if ( m_lastCalibrationPath.length() == 0 )
	{
		std::string _path;
		m_xmlConfiguration.handleData(false,_path, std::string(CALIBRATION_DIRECTORY_TAG));
		m_lastCalibrationPath = QString(_path.c_str());
	}

	QString fileName = fileDialog->getOpenFileName(0, "Select file ...", m_lastCalibrationPath, "Configuration files (*.xml)");
	std::string _temp = fileName. toLatin1().data();
	std::string _path = _temp.substr( 0, _temp.find_last_of( '/' ) +1 );
	m_guiConfigurationReady = !fileName.isEmpty();
	if(!fileName.isEmpty())
	{
		// write the directory used last time to the configuration file to remember it for the next time
		m_xmlConfiguration.handleData(true,_path, std::string(CALIBRATION_DIRECTORY_TAG));
		m_lastCalibrationPath = QString(_path.c_str());
		ui.configurationFileEdit->setText(fileName);
		m_configurationFileName = fileName;

		// store the filename in the XML, so we can use it next time
		m_xmlConfiguration.handleData(true, _temp, std::string(CALIBRATION_FILENAME_TAG));
	}
}

void Compounding::selectTrackingFile()
{
	// initialisation of the output values
	QStringList filterList;
	filterList << "*.txt";

	// creation of a new file dialog
	QFileDialog* fileDialog = new QFileDialog(this);
	fileDialog->setWindowTitle("Loading Data...");
	fileDialog->setNameFilters(filterList);


	// check if we previously loaded some data from a specific directory, if not use the one from the XML file
	if ( m_lastPath.length() == 0 )
	{
		std::string _path;
		m_xmlConfiguration.handleData(false,_path, std::string(DATA_DIRECTORY_TAG));
		m_lastPath = QString(_path.c_str());
	}


	QString fileName = fileDialog->getOpenFileName(0,"Select file ...", m_lastPath,"Tracking data (*.txt)");
	std::string _temp = fileName. toLatin1().data();
	std::string _path = _temp.substr( 0, _temp.find_last_of( '/' ) +1 );
	m_guiTrackingReady = !fileName.isEmpty();
	if(!fileName.isEmpty())
	{
		// write the directory used last time to the configuration file to remember it for the next time
		m_xmlConfiguration.handleData(true,_path, std::string(DATA_DIRECTORY_TAG));
		m_lastPath = QString(_path.c_str());
		ui.trackingFileEdit->setText(fileName);
		m_trackingFileName = fileName;

		// store the filename in the XML, so we can use it next time
		m_xmlConfiguration.handleData(true, _temp, std::string(TRACKING_FILENAME_TAG));
	}
}

void Compounding::selectVideoFile()
{
	// initialisation of the output values
	//QString movieFilter  = tr("Video files (*.avi)");		
	QStringList filterList;
	filterList << "*.avi";

	// creation of a new file dialog
	QFileDialog* fileDialog = new QFileDialog(this);
	fileDialog->setWindowTitle("Loading Data...");
	fileDialog->setNameFilters(filterList);

	// check if we previously loaded some data from a specific directory, if not use the one from the XML file
	if ( m_lastPath.length() == 0 )
	{
		std::string _path;
		m_xmlConfiguration.handleData(false,_path, std::string(DATA_DIRECTORY_TAG));
		m_lastPath = QString(_path.c_str());
	}

	QString fileName = fileDialog->getOpenFileName(0,"Select file ...", m_lastPath ,"Video File (*.avi)");
	std::string _temp = fileName. toLatin1().data();
	std::string _path = _temp.substr( 0, _temp.find_last_of( '/' ) +1 );
	m_guiVideoReady = !fileName.isEmpty();
	if(!fileName.isEmpty())
	{
		// write the directory used last time to the configuration file to remember it for the next time
		m_xmlConfiguration.handleData(true,_path, std::string(DATA_DIRECTORY_TAG));
		m_lastPath = QString(_path.c_str());
		ui.videoFileEdit->setText(fileName);
		m_videoFileName = fileName;

		// store the filename in the XML, so we can use it next time
		m_xmlConfiguration.handleData(true, _temp, std::string(VIDEO_FILENAME_TAG));
	}
}

void Compounding::showFrame(int frame)
{
	if ( m_imageCounter > 0 )
	{
		HBITMAP bmp;
		GLuint tex;

		int _width, _height;
		//m_videoTools->getFrame(frame, bmp, _width, _height);

		IplImage* _img;
		//m_videoTools->getFrame(frame, _bmp, _width, _height);
		int _index1 = (int)floor((float)frame/NUMBER_OF_SLICES); 
		int _index2 = frame % (int)(NUMBER_OF_SLICES);
		_img = m_imageVector[_index1][_index2];

		IplImage* _colImage = cvCreateImage(cvSize(_img->width, _img->height), IPL_DEPTH_8U, 3);
		cvCvtColor(_img, _colImage, CV_GRAY2RGB);
		bmp = m_videoTools->IplImage2HBITMAP(_colImage);
		cvReleaseImage(&_colImage);


		/*	int _minX, _minY, _maxX, _maxY;
		m_videoTools->getImageMinMax(bmp, _minX, _maxX, _minY, _maxY);


		*/
		m_graphicTools.Bitmap2GLTexture(bmp, tex);

		if ( m_imageForVisualization != 0)
			DeleteObject(m_imageForVisualization);

		m_imageForVisualization = bmp;

		/*
		VIDEO IMAGE VISUALIZATION

		*/

		// SLICES VISUALIZATION
		std::stringstream ss; 
		ss  << "Frame: "<< frame << " / " << m_imageCounter;
		ui.labelFrame->setText( ss.str().c_str() );

		glWidget->stopTimer();
		m_slicesVisualization->highlightSlice(true, frame, tex, bmp, _width, _height);
		glWidget->restartTimer();
	}

	if ( m_rfFrameCounter > 0 )
	{
		// SLICES VISUALIZATION
		std::stringstream ss; 
		ss  << "Frame: "<< frame << " / " << m_imageCounter;
		ui.labelFrame->setText( ss.str().c_str() );

		glWidget->stopTimer();
		m_slicesVisualization->highlightSlice(true, frame);
		FrameScanlineRange _fsr = m_frameScanlineRangeVector[frame];

		//m_slicesVisualization->drawScanlines(true, (int)_fsr.start , (int)_fsr.end, &m_scanlineDataVector);
		glWidget->restartTimer();
	}
}

bool Compounding::determineImageSize()
{
	HBITMAP bmp, seg_bmp;

	// determine the name of the segmentation movie from the video file name

	//m_videoSegmentationTools.loadFromFile(m_videoFileName. toLatin1().data());
	std::string _videoFileName = m_videoFileName. toLatin1().data();

	// load the video images from the specified file
	if ( m_videoTools->loadFromFile(m_videoFileName. toLatin1().data()) && m_videoTools->getFrame(0, bmp, m_imageWidth, m_imageHeight)  )
	{
		DeleteObject(bmp);
		return true;
	}
	else
	{
		return false;
	}
}

bool Compounding::loadVideoFromFile(float sigma, VIDEO_MODE mode, bool attach, bool auxiliary, bool applyNoise)
{
	{
		std::stringstream _ss;

		if ( !attach)
			_ss << "Loading video ... ";
		else
			_ss << "Attaching video ... ";
		ui.textEdit->append(QString(_ss.str().c_str()));
	}


	if ( !attach )
	{
		// reset the image counter, to store it at the right place in the vector
		m_imageCounter = 0;
		// if vector is not empty, empty it
		for (int j=0;j < m_imageVector.size();j++)
			for (int i=0;i < m_imageVector[j].size();i++)
			{
				IplImage* _img = m_imageVector[j][i];

				if ( _img )
					cvReleaseImage(&_img);
			}

#ifdef IMAGE_FILTER
			for (int j=0;j < m_filterVector.size();j++)
			for (int i=0;i < m_filterVector[j].size();i++)
			{
				IplImage* _img = m_filterVector[j][i];

				if ( _img )
					cvReleaseImage(&_img);
			}
#endif
			// clear the auxilary image vector if it contains stuff
			for (int j=0;j < m_auxImageVector.size();j++)
			{
				for (int i=0;i < m_auxImageVector[j].size();i++)
				{
					IplImage* _img = m_auxImageVector[j][i];

					if ( _img )
						cvReleaseImage(&_img);
				}
			}
	}
	else
	{
	}

	HBITMAP bmp;

	// determine the name of the segmentation movie from the video file name

	//m_videoSegmentationTools.loadFromFile(m_videoFileName. toLatin1().data());
	std::string _videoFileName = m_videoFileName. toLatin1().data();

	// auxilary file for Doppler Flow reconstruction is the video file containing the segmentation

	std::string::size_type i = _videoFileName.find(".avi");
	std::string _temp = "Segmentation.avi";
	std::string _segmentationFile = _videoFileName.replace(i, _temp.length(), _temp);


	VideoTools *_auxilaryVideoTools = 0;

	// if we really have auxilary video data, load the first frame, otherwise ignore it
	if ( auxiliary )
	{
		_auxilaryVideoTools = new VideoTools();
		_auxilaryVideoTools->loadFromFile(_segmentationFile.c_str());
	}

#ifdef IMAGE_FILTER
	
		
		_temp = "_m.avi";
		std::string _filteredFile = _videoFileName.replace(i,_videoFileName.length()-i, _temp);


		VideoTools *_filteredVideoTools = 0;
		_filteredVideoTools = new VideoTools();
		_filteredVideoTools->loadFromFile(_filteredFile.c_str());
	
#endif



	int _width, _height;



	// load the video images from the specified file
	if ( m_videoTools->loadFromFile(m_videoFileName. toLatin1().data()) && m_videoTools->getFrame(0, bmp,  _width, _height)  )
	{

		// set the range that it equals the number of frames
		ui.progressBar->setRange(0,m_videoTools->getNumberOfFrames()-1);
		ui.progressBar->setValue(0);


		/*GLuint tex;

		// generates textures for openGL display purposes
		m_graphicTools.Bitmap2GLTexture(bmp, tex);

		if ( m_imageForVisualization != 0)
		DeleteObject(m_imageForVisualization);

		m_imageForVisualization = bmp;
		*/
		imageViewerUpdate(0);

		//glWidget->setTexture(tex);


		// define the ROI for copying the interesting parts of the video image

		m_roiRect.x = m_roiMinX;
		m_roiRect.y = m_roiMinY;
		m_roiRect.width = m_roiMaxX - m_roiMinX;
		m_roiRect.height = m_roiMaxY - m_roiMinY;

		// create a new instance of the random number generator, if we require one
		RandomNumbers* _randomNumbers = 0;

		if ( applyNoise)
			_randomNumbers = new RandomNumbers();

		for ( int c=0; c<m_videoTools->getNumberOfFrames();c++)
		{
			ui.progressBar->setValue(c);
			// update the statistics for mapping the images into the right bin
			int _index = (int)floor((float)m_imageCounter++/NUMBER_OF_SLICES);

			switch(mode)
			{
			case STANDARD:
				{
					// getting the video image
					{
						HBITMAP bmp;
						// grab the frame from the video
						int _width, _height;
						m_videoTools->getFrame(c, bmp, _width, _height);

						// transform it from BMP to OpenCV IplImage for processing AND make a gray-scale image
						IplImage *_tmpImage = m_videoTools->bitmap2IplImage(bmp, true);
						// create an empty image with the size of the bounding box of the bounding box of the cropping mask
						IplImage *_smallImage =  cvCreateImage(cvSize (m_roiRect.width, m_roiRect.height), _tmpImage->depth, _tmpImage->nChannels); 

						// extract the area of interest of the image (as determined by the bounding box of the cropping mask)
						cvSetImageROI(_tmpImage, m_roiRect);
						cvCopy(_tmpImage, _smallImage, 0);
						cvResetImageROI(_tmpImage);
						cvReleaseImage(&_tmpImage);
						DeleteObject(bmp);

#if USE_DIRECTX
						// directX flips the images, so we need to unflip
						cvFlip(_smallImage,_smallImage,0);
#endif
						//if ( c == m_videoTools->getNumberOfFrames() / 2)
						//char buffer[255];
						//sprintf(buffer,".//temp//cropped%d.jpg",c);
						//cvSaveImage(buffer,_smallImage);

						// store the cropped image part in the vector
						m_imageVector[_index].push_back( _smallImage);
#ifdef IMAGE_FILTER
						{
							_filteredVideoTools->getFrame(c, bmp, _width, _height);
							// transform it from BMP to OpenCV IplImage for processing AND make a gray-scale image
						IplImage *_tmpImage = m_videoTools->bitmap2IplImage(bmp, true);
						// create an empty image with the size of the bounding box of the bounding box of the cropping mask
						IplImage *_smallImage =  cvCreateImage(cvSize (m_roiRect.width, m_roiRect.height), _tmpImage->depth, _tmpImage->nChannels); 

						// extract the area of interest of the image (as determined by the bounding box of the cropping mask)
						cvSetImageROI(_tmpImage, m_roiRect);
						cvCopy(_tmpImage, _smallImage, 0);
						cvResetImageROI(_tmpImage);
						cvReleaseImage(&_tmpImage);
						DeleteObject(bmp);
						m_filterVector[_index].push_back( _smallImage);
						}
#endif
					}

					// getting auxilary video image (if desired)
					{
					}
					break;
				}
			case DOPPLER:
				{

					// getting the video image and computing the velocity from it
					{
						Eigen::MatrixXf _velocityMatrix;
						_velocityMatrix = MatrixXf::Zero(m_roiRect.height, m_roiRect.width );

						HBITMAP bmp;
						// grab the frame from the video
						int _width, _height;
						m_videoTools->getFrame(c, bmp, _width, _height);

						// transform it from BMP to OpenCV IplImage for processing
						// working for intensity reconstruct: IplImage *_tmpImage = m_videoTools->bitmap2IplImage(bmp, true);
						IplImage *_tmpImage = m_videoTools->bitmap2IplImage(bmp, false);
						// create an empty image with the size of the bounding box of the bounding box of the cropping mask
						IplImage *_smallImage =  cvCreateImage(cvSize (m_roiRect.width, m_roiRect.height), _tmpImage->depth, _tmpImage->nChannels); 

						// extract the area of interest of the image (as determined by the bounding box of the cropping mask)
						cvSetImageROI(_tmpImage, m_roiRect);
						cvCopy(_tmpImage, _smallImage, 0);
						cvResetImageROI(_tmpImage);
						cvReleaseImage(&_tmpImage);
						DeleteObject(bmp);

#if USE_DIRECTX
						// directX flips the images, so we need to unflip
						cvFlip(_smallImage,_smallImage,0);
#endif



						//---------------

						std::vector<Eigen::Vector3f> _colorVec;

						// this image is going to contain the velocities - they are mapped from RGB (video) using a color LUT
						IplImage* _velocityImage = cvCreateImage(cvSize(_smallImage->width,_smallImage->height),IPL_DEPTH_32F,1);


						float *_velocityPointer = (float*)_velocityImage->imageData;
						unsigned char *_imgPointer= (unsigned char*)(_smallImage->imageData); // the R value of the RGB image
						for(int y=0;y<_smallImage->height;y++)
							for(int x=0;x<_smallImage->width;x++)
							{

								//unsigned char *_imgPointer= &(((unsigned char*)(_smallImage->imageData+ y*_smallImage->widthStep))[x*_smallImage->nChannels]); // the R value of the RGB image
								float _velocity;
								// BUG: RGB is always 0,0,0
								// get the intensity of the Doppler ultrasound image (RGB)
								int _intensityR = static_cast<int>(*_imgPointer++);
								int _intensityG = static_cast<int>(*_imgPointer++);
								int _intensityB = static_cast<int>(*_imgPointer++);


								// now the pointer points to the next pixel

								m_colLUT->colorToVelocity(_intensityR, _intensityG, _intensityB, _velocity);


								// check if we should apply some nice
								if ( applyNoise )
								{
									float _rndVal;
									_rndVal = _randomNumbers->univariateGaussian(_rndVal,sigma);

									// we now assume additive noise

									_velocity = _velocity + _rndVal;
								}

								_velocityMatrix(y,x)=_velocity;


								*_velocityPointer++=_velocity;

								Eigen::Vector3f _col((float)_intensityR,(float)_intensityG, (float)_intensityB);


								_colorVec.push_back(_col);

								/*if( y == 16 && x == 57 )
								{
								std::cout << "RGB: " << _intensityR << ", " << _intensityG << ", " << _intensityB << std::endl;
								m_colLUT->colorToVelocity(_intensityR, _intensityG, _intensityB, _velocity);
								std::cout << "Velocity: " << _velocity << std::endl;
								}*/
							}
							/*char buffer[250];
							sprintf(buffer,"Color%d.mat",c);
							saveMatlabMatrix3Dim(_colorVec,_smallImage->width, _smallImage->height, "ColorPlane", buffer);
							sprintf(buffer,"Velocity%d.mat",c);
							saveMatlabMatrix1Dim(_velocityMatrix, m_roiRect.width, m_roiRect.height, "VelocityPlane",buffer ); 

							sprintf(buffer,"C:\\cropped%d.jpg",c);
							cvSaveImage(buffer,_smallImage);

							*/	//---------------

							// store the cropped image part in the vector
							//m_imageVector.push_back( _smallImage);


							/*_velocityPointer = (float*)_velocityImage->imageData;
							for(int y=0;y<_smallImage->height;y++)
							for(int x=0;x<_smallImage->width;x++)
							{
							_velocityMatrix(y,x)=*_velocityPointer++;
							}


							char buffer[250];
							sprintf(buffer,"Velocity%d.mat",c);
							saveMatlabMatrix1Dim(_velocityMatrix, m_roiRect.width, m_roiRect.height, "VelocityPlane",buffer ); 
							*/
							cvReleaseImage(&_smallImage);
							// store the cropped image and into velocity transformed image in the vector
							m_imageVector[_index].push_back(_velocityImage);
					}

					// getting auxilary video image (if desired)
					if ( auxiliary )
					{
						HBITMAP bmp;
						// grab the frame from the video
						int _width, _height;
						_auxilaryVideoTools->getFrame(c, bmp, _width, _height);

						// transform it from BMP to OpenCV IplImage for processing AND make a gray-scale image
						IplImage *_tmpImage = _auxilaryVideoTools->bitmap2IplImage(bmp, true);
						// create an empty image with the size of the bounding box of the bounding box of the cropping mask
						IplImage *_smallImage =  cvCreateImage(cvSize (m_roiRect.width, m_roiRect.height), _tmpImage->depth, _tmpImage->nChannels); 

						// extract the area of interest of the image (as determined by the bounding box of the cropping mask)
						cvSetImageROI(_tmpImage, m_roiRect);
						cvCopy(_tmpImage, _smallImage, 0);
						cvResetImageROI(_tmpImage);
						cvReleaseImage(&_tmpImage);
						DeleteObject(bmp);

#if USE_DIRECTX
						// directX flips the images, so we need to unflip
						cvFlip(_smallImage,_smallImage,0);
#endif


						// store the cropped image part in the auxiliary image vector
						m_auxImageVector[_index].push_back( _smallImage);

					}

					break;
				}
			} // end of switch

		}


		// free the memory of the random number generator in case it was instantiated
		if ( _randomNumbers )
			delete _randomNumbers;

		// free memory if there was reserved something
		if ( _auxilaryVideoTools )
			delete _auxilaryVideoTools;

#ifdef IMAGE_FILTER
		if ( _filteredVideoTools )
			delete _filteredVideoTools;
#endif

		// update the GUI a little bit
		ui.dataScrollBar->setDisabled(false);
		ui.dataScrollBar->setMinimum(0);
		ui.dataScrollBar->setMaximum(m_imageCounter-1);

		ui.ImageLabelScrollBar->setDisabled(false);
		ui.ImageLabelScrollBar->setMinimum(0);
		ui.ImageLabelScrollBar->setMaximum(m_imageCounter-1);
		ui.ImageLabelScrollBar->setValue(0);
		// now we allow the editing of the number
		ui.MasterFrameEdit2->setEnabled(true);
		ui.bookmarkButton->setEnabled(true);

		std::stringstream ss; 
		ss  << "Frame: 0" << " / " <<m_imageCounter;
		ui.labelFrame->setText( ss.str().c_str() );


		return true;
	}
	else
	{
		// free memory if there was reserved something
		if ( _auxilaryVideoTools )
			delete _auxilaryVideoTools;

		return false;
	}
}

bool Compounding::loadTimeStampsFromFile(QString filename, bool RF, bool attach)
{
	std::string _filename;
	// load the time stamps for the video images from the specified file
	if (!RF )
		_filename = (filename+".ts"). toLatin1().data();
	else
		_filename = filename. toLatin1().data();

	if ( m_dataTools.getVideoTimestampsFromFile(_filename,m_videoTimeStampVector) )
	{
		std::stringstream ss;

		ss << "Time Stamps: " << m_videoTimeStampVector.size();

		// update the GUI, depending on whether it is RF or B-MODE take the appropriate panel
		if ( !RF )
			ui.labelTimestamps->setText(ss.str().c_str());
		else
			ui.labelTimestampsRF->setText(ss.str().c_str());

		int _start = 0;

		// if we attach, then we have to set the last endpoint as startpoint 
		if ( attach )
			_start = m_toVideoTimeStamp;
		if ( m_temporalOffset > 0.0 )
		{
			// substract the temporal offset from the video frames
			for (int i=_start;i<m_videoTimeStampVector.size();i++)
			{
				m_videoTimeStampVector[i] -= m_temporalOffset;
			}
		}
		// store the indices of the time stamps, in case we would like to attach data
		m_fromVideoTimeStamp = _start;
		m_toVideoTimeStamp = m_videoTimeStampVector.size();

		return true;
	}
	else
		return false;
}

bool Compounding::loadSettingFromXMLFile(QString filename)
{
	// load the ultrasound calibration form the XML file
	if ( m_xmlTools->loadFile(filename. toLatin1().data()) )
	{
		if ( m_xmlTools->handleData(m_matCalibration,"Matrix") )
			m_slicesVisualization->setCalibrationMatrix(m_matCalibration);
		else
		{
			//report some error of not existing calibration matrix
		}

		if ( m_xmlTools->handleData(m_scaleX, "Scale", "ScaleX") )
		{
		}
		else
		{
			//report some error of not existing scale-X value
		}

		if ( m_xmlTools->handleData(m_scaleY, "Scale", "ScaleY") )
		{
		}
		else
		{
			//report some error of not existing scale-X value
		}

		if ( m_xmlTools->handleData(m_apexX, "Geometry", "ApexX") )
		{
		}
		else
		{
			//report some error of not existing Apex-X value
		}

		if ( m_xmlTools->handleData(m_apexY, "Geometry", "ApexY") )
		{
		}
		else
		{
			//report some error of not existing Apex-X value
		}

		if ( m_xmlTools->handleData(m_innerRadius, "Geometry", "InnerRadius") )
		{
		}
		else
		{
			//report some error of not existing Apex-X value
		}

		if ( m_xmlTools->handleData(m_outerRadius, "Geometry", "OuterRadius") )
		{
		}
		else
		{
			//report some error of not existing Apex-X value
		}

		if ( m_xmlTools->handleData(m_angle, "Geometry", "Angle") )
		{
		}
		else
		{
			//report some error of not existing Apex-X value
		}


		if ( m_xmlTools->handleData(m_temporalOffset, "Temporal", "Offset") )
		{
		}
		else
		{
			//report some error of not existing Apex-X value
		}



		return true;
	}
	else
		return false;
}


bool Compounding::loadSettingFromXMLFile(QString filename, UltrasoundSettings &ultrasoundSettings)
{
	// load the ultrasound calibration form the XML file
	if ( m_xmlTools->loadFile(filename. toLatin1().data()) )
	{
		if ( m_xmlTools->handleData(ultrasoundSettings.calibrationMatrix,"Matrix") )
			m_slicesVisualization->setCalibrationMatrix(ultrasoundSettings.calibrationMatrix);
		else
		{
			//report some error of not existing calibration matrix
		}


		if ( m_xmlTools->handleData(ultrasoundSettings.apexX, "Geometry", "ApexX") )
		{
		}
		else
		{
			//report some error of not existing Apex-X value
		}

		if ( m_xmlTools->handleData(ultrasoundSettings.apexY, "Geometry", "ApexY") )
		{
		}
		else
		{
			//report some error of not existing Apex-X value
		}

		if ( m_xmlTools->handleData(ultrasoundSettings.scaleX, "Scale", "ScaleX") )
		{
		}
		else
		{
			//report some error of not existing scale-X value
		}

		if ( m_xmlTools->handleData(ultrasoundSettings.scaleY, "Scale", "ScaleY") )
		{
		}
		else
		{
			//report some error of not existing scale-X value
		}



		if ( m_xmlTools->handleData(ultrasoundSettings.angle, "Geometry", "Angle") )
		{
		}
		else
		{
			//report some error of not existing Apex-X value
		}


		if ( m_xmlTools->handleData(ultrasoundSettings.probeRadius, "Probe", "Radius") )
		{
		}
		else
		{
			//report some error of not existing Apex-X value
			ultrasoundSettings.probeRadius = 0; // by default assume phased array
		}

		if ( m_xmlTools->handleData(ultrasoundSettings.probeWidth, "Probe", "Width") )
		{
		}
		else
		{
			//report some error of not existing Apex-X value
			ultrasoundSettings.probeWidth = 0; // by default assume phased array, curvilinear etc.
		}


		if ( m_xmlTools->handleData(ultrasoundSettings.probeAngle, "Probe", "Angle") )
		{
		}
		else
		{
			//report some error of not existing Apex-X value
			ultrasoundSettings.probeAngle = 90.; // by default assume phased array
		}



		if ( m_xmlTools->handleData(m_temporalOffset, "Temporal", "Offset") )
		{
		}
		else
		{
			//report some error of not existing Apex-X value
		}

		return true;
	}
	else
		return false;
}

bool Compounding::cropMask()
{
	// now we will have to crop the cropping mask to the size of the bounding box


	// the new image has the size of the bounding box of the cropping mask
	IplImage *_smallImage =  cvCreateImage(cvSize (m_roiRect.width, m_roiRect.height), m_croppingMask->depth, m_croppingMask->nChannels); 

	// copy just the relevant part within the bounding box
	cvSetImageROI(m_croppingMask, m_roiRect);
	cvCopy(m_croppingMask, _smallImage, 0);
	cvResetImageROI(m_croppingMask);
	cvReleaseImage(&m_croppingMask);
	// now copy the pointer
	m_croppingMask = _smallImage;


	cvSaveImage("croppingMask.jpg",m_croppingMask);

	return true;
}


bool Compounding::determineCroppingArea()
{
	// generate the cropping mask (NOTE: this has to be done within a valid OpenGL context!!)
	if ( m_croppingMask ) 
		cvReleaseImage(&m_croppingMask);

	m_croppingMask = m_offlineVisualization->renderUSCroppingMask(m_imageWidth, m_imageHeight, m_apexX, m_apexY, m_innerRadius, m_outerRadius, m_angle);

	// now we determine the bounding box of the cropping mask (thus we can later on copy only the relevants parts of the ultrasound images)
	bool _found = m_videoTools->getImageMinMax(m_croppingMask, m_roiMinX, m_roiMaxX, m_roiMinY, m_roiMaxY);

	// check if we determined a valid ROI
	if ( _found )
	{
		m_roiHeight = m_roiMaxY - m_roiMinY;
		m_roiWidth = m_roiMaxX - m_roiMinX;
	}
	else // if we have no valid ROI, we set the whole image as ROI (this is a hack for invalid cropping regions such as very small images generated for doppler testing)
	{
		m_roiHeight = m_imageHeight;
		m_roiWidth = m_imageWidth;

		m_roiMinX = 0;
		m_roiMinY = 0;

		m_roiMaxX = m_imageWidth;
		m_roiMaxY = m_imageHeight;
	}

	return true;
}

bool Compounding::loadTrackingDataFromFile(QString filename, bool RF, bool useReferenceTarget, bool attach)
{
	// in case we don't attach, clear the lists containing the indices indicating the beginning of sequences
	if ( !attach) 
		m_indexVector.clear();

	// load the tracking data from the specified file and store it in a vector

	bool _loadedData = false;

	// use "ReferenceTarget" as the moving reference coordinate system, to cope with patient motion 
	if ( useReferenceTarget )
	{
		// update the index, from where data is being added
		if ( !attach )
			m_trackingIndex = 0;
		else
			m_trackingIndex = m_vecTracking.size()-1; // staring index is 0 (computer scientists joy)

		// get the probe name from the appropriate input
		if (!RF)
			_loadedData = m_dataTools.getTrackingDataFromFileRelative(filename. toLatin1().data(),ui.probeNameComboBox->currentText(). toLatin1().data(), ui.referenceNameComboBox->currentText(). toLatin1().data(), m_vecTracking);
		else
			_loadedData = m_dataTools.getTrackingDataFromFileRelative(filename. toLatin1().data(),ui.probeNameRFComboBox->currentText(). toLatin1().data(), ui.referenceNameRFComboBox->currentText(). toLatin1().data(), m_vecTracking);

	}
	else
	{

		// update the index, from where data is being added
		if ( !attach )
			m_trackingIndex = 0;
		else
			m_trackingIndex = m_vecTracking.size()-1; // staring index is 0 (computer scientists joy)

		// work in tracking system coordinate frame
		// get the probe name from the appropriate input
		if (!RF)
			_loadedData = m_dataTools.getTrackingDataFromFile(filename. toLatin1().data(), ui.probeNameComboBox->currentText(). toLatin1().data(), m_vecTracking);
		else
			_loadedData = m_dataTools.getTrackingDataFromFile(filename. toLatin1().data(), ui.probeNameRFComboBox->currentText(). toLatin1().data(), m_vecTracking);
	}



	if ( _loadedData )
	{

		/*
		if ( m_dataTools.getTrackingDataFromFile(m_trackingFileName. toLatin1().data(), ui.probeNameComboBox->currentText(). toLatin1().data(), m_vecTracking) )
		{*/
		std::stringstream ss;

		ss << "Tracking data: " << m_vecTracking.size();
		// choose the right panel for display RF or B-MODE
		if (!RF)
			ui.labelTracking->setText(ss.str().c_str());
		else
			ui.labelTrackingRF->setText(ss.str().c_str());

		// experimental
		Eigen::Matrix4f _H ;
		// --

		// there is no temporal offset!!
		/*
		for(int i=0;i<vecTracking.size();i++)
		{
		TrackingData  td = vecTracking[i];

		Eigen::Transform3f _tmpTrackingMatrix;
		Eigen::Matrix3f _tmpRotationMatrix;
		Eigen::Vector3f _tmpTranslationVector;
		_tmpRotationMatrix = Eigen::Quaterniond(td.quaternion[0], td.quaternion[1], td.quaternion[2], td.quaternion[3]);
		_tmpTranslationVector = Eigen::Vector3f(td.translation[0], td.translation[1], td.translation[2]);


		//here must be the bug
		_tmpTrackingMatrix.setIdentity();
		_tmpTrackingMatrix.rotate( _tmpRotationMatrix);
		_tmpTrackingMatrix.translation() = _tmpTranslationVector;

		//	std::cout << _tmpTrackingMatrix.matrix() << std::endl;


		Eigen::Matrix4f _tempData;
		_tempData = _tmpTrackingMatrix * _matCalibration;
		_matrixDataTracking.push_back(_tempData);

		_trackingData.push_back(_tmpTrackingMatrix);
		}
		*/

		int _futureIndex = 0;
		int _pastIndex = 0;

		if ( attach )
		{
			_futureIndex = m_trackingIndex+1;
			_pastIndex = m_trackingIndex+1;
		}

		float _minFutureDist, _minPastDist;
		// interpolation between the 
		for(int i=m_fromVideoTimeStamp;i<m_toVideoTimeStamp;i++,_minFutureDist = std::numeric_limits<float>::max(), _minPastDist = std::numeric_limits<float>::max())
		{
			//TrackingData _tempTD = vecTracking[i];
			// find closest tracking data

			float _tempStamp = m_videoTimeStampVector[i];
			float _tempTrStamp0 = m_vecTracking[_pastIndex].timestamp;
			float _tempTrStamp1 = m_vecTracking[_pastIndex+1].timestamp;


			while(_pastIndex+1 < m_vecTracking.size() &&  m_vecTracking[_pastIndex+1].timestamp < m_videoTimeStampVector[i] )
			{
				_pastIndex++;
				//_minPastDist = fabs(vecTracking[_pastIndex].timestamp-vecTimeStamp[i]);
			};

			if ( _futureIndex+1 < m_vecTracking.size() )
				_futureIndex = _pastIndex+1;


			float _t = 0;

			float _first  = m_vecTracking[_pastIndex].timestamp;
			float _second = m_vecTracking[_futureIndex].timestamp;
			float _video = m_videoTimeStampVector[i];
			float _diff = _second - _first;
			Eigen::Vector3f _transFirst;
			_transFirst << m_vecTracking[_pastIndex].translation[0], m_vecTracking[_pastIndex].translation[1], m_vecTracking[_pastIndex].translation[2];
			Eigen::Vector3f _transSecond;
			_transSecond << m_vecTracking[_futureIndex].translation[0], m_vecTracking[_futureIndex].translation[1], m_vecTracking[_futureIndex].translation[2];

			Eigen::Quaternionf _firstQuat( m_vecTracking[_pastIndex].quaternion[0], m_vecTracking[_pastIndex].quaternion[1], m_vecTracking[_pastIndex].quaternion[2], m_vecTracking[_pastIndex].quaternion[3]);
			Eigen::Quaternionf _secondQuat( m_vecTracking[_futureIndex].quaternion[0], m_vecTracking[_futureIndex].quaternion[1], m_vecTracking[_futureIndex].quaternion[2], m_vecTracking[_futureIndex].quaternion[3]);
			Eigen::Quaternionf _tempQuat;
			Eigen::Vector3f _tempTranslation;

			if ( _pastIndex != _futureIndex && m_temporalOffset > 0.0 )
			{

				_t =  (_video - _first) / (_second - _first);
				// make sure we interpolate along the shortest path
				// quaternions with opposite sign represent the same rotation
				// however, their differing components induce different interpolation ways
				// that is short vs. long way, therefore we have to check their dot product
				// and change sign on ALL components if necessary
				if ( _firstQuat.dot(_secondQuat)< 0.0)
				{
					_secondQuat.w() *=-1.0; 
					_secondQuat.x() *=-1.0;
					_secondQuat.y() *=-1.0;
					_secondQuat.z() *=-1.0;
				}
				_tempQuat = _firstQuat.slerp(_t, _secondQuat);
				_tempTranslation = _t*(_transSecond-_transFirst)+_transFirst;
			}
			else
			{
				_tempQuat = _firstQuat;
				_tempTranslation = _transFirst;
			}

			Eigen::Transform3f _tmpTrackingMatrix;
			_tmpTrackingMatrix.setIdentity();
			_tmpTrackingMatrix.rotate( _tempQuat.toRotationMatrix());
			_tmpTrackingMatrix.translation() = _tempTranslation;


			Eigen::Matrix4f _tempData;
			_tempData = _tmpTrackingMatrix.matrix() * m_calibrationMatrixVector[m_indexVector.size()];//m_matCalibration;
			m_matrixDataTracking.push_back(_tempData);



			m_trackingData.push_back(_tmpTrackingMatrix);


	//		std::cout << _tmpTrackingMatrix.matrix() << std::endl;

			//std::cout << "Index: " << _futureIndex << " vs. " << i << "    Second: "<< _pastIndex  << "  Interpolation: " << _t<< std::endl;

		}

		// VERY IMPORTANT: this provides the index to the corresponding calibration data and ultrasound settings!!!!
		m_indexVector.push_back(m_trackingData.size());

		// DEBUG CODE
		/*
		std::cout << (m_trackingData.at(367).matrix()) << std::endl;

		std::cout << "-----------------------------" << std::endl;

		std::cout << (m_trackingData.at(800).matrix()) << std::endl;
		*/
		// DEBUG CODE - END

		return true;
	}
	else
		return false;
}

void Compounding::exportTrackingData()
{
	DataTools _dt;
	_dt.writeTrackingDataToFile(m_trackingData,"TransducerInPhantomCoordinates","TrackingData.txt");
}

void  Compounding::clickedRFRadio()
{
	// update the transformation

	if ( ui.MasterFrameRFRadio->isChecked()) 
	{
		determineMinimumBoundingBoxfromMasterFrameRF();
	}

	else if ( ui.AutomaticRFRadio->isChecked())
	{
		determineMinimumBoundingBoxfromRF();
		m_masterFrameXML = -1; 
	}

	else if ( ui.ManualRFRadio->isChecked())
	{
		determineManualBoundingBoxfromRF();
		m_masterFrameXML = -1; 
	}
}

void  Compounding::clickedBModeRadio()
{
	// update the transformation

	if ( ui.MasterFrameBModeRadio->isChecked())
		determineBoundingBoxFromMasterFrame();
	else if ( ui.ManualBModeRadio->isChecked() )
		determineMinimumBoundingBoxManualFromBMode();
	else
		determineMinimumBoundingBoxFromBMode();
}

bool Compounding::determineBoundingBoxFromMasterFrame()
{


	int _masterFrameIndex = atoi(ui.MasterFrameEdit->text(). toLatin1().data());
	ui.ImageLabelScrollBar->setValue(_masterFrameIndex);

	// correct for wrong user input
	if (_masterFrameIndex > m_imageCounter  )
	{
		_masterFrameIndex = m_imageCounter;
		std::stringstream _ss;
		_ss << m_imageCounter;
		ui.MasterFrameEdit->setText(_ss.str().c_str());
		ui.MasterFrameEdit2->setText(_ss.str().c_str());
	}

	if ( _masterFrameIndex < 0 )
	{
		_masterFrameIndex = 0;
		std::stringstream _ss;
		_ss << "0";
		ui.MasterFrameEdit->setText(_ss.str().c_str());
		ui.MasterFrameEdit2->setText(_ss.str().c_str());
	}



	// determine the spatial bounding box for the US image

	float _bbMinX, _bbMinY, _bbMinZ;

	_bbMinX = std::numeric_limits<float>::max();
	_bbMinY = std::numeric_limits<float>::max();
	_bbMinZ = std::numeric_limits<float>::max();

	float _bbMaxX, _bbMaxY, _bbMaxZ;

	_bbMaxX = -std::numeric_limits<float>::max();
	_bbMaxY = -std::numeric_limits<float>::max();
	_bbMaxZ = -std::numeric_limits<float>::max();

	std::vector<Eigen::Vector3f> _dataPoints;


	Eigen::Matrix4f _referenceINV = (m_trackingData[_masterFrameIndex].matrix() * m_matCalibration ).inverse();


	for(int i=0; i<m_trackingData.size();i++)
	{
		Eigen::Transform3f _trans = m_trackingData[i];

		HBITMAP bmp;
		//m_videoSegmentationTools.getFrame(i, bmp);

		int _minX=m_imageMinX, _minY=m_imageMinY, _maxX=m_imageMaxX, _maxY=m_imageMaxY;

		//if ( m_videoSegmentationTools.getImageMinMax(bmp, _minX, _maxX, _minY, _maxY) )
		{

			Eigen::Vector4f _exVec[4];
			_exVec[0] << m_roiMinX * m_scaleX, m_roiMinY * m_scaleY, 0, 1;
			_exVec[1] << m_roiMaxX * m_scaleX, m_roiMaxY * m_scaleY, 0, 1;
			_exVec[2] << m_roiMaxX * m_scaleX, m_roiMinY * m_scaleY, 0, 1;
			_exVec[3] << m_roiMinX * m_scaleX, m_roiMaxY * m_scaleY, 0, 1;


			for (int j=0;j<4;j++)
			{
				// transform into the coordinate system of the master frame
				_exVec[j] = ( (_trans* m_matCalibration)) * _exVec[j];


				Eigen::Vector3f _tempVec;

				_tempVec << _exVec[j](0), _exVec[j](1), _exVec[j](2);

				_dataPoints.push_back(_tempVec);

				if ( _bbMinX > _exVec[j](0) )
					_bbMinX = _exVec[j](0);
				if ( _bbMinY > _exVec[j](1) )
					_bbMinY = _exVec[j](1);
				if ( _bbMinZ > _exVec[j](2) )
					_bbMinZ = _exVec[j](2);

				if ( _bbMaxX < _exVec[j](0) )
					_bbMaxX = _exVec[j](0);
				if ( _bbMaxY < _exVec[j](1) )
					_bbMaxY = _exVec[j](1);
				if ( _bbMaxZ < _exVec[j](2) )
					_bbMaxZ = _exVec[j](2);
			}
			//char buffer[200];
			//sprintf(buffer,"%f %f %f / %f %f %f",_minVec(0),_minVec(1),_minVec(2),_maxVec(0),_maxVec(1),_maxVec(2));
			//sprintf(buffer,"%d / %d %d / %d %d / %f %f %f / %f %f %f",i,_minX,_minY,_maxX,_maxY,_bbMinX,_bbMinY,_bbMinZ,_bbMaxX,_bbMaxY,_bbMaxZ);
			//ui.textEdit->append(QString(buffer));


			qApp->processEvents();
		}
	}

	Eigen::Vector3f _originInWorld;
	Eigen::Vector3f _centroid;
	Eigen::Matrix3f _tmpRotationMatrix;
	_tmpRotationMatrix.setIdentity();

	//m_boundingBoxTransformation)
	GeometricTools _geometricTools;
	Eigen::Transform3f _tempTrans(_referenceINV);
	_geometricTools.computeCentroid(_dataPoints, _tempTrans.rotation(), _centroid, _originInWorld, m_extent);

	m_boundingBoxTransformation.setIdentity();
	m_boundingBoxTransformation.rotate((m_trackingData[_masterFrameIndex]* m_matCalibration).rotation());
	m_boundingBoxTransformation.translate(1.0*_originInWorld);

	std::stringstream _ss2;
	_ss2 << m_extent.dimX << " / " << m_extent.dimY << " / " << m_extent.dimZ << std::endl;
	_ss2 << m_boundingBoxTransformation.matrix() << std::endl;

	ui.textEdit->append(_ss2.str().c_str());

	return true;
}

bool Compounding::determineManualBoundingBoxfromRF()
{
	// determine the spatial bounding box for the US image

	float _bbMinX, _bbMinY, _bbMinZ;

	_bbMinX = std::numeric_limits<float>::max();
	_bbMinY = std::numeric_limits<float>::max();
	_bbMinZ = std::numeric_limits<float>::max();

	float _bbMaxX, _bbMaxY, _bbMaxZ;

	_bbMaxX = -std::numeric_limits<float>::max();
	_bbMaxY = -std::numeric_limits<float>::max();
	_bbMaxZ = -std::numeric_limits<float>::max();

	std::vector<Eigen::Vector3f> _dataPoints;

	{
		std::stringstream _ss2;
		_ss2 << m_trackingData.size();

		ui.textEdit->append(_ss2.str().c_str());
	}
	int _size = m_trackingData.size();
	int _dataIndex = 0;
	Eigen::Matrix4f _calibrationMatrix = m_calibrationMatrixVector[_dataIndex];
	float _penetrationDepth = m_ultrasoundSettings[_dataIndex].penetrationDepth;
	float _angle = m_ultrasoundSettings[_dataIndex].probeAngle;
	float _apexX = m_ultrasoundSettings[_dataIndex].apexX;
	float _apexY = m_ultrasoundSettings[_dataIndex].apexY;
	float _scaleX = m_ultrasoundSettings[_dataIndex].scaleX;
	float _scaleY = m_ultrasoundSettings[_dataIndex].scaleY;
	float _radius =  (float)(m_ultrasoundSettings[_dataIndex].probeRadius);
	float _width = m_ultrasoundSettings[_dataIndex].probeWidth;

	Eigen::Matrix4f _corner_T_apex;
	_corner_T_apex << 1.0, 0.0, 0.0, _apexX*_scaleX,  0.0, 1.0, 0.0, _apexY*_scaleY,  0.0, 0.0, 1.0, 0.0,  0.0, 0.0, 0.0, 1.0;
	for(int i=0; i<_size;i++)
	{
		float _theta = -_angle/2.0*M_PI/180.0;
		// check if we have to jump to the next index
		if ( i == m_indexVector[_dataIndex])
		{
			_dataIndex++;
			_calibrationMatrix = m_calibrationMatrixVector[_dataIndex];
			_penetrationDepth = m_ultrasoundSettings[_dataIndex].penetrationDepth;
			_angle = m_ultrasoundSettings[_dataIndex].probeAngle;
			_radius =  m_ultrasoundSettings[_dataIndex].probeRadius;
			_apexX = m_ultrasoundSettings[_dataIndex].apexX;
			_apexY = m_ultrasoundSettings[_dataIndex].apexY;
			_scaleX = m_ultrasoundSettings[_dataIndex].scaleX;
			_scaleY = m_ultrasoundSettings[_dataIndex].scaleY;
			_width = m_ultrasoundSettings[_dataIndex].probeWidth;
			if ( _width < 1.0 ) // curvilinear probe
				_corner_T_apex << 1.0, 0.0, 0.0, _apexX*_scaleX,  0.0, 1.0, 0.0, _apexY*_scaleY,  0.0, 0.0, 1.0, 0.0,  0.0, 0.0, 0.0, 1.0;
			else
			_corner_T_apex << 1.0, 0.0, 0.0, 0.0,  0.0, 1.0, 0.0, 0.0,  0.0, 0.0, 1.0, 0.0,  0.0, 0.0, 0.0, 1.0;
		}

		Eigen::Transform3f _trans = m_trackingData[i];

		HBITMAP bmp;
		//m_videoSegmentationTools.getFrame(i, bmp);


		//if ( m_videoSegmentationTools.getImageMinMax(bmp, _minX, _maxX, _minY, _maxY) )
		{


			if ( _width < 1.0 ) // curvilinear or phased array probe
			{

			// now define the extreme points, which are the apex, the left-most, mid and right-most scanline
			Eigen::Vector4f _tmpPnts[4];
			Eigen::Vector4f _exVec[4];

			// those points are the extreme points of the sector
			// phased array assumption, scanline begings in the apex
			if ( _radius == 0) {
				_tmpPnts[0] = _corner_T_apex * Eigen::Vector4f(0.0, 0.0, 0.0, 1.0);
			}
			else // curvilinear probe, offset = inner radius of the probe
			{
				Eigen::Vector4f _tempVec = Eigen::Vector4f(_radius * sin(_theta), _radius * cos(_theta),  0.0, 1.0);
				_tmpPnts[0] = _corner_T_apex * _tempVec;
			}
			//_tmpPnts[0] = _corner_T_apex * Eigen::Vector4f(0.0, 0.0, 0.0, 1.0); // apex


			_tmpPnts[1] =  _corner_T_apex * Eigen::Vector4f( ( _radius + _penetrationDepth) * sin(_theta), ( _radius + _penetrationDepth ) * cos(_theta),  0.0, 1.0); // left-most scanline

			_theta = 0.0;
			_tmpPnts[2] = _corner_T_apex *  Eigen::Vector4f( (_radius + _penetrationDepth ) * sin(_theta), ( _radius + _penetrationDepth ) * cos(_theta),  0.0, 1.0); // central scanline

			_theta = _angle/2.0*M_PI/180.0;
			_tmpPnts[3] = _corner_T_apex *  Eigen::Vector4f( (_radius + _penetrationDepth ) * sin(_theta), ( _radius + _penetrationDepth ) * cos(_theta),  0.0, 1.0); // right-most scanline

			// now we will define the rectangular ROI containing the sector

			// left, top
			_exVec[0] << _tmpPnts[1].x(), _tmpPnts[0].y(), 0.0, 1.0;
			// left, bottom
			_exVec[1] << _tmpPnts[1].x(), _tmpPnts[2].y(), 0.0, 1.0;
			// right, bottom
			_exVec[2] << _tmpPnts[3].x(), _tmpPnts[2].y(), 0.0, 1.0;
			// right, top
			_exVec[3] << _tmpPnts[3].x(), _tmpPnts[0].y(), 0.0, 1.0;


			// now we fake the ROI values given the scanlines
			m_roiMinX = _exVec[1].x();
			m_roiMaxX = _exVec[3].x();

			m_roiMinY = _exVec[0].y();
			m_roiMaxY = _exVec[2].y();

			// now we fake scaling in X and Y diretion assuming some resolution of the image
			m_scaleX =150.0 / _penetrationDepth;
			m_scaleY = m_scaleX;

			for (int j=0;j<4;j++)
			{
				Eigen::Vector4f _temp = _trans * _calibrationMatrix *_exVec[j];


				Eigen::Vector3f _tempVec;

				_tempVec << _temp(0), _temp(1), _temp(2);

				_dataPoints.push_back(_tempVec);

				if ( _bbMinX > _tempVec(0) )
					_bbMinX = _tempVec(0);
				if ( _bbMinY > _tempVec(1) )
					_bbMinY = _tempVec(1);
				if ( _bbMinZ > _tempVec(2) )
					_bbMinZ = _tempVec(2);

				if ( _bbMaxX < _tempVec(0) )
					_bbMaxX = _tempVec(0);
				if ( _bbMaxY < _tempVec(1) )
					_bbMaxY = _tempVec(1);
				if ( _bbMaxZ < _tempVec(2) )
					_bbMaxZ = _tempVec(2);
			}
			}
			else // linear probe
			{

				// now define the extreme points, which are the apex, the left-most, mid and right-most scanline
			Eigen::Vector4f _tmpPnts[4];
			Eigen::Vector4f _exVec[4];

			_tmpPnts[0] = _corner_T_apex * Eigen::Vector4f(0.0, 0.0, 0.0, 1.0);
			

			_tmpPnts[1] =  _corner_T_apex * Eigen::Vector4f(0.0, _penetrationDepth,  0.0, 1.0); // left-most scanline

			_tmpPnts[2] = _corner_T_apex *  Eigen::Vector4f( _width/2.0, _penetrationDepth ,  0.0, 1.0); // central scanline

			_tmpPnts[3] = _corner_T_apex *  Eigen::Vector4f(_width,  _penetrationDepth ,  0.0, 1.0); // right-most scanline

			// now we will define the rectangular ROI containing the sector

			// left, top
			_exVec[0] << _tmpPnts[1].x(), _tmpPnts[0].y(), 0.0, 1.0;
			// left, bottom
			_exVec[1] << _tmpPnts[1].x(), _tmpPnts[2].y(), 0.0, 1.0;
			// right, bottom
			_exVec[2] << _tmpPnts[3].x(), _tmpPnts[2].y(), 0.0, 1.0;
			// right, top
			_exVec[3] << _tmpPnts[3].x(), _tmpPnts[0].y(), 0.0, 1.0;


			// now we fake the ROI values given the scanlines
			m_roiMinX = _exVec[1].x();
			m_roiMaxX = _exVec[3].x();

			m_roiMinY = _exVec[0].y();
			m_roiMaxY = _exVec[2].y();

			// now we fake scaling in X and Y diretion assuming some resolution of the image
			m_scaleX =150.0 / _penetrationDepth;
			m_scaleY = m_scaleX;

			for (int j=0;j<4;j++)
			{
				Eigen::Vector4f _temp = _trans * _calibrationMatrix *_exVec[j];


				Eigen::Vector3f _tempVec;

				_tempVec << _temp(0), _temp(1), _temp(2);

				_dataPoints.push_back(_tempVec);

				if ( _bbMinX > _tempVec(0) )
					_bbMinX = _tempVec(0);
				if ( _bbMinY > _tempVec(1) )
					_bbMinY = _tempVec(1);
				if ( _bbMinZ > _tempVec(2) )
					_bbMinZ = _tempVec(2);

				if ( _bbMaxX < _tempVec(0) )
					_bbMaxX = _tempVec(0);
				if ( _bbMaxY < _tempVec(1) )
					_bbMaxY = _tempVec(1);
				if ( _bbMaxZ < _tempVec(2) )
					_bbMaxZ = _tempVec(2);
			}



			}
			//char buffer[200];
			//sprintf(buffer,"%f %f %f / %f %f %f",_minVec(0),_minVec(1),_minVec(2),_maxVec(0),_maxVec(1),_maxVec(2));
			//sprintf(buffer,"%d / %d %d / %d %d / %f %f %f / %f %f %f",i,_minX,_minY,_maxX,_maxY,_bbMinX,_bbMinY,_bbMinZ,_bbMaxX,_bbMaxY,_bbMaxZ);
			//ui.textEdit->append(QString(buffer));


			qApp->processEvents();
		}
	}


	GeometricTools _geometricTools;

	Eigen::Vector3f _centroid;
	Eigen::Quaternionf _rotation;

	// get the transformation matrix from the text input field
	{
		std::string _tmpStr = ui.manualMatrixRF->toPlainText (). toLatin1().data();
		Eigen::Matrix4f _tmpMatrix;
		XMLTools _xmlTools;
		_xmlTools.string2matrix(_tmpMatrix, _tmpStr);

		m_boundingBoxTransformation.matrix() = _tmpMatrix;

		m_extent.dimX = atoi(ui.dimX->text().toLatin1().data());
		m_extent.dimY = atoi(ui.dimY->text().toLatin1().data());
		m_extent.dimZ = atoi(ui.dimZ->text().toLatin1().data());
	}
	//_geometricTools.determineExtent(_dataPoints, m_boundingBoxTransformation.matrix(), m_extent);


	//std::cout << _extent.dimX << " / " << _extent.dimY << " / " << _extent.dimZ << std::endl;

	std::stringstream _ss2;
	_ss2 << m_extent.dimX << " / " << m_extent.dimY << " / " << m_extent.dimZ << std::endl;
	_ss2 << m_boundingBoxTransformation.matrix() << std::endl;

	ui.textEdit->append(_ss2.str().c_str());
	return true;
}



bool Compounding::determineMinimumBoundingBoxfromRF()
{
	// determine the spatial bounding box for the US image

	float _bbMinX, _bbMinY, _bbMinZ;

	_bbMinX = std::numeric_limits<float>::max();
	_bbMinY = std::numeric_limits<float>::max();
	_bbMinZ = std::numeric_limits<float>::max();

	float _bbMaxX, _bbMaxY, _bbMaxZ;

	_bbMaxX = -std::numeric_limits<float>::max();
	_bbMaxY = -std::numeric_limits<float>::max();
	_bbMaxZ = -std::numeric_limits<float>::max();

	std::vector<Eigen::Vector3f> _dataPoints;

	{
		std::stringstream _ss2;
		_ss2 << m_trackingData.size();

		ui.textEdit->append(_ss2.str().c_str());
	}
	int _size = m_trackingData.size();
	int _dataIndex = 0;
	Eigen::Matrix4f _calibrationMatrix = m_calibrationMatrixVector[_dataIndex];
	float _penetrationDepth = m_ultrasoundSettings[_dataIndex].penetrationDepth;
	float _angle = m_ultrasoundSettings[_dataIndex].probeAngle;
	float _apexX = m_ultrasoundSettings[_dataIndex].apexX;
	float _apexY = m_ultrasoundSettings[_dataIndex].apexY;
	float _scaleX = m_ultrasoundSettings[_dataIndex].scaleX;
	float _scaleY = m_ultrasoundSettings[_dataIndex].scaleY;
	float _radius =  (float)(m_ultrasoundSettings[_dataIndex].probeRadius);
	float _width = m_ultrasoundSettings[_dataIndex].probeWidth;
	Eigen::Matrix4f _corner_T_apex;

	if ( _width < 1.0 ) // curvilinear or phased array probe
	{
	_corner_T_apex << 1.0, 0.0, 0.0, _apexX*_scaleX,  0.0, 1.0, 0.0, _apexY*_scaleY,  0.0, 0.0, 1.0, 0.0,  0.0, 0.0, 0.0, 1.0;
	for(int i=0; i<_size;i++)
	{
		float _theta = -_angle/2.0*M_PI/180.0;
		// check if we have to jump to the next index
		if ( i == m_indexVector[_dataIndex])
		{
			_dataIndex++;
			_calibrationMatrix = m_calibrationMatrixVector[_dataIndex];
			_penetrationDepth = m_ultrasoundSettings[_dataIndex].penetrationDepth;
			_angle = m_ultrasoundSettings[_dataIndex].probeAngle;
			_radius =  m_ultrasoundSettings[_dataIndex].probeRadius;
			_apexX = m_ultrasoundSettings[_dataIndex].apexX;
			_apexY = m_ultrasoundSettings[_dataIndex].apexY;
			_scaleX = m_ultrasoundSettings[_dataIndex].scaleX;
			_scaleY = m_ultrasoundSettings[_dataIndex].scaleY;
			_width = m_ultrasoundSettings[_dataIndex].probeWidth;
			_corner_T_apex << 1.0, 0.0, 0.0, _apexX*_scaleX,  0.0, 1.0, 0.0, _apexY*_scaleY,  0.0, 0.0, 1.0, 0.0,  0.0, 0.0, 0.0, 1.0;
		}

		Eigen::Transform3f _trans = m_trackingData[i];

		HBITMAP bmp;
		//m_videoSegmentationTools.getFrame(i, bmp);


		//if ( m_videoSegmentationTools.getImageMinMax(bmp, _minX, _maxX, _minY, _maxY) )
		{



			// now define the extreme points, which are the apex, the left-most, mid and right-most scanline
			Eigen::Vector4f _tmpPnts[4];
			Eigen::Vector4f _exVec[4];

			// those points are the extreme points of the sector
			// phased array assumption, scanline begings in the apex
			if ( _radius == 0) {
				_tmpPnts[0] = _corner_T_apex * Eigen::Vector4f(0.0, 0.0, 0.0, 1.0);
			}
			else // curvilinear probe, offset = inner radius of the probe
			{
				Eigen::Vector4f _tempVec = Eigen::Vector4f(_radius * sin(_theta), _radius * cos(_theta),  0.0, 1.0);
				_tmpPnts[0] = _corner_T_apex * _tempVec;
			}
			//_tmpPnts[0] = _corner_T_apex * Eigen::Vector4f(0.0, 0.0, 0.0, 1.0); // apex


			_tmpPnts[1] =  _corner_T_apex * Eigen::Vector4f( ( _radius + _penetrationDepth) * sin(_theta), ( _radius + _penetrationDepth ) * cos(_theta),  0.0, 1.0); // left-most scanline

			_theta = 0.0;
			_tmpPnts[2] = _corner_T_apex *  Eigen::Vector4f( (_radius + _penetrationDepth ) * sin(_theta), ( _radius + _penetrationDepth ) * cos(_theta),  0.0, 1.0); // central scanline

			_theta = _angle/2.0*M_PI/180.0;
			_tmpPnts[3] = _corner_T_apex *  Eigen::Vector4f( (_radius + _penetrationDepth ) * sin(_theta), ( _radius + _penetrationDepth ) * cos(_theta),  0.0, 1.0); // right-most scanline

			// now we will define the rectangular ROI containing the sector

			// left, top
			_exVec[0] << _tmpPnts[1].x(), _tmpPnts[0].y(), 0.0, 1.0;
			// left, bottom
			_exVec[1] << _tmpPnts[1].x(), _tmpPnts[2].y(), 0.0, 1.0;
			// right, bottom
			_exVec[2] << _tmpPnts[3].x(), _tmpPnts[2].y(), 0.0, 1.0;
			// right, top
			_exVec[3] << _tmpPnts[3].x(), _tmpPnts[0].y(), 0.0, 1.0;


			// now we fake the ROI values given the scanlines
			m_roiMinX = _exVec[1].x();
			m_roiMaxX = _exVec[3].x();

			m_roiMinY = _exVec[0].y();
			m_roiMaxY = _exVec[2].y();

			// now we fake scaling in X and Y diretion assuming some resolution of the image
			m_scaleX =150.0 / _penetrationDepth;
			m_scaleY = m_scaleX;

			for (int j=0;j<4;j++)
			{
				Eigen::Vector4f _temp = _trans * _calibrationMatrix *_exVec[j];


				Eigen::Vector3f _tempVec;

				_tempVec << _temp(0), _temp(1), _temp(2);

				_dataPoints.push_back(_tempVec);

				if ( _bbMinX > _tempVec(0) )
					_bbMinX = _tempVec(0);
				if ( _bbMinY > _tempVec(1) )
					_bbMinY = _tempVec(1);
				if ( _bbMinZ > _tempVec(2) )
					_bbMinZ = _tempVec(2);

				if ( _bbMaxX < _tempVec(0) )
					_bbMaxX = _tempVec(0);
				if ( _bbMaxY < _tempVec(1) )
					_bbMaxY = _tempVec(1);
				if ( _bbMaxZ < _tempVec(2) )
					_bbMaxZ = _tempVec(2);
			}
			//char buffer[200];
			//sprintf(buffer,"%f %f %f / %f %f %f",_minVec(0),_minVec(1),_minVec(2),_maxVec(0),_maxVec(1),_maxVec(2));
			//sprintf(buffer,"%d / %d %d / %d %d / %f %f %f / %f %f %f",i,_minX,_minY,_maxX,_maxY,_bbMinX,_bbMinY,_bbMinZ,_bbMaxX,_bbMaxY,_bbMaxZ);
			//ui.textEdit->append(QString(buffer));


			qApp->processEvents();
		}
	}
	}
	else // linear probe
	{
		_corner_T_apex << 1.0, 0.0, 0.0, 0.0,  0.0, 1.0, 0.0, 0.0,  0.0, 0.0, 1.0, 0.0,  0.0, 0.0, 0.0, 1.0;
	for(int i=0; i<_size;i++)
	{
		
		// check if we have to jump to the next index
		if ( i == m_indexVector[_dataIndex])
		{
			_dataIndex++;
			_calibrationMatrix = m_calibrationMatrixVector[_dataIndex];
			_penetrationDepth = m_ultrasoundSettings[_dataIndex].penetrationDepth;
			_angle = m_ultrasoundSettings[_dataIndex].probeAngle;
			_radius =  m_ultrasoundSettings[_dataIndex].probeRadius;
			_apexX = m_ultrasoundSettings[_dataIndex].apexX;
			_apexY = m_ultrasoundSettings[_dataIndex].apexY;
			_scaleX = m_ultrasoundSettings[_dataIndex].scaleX;
			_scaleY = m_ultrasoundSettings[_dataIndex].scaleY;
			_width = m_ultrasoundSettings[_dataIndex].probeWidth;
			_corner_T_apex << 1.0, 0.0, 0.0, 0.0,  0.0, 1.0, 0.0, 0.0,  0.0, 0.0, 1.0, 0.0,  0.0, 0.0, 0.0, 1.0;
		}

		Eigen::Transform3f _trans = m_trackingData[i];

		HBITMAP bmp;
		//m_videoSegmentationTools.getFrame(i, bmp);


		//if ( m_videoSegmentationTools.getImageMinMax(bmp, _minX, _maxX, _minY, _maxY) )
		{
		
			// now define the extreme points, which are the apex, the left-most, mid and right-most scanline
			Eigen::Vector4f _tmpPnts[4];
			Eigen::Vector4f _exVec[4];
		_tmpPnts[0] = _corner_T_apex * Eigen::Vector4f(0.0, 0.0, 0.0, 1.0);
			

			_tmpPnts[1] =  _corner_T_apex * Eigen::Vector4f(0.0, _penetrationDepth,  0.0, 1.0); // left-most scanline

			_tmpPnts[2] = _corner_T_apex *  Eigen::Vector4f( _width/2.0, _penetrationDepth ,  0.0, 1.0); // central scanline

			_tmpPnts[3] = _corner_T_apex *  Eigen::Vector4f(_width,  _penetrationDepth ,  0.0, 1.0); // right-most scanline

			// now we will define the rectangular ROI containing the sector

			// left, top
			_exVec[0] << _tmpPnts[1].x(), _tmpPnts[0].y(), 0.0, 1.0;
			// left, bottom
			_exVec[1] << _tmpPnts[1].x(), _tmpPnts[2].y(), 0.0, 1.0;
			// right, bottom
			_exVec[2] << _tmpPnts[3].x(), _tmpPnts[2].y(), 0.0, 1.0;
			// right, top
			_exVec[3] << _tmpPnts[3].x(), _tmpPnts[0].y(), 0.0, 1.0;


			// now we fake the ROI values given the scanlines
			m_roiMinX = _exVec[1].x();
			m_roiMaxX = _exVec[3].x();

			m_roiMinY = _exVec[0].y();
			m_roiMaxY = _exVec[2].y();

			// now we fake scaling in X and Y diretion assuming some resolution of the image
			m_scaleX =150.0 / _penetrationDepth;
			m_scaleY = m_scaleX;

			for (int j=0;j<4;j++)
			{
				Eigen::Vector4f _temp = _trans * _calibrationMatrix *_exVec[j];


				Eigen::Vector3f _tempVec;

				_tempVec << _temp(0), _temp(1), _temp(2);

				_dataPoints.push_back(_tempVec);

				if ( _bbMinX > _tempVec(0) )
					_bbMinX = _tempVec(0);
				if ( _bbMinY > _tempVec(1) )
					_bbMinY = _tempVec(1);
				if ( _bbMinZ > _tempVec(2) )
					_bbMinZ = _tempVec(2);

				if ( _bbMaxX < _tempVec(0) )
					_bbMaxX = _tempVec(0);
				if ( _bbMaxY < _tempVec(1) )
					_bbMaxY = _tempVec(1);
				if ( _bbMaxZ < _tempVec(2) )
					_bbMaxZ = _tempVec(2);
			}
			//char buffer[200];
			//sprintf(buffer,"%f %f %f / %f %f %f",_minVec(0),_minVec(1),_minVec(2),_maxVec(0),_maxVec(1),_maxVec(2));
			//sprintf(buffer,"%d / %d %d / %d %d / %f %f %f / %f %f %f",i,_minX,_minY,_maxX,_maxY,_bbMinX,_bbMinY,_bbMinZ,_bbMaxX,_bbMaxY,_bbMaxZ);
			//ui.textEdit->append(QString(buffer));


			qApp->processEvents();
		}
	}

}



	GeometricTools _geometricTools;

	Eigen::Vector3f _centroid;
	Eigen::Quaternionf _rotation;

	_geometricTools.computeOrientedBoundingBox(_dataPoints, m_extent, m_boundingBoxTransformation);


	//std::cout << _extent.dimX << " / " << _extent.dimY << " / " << _extent.dimZ << std::endl;

	std::stringstream _ss2;
	_ss2 << m_extent.dimX << " / " << m_extent.dimY << " / " << m_extent.dimZ << std::endl;
	_ss2 << m_boundingBoxTransformation.matrix() << std::endl;

	ui.textEdit->append(_ss2.str().c_str());


	//CAMP::Matrix4<float> _image2worldCAMP;
	//TJK_27_07_09: m_dataTools.eigen2CAMP(m_boundingBoxTransformation.matrix(), m_boundingBoxTransformationCAMP);


	/*Eigen::Vector3f _minVec, _maxVec;

	_minVec << _bbMinX, _bbMinY, _bbMinZ;
	_maxVec << _bbMaxX, _bbMaxY, _bbMaxZ;*/

	return true;
}


bool Compounding::determineMinimumBoundingBoxfromMasterFrameRF()
{

	int _masterFrameIndex = atoi(ui.RFMasterFrameEdit2->text(). toLatin1().data());
	m_masterFrameXML = _masterFrameIndex; 
	ui.ImageLabelScrollBar->setValue(_masterFrameIndex);


	// determine the spatial bounding box for the US image

	float _bbMinX, _bbMinY, _bbMinZ;

	_bbMinX = std::numeric_limits<float>::max();
	_bbMinY = std::numeric_limits<float>::max();
	_bbMinZ = std::numeric_limits<float>::max();

	float _bbMaxX, _bbMaxY, _bbMaxZ;

	_bbMaxX = -std::numeric_limits<float>::max();
	_bbMaxY = -std::numeric_limits<float>::max();
	_bbMaxZ = -std::numeric_limits<float>::max();

	std::vector<Eigen::Vector3f> _dataPoints;

	{
		std::stringstream _ss2;
		_ss2 << m_trackingData.size();

		ui.textEdit->append(_ss2.str().c_str());
	}
	int _size = m_trackingData.size();
	int _dataIndex = 0;
	Eigen::Matrix4f _calibrationMatrix = m_calibrationMatrixVector[_dataIndex];
	Eigen::Matrix3f _TMPtrans;
	float _penetrationDepth = m_ultrasoundSettings[_dataIndex].penetrationDepth;
	float _angle = m_ultrasoundSettings[_dataIndex].probeAngle;
	float _apexX = m_ultrasoundSettings[_dataIndex].apexX;
	float _apexY = m_ultrasoundSettings[_dataIndex].apexY;
	float _scaleX = m_ultrasoundSettings[_dataIndex].scaleX;
	float _scaleY = m_ultrasoundSettings[_dataIndex].scaleY;
	float _radius =  (float)(m_ultrasoundSettings[_dataIndex].probeRadius);
	float _width = m_ultrasoundSettings[_dataIndex].probeWidth;
	Eigen::Matrix4f _corner_T_apex;


	_corner_T_apex << 1.0, 0.0, 0.0, _apexX*_scaleX,  0.0, 1.0, 0.0, _apexY*_scaleY,  0.0, 0.0, 1.0, 0.0,  0.0, 0.0, 0.0, 1.0;
	for(int i=0; i<_size;i++)
	{
		float _theta = -_angle/2.0*M_PI/180.0;
		// check if we have to jump to the next index
		if ( i == m_indexVector[_dataIndex])
		{
			_dataIndex++;
			_calibrationMatrix = m_calibrationMatrixVector[_dataIndex];
			_penetrationDepth = m_ultrasoundSettings[_dataIndex].penetrationDepth;
			_angle = m_ultrasoundSettings[_dataIndex].probeAngle;
			_radius =  m_ultrasoundSettings[_dataIndex].probeRadius;
			_apexX = m_ultrasoundSettings[_dataIndex].apexX;
			_apexY = m_ultrasoundSettings[_dataIndex].apexY;
			_scaleX = m_ultrasoundSettings[_dataIndex].scaleX;
			_scaleY = m_ultrasoundSettings[_dataIndex].scaleY;
			_width = m_ultrasoundSettings[_dataIndex].probeWidth;
			_corner_T_apex << 1.0, 0.0, 0.0, _apexX*_scaleX,  0.0, 1.0, 0.0, _apexY*_scaleY,  0.0, 0.0, 1.0, 0.0,  0.0, 0.0, 0.0, 1.0;
		}

		if ( i == _masterFrameIndex )
			_TMPtrans = (m_trackingData[_masterFrameIndex]* _calibrationMatrix).rotation(); 

		Eigen::Transform3f _trans = m_trackingData[i];

		HBITMAP bmp;
		//m_videoSegmentationTools.getFrame(i, bmp);


		//if ( m_videoSegmentationTools.getImageMinMax(bmp, _minX, _maxX, _minY, _maxY) )
		{



			// now define the extreme points, which are the apex, the left-most, mid and right-most scanline
			Eigen::Vector4f _tmpPnts[4];
			Eigen::Vector4f _exVec[4];

				if ( _width < 1.0 ) // curvilinear or phased array probe
				{

			// those points are the extreme points of the sector
			// phased array assumption, scanline begings in the apex
			if ( _radius == 0) {
				_tmpPnts[0] = _corner_T_apex * Eigen::Vector4f(0.0, 0.0, 0.0, 1.0);
			}
			else // curvilinear probe, offset = inner radius of the probe
			{
				Eigen::Vector4f _tempVec = Eigen::Vector4f(_radius * sin(_theta), _radius * cos(_theta),  0.0, 1.0);
				_tmpPnts[0] = _corner_T_apex * _tempVec;
			}
			//_tmpPnts[0] = _corner_T_apex * Eigen::Vector4f(0.0, 0.0, 0.0, 1.0); // apex


			_tmpPnts[1] =  _corner_T_apex * Eigen::Vector4f( ( _radius + _penetrationDepth) * sin(_theta), ( _radius + _penetrationDepth ) * cos(_theta),  0.0, 1.0); // left-most scanline

			_theta = 0.0;
			_tmpPnts[2] = _corner_T_apex *  Eigen::Vector4f( (_radius + _penetrationDepth ) * sin(_theta), ( _radius + _penetrationDepth ) * cos(_theta),  0.0, 1.0); // central scanline

			_theta = _angle/2.0*M_PI/180.0;
			_tmpPnts[3] = _corner_T_apex *  Eigen::Vector4f( (_radius + _penetrationDepth ) * sin(_theta), ( _radius + _penetrationDepth ) * cos(_theta),  0.0, 1.0); // right-most scanline

				}
				else //linear probe
				{
					_tmpPnts[0] = _corner_T_apex * Eigen::Vector4f(0.0, 0.0, 0.0, 1.0);
			

			_tmpPnts[1] =  _corner_T_apex * Eigen::Vector4f(0.0, _penetrationDepth,  0.0, 1.0); // left-most scanline

			_tmpPnts[2] = _corner_T_apex *  Eigen::Vector4f( _width/2.0, _penetrationDepth ,  0.0, 1.0); // central scanline

			_tmpPnts[3] = _corner_T_apex *  Eigen::Vector4f(_width,  _penetrationDepth ,  0.0, 1.0); // right-most scanline
				}
			// now we will define the rectangular ROI containing the sector

			// left, top
			_exVec[0] << _tmpPnts[1].x(), _tmpPnts[0].y(), 0.0, 1.0;
			// left, bottom
			_exVec[1] << _tmpPnts[1].x(), _tmpPnts[2].y(), 0.0, 1.0;
			// right, bottom
			_exVec[2] << _tmpPnts[3].x(), _tmpPnts[2].y(), 0.0, 1.0;
			// right, top
			_exVec[3] << _tmpPnts[3].x(), _tmpPnts[0].y(), 0.0, 1.0;


			// now we fake the ROI values given the scanlines
			m_roiMinX = _exVec[1].x();
			m_roiMaxX = _exVec[3].x();

			m_roiMinY = _exVec[0].y();
			m_roiMaxY = _exVec[2].y();

			// now we fake scaling in X and Y diretion assuming some resolution of the image
			m_scaleX =150.0 / _penetrationDepth;
			m_scaleY = m_scaleX;

			for (int j=0;j<4;j++)
			{
				Eigen::Vector4f _temp = _trans * _calibrationMatrix *_exVec[j];


				Eigen::Vector3f _tempVec;

				_tempVec << _temp(0), _temp(1), _temp(2);

				_dataPoints.push_back(_tempVec);

				if ( _bbMinX > _tempVec(0) )
					_bbMinX = _tempVec(0);
				if ( _bbMinY > _tempVec(1) )
					_bbMinY = _tempVec(1);
				if ( _bbMinZ > _tempVec(2) )
					_bbMinZ = _tempVec(2);

				if ( _bbMaxX < _tempVec(0) )
					_bbMaxX = _tempVec(0);
				if ( _bbMaxY < _tempVec(1) )
					_bbMaxY = _tempVec(1);
				if ( _bbMaxZ < _tempVec(2) )
					_bbMaxZ = _tempVec(2);
			}
			//char buffer[200];
			//sprintf(buffer,"%f %f %f / %f %f %f",_minVec(0),_minVec(1),_minVec(2),_maxVec(0),_maxVec(1),_maxVec(2));
			//sprintf(buffer,"%d / %d %d / %d %d / %f %f %f / %f %f %f",i,_minX,_minY,_maxX,_maxY,_bbMinX,_bbMinY,_bbMinZ,_bbMaxX,_bbMaxY,_bbMaxZ);
			//ui.textEdit->append(QString(buffer));


			qApp->processEvents();
		}
	}


	GeometricTools _geometricTools;

	Eigen::Vector3f _centroid;
	Eigen::Quaternionf _rotation;
	Eigen::Matrix3f _rot = _TMPtrans;
	//_rot.setIdentity();

	_geometricTools.computeManualBoundingBox(_rot, _dataPoints, m_extent, m_boundingBoxTransformation);


	//std::cout << _extent.dimX << " / " << _extent.dimY << " / " << _extent.dimZ << std::endl;

	std::stringstream _ss2;
	_ss2 << m_extent.dimX << " / " << m_extent.dimY << " / " << m_extent.dimZ << std::endl;
	_ss2 << m_boundingBoxTransformation.matrix() << std::endl;

	ui.textEdit->append(_ss2.str().c_str());


	//CAMP::Matrix4<float> _image2worldCAMP;
	//TJK_27_07_09: m_dataTools.eigen2CAMP(m_boundingBoxTransformation.matrix(), m_boundingBoxTransformationCAMP);


	/*Eigen::Vector3f _minVec, _maxVec;

	_minVec << _bbMinX, _bbMinY, _bbMinZ;
	_maxVec << _bbMaxX, _bbMaxY, _bbMaxZ;*/

	return true;



	/*
	int _masterFrameIndex = atoi(ui.RFMasterFrameEdit2->text(). toLatin1().data());
	ui.ImageLabelScrollBar->setValue(_masterFrameIndex);



	// determine the spatial bounding box for the US image

	float _bbMinX, _bbMinY, _bbMinZ;

	_bbMinX = std::numeric_limits<float>::max();
	_bbMinY = std::numeric_limits<float>::max();
	_bbMinZ = std::numeric_limits<float>::max();

	float _bbMaxX, _bbMaxY, _bbMaxZ;

	_bbMaxX = -std::numeric_limits<float>::max();
	_bbMaxY = -std::numeric_limits<float>::max();
	_bbMaxZ = -std::numeric_limits<float>::max();

	std::vector<Eigen::Vector3f> _dataPoints;

	{
		std::stringstream _ss2;
		_ss2 << m_trackingData.size();

		ui.textEdit->append(_ss2.str().c_str());
	}
	int _size = m_trackingData.size();
	int _dataIndex = 0;
	Eigen::Matrix4f _calibrationMatrix = m_calibrationMatrixVector[_dataIndex];
	Eigen::Matrix4f _TMPcalibrationMatrix = m_calibrationMatrixVector[_dataIndex];
	float _penetrationDepth = m_ultrasoundSettings[_dataIndex].penetrationDepth;
	float _angle = m_ultrasoundSettings[_dataIndex].probeAngle;
	float _apexX = m_ultrasoundSettings[_dataIndex].apexX;
	float _apexY = m_ultrasoundSettings[_dataIndex].apexY;
	float _scaleX = m_ultrasoundSettings[_dataIndex].scaleX;
	float _scaleY = m_ultrasoundSettings[_dataIndex].scaleY;
	float _radius =  (float)(m_ultrasoundSettings[_dataIndex].probeRadius);
	Eigen::Matrix4f _corner_T_apex;
	_corner_T_apex << 1.0, 0.0, 0.0, _apexX*_scaleX,  0.0, 1.0, 0.0, _apexY*_scaleY,  0.0, 0.0, 1.0, 0.0,  0.0, 0.0, 0.0, 1.0;

	_TMPcalibrationMatrix = _TMPcalibrationMatrix *_corner_T_apex;

	for(int i=0; i<_masterFrameIndex;i++)
	{
		float _theta = -_angle/2.0*M_PI/180.0;
		// check if we have to jump to the next index
		if ( i == m_indexVector[_dataIndex])
		{
			_dataIndex++;
			_TMPcalibrationMatrix = m_calibrationMatrixVector[_dataIndex];
			_TMPcalibrationMatrix = _TMPcalibrationMatrix *_corner_T_apex;
			
		}
	}

	_TMPcalibrationMatrix.Identity();

	 _dataIndex = 0;

	Eigen::Matrix4f _referenceINV = (m_trackingData[_masterFrameIndex].matrix() * _TMPcalibrationMatrix ).inverse();
	for(int i=0; i<_size;i++)
	{
		float _theta = -_angle/2.0*M_PI/180.0;
		// check if we have to jump to the next index
		if ( i == m_indexVector[_dataIndex])
		{
			_dataIndex++;
			_calibrationMatrix = m_calibrationMatrixVector[_dataIndex];
			_penetrationDepth = m_ultrasoundSettings[_dataIndex].penetrationDepth;
			_angle = m_ultrasoundSettings[_dataIndex].probeAngle;
			_radius =  m_ultrasoundSettings[_dataIndex].probeRadius;
			_apexX = m_ultrasoundSettings[_dataIndex].apexX;
			_apexY = m_ultrasoundSettings[_dataIndex].apexY;
			_scaleX = m_ultrasoundSettings[_dataIndex].scaleX;
			_scaleY = m_ultrasoundSettings[_dataIndex].scaleY;
			_corner_T_apex << 1.0, 0.0, 0.0, _apexX*_scaleX,  0.0, 1.0, 0.0, _apexY*_scaleY,  0.0, 0.0, 1.0, 0.0,  0.0, 0.0, 0.0, 1.0;
		}

		Eigen::Transform3f _trans = m_trackingData[i];

		HBITMAP bmp;
		//m_videoSegmentationTools.getFrame(i, bmp);


		//if ( m_videoSegmentationTools.getImageMinMax(bmp, _minX, _maxX, _minY, _maxY) )
		{



			// now define the extreme points, which are the apex, the left-most, mid and right-most scanline
			Eigen::Vector4f _tmpPnts[4];
			Eigen::Vector4f _exVec[4];

			// those points are the extreme points of the sector
			// phased array assumption, scanline begings in the apex
			if ( _radius == 0) {
				_tmpPnts[0] = _corner_T_apex * Eigen::Vector4f(0.0, 0.0, 0.0, 1.0);
			}
			else // curvilinear probe, offset = inner radius of the probe
			{
				Eigen::Vector4f _tempVec = Eigen::Vector4f(_radius * sin(_theta), _radius * cos(_theta),  0.0, 1.0);
				_tmpPnts[0] = _corner_T_apex * _tempVec;
			}
			//_tmpPnts[0] = _corner_T_apex * Eigen::Vector4f(0.0, 0.0, 0.0, 1.0); // apex


			_tmpPnts[1] =  _corner_T_apex * Eigen::Vector4f( ( _radius + _penetrationDepth) * sin(_theta), ( _radius + _penetrationDepth ) * cos(_theta),  0.0, 1.0); // left-most scanline

			_theta = 0.0;
			_tmpPnts[2] = _corner_T_apex *  Eigen::Vector4f( (_radius + _penetrationDepth ) * sin(_theta), ( _radius + _penetrationDepth ) * cos(_theta),  0.0, 1.0); // central scanline

			_theta = _angle/2.0*M_PI/180.0;
			_tmpPnts[3] = _corner_T_apex *  Eigen::Vector4f( (_radius + _penetrationDepth ) * sin(_theta), ( _radius + _penetrationDepth ) * cos(_theta),  0.0, 1.0); // right-most scanline

			// now we will define the rectangular ROI containing the sector

			// left, top
			_exVec[0] << _tmpPnts[1].x(), _tmpPnts[0].y(), 0.0, 1.0;
			// left, bottom
			_exVec[1] << _tmpPnts[1].x(), _tmpPnts[2].y(), 0.0, 1.0;
			// right, bottom
			_exVec[2] << _tmpPnts[3].x(), _tmpPnts[2].y(), 0.0, 1.0;
			// right, top
			_exVec[3] << _tmpPnts[3].x(), _tmpPnts[0].y(), 0.0, 1.0;


			// now we fake the ROI values given the scanlines
			m_roiMinX = _exVec[1].x();
			m_roiMaxX = _exVec[3].x();

			m_roiMinY = _exVec[0].y();
			m_roiMaxY = _exVec[2].y();

			// now we fake scaling in X and Y diretion assuming some resolution of the image
			m_scaleX = 800.0 / _penetrationDepth;
			m_scaleY = m_scaleX;

			for (int j=0;j<4;j++)
			{
				Eigen::Vector4f _temp = _trans * _calibrationMatrix *_exVec[j];


				Eigen::Vector3f _tempVec;

				_tempVec << _temp(0), _temp(1), _temp(2);

				_dataPoints.push_back(_tempVec);

				if ( _bbMinX > _tempVec(0) )
					_bbMinX = _tempVec(0);
				if ( _bbMinY > _tempVec(1) )
					_bbMinY = _tempVec(1);
				if ( _bbMinZ > _tempVec(2) )
					_bbMinZ = _tempVec(2);

				if ( _bbMaxX < _tempVec(0) )
					_bbMaxX = _tempVec(0);
				if ( _bbMaxY < _tempVec(1) )
					_bbMaxY = _tempVec(1);
				if ( _bbMaxZ < _tempVec(2) )
					_bbMaxZ = _tempVec(2);
			}
			//char buffer[200];
			//sprintf(buffer,"%f %f %f / %f %f %f",_minVec(0),_minVec(1),_minVec(2),_maxVec(0),_maxVec(1),_maxVec(2));
			//sprintf(buffer,"%d / %d %d / %d %d / %f %f %f / %f %f %f",i,_minX,_minY,_maxX,_maxY,_bbMinX,_bbMinY,_bbMinZ,_bbMaxX,_bbMaxY,_bbMaxZ);
			//ui.textEdit->append(QString(buffer));


			qApp->processEvents();
		}
	}


	Eigen::Vector3f _originInWorld;
	Eigen::Vector3f _centroid;
	Eigen::Matrix3f _tmpRotationMatrix;
	_tmpRotationMatrix.setIdentity();

	//m_boundingBoxTransformation)
	GeometricTools _geometricTools;
	Eigen::Transform3f _tempTrans(_referenceINV);
	_geometricTools.computeCentroid(_dataPoints, _tempTrans.rotation(), _centroid, _originInWorld, m_extent);

	m_boundingBoxTransformation.setIdentity();
	m_boundingBoxTransformation.rotate((m_trackingData[_masterFrameIndex]* _TMPcalibrationMatrix).rotation());
	m_boundingBoxTransformation.translate(1.0*_originInWorld);

	std::stringstream _ss2;
	_ss2 << m_extent.dimX << " / " << m_extent.dimY << " / " << m_extent.dimZ << std::endl;
	_ss2 << m_boundingBoxTransformation.matrix() << std::endl;

	ui.textEdit->append(_ss2.str().c_str());


	return true;
	*/

/*
int _masterFrameIndex = atoi(ui.RFMasterFrameEdit2->text(). toLatin1().data());
	ui.ImageLabelScrollBar->setValue(_masterFrameIndex);
// determine the spatial bounding box for the US image

	float _bbMinX, _bbMinY, _bbMinZ;

	_bbMinX = std::numeric_limits<float>::max();
	_bbMinY = std::numeric_limits<float>::max();
	_bbMinZ = std::numeric_limits<float>::max();

	float _bbMaxX, _bbMaxY, _bbMaxZ;

	_bbMaxX = -std::numeric_limits<float>::max();
	_bbMaxY = -std::numeric_limits<float>::max();
	_bbMaxZ = -std::numeric_limits<float>::max();

	std::vector<Eigen::Vector3f> _dataPoints;

	{
		std::stringstream _ss2;
		_ss2 << m_trackingData.size();

		ui.textEdit->append(_ss2.str().c_str());
	}
	int _size = m_trackingData.size();
	int _dataIndex = 0;
	Eigen::Matrix4f _calibrationMatrix = m_calibrationMatrixVector[_dataIndex];
	Eigen::Matrix4f _TMPcalibrationMatrix = m_calibrationMatrixVector[_dataIndex];
	float _penetrationDepth = m_ultrasoundSettings[_dataIndex].penetrationDepth;
	float _angle = m_ultrasoundSettings[_dataIndex].probeAngle;
	float _apexX = m_ultrasoundSettings[_dataIndex].apexX;
	float _apexY = m_ultrasoundSettings[_dataIndex].apexY;
	float _scaleX = m_ultrasoundSettings[_dataIndex].scaleX;
	float _scaleY = m_ultrasoundSettings[_dataIndex].scaleY;
	float _radius =  (float)(m_ultrasoundSettings[_dataIndex].probeRadius);


	


	Eigen::Matrix4f _corner_T_apex;
	_corner_T_apex << 1.0, 0.0, 0.0, _apexX*_scaleX,  0.0, 1.0, 0.0, _apexY*_scaleY,  0.0, 0.0, 1.0, 0.0,  0.0, 0.0, 0.0, 1.0;




	for(int i=0; i<_masterFrameIndex;i++)
	{
		
		// check if we have to jump to the next index
		if ( i == m_indexVector[_dataIndex])
		{
			_dataIndex++;
			_TMPcalibrationMatrix = m_calibrationMatrixVector[_dataIndex];
			_TMPcalibrationMatrix = _TMPcalibrationMatrix;// *_corner_T_apex;
			
		}
	}

	Eigen::Matrix4f _alignTrans;

	/*
	_alignTrans << 0,0,1,0,
				  0,1,0,0,
				  1,0,0,0,
				  0,0,0,1;
				  // crash
	*/

	/*_alignTrans << 0,0,1,0,
				  0,1,0,0,
				  1,0,0,0,
				  0,0,0,1;

	// nonsense

	*/

	/*_alignTrans << 0,1,0,0,
				  1,0,0,0,
				  0,0,1,0,
				  0,0,0,1;
				  
				  // crash
				  */

	/*_alignTrans << 1,0,0,0,
				  0,1,0,0,
				  0,0,1,0,
				  0,0,0,1;*/


	/*_alignTrans << 1,0,0,0,
				  0,0,1,0,
				  0,1,0,0,
				  0,0,0,1;
				  // post: nonsense
				  */

	/*_alignTrans << 0,0,1,0,
				  1,0,0,0,
				  0,1,0,0,
				  0,0,0,1;
	// post: unten rechts okay
	// pre: schlecht
	*/


	/*_alignTrans << 0,1,0,0,
				  0,0,1,0,
				  1,0,0,0,
				  0,0,0,1;

				 // post: oben rechts zur hälfte okay
				 // pre: oben rechts okay
				  */
	


	

	Eigen::Matrix4f _referenceINV = ( m_trackingData[_masterFrameIndex].matrix() ).inverse();
	Eigen::Matrix4f _reference = ( m_trackingData[_masterFrameIndex].matrix()  );
	/*_referenceINV.setIdentity();
	_referenceINV << 0,0,1,0,
					1,0,0,0,
					0,1,0,0,
					0,0,0,1;*/
	_dataIndex = 0;



	for(int i=0; i<_size;i++)
	{
		float _theta = -_angle/2.0*M_PI/180.0;
		// check if we have to jump to the next index
		if ( i == m_indexVector[_dataIndex])
		{
			_dataIndex++;
			_calibrationMatrix = m_calibrationMatrixVector[_dataIndex];
			_penetrationDepth = m_ultrasoundSettings[_dataIndex].penetrationDepth;
			_angle = m_ultrasoundSettings[_dataIndex].probeAngle;
			_radius =  m_ultrasoundSettings[_dataIndex].probeRadius;
			_apexX = m_ultrasoundSettings[_dataIndex].apexX;
			_apexY = m_ultrasoundSettings[_dataIndex].apexY;
			_scaleX = m_ultrasoundSettings[_dataIndex].scaleX;
			_scaleY = m_ultrasoundSettings[_dataIndex].scaleY;
			_corner_T_apex << 1.0, 0.0, 0.0, _apexX*_scaleX,  0.0, 1.0, 0.0, _apexY*_scaleY,  0.0, 0.0, 1.0, 0.0,  0.0, 0.0, 0.0, 1.0;
		}

		Eigen::Transform3f _trans = m_trackingData[i];

		HBITMAP bmp;
		//m_videoSegmentationTools.getFrame(i, bmp);


		//if ( m_videoSegmentationTools.getImageMinMax(bmp, _minX, _maxX, _minY, _maxY) )
		{



			// now define the extreme points, which are the apex, the left-most, mid and right-most scanline
			Eigen::Vector4f _tmpPnts[4];
			Eigen::Vector4f _exVec[4];

			// those points are the extreme points of the sector
			// phased array assumption, scanline begings in the apex
			if ( _radius == 0) {
				_tmpPnts[0] = _corner_T_apex * Eigen::Vector4f(0.0, 0.0, 0.0, 1.0);
			}
			else // curvilinear probe, offset = inner radius of the probe
			{
				Eigen::Vector4f _tempVec = Eigen::Vector4f(_radius * sin(_theta), _radius * cos(_theta),  0.0, 1.0);
				_tmpPnts[0] = _corner_T_apex * _tempVec;
			}
			//_tmpPnts[0] = _corner_T_apex * Eigen::Vector4f(0.0, 0.0, 0.0, 1.0); // apex


			_tmpPnts[1] =  _corner_T_apex * Eigen::Vector4f( ( _radius + _penetrationDepth) * sin(_theta), ( _radius + _penetrationDepth ) * cos(_theta),  0.0, 1.0); // left-most scanline

			_theta = 0.0;
			_tmpPnts[2] = _corner_T_apex *  Eigen::Vector4f( (_radius + _penetrationDepth ) * sin(_theta), ( _radius + _penetrationDepth ) * cos(_theta),  0.0, 1.0); // central scanline

			_theta = _angle/2.0*M_PI/180.0;
			_tmpPnts[3] = _corner_T_apex *  Eigen::Vector4f( (_radius + _penetrationDepth ) * sin(_theta), ( _radius + _penetrationDepth ) * cos(_theta),  0.0, 1.0); // right-most scanline

			// now we will define the rectangular ROI containing the sector

			// left, top
			_exVec[0] << _tmpPnts[1].x(), _tmpPnts[0].y(), 0.0, 1.0;
			// left, bottom
			_exVec[1] << _tmpPnts[1].x(), _tmpPnts[2].y(), 0.0, 1.0;
			// right, bottom
			_exVec[2] << _tmpPnts[3].x(), _tmpPnts[2].y(), 0.0, 1.0;
			// right, top
			_exVec[3] << _tmpPnts[3].x(), _tmpPnts[0].y(), 0.0, 1.0;


			// now we fake the ROI values given the scanlines
			m_roiMinX = _exVec[1].x();
			m_roiMaxX = _exVec[3].x();

			m_roiMinY = _exVec[0].y();
			m_roiMaxY = _exVec[2].y();

			// now we fake scaling in X and Y diretion assuming some resolution of the image
			m_scaleX =150.0 / _penetrationDepth;
			m_scaleY = m_scaleX;

			for (int j=0;j<4;j++)
			{
				Eigen::Vector4f _temp = _trans * _calibrationMatrix *_exVec[j];


				Eigen::Vector3f _tempVec;

				_tempVec << _temp(0), _temp(1), _temp(2);

				_dataPoints.push_back(_tempVec);

				if ( _bbMinX > _tempVec(0) )
					_bbMinX = _tempVec(0);
				if ( _bbMinY > _tempVec(1) )
					_bbMinY = _tempVec(1);
				if ( _bbMinZ > _tempVec(2) )
					_bbMinZ = _tempVec(2);

				if ( _bbMaxX < _tempVec(0) )
					_bbMaxX = _tempVec(0);
				if ( _bbMaxY < _tempVec(1) )
					_bbMaxY = _tempVec(1);
				if ( _bbMaxZ < _tempVec(2) )
					_bbMaxZ = _tempVec(2);
			}
			//char buffer[200];
			//sprintf(buffer,"%f %f %f / %f %f %f",_minVec(0),_minVec(1),_minVec(2),_maxVec(0),_maxVec(1),_maxVec(2));
			//sprintf(buffer,"%d / %d %d / %d %d / %f %f %f / %f %f %f",i,_minX,_minY,_maxX,_maxY,_bbMinX,_bbMinY,_bbMinZ,_bbMaxX,_bbMaxY,_bbMaxZ);
			//ui.textEdit->append(QString(buffer));


			qApp->processEvents();
		}
	}

	/*
	GeometricTools _geometricTools;

	Eigen::Vector3f _centroid;
	Eigen::Quaternionf _rotation;

	_geometricTools.computeOrientedBoundingBox(_dataPoints, m_extent, m_boundingBoxTransformation);



	//std::cout << _extent.dimX << " / " << _extent.dimY << " / " << _extent.dimZ << std::endl;

	std::stringstream _ss2;
	_ss2 << m_extent.dimX << " / " << m_extent.dimY << " / " << m_extent.dimZ << std::endl;
	_ss2 << m_boundingBoxTransformation.matrix() << std::endl;

	ui.textEdit->append(_ss2.str().c_str());


	return true;*/
	/*
	Eigen::Vector3f _originInWorld;
	Eigen::Vector3f _centroid;
	Eigen::Matrix3f _tmpRotationMatrix;
	_tmpRotationMatrix.setIdentity();

	//m_boundingBoxTransformation)
	GeometricTools _geometricTools;
	Eigen::Transform3f _tempTrans(_reference);
	_geometricTools.computeCentroid(_dataPoints, _tempTrans.rotation(), _centroid, _originInWorld, m_extent);

	m_boundingBoxTransformation.setIdentity();
	Eigen::Matrix3f _rot = _referenceINV.block(0,0,3,3);
	m_boundingBoxTransformation.rotate(_rot);
	m_boundingBoxTransformation.translate(1.0*_originInWorld);


	// origin of trans is the lowest point
	// now we have to transform by transforming to another corner by reversing the Y-direction

	Eigen::Matrix4f _translate;
	Eigen::Matrix4f _mirror;

	_translate << 1.0, 0.0, 0.0, 0.0,		0.0, 1.0, 0.0, m_extent.dimY,	0.0, 0.0, 1.0, 0.0,		0.0, 0.0, 0.0, 1.0;
	_mirror << 1.0, 0.0, 0.0, 0.0,		0.0, -1.0, 0.0, 0.0,	0.0, 0.0, 1.0, 0.0,		0.0, 0.0, 0.0, 1.0;
	

	// bug inside :(
	m_boundingBoxTransformation =   m_boundingBoxTransformation.matrix() * _translate * _mirror;


	std::stringstream _ss2;
	_ss2 << m_extent.dimX << " / " << m_extent.dimY << " / " << m_extent.dimZ << std::endl;
	_ss2 << m_boundingBoxTransformation.matrix() << std::endl;

	ui.textEdit->append(_ss2.str().c_str());


	return true;*/
}

bool Compounding::determineMinimumBoundingBoxFromBMode()
{
	// determine the spatial bounding box for the US image

	float _bbMinX, _bbMinY, _bbMinZ;

	_bbMinX = std::numeric_limits<float>::max();
	_bbMinY = std::numeric_limits<float>::max();
	_bbMinZ = std::numeric_limits<float>::max();

	float _bbMaxX, _bbMaxY, _bbMaxZ;

	_bbMaxX = -std::numeric_limits<float>::max();
	_bbMaxY = -std::numeric_limits<float>::max();
	_bbMaxZ = -std::numeric_limits<float>::max();

	std::vector<Eigen::Vector3f> _dataPoints;

	{
		std::stringstream _ss2;
		_ss2 << m_trackingData.size();

		ui.textEdit->append(_ss2.str().c_str());
	}
	int _size = m_trackingData.size();
	for(int i=0; i<_size;i++)
	{
		Eigen::Transform3f _trans = m_trackingData[i];

		HBITMAP bmp;
		//m_videoSegmentationTools.getFrame(i, bmp);

		int _minX=m_imageMinX, _minY=m_imageMinY, _maxX=m_imageMaxX, _maxY=m_imageMaxY;

		//if ( m_videoSegmentationTools.getImageMinMax(bmp, _minX, _maxX, _minY, _maxY) )
		{

			Eigen::Vector4f _exVec[4];
			_exVec[0] << m_roiMinX * m_scaleX, m_roiMinY * m_scaleY, 0, 1;
			_exVec[1] << m_roiMaxX * m_scaleX, m_roiMaxY * m_scaleY, 0, 1;
			_exVec[2] << m_roiMaxX * m_scaleX, m_roiMinY * m_scaleY, 0, 1;
			_exVec[3] << m_roiMinX * m_scaleX, m_roiMaxY * m_scaleY, 0, 1;


			for (int j=0;j<4;j++)
			{
				_exVec[j] = _trans * m_matCalibration * _exVec[j];


				Eigen::Vector3f _tempVec;

				_tempVec << _exVec[j](0), _exVec[j](1), _exVec[j](2);

				_dataPoints.push_back(_tempVec);

				if ( _bbMinX > _exVec[j](0) )
					_bbMinX = _exVec[j](0);
				if ( _bbMinY > _exVec[j](1) )
					_bbMinY = _exVec[j](1);
				if ( _bbMinZ > _exVec[j](2) )
					_bbMinZ = _exVec[j](2);

				if ( _bbMaxX < _exVec[j](0) )
					_bbMaxX = _exVec[j](0);
				if ( _bbMaxY < _exVec[j](1) )
					_bbMaxY = _exVec[j](1);
				if ( _bbMaxZ < _exVec[j](2) )
					_bbMaxZ = _exVec[j](2);
			}
			//char buffer[200];
			//sprintf(buffer,"%f %f %f / %f %f %f",_minVec(0),_minVec(1),_minVec(2),_maxVec(0),_maxVec(1),_maxVec(2));
			//sprintf(buffer,"%d / %d %d / %d %d / %f %f %f / %f %f %f",i,_minX,_minY,_maxX,_maxY,_bbMinX,_bbMinY,_bbMinZ,_bbMaxX,_bbMaxY,_bbMaxZ);
			//ui.textEdit->append(QString(buffer));


			qApp->processEvents();
		}
	}


	GeometricTools _geometricTools;

	Eigen::Vector3f _centroid;
	Eigen::Quaternionf _rotation;

	_geometricTools.computeOrientedBoundingBox(_dataPoints, m_extent, m_boundingBoxTransformation);

	//std::cout << _extent.dimX << " / " << _extent.dimY << " / " << _extent.dimZ << std::endl;

	std::stringstream _ss2;
	_ss2 << m_extent.dimX << " / " << m_extent.dimY << " / " << m_extent.dimZ << std::endl;
	_ss2 << m_boundingBoxTransformation.matrix() << std::endl;

	ui.textEdit->append(_ss2.str().c_str());


	//CAMP::Matrix4<float> _image2worldCAMP;
	//TJK_27_07_09: m_dataTools.eigen2CAMP(m_boundingBoxTransformation.matrix(), m_boundingBoxTransformationCAMP);


	/*Eigen::Vector3f _minVec, _maxVec;

	_minVec << _bbMinX, _bbMinY, _bbMinZ;
	_maxVec << _bbMaxX, _bbMaxY, _bbMaxZ;*/

	return true;
}



bool Compounding::determineMinimumBoundingBoxManualFromBMode()
{
	// determine the spatial bounding box for the US image

	float _bbMinX, _bbMinY, _bbMinZ;

	_bbMinX = std::numeric_limits<float>::max();
	_bbMinY = std::numeric_limits<float>::max();
	_bbMinZ = std::numeric_limits<float>::max();

	float _bbMaxX, _bbMaxY, _bbMaxZ;

	_bbMaxX = -std::numeric_limits<float>::max();
	_bbMaxY = -std::numeric_limits<float>::max();
	_bbMaxZ = -std::numeric_limits<float>::max();

	std::vector<Eigen::Vector3f> _dataPoints;

	{
		std::stringstream _ss2;
		_ss2 << m_trackingData.size();

		ui.textEdit->append(_ss2.str().c_str());
	}
	int _size = m_trackingData.size();
	for(int i=0; i<_size;i++)
	{
		Eigen::Transform3f _trans = m_trackingData[i];

		HBITMAP bmp;
		//m_videoSegmentationTools.getFrame(i, bmp);

		int _minX=m_imageMinX, _minY=m_imageMinY, _maxX=m_imageMaxX, _maxY=m_imageMaxY;

		//if ( m_videoSegmentationTools.getImageMinMax(bmp, _minX, _maxX, _minY, _maxY) )
		{

			Eigen::Vector4f _exVec[4];
			_exVec[0] << m_roiMinX * m_scaleX, m_roiMinY * m_scaleY, 0, 1;
			_exVec[1] << m_roiMaxX * m_scaleX, m_roiMaxY * m_scaleY, 0, 1;
			_exVec[2] << m_roiMaxX * m_scaleX, m_roiMinY * m_scaleY, 0, 1;
			_exVec[3] << m_roiMinX * m_scaleX, m_roiMaxY * m_scaleY, 0, 1;


			for (int j=0;j<4;j++)
			{
				_exVec[j] = _trans * m_matCalibration * _exVec[j];


				Eigen::Vector3f _tempVec;

				_tempVec << _exVec[j](0), _exVec[j](1), _exVec[j](2);

				_dataPoints.push_back(_tempVec);

				if ( _bbMinX > _exVec[j](0) )
					_bbMinX = _exVec[j](0);
				if ( _bbMinY > _exVec[j](1) )
					_bbMinY = _exVec[j](1);
				if ( _bbMinZ > _exVec[j](2) )
					_bbMinZ = _exVec[j](2);

				if ( _bbMaxX < _exVec[j](0) )
					_bbMaxX = _exVec[j](0);
				if ( _bbMaxY < _exVec[j](1) )
					_bbMaxY = _exVec[j](1);
				if ( _bbMaxZ < _exVec[j](2) )
					_bbMaxZ = _exVec[j](2);
			}
			//char buffer[200];
			//sprintf(buffer,"%f %f %f / %f %f %f",_minVec(0),_minVec(1),_minVec(2),_maxVec(0),_maxVec(1),_maxVec(2));
			//sprintf(buffer,"%d / %d %d / %d %d / %f %f %f / %f %f %f",i,_minX,_minY,_maxX,_maxY,_bbMinX,_bbMinY,_bbMinZ,_bbMaxX,_bbMaxY,_bbMaxZ);
			//ui.textEdit->append(QString(buffer));


			qApp->processEvents();
		}
	}


	GeometricTools _geometricTools;

	Eigen::Vector3f _centroid;
	Eigen::Quaternionf _rotation;

	//_geometricTools.computeOrientedBoundingBox(_dataPoints, m_extent, m_boundingBoxTransformation);

		// get the transformation matrix from the text input field
	{
		std::string _tmpStr = ui.manualMatrixBMode->toPlainText (). toLatin1().data();
		Eigen::Matrix4f _tmpMatrix;
		XMLTools _xmlTools;
		_xmlTools.string2matrix(_tmpMatrix, _tmpStr);

		m_boundingBoxTransformation.matrix() = _tmpMatrix;
		m_extent.dimX = atoi(ui.dimX_2->text().toLatin1().data());
		m_extent.dimY = atoi(ui.dimY_2->text().toLatin1().data());
		m_extent.dimZ = atoi(ui.dimZ_2->text().toLatin1().data());
	}
//	_geometricTools.determineExtent(_dataPoints, m_boundingBoxTransformation.matrix(), m_extent);



	//std::cout << _extent.dimX << " / " << _extent.dimY << " / " << _extent.dimZ << std::endl;

	std::stringstream _ss2;
	_ss2 << m_extent.dimX << " / " << m_extent.dimY << " / " << m_extent.dimZ << std::endl;
	_ss2 << m_boundingBoxTransformation.matrix() << std::endl;

	ui.textEdit->append(_ss2.str().c_str());


	//CAMP::Matrix4<float> _image2worldCAMP;
	//TJK_27_07_09: m_dataTools.eigen2CAMP(m_boundingBoxTransformation.matrix(), m_boundingBoxTransformationCAMP);


	/*Eigen::Vector3f _minVec, _maxVec;

	_minVec << _bbMinX, _bbMinY, _bbMinZ;
	_maxVec << _bbMaxX, _bbMaxY, _bbMaxZ;*/

	return true;
}


bool Compounding::determineVoxelSize()
{
	// isotrophic voxel size 0.5 mm per voxel
	m_physicalX = m_voxelSize; 
	m_physicalY = m_voxelSize; 
	m_physicalZ = m_voxelSize;



	// determine the number of voxels
	m_voxelsX = static_cast<int>(ceil(m_extent.dimX/m_physicalX));
	m_voxelsY = static_cast<int>(ceil(m_extent.dimY/m_physicalY));
	m_voxelsZ = static_cast<int>(ceil(m_extent.dimZ/m_physicalZ));

	std::stringstream ss;

	ss << "Voxel Size: " << m_voxelSize;
	ui.voxelSizeLabel->setText(ss.str().c_str());


	std::stringstream ss2;
	ss2 << "Voxels: " << m_voxelsX*m_voxelsY*m_voxelsZ << " / " << m_voxelsX << " x " << m_voxelsY << " x " << m_voxelsZ;
	ui.VoxelsNumLabel->setText(ss2.str().c_str());

	return true;
}

bool Compounding::determineVoxelSizeRF()
{
	

	if ( ui.ManualRFRadio->isChecked()) // if manual is selected
	{
	m_physicalX = atof(ui.phyX->text().toLatin1().data());
		m_physicalY = atof(ui.phyY->text().toLatin1().data());
		m_physicalZ = atof(ui.phyZ->text().toLatin1().data());
	}
	else // automatically determined
	{
		// isotrophic voxel size 0.5 mm per voxel
	m_physicalX = m_voxelSize; 
	m_physicalY = m_voxelSize; 
	m_physicalZ = m_voxelSize;
	}



	// determine the number of voxels
	m_voxelsX = static_cast<int>(ceil(m_extent.dimX/m_physicalX));
	m_voxelsY = static_cast<int>(ceil(m_extent.dimY/m_physicalY));
	m_voxelsZ = static_cast<int>(ceil(m_extent.dimZ/m_physicalZ));

	// now determine if the voxelSize can be divided by 2 without remainder, otherwise adapt

	if ( m_voxelsX%2 != 0 ) m_voxelsX++;
	if ( m_voxelsY%2 != 0 ) m_voxelsY++;
	if ( m_voxelsZ%2 != 0 ) m_voxelsZ++;

	m_extent.dimX = ((float)m_voxelsX)*m_physicalX;
	m_extent.dimY = ((float)m_voxelsY)*m_physicalY;
	m_extent.dimZ = ((float)m_voxelsZ)*m_physicalZ;


	std::stringstream ss;

	ss << "Voxel Size: " << m_voxelSize;
	ui.voxelSizeLabel->setText(ss.str().c_str());


	std::stringstream ss2;
	ss2 << "Voxels: " << m_voxelsX*m_voxelsY*m_voxelsZ;
	ui.VoxelsNumLabel->setText(ss2.str().c_str());

	return true;
}

bool Compounding::compoundVolumeVelocity(float _tau, float _lambda1, float _lambda2, int NUMBER_OF_ITERATIORS)
{
	//int NUMBER_OF_ITERATIORS = 200;
	int _solution = 0;
	int _sumSol = 0;
	int _traversalCount = 0;
	clock_t _startTimer = clock();
	ui.textEdit->append(QString("Compounding started ..."));
	qApp->processEvents();


	int _maxVoxels =  m_voxelsX*m_voxelsY*m_voxelsZ*NUMBER_OF_ITERATIORS;
	int _1percent = (int)_maxVoxels/100.0;
	int _percentCounter = 0;
	int _percent = 0;
	ui.progressBar->setValue(0);
	qApp->processEvents();

	float _maxError =0.0f;
	int _voxelsConsidered = 0;
	// Nonlinear optimization begin
	/*float _lambda1 = 1.0;
	float _lambda2 = 0.1;
	float _tau = 0.0001;*/
	for (int c=0;c<NUMBER_OF_ITERATIORS;c++)
	{


		float _completeError = 0.0;
		float _completeTemporaryError = 0.0;

		//_traversalScheme.init();

		float _maxDeltaNorm = -999999999.9f;
		Eigen::Vector3f _maxDeltaVector;

		int _index = -1;

		//  #pragma omp parallel for
		for(int _z=0;_z<m_voxelsZ;_z++)
			for(int _y=0;_y<m_voxelsY;_y++)
				for(int _x=0;_x<m_voxelsX;_x++)
					//while(_traversalScheme.nextPosition(_x,_y,_z))


				{


					_percentCounter++;

					if ( _percentCounter == _1percent)
					{
						_percent++;
						_percentCounter = 0;
						ui.progressBar->setValue(_percent);
						qApp->processEvents();
					}


					_index++;

					Eigen::VectorXf _velocityError;

					//_deltaVector =  VectorXf::Zero(3);

					Eigen::VectorXf _deltaVectorX;
					if ( m_recon->getProjection(_index, _deltaVectorX)  &&  m_recon->getVelocityFieldValidity(_x,_y,_z) )

					{

						// just some statistics
						if ( c == 0 )
						{
							_voxelsConsidered++;
						}


						m_recon->getApproximationError(_index,_velocityError);
						Eigen::Vector3f _deltaVector(_deltaVectorX[0], _deltaVectorX[1], _deltaVectorX[2]);

						for (int k=0;k<_velocityError.rows();k++)
						{
							_completeError += sqrt(pow(_velocityError.row(k)[0],2.0f));
						}
						//std::cout << "Temp: " << _deltaVector << std::endl;



						const Eigen::Vector3f _laplacian = m_recon->computeLaplacian(_x,_y,_z, m_physicalX, m_physicalY, m_physicalZ);

						const Eigen::Vector3f _divergence = m_recon->computeDivergence(_x,_y,_z, m_physicalX, m_physicalY, m_physicalZ);

						const Eigen::Vector3f _deltaVelocity = _tau*(_deltaVector - _lambda1 * _laplacian - _lambda2 * _divergence);




						/*char c;
						*/
						/*
						std::cout << "_deltaVector: " << _deltaVector.transpose() << std::endl;
						std::cout << "_laplacian: " << _laplacian.transpose() << std::endl;
						std::cout << "_divergence: " << _divergence.transpose() << std::endl;
						/*
						/*
						std::cin >> c;*/

						//std::cout << "Delta: " << _deltaVelocity.norm() << " | " << _deltaVelocity.transpose() << std::endl;


						if ( _maxDeltaNorm < _deltaVelocity.norm() )
						{
							_maxDeltaNorm = _deltaVelocity.norm();
							_maxDeltaVector = _deltaVelocity;
						}

						m_recon->updateDeltaVelocity(_index,-1.0*_deltaVelocity);


						m_recon->getApproximationTemporaryError(_index,_velocityError);

						for (int k=0;k<_velocityError.rows();k++)
						{
							_completeTemporaryError += sqrt(pow(_velocityError.row(k)[0],2.0f));
						}

					}

				}


				std::cout << "Max Delta: " << _maxDeltaNorm << " | " << _maxDeltaVector.transpose() << std::endl;
				if ( fabs(_completeTemporaryError)/*.norm()*/ >  fabs(_completeError)/*.norm()*/ )
				{
					_tau = _tau*0.5f;
					std::cout << "New tau: " << _tau << std::endl;
					std::cout << "Error: " << _completeError/*.norm()*/ << std::endl;
				}
				else
				{
					std::cout << "Update vector field!" << std::endl;
					m_recon->performVectorFieldUpdate();
					std::cout << "Error: " << _completeTemporaryError/*.norm()*/ << std::endl;
				}


	}

	// Nonlinear optimization end


	std::cout << "Max Error: " << _maxError << std::endl;


	Eigen::Vector3f* _velocityField = m_recon->getVelocityField();

	m_velocityVisualization->setVelocityField(_velocityField, m_voxelsX, m_voxelsY, m_voxelsZ, m_physicalX, m_physicalY, m_physicalZ);

	//m_dataTools.saveVectorFieldVTK("test.vti",_velocityField, m_voxelsX, m_voxelsY, m_voxelsZ, m_physicalX, m_physicalY, m_physicalZ);
	//m_velocityVisualization->setPlane(_recon->getNormalPlane(), m_imageWidth, m_imageHeight, m_scaleX, m_scaleY);
	std::cout << _solution << " out of " << _sumSol << " have a valid solution." << std::endl;


	// free the memory


	double _sum = 0;
	double _sumNZ = 0;
	double _NZC = 0;
	for(int g=0;g<m_recon->m_measurementField.size();g++)
	{
		VelocityData *_data = m_recon->m_measurementField[g];

		_sum += static_cast<double>(_data->normalMatrix.rows() );
		if ( _data->normalMatrix.rows()  > 0 )
		{
			_sumNZ += static_cast<double>(_data->normalMatrix.rows() );
			_NZC++;
		}		
		//	std::cout << g << ".   " << _data->normalMatrix.rows() << std::endl;
	}

	_sum /= static_cast<double>(m_recon->m_measurementField.size());

	std::cout << "Avg. number of elements per voxel: " << _sum << std::endl;
	std::cout << "Avg. number of elements per non-zero voxel: " << _sumNZ / _NZC << std::endl;



	{
		std::stringstream _ss;
		_ss << "Voxels updated per iteration: " << _voxelsConsidered;

		ui.textEdit->append(_ss.str().c_str());

		_ss.str("");
	}

	// now compute the time it took for compounding the ultrasound volume

	// the last percentage update might have got lost
	//if ( m_percentage < 100.0)
	//	updatePercentage();

	{
		std::stringstream _ss;
		_ss << "Compounding finished after " << m_dataTools.computeElapsedTime(_startTimer, clock());

		ui.textEdit->append(_ss.str().c_str());


		_ss.str("");
	}

	//_ss << "Voxels visited: " << _traversalCount << " of " << m_voxelsX * m_voxelsY * m_voxelsZ;
	//ui.textEdit->append(_ss.str().c_str());

	qApp->processEvents();

	return true;
}


bool Compounding::computeMeasurementCollection()
{

	clock_t _startTimer = clock();
	ui.textEdit->append(QString("Collecting measurement data started ..."));
	qApp->processEvents();
	// set up the data structure that is going to store the reconstruction, measurement data and so on
	if ( m_recon != NULL ) delete m_recon;

	m_recon = new Reconstruction();

	// store the calibration values (ultrasound calibration)
	m_recon->setCalibration(m_matCalibration, m_scaleX, m_scaleY);
	// generate the vector field to hold the result
	m_recon->generateVelocityField(m_voxelsX, m_voxelsY, m_voxelsZ, m_physicalX, m_physicalY, m_physicalZ, m_boundingBoxTransformation.matrix());
	// compute the ultrasoud planes, where each pixel contains the normal with respect to the apex as well as the local space position
	m_recon->preComputePlanes(m_apexX, m_apexY, m_roiMinX, m_roiMinY, m_roiMaxX, m_roiMaxY);

	m_slicesVisualization->setSpacePlane(m_recon->getSpacePlane(), m_recon->getNormalPlane());


	// Create the MHD image and set its parameters
	/*CAMP::Image *compounding = new CAMP::ImageBase<unsigned char>;
	compounding->init( m_voxelsX, m_voxelsY, m_voxelsZ, 1,1);
	compounding->setPhysicalSize(m_physicalX, m_physicalY, m_physicalZ, 1.0);
	*/
	// determine the maximum dimension of the voxel spacing
	float _maxDim = m_voxelSize;

	// initialize the rotation queue
	RotationQueue<USFrame> _rotQueue(_maxDim,m_extent.dimX+m_extent.dimY+m_extent.dimZ, m_imageCounter);

	Eigen::Matrix4f _scalingMatrix;
	_scalingMatrix <<	m_scaleX,0,0,0, 
		0,m_scaleY,0,0, 
		0,0,1,0, 
		0,0,0,1;

	// now we set up a data structure containing information about the position of the ultrasound planes in space
	for (int i=0; i<m_matrixDataTracking.size();i++)
	{
		USFrame* _tempFrame = new USFrame();
		// the tracking data specifies where the plane lies in space, here the origin is assumed in the upper left corner of the image
		Plane *_plane = new Plane(m_matrixDataTracking[i] /*TJK_22_02_10* m_matCalibration*/);
		_tempFrame->distance = (float)0.0;
		// _matrixDataTracking.at(x) transforms from unit [mm] on the image to world [mm]
		// but we typically need it for pixels, thus we multiply it with the scaling factors

		//TJK_27_07_09: m_dataTools.eigen2CAMP(m_matrixDataTracking[i] * _scalingMatrix, _tempFrame->image2World);
		_tempFrame->image2World = (m_matrixDataTracking[i] /*TJK_22_02_10* m_matCalibration*/ * _scalingMatrix);
		_tempFrame->local2World =  /*TJK_22_02_10m_matrixDataTracking[i]*/ m_trackingData[i];
		// compute all the relevant information describing the ultrasound plane
		//TJK_27_07_09: _tempFrame->world2Image = _tempFrame->image2World.getInverse();
		_tempFrame->world2Image = _tempFrame->image2World.inverse();
		//_tempFrame->image2WorldEigen = _matrixDataTracking.at(i) * _scalingMatrix;
		//_tempFrame->world2ImageEigen = _tempFrame->image2WorldEigen.inverse();

		_tempFrame->vectorIndex = (int)floor((float)i/NUMBER_OF_SLICES); 
		_tempFrame->videoFrame = i % (int)(NUMBER_OF_SLICES);
		_tempFrame->plane = _plane;

		// set the ultrasound frame at initial distance 0.0
		_rotQueue.push(_tempFrame,0.0);

		// DEBUG CODE
		/*if ( i==2 )
		m_debugFrame = _tempFrame;
		*/
	}


	RadiusSelect _radiusSelect(m_voxelSize, std::max(m_scaleX, m_scaleY), 50);

	float _threshold = 2.0*_maxDim;
	float _thresholdSQR = pow(_threshold,2.0f);
	stdext::hash_map<int, std::vector<Segment>*> _hashMap;
	stdext::hash_map<int, std::vector<Segment>*>::iterator _hashMapIterator;


	TraversalScheme _traversalScheme(m_voxelsX, m_voxelsY, m_voxelsZ);
	int _x,_y,_z;
	int _traversalCount = 0;

	// Gaussian backward
	float _sigmaSQR = pow(1.0,2.0);

	Eigen::Vector3f _worldIntersection;







	/*
	_recon->saveNormalPlaneAsMAT();

	std::vector<Eigen::Vector3f> _spacePlane;
	for(int y=0;y<m_roiHeight;y++)
	for(int x=0;x<m_roiWidth;x++)
	{
	const Eigen::Vector3f _temp(x, y, 0.0);
	const Eigen::Vector3f _worldPoint = (m_debugFrame->local2World * _recon->m_normalPlane[y*m_roiWidth+x]);

	_spacePlane.push_back(_worldPoint);
	}

	saveMatlabMatrix3Dim(_spacePlane,m_roiWidth, m_roiHeight, "NormalPlane", "NormalPlane.mat");
	*/
	GraphicTools _graphicTools;

	// traverse through the volume

	ui.progressBar->setRange(0,100);
	// as to avoid too frequen updates of the GUI we do only updates when values change
	int _maxVoxels =  m_voxelsX*m_voxelsY*m_voxelsZ;
	int _1percent = (int)_maxVoxels/100.0;
	int _percentCounter = 0;
	int _percent = 0;
	ui.progressBar->setValue(_percent);
	qApp->processEvents();

	//#pragma omp parallel for

	//while(_traversalScheme.nextPosition(_x,_y,_z)) // now we will collect the measurments for each voxel

	for(int _z=0;_z<m_voxelsZ;_z++)
		for(int _y=0;_y<m_voxelsY;_y++)
			for(int _x=0;_x<m_voxelsX;_x++)
			{
				int _lastIndex = -1;
				_percentCounter++;
				_traversalCount++;

				if ( _percentCounter == _1percent)
				{
					_percent++;
					_percentCounter = 0;
					ui.progressBar->setValue(_percent);
					qApp->processEvents();
				}


				// compute the world coordinates from the voxel coordinates
				const Eigen::Vector3f _local(static_cast<float>(_x)*m_physicalX+m_physicalX/2.0,static_cast<float>(_y)*m_physicalY+m_physicalY/2.0,static_cast<float>(_z)*m_physicalZ+m_physicalZ/2.0);

				const Eigen::Vector3f _world = (m_boundingBoxTransformation*_local);



				// initalize the values which determine the color of the voxel (gaussian backward warping)
				float _intensityNominator = 0.0;
				float _intensityDenominator = 0.0;


				// pop the minimum distance bin vector from the rotation queue
				std::vector<USFrame*>* _usFrameQueue = _rotQueue.pop();
				int _oldFrameSize = _usFrameQueue->size();

				// for each slice within the vector returned form the queue do a backward warping
				while( _usFrameQueue->size() > 0 )
				{
					USFrame *_frame = _usFrameQueue->back();

					_usFrameQueue->pop_back();

					Plane *_plane = _frame->plane;


					// determine the orthogonal distance from the voxel to the ultrasound plane
					const float _distance = _plane->getOrthogonalDistance (_world, _worldIntersection );

					// determine the point of minimum distance to the voxel on the plane (in pixel coordinates)
					const Eigen::Vector3f _localIntersection = (_frame->world2Image*_worldIntersection);


					// point is relevant for the voxel, its intensity has to be considered / it has to be within the image
					// NOTE: the _roiMin/Max is because of the use of the cropped image area within the video frame
					
					if ( _distance < _threshold && _localIntersection[0] >= m_roiMinX && _localIntersection[0] <= m_roiMaxX && _localIntersection[1] >= m_roiMinY && _localIntersection[1] <= m_roiMaxY)
					{

						// update the distance
						_frame->distance = _distance;

						// and push it into the temporary queue, as the point is not to be considered again for the same voxel
						_rotQueue.push(_frame, _distance);




						// grab the video image

						// very slow code
						//HBITMAP bmp;
						//m_videoSegmentationTools.getFrame(_frame.videoFrame, bmp);
						//float _tempIntensity = m_videoSegmentationTools.getIntensity(bmp, (int)_localIntersection(0), (int)_localIntersection(1));




						// now compute the radius in the plane that is within distance specified by the threshold
						// that's simply pythagoras 
						//int _radius = static_cast<int>(sqrt ( _distance * _distance + _thresholdSQR ) / _scaleX);
						// instead of pythagoras a faster version is to take the precomputed distance values
						const int _radius = _radiusSelect.findRadius(_distance);

						std::vector<Segment> *_SegmentVector;

						// check the hash map if for the radius there exist already Segment information
						// the Segment contains the discretization of the circle for a given radius
						// e.g. K pixels above the center, move from -M to +M (relative position)
						// it is like the active edge table in the rasterizer
						_hashMapIterator = _hashMap.find(_radius);


						// in case there is no Segment information available, we have to generate it
						if ( _hashMapIterator == _hashMap.end() )
						{
							_SegmentVector = new std::vector<Segment>;

							// compute the relative position Segments of the circle
							// rasterization table for a circle of given radius (relative coordinates => 0,0 as first parameters)
							_graphicTools.bresenhamCircle(0,0, _radius, _SegmentVector);


							// place the information in the hashmap
							_hashMap[_radius]=_SegmentVector;
						} // Segment information is available, grab it
						else
						{
							// discretization was computed previously, so we simply have to acquire the pointer
							_SegmentVector = _hashMapIterator->second;
						}

						// now we look up follow the information contained in the vector (visit all the pixels, which are in range)
						// cache the size of the vector to save computations
						std::vector<Segment>::iterator _itBegin = _SegmentVector->begin();
						std::vector<Segment>::iterator _itEnd = _SegmentVector->end();

						for(std::vector<Segment>::iterator it = _itBegin; it != _itEnd;++it )
						{
							Segment _Segment = *it;

							int _originalBeginX = _Segment.beginX;
							int _originalEndX = _Segment.endX;
							int _originalY = _Segment.Y;


							// as we have only relative positional information from the mid-point/bresenham
							// we have to transform them in absolute positions on the image
							// we have to adapt the coordinates, as we have cropped the image in memory (that's why using _roiMin*)
							_Segment.Y			+= static_cast<int>(_localIntersection[1])- m_roiMinY;
							_Segment.beginX	+= static_cast<int>(_localIntersection[0])- m_roiMinX;
							_Segment.endX		+= static_cast<int>(_localIntersection[0])- m_roiMinX;

							// only consider points which are within the image
							if ( _Segment.Y < m_roiHeight && _Segment.Y >= 0 )
							{

								int _beginX, _endX;
								bool _draw = true;

								// if the Segment is outside the window, forget about it
								if ( _Segment.beginX >= m_roiWidth )
									_draw = false;
								if ( _Segment.endX < 0 )
									_draw = false;

								// if Segment-range is within the area of interest backward-warp it
								if ( _draw && ((_Segment.beginX	 >= 0 && _Segment.beginX <= m_roiWidth-1) ||( _Segment.endX >= 0 && _Segment.endX <= m_roiWidth-1)) )
								{
									// clip to the boundary, define the rasterization limits (the bounding box of the clipping mask)
									_beginX = std::max( std::min( _Segment.beginX, m_roiWidth-1), 0);
									_endX = std::min(std::max( _Segment.endX, 0 ), m_roiWidth-1 );



									IplImage *_img = m_imageVector[_frame->vectorIndex][(_frame->videoFrame)];
									// get the pointers to the ultrasound image and the cropping mask
									float *_imgPointer= &(((float*)(_img->imageData + _Segment.Y*_img->widthStep))[_beginX]); // the R value of the RGB image


									// Using the auxiliary video data (pre-segmentation)
									IplImage *_auxImg = m_auxImageVector[_frame->vectorIndex][(_frame->videoFrame)];
									unsigned char *_auxImgPointer = &(((unsigned char*)(_auxImg->imageData + _Segment.Y*_auxImg->widthStep))[_beginX]);  

									// now pixel by pixel backward warping (given the radius)
									unsigned char *_croppingMaskPointer = &(((unsigned char*)(m_croppingMask->imageData + _Segment.Y*m_croppingMask->widthStep))[_beginX]); 

									for (int x_pos=_beginX;x_pos <= _endX; x_pos++)
									{

										// now apply the cropping
										// get the intensity of the cropping mask
										const float _maskIntensity = (float)*_croppingMaskPointer++;

										// if the intensity equals white, then we are in non-cropped area
										//working for intensity reconstruction:					if ( _maskIntensity > 0.0)
										// using the auxiliary image data (pre-segmenation video)
										// the video is white at ROI, black otherwise
										const float _auxIntensity = (float)*_auxImgPointer++;
										if( _auxIntensity > 0.0 )
										{
											Eigen::Vector3f _normalVec, _spaceVec;

											const float _velocity = *_imgPointer++;
											// BUG: RGB is always 0,0,0
											// get the intensity of the Doppler ultrasound image (RGB)
											/*const int _intensityR = static_cast<int>(*_imgPointer++);
											const int _intensityG = static_cast<int>(*_imgPointer++);
											const int _intensityB = static_cast<int>(*_imgPointer++);

											// now the pointer points to the next pixel

											m_colLUT->colorToVelocity(_intensityR, _intensityG, _intensityG, _velocity);
											*/

											//cvSaveImage("c:\\Resultant.jpg",m_imageVector.at(_frame->videoFrame));
											// compute the distance for the point => _distance is just the orthogonal distance to the plane
											// Gaussian backward warping

											// determine the position of the pixel in world coordinates (mm)
											const Eigen::Vector3f _temp(x_pos+m_roiMinX, _Segment.Y+m_roiMinY, 0.0);
											const Eigen::Vector3f _worldPoint = (_frame->image2World * _temp);

											// determine squared euclidean distance of the point to the pixel in world coordinates
											const float _euclideanDistance =  static_cast<float>((_world - _worldPoint).norm());


											m_recon->getNormal(x_pos, _Segment.Y, _normalVec);


											/*	if ( x_pos == 65 && _Segment.Y == 13 && _frame->videoFrame == 2)
											{
											Eigen::Vector3f _tmp = _frame->local2World*_normalVec;
											std::cout << "Velocity: " << _velocity << std::endl;
											std::cout << "Normal: " << _normalVec << std::endl;
											std::cout << "Normal Space: " << _tmp << std::endl;
											//										std::cout << "RGB: " << _intensityR << ", " << _intensityG << ", " << _intensityB << std::endl;
											}*/
											//std::cout << _temp[0] << ", " << _temp[1] << std::endl;
											//if (fabs(_velocity) > 0.001 )
											{
												// TJK 29.03.2010: ROTATION BUG!!! [fixed]
												Eigen::Vector3f _tmp = (_frame->local2World.rotation()*_normalVec);//.normalized();
												int _index = _z*m_voxelsX*m_voxelsY+_y*m_voxelsX+_x;

												/*if (!_recon->getIndex(_worldPoint, _index))
												std::cout << "There must be some bug!" << std::endl;


												*/
												/*if ( _index == 3 )
												{
												std::cout << "Y: " << _Segment.Y+m_roiMinY << "  X: " << x_pos+m_roiMinX << std::endl;
												std::cout << "Original: " << _originalY << " | " << _originalBeginX << " to " << _originalEndX << "  Clipped: " << _beginX << " to " << _endX << " Radius: " << _radius << std::endl;
												std::cout << _worldPoint  << std::endl << std::endl << _index << std::endl << std::endl << _velocity << std::endl << std::endl << _temp << std::endl << "X: " << _x << "  Y: " << _y << "  Z: " << _z << "  Frame: " << (_frame->videoFrame) << std::endl << "----------------------" << std::endl;
												}*/


												/*if ( _x == m_voxelsX/2 && _y == m_voxelsY/2 && _z == m_voxelsZ/2 )
												{

												std::cout << "Velocity: " << _velocity << std::endl;
												std::cout << "Normal Space: " << _tmp << std::endl;
												}*/
												m_recon->pushMeasurement(_tmp, _velocity, _index, _frame->videoFrame); 
												/*if ( _lastIndex != -1 && _lastIndex != _frame->videoFrame)
												std::cout << std::endl;
												*/
												_lastIndex = _frame->videoFrame;
											}
										}

									}
								}
							}

						}

					}
					else // point too far away or not within the region of interest, but distance was updated so we have to update it in the rotation queue
					{
						// update the distance
						_frame->distance = _distance;

						// and push it into the temporary queue, as the point is not to be considered again for the same voxel
						_rotQueue.push(_frame, _distance);
					}

				}
				// instead of delete put empty vectors back to rotation queue, so we save new operations
				_usFrameQueue->reserve(_oldFrameSize);
				_rotQueue.recycle(_usFrameQueue);



			} // end of while


			// free memory
			stdext::hash_map<int, std::vector<Segment>*>::iterator _hashBegin = _hashMap.begin();
			stdext::hash_map<int, std::vector<Segment>*>::iterator _hashEnd = _hashMap.end();

			for(_hashMapIterator = _hashBegin; _hashMapIterator!= _hashEnd; ++_hashMapIterator)
			{
				if ( _hashMapIterator->second )
					delete _hashMapIterator->second;
			}





			std::stringstream _ss;
			_ss << "Measurement collection finished after " << m_dataTools.computeElapsedTime(_startTimer, clock());
			ui.textEdit->append(_ss.str().c_str());

			return true;
}

void Compounding::batchIntensityCompound(std::string filename,  COMPOUNDING_MODE compoundingMode, BATCH_DISTANCE_TYPE type, BATCH_DISTANCE distance, int logCompression)
{
	ui.postProcessingTargetFilenameEdit->setText(QString(filename.c_str()));

	// max distance
	ui.maxDistanceScalarCombo->setCurrentIndex(distance);

	// voxel extent or diameter
	ui.maxDistanceUnitCombo->setCurrentIndex(type);

	ui.checkBoxLogCompression->setChecked(logCompression == 1);


	ui.reconstructionModeComboBox->setCurrentIndex((int)compoundingMode);
	
	qApp->processEvents();

	postProcessRF();
}

void Compounding::batchButterwothEnvelope(std::string filename, int lowCutOff, int highCutOff, int order, int frequency)
{
	ui.postProcessingTargetFilenameEdit->setText(QString(filename.c_str()));

	{
		std::stringstream _ss;
		_ss << lowCutOff;
		ui.butterworthLow->setText(QString(_ss.str().c_str()));
	}

	{
		std::stringstream _ss;
		_ss << highCutOff;
		ui.butterworthHigh->setText(QString(_ss.str().c_str()));
	}

	{
		std::stringstream _ss;
		_ss << order;
		ui.butterworthOrder->setText(QString(_ss.str().c_str()));
	}

	{
		std::stringstream _ss;
		_ss << frequency;
		ui.butterworthSamplingFrequency->setText(QString(_ss.str().c_str()));
	}

	preProcessRF();
}

void Compounding::batchRun(BATCH_MODE mode, std::string calibrationFile, std::string configurationFile, std::string trackingFile, float voxelSize, BATCH_DISTANCE_TYPE type, BATCH_DISTANCE distance, COMPOUNDING_MODE compoundingMode, std::string targetName, std::string referenceName)
{
	// max distance
	ui.maxDistanceScalarCombo->setCurrentIndex(distance);

	// voxel extent or diameter
	ui.maxDistanceUnitCombo->setCurrentIndex(type);


	if ( mode == BATCH_MODE::RF ) {
		std::cout << "RF MODE" << std::endl;

		ui.calibrationRFFileEdit->setText(QString(calibrationFile.c_str()));
		ui.trackingRFFileEdit->setText(QString(trackingFile.c_str()));
		m_trackingFileName = QString(trackingFile.c_str());
		ui.RFXMLFileEdit->setText(QString(configurationFile.c_str()));
		
		m_rfDataFileName = QString(configurationFile.c_str());
		m_rfTrackingFileName = QString(trackingFile.c_str());
		m_rfCalibrationFileName = QString(calibrationFile.c_str());

		//output filename
		std::string _temp = configurationFile;
		int _idx1 =  _temp.find_last_of( '/' );
		int _idx2 = _temp.find_first_of( '.' );
		std::string _name = _temp.substr( _idx1+1, (_idx2-_idx1)  )+"RFC.xml";
		ui.outputRFFileName->setText(QString(_name.c_str()));

		loadRFData();
		if (referenceName.length() > 0 )
		{
			ui.referenceTargetRFCheckBox->setChecked(true);
			for(int i=0;i<ui.referenceNameRFComboBox->count();i++) {
				std::stringstream _ss;
				_ss << ui.referenceNameRFComboBox->itemText(i). toLatin1().data();
				if ( _ss.str() == referenceName)
					ui.referenceNameRFComboBox->setCurrentIndex(i);
			}
		}
		else
			ui.referenceTargetRFCheckBox->setChecked(false);

		int c= ui.probeNameRFComboBox->count();
		for(int i=0;i<ui.probeNameRFComboBox->count();i++) {
				std::stringstream _ss;
				_ss << ui.probeNameRFComboBox->itemText(i). toLatin1().data();
				if ( _ss.str() == targetName)
					ui.probeNameRFComboBox->setCurrentIndex(i);
			}

		initializeRFData();

		m_voxelSize = voxelSize;

		ui.reconstructionModeComboBox->setCurrentIndex((int)compoundingMode);
		compoundRF();

	}
	else
	{
		std::cout << "B-MODE" << std::endl;
	}
}

void Compounding::batchRunBMODE(BATCH_MODE mode, std::string calibrationFile, std::string outputFilename, std::vector< std::string> videoFile, std::vector<std::string> trackingFile, float voxelSize, BATCH_DISTANCE_TYPE type, BATCH_DISTANCE distance, COMPOUNDING_MODE compoundingMode, bool logCompress, std::string targetName, std::string referenceName)
{
	ui.checkBoxLogCompression->setChecked(logCompress);
	// max distance
	ui.maxDistanceScalarCombo->setCurrentIndex(distance);

	// voxel extent or diameter
	ui.maxDistanceUnitCombo->setCurrentIndex(type);


	if ( mode == BATCH_MODE::BMODE ) {
		std::cout << "B MODE" << std::endl;

		ui.configurationFileEdit->setText(QString(calibrationFile.c_str()));
		ui.trackingFileEdit->setText(QString(trackingFile[0].c_str()));
		m_trackingFileName = QString(trackingFile[0].c_str());
		ui.videoFileEdit->setText(QString(videoFile[0].c_str()));
		
		m_videoFileName = QString(videoFile[0].c_str());
		//m_rfDataFileName = QString(configurationFile.c_str());
		m_trackingFileName = QString(trackingFile[0].c_str());
		//m_rfCalibrationFileName = QString(calibrationFile.c_str());
		m_configurationFileName = QString(calibrationFile.c_str());

		//output filename
		std::string _temp = outputFilename;
		int _idx1 =  _temp.find_last_of( '/' );
		int _idx2 = _temp.find_first_of( '.' );
		std::string _name = _temp.substr( _idx1+1, (_idx2-_idx1)  )+"RFC.xml";
		ui.outputFileName->setText(QString(_name.c_str()));

		loadData();
		if (referenceName.length() > 0 )
		{
			ui.referenceTargetCheckBox->setChecked(true);
			for(int i=0;i<ui.referenceNameComboBox->count();i++) {
				std::stringstream _ss;
				_ss << ui.referenceNameComboBox->itemText(i). toLatin1().data();
				if ( _ss.str() == referenceName)
					ui.referenceNameComboBox->setCurrentIndex(i);
			}
		}
		else
			ui.referenceTargetCheckBox->setChecked(false);

		int c= ui.probeNameComboBox->count();
		for(int i=0;i<ui.probeNameComboBox->count();i++) {
				std::stringstream _ss;
				_ss << ui.probeNameComboBox->itemText(i). toLatin1().data();
				if ( _ss.str() == targetName)
					ui.probeNameComboBox->setCurrentIndex(i);
			}

		initializeData();

		

		for(int i=1;i<trackingFile.size();i++)
		{
			ui.trackingFileEdit->setText(QString(trackingFile[i].c_str()));
		m_trackingFileName = QString(trackingFile[i].c_str());
		ui.videoFileEdit->setText(QString(videoFile[i].c_str()));
		
		m_videoFileName = QString(videoFile[i].c_str());
		//m_rfDataFileName = QString(configurationFile.c_str());

		attachData();
	
		}


		m_voxelSize = voxelSize;

		ui.reconstructionModeComboBox->setCurrentIndex((int)compoundingMode);

		compound();

	}
	else
	{
		std::cout << "RF-MODE" << std::endl;
	}
}
	


bool Compounding::compoundRFVolume(COMPOUNDING_MODE _compoundingMode)
{

	glWidget->stopTimer();
	m_slicesVisualization->stopTimer();
	clock_t _tmpTimer = 0;
	clock_t _tmpTimer2 = 0;
	clock_t _startTimer = clock();

	m_rfRFDfilename = ui.outputRFFileName->text(). toLatin1().data();

	// first step is to determine the bounding box
	if ( ui.MasterFrameRFRadio->isChecked()) 
	{
		determineMinimumBoundingBoxfromMasterFrameRF();
	}

	else if ( ui.AutomaticRFRadio->isChecked())
	{
		determineMinimumBoundingBoxfromRF();
		m_masterFrameXML = -1; 
	}

	else if ( ui.ManualRFRadio->isChecked())
	{
		determineManualBoundingBoxfromRF();
		m_masterFrameXML = -1; 
	}


	// determine the maximum dimension of the voxel spacing
	float _maxDim = m_voxelSize;
	float _voxelDiameter = sqrt(m_physicalX*m_physicalX+m_physicalY*m_physicalY+m_physicalZ*m_physicalZ)/2.0f;

	// choose the distance threshold, that is the maximum distance a US slice can have
	// in order to be considered for contributing intensity to a voxel

	float _distanceScalar = 0.0f;

	switch(ui.maxDistanceScalarCombo->currentIndex())
	{
	case 0: { _distanceScalar = DISTANCE_SCALAR; break;}
	case 1: { _distanceScalar = 1.0; break;}
	case 2: { _distanceScalar = 1.1; break;}
	case 3: { _distanceScalar = 1.2; break;}
	case 4: { _distanceScalar = 1.3; break;}
	case 5: { _distanceScalar = 1.4; break;}
	case 6: { _distanceScalar = 1.5; break;}
	};

	float _threshold = 0.0f;
	switch(ui.maxDistanceUnitCombo->currentIndex())
	{
	case 0: { _threshold = _maxDim; break;}
	case 1: { _threshold = _voxelDiameter; break;}
	};

	// now we compute the actual threshold value from the UI selection
	_threshold = _threshold * _distanceScalar;


	// set up the RF file processing
	// TODO generate the pointer: m_rfFileStream

	float _physicalOffsetX = 0;
	float _physicalOffsetY = 0;
	float _physicalOffsetZ = 0;

	// subdivide into 8 sub-volumes
	float _extentX = m_extent.dimX/2.0;
	float _extentY = m_extent.dimY/2.0;
	float _extentZ = m_extent.dimZ/2.0;

	float _margin = 8.0;

	// the bounding box is a little bit extended for the intersecting scanlines
	float _extentXBB = m_extent.dimX/2.0 + 4.0*_threshold;
	float _extentYBB = m_extent.dimY/2.0 + 4.0*_threshold;
	float _extentZBB = m_extent.dimZ/2.0 + 4.0*_threshold;

	Eigen::Transform3f _subVolumeTransformation;
	Eigen::Transform3f _subVolumeTransformationBB;
	BBExtent _subVolumeExtent;


	std::string _filename;
	std::vector<std::string> _subVolumeFileNames;


	SYSTEM_INFO _info;
	GetSystemInfo(&_info );
	{
	std::stringstream _ss;
	_ss << _info.dwNumberOfProcessors << " cores detected!";
	ui.textEdit->append(QString(_ss.str().c_str()));
				qApp->processEvents();
	}
	int _numOfThreads = ceil(0.5 *  _info.dwNumberOfProcessors);

	omp_set_num_threads(_numOfThreads);
		




	for (int i=0;i<	8;i++)
	{
		
		// This will just contain existing pointers and therefore the pointer will not have to deleted
#ifdef USE_USHORT16
		std::vector<Scanline<unsigned short>*> _intersectingScanlines;
#endif

#ifdef USE_FLOAT32 || USE_FLOAT5X32
		std::vector<Scanline<float>*> _intersectingScanlines;
#endif

		{
			std::stringstream _ss;
			_ss << "Reconstructing part " << i+1 << "/ 8";
			ui.textEdit->append(_ss.str().data());
		}

		switch(i)
		{
			case 0: // origin 
			{
				_subVolumeTransformation = m_boundingBoxTransformation;
				Eigen::Transform3f _temp;
				_temp.setIdentity();
				_temp.translation() =  Eigen::Vector3f(-_margin/2.0f*_threshold, -_margin/2.0f*_threshold, -_margin/2.0f*_threshold);
				_subVolumeTransformationBB = m_boundingBoxTransformation * _temp;
				_physicalOffsetX = 0;
				_physicalOffsetY = 0;
				_physicalOffsetZ = 0;
				_filename = "subvolume0.tmp";
				_subVolumeExtent.dimX = _extentXBB;
				_subVolumeExtent.dimY = _extentYBB;
				_subVolumeExtent.dimZ = _extentZBB;
				break;
			}
		case 1:
			{
				Eigen::Transform3f _temp;
				_temp.setIdentity();
				_temp.translation() = Eigen::Vector3f(_extentX, 0.0, 0.0);
				_subVolumeTransformation = m_boundingBoxTransformation * _temp;
				_temp.translation() = Eigen::Vector3f( _extentX - _margin/2.0f*_threshold, -_margin/2.0f*_threshold, -_margin/2.0f*_threshold);
				_subVolumeTransformationBB = m_boundingBoxTransformation * _temp;
				_physicalOffsetX = (m_voxelsX/2.0)*m_physicalX;
				_physicalOffsetY = 0;
				_physicalOffsetZ = 0;
				_filename = "subvolume1.tmp";
				_subVolumeExtent.dimX = _extentXBB;
				_subVolumeExtent.dimY = _extentYBB;
				_subVolumeExtent.dimZ = _extentZBB;
				break;
			}

		case 4: // updated (fomerly 2)
			{
				Eigen::Transform3f _temp;
				_temp.setIdentity();
				_temp.translation() = Eigen::Vector3f(0.0, 0.0, _extentZ);
				_subVolumeTransformation = m_boundingBoxTransformation * _temp;
				_temp.translation() = Eigen::Vector3f(-_margin/2.0f*_threshold, -_margin/2.0f*_threshold, _extentZ -_margin/2.0f* _threshold);
				_subVolumeTransformationBB = m_boundingBoxTransformation * _temp;
				_physicalOffsetX = 0;
				_physicalOffsetY = 0;
				_physicalOffsetZ = (m_voxelsZ/2.0)*m_physicalZ;
				_filename = "subvolume4.tmp";
				_subVolumeExtent.dimX = _extentXBB;
				_subVolumeExtent.dimY = _extentYBB;
				_subVolumeExtent.dimZ = _extentZBB;
				break;
			}

		case 5: // updated (fomerly 4)
			{
				Eigen::Transform3f _temp;
				_temp.setIdentity();
				_temp.translation() = Eigen::Vector3f(_extentX, 0.0, _extentZ);
				_subVolumeTransformation = m_boundingBoxTransformation * _temp;
				_temp.translation() = Eigen::Vector3f(_extentX - _margin/2.0f*_threshold, -_margin/2.0f*_threshold, _extentZ - _margin/2.0f*_threshold);
				_subVolumeTransformationBB = m_boundingBoxTransformation * _temp;
				_physicalOffsetX = (m_voxelsX/2.0)*m_physicalX;
				_physicalOffsetY = 0;
				_physicalOffsetZ = (m_voxelsZ/2.0)*m_physicalZ;
				_filename = "subvolume5.tmp";
				_subVolumeExtent.dimX = _extentXBB;
				_subVolumeExtent.dimY = _extentYBB;
				_subVolumeExtent.dimZ = _extentZBB;
				break;
			}

		case 2: // updated (formerly 4)
			{
				Eigen::Transform3f _temp;
				_temp.setIdentity();
				_temp.translation() = Eigen::Vector3f(0.0, _extentY, 0.0);
				_subVolumeTransformation = m_boundingBoxTransformation * _temp;
				_temp.translation() = Eigen::Vector3f(-_margin/2.0f*_threshold, _extentY -_margin/2.0f*_threshold, -_margin/2.0f*_threshold);
				_subVolumeTransformationBB = m_boundingBoxTransformation * _temp;
				_physicalOffsetX = 0;
				_physicalOffsetY = (m_voxelsY/2.0)*m_physicalY;
				_physicalOffsetZ = 0;
				_filename = "subvolume2.tmp";
				_subVolumeExtent.dimX = _extentXBB;
				_subVolumeExtent.dimY = _extentYBB;
				_subVolumeExtent.dimZ = _extentZBB;
				break;
			}

		case 3: // updated (formerly 5)
			{
				Eigen::Transform3f _temp;
				_temp.setIdentity();
				_temp.translation() = Eigen::Vector3f(_extentX, _extentY, 0.0);
				_subVolumeTransformation = m_boundingBoxTransformation * _temp;
				_temp.translation() = Eigen::Vector3f( _extentX -_margin/2.0f* _threshold, _extentY -_margin/2.0f* _threshold, -_margin/2.0f*_threshold);
				_subVolumeTransformationBB = m_boundingBoxTransformation * _temp;
				_physicalOffsetX = (m_voxelsX/2.0)*m_physicalX;
				_physicalOffsetY = (m_voxelsY/2.0)*m_physicalY;
				_physicalOffsetZ = 0;
				_filename = "subvolume3.tmp";
				_subVolumeExtent.dimX = _extentXBB;
				_subVolumeExtent.dimY = _extentYBB;
				_subVolumeExtent.dimZ = _extentZBB;
				break;
			}
		case 6: // updated
			{
				Eigen::Transform3f _temp;
				_temp.setIdentity();
				_temp.translation() = Eigen::Vector3f(0.0, _extentY, _extentZ);
				_subVolumeTransformation = m_boundingBoxTransformation * _temp;
				_temp.translation() = Eigen::Vector3f(-_margin/2.0f*_threshold, _extentY -_margin/2.0f* _threshold, _extentZ -_margin/2.0f* _threshold);
				_subVolumeTransformationBB = m_boundingBoxTransformation * _temp;
				_physicalOffsetX = 0;
				_physicalOffsetY = (m_voxelsY/2.0)*m_physicalY;
				_physicalOffsetZ = (m_voxelsZ/2.0)*m_physicalZ;
				_filename = "subvolume6.tmp";
				_subVolumeExtent.dimX = _extentXBB;
				_subVolumeExtent.dimY = _extentYBB;
				_subVolumeExtent.dimZ = _extentZBB;
				break;
			}

		case 7:
			{
				Eigen::Transform3f _temp;
				_temp.setIdentity();
				_temp.translation() = Eigen::Vector3f(_extentX, _extentY, _extentZ);
				_subVolumeTransformation = m_boundingBoxTransformation * _temp;
				_temp.translation() = Eigen::Vector3f(_extentX -_margin/2.0f* _threshold, _extentY - _margin/2.0f*_threshold, _extentZ -_margin/2.0f* _threshold);
				_subVolumeTransformationBB = m_boundingBoxTransformation * _temp;
				_physicalOffsetX = (m_voxelsX/2.0)*m_physicalX;
				_physicalOffsetY = (m_voxelsY/2.0)*m_physicalY;
				_physicalOffsetZ = (m_voxelsZ/2.0)*m_physicalZ;
				_filename = "subvolume7.tmp";
				_subVolumeExtent.dimX = _extentXBB;
				_subVolumeExtent.dimY = _extentYBB;
				_subVolumeExtent.dimZ = _extentZBB;
				break;
			}
	
		}

		// putting the filenames into vector not needed anymore
		/*{
			// get the path to the the RF data file
			std::string _temp = m_rfDataFileName. toLatin1().data();
			std::string _path = _temp.substr( 0, _temp.find_last_of( '/' ) +1 );
			_filename = _path + _filename;
			// now add the subvolume file to the vector
			_subVolumeFileNames.push_back(_filename);
		}*/

	//	if ( i < 2 )
		//	continue;
		OrientedBoundingBox *_obb = new OrientedBoundingBox(_subVolumeTransformationBB, _subVolumeExtent);

		// now check which scanlines intersect the subvolume

		{
			std::stringstream _ss;
			_ss << "Computing intersecting scanlines ... ";
			ui.textEdit->append(QString(_ss.str().c_str()));
			qApp->processEvents();
		}



		if ( true ) // compounding
		{

			int _scanlines = m_scanlineDataVector.size();

			for(int j=0; j<_scanlines;j++)
			{
#ifdef USE_USHORT16
				Scanline<unsigned short>* _scl = m_scanlineDataVector[j];
#endif

#ifdef USE_FLOAT32 || USE_FLOAT5X32
				Scanline<float>* _scl = m_scanlineDataVector[j];
#endif

				// check if it intersect, if 'YES' then add it to the vector
				if ( _obb->lineSegmentIntersection(_scl->getLine() ) )
				{
					_intersectingScanlines.push_back(_scl);
				}
			}

			// now sort the scanlines according to image and then scanline index, which will then allow efficient gathering of the data
			{
				std::stringstream _ss;
				_ss << "Intersecting scanlines: " << _intersectingScanlines.size() << "  ( " << 100.0f* (float)_intersectingScanlines.size()/(float)_scanlines << " % )";
				ui.textEdit->append(_ss.str().data());
			}
#ifdef USE_USHORT16
			std::sort(_intersectingScanlines.begin(), _intersectingScanlines.end(), ScanlineComparator<unsigned short>(true));
#endif

#ifdef USE_FLOAT32 || USE_FLOAT5X32
			std::sort(_intersectingScanlines.begin(), _intersectingScanlines.end(), ScanlineComparator<float>(true));
#endif

			// in the first iteration make a time measurement 
			if ( i == 0 ) {
				_tmpTimer = clock(); }
			else
			{
				// we put existing writer threads to sleep (20 sec), in order to speed up reading data of the scanline

				for (int pro=0; pro < m_subvolumeWriting.size();pro++)
				{
					CompoundRFSubvolumeWritingThread* _thread = m_subvolumeWriting[pro];

					/*if (! _thread->isFinished() )
						_thread->wait(2.0 * _tmpTimer/CLOCKS_PER_SEC * 1000);*/
				}

				_tmpTimer2 = clock(); 
			}

			// now the scanlines are sorted and we can linearly load the data (but don't forget to unload it at the end)
			// BUG: no intersecting scanlines leads to exception
		//	int _startIndex = _intersectingScanlines[0]->getRfImageIndex();
		//	int _endIndex = _intersectingScanlines[_intersectingScanlines.size()-1]->getRfImageIndex();
			// data conversion
#ifdef USE_USHORT16
			RFData<unsigned short> * _rfData = NULL;
#endif

#ifdef USE_FLOAT32 || USE_FLOAT5X32
			RFData<float> * _rfData = NULL;
#endif

			int _lastRF = -1;
			for(int j=0; j < _intersectingScanlines.size(); j++ )
			{
				int _currentRFidx = _intersectingScanlines[j]->getRfImageIndex();
				int _currentBeamIdx =_intersectingScanlines[j]->getBeamIndex();
				UltrasoundSettings *_ultrasoundSettings = _intersectingScanlines[j]->getUltrasoundSettings();


				// as we only copy entire RF frames, we need to copy only when the image index changes
				if (_lastRF != _currentRFidx )
				{
					_lastRF = _currentRFidx;


					// if the RF frame is changing, dispose the old before building a new one
					if ( _rfData ) delete _rfData;

#ifdef USE_USHORT16
					_rfData = new RFData<unsigned short>();
#endif
#ifdef USE_FLOAT32 || USE_FLOAT5X32
					_rfData = new RFData<float>();
#endif
					_rfData->setNumberOfSamples(_ultrasoundSettings->samplesAlongScanline);
					_rfData->setNumberOfScanlines(_ultrasoundSettings->scanlines);
					// now copy the frame from the file to pointer
					// data conversion
#ifdef USE_USHORT16
					ushort16 *_rfPtr = reinterpret_cast<ushort16*>(m_rfFileStreamVector[_intersectingScanlines[j]->getFileIndex()]->getFrameByIndex(_currentRFidx));
					// have to copy the ENTIRE RF frame into memory, no FILE pointer
					unsigned short * _memPtr = new unsigned short[_ultrasoundSettings->rfFrameSizeElements];
#endif

#ifdef USE_FLOAT32 || USE_FLOAT5X32
					float32 *_rfPtr = reinterpret_cast<float32*>(m_rfFileStreamVector[_intersectingScanlines[j]->getFileIndex()]->getFrameByIndex(_currentRFidx));
					// have to copy the ENTIRE RF frame into memory, no FILE pointer
					float * _memPtr = new float[_ultrasoundSettings->rfFrameSizeElements];
#endif

					memcpy(_memPtr,&_rfPtr->byte.c[0],_ultrasoundSettings->rfFrameSizeBytes);

					_rfData->setData(_memPtr);



					//std::cout << "Frame: "<< _currentRFidx << " found: "<< *_rfData->getData()-2 << std::endl;

					// now free the pointer from file
					delete[] _rfPtr;


				}


				// now copy the scanline (which of couse has to be deleted at the end)
				// we might not need all scanlines, therefore copy on a scanline to scanline basis to be memory efficient
				// data conversion
#ifdef USE_USHORT16

				unsigned short *_tmpScanline = _rfData->getScanline(_currentBeamIdx);
				unsigned short *_tmpData = new unsigned short[_ultrasoundSettings->samplesAlongScanline];
				memcpy(_tmpData,_tmpScanline, sizeof(unsigned short)*_ultrasoundSettings->samplesAlongScanline);

#endif

#ifdef USE_FLOAT32

				float *_tmpScanline = _rfData->getScanline(_currentBeamIdx);
				float *_tmpData = new float[_ultrasoundSettings->samplesAlongScanline];
				memcpy(_tmpData,_tmpScanline, sizeof(float)*_ultrasoundSettings->samplesAlongScanline);
#endif

#ifdef USE_FLOAT5X32

				float *_tmpScanline = _rfData->getScanline(_currentBeamIdx);
				float *_tmpData = new float[_ultrasoundSettings->samplesAlongScanline*5];
				memcpy(_tmpData,_tmpScanline, sizeof(float)*5*_ultrasoundSettings->samplesAlongScanline);
#endif

				/*	if ( _currentRFidx == 366 && _currentBeamIdx == 80 )
				{
				float *_ptr = _tmpData;
				float _last = -1;
				for (int i=0; i<3072;i++)
				{
				if ( _last != *_ptr)
				{
				std::cout << i << " => " << *_ptr << std::endl;
				_last = *_ptr;
				}
				_ptr++;
				}
				int _debug = 5;
				}*/


				// now set the data in the scanline
				_intersectingScanlines[j]->setData(_tmpData);
			}


			if ( _rfData) delete _rfData;

			if ( i== 0) {
				_tmpTimer = clock() - _tmpTimer; }
			else
			{
				// update the timer to the max
				_tmpTimer = std::max(_tmpTimer, clock() - _tmpTimer2);
			}


			// now we are ready to compound a subvolume
			StandardMeasurementCollectionManager<IntensityAndDistance> * _collectionManager = new StandardMeasurementCollectionManager<IntensityAndDistance>((m_voxelsX/2.0),(m_voxelsY/2.0), (m_voxelsZ/2.0)); // Potential source of error!!!



			// next we subdivide the volume in sub-volumes, for each subvolume we determine the scanlines that intersect
			// followed by loading the data into the memory

			

			/*_info.dwNumberOfProcessors = 1;
			m_threadingFactor = 1;
			*/
			std::stringstream _ss;



			clock_t _startTimer = clock();
			ui.textEdit->append(QString("Compounding started"));
			qApp->processEvents();

			

			// here we count the bytes for each subvolume, to speed up later on file writing
			__int64 *_byteCounter = new __int64[(int)_numOfThreads];
			__int64 _fileSize = 0;

			

			int _share = (int)floor((float)(m_voxelsZ/2)/(float)_numOfThreads);
			int _startPos = 0;
			int _endPos = _share;

			ui.progressBar->setRange(0,100);
			// as to avoid too frequen updates of the GUI we do only updates when values change
			int _maxVoxels =  ceil(m_voxelsX/2.0)*ceil(m_voxelsY/2.0)*ceil(m_voxelsZ/2.0); // POTENTIAL SOURCE OF ERROR
			int _1percent = (int)((float)_maxVoxels/100.0);
			int _percentCounter = 0;
			int _percent = 0;
			ui.progressBar->setValue(_percent);
			qApp->processEvents();

			{
				std::stringstream _ss;
				_ss << "Start " << _numOfThreads << " compounding threads ...";
				ui.textEdit->append(QString(_ss.str().c_str()));
				qApp->processEvents();
			}
				m_crfvThread = new CompoundRFVolumeThread[(int)_numOfThreads];
			for (int pro=0; pro < _numOfThreads;pro++)
			{
				_byteCounter[pro] = 0;

				/*
				std::string _tmp = _filename.substr( 0, _filename.find_last_of( '.' )  );
				std::stringstream _ss;
				_ss << _tmp << "_" << pro << ".tmp";
				*/
				// set the UI input 

				m_crfvThread[pro].setUserConfiguration(ui.maxDistanceScalarCombo->currentIndex(), ui.maxDistanceUnitCombo->currentIndex(), _1percent);
				m_crfvThread[pro].setDataSize(m_voxelsX/2, m_voxelsY/2, _endPos, m_physicalX, m_physicalY, m_physicalZ, _subVolumeExtent, 0, 0, _startPos);
				m_crfvThread[pro].setDataTarget(_collectionManager);	
				//m_crfvThread[pro].setFilename(_ss.str());
				m_crfvThread[pro].setTransformation(_subVolumeTransformation);
				m_crfvThread[pro].setData(&_intersectingScanlines);
				m_crfvThread[pro].setByteCounter(&_byteCounter[pro]);
				m_crfvThread[pro].setIndex(pro);

				// calculate voxel section for next thread
				_startPos = _endPos;
				_endPos = _startPos+_share;
				// the last processor will get all the remainaining slices, so make sure we compound to the last voxel
				if ( pro+2 == _numOfThreads)
					_endPos = (int)floor((float)(m_voxelsZ/2));
			}


			for (int pro=0; pro < _numOfThreads;pro++)
			{
				//if ( i == 0)
				{
				connect(&m_crfvThread[pro], SIGNAL(percentEvent()),
					this, SLOT(updatePercentage()), Qt::QueuedConnection);
				}
				
					
				m_crfvThread[pro].start();
				//m_crfvThread[pro].setPriority(QThread::Priority::HighestPriority);
			}



			{
				std::stringstream _ss;
				_ss << "Waiting for compounding threads to finish ...";
				ui.textEdit->append(QString(_ss.str().c_str()));
				qApp->processEvents();
			}

			bool _finished = false;
			// we have to run the thread for a certain amount of time, then we update the gui
			// if the threads have not finished with wait again for a certain time and so on...until all the threads are finished
			// this way we can update UI elements
			int _waits = 0;
			
			int _oldThreadSum = -1;
			do 
			{

				/*for (int pro=0; pro < _numOfThreads;pro++)
				{
				m_crfvThread[pro].wait(150);
				}*/
				//qApp->processEvents();

				std::string _tmpStr;

				_finished = true;
				int _threadSum = 0;
				for (int pro=0; pro < _numOfThreads;pro++)
				{
					if (!m_crfvThread[pro].isFinished()) {
						_finished = false;
						std::stringstream _ss;
						_ss << pro;
						if ( _tmpStr.length() > 0 )
							_tmpStr= _tmpStr+","+_ss.str();
						else
							_tmpStr = _ss.str();
						m_crfvThread[pro].wait(1500);
						_threadSum+=pro;
					}
				}
				if ( !_finished) {
					if ( _threadSum == 11 )
						_waits++;
					if ( _waits > 25 )
						int _debug = 5;
					if (_threadSum != _oldThreadSum )
					{
					std::stringstream _ss;
					_ss << "Waiting for compounding thread(s) " << _tmpStr << " to finish ...";
					ui.textEdit->append(QString(_ss.str().c_str()));
					qApp->processEvents();
					_oldThreadSum =_threadSum;

	
					} 
				}
			}

			while(!_finished);

			{
				std::stringstream _ss;
				_ss << "Compounding threads finished ...";
				ui.textEdit->append(QString(_ss.str().c_str()));
				qApp->processEvents();
			}
			for (int pro=0; pro < _numOfThreads;pro++)
				{
			disconnect(&m_crfvThread[pro], SIGNAL(percentEvent()),
					this, SLOT(updatePercentage()));
			}
			// remove the compounding threads
			delete[] m_crfvThread;

			_fileSize = 0;
			for (int pro=0; pro < _numOfThreads;pro++)
			{
				_fileSize += _byteCounter[pro];
			}

			/*{
			remove(_filename.c_str());
	
			_ss << "contig -n " << _filename << " " << _fileSize;
			std::cout << _ss.str() << std::endl;
			system(_ss.str().c_str());
			}
		
			HANDLE _hFile = CreateFileA( _filename.c_str(), GENERIC_WRITE ,  0, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL,  NULL);
			if ( _hFile ==  INVALID_HANDLE_VALUE )
				std::cout << "error opening file!" << std::endl;
				
			{
				std::stringstream _ss;
				_ss << "Merging sub-volumes...";
				ui.textEdit->append(QString(_ss.str().c_str()));
				qApp->processEvents();
			}
	


			// now we compute the file size from the individual components
			for (int pro=0; pro < _numOfThreads;pro++)
			{
				
				std::string _tmp = _filename.substr( 0, _filename.find_last_of( '.' )  );
				std::stringstream _ss;
				_ss << _tmp << "_" << pro << ".tmp";
				HANDLE _tmpFile = CreateFileA( _ss.str().c_str(), GENERIC_READ ,  0, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL,  NULL);
				int _size = (int)GetFileSize(_tmpFile, NULL);
				char *_data = new char[ _size];
				DWORD _readBytes;
				DWORD _writeBytes;
				ReadFile(_tmpFile,_data,_size, &_readBytes, NULL);
				WriteFile(_hFile,_data,_size, &_writeBytes,NULL);
				CloseHandle(_tmpFile);
				//remove( _ss.str().c_str());
				delete[] _data;
				//std::cout << "File size: " << _byteCounter[pro] << std::endl;
			}
			CloseHandle(_hFile);

			{
				std::stringstream _ss;
				_ss << "Merge finished!";
				ui.textEdit->append(QString(_ss.str().c_str()));
				qApp->processEvents();
			}*/

			delete[] _byteCounter;





			// in case we are in the last volume compounding or an old one has not finished writing we need to wait
			/*if (m_subvolumeWriting != NULL )
			{
			_finished = false;
			//m_subvolumeWriting->setPriority(QThread::HighestPriority);
			do
			{
			//m_subvolumeWriting->wait(1450);
			qApp->processEvents();

			_finished = m_subvolumeWriting->isFinished();

			if ( !_finished )
			Sleep(1500);
			}
			while(!_finished);

			delete m_subvolumeWriting;
			m_subvolumeWriting = NULL;

			}*/



			// now we may create a new working thread
			/*if ( i < 7 )
			{*/
			// now we create the thread that write the data to temporary file and removes the data from memory
			{
				std::string _temp = m_rfDataFileName. toLatin1().data();
				std::string _path = _temp.substr( 0, _temp.find_last_of( '/' ) +1 );

				CompoundRFSubvolumeWritingThread* _subvolumeWriting = new CompoundRFSubvolumeWritingThread(_filename, (m_voxelsX/2.0),  (m_voxelsY/2.0),  (m_voxelsZ/2.0), _maxVoxels, _fileSize);
				_subvolumeWriting->setDataTarget(_collectionManager);
				_subvolumeWriting->run();
				
				//_subvolumeWriting->setPriority(QThread::Priority::LowestPriority);

				//m_subvolumeWriting.push_back(_subvolumeWriting);
			}
			//_subvolumeWriting->setPriority(QThread::HighestPriority);
			/*}
			else
			{

			m_subvolumeWriting = new CompoundRFSubvolumeWritingThread(_filename, (m_voxelsX/2.0),  (m_voxelsY/2.0),  (m_voxelsZ/2.0), _maxVoxels);
			m_subvolumeWriting->setDataTarget(_collectionManager);
			m_subvolumeWriting->start();
			_finished = false;
			m_subvolumeWriting->setPriority(QThread::HighestPriority);
			do
			{
			//m_subvolumeWriting->wait(1500);
			qApp->processEvents();

			_finished = m_subvolumeWriting->isFinished();
			}
			while(!_finished);

			delete m_subvolumeWriting;
			m_subvolumeWriting = NULL;

			}*/



			//if ( _numOfThreads  > 1 )
			





			// finished compounding subvolume, free the memory!

			// the last thing we have to do is to remove all the scanline RF samples (we could also recycle them for the next subvolume though)
			#pragma omp parallel for
			for(int j=0;j<_intersectingScanlines.size();j++)
			{
#ifdef USE_USHORT16
				unsigned short *_ptr = _intersectingScanlines[j]->getData();
#endif

#ifdef USE_FLOAT32 || USE_FLOAT5X32
				float *_ptr = _intersectingScanlines[j]->getData();
#endif

				if (_ptr) delete[] _ptr;
			}


			// update the progress bar
			qApp->processEvents();
			m_percentage = -1;
			updatePercentage();
		}
	} // compounding

	
	{
		std::stringstream _ss;
		_ss << "Waiting for writing threads to finish ...";
		ui.textEdit->append(QString(_ss.str().c_str()));
		qApp->processEvents();
	}
	
	bool _finished = true;
	//m_subvolumeWriting->setPriority(QThread::HighestPriority);
	/*int _oldThreadSum = 0;
	do
	{
		_finished = true;
		std::string _tmpStr;
		//m_subvolumeWriting->wait(1450);
		qApp->processEvents();
		int _threadSum = 0;
		for (int pro=0; pro < m_subvolumeWriting.size();pro++)
		{
			CompoundRFSubvolumeWritingThread* _thread = m_subvolumeWriting[pro];

			if (! _thread->isFinished() )
			{
				_thread->setPriority(QThread::Priority::NormalPriority);
				_finished =false;

				std::stringstream _ss;
				_ss << pro;
				if ( _tmpStr.length() > 0 )
					_tmpStr= _tmpStr+","+_ss.str();
				else
					_tmpStr = _ss.str();
				_thread->wait(1500);
				_threadSum +=pro;
			}
		}
		if ( !_finished) {
			if ( _threadSum != _oldThreadSum ) {
			std::stringstream _ss;
			_ss << "Waiting for writing thread(s) " << _tmpStr << " to finish ...";
			ui.textEdit->append(QString(_ss.str().c_str()));
			qApp->processEvents();
			_oldThreadSum = _threadSum;
			}
		}
		else
		{
			int debug=5;
		}
	}
	while(!_finished);*/

	// merge the files together
	/*for(int i=0;i<8;i++)
		{
			std::string _temp = m_rfDataFileName. toLatin1().data();
			std::string _path = _temp.substr( 0, _temp.find_last_of( '/' ) +1 );

			std::stringstream _ss0, _ss1, _ss2, _ss3, _ss4;
			std::string _filename;
			_ss1 << "subvolume" << i << "1.tmp";
			_ss2 << "subvolume" << i << "2.tmp";
			_ss3 << "subvolume" << i << "3.tmp";
			_ss4 << "subvolume" << i << "4.tmp";
			_ss0 << "subvolume" << i << ".tmp";
		
			std::fstream _file0, _file1, _file2, _file3, _file4;
		
			// load the  file and check size
			_file1.open(_ss1.str(), std::ios::in  | std::ios::binary);
			// check file length and calculate the number of elements from the size
			_file1.seekg (0, std::ios::end);
			int _fileLength1 = _file1.tellg();
			_file1.seekg (0, std::ios::beg);

			char *_data1 = new char[_fileLength1];
			_file1.read(_data1,_fileLength1);
			_file1.close();
			// load the  file and check size
			_file2.open(_ss2.str(), std::ios::in  | std::ios::binary);
			// check file length and calculate the number of elements from the size
			_file2.seekg (0, std::ios::end);
			int _fileLength2 = _file2.tellg();
			_file2.seekg (0, std::ios::beg);

			char *_data2 = new char[_fileLength2];
			_file2.read(_data2,_fileLength2);
			_file2.close();

			// load the  file and check size
			_file3.open(_ss3.str(), std::ios::in  | std::ios::binary);
			// check file length and calculate the number of elements from the size
			_file3.seekg (0, std::ios::end);
			int _fileLength3 = _file3.tellg();
			_file3.seekg (0, std::ios::beg);

			char *_data3 = new char[_fileLength3];
			_file3.read(_data3,_fileLength3);
			_file3.close();

			// load the  file and check size
			_file4.open(_ss4.str(), std::ios::in  | std::ios::binary);
			// check file length and calculate the number of elements from the size
			_file4.seekg (0, std::ios::end);
			int _fileLength4 = _file4.tellg();
			_file4.seekg (0, std::ios::beg);

			char *_data4 = new char[_fileLength4];
			_file4.read(_data4,_fileLength4);
			_file4.close();



			_file0.open(_ss0.str(),std::ios::out  | std::ios::binary);

			_file0.write(_data1,_fileLength1);
			_file0.write(_data2,_fileLength2);
			_file0.write(_data3,_fileLength3);
			_file0.write(_data4,_fileLength4);
			_file0.close();
			delete[] _data1, _data2, _data3, _data4;
		}*/

	// add the filenames to the name vector

	for(int i=0;i<8;i++)
		{
			std::stringstream _ss;
			std::string _filename;
			_ss << "subvolume" << i << ".tmp";
			// get the path to the the RF data file
			std::string _temp = m_rfDataFileName. toLatin1().data();
			std::string _path = _temp.substr( 0, _temp.find_last_of( '/' ) +1 );
			_filename = /*_path +*/ _ss.str();
			// now add the subvolume file to the vector
			_subVolumeFileNames.push_back(_filename);
		}


	m_percentage = -1;
	updatePercentage();

	std::stringstream _ss;
	_ss << "Finished writing sub-volumes ... ";
	ui.textEdit->append(QString(_ss.str().c_str()));
	qApp->processEvents();

	// now we can merge together the individual subvolumes

	std::fstream _rfSubDataFile[8];
	std::fstream _mergedFile;
	std::fstream _indexFile;

	for(int i=0;i<8;i++)
	{
		_rfSubDataFile[i].open(_subVolumeFileNames[i], std::ios::in  | std::ios::binary);
		_rfSubDataFile[i].seekg (0, std::ios::beg);
	}


	XMLTools *_xmlTools = new XMLTools();

	

	_xmlTools->createFile("Reconstruction");
	{
		// dirty hack! we are taking the first value
		std::stringstream _ss;
		_ss << m_ultrasoundSettings[0].penetrationDepth;
		_xmlTools->handleData(true, _ss.str(), std::string(PENETRATION_DEPTH_TAG));
	}

	{
		std::stringstream _ss;
		_ss << m_voxelsX;
		_xmlTools->handleData(true, _ss.str(), "DataSize","VoxelsX");
	}

	{
		std::stringstream _ss;
		_ss << m_voxelsY;
		_xmlTools->handleData(true, _ss.str(), "DataSize","VoxelsY");
	}

	{
		std::stringstream _ss;
		_ss << m_voxelsZ;
		_xmlTools->handleData(true, _ss.str(), "DataSize","VoxelsZ");
	}


	{
		std::stringstream _ss;
		_ss << m_physicalX;
		_xmlTools->handleData(true, _ss.str(), "DataSize","PhysicalX");
	}

	{
		std::stringstream _ss;
		_ss << m_physicalY;
		_xmlTools->handleData(true, _ss.str(), "DataSize","PhysicalY");
	}

	{
		std::stringstream _ss;
		_ss << m_physicalZ;
		_xmlTools->handleData(true, _ss.str(), "DataSize","PhysicalZ");
	}

	{
		std::stringstream _ss;
		_ss << m_calibratinMatrixXML;
		_xmlTools->handleData(true,_ss.str(),"Calibration");
	}

	{
		std::stringstream _ss;
		_ss << m_masterFrameXML;
		_xmlTools->handleData(true,_ss.str(), "ReconstructionInfo","MasterFrame");
	}

	{
		std::stringstream _ss;
		for (int i=0;i<m_imageFilesXML.size();i++)
		_ss << m_imageFilesXML.at(i) << "\r\n";
		_xmlTools->handleData(true,_ss.str(), "ReconstructionFiles");
	}

	{
		std::stringstream _ss;

		QString _tmp = (QLocale::system().toString(QDate::currentDate(), QLocale::LongFormat));
		_ss <<  _tmp.toUtf8().constData();
		_xmlTools->handleData(true,_ss.str(), "ReconstructionInfo","Date");
	}

	{
		std::stringstream _ss;

		QString _tmp = QLocale::system().toString(QTime::currentTime(), QLocale::LongFormat);
		_ss <<  _tmp.toUtf8().constData();
		_xmlTools->handleData(true,_ss.str(), "ReconstructionInfo","Time");
	}


	// read in the output filename from the gui
	std::string _temp =  ui.outputRFFileName->text(). toLatin1().data();//m_rfDataFileName. toLatin1().data();
	std::string _path = _temp.substr( 0, _temp.find_last_of( '/' ) +1 );

	_temp = ui.outputRFFileName->text(). toLatin1().data();
	_filename = _temp.substr( 0, _temp.find_last_of( '.xml' )-3 ).c_str(),
		_xmlTools->handleData(true,_filename+".3drf",RF_RECONSTRUCTIONFILE_TAG);

	_xmlTools->handleData(true,_filename+".3dix",RF_INDEXFILE_TAG);
	_xmlTools->handleData(true,_filename+".3ddv",RF_DIRECTIONFILE_TAG);

	// rename the tempory filename to the final filename
	rename((_path+"temporaryDirectionData.3ddv").c_str(),(_path+_filename+".3ddv").c_str());

	_xmlTools->saveToFile(_path+_filename+".xml");
	_xmlTools->closeHandle();
	delete _xmlTools;

	// now write out the transformation file
	{
		std::ofstream _transfile;
		_transfile.open(_path+_filename+".txt",std::ios::out);
		_transfile << m_boundingBoxTransformation.matrix();
		_transfile.close();
	}

	_mergedFile.open(_path+_filename+".3drf",std::ios::out  | std::ios::binary);
	_indexFile.open(_path+_filename+".3dix",std::ios::out  | std::ios::binary);

	int _volumeIndexStart[] = {0,2,4,6};
	int _volWriteCounter[] = {0,0,0,0,0,0,0,0};
	int _volumeIndexEnd[] = {1,3,5,7};

	ulong64 _byteCounter;
	_byteCounter.value = (unsigned long long)0;


	int _percentCounter = 0;
	int _voxels = m_voxelsZ * m_voxelsX * m_voxelsY;
	int _1percent = (int)_voxels/100.0;

	{
		m_percentage = -1;
		updatePercentage();	
		std::stringstream _ss;
		_ss << "Merging sub volumes ... ";
		ui.textEdit->append(QString(_ss.str().c_str()));
		qApp->processEvents();
	}

#if DEBUG_TXT_OUTPUT
	std::ofstream _debugFile;
	_debugFile.open(_filename+".mergeinfo.txt");
#endif

	int _elementsCounter = 0;
	for(int _index = 0; _index <= 2; _index+=2)

	{


		for (int c=0; c < ceil(m_voxelsZ/2.);c++)
		{

			// lower part


			for (int _low = 0; _low < ceil(m_voxelsY/2.0);_low++)
			{

				for (int _vol=_volumeIndexStart[_index]; _vol <= _volumeIndexEnd[_index]; _vol++) // 0,1 int first part, then 4,5
				{


					// read rows in sub volume 0 and 1

					for (int _row = 0; _row < ceil(m_voxelsX/2.0); _row++)
					{

						_elementsCounter++;

						_percentCounter++;

						if ( _percentCounter == _1percent)
						{
							updatePercentage();
							_percentCounter = 0;
						}

						// first value indicating the number of elements
						ushort16 _ushort16;
						_rfSubDataFile[_vol].read(&_ushort16.byte.c[0], sizeof(ushort16));
						_mergedFile.write(&_ushort16.byte.c[0], sizeof(ushort16));

						_byteCounter.value += sizeof(ushort16);

					


#if DEBUG_TXT_OUTPUT
						_debugFile << _elementsCounter << "\t" << _ushort16.value  << "\t" << "[ "<< _vol <<" ]" << "\t" << _volWriteCounter[_vol]++ << "\t"; // << std::endl;
#endif

						// now read intensity and distance
						for(int i=0;i<_ushort16.value;i++)
						{



							ushort16 _ushort;
							_rfSubDataFile[_vol].read(&_ushort.byte.c[0], sizeof(ushort16));
							_mergedFile.write(&_ushort.byte.c[0], sizeof(ushort16));

							_byteCounter.value += sizeof(ushort16);

							#if DEBUG_TXT_OUTPUT
							_debugFile << _ushort.value << "\t" ;
#endif

							// handling the direction index
							{
								uint32 _uint32;
								_rfSubDataFile[_vol].read(&_uint32.byte.c[0], sizeof(uint32));
								_mergedFile.write(&_uint32.byte.c[0], sizeof(uint32));
								_byteCounter.value += sizeof(uint32);
							}

							for (int j=0; j<_ushort.value;j++)
							{
								IntensityAndDistance _iad;
								//ushort16 _ushort16;
								
								ushort16 _scanlineIndex;

								/*_rfSubDataFile[_vol].read(&_ushort16.byte.c[0], sizeof(ushort16));
								_mergedFile.write(&_ushort16.byte.c[0], sizeof(ushort16));
								_byteCounter.value += sizeof(ushort16);*/
								#ifdef SAVE_DISTANCE_SCANLINE
#ifdef  QUANTIZE_DISTANCE
								char _distance;
								_rfSubDataFile[_vol].read(&_distance, sizeof(char));
								_mergedFile.write(&_distance, sizeof(char));
								_byteCounter.value += sizeof(char);
#else

								float32 _distance;
								_rfSubDataFile[_vol].read(&_distance.byte.c[0], sizeof(float32));
								_mergedFile.write(&_distance.byte.c[0], sizeof(float32));
								_byteCounter.value += sizeof(float32);
#endif
#endif

#ifdef USE_USHORT16
								ushort16 _intensity16;

								_rfSubDataFile[_vol].read(&_intensity16.byte.c[0], sizeof(ushort16));
								_mergedFile.write(&_intensity16.byte.c[0], sizeof(ushort16));
								_byteCounter.value += sizeof(ushort16);
#endif

#ifdef USE_FLOAT32
								float32 _intensity32;
								_rfSubDataFile[_vol].read(&_intensity32.byte.c[0], sizeof(float32));
								_mergedFile.write(&_intensity32.byte.c[0], sizeof(float32));
								_byteCounter.value += sizeof(float32);
#endif

#ifdef USE_FLOAT5X32
								for(int k=0;k<5;k++) {
								float32 _intensity32;
								_rfSubDataFile[_vol].read(&_intensity32.byte.c[0], sizeof(float32));
								_mergedFile.write(&_intensity32.byte.c[0], sizeof(float32));
								_byteCounter.value += sizeof(float32);
								}
#endif

								#ifdef SAVE_DISTANCE_SCANLINE
								_rfSubDataFile[_vol].read(&_scanlineIndex.byte.c[0], sizeof(ushort16));
								_mergedFile.write(&_scanlineIndex.byte.c[0], sizeof(ushort16));
								_byteCounter.value += sizeof(ushort16);
#endif
							}

						}
						#if DEBUG_TXT_OUTPUT
							_debugFile << std::endl;
#endif
						// write index
						_indexFile.write(&_byteCounter.byte.c[0], sizeof(ulong64));
						__int64 _tmpPos = _mergedFile.tellg();
						int _debug = 5;
					}
				}

			}


			// upper part

			for (int _hi = 0; _hi < ceil(m_voxelsY/2.0);_hi++)
			{

				for (int _vol=_volumeIndexStart[_index+1]; _vol <= _volumeIndexEnd[_index+1]; _vol++) // 2,3 in first and 6,7 second
				{
					// read rows in sub volume 0 and 1
					for (int _row = 0; _row < ceil(m_voxelsX/2.0); _row++)
					{
						_elementsCounter++;
						_percentCounter++;

						if ( _percentCounter == _1percent)
						{
							updatePercentage();
							_percentCounter = 0;
						}

						// first value indicating the number of elements
						ushort16 _ushort16;
						_rfSubDataFile[_vol].read(&_ushort16.byte.c[0], sizeof(ushort16));
						_mergedFile.write(&_ushort16.byte.c[0], sizeof(ushort16));
						_byteCounter.value += sizeof(ushort16);

				

#if DEBUG_TXT_OUTPUT
						_debugFile << _elementsCounter << "\t" << _ushort16.value  << "\t" << "[ "<< _vol <<" ]" << "\t" << _volWriteCounter[_vol]++ << "\t"; // << std::endl;
#endif

						// now read intensity and distance
						for(int i=0;i<_ushort16.value;i++)
						{

							if ( i+1 == 3464 && _elementsCounter == 4701076 && _vol == 2 )
								std::cout << "debug" << std::endl;
							
							ushort16 _ushort;
							_rfSubDataFile[_vol].read(&_ushort.byte.c[0], sizeof(ushort16));
							_mergedFile.write(&_ushort.byte.c[0], sizeof(ushort16));
							_byteCounter.value += sizeof(ushort16);

							#if DEBUG_TXT_OUTPUT
							_debugFile << _ushort.value << "\t" ;
#endif


							// handling the direction index

							{
								uint32 _uint32;
								_rfSubDataFile[_vol].read(&_uint32.byte.c[0], sizeof(uint32));
								_mergedFile.write(&_uint32.byte.c[0], sizeof(uint32));
								_byteCounter.value += sizeof(uint32);
							}

							for (int j=0; j<_ushort.value;j++)
							{
								IntensityAndDistance _iad;
								//ushort16 _ushort16;
								
								//ushort16 _intensity16;
								ushort16 _scanlineIndex;

								/*_rfSubDataFile[_vol].read(&_ushort16.byte.c[0], sizeof(ushort16));
								_mergedFile.write(&_ushort16.byte.c[0], sizeof(ushort16));
								_byteCounter.value += sizeof(ushort16);*/
								#ifdef SAVE_DISTANCE_SCANLINE
#ifdef  QUANTIZE_DISTANCE
								char _distance;
								_rfSubDataFile[_vol].read(&_distance, sizeof(char));
								_mergedFile.write(&_distance, sizeof(char));
								_byteCounter.value += sizeof(char);
#else
								float32 _distance;
								_rfSubDataFile[_vol].read(&_distance.byte.c[0], sizeof(float32));
								_mergedFile.write(&_distance.byte.c[0], sizeof(float32));
								_byteCounter.value += sizeof(float32);
#endif
#endif


								#ifdef USE_USHORT16
								ushort16 _intensity16;

								_rfSubDataFile[_vol].read(&_intensity16.byte.c[0], sizeof(ushort16));
								_mergedFile.write(&_intensity16.byte.c[0], sizeof(ushort16));
								_byteCounter.value += sizeof(ushort16);
#endif

#ifdef USE_FLOAT32
								float32 _intensity32;
								_rfSubDataFile[_vol].read(&_intensity32.byte.c[0], sizeof(float32));
								_mergedFile.write(&_intensity32.byte.c[0], sizeof(float32));
								_byteCounter.value += sizeof(float32);
#endif

#ifdef USE_FLOAT5X32
								for(int k=0;k<5;k++) {
								float32 _intensity32;
								_rfSubDataFile[_vol].read(&_intensity32.byte.c[0], sizeof(float32));
								_mergedFile.write(&_intensity32.byte.c[0], sizeof(float32));
								_byteCounter.value += sizeof(float32);
								}
#endif

								#ifdef SAVE_DISTANCE_SCANLINE
								_rfSubDataFile[_vol].read(&_scanlineIndex.byte.c[0], sizeof(ushort16));
								_mergedFile.write(&_scanlineIndex.byte.c[0], sizeof(ushort16));
								_byteCounter.value += sizeof(ushort16);
#endif
							}
						}
						#if DEBUG_TXT_OUTPUT
							_debugFile << std::endl;
#endif
						// write index
						_indexFile.write(&_byteCounter.byte.c[0], sizeof(ulong64));
						__int64 _tmpPos = _mergedFile.tellg();
						int debug = 5;
					}
				}
			}



		}


	}

#if DEBUG_TXT_OUTPUT
	_debugFile.close();
#endif
	_mergedFile.close();
	_indexFile.close();

	for(int i=0;i<8;i++)
	{
		_rfSubDataFile[i].close();
	}

	_percentCounter = -1;
	updatePercentage();



	{
		clock_t _stopTimer = clock();
		std::stringstream _ss;

		int _iHours =(int)floor((float)(_stopTimer - _startTimer)/(CLOCKS_PER_SEC*60*60));
		int _iMinutes = (int)floor((float)(_stopTimer - _startTimer))/(CLOCKS_PER_SEC*60)-_iHours*60;
		int _iSeconds = ((float)(_stopTimer - _startTimer))/CLOCKS_PER_SEC - _iMinutes*60 - _iHours*60*60;

		char _cHours[3], _cMinutes[3], _cSeconds[3];
		sprintf(_cHours,"%d",_iHours);
		if ( strlen(_cHours) < 2)
			sprintf(_cHours,"0%d",_iHours);
		sprintf(_cMinutes,"%d",_iMinutes);
		if ( strlen(_cMinutes) < 2)
			sprintf(_cMinutes,"0%d",_iMinutes);
		sprintf(_cSeconds,"%d",_iSeconds);
		if ( strlen(_cSeconds) < 2)
			sprintf(_cSeconds,"0%d",_iSeconds);

		// the last percentage update might have got lost
		if ( m_percentage < 100.0)
			updatePercentage();

		_ss.str("");
		_ss << "Compounding finished after " << _cHours << ":" << _cMinutes << ":" << _cSeconds;
		ui.textEdit->append(_ss.str().c_str());
	}

	// now remove the temporary files (subvolumes)
	for(int i=0;i<8;i++)
	{
		remove(_subVolumeFileNames[i].c_str());

	}

	glWidget->restartTimer();
	m_slicesVisualization->restartTimer();

	// free the XML data handler
	if ( m_xmlTools ) 
	{
		delete m_xmlTools;
		m_xmlTools = 0;
	}

	return true;

}

void Compounding::computeScanlines(unsigned char settingsIndex, unsigned char fileIndex)
{
	// 
	//std::string _test = ui.outputRFFileName. toLatin1().data();
	std::string _temp = ui.outputRFFileName->text(). toLatin1().data();
	std::string _filename = _temp.substr( 0, _temp.find_last_of( '.xml' )-3 ).c_str();
	std::string _path = _temp.substr( 0, _temp.find_last_of( '/' ) +1 );
	_filename = "temporaryDirectionData.3ddv";

	std::fstream m_rfDirectionFile;
	if ( (int)fileIndex == 0 ) // append data to the file
		{
			remove((_path+_filename).c_str());
		m_rfDirectionFile.open(_path+_filename, std::ios::out  | std::ios::binary );
		}
	else // overwrite
		m_rfDirectionFile.open(_path+_filename, std::ios::out  | std::ios::binary | std::ios::app);


	int _rfFrames = m_rfDataVector.size();

	// don't compute old scanlines again, add them if there is already previous data, so adapt offset
	int _offset = 0;
	if ( m_indexVector.size() > 1 )
		_offset = m_indexVector[m_indexVector.size()-2];


	//std::remove("scanlines.txt");

	for(int i=_offset; i<_rfFrames; i++)
	{
		// data conversion
#ifdef USE_USHORT16
		StandardRF<unsigned short> *_rf = reinterpret_cast<StandardRF<unsigned short>*>(m_rfDataVector[i]);
#endif
#ifdef USE_FLOAT32
		StandardRF<float> *_rf = reinterpret_cast<StandardRF<float>*>(m_rfDataVector[i]);
#endif
#ifdef USE_FLOAT5X32
		CompositeRF<float> *_rf = reinterpret_cast<CompositeRF<float>*>(m_rfDataVector[i]);
#endif

		_rf->computeScanlines(m_rfDirectionFile, m_scanlineDataVector, m_frameScanlineRangeVector, fileIndex, m_rfScanlineIndex, i - _offset, &m_ultrasoundSettings, settingsIndex);
		//m_scanlineDataVector;
	}

	m_rfDirectionFile.close();
}

bool Compounding::compoundVolume(COMPOUNDING_MODE _compoundingMode)
{
	// future use, not yet implemented
	bool _adaptiveRange = (ui.adaptiveRangeCheckbox->checkState() == Qt::Checked);

	// Create the MHD image and set its parameters

	// TJK - BEGIN: 17-03-10 change from CAMP to VTK
	// 
	// CAMP::Image *compounding = new CAMP::ImageBase<unsigned char>;
	// compounding->init( m_voxelsX, m_voxelsY, m_voxelsZ, 1,1);
	// compounding->setPhysicalSize(m_physicalX, m_physicalY, m_physicalZ, 1.0);
	// TJK - END

	// legacy code VTK 5.0
	/*
	vtkImageData *_compoundingOutput = vtkImageData::New();
	_compoundingOutput->SetScalarTypeToUnsignedChar(); 
	_compoundingOutput->SetDimensions(m_voxelsX, m_voxelsY, m_voxelsZ); 
	_compoundingOutput->SetOrigin(0.0, 0.0, 0.0); 
	_compoundingOutput->SetSpacing(m_physicalX, m_physicalY, m_physicalZ); 
	_compoundingOutput->Update(); 
*/

	// VTK 6.0
	vtkSmartPointer<vtkImageData> _compoundingOutput = vtkSmartPointer<vtkImageData>::New();
	// Setup the image
	_compoundingOutput->AllocateScalars(VTK_UNSIGNED_CHAR, 1);
	_compoundingOutput->SetDimensions(m_voxelsX, m_voxelsY, m_voxelsZ);
	_compoundingOutput->SetOrigin(0.0, 0.0, 0.0);
	_compoundingOutput->SetSpacing(m_physicalX, m_physicalY, m_physicalZ);
	//_compoundingOutput->Update();

	// clear the memory
	COMPOUNDING_DATA_TYPE *_tmpIntensityPtr = static_cast<COMPOUNDING_DATA_TYPE*>(_compoundingOutput->GetScalarPointer(0,0,0));

	for(int _z=0;_z<m_voxelsZ;_z++)
		for(int _y=0;_y<m_voxelsY;_y++)
			for(int _x=0;_x<m_voxelsX;_x++)
			{
				*_tmpIntensityPtr = 0;
				_tmpIntensityPtr++;
			}


			/*
			// determine the maximum dimension of the voxel spacing
			float _maxDim = m_voxelSize;
			float _voxelDiameter = sqrt(m_physicalX*m_physicalX+m_physicalY*m_physicalY+m_physicalZ*m_physicalZ)/2.0f;

			// initialize the rotation queue
			RotationQueue _rotQueue(_maxDim,m_extent.dimX+m_extent.dimY+m_extent.dimZ, m_videoTools->getNumberOfFrames());

			Eigen::Matrix4f _scalingMatrix;
			_scalingMatrix <<	m_scaleX,0,0,0, 
			0,m_scaleY,0,0, 
			0,0,1,0, 
			0,0,0,1;

			// now we set up a data structure containing information about the position of the ultrasound planes in space
			for (int i=0; i<m_matrixDataTracking.size();i++)
			{
			USFrame* _tempFrame = new USFrame();
			// the tracking data specifies where the plane lies in space, here the origin is assumed in the upper left corner of the image
			Plane *_plane = new Plane(m_matrixDataTracking[i]);
			_tempFrame->distance = (float)0.0;
			// _matrixDataTracking.at(x) transforms from unit [mm] on the image to world [mm]
			// but we typically need it for pixels, thus we multiply it with the scaling factors

			//TJK_27_07_09: m_dataTools.eigen2CAMP(m_matrixDataTracking[i] * _scalingMatrix, _tempFrame->image2World);
			_tempFrame->image2World = (m_matrixDataTracking[i] * _scalingMatrix);

			// compute all the relevant information describing the ultrasound plane
			//TJK_27_07_09: _tempFrame->world2Image = _tempFrame->image2World.getInverse();
			_tempFrame->world2Image = _tempFrame->image2World.inverse();
			//_tempFrame->image2WorldEigen = _matrixDataTracking.at(i) * _scalingMatrix;
			//_tempFrame->world2ImageEigen = _tempFrame->image2WorldEigen.inverse();
			_tempFrame->vectorIndex = (int)floor((float)i/NUMBER_OF_SLICES); 
			_tempFrame->videoFrame = i % (int)(NUMBER_OF_SLICES);
			_tempFrame->plane = _plane;

			// set the ultrasound frame at initial distance 0.0
			_rotQueue.push(_tempFrame,0.0);
			}

			int test = m_matrixDataTracking.size();

			RadiusSelect _radiusSelect(m_voxelSize, std::max(m_scaleX, m_scaleY), 50);


			// choose the distance threshold, that is the maximum distance a US slice can have
			// in order to be considered for contributing intensity to a voxel

			float _distanceScalar = 0.0f;
			switch(ui.maxDistanceScalarCombo->currentIndex())
			{
			case 0: { _distanceScalar = 0.5f; break;}
			default: { _distanceScalar = (float)(ui.maxDistanceScalarCombo->currentIndex()); }
			};

			float _threshold = 0.0f;
			switch(ui.maxDistanceUnitCombo->currentIndex())
			{
			case 0: { _threshold = _maxDim; break;}
			case 1: { _threshold = _voxelDiameter; break;}
			};

			// now we compute the actual threshold value from the UI selection
			_threshold = _threshold * _distanceScalar;

			// precompute the squared distance to speed-up the computation
			//float _thresholdSQR = pow(_threshold,2.0f);
			stdext::hash_map<int, std::vector<Segment>*> _hashMap;
			stdext::hash_map<int, std::vector<Segment>*>::iterator _hashMapIterator;


			clock_t _startTimer = clock();
			ui.textEdit->append(QString("Compounding started"));
			qApp->processEvents();

			TraversalScheme _traversalScheme(m_voxelsX, m_voxelsY, m_voxelsZ);
			int _x,_y,_z;
			int _traversalCount = 0;

			// Gaussian backward
			// sigma is only needed as squared, so precalculate to save computation time
			float _sigmaSQR = pow(atof(ui.gaussianSigmaEdit->text(). toLatin1().data()),2.0);

			// Inverse distance weighted smoothness factor
			float _inverseWeightedSmoothness = atof(ui.inverseDistanceFactorEdit->text(). toLatin1().data());

			Eigen::Vector3f _worldIntersection;


			GraphicTools _graphicTools;

			// traverse through the volume

			ui.progressBar->setRange(0,100);
			// as to avoid too frequen updates of the GUI we do only updates when values change
			int _maxVoxels =  m_voxelsX*m_voxelsY*m_voxelsZ;
			int _1percent = (int)_maxVoxels/100.0;
			int _percentCounter = 0;
			int _percent = 0;
			ui.progressBar->setValue(_percent);
			qApp->processEvents();
			while(_traversalScheme.nextPosition(_x,_y,_z))
			{
			_percentCounter++;
			_traversalCount++;

			//  get the memory addess where the intensity to be set in the MHD data is located
			COMPOUNDING_DATA_TYPE *_outputImageIntensityPtr = static_cast<COMPOUNDING_DATA_TYPE*>(_compoundingOutput->GetScalarPointer(_x, _y, _z));

			if ( _percentCounter == _1percent)
			{
			_percent++;
			_percentCounter = 0;
			ui.progressBar->setValue(_percent);
			char buffer[50];
			sprintf(buffer,"TASScomp - %d%%",_percent);
			this->setWindowTitle(QString(buffer));
			qApp->processEvents();
			}


			// compute the world coordinates from the voxel coordinates
			const Eigen::Vector3f _local(static_cast<float>(_x)*m_physicalX+m_physicalX/2.0,static_cast<float>(_y)*m_physicalY+m_physicalY/2.0,static_cast<float>(_z)*m_physicalZ+m_physicalZ/2.0);

			const Eigen::Vector3f _world = (m_boundingBoxTransformation*_local);



			// initalize the values which determine the color of the voxel (gaussian backward warping, inverse distance weighting)
			float _intensityNominator = 0.0;
			float _intensityDenominator = 0.0;


			ADAPTIVE_REPEAT:

			int _numberOfVoxelContributors = 0;

			// pop the minimum distance bin vector from the rotation queue
			std::vector<USFrame*>* _usFrameQueue = _rotQueue.pop();
			int _oldFrameSize = _usFrameQueue->size();

			// this vector is needed for determining the mean color
			std::vector<float> _intensityVector;

			std::vector<WeightedMedianData> _weightedMedianVector;
			float _weightedIntensitySum = 0.0f;

			// this variables are needed for nearest neighbor
			float _minDistance = std::numeric_limits<float>::max();
			float _minDistanceIntensity;

			// for each slice within the vector returned form the queue do a backward warping
			while( _usFrameQueue->size() > 0 )
			{

			USFrame *_frame = _usFrameQueue->back();

			_usFrameQueue->pop_back();

			Plane *_plane = _frame->plane;


			// determine the orthogonal distance from the voxel to the ultrasound plane
			const float _distance = _plane->getOrthogonalDistance (_world, _worldIntersection );

			// determine the point of minimum distance to the voxel on the plane (in pixel coordinates)
			const Eigen::Vector3f _localIntersection = (_frame->world2Image*_worldIntersection);

			// point is relevant for the voxel, its intensity has to be considered / it has to be within the image
			// NOTE: the _roiMin/Max is because of the use of the cropped image area within the video frame
			if ( _distance < _threshold && _localIntersection[0] >= m_roiMinX && _localIntersection[0] <= m_roiMaxX && _localIntersection[1] >= m_roiMinY && _localIntersection[1] <= m_roiMaxY)
			{

			// update the distance
			_frame->distance = _distance;

			// and push it into the temporary queue, as the point is not to be considered again for the same voxel
			_rotQueue.push(_frame, _distance);




			// grab the video image

			// very slow code
			//HBITMAP bmp;
			//m_videoSegmentationTools.getFrame(_frame.videoFrame, bmp);
			//float _tempIntensity = m_videoSegmentationTools.getIntensity(bmp, (int)_localIntersection(0), (int)_localIntersection(1));




			// now compute the radius in the plane that is within distance specified by the threshold
			// that's simply pythagoras 
			//int _radius = static_cast<int>(sqrt ( _distance * _distance + _thresholdSQR ) / _scaleX);
			// instead of pythagoras a faster version is to take the precomputed distance values

			int _radius = 0; 
			switch(_compoundingMode )
			{
			case COMPOUNDING_MODE::GAUSSIAN: 
			{
			// determine the pixel area contributing to backward-warping
			_radius = _radiusSelect.findRadius(_distance);
			break;
			}
			case COMPOUNDING_MODE::INVERSE_DISTANCE:
			{
			// determine the pixel area contributing to backward-warping
			_radius = _radiusSelect.findRadius(_distance);
			break;
			}
			case COMPOUNDING_MODE::NEAREST_NEIGHBOR:
			{
			// assume _radius = 0, so on the slice, only the pixel with min.
			// orthogonal distance is selected
			break;
			}
			case COMPOUNDING_MODE::MEDIAN:
			{
			// assume _radius = 0, so on the slice, only the pixel with min.
			// orthogonal distance is selected
			break;
			}
			case COMPOUNDING_MODE::WEIGHTED_MEDIAN:
			{
			// assume _radius = 0, so on the slice, only the pixel with min.
			// orthogonal distance is selected
			break;
			}
			};

			std::vector<Segment> *_SegmentVector;

			// check the hash map if for the radius there exist already Segment information
			// the Segment contains the discretization of the circle for a given radius
			// e.g. K pixels above the center, move from -M to +M (relative position)
			// it is like the active edge table in the rasterizer
			_hashMapIterator = _hashMap.find(_radius);


			// in case there is no Segment information available, we have to generate it
			if ( _hashMapIterator == _hashMap.end() )
			{
			_SegmentVector = new std::vector<Segment>;

			// compute the relative position Segments of the circle
			// rasterization table for a circle of given radius (relative coordinates => 0,0 as first parameters)
			_graphicTools.bresenhamCircle(0,0, _radius, _SegmentVector);

			// place the information in the hashmap
			_hashMap[_radius]=_SegmentVector;
			} // Segment information is available, grab it
			else
			{
			// discretization was computed previously, so we simply have to acquire the pointer
			_SegmentVector = _hashMapIterator->second;
			}

			// now we look up follow the information contained in the vector (visit all the pixels, which are in range)
			// cache the size of the vector to save computations
			std::vector<Segment>::iterator _itBegin = _SegmentVector->begin();
			std::vector<Segment>::iterator _itEnd = _SegmentVector->end();

			for(std::vector<Segment>::iterator it = _itBegin; it != _itEnd;++it )
			{
			Segment _Segment = *it;

			// as we have only relative positional information from the mid-point/bresenham
			// we have to transform them in absolute positions on the image
			// we have to adapt the coordinates, as we have cropped the image in memory (that's why using _roiMin*)
			_Segment.Y			+= static_cast<int>(_localIntersection[1])- m_roiMinY;
			_Segment.beginX	+= static_cast<int>(_localIntersection[0])- m_roiMinX;
			_Segment.endX		+= static_cast<int>(_localIntersection[0])- m_roiMinX;
			// only consider points which are within the image
			if ( _Segment.Y < m_roiHeight && _Segment.Y >= 0 )
			{

			int _beginX, _endX;
			bool _draw = true;

			// if the Segment is outside the window, forget about it
			if ( _Segment.beginX >= m_roiWidth )
			_draw = false;
			if ( _Segment.endX < 0 )
			_draw = false;

			// if Segment-range is within the area of interest backward-warp it
			if ( _draw && ((_Segment.beginX	 >= 0 && _Segment.beginX <= m_roiWidth-1) ||( _Segment.endX >= 0 && _Segment.endX <= m_roiWidth-1)) )
			{
			// clip to the boundary, define the rasterization limits (the bounding box of the clipping mask)
			_beginX = std::max( std::min( _Segment.beginX, m_roiWidth-1), 0);
			_endX = std::min(std::max( _Segment.endX, 0 ), m_roiWidth-1 );


			IplImage *_img = m_imageVector[_frame->vectorIndex][(_frame->videoFrame)];
			// get the pointers to the ultrasound image and the cropping mask
			unsigned char *_imgPointer= &(((unsigned char*)(_img->imageData + _Segment.Y*_img->widthStep))[_beginX]);
			unsigned char *_croppingMaskPointer = &(((unsigned char*)(m_croppingMask->imageData + _Segment.Y*m_croppingMask->widthStep))[_beginX]);
			// now pixel by pixel backward warping (given the radius)

			for (int x_pos=_beginX;x_pos <= _endX; x_pos++)
			{

			// now apply the cropping
			// get the intensity of the cropping mask
			const float _maskIntensity = (float)*_croppingMaskPointer++;

			// if the intensity equals white, then we are in non-cropped area
			if ( _maskIntensity > 0.0)
			{
			// here we go: we have a pixel, whose itensity is contributing to the voxel intensity!
			_numberOfVoxelContributors ++;
			// get the intensity of the ultrasound image
			const float _tempIntensity = (float)*_imgPointer++;			

			//cvSaveImage("c:\\Resultant.jpg",m_imageVector.at(_frame->videoFrame));
			// compute the distance for the point => _distance is just the orthogonal distance to the plane


			switch(_compoundingMode )
			{
			case COMPOUNDING_MODE::GAUSSIAN: 
			{
			// Gaussian backward warping

			// determine the position of the pixel in world coordinates (mm)
			const Eigen::Vector3f _temp(x_pos+m_roiMinX, _Segment.Y+m_roiMinY, 0.0);
			const Eigen::Vector3f _worldPoint = (_frame->image2World * _temp);

			// determine squared euclidean distance of the point to the pixel in world coordinates
			const float _euclideanDistanceSQR =  pow(static_cast<float>((_world - _worldPoint).norm()),2.0f);

			// determine the interpolation ratio
			const float _tempRatio = exp(-(_euclideanDistanceSQR)/(_sigmaSQR));

			// add the intensity contribution by weighting distance to the nominator and denominator
			// which after considering all relevant pixels will by division supply the gaussian interpolated intensity
			_intensityDenominator += _tempRatio;
			_intensityNominator += _tempIntensity * _tempRatio; 
			break;
			}
			case COMPOUNDING_MODE::INVERSE_DISTANCE: 
			{
			// Inverse distance backward warping

			// determine the position of the pixel in world coordinates (mm)
			const Eigen::Vector3f _temp(x_pos+m_roiMinX, _Segment.Y+m_roiMinY, 0.0);
			const Eigen::Vector3f _worldPoint = (_frame->image2World * _temp);

			// determine euclidean distance by power of smoothness factor of the point to the pixel in world coordinates
			const float _euclideanDistanceSmoothed =  pow(static_cast<float>((_world - _worldPoint).norm()),2.0f);


			// add the intensity contribution by weighting distance to the nominator and denominator
			// which after considering all relevant pixels will by division (and denominator * N) give the inverse distance measure
			_intensityDenominator += _euclideanDistanceSmoothed;
			_intensityNominator += _tempIntensity * _euclideanDistanceSmoothed; 
			break;
			}
			case COMPOUNDING_MODE::NEAREST_NEIGHBOR:
			{
			if ( _distance < _minDistance )
			{
			_minDistance = _distance;
			_minDistanceIntensity = _tempIntensity;
			}
			break;
			}
			case COMPOUNDING_MODE::MEDIAN:
			{
			_intensityVector.push_back(_tempIntensity);
			break;
			}
			case COMPOUNDING_MODE::WEIGHTED_MEDIAN:
			{
			WeightedMedianData _wmd;

			_wmd.weight = 1.0 - _distance / _threshold;
			_wmd.intensity = _tempIntensity;
			_weightedMedianVector.push_back(_wmd);
			_weightedIntensitySum += _wmd.weight;
			break;
			}
			}
			}
			}
			}
			}

			}

			}
			else // point too far away or not within the region of interest, but distance was updated so we have to update it in the rotation queue
			{
			// update the distance
			_frame->distance = _distance;

			// and push it into the temporary queue, as the point is not to be considered again for the same voxel
			_rotQueue.push(_frame, _distance);
			}

			}
			// instead of delete put empty vectors back to rotation queue, so we save new operations
			_usFrameQueue->reserve(_oldFrameSize);
			_rotQueue.recycle(_usFrameQueue);

			// now set the intensity of the voxel
			// small box for showing the origin and the orientation
			//if ( _x <=25 && _y <=12 && _z <=5 )
			//{
			//	_intensityDenominator = 1.0;
			//	_intensityNominator = 255.0;
			//}



			// calculate the intensity of the voxel given the specified method

			switch(_compoundingMode )
			{
			case COMPOUNDING_MODE::GAUSSIAN: 
			{
			// avoid division by zero
			if ( _intensityDenominator > 0.0 && _intensityNominator > 0.0 ) // avoid division by ZERO
			{
			// CAMP Image replaced
			//compounding->setVoxelFloat( _intensityNominator / _intensityDenominator,_x+_y*m_voxelsX+_z*m_voxelsX*m_voxelsY);
			*_outputImageIntensityPtr = static_cast<COMPOUNDING_DATA_TYPE>( _intensityNominator / _intensityDenominator );
			}
			break;
			}

			case COMPOUNDING_MODE::INVERSE_DISTANCE: 
			{
			// avoid division by zero
			if ( _intensityDenominator > 0.0 && _intensityNominator > 0.0 ) // avoid division by ZERO
			{
			// CAMP Image replaced
			//compounding->setVoxelFloat( _intensityNominator / ((float)(_numberOfVoxelContributors)*_intensityDenominator),_x+_y*m_voxelsX+_z*m_voxelsX*m_voxelsY);
			*_outputImageIntensityPtr = static_cast<COMPOUNDING_DATA_TYPE>( _intensityNominator / ((float)(_numberOfVoxelContributors)*_intensityDenominator) );
			}
			break;
			}

			case COMPOUNDING_MODE::NEAREST_NEIGHBOR:
			{
			// check if there was a slice within range
			if ( _minDistance != std::numeric_limits<float>::max() )
			{
			// CAMP Image replaced
			//compounding->setVoxelFloat( _minDistanceIntensity,_x+_y*m_voxelsX+_z*m_voxelsX*m_voxelsY);
			*_outputImageIntensityPtr = static_cast<COMPOUNDING_DATA_TYPE>( _minDistanceIntensity );
			// reset the minimum distance for the next run
			_minDistance = std::numeric_limits<float>::max();
			}
			break;
			}
			case COMPOUNDING_MODE::MEDIAN:
			{
			int _vecSize = _intensityVector.size();
			if ( _vecSize> 1 )
			{
			// sort the intensities and take the mean
			sort (_intensityVector.begin(), _intensityVector.end());
			// CAMP Image replaced
			//compounding->setVoxelFloat( _intensityVector[round((float)_intensityVector.size()/2.0,2)],_x+_y*m_voxelsX+_z*m_voxelsX*m_voxelsY);
			*_outputImageIntensityPtr = static_cast<COMPOUNDING_DATA_TYPE>( _intensityVector[round((float)_intensityVector.size()/2.0,2)] );
			}
			else if ( _vecSize == 1) 
			{
			// we only have one intensity...no need to sort, just take it
			// CAMP Image replaced
			//compounding->setVoxelFloat( _intensityVector[0],_x+_y*m_voxelsX+_z*m_voxelsX*m_voxelsY);
			*_outputImageIntensityPtr = static_cast<COMPOUNDING_DATA_TYPE>( _intensityVector[0] );
			}
			else
			{
			// do nothing, leave the voxel black
			}
			// reset the intensity vector for the next run
			_intensityVector.clear();
			break;
			}
			case COMPOUNDING_MODE::WEIGHTED_MEDIAN:
			{
			int _vecSize = _weightedMedianVector.size();
			if ( _vecSize > 1 )
			{
			float _tempLast = 0.0f;
			float _half = _weightedIntensitySum / 2.0f;

			for(int i=1;i<_vecSize;i++)
			{
			float _w1 = _weightedMedianVector[i-1].weight+_tempLast;
			float _w2 = _weightedMedianVector[i].weight+_w1;
			_tempLast = _w1;

			if ( _w1 <=  _half && _half <= _w2)
			{
			// CAMP Image replaced
			//compounding->setVoxelFloat( _weightedMedianVector[i].intensity,_x+_y*m_voxelsX+_z*m_voxelsX*m_voxelsY);
			*_outputImageIntensityPtr = static_cast<COMPOUNDING_DATA_TYPE>(  _weightedMedianVector[i].intensity );
			break;
			}
			}
			}
			else if ( _vecSize == 1)
			{
			// CAMP Image replaced
			//compounding->setVoxelFloat( _weightedMedianVector[0].intensity,_x+_y*m_voxelsX+_z*m_voxelsX*m_voxelsY);
			*_outputImageIntensityPtr = static_cast<COMPOUNDING_DATA_TYPE>( _weightedMedianVector[0].intensity );
			break;
			}
			else
			{
			// leave it black
			}
			// reset the intensity vector for the next run
			_weightedMedianVector.clear();
			break;
			}
			}



			} // end of while

			*/

			SYSTEM_INFO _info;
			GetSystemInfo(&_info );

			std::stringstream _ss;



			clock_t _startTimer = clock();
			ui.textEdit->append(QString("Compounding started"));
			qApp->processEvents();

			_ss << _info.dwNumberOfProcessors << " cores detected!";

			int _numOfThreads = ceil(m_threadingFactor *  _info.dwNumberOfProcessors);

			ui.textEdit->append(QString(_ss.str().c_str()));


			m_cvThread = new CompoundVolumeThread[_numOfThreads];

			int _share = (int)floor((float)m_voxelsZ/(float)_numOfThreads);
			int _startPos = 0;
			int _endPos = _share;

			ui.progressBar->setRange(0,100);
			// as to avoid too frequen updates of the GUI we do only updates when values change
			int _maxVoxels =  m_voxelsX*m_voxelsY*m_voxelsZ;
			int _1percent = (int)_maxVoxels/100.0;
			int _percentCounter = 0;
			int _percent = 0;
			ui.progressBar->setValue(_percent);
			qApp->processEvents();

			for (int pro=0; pro < _numOfThreads;pro++)
			{
				_ss.str("");
				_ss << "Start Compounding Thread No. " << pro ;
				ui.textEdit->append(QString(_ss.str().c_str()));
				qApp->processEvents();

				// set the UI input 
				m_cvThread[pro].setUserConfiguration(ui.maxDistanceScalarCombo->currentIndex(), ui.maxDistanceUnitCombo->currentIndex(), atof(ui.gaussianSigmaEdit->text(). toLatin1().data()), atof(ui.inverseDistanceFactorEdit->text(). toLatin1().data()), _1percent);

				m_cvThread[pro].setROI(m_roiMaxX, m_roiMaxY, m_roiMinX, m_roiMinY, m_roiHeight, m_roiWidth);

#ifdef  IMAGE_FILTER
				m_cvThread[pro].setDataPointers( &m_imageVector, &m_filterVector, &m_matrixDataTracking, &m_trackingData);
#else
				m_cvThread[pro].setDataPointers( &m_imageVector, NULL, &m_matrixDataTracking, &m_trackingData);
#endif
				m_cvThread[pro].setAuxiliaryPointers(m_videoTools, m_croppingMask);
				m_cvThread[pro].setTransformation(m_boundingBoxTransformation);
				m_cvThread[pro].setCompoundingMode(_compoundingMode);
				m_cvThread[pro].setTargetImage(_compoundingOutput);
				m_cvThread[pro].setDataSize(m_voxelsX, m_voxelsY, _endPos, m_physicalX, m_physicalY, m_physicalZ, m_extent, 0, 0, _startPos, m_imageCounter);
				m_cvThread[pro].setUSResolution(m_scaleX, m_scaleY);
				//m_cvThread[pro].setPriority(QThread::HighestPriority);

				// calculate voxel section for next thread
				_startPos = _endPos;
				_endPos = _startPos+_share;
				// the last processor will get all the remainaining slices, so make sure we compound to the last voxel
				if ( pro+2 == _numOfThreads)
					_endPos = m_voxelsZ;
			}


			for (int pro=0; pro < _numOfThreads;pro++)
			{
				connect(&m_cvThread[pro], SIGNAL(percentEvent()),
					this, SLOT(updatePercentage()), Qt::QueuedConnection);

				m_cvThread[pro].start();
			}

			bool _finished = false;
			// we have to run the thread for a certain amount of time, then we update the gui
			// if the threads have not finished with wait again for a certain time and so on...until all the threads are finished
			// this way we can update UI elements
			do 
			{

				for (int pro=0; pro < _numOfThreads;pro++)
				{
					m_cvThread[pro].wait(350);
				}
				qApp->processEvents();

				_finished = true;
				for (int pro=0; pro < _numOfThreads;pro++)
				{
					if (!m_cvThread[pro].isFinished())
						_finished = false;
				}
			}
			while(!_finished);
			if ( _numOfThreads  > 1 )
				delete[] m_cvThread;
			else
				delete m_cvThread;


			//TJK_removed_CAMP: CAMP::Matrix4<float> _patient;
			//TJK_removed_CAMP: m_dataTools.eigen2CAMP(m_boundingBoxTransformation.matrix(), _patient);

			// Removed CAMP Image dependency
			//compounding->setImagePosPatient(_patient.getTranslation());
			//compounding->setImageOriPatient(_patient.getRotation());

			//TJK_removed_CAMP: std::cout << "Rotation: " << std::endl << _patient.getRotation() << std::endl;

			//TJK_removed_CAMP: std::cout << "Translation: " << std::endl << _patient.getTranslation() << std::endl;


			// extract the filename without extension
			std::string _temp = m_mhdFileName. toLatin1().data();

			//CAMP::Tools::saveMetaImage(_temp.substr( 0, _temp.find_last_of( '.mhd' )-3 ).c_str(), compounding );

			// save the file


			vtkSmartPointer<vtkMetaImageWriter> _writer = vtkSmartPointer<vtkMetaImageWriter>::New();
			
			vtkSmartPointer<vtkImageCast> _castFilter =
				vtkSmartPointer<vtkImageCast>::New();
			_castFilter->SetOutputScalarTypeToUnsignedChar();
			_castFilter->SetInputData(_compoundingOutput);

			_writer->SetInputConnection(_castFilter->GetOutputPort());
			_writer->SetFileName(ui.outputFileName->text().toLatin1().data());
			_writer->SetCompression(false);
			_writer->Write();

			_writer->Delete();

			



			//delete compounding;

			// now compute the time it took for compounding the ultrasound volume

			clock_t _stopTimer = clock();

			int _iHours =(int)floor((float)(_stopTimer - _startTimer)/(CLOCKS_PER_SEC*60*60));
			int _iMinutes = (int)floor((float)(_stopTimer - _startTimer))/(CLOCKS_PER_SEC*60)-_iHours*60;
			int _iSeconds = ((float)(_stopTimer - _startTimer))/CLOCKS_PER_SEC - _iMinutes*60 - _iHours*60*60;

			char _cHours[3], _cMinutes[3], _cSeconds[3];
			sprintf(_cHours,"%d",_iHours);
			if ( strlen(_cHours) < 2)
				sprintf(_cHours,"0%d",_iHours);
			sprintf(_cMinutes,"%d",_iMinutes);
			if ( strlen(_cMinutes) < 2)
				sprintf(_cMinutes,"0%d",_iMinutes);
			sprintf(_cSeconds,"%d",_iSeconds);
			if ( strlen(_cSeconds) < 2)
				sprintf(_cSeconds,"0%d",_iSeconds);

			// the last percentage update might have got lost
			if ( m_percentage < 100.0)
				updatePercentage();

			_ss.str("");
			_ss << "Compounding finished after " << _cHours << ":" << _cMinutes << ":" << _cSeconds;

			ui.textEdit->append(_ss.str().c_str());

			_ss.str("");

			/*_ss << "Voxels visited: " << _traversalCount << " of " << m_voxelsX * m_voxelsY * m_voxelsZ;
			ui.textEdit->append(_ss.str().c_str());
			*/
			qApp->processEvents();

			return true;
}


bool Compounding::suggestVoxelSize()
{
	float _niceNumberOfVoxels = 15000000.0;

	m_voxelSize=ceilf(10.0f*pow(((m_extent.dimX*m_extent.dimY*m_extent.dimZ) / _niceNumberOfVoxels),1.0f/3.0f))/10.0f;

	return true;
}

void Compounding::attachData()
{
	bool _doppler = false;

	// identify the sender
	QPushButton *_tempButton = (QPushButton *)sender();
	if (_tempButton!=NULL)
	if ( strcmp(_tempButton->objectName(). toLatin1().data(),"AttachDopplerButton") == 0 )
	{
		// we want to reconstruct Doppler data, this need to be handled as a special case
		_doppler = true;
	}

	loadTimeStampsFromFile(m_videoFileName, false, true);

	// LEGACY SOLUTION
	m_calibrationMatrixVector.push_back(m_matCalibration);

	float _velocityVariance;
	//Eigen::Matrix3f _velocityNoise;
	// get the noise covariance matrix, when noise
	if ( ui.VelocityNoiseCheckBox->isChecked())
	{


		/*{
		std::stringstream _ss(ui.VelocityCovarianceTextEdit->toPlainText(). toLatin1().data());

		float _m[9];
		_ss >> _m[0] >> _m[1] >> _m[2] >> _m[3] >> _m[4] >> _m[5] >> _m[6] >> _m[7] >> _m[8]; 


		_velocityNoise << _m[0],_m[1],_m[2],_m[3], _m[4],_m[5],_m[6],_m[7],_m[8];
		}
		{
		std::stringstream _ss;
		_ss << _velocityNoise << "\n applying velocity noise... ";
		ui.textEdit->append(QString(_ss.str().c_str()));
		}*/

		_velocityVariance = atof(ui.VelocityNoiseSigma->text(). toLatin1().data());
	}

	// Doppler reconstruction again needs a special handling of the video data
	if (!_doppler)
		loadVideoFromFile(_velocityVariance, STANDARD,  true);
	else // here we use Doppler data
		loadVideoFromFile(_velocityVariance, DOPPLER,  true, true, ui.VelocityNoiseCheckBox->isChecked()); // also load the auxiliary file containing the segmentation

	loadTrackingDataFromFile(m_trackingFileName, false, ui.referenceTargetCheckBox->isChecked(), true);
	m_slicesVisualization->setTrackingData(m_trackingData);


	{
		std::stringstream _ss;
		_ss << m_videoFileName. toLatin1().data() << " attached successfully ... ";
		ui.textEdit->append(QString(_ss.str().c_str()));
	}

	if ( ui.MasterFrameBModeRadio->isChecked())
		determineBoundingBoxFromMasterFrame();
	else if ( ui.ManualBModeRadio->isChecked() )
		determineMinimumBoundingBoxManualFromBMode();
	else
		determineMinimumBoundingBoxFromBMode();

	m_slicesVisualization->defineOrientedBoundingBox(m_extent, m_boundingBoxTransformation);



	suggestVoxelSize();
	ui.voxelSizeSlider->setEnabled(true);
	ui.voxelSizeSlider->setValue(m_voxelSize/0.025);
	determineVoxelSize();

	// make the image slider fit the new data size
	ui.ImageLabelScrollBar->setMaximum(m_imageCounter-1);
	ui.dataScrollBar->setMaximum(m_imageCounter-1);

	// now we allow the editing of the number
	ui.MasterFrameEdit2->setEnabled(true);
}

void Compounding::preprocessDoppler()
{
	computeMeasurementCollection();
	computeInitialVectorField();
	ui.DopplerCompoundButton->setEnabled(true);
	ui.DopplerRestoreButton->setEnabled(true);
	ui.DopplerAddNoiseButton->setEnabled(true);
	m_recon->createFieldBackup();
}

void Compounding::loadData()
{
	if ( m_xmlTools ) delete m_xmlTools;
	m_xmlTools = new XMLTools();

	// clear the name items of the tracking target combo box
	ui.probeNameComboBox->clear();
	ui.referenceNameComboBox->clear();


	std::vector<std::string> _targetNames;
	// use "ReferenceTarget" as the moving reference coordinate system, to cope with patient motion 
	bool _foundReference = false;
	bool _foundProbe = false;
	if (m_dataTools.getTargetNamesFromFile(m_trackingFileName. toLatin1().data(), _targetNames) )
	{
		for(int i=0;i<_targetNames.size();i++)
		{
			ui.probeNameComboBox->addItem(QString(_targetNames[i].c_str()));
			if (_targetNames[i].find("Probe") != std::string::npos  || _targetNames[i].find("Transducer") != std::string::npos )
			{
				_foundProbe = true;
				ui.probeNameComboBox->setCurrentIndex(ui.probeNameComboBox->count()-1);
			}
			ui.referenceNameComboBox->addItem(QString(_targetNames[i].c_str()));
			if (_targetNames[i].find("Reference") != std::string::npos  || _targetNames[i].find("reference") != std::string::npos )
			{
				_foundReference = true;
				ui.referenceNameComboBox->setCurrentIndex(ui.referenceNameComboBox->count()-1);
			}
		}
	}

	// if we found a reference target ID, let's assume we use it
	if ( _foundReference == true)
		ui.referenceTargetCheckBox->setChecked(true);
	else
		ui.referenceTargetCheckBox->setChecked(false);


	bool _doppler = false;

	// identify the sender
	QPushButton *_tempButton = (QPushButton *)sender();
	if ( _tempButton != NULL )
	if ( strcmp(_tempButton->objectName(). toLatin1().data(),"LoadDopplerButton") == 0 )
	{
		// we want to reconstruct Doppler data, this need to be handled as a special case
		_doppler = true;
	}

	if ( !_doppler )
		ui.initializeButton->setDisabled(false);
	else
		ui.InitializeDopplerButton->setDisabled(false);
}

void Compounding::setoutputFileName()
{
	// initialisation of the output values
	//QString movieFilter  = tr("MHD files (*.mhd)");		
	QStringList filterList;// = ((QStringList)movieFilter);
	filterList << "*.mhd";

	// creation of a new file dialog
	QFileDialog* fileDialog = new QFileDialog(this);
	fileDialog->setWindowTitle("Save Compounding Data...");
	fileDialog->setNameFilters(filterList);

	// check if we previously loaded some data from a specific directory, if not use the one from the XML file
	if ( m_lastPath.length() == 0 )
	{
		std::string _path;
		m_xmlConfiguration.handleData(false,_path, std::string(TARGET_DIRECTORY_TAG));
		m_lastTargetPath = QString(_path.c_str());
	}

	QString fileName = fileDialog->getSaveFileName(0,"Save file ...", m_lastTargetPath ,"MHD (*.mhd)");
	std::string _temp = fileName. toLatin1().data();
	std::string _path = _temp.substr( 0, _temp.find_last_of( '/' ) +1 );
	m_guiVideoReady = !fileName.isEmpty();
	if(!fileName.isEmpty())
	{
		// write the directory used last time to the configuration file to remember it for the next time
		m_xmlConfiguration.handleData(true,_path, std::string(TARGET_DIRECTORY_TAG));
		m_lastTargetPath = QString(_path.c_str());
		ui.outputFileName->setText(fileName);
		m_mhdFileName = fileName;
	}
}

void Compounding::useReferenceFrame(int state)
{
	if ( ui.referenceTargetCheckBox->isChecked())
	{
		ui.referenceNameComboBox->setEnabled(true);
	}
	else
	{
		ui.referenceNameComboBox->setEnabled(false);
	}
}

void Compounding::computeInitialVectorField()
{

	// this vector is going to contain the voxel coordinates, of voxels
	// for which no velocity could be reconstructed
	std::vector<Coordinate> _interpolateData;
	// the number affecting the acceptance of reconstruction results
	float _conditionNumberThreshold = 5.0f;

	ui.textEdit->append(QString("Compute initial vector field ..."));
	qApp->processEvents();

	clock_t _startTimer = clock();
	TraversalScheme _traversalScheme(m_voxelsX, m_voxelsY, m_voxelsZ);

	_traversalScheme.init();

	int _solution = 0;
	int _sumSol = 0;
	int _traversalCount = 0;

	srand((unsigned)time(0)); 

	ui.progressBar->setValue(0);
	int _voxelsRemoved = 0;
	int _maxVoxels =  m_voxelsX*m_voxelsY*m_voxelsZ;
	int _1percent = (int)_maxVoxels/100.0;
	int _percentCounter = 0;
	int _percent = 0;
	int _x, _y, _z;
	float _error =0.0f, _maxError =0.0f;
	float _conditionNumber = 0.0f, _maxConditionNumber =0.0, _minConditionNumber = 9999999999999.0f;
	//while(_traversalScheme.nextPosition(_x,_y,_z))
	//#pragma omp parallel for
	for(int _z=0;_z<m_voxelsZ;_z++)
		for(int _y=0;_y<m_voxelsY;_y++)
			for(int _x=0;_x<m_voxelsX;_x++)
			{
				_percentCounter++;
				_traversalCount++;

				if ( _percentCounter == _1percent)
				{
					_percent++;
					_percentCounter = 0;
					ui.progressBar->setValue(_percent);
					qApp->processEvents();
				}
				Eigen::VectorXf _velocity;
				int _measurements = 0;
				if ( m_recon->approximateVelocityVector(_x,_y,_z, _velocity, _error, _conditionNumber, _measurements) )
				{
					_solution+=1;
					_maxError = std::max(_maxError, _error);
					_minConditionNumber = std::min(_minConditionNumber, _conditionNumber);
					_maxConditionNumber = std::max(_maxConditionNumber, _conditionNumber);

					if ( /*_error > 10.0 ||*/ _conditionNumber >= _conditionNumberThreshold )
					{
						_velocity = VectorXf::Zero(3);

						// store the position of the voxel, such that we can interpolate
						Coordinate _coordinate;
						_coordinate.X = _x;
						_coordinate.Y = _y;
						_coordinate.Z = _z;

						// this vector will contain all the voxels that will be interpolated
						_interpolateData.push_back(_coordinate);
					}
				}
				else // no solution, e.g. no measurement data available 
				{

					_velocity = VectorXf::Zero(3);

					// store the position of the voxel, such that we can interpolate
					Coordinate _coordinate;
					_coordinate.X = _x;
					_coordinate.Y = _y;
					_coordinate.Z = _z;

					// this vector will contain all the voxels that will be interpolated
					if ( _measurements > 0 ) // only interpolate if we have any measurements
						_interpolateData.push_back(_coordinate);

				}

				// int random_integer = (rand()%10)+1; 

				// add some noise
				// if ( random_integer == 5 )
				// now we set the velocity in the velocity volume
				//m_recon->setVelocity(_x,_y,_z,_velocity+_velocity.norm()*0.1*Eigen::Vector3f((float)rand()/((float)RAND_MAX + 1.0),(float)rand()/((float)RAND_MAX + 1.0),/*(float)rand()/((float)RAND_MAX + 1.0)*/0.0), (_measurements > 0));
				//else

				// if we have more than one measurement (pre-segmented), then we can interpolate here

				float _avgVelocity = 0.0f;
				for (int j=0;j<_velocity.rows();j++)
				{
					_avgVelocity += fabs(_velocity[j]);
				}
				_avgVelocity /= (float)_velocity.rows();

				bool _validVoxel = (_measurements > 0 && _avgVelocity > 0.5);
				m_recon->setVelocity(_x,_y,_z,_velocity, _validVoxel);

				if ( !_validVoxel )
					_voxelsRemoved++;


				//std::cout << "Laplacian: " <<_recon->computeLaplacian(_x,_y,_z, m_physicalX, m_physicalY, m_physicalZ).transpose() << std::endl;

				_sumSol+=1;
			}

			std::cout << "Condition numbers (max, min): " << _maxConditionNumber << " / " << _minConditionNumber << std::endl;


			{
				std::stringstream _ss;
				_ss << "Voxels removed due to inconsitency: " << _voxelsRemoved;

				ui.textEdit->append(_ss.str().c_str());

				_ss.str("");
			}


			// now we interpolate trilinearly with the 8 neighbors of the voxels stores in the interpolation vector


			{
				std::stringstream _ss;
				_ss << "Trilinear interpolation (" << _interpolateData.size() << " Elements) ...";
				ui.textEdit->append(_ss.str().c_str());
			}

			_1percent = (int)_interpolateData.size()/100.0;

			_percentCounter = 0;
			_percent = 0;
			for (int i=0;i < _interpolateData.size();i++)
			{
				// update the progress bar
				if ( _percentCounter == _1percent)
				{
					_percent++;
					_percentCounter = 0;
					ui.progressBar->setValue(_percent);
					qApp->processEvents();
				}

				Coordinate _coordinate = _interpolateData[i];
				int _index = 0;

				Eigen::Vector3f _i11 = m_recon->getVelocityVectorByCoordinates(_coordinate.X, _coordinate.Y, _coordinate.Z);
				Eigen::Vector3f _i12 = VectorXf::Zero(3);
				// check if we are still in range, otherwise keep on using the zero vector
				if ( m_recon->getIndex(_coordinate.X, _coordinate.Y, _coordinate.Z+1, _index) )
				{
					// there is valid data, USE it
					_i12 = m_recon->getVelocityVectorByIndex(_index);
				}

				Eigen::Vector3f _i1 = _i11 * 0.5 + _i12 * 0.5;

				//---------

				Eigen::Vector3f _i21 = VectorXf::Zero(3);
				// check if we are still in range, otherwise keep on using the zero vector
				if ( m_recon->getIndex(_coordinate.X, _coordinate.Y+1, _coordinate.Z, _index) )
				{
					// there is valid data, USE it
					_i21 = m_recon->getVelocityVectorByIndex(_index);
				}

				Eigen::Vector3f _i22 = VectorXf::Zero(3);
				// check if we are still in range, otherwise keep on using the zero vector
				if ( m_recon->getIndex(_coordinate.X, _coordinate.Y+1, _coordinate.Z+1, _index) )
				{
					// there is valid data, USE it
					_i22 = m_recon->getVelocityVectorByIndex(_index);
				}

				Eigen::Vector3f _i2 = _i21 * 0.5 + _i22 * 0.5;

				//---------

				Eigen::Vector3f _j11 = VectorXf::Zero(3);
				// check if we are still in range, otherwise keep on using the zero vector
				if ( m_recon->getIndex(_coordinate.X+1, _coordinate.Y, _coordinate.Z, _index) )
				{
					// there is valid data, USE it
					_j11 = m_recon->getVelocityVectorByIndex(_index);
				}

				Eigen::Vector3f _j12 = VectorXf::Zero(3);
				// check if we are still in range, otherwise keep on using the zero vector
				if ( m_recon->getIndex(_coordinate.X+1, _coordinate.Y, _coordinate.Z+1, _index) )
				{
					// there is valid data, USE it
					_j12 = m_recon->getVelocityVectorByIndex(_index);
				}

				Eigen::Vector3f _j1 = _j11 * 0.5 + _j12 * 0.5;

				//-----------


				Eigen::Vector3f _j21 = VectorXf::Zero(3);
				// check if we are still in range, otherwise keep on using the zero vector
				if ( m_recon->getIndex(_coordinate.X+1, _coordinate.Y+1, _coordinate.Z, _index) )
				{
					// there is valid data, USE it
					_j21 = m_recon->getVelocityVectorByIndex(_index);
				}

				Eigen::Vector3f _j22 = VectorXf::Zero(3);
				// check if we are still in range, otherwise keep on using the zero vector
				if ( m_recon->getIndex(_coordinate.X+1, _coordinate.Y+1, _coordinate.Z+1, _index) )
				{
					// there is valid data, USE it
					_j22 = m_recon->getVelocityVectorByIndex(_index);
				}

				Eigen::Vector3f _j2 = _j21 * 0.5 + _j22 * 0.5;

				Eigen::Vector3f _w1 = _i1 * 0.5 + _i2 * 0.5;
				Eigen::Vector3f _w2 = _j1 * 0.5 + _j2 * 0.5;

				Eigen::Vector3f _interpolationResult = _w1 * 0.5 + _w2 * 0.5;


				// store the new velocity, because data will be set in the update step
				_interpolateData[i].Velocity = _interpolationResult;



				// look at the neighbors and weight them equally, if border assume zero velocity vector)

				_percentCounter++;
			}


			// now we finished trilinear interpolation, so we can update the velocities


			{
				std::stringstream _ss;
				_ss << "Updating the velocity field (" << _interpolateData.size() << " Elements) ..." ;
				ui.textEdit->append(_ss.str().c_str());
			}

			_percentCounter = 0;
			_percent = 0;
			for (int i=0;i < _interpolateData.size();i++)
			{
				// update the progress bar
				if ( _percentCounter == _1percent)
				{
					_percent++;
					_percentCounter = 0;
					ui.progressBar->setValue(_percent);
					qApp->processEvents();
				}

				Coordinate _coordinate = _interpolateData[i]; 


				m_recon->setVelocity(_coordinate.X, _coordinate.Y, _coordinate.Z, _interpolateData[i].Velocity, true); 

				_percentCounter++;
			}


			{
				std::stringstream _ss;
				_ss << "Initial vector field computation finished after " << m_dataTools.computeElapsedTime(_startTimer, clock());
				ui.textEdit->append(_ss.str().c_str());
			}

			m_recon->setMaxError(_maxError);
			Eigen::Vector3f* _velocityField = m_recon->getVelocityField();

			m_velocityVisualization->setVelocityField(_velocityField, m_voxelsX, m_voxelsY, m_voxelsZ, m_physicalX, m_physicalY, m_physicalZ);


			// determine which voxels contain real velocity information (for debugging purposes)
			m_slicesVisualization->setReconstruction(m_recon, m_voxelsX, m_voxelsY, m_voxelsZ, m_physicalX, m_physicalY, m_physicalZ);


			m_slicesVisualization->defineOrientedBoundingBox(m_extent, m_boundingBoxTransformation);
}

void Compounding::restoreVectorField()
{
	m_recon->restoreFieldBackup();
	Eigen::Vector3f* _velocityField = m_recon->getVelocityField();

	m_velocityVisualization->setVelocityField(_velocityField, m_voxelsX, m_voxelsY, m_voxelsZ, m_physicalX, m_physicalY, m_physicalZ);

	std::stringstream _ss;

	_ss << "Restored vector field ... ";
	ui.textEdit->append(QString(_ss.str().c_str()));
}

void Compounding::setLastSettings()
{
	// store the filename in the XML, so we can use it next time
	std::string _temp;
	m_xmlConfiguration.handleData(false, _temp, std::string(CALIBRATION_FILENAME_TAG));
	m_configurationFileName = QString(_temp.c_str());

	m_xmlConfiguration.handleData(false, _temp, std::string(TRACKING_FILENAME_TAG));
	m_trackingFileName = QString(_temp.c_str());

	m_xmlConfiguration.handleData(false, _temp, std::string(VIDEO_FILENAME_TAG));
	m_videoFileName = QString(_temp.c_str());

	// update the GUI
	ui.videoFileEdit->setText( m_videoFileName);
	ui.trackingFileEdit->setText( m_trackingFileName);
	ui.configurationFileEdit->setText( m_configurationFileName);

	// indicate that data is ready OR not (and update the GUI accordingly)
	if ( QFile(m_configurationFileName).exists() && QFile(m_trackingFileName).exists() && QFile(m_videoFileName).exists() )
	{

		m_guiVideoReady = true;
		m_guiTrackingReady = true;
		m_guiConfigurationReady = true;
	}
	else
	{
		m_guiVideoReady = false;
		m_guiTrackingReady = false;
		m_guiConfigurationReady = false;
	}

	guiUpdated();

	/*
	if ( System::IO::File::Exists(m_configurationFileName) )
	{
	}*/

}

void Compounding::addNoise()
{
	double _mean = atof(ui.MeanEdit->text(). toLatin1().data());
	double _std = sqrt(atof(ui.VarianceEdit->text(). toLatin1().data()));
	double _rate = atof(ui.RateEdit->text(). toLatin1().data());;
	for (int i=0;i<m_voxelsX*m_voxelsY*m_voxelsZ;i++)
	{
		if ( m_recon->getVelocityFieldValidity(i) &&  m_numericTools->generateUniformRand() <= _rate)
		{
			Eigen::Vector3f _vec = m_recon->getVelocityVectorByIndex(i);

			float _x = m_numericTools->generateStandardDistributionRand(_mean, _std);
			float _y = m_numericTools->generateStandardDistributionRand(_mean, _std);
			float _z = m_numericTools->generateStandardDistributionRand(_mean, _std);

			_vec[0]+=_x;
			_vec[1]+=_y;
			_vec[2]+=_z;

			m_recon->setVelocityVectorByIndex(i, _vec);
		}
	}

	{
		std::stringstream _ss;
		_ss << "Noise added ... ";

		ui.textEdit->append(_ss.str().c_str());

		_ss.str("");
	}

	Eigen::Vector3f* _velocityField = m_recon->getVelocityField();

	m_velocityVisualization->setVelocityField(_velocityField, m_voxelsX, m_voxelsY, m_voxelsZ, m_physicalX, m_physicalY, m_physicalZ);
}

void Compounding::initializeData()
{
	bool _doppler = false;

	// identify the sender
	QPushButton *_tempButton = (QPushButton *)sender();
	if ( _tempButton != NULL )
	if ( strcmp(_tempButton->objectName(). toLatin1().data(),"InitializeDopplerButton") == 0 )
	{
		// we want to reconstruct Doppler data, this need to be handled as a special case
		_doppler = true;
	}

	loadSettingFromXMLFile(m_configurationFileName);

	// LEGACY SOLUTION
	m_calibrationMatrixVector.push_back(m_matCalibration);

	// for RF data processing we require calibration matrix from apex
	//Eigen::Matrix4f _corner_T_apex;
	//_corner_T_apex << 1.0, 0.0, 0.0, m_apexX,  0.0, 1.0, 0.0, m_apexY,  0.0, 0.0, 1.0, 0.0,  0.0, 0.0, 0.0, 1.0;
	//m_matCalibrationApex = m_matCalibration * _corner_T_apex;

	determineImageSize();
	determineCroppingArea();

	// set the parameters of the region of interest within the visualization module, such that it can be drawn properly
	m_slicesVisualization->setImageSize(m_roiMinX, m_roiMinY, m_roiMaxX, m_roiMaxY, m_scaleX, m_scaleY);


	loadTimeStampsFromFile(m_videoFileName, false);

	float _velocityVariance;
	//Eigen::Matrix3f _velocityNoise;
	// get the noise covariance matrix, when noise
	if ( ui.VelocityNoiseCheckBox->isChecked())
	{


		/*{
		std::stringstream _ss(ui.VelocityCovarianceTextEdit->toPlainText(). toLatin1().data());

		float _m[9];
		_ss >> _m[0] >> _m[1] >> _m[2] >> _m[3] >> _m[4] >> _m[5] >> _m[6] >> _m[7] >> _m[8]; 


		_velocityNoise << _m[0],_m[1],_m[2],_m[3], _m[4],_m[5],_m[6],_m[7],_m[8];
		}
		{
		std::stringstream _ss;
		_ss << _velocityNoise << "\n applying velocity noise... ";
		ui.textEdit->append(QString(_ss.str().c_str()));
		}*/

		_velocityVariance = atof(ui.VelocityNoiseSigma->text(). toLatin1().data());
	}


	// Doppler videos are handled as a special case, as we images are transformed there 
	if ( !_doppler )
		loadVideoFromFile(_velocityVariance, STANDARD);
	else // now in case of Doppler
		loadVideoFromFile(_velocityVariance, DOPPLER, false, true, ui.VelocityNoiseCheckBox->isChecked());


	{
		std::stringstream _ss;
		_ss << m_videoFileName. toLatin1().data() << "  loaded successfully ... ";
		ui.textEdit->append(QString(_ss.str().c_str()));
	}

	loadTrackingDataFromFile(m_trackingFileName, false, ui.referenceTargetCheckBox->isChecked());

	/*
	Eigen::Vector3f _diffVec = m_trackingData[1].translation();;
	float _maxDist = 0.0;
	int _maxDex = 0;
	for (int i=0;i<m_trackingData.size();i++)
	{
	std::stringstream _ss;



	// conversion to Euler angles

	Eigen::Quaternionf _tempQuat(m_trackingData[i].rotation());
	float angleX, angleY, angleZ;
	m_geometricTools.getEulerAngles(_tempQuat,angleX,angleY,angleZ,false);

	//_ss << "Tracking: " << i << " \ " << angleX << ", " << angleY << ", "  << angleZ << ", "<<  m_trackingData[i].translation().transpose() << " \ " << (m_trackingData[i].translation()-_diffVec).norm();
	_ss << "Tracking: " << i << " \ " << _tempQuat.w() << ", " << _tempQuat.x() << ", " << _tempQuat.y() << ", " << _tempQuat.z();
	if ( (m_trackingData[i].translation()-_diffVec).norm() > _maxDist  )
	{
	_maxDist = (m_trackingData[i].translation()-_diffVec).norm();
	_maxDex = i;
	}
	_diffVec = m_trackingData[i].translation();

	ui.textEdit->append(QString(_ss.str().c_str()));

	}

	std::stringstream _ss;

	_ss << "Max Distance: " << _maxDist << " \t Index: " << _maxDex;
	ui.textEdit->append(QString(_ss.str().c_str()));*/
	//char buffer[200];
	//sprintf(buffer,"%f %f %f / %f %f %f",_minVec(0),_minVec(1),_minVec(2),_maxVec(0),_maxVec(1),_maxVec(2));
	//sprintf(buffer,"%d / %d %d / %d %d / %f %f %f / %f %f %f",i,_minX,_minY,_maxX,_maxY,_bbMinX,_bbMinY,_bbMinZ,_bbMaxX,_bbMaxY,_bbMaxZ);
	//ui.textEdit->append(QString(buffer));


	m_slicesVisualization->setTrackingData(m_trackingData);



	cropMask();



	if ( ui.MasterFrameBModeRadio->isChecked())
		determineBoundingBoxFromMasterFrame();
	else if ( ui.ManualBModeRadio->isChecked() )
		determineMinimumBoundingBoxManualFromBMode();
	else
		determineMinimumBoundingBoxFromBMode();

	m_slicesVisualization->defineOrientedBoundingBox(m_extent, m_boundingBoxTransformation);


	suggestVoxelSize();
	ui.voxelSizeSlider->setEnabled(true);
	ui.voxelSizeSlider->setValue(m_voxelSize/0.025);
	determineVoxelSize();




	if (!_doppler)
	{
		ui.attachButton->setDisabled(false);
		ui.compoundButton->setEnabled(true);
	}
	else
	{
		ui.AttachDopplerButton->setDisabled(false);
		ui.PreProcessDopplerButton->setEnabled(true);
	}
}



void Compounding::initializeRFData()
{
	bool _doppler = false;

	// identify the sender
	QPushButton *_tempButton = (QPushButton *)sender();
	if ( _tempButton != NULL && strcmp(_tempButton->objectName(). toLatin1().data(),"InitializeDopplerRFButton") == 0 )
	{
		// we want to reconstruct Doppler data, this need to be handled as a special case
		_doppler = true;
	}

	UltrasoundSettings _ultrasoundSettings;

	loadSettingFromXMLFile(m_rfCalibrationFileName, _ultrasoundSettings);

	// LEGACY SOLUTION
	m_calibrationMatrixVector.push_back(_ultrasoundSettings.calibrationMatrix);

	// data to be stored in the XML
	
	m_imageFilesXML.clear();
	// legacy: this assumes all matrices to be identical
	m_calibratinMatrixXML = _ultrasoundSettings.calibrationMatrix;
	m_masterFrameXML = -1;
	m_imageFilesXML.push_back(m_rfDataFileName. toLatin1().data());

	// for RF data processing we require calibration matrix from apex
	//Eigen::Matrix4f _corner_T_apex;
	//_corner_T_apex << 1.0, 0.0, 0.0, _ultrasoundSettings.apexX,  0.0, 1.0, 0.0, _ultrasoundSettings.apexY,  0.0, 0.0, 1.0, 0.0,  0.0, 0.0, 0.0, 1.0;
	//m_matCalibrationApex = _ultrasoundSettings.calibrationMatrix * _corner_T_apex; // NEEDS TO BE PUT INTO ULTRASOUND SETTINGS




	// TODO: obtain m_rfRFDfilename from the xml file!

	// Now load all the RF specific variables from the XML file

	//if ( m_rfFileStream) delete m_rfFileStream;

	RFFileStream *_rfFileStream = new RFFileStream();
	m_rfFileStreamVector.push_back(_rfFileStream);

#ifdef USE_USHORT16
	_rfFileStream->openFile(m_rfDataFileName. toLatin1().data(), DATA_TYPE::USHORT16);
#endif

#ifdef USE_FLOAT32
	_rfFileStream->openFile(m_rfDataFileName. toLatin1().data(), DATA_TYPE::FLOAT32);
#endif

#ifdef USE_FLOAT5X32
	_rfFileStream->openFile(m_rfDataFileName. toLatin1().data(), DATA_TYPE::FLOAT5X32);
#endif

	std::string _temp = m_rfDataFileName. toLatin1().data();
	std::string _path = _temp.substr( 0, _temp.find_last_of( '/' ) +1 );

	// update the GUI 
	int _rfFrames = 0;
	for (int i=0; i < m_rfFileStreamVector.size(); i++)
	{
		_rfFrames += m_rfFileStreamVector.at(i)->getRecordingInfo().frames;
	}

	ui.RFImageLabelScrollBar->setMaximum(_rfFrames-1);
	{
		std::stringstream ss; 
		ss  << "Frame: 0" << " / " << _rfFrames;
		ui.labelFrameRF->setText( ss.str().c_str() );
	}

	// set variables
	_ultrasoundSettings.penetrationDepth = _rfFileStream->getRecordingInfo().penetrationDepth;
	_ultrasoundSettings.scanlines = _rfFileStream->getRecordingInfo().width;
	_ultrasoundSettings.samplesAlongScanline = _rfFileStream->getRecordingInfo().height;
	_ultrasoundSettings.rfFrameSizeElements = _ultrasoundSettings.scanlines*_ultrasoundSettings.samplesAlongScanline; // NEEDS TO BE PUT INTO ULTRASOUND SETTINGS
#ifdef USE_USHORT16

	_ultrasoundSettings.rfFrameSizeBytes = sizeof(unsigned short)*_ultrasoundSettings.scanlines*_ultrasoundSettings.samplesAlongScanline; // NEEDS TO BE PUT INTO ULTRASOUND SETTINGS
#endif

	#ifdef USE_FLOAT32

	_ultrasoundSettings.rfFrameSizeBytes = sizeof(float)*_ultrasoundSettings.scanlines*_ultrasoundSettings.samplesAlongScanline; // NEEDS TO BE PUT INTO ULTRASOUND SETTINGS
#endif




	_ultrasoundSettings.mmPerSample = (float)_ultrasoundSettings.penetrationDepth / (float)_ultrasoundSettings.samplesAlongScanline;
	_ultrasoundSettings.samplesPerMM = (float)_ultrasoundSettings.samplesAlongScanline / (float)_ultrasoundSettings.penetrationDepth;
	// for phased array the ray begins in the apex, therefore the offset is zero
	//_ultrasoundSettings.radius = 0;
	// for curved array the ray begins (offset) mm from the apex
	loadTimeStampsFromFile(QString((_path+_rfFileStream->getTimeStampFileName()).c_str()), true);
	loadTrackingDataFromFile(m_rfTrackingFileName, true, ui.referenceTargetRFCheckBox->isChecked());
	m_ultrasoundSettings.push_back(_ultrasoundSettings);
	// TODO: here must be a bug, wrong extent and wrong ROI etc
	if ( ui.MasterFrameRFRadio->isChecked()) 
	{
		determineMinimumBoundingBoxfromMasterFrameRF();
	}

	else if ( ui.AutomaticRFRadio->isChecked())
	{
		determineMinimumBoundingBoxfromRF();
		m_masterFrameXML = -1; 
	}

	else if ( ui.ManualRFRadio->isChecked())
	{
		determineManualBoundingBoxfromRF();
		m_masterFrameXML = -1; 
	}
		

	// set the parameters of the region of interest within the visualization module, such that it can be drawn properly
	m_slicesVisualization->setImageSize(m_roiMinX, m_roiMinY, m_roiMaxX, m_roiMaxY, m_scaleX, m_scaleY);



	m_slicesVisualization->setTrackingData(m_trackingData);


	m_slicesVisualization->defineOrientedBoundingBox(m_extent, m_boundingBoxTransformation);

	// compute the scanlines
	for (int i=0;i < m_matrixDataTracking.size(); i++)
	{
#ifdef USE_USHORT16
		StandardRF<unsigned short> *_data = new StandardRF<unsigned short>();
#endif

#ifdef USE_FLOAT32 || USE_FLOAT5X32
		StandardRF<float> *_data = new StandardRF<float>();
#endif
		_data->setTransformation(m_matrixDataTracking[i]);
		_data->setNumberOfSamples(_rfFileStream->getRecordingInfo().height);
		_data->setNumberOfScanlines(_rfFileStream->getRecordingInfo().width);
		m_rfDataVector.push_back(_data);
	}

	computeScanlines(0, 0);

	suggestVoxelSize();
	ui.voxelSizeSlider->setEnabled(true);
	ui.voxelSizeSlider->setValue(m_voxelSize/0.025);
	determineVoxelSizeRF();

	ui.compoundRFButton->setEnabled(true);
	ui.attachRFButton->setEnabled(true);
	ui.exportTrackingButton->setEnabled(true);


	// update the GUI a little bit
	ui.dataScrollBar->setDisabled(false);
	ui.dataScrollBar->setMinimum(0);
	ui.dataScrollBar->setMaximum(_rfFrames-1);

	ui.ImageLabelScrollBar->setDisabled(false);
	ui.ImageLabelScrollBar->setMinimum(0);
	ui.ImageLabelScrollBar->setMaximum(_rfFrames-1);
	ui.ImageLabelScrollBar->setValue(0);
	// now we allow the editing of the number
	ui.MasterFrameEdit2->setEnabled(true);
	ui.bookmarkButton->setEnabled(true);

	m_rfFrameCounter = _rfFrames-1;
}

void Compounding::imageViewerUpdate(int frame)
{

	if ( m_imageCounter > 0 )
	{
		showFrame(frame);

		HBITMAP _bmp;

		int _width, _height;
		IplImage* _img;
		//m_videoTools->getFrame(frame, _bmp, _width, _height);
		int _index1 = (int)floor((float)frame/NUMBER_OF_SLICES); 
		int _index2 = frame % (int)(NUMBER_OF_SLICES);
		_img = m_imageVector[_index1][_index2];

		IplImage* _colImage = cvCreateImage(cvSize(_img->width, _img->height), IPL_DEPTH_8U, 3);
		cvCvtColor(_img, _colImage, CV_GRAY2RGB);
		_bmp = m_videoTools->IplImage2HBITMAP(_colImage);
	//	cvSaveImage("color.jpg",_colImage);
		cvReleaseImage(&_colImage);

		_width = _img->width;
		_height = _img->height;
	
		QPixmap _pixmap = QtWin::fromHBITMAP(_bmp, QtWin::HBitmapNoAlpha); 


		// some given arbitary video image size
		float _targetSize = 600.0;
		float _scaleFactor = _targetSize / (float) std::max(_width,_height);
		ui.ImageLabel->setPixmap(_pixmap.scaled(_scaleFactor*_width, _scaleFactor*_height));

		//ui.ImageLabel->resize(_scaleFactor *ui.ImageLabel->pixmap()->size());


		if ( m_masterImage != 0)
			DeleteObject(m_masterImage);

		m_masterImage = _bmp;

		// update the edit boxes
		{
			std::stringstream _ss;

			_ss << frame;

			ui.MasterFrameEdit->setText(_ss.str().c_str());
			ui.MasterFrameEdit2->setText(_ss.str().c_str());
		}
	}
}

void Compounding::masterFrame(int changed)
{
	ui.MasterFrameEdit->setEnabled(ui.MasterFrameBModeRadio->isChecked());

	if ( ui.MasterFrameBModeRadio->isChecked())
		determineBoundingBoxFromMasterFrame();
	/*else
	determineMinimumBoundingBox();
	*/
	m_slicesVisualization->defineOrientedBoundingBox(m_extent, m_boundingBoxTransformation);


	//suggestVoxelSize();
	ui.voxelSizeSlider->setEnabled(true);
	ui.voxelSizeSlider->setValue(m_voxelSize/0.025);
	determineVoxelSize();
}

void Compounding::masterFrameRF(int changed)
{
	ui.RFMasterFrameEdit2->setEnabled(ui.MasterFrameRFRadio->isChecked());

	if ( ui.MasterFrameRFRadio->isChecked()) 
	{
		determineMinimumBoundingBoxfromMasterFrameRF();
	}

	else if ( ui.AutomaticRFRadio->isChecked())
	{
		determineMinimumBoundingBoxfromRF();
		m_masterFrameXML = -1; 
	}

	else if ( ui.ManualRFRadio->isChecked())
	{
		determineManualBoundingBoxfromRF();
		m_masterFrameXML = -1; 
	}
	
	m_slicesVisualization->defineOrientedBoundingBox(m_extent, m_boundingBoxTransformation);


	//suggestVoxelSize();
	ui.voxelSizeSlider->setEnabled(true);
	ui.voxelSizeSlider->setValue(m_voxelSize/0.025);
	determineVoxelSize();
}


void Compounding::masterSliceChanged(QString str)
{

	if ( ui.MasterFrameBModeRadio->isChecked())
		determineBoundingBoxFromMasterFrame();
	else // obviously the text of the edit field was changed, so update the slicer viewer
	{
		int _masterFrameIndex = atoi(ui.MasterFrameEdit2->text(). toLatin1().data());
		ui.ImageLabelScrollBar->setValue(_masterFrameIndex);
	}

	/*else
	determineMinimumBoundingBox();
	*/
	m_slicesVisualization->defineOrientedBoundingBox(m_extent, m_boundingBoxTransformation);


	//suggestVoxelSize();
	ui.voxelSizeSlider->setEnabled(true);
	ui.voxelSizeSlider->setValue(m_voxelSize/0.025);
	determineVoxelSize();
}


void Compounding::finishedOutputFilenameEditing()

{
	m_mhdFileName = ui.outputFileName->text();
}

void Compounding::addBookmark()
{
	ui.ImageLabelScrollBar->value();
	std::stringstream _ss;
	_ss << ui.ImageLabelScrollBar->value();
	ui.bookmarkListView->addItem(_ss.str().c_str());
	ui.bookmarkListView->update();
}

void Compounding::doubleClickBookmark(QListWidgetItem *item)
{
	ui.ImageLabelScrollBar->setValue(atoi(item->text(). toLatin1().data()));
}


Eigen::VectorXf Compounding::jacobianNormalUncertainty(Eigen::Vector3f &target, Eigen::Vector3f &apex, Eigen::Vector3f &velocity)
{

	float _res[7];
	{

		// diff J_ax

		//["cg", "cg1", "cg3", "cg5", "cg7", "cg9"] = ["tx~", "ty~", "tz~", "vx~", "vy~", "vz~"]
		float t1 = (int) (apex.x() - target.x());
		float t2 = abs(t1);
		float t4 = fabs((double) t1) / (double) t1;
		float t5 = t4 * (double) t2 * velocity.x();
		float t6 = (apex.x()>0 ? 1 : (apex.x()<0 ? -1 : 0));
		float t7 = (double) t6 * apex.x();
		float t10 = apex.x(); //conjugate(apex.x());
		float t13 = t10 * velocity.x();
		float t14 = t2 * t2;
		float t18 = abs(apex.y() - target.y());
		float t19 = t18 * t18;
		float t23 = abs(apex.z() - target.z());
		float t24 = t23 * t23;
		float t27 = fabs(apex.x()) / apex.x();
		float t28 = t27 * velocity.x();
		float t39 = t4 * (double) t2 * velocity.y();
		float t42 = apex.y(); //conjugate(apex.y());
		float t46 = t4 * (double) t2 * velocity.z();
		float t49 = apex.z();//;conjugate(apex.z());
		float t52 = -target.x() * t7 * t5 + t10 * t7 * t5 + (double) t14 * (double) t6 * t13 + (double) t19 * (double) t6 * t13 + (double) t24 * (double) t6 * t13 - 0.2e1 * (double) t14 * apex.x() * t28 - 0.2e1 * (double) t19 * apex.x() * t28 - 0.2e1 * (double) t24 * apex.x() * t28 - (double) target.y() * t7 * t39 + t42 * t7 * t39 - (double) target.z() * t7 * t46 + t49 * t7 * t46;
		float t53 = t14 + t19 + t24;
		float t54 = sqrt((double) t53);
		_res[0] = -0.1e1 / (double) t6 / apex.x() / t54 / (double) t53 * t52;
	}


	{
		// diff J_ay

		// ["cg", "cg1", "cg3", "cg5", "cg7", "cg9"] = ["tx~", "ty~", "tz~", "vx~", "vy~", "vz~"]
		float t1 = (int) (apex.y() - target.y());
		float t2 = abs(t1);
		float t4 = fabs((double) t1) / (double) t1;
		float t5 = t4 * (double) t2 * velocity.x();
		float t6 = (apex.y()>0 ? 1 : (apex.y()<0 ? -1 : 0));
		float t7 = (double) t6 * apex.y();
		float t10 = apex.x(); //conjugate(apex.x());
		float t14 = t4 * (double) t2 * velocity.y();
		float t17 = apex.y(); // conjugate(apex.y());
		float t20 = t17 * velocity.y();
		float t22 = abs(apex.x() - target.x());
		float t23 = t22 * t22;
		float t26 = t2 * t2;
		float t30 = abs(apex.z() - target.z());
		float t31 = t30 * t30;
		float t34 = fabs(apex.y()) / apex.y();
		float t35 = t34 * velocity.y();
		float t46 = t4 * (double) t2 * velocity.z();
		float t49 = apex.z(); // conjugate(apex.z());
		float t52 = -(double) target.x() * t7 * t5 + t10 * t7 * t5 - target.y() * t7 * t14 + t17 * t7 * t14 + (double) t23 * (double) t6 * t20 + (double) t26 * (double) t6 * t20 + (double) t31 * (double) t6 * t20 - 0.2e1 * (double) t23 * apex.y() * t35 - 0.2e1 * (double) t26 * apex.y() * t35 - 0.2e1 * (double) t31 * apex.y() * t35 - (double) target.z() * t7 * t46 + t49 * t7 * t46;
		float t53 = t23 + t26 + t31;
		float t54 = sqrt((double) t53);
		_res[1] = -0.1e1 / (double) t6 / apex.y() / t54 / (double) t53 * t52;
	}

	{
		// diff J_az

		//["cg", "cg1", "cg3", "cg5", "cg7", "cg9"] = ["tx~", "ty~", "tz~", "vx~", "vy~", "vz~"]
		float t1 = (int) (apex.z() - target.z());
		float t2 = abs(t1);
		float t4 = fabs((double) t1) / (double) t1;
		float t5 = t4 * (double) t2 * velocity.x();
		float t6 = (apex.z()>0 ? 1 : (apex.z()<0 ? -1 : 0));
		float t7 = (double) t6 * apex.z();
		float t10 = apex.x(); // conjugate(apex.x());
		float t14 = t4 * (double) t2 * velocity.y();
		float t17 = apex.y(); // conjugate(apex.y());
		float t21 = t4 * (double) t2 * velocity.z();
		float t24 = apex.z(); // conjugate(apex.z());
		float t27 = t24 * velocity.z();
		float t29 = abs(apex.x() - target.x());
		float t30 = t29 * t29;
		float t34 = abs(apex.y() - target.y());
		float t35 = t34 * t34;
		float t38 = t2 * t2;
		float t41 = fabs(apex.z()) / apex.z();
		float t42 = t41 * velocity.z();
		float t52 = -(double) target.x() * t7 * t5 + t10 * t7 * t5 - (double) target.y() * t7 * t14 + t17 * t7 * t14 - target.z() * t7 * t21 + t24 * t7 * t21 + (double) t30 * (double) t6 * t27 + (double) t35 * (double) t6 * t27 + (double) t38 * (double) t6 * t27 - 0.2e1 * (double) t30 * apex.z() * t42 - 0.2e1 * (double) t35 * apex.z() * t42 - 0.2e1 * (double) t38 * apex.z() * t42;
		float t53 = t30 + t35 + t38;
		float t54 = sqrt((double) t53);
		_res[2] = -0.1e1 / (double) t6 / apex.z() / t54 / (double) t53 * t52;
	}


	{
		// diff J_tx

		//["cg", "cg1", "cg3", "cg5", "cg7", "cg9"] = ["tx~", "ty~", "tz~", "vx~", "vy~", "vz~"]
		float t1 = (int) (apex.x() - target.x());
		float t2 = abs(t1);
		float t3 = t2 * velocity.x();
		float t4 = fabs((double) t1) / (double) t1;
		float t7 = apex.x(); // conjugate(apex.x());
		float t10 = t2 * t2;
		float t13 = abs(apex.y() - target.y());
		float t14 = t13 * t13;
		float t17 = abs(apex.z() - target.z());
		float t18 = t17 * t17;
		float t20 = t2 * velocity.y();
		float t23 = apex.y(); // conjugate(apex.y());
		float t26 = t2 * velocity.z();
		float t29 = apex.z(); //conjugate(apex.z());
		float t33 = t10 + t14 + t18;
		float t34 = sqrt((double) t33);
		_res[3] = 0.1e1 / t34 / (double) t33 * (-target.x() * t4 * (double) t3 + t7 * t4 * (double) t3 - (double) (t10 * velocity.x()) - (double) (t14 * velocity.x()) - (double) (t18 * velocity.x()) - (double) target.y() * t4 * (double) t20 + t23 * t4 * (double) t20 - (double) target.z() * t4 * (double) t26 + t29 * t4 * (double) t26);

	}

	{

		// diff J_ty

		//["cg", "cg1", "cg3", "cg5", "cg7", "cg9"] = ["tx~", "ty~", "tz~", "vx~", "vy~", "vz~"]
		float t1 = (int) (apex.y() - target.y());
		float t2 = abs(t1);
		float t3 = t2 * velocity.x();
		float t4 = fabs((double) t1) / (double) t1;
		float t7 = apex.x(); // conjugate(apex.x());
		float t10 = t2 * velocity.y();
		float t13 = apex.y(); // conjugate(apex.y());
		float t17 = abs(apex.x() - target.x());
		float t18 = t17 * t17;
		float t20 = t2 * t2;
		float t23 = abs(apex.z() - target.z());
		float t24 = t23 * t23;
		float t26 = t2 * velocity.z();
		float t29 = apex.z(); // conjugate(apex.z());
		float t33 = t18 + t20 + t24;
		float t34 = sqrt((double) t33);
		_res[4] = 0.1e1 / t34 / (double) t33 * (-(double) target.x() * t4 * (double) t3 + t7 * t4 * (double) t3 - target.y() * t4 * (double) t10 + t13 * t4 * (double) t10 - (double) (t18 * velocity.y()) - (double) (t20 * velocity.y()) - (double) (t24 * velocity.y()) - (double) target.z() * t4 * (double) t26 + t29 * t4 * (double) t26);

	}


	{

		// diff J_tz

		//["cg", "cg1", "cg3", "cg5", "cg7", "cg9"] = ["tx~", "ty~", "tz~", "vx~", "vy~", "vz~"]
		float t1 = (int) (apex.z() - target.z());
		float t2 = abs(t1);
		float t3 = t2 * velocity.x();
		float t4 = fabs((double) t1) / (double) t1;
		float t7 = apex.x(); // conjugate(apex.x());
		float t10 = t2 * velocity.y();
		float t13 = apex.y(); // conjugate(apex.y());
		float t16 = t2 * velocity.z();
		float t19 = apex.z(); // conjugate(apex.z());
		float t23 = abs(apex.x() - target.x());
		float t24 = t23 * t23;
		float t27 = abs(apex.y() - target.y());
		float t28 = t27 * t27;
		float t30 = t2 * t2;
		float t33 = t24 + t28 + t30;
		float t34 = sqrt((double) t33);
		_res[5] = 0.1e1 / t34 / (double) t33 * (-(double) target.x() * t4 * (double) t3 + t7 * t4 * (double) t3 - (double) target.y() * t4 * (double) t10 + t13 * t4 * (double) t10 - target.z() * t4 * (double) t16 + t19 * t4 * (double) t16 - (double) (t24 * velocity.z()) - (double) (t28 * velocity.z()) - (double) (t30 * velocity.z()));

	}
	Eigen::VectorXf _result = VectorXf::Zero(6);

	_result << _res[0], _res[1], _res[2], _res[3], _res[4], _res[5];

	return _result;
}

void CompoundRFSubvolumeWritingThread::setDataTarget(StandardMeasurementCollectionManager<IntensityAndDistance> *ptr)
{
	m_measurementsCollectionManager = ptr;
}


CompoundRFSubvolumeWritingThread::CompoundRFSubvolumeWritingThread(std::string filename, int voxelsX, int voxelsY, int voxelsZ, int maxElements, __int64 size)
{
	/*std::cout << filename.c_str() << std::endl;
	int _fh;
	if ( _sopen_s( &_fh, filename.c_str(), _O_RDWR | _O_CREAT, _SH_DENYNO, _S_IREAD | _S_IWRITE ) == 0)
	std::cout << "file created!" << std::endl;
	else
	std::cout << "error while creating file!" << std::endl;

	if ( _chsize(_fh, size) == 0 )
	std::cout << "file size changed!" << std::endl;
	else
	std::cout << "error while changing file size!" << std::endl;

	_close(_fh);
	m_rfDataFile.open(filename.c_str(), std::ios::binary | std::ios::ate | std::ios::in | std::ios::out);
	if ( m_rfDataFile.fail() )
	std::cout << "failure loading file!" << std::endl;

	m_rfDataFile.seekg(0, std::ios::beg);
	m_rfDataFile.seekp(0, std::ios::beg);
	*/
	//m_rfDataFile.open(filename.c_str(), std::ios::binary | std::ios::out);
	
	std::stringstream _ss;
	
	remove(filename.c_str());
	
	_ss << "contig -n " << filename << " " << size;
	std::cout << _ss.str() << std::endl;
	system(_ss.str().c_str());
	m_rfDataFile = CreateFileA( filename.c_str(), GENERIC_WRITE ,  0, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL,  NULL);
	if ( m_rfDataFile ==  INVALID_HANDLE_VALUE )
		std::cout << "error opening file!" << std::endl;
	m_maxElements = maxElements;
	m_voxelsX = voxelsX;
	m_voxelsY = voxelsY;
	m_voxelsZ = voxelsZ;
	m_filename = filename;
	m_fileSize = size;
}

CompoundRFSubvolumeWritingThread::~CompoundRFSubvolumeWritingThread()
{
}

void CompoundRFSubvolumeWritingThread::run()
{
	std::cout << "Started writing " << m_filename << std::endl;
#if DEBUG_MHD_OUTPUT

	vtkImageData *_compoundingOutput = vtkImageData::New();
	// set the parameters of the image

	//_compoundingOutput->SetScalarTypeToUnsignedChar(); 
	
	_compoundingOutput->SetDimensions(m_voxelsX, m_voxelsY, m_voxelsZ); 
	_compoundingOutput->SetOrigin(0.0, 0.0, 0.0);  
	_compoundingOutput->SetSpacing(0.4, 0.4, 0.4); 
	_compoundingOutput->AllocateScalars(VTK_DOUBLE, 1);

	/*vtkImageData *_compoundingOutput = vtkImageData::New();
	// set the parameters of the image

	//_compoundingOutput->SetScalarTypeToUnsignedChar(); 
	_compoundingOutput->SetScalarTypeToFloat();
	_compoundingOutput->SetDimensions(m_voxelsX, m_voxelsY, m_voxelsZ); 
	_compoundingOutput->SetOrigin(0.0, 0.0, 0.0); 
	_compoundingOutput->SetSpacing(0.4, 0.4, 0.4); 
	_compoundingOutput->Update(); 
*/
	// clear the memory
	float *_tmpIntensityPtr = static_cast<float*>(_compoundingOutput->GetScalarPointer(0,0,0));
	
#endif

#if DEBUG_TXT_OUTPUT
	std::ofstream _debugFile;
	_debugFile.open(m_filename+".txt");
#endif

	float _1percent = (float)m_fileSize/100.0;
	long _localCount = 0;
	int _percent = 0;
	DataTools _dataTools;

	clock_t _startTimer = clock();
	
	char *_dataChunk = new char[1024*1024*30]; 
	int _offset = 0;
	int _counter[3]={0,0,0};

	std::vector<StandardMeasurementCollection<IntensityAndDistance>*>  *_ptr = m_measurementsCollectionManager->getData(0);

	int _intensityCounter = 0;
	for (int k=0; k< m_maxElements;k++)
	{
		_counter[0]++;
		//std::cout << "Thread No. " << m_filename << "/ Comp. . "<< k << " / "<< m_maxElements<< std::endl;
		float _maxVal = 0.0;		

		int _val = _ptr->size();

		ushort16 _elements;
		_elements.value = (unsigned short)_ptr->size();

		//m_rfDataFile.write(&_elements.byte.c[0], sizeof(ushort16));
		DWORD _writtenBytes;
		//WriteFile(m_rfDataFile, &_elements.byte.c[0], sizeof(ushort16), &_writtenBytes, NULL);
		memcpy(&_dataChunk[_offset],&_elements.byte.c[0], sizeof(ushort16));
		_offset+=sizeof(ushort16);

#if DEBUG_TXT_OUTPUT
		_debugFile << k << ". {" << _elements.value << "}\t";
#endif

			int _sumElements = 0;
		// debug dump
		//	std::ofstream _dumpFile;
//bool _dump = false;
			// debug
	//		if ( strcmp(m_filename.c_str(),"subvolume2.tmp") == 0 && k == 1173114)
		//		{
			//		_dumpFile.open("dump.txt");
			//		_dump = true;
			//	}
		for(int j=0;j<_val;j++)
		{
		
			_counter[1]++;
			StandardMeasurementCollection<IntensityAndDistance>* _smc = (*_ptr)[j];//_ptr->back();

		

			int _elements = _smc->count; //.size();

			int _size = sizeof(ushort16)+sizeof(uint32);
			
#ifdef USE_FLOAT32

			_size += sizeof(float32)*_elements*2+sizeof(ushort16)*_elements;
#endif

#ifdef USE_USHORT16
			_size += sizeof(float32)*_elements+sizeof(ushort16)*_elements*2;
#endif

#ifdef USE_FLOAT5X32
			_size += sizeof(float32)*_elements+sizeof(ushort16)*_elements+5*sizeof(float32)*_elements;
#endif

			//char *_tmp = new char[_size];
			//int _offset = 0;
			ushort16 _ushort16;
			_ushort16.value = (unsigned short)_smc->count;//data.size();
			//m_rfDataFile.write(&_ushort16.byte.c[0], sizeof(ushort16));
			//if(!WriteFile(m_rfDataFile, &_ushort16.byte.c[0], sizeof(ushort16), &_writtenBytes, NULL))
			//	std::cout << "Error while  writing data! " << std::endl;
			memcpy(&_dataChunk[_offset],&_ushort16.byte.c[0], sizeof(ushort16));
			_offset+=sizeof(ushort16);

			// handle the direction index
			{
				uint32 _uint32;
				_uint32.value = (unsigned int)_smc->uniqueID;
				//m_rfDataFile.write(&_uint32.byte.c[0], sizeof(uint32));
				//if(!WriteFile(m_rfDataFile, &_uint32.byte.c[0], sizeof(uint32), &_writtenBytes, NULL))
				//	std::cout << "Error while  writing data! " << std::endl;
				memcpy(&_dataChunk[_offset],&_uint32.byte.c[0], sizeof(uint32));
				_offset+=sizeof(uint32);
	

#if DEBUG_TXT_OUTPUT
			_debugFile << _elements << "\t";
			if ( _dump )
			{
				_dumpFile << "DirectionIndex: " << _uint32.value << std::endl;
				_dumpFile << "Elements: " << _elements << std::endl;
			}
#endif
					}
	
		
			for(int i=0;i<_elements;i++)
			{
				_sumElements++;
				_counter[2]++;
				IntensityAndDistance _iad = _smc->data[i];// _smc->data.back();
				//_smc->data.pop_back();
				//ushort16 _ushort16Distance;
				//_ushort16Distance.value = (unsigned short)_iad.distance;
				
				
				#ifdef SAVE_DISTANCE_SCANLINE
#ifdef  QUANTIZE_DISTANCE
				char _distance;
				_distance = (char)_iad.distance;
				ushort16 _scanlineIndex;
				_scanlineIndex.value = /*(float)*/_iad.scanlineIndex;
				memcpy(&_dataChunk[_offset],&_distance, sizeof(char));
				_offset+=sizeof(char);


				#if DEBUG_TXT_OUTPUT
					if ( _dump )
						_dumpFile << "{ " <<_scanlineIndex.value  << "\t";
				#endif


#else
				float32 _distance;
				_distance.value = (float)_iad.distance;
				ushort16 _scanlineIndex;
				_scanlineIndex.value = /*(float)*/_iad.scanlineIndex;
				memcpy(&_dataChunk[_offset],&_distance.byte.c[0], sizeof(float32));
				_offset+=sizeof(float32);
#endif
				#endif
				
				
				

				//m_rfDataFile.write(&_float32.byte.c[0], sizeof(float32));
				//if(!WriteFile(m_rfDataFile, &_float32.byte.c[0], sizeof(float32), &_writtenBytes, NULL))
				//	std::cout << "Error while  writing data! " << std::endl;
#ifdef USE_USHORT16
				ushort16 _intensity16;
				_intensity16.value = _iad.intensity;
				memcpy(&_dataChunk[_offset],&_intensity16.byte.c[0], sizeof(ushort16));
				_offset+=sizeof(ushort16);


				#if DEBUG_TXT_OUTPUT
					if ( _dump )
						_dumpFile << _intensity16.value  << "\t";
				#endif
#endif

#ifdef USE_FLOAT32

				float32 _intensity32;
				_intensity32.value = _iad.intensity;
				memcpy(&_dataChunk[_offset],&_intensity32.byte.c[0], sizeof(float32));
				_offset+=sizeof(float32);
#endif

#ifdef USE_FLOAT5X32
				for(int k=0;k<5;k++) {
				float32 _intensity32;
				_intensity32.value = _iad.mrf[k];
				memcpy(&_dataChunk[_offset],&_intensity32.byte.c[0], sizeof(float32));
				_offset+=sizeof(float32);
				}
#endif
				//m_rfDataFile.write(&_scanlineIndex.byte.c[0], sizeof(float32));
				//if(!WriteFile(m_rfDataFile, &_scanlineIndex.byte.c[0], sizeof(float32), &_writtenBytes, NULL))
				//	std::cout << "Error while  writing data! " << std::endl;
				#ifdef SAVE_DISTANCE_SCANLINE
				memcpy(&_dataChunk[_offset],&_scanlineIndex.byte.c[0], sizeof(ushort16));
				_offset+=sizeof(ushort16);


				#if DEBUG_TXT_OUTPUT
					if ( _dump )
						_dumpFile << _scanlineIndex.value << "}" << "\t";
				#endif
#endif
				if ( _maxVal < _iad.intensity )
				{
				_maxVal = _iad.intensity;
				_intensityCounter++;
				}/**/	

				
			}

			#if DEBUG_TXT_OUTPUT
					if ( _dump )
						_dumpFile << std::endl;
				#endif

		
			//WriteFile(m_rfDataFile, &_tmp, _size, &_writtenBytes, NULL);
			//std::cout << "Wrote: " << _size << " / " << _writtenBytes << std::endl;
			//delete[] _tmp;

			_localCount += _size;
			
/*
					#if DEBUG_TXT_OUTPUT
			_debugFile << _sumElements << "\t";
#endif*/
			/*if ( _smc)
				delete _smc;
			*/
			//_junkVec.push_back(_smc);
			//_ptr->pop_back();
			
			
			//std::cout << "Finished Scanline Thread No. " << m_filename << std::endl;

			

		}

		#if DEBUG_TXT_OUTPUT
			
			if ( _dump )
			{
				_dumpFile.close();
			}
#endif


	

		if  ( _offset >= 1024*1024*5 || k == m_maxElements-1)
		{
			WriteFile(m_rfDataFile, &_dataChunk[0], _offset, &_writtenBytes, NULL);
		//std::cout << _writtenBytes << " bytes written!" << std::endl;
			_offset = 0;
		}

		if ( (float)_localCount >= _1percent )
			{
				//qApp->processEvents();
				std::cout << "\rFile writing: " << ++_percent << " %  - time: " << _dataTools.computeElapsedTime(_startTimer, clock())<<"  " <<_counter[0]<< " / " << _counter[1] << " / " << _counter[2] << "                ";
				_localCount = 0;
				_counter[0]=0;
				_counter[1]=0;
				_counter[2]=0;
				_startTimer = clock();
			}
			/*if ( _percent%10 == 0 )
				if ( !FlushFileBuffers(m_rfDataFile) )
					std::cout << "Error while  flushing data! " << std::endl;

					*/
#if DEBUG_TXT_OUTPUT
		//_debugFile << "["<< _maxVal << "]" << std::endl;
		_debugFile << "["<< _sumElements << "]" << std::endl;
#endif
		_ptr++;
#if DEBUG_MHD_OUTPUT
		(*_tmpIntensityPtr) = _maxVal;
		_tmpIntensityPtr++;
#endif

	}
	delete []_dataChunk;
	//m_rfDataFile.flush();
	//m_rfDataFile.close();
	CloseHandle(m_rfDataFile);
	std::cout << std::endl <<"Finished writing " << m_filename << std::endl;


	//setPriority(QThread::Priority::IdlePriority);


	_startTimer = clock();
	std::cout << "Deleting temporary data ...!" << std::endl;
	delete m_measurementsCollectionManager;
	std::cout << "Delete finished after: " << _dataTools.computeElapsedTime(_startTimer, clock()) << std::endl;
#if DEBUG_TXT_OUTPUT
	_debugFile.close();
#endif

#if DEBUG_MHD_OUTPUT

	vtkSmartPointer<vtkMetaImageWriter> _writer = vtkSmartPointer<vtkMetaImageWriter>::New();
			
			vtkSmartPointer<vtkImageCast> _castFilter =
				vtkSmartPointer<vtkImageCast>::New();
			_castFilter->SetOutputScalarTypeToUnsignedChar();
			_castFilter->SetInputData(_compoundingOutput);

			_writer->SetInputConnection(_castFilter->GetOutputPort());
			_writer->SetFileName((m_filename+".mhd").c_str()); 
			_writer->SetCompression(false);
			_writer->Write();
	/*vtkMetaImageWriter *_writer = vtkMetaImageWriter::New(); 
	_writer->SetInput(_compoundingOutput); 
	_writer->SetFileName((m_filename+".mhd").c_str()); 
	_writer->SetCompression( false );
	_writer->Write();

	_writer->Delete();*/
#endif
}


CompoundRFVolumeThread::CompoundRFVolumeThread() 
{
}

void CompoundRFVolumeThread::run()
{
	NumericTools _nr;
	// we have to create our own heap now in order to avoid the global lock of the
	// process' own heap
	HANDLE _heap = HeapCreate(HEAP_NO_SERIALIZE, 0x100000, 0);
	if ( _heap == NULL )
		std::cout << "Failed to create heap! " << std::endl;
	
	/*	HANDLE _rfDataFile = CreateFileA( m_filename.c_str(), GENERIC_WRITE ,  0, NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL,  NULL);
	if ( _rfDataFile ==  INVALID_HANDLE_VALUE )
		std::cout << "error opening file!" << std::endl;*/

	std::cout << "Launched Thread No. : " << m_index /*<< " / " << m_filename */<< std::endl;
	// determine the maximum dimension of the voxel spacing
	float _maxDim = m_voxelSize;
	float _voxelDiameter = sqrt(m_physicalX*m_physicalX+m_physicalY*m_physicalY+m_physicalZ*m_physicalZ)/2.0f;

	// choose the distance threshold, that is the maximum distance a US slice can have
	// in order to be considered for contributing intensity to a voxel

	float _distanceScalar = 0.0f;
	switch(m_maxDistanceScalar)
	{
	case 0: { _distanceScalar = DISTANCE_SCALAR; break;}
	case 1: { _distanceScalar = 1.0; break;}
	case 2: { _distanceScalar = 1.1; break;}
	case 3: { _distanceScalar = 1.2; break;}
	case 4: { _distanceScalar = 1.3; break;}
	case 5: { _distanceScalar = 1.4; break;}
	case 6: { _distanceScalar = 1.5; break;}
	};

	float _threshold = 0.0f;
	switch(m_maxDistanceUnit)
	{
	case 0: { _threshold = _maxDim; break;}
	case 1: { _threshold = _voxelDiameter; break;}
	};

	// now we compute the actual threshold value from the UI selection
	_threshold = _threshold * _distanceScalar;
	
	// now perform a test reconstruction in order to determine the number of structure elements needed
	// that is necessary because deallocating a large number of pointers is very slow
	// we want to allocate one contiguous pointer block
	int _smcCounter = 0; // this is what we determine
	{ // BEGIN TEST

		// initialize the rotation queue
#ifdef USE_USHORT16
	RotationQueue<Scanline<unsigned short>> _rotQueue(_maxDim,m_extent.dimX+m_extent.dimY+m_extent.dimZ, (int)((float)m_scanlineVector->size()/10.0));
#endif

#ifdef USE_FLOAT32 || USE_FLOAT5X32
	RotationQueue<Scanline<float>> _rotQueue(_maxDim,m_extent.dimX+m_extent.dimY+m_extent.dimZ, (int)((float)m_scanlineVector->size()/10.0));
#endif
	// now fill the rotation queue with the scanlines
	// now we set up a data structure containing information about the position of the ultrasound planes in space
	for (int i=0; i<m_scanlineVector->size();i++)
	{
		// set the ultrasound scanline at initial distance 0.0
		_rotQueue.push((*m_scanlineVector)[i], 0.0);
	}


	// now perform the actual (pseudo) compounding


	TraversalScheme _traversalScheme(m_voxelsX, m_voxelsY, m_voxelsZ, m_offsetX, m_offsetY, m_offsetZ);
	_traversalScheme.init();


	int _x,_y,_z;


	while(_traversalScheme.nextPosition(_x,_y,_z))
	{
		int _intersectionCounter = 0;

		// compute the world coordinates from the voxel coordinates
		const Eigen::Vector3f _local(static_cast<float>(_x)*m_physicalX+m_physicalX/2.0,static_cast<float>(_y)*m_physicalY+m_physicalY/2.0,static_cast<float>(_z)*m_physicalZ+m_physicalZ/2.0);


		const Eigen::Vector3f _world = ( m_boundingBoxTransformation *_local );


		// now we can look for the scanlines in vicinity

		// pop the minimum distance bin vector from the rotation queue
#ifdef USE_USHORT16
		std::vector<Scanline<unsigned short>*>* _usScanlineQueue = _rotQueue.pop();
#endif
#ifdef USE_FLOAT32 || USE_FLOAT5X32
		std::vector<Scanline<float>*>* _usScanlineQueue = _rotQueue.pop();
#endif
		int _oldNumberOfScanlinesInQueue = _usScanlineQueue->size();

		// for each slice within the vector returned form the queue do a backward warping
		while( _usScanlineQueue->size() > 0 )
		{

#ifdef USE_USHORT16
			Scanline<unsigned short> *_scanline = _usScanlineQueue->back();
#endif

#ifdef USE_FLOAT32 || USE_FLOAT5X32
			Scanline<float> *_scanline = _usScanlineQueue->back();
#endif

			_usScanlineQueue->pop_back();


			// determine the orthogonal distance from the voxel to the ultrasound plane
			float _lineParam = 0.0f;
			const float _distance = _scanline->distancePointToScanline(_world, _lineParam);

			// HACK to avoid overflow
			// the CONDITION: _intersectionCounter >= std::numeric_limits<unsigned short>::max()-1
			if ( _distance <= _threshold && _lineParam >= 0 && _lineParam <= 1.0 && _intersectionCounter <= std::numeric_limits<unsigned short>::max()-1)
			{
				//_smc = new StandardMeasurementCollection<IntensityAndDistance>();

				// and push it into the temporary queue, as the point is not to be considered again for the same voxel
				_rotQueue.push(_scanline, _distance);

				_smcCounter++;
				_intersectionCounter++;
		
				

			}
			else // too far away
			{

				// and push it into the temporary queue, as the point is not to be considered again for the same voxel
				_rotQueue.push(_scanline, _distance);
			}



			
		}
	


		// instead of delete put empty vectors back to rotation queue, so we save new operations
		_usScanlineQueue->reserve(_oldNumberOfScanlinesInQueue);
		_rotQueue.recycle(_usScanlineQueue);

	}

	} // END TEST

	// initialize the rotation queue
#ifdef USE_USHORT16
	RotationQueue<Scanline<unsigned short>> _rotQueue(_maxDim,m_extent.dimX+m_extent.dimY+m_extent.dimZ, (int)((float)m_scanlineVector->size()/10.0));
#endif

#ifdef USE_FLOAT32 || USE_FLOAT5X32
	RotationQueue<Scanline<float>> _rotQueue(_maxDim,m_extent.dimX+m_extent.dimY+m_extent.dimZ, (int)((float)m_scanlineVector->size()/10.0));
#endif
	// now fill the rotation queue with the scanlines
	// now we set up a data structure containing information about the position of the ultrasound planes in space
	for (int i=0; i<m_scanlineVector->size();i++)
	{
		// set the ultrasound scanline at initial distance 0.0
		_rotQueue.push((*m_scanlineVector)[i], 0.0);
	}


	// now perform the actual compounding


	TraversalScheme _traversalScheme(m_voxelsX, m_voxelsY, m_voxelsZ, m_offsetX, m_offsetY, m_offsetZ);
	_traversalScheme.init();


	int _x,_y,_z;
	int _traversalCount = 0;
	int _percentCounter = 0;

	/*char *_dataChunk = new char[1024*1024*30]; 
	int _offset = 0;
	int _offsetIndex = 0;
	DWORD _writtenBytes = 0;*/
	// as we know exactly how many elements we need, we have a shrink to fit contiguous data block
	// which allows for fast deallocation!

	
	StandardMeasurementCollection<IntensityAndDistance> *_smcBlock = (StandardMeasurementCollection<IntensityAndDistance> *)HeapAlloc(_heap, HEAP_NO_SERIALIZE, sizeof(StandardMeasurementCollection<IntensityAndDistance>)*_smcCounter);
	if ( _smcBlock == NULL)
		std::cout << "Memory allocation failed on the heap! " << _smcCounter << std::endl;
	else
		std::cout << "SMC: " << _smcCounter << std::endl;

	//StandardMeasurementCollection<IntensityAndDistance> *_smcBlock = new StandardMeasurementCollection<IntensityAndDistance>[_smcCounter];
	m_measurementsCollectionManager->addBlockPointer(_smcBlock, _heap); // this stores the actual data
	int _maxC = _smcCounter;
	_smcCounter = 0; // now we use the counter to keep track of the data
	while(_traversalScheme.nextPosition(_x,_y,_z))
	{
		// count the number of scanlines intersecting
		int _intersectionCounter = 0;

		//std::cout << "Thread No. " << m_index << "/ Comp. "<< _traversalCount << std::endl;
		//int _intersectingSlices = 0;
		//std::cout << "Begin: " << _traversalCount << std::endl;
		_percentCounter++;
		_traversalCount++;

		// this is the index where the number of intersecting slices is saved at
		//_offsetIndex = _offset;
		//_offset += sizeof(ushort16); // jump to the next position beyond the one that stores the number of slices

		if ( _percentCounter == m_1percent)
		{
			//std::cout << "1%: " << m_1percent << std::endl;
			emit percentEvent();
			QCoreApplication::processEvents();
			_percentCounter = 0;
		}

		// number of elements
		*m_byteCounter += sizeof(unsigned short);


		// compute the world coordinates from the voxel coordinates
		const Eigen::Vector3f _local(static_cast<float>(_x)*m_physicalX+m_physicalX/2.0,static_cast<float>(_y)*m_physicalY+m_physicalY/2.0,static_cast<float>(_z)*m_physicalZ+m_physicalZ/2.0);


		const Eigen::Vector3f _world = ( m_boundingBoxTransformation *_local );


		// now we can look for the scanlines in vicinity

		// pop the minimum distance bin vector from the rotation queue

#ifdef USE_USHORT16
		std::vector<Scanline<unsigned short>*>* _usScanlineQueue = _rotQueue.pop();
#endif
#ifdef USE_FLOAT32 || USE_FLOAT5X32
		std::vector<Scanline<float>*>* _usScanlineQueue = _rotQueue.pop();
#endif
		int _oldNumberOfScanlinesInQueue = _usScanlineQueue->size();

		// for each slice within the vector returned form the queue do a backward warping
		while( _usScanlineQueue->size() > 0 )
		{
			StandardMeasurementCollection<IntensityAndDistance>* _smc = NULL;

			#ifdef USE_USHORT16
				Scanline<unsigned short> *_scanline = _usScanlineQueue->back();
			#endif

			#ifdef USE_FLOAT32 || USE_FLOAT5X32
				Scanline<float> *_scanline = _usScanlineQueue->back();
			#endif


			_usScanlineQueue->pop_back();

			/*if (_scanline->getUniqueIndex() == 1968754124)
			int _debug = 5;
			*/


			// determine the orthogonal distance from the voxel to the ultrasound plane
			float _lineParam = 0.0f;
			const float _distance = _scanline->distancePointToScanline(_world, _lineParam);
			
	//		if ( _distance > 0.25 && _distance < 0.4 && _lineParam >= 0 && _lineParam <= 1.0)
	//					std::cout << "DISTANCE: " << _distance << std::endl;
			// HACK to avoid overflow
			// the CONDITION: _intersectionCounter >= std::numeric_limits<unsigned short>::max()-2
			if ( _distance <= _threshold && _lineParam >= 0 && _lineParam <= 1.0 && _intersectionCounter <= std::numeric_limits<unsigned short>::max()-2)
			{
				_intersectionCounter++;
			
				//_smc = new StandardMeasurementCollection<IntensityAndDistance>();
				_smc = &_smcBlock[_smcCounter++];
				_smc->count = 0;
				// updating the byte counter (number of elements, direction index of scanline)
				*m_byteCounter += sizeof(unsigned int) + sizeof(unsigned short);

				// set the direction of the scanline to the measurement collection

				ScanlineIndex _di;
				_di.scanlineID = _scanline->getBeamIndex();
				_di.imageID = _scanline->getRfImageIndex();
				_di.fileID = _scanline->getFileIndex();


				

				_smc->scanlineIndex = (_di);
				_smc->uniqueID = _scanline->getUniqueIndex();

				// and push it into the temporary queue, as the point is not to be considered again for the same voxel
				_rotQueue.push(_scanline, _distance);

				// now collect the data points along the scanline

				// at first we have to find the limiting points (within distance)

				// compute the maximum distance we can go from the intersection point still being in the distance range

				float _maxDistanceOnScanline = sqrt(_threshold*_threshold - _distance*_distance);
				int _samplesRange = static_cast<int>(_scanline->getUltrasoundSettings()->samplesPerMM * _maxDistanceOnScanline);
				int _intersectionIndex = static_cast<int>(_lineParam * _scanline->getUltrasoundSettings()->samplesAlongScanline);


				int _minDex = std::max(0, _intersectionIndex-_samplesRange);
				int _maxDex = std::min(_scanline->getUltrasoundSettings()->samplesAlongScanline-1, _intersectionIndex+_samplesRange);

				// compute the range of the line parameter (of the scanline ray)
				// the minimum line parameter
				float _deltaLineParam1 = (_intersectionIndex-_minDex)*_scanline->getUltrasoundSettings()->mmPerSample/static_cast<float>(_scanline->getLength());
				// the maximum line parameter
				float _deltaLineParam2 = (_maxDex-_intersectionIndex)*_scanline->getUltrasoundSettings()->mmPerSample/static_cast<float>(_scanline->getLength());

				// now grab the values!
				// data conversion

#ifdef USE_USHORT16
				unsigned short *_ptr = &_scanline->getData()[0];
				_ptr += _minDex;
#endif

#ifdef USE_FLOAT32
				float *_ptr = &_scanline->getData()[0];
				_ptr += _minDex;
#endif

#ifdef USE_FLOAT5X32
				float *_ptr = &_scanline->getData()[0];
				_ptr += 5*_minDex;
#endif

				// from upper the beam to the intersection
				int _index = _minDex;
				// now divide the line parameter in the number of steps (measurements)
				float _deltaLineParamPerIndex = 0.0;
				if (_intersectionIndex-_index-1 > 0)
					_deltaLineParamPerIndex = _deltaLineParam1/static_cast<float>(_intersectionIndex-_index-1);
				else
					_deltaLineParamPerIndex = _deltaLineParam1;

				float _tmpLineParam = _lineParam - _deltaLineParam1 ;
				while(_index <= _intersectionIndex)
				{
					// compute the distance from the pixel to the voxel center
					float _tmpDistance = sqrt(pow((_intersectionIndex - _index) * _scanline->getUltrasoundSettings()->mmPerSample,2.0f)+pow(_distance,2.0f));

					// set the measured datata
					IntensityAndDistance _iad;

#ifdef USE_USHORT16
					_iad.intensity = *_ptr;
					// now go to the next;
					_ptr++;
					*m_byteCounter += sizeof(unsigned short);
#endif

#ifdef USE_FLOAT32
					_iad.intensity = *_ptr;
					// now go to the next;
					_ptr++;
					*m_byteCounter += sizeof(float);
#endif

#ifdef USE_FLOAT5X32
					_iad.mrf[0] = *_ptr; _ptr++;
					_iad.mrf[1] = *_ptr; _ptr++;
					_iad.mrf[2] = *_ptr; _ptr++;
					_iad.mrf[3] = *_ptr; _ptr++;
					_iad.mrf[4] = *_ptr; _ptr++;
					*m_byteCounter += sizeof(float)*5;
#endif
					#ifdef SAVE_DISTANCE_SCANLINE
					
#ifdef QUANTIZE_DISTANCE
					_iad.distance = _nr.quantize8BIT(_tmpDistance, m_voxelSize);
					*m_byteCounter += sizeof(char); // size of the distance element

#else
					_iad.distance = _tmpDistance;//compactFloat(_tmpDistance);
					*m_byteCounter += sizeof(float); // size of the distance element
#endif	
					// data conversion
					_iad.scanlineIndex = _index;
					*m_byteCounter += sizeof(unsigned short);
					//_iad.scanlineIndex = static_cast<float>(_index) * _scanline->getUltrasoundSettings()->mmPerSample; 
					#endif
					// distance  + scanlineIndex (each float)
					//*m_byteCounter += sizeof(float);
					

					// push the measurement data into the measurement collection
					//_smc->data.push_back(_iad);
//std::cout << _smc->count << std::endl;
					_smc->data[_smc->count++] = _iad;


					// now go to the next;
					//_ptr++;
					_index++;

					// now update the line parameter
					_tmpLineParam += _deltaLineParamPerIndex;
				}

				// from below the intersection up to the intersection
				_ptr = &_scanline->getData()[0];
				_ptr += _maxDex;
				_index = _maxDex;
				// now divide the line parameter in the number of steps (measurements)
				if ( _index-_intersectionIndex-1 > 0)
					_deltaLineParamPerIndex = _deltaLineParam2/static_cast<float>(_index-_intersectionIndex-1);
				else
					_deltaLineParamPerIndex = _deltaLineParam2;

				_tmpLineParam = _lineParam + _deltaLineParam2;
				while(_index > _intersectionIndex)
				{
					// compute the distance from the pixel to the voxel center
					float _tmpDistance = sqrt(pow((_index - _intersectionIndex) * _scanline->getUltrasoundSettings()->mmPerSample,2.0f)+pow(_distance,2.0f));

					// set the measured datata
					IntensityAndDistance _iad;

#ifdef USE_USHORT16
					_iad.intensity = *_ptr;
					*m_byteCounter += sizeof(unsigned short);
					_ptr--;
#endif

#ifdef USE_FLOAT32
					_iad.intensity = *_ptr;
					*m_byteCounter += sizeof(float);
					_ptr--;
#endif

#ifdef USE_FLOAT5X32
					_iad.mrf[0] = *_ptr; _ptr--;
					_iad.mrf[1] = *_ptr; _ptr--;
					_iad.mrf[2] = *_ptr; _ptr--;
					_iad.mrf[3] = *_ptr; _ptr--;
					_iad.mrf[4] = *_ptr; _ptr--;
					*m_byteCounter += sizeof(float)*5;
#endif

					#ifdef SAVE_DISTANCE_SCANLINE

#ifdef QUANTIZE_DISTANCE
					_iad.distance = _nr.quantize8BIT(_tmpDistance, m_voxelSize);
					*m_byteCounter += sizeof(char); // size of the distance element

#else
					_iad.distance = _tmpDistance;//compactFloat(_tmpDistance);
					*m_byteCounter += sizeof(float); // size of the distance element
#endif					
					// TESTWISE scanlineIndex will containe the line parameter
					//_iad.scanlineIndex = _tmpLineParam;
					_iad.scanlineIndex = _index;
					//_iad.scanlineIndex =  static_cast<float>(_index) * _scanline->getUltrasoundSettings()->mmPerSample; 
					*m_byteCounter += sizeof(unsigned short);
					#endif
					// distance + intensity + scanlineIndex (each float)
					//*m_byteCounter += sizeof(float) * 3;

					// push the measurement data into the measurement collection
					//_smc->data.push_back(_iad);
					_smc->data[_smc->count++] = _iad;


		// now go to the next;
					_index--;

						// now update the line parameter
					_tmpLineParam -= _deltaLineParamPerIndex;
				}

			}
			else // too far away
			{

				// and push it into the temporary queue, as the point is not to be considered again for the same voxel
				_rotQueue.push(_scanline, _distance);
			}

			// now push the measurement collection to the vector

			
			
			// we don't need NO null pointer
			if ( _smc != NULL )
			//if ( _smc.data.size() > 0 ) 
			{
				// prior to saving get rid of the excess capacity to the swap trick
			// shrink to fit the std::vector
			//std::vector<IntensityAndDistance>(_smc->data).swap(_smc->data);

				
				m_measurementsCollectionManager->setData(_smc,_x, _y, _z);
				//_intersectingSlices++;
				}



			// TJKWRITE
			

			
			
		}

		// check if there is an element overflow, if yes remove elements that are too far away
		//if (_intersectionCounter >= std::numeric_limits<unsigned short>::max()-1)
	//	{
		//	m_measurementsCollectionManager->getData(_x,_y,_z);
	//	}
//m_measurementsCollectionManager->trim(_x, _y, _z);
	
		// instead of delete put empty vectors back to rotation queue, so we save new operations
		_usScanlineQueue->reserve(_oldNumberOfScanlinesInQueue);
		_rotQueue.recycle(_usScanlineQueue);

		//std::cout << "CMP: " << _traversalCount;


		//std::cout << "End: " << _traversalCount-1 << std::endl;
		//std::cout << "Finisheded Thread No. " << m_index << "/ Comp. "<< _traversalCount << std::endl;
	}

	
	/*std::cout << "MIN: " << _minX << "\t" << _minY << "\t" << _minZ << std::endl;
	std::cout << "MAX: " << _maxX << "\t" << _maxY << "\t" << _maxZ << std::endl;
	std::cout << "MIN DIST: " << _minDist << "\t" << _minBeam.imageID << "\t" << _minBeam.scanlineID << std::endl << std::endl;
	std::cout << m_boundingBoxTransformation.matrix() << std::endl;
	std::cout << std::endl;*/
	std::cout << "Finished Thread No. : " << m_index << std::endl; 

}


void CompoundRFVolumeThread::setFilename(std::string filename)
{
	m_filename = filename;
}

void CompoundRFVolumeThread::setDataTarget(StandardMeasurementCollectionManager<IntensityAndDistance> *ptr)
{
	m_measurementsCollectionManager = ptr;
}

void CompoundRFVolumeThread::setPhysicalOffset(float offsetX, float offsetY, float offsetZ)
{
	m_physicalOffsetX = offsetX; 
	m_physicalOffsetY = offsetY;
	m_physicalOffsetZ = offsetZ;
}

void CompoundRFVolumeThread::setIndex(int index)
{
	m_index = index;
}

void CompoundRFVolumeThread::setByteCounter(__int64 *byteCounter)
{
	m_byteCounter = byteCounter;
}


void CompoundRFVolumeThread::setDataSize(const int dimX, const int dimY, const int dimZ, const float physicalX, const float physicalY, const float physicalZ, BBExtent extent, int offsetX, int offsetY, int offsetZ)
{
	m_voxelsX = dimX;
	m_voxelsY = dimY;
	m_voxelsZ = dimZ;

	m_physicalX = physicalX;
	m_physicalY = physicalY;
	m_physicalZ = physicalZ;

	m_voxelSize = std::min(std::min(physicalX,physicalY),physicalZ);
	m_extent = extent;

	m_offsetX = offsetX;
	m_offsetY = offsetY;
	m_offsetZ = offsetZ;
}


void CompoundRFVolumeThread::setUserConfiguration(int maxDistanceScalar, int maxDistanceUnit, int _1percent )
{
	m_maxDistanceScalar = maxDistanceScalar;
	m_maxDistanceUnit = maxDistanceUnit;
	m_1percent = _1percent;
}

void CompoundRFVolumeThread::setData(std::vector<Scanline<unsigned short>*>* scanlineVector)
{
#ifdef USE_USHORT16
	m_scanlineVector = scanlineVector;
#endif
}

void CompoundRFVolumeThread::setData(std::vector<Scanline<float>*>* scanlineVector)
{
#ifdef USE_FLOAT32 || USE_FLOAT5X32
	m_scanlineVector = scanlineVector;
#endif
}


void CompoundRFVolumeThread::setTransformation(Eigen::Transform3f& boundingBoxTransformation)
{
	m_boundingBoxTransformation = boundingBoxTransformation;
}


// UI RF

void Compounding::loadRFData()
{
	

	if ( m_xmlTools ) delete m_xmlTools;
	m_xmlTools = new XMLTools();
	// clear the name items of the tracking target combo box
	ui.probeNameRFComboBox->clear();
	ui.referenceNameRFComboBox->clear();


	std::vector<std::string> _targetNames;
	// use "ReferenceTarget" as the moving reference coordinate system, to cope with patient motion 
	bool _foundReference = false;
	bool _foundProbe = false;
	if (m_dataTools.getTargetNamesFromFile(m_rfTrackingFileName. toLatin1().data(), _targetNames) )
	{
		for(int i=0;i<_targetNames.size();i++)
		{
			ui.probeNameRFComboBox->addItem(QString(_targetNames[i].c_str()));
			if (_targetNames[i].find("Probe") != std::string::npos  || _targetNames[i].find("Transducer") != std::string::npos )
			{
				_foundProbe = true;
				ui.probeNameRFComboBox->setCurrentIndex(ui.probeNameRFComboBox->count()-1);
			}
			ui.referenceNameRFComboBox->addItem(QString(_targetNames[i].c_str()));
			if (_targetNames[i].find("Reference") != std::string::npos  || _targetNames[i].find("reference") != std::string::npos )
			{
				_foundReference = true;
				ui.referenceNameRFComboBox->setCurrentIndex(ui.referenceNameRFComboBox->count()-1);
			}
		}
	}

	// if we found a reference target ID, let's assume we use it
	if ( _foundReference == true)
		ui.referenceTargetRFCheckBox->setChecked(true);
	else
		ui.referenceTargetRFCheckBox->setChecked(false);


	/*	bool _doppler = false;

	// identify the sender
	QPushButton *_tempButton = (QPushButton *)sender();
	if ( strcmp(_tempButton->objectName(). toLatin1().data(),"LoadDopplerButton") == 0 )
	{
	// we want to reconstruct Doppler data, this need to be handled as a special case
	_doppler = true;
	}*/

	//if ( !_doppler )
	ui.initializeRFButton->setDisabled(false);
	//else
	//ui.InitializeDopplerButton->setDisabled(false);
}

void Compounding::selectRFFile()
{
	// initialisation of the output values
	//QString configurationFilter  = tr("Configuration files (*.xml)");		
	QStringList filterList; // = ((QStringList)configurationFilter);
	filterList << ".xml";

	// creation of a new file dialog
	QFileDialog* fileDialog = new QFileDialog(this);
	fileDialog->setWindowTitle("Loading Data...");
	fileDialog->setNameFilters(filterList);

	// check if we previously loaded some data from a specific directory, if not use the one from the XML file
	if ( m_lastDataRFPath.length() == 0 )
	{
		std::string _path;
		m_xmlConfiguration.handleData(false,_path, std::string(RF_DATA_DIRECTORY_TAG));
		m_lastDataRFPath = QString(_path.c_str());
	}

	QString fileName = fileDialog->getOpenFileName(0,"Select file ...", m_lastDataRFPath,"Configuration files (*.xml)");
	std::string _temp = fileName. toLatin1().data();
	std::string _path = _temp.substr( 0, _temp.find_last_of( '/' ) +1 );
	m_guiRFDataReady = !fileName.isEmpty();
	if(!fileName.isEmpty())
	{
		// write the directory used last time to the configuration file to remember it for the next time
		m_xmlConfiguration.handleData(true,_path, std::string(RF_DATA_DIRECTORY_TAG));
		m_lastDataRFPath = QString(_path.c_str());
		ui.RFXMLFileEdit->setText(fileName);
		m_rfDataFileName = fileName;

		// store the filename in the XML, so we can use it next time
		m_xmlConfiguration.handleData(true, _temp, std::string(RF_FILENAME_TAG));

		std::string _temp = fileName. toLatin1().data();
		int _idx1 =  _temp.find_last_of( '/' );
		int _idx2 = _temp.find_first_of( '.' );
		std::string _name = _temp.substr( _idx1+1, (_idx2-_idx1)  )+"RFC.xml";
		ui.outputRFFileName->setText(QString(_name.c_str()));
	}

	guiUpdated();
}

void Compounding::compoundRF()
{
	determineVoxelSizeRF();

	// refresh the widget once every minute (to speed up computation)
	glWidget->stopTimer();

	// special handling of Doppler compounding

	// compound the volume with specified compounding mode such as gaussian distance weighting, nearest neighbor, mean
	compoundRFVolume((COMPOUNDING_MODE)ui.reconstructionModeComboBox->currentIndex ());
}


void Compounding::attachRFData()
{
	bool _doppler = false;

	// identify the sender
	QPushButton *_tempButton = (QPushButton *)sender();
	if ( strcmp(_tempButton->objectName(). toLatin1().data(),"InitializeDopplerRFButton") == 0 )
	{
		// we want to reconstruct Doppler data, this need to be handled as a special case
		_doppler = true;
	}

	UltrasoundSettings _ultrasoundSettings;

	loadSettingFromXMLFile(m_rfCalibrationFileName, _ultrasoundSettings);

	// LEGACY SOLUTION
	m_calibrationMatrixVector.push_back(_ultrasoundSettings.calibrationMatrix);

	// for RF data processing we require calibration matrix from apex
	Eigen::Matrix4f _corner_T_apex;
	_corner_T_apex << 1.0, 0.0, 0.0, _ultrasoundSettings.apexX * _ultrasoundSettings.scaleX,  0.0, 1.0, 0.0, _ultrasoundSettings.apexY * _ultrasoundSettings.scaleY,  0.0, 0.0, 1.0, 0.0,  0.0, 0.0, 0.0, 1.0;
	//m_matCalibrationApex = _ultrasoundSettings.calibrationMatrix * _corner_T_apex; // NEEDS TO BE PUT INTO ULTRASOUND SETTINGS




	// TODO: obtain m_rfRFDfilename from the xml file!

	// Now load all the RF specific variables from the XML file

	//if ( m_rfFileStream) delete m_rfFileStream;
	// count the number of frames before adding new ones
	int _oldRFFrames = 0;
	for (int j=0;j<m_rfFileStreamVector.size();j++)
	{
		_oldRFFrames += m_rfFileStreamVector.at(j)->getRecordingInfo().frames;
	}

/*	RFFileStream *_rfFileStream = new RFFileStream();
	m_rfFileStreamVector.push_back(_rfFileStream);

	_rfFileStream->openFile(m_rfDataFileName. toLatin1().data());
	*/

	RFFileStream *_rfFileStream = new RFFileStream();
	m_rfFileStreamVector.push_back(_rfFileStream);

	_rfFileStream->openFile(m_rfDataFileName. toLatin1().data(), DATA_TYPE::USHORT16);

	// XML file info
	m_imageFilesXML.push_back(m_rfDataFileName. toLatin1().data());

	std::string _temp = m_rfDataFileName. toLatin1().data();
	std::string _path = _temp.substr( 0, _temp.find_last_of( '/' ) +1 );


	// update the GUI 
	int _rfFrames = 0;
	for (int j=0;j<m_rfFileStreamVector.size();j++)
	{
		_rfFrames += m_rfFileStreamVector.at(j)->getRecordingInfo().frames;
	}
	ui.RFImageLabelScrollBar->setMaximum(_rfFrames-1);
	{
		std::stringstream ss; 
		ss  << "Frame: 0" << " / " << _rfFrames;
		ui.labelFrameRF->setText( ss.str().c_str() );
	}

	// set variables
	_ultrasoundSettings.penetrationDepth = _rfFileStream->getRecordingInfo().penetrationDepth;
	_ultrasoundSettings.scanlines = _rfFileStream->getRecordingInfo().width;
	_ultrasoundSettings.samplesAlongScanline = _rfFileStream->getRecordingInfo().height;
	_ultrasoundSettings.rfFrameSizeElements = _ultrasoundSettings.scanlines*_ultrasoundSettings.samplesAlongScanline; // NEEDS TO BE PUT INTO ULTRASOUND SETTINGS
	#ifdef USE_USHORT16

	_ultrasoundSettings.rfFrameSizeBytes = sizeof(unsigned short)*_ultrasoundSettings.scanlines*_ultrasoundSettings.samplesAlongScanline; // NEEDS TO BE PUT INTO ULTRASOUND SETTINGS
#endif

	#ifdef USE_FLOAT32

	_ultrasoundSettings.rfFrameSizeBytes = sizeof(float)*_ultrasoundSettings.scanlines*_ultrasoundSettings.samplesAlongScanline; // NEEDS TO BE PUT INTO ULTRASOUND SETTINGS
#endif


	_ultrasoundSettings.mmPerSample = (float)_ultrasoundSettings.penetrationDepth / (float)_ultrasoundSettings.samplesAlongScanline;
	_ultrasoundSettings.samplesPerMM = (float)_ultrasoundSettings.samplesAlongScanline / (float)_ultrasoundSettings.penetrationDepth;
	loadTimeStampsFromFile(QString((_path+_rfFileStream->getTimeStampFileName()).c_str()), true, true);
	loadTrackingDataFromFile(m_rfTrackingFileName, true, ui.referenceTargetRFCheckBox->isChecked(), true);
	m_ultrasoundSettings.push_back(_ultrasoundSettings);
	// TODO: here must be a bug, wrong extent and wrong ROI etc
	if ( ui.MasterFrameRFRadio->isChecked()) 
	{
		determineMinimumBoundingBoxfromMasterFrameRF();
	}

	else if ( ui.AutomaticRFRadio->isChecked())
	{
		determineMinimumBoundingBoxfromRF();
		m_masterFrameXML = -1; 
	}

	else if ( ui.ManualRFRadio->isChecked())
	{
		determineManualBoundingBoxfromRF();
		m_masterFrameXML = -1; 
	}

	// set the parameters of the region of interest within the visualization module, such that it can be drawn properly
	m_slicesVisualization->setImageSize(m_roiMinX, m_roiMinY, m_roiMaxX, m_roiMaxY, m_scaleX, m_scaleY);



	m_slicesVisualization->setTrackingData(m_trackingData);


	m_slicesVisualization->defineOrientedBoundingBox(m_extent, m_boundingBoxTransformation);

	// compute the scanlines
	for (int i=_oldRFFrames;i < m_matrixDataTracking.size(); i++)
	{
#ifdef USE_USHORT16
		StandardRF<unsigned short> *_data = new StandardRF<unsigned short>();
#endif

#ifdef USE_FLOAT32 || USE_FLOAT5X32
		StandardRF<float> *_data = new StandardRF<float>();
#endif
		_data->setTransformation(m_matrixDataTracking[i]);
		_data->setNumberOfSamples(_rfFileStream->getRecordingInfo().height);
		_data->setNumberOfScanlines(_rfFileStream->getRecordingInfo().width);
		m_rfDataVector.push_back(_data);
	}
	m_fileIndex = (int)m_fileIndex +1;
	m_settingsIndex = (int)m_settingsIndex+1;
	computeScanlines(m_settingsIndex, m_fileIndex );

	suggestVoxelSize();
	ui.voxelSizeSlider->setEnabled(true);
	ui.voxelSizeSlider->setValue(m_voxelSize/0.025);
	determineVoxelSizeRF();

	ui.compoundRFButton->setEnabled(true);

}

void Compounding::selectConfigurationFileRF()
{
	// initialisation of the output values
	//QString configurationFilter  = tr("Configuration files (*.xml)");		
	QStringList filterList;// = ((QStringList)configurationFilter);
	filterList << "*.xml";

	// creation of a new file dialog
	QFileDialog* fileDialog = new QFileDialog(this);
	fileDialog->setWindowTitle("Loading Data...");
	fileDialog->setNameFilters(filterList);

	// check if we previously loaded some data from a specific directory, if not use the one from the XML file
	if ( m_lastCalibrationRFPath.length() == 0 )
	{
		std::string _path;
		m_xmlConfiguration.handleData(false,_path, std::string(RF_CALIBRATION_DIRECTORY_TAG));
		m_lastCalibrationRFPath = QString(_path.c_str());
	}

	QString fileName = fileDialog->getOpenFileName(0,"Select file ...", m_lastCalibrationRFPath,"Configuration files (*.xml)");
	std::string _temp = fileName. toLatin1().data();
	std::string _path = _temp.substr( 0, _temp.find_last_of( '/' ) +1 );
	m_guiRFCalibrationReady = !fileName.isEmpty();
	if(!fileName.isEmpty())
	{
		// write the directory used last time to the configuration file to remember it for the next time
		m_xmlConfiguration.handleData(true,_path, std::string(RF_CALIBRATION_DIRECTORY_TAG));
		m_lastCalibrationRFPath = QString(_path.c_str());
		ui.calibrationRFFileEdit->setText(fileName);
		m_rfCalibrationFileName = fileName;

		// store the filename in the XML, so we can use it next time
		m_xmlConfiguration.handleData(true, _temp, std::string(RF_CALIBRATION_FILENAME_TAG));
	}

	guiUpdated();
}

void Compounding::selectTrackingFileRF()
{
	// initialisation of the output values
	//QString trackingFilter  = tr("Tracking files (*.txt)");		
	QStringList filterList;// = ((QStringList)trackingFilter);
	filterList << "*.txt";

	// creation of a new file dialog
	QFileDialog* fileDialog = new QFileDialog(this);
	fileDialog->setWindowTitle("Loading Data...");
	fileDialog->setNameFilters(filterList);


	// check if we previously loaded some data from a specific directory, if not use the one from the XML file
	if ( m_lastTrackingRFPath.length() == 0 )
	{
		std::string _path;
		m_xmlConfiguration.handleData(false,_path, std::string(RF_TRACKING_DIRECTORY_TAG));
		m_lastTrackingRFPath = QString(_path.c_str());
	}


	QString fileName = fileDialog->getOpenFileName(0,"Select file ...", m_lastTrackingRFPath,"Tracking files (*.txt)");
	std::string _temp = fileName. toLatin1().data();
	std::string _path = _temp.substr( 0, _temp.find_last_of( '/' ) +1 );
	m_guiRFTrackingReady = !fileName.isEmpty();
	if(!fileName.isEmpty())
	{
		// write the directory used last time to the configuration file to remember it for the next time
		m_xmlConfiguration.handleData(true,_path, std::string(RF_TRACKING_DIRECTORY_TAG));
		m_lastTrackingRFPath = QString(_path.c_str());
		ui.trackingRFFileEdit->setText(fileName);
		m_rfTrackingFileName = fileName;

		// store the filename in the XML, so we can use it next time
		m_xmlConfiguration.handleData(true, _temp, std::string(RF_TRACKING_FILENAME_TAG));
	}

	guiUpdated();
}

void Compounding::setLastSettingsRF()
{
	// store the filename in the XML, so we can use it next time
	std::string _temp;
	m_xmlConfiguration.handleData(false, _temp, std::string(RF_CALIBRATION_FILENAME_TAG));
	m_rfCalibrationFileName = QString(_temp.c_str());

	m_xmlConfiguration.handleData(false, _temp, std::string(RF_TRACKING_FILENAME_TAG));
	m_rfTrackingFileName = QString(_temp.c_str());

	m_xmlConfiguration.handleData(false, _temp, std::string(RF_FILENAME_TAG));
	m_rfDataFileName = QString(_temp.c_str());

	m_xmlConfiguration.handleData(false, _temp, std::string(RF_TARGET_FILENAME_TAG));
	m_rfRFDfilename = QString(_temp.c_str());



	// update the GUI
	ui.RFXMLFileEdit->setText( m_rfDataFileName);
	ui.trackingRFFileEdit->setText( m_rfTrackingFileName);
	ui.calibrationRFFileEdit->setText( m_rfCalibrationFileName);

	{
		std::string _temp = m_rfDataFileName. toLatin1().data();
		int _idx1 =  _temp.find_last_of( '/' );
		int _idx2 = _temp.find_first_of( '.' );
		std::string _name = _temp.substr( _idx1+1, (_idx2-_idx1)  )+"RFC.xml";
		ui.outputRFFileName->setText(QString(_name.c_str()));
	}
	// indicate that data is ready OR not (and update the GUI accordingly)
	if ( QFile(m_rfCalibrationFileName).exists() && QFile(m_rfTrackingFileName).exists() && QFile(m_rfDataFileName).exists() )
	{
		m_guiRFCalibrationReady = true;
		m_guiRFDataReady = true;
		m_guiRFTrackingReady = true;
	}
	else
	{
		m_guiRFCalibrationReady = false;
		m_guiRFDataReady = false;
		m_guiRFTrackingReady = false;
	}

	guiUpdated();
}

void Compounding::selectOutputRF()
{
	// initialisation of the output values
	//QString configurationFilter  = tr("Output file (*.xml)");		
	QStringList filterList;// = ((QStringList)configurationFilter);
	filterList << "*.xml";

	// creation of a new file dialog
	QFileDialog* fileDialog = new QFileDialog(this);
	fileDialog->setWindowTitle("Save Data...");
	fileDialog->setNameFilters(filterList);

	// check if we previously loaded some data from a specific directory, if not use the one from the XML file
	if ( m_lastDataRFPath.length() == 0 )
	{
		std::string _path;
		m_xmlConfiguration.handleData(false,_path, std::string(RF_TARGET_DIRECTORY_TAG));
		m_lastDataRFPath = QString(_path.c_str());
	}

	QString fileName = fileDialog->getSaveFileName(0,"Select file ...", m_lastTargetRFPath,"Ouput file (*.xml)");
	std::string _temp = fileName. toLatin1().data();
	std::string _path = _temp.substr( 0, _temp.find_last_of( '/' ) +1 );
	m_guiRFDataReady = !fileName.isEmpty();
	if(!fileName.isEmpty())
	{
		// write the directory used last time to the configuration file to remember it for the next time
		m_xmlConfiguration.handleData(true,_path, std::string(RF_TARGET_DIRECTORY_TAG));
		m_lastTargetRFPath = QString(_path.c_str());
		ui.outputRFFileName->setText(fileName);
		m_rfRFDfilename = fileName;

		// store the filename in the XML, so we can use it next time
		m_xmlConfiguration.handleData(true, _temp, std::string(RF_TARGET_FILENAME_TAG));
	}

	guiUpdated();
}

void Compounding::postProcessRF(std::string postfix)
{
	NumericTools _nr;
	if  (ui.checkBoxLogCompression->checkState() == Qt::Checked)
		postfix = postfix+"-LOG";

	IntensityCompounding *_intensityCompounding;

	switch((COMPOUNDING_MODE)ui.reconstructionModeComboBox->currentIndex () )
	{
	case COMPOUNDING_MODE::GAUSSIAN: 
		{
			_intensityCompounding = new GaussianWeightedCompounding(atof(ui.gaussianSigmaEdit->text(). toLatin1().data()));
			break;
		}

	case COMPOUNDING_MODE::MEAN: 
		{
			_intensityCompounding = new MeanCompounding();
			break;
		}
	}

	std::cout << "Started post processing ... " << std::endl;

	XMLTools *_xmlTools = new XMLTools();

	_xmlTools->loadFile(ui.postProcessingTargetFilenameEdit->text(). toLatin1().data());

	std::string filename = ui.postProcessingTargetFilenameEdit->text(). toLatin1().data();

	std::string _path = filename.substr( 0, filename.find_last_of( '/' ) +1 );

	std::string _reconstructionFileName;
	_xmlTools->handleData(false,_reconstructionFileName,"ReconstructionFile");

	std::string _voxelsX, _voxelsY, _voxelsZ;
	_xmlTools->handleData(false, _voxelsX, "DataSize","VoxelsX");
	_xmlTools->handleData(false, _voxelsY, "DataSize","VoxelsY");
	_xmlTools->handleData(false, _voxelsZ, "DataSize","VoxelsZ");


	std::string _physicalX, _physicalY, _physicalZ;
	_xmlTools->handleData(false, _physicalX, "DataSize","PhysicalX");
	_xmlTools->handleData(false, _physicalY, "DataSize","PhysicalY");
	_xmlTools->handleData(false, _physicalZ, "DataSize","PhysicalZ");

	// we assume isotropic voxels

	float _voxelSize = atof(_physicalX.c_str()); 

	delete _xmlTools;

	std::fstream _reconstructionFile;

	_reconstructionFile.open(_path+_reconstructionFileName, std::ios::in  | std::ios::binary);

	int _volumeIndexStart[] = {0,2,4,6};
	int _volumeIndexEnd[] = {1,3,5,7};

#if DEBUG_TXT_OUTPUT
	std::ofstream _debugFile;
	_debugFile.open(_reconstructionFileName+".debug.txt");
#endif

	float _maxVal = -1000000000.0;
	float _minVal =  1000000000.0;
	
#ifdef  USE_USHORT16
	vtkSmartPointer<vtkImageData> _compoundingOutput = vtkSmartPointer<vtkImageData>::New();
	// set the parameters of the image

	//_compoundingOutput->SetScalarTypeToUnsignedChar(); 
	
	//_compoundingOutput->SetDimensions( atoi(_voxelsX.c_str()),  atoi(_voxelsY.c_str()),  atoi(_voxelsZ.c_str())); 
	_compoundingOutput->SetExtent(0,atoi(_voxelsX.c_str())-1,  0, atoi(_voxelsY.c_str())-1,  0, atoi(_voxelsZ.c_str())-1); 
	_compoundingOutput->SetOrigin(0.0, 0.0, 0.0);  
	_compoundingOutput->SetSpacing( atof(_physicalX.c_str()),  atof(_physicalY.c_str()),  atof(_physicalZ.c_str())); 
	_compoundingOutput->AllocateScalars(VTK_DOUBLE, 1);
	//_compoundingOutput->Update(); 
	// clear the memory
	//double *_tmpIntensityPtr = static_cast<double*>(_compoundingOutput->GetScalarPointer(0,0,0));
	//double *_tmpIntensityPtr = static_cast<double*>(_compoundingOutput->GetScalarPointer());
#endif

#ifdef USE_FLOAT32
	vtkImageData *_compoundingOutput = vtkImageData::New();
	// set the parameters of the image

	//_compoundingOutput->SetScalarTypeToUnsignedChar(); 
	_compoundingOutput->SetScalarTypeToFloat();
	_compoundingOutput->SetDimensions( atoi(_voxelsX.c_str()),  atoi(_voxelsY.c_str()),  atoi(_voxelsZ.c_str())); 
	_compoundingOutput->SetOrigin(0.0, 0.0, 0.0);  
	_compoundingOutput->SetSpacing( atof(_physicalX.c_str()),  atof(_physicalY.c_str()),  atof(_physicalZ.c_str())); 
	_compoundingOutput->Update(); 
	// clear the memory
	float *_tmpIntensityPtr = static_cast<float*>(_compoundingOutput->GetScalarPointer(0,0,0));
#endif

#ifdef USE_FLOAT5X32
	vtkImageData *_compoundingOutput[5]
	float *_tmpIntensityPtr[5];
	
	for(int k=0;k<5;k++){
	_compoundingOutput[k] = vtkImageData::New();
	// set the parameters of the image

	//_compoundingOutput->SetScalarTypeToUnsignedChar(); 
	_compoundingOutput->SetScalarTypeToFloat();
	_compoundingOutput->SetDimensions( atoi(_voxelsX.c_str()),  atoi(_voxelsY.c_str()),  atoi(_voxelsZ.c_str())); 
	_compoundingOutput->SetOrigin(0.0, 0.0, 0.0);  
	_compoundingOutput->SetSpacing( atof(_physicalX.c_str()),  atof(_physicalY.c_str()),  atof(_physicalZ.c_str())); 
	_compoundingOutput->Update(); 
	// clear the memory
	_tmpIntensityPtr[k] = static_cast<float*>(_compoundingOutput->GetScalarPointer(0,0,0));
	}
#endif
	int _percentCounter = 0;
	int _voxels = atoi(_voxelsZ.c_str()) * atoi(_voxelsX.c_str()) * atoi(_voxelsY.c_str());
	int _1percent = (int)_voxels/100.0;
	//for (int c=0; c < _voxels;c++)
	std::cout << "Writing "<< _voxels << " | " << atoi(_voxelsZ.c_str()) << " x " << atoi(_voxelsY.c_str()) << " x " << atoi(_voxelsX.c_str()) << std::endl;
	//double _c = 0.0;
	for (int _z=0; _z <  atoi(_voxelsZ.c_str()); _z++) {
	for (int _y=0; _y <  atoi(_voxelsY.c_str()); _y++) {
	for (int _x=0; _x <  atoi(_voxelsX.c_str()); _x++)
	{
		//_c = _c +1.0;
		_percentCounter++;

		if ( _percentCounter == _1percent)
		{
			updatePercentage();
			_percentCounter = 0;
		}

		// first value indicating the number of elements
		ushort16 _ushort16;
		_reconstructionFile.read(&_ushort16.byte.c[0], sizeof(ushort16));

		//float _maxVal = 0.0;


		std::vector<std::vector<float>> _intensities;
#ifdef USE_USHORT16
		std::vector<float> _tmp;
		_intensities.push_back(_tmp);
#endif
#ifdef USE_FLOAT32
		std::vector<float> _tmp;
		_intensities.push_back(_tmp);
#endif
#ifdef USE_FLOAT5X32
		for(int k=0;k<5;k++){
		std::vector<float> _tmp;
		_intensities.push_back(_tmp);}
#endif


		std::vector<float> _distances;

		// now read intensity and distance
		for(int i=0;i<_ushort16.value;i++)
		{
			ushort16 _ushort;
			_reconstructionFile.read(&_ushort.byte.c[0], sizeof(ushort16));

			// handling the direction index
			{
				uint32 _uint32;
				_reconstructionFile.read(&_uint32.byte.c[0], sizeof(uint32));

			}


			for (int j=0; j<_ushort.value;j++)
			{
				IntensityAndDistance _iad;
				//ushort16 _ushort16;
				
				ushort16 _scanlineIndex;

				/*_rfSubDataFile[_vol].read(&_ushort16.byte.c[0], sizeof(ushort16));
				_mergedFile.write(&_ushort16.byte.c[0], sizeof(ushort16));
				_byteCounter.value += sizeof(ushort16);*/
				float _distance;
#ifdef QUANTIZE_DISTANCE
				char _tmpdistance;
				_reconstructionFile.read(&_tmpdistance, sizeof(char));
				_distance = _nr.dequantize8BIT((unsigned char)_tmpdistance, _voxelSize);
#else
				float32 _tmpdistance;
				_reconstructionFile.read(&_tmpdistance.byte.c[0], sizeof(float32));
				_distance = _tmpdistance.value;
#endif

#ifdef USE_USHORT16
				ushort16 _intensity16;
				_reconstructionFile.read(&_intensity16.byte.c[0], sizeof(ushort16));
				_intensities[0].push_back(fabs((float)_intensity16.value));
#endif

#ifdef USE_FLOAT32
				float32 _intensity32;
				_reconstructionFile.read(&_intensity32.byte.c[0], sizeof(float32));
				_intensities[0].push_back(fabs((float)_intensity32.value));
#endif

#ifdef USE_FLOAT5X32
				for(int k=0;k<5;k++) {
				float32 _intensity32;
				_reconstructionFile.read(&_intensity32.byte.c[0], sizeof(float32));
				_intensities[k].push_back(fabs((float)_intensity32.value));
				}
#endif

				_reconstructionFile.read(&_scanlineIndex.byte.c[0], sizeof(ushort16));

				//if  (_float32.value == 255.0)
				//	_maxVal = 255.0;

				_distances.push_back(_distance);

			}

		}

		//if ( _maxVal >= 255.0 )
		///	std::cout << "Max: " << _maxVal << " \ " << _intensityCompounding->evaluate(&_intensities, &_distances) << std::endl;
#ifdef USE_USHORT16 
		// perform log compression ?
		if (m_logCompress) {
			double *_tmpIntensityPtr = static_cast<double*>(_compoundingOutput->GetScalarPointer(_x,_y,_z));
			_tmpIntensityPtr[0] = log10(1.0+_intensityCompounding->evaluate(&_intensities[0], &_distances));
		}
		else
		{
			double *_tmpIntensityPtr = static_cast<double*>(_compoundingOutput->GetScalarPointer(_x,_y,_z));
			_tmpIntensityPtr[0] = (1.0+_intensityCompounding->evaluate(&_intensities[0], &_distances));
		}
		//_tmpIntensityPtr++;
#endif

		#ifdef USE_FLOAT32
		// perform log compression ?
		if (m_logCompress) {
			*_tmpIntensityPtr = /*_maxVal; //*/(float)log10(1.0+_intensityCompounding->evaluate(&_intensities[0], &_distances));
		}
		else
		{
			*_tmpIntensityPtr = /*_maxVal; //*/(float)(1.0+_intensityCompounding->evaluate(&_intensities[0], &_distances));
		}
		_tmpIntensityPtr++;
#endif

#ifdef USE_FLOAT5X32
		for(int k=0;k<5;k++){
		// perform log compression ?
		if (m_logCompress) {
			*_tmpIntensityPtr[k] = /*_maxVal; //*/(float)log10(1.0+_intensityCompounding->evaluate(&_intensities[k], &_distances));
		}
		else
		{
			*_tmpIntensityPtr[k] = /*_maxVal; //*/(float)(1.0+_intensityCompounding->evaluate(&_intensities[k], &_distances));
		}
		_tmpIntensityPtr[k]++;
		}
#endif

		/*if ( *_tmpIntensityPtr  > _maxVal )
			_maxVal = *_tmpIntensityPtr;

		if ( *_tmpIntensityPtr < _minVal )
			_minVal = *_tmpIntensityPtr;
		*/
		

	}
	}}

	_reconstructionFile.close();

	std::cout << "Min: " << _minVal << " / Max: " << _maxVal << std::endl;

#if DEBUG_TXT_OUTPUT
	_debugFile.close();
#endif
	#ifdef USE_USHORT16


	

	
	



	vtkSmartPointer<vtkMetaImageWriter> _writer = vtkSmartPointer<vtkMetaImageWriter>::New();
	//vtkSmartPointer<vtkImageCast> _castFilter =
	//	vtkSmartPointer<vtkImageCast>::New();
	//_castFilter->SetInputData(_compoundingOutput);
	//_writer->SetInputConnection(_castFilter->GetOutputPort());
	_writer->SetInputData(_compoundingOutput);

	std::string _temp = ui.postProcessingTargetFilenameEdit->text(). toLatin1().data();
	std::string _xmlFilename = _temp.substr( 0, _temp.find_last_of( ".xml" )-3 ).c_str() + postfix;

	_writer->SetFileName((_xmlFilename+".mhd").c_str()); 

	_writer->SetCompression( false );
	_writer->Write();

	//_writer->Delete();
#endif

		#ifdef USE_FLOAT32
	vtkMetaImageWriter *_writer = vtkMetaImageWriter::New(); 
	_writer->SetInput(_compoundingOutput); 


	std::string _temp = ui.postProcessingTargetFilenameEdit->text(). toLatin1().data();
	std::string _xmlFilename = _temp.substr( 0, _temp.find_last_of( ".xml" )-3 ).c_str() + postfix;

	_writer->SetFileName((_xmlFilename+".mhd").c_str()); 

	_writer->SetCompression( false );
	_writer->Write();

	_writer->Delete();
#endif

#ifdef USE_FLOAT5X32
	for(int k=0;k<5;k++){
		std::stringstream _ss;
		_ss << _xmlFilename << i << ".mhd";
	vtkMetaImageWriter *_writer = vtkMetaImageWriter::New(); 
	_writer->SetInput(_compoundingOutput); 


	std::string _temp = ui.postProcessingTargetFilenameEdit->text(). toLatin1().data();
	std::string _xmlFilename = _temp.substr( 0, _temp.find_last_of( ".xml" )-3 ).c_str() + postfix;

	_writer->SetFileName(_ss.c_str()); 

	_writer->SetCompression( false );
	_writer->Write();

	_writer->Delete();
	}
#endif

	m_percentage = -1;
	updatePercentage();

	std::cout << "Finished post processing ... " << std::endl;

	delete _intensityCompounding;

}


void Compounding::selectPostProcessingFile()
{
	// initialisation of the output values
	//QString configurationFilter  = tr("Configuration files (*.xml)");		
	QStringList filterList;// = ((QStringList)configurationFilter);
	filterList << ".xml";

	// creation of a new file dialog
	QFileDialog* fileDialog = new QFileDialog(this);
	fileDialog->setWindowTitle("Loading Data...");
	fileDialog->setNameFilters(filterList);

	QString fileName = fileDialog->getOpenFileName(0,"Select file ...", "","Configuration files (*.xml)");
	std::string _temp = fileName. toLatin1().data();
	std::string _path = _temp.substr( 0, _temp.find_last_of( '/' ) +1 );

	if(!fileName.isEmpty())
	{
		ui.postProcessingTargetFilenameEdit->setText(fileName);
		ui.PostProcess->setEnabled(true);
	}
}

void Compounding::extractDataRF()
{
	// the number of features to be saved
	int _extractionFeatures = 4;

	if ( !ui.checkBoxExtractDepth->isChecked())
		_extractionFeatures--;
	if (!ui.checkBoxExtractDistance->isChecked())
		_extractionFeatures--;
	if (!ui.checkBoxExtractIntensity->isChecked())
		_extractionFeatures--;
	if (!ui.checkBoxExtractScanline->isChecked())
		_extractionFeatures--;

	NumericTools _nr;
	srand((unsigned)time(0)); 

	int _rangeX[] = {ui.ComboFromX->currentIndex()+1, ui.ComboToX->currentIndex()+1};
	int _rangeY[] = {ui.ComboFromY->currentIndex()+1, ui.ComboToY->currentIndex()+1};
	int _rangeZ[] = {ui.ComboFromZ->currentIndex()+1, ui.ComboToZ->currentIndex()+1};

	XMLTools *_xmlTools = new XMLTools();
	int _elements;


	std::string _temp = ui.postProcessingTargetFilenameEdit->text(). toLatin1().data();
	std::string _xmlFilename = _temp.substr( 0, _temp.find_last_of( ".xml" )-3 ).c_str();
	_xmlTools->loadFile(ui.postProcessingTargetFilenameEdit->text(). toLatin1().data());
	std::string _filename = (ui.postProcessingTargetFilenameEdit->text(). toLatin1().data());
	std::string _path = _filename.substr( 0, _filename.find_last_of( '/' ) +1 );

	std::string _reconstructionFilename, _indexFilename, _directionFilename;
	std::string _voxelsX, _voxelsY, _voxelsZ;
	_xmlTools->handleData(false, _voxelsX, "DataSize","VoxelsX");
	_xmlTools->handleData(false, _voxelsY, "DataSize","VoxelsY");
	_xmlTools->handleData(false, _voxelsZ, "DataSize","VoxelsZ");

	std::string _physicalX, _physicalY, _physicalZ;
	_xmlTools->handleData(false, _physicalX, "DataSize","PhysicalX");
	_xmlTools->handleData(false, _physicalY, "DataSize","PhysicalY");
	_xmlTools->handleData(false, _physicalZ, "DataSize","PhysicalZ");


	std::string _depthStr;
	_xmlTools->handleData(false, _depthStr, std::string(PENETRATION_DEPTH_TAG));

	// we assume isotropic voxels

	float _voxelSize = atof(_physicalX.c_str()); 


	int _voxels[] = {atoi(_voxelsX.c_str()), atoi(_voxelsY.c_str()),atoi(_voxelsZ.c_str())};

	_elements = atoi(_voxelsX.c_str()) * atoi(_voxelsY.c_str()) * atoi(_voxelsZ.c_str());

	_xmlTools->handleData(false,_reconstructionFilename, RF_RECONSTRUCTIONFILE_TAG);
	_xmlTools->handleData(false,_indexFilename, RF_INDEXFILE_TAG);
	_xmlTools->handleData(false,_directionFilename, RF_DIRECTIONFILE_TAG);

	delete _xmlTools;
	// TODO: need to save the number of Scanlines & RF images

	std::fstream _directionFile;
	//std::fstream _reconstructionFile;
	HANDLE _reconstructionFile;
	std::fstream _indexFile;

	// load the direction file, and insert all the directions into a vector for later-on fast handling
	_directionFile.open(_path+_directionFilename, std::ios::in  | std::ios::binary);

	// check file length and calculate the number of elements from the size
	_directionFile.seekg (0, std::ios::end);
	int _directionFileLength = _directionFile.tellg();
	_directionFile.seekg (0, std::ios::beg);
	int _directionElements = _directionFileLength/(6.0*sizeof(float));


	matlab::mxArray* _startPoints = matlab::mxCreateDoubleMatrix(3, _directionElements, matlab::mxREAL);
	matlab::mxArray* _endPoints = matlab::mxCreateDoubleMatrix(3, _directionElements, matlab::mxREAL);
	double *_dataStartPtr = mxGetPr(_startPoints);
	double *_dataEndPtr = mxGetPr(_endPoints);

	std::vector<Eigen::Vector3f> _directionVectorStart;
	std::vector<Eigen::Vector3f> _directionVectorEnd;
	for (int i=0; i<_directionElements;i++)
	{
		Eigen::Vector3f _tempVec;
		float32 _x, _y, _z;

		_directionFile.read(&_x.byte.c[0],sizeof(float32));
		_directionFile.read(&_y.byte.c[0],sizeof(float32));
		_directionFile.read(&_z.byte.c[0],sizeof(float32));

		_tempVec << _x.value, _y.value, _z.value;
		*_dataStartPtr = (double)(_x.value);
		_dataStartPtr++;
		*_dataStartPtr = (double)(_y.value);
		_dataStartPtr++;
		*_dataStartPtr = (double)(_z.value);
		_dataStartPtr++;


		_directionVectorStart.push_back(_tempVec);

		_directionFile.read(&_x.byte.c[0],sizeof(float32));
		_directionFile.read(&_y.byte.c[0],sizeof(float32));
		_directionFile.read(&_z.byte.c[0],sizeof(float32));

		_tempVec << _x.value, _y.value, _z.value;

		*_dataEndPtr = (double)(_x.value);
		_dataEndPtr++;
		*_dataEndPtr = (double)(_y.value);
		_dataEndPtr++;
		*_dataEndPtr = (double)(_z.value);
		_dataEndPtr++;

		_directionVectorEnd.push_back(_tempVec);
	}
	_directionFile.close();

	matlab::mxArray* _offset = matlab::mxCreateDoubleMatrix(3, 1, matlab::mxREAL);
	double *_offsetPtr = mxGetPr(_offset);
	*_offsetPtr = _rangeX[0];
	_offsetPtr++;
	*_offsetPtr = _rangeY[0];
	_offsetPtr++;
	*_offsetPtr = _rangeZ[0];
	_offsetPtr++;


	matlab::mxArray* _physical = matlab::mxCreateDoubleMatrix(3, 1, matlab::mxREAL);
	double *_physicalPtr = mxGetPr(_physical);
	*_physicalPtr = atof(_physicalX.c_str()); 
	_physicalPtr++;
	*_physicalPtr = atof(_physicalY.c_str()); 
	_physicalPtr++;
	*_physicalPtr = atof(_physicalZ.c_str()); 
	_physicalPtr++;


	matlab::mxArray* _depth = matlab::mxCreateDoubleMatrix(1, 1, matlab::mxREAL);
	double *_depthPtr = mxGetPr(_depth);
	*_depthPtr = atof(_depthStr.c_str());

	{
		matlab::MATFile *pmatInfo;
		pmatInfo = matlab::matOpen("extractionInfo.mat", "w");

		int status = matPutVariable(pmatInfo, "ExtractionOffset", _offset);
		if ( status != 0 )
			std::cout << "Failed writing offset data!" << std::endl;
		status = matPutVariable(pmatInfo, "VoxelSize", _physical);
		if ( status != 0 )
			std::cout << "Failed writing voxel size data!" << std::endl;
		status = matPutVariable(pmatInfo, "PenetrationDepth", _depth);
		if ( status != 0 )
			std::cout << "Failed writing penetration depth data!" << std::endl;
		status = matPutVariable(pmatInfo, "ScanlineStart", _startPoints);
		if ( status != 0 )
			std::cout << "Failed writing scanline start data!" << std::endl;
		status = matPutVariable(pmatInfo, "ScanlineEnd", _endPoints);
		if ( status != 0 )
			std::cout << "Failed writing scanline end data!" << std::endl;
		matlab::mxDestroyArray(_physical);
		matlab::mxDestroyArray(_offset);
		matlab::mxDestroyArray(_depth);
		matlab::mxDestroyArray(_startPoints);
		matlab::mxDestroyArray(_endPoints);

	}


	std::vector<int> _indicesVector;
	std::vector<Eigen::Vector3i> _posVector;
	// now compute the indices
	for(int z=_rangeZ[0]; z<=_rangeZ[1];z++)
	{
		for(int y=_rangeY[0]; y<=_rangeY[1];y++)
		{
			for(int x=_rangeX[0]; x<=_rangeX[1];x++)
			{
				int _index = z*(_voxels[0] * _voxels[1]) + y * _voxels[0] + x;
				_indicesVector.push_back(_index);
				_posVector.push_back(Eigen::Vector3i(x,y,z));
			}
		}
	}

	

	// now read the byte offsets
	std::vector<__int64> _offsetVector;
	_indexFile.open(_path+_indexFilename, std::ios::in  | std::ios::binary);

	if ( _indexFile.eof() || _indexFile.bad() || _indexFile.fail() )
		std::cout << "Error while opening file!" << std::endl;

	_offsetVector.push_back(0);
	for (int i=0;i<_elements-1;i++)
	{
		ulong64 _tmp;
		_indexFile.read(&_tmp.byte.c[0],sizeof(ulong64));
		_offsetVector.push_back(_tmp.value);
	}
	_indexFile.close();

	matlab::mwSize _cellSize[]={_rangeX[1]-_rangeX[0]+1,_rangeY[1]-_rangeY[0]+1,_rangeZ[1]-_rangeZ[0]+1};
	matlab::mxArray* _cellDataIntensity  = matlab::mxCreateCellArray(3, _cellSize);
	matlab::mxArray* _cellDataScanline  = matlab::mxCreateCellArray(3, _cellSize);
	matlab::mxArray* _cellDataDistance  = matlab::mxCreateCellArray(3, _cellSize);
	matlab::mxArray* _cellDataDirection  = matlab::mxCreateCellArray(3, _cellSize);
 


	/*{
	std::ofstream _outStream;
	_outStream.open("Extraction.txt");
	_outStream << "test!!" << std::endl;
	_outStream.close();
	}*/
	std::ofstream _outStream;
	_outStream.open(_xmlFilename+".Extraction.txt");
	if ( _outStream.fail() )
		std::cout << "Error while opening file!" << std::endl;
	
	//_reconstructionFile.open(_path+_reconstructionFilename, std::ios::in  | std::ios::binary);

	_reconstructionFile = CreateFileA( (_path+_reconstructionFilename).c_str(), GENERIC_READ ,  FILE_SHARE_READ, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL,  NULL);

	__int64 _lastOffset = 0;


	int _lineCounter = 0;
	// progress bar
	m_percentage = -1;
	updatePercentage();
	int _1percent = static_cast<int>(_indicesVector.size())/100.0;
	int _percentCounter = 0;

	bool _voxelSkip = ui.VoxelSubsamplingCheckBox->isChecked();
	int _voxelPercentage = 0;
	int _voxelKeep = 0;
	if ( _voxelSkip )
	{
		std::string _temp = ui.VoxelSubsamplingComboBox->currentText(). toLatin1().data();

		if ( _temp.compare("2nd") == 0)
		{
			_voxelKeep = 2;
		}
		else if ( _temp.compare("3rd") == 0)
		{
			_voxelKeep = 3;
		}
		else if ( _temp.compare("4th") == 0)
		{
			_voxelKeep = 4;
		}
		else if ( _temp.compare("5th") == 0)
		{
			_voxelKeep = 5;
		}
		else if ( _temp.compare("6th") == 0)
		{
			_voxelKeep = 6;
		}
		else if ( _temp.compare("7th") == 0)
		{
			_voxelKeep = 7;
		}
		else if ( _temp.compare("8th") == 0)
		{
			_voxelKeep = 8;
		}
		else if ( _temp.compare("9th") == 0)
		{
			_voxelKeep = 9;
		}
		else if ( _temp.compare("10th") == 0)
		{
			_voxelKeep = 10;
		}

	}
	bool _intensitySkip = ui.WithinVoxelSubsamplingCheckbox->isChecked();
	if ( _intensitySkip )
	{
		std::string _temp = ui.WithinVoxelSubsamplingComboBox->currentText(). toLatin1().data();

		if ( _temp.compare("90%") == 0)
		{
			_voxelPercentage = 9;
		}
		else if ( _temp.compare("80%") == 0)
		{
			_voxelPercentage = 8;
		}
		else if ( _temp.compare("70%") == 0)
		{
			_voxelPercentage = 7;
		}
		else if ( _temp.compare("60%") == 0)
		{
			_voxelPercentage = 6;
		}
		else if ( _temp.compare("50%") == 0)
		{
			_voxelPercentage = 5;
		}
		else if ( _temp.compare("40%") == 0)
		{
			_voxelPercentage = 4;
		}
		else if ( _temp.compare("20%") == 0)
		{
			_voxelPercentage = 2;
		}
		else if ( _temp.compare("10%") == 0)
		{
			_voxelPercentage = 1;
		}
	}

	bool _write = true;
	int _voxelCount = 0;

	for(int k=0; k<_indicesVector.size();k++ )
	{

		
		_CrtDumpMemoryLeaks();

		std::vector<unsigned short> _intensityVector;
		std::vector<unsigned char> _distanceVector;
		std::vector<unsigned short> _scanlineindexVector;
		std::vector<unsigned int> _directionindexVector;


		if ( _voxelSkip )
		{
			_voxelCount++;
			if ( _voxelCount == _voxelKeep )
			{
				_write = true;
				_voxelCount = 0;
			}
			else
				_write = false;
		}

		_percentCounter++;
		if ( _percentCounter == _1percent)
		{
			_percentCounter = 0;
			updatePercentage();
		}

		int _seekIndex = _indicesVector[k];
		Eigen::Vector3i _coordinates = _posVector[k];
		__int64 _seekPos = _offsetVector[_seekIndex] - _lastOffset;
		
		//_reconstructionFile.seekg( _seekPos, std::ios_base::cur);
		
		LARGE_INTEGER _pint;
		LARGE_INTEGER _filePos;
		_filePos.QuadPart = _seekPos;
		//if ( _seekPos > 0)
			if ( SetFilePointerEx(_reconstructionFile, _filePos, &_pint, FILE_CURRENT) == 0)
				std::cout << "Error seeking file!" << std::endl;


		// now read the data
		DWORD  dwBytesRead = 0;

		// first value indicating the number of elements
		ushort16 _ushort16;

		//_reconstructionFile.read(&_ushort16.byte.c[0], sizeof(ushort16));
		ReadFile(_reconstructionFile,&_ushort16.byte.c[0],sizeof(ushort16),&dwBytesRead,NULL);
		
		
		float _maxVal = 0.0;

		// now read intensity and distance
		for(int i=0;i<_ushort16.value;i++)
		{
			ushort16 _ushort;
			//_reconstructionFile.read(&_ushort.byte.c[0], sizeof(ushort16));
			ReadFile(_reconstructionFile,&_ushort.byte.c[0],sizeof(ushort16),&dwBytesRead,NULL);

			// handling the direction index

			uint32 _uint32;
			//_reconstructionFile.read(&_uint32.byte.c[0], sizeof(uint32));
			ReadFile(_reconstructionFile,&_uint32.byte.c[0], sizeof(uint32), &dwBytesRead,NULL);


			// save the scanline index and the number of occurences
			_directionindexVector.push_back(_uint32.value);
			_directionindexVector.push_back((unsigned int)(_ushort.value));

			for (int j=0; j<_ushort.value;j++)
			{
				IntensityAndDistance _iad;
				//ushort16 _ushort16;
				float _distance;
				
				ushort16 _scanlineIndex;

				/*_rfSubDataFile[_vol].read(&_ushort16.byte.c[0], sizeof(ushort16));
				_mergedFile.write(&_ushort16.byte.c[0], sizeof(ushort16));
				_byteCounter.value += sizeof(ushort16);*/
				
				//_reconstructionFile.read(&_distance.byte.c[0], sizeof(float32));
#ifdef  QUANTIZE_DISTANCE
				char _distanceTMP;
				ReadFile(_reconstructionFile,&_distanceTMP, sizeof(char), &dwBytesRead,NULL);
				_distance = _nr.dequantize8BIT((unsigned char)_distanceTMP, _voxelSize);
#else // 32-bit float distance
				float32 _distanceTMP;
				ReadFile(_reconstructionFile,&_distanceTMP.byte.c[0], sizeof(float32), &dwBytesRead,NULL);
				_distance = _distanceTMP.value;
#endif

#ifdef USE_USHORT16
				ushort16 _intensity;
				//_reconstructionFile.read(&_float32.byte.c[0], sizeof(float32));
				ReadFile(_reconstructionFile,&_intensity.byte.c[0], sizeof(ushort16), &dwBytesRead,NULL);
#endif

#ifdef USE_FLOAT32
				float32 _intensity;
				//_reconstructionFile.read(&_float32.byte.c[0], sizeof(float32));
				ReadFile(_reconstructionFile,&_intensity.byte.c[0], sizeof(float32), &dwBytesRead,NULL);
#endif

				//_reconstructionFile.read(&_scanlineIndex.byte.c[0], sizeof(float32));
				ReadFile(_reconstructionFile,&_scanlineIndex.byte.c[0], sizeof(ushort16), &dwBytesRead,NULL);

				//if ( _float32.value == 255.0 )
				//	std::cout << "HALLO" << std::endl;

				if ( (_write && _intensitySkip && (rand() % 11  ) < _voxelPercentage) || (_write && !_intensitySkip )) {


					_intensityVector.push_back(_intensity.value);
					_distanceVector.push_back(_distanceTMP);
					_scanlineindexVector.push_back(_scanlineIndex.value);
					//_directionindexVector.push_back(_uint32.value);

				//	_outStream << (float)_float32.value << "\t" << (float)_distance.value << "\t" << (float)_scanlineIndex.value << "\t" << _directionVectorStart[_uint32.value].x() << "\t" << _directionVectorStart[_uint32.value].y() << "\t" << _directionVectorStart[_uint32.value].z() << "\t" << _directionVectorEnd[_uint32.value].x() << "\t" << _directionVectorEnd[_uint32.value].y() << "\t" << _directionVectorEnd[_uint32.value].z() << std::endl;
					//_outStream << k << "\t" << _intensity.value << "\t" << _distance << "\t" << _scanlineIndex.value << "\t" << _directionVectorStart[_uint32.value].x() << "\t" << _directionVectorStart[_uint32.value].y() << "\t" << _directionVectorStart[_uint32.value].z() << "\t" << _directionVectorEnd[_uint32.value].x() << "\t" << _directionVectorEnd[_uint32.value].y() << "\t" << _directionVectorEnd[_uint32.value].z() << "\t" << _coordinates.x() << "\t" << _coordinates.y() << "\t" << _coordinates.z() << std::endl;

					_lineCounter++;
				}
			}

		}

		_outStream.flush();

		

		matlab::mwSize _arraySize[]={ _intensityVector.size()};
		matlab::mwSize _arraySizeDirectionVector[]={ _directionindexVector.size()};
	
		
		matlab::mxArray* _dataMatrixIntensity;
		unsigned short *_dataIntensityPtr;
		if (ui.checkBoxExtractIntensity->isChecked()) 
		{
			_dataMatrixIntensity = matlab::mxCreateNumericArray(1,_arraySize, matlab::mxUINT16_CLASS, matlab::mxREAL);
			_dataIntensityPtr = (unsigned short*) mxGetData(_dataMatrixIntensity); 
		}

		matlab::mxArray* _dataMatrixScanline;
		unsigned short *_dataScanlinePtr;
		if ( ui.checkBoxExtractDepth->isChecked()) 
		{
			_dataMatrixScanline = matlab::mxCreateNumericArray(1,_arraySize, matlab::mxUINT16_CLASS, matlab::mxREAL);
			_dataScanlinePtr = (unsigned short*) mxGetData(_dataMatrixScanline);
		}

		matlab::mxArray* _dataMatrixDistance;
		unsigned char *_dataDistancePtr;
		if (ui.checkBoxExtractDistance->isChecked()) 
		{
			_dataMatrixDistance = matlab::mxCreateNumericArray(1, _arraySize, matlab::mxUINT8_CLASS, matlab::mxREAL); 
			_dataDistancePtr = (unsigned char*) mxGetData(_dataMatrixDistance);
		}

		
		matlab::mxArray* _dataMatrixDirection;
		int *_dataDirectionPtr;
		if (ui.checkBoxExtractScanline->isChecked())
		{
			
			_dataMatrixDirection = matlab::mxCreateNumericArray(1, _arraySizeDirectionVector,  matlab::mxUINT32_CLASS, matlab::mxREAL);
			
			_dataDirectionPtr = (int*) mxGetData(_dataMatrixDirection);
			
		}
	
		 
		
		
		
		
		for (int u=0;u<_intensityVector.size();u++)
		{
			if (ui.checkBoxExtractIntensity->isChecked()) {
				*_dataIntensityPtr = _intensityVector[u];
				_dataIntensityPtr++; }

			if (ui.checkBoxExtractDistance->isChecked()) {
				*_dataDistancePtr = _distanceVector[u];
				_dataDistancePtr++; }

			if ( ui.checkBoxExtractDepth->isChecked())  {
				*_dataScanlinePtr = _scanlineindexVector.at(u);
				_dataScanlinePtr++; }

			}

		
		for (int u=0;u<_directionindexVector.size();u++)
{
			if (ui.checkBoxExtractScanline->isChecked()) {
				*_dataDirectionPtr =  _directionindexVector[u];
				_dataDirectionPtr++; }

	}	

		
		if (ui.checkBoxExtractIntensity->isChecked()) 
		{
			mxSetCell(_cellDataIntensity,k,_dataMatrixIntensity); 
		}
		//else matlab::mxDestroyArray(_cellDataIntensity);

		if ( ui.checkBoxExtractDepth->isChecked())  
		{
			mxSetCell(_cellDataScanline,k,_dataMatrixScanline); 
		}
		//else matlab::mxDestroyArray(_cellDataScanline);

		if (ui.checkBoxExtractDistance->isChecked()) 
		{
			mxSetCell(_cellDataDistance,k,_dataMatrixDistance); 
		}
		//else matlab::mxDestroyArray(_cellDataDistance);
			
		if (ui.checkBoxExtractScanline->isChecked()) 
		{
			mxSetCell(_cellDataDirection,k,_dataMatrixDirection); 
		}

		
		//else matlab::mxDestroyArray(_cellDataDirection);
		// now update the lastOffset to the current Index / note special case the last position

		if ( _seekIndex+1 < _offsetVector.size() )
			_lastOffset = _offsetVector[_seekIndex+1];
	}
	//_reconstructionFile.close();
	CloseHandle(_reconstructionFile);
	_outStream.close();


	
	if (ui.checkBoxExtractIntensity->isChecked())
	{
		matlab::MATFile *pmatIntensity;
		pmatIntensity = matlab::matOpen("extractionIntensity.mat", "w");
		int status = matPutVariable(pmatIntensity, "Intensity", _cellDataIntensity);
		if ( status != 0 )
			std::cout << "Failed writing intensity data!" << std::endl;
		matlab::matClose(pmatIntensity);
		matlab::mxDestroyArray(_cellDataIntensity);
	}

	if (ui.checkBoxExtractDistance->isChecked()) 
	{
		matlab::MATFile *pmatDistance;
		pmatDistance = matlab::matOpen("extractionDistance.mat", "w");
		int status = matPutVariable(pmatDistance, "Distance", _cellDataDistance);
		if ( status != 0 )
			std::cout << "Failed writing distance data!" << std::endl;
		matlab::matClose(pmatDistance);
		matlab::mxDestroyArray(_cellDataDistance);
	}

	if ( ui.checkBoxExtractDepth->isChecked())
	{
		matlab::MATFile *pmatScanline;
		pmatScanline = matlab::matOpen("extractionScanline.mat", "w");
		int status = matPutVariable(pmatScanline, "Scanline", _cellDataScanline);
		if ( status != 0 )
			std::cout << "Failed writing scanline/depth data!" << std::endl;
		matlab::matClose(pmatScanline);
		matlab::mxDestroyArray(_cellDataScanline);

	}

	if (ui.checkBoxExtractScanline->isChecked())
	{
		matlab::MATFile *pmatDirection;
		pmatDirection = matlab::matOpen("extractionDirection.mat", "w");
		int status = matPutVariable(pmatDirection, "Direction", _cellDataDirection);
		if ( status != 0 )
			std::cout << "Failed writing direction data!" << std::endl;
		matlab::matClose(pmatDirection);
		matlab::mxDestroyArray(_cellDataDirection);
	}
	
	
	
	std::cout << _lineCounter << " lines written to file, " << _indicesVector.size() << " voxels" << std::endl;
	// progress bar
	m_percentage = -1;
	updatePercentage();
}

void Compounding::loadFileForPostProcessingRF()
{
	XMLTools *_xmlTools = new XMLTools();
	int _elements;
	_xmlTools->loadFile(ui.postProcessingTargetFilenameEdit->text(). toLatin1().data());

	std::string _reconstructionFilename, _indexFilename, _directionFilename;
	std::string _voxelsX, _voxelsY, _voxelsZ;
	_xmlTools->handleData(false, _voxelsX, "DataSize","VoxelsX");
	_xmlTools->handleData(false, _voxelsY, "DataSize","VoxelsY");
	_xmlTools->handleData(false, _voxelsZ, "DataSize","VoxelsZ");
	int _voxels[] = {atoi(_voxelsX.c_str()), atoi(_voxelsY.c_str()),atoi(_voxelsZ.c_str())};

	// clear the combo boxes from previous values
	ui.ComboFromX->clear();
	ui.ComboToX->clear();

	ui.ComboFromY->clear();
	ui.ComboToY->clear();

	ui.ComboFromZ->clear();
	ui.ComboToZ->clear();



	// fill in the voxel entries into the combo box
	for(int x=1; x<=atoi(_voxelsX.c_str());x++)
	{
		std::stringstream _ss;
		_ss << x;
		ui.ComboFromX->addItem(_ss.str().c_str());
		ui.ComboToX->addItem(_ss.str().c_str());
	}

	for(int y=1; y<=atoi(_voxelsY.c_str());y++)
	{
		std::stringstream _ss;
		_ss << y;
		ui.ComboFromY->addItem(_ss.str().c_str());
		ui.ComboToY->addItem(_ss.str().c_str());
	}

	for(int z=1; z<=atoi(_voxelsZ.c_str());z++)
	{
		std::stringstream _ss;
		_ss << z;
		ui.ComboFromZ->addItem(_ss.str().c_str());
		ui.ComboToZ->addItem(_ss.str().c_str());
	}

	// now set the TO combo box to the last element
	ui.ComboToX->setCurrentIndex(ui.ComboToX->count()-1);
	ui.ComboToY->setCurrentIndex(ui.ComboToY->count()-1);
	ui.ComboToZ->setCurrentIndex(ui.ComboToZ->count()-1);

	delete _xmlTools;
}

void Compounding::preProcessRF()
{


	if ( ui.filterList->currentIndex() == 0 ) // Butterworth & Envelope detection
	{
		std::string _temp = ui.postProcessingTargetFilenameEdit->text(). toLatin1().data();
		std::string _xmlFilename = _temp.substr( 0, _temp.find_last_of( ".xml" )-3 ).c_str();
		RFFileStream *_rfs = new RFFileStream((STREAMING_MODE)ui.comboBoxRFStreaming->currentIndex ());


		_rfs->openFile(ui.postProcessingTargetFilenameEdit->text(). toLatin1().data(), DATA_TYPE::USHORT16);

		std::string filename = ui.postProcessingTargetFilenameEdit->text(). toLatin1().data();
		std::string _path = filename.substr( 0, filename.find_last_of( '/' ) +1 );

		_temp = _rfs->getDataFileName();
		std::string _filename = _temp.substr( 0, _temp.find_last_of( ".3drf" )-3 ).c_str();

		// create new XML file containing information about the filtered data
		XMLTools *_xmlTools = new XMLTools();
		_xmlTools->createFile("Reconstruction");

		_xmlTools->handleData(true, _filename+".Envelope.3drf", std::string(IMAGEDATAFILENAME_TAG));
		_xmlTools->handleData(true, _rfs->getTimeStampFileName(), std::string(TIMESTAMPSFILENAME_TAG));
		{
			std::stringstream _ss, _ss2;
			_ss << "Frames: " <<_rfs->getRecordingInfo().frames;
			_ss2 << _rfs->getRecordingInfo().frames;
			_xmlTools->handleData(true, _ss2.str(), std::string(NUMBER_OF_FRAMES_TAG));
			ui.textEdit->append(QString(_ss.str().c_str()));
		}
		{
			std::stringstream _ss, _ss2;
			_ss << "Samples per Scanline: " <<_rfs->getRecordingInfo().height;
			_ss2 << _rfs->getRecordingInfo().height;
			_xmlTools->handleData(true,  _ss2.str(), std::string(SAMPLES_SCANLINE_TAG));
			ui.textEdit->append(QString(_ss.str().c_str()));
		}
		{
			std::stringstream _ss, _ss2;
			_ss << "Scanlines: " << _rfs->getRecordingInfo().width;
			_ss2 << _rfs->getRecordingInfo().width;
			_xmlTools->handleData(true, _ss2.str(), std::string(NUMBER_OF_SCANLINES_TAG));
			ui.textEdit->append(QString(_ss.str().c_str()));
		}
		{
			std::stringstream _ss, _ss2;
			_ss << "Sampling Frequency: " << _rfs->getRecordingInfo().samplingFrequency;
			_ss2 << _rfs->getRecordingInfo().samplingFrequency;
			_xmlTools->handleData(true, _ss2.str(), std::string(SAMPLING_FREQUENCY_TAG));
			ui.textEdit->append(QString(_ss.str().c_str()));
		}
		{
			std::stringstream _ss, _ss2;
			_ss << "Penetration Depth: " << _rfs->getRecordingInfo().penetrationDepth;
			_ss2 << _rfs->getRecordingInfo().penetrationDepth;
			_xmlTools->handleData(true, _ss2.str(), std::string(PENETRATION_DEPTH_TAG));
			ui.textEdit->append(QString(_ss.str().c_str()));
		}
		{
			std::stringstream _ss, _ss2;
			_ss << "Recording Frequency: " <<_rfs->getRecordingInfo().recordingFrequency;
			_ss2 << _rfs->getRecordingInfo().recordingFrequency;
			_xmlTools->handleData(true,  _ss2.str(), std::string(RECORDING_FREQUENCY_TAG));
			ui.textEdit->append(QString(_ss.str().c_str()));
		}
		{
			std::stringstream  _ss2;
		
			_ss2 << _rfs->getRecordingInfo().probeName;
			_xmlTools->handleData(true,  _ss2.str(), std::string(PROBE_NAME_TAG));
			
		}




		_xmlTools->handleData(true, std::string(FILTER_BANDPASS), std::string(FILTER_OPERATION_TAG), std::string(FILTER_TYPE_ATTRIBUTE));
		_xmlTools->handleData(true, std::string(FILTER_ENVELOPE), std::string(FILTER_OPERATION_TAG), std::string(FILTER_TYPE_ATTRIBUTE));

		_xmlTools->saveToFile((_xmlFilename+".Envelope.xml"));
		_xmlTools->closeHandle();
		delete _xmlTools;

		std::fstream _filteredFile;

		_filteredFile.open(_path+_filename+".Envelope.3drf", std::ios::out  | std::ios::binary);

		// progress bar
		m_percentage = -1;
		updatePercentage();
		int _1percent = static_cast<int>(_rfs->getRecordingInfo().frames)/100.0;
		int _percentCounter = 0;

		// data conversion
		RFProcessing<float> *_rfp = new RFProcessing<float>(_rfs->getRecordingInfo().width, _rfs->getRecordingInfo().height, atof(ui.butterworthLow->text(). toLatin1().data()), atof(ui.butterworthHigh->text(). toLatin1().data()), atof(ui.butterworthSamplingFrequency->text(). toLatin1().data()), atoi(ui.butterworthOrder->text(). toLatin1().data()));
		for(int i=0;i<_rfs->getRecordingInfo().frames;i++)
		{
			_percentCounter++;
			if ( _percentCounter == _1percent)
			{
				_percentCounter = 0;
				updatePercentage();
			}
			int16 *_data = reinterpret_cast<int16*>(_rfs->getFrameByIndex(i));
			// have to copy the ENTIRE RF frame into memory, no FILE pointer
			short * _memPtr = new short[_rfs->getRecordingInfo().width * _rfs->getRecordingInfo().height];

			memcpy(_memPtr,&_data->byte.low,_rfs->getFrameSize());
			/*short *_mP = &_memPtr[0];
			int16 *_mP2 = &_data[0];
			for(int j=0;j<_rfs->getRecordingInfo().width * _rfs->getRecordingInfo().height;j++)
			{
			*_mP = (*_mP2).value;
			_mP2++;
			_mP++;
			}*/

			// data conersion
			RFData<float> *_rfRawData = _rfp->convert(&_memPtr[0], _rfs->getFrameSize()); 
			RFData<float> *_rfButterworthFilteredData = _rfp->bandpassFilter(_rfRawData);
			float *_rfEnvelopeFilteredData = _rfp->envelopeDetection(_rfButterworthFilteredData->getData());

			// Image Preview - BEGIN

			{
				//IplImage* _tmpImage = cvCreateImage( cvSize(_rfs->getRecordingInfo().width, _rfs->getRecordingInfo().height), 8, 3 );



				float *_envPtr = &_rfEnvelopeFilteredData[0];
				//char *_chPtr = &_tmpImage->imageData[0];

				char *_dataTMP = new char [_rfs->getRecordingInfo().width * _rfs->getRecordingInfo().height];
				char *_dataPtr = &_dataTMP[0];
				int _maxVal = 0;

				/*std::fstream _rfDataFile;
				std::stringstream _ss;
				_ss << "gray" << i << ".rf";
				_rfDataFile.open(_ss.str(),std::ios::out   | std::ios::binary);*/

				for(int y=0;y<_rfs->getRecordingInfo().width * _rfs->getRecordingInfo().height;y++)
				{
					char _val = std::min(20*(int)log(*_envPtr),255);

					*_dataPtr++ = _val;

					//_rfDataFile.write(&_val, sizeof(char));
					//if (*_chPtr > _maxVal )
					//	_maxVal = *_chPtr;

					/**_chPtr = _val; // R
					_chPtr++;

					*_chPtr = _val; // G
					_chPtr++;

					*_chPtr = _val; // B
					_chPtr++;*/

					_envPtr++;

				}
				//_rfDataFile.close();



				VideoTools _videoTools;

				//HBITMAP _bmp = m_videoTools->IplImage2HBITMAP(_tmpImage);
				HBITMAP _bmp = m_videoTools->uchar2HBITMAP(_rfs->getRecordingInfo().height, _rfs->getRecordingInfo().width, (unsigned char*)_dataTMP);
				delete[] _dataTMP;

				/*{
				std::stringstream _ss;
				_ss << "color" << i << ".bmp";
				HDC hnewDC;
				hnewDC = ::CreateCompatibleDC( NULL );
				GraphicTools _gt;
				PBITMAPINFO pBitmapInfo = _gt.CreateBitmapInfoStruct(_bmp);

				_gt.CreateBMPFile(_ss.str(), pBitmapInfo,_bmp,hnewDC);

				::DeleteDC(hnewDC);
				}*/
				//std::stringstream _ss;
				//_ss << "color" << i << ".jpg";
				//cvSaveImage(_ss.str().c_str(),_tmpImage);
				//cvReleaseImage(&_colImage);


				if ( m_rfPreviewImage != 0)
					DeleteObject(m_rfPreviewImage);

				m_rfPreviewImage = _bmp;


				//QPixmap::HBitmapFormat _format = QPixmap::NoAlpha;
				QPixmap _pixmap = QtWin::fromHBITMAP(_bmp, QtWin::HBitmapNoAlpha); 

				float _targetSize = 600.0;
				//float _scaleFactor = _targetSize / (float) std::max(_width,_height);
				ui.RFImagePreviewLabel->setPixmap(_pixmap.scaled(_targetSize, _targetSize));

				//cvReleaseImage(&_tmpImage);

				//delete[] _tmpImg;


			}

			// Image Preview - END


			/*
			std::fstream _rfDataFile;
			std::stringstream _ss;
			_ss << "envelope" << i << ".rf";
			_rfDataFile.open(_ss.str(),std::ios::out   | std::ios::binary);
			*/
			float *_envPtr = &_rfEnvelopeFilteredData[0];
			//float *_envPtr = &_rfRawData->getData()[0];
			//float *_envPtr = &_rfButterworthFilteredData->getData()[0];
			//short* _envPtr = &_memPtr[0];
			// now write the data to file
			for (int j=0;j<_rfs->getRecordingInfo().width*_rfs->getRecordingInfo().height;j++)
			{
				// data conversion
				//float32 _float32;
				//_float32.value = fabs(*_envPtr);

#ifdef USE_FLOAT32
				float32 _float32;
				_float32.value = fabs(*_envPtr);
				_filteredFile.write(&_float32.byte.c[0],sizeof(float32));
#endif

#ifdef USE_USHORT16
				ushort16 _ushort16;
				_ushort16.value = roundf(fabs(*_envPtr));
				_filteredFile.write(&_ushort16.byte.c[0],sizeof(ushort16));
#endif
				_envPtr++;
				/*short16 _short16;
				_short16.value = *_envPtr;
				_filteredFile.write(&_short16.byte.c[0],sizeof(short16));
				_envPtr++;*/
				// DEBUG
				//_rfDataFile.write(&_float32.byte.c[0],sizeof(float32));




			}
			//_rfDataFile.close();


			delete _rfRawData;
			delete _rfButterworthFilteredData; 
			delete[] _memPtr;
			delete[] _rfEnvelopeFilteredData;
		}

		delete _rfp;

		_filteredFile.close();
		//float *data = reinterpret_cast<float*>(_rfs->getFrameByIndex(0));


		// progress bar
		m_percentage = -1;
		updatePercentage();

		delete _rfs;
	}
	else
	{
		XMLTools *_xmlTools = new XMLTools();
		int _elements;


		std::string _temp = ui.postProcessingTargetFilenameEdit->text(). toLatin1().data();
		std::string _xmlFilename = _temp.substr( 0, _temp.find_last_of( ".xml" )-3 ).c_str();
		_xmlTools->loadFile(ui.postProcessingTargetFilenameEdit->text(). toLatin1().data());
		std::string _filename = (ui.postProcessingTargetFilenameEdit->text(). toLatin1().data());
		std::string _path = _filename.substr( 0, _filename.find_last_of( '/' ) +1 );

		std::string _reconstructionFilename, _indexFilename, _directionFilename;
		std::string _voxelsX, _voxelsY, _voxelsZ;
		_xmlTools->handleData(false, _voxelsX, "DataSize","VoxelsX");
		_xmlTools->handleData(false, _voxelsY, "DataSize","VoxelsY");
		_xmlTools->handleData(false, _voxelsZ, "DataSize","VoxelsZ");
		int _voxels[] = {atoi(_voxelsX.c_str()), atoi(_voxelsY.c_str()),atoi(_voxelsZ.c_str())};

		_elements = atoi(_voxelsX.c_str()) * atoi(_voxelsY.c_str()) * atoi(_voxelsZ.c_str());

		_xmlTools->handleData(false,_reconstructionFilename, RF_RECONSTRUCTIONFILE_TAG);
		_xmlTools->handleData(false,_indexFilename, RF_INDEXFILE_TAG);
		_xmlTools->handleData(false,_directionFilename, RF_DIRECTIONFILE_TAG);

		delete _xmlTools;
	}
	
}

void Compounding::rfFrameChanged(int)
{
	//ui.RFImageLabelScrollBar

	int _rfFrames = 0;
	int _index = 0;
	int *_index2 = new int[m_rfFileStreamVector.size()];
	_index2[0] = 0;
	for (int i=0; i < m_rfFileStreamVector.size(); i++)
	{
		if ( i > 0 )
			_index2[i]=_index2[i-1]+ m_rfFileStreamVector.at(i-1)->getRecordingInfo().frames;

		_index = i;
		_rfFrames += m_rfFileStreamVector.at(i)->getRecordingInfo().frames;
		if ( ui.RFImageLabelScrollBar->value() <= _rfFrames )
			break;
	}



	std::stringstream _ss;

	_ss << ui.RFImageLabelScrollBar->value();
	ui.RFMasterFrameEdit->setText(QString(_ss.str().c_str()));
	ui.RFMasterFrameEdit2->setText(QString(_ss.str().c_str()));

	int _frame = /*_rfFrames - */ui.RFImageLabelScrollBar->value();

	//std::cout << "Frame: " << _frame << std::endl;

#ifdef USE_USHORT16
	ushort16 *_rfPtr = reinterpret_cast<ushort16*>(m_rfFileStreamVector[_index]->getFrameByIndex(_frame - _index2[_index]));

	// have to copy the ENTIRE RF frame into memory, no FILE pointer
	unsigned short * _memPtr = new unsigned short[m_ultrasoundSettings[_index].rfFrameSizeElements];
	memcpy(_memPtr,&_rfPtr->byte.c[0],m_ultrasoundSettings[_index].rfFrameSizeBytes);

	//std::cout << "Bytes: " << m_ultrasoundSettings[_index].rfFrameSizeBytes << std::endl;
	//std::cout << "Elements: " << m_ultrasoundSettings[_index].rfFrameSizeElements << std::endl;


	//IplImage* _tmpImage = cvCreateImage( cvSize(m_ultrasoundSettings[_index].samplesAlongScanline, m_ultrasoundSettings[_index].scanlines), 8, 3 );

	//unsigned char *_tmpImg = new unsigned char[m_ultrasoundSettings[_index].rfFrameSizeElements*3];
	unsigned short *_envPtr = /*&_rfRawData->getData()[0];*/&_memPtr[0];
	
	//char *_chPtr = &_tmpImage->imageData[0];

	char *_tmpImg = new char [m_ultrasoundSettings[_index].rfFrameSizeElements];
	char *_chPtr = &_tmpImg[0];

#endif

#ifdef USE_FLOAT32

	float32 *_rfPtr = reinterpret_cast<float32*>(m_rfFileStreamVector[_index]->getFrameByIndex(_frame - _index2[_index]));

	// have to copy the ENTIRE RF frame into memory, no FILE pointer
	float * _memPtr = new float[m_ultrasoundSettings[_index].rfFrameSizeElements];
	memcpy(_memPtr,&_rfPtr->byte.c[0],m_ultrasoundSettings[_index].rfFrameSizeBytes);

	//std::cout << "Bytes: " << m_ultrasoundSettings[_index].rfFrameSizeBytes << std::endl;
	//std::cout << "Elements: " << m_ultrasoundSettings[_index].rfFrameSizeElements << std::endl;


	//IplImage* _tmpImage = cvCreateImage( cvSize(m_ultrasoundSettings[_index].samplesAlongScanline, m_ultrasoundSettings[_index].scanlines), 8, 3 );

	//unsigned char *_tmpImg = new unsigned char[m_ultrasoundSettings[_index].rfFrameSizeElements*3];
	float *_envPtr = /*&_rfRawData->getData()[0];*/&_memPtr[0];
	
	//char *_chPtr = &_tmpImage->imageData[0];

	char *_tmpImg = new char [m_ultrasoundSettings[_index].rfFrameSizeElements];
	char *_chPtr = &_tmpImg[0];
#endif
	int _maxVal = 0;


	for(int y=0;y<m_ultrasoundSettings[_index].rfFrameSizeElements;y++)
	{
		char _val = std::min(20*(int)log((float)(*_envPtr)),255);

		//if (*_chPtr > _maxVal )
		//	_maxVal = *_chPtr;

		*_chPtr = _val; // R
		_chPtr++;

		/**_chPtr = _val; // G
		_chPtr++;

		*_chPtr = _val; // B
		_chPtr++;*/

		_envPtr++;

	}


	VideoTools _videoTools;
	//IplImage *_tmpImage=  _videoTools.uchar2IplImage(_tmpImg, m_ultrasoundSettings[_index].samplesAlongScanline, m_ultrasoundSettings[_index].scanlines ,3,8);

	//QImage m_image = QImage((uchar *) _tmpImg, _rfs->getRecordingInfo().width, _rfs->getRecordingInfo().height, QImage::Format_Mono);

	//ui.RFImageLabel->setPixmap(QPixmap::fromImage(m_image));

	//IplImage* _colImage = cvCreateImage(cvSize( m_ultrasoundSettings[_index].samplesAlongScanline, m_ultrasoundSettings[_index].scanlines), IPL_DEPTH_8U, 3);
	//cvCvtColor(_tmpImage, _colImage, CV_GRAY2RGB);
	
	//HBITMAP _bmp = m_videoTools->IplImage2HBITMAP(_tmpImage);
	HBITMAP _bmp = m_videoTools->uchar2HBITMAP(m_ultrasoundSettings[_index].samplesAlongScanline, m_ultrasoundSettings[_index].scanlines, (unsigned char*)_tmpImg);
	delete[] _tmpImg;
	//cvSaveImage("color.jpg",_tmpImage);
	//cvReleaseImage(&_colImage);

	/*{
				std::stringstream _ss;
				_ss << "RFimage" << _frame << ".bmp";
				HDC hnewDC;
				hnewDC = ::CreateCompatibleDC( NULL );
				GraphicTools _gt;
				PBITMAPINFO pBitmapInfo = _gt.CreateBitmapInfoStruct(_bmp);

				_gt.CreateBMPFile(_ss.str(), pBitmapInfo,_bmp,hnewDC);

				::DeleteDC(hnewDC);
				}*/

	if ( m_rfImage != 0)
		DeleteObject(m_rfImage);

	m_rfImage = _bmp;


	//QPixmap::HBitmapFormat _format = QPixmap::NoAlpha;
	QPixmap _pixmap = QtWin::fromHBITMAP(_bmp, QtWin::HBitmapNoAlpha); 

	float _targetSize = 600.0;
	//float _scaleFactor = _targetSize / (float) std::max(_width,_height);
	ui.RFImageLabel->setPixmap(_pixmap.scaled(_targetSize, _targetSize));

	//cvReleaseImage(&_tmpImage);

	//delete[] _tmpImg;

	delete[] _index2;
	delete[] _rfPtr;
}

void Compounding::voxelSubsampling(int state)
{

	ui.VoxelSubsamplingComboBox->setEnabled(ui.VoxelSubsamplingCheckBox->checkState() == Qt::Checked);

}

void Compounding::withinVoxelSubsampling(int state)
{
	ui.WithinVoxelSubsamplingComboBox->setEnabled(ui.WithinVoxelSubsamplingCheckbox->checkState() == Qt::Checked);
}

void Compounding::RFLogCompression(int state)
{
	m_logCompress = (ui.checkBoxLogCompression->checkState() == Qt::Checked);

}