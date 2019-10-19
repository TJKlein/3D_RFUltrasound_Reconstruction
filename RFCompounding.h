#ifndef RFCOMPOUNDING_H
#define RFCOMPOUNDING_H

#include <Windows.h>


/**
 * \brief 3D Reconstruction software of ultrasound data using RF-data
 *
 * \authors Tassilo Klein
 *          <br>
 *          kleint\@cs.tum.edu
 * \ingroup Ultrasound
 * \version 1.0a
 * \date 07-06-2010
 *
 * \par License:
 * Copyright (c) 2005, 2006, 2007, 2008, 2009
 * Chair for Computer Aided Medical Procedures & Augmented Reality (CAMPAR) I-16,
 * Technische Universit&auml;t M&uuml;nchen
 * <br>
 * <br>
 * All rights reserved.
 * <br>
 * <br>
 * See <a href="COPYRIGHT.txt">COPYRIGHT.txt</a> for details.
 * <br>
 * <br>
 * This software is distributed WITHOUT ANY WARRANTY; without even 
 * <br>
 * the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
 * <br>
 * PURPOSE.  See the <a href="COPYRIGHT.txt">COPYRIGHT.txt</a> notice
 * for more information.
 *
 */

//#define USE_FLOAT32 0
#define USE_USHORT16 1
//#define USE_FLOAT5X32 0

#define QUANTIZE_DISTANCE 1

enum STREAMING_MODE { NORMAL, MISSING4BYTE, UNNECESSARY4BYTEHEADER};

struct FrameScanlineRange
{
	short start;
	short end;
};

union ulong64
{
	struct {
		char c[8];
	} byte;
	unsigned __int64 value;
};


union ushort16
{
	struct
	{
		char c[2];
	} byte;

	unsigned short value;
};

union short16
{
	struct
	{
		char c[2];
	} byte;

	short value;
};

union ulong32
{
	struct
	{
		char c[4];
	} byte;

	unsigned long value;
};

union uint32
{
	struct
	{
		char c[4];
	}byte;
	unsigned int value;
};

union float32
{
	struct
	{
		char c[4];
	}byte;
	float value;
};


union int16
{
	struct
	{
		 char low;
		 char high;
	} byte;
	short value;
} ;

// define the tags for the XML configration file
const char VERSION_TAG[] = {"Version"};
const char IMAGEDATAFILENAME_TAG[] = {"ImageDataFile"};
const char TIMESTAMPSFILENAME_TAG[] = {"TimeStampsFile"};
const char NUMBER_OF_FRAMES_TAG[] = {"Frames"};
const char SAMPLES_SCANLINE_TAG[] = {"ScanlineSamples"};
const char NUMBER_OF_SCANLINES_TAG[] = {"Scanlines"};
const char SAMPLING_FREQUENCY_TAG[] = {"SamplingFrequency"};
const char PENETRATION_DEPTH_TAG[] = {"PenetrationDepth"};
const char RECORDING_FREQUENCY_TAG[] = {"RecordingFrequency"};
const char FILTER_OPERATION_TAG[] = {"AppliedFilter"};
const char FILTER_TYPE_ATTRIBUTE[] = {"type"};
const char FILTER_ENVELOPE[] = {"Envelope"};
const char FILTER_BANDPASS[] = {"Butterworth"};
const char PROBE_NAME_TAG[] = {"ProbeName"};

enum DATA_TYPE {USHORT16, FLOAT32, FLOAT5X32};



#include <Windows.h>
#include <Eigen/StdVector>
#include <Eigen/Core>
#include <iostream>
#include <fstream>
#include "GeometricTools.h"
#include "NumericTools.h"
#include "iir.h"
#include "XMLTools.h"

#define roundf(x) (x<0?ceil((x)-0.5):floor((x)+0.5))

const int compact_range = 1024;
const int MAX_VAL = 5;

static short compactFloat(float input) {
    return roundf(input * compact_range / MAX_VAL);
}
static double expandToFloat(short input) {
    return ((double)input) * MAX_VAL / compact_range;
}


/*union ScanlineIndex
{
	struct { // size and order matters
		unsigned short scanlineID; // least significant bytes
		unsigned short imageID; // most significant (to later one, file ID will be MSB, with lower part the image ID)
	} long32;
	unsigned long value;
};*/

struct ScanlineIndex
{
	unsigned short fileID; // when merging RF files, this indicates the index (integer) to the file
	unsigned short scanlineID; // least significant bytes
	unsigned short imageID; // most significant (to later one, file ID will be MSB, with lower part the image ID)
};

struct UltrasoundSettings
{
	// RF specific
	float probeAngle;
	int probeRadius; // mm
	float probeWidth; // mm (if it is a linear transducer we don't have an angle!!)

	unsigned short samplesAlongScanline;
	unsigned short penetrationDepth;
	unsigned short scanlines;

	float scaleX;
	float scaleY;
	float apexX;
	float apexY;
	float angle;
	float mmPerSample;
	float samplesPerMM;
	Eigen::Matrix4f calibrationMatrix;
	int rfFrameSizeElements;
	int rfFrameSizeBytes;
};

template< typename T >
class Scanline
{
private:

	ScanlineIndex m_scanlineIndex;

	unsigned int m_uniqueIndex;

	std::vector<UltrasoundSettings> *m_ultrasoundSettings;
	unsigned char m_settingsIndex;


	/**
	\brief Specifies the line going from the apex to the end of beam
	**/
	Line *m_line;

	T *m_data;

public:


	void setData(T* data) { m_data = data;}
	T* getData() { return m_data;}

	Line* getLine() { return m_line; }

	/**
	\brief Constructor of the scanline class, that contains information about the spatial location of the scanline as well as its indices, to find the corresponding image data
	\param apexPosition Position of the ultrasound apex of the RF image in world space coordinates
	\param beamEnd End position of the sanline in world space coordinates
	\param numberOfSamples Number of RF samples along scanline
	\param length Length in mm of the scanline (penetration depth)
	\param uniqueIndex Each beam direction will be given an unique index, which will be saved to the reconstruction file later on
	\param fileIndex Index indicating the file where data is located (when merging multiple RF data sets)
	\param rfImageIndex Index of the RF image (in the data set)
	\param beamIndex Index of the scanline on the RF image
	**/
	Scanline(Eigen::Vector3f &apexPosition, Eigen::Vector3f &beamEnd, std::vector<UltrasoundSettings> *ultrasoundSettings, unsigned char settingsIndex, unsigned int uniqueIndex, unsigned char fileIndex, unsigned short rfImageIndex, unsigned short beamIndex = 0);


	/**
	\brief Computes the distance from a point in world coordinate space to the scaline
	\param point A point in world coordinate space
	\param t Given the line paramater t from: L(t)=apex + t * (beamEnd - apex)
	**/
	float distancePointToScanline(const Eigen::Vector3f &point, float &t);

	/**
	\brief Computes the distance from a point in world coordinate space to the scaline maximum norm
	\param point A point in world coordinate space
	\param t Given the line paramater t from: L(t)=apex + t * (beamEnd - apex)
	**/
	float distancePointToScanlineMAX(const Eigen::Vector3f &point, float &t);
	
	int getFileIndex() { return static_cast<int>(m_scanlineIndex.fileID); }
	/**
	\brief Access to the beam position of in the RF image
	\return Scanline position in the corresponding RF image
	**/
	unsigned short getBeamIndex() { return m_scanlineIndex.scanlineID; }

	UltrasoundSettings* getUltrasoundSettings() { return &((*m_ultrasoundSettings)[static_cast<int>(m_settingsIndex)]); }

	/**
	\brief Access to the RF image index
	\return RF image index
	**/
	unsigned short getRfImageIndex() { return m_scanlineIndex.imageID; }


	//unsigned long getIndex() { return m_scanlineIndex.value; }

	unsigned int getUniqueIndex() { return m_uniqueIndex; }

	unsigned short getNumberOfSamples() { return (*m_ultrasoundSettings)[static_cast<int>(m_settingsIndex)].samplesAlongScanline; }

	/**
	\return Penetration depth of the ultrasound beam (mm)
	**/
	unsigned short getLength() { return  (*m_ultrasoundSettings)[static_cast<int>(m_settingsIndex)].penetrationDepth; }

	~Scanline()
	{
		delete m_line;
	}
};

#ifdef USE_USHORT16
template class Scanline<unsigned short>;
#endif

#ifdef USE_FLOAT32 || USEFLOAT5X32
template class Scanline<float>;
#endif

template< typename T >
class ScanlineComparator
{
	bool reverse;

public:
	ScanlineComparator(const bool& revparam=false)	{ reverse=revparam; }
	bool operator() (Scanline<T>* lhs, Scanline<T>* rhs) const
	{
		unsigned long _ldx =  lhs->getUniqueIndex();// lhs->getIndex();
		unsigned long _rdx =  rhs->getUniqueIndex();//rhs->getIndex();
		if (reverse) 
			return (_ldx < _rdx );
		else 
			return ( _ldx > _rdx );
			
	}
};

template< typename T >
class RFData
{
protected:
	/**
	\brief Pointer to the RF data
	**/
	T* m_data;
	/**
	\brief Number of beams
	**/
	int m_width;
	/**
	\brief Number of samples along beam
	**/
	int m_height;
public:
	RFData();
	~RFData() { if (m_data) 
				delete[] m_data; } // !! make sure data was created with new operator
	void setData(T *data) { m_data = data; }
	void setData(T *data, int width, int height) { setData(data); m_width = width; m_height = height; }
	T* getData() { return m_data; }
	T* getScanline(int index) { return &m_data[index*m_height]; }

	void setNumberOfSamples(int size) { m_height = size; }
	void setNumberOfScanlines(int size) { m_width = size; }
	int numberOfSamples() { return m_height; }
	int numberOfScanlines() { return m_width; }

};

template< typename T >
class StandardRF : public RFData<T>
{
private:
	/**
	\brief This transformation describes the maping from the image plane (mm) to global coordinates, that is tracking_transform * calibration matrix
	**/
	Eigen::Matrix4f m_transformation;

	std::vector<Scanline<T>> m_scanlines;

public:
	StandardRF();

	void setTransformation(Eigen::Matrix4f &transform) { m_transformation = transform; }
	/**
	\brief Given the RF image, transformation matrix and the position of the fan during calibration
	\param scanlineVector Reference to the vector, where the individual scanlines are pushed in. As beam index serves the position in the vector.
	\param scanlineIndex Every scanline has an unique number identification which will be stored in the 3d reconstruction file later on
	\param imageIndex Index of RF image
	\param imageWidth Width of the image on which calibration was performed
	\param imageHeight Height of the image on which calibration was performed
	**/
	void computeScanlines(std::fstream &file, std::vector<Scanline<T>*> &scanlineVector, std::vector<FrameScanlineRange> &fsr, unsigned char fileIndex, unsigned int &scanlineIndex, int imageIndex, std::vector<UltrasoundSettings> *ultrasoundSettings, unsigned char settingsIndex/*, unsigned short samplesPerScanline, unsigned short penetrationDepth, int scanlines, float apexX, float apexY, float angle*/);
};

#ifdef USE_USHORT16
template class StandardRF<unsigned short>;
#endif

#ifdef USE_FLOAT32 || USEFLOAT5X32
template class StandardRF<float>;
#endif
/**
\brief This class could be used for incorporating pre-processed RF data into the reconstruction process
**/
template< typename T >
class CompositeRF : public StandardRF<T>
{
private:
	float* m_compositeData;
public:
	CompositeRF();
	float* getCompositeData() { return m_compositeData; }
};

template< typename T >
class ColorRF : public RFData<T>
{
public:
	ColorRF();
};

struct RF_RECORDING_INFO
{
	int frames;
	int height;
	int width;
	int samplingFrequency;
	int recordingFrequency;
	int penetrationDepth;
	std::string probeName;
};

class RFFileStream
{
private:
	std::vector<float> m_timeStampVector;
	//std::ifstream m_rfDataFile;
	HANDLE m_rfDataFile;
	DWORD m_rfDataFileLength;
	//int m_rfDataFileLength;
	XMLTools m_xmlTools;
	std::string m_dataFileName;
	std::string m_timeStampFileName;
	std::string m_probeName;
	RF_RECORDING_INFO m_recordingInfo;
	float m_version;
	/*int m_frames;
	int m_width;
	int m_height;
	int m_samplingFrequency;
	int m_penetrationDepth;
	int m_recordingFrequency;*/
	int m_rfFrameSize;
	bool m_evenlopeSignal;
	int m_currentFrame;
	int m_rfFilePointer;
	STREAMING_MODE m_mode;

	//std::ofstream m_file;

public:
	~RFFileStream();
	RFFileStream(STREAMING_MODE mode = STREAMING_MODE::NORMAL);
	/**
	\brief Opens a xml file describing the data set
	\param filename Name of the .xml file describing the recorded data
	**/
	bool openFile(std::string filename, DATA_TYPE type);

	/**
	\brief Write the file from the buffer to the disk
	**/
	bool writeFile();

	/**
	\brief Close the file handle
	**/
	void closeFile();

	/**
	\brief Returns the index of the current file frame position
	\return Index of the current frame
	**/
	int currentFrameIndex();
	
	/**
	\brief Returns pointer given the frame index (starting with 0).
	\return Either it is an array of float (filtered signal) or array of short (raw RF data)
	**/
	void* getFrameByIndex(int index);

	/**
	\brief Returns the next frame
	\return Either it is an array of float (filtered signal) or array of short (raw RF data)
	**/
	void* getNextFrame();

	/**
	\brief Returns the number of bytes of a RF frame
	**/
	int getFrameSize() { return m_rfFrameSize; }

	RF_RECORDING_INFO getRecordingInfo() { return m_recordingInfo; }

	bool isEnvelopeSignal() { return m_evenlopeSignal; }

	std::string getDataFileName() { return m_dataFileName; }
	std::string getTimeStampFileName() { return m_timeStampFileName; }
};

template< typename T >
class RFProcessing
{
private:
	int m_width;
	int m_height;
	double m_scalingFactor;        // scaling factor
    double *m_coefficientsD;     // d coefficients
    int *m_coefficientsC; // c coefficients
	int m_filterOrder;
public:
	/**
	\brief Constructor of the RF Processing
	\param width The number of beams / scanlines
	\param height The number of samples along a beam
	\param lowFreq Low frequency cut-off
	\param highFreq High frequency cut-off
	\param samplingFreq Sampling frequency of the signal
	\param filterOrder Order of the filter
	**/
	RFProcessing(int width, int height, float lowFreq, float highFreq, float samplingFreq, int filterOrder);

	~RFProcessing();



	/**
	\brief Perfoms a butterworth bandpass filter on a scanline per scanline basis (not in place)
	\param rf Pointer to the signal to be filtered
	\return New data array that has to be deleted when not needed anymore
	**/
	RFData<T>* bandpassFilter(RFData<T> *rf);


	/**
	\brief Performs an envelope detection by computing the absolute value of the hilbert transformed data set (not in place)
	\param data Pointer to the array holding the RF data values
	\return New data array that has to be deleted when not needed anymore
	**/
	float* envelopeDetection(float *data);

	/**
	\brief Converts a block of signed short (16-bits) into a array of float. The data is beam by beam. If beam length is N, data[N], will be the first element of the second beam
	\param data The block of beam data
	\param size The number of bytes. If we have K beams of length N with each element 2 bytes - size = K*N*2 (bytes)
	**/
	RFData<T>* convert(char *data, int size);


	/**
	\brief Converts a block of signed short (16-bits) into a array of float. The data is beam by beam. If beam length is N, data[N], will be the first element of the second beam
	\param data The block of beam data
	\param size The number of bytes. If we have K beams of length N with each element 2 bytes - size = K*N*2 (bytes)
	**/
	RFData<T>* convert(short *data, int size);

private:
		/**
	\brief Performs a butterwoth filter operation on a position of the beam
	\param output Filter output (results will be put here)
	\param input Filter input (unfiltered data)
	**/
	void filterBeam(double *output, double *input, int pos); 

};

template class RFProcessing<float>;


class Measurement
{
};


/**
\brief This class saves for each measurement all data such as value, distance and ScanlineIndex
**/
class StandardMeasurement : public Measurement
{
private:
	float m_data;
	float m_distance;
	ScanlineIndex m_directionIndex;
	
public:
	float getData() { return m_data; }
	float getDistance() { return m_distance; }
	ScanlineIndex getScanlineIndex() { return m_directionIndex; }

	void setDirectionIndex(ScanlineIndex directionIndex) { m_directionIndex = directionIndex; }
	void setData(float data) { m_data = data; }
	void setDistance(float distance) { m_distance = distance; }
};

#define SAVE_DISTANCE_SCANLINE 1

#ifdef USE_USHORT16
struct IntensityAndDistance
{
	unsigned short intensity;
	//short distance;
#ifdef SAVE_DISTANCE_SCANLINE
	#ifdef QUANTIZE_DISTANCE
	unsigned char distance;
#else
	float distance;
#endif
	unsigned short scanlineIndex;
#endif
};
#endif

#ifdef USE_FLOAT32
struct IntensityAndDistance
{
	float intensity;
	//short distance;
#ifdef SAVE_DISTANCE_SCANLINE
#ifdef QUANTIZE_DISTANCE
	unsigned char distance;
#else
	float distance;
#endif
	unsigned short scanlineIndex;
#endif
};
#endif

#ifdef USE_FLOAT5X32
struct IntensityAndDistance
{
	float mrf[5];
	//short distance;
#ifdef SAVE_DISTANCE_SCANLINE
	#ifdef QUANTIZE_DISTANCE
	unsigned char distance;
#else
	float distance;
#endif
	unsigned short scanlineIndex;
#endif
};
#endif



/**
\brief This class is especially suitable when different measurements are taken from the same direction. Direction is saved once,
whereas distance and value are stored in a vector
**/
template <typename T>
struct StandardMeasurementCollection
{
	//std::vector<T> data;
	T data[120]; // was 45 before
	short count; // was short before
	ScanlineIndex scanlineIndex;
	unsigned int uniqueID;
};
/*class StandardMeasurementCollection : public Measurement
{
	private:
	std::vector<T> m_data;
	ScanlineIndex m_directionIndex;
	
public:
	T& getData() { return m_data; }
	ScanlineIndex getDirectionIndex() { return m_directionIndex; }

	void setScanlineIndex(ScanlineIndex directionIndex) { m_directionIndex = directionIndex; }
	void push(T data) { m_data.push_back(data); }
	T pop() { return m_data.pop_back(); }

	int measurements() { return m_data.size(); }
	void setData(T data) { m_data = data; }
};*/


template <typename TX> struct DataAndHeap
{
	HANDLE heap;
	StandardMeasurementCollection<TX>* data;
};	

/**
\brief Class for handling access to a standard measurement collection
**/
template <typename TX>
class StandardMeasurementCollectionManager
{
	private:
		std::vector<StandardMeasurementCollection<TX>*>* m_measurementsCollection;
		int m_elementsPerSlice; // X-Y plane
		int m_elementsPerStride; // X vector
		int m_maxX;
		int m_maxY;
		int m_maxZ;
		std::vector<DataAndHeap<TX>> m_smcBlocks;

public:
	StandardMeasurementCollectionManager(int maxX, int maxY, int maxZ);
	~StandardMeasurementCollectionManager();

	void addBlockPointer(StandardMeasurementCollection<TX> *_smcBlock, HANDLE _heap);
	void setData(StandardMeasurementCollection<TX> *val, int posX, int posY, int posZ);
	void trim(int posX, int posY, int posZ);
	std::vector<StandardMeasurementCollection<TX>*>* getData(int posX, int posY, int posZ);
	std::vector<StandardMeasurementCollection<TX>*>* getData(int index);
	
};

template  <typename TX> void StandardMeasurementCollectionManager<TX>::addBlockPointer(StandardMeasurementCollection<TX> *_smcBlock, HANDLE _heap)
{
	DataAndHeap<TX> _data;
	_data.data = _smcBlock;
	_data.heap = _heap;
	
	m_smcBlocks.push_back(_data);
}

template  <typename TX> StandardMeasurementCollectionManager<TX>::StandardMeasurementCollectionManager(int maxX, int maxY, int maxZ) 
{
	m_measurementsCollection = new std::vector<StandardMeasurementCollection<TX>*>[maxX*maxY*maxZ];

	m_elementsPerSlice = maxX*maxY;
	m_elementsPerStride = maxX;

	m_maxX = maxX;
	m_maxY = maxY;
	m_maxZ = maxZ;

	m_smcBlocks.reserve(maxX*maxY*maxZ);

	//hMutex = CreateMutex(NULL, FALSE, NULL);
}


template  <typename TX> void StandardMeasurementCollectionManager<TX>::setData(StandardMeasurementCollection<TX> *val, int posX, int posY, int posZ)
{
	int _index = posZ*m_elementsPerSlice+posY*m_elementsPerStride+posX;

//	if ( _index == 1198899)
//		std::cout << "BUG" << std::endl;


	m_measurementsCollection[_index].push_back( val);


}

template  <typename TX> void StandardMeasurementCollectionManager<TX>::trim(int posX, int posY, int posZ)
{
	int _index = posZ*m_elementsPerSlice+posY*m_elementsPerStride+posX;

	//m_measurementsCollection[_index].push_back( val);

	// prior to saving get rid of the excess capacity to the swap trick
			// shrink to fit the std::vector
	std::vector<StandardMeasurementCollection<TX>*>(m_measurementsCollection[_index]).swap(m_measurementsCollection[_index]);
}


template  <typename TX> std::vector<StandardMeasurementCollection<TX>*>* StandardMeasurementCollectionManager<TX>::getData(int posX, int posY, int posZ)
{
	int _index = posZ*m_elementsPerSlice+posY*m_elementsPerStride+posX;

	
	std::vector<StandardMeasurementCollection<TX>*>* _ptr = &m_measurementsCollection[_index];
	
	return _ptr;
}

template  <typename TX> std::vector<StandardMeasurementCollection<TX>*>* StandardMeasurementCollectionManager<TX>::getData(int index)
{
	
	std::vector<StandardMeasurementCollection<TX>*>* _ptr = &m_measurementsCollection[index];
	
	return _ptr;
}

template  <typename TX> StandardMeasurementCollectionManager<TX>::~StandardMeasurementCollectionManager()
{
	int _max = m_maxX * m_maxY * m_maxZ;

	for(int i=0;i<m_smcBlocks.size();i++)
	{
		//std::cout << "Deleting block " << i+1 << " of " << m_smcBlocks.size() << std::endl;
		DataAndHeap<TX> _data =m_smcBlocks[i];
		HeapFree( _data.heap, HEAP_NO_SERIALIZE, _data.data );

		// now we have to kill the heap again
		HeapDestroy(_data.heap);
	}
	//int _1percent = (int)_max/100.0;
	std::vector<StandardMeasurementCollection<IntensityAndDistance>*>*  _vec = &m_measurementsCollection[0];
	/*int _counter = 0;
	int _percent = 0;*/
	/*for (int i=0; i<_max;i++)
	{
		
		int _elements = _vec->size();

		
		
		_vec++;
	}*/
	std::cout << "Deleting Collection" << std::endl;
	if ( m_measurementsCollection ) delete[] m_measurementsCollection;
}


class VoxelData
{
private:
	std::vector<Measurement> m_measurementCollection;
public:
	std::vector<Measurement>& getMeasurementCollection() { return m_measurementCollection; }
};

class DopplerVoxel : public VoxelData
{
private:
	Eigen::Vector3f m_data;
public:
	Eigen::Vector3f& getData() { return m_data; } 
	
};

class IntensityVoxel : public VoxelData
{
private:
		float m_data;
public:
	float& getData() { return m_data; }
};

class VoxelDataSet
{
private:
	std::vector<VoxelData> m_data;
	int m_width;
	int m_height;
	int m_depth;
public:
	VoxelDataSet();
	void setVoxelDimensions(int width, int height, int depth);
	VoxelDataSet(int width, int height, int depth);
};



#endif