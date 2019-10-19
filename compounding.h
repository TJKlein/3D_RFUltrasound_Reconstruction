#ifndef DOPPLERCOMPOUNDING_H
#define DOPPLERCOMPOUNDING_H

#include <Eigen/StdVector>
#include <QtWidgets/QMainWindow>
#include "ui_compounding.h"
#include "GraphicTools.h"
#include "glwidget.h"
#include "SlicesVisualization.h"
#include "VelocityVisualization.h"
#include "GraphicTools.h"
#include "VideoTools.h"
#include "DataTools.h"
#include "XMLTools.h"
#include "GeometricTools.h"
#include "Reconstruction.h"
#include "NumericTools.h"

#include "RFCompounding.h"

#include <vtkImageData.h>
#include <vtkSmartPointer.h> 
#include <vtkMetaImageWriter.h> 
#include <vtkPointData.h>
#include <vtkImageCast.h>
#include <qt_windows.h>
#include <QtWinExtras/QWinFunctions.h>
#include <QPixmap>

//#include <amplio.h>
//#include <pando.h>

//#include <vld.h>

// XML name tags for configuration
const char COMPOUNDING_CONFIGURATION[] = {"Configuration.xml"};
const char DATA_DIRECTORY_TAG[] = {"DataDirectory"};
const char CALIBRATION_DIRECTORY_TAG[] = {"CalibrationDirectory"};
const char TARGET_DIRECTORY_TAG[]={"TargetDirectory"};


const char RF_DATA_DIRECTORY_TAG[] = {"RFDataDirectory"};
const char RF_CALIBRATION_DIRECTORY_TAG[] = {"RFCalibrationDirectory"};
const char RF_TRACKING_DIRECTORY_TAG[] = {"RFTrackingDirectory"};
const char RF_TARGET_DIRECTORY_TAG[]={"RFTargetDirectory"};

const char TRACKING_FILENAME_TAG[] = {"TrackingFilename"};
const char CALIBRATION_FILENAME_TAG[] = {"CalibrationFilename"};
const char VIDEO_FILENAME_TAG[] = {"VideoFilename"};

const char RF_TARGET_FILENAME_TAG[] = {"RFTargetFilename"};
const char RF_FILENAME_TAG[] = {"RFFilename"};
const char RF_TRACKING_FILENAME_TAG[] = {"RFTrackingFilename"};
const char RF_CALIBRATION_FILENAME_TAG[] = {"RFCalibrationFilename"};

const char RF_INDEXFILE_TAG[] = {"IndexFile"};
const char RF_DIRECTIONFILE_TAG[] = {"DirectionFile"};
const char RF_RECONSTRUCTIONFILE_TAG[] = {"ReconstructionFile"};


typedef unsigned char COMPOUNDING_DATA_TYPE;


//#define IMAGE_FILTER 0

// using DirectX requires Windows SDK (Platform SDK) and DirectX 9.1
// the following lib files have then to be included: comdlg32.lib Strmbase.lib
#define USE_DIRECTX 0

#define _HAS_ITERATOR_DEBUGGING 0
#define _NOTHREADS

#ifdef _WIN32
#pragma warning ( disable : 4251 4311 4312 4995 )
#endif
#include <QtWidgets/QMenu.h>
#include <QtWidgets/QLabel.h>
#include <QObject.h>
#include <QtWidgets/QMenuBar.h>
#include <QtWidgets/QFileDialog.h>
#include <QtWidgets/QGridLayout.h>
#include <QtWidgets/QTextEdit.h>
#ifdef _WIN32
#pragma warning ( default : 4251 4311 4312 4995 )
#endif


/**
 * \brief 3D Reconstruction software of ultrasound data
 *
 * \authors Tassilo Klein
 *          <br>
 *          kleint\@cs.tum.edu
 * \ingroup Ultrasound
 * \version 1.0a
 * \date 02-03-2009
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

// not the entire video data may be stored into one contiguous block of memory
// therefore it is broken-up into BLOCKS which can contain a limited amount of SLICES
#define NUMBER_OF_BLOCKS 25
#define NUMBER_OF_SLICES 750.0f
#define DISTANCE_SCALAR 0.7f // 0.55f


static HANDLE hMutex;

class CompoundVolumeThread;

struct Coordinate
{
	int X;
	int Y;
	int Z;
	Eigen::Vector3f Velocity;
};

enum COMPOUNDING_MODE { GAUSSIAN, NEAREST_NEIGHBOR, MEDIAN, WEIGHTED_MEDIAN, INVERSE_DISTANCE, MEAN, SIGMOID };

class CompoundRFSubvolumeWritingThread // :public QThread
{
	//Q_OBJECT

private:
	//std::fstream m_rfDataFile;
	HANDLE m_rfDataFile;
	int m_maxElements;
	int m_voxelsX, m_voxelsY, m_voxelsZ;
	__int64 m_fileSize;
	std::string m_filename;
	StandardMeasurementCollectionManager<IntensityAndDistance> *m_measurementsCollectionManager;
public:
	CompoundRFSubvolumeWritingThread(std::string filename,int voxelsX, int voxelsY, int voxelsZ, int maxElements, __int64 size);
	~CompoundRFSubvolumeWritingThread();

	void run();
	void setDataTarget(StandardMeasurementCollectionManager<IntensityAndDistance> *ptr);
};

class CompoundRFVolumeThread : public QThread
{
	Q_OBJECT

	signals:
		void percentEvent();
	private:
	
		float m_physicalOffsetX, m_physicalOffsetY, m_physicalOffsetZ;
		int m_offsetX, m_offsetY, m_offsetZ;
		int m_voxelsX;
		int m_voxelsY;
		int m_voxelsZ;
		float m_physicalX;
		float m_physicalY;
		float m_physicalZ;
		int m_numberOfScanlines;
		BBExtent m_extent;
		float m_voxelSize;
		int m_1percent;
#ifdef USE_USHORT16
		std::vector<Scanline<unsigned short>*>* m_scanlineVector;
#endif

#ifdef USE_FLOAT32
		std::vector<Scanline<float>*>* m_scanlineVector;
#endif
		Eigen::Transform3f m_boundingBoxTransformation;
		int *m_finishedIndex;

		 // UI input data
		int m_maxDistanceScalar;
		int m_maxDistanceUnit;
		int m_index;
		__int64 *m_byteCounter;
		std::string m_filename;

		// here all the measurements are stored
		StandardMeasurementCollectionManager<IntensityAndDistance> *m_measurementsCollectionManager;
		
	public:
		
		CompoundRFVolumeThread();
		void run();
		
		void setFilename(std::string filename);
		void setData(std::vector<Scanline<unsigned short>*>* scanlineVector);
		void setData(std::vector<Scanline<float>*>* scanlineVector);
		void setDataTarget(StandardMeasurementCollectionManager<IntensityAndDistance> *ptr);
		void setPhysicalOffset(float offsetX, float offsetY, float offsetZ);
		void setIndex(int index);
		void setByteCounter(__int64 *byteCounter);
		void setDataSize(const int dimX, const int dimY, const int dimZ, const float physicalX, const float physicalY, const float physicalZ, BBExtent extent, int offsetX, int offsetY, int offsetZ);
		void setUserConfiguration(int maxDistanceScalar, int maxDistanceUnit, int _1percent);
		void setTransformation(Eigen::Transform3f &boundingBoxTransformation);
		
};


class CompoundVolumeThread : public QThread
 {
	 Q_OBJECT

	 signals:
		void percentEvent();

 private:
	 COMPOUNDING_MODE m_compoundingMode;
	 vtkImageData *m_targetImage;
	 int m_voxelsX;
	 int m_voxelsY;
	 int m_voxelsZ;
	 int m_roiMaxX, m_roiMaxY, m_roiMinX, m_roiMinY, m_roiHeight, m_roiWidth;
	 int m_offsetX, m_offsetY, m_offsetZ;
	 float m_physicalX;
	 float m_physicalY;
	 float m_physicalZ;
	 float m_scaleX, m_scaleY;
	 VideoTools *m_videoTools;
	 std::vector< Eigen::Matrix4f, Eigen::aligned_allocator<Eigen::Matrix4f> > *m_matrixDataTracking;
	 std::vector<Eigen::Transform3f, Eigen::aligned_allocator<Eigen::Transform3f> > *m_trackingData;
	 std::vector<std::vector<IplImage*>> *m_imageVector;

#ifdef IMAGE_FILTER
	 std::vector<std::vector<IplImage*>> *m_filterVector;
#endif
	 BBExtent m_extent;
	 float m_voxelSize;
	 Eigen::Transform3f m_boundingBoxTransformation;
	 IplImage* m_croppingMask;
	 // UI input data
     int m_maxDistanceScalar;
	 int m_maxDistanceUnit;
	 float m_gaussianSigma;
	 float m_inverseDistanceFactor;
	 int m_1percent;
	 int m_imageCounter;
 public:
		void setUserConfiguration(int maxDistanceScalar, int maxDistanceUnit, float gaussianSigma, float inverseDistanceFactor, int _1percent);
		void setROI(int roiMaxX, int roiMaxY, int roiMinX, int roiMinY, int roiHeight, int roiWidth);
		void setDataPointers( std::vector<std::vector<IplImage*>> *imageVector, std::vector<std::vector<IplImage*>> *filterVector, std::vector< Eigen::Matrix4f, Eigen::aligned_allocator<Eigen::Matrix4f> > *matrixDataTracking, std::vector<Eigen::Transform3f, Eigen::aligned_allocator<Eigen::Transform3f> > *trackingData);
		void setAuxiliaryPointers(VideoTools *videoTools, IplImage* croppingMask);
		void setTransformation(Eigen::Transform3f &boundingBoxTransformation);
		void setCompoundingMode(COMPOUNDING_MODE mode);
		void setTargetImage(vtkImageData* image);
		void setDataSize(const int dimX, const int dimY, const int dimZ, const float physicalX, const float physicalY, const float physicalZ, BBExtent extent, int offsetX, int offsetY, int offsetZ, int imageCounter);
		void setUSResolution(float scaleX, float scaleY);
        void run();
 };



// for standard reconstruction VIDEO_MODE::STANDARD, for Doppler reconstruction the RGB is mapped to actual  velocity employing the color LUT
enum VIDEO_MODE { STANDARD, DOPPLER };

enum BATCH_MODE { RF, BMODE};
enum BATCH_DISTANCE_TYPE {VOXEL_EXTENT, VOXEL_DIAMETER};
enum BATCH_DISTANCE { point5, one, onePoint1, onePoint2, onePoint3, onePoint4, onePoint5 };
class Compounding : public QMainWindow
{
	Q_OBJECT

public:
	Compounding(QWidget *parent = 0, Qt::WindowFlags flags = 0);
	void batchIntensityCompound(std::string filename,  COMPOUNDING_MODE compoundingMode, BATCH_DISTANCE_TYPE type, BATCH_DISTANCE distance, int logCompression);
	void batchRun(BATCH_MODE mode, std::string calibrationFile, std::string configurationFile, std::string trackingFile, float voxelSize, BATCH_DISTANCE_TYPE type, BATCH_DISTANCE distance, COMPOUNDING_MODE compoundingMode, std::string targetName, std::string referenceName);
	void batchRunBMODE(BATCH_MODE mode, std::string calibrationFile, std::string outputFileName, std::vector<std::string> videoFile, std::vector<std::string> trackingFile, float voxelSize, BATCH_DISTANCE_TYPE type, BATCH_DISTANCE distance, COMPOUNDING_MODE compoundingMode, bool logCompress, std::string targetName, std::string referenceName);
	
	void batchButterwothEnvelope(std::string filename, int lowCutOff, int highCutOff, int order, int frequency);
	~Compounding();

private:

	 // information to be stored in the XML file
	 std::vector<std::string> m_imageFilesXML;
	 Eigen::Matrix4f m_calibratinMatrixXML;
	 int m_masterFrameXML;

	// Ultrasonix stuff
	//pando *m_pando;
	//amplio *m_amplio;

	// working thread that writes the sub-volume files to disk
	std::vector<CompoundRFSubvolumeWritingThread*>	m_subvolumeWriting; 
	//--------
	// RF data related variables and methods
	//Eigen::Matrix4f m_matCalibrationApex;
	//int m_scanlines;
	//int m_samplesAlongScanline;
	//int m_rfFrameSizeElements;
	//int m_rfFrameSizeBytes;
	unsigned int m_rfScanlineIndex; // each beam direction is indexed with this variable, which is later saved in the 3d reconstruction file

	//float m_penetrationDepth; // mm
	// data conversion
#ifdef USE_USHORT16
	std::vector<RFData<unsigned short>*> m_rfDataVector; // pointer to all RF images (which contain at least the geometrical position, data itself can be zero)
	std::vector<Scanline<unsigned short>*> m_scanlineDataVector; // vector containing all the scanlines (which contain at least the geometrical position, data itself can be zero)
#endif

#ifdef USE_FLOAT32 || USE_FLOAT5X32
	std::vector<RFData<float>*> m_rfDataVector; // pointer to all RF images (which contain at least the geometrical position, data itself can be zero)
	std::vector<Scanline<float>*> m_scanlineDataVector; // vector containing all the scanlines (which contain at least the geometrical position, data itself can be zero)
#endif
	/**
	\brief This function computes the positions of all scanlines in space and puts it into a vector
	**/
	void computeScanlines(unsigned char settingsIndex, unsigned char fileIndex);

	//-------
	/**
	\brief Loads video data into the corresponding video data vector
	\param VIDEO_MODE Determines the color of the video, STANARD is GRAY, DOPPLER is color
	\param attach Determines if the video data will be attached to the current video vector (default false)
	\param auxiliaryFile Determines whether an auxilliary file will be loaded (default false)
	**/
	bool loadVideoFromFile(float sigma, VIDEO_MODE mode = STANDARD, bool attach = false, bool auxiliaryFile = false, bool applyDopplerNoise = false);
	bool loadTrackingDataFromFile(QString file, bool RF, bool useReferenceTarget, bool attach = false);
	bool loadTimeStampsFromFile(QString filename, bool RF, bool attach = false);
	bool loadSettingFromXMLFile(QString filename);
	bool loadSettingFromXMLFile(QString filename, UltrasoundSettings &ultrasoundSettings);
	void computeInitialVectorField();



	/**
	\brief Given a target and apex position on the US plane as well as velocity direction in world coordinate system, determine the jacobian for propagating the uncertainty from tracking
	**/
	Eigen::VectorXf jacobianNormalUncertainty(Eigen::Vector3f &target, Eigen::Vector3f &apex, Eigen::Vector3f &velocity);

	/**
	\brief Crops the mask to the minimum size (region of interest)
	**/
	bool cropMask();

	/**
	\brief Determines the area of interest within the video image, from which later on reconstruction will be performed
	**/
	bool determineCroppingArea();

	/**
	\brief Determines the size of the video image
	**/
	bool determineImageSize();

	/**
	\brief Determines the minimum bounding box of the 2D slice data
	**/
	bool determineMinimumBoundingBoxFromBMode();

	/**
	\brief Determines the minimum bounding box given a transformation matrix
	**/
	bool determineMinimumBoundingBoxManualFromBMode();

	/**
	\brief Determine the minimum bounding box of RF data
	**/
	bool determineMinimumBoundingBoxfromRF();

	/**
	\brief Determine the minimum bounding box of RF data from manual transformation
	**/
	bool determineManualBoundingBoxfromRF();

	/**
	\brief Determines a bounding box, by setting a specific slice as parallel to coordinate plane
	**/
	bool determineBoundingBoxFromMasterFrame();

	/**
	\brief Determines a bounding box, by setting a specific slice as parallel to coordinate plane
	**/
	bool determineMinimumBoundingBoxfromMasterFrameRF();

	bool determineVoxelSize();

	bool determineVoxelSizeRF();

	/**
	\brief Reconstruction of the velocity field from 2D US Doppler images
	**/
	bool compoundVolumeVelocity(float tau, float lambda1, float lambda2, int iterations);

	/**
	\brief Reconstruction of the 3D volume from 2D ultrasound slices
	**/
	bool compoundVolume(COMPOUNDING_MODE _mode);


	/**
	\brief Reconstruction of the 3D volume from RD data
	**/
	bool compoundRFVolume(COMPOUNDING_MODE _mode);

	bool computeMeasurementCollection();

	bool suggestVoxelSize();
	// the vector containing is subdivided into memory chunks in order to avoid 
	// the allocation of one giant contiguous block of memory (which might not even be possible for long vidos..except if you got plenty of RAM)
	std::vector<std::vector<IplImage*>> m_imageVector;
	std::vector<std::vector<IplImage*>> m_auxImageVector;
	#ifdef IMAGE_FILTER
	 std::vector<std::vector<IplImage*>> m_filterVector;
	#endif
	Ui::CompoundingClass ui;
	GLWidget *glWidget;
	SlicesVisualization *m_offlineVisualization;
	SlicesVisualization *m_slicesVisualization;
	VelocityVisualization *m_velocityVisualization; // Visualization of the reconstructed flow from Doppler data


	std::vector<FrameScanlineRange> m_frameScanlineRangeVector;

	// compounding related variables
	int m_imageWidth;
	int m_imageHeight;
	int m_imageMinX;
	int m_imageMaxX;
	int m_imageMinY;
	int m_imageMaxY;

	bool m_logCompress;

	// pixel to mm
	float m_scaleX;
	float m_scaleY;

	float m_physicalX, m_physicalY, m_physicalZ;

	// US fan geometry
	float m_apexX, m_apexY;
	float m_innerRadius, m_outerRadius;
	float m_angle;

	// Temporal offset
	float m_temporalOffset;
	int m_roiMaxX, m_roiMaxY, m_roiMinX, m_roiMinY, m_roiHeight, m_roiWidth;
	CvRect m_roiRect;
	BBExtent m_extent;
	Eigen::Transform3f m_boundingBoxTransformation;
	//TJK_removed_CAMP: CAMP::Matrix4<float> m_boundingBoxTransformationCAMP;

	int m_voxelsX, m_voxelsY, m_voxelsZ;
	float m_voxelSize;
	std::vector<TrackingData> m_vecTracking;
	/*std::vector<Eigen::Matrix4f> m_matrixDataTracking;
	std::vector<Eigen::Transform3f> m_trackingData;*/

	std::vector< Eigen::Matrix4f, Eigen::aligned_allocator<Eigen::Matrix4f> > m_matrixDataTracking; // with calibration matrix post-multiplied
	std::vector<Eigen::Transform3f, Eigen::aligned_allocator<Eigen::Transform3f> > m_trackingData; // only the tracking data!!

	std::vector<float> m_videoTimeStampVector;

	Eigen::Matrix4f m_matCalibration;
	
	IplImage* m_croppingMask;
	HBITMAP m_imageForVisualization;
	HBITMAP m_masterImage;
	HBITMAP m_rfImage;
	HBITMAP m_rfPreviewImage;

	// Doppler specific variables
	ColorLUT* m_colLUT; 
	Reconstruction *m_recon;

	// GUI related variables

	bool	m_guiVideoReady;
	bool	m_guiTrackingReady;
	bool	m_guiDataLoaded;
	bool	m_guiConfigurationReady;

	// RF GUI related variables
	bool m_guiRFCalibrationReady;
	bool m_guiRFDataReady;
	bool m_guiRFTrackingReady;
	
	// RF data related variables
	QString m_rfCalibrationFileName;
	QString m_rfDataFileName;
	QString m_rfTrackingFileName;
	QString m_rfRFDfilename;

	// RF processing related variables
	std::vector<RFFileStream*> m_rfFileStreamVector;

	// Data related variables

	QString	m_videoFileName;
	QString m_trackingFileName;
	QString m_configurationFileName;
	QString m_mhdFileName;

	// Tool libraries
	NumericTools *m_numericTools;
	GeometricTools m_geometricTools;
	VideoTools	*m_videoTools;
	VideoTools	m_videoSegmentationTools;
	GraphicTools	m_graphicTools;
	DataTools	m_dataTools;
	XMLTools	*m_xmlTools;
	XMLTools	m_xmlConfiguration;
	QString m_lastPath;
	QString m_lastCalibrationPath;
	QString m_lastTargetPath;

	// RF specific data
	QString m_lastCalibrationRFPath;
	QString m_lastTargetRFPath;
	QString m_lastDataRFPath;
	QString m_lastTrackingRFPath;

	double m_threadingFactor;
	CompoundVolumeThread *m_cvThread;
	CompoundRFVolumeThread *m_crfvThread;
	int m_percentage;

	unsigned char m_settingsIndex;
	unsigned char m_fileIndex;

	// variables required for attaching video and tracking data
	int m_fromVideoTimeStamp;
	int m_toVideoTimeStamp;
	int m_imageCounter;

	int m_rfFrameCounter;

	// this variable keeps the starting index in the vector of tracking data, where data was added from
	int m_trackingIndex;

	std::vector<UltrasoundSettings> m_ultrasoundSettings;
	std::vector<Eigen::Matrix4f> m_calibrationMatrixVector;
	std::vector<int> m_indexVector; // the index is the index of the corresponding ultrasoundsettings as well as calibration matrix

	public Q_SLOTS:
		/**
		@Brief Notification when the mode for computing the transformation matrix ix changed (automatic, from slice, manual)
		**/
		void clickedRFRadio();
		void clickedBModeRadio();
		void updatePercentage();
		void loadData();
		void selectVideoFile();
		void selectTrackingFile();
		void selectConfigurationFile();
		void guiUpdated();
		void showFrame(int);
		void compound();
		void initializeData();
		void updateVoxelDimension(int);
		void useReferenceFrame(int);
		void attachData();
		void setoutputFileName();
		void preprocessDoppler();
		void restoreVectorField();
		void setLastSettings();
		void addNoise();
		void imageViewerUpdate(int);
		void masterFrame(int);
		void masterFrameRF(int);
		void masterSliceChanged(QString);
		void finishedOutputFilenameEditing();
		void addBookmark();
		void doubleClickBookmark(QListWidgetItem *item);

		// RF specific signals
		void loadRFData();
		void selectRFFile();
		void compoundRF();
		void initializeRFData(); //
		void attachRFData();
		void selectConfigurationFileRF();
		void setLastSettingsRF();
		void selectOutputRF();
		void selectTrackingFileRF();
		void postProcessRF(std::string postfix = "");
		void preProcessRF();

		void selectPostProcessingFile();
		void extractDataRF();
		void loadFileForPostProcessingRF();
		void rfFrameChanged(int);
		void voxelSubsampling(int);
		void withinVoxelSubsampling(int);
		void RFLogCompression(int);
		void exportTrackingData();
};


#endif // COMPOUNDING_H
