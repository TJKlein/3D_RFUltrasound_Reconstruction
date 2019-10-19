#ifndef SLICESVISUALIZATION_H
#define SLICESVISUALIZATION_H
#include <GL/glew.h>
#include <Eigen/StdVector>
#include <QGLWidget>
#include <QtOpenGL>
#include <opencv/cv.h>
#include <opencv/highgui.h>

#include <Eigen\Core>
#include <Eigen/Geometry>
#include <Eigen/Array>
#include <Eigen/QR>
#include <Eigen/LU>

#include "GeometricTools.h"
#include "GraphicTools.h"
#include "Reconstruction.h"
#include "RFCompounding.h"

USING_PART_OF_NAMESPACE_EIGEN

class SlicesVisualization : public QGLWidget
{
public:
	 EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    SlicesVisualization();
	~SlicesVisualization();

	void stopTimer();
	void restartTimer();
	void setTrackingData(std::vector<Eigen::Transform3f, Eigen::aligned_allocator<Eigen::Transform3f> > m_trackingData);
	void setCalibrationMatrix(Eigen::Matrix4f &calibrationMatrix);
	void setImageSize(const int &minX, const int &minY, const int &maxX, const int &maxY, const double &scaleX, const double &scaleY);
	void setLine(Eigen::Vector3f start, Eigen::Vector3f stop);

    double rotX, rotY, rotZ; //i want access these variable so public
    void resizeGL(int width, int height);
    short int col;

	void drawScanlines(bool draw, int startIndex, int endIndex, std::vector<Scanline<unsigned short>*> *sclVector);
	void highlightSlice(bool highlight = false, int sliceId = 0, GLuint texture = 0, HBITMAP bmp = 0, int with = 0, int height = 0);
	void defineAxisAlignedBoundingBox(Eigen::Vector3f minPoint, Eigen::Vector3f maxPoint);
	void defineOrientedBoundingBox(BBExtent extent, Eigen::Transform3f &trans);

	void setSpacePlane(std::vector<Eigen::Vector3f> *spacePoint, std::vector<Eigen::Vector3f> *normal);
	IplImage* renderUSCroppingMask(const int &imageWidth, const int &imageHeight, const double &apexPosX, const double &apexPosY, const double &innerRadius, const double &outerRadius, const double &angle);
	void setReconstruction(Reconstruction *recon, const int &dimX, const int &dimY, const int &dimZ, const float &scaleX, const float &scaleY, const float &scaleZ);
protected:
    void initializeGL();
    void paintGL();
    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
	void wheelEvent(QWheelEvent *event);
    void timerEvent(QTimerEvent *event);
	void keyPressEvent ( QKeyEvent * event );
private:
	void drawAxisAlignedBoundingBox();
	void drawOrientedBoundingBox();
	void drawBox(Eigen::Vector3f &minPos, Eigen::Vector3f &maxPos);

private:
	// variables related to scanline visualization
	bool m_drawScanlines;
	int  m_scanlineStartIndex;
	int  m_scanlineEndIndex;
	std::vector<Scanline<unsigned short>*> *m_scanlineVector;

	GraphicTools m_graphicTools;
	Reconstruction* m_reconstruction;
	Eigen::Vector3f m_lineStart,m_lineStop;
	bool m_drawLine;
    BBExtent m_obbExtent;
	Eigen::Transform3f m_obbTransformation;
	Eigen::Vector3f m_minPos, m_maxPos;
	Eigen::Vector3f m_imageMin, m_imageMax;
	bool m_axisAlignedBoundingBox;
	bool m_orientedBoundingBox;
	GLuint m_sliceTexture;
	GLuint m_usPlanes;
	GLuint m_normalPlane;
	GLuint m_rays;
	GLuint m_velocityVoxels;
	HBITMAP m_sliceBMP;
	int m_highlightSliceId;
	bool m_highlightSlice;
	bool m_imageSizeSet;
	bool m_trackingDataSet;
	bool m_calibrationMatrixSet;
	std::vector<Eigen::Transform3f, Eigen::aligned_allocator<Eigen::Transform3f> > m_trackingData;
	Eigen::Transform3f m_calibrationMatrix;
	std::vector< Eigen::Vector4f, Eigen::aligned_allocator<Eigen::Vector4f> > m_imageCorners;
	Eigen::Vector4f m_centroid;

	// mouse interaction variables
	double m_rotationX;
	double m_rotationY;
	double m_rotationZ;
	double m_scale;
	double m_translationZ;
	QPointF m_lastPosition;
};

#endif // SlicesVisualization_H
