#ifndef VELOCITYVISUALIZATION_H
#define VELOCITYVISUALIZATION_H
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

USING_PART_OF_NAMESPACE_EIGEN

class VelocityVisualization : public QGLWidget
{
public:
	 EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    VelocityVisualization();
	~VelocityVisualization();

    double rotX, rotY, rotZ;
    void resizeGL(int width, int height);
    short int col;


	void setPlanes(std::vector<Eigen::Transform3f, Eigen::aligned_allocator<Eigen::Transform3f> > trackingData, Eigen::Matrix4f &calibrationMatrix, int maxX, int maxY);



	void setVelocityField(Eigen::Vector3f* velocityField, const int &maxX, const int &maxY, const int &maxZ, const float &physicalX, const float &physicalY, const float &physicalZ);

	void setPlane(std::vector<Eigen::Vector3f> *plane, const int &maxX, const int &maxY, const float &scaleX, const float &scaleY);
protected:
    void initializeGL();
    void paintGL();
    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
	void wheelEvent(QWheelEvent *event);
    void timerEvent(QTimerEvent *event);
	void keyPressEvent ( QKeyEvent * event );

private:
	GraphicTools m_graphicTools;
	unsigned int m_velocityField;
	unsigned int m_plane;
	Eigen::Vector4f m_centroid;
	// mouse interaction variables
	double m_rotationX;
	double m_rotationY;
	double m_rotationZ;
	double m_scale;
	double m_translationZ;
	QPointF m_lastPosition;
};

#endif // VelocityVisualization_H
