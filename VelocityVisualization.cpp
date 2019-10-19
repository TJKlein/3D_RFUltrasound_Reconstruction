#include "VelocityVisualization.h"
#include <stdio.h>
#include <iostream>

#include "VideoTools.h"

// code removed because of 64-bit ISSUE: ISSUE_X64

// ISSUE_X64
/* sorry GLUT was removed because of 64 bit issues
#ifdef MAC
	#include <GLUT/glut.h>
#else
	#include <GL/glut.h>
#endif
*/
//#include <glut.h>

//#include <GLTools/GLUtils.h>
#include "GeometricTools.h"

/*ISSUE_X64
#define glError() { \
	
	GLenum err = glGetError(); \
	while (err != GL_NO_ERROR) { \
	fprintf(stderr, "glError: %s caught at %s:%u\n", (char *)gluErrorString(err), __FILE__, __LINE__); \
	} \
	
}*/

VelocityVisualization::VelocityVisualization() :
m_rotationX(0.0),
m_rotationY(0.0),
m_rotationZ(0.0),
m_scale(1.0),
m_translationZ(-100.0),
m_velocityField(0),
m_plane(0)
{
	startTimer( 1509999999 ); //64-65fps
	rotX = rotY = rotZ = 0.f;
	col = 0;
	m_centroid = Eigen::Vector4f::Zero();
}

VelocityVisualization::~VelocityVisualization()
{
	if ( m_velocityField )
		glDeleteLists(m_velocityField, 1);

	if ( m_plane )
		glDeleteLists(m_plane, 1);
}

void VelocityVisualization::initializeGL()
{
	
	glewInit();
	//initialization of OpenGL

	glClearColor(0.0f, 0.0f, 0.0f, 0.f);
	//resizeGL( 400 , 300 );

	glShadeModel( GL_SMOOTH );
	//glHint( GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST );

	

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	float _aspect = static_cast<float>(this->size().width()) / static_cast<float>(this->size().height());
	glOrtho(-_aspect*650.0f, _aspect*650.0f, -650.0f, 650.0f, -1000, 1000);
	glMatrixMode(GL_MODELVIEW);

	//glEnable( GL_TEXTURE_2D );
	//glEnable( GL_CULL_FACE );
	glEnable( GL_DEPTH_TEST );
}

void VelocityVisualization::setPlane(std::vector<Eigen::Vector3f> *plane, const int &maxX, const int &maxY, const float &scaleX, const float &scaleY)
{
	if ( m_plane )
		glDeleteLists(m_plane, 1);

	m_plane = glGenLists(1);
	glNewList(m_plane, GL_COMPILE);

	for (int y=0;y<maxY;y++)
	{
		for(int x=0;x<maxX;x++)
		{
			int _index = y*maxX+x;
			Eigen::Vector3f _normal =  plane->at(_index);
			//if (_velocity[0] != 0.0 || _velocity[1] != 0.0 || _velocity[2] != 0.0 )
			{
				glBegin(GL_LINES);
				glVertex3f(static_cast<float>(x)*scaleX, static_cast<float>(y)*scaleY, 0.0f);
				glVertex3f(static_cast<float>(x)*scaleX+_normal[0], static_cast<float>(y)*scaleY+_normal[1], 0.0f);
				glEnd();
			}

		}
	}

	glEndList();
}



void VelocityVisualization::setVelocityField(Eigen::Vector3f* velocityField, const int &maxX, const int &maxY, const int &maxZ, const float &physicalX, const float &physicalY, const float &physicalZ)
{
	makeCurrent();

	GLint _maxVerts;
	glGetIntegerv( GL_MAX_ELEMENTS_VERTICES, &_maxVerts );

	if ( m_velocityField )
		glDeleteLists(m_velocityField, 1);
	// create a display list to speed up the process
	std::cout << "Generating list..." << std::endl;
	m_velocityField = glGenLists(1);
	glNewList(m_velocityField, GL_COMPILE);
	glLineWidth(1.0);
	//CAMP::GLUtils::drawRectangles3D(CAMP::Vector3<float>(static_cast<float>(maxX)*physicalX, static_cast<float>(maxY)*physicalY, static_cast<float>(maxZ)*physicalZ),CAMP::Vector4<float>(1.0,0.0,0.0,1.0),true,false);

	//TJK_removed_CAMP: std::cout << CAMP::Vector3<float>(static_cast<float>(maxX)*physicalX, static_cast<float>(maxY)*physicalY, static_cast<float>(maxZ)*physicalZ) << std::endl;

	glColor3f(0.0,1.0,0.0);
	int _counter = 0;

	Eigen::Vector3f _meanVelocity(0.0f, 0.0f, 0.0f);
	Eigen::Vector3f _maxVelocity(-10000.f,-10000.f,-10000.f);
	Eigen::Vector3f _minVelocity(10000.f,10000.f,10000.f);
	GeometricTools _gt;
	glPushMatrix();
	glTranslatef(-maxX/2.0f, -maxY/2.0, -maxZ/2.0);
	for (int z=0;z<maxZ;z+=5)
	{
		for(int y=0;y<maxY;y+=5)
		{
			for(int x=0;x<maxX;x+=5)
			{

				int _index = z*maxX*maxY+y*maxX+x;
				Eigen::Vector3f _velocity =  velocityField[_index];
				//std::cout << _velocity << "\t\t\t";
				//if (_velocity[0] != 0.0 || _velocity[1] != 0.0 || _velocity[2] != 0.0 )
				_meanVelocity+=_velocity;
				if ( _velocity.norm() > _maxVelocity.norm() )
					_maxVelocity = _velocity;

				if ( _velocity.norm() < _minVelocity.norm() )
					_minVelocity = _velocity;

				if (  _velocity.norm() < 10.0 && _velocity.norm() > 1.0e-010)
				{
					float _velocityScaleFactor = 1.f;
					//if ( _velocity.norm() > 0.0 )
					//std::cout << _velocity.norm() << std::endl;
					_counter+=1;
					glBegin(GL_LINES);
					glColor3f(0.0,1.0,0.0);
					glVertex3f(static_cast<float>(x), static_cast<float>(y), static_cast<float>(z));
					glColor3f(1.0,0.0,0.0);
					glVertex3f(static_cast<float>(x)+(_velocity[0]*_velocityScaleFactor), static_cast<float>(y)+(_velocity[1]*_velocityScaleFactor), static_cast<float>(z)+(_velocity[2]*_velocityScaleFactor));
					glEnd();

					/*glPushMatrix();
					{

					Eigen::Matrix3f _rot;
					Eigen::Transform3f _trans;
					Eigen::Vector3f _tail(static_cast<float>(x)*physicalX, static_cast<float>(y)*physicalY, static_cast<float>(z)*physicalZ);
					Eigen::Vector3f _head(static_cast<float>(x)*physicalX+_velocity[0]*0.6, static_cast<float>(y)*physicalY+_velocity[1]*0.6, static_cast<float>(z)*physicalZ+_velocity[2]*0.6);
					_gt.getRotationMatrix(_tail,_head, _rot);
					_trans.setIdentity();
					_trans.rotate( _rot );

					glMultMatrixf(_trans.data() );

					CAMP::GLUtils::drawArrow((_head-_tail).norm());

					}
					glPopMatrix();*/
				}
				//else
				{

					//	CAMP::GLUtils::drawSphere(1.5f, static_cast<float>(x)*physicalX, static_cast<float>(y)*physicalY, static_cast<float>(z)*physicalZ,1.0,0.0,0.0, 1.0,3.0);

				}
			}
		}
	}
	glPopMatrix();
	glEndList();

	_meanVelocity/=maxX*maxY*maxZ;

	std::cout << "Max velocity: " << _maxVelocity << std::endl;
	std::cout << "Min velocity: " << _minVelocity << std::endl;
	std::cout << "Mean velocity: " << _meanVelocity << std::endl;
	std::cout << "Finished list - containing " << 2*_counter << " vertices ... " << std::endl;
	glLineWidth(1.0);

	if ( 2*_counter > _maxVerts )
		std::cout << "Max number of vertices EXCEEDED!!!" << std::endl;
	std::cout << "Maximum allowed number of vertices: " <<  _maxVerts << std::endl; 

	m_centroid << static_cast<float>(maxX)/2.0*physicalX, static_cast<float>(maxY)/2.0*physicalY, static_cast<float>(maxZ)/2.0*physicalZ, 1.0;

}



void VelocityVisualization::paintGL()
{
	
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();


	//if ()
	{

		glPushMatrix();
		{
			glTranslatef(0,0,m_translationZ);

			glRotatef(m_rotationX,1.0,0.0,0.0);
			glRotatef(m_rotationY,0.0,1.0,0.0);
			glRotatef(m_rotationZ,0.0,0.0,1.0);


			glScalef(m_scale, m_scale, m_scale);

			glTranslatef(-m_centroid(0), -m_centroid(1), -m_centroid(2));

			//ISSUE_X64: m_graphicTools.drawCoordinateSystem(20.0f, "Center");
			
			if ( m_velocityField )
			{
				glCallList(m_velocityField);
			}
	
			if ( m_plane )
				glCallList(m_plane);

		}
		glPopMatrix();
	}

}


void VelocityVisualization::resizeGL(int width, int height)
{
	//proces resize keep good aspect ratio for 3f scene
	glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);



	int side = qMin(width, height);
	glViewport(0,0,this->size().width(), this->size().height());

	//float _aspect = this->size().width() / this->size().height();
	// glViewport(0, 0, width, height);

	//glFrustum(-100,100,-100,100,0,100);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	float _aspect = static_cast<float>(this->size().width()) / static_cast<float>(this->size().height());
	glOrtho(-_aspect*650.0f, _aspect*650.0f, -650.0f, 650.0f, -1000, 1000);
	//gluPerspective(45.0, (GLdouble)width/(GLdouble)height, 0.01f, 1000.0f);
	glMatrixMode(GL_MODELVIEW);
}

void VelocityVisualization::mousePressEvent(QMouseEvent *event)
{
	//proces mouse events for rotate/move inside 3f scene

	m_lastPosition = event->pos();
}

void VelocityVisualization::mouseMoveEvent(QMouseEvent *event)
{
	//proces keyboard events
	GLdouble dx = (GLdouble)(event->x() - m_lastPosition.x()) / width();
	GLdouble dy = (GLdouble)(event->y() - m_lastPosition.y()) / height();

	if ( event->buttons() & Qt::MouseButton::LeftButton)
	{
		m_rotationX += 180.0 * dy;
		m_rotationY += 180.0 * dx;
	}
	else if ( event->buttons() & Qt::MouseButton::RightButton )
	{
		m_rotationX += 180.0 * dy;
		m_rotationZ += 180.0 * dx;
	}
	m_lastPosition = event->pos();
}

void VelocityVisualization::timerEvent(QTimerEvent *event)
{
	updateGL();
}

void VelocityVisualization::wheelEvent(QWheelEvent *event)
{
	event->delta() > 0 ? m_scale += m_scale*0.1f : m_scale -= m_scale*0.1f;
	updateGL();
}

void VelocityVisualization::keyPressEvent ( QKeyEvent * event )
{	
	switch (event->key())
	{
	case Qt::Key_Up : 
		{ 
			m_translationZ += 50.0; 
			break;
		}

	case Qt::Key_Down : 
		{ 
			m_translationZ -= 50.0; 
			break;
		}
	}

	/*switch(key)
	{
	}*/
}
