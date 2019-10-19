#include "SlicesVisualization.h"
#include <stdio.h>
#include <iostream>
#include <math.h>

#include "VideoTools.h"

/* ISSUE_X64
#ifdef MAC
	#include <GLUT/glut.h>
#else
	#include <GL/glut.h>
#endif
	*/

//#include <glut.h>

//#include <GLTools/GLUtils.h>


#define glError() { \
	GLenum err = glGetError(); \
	while (err != GL_NO_ERROR) { \
		fprintf(stderr, "glError: %s caught at %s:%u\n", (char *)gluErrorString(err), __FILE__, __LINE__); \
		err = glGetError(); \
	} \
}

SlicesVisualization::SlicesVisualization() :
m_imageSizeSet(false),
m_trackingDataSet(false),
m_calibrationMatrixSet(false),
m_rotationX(0.0),
m_rotationY(0.0),
m_rotationZ(0.0),
m_scale(1.0),
m_translationZ(-100.0),
m_highlightSliceId(0),
m_highlightSlice(false),
m_sliceTexture(0),
m_axisAlignedBoundingBox(false),
m_orientedBoundingBox(false),
m_drawLine(false),
m_normalPlane(0),
m_rays(0)
{
	startTimer( 150 ); //64-65fps
	rotX = rotY = rotZ = 0.f;
	col = 0;
	m_centroid = Eigen::Vector4f::Zero();
}

SlicesVisualization::~SlicesVisualization()
{
}

void SlicesVisualization::restartTimer()
{
	startTimer(150);
}

void SlicesVisualization::stopTimer()
{
	
		startTimer(1000000000000000000);
}


void SlicesVisualization::setReconstruction(Reconstruction *recon, const int &dimX, const int &dimY, const int &dimZ, const float &scaleX, const float &scaleY, const float &scaleZ)
{
	makeCurrent();
	m_reconstruction = recon;

	
	if ( glIsList(m_velocityVoxels) )
		glDeleteLists(m_velocityVoxels, 1);

	m_velocityVoxels = glGenLists(1);
	glNewList(m_velocityVoxels, GL_COMPILE);

	int _index = 0;
	int _voxelsWithVelocity = 0;

glEnable (GL_BLEND); 
	{
		glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	for(int _z=0;_z<dimZ;_z++)
			for(int _y=0;_y<dimY;_y++)
				for(int _x=0;_x<dimX;_x++)
				{
					VelocityData *_data = m_reconstruction->getMeasurementsByIndex(_index);
					_index++;

					bool _draw = false;
					for(int i=0; i< _data->velocityVector.rows();i++)
					{
						if ( abs(_data->velocityVector(i)) >= 0.01)
						{
							_draw = true;
							continue;
						}
					}

					if ( _draw )
					{
						_voxelsWithVelocity++;
						float _col = _data->approximationError/m_reconstruction->getMaxError();

						glColor4f(_col,1.0-_col,0.0,2.0);
						glPushMatrix();
						{
							glTranslatef(static_cast<float>(_x)*scaleX, static_cast<float>(_y)*scaleY, static_cast<float>(_z)*scaleZ);
							// ISSUE_X64: glutWireCube(scaleX);
						}
						glPopMatrix();
					}
				}
	}
				glDisable(GL_BLEND);
	glEndList();

	std::cout << _voxelsWithVelocity << " out of " << dimX*dimY*dimZ << " voxels contain color data!" << std::endl;
}

void SlicesVisualization::setLine(Eigen::Vector3f start, Eigen::Vector3f stop)
{
	m_drawLine = true;

	m_lineStart = start;
	m_lineStop = stop;
}

void SlicesVisualization::drawOrientedBoundingBox()
{
	double _maxX = m_obbExtent.dimX;
	double _minX = 0.0;

	double _maxY = m_obbExtent.dimY;
	double _minY = 0.0;

	double _maxZ = m_obbExtent.dimZ;
	double _minZ = 0.0;

	Eigen::Vector3f _minPos, _maxPos;

	_minPos << _minX,_minY,_minZ;
	_maxPos << _maxX/2.0,_maxY/2.0,_maxZ/2.0;

	
	
}


void SlicesVisualization::drawBox(Eigen::Vector3f &_minPos, Eigen::Vector3f &_maxPos)
{
	glBegin( GL_QUADS);
	{
		glVertex3f(_maxPos(0), _maxPos(1), _minPos(2));
		glVertex3f(_maxPos(0), _maxPos(1), _maxPos(2));
		glVertex3f(_maxPos(0), _minPos(1), _maxPos(2));
		glVertex3f(_maxPos(0), _minPos(1), _minPos(2));

		glVertex3f(_minPos(0), _maxPos(1), _minPos(2));
		glVertex3f(_minPos(0), _maxPos(1), _maxPos(2));
		glVertex3f(_minPos(0), _minPos(1), _maxPos(2));
		glVertex3f(_minPos(0), _minPos(1), _minPos(2));

		//-----

		glVertex3f(_minPos(0), _maxPos(1), _minPos(2));
		glVertex3f(_maxPos(0), _maxPos(1), _minPos(2));
		glVertex3f(_maxPos(0), _minPos(1), _minPos(2));
		glVertex3f(_minPos(0), _minPos(1), _minPos(2));

		glVertex3f(_minPos(0), _maxPos(1), _maxPos(2));
		glVertex3f(_maxPos(0), _maxPos(1), _maxPos(2));
		glVertex3f(_maxPos(0), _minPos(1), _maxPos(2));
		glVertex3f(_minPos(0), _minPos(1), _maxPos(2));

		//-----

		glVertex3f(_minPos(0), _minPos(1), _minPos(2));
		glVertex3f(_minPos(0), _minPos(1), _maxPos(2));
		glVertex3f(_maxPos(0), _minPos(1), _maxPos(2));
		glVertex3f(_maxPos(0), _minPos(1), _minPos(2));

		glVertex3f(_minPos(0), _maxPos(1), _minPos(2));
		glVertex3f(_minPos(0), _maxPos(1), _maxPos(2));
		glVertex3f(_maxPos(0), _maxPos(1), _maxPos(2));
		glVertex3f(_maxPos(0), _maxPos(1), _minPos(2));
	}
	glEnd();
}

void SlicesVisualization::drawAxisAlignedBoundingBox()
{
	glColor4f(1.0, 1.0, 1.0,0.2);
	glBegin( GL_QUADS );
	{
		glVertex3f(m_maxPos(0), m_maxPos(1), m_minPos(2));
		glVertex3f(m_maxPos(0), m_maxPos(1), m_maxPos(2));
		glVertex3f(m_maxPos(0), m_minPos(1), m_maxPos(2));
		glVertex3f(m_maxPos(0), m_minPos(1), m_minPos(2));

		glVertex3f(m_minPos(0), m_maxPos(1), m_minPos(2));
		glVertex3f(m_minPos(0), m_maxPos(1), m_maxPos(2));
		glVertex3f(m_minPos(0), m_minPos(1), m_maxPos(2));
		glVertex3f(m_minPos(0), m_minPos(1), m_minPos(2));

		//-----

		glVertex3f(m_minPos(0), m_maxPos(1), m_minPos(2));
		glVertex3f(m_maxPos(0), m_maxPos(1), m_minPos(2));
		glVertex3f(m_maxPos(0), m_minPos(1), m_minPos(2));
		glVertex3f(m_minPos(0), m_minPos(1), m_minPos(2));

		glVertex3f(m_minPos(0), m_maxPos(1), m_maxPos(2));
		glVertex3f(m_maxPos(0), m_maxPos(1), m_maxPos(2));
		glVertex3f(m_maxPos(0), m_minPos(1), m_maxPos(2));
		glVertex3f(m_minPos(0), m_minPos(1), m_maxPos(2));

		//-----

		glVertex3f(m_minPos(0), m_minPos(1), m_minPos(2));
		glVertex3f(m_minPos(0), m_minPos(1), m_maxPos(2));
		glVertex3f(m_maxPos(0), m_minPos(1), m_maxPos(2));
		glVertex3f(m_maxPos(0), m_minPos(1), m_minPos(2));

		glVertex3f(m_minPos(0), m_maxPos(1), m_minPos(2));
		glVertex3f(m_minPos(0), m_maxPos(1), m_maxPos(2));
		glVertex3f(m_maxPos(0), m_maxPos(1), m_maxPos(2));
		glVertex3f(m_maxPos(0), m_maxPos(1), m_minPos(2));
	}
	glEnd();
}

void SlicesVisualization::defineAxisAlignedBoundingBox(Eigen::Vector3f minPoint, Eigen::Vector3f maxPoint)
{
	m_minPos = minPoint;
	m_maxPos = maxPoint;

	m_axisAlignedBoundingBox = true;
}

void SlicesVisualization::defineOrientedBoundingBox(BBExtent extent, Eigen::Transform3f &trans)
{
	m_obbExtent = extent;
	m_obbTransformation = trans;
	
	m_orientedBoundingBox = true;
}

void SlicesVisualization::highlightSlice(bool highlight, int sliceId, GLuint texture, HBITMAP bmp, int width, int height)
{
	m_highlightSlice = highlight;
	m_highlightSliceId = sliceId;
	m_sliceTexture = texture;
	m_sliceBMP = bmp;
	VideoTools _vtools;

	/*int _minX, _maxX, _minY, _maxY;
	_vtools.getImageMinMax(bmp, _minX, _maxX, _minY, _maxY);
	m_imageMin << _minX,_minY,0;
	m_imageMax << _maxX,_maxY,0;*/

	

	m_imageMin << 0,0,0;
	m_imageMax << width,height,0;
}

void SlicesVisualization::setImageSize(const int &minX, const int &minY, const int &maxX, const int &maxY, const double &scaleX, const double &scaleY)
{
	m_imageCorners.clear();

	Eigen::Vector4f _tempPoint;

	_tempPoint << minX*scaleX, maxY*scaleY,0,1;
	m_imageCorners.push_back(_tempPoint);

	_tempPoint << minX*scaleX, minY*scaleY,0,1;
	m_imageCorners.push_back(_tempPoint);

	_tempPoint << maxX*scaleX, minY*scaleY,0,1;
	m_imageCorners.push_back(_tempPoint);

	_tempPoint << maxX*scaleX, maxY*scaleY,0,1;
	m_imageCorners.push_back(_tempPoint);

	_tempPoint << minX*scaleX, maxY*scaleY,0,1;
	m_imageCorners.push_back(_tempPoint);

	m_imageSizeSet = true;
}

void SlicesVisualization::setTrackingData(std::vector<Eigen::Transform3f, Eigen::aligned_allocator<Eigen::Transform3f> > trackingData)
{

	m_trackingDataSet = true;


	m_trackingData = trackingData;

	for(int i=0;i<trackingData.size();i++)
	{
		for(int j=0;j<m_imageCorners.size()-1;j++)
		{
			m_centroid += trackingData[i]*m_calibrationMatrix*m_imageCorners[j];

		}
		//std::cout << (trackingData[i]*m_calibrationMatrix).matrix() << std::endl;

	}

	m_centroid /= m_centroid(3);
}

void SlicesVisualization::setCalibrationMatrix(Eigen::Matrix4f &calibrationMatrix)
{
	m_calibrationMatrix = calibrationMatrix;

	m_calibrationMatrixSet = true;
}

void SlicesVisualization::initializeGL()
{
	glewInit();
	//initialization of OpenGL

	glClearColor(0.85f, 0.85f, 0.85f, 0.f);
	//resizeGL( 400 , 300 );

	glShadeModel( GL_SMOOTH );
	//glHint( GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST );


	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-500, 500, -500, 500, -1000, 1000);
	glMatrixMode(GL_MODELVIEW);


	//glEnable( GL_TEXTURE_2D );
	//glEnable( GL_CULL_FACE );
	glEnable( GL_DEPTH_TEST );

	

}

IplImage* SlicesVisualization::renderUSCroppingMask(const int &imageWidth, const int &imageHeight, const double &apexPosX, const double &apexPosY, const double &innerRadius, const double &outerRadius, const double &angle)
{
	makeCurrent();
	/*GraphicTools _gt;
	return _gt.createUSCroppingMask(imageWidth,imageHeight,apexPosX, apexPosY, innerRadius, outerRadius, angle);

	*/

// NOTE: must be called within a valid GL context!!

	//glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);

	/*double _apexPosX = 500.0, _apexPosY = 0.0;
	double _innerRadius = 100.0, _outerRadius=500.0;
	double _angle=30.0;
	double _scaleX=1.0, _scaleY=1.0;*/
		IplImage* _grayImage = 0;

	glewInit();





	// offscreen rendering
	GLuint _fbo = 0;
	GLuint _depthbuffer;
	GLuint _offlineTexture;

	// We define our frame buffer object for offline rendering and then storing it into bitmap

	// Setup our FBO 1
	glGenFramebuffersEXT(1, &_fbo);
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, _fbo);

	// Create the render buffer for depth	
	glGenRenderbuffersEXT(1, &_depthbuffer);
	glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, _depthbuffer);
	glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, GL_DEPTH_COMPONENT, imageWidth, imageHeight);


	GLint _maxbuffers;
	glGetIntegerv(GL_MAX_COLOR_ATTACHMENTS_EXT, &_maxbuffers);


	// Now setup a texture to render to
	glGenTextures(1, &_offlineTexture);
	glBindTexture(GL_TEXTURE_2D, _offlineTexture);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8,  imageWidth, imageHeight, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

	//  The following 3 lines enable mipmap filtering and generate the mipmap data so rendering works
	//	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	//	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
	//	glGenerateMipmapEXT(GL_TEXTURE_2D);

	// And attach it to the FBO so we can render to it
	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_2D, _offlineTexture, 0);

	GLenum status = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
	if(status != GL_FRAMEBUFFER_COMPLETE_EXT)
	{
		// do some error handling here
	}

	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);	// Unbind the FBO for now


	// now draw the US fan into the FBO, which will then later on be used for masking


	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, _fbo);
	{

		glPushAttrib(GL_ALL_ATTRIB_BITS);
		{
			glViewport(0, 0, imageWidth, imageHeight);

			glMatrixMode(GL_PROJECTION);
			glPushMatrix();
			{
				glLoadIdentity();

				//ISSUE_X64
				glOrtho(0.0,  imageWidth, imageHeight,0.0, -1.0, 1.0 );

				glMatrixMode(GL_MODELVIEW);
				glPushMatrix();
				{
					glLoadIdentity();
					
					// viewport has to correspond to FBO texture size!!!!
					glViewport(0, 0, imageWidth, imageHeight );
					
					//glClearColor(0.0,0.0,0.0,0.0);
					glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

					glDrawBuffer(GL_COLOR_ATTACHMENT0_EXT);
					

					glPushMatrix();
					{

						glColor3f(0.0,0.0,0.0);
						glBegin(GL_QUADS);
						glVertex3f(0.0,0,0);
						glVertex3f(imageWidth,0,0);
						glVertex3f(imageWidth,imageHeight,0);
						glVertex3f(0,imageHeight,0);
						glEnd();

						glColor4f(1.0f, 1.0f, 1.0f, 1.0);
						/*CAMP::UltrasoundGeometry probe;
						m_croppingArea->copyToGeometry(probe);*/

						// draw circle segment
						//ISSUE_X64: GLUquadricObj *quadObject = NULL;
						//ISSUE_X64: if (! quadObject)	quadObject = gluNewQuadric();
	
						//glScalef(_scaleX, _scaleY, 1.0);
						//glLoadIdentity();
						//ISSUE_X64: gluQuadricDrawStyle(quadObject, GLU_FILL);


						glTranslatef(apexPosX, apexPosY, 0.0);



						m_graphicTools.drawPartialDisc(innerRadius, outerRadius, 360.0, angle, 0.1);

						//ISSUE_X64: gluPartialDisk( quadObject, innerRadius,  outerRadius, 32, 1, 360-angle, angle*2.0);

						//ISSUE_X64: gluDeleteQuadric(quadObject);
					}
					glPopMatrix();

				}
				glPopMatrix();
			}
			glPopMatrix();


			// now save the FBO to a file

			GLubyte* image  =0;
			image = new GLubyte[imageWidth * imageHeight * 3];
			glReadBuffer(GL_COLOR_ATTACHMENT0_EXT); 
			glReadPixels(0, 0, imageWidth, imageHeight, GL_RGB, GL_UNSIGNED_BYTE, image); 

//			writeStandardBitmap("c:\\test.bmp", image, imageWidth, imageHeight);
			
		

			VideoTools _videoTools;
			IplImage *_tmpImage=  _videoTools.uchar2IplImage(image, imageWidth, imageHeight,3,8);

			//cvSaveImage("renderedCroppingMaskColor.bmp",_tmpImage);
			_grayImage = cvCreateImage(cvSize(imageWidth, imageHeight), IPL_DEPTH_8U, 1);
			cvCvtColor(_tmpImage, _grayImage, CV_RGB2GRAY);
			cvReleaseImage(&_tmpImage);
			cvSaveImage("renderedCroppingMask.jpg",_grayImage);
			//cvShowImage("test",_grayImage);
			//cvReleaseImage(&_img);
			
			delete image;
		}
		glPopAttrib();

	}
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);


	// free the memory of the FBO and the textures

	glDeleteFramebuffersEXT(1, &_fbo);
	glDeleteRenderbuffersEXT(1, &_depthbuffer);
	glDeleteTextures(1,&_offlineTexture);


	return _grayImage;
	
}

void SlicesVisualization::paintGL()
{
	

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();


	if ( m_imageSizeSet && m_trackingDataSet && m_calibrationMatrixSet )
	{

		glPushMatrix();
		{
			glTranslatef(0,0,m_translationZ);

			glRotatef(m_rotationX,1.0,0.0,0.0);
			glRotatef(m_rotationY,0.0,1.0,0.0);
			glRotatef(m_rotationZ,0.0,0.0,1.0);



			glScalef(m_scale, m_scale, m_scale);


			
if ( m_drawScanlines )
	{
		/*for (int i=m_scanlineStartIndex; i<m_scanlineEndIndex;i++)
		{
			glBegin(GL_LINES);
			{
				Scanline* _scl = (m_scanlineVector->at(i));
				glColor3f(0.0,1.0,0.0);
				glVertex3f(_scl->getLine()->getStartPoint().x(),_scl->getLine()->getStartPoint().y(),_scl->getLine()->getStartPoint().z());
				glColor3f(1.0,0.0,0.0);
				glVertex3f(_scl->getLine()->getEndPoint().x(),_scl->getLine()->getEndPoint().y(),_scl->getLine()->getEndPoint().z());

				//glVertex3f
			}
			glEnd();
		}*/
		/*
		if ( m_rays )
				glCallList(m_rays);
	*/
	}


glTranslatef(-m_centroid(0), -m_centroid(1), -m_centroid(2));
			
			/*glPushMatrix();
			{
				glTranslatef(m_minPos(0),m_minPos(1),m_minPos(2));
				glutSolidSphere(1.0,10,10);
			}
			glPopMatrix();
*/


			if ( m_drawLine)
			{
			glPushMatrix();
			{
				glLineWidth(2.0);
				glColor4f(1.0,0.0,0.7,1.0);
				
				glBegin(GL_LINES);
				{
					glVertex3f(m_lineStart(0), m_lineStart(1), m_lineStart(2));
					glVertex3f(m_lineStop(0), m_lineStop(1), m_lineStop(2));
				}
				glEnd();

				glLineWidth(1.0);
			}
			glPopMatrix();
			}
/*
			glPushMatrix();
			{
				glTranslatef(m_maxPos(0),m_maxPos(1),m_maxPos(2));
				glutSolidSphere(1.0,10,10);
			}
			glPopMatrix();
*/

			glPushMatrix();
			{


			for(int i=0;i<m_trackingData.size();i++)
			{
				if ( m_highlightSlice && i == m_highlightSliceId )
				{
					glColor3f(0.0, 0.0, 0.0);
					glLineWidth(4.0);
					if ( m_sliceTexture != 0 )
					{


						glPushMatrix();
						{
							glMultMatrixf(m_trackingData[i].data());
							glMultMatrixf(m_calibrationMatrix.data());

							glEnable( GL_TEXTURE_2D );
							//glBindTexture(GL_TEXTURE_2D, m_sliceTexture);
							glLineWidth(5.0);
							glBegin(GL_QUADS);
							{
								for(int j=0;j<m_imageCorners.size()-1;j++)
								{
									// if highlight of slice is set and texture is available, draw it
									if ( m_highlightSlice && m_sliceTexture != 0)
									{
										switch(j)
										{
										case 2: glTexCoord2f(0.0f, 0.0f); break;
										case 3: glTexCoord2f(1.0f, 0.0f); break;
										case 0: glTexCoord2f(1.0f, 1.0f); break;
										case 1: glTexCoord2f(0.0f, 1.0f); break;
										default: break;
										}
									}
									Eigen::Vector4f _tempPoint = m_imageCorners[j];
									glVertex3f(_tempPoint(0), _tempPoint(1), _tempPoint(2) );
								}
							}
							glEnd();
							glDisable(GL_TEXTURE_2D);



							/*glColor3f(1.0,0.0,1.0);
							glLineWidth(5.0);

							glBegin(GL_QUADS);
							{
								glVertex3f(m_imageMin(0), m_imageMin(1), m_imageMin(2)-1);
								glVertex3f(m_imageMax(0), m_imageMin(1), m_imageMin(2)-1);
								glVertex3f(m_imageMax(0), m_imageMax(1), m_imageMin(2)-1);
								glVertex3f(m_imageMin(0), m_imageMax(1), m_imageMin(2)-1);

								glVertex3f(m_imageMin(0), m_imageMin(1), m_imageMin(2)+1);
								glVertex3f(m_imageMax(0), m_imageMin(1), m_imageMin(2)+1);
								glVertex3f(m_imageMax(0), m_imageMax(1), m_imageMin(2)+1);
								glVertex3f(m_imageMin(0), m_imageMax(1), m_imageMin(2)+1);
							}
							glEnd();

							glLineWidth(1.0);*/
						}
						glPopMatrix();


					}
					else // just draw the frame
					{
						glPushMatrix();
						{
							glMultMatrixf(m_trackingData[i].data());
							glMultMatrixf(m_calibrationMatrix.data());
							glBegin(GL_QUADS);
							{
								for(int j=0;j<m_imageCorners.size();j++)
								{
									Eigen::Vector4f _tempPoint = m_imageCorners[j];
									glVertex3f(_tempPoint(0), _tempPoint(1), _tempPoint(2) );
								}
							}
							glEnd();
						}
						glPopMatrix();
					}
				}
				else
				{
					glEnable (GL_BLEND); 
					{
						glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
						glColor4f(0.0, 0.5, i*1.0/m_trackingData.size(),0.2);
						glLineWidth(1.0);

						glPushMatrix();
						{
							glMultMatrixf(m_trackingData[i].data());
							glMultMatrixf(m_calibrationMatrix.data());
							glBegin(GL_LINE_STRIP);
							{
								for(int j=0;j<m_imageCorners.size();j++)
								{
									Eigen::Vector4f _tempPoint = m_imageCorners[j];
									glVertex3f(_tempPoint(0), _tempPoint(1), _tempPoint(2) );
								}
							}
							glEnd();
						}
						glPopMatrix();

					}
					glDisable(GL_BLEND);

				}

			}

			if ( m_axisAlignedBoundingBox )
			{
				glEnable (GL_BLEND); 
				{
					glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
					drawAxisAlignedBoundingBox();
				}
				glDisable(GL_BLEND);
			}

			}
			glPopMatrix();

			if ( m_orientedBoundingBox )
			{
				glEnable (GL_BLEND); 
				{
					glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
					drawOrientedBoundingBox();
				}
				glDisable(GL_BLEND);
			}

		}
		glPopMatrix();
	}
	
}

void SlicesVisualization::drawScanlines(bool draw, int startIndex, int endIndex, std::vector<Scanline<unsigned short>*> *sclVector)
{

	makeCurrent();

	if ( glIsList(m_rays) )
		glDeleteLists(m_rays, 1);

	m_rays = glGenLists(1);
	glNewList(m_rays, GL_COMPILE);

	m_drawScanlines = draw;
	m_scanlineStartIndex = startIndex;
	m_scanlineEndIndex = endIndex;
	m_scanlineVector = sclVector;

	for (int i=m_scanlineStartIndex; i<m_scanlineEndIndex;i++)
		{
			glBegin(GL_LINES);
			{
				Scanline<unsigned short>* _scl = (m_scanlineVector->at(i));
				glColor3f(0.0,1.0,0.0);
				glVertex3f(_scl->getLine()->getStartPoint().x(),_scl->getLine()->getStartPoint().y(),_scl->getLine()->getStartPoint().z());
				glColor3f(1.0,0.0,0.0);
				glVertex3f(_scl->getLine()->getEndPoint().x(),_scl->getLine()->getEndPoint().y(),_scl->getLine()->getEndPoint().z());

				//glVertex3f
			}
			glEnd();
		}

	glEndList();
}


void SlicesVisualization::resizeGL(int width, int height)
{
	//proces resize keep good aspect ratio for 3f scene

	glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	int side = qMin(width, height);
	glViewport(0,0,this->size().width(), this->size().height());

	// glViewport(0, 0, width, height);

	//glFrustum(-1000,1000,1000,-1000,-100000,1000000);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	glOrtho(-500, 500, -500, 500, -1000, 1000);
	//gluPerspective(45.0, (GLdouble)width/(GLdouble)height, 0.01f, 1000.0f);
	glMatrixMode(GL_MODELVIEW);
}

void SlicesVisualization::mousePressEvent(QMouseEvent *event)
{
	//proces mouse events for rotate/move inside 3f scene

	m_lastPosition = event->pos();
}

void SlicesVisualization::mouseMoveEvent(QMouseEvent *event)
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

void SlicesVisualization::timerEvent(QTimerEvent *event)
{
	updateGL();
}

void SlicesVisualization::wheelEvent(QWheelEvent *event)
{
	event->delta() > 0 ? m_scale += m_scale*0.1f : m_scale -= m_scale*0.1f;
	updateGL();
}

void SlicesVisualization::setSpacePlane(std::vector<Eigen::Vector3f> *spacePoint, std::vector<Eigen::Vector3f> *normal)
{
	makeCurrent();

	if ( glIsList(m_normalPlane) )
		glDeleteLists(m_normalPlane, 1);

	m_normalPlane = glGenLists(1);
	glNewList(m_normalPlane, GL_COMPILE);
	glLineWidth(1.0);
	for (int i=0;i<spacePoint->size();i+=3)
	{
		Eigen::Vector3f _startPoint = spacePoint->at(i);
		Eigen::Vector3f _endPoint = normal->at(i);

		// draw the vector field of the normals (for each pixel of the texture)
		glPushMatrix();
		{
			glTranslatef(_startPoint.x(), _startPoint.y(),_startPoint.z());
			//glutSolidSphere(0.1,10,10);

			// draw the vector field, one set of arrows in front, and one behind the US plane
			glBegin(GL_LINES);
			glColor3f(0.0,1.0,0.0);
			glVertex3f(0.01, 0.01, 0.01);
			glColor3f(1.0,0.0,0.0);
			glVertex3f(_endPoint.x()+0.01, _endPoint.y()+0.01, _endPoint.z()+0.01);
			glEnd();

			glBegin(GL_LINES);
			glColor3f(0.0,1.0,0.0);
			glVertex3f(-0.01, -0.01, -0.01);
			glColor3f(1.0,0.0,0.0);
			glVertex3f(_endPoint.x()-0.01, _endPoint.y()-0.01, _endPoint.z()-0.01);
			glEnd();
		}
		glPopMatrix();
	}

	glEndList();
}

void SlicesVisualization::keyPressEvent ( QKeyEvent * event )
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
