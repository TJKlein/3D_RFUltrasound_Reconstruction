#include "glwidget.h"
#include <stdio.h>
#include <iostream>

GLWidget::GLWidget() :
m_texture(false),
m_glTexture(0),
m_timerID(0),
m_lastInterval(100)
{
    m_timerID = startTimer( 100 ); //64-65fps
    rotX = rotY = rotZ = 0.f;
    col = 0;
}

void GLWidget::refreshSceneInterval(int ms)
{
	if (m_timerID)
		killTimer(m_timerID);

	m_lastInterval = ms;
	m_timerID = startTimer(ms);
}

void GLWidget::restartTimer()
{
	if (m_timerID)
		killTimer(m_timerID);

	m_timerID = startTimer(m_lastInterval);
}

void GLWidget::stopTimer()
{
	if (m_timerID)
		killTimer(m_timerID);
}

GLWidget::~GLWidget()
{
	if ( m_glTexture) glDeleteTextures(1, &m_glTexture);
}

void GLWidget::setTexture(GLuint tex)
{
	m_texture = true;

	m_glTexture = tex;
}

void GLWidget::initializeGL()
{
    //initialization of OpenGL

    glClearColor(0.85f, 0.85f, 0.85f, 0.f);
    //resizeGL( 400 , 300 );

    glShadeModel( GL_SMOOTH );
    glHint( GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST );
    glEnable( GL_TEXTURE_2D );
    glEnable( GL_CULL_FACE );
    glEnable( GL_DEPTH_TEST );

}

void GLWidget::paintGL()
{
	
    //draw scene here

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
	if ( m_texture )
	{
		glEnable( GL_TEXTURE_2D );
	glBindTexture(GL_TEXTURE_2D, m_glTexture);
	glBegin(GL_QUADS);
	{
		// Front Face
		glTexCoord2f(0.0f, 0.0f); glVertex3d(-1.0f, -1.0f,  -1.0f);	// Bottom Left Of The Texture and Quad
		glTexCoord2f(1.0f, 0.0f); glVertex3d( 1.0f, -1.0f,  -1.0f);	// Bottom Right Of The Texture and Quad
		glTexCoord2f(1.0f, 1.0f); glVertex3d( 1.0f,  1.0f,  -1.0f);	// Top Right Of The Texture and Quad
		glTexCoord2f(0.0f, 1.0f); glVertex3d(-1.0f,  1.0f,  -1.0f);	// Top Left Of The Texture and Quad
	}
	glEnd();
	glDisable(GL_TEXTURE_2D);
	}

    glTranslated(0.0, 0.0, -1.0);

    glRotatef( rotX, 1, 0, 0 );
    glRotatef( rotY, 0, 1, 0 );
    glRotatef( rotZ, 0, 0, 1 );

      static const int coords[6][4][3] = {
        { { +1, -1, -1 }, { -1, -1, -1 }, { -1, +1, -1 }, { +1, +1, -1 } },
        { { +1, +1, -1 }, { -1, +1, -1 }, { -1, +1, +1 }, { +1, +1, +1 } },
        { { +1, -1, +1 }, { +1, -1, -1 }, { +1, +1, -1 }, { +1, +1, +1 } },
        { { -1, -1, -1 }, { -1, -1, +1 }, { -1, +1, +1 }, { -1, +1, -1 } },
        { { +1, -1, +1 }, { -1, -1, +1 }, { -1, -1, -1 }, { +1, -1, -1 } },
        { { -1, -1, +1 }, { +1, -1, +1 }, { +1, +1, +1 }, { -1, +1, +1 } }
    };

    for (int i = 0; i < 6; ++i) {
        glColor3ub( i*20+col, 100+i*10+col, i*42+col );
        glBegin(GL_QUADS);
        for (int j = 0; j < 4; ++j) {
            glVertex3d(0.2 * coords[i][j][0], 0.2 * coords[i][j][1], 0.2 * coords[i][j][2]);
        }
        glEnd();
    }
if(col > 200) col = 0;
col++;

    rotX += 1.1;
    rotY += 1.4;
    rotZ += 1.6;

}

void GLWidget::resizeGL(int width, int height)
{
    //proces resize keep good aspect ratio for 3D scene

    glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    int side = qMin(width, height);
    glViewport((width - side) / 2, (height - side) / 2, side, side);

   // glViewport(0, 0, width, height);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    //ISSUE_X64: gluPerspective(45.0, (GLdouble)width/(GLdouble)height, 0.01f, 1000.0f);
	// as we don't want to use gluPerspective, use glFrustrum instead with the values
	// matching those of gluPerspective
	float _fovy = 45.0f;
	float _aspect = (double)width/(double)height;
	float _zNear = 0.01f;
	float _zFar = 1000.0f;
	float _top= tan(_fovy*3.14159/360.0) * _zNear;
	float _bottom = -_top;
	float _left = _aspect * _bottom;
	float _right = _aspect * _top;
	glFrustum(_left, _right, _bottom, _top, _zNear, _zFar);
    glMatrixMode(GL_MODELVIEW);
}

void GLWidget::mousePressEvent(QMouseEvent *event)
{
  //proces mouse events for rotate/move inside 3D scene
}

void GLWidget::mouseMoveEvent(QMouseEvent *event)
{
 //proces keyboard events
}

void GLWidget::timerEvent(QTimerEvent *event)
{
    updateGL();
}
