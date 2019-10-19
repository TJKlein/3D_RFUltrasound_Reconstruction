#ifndef GLWIDGET_H
#define GLWIDGET_H

#include <QGLWidget>
#include <QtOpenGL>

class GLWidget : public QGLWidget
{
public:
    GLWidget();
	~GLWidget();

    double rotX, rotY, rotZ; //i want access these variable so public
    void resizeGL(int width, int height);
    short int col;
	void setTexture(GLuint test);
	void refreshSceneInterval(int ms);
	void stopTimer();
	void restartTimer();

protected:
    void initializeGL();
    void paintGL();
    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    void timerEvent(QTimerEvent *event);

private:
	bool m_texture;
	GLuint m_glTexture;
	int m_timerID;
	int m_lastInterval;

};

#endif // GLWIDGET_H
