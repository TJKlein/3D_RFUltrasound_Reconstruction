#ifndef GRAPHICTOOLS_H
#define GRAPHICTOOLS_H

#include <GL/glew.h>
//#include <glut.h>
#include <windows.h>
#include <Eigen/StdVector>
#include <opencv/cv.h>
#include <opencv/highgui.h>

#define PI 3.1415926536

class Segment
{
public:
	int beginX, endX;
	int Y;
};


class GraphicTools
{

	private:
		TCHAR* StringToTCHAR(std::string& s);
	bool smartInsert(std::vector<Segment> *SegmentVector, const int &y, const int &beginX, const int &endX, const int &radius);
	public:
		GraphicTools(){};
		~GraphicTools(){};
		bool Bitmap2GLTexture(HBITMAP bmp,GLuint &textureID);
		IplImage* createUSCroppingMask(const int &imageWidth, const int &imageHeight, const double &apexPosX, const double &apexPosY, const double &innerRadius, const double &outerRadius, const double &angle);
		void GraphicTools::CreateBMPFile(std::string filename, PBITMAPINFO pbi, HBITMAP hBMP, HDC hDC);
		PBITMAPINFO CreateBitmapInfoStruct(HBITMAP hBmp);
		void writeStandardBitmap(const char *filename, GLubyte *image, int width, int height);
		void bresenhamCircle(const int &x0, const int &y0, const int &radius, std::vector<Segment> *SegmentVector);
		void draw3DText(float x, float y, float z, void *font, const char *string);
		//ISSUE_X64: void drawCoordinateSystem(float size, std::string label);
		//ISSUE_X64: void drawCoordinateSystem(float size);
		void drawArrow(float length);
		void drawPartialDisc(float innerRadius, float outerRadius, float angle, float startAngle, float stepSize);
}; 

#endif