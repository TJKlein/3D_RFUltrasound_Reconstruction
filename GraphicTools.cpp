#include "GraphicTools.h"
#include "VideoTools.h"
//#include "Common/ImageBase.h"
#include <stdio.h>
#include <iostream>
#include <fstream> 

 

bool GraphicTools::smartInsert(std::vector<Segment> *SegmentVector, const int &y, const int &beginX, const int &endX, const int &radius)
{
	Segment *_seg = &SegmentVector->at(radius-y);

	if ( _seg->beginX > beginX )
		_seg->beginX = beginX;

	if ( _seg->endX < endX )
		_seg->endX = endX;

	return true;
}

/**

\brief Bresenham (Mid-point) circle rasterization algorithms. Computes the outer points of a circle, avoiding complex computation following an iterative approach.

\param x0 Origin x coordinate of the circle
\param y0 Origin y coordinate of the circle
\param radius Radius of the circle
\param SegmentVector Vector containing the circle's Segment segments

**/
void GraphicTools::bresenhamCircle(const int &x0, const int& y0, const int &radius, std::vector<Segment> *SegmentVector)
{
	Segment _Segment;
	// if radius is zero, simply return the mid-point Segment
	// for nearest-neighbor
	if ( radius == 0 )
	{
		_Segment.beginX = x0;
		_Segment.endX = x0;
		_Segment.Y = y0;

		SegmentVector->push_back(_Segment);
		return;
	}
	int f = 1 - radius;
	int ddF_x = 1;
	int ddF_y = -2 * radius;
	int x = 0;
	int y = radius;


	for (int c=0;c<2*radius+1;c++)
	{
		Segment _empty;
		_empty.beginX = 1000000;
		_empty.endX = -1000000;
		_empty.Y = radius-c;
		SegmentVector->push_back(_empty);
	}


	// top pixel of the circle
	_Segment.beginX = x0;
	_Segment.endX = x0;
	_Segment.Y = y0 + radius;

	//SegmentVector->push_back(_Segment);
	smartInsert(SegmentVector,_Segment.Y, _Segment.beginX, _Segment.endX, radius);


	// bottom pixel of the circle
	_Segment.Y = y0 - radius;

	//SegmentVector->push_back(_Segment);
	smartInsert(SegmentVector,_Segment.Y, _Segment.beginX, _Segment.endX, radius);

	_Segment.Y = y0;

	// right-most pixel of the circle
	_Segment.beginX = x0 + radius;
	_Segment.endX = x0 + radius;

	//SegmentVector->push_back(_Segment);
	smartInsert(SegmentVector,_Segment.Y, _Segment.beginX, _Segment.endX, radius);

	// left-most pixel of the circle
	_Segment.beginX = x0 - radius;
	_Segment.endX = x0 - radius;

	//SegmentVector->push_back(_Segment);
	smartInsert(SegmentVector,_Segment.Y, _Segment.beginX, _Segment.endX, radius);



	while(x < y)
	{
		// ddF_x == 2 * x + 1;
		// ddF_y == -2 * y;
		// f == x*x + y*y - radius*radius + 2*x - y + 1;
		if(f >= 0) 
		{
			y--;
			ddF_y += 2;
			f += ddF_y;
		}
		x++;
		ddF_x += 2;
		f += ddF_x;    

		// define the circle Segments, they are equal by mirroring in 4 quadrants

		_Segment.beginX = x0 - x;
		_Segment.endX = x0 + x;
		_Segment.Y = y0 + y;

		//SegmentVector->push_back(_Segment);
		smartInsert(SegmentVector,_Segment.Y, _Segment.beginX, _Segment.endX, radius);

		_Segment.Y = y0 - y;

		//SegmentVector->push_back(_Segment);
		smartInsert(SegmentVector,_Segment.Y, _Segment.beginX, _Segment.endX, radius);

		_Segment.beginX = x0 - y;
		_Segment.endX = x0 + y;
		_Segment.Y = y0 + x;

		//SegmentVector->push_back(_Segment);
		smartInsert(SegmentVector,_Segment.Y, _Segment.beginX, _Segment.endX, radius);

		_Segment.beginX = x0 - y;
		_Segment.endX = x0 + y;
		_Segment.Y = y0 - x;

		//SegmentVector->push_back(_Segment);
		smartInsert(SegmentVector,_Segment.Y, _Segment.beginX, _Segment.endX, radius);
	}
}



IplImage* GraphicTools::createUSCroppingMask(const int &imageWidth, const int &imageHeight, const double &apexPosX, const double &apexPosY, const double &innerRadius, const double &outerRadius, const double &angle)
{

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

						glColor3d(0.0,0.0,0.0);
						glBegin(GL_QUADS);
						glVertex3d(0.0,0,0);
						glVertex3d(imageWidth,0,0);
						glVertex3d(imageWidth,imageHeight,0);
						glVertex3d(0,imageHeight,0);
						glEnd();

						glColor4d(1.0f, 1.0f, 1.0f, 1.0);
						/*CAMP::UltrasoundGeometry probe;
						m_croppingArea->copyToGeometry(probe);*/

						// draw circle segment
						//ISSUE_X64: GLUquadricObj *quadObject = NULL;
						//ISSUE_X64: if (! quadObject)	quadObject = gluNewQuadric();
						glPushMatrix();
						//glScalef(_scaleX, _scaleY, 1.0);
						//glLoadIdentity();
						//ISSUE_X64: gluQuadricDrawStyle(quadObject, GLU_FILL);


						glTranslatef(apexPosX, apexPosY, 0.0);



						//ISSUE_X64: gluPartialDisk( quadObject, innerRadius,  outerRadius, 32, 1, 360-angle, angle*2.0);


						drawPartialDisc(innerRadius, outerRadius, 360.0, angle, 0.1);

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

			//writeStandardBitmap("c:\\test.bmp", image, imageWidth, imageHeight);



			VideoTools _videoTools;
			IplImage *_tmpImage=  _videoTools.uchar2IplImage(image, imageWidth, imageHeight,3,8);
			_grayImage = cvCreateImage(cvSize(imageWidth, imageHeight), IPL_DEPTH_8U, 1);
			cvCvtColor(_tmpImage, _grayImage, CV_RGB2GRAY);
			cvReleaseImage(&_tmpImage);
			//cvSaveImage("C:\\testIMG.jpg",_img);
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

void GraphicTools::drawPartialDisc(float innerRadius, float outerRadius, float startAngle, float angle, float stepSize)
{

	//std::ofstream out("test.txt");

	glColor3f(1.0,1.0,1.0);
	glBegin(GL_QUADS);
	{
		float x, y,z;
		for(float k = startAngle-angle; k < startAngle+angle; k=k+stepSize)
		{

			x = innerRadius*sin(PI/180.0*(k));
			y = innerRadius*cos(PI/180.0*(k));
			z = 0;
			glVertex3f(x,y,z);	

			x = outerRadius*sin(PI/180.0*(k));
			y = outerRadius*cos(PI/180.0*(k));
			z = 0;
			glVertex3f(x,y,z);

			//out << "(" << x << "," << y << ")" << std::endl;



			x = outerRadius*sin(PI/180.0*(k+stepSize));
			y = outerRadius*cos(PI/180.0*(k+stepSize));
			z = 0;
			glVertex3f(x,y,z);

			//out << "\t\t(" << x << "," << y << ")" << std::endl;


			x = innerRadius*sin(PI/180.0*(k+stepSize));
			y = innerRadius*cos(PI/180.0*(k+stepSize));
			z = 0;
			glVertex3f(x,y,z);	





		}
	}
	glEnd();

	//out.close();
}


void GraphicTools::writeStandardBitmap(const char *filename, GLubyte *image, int width, int height)
{
	FILE *file;

	if((file = fopen(filename, "wb"))==NULL)
		printf("Error while saving data into file %s.\n",filename);
	else 
	{
		// ************* write BMP-Header
		fprintf(file,"BM");               // all BMP-Files start with "BM"
		unsigned int header[3];
		header[0] = 54+width*height*3;		// filesize = 54 (header) + size * size *channels
		header[1] = 0;						// reserved = 0
		header[2] = 54;						// File offset to Raster Data
		fwrite(header, 4, 3, file);

		// ************* write BMP-Info-Header
		unsigned int infoHeader[10];
		infoHeader[0] = 40;	                // size of info header 
		infoHeader[1] = width;              // Bitmap Width
		infoHeader[2] = height;             // Bitmap Height
		infoHeader[3] = 1+65536*8*3;        // first 2 bytes=Number of Planes (=1)
		// next  2 bytes=BPP
		infoHeader[4] = 0;					// compression (0 = none)
		infoHeader[5] = 0;					// compressed file size (0 if no compression)
		infoHeader[6] = 0;					// horizontal resolution: Pixels/meter
		infoHeader[7] = 0;					// vertical resolution: Pixels/meter
		infoHeader[8] = 0;					// Number of actually used colors
		infoHeader[9] = 0;					// Number of important colors  0 = all		
		fwrite(infoHeader, 4, 10, file);

		// for some reason the data in BMP is stored BGR, so convert RGB to BGR
		GLubyte *tempImage;
		tempImage = new GLubyte[width*height*3];

		for (int i = 0; i < width*height*3; i += 3)
		{
			tempImage[i]   = image[i+2];
			tempImage[i+1] = image[i+1];
			tempImage[i+2] = image[i];
		}

		// ************* write Data
		fwrite(tempImage, sizeof(GLubyte), width*height*3, file);

		fclose(file);

		delete[] tempImage;
	}
}


TCHAR* GraphicTools::StringToTCHAR(std::string& s)
{
  const char* all = s.c_str();
  int len = 1 + strlen(all);
  wchar_t* t = new wchar_t[len]; 
  if (NULL == t) throw std::bad_alloc();
  mbstowcs(t, all, len);
  return (TCHAR*)t;
}


void GraphicTools::CreateBMPFile(std::string filename, PBITMAPINFO pbi, HBITMAP hBMP, HDC hDC) 
{ 
	LPTSTR pszFile = (LPTSTR)StringToTCHAR(filename);

 HANDLE hf;                 // file handle 
 BITMAPFILEHEADER hdr;       // bitmap file-header 
 PBITMAPINFOHEADER pbih;     // bitmap info-header 
 LPBYTE lpBits;              // memory pointer 
 DWORD dwTotal;              // total count of bytes 
 DWORD cb;                   // incremental count of bytes 
 BYTE *hp;                   // byte pointer 
 DWORD dwTmp; 
 
 pbih = (PBITMAPINFOHEADER) pbi; 
 lpBits = (LPBYTE) GlobalAlloc(GMEM_FIXED, pbih->biSizeImage);
 
 if (!lpBits) 
 return; 
 
 // Retrieve the color table (RGBQUAD array) and the bits 
 // (array of palette indices) from the DIB. 
 if (!GetDIBits(hDC, hBMP, 0, (WORD) pbih->biHeight, lpBits, pbi, 
 DIB_RGB_COLORS)) 
 {
 return;
 }
 
 // Create the .BMP file. 
 hf = CreateFile(pszFile, 
 GENERIC_READ | GENERIC_WRITE, 
 (DWORD) 0, 
 NULL, 
 CREATE_ALWAYS, 
 FILE_ATTRIBUTE_NORMAL, 
 (HANDLE) NULL); 
 if (hf == INVALID_HANDLE_VALUE) 
 return; 
 hdr.bfType = 0x4d42;        // 0x42 = "B" 0x4d = "M" 
 // Compute the size of the entire file. 
 hdr.bfSize = (DWORD) (sizeof(BITMAPFILEHEADER) + 
 pbih->biSize + pbih->biClrUsed 
 * sizeof(RGBQUAD) + pbih->biSizeImage); 
 hdr.bfReserved1 = 0; 
 hdr.bfReserved2 = 0; 
 
 // Compute the offset to the array of color indices. 
 hdr.bfOffBits = (DWORD) sizeof(BITMAPFILEHEADER) + 
 pbih->biSize + pbih->biClrUsed 
 * sizeof (RGBQUAD); 
 
 // Copy the BITMAPFILEHEADER into the .BMP file. 
 if (!WriteFile(hf, (LPVOID) &hdr, sizeof(BITMAPFILEHEADER), 
 (LPDWORD) &dwTmp,  NULL)) 
 {
 return; 
 }
 
 // Copy the BITMAPINFOHEADER and RGBQUAD array into the file. 
 if (!WriteFile(hf, (LPVOID) pbih, sizeof(BITMAPINFOHEADER) 
 + pbih->biClrUsed * sizeof (RGBQUAD), 
 (LPDWORD) &dwTmp, ( NULL))) 
 return; 
 
 // Copy the array of color indices into the .BMP file. 
 dwTotal = cb = pbih->biSizeImage; 
 hp = lpBits; 
 if (!WriteFile(hf, (LPSTR) hp, (int) cb, (LPDWORD) &dwTmp,NULL)) 
 return; 
 
 // Close the .BMP file. 
 if (!CloseHandle(hf)) 
 return; 
 
 // Free memory. 
 GlobalFree((HGLOBAL)lpBits);
}

PBITMAPINFO GraphicTools::CreateBitmapInfoStruct(HBITMAP hBmp)
{ 
 BITMAP bmp; 
 PBITMAPINFO pbmi; 
 WORD    cClrBits;
 // Retrieve the bitmap color format, width, and height. 
 if (!GetObject(hBmp, sizeof(BITMAP), (LPSTR)&bmp)) 
 return NULL;
 
 // Convert the color format to a count of bits. 
 cClrBits = (WORD)(bmp.bmPlanes * bmp.bmBitsPixel); 
 if (cClrBits == 1) 
 cClrBits = 1; 
 else if (cClrBits <= 4) 
 cClrBits = 4; 
 else if (cClrBits <= 8) 
 cClrBits = 8; 
 else if (cClrBits <= 16) 
 cClrBits = 16; 
 else if (cClrBits <= 24) 
 cClrBits = 24; 
 else cClrBits = 32; 
 
 // Allocate memory for the BITMAPINFO structure. (This structure 
 // contains a BITMAPINFOHEADER structure and an array of RGBQUAD 
 // data structures.) 
 
 if (cClrBits != 24) 
 pbmi = (PBITMAPINFO) LocalAlloc(LPTR, 
 sizeof(BITMAPINFOHEADER) + 
 sizeof(RGBQUAD) * (1<< cClrBits)); 
 
 // There is no RGBQUAD array for the 24-bit-per-pixel format. 
 
 else 
 pbmi = (PBITMAPINFO) LocalAlloc(LPTR, 
 sizeof(BITMAPINFOHEADER)); 
 
 // Initialize the fields in the BITMAPINFO structure. 
 
 pbmi->bmiHeader.biSize = sizeof(BITMAPINFOHEADER); 
 pbmi->bmiHeader.biWidth = bmp.bmWidth; 
 pbmi->bmiHeader.biHeight = bmp.bmHeight; 
 pbmi->bmiHeader.biPlanes = bmp.bmPlanes; 
 pbmi->bmiHeader.biBitCount = bmp.bmBitsPixel; 
 if (cClrBits < 24) 
 pbmi->bmiHeader.biClrUsed = (1<<cClrBits); 
 
 // If the bitmap is not compressed, set the BI_RGB flag. 
 pbmi->bmiHeader.biCompression = BI_RGB; 
 
 // Compute the number of bytes in the array of color 
 // indices and store the result in biSizeImage. 
 // For Windows NT, the width must be DWORD aligned unless 
 // the bitmap is RLE compressed. This example shows this. 
 // For Windows 95/98/Me, the width must be WORD aligned unless the 
 // bitmap is RLE compressed.
 pbmi->bmiHeader.biSizeImage = ((pbmi->bmiHeader.biWidth * cClrBits +31) & ~31) /8
 * pbmi->bmiHeader.biHeight; 
 // Set biClrImportant to 0, indicating that all of the 
 // device colors are important. 
 pbmi->bmiHeader.biClrImportant = 0; 
 return pbmi; 
} 



bool GraphicTools::Bitmap2GLTexture(HBITMAP bmp,GLuint &textureID)
{
	// empty image ?
	if ( bmp == NULL ) return false;

	// acquire bitmap information
	BITMAP BM;
	::GetObject(bmp, sizeof(BM), &BM);

	// bitmap has 24 bits per pixel ?
	if ( BM.bmBitsPixel != 24 )	return false;

	// texture already existing ?
	if (!textureID)
		glGenTextures(1, &textureID);
	else glDeleteTextures(1, &textureID);


	glBindTexture(GL_TEXTURE_2D, textureID);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 4);
	glPixelStorei(GL_UNPACK_ROW_LENGTH, 0);
	glPixelStorei(GL_UNPACK_SKIP_ROWS, 0);
	glPixelStorei(GL_UNPACK_SKIP_PIXELS, 0);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);

	//ISSUE_X64: gluBuild2DMipmaps(GL_TEXTURE_2D, 3, BM.bmWidth, BM.bmHeight, GL_BGR_EXT, GL_UNSIGNED_BYTE, BM.bmBits);

	glTexImage2D(GL_TEXTURE_2D,0,3,BM.bmWidth, BM.bmHeight,0, GL_RGB, GL_UNSIGNED_BYTE, BM.bmBits);

	glBindTexture(GL_TEXTURE_2D, 0);

	return true;
}

/* ISSUE_X64
void GraphicTools::drawCoordinateSystem(float size)
{
// Store GL state.
glPushAttrib(GL_ALL_ATTRIB_BITS);

glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);


// Bitmap representations of an "X", an "Y", and a "Z" for the axis cross.
GLubyte xbmp[] = { 0x11,0x11,0x0a,0x04,0x0a,0x11,0x11 };
GLubyte ybmp[] = { 0x04,0x04,0x04,0x04,0x0a,0x11,0x11 };
GLubyte zbmp[] = { 0x1f,0x10,0x08,0x04,0x02,0x01,0x1f };

double xpos[3];
double ypos[3];
double zpos[3];

GLint viewport[4];								// Where The Viewport Values Will Be Stored
glGetIntegerv(GL_VIEWPORT, viewport);			// Retrieves The Viewport Values (X, Y, Width, Height)

GLdouble modelview[16];							// Where The 16 Doubles Of The Modelview Matrix Are To Be Stored
glGetDoublev(GL_MODELVIEW_MATRIX, modelview);	// Retrieve The Modelview Matrix

GLdouble projection[16];						// Where The 16 Doubles Of The Projection Matrix Are To Be Stored
glGetDoublev(GL_PROJECTION_MATRIX, projection);	// Retrieve The Projection Matrix

gluProject(size+(size/50.0f), 0, 0, modelview, projection, viewport, &xpos[0], &xpos[1], &xpos[2]);
gluProject(0, size+(size/50.0f), 0, modelview, projection, viewport, &ypos[0], &ypos[1], &ypos[2]);
gluProject(0, 0, size+(size/50.0f), modelview, projection, viewport, &zpos[0], &zpos[1], &zpos[2]);


// Render the cross.
glLineWidth(2.0);

enum { XAXIS, YAXIS, ZAXIS };
int idx[3] = { XAXIS, YAXIS, ZAXIS };

for (int i=0; i < 3; i++)
{
glPushMatrix();
if (idx[i] == XAXIS) {                       // X axis.
glColor3f(0.500f, 0.125f, 0.125f);
} else if (idx[i] == YAXIS) {                // Y axis.
glRotatef(90, 0, 0, 1);
glColor3f(0.125f, 0.500f, 0.125f);
} else {                                     // Z axis.
glRotatef(-90, 0, 1, 0);
glColor3f(0.125f, 0.125f, 0.500f);
}
drawArrow(size);
glPopMatrix();
}

// Render axis notation letters ("X", "Y", "Z").
glMatrixMode(GL_PROJECTION);
glPushMatrix();			/// push PROJECTION matrix
glLoadIdentity();
glOrtho(0, viewport[2], 0, viewport[3], -1, 1);

glMatrixMode(GL_MODELVIEW);
glPushMatrix();			/// push MODELVIEW matrix
glLoadIdentity();

GLint unpack;
glGetIntegerv(GL_UNPACK_ALIGNMENT, &unpack);
glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

//glColor3f(0.8f, 0.8f, 0.0f);
glColor3f(1.0f, 1.0f, 1.0f);

glRasterPos2d(xpos[0]-viewport[0], xpos[1]);
glBitmap(8, 7, 0, 0, 0, 0, xbmp);
glRasterPos2d(ypos[0]-viewport[0], ypos[1]);
glBitmap(8, 7, 0, 0, 0, 0, ybmp);
glRasterPos2d(zpos[0]-viewport[0], zpos[1]);
glBitmap(8, 7, 0, 0, 0, 0, zbmp);

glPixelStorei(GL_UNPACK_ALIGNMENT, unpack);


glPopMatrix();			/// pop MODELVIEW matrix
glMatrixMode(GL_PROJECTION);
glPopMatrix();			/// pop PROJECTION matrix
glMatrixMode(GL_MODELVIEW);	/// go back to GL_MODELVIEW mode

// Restore GL state
glPopAttrib();
}
*/

/* ISSUE_X64
void GraphicTools::drawCoordinateSystem(float size, std::string label)
{
drawCoordinateSystem(size);
//draw3DText(0.0,10.0, 0.0, GLUT_BITMAP_9_BY_15, label.c_str());
}
*/

void GraphicTools::draw3DText(float x, float y, float z, void *font, const char *string)
{  
	glRasterPos3f(x,y,z);
	for (; *string != '\0'; string++) {
		// ISSUE_X64: glutBitmapCharacter(font, *string);
	}
}

void GraphicTools::drawArrow(float length)
{
	glBegin(GL_LINES);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(length, 0.0f, 0.0f);
	glEnd();
	glDisable(GL_CULL_FACE);
	glBegin(GL_TRIANGLES);
	glVertex3f(length, 0.0f, 0.0f);
	glVertex3f(length - length / (length*3.0f/10.0f), +(length/2.0f) / (length*2.0f/5.0f), 0.0f);
	glVertex3f(length - length / (length*3.0f/10.0f), -(length/2.0f) / (length*2.0f/5.0f), 0.0f);
	glVertex3f(length, 0.0f, 0.0f);
	glVertex3f(length - length / (length*3.0f/10.0f), 0.0f, +(length/2.0f) / (length*2.0f/5.0f));
	glVertex3f(length - length / (length*3.0f/10.0f), 0.0f, -(length/2.0f) / (length*2.0f/5.0f));
	glEnd();
	glBegin(GL_QUADS);
	glVertex3f(length - length / (length*3.0f/10.0f), +(length/2.0f) / (length*2.0f/5.0f), 0.0f);
	glVertex3f(length - length / (length*3.0f/10.0f), 0.0f, +(length/2.0f) / (length*2.0f/5.0f));
	glVertex3f(length - length / (length*3.0f/10.0f), -(length/2.0f) / (length*2.0f/5.0f), 0.0f);
	glVertex3f(length - length / (length*3.0f/10.0f), 0.0f, -(length/2.0f) / (length*2.0f/5.0f));
	glEnd();
}