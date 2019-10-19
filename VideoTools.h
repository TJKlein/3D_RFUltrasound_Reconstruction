#ifndef VIDEOTOOLS_H
#define VIDEOTOOLS_H

#define USE_DIRECTX 0

#include <windows.h>
#if USE_DIRECTX
#include "AtlBase.h"
#include "dShow.h"
#include "Qedit.h"
#endif


#include <opencv/cv.h>
#include <opencv/highgui.h>
//#include <opencv2/core/core.hpp>
//#include <opencv2/imgproc/imgproc.hpp>
//#include <opencv2/highgui/highgui.hpp>

class VideoTools
{
	public:
		VideoTools();
		~VideoTools();
		bool	loadFromFile(std::string fileName);
		int		getNumberOfFrames();
		bool	getFrame(int frameNumber, HBITMAP& bmp, int &Width, int &Height);
		bool	getImageMinMax(HBITMAP hBmp, int &minX, int &maxX, int &minY, int &maxY);
		bool	getImageMinMax(IplImage *img, int &minX, int &maxX, int &minY, int &maxY);
		float	getIntensity(HBITMAP hBmp, int posX, int posY);
		HBITMAP VideoTools::IplImage2HBITMAP(const IplImage *Image);
		IplImage* bitmap2IplImage(HBITMAP hBmp, bool singleChannel = false);
		IplImage* uchar2IplImage(unsigned char* data, const int &width, const int &height, const int &channels, const int & depth);
		HBITMAP uchar2HBITMAP(int width,int height, unsigned char *data);
		void	getRGB(IplImage* image, const int &x, const int &y, int &R, int &G, int &B);
		void	getIntensity(IplImage* image, const int &x, const int &y, float &intensity);
		void	getIntensitySingleChannel(IplImage* image, const int &x, const int &y, float &intensity);
	private:
		// variables necessary for accessing video data
		
		int		m_frames;
#if USE_DIRECTX
		CComPtr<IGraphBuilder>	m_pGraphBuilder;
		CComQIPtr<IMediaSeeking> m_pSeek;
		CComPtr<IMediaControl>	m_pMediaControl;
		CComPtr<IMediaEvent>	m_pMediaEventEx;
		CComPtr<ISampleGrabber> m_pSampleGrabber;
#else
		CvCapture* m_videoCapture;
#endif
		std::wstring string2widestring(const std::string& s);
}; 

#endif