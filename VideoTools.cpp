#include <stdio.h>
#include <iostream>
#include <limits>
#include <windows.h>
#include "VideoTools.h"


VideoTools::VideoTools()
{
	m_frames = 0;
	m_videoCapture = 0;
}
void VideoTools::getIntensitySingleChannel(IplImage* img, const int &x, const int &y, float &intensity)
{
	intensity = (float)(((unsigned char *)(img->imageData + y*img->widthStep))[x]);
}

void VideoTools::getIntensity(IplImage* img, const int &x, const int &y, float &intensity)
{
	if ( img->nChannels == 3 )
	{
	int R,G,B;
	getRGB(img, x, y, R, G, B);
/*
	CvScalar s;
	s=cvGet2D(img,y,x);

	R = s.val[0];
	G = s.val[1];
	B = s.val[2];
*/
	intensity = 0.3*static_cast<float>(R) + 0.59*static_cast<float>(G) + 0.11*static_cast<float>(B);
	}
	else

	{
		intensity = (float)(((unsigned char *)(img->imageData + y*img->widthStep))[x]);
	}
}

void VideoTools::getRGB(IplImage* image, const int &x, const int &y, int &R, int &G, int &B)
{

	
	//unsigned char* data    = (unsigned char *)image->imageData;
R= ((unsigned char *)(image->imageData + y*image->widthStep))[x*image->nChannels + 0]; // B
G=((unsigned char *)(image->imageData + y*image->widthStep))[x*image->nChannels + 1]; // G
B=((unsigned char *)(image->imageData + y*image->widthStep))[x*image->nChannels + 2]; // R
}


IplImage* VideoTools::uchar2IplImage(unsigned char* data, const int &width, const int &height, const int &channels, const int & depth)
{
	// NOTE: memory must be release using ReleaseImageHeader otherwise there will be a memory leak!!
	IplImage* _img = cvCreateImage( cvSize(width, height), depth, channels );
	//_img->imageData =(char*)malloc(height*width*channels*sizeof(char));
	memcpy(_img->imageData,(char*)(data),height*width*channels);
	cvFlip(_img,_img,0);
	return _img;
}

HBITMAP VideoTools::uchar2HBITMAP(int width,int height, unsigned char *data)
{
  HDC hDC=GetDC(NULL);
  BITMAPINFO bmi;
  HBITMAP ret;
  unsigned char *pImage;
  RGBTRIPLE *pRGBData;
  int WidthInBytes=4*((width*3+3)/4);

  bmi.bmiHeader.biSize=sizeof(BITMAPINFOHEADER);
  bmi.bmiHeader.biWidth=width;
  bmi.bmiHeader.biHeight=height;
  bmi.bmiHeader.biPlanes=1;
  bmi.bmiHeader.biBitCount=24;
  bmi.bmiHeader.biCompression=BI_RGB;
  bmi.bmiHeader.biSizeImage=height*WidthInBytes;
  bmi.bmiHeader.biXPelsPerMeter=0;
  bmi.bmiHeader.biYPelsPerMeter=0;
  bmi.bmiHeader.biClrUsed=0;
  bmi.bmiHeader.biClrImportant=0;

  pImage=new unsigned char[bmi.bmiHeader.biSizeImage];
  
  for(int y=0;y<height;y++){
    pRGBData=(RGBTRIPLE *)&pImage[y*WidthInBytes];
    for(int x=0;x<width;x++){      
      pRGBData[x].rgbtRed=data[y*width+x];
      pRGBData[x].rgbtGreen=data[y*width+x];
      pRGBData[x].rgbtBlue=data[y*width+x];
    }
  }
  ret=CreateDIBitmap(hDC,&bmi.bmiHeader,CBM_INIT,pImage,&bmi,DIB_RGB_COLORS);
  delete[] pImage;
  return ret;
}

IplImage* VideoTools::bitmap2IplImage(HBITMAP hBmp, bool singleChannel)
{
	BITMAP _bmp;

	::GetObject(hBmp,sizeof(BITMAP),&_bmp);

	int _nChannels	= _bmp.bmBitsPixel == 1 ? 1 : _bmp.bmBitsPixel/8 ;
	int _depth		= _bmp.bmBitsPixel == 1 ? IPL_DEPTH_1U : IPL_DEPTH_8U;
	/*IplImage* _img = cvCreateImageHeader( cvSize(_bmp.bmWidth, _bmp.bmHeight), _depth, _nChannels );
	_img->imageData =(char*)malloc(_bmp.bmHeight*_bmp.bmWidth*_nChannels*sizeof(char));
	*/
	IplImage* _img = cvCreateImage(cvSize(_bmp.bmWidth, _bmp.bmHeight), _depth, _nChannels );
	memcpy(_img->imageData,( char*)(_bmp.bmBits),_bmp.bmHeight*_bmp.bmWidth*_nChannels);

	if (_nChannels == 3 && singleChannel )
	{
		IplImage* grayImage = cvCreateImage(cvSize(_bmp.bmWidth, _bmp.bmHeight), IPL_DEPTH_8U, 1);
		cvCvtColor(_img, grayImage, CV_RGB2GRAY);
		cvReleaseImage(&_img);
		
		return grayImage;
	}
	return _img;
} 

bool VideoTools::loadFromFile(std::string fileName)
{
#if USE_DIRECTX
	// The code below requires Windows Platform SDK and DirectX 9.1
	// Quite nasty, therefore it was replaced by OpenCV
	
	if ( m_pGraphBuilder ) m_pGraphBuilder.Release();
	if ( m_pSeek ) m_pSeek.Release();
	if ( m_pMediaControl ) m_pMediaControl.Release();
	if ( m_pMediaEventEx ) m_pMediaEventEx.Release();
	if ( m_pSampleGrabber ) m_pSampleGrabber.Release();

	// Create the graph builder
	HRESULT hr = ::CoCreateInstance(CLSID_FilterGraph, NULL, CLSCTX_INPROC_SERVER, 
								    IID_IGraphBuilder, (void**)&m_pGraphBuilder);
    if (FAILED(hr)) return false;
	//ASSERT(pGraphBuilder != NULL);

	// Create the "Grabber filter"
	CComPtr<IBaseFilter>    pGrabberBaseFilter;
	AM_MEDIA_TYPE   mt;
	hr = ::CoCreateInstance(CLSID_SampleGrabber, NULL, CLSCTX_INPROC_SERVER, 
						    IID_IBaseFilter, (LPVOID *)&pGrabberBaseFilter);

	if (FAILED(hr)) return false;
	pGrabberBaseFilter->QueryInterface(IID_ISampleGrabber, (void**)&m_pSampleGrabber);
	if (m_pSampleGrabber == NULL) return false; //return E_NOINTERFACE;
	
	hr = m_pGraphBuilder->AddFilter(pGrabberBaseFilter,L"Grabber");
	if (FAILED(hr)) return false; //return hr;


	ZeroMemory(&mt, sizeof(AM_MEDIA_TYPE));
	mt.majortype = MEDIATYPE_Video;
	mt.subtype = MEDIASUBTYPE_RGB24;
	mt.formattype = FORMAT_VideoInfo; 
	hr = m_pSampleGrabber->SetMediaType(&mt);        
	if (FAILED(hr)) return false;
		//return hr;
	hr = m_pGraphBuilder->RenderFile(string2widestring(fileName).c_str(),NULL);
	if (FAILED(hr)) return false;
		//return hr;
		
	// QueryInterface for some basic interfaces
    m_pGraphBuilder->QueryInterface(IID_IMediaControl, (void **)&m_pMediaControl);
    m_pGraphBuilder->QueryInterface(IID_IMediaEvent, (void **)&m_pMediaEventEx);

	if (m_pMediaControl == NULL || m_pMediaEventEx == NULL) return false;
		//return E_NOINTERFACE;

	// Set up one-shot mode.
	hr = m_pSampleGrabber->SetBufferSamples(TRUE);
	if (FAILED(hr)) return false;
		//return hr;

	hr = m_pSampleGrabber->SetOneShot(TRUE);
	if (FAILED(hr)) return false;
		//return hr;

	m_pSeek = m_pMediaControl;
	if (m_pSeek == NULL) return false;
		//return E_NOINTERFACE;


	double rate;
	hr = m_pSeek->GetRate(&rate);
	if (FAILED(hr)) return false;
		//return hr;

	const GUID timeFormat = TIME_FORMAT_FRAME;
	hr = m_pSeek->SetTimeFormat(&timeFormat);


	LONGLONG Duration;
	hr = m_pSeek->GetDuration(&Duration);
	if (FAILED(hr)) return false;
		//return hr;
	m_frames = int(Duration);

	return true;
#else
	 if ( m_videoCapture )
		cvReleaseCapture(&m_videoCapture);


	m_videoCapture = cvCaptureFromAVI(fileName.c_str());

	m_frames = (int) cvGetCaptureProperty(m_videoCapture,  CV_CAP_PROP_FRAME_COUNT);


	if ( m_videoCapture && m_frames > 0)
		return true;
	else
		return false;
#endif
}

VideoTools::~VideoTools()
{

#if USE_DIRECTX
	if ( m_pGraphBuilder ) m_pGraphBuilder.Release();
	if ( m_pSeek ) m_pSeek.Release();
	if ( m_pMediaControl ) m_pMediaControl.Release();
	if ( m_pMediaEventEx ) m_pMediaEventEx.Release();
	if ( m_pSampleGrabber ) m_pSampleGrabber.Release();
#else
	if ( m_videoCapture )
		cvReleaseCapture(&m_videoCapture);
#endif
}
bool VideoTools::getFrame(int frameNumber, HBITMAP& bmp, int &Width, int &Height)
{
#if USE_DIRECTX
	HRESULT hr;

	if ( frameNumber < 0 || frameNumber > m_frames ) return false;

	REFERENCE_TIME rtStart = frameNumber;
	

	REFERENCE_TIME rtStop = rtStart; 
			
	hr = m_pSeek->SetPositions(&rtStart, AM_SEEKING_AbsolutePositioning, 
							 &rtStop, AM_SEEKING_AbsolutePositioning);
	if (FAILED(hr)) return false;
		//return hr;

	CComQIPtr<IVideoWindow> pVideoWindow = m_pGraphBuilder;

	hr = pVideoWindow->put_AutoShow(OAFALSE);
	if (FAILED(hr)) return false;
		//return hr;

	// Run the graph and wait for completion.
	hr = m_pMediaControl->Run();
	if (FAILED(hr)) return false;
		//return hr;

	long evCode;
	hr = m_pMediaEventEx->WaitForCompletion(INFINITE, &evCode);
	if (FAILED(hr)) return false;
		//return hr;
				
	AM_MEDIA_TYPE MediaType;
	ZeroMemory(&MediaType,sizeof(MediaType));
	hr = m_pSampleGrabber->GetConnectedMediaType(&MediaType); 
	if (FAILED(hr)) return false;
		//return hr;

	// Get a pointer to the video header.
	VIDEOINFOHEADER *pVideoHeader = (VIDEOINFOHEADER*)MediaType.pbFormat;
	if (pVideoHeader == NULL) return false;
		//return E_FAIL;

	// The video header contains the bitmap information. 
	// Copy it into a BITMAPINFO structure.
	BITMAPINFO BitmapInfo;
	ZeroMemory(&BitmapInfo, sizeof(BitmapInfo));
	CopyMemory(&BitmapInfo.bmiHeader, &(pVideoHeader->bmiHeader), sizeof(BITMAPINFOHEADER));

	// Create a DIB from the bitmap header, and get a pointer to the buffer.
	void *buffer = NULL;
	HBITMAP hBitmap = ::CreateDIBSection(0, &BitmapInfo, DIB_RGB_COLORS, &buffer, NULL, 0);
	GdiFlush();
	// Copy the image into the buffer.
	long size = BitmapInfo.bmiHeader.biSizeImage;
	hr = m_pSampleGrabber->GetCurrentBuffer(&size, (long *)buffer);
	if (FAILED(hr)) return false;
		//return hr;
	Width = pVideoHeader->bmiHeader.biWidth;
	Height = pVideoHeader->bmiHeader.biHeight;
	
	bmp = hBitmap;
	return true;
#else

	int ret = cvSetCaptureProperty(m_videoCapture, CV_CAP_PROP_POS_FRAMES, frameNumber);
	IplImage* _image = cvQueryFrame(m_videoCapture);
	
	bmp = IplImage2HBITMAP(_image);
	Width = _image->width;
	Height = _image->height;
	//cvReleaseImage(&_image);
	return true;
#endif
}

int VideoTools::getNumberOfFrames()
{
	return m_frames;
}

std::wstring VideoTools::string2widestring(const std::string& s)
{
    int len;
    int slength = (int)s.length() + 1;
    len = MultiByteToWideChar(CP_ACP, 0, s.c_str(), slength, 0, 0); 
    wchar_t* buf = new wchar_t[len];
    MultiByteToWideChar(CP_ACP, 0, s.c_str(), slength, buf, len);
    std::wstring r(buf);
    delete[] buf;
    return r;
}

float VideoTools::getIntensity(HBITMAP hBmp, int posX, int posY)
{
	BITMAP    bm;
	HDC       hDC, hMemDC;

	hDC = GetDC(NULL);
	hMemDC = CreateCompatibleDC(hDC);
	GetObject(hBmp, sizeof(BITMAP), (LPSTR)&bm);
	SelectObject(hMemDC, hBmp);

	// determine the color of each pixel
	COLORREF colorref = GetPixel(hMemDC, posX, posY);

	int R = GetRValue(colorref);
	int G = GetGValue(colorref);
	int B = GetBValue(colorref);

	return 0.3*R + 0.59*G + 0.11*B;

	DeleteDC(hMemDC);
	ReleaseDC(NULL, hDC);
}

bool VideoTools::getImageMinMax(IplImage *img, int &minX, int &maxX, int &minY, int &maxY)
{
	bool _found = false;
	minX = std::numeric_limits<int>::max();
	minY = std::numeric_limits<int>::max();

	maxX = -std::numeric_limits<int>::max();
	maxY = -std::numeric_limits<int>::max();
	IplImage *_img = 0;


	if (img->nChannels == 3)
	{
		_img = cvCreateImage(cvSize(_img->width, _img->height), IPL_DEPTH_8U, 1);
		cvCvtColor(img, _img, CV_RGB2GRAY);
	}
	else
	{
		_img = img;
	}


	for(int y = 0; y < _img->height; y++)
	{
		for(int x = 0; x < _img->width; x++)
		{
			

			// if the pixel is segmented, consider its position
			if ( (float)(((unsigned char *)(_img->imageData + y*_img->widthStep))[x]) > 0.0)
			{
				_found = true;
				if ( x > maxX )
					maxX = x;
				if ( x < minX )
					minX = x;

				if ( y > maxY )
					maxY = y;
				if ( y < minY )
					minY = y;
			}
		}
	}



	if (img->nChannels == 3 )
		cvReleaseImage(&_img);


	if (!_found)
		return false;

	if ( minY >= 0 && minX >= 0 && maxX < _img->width && maxY < _img->height)
		return true;
	else
		return false;
}

bool VideoTools::getImageMinMax(HBITMAP hBmp, int &minX, int &maxX, int &minY, int &maxY)
{
	bool _found = false;
	minX = std::numeric_limits<int>::max();
	minY = std::numeric_limits<int>::max();

	maxX = -std::numeric_limits<int>::max();
	maxY = -std::numeric_limits<int>::max();

	BITMAP    bm;
	HDC       hDC, hMemDC;
	COLORREF  colorref;
	
	hDC = GetDC(NULL);
	hMemDC = CreateCompatibleDC(hDC);
	GetObject(hBmp, sizeof(BITMAP), (LPSTR)&bm);
	SelectObject(hMemDC, hBmp);

	for(int y = 0; y < bm.bmHeight; y++)
	{
		for(int x = 0; x < bm.bmWidth; x++)
		{
			// determine the color of each pixel
			colorref = GetPixel(hMemDC, x, y);

			int colR = GetRValue(colorref);
			int colG = GetGValue(colorref);
			int colB = GetBValue(colorref);

			// if the pixel is segmented, consider its position
			if ( colorref != RGB(0,0,0))
			{
				_found = true;
				if ( x > maxX )
					maxX = x;
				if ( x < minX )
					minX = x;

				if ( y > maxY )
					maxY = y;
				if ( y < minY )
					minY = y;
			}
		}
	}
DeleteDC(hMemDC);
	ReleaseDC(NULL, hDC);

	if (!_found)
		return false;

	if ( minY >= 0 && minX >= 0 && maxX < bm.bmWidth && maxY < bm.bmHeight)
		return true;
	else
		return false;
}



HBITMAP VideoTools::IplImage2HBITMAP(const IplImage *Image)
{
        int bpp = Image->nChannels * 8;
        assert(Image->width >= 0 && Image->height >= 0 &&
                (bpp == 8 || bpp == 24 || bpp == 32));
        CvMat dst;
        void* dst_ptr = 0;
        HBITMAP hbmp = NULL;
        unsigned char buffer[sizeof(BITMAPINFO) + 255*sizeof(RGBQUAD)];
        BITMAPINFO* bmi = (BITMAPINFO*)buffer;
        BITMAPINFOHEADER* bmih = &(bmi->bmiHeader);
        
        ZeroMemory(bmih, sizeof(BITMAPINFOHEADER));
        bmih->biSize = sizeof(BITMAPINFOHEADER);
        bmih->biWidth = Image->width;
        bmih->biHeight = Image->origin ? abs(Image->height) :
        -abs(Image->height);
        bmih->biPlanes = 1;
        bmih->biBitCount = bpp;
        bmih->biCompression = BI_RGB;
        
        if (bpp == 8) {
                RGBQUAD* palette = bmi->bmiColors;
                int i;
                for (i = 0; i < 256; i++) {
                        palette[i].rgbRed = palette[i].rgbGreen = palette[i].rgbBlue
                                = (BYTE)i;
                        palette[i].rgbReserved = 0;
                }
        }
        
        hbmp = CreateDIBSection(NULL, bmi, DIB_RGB_COLORS, &dst_ptr, 0, 0);
        cvInitMatHeader(&dst, Image->height, Image->width, CV_8UC3,
                dst_ptr, (Image->width * Image->nChannels + 3) & -4);
        cvConvertImage(Image, &dst, Image->origin ? CV_CVTIMG_FLIP : 0);
        
        return hbmp;
}