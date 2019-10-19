#include "RFCompounding.h"
#include <iostream>
#include <fftw3.h>

std::wstring s2ws(const std::string& s)
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


	RFFileStream::RFFileStream(STREAMING_MODE mode) :
		m_evenlopeSignal(0), // assume unfiltered signal for now
		m_currentFrame(0),
		m_rfDataFile(0)
	{
		m_mode = mode;
		//m_file.open("info.txt");
	}

		RFFileStream::~RFFileStream()
		{
			//m_file.close();
		}





	bool RFFileStream::openFile(std::string filename, DATA_TYPE dataType)
	{	
		std::string _path = filename.substr( 0, filename.find_last_of( '/' ) +1 );

		// open the configuration file and read out the information specifying the image recording settings
		m_xmlTools.loadFile(filename);

		std::string _tmpString;


		m_xmlTools.handleData(false, _tmpString, std::string(VERSION_TAG));
		m_version = atof(_tmpString.c_str());
		m_xmlTools.handleData(false, m_dataFileName, std::string(IMAGEDATAFILENAME_TAG));
		m_xmlTools.handleData(false, m_timeStampFileName, std::string(TIMESTAMPSFILENAME_TAG));
		m_xmlTools.handleData(false, _tmpString, std::string(NUMBER_OF_FRAMES_TAG));
		m_recordingInfo.frames = atoi(_tmpString.c_str());
		m_xmlTools.handleData(false, _tmpString, std::string(SAMPLES_SCANLINE_TAG));
		m_recordingInfo.height = atoi(_tmpString.c_str());
		m_xmlTools.handleData(false, _tmpString, std::string(NUMBER_OF_SCANLINES_TAG));
		m_recordingInfo.width = atoi(_tmpString.c_str());

		m_xmlTools.handleData(false, _tmpString, std::string(SAMPLING_FREQUENCY_TAG));
		m_recordingInfo.samplingFrequency = atoi(_tmpString.c_str());
		m_xmlTools.handleData(false, _tmpString, std::string(PENETRATION_DEPTH_TAG));
		m_recordingInfo.penetrationDepth = atoi(_tmpString.c_str());
		m_xmlTools.handleData(false, _tmpString, std::string(RECORDING_FREQUENCY_TAG));
		m_recordingInfo.recordingFrequency = atoi(_tmpString.c_str());
		m_xmlTools.handleData(false, _tmpString, std::string(PROBE_NAME_TAG));
		m_recordingInfo.probeName = _tmpString;
		//m_recordingInfo.pro = atoi(_tmpString.c_str());


		// determine if the signal was envelope detected
		m_evenlopeSignal = false;
		int _filterOperations = m_xmlTools.numberOfElements(std::string(FILTER_OPERATION_TAG));
		for (int i=0;i<_filterOperations;i++)
		{
			std::string _type;
			m_xmlTools.handleData(false, _type, std::string(FILTER_OPERATION_TAG), std::string(FILTER_TYPE_ATTRIBUTE), i);
			if ( _type.compare(std::string(FILTER_ENVELOPE)) == 0 )
			{
				m_evenlopeSignal = true;
			}
		}

		// now comppute the size of a frame, then we can jump from frame to frame (when the signal is filtered it is 32 bit float, otherwise 16 bit signed short
		switch(dataType)
		{
		case DATA_TYPE::FLOAT32: {
			m_rfFrameSize = sizeof(float)*m_recordingInfo.width*m_recordingInfo.height; break; }
		case DATA_TYPE::USHORT16: {
			m_rfFrameSize = sizeof(char)*2*m_recordingInfo.width*m_recordingInfo.height; break; }
		case DATA_TYPE::FLOAT5X32: {
			m_rfFrameSize = sizeof(float)*5*m_recordingInfo.width*m_recordingInfo.height; break; }
		}
		if ( m_rfDataFile == NULL )
		{
			std::string _filename = (_path+m_dataFileName);
			
			SetLastError(0);
			m_rfDataFile = CreateFileA( _filename.c_str(), GENERIC_READ ,  FILE_SHARE_READ, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL,  NULL);

			 if (m_rfDataFile == INVALID_HANDLE_VALUE) 
			{ 
				std::cout << "Error opening file: " << _filename << std::endl;

				return false;
			 }
 
			//m_rfDataFile.open(_path+m_dataFileName,std::ios::binary);
			/*m_rfDataFile.seekg (0, std::ios::end);
			m_rfDataFileLength = m_rfDataFile.tellg();
			m_rfDataFile.seekg (0, std::ios::beg);
			*/
			 m_rfDataFileLength = (int)GetFileSize(m_rfDataFile, NULL);

			 std::cout << "Size: " << m_rfDataFileLength << " / " << GetLastError() << std::endl;
			/*if ( m_rfDataFile.fail() )
				std::cout << "Seek error!" << std::endl;
				*/
			// now check if the file is okay, or if it is corrupted
			if ( (m_rfDataFileLength/m_rfFrameSize) != m_recordingInfo.frames )
				return false;
			/*	*/
			m_currentFrame = 0;
			return true;
		}
		else
			return false;
	}


	bool RFFileStream::writeFile()
	{
		if ( m_rfDataFile )
		{
			return true;
		}
		else
			return false;
	}

	void RFFileStream::closeFile()
	{
		//m_rfDataFile.close();
		CloseHandle(m_rfDataFile);
	}

	
	int RFFileStream::currentFrameIndex()
	{
		return m_currentFrame;
	}
	

	void* RFFileStream::getFrameByIndex(int index)
	{
		LARGE_INTEGER _pint;
		
		char *_data = new char[m_rfFrameSize];
		if ( index == m_currentFrame)
		{
			// no need to seek
			//std::cout << "No need to seek!" << std::endl;
		}
	/*	if ( index == m_currentFrame+1)
		{
			//m_rfDataFile.seekg(m_rfFrameSize, std::ios_base::cur);
		}
		*/else //if ( index > m_currentFrame+1)
		{
			//m_rfDataFile.seekg( (index-m_currentFrame)*m_rfFrameSize, std::ios_base::cur);
			
			LARGE_INTEGER _filePos;
			_filePos.QuadPart =  static_cast<__int64>(index-m_currentFrame)*static_cast<__int64>(m_rfFrameSize);
			if ( SetFilePointerEx(m_rfDataFile, _filePos, &_pint, FILE_CURRENT) == 0)
				std::cout << "Error seeking file!" << std::endl;

			//std::cout << "Shift: " << _filePos.QuadPart << " / " << GetLastError()  << " / " << _pint.QuadPart << std::endl;
		}
	/*	else // slow
		{
			
			//m_rfDataFile.seekg(m_rfFrameSize*index, std::ios_base::beg);
			
			

			LARGE_INTEGER _filePos;
			_filePos.QuadPart =  static_cast<__int64>(m_rfFrameSize*index)*static_cast<__int64>(m_rfFrameSize);
			if ( SetFilePointerEx(m_rfDataFile, _filePos, NULL, FILE_BEGIN) == 0)
				std::cout << "Error seeking file!" << std::endl;
		}*/

		//if ( m_rfDataFile.fail() )
		//	std::cout << "Seek error!" << std::endl;

		//m_file << m_rfDataFile.tellg() << std::endl;
		//m_file.flush();

		switch(m_mode)
		{
		case STREAMING_MODE::NORMAL: {
					//m_rfDataFile.read(&_data[0], m_rfFrameSize);
					DWORD  dwBytesRead = 0;

					if (! ReadFile(m_rfDataFile, &_data[0], m_rfFrameSize, &dwBytesRead, NULL)  )

						std::cout << "Error reading data! " << GetLastError() << std::endl;

					

					 /*std::cout << "Read: "  <<  dwBytesRead << " of " << m_rfFrameSize << std::endl;
					 if ( dwBytesRead == 0)
					 {

						 LPVOID lpMsgBuf;
  FormatMessage(
      FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM,
      0,
      GetLastError(),
      MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), // Default language
      (LPTSTR) &lpMsgBuf,
      0,
      0
  );
  //Transform to std::string
  const std::string lastError = (char*)lpMsgBuf; //A dirty cast!
  // Free the buffer.
  LocalFree( lpMsgBuf );


						 std::cout << "Last Error: " << GetLastError() << std::endl;
					 }*/
					break; }
									
		case STREAMING_MODE::MISSING4BYTE: {
			/*char *_ptr = &_data[4];
			_data[0] = 0; // first four bytes set to zero
			_data[1] = 0;
			_data[2] = 0;
			_data[3] = 0;
			m_rfDataFile.read(_ptr, m_rfFrameSize-4*sizeof(char));
			// now read the remaining 4 bytes to stay in line
			char _dummy[4];
			m_rfDataFile.read(&_dummy[0], 4*sizeof(char));*/
			
			break; }
		case STREAMING_MODE::UNNECESSARY4BYTEHEADER: { 
			// now read the unnecessary leading 4 bytes to stay in line
			/*char _dummy[4];
			m_rfDataFile.read(&_dummy[0], 4*sizeof(char));
			m_rfDataFile.read(_data, m_rfFrameSize);
			break;*/
													 }
		}
		/*std::fstream _rfDataFile;
		std::stringstream _ss;
		_ss << "crap" << index << ".rf";
		_rfDataFile.open(_ss.str(),std::ios::out   | std::ios::binary);
		_rfDataFile.write(&_data[0],m_rfFrameSize);
		_rfDataFile.close();
		
			*/
		m_currentFrame = index+1;

		/*if ( m_rfDataFile.fail() )
			std::cout << "Read error!" << std::endl;

		if ( m_rfDataFile.eof() )
		{
			// now go back to the beginning
			m_rfDataFile.clear();
			m_currentFrame = 0;
			m_rfDataFile.seekg(0,  std::ios_base::beg);
			std::cout << "End of file reached!!" << std::endl;
		}*/
		return reinterpret_cast<void*>(_data);
	}

	
	void* RFFileStream::getNextFrame()
	{
		char *_data = new char[m_rfFrameSize];

		m_currentFrame++;
		//m_rfDataFile.seekg(m_rfFrameSize, std::ios_base::cur);
		switch(m_mode)
		{
		case STREAMING_MODE::NORMAL: {
					//m_rfDataFile.read(_data, m_rfFrameSize);
				DWORD  dwBytesRead = 0;	
					if ( ReadFile(m_rfDataFile, &_data[0], m_rfFrameSize, &dwBytesRead, NULL) == 0)
						std::cout << "Error reading data!" << std::endl;
					break; }
		case STREAMING_MODE::MISSING4BYTE: {
			/*char *_ptr = &_data[4];
			_data[0] = 0; // first four bytes set to zero
			_data[1] = 0;
			_data[2] = 0;
			_data[3] = 0;
			m_rfDataFile.read(_ptr, m_rfFrameSize-4*sizeof(char));
			// now read the remaining 4 bytes to stay in line
			char _dummy[4];
			m_rfDataFile.read(&_dummy[0], 4*sizeof(char));
			break;*/ }
		case STREAMING_MODE::UNNECESSARY4BYTEHEADER: { 
			// now read the unnecessary leading 4 bytes to stay in line
			/*char _dummy[4];
			m_rfDataFile.read(&_dummy[0], 4*sizeof(char));
			m_rfDataFile.read(_data, m_rfFrameSize);
			break;*/}
		}

		/*if ( m_rfDataFile.fail() )
			std::cout << "Read error!" << std::endl;
			*/
		return reinterpret_cast<void*>(_data);
	}

template< typename T >
RFProcessing<T>::RFProcessing(int width, int height, float lowFreq, float highFreq, float samplingFreq, int filterOrder )
{
	m_width = width;
	m_height = height;

	// normalize the frequencies [0,1]
	float _normalizedHighFreq = highFreq / (samplingFreq / 2.0 );
	float _normalizedLowFreq = lowFreq / (samplingFreq / 2.0 );

	// now compute the butterworth filter coefficients

	/* calculate the d coefficients */
    m_coefficientsD = dcof_bwbp( filterOrder, _normalizedLowFreq, _normalizedHighFreq );
	m_coefficientsC = ccof_bwbp( filterOrder );
	m_scalingFactor = sf_bwbp( filterOrder, _normalizedLowFreq, _normalizedHighFreq ); 

	m_filterOrder = filterOrder;
}

template< typename T >
RFData<T>* RFProcessing<T>::convert(char* data, int size)
{
	int16 *_data16 = reinterpret_cast<int16*>(data);

	int _currentPos = 0 ;
	int _elements = size/(sizeof(char)*2);

	float *_rfData = new float[_elements];
	float *_fp = &_rfData[0];


	while(_currentPos < _elements)
	{
		short _value = _data16->value;
		//std::cout << _value << std::endl;;
		(*_fp++) = static_cast<float>(_value);
		_currentPos++;
		_data16++;
	};

	RFData<T> *_result = new RFData<T>();
	_result->setData(_rfData, m_width, m_height);
	return _result;
}

template< typename T >
RFData<T>* RFProcessing<T>::convert( short *data, int size)
{
	//int16 *_data16 = reinterpret_cast<int16*>(data);
	short *_data16 = &data[0];
	int _currentPos = 0 ;
	int _elements = size/(sizeof(char)*2);

	float *_rfData = new float[_elements];
	float *_fp = &_rfData[0];


	while(_currentPos < _elements)
	{
		short _value = *_data16; // _data16->value;
		//std::cout << _value << std::endl;;
		(*_fp++) = static_cast<float>(_value);
		_currentPos++;
		_data16++;
	};

	RFData<T> *_result = new RFData<T>();
	_result->setData(_rfData, m_width, m_height);
	return _result;
}

template< typename T >
RFData<T>* RFProcessing<T>::bandpassFilter(RFData<T> *rf)
{
	// this is where the result goes
	double *_data = new double[m_width * m_height];
	float *_dataF = new float[m_width * m_height];
	memset(_data,0.0,sizeof(double)*m_height*m_width);
	for(int i=0; i < m_width; i++)
	{
		float *_rawBeam = &rf->getData()[i*m_height];
		double *_rawBeamD = new double[m_height];
		for (int j=0; j < m_height; j++)
		{
			_rawBeamD[j] = static_cast<double>(_rawBeam[j]);
		}
		double *_filteredBeam  = &_data[i*m_height];

		for (int j=0; j < m_height; j++)
		{
			filterBeam(_filteredBeam, _rawBeamD, j);
		}
		delete[] _rawBeamD;
	}
	// now convert back to float
		for(int i=0; i < m_width * m_height; i++)
		{
			_dataF[i] = static_cast<float>(_data[i]);
		}
		delete _data;
	RFData<T> *_result = new RFData<T>();
	_result->setData(_dataF, m_width, m_height);
	return _result;
}

template< typename T >
void RFProcessing<T>::filterBeam(double *output, double *input, int pos)
{
	double _resultA = 0.0f;
	double _resultB = 0.0f;
	for(int i=0;i<=2*m_filterOrder;i++)
	{
		int _currentPos = pos - i;

		// positions outside the data range are assumed to be zero
		double _x = 0.0f;
		
		if ( _currentPos >= 0 )
		{
			_x = input[_currentPos];
		}
		//else _x = input[-(_currentPos)];

		// mix the current signal with the previous signal in combinations with the filter coefficients
		_resultB += _x * (m_coefficientsC[i] *m_scalingFactor);
		
		//std::cout << i << ". " << (m_coefficientsC[i] *m_scalingFactor) << std::endl;
	}

	for(int i=1;i<=2*m_filterOrder;i++)
	{
		double _y = 0.0f;
		int _currentPos = pos - i;

		if ( _currentPos >= 0 )
		{
			_y = output[_currentPos];
		}
		
		
		_resultA += _y *  (m_coefficientsD[i]);
	//	std::cout << i << ". " << (m_coefficientsD[i]) << std::endl;
		
	}

	//std::cout << "Out: " << _resultA  << "\t" << _resultB << std::endl;
	output[pos] = _resultB - _resultA;
}

template< typename T >
float* RFProcessing<T>::envelopeDetection(float *data)
{
	// each beam is FFT transformed
	fftwf_complex *_out;
	_out = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * m_height);

	fftwf_complex *_hilbert;
	_hilbert = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * m_height);

	float *_coefficients = new float[m_height];
	float *_result = new float[m_height*m_width];


	// N/2 equals: m_height/2+1
	
	//the first and the N/2+1 elements are filled with one
	_coefficients[0] = 1.0;
	_coefficients[m_height/2] = 1.0;

	// fill starting from 2nd to the N/2-1 with 2
	for(int i=1;i<=m_height/2-1;i++)
		_coefficients[i]=2.0;
	// fill with zero N+1 to last element
	memset(&_coefficients[m_height/2+1],0,sizeof(float)*(m_height/2-1));

	// normalizing coefficient as IFFT is unnormalized (unfortunately!!)
	float _coeff = 1.0/(float)m_height;
	float *_tmpPtr = &_result[0];
	for (int i=0; i<m_width;i++)
	{
		float *_beam = &data[i*m_height];
		float *_resultPtr = &_result[i*m_height];

		// Forward DFT
		fftwf_plan _planForward = fftwf_plan_dft_r2c_1d(m_height, _beam, _out, FFTW_ESTIMATE);
		fftwf_execute(_planForward); 
        fftwf_destroy_plan(_planForward);

		{
		// now we have to multiply the 'complex' fourier coefficients with the _coeffcients (assuming imaginary part 0)
		float *_coeffPtr = _out[0];
		for (int j=0; j<m_height;j++)
		{
			(*_coeffPtr++)*=_coefficients[j];
			(*_coeffPtr++)*=_coefficients[j];
		}
		}

		// Backward (inverse) FFT

		fftwf_plan _planBackward = fftwf_plan_dft_1d(m_height, _out, _hilbert, 1, FFTW_ESTIMATE);
		fftwf_execute(_planBackward); 
        fftwf_destroy_plan(_planBackward);
		{
		// take the absolute value of hilbert to become REAL again (puuhh...computer scientists don't know how to handle complex values...abs comes to rescue ;-) )
		float *_coeffPtr = _hilbert[0];
		for (int j=0; j<m_height;j++)
		{
			// now normalize the result by 1/N 
			float _a = pow((*_coeffPtr++)*_coeff,2.0f);
			float _b = pow((*_coeffPtr++)*_coeff,2.0f);

			// back to the safe haven of real values!
			(*_tmpPtr++) = sqrt(_a+_b);
		}
		}
	}

	fftwf_complex *_test;
	_test = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * m_height);


	 fftwf_free(_out);

	 return _result;
}

template< typename T >
RFProcessing<T>::~RFProcessing()
{
	// free the butterworth filter coefficients
	if ( m_coefficientsD != NULL ) free( m_coefficientsD );
    if ( m_coefficientsC != NULL ) free( m_coefficientsC );

}


template< typename T >
ColorRF<T>::ColorRF()
{
}

template< typename T >
RFData<T>::RFData() : 
m_data(0)
{
}

template< typename T >
StandardRF<T>::StandardRF()
{
}

template< typename T >
void StandardRF<T>::computeScanlines(std::fstream &file, std::vector<Scanline<T>*> &scanlineVector,std::vector<FrameScanlineRange> &fsr, unsigned char fileIndex, unsigned int &scanlineIndex, int imageIndex, std::vector<UltrasoundSettings> *ultrasoundSettings, unsigned char settingsIndex)
{
	Eigen::Matrix4f _corner_T_apex;

	_corner_T_apex << 1.0, 0.0, 0.0, (*ultrasoundSettings)[static_cast<int>(settingsIndex)].apexX * (*ultrasoundSettings)[static_cast<int>(settingsIndex)].scaleX,  0.0, 1.0, 0.0, (*ultrasoundSettings)[static_cast<int>(settingsIndex)].apexY * (*ultrasoundSettings)[static_cast<int>(settingsIndex)].scaleY,  0.0, 0.0, 1.0, 0.0,  0.0, 0.0, 0.0, 1.0;
	
	// m_transformation = tracking data * calibration matrix
	Eigen::Matrix4f _world_T_apex = m_transformation * _corner_T_apex;

	Eigen::Vector4f _world_P_apex = _world_T_apex * Eigen::Vector4f(0.0 ,0.0, 0.0, 1.0);

	if ((*ultrasoundSettings)[static_cast<int>(settingsIndex)].probeWidth > 0 ) // linear probe
	{

		float _widthIncrement = (float)(*ultrasoundSettings)[static_cast<int>(settingsIndex)].probeWidth/static_cast<float>((*ultrasoundSettings)[static_cast<int>(settingsIndex)].scanlines-1);

		FrameScanlineRange _fsr;
	_fsr.start = scanlineVector.size()-1;
	_fsr.end = _fsr.start + (*ultrasoundSettings)[static_cast<int>(settingsIndex)].scanlines;
	
	// scanlines are supposed to be a power of 2
	fsr.push_back(_fsr);

	UltrasoundSettings _tmp = (*ultrasoundSettings)[static_cast<int>(settingsIndex)];

	float _radius = (float)(*ultrasoundSettings)[static_cast<int>(settingsIndex)].probeRadius;
	float _depth = (float)(*ultrasoundSettings)[static_cast<int>(settingsIndex)].penetrationDepth;

	for (int i=0; i<(*ultrasoundSettings)[static_cast<int>(settingsIndex)].scanlines;i++)
	{
		float _posX = i*_widthIncrement;

		
		Eigen::Vector4f _scanlineStartPos;
		
		
			Eigen::Vector4f _tempVec = Eigen::Vector4f( _posX, 0.0,  0.0, 1.0);
			_scanlineStartPos = _world_T_apex * _tempVec;
		

		
		 _tempVec = Eigen::Vector4f(_posX, _depth,  0.0, 1.0);
		Eigen::Vector4f _scanlineEndPos = _world_T_apex * _tempVec;

		Eigen::Vector4f _beamDirection = _scanlineEndPos - _scanlineStartPos;

		Eigen::Vector4f _beamDirN = _beamDirection.normalized();

		// now insert the beam direction into the beam collection
#ifdef USE_USHORT16

		Scanline<unsigned short> *_scanline = new Scanline<unsigned short>(Eigen::Vector3f(_scanlineStartPos.x(), _scanlineStartPos.y(), _scanlineStartPos.z()), Eigen::Vector3f(_scanlineEndPos.x(), _scanlineEndPos.y(), _scanlineEndPos.z()),ultrasoundSettings, settingsIndex, scanlineIndex++, fileIndex, imageIndex, i);
		
#endif

#ifdef USE_FLOAT32 || USE_FLOAT5X32
		Scanline<float> *_scanline = new Scanline<float>(Eigen::Vector3f(_scanlineStartPos.x(), _scanlineStartPos.y(), _scanlineStartPos.z()), Eigen::Vector3f(_scanlineEndPos.x(), _scanlineEndPos.y(), _scanlineEndPos.z()),ultrasoundSettings, settingsIndex, scanlineIndex++, fileIndex, imageIndex, i);
#endif
		scanlineVector.push_back(_scanline);

		// now write the direction to the file
		float32 _x, _y, _z;
		_x.value = (_scanlineStartPos/*+_beamDirN*/).x();
		_y.value = (_scanlineStartPos/*+_beamDirN*/).y();
		_z.value = (_scanlineStartPos/*+_beamDirN*/).z();
		file.write(&_x.byte.c[0], sizeof(float32));
		file.write(&_y.byte.c[0], sizeof(float32));
		file.write(&_z.byte.c[0], sizeof(float32));

		_x.value = (_scanlineEndPos/*+_beamDirN*/).x();
		_y.value = (_scanlineEndPos/*+_beamDirN*/).y();
		_z.value = (_scanlineEndPos/*+_beamDirN*/).z();
		file.write(&_x.byte.c[0], sizeof(float32));
		file.write(&_y.byte.c[0], sizeof(float32));
		file.write(&_z.byte.c[0], sizeof(float32));

		//_debugFile << /*scanlineVector.size() << ". " <<*/ _scanlineStartPos.x() << "\t" << _scanlineStartPos.y() << "\t" << _scanlineStartPos.z() << "\t" << _scanlineEndPos.x() << "\t" << _scanlineEndPos.y() << "\t" << _scanlineEndPos.z() << std::endl;
	}





	}
	else //curvilinear or phased-array probe
	{
	float _angleIncrement = M_PI/180.0*(*ultrasoundSettings)[static_cast<int>(settingsIndex)].probeAngle/((float) (*ultrasoundSettings)[static_cast<int>(settingsIndex)].scanlines - 1.0);

	//std::ofstream _debugFile;
	//_debugFile.open("scanlines.txt", std::ios_base::app);

	
	FrameScanlineRange _fsr;
	_fsr.start = scanlineVector.size()-1;
	_fsr.end = _fsr.start + (*ultrasoundSettings)[static_cast<int>(settingsIndex)].scanlines;
	
	// scanlines are supposed to be a power of 2
	fsr.push_back(_fsr);

	UltrasoundSettings _tmp = (*ultrasoundSettings)[static_cast<int>(settingsIndex)];

	float _radius = (float)(*ultrasoundSettings)[static_cast<int>(settingsIndex)].probeRadius;
	float _depth = (float)(*ultrasoundSettings)[static_cast<int>(settingsIndex)].penetrationDepth;

	for (int i=0; i<(*ultrasoundSettings)[static_cast<int>(settingsIndex)].scanlines;i++)
	{
		float _theta = -(*ultrasoundSettings)[static_cast<int>(settingsIndex)].probeAngle/2.0*M_PI/180.0 + (float)i * _angleIncrement;

		
		Eigen::Vector4f _scanlineStartPos;
		
		// phased array assumption, scanline begings in the apex
		if ( (*ultrasoundSettings)[static_cast<int>(settingsIndex)].probeRadius == 0) {
			_scanlineStartPos = _world_T_apex * Eigen::Vector4f(0.0, 0.0, 0.0, 1.0);
		}
		else // curvilinear probe, offset = inner radius of the probe
		{
			Eigen::Vector4f _tempVec = Eigen::Vector4f( _radius * sin(_theta), _radius * cos(_theta),  0.0, 1.0);
			_scanlineStartPos = _world_T_apex * _tempVec;
		}

		
		Eigen::Vector4f _tempVec = Eigen::Vector4f((_depth + _radius) * sin(_theta), (_depth + _radius) * cos(_theta),  0.0, 1.0);
		Eigen::Vector4f _scanlineEndPos = _world_T_apex * _tempVec;

		Eigen::Vector4f _beamDirection = _scanlineEndPos - _scanlineStartPos;

		Eigen::Vector4f _beamDirN = _beamDirection.normalized();

		// now insert the beam direction into the beam collection
#ifdef USE_USHORT16

		Scanline<unsigned short> *_scanline = new Scanline<unsigned short>(Eigen::Vector3f(_scanlineStartPos.x(), _scanlineStartPos.y(), _scanlineStartPos.z()), Eigen::Vector3f(_scanlineEndPos.x(), _scanlineEndPos.y(), _scanlineEndPos.z()),ultrasoundSettings, settingsIndex, scanlineIndex++, fileIndex, imageIndex, i);
		
#endif

#ifdef USE_FLOAT32 || USE_FLOAT5X32
		Scanline<float> *_scanline = new Scanline<float>(Eigen::Vector3f(_scanlineStartPos.x(), _scanlineStartPos.y(), _scanlineStartPos.z()), Eigen::Vector3f(_scanlineEndPos.x(), _scanlineEndPos.y(), _scanlineEndPos.z()),ultrasoundSettings, settingsIndex, scanlineIndex++, fileIndex, imageIndex, i);
#endif
		scanlineVector.push_back(_scanline);

		// now write the direction to the file
		float32 _x, _y, _z;
		_x.value = (_scanlineStartPos/*+_beamDirN*/).x();
		_y.value = (_scanlineStartPos/*+_beamDirN*/).y();
		_z.value = (_scanlineStartPos/*+_beamDirN*/).z();
		file.write(&_x.byte.c[0], sizeof(float32));
		file.write(&_y.byte.c[0], sizeof(float32));
		file.write(&_z.byte.c[0], sizeof(float32));

		_x.value = (_scanlineEndPos/*+_beamDirN*/).x();
		_y.value = (_scanlineEndPos/*+_beamDirN*/).y();
		_z.value = (_scanlineEndPos/*+_beamDirN*/).z();
		file.write(&_x.byte.c[0], sizeof(float32));
		file.write(&_y.byte.c[0], sizeof(float32));
		file.write(&_z.byte.c[0], sizeof(float32));

		//_debugFile << /*scanlineVector.size() << ". " <<*/ _scanlineStartPos.x() << "\t" << _scanlineStartPos.y() << "\t" << _scanlineStartPos.z() << "\t" << _scanlineEndPos.x() << "\t" << _scanlineEndPos.y() << "\t" << _scanlineEndPos.z() << std::endl;
	}
	}

	//_debugFile.close();
}

VoxelDataSet::VoxelDataSet()
{
}


VoxelDataSet::VoxelDataSet(int width, int height, int depth)
{
	setVoxelDimensions(width,height, depth);
}

void VoxelDataSet::setVoxelDimensions(int width, int height, int depth)
{
	m_width = width;
	m_height = height;
	m_depth = depth;
}

template< typename T >
Scanline<T>::Scanline(Eigen::Vector3f &apexPosition, Eigen::Vector3f &beamEnd, std::vector<UltrasoundSettings> *ultrasoundSettings, unsigned char settingsIndex, unsigned int uniqueIndex, unsigned char fileIndex, unsigned short rfImageIndex, unsigned short beamIndex)
	{
		m_line = new Line(apexPosition, beamEnd);
		m_scanlineIndex.fileID = fileIndex;
		m_scanlineIndex.scanlineID = beamIndex;
		m_scanlineIndex.imageID = rfImageIndex;
		m_uniqueIndex = uniqueIndex;
		m_ultrasoundSettings = ultrasoundSettings;
		m_settingsIndex = settingsIndex;

		//m_data = NULL;
	}

template< typename T >
float Scanline<T>::distancePointToScanline(const Eigen::Vector3f &point, float &t)
{
	float _distance;
	Eigen::Vector3f _intersectionPoint;
	m_line->distanceToPoint(point, _distance, t, _intersectionPoint); 
	return _distance;
}
template< typename T >
float Scanline<T>::distancePointToScanlineMAX(const Eigen::Vector3f &point, float &t)
{
	float _distance;
	Eigen::Vector3f _intersectionPoint;
	m_line->distanceToPoint(point, _distance, t, _intersectionPoint); 
	return _distance;
}

