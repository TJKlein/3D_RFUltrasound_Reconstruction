#include "compounding.h"
//#include <QtWidgets/QApplication>
#include <iostream>
#include <sstream>
#include "Reconstruction.h"

#include <QtWidgets\QApplication>
#include <QtDebug>
#include <QtGlobal>

#include <stdio.h>
#include <stdlib.h>
#include <Windows.h>

//#include "eigen_plugin.h"

/*
#include "Timer.h"

#include "Common/DataTypes.h"

#include "Eigen/Core"

void camp_runner(const CAMP::Matrix4<float>& m, const CAMP::Vector4<float>& v)
{
  CAMP::Vector4<float> n = m*n;
}

void eigen_runner(const Eigen::Matrix4f& m, const Eigen::Vector4f& v)
{
  Eigen::Vector4f n = (m*n).lazy();
}

int main(int argc, char* argv[])
{
	if ( AllocConsole() == 0 )
      {
      ::MessageBox( 0, L"Can't allocate console!\n",
                       L"** PROGRAM NAME ***",
                       MB_OK | MB_ICONSTOP );

      return FALSE;
      }

   // We assume this stuff never fails
   //(!) add better error handling code later...   freopen( "CON", "w", stdout);

   freopen( "CON", "r", stdin );
   freopen( "CON", "w", stderr );
   freopen( "CON", "w", stdout );

   std::ios::sync_with_stdio();
  
  const int num_iter = 10000;
  Timer timer;

  double start = timer.getTimeStamp();
  CAMP::Matrix4<float> m(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16);  
  for (int i=0; i<num_iter; ++i)
  {
    camp_runner(m, CAMP::Vector4<float>(1,2,3,4));
  }
  const double elapsed_camp = (timer.getTimeStamp() - start)*1000.0;

  start = timer.getTimeStamp();
  Eigen::Matrix4f em; em << 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16;
  for (int i=0; i<num_iter; ++i)
  {
    eigen_runner(em, Eigen::Vector4f(1,2,3,4));
  }
  const double elapsed_eigen = (timer.getTimeStamp() - start)*1000.0;

  std::cout << "CAMP: " << elapsed_camp << " ms" << std::endl;
  std::cout << "Eigen: " << elapsed_eigen << " ms" << std::endl;

  char c;
  std::cin >>c;
}
*/
std::string computeElapsedTime(const clock_t &_startTimer, const clock_t &_stopTimer)
{
	int _iHours =(int)floor((float)(_stopTimer - _startTimer)/(CLOCKS_PER_SEC*60*60));
	int _iMinutes = (int)floor((float)(_stopTimer - _startTimer))/(CLOCKS_PER_SEC*60)-_iHours*60;
	int _iSeconds = ((float)(_stopTimer - _startTimer))/CLOCKS_PER_SEC - _iMinutes*60 - _iHours*60*60;
		char _cHours[3], _cMinutes[3], _cSeconds[3];
		sprintf(_cHours,"%d",_iHours);
		if ( strlen(_cHours) < 2)
			sprintf(_cHours,"0%d",_iHours);
		sprintf(_cMinutes,"%d",_iMinutes);
		if ( strlen(_cMinutes) < 2)
			sprintf(_cMinutes,"0%d",_iMinutes);
		sprintf(_cSeconds,"%d",_iSeconds);
		if ( strlen(_cSeconds) < 2)
			sprintf(_cSeconds,"0%d",_iSeconds);

	std::stringstream _ss;
	_ss <<  _cHours << ":" << _cMinutes << ":" << _cSeconds;

	return _ss.str();

}

int main(int argc, char *argv[])
{

	// Create the application's console
  // ISSUE_X64
	#if NDEBUG
    {
	 

  if ( AllocConsole() == 0 )
      {
      ::MessageBox( 0, L"Can't allocate console!\n",
                       L"** PROGRAM NAME ***",
                       MB_OK | MB_ICONSTOP );

      return FALSE;
      }

   // We assume this stuff never fails
   //(!) add better error handling code later...   freopen( "CON", "w", stdout);

   freopen( "CON", "r", stdin );
   freopen( "CON", "w", stderr );
   freopen( "CON", "w", stdout );

   std::ios::sync_with_stdio();
   }
#endif
   
/*

  Reconstruction *_recon = new Reconstruction();

  Eigen::Matrix4f _calib =  Matrix4f::Identity();

  int dim=25;
  int scale=1.0;

  _recon->generateVelocityField(dim, dim, dim, scale, scale, scale, _calib);

    Eigen::Vector3f _position(1.0f, 2.0f, 3.0f);
   Eigen::Vector3f _normal(4.0f, 1.0f, 2.0f);

   float _velocity = 10.0f;

   _recon->pushMeasurement(_normal, _velocity, _position);

   _normal << 1.0f, 1.0f, 6.0f;
   _velocity = 11.0f;

    _recon->pushMeasurement(_normal, _velocity, _position);

	_normal << 9.0f, 4.0f, 1.0f;
	_velocity = 12.0f;

	_recon->pushMeasurement(_normal, _velocity, _position);

	_normal << 7.0f, 1.0f, 9.0f;
	_velocity = 13.0f;


	_recon->pushMeasurement(_normal, _velocity, _position);





   Eigen::MatrixXf _mat;
   Eigen::VectorXf _vec;
   if ( _recon->getMeasurementData(_position, _mat, _vec) )
   {
	   std::cout << _mat << std::endl << std::endl << _vec << std::endl;


	   VectorXf x;

	   _recon->approximateVelocityVector(_position, x);

		std::cout << _mat << std::endl << _vec << std::endl << x << std::endl;
   }


   delete _recon;

 */ 

	/*struct Test
	{
		std::vector<float> data;
		float val;
	};

	std::vector<Test*> _vec;
	clock_t _timer[5];

	_timer[0] = clock();
	for (int i=0;i<500000;i++)
	{
		Test *_test = new Test;
		_test->data.reserve(50);
		_vec.push_back(_test);
	}
	_timer[1] = clock();
	for (int i=0;i<500000;i++)
	{
		
		delete _vec[i];
	}
	_timer[2] = clock();

	std::cout << "Create: " << computeElapsedTime(_timer[0], _timer[1]) << std::endl;
	std::cout << "Delete: " << computeElapsedTime(_timer[1], _timer[2]) << std::endl;
*/
	QApplication a(argc, argv);

	//std::cout << argv[0] << " " << argv[1] << std::endl;
	Compounding w;
	w.show();
	
	if ( argc > 1 )
	{
		
		if ( strcmp(argv[1],"RF")==0 ) 
		{
			// MODE calibrationfile configurationfile trackingfile voxelSize voxelExtent=0,voxelDiameter=1 distance(0=0.5, 1=1,...) compoundingMode(0=GAUSSIAN, 1=NEAREST_NEIGHBOR, 2=MEDIAN, 3=WEIGHTED_MEDIAN, 4=INVERSE_DISTANCE, 5=MEAN) targetName referenceName(optional)
			BATCH_MODE mode;
			mode = BATCH_MODE::RF;

			std::string _reference = "";
			if (argc == 11 )
				_reference = std::string(argv[10]);

		w.batchRun(mode, std::string(argv[2]), std::string(argv[3]), std::string(argv[4]), atof(argv[5]),  (BATCH_DISTANCE_TYPE)atoi(argv[6]), (BATCH_DISTANCE)atoi(argv[7]),  (COMPOUNDING_MODE)atoi(argv[8]), std::string(argv[9]), _reference);

		}
		else if ( strcmp(argv[1],"BMODE")==0 )
		{

			// BMODE <Number of data sets> ouputFilename calibrationFilename video1 tracking1 ...
			BATCH_MODE mode;
			mode = BATCH_MODE::BMODE;

			

			std::vector<std::string> _trackingFiles;
			std::vector<std::string> _videoFiles;
		for(int k=0;k<atoi(argv[2])*2;k+=2)
		{
			_trackingFiles.push_back(std::string(argv[5+k+1]));
			_videoFiles.push_back(std::string(argv[5+k]));
		}


		int _offset = atoi(argv[2])*2 +5;

		std::string _reference = "";
			if (argc == _offset+8-1 )
				_reference = std::string(argv[_offset+6]);

			w.batchRunBMODE(mode, std::string(argv[4]), std::string(argv[3]), _videoFiles, _trackingFiles, atof(argv[_offset+0]),  (BATCH_DISTANCE_TYPE)atoi(argv[_offset+1]), (BATCH_DISTANCE)atoi(argv[_offset+2]),  (COMPOUNDING_MODE)atoi(argv[_offset+3]),(bool)atoi(argv[_offset+4]), std::string(argv[_offset+5]),  _reference);

		}
		else if ( strcmp(argv[1],"BUTTERWORTHENVELOPE")==0 )
		{
			// filename lowCut highCut order frequency
			w.batchButterwothEnvelope(std::string(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), atoi(argv[6]));
		}
		else if ( strcmp(argv[1],"INTENSITYCOMPOUND")==0 )
		{
			// filename voxelExtent=0,voxelDiameter=1 distance(0=0.5, 1=1,...) logCompress(0=false, 1=true)
			w.batchIntensityCompound(std::string(argv[2]),  (COMPOUNDING_MODE)atoi(argv[3]), (BATCH_DISTANCE_TYPE)atoi(argv[4]), (BATCH_DISTANCE)atoi(argv[5]), atoi(argv[6]));

		}
		
		return 0;
	}
	return a.exec();
}
