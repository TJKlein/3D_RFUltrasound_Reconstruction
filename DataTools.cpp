#include "DataTools.h"
#include <boost/tokenizer.hpp> 


DataTools::DataTools()
{}

DataTools::~DataTools()
{
}


std::string DataTools::computeElapsedTime(const clock_t &_startTimer, const clock_t &_stopTimer)
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

bool DataTools::getTrackingDataFromFile(std::string filename, std::string target, std::vector<TrackingData> &trackingVector)
{
	std::fstream _file;
	std::string _header;
	std::streamoff _headerOffset;
	int _fileLength;

	_file.open(filename.c_str() , std::ios::in);

	_file.seekg (0, std::ios::end);
	_fileLength = _file.tellg();
	_file.seekg (0, std::ios::beg);


	// parse the header 
	bool isHeader = true;
	bool lastStar = false;
	while (_file.good() && isHeader)
	{
		char ch;
		ch = (char) _file.get();
		if (ch == '*')
		{
			lastStar = true;
		} else
			if (ch == '/' && lastStar)
			{
				isHeader = false;
			} else
			{
				lastStar = false;
			}
			_header.append(1,ch);
	}

	// after header is succesfully parsed calculate header offset
	_headerOffset = _file.tellg();

	while(! ( static_cast<int>(_file.tellg()) == -1  || static_cast<int>(_file.tellg()) >= _fileLength ))
	{
		float tmpTimeStamp;
		int nrTargets, nrPoints;


		_file >> tmpTimeStamp;
		//m_allTimeStamps.push_back(tmpTimeStamp);
		_file >> nrTargets;
		_file >> nrPoints;
		/*
		std::cout << "Time Stamp: " << tmpTimeStamp << std::endl;
		std::cout << "Number of Targets: " << nrTargets << std::endl;
		std::cout << "Number of Points: " << nrPoints << std::endl;
		*/

		for ( int i= 0; i < nrTargets; i++)
		{
			TrackingData td;

			//m_file.flags( std::ios_base::scientific);	
			_file >> td.targetID	 
				>> td.translation[0]			>> td.translation[1]		>> td.translation[2] 
			>> td.quaternion[0]			>> td.quaternion[1]		>>  td.quaternion[2]  >> td.quaternion[3]
			>> td.scale[0]			>> td.scale[1]		>> td.scale[2] 
			>> td.error		>> td.isValid;

			td.timestamp = tmpTimeStamp;
			if ( td.targetID.compare( target) == 0 && td.isValid  && !_file.fail()) // avoid the empty line in file problem
			{
				trackingVector.push_back(td);
			}
		}
	}

	_file.clear();

	return true;
}


bool DataTools::getTargetNamesFromFile(std::string filename, std::vector<std::string> &targetNameVector)
{
	std::fstream _file;
	std::string _header;
	std::streamoff _headerOffset;
	int _fileLength;

	_file.open(filename.c_str() , std::ios::in);

	_file.seekg (0, std::ios::end);
	_fileLength = _file.tellg();
	_file.seekg (0, std::ios::beg);


	// parse the header 
	bool isHeader = true;
	bool lastStar = false;
	while (_file.good() && isHeader)
	{
		char ch;
		ch = (char) _file.get();
		if (ch == '*')
		{
			lastStar = true;
		} else
			if (ch == '/' && lastStar)
			{
				isHeader = false;
			} else
			{
				lastStar = false;
			}
			_header.append(1,ch);
	}

	// after header is succesfully parsed calculate header offset
	_headerOffset = _file.tellg();
	
	// only load up to _maxCOUNT tracking data to check for the target IDs
	int _counter = 0;
	
	int _maxCOUNT = 100;

	while(! ( static_cast<int>(_file.tellg()) == -1  || static_cast<int>(_file.tellg()) >= _fileLength ) && _counter <= _maxCOUNT)
	{
		_counter++;

		float tmpTimeStamp;
		int nrTargets, nrPoints;


		_file >> tmpTimeStamp;
		//m_allTimeStamps.push_back(tmpTimeStamp);
		_file >> nrTargets;
		_file >> nrPoints;
		/*
		std::cout << "Time Stamp: " << tmpTimeStamp << std::endl;
		std::cout << "Number of Targets: " << nrTargets << std::endl;
		std::cout << "Number of Points: " << nrPoints << std::endl;
		*/

		for ( int i= 0; i < nrTargets; i++)
		{
			TrackingData td;

			//m_file.flags( std::ios_base::scientific);	
			_file >> td.targetID	 
				>> td.translation[0]			>> td.translation[1]		>> td.translation[2] 
			>> td.quaternion[0]			>> td.quaternion[1]		>>  td.quaternion[2]  >> td.quaternion[3]
			>> td.scale[0]			>> td.scale[1]		>> td.scale[2] 
			>> td.error		>> td.isValid;

			
			std::vector<std::string>::iterator it = std::find(targetNameVector.begin(), targetNameVector.end(), td.targetID);
			// not found, so we can add it to the target names
			if (it == targetNameVector.end())
			{
				targetNameVector.push_back(td.targetID);
			}
		}
	}

	_file.clear();

	return true;

}

bool DataTools::writeTrackingDataToFile(std::vector<Eigen::Transform3f, Eigen::aligned_allocator<Eigen::Transform3f>> &_data, std::string targetName, std::string filename)
{
	std::ofstream _trackingFile;
	_trackingFile.open(filename);

	// write mandatory header
	_trackingFile << " Tracking data: timestamp <NumOfTargets> <NumOfPoints>  <TargetID> x z z qw qx qy qz sx sy sz error isValid */" << std::endl;

	for (int i=0; i < _data.size();i++)
	{
		_trackingFile << i << "\t1 \t0" << std::endl;

		Eigen::Transform3f _trans  = _data[i];
		Eigen::Quaternionf _quat(_trans.rotation());

		_trackingFile << targetName << "\t" << _trans.translation().x() << "\t" << _trans.translation().y() << "\t" << _trans.translation().z() << "\t" << _quat.w() << "\t" << _quat.x() << "\t" << _quat.y() << "\t" << _quat.z() << "\t1.0\t1.0\t1.0\t0\t1" << std::endl;
	} 

	_trackingFile.close();

	return true;
}

bool DataTools::getTrackingDataFromFileRelative(std::string filename, std::string target, std::string referenceTarget, std::vector<TrackingData> &trackingVector)
{
	std::fstream _file;
	std::string _header;
	std::streamoff _headerOffset;
	int _fileLength;

	_file.open(filename.c_str() , std::ios::in);

	_file.seekg (0, std::ios::end);
	_fileLength = _file.tellg();
	_file.seekg (0, std::ios::beg);


	// parse the header 
	bool isHeader = true;
	bool lastStar = false;
	while (_file.good() && isHeader)
	{
		char ch;
		ch = (char) _file.get();
		if (ch == '*')
		{
			lastStar = true;
		} else
			if (ch == '/' && lastStar)
			{
				isHeader = false;
			} else
			{
				lastStar = false;
			}
			_header.append(1,ch);
	}

	// after header is succesfully parsed calculate header offset
	_headerOffset = _file.tellg();

	int counter = 0;

	while(! ( static_cast<int>(_file.tellg()) == -1  || static_cast<int>(_file.tellg()) >= _fileLength ))
	{
		float tmpTimeStamp;
		int nrTargets, nrPoints;


		_file >> tmpTimeStamp;
		//m_allTimeStamps.push_back(tmpTimeStamp);
		_file >> nrTargets;
		_file >> nrPoints;

		/*
		std::cout << "Time Stamp: " << tmpTimeStamp << std::endl;
		std::cout << "Number of Targets: " << nrTargets << std::endl;
		std::cout << "Number of Points: " << nrPoints << std::endl;
		*/

		TrackingData _referenceTarget;
		TrackingData _toolTarget;
		bool _validReferenceTarget = false;
		bool _validToolTarget = false;


		for ( int i= 0; i < nrTargets; i++)
		{
			TrackingData td;

			//m_file.flags( std::ios_base::scientific);	
			_file >> td.targetID	 
				>> td.translation[0]			>> td.translation[1]		>> td.translation[2] 
			>> td.quaternion[0]			>> td.quaternion[1]		>>  td.quaternion[2]  >> td.quaternion[3]
			>> td.scale[0]			>> td.scale[1]		>> td.scale[2] 
			>> td.error		>> td.isValid;

			td.timestamp = tmpTimeStamp;

			// get the tracking data of the tool target
			if ( td.targetID.compare( target) == 0 && !_file.fail() && td.isValid) // avoid the empty line in file problem
			{
				_validToolTarget = true;
				_toolTarget = td;
				//trackingVector.push_back(td);
			} // get the tracking data of the reference target
			else if ( td.targetID.compare( referenceTarget) == 0 && !_file.fail() && td.isValid) // avoid the empty line in file problem
			{
				_validReferenceTarget = true;
				_referenceTarget = td;
				//trackingVector.push_back(td);
			}
		}

		// if we have acquired both, the tool and the reference position, we got enough data
		// to provide the tool position in the reference target coordinate system
		if ( _validReferenceTarget && _validToolTarget )
		{

			

			Eigen::Quaternionf _toolQuat( _toolTarget.quaternion[0], _toolTarget.quaternion[1], _toolTarget.quaternion[2], _toolTarget.quaternion[3]);
			Eigen::Vector3f _transTool;
			_transTool << _toolTarget.translation[0], _toolTarget.translation[1], _toolTarget.translation[2];

			Eigen::Quaternionf _referenceQuat( _referenceTarget.quaternion[0], _referenceTarget.quaternion[1], _referenceTarget.quaternion[2], _referenceTarget.quaternion[3]);
			Eigen::Vector3f _transReference;
			_transReference << _referenceTarget.translation[0], _referenceTarget.translation[1], _referenceTarget.translation[2];


			// build the 4x4 transformation matrix from the quaternion and 3-vector (tool target)
			Eigen::Transform3f _toolTrackingMatrix;
			_toolTrackingMatrix.setIdentity();
			_toolTrackingMatrix.rotate( _toolQuat.toRotationMatrix());
			_toolTrackingMatrix.translation() = _transTool;
			
			// build the 4x4 transformation matrix from the quaterion and 3-vector (reference target)
			Eigen::Transform3f _referenceTrackingMatrix;
			_referenceTrackingMatrix.setIdentity();
			_referenceTrackingMatrix.rotate( _referenceQuat.toRotationMatrix());
			_referenceTrackingMatrix.translation() = _transReference;



		
			Eigen::Transform3f  _newToolCoordinates = _referenceTrackingMatrix.inverse() * _toolTrackingMatrix;
			Eigen::Quaternionf _newToolQuat(_newToolCoordinates.rotation());
			Eigen::Vector3f _newToolTrans(_newToolCoordinates.translation());


			/*if ( counter == 889 || counter == 890 )
			{
				std::cout << "Reference: " << std::endl;
				std::cout << _referenceTrackingMatrix.matrix() << std::endl << std::endl;
				std::cout << _referenceQuat.w() << ", " << _referenceQuat.x() << ", " << _referenceQuat.y() << ", " << _referenceQuat.z() << std::endl;

				std::cout << "USProbe: " << std::endl;
				std::cout << _toolTrackingMatrix.matrix() << std::endl << std::endl;
				std::cout << _toolQuat.w() << ", " << _toolQuat.x() << ", " << _toolQuat.y() << ", " << _toolQuat.z() << std::endl;
			
				std::cout << "Trans: " << std::endl;
				std::cout << _newToolCoordinates.matrix() << std::endl;
				std::cout << std::endl;

				std::cout << "New Quat: "<< std::endl;
				std::cout << _newToolQuat.w() << ", " << _newToolQuat.x() << ", " << _newToolQuat.y() << ", " << _newToolQuat.z() << std::endl;
			

				std::cout << "Back to Rotation: "<< std::endl;
				Eigen::Matrix3f _testMat = _newToolQuat.toRotationMatrix();
				std::cout << _testMat << std::endl;

					std::cout << std::endl;
			}*/

			TrackingData _td;

			_td = _toolTarget;
			_td.quaternion[0] = _newToolQuat.w();
			_td.quaternion[1] = _newToolQuat.x();
			_td.quaternion[2] = _newToolQuat.y();
			_td.quaternion[3] = _newToolQuat.z();

			_td.translation[0] = _newToolTrans.x();
			_td.translation[1] = _newToolTrans.y();
			_td.translation[2] = _newToolTrans.z();

			trackingVector.push_back(_td);
		}
		counter++;
	}

	_file.clear();

	return true;
}

bool DataTools::getVideoTimestampsFromFile(std::string filename, std::vector<float> &timestampVector)
{
	std::fstream _file;
	int _fileLength;

	_file.open(filename.c_str() , std::ios::in);

	_file.seekg (0, std::ios::end);
	_fileLength = _file.tellg();
	_file.seekg (0, std::ios::beg);

	int counter=0;
	while(! ( static_cast<int>(_file.tellg()) == -1  || static_cast<int>(_file.tellg()) >= _fileLength ))
	{
		counter++;
		float tmpTimeStamp;
		
	
		_file >> tmpTimeStamp;

		
		if( !_file.fail() ) // avoid the empty line at end of file problem
		timestampVector.push_back(tmpTimeStamp);
			
	}

	
	return true;
}

/* TJK_removed_CAMP:
void DataTools::eigen2CAMP(const Eigen::Matrix4f &eigMat, CAMP::Matrix4<float> &campMat)
{
	for(int i=0;i<16;i++)
	{
		campMat.c_array()[i] = eigMat(i);
	}
	campMat = campMat.transpose();
}

void DataTools::eigen2CAMP(const Eigen::Matrix3f &eigMat, CAMP::Matrix3<float> &campMat)
{
	for(int i=0;i<9;i++)
	{
		campMat.c_array()[i] = eigMat(i);
	}
	campMat = campMat.transpose();
}
*/
/*void DataTools::eigen2CAMP(const Eigen::Matrix4f &eigMat, CAMP::Matrix4<float> &campMat)
{
	for(int i=0;i<16;i++)
	{
		campMat.c_array()[i] = eigMat(i);
	}
	campMat = campMat.transpose();
}

void DataTools::eigen2CAMP(const Eigen::Matrix3f &eigMat, CAMP::Matrix3<float> &campMat)
{
	for(int i=0;i<9;i++)
	{
		campMat.c_array()[i] = eigMat(i);
	}
	campMat = campMat.transpose();
}*/

/*
void DataTools::eigen2CAMP(const Eigen::Matrix3f &eigMat, CAMP::Matrix3<float> &campMat)
{
}

void DataTools::eigen2CAMP(const Eigen::Matrix3f &eigMat, CAMP::Matrix3<float> &ampMat)
{
}*/