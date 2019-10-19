#include "Reconstruction.h"
#include <iostream>
#include <limits>

#if MATLAB_SUPPORT

#include "mat.h"
#include "matrix.h"

#endif



ColorLUT::~ColorLUT()
{
}


bool ColorLUT::colorToVelocity(const int &R, const int &G, const int &B, float &velocity)
{
	stdext::hash_map<DWORD, float>::iterator _hashMapIterator;

	_hashMapIterator = m_hashMapVelocity.find(RGB(R,G,B));

	if ( _hashMapIterator == m_hashMapVelocity.end() )
	{
		_hashMapIterator = m_hashMapVelocity.begin();

		float distance = 99999999.0;
		float _foundV;
		int _foundR, _foundG, _foundB;
		COLORREF _foundCol;
		while(_hashMapIterator != m_hashMapVelocity.end())
		{
			COLORREF _col = _hashMapIterator->first;
			int _R = (int)GetRValue(_col);
			int _G = (int)GetGValue(_col);
			int _B = (int)GetBValue(_col);

			Eigen::Vector3f _temp((float)(_R-R),(float)(_G-G),(float)(_B-B));

			if ( _temp.norm() < distance)
			{
				distance = _temp.norm();
				_foundR = _R;
				_foundG = _G;
				_foundB = _B;
				_foundV = _hashMapIterator->second;
				_foundCol = _col;
			}
			_hashMapIterator++;
		}
			
			m_hashMapVelocity[_foundCol] = _foundV;
			velocity = _foundV;
		// find the closest
		return false;
	}
	
	velocity = _hashMapIterator->second;

	return true;
}

ColorLUT::ColorLUT(std::vector<int> interpolation, std::vector<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f>> color, const float &maxVelocity) :
m_LUT(0)
{

	m_maxVelocity = maxVelocity;

	int _length = 0;
	for(int i=0;i<interpolation.size();i++)
	{
		_length += interpolation.at(i);
	}


	float _velocityDiff = 2.0f*m_maxVelocity / static_cast<float>(_length);

	int counter = 0;
	for(int j=0;j<interpolation.size();j++)
	{
		Eigen::Vector3f _colorDiff;
		_colorDiff.noalias() = (color[j+1] - color[j])/static_cast<float>(interpolation[j]);
		for(int i=0;i<interpolation[j];i++)
		{
			m_LUT.push_back(_colorDiff*static_cast<float>(i)+color[j]);
			m_velocity.push_back(m_maxVelocity-static_cast<float>(counter)*_velocityDiff);

			COLORREF _tempCol = RGB(static_cast<int>(m_LUT[counter][0]), static_cast<int>(m_LUT[counter][1]), static_cast<int>(m_LUT[counter][2]));
			
			m_hashMapVelocity[_tempCol] = m_velocity[counter];
			counter++;
		}
	}
	
	Eigen::Vector3f _lastColor =color[interpolation.size()];
	COLORREF _colRef = RGB(static_cast<int>(_lastColor(0)), static_cast<int>(_lastColor(1)), static_cast<int>(_lastColor(2)));

	m_LUT.push_back(_lastColor);
	m_velocity.push_back(-m_maxVelocity);
	m_hashMapVelocity[_colRef] = -m_maxVelocity;

}




Reconstruction::Reconstruction() :
m_calibrationMatrixSet(false),
m_velocityField(0),
m_velocityFieldUpdated(0),
m_velocityFieldValid(0)
{

}

Reconstruction::~Reconstruction()
{
	if ( m_velocityField) delete[] m_velocityField;
	if ( m_velocityFieldUpdated) delete[] m_velocityFieldUpdated;
	if ( m_velocityFieldValid) delete[] m_velocityFieldValid;


	if ( m_velocityFieldBackup) delete[] m_velocityFieldBackup;
	if ( m_velocityFieldValidBackup) delete[] m_velocityFieldValidBackup;
	
	for (int i=0;i<m_measurementField.size();i++)
	{
		delete m_measurementField[i];
	}
}

void Reconstruction::createFieldBackup()
{
	for (int i=0; i<m_dimX*m_dimY*m_dimZ;i++)
	{
		m_velocityFieldBackup[i] = m_velocityField[i];
		m_velocityFieldValidBackup[i] = m_velocityFieldValid[i];
	}
}

void Reconstruction::restoreFieldBackup()
{
	for (int i=0; i<m_dimX*m_dimY*m_dimZ;i++)
	{
		m_velocityField[i] = m_velocityFieldBackup[i];
		m_velocityFieldValid[i] = m_velocityFieldValidBackup[i];
	}
}

bool Reconstruction::getVelocityFieldValidity(const int &x, const int &y, const int &z)
{
	int _index;
	getIndex(x,y,z,_index);
	return m_velocityFieldValid[_index];
}  

bool Reconstruction::getVelocityFieldValidity(const int &index)
{
	return m_velocityFieldValid[index];
}  

bool Reconstruction::setVelocityFieldValidity(const int &x, const int &y, const int &z, const bool &validity)
{
	int _index;
	getIndex(x,y,z,_index);
	m_velocityFieldValid[_index] = validity;

	return true;
}

void Reconstruction::setMaxError(const float &error)
{
	m_maxError = error;
}

float Reconstruction::getMaxError()
{
	return m_maxError;
}

bool Reconstruction::setVelocity(const int &x, const int &y, const int &z, const Eigen::Vector3f &velocity, bool validity)
{
	int _index = z*m_dimX*m_dimY+y*m_dimX+x;

	if ( _index < m_dimX*m_dimY*m_dimZ )
	{
		m_velocityField[_index] = velocity;
		m_velocityFieldValid[_index] = validity;
		return true;
	}
	else
		return false;
}


bool Reconstruction::getProjection(const int &_index, Eigen::VectorXf &solutionVector)
{


	if (  m_velocityFieldValid[_index] )
	{
		VelocityData *_data = m_measurementField[_index];


		Eigen::VectorXf _error;
		_error = Eigen::VectorXf::Zero(_data->normalMatrix.rows()/*, _data->normalMatrix.cols()*/);

		// check if we have acquired data for that field, otherwise the matrix will be empty
		if ( !(_data->normalMatrix.rows() > 0) || !(_data->normalMatrix.cols() > 0) )
			return false;

		/*std::cout << "1.\t" <<_data->normalMatrix << std::endl << std::endl;
		std::cout << "2.\t" << m_velocityField[_index] << std::endl << std::endl;
		std::cout << "3.\t" << _data->velocityVector << std::endl;
		*/
		_error = (_data->normalMatrix * m_velocityField[_index] - _data->velocityVector);
		

		//std::cout << "Normals: " << _data->normalMatrix << std::endl;

		//std::cout <<  "Error: " << _error.transpose() <<std::endl;

		float _sumError = 0.0f;
		Eigen::VectorXf _sumNormal;
		_sumNormal = Eigen::VectorXf::Zero(3);
		

		for(int i=0;i<_data->normalMatrix.rows();i++)
		{
			// actually this should not be norm but squared only
			_sumNormal += _data->normalMatrix.row(i).transpose() * _error.row(i)[0];

			float _temp = _error.row(i)[0];
			//std::cout <<"Error: " << _temp << std::endl; 
			//float _temp = _error.row(i)[0];
			_sumError += _temp;
		}

		solutionVector = _sumNormal;/* * _sumError;*/
		//std::cout << "Error Matrix: " << errorMatrix << std::endl;


		return true;
	}
	return false;

}

Eigen::Vector3f Reconstruction::computeDivergence(const int &x, const int &y, const int &z, const float &physicalX, const float &physicalY, const float &physicalZ)
{
	float _temp[] = {0.0, 0.0, 0.0};
	_temp[0] += partialDerivativeComponent(x,y,z,0,0,0,physicalX,physicalX);
	_temp[0] += partialDerivativeComponent(x,y,z,0,1,1,physicalX,physicalY);
	_temp[0] += partialDerivativeComponent(x,y,z,0,2,2,physicalX,physicalZ);


	_temp[1] += partialDerivativeComponent(x,y,z,1,0,0,physicalY,physicalX);
	_temp[1] += partialDerivativeComponent(x,y,z,1,1,1,physicalY,physicalY);
	_temp[1] += partialDerivativeComponent(x,y,z,1,2,2,physicalY,physicalZ);

	_temp[2] += partialDerivativeComponent(x,y,z,2,0,0,physicalZ,physicalX);
	_temp[2] += partialDerivativeComponent(x,y,z,2,1,1,physicalZ,physicalY);
	_temp[2] += partialDerivativeComponent(x,y,z,2,2,2,physicalZ,physicalZ);


	return Eigen::Vector3f(_temp[0], _temp[1], _temp[2]);
}


float Reconstruction::partialDerivativeComponent(int x, int y, int z, const int &velocityComponent, const int &x1, const int &x2, const float &hx, const float &hy)
{
	int _index[4];

	int _incX1[3];
	int _incX2[3];

	switch(x1)
	{
	case 0:
		{
			_incX1[0]=1;
			_incX1[1]=0;
			_incX1[2]=0;
			break;
		}
	case 1:
		{
			_incX1[0]=0;
			_incX1[1]=1;
			_incX1[2]=0;
			break;
		}
	case 2:
		{
			_incX1[0]=0;
			_incX1[1]=0;
			_incX1[2]=1;
			break;
		}
	}

	switch(x2)
	{
	case 0:
		{
			_incX2[0]=1;
			_incX2[1]=0;
			_incX2[2]=0;
			break;
		}
	case 1:
		{
			_incX2[0]=0;
			_incX2[1]=1;
			_incX2[2]=0;
			break;
		}
	case 2:
		{
			_incX2[0]=0;
			_incX2[1]=0;
			_incX2[2]=1;
			break;
		}
	}

	if ( getIndex(x,y,z,_index[0]) && m_velocityFieldValid[_index[0]])
	{

		// check for border proximity an move inwards (somewhat dirty)
		if ( x == 0 )
			x++;
		else if ( x == m_dimX-1 )
			x--;
		if ( y == 0 )
			y++;
		else if ( y == m_dimY-1 )
			y--;
		if ( z == 0)
			z++;
		else if ( z == m_dimZ-1 )
			z--;

		if ( x1 != x2 )
		{
			
				float _velocity[4];

				
				
				getIndex(x+_incX1[0]+_incX2[0],y+_incX1[1]+_incX2[1],z+_incX1[2]+_incX2[2],_index[0]);
				getIndex(x+_incX1[0]-_incX2[0],y+_incX1[1]-_incX2[1],z+_incX1[2]-_incX2[2],_index[1]);
				getIndex(x-_incX1[0]+_incX2[0],y-_incX1[1]+_incX2[1],z-_incX1[2]+_incX2[2],_index[2]);
				getIndex(x-_incX1[0]-_incX2[0],y-_incX1[1]-_incX2[1],z-_incX1[2]-_incX2[2],_index[3]);
			

				_velocity[0] = m_velocityField[_index[0]][velocityComponent];
				_velocity[1] = m_velocityField[_index[1]][velocityComponent];
				_velocity[2] = m_velocityField[_index[2]][velocityComponent];
				_velocity[3] = m_velocityField[_index[3]][velocityComponent];

				return (_velocity[0]-_velocity[1])-(_velocity[2]-_velocity[3])/4.0*hx*hy;

		}
		else
		{
			float _velocity[3];

			getIndex(x+_incX1[0],y+_incX1[1],z+_incX1[2],_index[0]);
			getIndex(x-_incX1[0],y-_incX1[1],z-_incX1[2],_index[1]);
			getIndex(x,y,z,_index[2]);

			_velocity[0] = m_velocityField[_index[0]][velocityComponent];
			_velocity[1] = m_velocityField[_index[1]][velocityComponent];
			_velocity[2] = m_velocityField[_index[2]][velocityComponent];

			return (_velocity[0]+_velocity[1]-2.0*_velocity[3])/(hx*hx);
		}
	}
}

bool Reconstruction::getProjection(const int &x, const int &y, const int &z, Eigen::VectorXf &solutionVector)
{
	int _index ;

	if ( getIndex(x,y,z,_index) && m_velocityFieldValid[_index] )
	{
		VelocityData *_data = m_measurementField[_index];


		Eigen::VectorXf _error;
		_error = Eigen::VectorXf::Zero(_data->normalMatrix.rows()/*, _data->normalMatrix.cols()*/);
		_error = (_data->normalMatrix * m_velocityField[_index] - _data->velocityVector);
		

		//std::cout << "Normals: " << _data->normalMatrix << std::endl;

		//std::cout <<  "Error: " << _error.transpose() <<std::endl;

		float _sumError = 0.0f;
		Eigen::VectorXf _sumNormal;
		_sumNormal = Eigen::VectorXf::Zero(3);
		

		for(int i=0;i<_data->normalMatrix.rows();i++)
		{
			// actually this should not be norm but squared only
			_sumNormal += _data->normalMatrix.row(i).transpose();

			float _temp = _error.row(i)[0];
			//std::cout <<"Error: " << _temp << std::endl; 
			//float _temp = _error.row(i)[0];
			_sumError += _temp;
		}

		solutionVector = _sumNormal * _sumError;
		//std::cout << "Error Matrix: " << errorMatrix << std::endl;


		return true;
	}
	return false;
}


bool Reconstruction::getApproximationError(const int &_index, Eigen::VectorXf  &errorMatrix)
{

	if ( m_velocityFieldValid[_index] )
	{
		VelocityData *_data = m_measurementField[_index];

		errorMatrix = Eigen::VectorXf::Zero(_data->normalMatrix.rows()/*, _data->normalMatrix.cols()*/);

		//errorMatrix = Eigen::VectorXf::Zero(_data->normalMatrix.rows());
		errorMatrix = (_data->normalMatrix * m_velocityField[_index] - _data->velocityVector);
		

		//std::cout << "Normals: " << _data->normalMatrix << std::endl;

		//std::cout <<  "Error: " << _error.transpose() <<std::endl;
		

		/*for(int i=0;i<_data->normalMatrix.rows();i++)
		{
			// actually this should not be norm but squared only
			errorMatrix.row(i) = (_data->normalMatrix.row(i) * _error(i));
		}*/

		//std::cout << "Error Matrix: " << errorMatrix << std::endl;


		return true;
	}
	return false;
}



bool Reconstruction::getApproximationError(const int &x, const int &y, const int &z, Eigen::VectorXf  &errorMatrix)
{

	int _index ;

	if ( getIndex(x,y,z,_index) && m_velocityFieldValid[_index] )
	{
		VelocityData *_data = m_measurementField[_index];

		errorMatrix = Eigen::VectorXf::Zero(_data->normalMatrix.rows()/*, _data->normalMatrix.cols()*/);

		//errorMatrix = Eigen::VectorXf::Zero(_data->normalMatrix.rows());
		errorMatrix = (_data->normalMatrix * m_velocityField[_index] - _data->velocityVector);
		

		//std::cout << "Normals: " << _data->normalMatrix << std::endl;

		//std::cout <<  "Error: " << _error.transpose() <<std::endl;
		

		/*for(int i=0;i<_data->normalMatrix.rows();i++)
		{
			// actually this should not be norm but squared only
			errorMatrix.row(i) = (_data->normalMatrix.row(i) * _error(i));
		}*/

		//std::cout << "Error Matrix: " << errorMatrix << std::endl;


		return true;
	}
	return false;
}


bool Reconstruction::getApproximationTemporaryError(const int &_index, Eigen::VectorXf  &errorMatrix)
{


	if (  m_velocityFieldValid[_index] )
	{
		VelocityData *_data = m_measurementField[_index];

		errorMatrix = Eigen::VectorXf::Zero(_data->normalMatrix.rows()/*, _data->normalMatrix.cols()*/);

		//errorVector = Eigen::VectorXf::Zero(_data->normalMatrix.rows());
		errorMatrix = (_data->normalMatrix * m_velocityFieldUpdated[_index] - _data->velocityVector);
		

		//std::cout << "Normals: " << _data->normalMatrix << std::endl;

		//std::cout <<  "Error: " << _error.transpose() <<std::endl;
		

		/*for(int i=0;i<_data->normalMatrix.rows();i++)
		{
			// actually this should not be norm but squared only
			errorMatrix.row(i) = (_data->normalMatrix.row(i) * _error(i));
		}*/

		//std::cout << "Error Matrix: " << errorMatrix << std::endl;


		return true;
	}
	return false;

/*
	int _index ;

	if ( getIndex(x,y,z,_index)  )
	{
		VelocityData *_data = m_measurementField[_index];

		errorMatrix = Eigen::MatrixXf::Zero(_data->normalMatrix.rows(), _data->normalMatrix.cols());

		Eigen::VectorXf _error = (_data->normalMatrix * m_velocityFieldUpdated[_index] - _data->velocityVector);
		

		//std::cout << "Normals: " << _data->normalMatrix << std::endl;

		//std::cout <<  "Error: " << _error.transpose() <<std::endl;
		
		for(int i=0;i<_data->normalMatrix.rows();i++)
		{
			// actually this should not be norm but squared only
			errorMatrix.row(i) = (_data->normalMatrix.row(i) * _error(i));
		}

		//std::cout << "Error Matrix: " << errorMatrix << std::endl;


		return true;
	}
	return false;*/
}


bool Reconstruction::getApproximationTemporaryError(const int &x, const int &y, const int &z, Eigen::VectorXf  &errorMatrix)
{

	
	int _index ;

	if ( getIndex(x,y,z,_index) && m_velocityFieldValid[_index] )
	{
		VelocityData *_data = m_measurementField[_index];

		errorMatrix = Eigen::VectorXf::Zero(_data->normalMatrix.rows()/*, _data->normalMatrix.cols()*/);

		//errorVector = Eigen::VectorXf::Zero(_data->normalMatrix.rows());
		errorMatrix = (_data->normalMatrix * m_velocityFieldUpdated[_index] - _data->velocityVector);
		

		//std::cout << "Normals: " << _data->normalMatrix << std::endl;

		//std::cout <<  "Error: " << _error.transpose() <<std::endl;
		

		/*for(int i=0;i<_data->normalMatrix.rows();i++)
		{
			// actually this should not be norm but squared only
			errorMatrix.row(i) = (_data->normalMatrix.row(i) * _error(i));
		}*/

		//std::cout << "Error Matrix: " << errorMatrix << std::endl;


		return true;
	}
	return false;

/*
	int _index ;

	if ( getIndex(x,y,z,_index)  )
	{
		VelocityData *_data = m_measurementField[_index];

		errorMatrix = Eigen::MatrixXf::Zero(_data->normalMatrix.rows(), _data->normalMatrix.cols());

		Eigen::VectorXf _error = (_data->normalMatrix * m_velocityFieldUpdated[_index] - _data->velocityVector);
		

		//std::cout << "Normals: " << _data->normalMatrix << std::endl;

		//std::cout <<  "Error: " << _error.transpose() <<std::endl;
		
		for(int i=0;i<_data->normalMatrix.rows();i++)
		{
			// actually this should not be norm but squared only
			errorMatrix.row(i) = (_data->normalMatrix.row(i) * _error(i));
		}

		//std::cout << "Error Matrix: " << errorMatrix << std::endl;


		return true;
	}
	return false;*/
}

bool Reconstruction::updateVelocity(const int &x, const int &y, const int &z, const Eigen::Vector3f &velocity)
{
	int _index;

	if ( getIndex(x,y,z, _index) )
	{
		m_velocityFieldUpdated[_index] = velocity;
		return true;
	}
	else
		return false;
}


bool Reconstruction::updateDeltaVelocity(const int &_index, const Eigen::Vector3f &velocity)
{

		m_velocityFieldUpdated[_index] = m_velocityField[_index] + velocity;
		return true;
}

bool Reconstruction::updateDeltaVelocity(const int &x, const int &y, const int &z, const Eigen::Vector3f &velocity)
{
	int _index;

	if ( getIndex(x,y,z, _index) )
	{
		m_velocityFieldUpdated[_index] = m_velocityField[_index] + velocity;
		return true;
	}
	else
		return false;
}


void Reconstruction::performVectorFieldUpdate()
{
	Eigen::Vector3f *temp;
	temp = m_velocityField;
	m_velocityField = m_velocityFieldUpdated;
	m_velocityFieldUpdated = temp;
}
						 
bool Reconstruction::getNormal(int x, int y, Eigen::Vector3f &normal)
{
	if ( y*m_roiWidth+x < m_roiHeight*m_roiWidth )
	{
		normal = m_normalPlane[y*m_roiWidth+x];
		return true;
	}
	else return false;
}

bool Reconstruction::getSpace(int x, int y, Eigen::Vector3f &space)
{
	if ( y*m_roiWidth+x < m_roiHeight*m_roiWidth )
	{
		space = m_spacePlane[y*m_roiWidth+x];
		return true;
	}
	else return false;
}

bool Reconstruction::preComputePlanes(const float &apexX, const float &apexY, const int &minRoiX, const int &minRoiY, const int &maxRoiX, const int &maxRoiY)
{
	m_roiWidth = (maxRoiX - minRoiX);
	m_roiHeight = (maxRoiY - minRoiY);



	m_normalPlane.clear();
	m_spacePlane.clear();

	for (int y=minRoiY;y<maxRoiY;y++)
	{
		for (int x=minRoiX;x<maxRoiX;x++)
		{
			Eigen::Vector3f _temp(static_cast<float>(x*m_usScaleX),static_cast<float>(y*m_usScaleY),0.0);

			Eigen::Vector3f _pointWorld = (m_calibrationMatrix * _temp);


			Eigen::Vector3f _apexWorld = (m_calibrationMatrix * Eigen::Vector3f(apexX, apexY,0.0));
			Eigen::Vector3f _normal = _pointWorld - _apexWorld;
			_normal = _normal.normalized();

		//	std::cout << _normal << std::endl;

			//char c;
		//	std::cin >> c;

			//std::cout << y*m_roiWidth+x << " vs. " << m_spacePlane.size() << std::endl; 
			m_spacePlane.push_back(_pointWorld);
			m_normalPlane.push_back(_normal);
		}
	}

	return true;
}


void Reconstruction::setCalibration(Eigen::Matrix4f &calibrationMatrix, const float &scaleX, const float &scaleY)
{
	m_usScaleX = scaleX;
	m_usScaleY = scaleY;

	m_calibrationMatrix = calibrationMatrix;

	m_calibrationMatrixSet = true;
}

#if MATLAB_SUPPORT

void Reconstruction::saveNormalPlaneAsMAT()
{
	MATFile *pmat;
	mxArray *pa1;

	int status; 

	pmat = matOpen("Normal.mat", "w");

	const mwSize dims[]={m_roiHeight,m_roiWidth, 3};
	pa1 = mxCreateNumericArray (3,dims,	mxDOUBLE_CLASS, mxREAL);
	char *pdata = reinterpret_cast<char *>(mxGetData(pa1));

	

	double counter = 0.0;
	/*for (int y=0;y<10;y++)
	{
		for(int x=0;x<10;x++)
		{
			for(int z=0;z<3;z++)
			{
				counter += 1.0;

				std::cout << counter << std::endl;
				memcpy(pdata,&counter,mxGetElementSize(pa1));
				pdata+=mxGetElementSize(pa1);
			}
		}
	}*/
	for(int x=0;x<m_roiWidth;x++)
	{
		for (int y=0;y<m_roiHeight;y++)
		{
			const int _index = y*m_roiWidth+x;
			const Eigen::Vector3d normal = m_normalPlane[_index].cast<double>();

			//std::cout << normal << std::endl;

			memcpy(pdata,&normal[0],mxGetElementSize(pa1));
			memcpy(pdata+m_roiHeight*m_roiWidth*mxGetElementSize(pa1),&normal[1],mxGetElementSize(pa1));
			memcpy(pdata+2*m_roiHeight*m_roiWidth*mxGetElementSize(pa1),&normal[2],mxGetElementSize(pa1));

			pdata+=mxGetElementSize(pa1);
		}
	}
status = matPutVariable(pmat, "NormalPlane", pa1);
	mxDestroyArray(pa1);

	matClose(pmat);
}



void Reconstruction::saveMatlabMatrix1Dim(const Eigen::MatrixXf &data, const int &width, const int &height, char* name, char* filename)
{
	MATFile *pmat;
	mxArray *pa1;

	int status; 

	pmat = matOpen(filename, "w");

	const mwSize dims[]={height,width, 1};
	pa1 = mxCreateNumericArray (2,dims,	mxDOUBLE_CLASS, mxREAL);
	char *pdata = reinterpret_cast<char *>(mxGetData(pa1));


	for(int x=0;x<width;x++)
	{
		for (int y=0;y<height;y++)
		{
			double _tempData = static_cast<double>(data(y,x));

			memcpy(pdata,&_tempData,mxGetElementSize(pa1));

			pdata+=mxGetElementSize(pa1);
		}
	}
	status = matPutVariable(pmat, name, pa1);
	mxDestroyArray(pa1);

	matClose(pmat);
}

void Reconstruction::saveMatlabMatrix3Dim(std::vector<Eigen::Vector3f> &data, const int &width, const int &height, char* name, char* filename)
{
	MATFile *pmat;
	mxArray *pa1;

	int status; 

	pmat = matOpen(filename, "w");

	const mwSize dims[]={height,width, 3};
	pa1 = mxCreateNumericArray (3,dims,	mxDOUBLE_CLASS, mxREAL);
	char *pdata = reinterpret_cast<char *>(mxGetData(pa1));



	double counter = 0.0;

	for(int x=0;x<width;x++)
	{
		for (int y=0;y<height;y++)
		{
			const int _index = y*width+x;
			const Eigen::Vector3d _dataElement = data[_index].cast<double>();

			//std::cout << normal << std::endl;

			memcpy(pdata,&_dataElement[0],mxGetElementSize(pa1));
			memcpy(pdata+height*width*mxGetElementSize(pa1),&_dataElement[1],mxGetElementSize(pa1));
			memcpy(pdata+2*height*width*mxGetElementSize(pa1),&_dataElement[2],mxGetElementSize(pa1));

			pdata+=mxGetElementSize(pa1);
		}
	}
	status = matPutVariable(pmat, name, pa1);
	mxDestroyArray(pa1);

	matClose(pmat);
}


#endif


void Reconstruction::generateVelocityField(const int &dimX, const int &dimY, const int &dimZ, const float &scaleX, const float &scaleY, const float &scaleZ, const Eigen::Matrix4f &origin)
{

	m_dimX = dimX;
	m_dimY = dimY;
	m_dimZ = dimZ;

	m_scaleX = scaleX;
	m_scaleY = scaleY;
	m_scaleZ = scaleZ;

	if ( m_velocityField) delete[] m_velocityField;
	if ( m_velocityFieldUpdated) delete[] m_velocityFieldUpdated;
	if ( m_velocityFieldValid) delete[] m_velocityFieldValid;
	//if ( m_measurementField ) delete[] m_measurementField;

	// the following data will represent the vector field
	m_velocityField = new Eigen::Vector3f[dimX*dimY*dimZ];
	m_velocityFieldUpdated = new Eigen::Vector3f[dimX*dimY*dimZ];
	m_velocityFieldValid = new bool[dimX*dimY*dimZ];

	// backup of the vector field representation
	m_velocityFieldBackup = new Eigen::Vector3f[dimX*dimY*dimZ];
	m_velocityFieldValidBackup = new bool[dimX*dimY*dimZ];
	//m_measurementField = new std::vector<VelocityData*>[dimX*dimY*dimZ];

	m_measurementField.reserve(dimX*dimY*dimZ);

	for (int i=0;i<dimX*dimY*dimZ;i++)
	{
		VelocityData *_data = new VelocityData();
		m_measurementField.push_back(_data);
	}

	m_velocityField2World = origin;
	m_world2velocityField.matrix() = origin.inverse();
}

bool Reconstruction::getIndex(const Eigen::Vector3f& worldPosition, int &index)
{
	Eigen::Vector3f _local;
	_local.noalias() = (m_world2velocityField*worldPosition);

	if ( _local[0] < 0.0 || _local[1] < 0.0 || _local[2] < 0.0 )
		return false; // not inside the volume

	int _x = floor(_local[0]/m_scaleX ), _y = floor(_local[1] / m_scaleY ), _z = floor(_local[2] / m_scaleZ);

	if ( _x >= m_dimX || _y >= m_dimY || _z >= m_dimZ )
		return false; // not inside the volume

	index = _z*m_dimX*m_dimY+_y*m_dimX+_x;

	return true;
}

Eigen::Vector3f Reconstruction::getVelocityVectorByCoordinates(const int &_x, const int &_y, const int &_z)
{
	int index = _z*m_dimX*m_dimY+_y*m_dimX+_x;

	return m_velocityField[index];
}

Eigen::Vector3f Reconstruction::getVelocityVectorByIndex(const int &index)
{

	return m_velocityField[index];
}

void Reconstruction::setVelocityVectorByIndex(const int &index, Eigen::Vector3f &_vec)
{

	m_velocityField[index] = _vec;
}

VelocityData* Reconstruction::getMeasurementsByCoordinates(const int &_x, const int &_y, const int &_z)
{
	int index = _z*m_dimX*m_dimY+_y*m_dimX+_x;

	return m_measurementField[index];
}

VelocityData* Reconstruction::getMeasurementsByIndex(int index)
{
	return m_measurementField[index];
}
bool Reconstruction::getClosestVelocity(const Eigen::Vector3f &worldPosition, Eigen::Vector3f &velocity)
{

	int _index;
	if (  getIndex(worldPosition, _index) )
	{
		// position is within the volume, grab the right one
		velocity = m_velocityField[_index];
		return true;
	}
	else
		return false;
}

bool Reconstruction::getClosestMeasurements(const Eigen::Vector3f &worldPosition, VelocityData * velocityData)
{
	int _index;

	if (  getIndex(worldPosition, _index) )
	{
		// position is within the volume, grab the right one
		velocityData = m_measurementField[_index];
		return true;
	}
	else
		return false;
}

bool Reconstruction::approximateVelocityVector(const Eigen::Vector3f &worldPosition, Eigen::VectorXf &velocity)
{
	int _index;

	if ( getIndex(worldPosition, _index) )
	{
		VelocityData *_data = m_measurementField[_index];

		/*Eigen::MatrixXf D = _data->normalMatrix;
		Eigen::MatrixXf A = D.transpose() * D;
		Eigen::VectorXf b = D.transpose() * _data->velocityVector;

		A.llt().solve(b,&velocity);   // using a LLT factorization
		A.ldlt().solve(b,&velocity);  // using a LDLT factorization
		*/

		velocity = _data->normalMatrix.svd().solve(_data->velocityVector);
		return true;
	}
	else
		return false;
}



bool Reconstruction::getIndex(const int &x, const int &y, const int &z, int &index)
{
	if ( x < 0 || y < 0 || z < 0 ||x >= m_dimX || y >= m_dimY || z >= m_dimZ )
		return false;
	index = z*m_dimX*m_dimY+y*m_dimX+x;

	return true;
}

void Reconstruction::getCoordinatesFromIndex(const int &index, int &x, int &y, int &z)
{
	z = (int)((float)index / (float)(m_dimX*m_dimY));
	y = (int)((float)(index-z*(m_dimX*m_dimY) ) / (float)(m_dimX));
	x = index - z*(m_dimX*m_dimY) - y * m_dimX;
}

Eigen::Vector3f Reconstruction::computeLaplacian(const int &x, const int &y, const int &z, const float &physicalX, const float &physicalY, const float &physicalZ)
{
	bool _borderProximity = false;
	int _index[3];
	Eigen::Vector3f _laplacian(0.0, 0.0, 0.0);
	Eigen::Vector3f _velocity[3];
	if ( x == 0 || y == 0 || z == 0 || x == m_dimX-1 || y == m_dimY-1 || z == m_dimZ-1)
		_borderProximity = true;

	if ( !_borderProximity)
	{

		getIndex(x-1,y,z,_index[0]);
		getIndex(x+1,y,z,_index[1]);
		getIndex(x,y,z,_index[2]);
		_laplacian += (m_velocityField[_index[0]]+m_velocityField[_index[1]]-2.0f*m_velocityField[_index[2]] ) / (4.0f*physicalX*physicalX);


		getIndex(x,y-1,z,_index[0]);
		getIndex(x,y+1,z,_index[1]);
		getIndex(x,y,z,_index[2]);
		_laplacian += (m_velocityField[_index[0]]+m_velocityField[_index[1]]-2.0f*m_velocityField[_index[2]] ) / (4.0f*physicalY*physicalY);


		getIndex(x,y,z-1,_index[0]);
		getIndex(x,y,z+1,_index[1]);
		getIndex(x,y,z,_index[2]);
		_laplacian += (m_velocityField[_index[0]]+m_velocityField[_index[1]]-2.0f*m_velocityField[_index[2]] ) / (4.0f*physicalZ*physicalZ);

	}
	else
	{
		if (getIndex(x-1,y,z,_index[0]) )
			_velocity[0] = m_velocityField[_index[0]];
		else
			_velocity[0] = Eigen::Vector3f::Zero(3,1);

		if ( getIndex(x+1,y,z,_index[1]) )
			_velocity[1] = m_velocityField[_index[1]];
		else
			_velocity[1] = Eigen::Vector3f::Zero(3,1);

		if ( getIndex(x,y,z,_index[2]) )
			_velocity[2] = m_velocityField[_index[2]];
		else
			_velocity[2] = Eigen::Vector3f::Zero(3,1);
		_laplacian += (_velocity[0]+_velocity[1]-2.0f*_velocity[2] ) / (4.0f*physicalX*physicalX);




		if (getIndex(x,y-1,z,_index[0]) )
			_velocity[0] = m_velocityField[_index[0]];
		else
			_velocity[0] = Eigen::Vector3f::Zero(3);

		if ( getIndex(x,y+1,z,_index[1]) )
			_velocity[1] = m_velocityField[_index[1]];
		else
			_velocity[1] = Eigen::Vector3f::Zero(3);

		if ( getIndex(x,y,z,_index[2]) )
			_velocity[2] = m_velocityField[_index[2]];
		else
			_velocity[2] = Eigen::Vector3f::Zero(3);
		_laplacian += (_velocity[0]+_velocity[1]-2.0f*_velocity[2] ) / (4.0f*physicalX*physicalX);



		if (getIndex(x,y,z-1,_index[0]) )
			_velocity[0] = m_velocityField[_index[0]];
		else
			_velocity[0] = Eigen::Vector3f::Zero(3);

		if ( getIndex(x,y,z+1,_index[1]) )
			_velocity[1] = m_velocityField[_index[1]];
		else
			_velocity[1] = Eigen::Vector3f::Zero(3);

		if ( getIndex(x,y,z,_index[2]) )
			_velocity[2] = m_velocityField[_index[2]];
		else
			_velocity[2] = Eigen::Vector3f::Zero(3);
		_laplacian += (_velocity[0]+_velocity[1]-2.0f*_velocity[2] ) / (4.0f*physicalZ*physicalZ);
	}

	return _laplacian;
}

bool Reconstruction::approximateVelocityVector(const int &x, const int &y, const int &z, Eigen::VectorXf &velocity, float &error, float &conditionNumber, int &measurments)
{
	int _index = z*m_dimX*m_dimY+y*m_dimX+x;

	if ( _index < m_measurementField.size() && m_measurementField[_index]->normalMatrix.rows() >= 3 )
	{

		

		VelocityData *_data = m_measurementField[_index];

		measurments = _data->normalMatrix.rows();

		//Eigen::LUDecomposition<Eigen::MatrixXf> lu(_data->normalMatrix);

		//_data->normalMatrix.lu().rank();
		Eigen::FullPivLU<Eigen::MatrixXf> _normalMatrixLU(_data->normalMatrix);
		// insufficient rank condition
		
		//std::cout << "Orientations: " << std::endl << _data->normalMatrix << std::endl << std::endl;

		/*std::cout << "Velocities: " << std::endl << _data->velocityVector << std::endl << std::endl;

		std::cout << "Frames: " << std::endl << _data->frameVector << std::endl << std::endl;
*/

		if ( _normalMatrixLU.rank() < 3 )
			return false;/**/

		/*Eigen::MatrixXf D = _data->normalMatrix;
		Eigen::MatrixXf A = D.transpose() * D;
		Eigen::VectorXf b = D.transpose() * _data->velocityVector;

		A.llt().solve(b,&velocity);   // using a LLT factorization
		A.ldlt().solve(b,&velocity);  // using a LDLT factorization
*/

		velocity = _data->normalMatrix.svd().solve(_data->velocityVector);
		Eigen::VectorXf _testSV = _data->normalMatrix.svd().singularValues();

		conditionNumber = _testSV(0) / _testSV(_testSV.rows()-1);


		//std::cout << _testSV.transpose() << std::endl;



		/*if ( x == m_dimX/2 && y == m_dimY/2 && z == m_dimZ/2)
		{
			std::cout << "A: " << std::endl << _data->normalMatrix << std::endl << std::endl;

			std::cout << "b: " << std::endl << _data->velocityVector << std::endl << std::endl;

			std::cout << "Solution: " << std::endl << velocity << std::endl;
		}*/

		error = (_data->normalMatrix * velocity - _data->velocityVector).norm();

		_data->approximationError = error;
		/*if ( error == error )
		{*/

		/*std::cout << "Frames: " << std::endl << _data->frameVector << std::endl << std::endl;

		
		std::cout << "Norm: " << (_data->normalMatrix * velocity - _data->velocityVector).norm() << std::endl;

		std::cout << "Solution: " << std::endl << velocity << std::endl << std::endl << std::endl ;
		*/
		//}
		return (error == error);
		
		//char c;
		//std::cin >> c;
		//return true;
	}
	else
		return false;
}


bool Reconstruction::getMeasurementData(const Eigen::Vector3f &worldPosition, Eigen::MatrixXf &normalMatrix, Eigen::VectorXf &velocityVector)
{
	int _index;

	if ( getIndex(worldPosition, _index) )
	{
		VelocityData *_data = m_measurementField[_index];
		

		normalMatrix = _data->normalMatrix;
		velocityVector = _data->velocityVector;
		
		return true;
	}
	return false;
}


bool Reconstruction::pushMeasurement( const Eigen::Vector3f &normal, const float &velocity, const int &_index, const int &frame)
{

	/*int _index;

	if ( getIndex(worldPosition, _index) )
	{*/
		VelocityData *_data = m_measurementField[_index];

		if ( _data->normalMatrix.rows() > 0 && _data->normalMatrix.cols() > 0 )
			_data->normalMatrix.conservativeResize(_data->normalMatrix.rows()+1, _data->normalMatrix.cols());
		else
			_data->normalMatrix = Eigen::MatrixXf::Zero(1,3);

		_data->normalMatrix(_data->normalMatrix.rows()-1,0) = normal(0);
		_data->normalMatrix(_data->normalMatrix.rows()-1,1) = normal(1);
		_data->normalMatrix(_data->normalMatrix.rows()-1,2) = normal(2);

		typedef Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> MatrixType;

		MatrixType m;
		m.setZero(1,3);
		m.conservativeResize(3,5);

		if ( _data->velocityVector.size() > 0 )
			_data->velocityVector.conservativeResize(_data->velocityVector.size()+1);
		else
			_data->velocityVector = Eigen::VectorXf(1);

		_data->velocityVector(_data->velocityVector.size()-1) = velocity;


		if ( _data->frameVector.size() > 0 )
			_data->frameVector.conservativeResize(_data->frameVector.size()+1);
		else
		_data->frameVector = Eigen::VectorXi(1);

		_data->frameVector(_data->frameVector.size()-1) = frame;
		//(data);
		return true;
	/*}
	else
		return false;*/
}