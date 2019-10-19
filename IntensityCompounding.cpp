#include "IntensityCompounding.h"
#include <stdio.h>
#include <iostream>

GaussianWeightedCompounding::GaussianWeightedCompounding(float sigma) 
{
	m_sigma = sigma;
	m_sigmaSQR = pow(sigma, 2.0f);
}

float GaussianWeightedCompounding::evaluate(std::vector<float> *intensityVector, std::vector<float> *distanceVector)
{
	int _size = intensityVector->size();

	float _intensityDenominator = 0.0f;
	float _intensityNominator = 0.0;


	 std::vector<float>::iterator _distance = distanceVector->begin();
	for(std::vector<float>::iterator _intensity = intensityVector->begin(); _intensity !=intensityVector->end();_intensity++, _distance++ )
	{
		float _euclideanDistanceSQR = pow((*_distance),2.0f);
		const float _tempRatio = exp(-(_euclideanDistanceSQR)/(m_sigmaSQR));


		// add the intensity contribution by weighting distance to the nominator and denominator
		// which after considering all relevant pixels will by division supply the gaussian interpolated intensity

		_intensityDenominator += _tempRatio;

		_intensityNominator += (*_intensity) * _tempRatio; 


	}

	if ( _intensityDenominator > 0.0 )
	 return _intensityNominator / _intensityDenominator;
	else
		return 0.0;
}

MeanCompounding::MeanCompounding()
{
	std::cout << "Mean Compounding ... " << std::endl;
}


float MeanCompounding::evaluate(std::vector<float> *intensityVector, std::vector<float> *distanceVector)
{
	int _size = intensityVector->size();

	float _intensityDenominator = 0.0f;
	float _intensityNominator = 0.0;


	 std::vector<float>::iterator _distance = distanceVector->begin();
	for(std::vector<float>::iterator _intensity = intensityVector->begin(); _intensity !=intensityVector->end();_intensity++, _distance++ )
	{
		


		// add the intensity contribution by weighting distance to the nominator and denominator
		// which after considering all relevant pixels will by division supply the gaussian interpolated intensity

		_intensityNominator += (*_intensity); 


	}

	if ( _size > 0 )
	 return _intensityNominator / static_cast<float>(_size);
	else
		return 0.0;
}