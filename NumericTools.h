#ifndef NUMERICTOOLS_H
#define NUMERICTOOLS_H

#include <vector>

class NumericTools
{
private:
public:

	NumericTools();

	/**
	\brief Low-quality random number generation following the standard distribution
	**/
	double generateStandardDistributionRand(double &mu, double &sigma);

	double generateUniformRand();

	void variance( std::vector<float> arr, float *var, float *avg);
	void variance( std::vector<float> arr, float *var, float *avg, float *maxV, float *minV);
	unsigned char quantize8BIT(float d, float max);
	float dequantize8BIT(unsigned char d, float max);
	float dequantize16BIT(unsigned short d, float max);
	unsigned short quantize16BIT(float d, float max);
};



#endif