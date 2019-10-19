#include "NumericTools.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

NumericTools::NumericTools()
{
	// initialize random number generator
	srand((unsigned)time(NULL));
}

unsigned char NumericTools::quantize8BIT(float d, float max)
{
    return (unsigned char)((d / max) * 255.0);
}

float NumericTools::dequantize8BIT(unsigned char d, float max)
{
	return ((float)d)*max/255.0;
}

float NumericTools::dequantize16BIT(unsigned short d, float max)
{
	return ((float)d)*max/65535.0;
}

unsigned short NumericTools::quantize16BIT(float d, float max)
{
    return (unsigned short)((d / max) * 65535.0);
}

double NumericTools::generateUniformRand()
{
	return ((double) rand() / (RAND_MAX+1));
}
double NumericTools::generateStandardDistributionRand(double &mu, double &sigma)
{
 double u, v, x, y, q;

	do {
			u = ((double) rand() / (RAND_MAX+1)) ;
			v = 1.7156*(((double) rand() / (RAND_MAX+1)) -0.5);
			x = u - 0.449871;
			y = fabs(v) + 0.386595;
			q = pow(x,2.0) + y*(0.19600*y-0.25472*x);
		} while (q > 0.27597&& (q > 0.27846 || pow(v,2.0) > -4.*log(u)*pow(u,2.0)));
		return mu + sigma*v/u;

}



void NumericTools::variance( std::vector<float> arr, float *var, float *avg)
{
    int i;
    float sum = 0.0, sum2 = 0.0, tavg;

    for (i = 0; i < arr.size(); i++)
	sum += arr[i];
    tavg = sum / (float)arr.size();;

	for (i = 0; i < arr.size(); i++)
	sum2 += (tavg - arr[i]) * (tavg - arr[i]);

    *avg = tavg;
	*var = sum2 / (float) (arr.size() - 1);
}

void NumericTools::variance( std::vector<float> arr, float *var, float *avg, float *maxV, float *minV)
{
	* minV = std::numeric_limits<float>::max();
	*maxV = std::numeric_limits<float>::min();
    int i;
    float sum = 0.0, sum2 = 0.0, tavg;

    for (i = 0; i < arr.size(); i++)
	{
		sum += arr[i];
		if ( arr[i] < *minV ) *minV = arr[i];
		if ( arr[i] > *maxV ) *maxV = arr[i];
	}
	
    tavg = sum / (float)arr.size();;

	for (i = 0; i < arr.size(); i++)
	sum2 += (tavg - arr[i]) * (tavg - arr[i]);

    *avg = tavg;
	*var = sum2 / (float) (arr.size() - 1);
}
