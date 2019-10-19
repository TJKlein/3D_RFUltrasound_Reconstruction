#include "RandomNumbers.h"
#include <iostream>

RandomNumbers::RandomNumbers() :
m_numericTools(0)
{
	m_numericTools = new NumericTools();
}

RandomNumbers::~RandomNumbers()
{
	if ( m_numericTools ) delete m_numericTools;
}


double RandomNumbers::univariateGaussian(double mu, double sigma)
{
	return m_numericTools->generateStandardDistributionRand(mu, sigma);
}

Eigen::Vector3f RandomNumbers::multivariateGaussian(Eigen::Vector3f &mu, Eigen::Matrix3f &covariance)
{
	// to generate multivariate gaussian random numbers given mean and covariance
	// 1. decompose the covariance with cholesky LU: cov = U'U  (size n x n)
	// 2. generate n normally distributed random numbers in a vector RND (size n x 1)
	// 3. generate the multivariate random number: mvrnd = mu + U' * RND

	Eigen::FullPivHouseholderQR<Eigen::Matrix3f> _qr(covariance);

	Eigen::LLT<Eigen::Matrix3f> _test;

	_test.compute(covariance);

	Eigen::Matrix3f _Uprime = Eigen::Matrix3f::Zero();
	_Uprime.block<3,3>(0,0).part<Eigen::Lower>() = _test.matrixL();

	//std::cout <<"L" << std::endl << _Uprime << std::endl;

	
	// random vector
	Eigen::Vector3f _randomVec = Eigen::Vector3f::Zero(3);

	double _stdVal1 = 0.0;
	double _stdVal2 = 1.0;
	_randomVec << m_numericTools->generateStandardDistributionRand(_stdVal1, _stdVal2), m_numericTools->generateStandardDistributionRand(_stdVal1, _stdVal2), m_numericTools->generateStandardDistributionRand(_stdVal1, _stdVal2);;
	return _Uprime*_randomVec+mu;	

	


	//Eigen::LU

}