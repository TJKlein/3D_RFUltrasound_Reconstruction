#ifndef RANDOMNUMBERS_H
#define RANDOMNUMBERS_H

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Array>
#include <Eigen/QR>
#include <Eigen/LU>
#include <Eigen/SVD>

#include "NumericTools.h"


class RandomNumbers
{
	private:
		NumericTools *m_numericTools;
	public:
		RandomNumbers();
		~RandomNumbers();
		Eigen::Vector3f multivariateGaussian(Eigen::Vector3f &mu, Eigen::Matrix3f &covariance);
		double univariateGaussian(double mu, double sigma);
};

#endif