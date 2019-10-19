#ifndef INTENSITYCOMPOUNDING_H
#define INTENSITYCOMPOUNDING_H


#include <vector>

/**
 * \brief 3D Reconstruction software of ultrasound data using RF-data
 *
 * \authors Tassilo Klein
 *          <br>
 *          kleint\@cs.tum.edu
 * \ingroup Ultrasound
 * \version 1.0a
 * \date 07-06-2010
 *
 * \par License:
 * Copyright (c) 2005, 2006, 2007, 2008, 2009
 * Chair for Computer Aided Medical Procedures & Augmented Reality (CAMPAR) I-16,
 * Technische Universit&auml;t M&uuml;nchen
 * <br>
 * <br>
 * All rights reserved.
 * <br>
 * <br>
 * See <a href="COPYRIGHT.txt">COPYRIGHT.txt</a> for details.
 * <br>
 * <br>
 * This software is distributed WITHOUT ANY WARRANTY; without even 
 * <br>
 * the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
 * <br>
 * PURPOSE.  See the <a href="COPYRIGHT.txt">COPYRIGHT.txt</a> notice
 * for more information.
 *
 */


class IntensityCompounding
{
public:
	IntensityCompounding() {}
	virtual float evaluate(std::vector<float> *intensityVector, std::vector<float> *distanceVector) { return 0; }
};

class GaussianWeightedCompounding : public IntensityCompounding
{
private:
	float m_sigma;
	float m_sigmaSQR;
public:
	GaussianWeightedCompounding(float sigma);
	virtual float evaluate(std::vector<float> *intensityVector, std::vector<float> *distanceVector);

};

class MeanCompounding : public IntensityCompounding
{
public:
	MeanCompounding();
	virtual float evaluate(std::vector<float> *intensityVector, std::vector<float> *distanceVector);
};


#endif