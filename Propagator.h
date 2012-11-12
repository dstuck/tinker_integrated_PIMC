/*
 * Propagator.h
 *
 *  Created on: May 9, 2012
 *      Author: dstuck
 */

#ifndef PROPAGATOR_H_
#define PROPAGATOR_H_

#include "Particle.h"

class Propagator {
public:
	Propagator();
	virtual ~Propagator();

	virtual double GetRho(vector<Particle>,vector<Particle>, double) = 0;
	virtual double ModifyPotential(vector<Particle>) = 0;
	virtual double Estimate(vector<Particle>,vector<Particle>, double, int) = 0;
	virtual vector<double> GetSpringLength(vector<Particle>, double) = 0;		//	TODO: Remove GetSpringLength
	virtual vector<double> GetLevyMean(Particle, Particle, double, double, int) = 0;
	virtual double GetLevySigma(double, double, int) = 0;
	virtual string GetType() = 0;
};

#endif /* PROPAGATOR_H_ */
