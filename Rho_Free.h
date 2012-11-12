/*
 * RhoFree.h
 *
 *  Created on: May 9, 2012
 *      Author: dstuck
 */

#ifndef RHOFREE_H_
#define RHOFREE_H_

#include <math.h>
#include <iostream>
#include "Propagator.h"
using namespace std;

class Rho_Free: public Propagator {
public:
	Rho_Free();
	virtual ~Rho_Free();

	double GetRho(vector<Particle>, vector<Particle>, double);
	double ModifyPotential(vector<Particle>);
	double Estimate(vector<Particle>, vector<Particle>, double, int);
	vector<double> GetSpringLength(vector<Particle>, double);
	vector<double> GetLevyMean(Particle, Particle, double, double, int);
	double GetLevySigma(double, double, int);
	string GetType();

};

#endif /* RHOFREE_H_ */
