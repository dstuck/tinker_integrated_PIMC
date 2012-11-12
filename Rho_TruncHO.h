/*
 * Rho_TruncHO.h
 *
 *  Created on: Jun 13, 2012
 *      Author: dstuck
 */

#ifndef RHO_TRUNCHO_H_
#define RHO_TRUNCHO_H_

#include <math.h>
#include <iostream>
#include "Propagator.h"
#include "debug.h"
using namespace std;

class Rho_TruncHO: public Propagator {
public:
	Rho_TruncHO(vector<double>);
	Rho_TruncHO(double, int);
	virtual ~Rho_TruncHO();

	double GetRho(vector<Particle>, vector<Particle>, double);
	double ModifyPotential(vector<Particle>);
	double Estimate(vector<Particle>, vector<Particle>, double, int);
	vector<double> GetSpringLength(vector<Particle>, double);
	vector<double> GetLevyMean(Particle, Particle, double, double, int);
	double GetLevySigma(double, double, int);
	string GetType();

	vector<double> omega;
};

#endif /* RHO_TRUNCHO_H_ */
