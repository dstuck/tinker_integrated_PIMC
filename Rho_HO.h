/*
 * RhoHO.h
 *
 *  Created on: May 17, 2012
 *      Author: dstuck
 */

#ifndef RHOHO_H_
#define RHOHO_H_

#include <math.h>
#include <iostream>
#include "Propagator.h"
#include "debug.h"
using namespace std;

class Rho_HO: public Propagator {
public:
	Rho_HO(vector<double>, int, int=0);
	Rho_HO(double, int);
	virtual ~Rho_HO();

	double GetRho(vector<Particle>, vector<Particle>, double);
	double ModifyPotential(vector<Particle>);
	double Estimate(vector<Particle>, vector<Particle>, double, int);
	vector<double> GetSpringLength(vector<Particle>, double);
	vector<double> GetLevyMean(Particle, Particle, double, double, int);
	double GetLevySigma(double, double, int);
	string GetType();

	vector<double> omega;
        int lowFrozModes;
        int highFrozModes;
};

#endif /* RHOHO_H_ */
