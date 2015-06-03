/*
 * RhoHO.h
 *
 *  Created on: May 21, 2015
 *      Author: dstuck
 */

#ifndef RHO_MIXED_H_
#define RHO_MIXED_H_

#include <math.h>
#include <iostream>
#include "Propagator.h"
#include "debug.h"
using namespace std;

class Rho_Mixed: public Propagator {
public:
	Rho_Mixed(vector<double>, double, int, int=0);
	virtual ~Rho_Mixed();

	double GetRho(vector<Particle>, vector<Particle>, double);
	double ModifyPotential(vector<Particle>);
	double Estimate(vector<Particle>, vector<Particle>, double, int);
	vector<double> GetSpringLength(vector<Particle>, double);
	vector<double> GetLevyMean(Particle, Particle, double, double, int);
	double GetLevySigma(double, double, int);
	string GetType();

        int lowFrozModes;
        int highFrozModes;
        int cutoffIndex;
        double cutoffVal;
        vector<double> mass;
	vector<double> omega;
};

#endif /* RHO_MIXED_H_ */
