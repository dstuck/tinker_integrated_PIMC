/*
 * VAQO.h
 *
 *  Created on: May 24, 2012
 *      Author: dstuck
 */

#ifndef VAQO_H_
#define VAQO_H_

#include "debug.h"
#include "CoordUtil.h"
#include "Potential.h"
#include "Propagator.h"
#include "Particle.h"
#include <iostream>
#include <vector>
using namespace std;

class V_AQO: public Potential {
public:
	V_AQO();
	V_AQO(double, double, int);
	V_AQO(vector<double>, vector<double>);
	virtual ~V_AQO();
	double GetV(vector<Particle>, Propagator *);
	string GetType();
        CoordUtil* GetCoordUtil();

	double V;
	vector<double> omega;
	vector<double> lambda;
};

#endif /* VAQO_H_ */
