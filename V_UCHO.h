/*
 * V_UCHO.h
 *
 *  Created on: May 9, 2012
 *      Author: dstuck
 */

#ifndef V_UCHO_H_
#define V_UCHO_H_

#include "debug.h"
#include "Potential.h"
#include "Propagator.h"
#include "Particle.h"
#include <iostream>
#include <vector>
using namespace std;


class V_UCHO: public Potential {
public:
	V_UCHO();
	V_UCHO(vector<double>);
	V_UCHO(double, int);
	virtual ~V_UCHO();
	double GetV(vector<Particle>, Propagator *);
	string GetType();

	double V;
	vector<double> omega;
};

#endif /* V_UCHO_H_ */
