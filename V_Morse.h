/*
 * V_Morse.h
 *
 *  Created on: May 24, 2012
 *      Author: dstuck
 */

#ifndef V_MORSE_H_
#define V_MORSE_H_

#include "debug.h"
#include "Potential.h"
#include "Propagator.h"
#include "Particle.h"
#include <math.h>
#include <iostream>
#include <vector>
using namespace std;

class V_Morse: public Potential {
public:
	V_Morse();
	V_Morse(double, double, int);
	V_Morse(vector<double>, vector<double>);
	virtual ~V_Morse();
	double GetV(vector<Particle>, Propagator *);
	string GetType();

	vector<double> de;
	vector<double> a;
};

#endif /* V_MORSE_H_ */
