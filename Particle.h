/*
 * Particles.h
 *
 *  Created on: May 7, 2012
 *      Author: dstuck
 */

#ifndef PARTICLES_H_
#define PARTICLES_H_
#include <string>
#include <vector>
using namespace std;

class Particle {
public:
	Particle();
	Particle(double, int);
	Particle(vector<double>);
	Particle(vector<double>, double);
	virtual ~Particle();

	vector<double> pos;
	string atomType;
	double mass;
	double charge;
};

#endif /* PARTICLES_H_ */
