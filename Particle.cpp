/*
 * Particles.cpp
 *
 *  Created on: May 7, 2012
 *      Author: dstuck
 */

#include "Particle.h"
#include <string.h>
#include <iostream>
using namespace std;

Particle::Particle() {}

Particle::Particle(double x, int dim) {
	for(int i=0; i<dim; i++) {
		pos.push_back(x);
	}
}

Particle::Particle(vector<double> x) {
	pos = x;
}

Particle::Particle(vector<double> x, double m) {
	pos = x;
	mass = m;
}

Particle::~Particle() {
	// TODO Auto-generated destructor stub
}

