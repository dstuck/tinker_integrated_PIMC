/*
 * Potential.h
 *
 *  Created on: May 9, 2012
 *      Author: dstuck
 */

#ifndef POTENTIAL_H_
#define POTENTIAL_H_

#include "Particle.h"
#include "Propagator.h"
#include <vector>

class CoordUtil;

using namespace std;

class Potential {
public:
	Potential();
	virtual ~Potential();

	virtual double GetV(vector<Particle>, Propagator *) = 0;
	virtual string GetType() = 0;
        virtual CoordUtil* GetCoordUtil() = 0;
};

#endif /* POTENTIAL_H_ */
