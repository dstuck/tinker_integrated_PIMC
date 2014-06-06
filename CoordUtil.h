/*
 * CoordUtil.h
 *
 *  Created on: Sep 16, 2012
 *      Author: dstuck
 */

#ifndef COORDUTIL_H_
#define COORDUTIL_H_

#include "debug.h"
#include "Particle.h"
#include <math.h>
#include <iostream>
#include <vector>
using namespace std;

class CoordUtil {
public:
	CoordUtil();
	CoordUtil(int, int, vector< vector< vector<double> > >, vector<double>, vector<double>, vector< vector<double> >, vector<string>, vector<int>, vector< vector<int> >, string, string,bool);
	virtual ~CoordUtil();

	vector< vector<double> > normalModeToCart(vector<Particle>);

        bool readOmega;
	int numModes;
	int numPart;
	int dim;
	std::string tinkerName;
	std::string prmName;
	vector<string> atomType;
	vector<int> paramType;
//        vector<double> mass;
        vector<double> reducedMass;
	vector<double> omega;
	vector< vector<int> > connectivity;
	vector< vector<double> > initCart;
	vector< vector< vector<double> > > normModes;
};

#endif /* COORDUTIL_H_ */
