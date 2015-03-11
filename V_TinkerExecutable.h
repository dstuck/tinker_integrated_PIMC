/*
 * V_TinkerExecutable.h
 *
 *  Created on: Sep 7, 2012
 *      Author: dstuck
 */

#ifndef V_TINKEREXECUTABLE_H_
#define V_TINKEREXECUTABLE_H_

#include "CoordUtil.h"
#include "debug.h"
#include "Potential.h"
#include "Propagator.h"
#include "Particle.h"
#include "Rho_Free.h"
#include <math.h>
#include <iostream>
#include <vector>
#include <fstream>
using namespace std;

class V_TinkerExecutable: public Potential {
public:
	V_TinkerExecutable(CoordUtil*);
	virtual ~V_TinkerExecutable();
	double GetV(vector<Particle>, Propagator *);
	double GetV(vector<Particle>);
	string GetType();
        CoordUtil* GetCoordUtil();
	void Tokenize(const string&, vector<string>&, const string& = " ");

	double vEquib;
	std::string tinkInFileName;
	std::string tinkOutFileName;
	std::string tinkPrmFileName;
	ofstream inFile;
	ifstream outFile;
	CoordUtil * coordKeeper;
};

#endif /* V_TINKEREXECUTABLE_H_ */
