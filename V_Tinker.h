/*
 * V_Tinker.h
 *
 *  Created on: Sep 7, 2012
 *      Author: dstuck
 */

#ifndef V_TINKER_H_
#define V_TINKER_H_

#include "CallTinker.h"
#include "CoordUtil.h"
#include "debug.h"
#include "Potential.h"
#include "Propagator.h"
#include "Particle.h"
#include "Rho_Free.h"
#include <math.h>
#include <iostream>
#include <vector>
using namespace std;

class V_Tinker: public Potential {
public:
	V_Tinker(CoordUtil*, double, double);
	virtual ~V_Tinker();
	double GetV(vector<Particle>, Propagator *);
	double GetV(vector<Particle>);
	string GetType();
        CoordUtil* GetCoordUtil();
	void Tokenize(const string&, vector<string>&, const string& = " ");

	double vEquib;        //in kcal/mole
	std::string tinkPrmFileName;
	CoordUtil* coordKeeper;
        vector<double> wTinker;
};

#endif /* V_TINKER_H_ */
