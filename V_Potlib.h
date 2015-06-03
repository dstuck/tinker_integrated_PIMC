/*
 * V_Potlib.h
 *
 *  Created on: Sep 7, 2012
 *      Author: dstuck
 */

#ifndef V_POTLIB_H_
#define V_POTLIB_H_

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

class V_Potlib: public Potential {
public:
	V_Potlib(CoordUtil*, double, double);
	virtual ~V_Potlib();
	double GetV(vector<Particle>, Propagator *);
	double GetV(vector<Particle>);
	string GetType();
        CoordUtil* GetCoordUtil();
	void Tokenize(const string&, vector<string>&, const string& = " ");

	double vEquib;        //in kcal/mole
	//std::string tinkPrmFileName;
	CoordUtil* coordKeeper;
        vector<double> wTinker;
};

#endif /* V_POTLIB_H_ */
