/*
 * ClassicalSys.h
 *
 *  Created on: May 7, 2012
 *      Author: dstuck
 */

#ifndef CLASSICALSYS_H_
#define CLASSICALSYS_H_

#include "debug.h"
#include "System.h"
#include "Potential.h"
#include "V_UCHO.h"
#include "V_AQO.h"
#include "V_Morse.h"
#include "V_TinkerExecutable.h"
#include "V_Tinker.h"
#include "Propagator.h"
#include "Rho_Free.h"
#include "Rho_HO.h"
#include "Rho_TruncHO.h"
#include "Random.h"
#include "CoordUtil.h"
#include "PhysicsUtil.h"
#include <math.h>
#include <iostream>
#include <vector>
#include <fstream>
using namespace std;

class TestSystem: public System {
public:
	TestSystem(int, double, CoordUtil*, PhysicsUtil *);
	virtual ~TestSystem();
	void CalcEnergy();
	void CalcPotential();
	void Move(vector<double>, int);
	void Forget();
	void Undo();
	void Reset();
	double GetWeight();
	double EstimatorV();
	double EstimatorE();
	double GetOldWeight();
	vector< vector<Particle> > GetParticle();
	string GetVType();
	string GetRhoType();
	PhysicsUtil* GetPhysics();
	double Debug();

	bool ECheckFlag;
	int P;
	int N;
	int coorDim;
	int levyLength;
	double eps;
	double weight;
	double oldWeight;
	double energy;		//TODO: Phase out?
	double oldEnergy;	//TODO: Phase out?
	double potE;		//Warning: this is the sum of V over all the beads (needs to be divided by P to calc avg)
	double oldPotE;
	double avgV;
	int numSteps;
	Potential * V;
	Propagator * rho;
        PhysicsUtil * physics;
	vector<bool> upToDate;
	vector<double> sliceV;
	vector<double> oldSliceV;
	vector< vector<Particle> > part;
	vector< vector<Particle> > oldPart;

//DES Temp
        int tempNum;     //TODO: Delete this
        ofstream qchemFile;  //TODO: Delete this
};

#endif /* CLASSICALSYS_H_ */
