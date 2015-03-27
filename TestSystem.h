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
#include "V_QChem.h"
#include "Propagator.h"
#include "Rho_Free.h"
#include "Rho_HO.h"
#include "Rho_TruncHO.h"
#include "Random.h"
#include "CoordUtil.h"
#include "PhysicsUtil.h"
#include <math.h>
#include <iostream>
#include <cstdio>
#include <vector>
#include <sstream>
#include <fstream>
using namespace std;

class TestSystem: public System {
public:
	TestSystem(int, double, CoordUtil*, PhysicsUtil *);
	virtual ~TestSystem();
	void CalcEnergy();
	void CalcPotential();
	void Move(vector<double> prob, int levyNum, int levyPart = 1);
	void Forget();
	void Undo();
	void Reset();
	double EstimatorV();
	double EstimatorE();
	double GetWeight();
	double GetOldWeight();
	double GetHarmonicE();
	vector< vector<Particle> > GetParticle();
	string GetVType();
	string GetRhoType();
	PhysicsUtil* GetPhysics();
	double Debug();
        string toString(int);

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
        double harmonicE;
	int numSteps;
	Potential * V;
	Potential * V2;
	Propagator * rho;
        PhysicsUtil * physics;
	vector<bool> upToDate;
	vector<double> sliceV;
	vector<double> oldSliceV;
	vector< vector<Particle> > part;
	vector< vector<Particle> > oldPart;
	vector< vector< vector<double> > > oldCarts;
	vector< vector< vector<double> > > newCarts;

//DES Temp
      ofstream vFile;      //TODO: Remove
      int tempNum;      //TODO: Remove
//        ofstream qchemFile;  //TODO: Delete this
};

#endif /* CLASSICALSYS_H_ */
