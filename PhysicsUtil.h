/*
 * PhysicsUtil.h
 *
 *  Created on: Sep 21, 2012
 *      Author: dstuck
 */

#ifndef PHYSICSUTIL_H_
#define PHYSICSUTIL_H_

#include "debug.h"
#include <iostream>
#include <vector>
#include <string>
using namespace std;

class PhysicsUtil {
public:
	PhysicsUtil();
	virtual ~PhysicsUtil();
        bool isDeltaAI();

//	TODO: Make this a struct
//	Propagator variables
   string rhoType;
   double freqCutoff;
//	Potential Variables
   string vType;
   bool deltaAbInit;
   int charge;
   int multiplicity;
   int numInit;
   int lowFrozModes;
   int highFrozModes;
   double morseDE;
   double morseAlpha;
   double morseMass;
   double lambdaTI;
};

#endif /* PHYSICSUTIL_H_ */
