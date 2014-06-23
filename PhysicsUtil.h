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
//	Potential Variables
	string vType;
	double morseDE;
	double morseAlpha;
	int numInit;
    int numFrozModes;
    double lambdaTI;
    bool deltaAbInit;
};

#endif /* PHYSICSUTIL_H_ */
