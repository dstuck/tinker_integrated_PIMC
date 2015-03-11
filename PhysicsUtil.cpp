/*
 * PhysicsUtil.cpp
 *
 *  Created on: Sep 21, 2012
 *      Author: dstuck
 */

#include "PhysicsUtil.h"

PhysicsUtil::PhysicsUtil() {
//	Makes sure it gets set
	numInit = -1;
   lowFrozModes = 0;
   highFrozModes = 0;
   lambdaTI = 1.0;
   deltaAbInit = false;
   charge = 0;
   multiplicity = 1;
}

PhysicsUtil::~PhysicsUtil() {
	// TODO Auto-generated destructor stub
}

bool PhysicsUtil::isDeltaAI() {
    return deltaAbInit;
}
