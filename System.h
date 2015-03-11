/*
 * System.h
 *
 *  Created on: May 7, 2012
 *      Author: dstuck
 */

#ifndef SYSTEM_H_
#define SYSTEM_H_
#include "Particle.h"
#include <vector>
class PhysicsUtil;

class System {
public:
	System();
	virtual ~System();
//	virtual void CalcEnergy() = 0;		//Used internally
//	virtual void CalcPotential() = 0;	//Used internally
	virtual void Move(vector<double>, int, int levyPart = 1) = 0;
	virtual void Forget() = 0;
	virtual void Undo() = 0;
    virtual void Reset() = 0;
	virtual double GetWeight() = 0;
	virtual double EstimatorE() = 0;
	virtual double EstimatorV() = 0;
	virtual double GetOldWeight() = 0;
	virtual vector< vector<Particle> > GetParticle() = 0;
	virtual double Debug() = 0;
	virtual string GetVType() = 0;
	virtual string GetRhoType() = 0;
    virtual PhysicsUtil* GetPhysics() = 0;

};

#endif /* SYSTEM_H_ */
