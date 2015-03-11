/*
 * V_QChem.h
 *
 *  Created on: June 23, 2014
 *      Author: dstuck
 */

#ifndef V_QCHEM_H_
#define V_QCHEM_H_

#include "CallTinker.h"
#include "CoordUtil.h"
#include "debug.h"
#include "Potential.h"
#include "Propagator.h"
#include "Particle.h"
#include "Rho_Free.h"
#include "V_UCHO.h"
#include <math.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
//#include <algorithm>
using namespace std;

class V_QChem: public Potential {
public:
   V_QChem(CoordUtil*, CoordUtil*, double, double, int, int);
   virtual ~V_QChem();
   double GetV(vector<Particle>, Propagator *);
   double GetV(vector<Particle>);
   string GetType();
   CoordUtil* GetCoordUtil();
   double GetVHO(vector<Particle>);
   void Tokenize(const string&, vector<string>&, const string& = " ");

   bool doubleTI;
   double vEquib;        //in kcal/mole
   CoordUtil* coordKeeper;
   CoordUtil* tinkerCoords;
   vector<double> wQChem;
   ofstream inFile;
   V_UCHO harmV;
   int charge;
   int spin;

//DES Temp
   int tempNum;     //TODO: Delete this
   ofstream qchemFile;  //TODO: Delete this
   ofstream qharmFile;  //TODO: Delete this
};

#endif /* V_QCHEM_H_ */
