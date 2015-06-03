/*
 * V_Potlib.cpp
 *
 *  Created on: Sep 7, 2012
 *      Author: dstuck
 */

#include "V_Potlib.h"


V_Potlib::V_Potlib(CoordUtil* coords, double eps, double beta) : coordKeeper(coords) {

// For now just read in frequency and normaal modes

//    Get vEquib
   double hartreeToKcal = 627.509469;						//From http://en.wikipedia.org/wiki/Hartree 8/17/2012
   vector<Particle> parts;
   Particle p0(0.0, 1);
   for(int i=0; i<coordKeeper->numModes; i++) {
      parts.push_back(p0);
   }
   vEquib = 0.0;
   vEquib = GetV(parts);
}

V_Potlib::~V_Potlib() {
   delete coordKeeper;
}

double V_Potlib::GetV(vector<Particle> part, Propagator * rho){

   double V = 0.0;
   V += GetV(part);
//	cout << "Unmodified V =\t" << V << endl;
   V += rho->ModifyPotential(part);
//	cout << "Modified V =\t" << V << endl;
   return V;
}

double V_Potlib::GetV(vector<Particle> part){

   //double hartreeToKcal = 627.509469;						//From http://en.wikipedia.org/wiki/Hartree 8/17/2012
   double V = 0.0;
//	Get cartesians from normal modes
   vector< vector<double> > cartPos = coordKeeper->normalModeToCart(part);

// Call Potlib
   double potEnergy = 0.0; 
   int nPart = coordKeeper->numPart;
   double xyzCoord[nPart*3];
   for(int i=0; i<nPart; i++) {
      //typeVec[i] = coordKeeper->connectivity[i][0];
      for(int k=0; k<3; k++) {
         xyzCoord[k+i*3] = cartPos[i][k];
      }
   }
   //dstuckenergy_(nPart, typeVec, xyzCoord, connMat, tinkPrm.c_str(), prmLen, tinkEnergy, init);
   pot_(potEnergy,xyzCoord);
   V = potEnergy;
//	cout << "Tinker V =\t" << V << endl;
   V -= vEquib;
//   V /= hartreeToKcal;
   return V;
}

string V_Potlib::GetType() {
   string name = "Potlib";
   return name;
}

CoordUtil* V_Potlib::GetCoordUtil() {
   return coordKeeper;
}
