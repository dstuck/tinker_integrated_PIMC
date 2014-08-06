/*
 * V_Tinker.cpp
 *
 *  Created on: Sep 7, 2012
 *      Author: dstuck
 */

#include "V_Tinker.h"


V_Tinker::V_Tinker(CoordUtil* coords, double eps, double beta) : coordKeeper(coords) {
   tinkPrmFileName = coordKeeper->prmName;
//	coordKeeper = coords;
//    Initialize forcefield
   vector<Particle> parts;
   Particle p0(0.0, 1);
   for(int i=0; i<coordKeeper->numModes; i++) {
      parts.push_back(p0);
   }
//	Get cartesians from normal modes
   vector< vector<double> > cartPos = coordKeeper->initCart;
// DES Temp for CN only!:
   //cartPos[0][2] -= -0.00846;
   //cartPos[1][2] += 0.00846;
//bool init = true;
   int init = 1;
//      cout << "inLen is: " << inLen << endl;
   int prmLen = tinkPrmFileName.size();
   double tinkEnergy = 0.0; 
   int nPart = coordKeeper->numPart;
   int typeVec[nPart];
   double xyzCoord[nPart*3];
   double connMat[nPart*8];
   for(int i=0; i<nPart*8 ; i++){
      connMat[i]=0.0;
   }
// TODO: Oops, the first column of connectivity contains tinker atomtypes...
   for(int i=0; i<nPart; i++) {
      typeVec[i] = coordKeeper->connectivity[i][0];
      for(int k=0; k<3; k++) {
         xyzCoord[k+i*3] = cartPos[i][k];
      }
      for(int j=1; j<coordKeeper->connectivity[i].size(); j++) {
         connMat[j-1+i*8] = coordKeeper->connectivity[i][j];
      }
   }
   string tinkPrm = tinkPrmFileName;
   dstuckenergy_(nPart, typeVec, xyzCoord, connMat, tinkPrm.c_str(), prmLen, tinkEnergy, init);
//    Set up tinker normal modes and frequencies!
   int N = coordKeeper->numModes;
   if(!coordKeeper->readOmega) {
      double tinkModes[N*coordKeeper->numPart*3];
      double tinkFreqs[N];
      double redMass[N];
      dstuckvibrate_(typeVec, xyzCoord, connMat, N, coordKeeper->numPart, redMass, tinkModes, tinkFreqs);
      coordKeeper->reducedMass.clear();
      for(int i=0; i<N; i++) {
         coordKeeper->reducedMass.push_back(redMass[i]*1822.8886);
//         cout << coordKeeper->reducedMass[i]/1822.8886 << endl;
      }
      coordKeeper->omega.clear();
      for(int i=0; i<N; i++) {
// Don't want frequencies == 0      TODO: double check this
         if(tinkFreqs[i] < 1.0) {
            cout << "Frequency of " << tinkFreqs[i] << " set to 0.09 cm-1" << endl;
            coordKeeper->omega.push_back(0.00000455633);
         }
         else {
            coordKeeper->omega.push_back(tinkFreqs[i]*0.00000455633);
         }
//         coordKeeper.omega[i] = tinkFreqs[i]*0.00000455633;
//         cout << tinkFreqs[i] << endl;
      }
      coordKeeper->normModes.clear();
      coordKeeper->normModes.resize(N);
      for(int i=0; i<N; i++) {
         coordKeeper->normModes[i].resize(coordKeeper->numPart);
         double vnorm = 0.0;
         double tempVal = 0.0;
         for(int j=0; j<coordKeeper->numPart; j++) {
            coordKeeper->normModes[i][j].resize(3,0.0);
            for(int k=0; k<3; k++) {
// TODO: remove vnorm
//                  tempVal = tinkModes[k+3*j+3*coordKeeper->numPart*i];
               coordKeeper->normModes[i][j][k] = tinkModes[k+3*j+3*coordKeeper->numPart*i];
//                  coordKeeper->normModes[i][j][k] = tempVal;
//cout << "DES Temp:" << tinkModes[k+3*j+3*coordKeeper->numPart*i] << endl;
            }
         }
//            for(int j=0; j<coordKeeper->numPart; j++) {
//               for(int k=0; k<3; k++) {
//                  coordKeeper->normModes[i][j][k] /= vnorm;
//               }
//            }
      }
   }
   else {
//cout << "DES Temp: Reading in frequencies from file" << endl;
   }

   double clHarmFull = 0.0;
   for(int i=0; i<N; i++) {
      clHarmFull += coordKeeper->omega[i]/2.0;
   }
   cout << "Harmonic Energy is: " << clHarmFull << endl;

//    Get vEquib
   double hartreeToKcal = 627.509469;						//From http://en.wikipedia.org/wiki/Hartree 8/17/2012
   vEquib = 0.0;
   vEquib = GetV(parts)*hartreeToKcal;
   /*
      for(int i=0; i<N; i++) {
      cout << coordKeeper->reducedMass[i]/1822.8886 << endl;
      cout << coordKeeper->omega[i]/0.00000455633 << endl;
      for(int j=0; j<coordKeeper->numPart; j++) {
      for(int k=0; k<3; k++) {
      cout << coordKeeper->normModes[i][j][k] << "\t";
      }
      cout << endl;
      }
      }
    */
}

V_Tinker::~V_Tinker() {
   delete coordKeeper;
}

double V_Tinker::GetV(vector<Particle> part, Propagator * rho){

   double V = 0.0;
   V += GetV(part);
//	cout << "Unmodified V =\t" << V << endl;
   V += rho->ModifyPotential(part);
//	cout << "Modified V =\t" << V << endl;
   return V;
}

double V_Tinker::GetV(vector<Particle> part){

   double hartreeToKcal = 627.509469;						//From http://en.wikipedia.org/wiki/Hartree 8/17/2012
   double V = 0.0;
//	Get cartesians from normal modes
   vector< vector<double> > cartPos = coordKeeper->normalModeToCart(part);
   /*
//DES Temp Temp:
   if(tempcounter==0||tempcounter==2) {
   tempcounter++;
   int N = coordKeeper->numPart;
   vector<double> normMode = coordKeeper->cartToNormalMode(cartPos);
   cout << "Normal modes" << endl;
   for(int i=0; i<part.size(); i++) {
   cout << part[i].pos[0] << "\t" << normMode[i] << endl;
   }
//      cout << "Energy 1 = " << GetV(part) << endl;
   cout << "Cartesians1" << endl;
   for(int i=0; i<N; i++) {
   cout << cartPos[i][0] << "\t" << cartPos[i][1] << "\t" << cartPos[i][2] << endl;
   }

   for(int i=0; i<part.size(); i++) {
   part[i].pos[0] = normMode[i];
   }
//      cout << "Energy 2 = " << GetV(part) << endl;
   cartPos = coordKeeper->normalModeToCart(part);
   cout << "Cartesians2" << endl;
   for(int i=0; i<N; i++) {
   cout << cartPos[i][0] << "\t" << cartPos[i][1] << "\t" << cartPos[i][2] << endl;
   }
   if(tempcounter==1) {
   tempcounter++;
   }
   else {
   exit(-1);
   }
   }
    */

// Call Tinker
//bool init = false;
   int init = 0;
   int prmLen = tinkPrmFileName.size();
   double tinkEnergy = 0.0; 
   int nPart = coordKeeper->numPart;
   int typeVec[nPart];
   double xyzCoord[nPart*3];
   double connMat[nPart*8];
   for(int i=0; i<nPart*8 ; i++){
      connMat[i]=0.0;
   }
// DES Temp for CN only!:
   //cartPos[0][2] += -0.00846;
   //cartPos[1][2] += 0.00846;
// TODO: Oops, the first column of connectivity contains tinker atomtypes...
   for(int i=0; i<nPart; i++) {
      typeVec[i] = coordKeeper->connectivity[i][0];
      for(int k=0; k<3; k++) {
         xyzCoord[k+i*3] = cartPos[i][k];
      }
      for(int j=1; j<coordKeeper->connectivity[i].size(); j++) {
         connMat[j-1+i*8] = coordKeeper->connectivity[i][j];
      }
   }
   string tinkPrm = tinkPrmFileName;
   dstuckenergy_(nPart, typeVec, xyzCoord, connMat, tinkPrm.c_str(), prmLen, tinkEnergy, init);
   V = tinkEnergy;
//	cout << "Tinker V =\t" << V << endl;
   V -= vEquib;
   V /= hartreeToKcal;
   return V;
}

string V_Tinker::GetType() {
   string name = "Tinker";
   return name;
}

CoordUtil* V_Tinker::GetCoordUtil() {
   return coordKeeper;
}
