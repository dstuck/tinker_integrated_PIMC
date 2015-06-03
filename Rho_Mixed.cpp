/*
 * RhoHO.cpp
 *
 *  Created on: May 21, 2015
 *      Author: dstuck
 */

#include "Rho_Mixed.h"

Rho_Mixed::Rho_Mixed(vector<double> w, double cutVal, int nFrozModes, int hFrozModes) : omega(w), cutoffVal(cutVal), lowFrozModes(nFrozModes), highFrozModes(hFrozModes) {
   cutoffIndex = omega.size();
   for(int i=0; i<omega.size(); i++) {
      if(omega[i] > cutoffVal) {
         cutoffIndex = i;
         break;
      }
   }
   cout << "DES: cutoffIndex = " << cutoffIndex << endl;
}

Rho_Mixed::~Rho_Mixed() {
// TODO Auto-generated destructor stub
}

double Rho_Mixed::GetRho(vector<Particle> slice1, vector<Particle> slice2, double eps) {
// 	This function is not called if using Levy Flights
   double rho = 0;
   if(omega.size() != slice1.size()) {
      cout << "Error! Number of particles not equal to number of frequencies! Exiting with -99999." << endl;
      return -99999;
   }
//Free particle part
   for(int j=0; j<cutoffIndex; j++) {
      rho += slice1[j].mass/2.0/eps/eps*(slice1[j].pos[0]-slice2[j].pos[0])*(slice1[j].pos[0]-slice2[j].pos[0]);
   }
//Harmonic Oscillator part
   for(int j=cutoffIndex; j<(int)slice1.size(); j++) {
      rho += slice1[j].mass*omega[j]/2.0/eps * ((slice1[j].pos[0] - slice2[j].pos[0])*(slice1[j].pos[0] - slice2[j].pos[0])/sinh(omega[j]*eps) + tanh(eps*omega[j]/2.0)*(slice1[j].pos[0]*slice1[j].pos[0]+slice2[j].pos[0]*slice2[j].pos[0]));
   }
   return rho;
}

double Rho_Mixed::ModifyPotential(vector<Particle> part) {
   double dV = 0.0;
   if(part[0].pos.size() != 1) {
      cout << "DES Error: Must be using normal modes with mixed propagator" << endl;
      LINE
      exit(-1);
   }
   for(int j=cutoffIndex; j<(int)part.size(); j++) {
         dV -= (part[j].pos[0])*(part[j].pos[0])/2.0*omega[j]*omega[j]*part[j].mass;
   }
//	cout << "\tV_HO =\t" << dV << endl;
//	cout << dV << endl;
   return dV;
}

double Rho_Mixed::Estimate(vector<Particle> slice1, vector<Particle> slice2, double eps, int P) {
   /*
    * 	Estimator taken from Whitfield and Martyna, J. Chem. Phys. 126, 074104 (2007) with missing factor of 1/2 added to term 2
    */
   if(omega.size() != slice1.size()) {
      cout << omega.size() << " frequencies for " << slice1.size() << " modes." << endl;
      cout << "Error! Number of particles not equal to number of frequencies! Exiting with -99999." << endl;
      return -99999;
   }

   double est = 0;
//Free particle part
   double delta;
   double kin;
   for(int j=0; j<cutoffIndex; j++) {
      delta = -(slice1[j].pos[0]-slice2[j].pos[0])*(slice1[j].pos[0]-slice2[j].pos[0]);
      kin = 1.0/2.0/eps;
      est += (slice1[j].mass/2.0/eps/eps*(delta)+kin);
   }
//Harmonic part
   for(int j=cutoffIndex; j< lowFrozModes; j++) {
      for(int k=0; k<(int)slice1[j].pos.size(); k++) {
         est += omega[j]/2.0/tanh(eps*(double)P*omega[j]);
      }
   }
   for(int j=(int)slice1.size()-highFrozModes; j<(int)slice1.size(); j++) {
      for(int k=0; k<(int)slice1[j].pos.size(); k++) {
         est += omega[j]/2.0/tanh(eps*(double)P*omega[j]);
      }
   }
   for(int j=max(lowFrozModes,cutoffIndex); j<(int)slice1.size()-highFrozModes; j++) {
      for(int k=0; k<(int)slice1[j].pos.size(); k++) {
         est += omega[j]/2.0/tanh(eps*omega[j]);
         est += -slice1[j].mass*omega[j]*omega[j]/2.0/tanh(eps*omega[j])/sinh(eps*omega[j])*(slice1[j].pos[k]-slice2[j].pos[k])*(slice1[j].pos[k]-slice2[j].pos[k]);
         est += slice1[j].mass*omega[j]*omega[j]/2.0/cosh(eps*omega[j]/2.0)/cosh(eps*omega[j]/2.0)*(slice1[j].pos[k]*slice1[j].pos[k]);
      }
   }
   est /= (double)P;
   return est;
}


vector<double> Rho_Mixed::GetSpringLength(vector<Particle> part, double eps) {
   cout << "Should not be calling this function!" << endl;
   exit(-1);
   vector<double> r0;
   for(int j=0; j<(int)part.size(); j++) {
      r0.push_back(sqrt(tanh(eps*omega[j])/part[j].mass/omega[j]));
   }
   return r0;
}


vector<double> Rho_Mixed::GetLevyMean(Particle partI, Particle partF, double delta, double eps, int partNum) {
   if(partI.pos.size()>1) {
      cout << "DES Error: Must be using normal modes with mixed propagator" << endl;
      LINE
      exit(-1);
   }
   vector<double> mean;
//Free Particle part
   if(partNum < cutoffIndex) {
      mean.push_back( (partI.pos[0]*delta + partF.pos[0]) / (delta+1.0) );
   }
   else {
//Harmonic part
      double tDelI = tanh(eps*omega[partNum]);
      double tDelF = tanh(delta*eps*omega[partNum]);
      double sDelI = sinh(eps*omega[partNum]);
      double sDelF = sinh(delta*eps*omega[partNum]);
      mean.push_back( (partI.pos[0]/sDelI + partF.pos[0]/sDelF) * tDelI*tDelF/(tDelI+tDelF) );
   }
   return mean;
}

double Rho_Mixed::GetLevySigma(double delta, double eps, int partNum) {
   double sigma;
//Free Particle part
   if(partNum < cutoffIndex) {
      sigma = sqrt(eps/(1.0+1.0/delta));
   }
   else {
//Harmonic part
      double tDelI = tanh(eps*omega[partNum]);
      double tDelF = tanh(delta*eps*omega[partNum]);
      sigma = sqrt(tDelI*tDelF/(tDelI+tDelF) / omega[partNum]);
   }
   return sigma;
}


string Rho_Mixed::GetType() {
//TODO: add cutoff here
   string name = "Mixed";
   return name;
}
