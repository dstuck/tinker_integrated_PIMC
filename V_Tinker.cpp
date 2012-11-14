/*
 * V_Tinker.cpp
 *
 *  Created on: Sep 7, 2012
 *      Author: dstuck
 */

#include "V_Tinker.h"


V_Tinker::V_Tinker(CoordUtil* coords, double eps, double beta) : coordKeeper(coords) {
	tinkInFileName = coords->tinkerName + ".xyz";
	tinkOutFileName = coords->tinkerName + ".outTinker";
	tinkPrmFileName = coords->prmName;
//	coordKeeper = coords;
//    Initialize forcefield
	vector<Particle> parts;
	Particle p0(0.0, 1);
	Propagator * rhoFree = new Rho_Free;
	for(int i=0; i<coordKeeper->numModes; i++) {
		parts.push_back(p0);
	}
	inFile.open(tinkInFileName.c_str());
	outFile.open(tinkOutFileName.c_str());
	int numPart = coordKeeper->numPart;
//	Get cartesians from normal modes
//	vector< vector<double> > cartPos = coordKeeper->normalModeToCart(parts);
	vector< vector<double> > cartPos = coordKeeper->initCart;
//	Write tinker inputfile
	inFile << numPart << endl;
	for(int i=0; i<numPart; i++) {
		inFile << i+1 << "\t" << coordKeeper->atomType[i] << "\t" << cartPos[i][0] << "\t" << cartPos[i][1] << "\t" << cartPos[i][2];
		for(int c=0; c<int(coordKeeper->connectivity[i].size()); c++) {
				inFile << "\t" << coordKeeper->connectivity[i][c];
		}
		inFile << endl;
	}
	inFile.close();
        bool init = true;
        int inLen = tinkInFileName.size();
//      cout << "inLen is: " << inLen << endl;
        int prmLen = tinkPrmFileName.size();
        double tinkEnergy = 0.0; 
        char * tinkIn = new char[120];
        strcpy(tinkIn, tinkInFileName.c_str());      //So fortran doesn't mess up the original
        string tinkPrm = tinkPrmFileName;
        dstuckenergy_(tinkIn, inLen, tinkPrm.c_str(), prmLen, tinkEnergy, init);

//    Set up tinker normal modes and frequencies!
        int N = coordKeeper->numModes;
        double tinkModes[N*coordKeeper->numPart*3];
        double tinkFreqs[N];
        double redMass[N];
        dstuckvibrate_(tinkIn, inLen, N, coordKeeper->numPart, redMass, tinkModes, tinkFreqs);
      delete[] tinkIn;
      coordKeeper->reducedMass.clear();
      for(int i=0; i<N; i++) {
         coordKeeper->reducedMass.push_back(redMass[i]*1822.8886);
      }
      coordKeeper->omega.clear();
      for(int i=0; i<N; i++) {
         coordKeeper->omega.push_back(tinkFreqs[i]*0.00000455633);
//         coords.omega[i] = tinkFreqs[i]*0.00000455633;
//         cout << tinkFreqs[i] << endl;
      }
      coordKeeper->normModes.clear();
      coordKeeper->normModes.resize(N);
      for(int i=0; i<N; i++) {
         coordKeeper->normModes[i].resize(coordKeeper->numPart);
         for(int j=0; j<coordKeeper->numPart; j++) {
            coordKeeper->normModes[i][j].resize(3,0.0);
            for(int k=0; k<3; k++) {
               coordKeeper->normModes[i][j][k] = tinkModes[k+3*j+3*coordKeeper->numPart*i];
            }
         }
      }
//      for(int i=0; i<N; i++) {
//         cout << tinkFreqs[i]*0.00000455633 << endl;
//         for(int j=0; j<coordKeeper->numPart; j++) {
//            for(int k=0; k<3; k++) {
//               cout << coordKeeper->normModes[i][j][k] << "\t";
//            }
//         cout << endl;
//         }
//      }
      double clHarmFull = 0.0;
      for(int i=0; i<N; i++) {
         clHarmFull += coordKeeper->omega[i]/2.0;
      }
      cout << "Harmonic Energy is: " << clHarmFull << endl;

//    Get vEquib
	double hartreeToKcal = 627.509469;						//From http://en.wikipedia.org/wiki/Hartree 8/17/2012
	vEquib = 0.0;
	vEquib = GetV(parts, rhoFree)*hartreeToKcal;
     
/* 
//    Trying to use the 5 point second derivative approximation from http://en.wikipedia.org/wiki/Five-point_stencil
   vector< vector<double> > derivPoints;
   vector<double> tempDeriv;
   vector<double> numStepSize;
   vector<double> numDeriv;
   double clHarmZPE = 0.0;
   double clHarmFull = 0.0;
   for(int i=0; i<coords.numModes; i++) {
      numStepSize.push_back(tanh(eps*w[i])/2.0/mass[i]/w[i]);
   }
   for(int i=0; i<coords.numModes; i++) {
      for(int p=0; p<5; p++) {
         parts[i].pos[0] = (double)(p-2)*numStepSize[i];
         double vtemp = GetV(parts, rhoFree);
         tempDeriv.push_back(vtemp);
      }
      parts[i].pos[0]=0.0;
      derivPoints.push_back(tempDeriv);
      tempDeriv.clear();
   }
// TODO: Make this one loop
   for(int i=0; i<coords.numModes; i++) {
      numDeriv.push_back(0.0);
      numDeriv[i] += -derivPoints[i][0];
      numDeriv[i] += 16.0*derivPoints[i][1];
      numDeriv[i] += -30.0*derivPoints[i][2];
      numDeriv[i] += 16.0*derivPoints[i][3];
      numDeriv[i] += -derivPoints[i][4];
      numDeriv[i] /= (12.0*numStepSize[i]*numStepSize[i]);
//      for(int p=0; p<5; p++) {
//         cout << derivPoints[i][p] << "\t";
//      }
//      cout << "\tStep: " << numStepSize[i] <<  "\tDerive is: " << numDeriv[i];
//      cout << endl;
      numDeriv[i] = sqrt(numDeriv[i]/mass[i]);
      clHarmZPE += numDeriv[i]/2.0;
      clHarmFull += numDeriv[i]/2.0/tanh(beta*numDeriv[i]/2.0);
   }
//	TODO: Remove this
//	cout << "Tinker Harmonic ZPE:\t" << clHarmZPE << endl;
	cout << "Tinker Harmonic Energy:\t" << clHarmFull << endl;
*/
}

V_Tinker::~V_Tinker() {
}

double V_Tinker::GetV(vector<Particle> part, Propagator * rho){

	double hartreeToKcal = 627.509469;						//From http://en.wikipedia.org/wiki/Hartree 8/17/2012
	inFile.open(tinkInFileName.c_str());
	outFile.open(tinkOutFileName.c_str());
//      cout << "inFileName: " << tinkInFileName.c_str() << endl;

	double V = 0;
	int N = coordKeeper->numPart;
	//	Get cartesians from normal modes
//	vector< vector<double> > cartPos = coordKeeper->normalModeToCart(part);
	vector< vector<double> > cartPos = coordKeeper->initCart;

//	Write tinker inputfile
	inFile << N << endl;
	for(int i=0; i<N; i++) {
		inFile << i+1 << "\t" << coordKeeper->atomType[i] << "\t" << cartPos[i][0] << "\t" << cartPos[i][1] << "\t" << cartPos[i][2];
		for(int c=0; c<int(coordKeeper->connectivity[i].size()); c++) {
				inFile << "\t" << coordKeeper->connectivity[i][c];
		}
		inFile << endl;
	}
	inFile.close();

// Call Tinker
      bool init = false;
      int inLen = tinkInFileName.size();
//   cout << "inLen is: " << inLen << endl;
      int prmLen = tinkPrmFileName.size();
      double tinkEnergy = 0.0; 
      char * tinkIn = new char[120];
      strcpy(tinkIn, tinkInFileName.c_str());      //So fortran doesn't mess up the original
      string tinkPrm = tinkPrmFileName;
//   cout <<  "Prm file is " << tinkPrm;
      dstuckenergy_(tinkIn, inLen, tinkPrm.c_str(), prmLen, tinkEnergy, init);
      V = tinkEnergy;
//   cout << "tinkEnergy is " << tinkEnergy << endl;
	V -= vEquib;
//	cout << "Unmodified V =\t" << V << endl;
      V /= hartreeToKcal;
	V += rho->ModifyPotential(part);
//	cout << "Modified V =\t" << V << endl;
	outFile.close();
      delete[] tinkIn;
	return V;
}

string V_Tinker::GetType() {
	string name = "Tinker";
	return name;
}
