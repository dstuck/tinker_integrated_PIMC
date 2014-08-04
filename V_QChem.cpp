/*
 * V_QChem.cpp
 *
 *  Created on: June 23, 2014
 *      Author: dstuck
 */

#include "V_QChem.h"


V_QChem::V_QChem(CoordUtil* coords, CoordUtil* tinkCoords,  double eps, double beta) : coordKeeper(coords), tinkerCoords(tinkCoords) {
// Set rems

// Set equib pos
   vector<Particle> parts;
   Particle p0(0.0, 1);
   for(int i=0; i<coordKeeper->numModes; i++) {
      parts.push_back(p0);
   }

   vector< vector<double> > cartPos = coordKeeper->initCart;
   int nPart = coordKeeper->numPart;
   /*
      double xyzCoord[nPart*3];
      for(int i=0; i<nPart; i++) {
      for(int k=0; k<3; k++) {
      xyzCoord[k+i*3] = cartPos[i][k];
      }
      }
    */

// Run equib frequency
// qchem_freq();
// DES: For now just used read in frequency
   tempNum=0;

// Save vEquib and wQChem and set normal modes

   double clHarmFull = 0.0;
   for(int i=0; i<coordKeeper->omega.size(); i++) {
      clHarmFull += coordKeeper->omega[i]/2.0;
   }
   cout << "QChem harmonic Energy is: " << clHarmFull << endl;

// Set up harmonic V for later
   qharmFile.open("qPotential.txt");
   harmV.omega = coordKeeper->omega;

}

V_QChem::~V_QChem() {
   delete coordKeeper;
   qharmFile.close();
}

double V_QChem::GetV(vector<Particle> part, Propagator * rho){

   double V = 0.0;
   V += GetV(part);
   V += rho->ModifyPotential(part);
   return V;
}

double V_QChem::GetV(vector<Particle> part){

   double hartreeToKcal = 627.509469;						//From http://en.wikipedia.org/wiki/Hartree 8/17/2012
   double V = 0.0;
//	Get cartesians from normal modes
   vector< vector<double> > cartPos = tinkerCoords->normalModeToCart(part);

// DES: running ab initio thermodynamic integration correction to molecular mechanics
//   int pickedSlice = 0;        //TODO: Make this random
//            Propagator * freerho = new Rho_Free;
//            double pickedV = V->GetV(part[pickedSlice], freerho);
//            delete freerho;
// DES Temp: Just print a .in file with geometry
   ostringstream convert;
   convert << tempNum;
   string qchemName = convert.str() + "_pimc.in";
//string qchemName = std::to_string(tempNum) + "_pimc.in";
   qchemFile.open(qchemName.c_str());
   int N = coordKeeper->numPart;
//	Get cartesians from normal modes
   //vector< vector<double> > cartPos = coordKeeper->normalModeToCart(part[pickedSlice]);
/*
//DES Temp Temp:
   cout << "QC cartPos in V_QChem" << endl;
   for(int j=0; j<cartPos.size(); j++) {
      for(int k=0; k<cartPos[j].size(); k++) {
         cout << (cartPos[j][k] - coordKeeper->initCart[j][k])/0.52918 << endl;
      }
   }

   cout <<"DES Temp: QC cton" << endl;
   vector<double> normMode = coordKeeper->cartToNormalMode(cartPos);
   cout << "QC Normal modes" << endl;
   for(int i=0; i<part.size(); i++) {
      cout << part[i].pos[0] << "\t" << normMode[i] << endl;
   }
   exit(-1);
   cout << "QC Cartesians1" << endl;
   for(int i=0; i<N; i++) {
   cout << cartPos[i][0] << "\t" << cartPos[i][1] << "\t" << cartPos[i][2] << endl;
   }

   for(int i=0; i<part.size(); i++) {
      part[i].pos[0] = normMode[i];
   }
   cartPos = coordKeeper->normalModeToCart(part);
   cout << "QC Cartesians2" << endl;
   for(int i=0; i<N; i++) {
   cout << cartPos[i][0] << "\t" << cartPos[i][1] << "\t" << cartPos[i][2] << endl;
   }
   exit(-1);
*/
//	Write qchem inputfile
   qchemFile << "$comments\nDES: deltaPIMC job\n$end\n" << endl;
   qchemFile << "$molecule\n 0 1" << endl;
   for(int i=0; i<N; i++) {
      qchemFile << coordKeeper->atomType[i] << "\t" << cartPos[i][0] << "\t" << cartPos[i][1] << "\t" << cartPos[i][2];
      qchemFile << endl;
   }
   qchemFile << "$end\n" << endl;
// B3LYP
   qchemFile << "$rem\n         JOBTYPE           sp\n         EXCHANGE          B3LYP\n         BASIS             cc-pVTZ\n         SCF_GUESS         SAD\n         SCF_ALGORITHM     DIIS\n         MAX_SCF_CYCLES    200\n         SYM_IGNORE        TRUE\n         SYMMETRY          FALSE\n         UNRESTRICTED      TRUE\n         SCF_CONVERGENCE   6\n         THRESH         9\n         MEM_STATIC       1000 \n         MEM_TOTAL       4000\n$end" << endl;
//            qchemFile << "$rem\n         JOBTYPE           sp\n         EXCHANGE          HF\n         CORRELATION          MP2\n         BASIS             cc-pVTZ\n   MP2_RESTART_NO_SCF   TRUE\n      purecart       2222\n         SCF_GUESS         READ\n         SCF_ALGORITHM     DIIS\n         THRESH_DIIS_SWITCH  4\n         MAX_SCF_CYCLES    200\n         SYM_IGNORE        TRUE\n         SYMMETRY          FALSE\n         UNRESTRICTED      TRUE\n         SCF_CONVERGENCE   8\n         DO_O2          0\n         THRESH         14\n         MEM_STATIC       1000 \n         MEM_TOTAL       4000\n$end" << endl;
   qchemFile.close();
 
   vector<double> modes = coordKeeper->cartToNormalMode(cartPos);
//   cout << "QC Normal modes" << endl;
//   for(int i=0; i<part.size(); i++) {
//      cout << part[i].pos[0] << "\t" << modes[i] << endl;
//   }
   vector<Particle> tempPart;
   for(int i=0; i<modes.size(); i++) {
      tempPart.push_back(Particle(modes[i],1));
      tempPart[i].mass = coordKeeper->reducedMass[i];
   }
   double potHarm = harmV.GetV(tempPart);
   qharmFile << tempNum << "\t" << potHarm << endl;

   tempNum++;

//    V = qchem_energy();
//      V -= vEquib;
//	return V;
   return 0.0;
}

string V_QChem::GetType() {
   string name = "QChem";
   return name;
}

CoordUtil* V_QChem::GetCoordUtil() {
   return coordKeeper;
}
