/*
 * CoordUtil.cpp
 *
 *  Created on: Sep 16, 2012
 *      Author: dstuck
 */

#include "CoordUtil.h"
#include <stdlib.h> 
//#include "optking/molecule.h"
#include "optking/globals.h"

namespace opt {
   void set_params(void);
}

CoordUtil::CoordUtil() {
}

//CoordUtil::CoordUtil(int nMode, int nPart, bool internalCoords, vector <vector < vector <double> > > modes, vector<double> freqs, vector<double> m, vector <vector <double> > initPos, vector<string> atomicSymbols, vector<int> params, vector< vector<int> > conn, string tinkName, string prmFile, bool readG, bool readW) : numModes(nMode), numPart(nPart), internals(internalCoords), normModes(modes), omega(freqs), reducedMass(m), initCart(initPos), atomType(atomicSymbols), paramType(params), connectivity(conn), outFileName(tinkName), prmName(prmFile), readGeom(readG), readOmega(readW)  {
CoordUtil::CoordUtil(int nMode, int nPart, bool internalCoords, vector <vector < vector <double> > > modes, vector<double> freqs, vector<double> m, vector <vector <double> > initPos, vector<string> atomicSymbols, vector<int> params, vector< vector<int> > conn, string tinkName, string prmFile, bool readG, bool readW) {
   numModes=nMode;
   numPart=nPart;
   internals=internalCoords;
   normModes=modes;
   omega=freqs;
   reducedMass=m;
   initCart=initPos;
   atomType=atomicSymbols;
   paramType=params;
   connectivity=conn;
   outFileName=tinkName;
   prmName=prmFile;
   readGeom=readG;
   readOmega=readW;

   useGuess = false;
   guessCarts = initCart;
   Particle part;
   for(int i=0; i<numModes; i++) {
      part = Particle(0,1);
      guessPart.push_back(part);
   }
//	normModes = modes;				//normModes[i][j][k] is kth dimension of jth atom for the ith normal mode
//        omega = freqs;
//	initCart = initPos;				//initPos[i][j] is jth dimension of ith atom
//	connectivity = conn;			//connectivity[i][c] is the number of the cth atom connected to atom i (not square)
//	outFileName = tinkName;
//	prmName = prmFile;
//	paramType = params;
//	atomType = atomicSymbols;
//	numModes = modes.size();
//	numPart = modes[0].size();
//	dim = modes[0][0].size();
        scaleGeom = false;
        dim = 3;
	if(dim != 3) {
		cout << "Warning: Not working in xyz coordinates!" << endl;
	}
      initMass();
/*
      atomMass = arma::vec(numPart);
      for(int i=0; i<numPart; i++) {
         atomMass(i) = atomToMass[atomType[i]];
      }
*/
/*
      armaNormModes = arma::mat(numPart*dim,numModes);
      for(int j=0; j<numPart; j++) {
         for(int k=0; k<dim; k++) {
            //armaCart(3*j+k) = (cart[j][k]-initCart[j][k])/0.52918;
            for(int i=0; i<numModes; i++) {
               armaNormModes(dim*j+k,i) = normModes[i][j][k];
            }
         }
      }
*/

//   initMass();
}

void CoordUtil::initInternals() {
//Internal Coordinates
   if(internals) {
      atomMass = arma::vec(numPart);
      for(int i=0; i<numPart; i++) {
         atomMass(i) = atomToMass[atomType[i]];
      }
      armaNormModes = arma::mat(numPart*dim,numModes);
      for(int j=0; j<numPart; j++) {
         for(int k=0; k<dim; k++) {
            //armaCart(3*j+k) = (cart[j][k]-initCart[j][k])/0.52918;
            for(int i=0; i<numModes; i++) {
               armaNormModes(dim*j+k,i) = normModes[i][j][k];
            }
         }
      }
      
      opt::set_params();
      workingMol = new opt::MOLECULE(numPart);
      int cnt=0;
      for (int f=0; f<workingMol->fragments.size(); ++f) {        //TODO: remove fragment loop here...
         double *Z = workingMol->fragments[f]->g_Z_pointer();
         double **geom = workingMol->fragments[f]->g_geom_pointer();

         for (int i=0; i<workingMol->fragments[f]->g_natom(); ++i) {
            Z[i] = atomStrToInt(atomType[cnt]);
            for (int xyz=0; xyz<3; ++xyz) {
               geom[i][xyz] = initCart[i][xyz];
            }
            ++cnt;
         }
         //workingMol->fragments[f]->set_default_masses();
//         for(int i=0; i<3; i++) {
//            cout << "DES: test mass = " << workingMol->fragments[f]->mass[i] << endl;
//         }
      }
      opt::outfile=stdout;
      //workingMol->update_connectivity_by_distances();
// DES: set connectivity ourselves!
      for(int i=0; i<numPart; i++) {
         for(int j=0; j<numPart; j++) {
            workingMol->fragments[0]->connectivity[i][j] = false;
         }
      }
      for(int i=0; i<connectivity.size(); i++) {
         for(int j=1; j<connectivity[i].size(); j++) {
            workingMol->fragments[0]->connectivity[i][connectivity[i][j]-1] = true;
         }
      }
//      for(int i=0; i<numPart; i++) {
//         for(int j=0; j<numPart; j++) {
//            cout << workingMol->fragments[0]->connectivity[i][j] << "\t";
//         }
//         cout << endl;
//      }
      workingMol->fragmentize();
// DES TODO: get rid of this if not using inertia_tensor
      for (int f=0; f<workingMol->fragments.size(); ++f) {
         workingMol->fragments[f]->set_default_masses();
         workingMol->fragments[f]->print_geom(stdout);
      }
// DES: fragment_mode=simple:
      //workingMol->add_intrafragment_simples_by_connectivity();
// DES: fragment_mode=multi:
      //workingMol->update_connectivity_by_distances();
      workingMol->add_intrafragment_simples_by_connectivity();
      //workingMol->add_interfragment();
      workingMol->print_connectivity(stdout);
      //workingMol->print_intcos(stdout);    //Fuller report
      //workingMol->print_intco_dat(stdout);
      if(workingMol->g_nfragment() > 1) {
         cout << "---Interfragment Internal---" << endl;
         cout << "Adding " << (workingMol->g_nfragment()-1)*3 << " fragment translations" << endl;
         cout << "Adding " << workingMol->g_nfragment()*3 << " fragment rotations" << endl;
      }
      //double ** bArray = workingMol->fragments[0]->compute_B();
      double ** bArray = workingMol->compute_B();
      int nPrimInt = workingMol->g_nintco();
      equibBMat = arma::mat(nPrimInt,numPart*3);
      for(int i=0; i< nPrimInt; i++) {
         for(int j=0; j<numPart*3; j++) {
            equibBMat(i,j) = bArray[i][j];
         }
      }
      opt::free_matrix(bArray);
//DES: Add fragment translation and rotation
// First need to get fragment info
// workingMol->g_nfragment()
// workingMol->fragments[f]->globalIndex[i] gets index of ith atom in fth frag

      arma::vec bCol;
      //equibBMat.print("preB");
      double alpha, sinT, cosT, massAB;
      //fragExtraB = arma::mat(0,nPart*3);
//Add cartessians 3*nfrag displacements
/* Non COM preserving
      for(int f=0; f<workingMol->g_nfragment(); f++) {
         alpha = 1/sqrt(double(workingMol->g_nfragment()));
         for(int j=0; j<3; j++) {
            bCol = arma::vec(3*numPart,arma::fill::zeros);
            for(int i=0; i<workingMol->fragments[f]->g_natom(); i++) {
               bCol(3*workingMol->fragments[f]->globalIndex[i] + j) += alpha;
            }
            fragExtraB = arma::join_cols(fragExtraB,bCol.t());
         }
      }
*/
// Add COM preserving translations
      double totalMass = 0.0;
      for(int i=0; i<numPart; i++) {
         totalMass += atomMass(i);
      }
//DES Temp prim
      for(int f=0; f<workingMol->g_nfragment()-1; f++) {
      //for(int f=0; f<workingMol->g_nfragment(); f++) {}
         massAB=0.0;
         for(int i=0; i<workingMol->fragments[f]->g_natom(); i++) {
            massAB += atomMass(workingMol->fragments[f]->globalIndex[i]);
         }
         alpha = totalMass - massAB;
         massAB /= alpha;
         alpha=1/sqrt(workingMol->fragments[f]->g_natom()+massAB*massAB*(numPart-workingMol->fragments[f]->g_natom()));
         for(int j=0; j<3; j++) {
            bCol = arma::vec(3*numPart,arma::fill::zeros);
            for(int i=0; i<numPart; i++) {
               bCol(3*i + j) = -alpha*massAB*atomMass(i);
               //bCol(3*i + j) = -alpha*massAB;
            }
            for(int i=0; i<workingMol->fragments[f]->g_natom(); i++) {
               bCol(3*workingMol->fragments[f]->globalIndex[i] + j) = alpha * atomMass(workingMol->fragments[f]->globalIndex[i]);
               //bCol(3*workingMol->fragments[f]->globalIndex[i] + j) = alpha;
            }
            fragExtraB = arma::join_cols(fragExtraB,bCol.t());
         }
      }
//Add 3*nfrag rotations
// alpha is set to norm of 
      if(workingMol->g_nfragment()>1) {
         double * fragCOM;
         for(int f=0; f<workingMol->g_nfragment(); f++) {
            fragCOM = workingMol->fragments[f]->com();
            //for(int j=0;j<3;j++) fragCOM[j] = 0.0;
            //for(int j=0;j<3;j++) cout << "DES: COM = " << fragCOM[j] << endl;
            for(int j=0; j<3; j++) {
               bCol = arma::vec(3*numPart,arma::fill::zeros);
/* Old angle based (rather than tangent)
               alpha=0;
               for(int k=0; k<2; k++) {
                  int tempk = (k+j)%3;
                  for(int i=0; i<workingMol->fragments[f]->g_natom(); i++) {
                     alpha += (initCart[workingMol->fragments[f]->globalIndex[i]][tempk]-fragCOM[tempk])*(initCart[workingMol->fragments[f]->globalIndex[i]][tempk]-fragCOM[tempk]);
                  }
               }
               cosT = 1.0/alpha;
               sinT = sqrt(2.0*alpha-1.0)/alpha;
               for(int i=0; i<workingMol->fragments[f]->g_natom(); i++) {
                  bCol(3*workingMol->fragments[f]->globalIndex[i] + (j+1)%3) = cosT*(initCart[workingMol->fragments[f]->globalIndex[i]][(j+1)%3]-fragCOM[(j+1)%3]) - sinT*(initCart[workingMol->fragments[f]->globalIndex[i]][(j+2)%3]-fragCOM[(j+2)%3]) * atomMass(workingMol->fragments[f]->globalIndex[i]);
                  bCol(3*workingMol->fragments[f]->globalIndex[i] + (j+2)%3) = sinT*(initCart[workingMol->fragments[f]->globalIndex[i]][(j+1)%3]-fragCOM[(j+1)%3]) + cosT*(initCart[workingMol->fragments[f]->globalIndex[i]][(j+2)%3]-fragCOM[(j+2)%3]) * atomMass(workingMol->fragments[f]->globalIndex[i]);
                  //bCol(3*workingMol->fragments[f]->globalIndex[i] + (j+1)%3) = cosT*(initCart[workingMol->fragments[f]->globalIndex[i]][(j+1)%3]-fragCOM[(j+1)%3]) - sinT*(initCart[workingMol->fragments[f]->globalIndex[i]][(j+2)%3]-fragCOM[(j+2)%3]) / atomMass(workingMol->fragments[f]->globalIndex[i]);
                  //bCol(3*workingMol->fragments[f]->globalIndex[i] + (j+2)%3) = sinT*(initCart[workingMol->fragments[f]->globalIndex[i]][(j+1)%3]-fragCOM[(j+1)%3]) + cosT*(initCart[workingMol->fragments[f]->globalIndex[i]][(j+2)%3]-fragCOM[(j+2)%3]) / atomMass(workingMol->fragments[f]->globalIndex[i]);
                  //bCol(3*workingMol->fragments[f]->globalIndex[i] + (j+1)%3) = cosT*(initCart[workingMol->fragments[f]->globalIndex[i]][(j+1)%3]-fragCOM[(j+1)%3]) - sinT*(initCart[workingMol->fragments[f]->globalIndex[i]][(j+2)%3]-fragCOM[(j+2)%3]);
                  //bCol(3*workingMol->fragments[f]->globalIndex[i] + (j+2)%3) = sinT*(initCart[workingMol->fragments[f]->globalIndex[i]][(j+1)%3]-fragCOM[(j+1)%3]) + cosT*(initCart[workingMol->fragments[f]->globalIndex[i]][(j+2)%3]-fragCOM[(j+2)%3]);
               }
               //fragExtraB = arma::join_cols(fragExtraB,bCol.t());
*/
               for(int i=0; i<workingMol->fragments[f]->g_natom(); i++) {
                  bCol(3*workingMol->fragments[f]->globalIndex[i] + (j+1)%3) = (initCart[workingMol->fragments[f]->globalIndex[i]][(j+2)%3]-fragCOM[(j+2)%3]) * atomMass(workingMol->fragments[f]->globalIndex[i]);
                  bCol(3*workingMol->fragments[f]->globalIndex[i] + (j+2)%3) = -(initCart[workingMol->fragments[f]->globalIndex[i]][(j+1)%3]-fragCOM[(j+1)%3]) * atomMass(workingMol->fragments[f]->globalIndex[i]);
               }
   //DES Temp prim:
               equibBMat = arma::join_cols(equibBMat,bCol.t());
               nPrimInt++;
            }
            opt::free_array(fragCOM);
         }
      }
      

//DES Temp prim:
      equibBMat = arma::join_cols(equibBMat,fragExtraB);
      nPrimInt += 3*(workingMol->g_nfragment()-1);
      //nPrimInt += 3*workingMol->g_nfragment();   
      //equibBMat.print("postB");

//DES: Project out global rotations/translations because of fragment motions
//TODO: This will crash on linear molecules!
      //arma::mat projector(3*numPart,3*numPart,arma::fill::zeros);
      arma::mat projector,projVecs,tempMat;
      arma::vec tempVec;
      arma::vec totalCOM(3,arma::fill::zeros);
      for(int i=0; i<numPart; i++) {
         for(int j=0; j<3; j++) {
            totalCOM(j) += atomMass(i) * initCart[i][j];
         }
      }
   
// global translations
      for(int n=0; n<3; n++) {
         tempVec = arma::vec(3*numPart,arma::fill::zeros);
         for(int i=0; i<numPart; i++) {
            tempVec(3*i + n) = atomMass(i);
         }
         tempVec = arma::normalise(tempVec);
         //tempVec.print("global translation");
         //projector += tempVec*tempVec.t();
         projVecs = join_rows(projVecs,tempVec);
      }
// global rotations
      if(numPart<3) {
         cout << "Error! Not set up to use internal coordinates with linear molecules:\n\t\tJust use cartesians anyway!" << endl;
         exit(-1);
      }
      for(int n=0; n<3; n++) {
         tempVec = arma::vec(3*numPart,arma::fill::zeros);
         for(int i=0; i<numPart; i++) {
            tempVec(3*i + (n+1)%3) = atomMass(i) * (initCart[i][(n+2)%3] - totalCOM((n+2)%3));
            tempVec(3*i + (n+2)%3) = -atomMass(i) * (initCart[i][(n+1)%3] - totalCOM((n+1)%3));
         }
         tempVec = arma::normalise(tempVec);
         //tempVec.print("global rotation");
         //projector += tempVec*tempVec.t();
         projVecs = join_rows(projVecs,tempVec);
      }
      //projVecs.print("projVecs");
      arma::qr(projector,tempMat,projVecs);
      //projector.print("Q");
      projector = projector.cols(0,5)*projector.cols(0,5).t();
      //projector.print("projector");
      //(projector*projector-projector).print("projector^2-projector");
      equibBMat -= equibBMat * projector; 
      //equibBMat.print("projected B");

      equibBMat = arma::normalise(equibBMat,2,1);
      //cout << "bMat row norm " <<  arma::dot(equibBMat.t(),equibBMat.t())/equibBMat.n_rows << endl;
      //equibBMat.print("BMat in primitives");
      //cout << "bMat fro norm " <<  norm(equibBMat,"fro") << endl;
   // Set up delocalized internals
      sqrtInvMass = arma::mat(numPart*3,numPart*3,arma::fill::zeros);
      for(int n=0; n<numPart; n++) {
         for(int k=0; k<3; k++) {
            sqrtInvMass(3*n+k,3*n+k) = 1.0/sqrt(atomMass(n))/1822.888;
            //sqrtInvMass(3*n+k,3*n+k) = 1.0/sqrt(atomToMass[atomType[n]])/1822.888;
         }
      }
      arma::mat tempV;
      arma::vec tempS;
      arma::svd(delocIntCos,tempS,tempV,equibBMat);
      //tempS.print("Singular vals");
      //delocIntCos.print("eVecs");
      delocIntCos.resize(nPrimInt,numModes);
//DES Temp prim:
      //delocIntCos = arma::eye(numModes,numModes);
      //armaNormModes.print("normal modes in cartesians");
      //(equibBMat*armaNormModes).print("normal modes in primitive internals");
      equibBMat = delocIntCos.t()*equibBMat;
      //(equibBMat*armaNormModes).print("normal modes in delocalized internals");
      //equibBInv = arma::normalise(sqrtInvMass*arma::pinv(equibBMat*sqrtInvMass));
      
      //equibBInv = sqrtInvMass*arma::pinv(equibBMat*sqrtInvMass);


      //equibBMat.print("BMat in deloc");
      //equibBInv.print("bInv in deloc");
      //(equibBMat*armaNormModes).print("normal modes in deloc internals");
      //(equibBMat*equibBInv).print("B*B-1");
      //(equibBInv*equibBMat).print("B^-1*B");
      //cout << "bMat row norm " <<  arma::dot(equibBMat.t(),equibBMat.t())/equibBMat.n_rows << endl;
      //cout << "bMat col norm " <<  norm(equibBMat,2) << endl;
      //cout << "bInv col norm " <<  arma::dot(equibBInv,equibBInv)/equibBInv.n_cols << endl;

      //delocIntCos.print("deloc coords");
      //equibBMat.print("bMat");
      //(armaNormModes).print("internal normalModes");
      //(equibBMat*armaNormModes).print("delocal internalNormalModes");
      //(delocIntCos*equibBMat*armaNormModes).print("primitive internalNormalModes");

      //exit(-1);

//    DES: Testing internals
/*
      vector< vector<double> > carts(initCart);
      //arma::vec step(nPrimInt);
      vector<double> step(numModes,0.0);
      for(int n=0; n<21; n++) {
      //for(int n=0; n<10; n++) {}
         //step[6] = 0.04*double(n+1)/0.52918;
         step[61] = 0.04*double(-10+n)/0.52918;
         //cout << "Step is " << step[0] << endl;
         carts = normalModeToCart(step);
         cout << numPart << endl;
         cout << endl;
         for(int k=0; k<numPart; k++) {
            cout << atomType[k]<<"  " << carts[k][0] << "\t" << carts[k][1] << "\t" << carts[k][2] << "\t" << endl;
         }
      }
      exit(-1);
*/
   }
   return;
}

CoordUtil::~CoordUtil() {
   atomToMass.clear();

	// TODO Auto-generated destructor stub
}

vector< vector<double> > CoordUtil::normalModeToCart(vector<double> modes) {
   vector<Particle> partVec;
   Particle part;
   for(int i=0; i< modes.size(); i++) {
      part = Particle(modes[i],1);
      partVec.push_back(part);
   }
   return normalModeToCart(partVec);
}

vector< vector<double> > CoordUtil::normalModeToCart(vector<Particle> part) {
// Really converts any part to cartesians (normal or already cartesian)
	vector< vector<double> > carts;
	//vector< vector<double> > carts(initCart);

//	For Normal modes into Cartesians in Angstroms
        if(part[0].pos.size()==1 && !internals) {
            arma::vec armaCart(numPart*dim);
            arma::vec armaVecs(numModes);
            for(int i=0; i<numModes; i++) {
               armaVecs(i) = part[i].pos[0];
            }
            if(armaNormModes.is_empty()) {
               armaNormModes = arma::mat(numPart*dim,numModes);
               for(int j=0; j<numPart; j++) {
                  for(int k=0; k<dim; k++) {
                     for(int i=0; i<numModes; i++) {
                        armaNormModes(3*j+k,i) = normModes[i][j][k];
                     }
                  }
               }
            }
            if(armaInitCart.is_empty()) {
               armaInitCart = arma::vec(numPart*dim);
               for(int j=0; j<numPart; j++) {
                  for(int k=0; k<dim; k++) {
                     armaInitCart(3*j+k) = initCart[j][k];
                  }
               }
            }
            armaCart = armaNormModes*armaVecs + armaInitCart;
            vector<double> tempVec;
            for(int j=0; j<numPart; j++) {
//               carts.push_back(arma::conv_to< vector<double> >::from(armaCart(armaInitCart(arma::span(3*j,3*j+2)))));
               for(int k=0; k<3; k++) {
                  tempVec.push_back(armaCart(3*j+k));
               }
               carts.push_back(tempVec);
               tempVec.clear();
            }
            return carts;
            
/* Old code spends all its time accessing std::vector elements
            carts = initCart;
            for(int j=0; j<numPart; j++) {
               for(int k=0; k<dim; k++) {
                  for(int i=0; i<numModes; i++) {
                     carts[j][k] += normModes[i][j][k]*part[i].pos[0] * 0.52918;				// Converts from Bohr radii to A, value from wikipedia 9/18/12
                  }
               }
            }
*/
         }
         else if(internals) {
//            for(int i=0; i<numModes; i++) {
//               cout << "DES: guessPart = " << guessPart[i].pos[0] << endl;
//            }
            if(useGuess) {
               carts = guessCarts;
            }
            else {
               carts = initCart;
            }
            //double maxStepSize = 0.005;
            double maxStepSize = 0.02;
            //double maxStepSize = 0.04;
            //arma::mat deltaCart;
            double modeNorm = 0.0;
            //double tempNorm = 0.0;

            if(useGuess) {
               for(int i=0; i<numModes; i++) {
                  modeNorm = max((part[i].pos[0]-guessPart[i].pos[0])*(part[i].pos[0]-guessPart[i].pos[0]),modeNorm);
               }
            }
            else {
               for(int i=0; i<numModes; i++) {
                  modeNorm = max(part[i].pos[0]*part[i].pos[0],modeNorm);
               }
            }

            modeNorm = sqrt(modeNorm) * 0.52918;      //In Angstrom
// DES Temp: stepCheck
//            double stepCheck  = 0.0;
//            cout << "DES: Step size is " << modeNorm << " angstrom" << endl;
            //cout << "DES: Max modeNorm = " << modeNorm << endl;
            int numSteps = int(modeNorm/maxStepSize)+1;
            //int numSteps = max(1, int(modeNorm/maxStepSize));
// DES Temp: For carts
            //numSteps = 1;
            //cout << "numSteps = " << numSteps <<  endl;
            //cout << "maxNorm  = " << modeNorm <<  endl;
// DES Temp prim:
///*            
            if(useGuess) {
               cout << "Shouldn't be in useGuess in CoordUtil.cpp" << endl;
               for(int j=0; j<numPart; j++) {
                  for(int k=0; k<dim; k++) {
                     for(int i=0; i<numModes; i++) {
                        carts[j][k] += normModes[i][j][k]*(part[i].pos[0]-guessPart[i].pos[0]) / double(numSteps) * 0.52918;				// Converts from Bohr to A radii, value from wikipedia 9/18/12
                        //carts[j][k] += normModes[i][j][k]*part[i].pos[0] / double(numSteps) * 0.52918;				// Converts from Bohr to A radii, value from wikipedia 9/18/12
                     }
                  }
               }
            }
            else {
               for(int j=0; j<numPart; j++) {
                  for(int k=0; k<dim; k++) {
                     for(int i=0; i<numModes; i++) {
                        carts[j][k] += normModes[i][j][k]*part[i].pos[0] / double(numSteps) * 0.52918;				// Converts from Bohr to A radii, value from wikipedia 9/18/12
                        //carts[j][k] += normModes[i][j][k]*part[i].pos[0] / double(numSteps) * 0.52918;				// Converts from Bohr to A radii, value from wikipedia 9/18/12
                     }
                  }
               }
            }
// DES Temp: stepCheck
//            stepCheck += modeNorm/double(numSteps);
//            cout << "DES: step 1 =" << modeNorm/double(numSteps) << endl;
//*/

            arma::vec step(numModes);
            if(numSteps>1) {
//               for (int j=0; j<numPart; ++j) {
//                  for (int k=0; k<3; ++k) {
//                     cout << carts[j][k] << "\t";
//                  }
//                  cout << endl;
//               }
               double geomArray[3*numPart];
               for (int j=0; j<numPart; ++j) {
                  for (int k=0; k<3; ++k) {
                     geomArray[3*j+k] = carts[j][k];
                  }
               }
               workingMol->set_geom_array(geomArray);

               for(int i=0; i<numModes; i++) {
                  if(useGuess) {
                     step(i) = (part[i].pos[0]-guessPart[i].pos[0]);
                  }
                  else {
                     step(i) = part[i].pos[0];
                  }
               }
//DES: Project out global rotations/translations because of fragment motions
//TODO: This will crash on linear fragments!
               arma::vec totalCOM(3,arma::fill::zeros);
               for(int i=0; i<numPart; i++) {
                  for(int j=0; j<3; j++) {
                     totalCOM(j) += atomMass(i) * carts[i][j];
                  }
               }
/*
               arma::mat projector,projVecs,tempMat;
               arma::vec tempVec;
// global translations
               for(int n=0; n<3; n++) {
                  tempVec = arma::vec(3*numPart,arma::fill::zeros);
                  for(int i=0; i<numPart; i++) {
                     tempVec(3*i + n) = atomMass(i);
                  }
                  tempVec = arma::normalise(tempVec);
                  //tempVec.print("global translation");
                  //projector += tempVec*tempVec.t();
                  projVecs = join_rows(projVecs,tempVec);
               }
// global rotations
               for(int n=0; n<3; n++) {
                  tempVec = arma::vec(3*numPart,arma::fill::zeros);
                  for(int i=0; i<numPart; i++) {
                     tempVec(3*i + (n+1)%3) = atomMass(i) * (carts[i][(n+2)%3] - totalCOM((n+2)%3));
                     tempVec(3*i + (n+2)%3) = -atomMass(i) * (carts[i][(n+1)%3] - totalCOM((n+1)%3));
                  }
                  tempVec = arma::normalise(tempVec);
                  //tempVec.print("global rotation");
                  //projector += tempVec*tempVec.t();
                  projVecs = join_rows(projVecs,tempVec);
               }
               arma::qr(projector,tempMat,projVecs);
               projector = projector.cols(0,5)*projector.cols(0,5).t();
*/
//DES Temp prim:
               //for(int n=0; n<numSteps; n++) 
               for(int n=1; n<numSteps; n++) {
                  arma::mat projector,projVecs,tempMat;
                  arma::vec tempVec;
   // global translations
                  for(int nn=0; nn<3; nn++) {
                     tempVec = arma::vec(3*numPart,arma::fill::zeros);
                     for(int i=0; i<numPart; i++) {
                        tempVec(3*i + nn) = atomMass(i);
                     }
                     //tempVec = arma::normalise(tempVec);
                     //tempVec.print("global translation");
                     //projector += tempVec*tempVec.t();
                     projVecs = join_rows(projVecs,tempVec);
                  }
   // global rotations
                  for(int nn=0; nn<3; nn++) {
                     tempVec = arma::vec(3*numPart,arma::fill::zeros);
                     for(int i=0; i<numPart; i++) {
                        tempVec(3*i + (nn+1)%3) = atomMass(i) * (carts[i][(nn+2)%3] - totalCOM((nn+2)%3));
                        tempVec(3*i + (nn+2)%3) = -atomMass(i) * (carts[i][(nn+1)%3] - totalCOM((nn+1)%3));
                     }
                     //tempVec = arma::normalise(tempVec);
                     //tempVec.print("global rotation");
                     //projector += tempVec*tempVec.t();
                     projVecs = join_rows(projVecs,tempVec);
                  }
                  arma::qr(projector,tempMat,projVecs);
                  projector = projector.cols(0,5)*projector.cols(0,5).t();

                  int nPrimInt = workingMol->g_nintco();
                  arma::mat bTemp(nPrimInt,numPart*3);
                  double ** bArray = workingMol->compute_B();
                  for(int i=0; i< nPrimInt; i++) {
                     for(int j=0; j<numPart*3; j++) {
                        bTemp(i,j) = bArray[i][j];
                     }
                  }
                  opt::free_matrix(bArray);
// Add rotations
                  if(workingMol->g_nfragment()>1) {
                     arma::vec bCol;
                     //double alpha,sinT,cosT;
                     double * fragCOM;
                     for(int f=0; f<workingMol->g_nfragment(); f++) {
                        fragCOM = workingMol->fragments[f]->com();
                        //for(int j=0;j<3;j++) cout << "DES: COM = " << fragCOM[j] << endl;;
                        for(int j=0; j<3; j++) {
                           //alpha=0;
                           bCol = arma::vec(3*numPart,arma::fill::zeros);
                           for(int i=0; i<workingMol->fragments[f]->g_natom(); i++) {
                              bCol(3*workingMol->fragments[f]->globalIndex[i] + (j+1)%3) = (carts[workingMol->fragments[f]->globalIndex[i]][(j+2)%3]-fragCOM[(j+2)%3]) * atomMass(workingMol->fragments[f]->globalIndex[i]);
                              bCol(3*workingMol->fragments[f]->globalIndex[i] + (j+2)%3) = -(carts[workingMol->fragments[f]->globalIndex[i]][(j+1)%3]-fragCOM[(j+1)%3]) * atomMass(workingMol->fragments[f]->globalIndex[i]);
                           }
   //DES Temp prim:
                           bTemp = arma::join_cols(bTemp,bCol.t());
                           nPrimInt++;
                        }
                        opt::free_array(fragCOM);
                     }
                  }
// TODO: recalculate translations here
//DES Temp prim:
                  bTemp = arma::join_cols(bTemp,fragExtraB);
                  nPrimInt += 3*(workingMol->g_nfragment()-1);
                  //nPrimInt += 3*(workingMol->g_nfragment());
                  //bTemp.print("bTemp in step");

                  bTemp -= bTemp * projector; 

                  bTemp = delocIntCos.t()*arma::normalise(bTemp,2,1);
                  //bTemp = arma::normalise(sqrtInvMass*arma::pinv(bTemp*sqrtInvMass));
                  bTemp = sqrtInvMass*arma::pinv(bTemp*sqrtInvMass);
   //Overwriting bTemp with deltaCart
//DES Temp prim:
                  bTemp = arma::normalise(bTemp*equibBMat*armaNormModes)*step / double(numSteps) * 0.52918;
                  //bTemp = arma::normalise(bTemp)*step / double(numSteps) * 0.52918;     //For primitives only
                  //bTemp = armaNormModes*step / double(numSteps) * 0.52918;
                  //bTemp = equibBInv*arma::diagmat(wMat)*bTemp*armaNormModes*step / double(numSteps) * 0.52918;
                  //bTemp = armaNormModes*arma::diagmat(wMat)*bTemp*equibBInv*step / double(numSteps) * 0.52918;
                  for(int j=0; j<numPart; j++) {
                     for(int k=0; k<dim; k++) {
                        //geom[j][k] += bTemp(3*j+k);
                        geomArray[j*3+k] += bTemp(3*j+k);
                        carts[j][k] += bTemp(3*j+k);
                     }
                  }
                  workingMol->set_geom_array(geomArray);

// DES Temp: stepCheck
//                  stepCheck += arma::norm(bTemp,"fro");
//                  cout << "DES: step " << n+1 << " =" << arma::norm(bTemp,"fro") << endl;
      
               }
// DES Temp: stepCheck
//               cout << "DES: stepCheck = " << stepCheck << endl;
            }
// DES: Save coords for reading in next time
            guessCarts = carts;
            guessPart = part;
            useGuess=false;
         }
//	For cartesian coordinates
        else if(part[0].pos.size()==3) {
           cout << "I shouldn't be here in CoordUtil" << endl;
           for(int j=0; j<numPart; j++) {
              for(int k=0; k<3; k++) {
                 carts[j][k] += part[j].pos[k] * 0.52918;				// Converts to A from Bohr radii, value from wikipedia 9/18/12
              }
           }
        }
	else {
		cout << "Error in CoordUtil: mistaken dimension" << endl;
		exit(-1);
	}
	return carts;
}

//DES TODO: normConst is just reduced mass! Don't need to do pinv!
vector<double> CoordUtil::cartToNormalMode(vector< vector<double> > cart) {
   vector<double> nmodes;
// *** Should only be called by V_QChem at the moment ***
//	For cartesians in Angstroms into normal modes
/* Old way needed normConstant and masses
   if(cart[0].size()==3) {
      double modeVal = 0.0;
      for(int i=0; i<numModes; i++) {
         for(int j=0; j<numPart; j++) {
            for(int k=0; k<dim; k++) {
	       modeVal += normModes[i][j][k] * atomToMass[atomType[j]]*normConst[i]*normConst[i]*(cart[j][k]-initCart[j][k])/0.52918;    // Converts to A from Bohr radii, value from wikipedia 9/18/12
	    }
         }
         nmodes.push_back(modeVal);
         modeVal = 0.0;
      }
   }
   else {
      cout << "Error in CoordUtil: calling normalModeToCart without dim=3" << endl;
      exit(-1);
   }
*/
//Can't conv_to cube from std::vector
   //arma::cube armaCart = arma::conv_to<arma::cube>::from(cart);
   //arma::cube armaModes = arma::conv_to<arma::cube>::from(normModes);
   arma::vec armaCart(numPart*dim);
   arma::vec armaVecs(numModes);
   if(armaNormModes.is_empty()) {
      armaNormModes = arma::mat(numPart*dim,numModes);
      for(int j=0; j<numPart; j++) {
         for(int k=0; k<dim; k++) {
            for(int i=0; i<numModes; i++) {
               armaNormModes(3*j+k,i) = normModes[i][j][k];
            }
         }
      }
   }
   if(armaInvModes.is_empty()) {
      armaInvModes = arma::pinv(armaNormModes);
   }
   for(int j=0; j<numPart; j++) {
      for(int k=0; k<dim; k++) {
         armaCart(3*j+k) = (cart[j][k]-initCart[j][k])/0.52918;
      }
   }
   armaVecs = armaInvModes*armaCart;
   nmodes = arma::conv_to< vector<double> >::from(armaVecs);

   return nmodes;
}

double CoordUtil::GetHarmV(vector< Particle > part) {
   double vHO = 0.0;
   for(int i=0; i<part.size(); i++) {
      vHO += (part[i].pos[0]*part[i].pos[0])/2.0*omega[i]*omega[i]*reducedMass[i];
   }
   return vHO;
}

void CoordUtil::MakeScalingVec (CoordUtil* otherCoords) {
   arma::mat alpha;
   if(armaInvModes.is_empty()) {
      armaInvModes = pinv(armaNormModes);
   }
   if(otherCoords->armaNormModes.is_empty()) {
      otherCoords->armaNormModes=arma::mat(dim*numPart,numModes);
      for(int j=0; j<numPart; j++) {
         for(int k=0; k<dim; k++) {
            for(int i=0; i<numModes; i++) {
               otherCoords->armaNormModes(dim*j+k,i) = otherCoords->normModes[i][j][k];
            }
         }
      }
   }

//   (armaNormModes.t()*armaNormModes).print("mode.T*modes");
//   (otherCoords->armaNormModes.t()*otherCoords->armaNormModes).print("mode.T*modes");
//   (armaNormModes.t()*sqrtInvMass*sqrtInvMass*armaNormModes).print("mode.T*mass-1*modes");
//   (otherCoords->armaNormModes.t()*otherCoords->sqrtInvMass*otherCoords->sqrtInvMass*otherCoords->armaNormModes).print("otherMode mode.T*mass-1*modes");
//   (armaNormModes*armaNormModes.t()).print("mode*modes.T");
//   exit(-1);

//DES: Project out global rotations/translations
//TODO: This will crash on linear fragments!
   //arma::mat projector(3*numPart,3*numPart,arma::fill::zeros);
   arma::mat projector,projVecs,tempMat;
   arma::vec tempVec;
   arma::vec totalCOM(3,arma::fill::zeros);
   for(int i=0; i<numPart; i++) {
      for(int j=0; j<3; j++) {
         totalCOM(j) += atomMass(i) * initCart[i][j];
      }
   }
// global translations
   for(int n=0; n<3; n++) {
      tempVec = arma::vec(3*numPart,arma::fill::zeros);
      for(int i=0; i<numPart; i++) {
         tempVec(3*i + n) = atomMass(i);
      }
      tempVec = arma::normalise(tempVec);
      //tempVec.print("global translation");
      //projector += tempVec*tempVec.t();
      projVecs = join_rows(projVecs,tempVec);
   }
// global rotations
   for(int n=0; n<3; n++) {
      tempVec = arma::vec(3*numPart,arma::fill::zeros);
      for(int i=0; i<numPart; i++) {
         tempVec(3*i + (n+1)%3) = atomMass(i) * (initCart[i][(n+2)%3] - totalCOM((n+2)%3));
         tempVec(3*i + (n+2)%3) = -atomMass(i) * (initCart[i][(n+1)%3] - totalCOM((n+1)%3));
      }
      tempVec = arma::normalise(tempVec);
      //tempVec.print("global rotation");
      //projector += tempVec*tempVec.t();
      projVecs = join_rows(projVecs,tempVec);
   }
   //projVecs.print("projVecs");
   //arma::qr(projector,tempMat,projVecs);
   //projector.print("Q");
   arma::vec sVals;
   arma::mat leftVecs,rightVecs;
// DES: doing svd here always since only done once anyway
   svd(leftVecs,sVals,rightVecs,projVecs);
   //sVals.print("Singular values");
   //rightVecs.print("Right vecs");
   //leftVecs.print("Left vecs");

   double sValThresh = 0.000001;
   int numVecs = 6;
   for(int i=0; i<6; i++) {
      if(sVals(i)<sValThresh) {
         numVecs--;
      }
   }

   //projector = projector.cols(0,5)*projector.cols(0,5).t();
   projector = leftVecs.cols(0,numVecs-1)*leftVecs.cols(0,numVecs-1).t();

   //armaInvModes.print("modes-1 pre projector");
   armaInvModes -= armaInvModes*projector;
   //armaInvModes.print("modes-1 post projector");


   //(armaInvModes*armaNormModes).print("modes-1*modes");
   alpha = armaInvModes * otherCoords->armaNormModes;
   //alpha.print("alphaMat");
   //(alpha.t()*alpha).print("aT.a");
   alpha = arma::square(alpha);
   arma::vec tempOmega(numModes,arma::fill::zeros);
   for(int i=0; i<numModes; i++) {
      tempOmega[i] = reducedMass[i] * omega[i] * omega[i];
   }
   alpha = tempOmega.t()*alpha;
   //alpha = alpha*tempOmega;
   for(int i=0; i<numModes; i++) {
      tempOmega[i] = otherCoords->reducedMass[i] * otherCoords->omega[i] * otherCoords->omega[i];
   }
   scalingVec = arma::sqrt(tempOmega/alpha.t());
   
   scalingVec.print("ScalingVec");
   //exit(-1);
}

void CoordUtil::MatchModes (CoordUtil* matchCoords) {
   arma::mat S(numModes,numModes,arma::fill::zeros);
   arma::mat absS(numModes,numModes);
   for(int i=0; i<numModes; i++) {
      for(int ip=0; ip<numModes; ip++) {
         for(int j=0; j<numPart; j++) {
            for(int k=0; k<3; k++) {
               S(i,ip) += normModes[i][j][k]*matchCoords->normModes[ip][j][k];
            }
         }
      }
   }
   //S.print("S");
   arma::uword row;
   arma::uword col;
   arma::vec overlap(numModes,arma::fill::zeros);
   arma::uvec order(numModes,arma::fill::zeros);
   for(int i=0; i<numModes; i++) {
      absS = arma::abs(S);
      absS.max(row,col);
      order(row) = col;
      overlap(col) = S(row,col);
      S.row(row).zeros();
      S.col(col).zeros();
      //S.print("S");
   }
   //overlap.print("overlap");
   //order.print("order");
// Reorder normModes, omega, reducedMass
   for ( int s = 1, d; s < order.size(); ++ s ) {
      for ( d = order[s]; d < s; d = order[d] ) ;
         if ( d == s ) {
            while ( d = order[d], d != s ) {
               swap( omega[s], omega[d] );
               swap( reducedMass[s], reducedMass[d] );
               swap( normModes[s], normModes[d] );
            }
         }
   }
   for(int i=0; i<numModes; i++) {
      if(fabs(overlap(i)) < 0.5) {
         cout << "Warning! Bad overlap between tinker and qchem normal modes: " << fabs(overlap(i)) << endl;
         cout << "omega values:\t" << omega[i] << ",\t" << matchCoords->omega[i] << endl;
      }
      if(overlap(i) < 0.0) {
         for(int j=0; j<numPart; j++) {
            for(int k=0; k<3; k++) {
               normModes[i][j][k] *= -1.0;
            }
         }
      }
   }
//   cout << "omega" << endl;
//   for(int i=0; i<numModes; i++) {
//      cout << omega[i] << endl;
//   }
   armaInvModes.reset();
}

CoordUtil* CoordUtil::Clone() {
   CoordUtil* newCoords = new CoordUtil(numModes,numPart,internals,normModes,omega,reducedMass,initCart,atomType,paramType,connectivity,outFileName,prmName,readGeom,readOmega);
   newCoords->scaleGeom = scaleGeom;
   return newCoords;
}

void CoordUtil::initMass() {
   atomToMass["H"] = 1.0079;
   atomToMass["He"] = 4.0026;                                                             
   atomToMass["Li"] = 6.941;                                                              
   atomToMass["Be"] = 9.0122;                                                             
   atomToMass["B"] = 10.811;                                                              
   atomToMass["C"] = 12.0107;                                                             
   atomToMass["N"] = 14.0067;                                                             
   atomToMass["O"] = 15.9994;                                                             
   atomToMass["F"] = 18.9984;                                                             
   atomToMass["Ne"] = 20.1797;                                                            
   atomToMass["Na"] = 22.9897;                                                            
   atomToMass["Mg"] = 24.305;                                                             
   atomToMass["Al"] = 26.9815;                                                            
   atomToMass["Si"] = 28.0855;                                                            
   atomToMass["P"] = 30.9738;                                                             
   atomToMass["S"] = 32.065;                                                              
   atomToMass["Cl"] = 35.453;                                                             
   atomToMass["K"] = 39.0983;                                                             
   atomToMass["Ar"] = 39.948;                                                             
   atomToMass["Ca"] = 40.078;                                                             
   atomToMass["Sc"] = 44.9559;                                                            
   atomToMass["Ti"] = 47.867;                                                             
   atomToMass["V"] = 50.9415;                                                             
   atomToMass["Cr"] = 51.9961;                                                            
   atomToMass["Mn"] = 54.938;                                                             
   atomToMass["Fe"] = 55.845;                                                             
   atomToMass["Ni"] = 58.6934;                                                            
   atomToMass["Co"] = 58.9332;                                                            
   atomToMass["Cu"] = 63.546;                                                             
   atomToMass["Zn"] = 65.39;                                                              
   atomToMass["Ga"] = 69.723;                                                             
   atomToMass["Ge"] = 72.64;                                                              
   atomToMass["As"] = 74.9216;                                                            
   atomToMass["Se"] = 78.96;                                                              
   atomToMass["Br"] = 79.904;                                                             
   atomToMass["Kr"] = 83.8;                                                               
   atomToMass["Rb"] = 85.4678;                                                            
   atomToMass["Sr"] = 87.62;                                                              
   atomToMass["Y"] = 88.9059;                                                             
   atomToMass["Zr"] = 91.224;                                                             
   atomToMass["Nb"] = 92.9064;                                                            
   atomToMass["Mo"] = 95.94;                                                              
   atomToMass["Tc"] = 98;                                                                 
   atomToMass["Ru"] = 101.07;                                                             
   atomToMass["Rh"] = 102.9055;                                                           
   atomToMass["Pd"] = 106.42;                                                             
   atomToMass["Ag"] = 107.8682;                                                           
   atomToMass["Cd"] = 112.411;                                                            
   atomToMass["In"] = 114.818;                                                            
   atomToMass["Sn"] = 118.71;                                                             
   atomToMass["Sb"] = 121.76;                                                             
   atomToMass["I"] = 126.9045;                                                            
   atomToMass["Te"] = 127.6;                                                              
   atomToMass["Xe"] = 131.293;                                                            
   atomToMass["Cs"] = 132.9055;
}

int CoordUtil::atomStrToInt(string atomStr) {
if ("H"==atomStr) return  1  ;
if ("He"==atomStr) return  2  ;
if ("Li"==atomStr) return  3  ;
if ("Be"==atomStr) return  4  ;
if ("B"==atomStr) return  5  ;
if ("C"==atomStr) return  6  ;
if ("N"==atomStr) return  7  ;
if ("O"==atomStr) return  8  ;
if ("F"==atomStr) return  9  ;
if ("Ne"==atomStr) return  10 ;
if ("Na"==atomStr) return  11 ;
if ("Mg"==atomStr) return  12 ;
if ("Al"==atomStr) return  13 ;
if ("Si"==atomStr) return  14 ;
if ("P"==atomStr) return  15 ;
if ("S"==atomStr) return  16 ;
if ("Cl"==atomStr) return  17 ;
if ("K"==atomStr) return  18 ;
if ("Ar"==atomStr) return  19 ;
if ("Ca"==atomStr) return  20 ;
if ("Sc"==atomStr) return  21 ;
if ("Ti"==atomStr) return  22 ;
if ("V"==atomStr) return  23 ;
if ("Cr"==atomStr) return  24 ;
if ("Mn"==atomStr) return  25 ;
if ("Fe"==atomStr) return  26 ;
if ("Ni"==atomStr) return  27 ;
if ("Co"==atomStr) return  28 ;
if ("Cu"==atomStr) return  29 ;
if ("Zn"==atomStr) return  30 ;
if ("Ga"==atomStr) return  31 ;
if ("Ge"==atomStr) return  32 ;
if ("As"==atomStr) return  33 ;
if ("Se"==atomStr) return  34 ;
if ("Br"==atomStr) return  35 ;
if ("Kr"==atomStr) return  36 ;
if ("Rb"==atomStr) return  37 ;
if ("Sr"==atomStr) return  38 ;
if ("Y"==atomStr) return  39 ;
if ("Zr"==atomStr) return  40 ;
if ("Nb"==atomStr) return  41 ;
if ("Mo"==atomStr) return  42 ;
if ("Tc"==atomStr) return  43 ;
if ("Ru"==atomStr) return  44 ;
if ("Rh"==atomStr) return  45 ;
if ("Pd"==atomStr) return  46 ;
if ("Ag"==atomStr) return  47 ;
if ("Cd"==atomStr) return  48 ;
if ("In"==atomStr) return  49 ;
if ("Sn"==atomStr) return  50 ;
if ("Sb"==atomStr) return  51 ;
if ("I"==atomStr) return  52 ;
if ("Te"==atomStr) return  53 ;
if ("Xe"==atomStr) return  54 ;
if ("Cs"==atomStr) return  55 ;
else return -1;
}
