/*
 * CoordUtil.cpp
 *
 *  Created on: Sep 16, 2012
 *      Author: dstuck
 */

#include "CoordUtil.h"
#include <stdlib.h> 

CoordUtil::CoordUtil() {
}

CoordUtil::CoordUtil(int nMode, int nPart, vector <vector < vector <double> > > modes, vector<double> freqs, vector<double> m, vector <vector <double> > initPos, vector<string> atomicSymbols, vector<int> params, vector< vector<int> > conn, string tinkName, string prmFile, bool readW) : numModes(nMode), numPart(nPart), normModes(modes), omega(freqs), reducedMass(m), initCart(initPos), atomType(atomicSymbols), paramType(params), connectivity(conn), tinkerName(tinkName), prmName(prmFile), readOmega(readW)  {
//	normModes = modes;				//normModes[i][j][k] is kth dimension of jth atom for the ith normal mode
//        omega = freqs;
//	initCart = initPos;				//initPos[i][j] is jth dimension of ith atom
//	connectivity = conn;			//connectivity[i][c] is the number of the cth atom connected to atom i (not square)
//	tinkerName = tinkName;
//	prmName = prmFile;
//	paramType = params;
//	atomType = atomicSymbols;
//	numModes = modes.size();
//	numPart = modes[0].size();
//	dim = modes[0][0].size();
        dim = 3;
	if(dim != 3) {
		cout << "Warning: Not working in xyz coordinates!" << endl;
	}
   initMass();
}

CoordUtil::~CoordUtil() {
	// TODO Auto-generated destructor stub
}

vector< vector<double> > CoordUtil::normalModeToCart(vector<Particle> part) {
// Really converts any part to cartesians (normal or already cartesian)
	vector< vector<double> > carts(initCart);
//	for(int i=0; i<(int)initCart.size(); i++) {
//            carts.push_back(initCart[i]);
//        }
//DES Temp:
/*
   cout << "DES: normConst" << endl;
   for(int i=0; i<numModes; i++) {
      cout << normConst[i] << endl;
   }
   vector<double> overlap;
   double dot = 0.0;
   for(int n=0; n<numModes; n++) {
      for(int m=0; m<numModes; m++) {
         for(int j=0; j<numPart; j++) {
            for(int k=0; k<dim; k++) {
               //dot += normModes[n][j][k]*normModes[m][j][k]*atomToMass[atomType[j]];
               dot += normModes[n][j][k]*normModes[m][j][k]*atomToMass[atomType[j]]*normConst[n]*normConst[m];
            }
         }
         overlap.push_back(dot);
         dot = 0.0;
      }
   }
   cout << "DES: Overlap matrix" << endl;
   for(int n=0; n<numModes; n++) {
      for(int m=0; m<numModes; m++) {
         cout << overlap[m+numModes*n] << "\t";
      }
      cout << endl;
   }
   cout << endl;
   for(int i=0; i<4; i++) {
      for(int j=0; j<numPart; j++) {
         cout << normModes[i][j][0] << "\t" << normModes[i][j][1] << "\t" << normModes[i][j][2] << endl;
      }
      cout << endl;
   }
   exit(-1);
*/


//	For Normal modes into Cartesians in Angstroms
	if(part[0].pos.size()==1) {
/*
            arma::vec armaCart(numPart*dim);
            arma::mat armaModes(numPart*dim,numModes);
            arma::vec armaVec(numModes);
            arma::vec armaInitCart(numPart*dim);
            for(int j=0; j<numPart; j++) {
               for(int k=0; k<dim; k++) {
                  armaInitCart(3*j+k) = carts[j][k];
                  for(int i=0; i<numModes; i++) {
                     armaModes(3*j+k,i) = normModes[i][j][k];
                  }
               }
            }
            //armaInitCart.print("armaInitCart");
            for(int i=0; i<numModes; i++) {
               armaVec(i) = part[i].pos[0];
            }
            armaVec.print("armaVec in ntoc");
            armaCart.print("armaCart in ntoc");
            armaCart = armaModes*armaVec * 0.52918 + armaInitCart;
            arma::mat armaInvModes = arma::pinv(armaModes);
            armaVec = armaInvModes*(armaCart-armaInitCart)/0.52918;
            //armaVec.print("reModes in ntoc");
            //armaModes.print("armaModes in ntoc");
*/
///*
		for(int j=0; j<numPart; j++) {
			for(int k=0; k<dim; k++) {
				for(int i=0; i<numModes; i++) {
					carts[j][k] += normModes[i][j][k]*part[i].pos[0] * 0.52918;				// Converts to A from Bohr radii, value from wikipedia 9/18/12
				}
			}
		}
//*/

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
   if(armaInvModes.is_empty()) {
      arma::mat armaModes(numPart*dim,numModes);
      for(int j=0; j<numPart; j++) {
         for(int k=0; k<dim; k++) {
            armaCart(3*j+k) = (cart[j][k]-initCart[j][k])/0.52918;
            for(int i=0; i<numModes; i++) {
               armaModes(3*j+k,i) = normModes[i][j][k];
            }
         }
      }
      armaInvModes = arma::pinv(armaModes);
   }
   else {
      for(int j=0; j<numPart; j++) {
         for(int k=0; k<dim; k++) {
            armaCart(3*j+k) = (cart[j][k]-initCart[j][k])/0.52918;
         }
      }
   }
   armaVecs = armaInvModes*armaCart;
   armaVecs = armaInvModes*armaCart;
   nmodes = arma::conv_to< vector<double> >::from(armaVecs);

   return nmodes;
}

CoordUtil* CoordUtil::Clone() {
   CoordUtil* newCoords = new CoordUtil(numModes,numPart,normModes,omega,reducedMass,initCart,atomType,paramType,connectivity,tinkerName,prmName,readOmega);
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
