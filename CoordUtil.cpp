/*
 * CoordUtil.cpp
 *
 *  Created on: Sep 16, 2012
 *      Author: dstuck
 */

#include "CoordUtil.h"

CoordUtil::CoordUtil() {
}

CoordUtil::CoordUtil(vector <vector < vector <double> > > modes, vector<double> freqs, vector <vector <double> > initPos, vector<string> atomicSymbols, vector<int> params, vector< vector<int> > conn, string tinkName, string prmFile) : normModes(modes) {
//	normModes = modes;				//normModes[i][j][k] is kth dimension of jth atom for the ith normal mode
//        omega = freqs;
         omega = vector<double>(3,0);
	initCart = initPos;				//initPos[i][j] is jth dimension of ith atom
	connectivity = conn;			//connectivity[i][c] is the number of the cth atom connected to atom i (not square)
	tinkerName = tinkName;
	prmName = prmFile;
	paramType = params;
	atomType = atomicSymbols;
	numModes = modes.size();
	numPart = modes[0].size();
	dim = modes[0][0].size();
	if(dim != 3) {
		cout << "Warning: Not working in xyz coordinates!" << endl;
	}
}


CoordUtil::~CoordUtil() {
	// TODO Auto-generated destructor stub
}

vector< vector<double> > CoordUtil::normalModeToCart(vector<Particle> part) {
	vector< vector<double> > carts;
	carts = initCart;
//	For Normal modes
	if(part[0].pos.size()==1) {
		for(int j=0; j<numPart; j++) {
			for(int k=0; k<dim; k++) {
				for(int i=0; i<numModes; i++) {
					carts[j][k] += normModes[i][j][k]*part[i].pos[0] * 0.52918;				// Converts to A from Bohr radii, value from wikipedia 9/18/12
				}
			}
		}
	}
//	For cartesian coordinates
	else if(part[0].pos.size()==3) {
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
