/*
 * V_TinkerExecutable.cpp
 *
 *  Created on: Sep 7, 2012
 *      Author: dstuck
 */

#include "V_TinkerExecutable.h"
#include <stdlib.h>

V_TinkerExecutable::V_TinkerExecutable(CoordUtil * coords) {
//	tinkInFileName = "pimc.xyz";
//	tinkOutFileName = "pimc.out";
	tinkInFileName = coords->outFileName + ".xyz";
	tinkOutFileName = coords->outFileName + ".outTinker";
	tinkPrmFileName = coords->prmName;
	coordKeeper = coords;

	vector<Particle> parts;
	Particle p0(0.0, 1);
	Propagator * rhoFree = new Rho_Free;
	for(int i=0; i<coordKeeper->numModes; i++) {
		parts.push_back(p0);
	}
	vEquib = 0.0;
	vEquib = GetV(parts, rhoFree);

//	cout << "vEquib = " << vEquib << endl;
}

V_TinkerExecutable::~V_TinkerExecutable() {
}

double V_TinkerExecutable::GetV(vector<Particle> part, Propagator * rho){
	double V = 0;
        V += GetV(part);
//	cout << V << endl;
	V += rho->ModifyPotential(part);
//	cout << "Modified V =\t" << V << endl;
	return V;
}

double V_TinkerExecutable::GetV(vector<Particle> part){

	double hartreeToKcal = 627.509469;						//From http://en.wikipedia.org/wiki/Hartree 8/17/2012

	inFile.open(tinkInFileName.c_str());
	outFile.open(tinkOutFileName.c_str());

	double V = 0;
	int N = coordKeeper->numPart;
	//	Get cartesians from normal modes
	vector< vector<double> > cartPos = coordKeeper->normalModeToCart(part);

//	Write tinker inputfile
	inFile << N << endl;
	for(int i=0; i<N; i++) {
		inFile << i+1 << "\t" << coordKeeper->atomType[i] << "\t" << cartPos[i][0] << "\t" << cartPos[i][1] << "\t" << cartPos[i][2];
		for(int c=0; c<int(coordKeeper->connectivity[i].size()); c++) {
				inFile << "\t" << coordKeeper->connectivity[i][c];
		}
		inFile << endl;
	}

//	Make analyze call
	string command = "/Users/dstuck/projects/PIMC/executablesTinker/analyze "+tinkInFileName+" "+tinkPrmFileName+" E > "+tinkOutFileName;
	system(command.c_str());

	string lineRead;
	bool found = false;
	while(!found && !(getline(outFile, lineRead).eof())) {
		vector<string> lineTokens;
		Tokenize(lineRead, lineTokens);
		if(!(lineTokens.empty())) {
			if(!(lineTokens[0].compare("Total"))&&!(lineTokens[1].compare("Potential"))) {			//if on the potential energy line
				V = atof(lineTokens[4].c_str()) / hartreeToKcal;
				found = true;
			}
		}
	}
	V -= vEquib;

	inFile.close();
	outFile.close();

	return V;
}

string V_TinkerExecutable::GetType() {
	string name = "Tinker";
	return name;
}

CoordUtil* V_TinkerExecutable::GetCoordUtil() {
        return coordKeeper;
}



//	Utility function stolen from http://www.oopweb.com/CPP/Documents/CPPHOWTO/Volume/C++Programming-HOWTO-7.html
void V_TinkerExecutable::Tokenize(const string& str, vector<string>& tokens, const string& delimiters)
{
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}
