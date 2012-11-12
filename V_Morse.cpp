/*
 * V_Morse.cpp
 *
 *  Created on: May 24, 2012
 *      Author: dstuck
 */

#include "V_Morse.h"

V_Morse::V_Morse() {
}

V_Morse::V_Morse(double deVar, double aVar, int N) {
	for(int j=0; j<N; j++) {
		de.push_back(deVar);
		a.push_back(aVar);
	}
}

V_Morse::V_Morse(vector<double> deVar, vector<double> aVar) {
	if(deVar.size()!=aVar.size()){
		cout << "!!!Error: De and a dimensions different in V_Morse initialization!!!" << endl;
	}
	de = deVar;
	a = aVar;
}

V_Morse::~V_Morse() {
	// TODO Auto-generated destructor stub
}

double V_Morse::GetV(vector<Particle> part, Propagator * rho) {
//	Really, should only have dimension 1
	double V = 0;
	int dim = part[0].pos.size();
	double expTemp = 0;
	for(int j=0; j<(int)part.size(); j++) {
		for(int k=0; k<dim; k++) {
			expTemp = exp(-a[j]*part[j].pos[k]);
			V += de[j]*(1.0 - expTemp)*(1.0-expTemp);
		}
	}
//	cout << "preV: " << V << endl;
	V += rho->ModifyPotential(part);
//	cout << "V: " << V << endl;

	return V;
}

string V_Morse::GetType() {
	string name = "Morse";
	return name;
}
