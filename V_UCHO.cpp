/*
 * V_UCHO.cpp
 *
 *  Created on: May 9, 2012
 *      Author: dstuck
 */

#include "V_UCHO.h"

V_UCHO::V_UCHO() {
	// TODO Auto-generated constructor stub
}

V_UCHO::V_UCHO(vector<double> w) {
	omega = w;
}

V_UCHO::V_UCHO(double w, int num) {
	for(int n=0; n<num; n++) {
		omega.push_back(w);
	}
}

V_UCHO::~V_UCHO() {
	// TODO Auto-generated destructor stub
}

double V_UCHO::GetV(vector<Particle> part, Propagator * rho) {
	V = 0;
//	cout <<"w = ";
//	for(int i = 0; i < int(omega.size()); i++) {
//		cout << omega[i] << ", ";
//	}
//	cout << endl;
//	cout <<"m = ";
//		for(int i = 0; i < int(part.size()); i++) {
//			cout << part[i].mass << ", ";
//		}
//	cout << endl;
	int dim = part[0].pos.size();
	for(int j=0; j<(int)part.size(); j++) {
		for(int k=0; k<dim; k++) {
			V += (part[j].pos[k])*(part[j].pos[k])/2.0*omega[j]*omega[j]*part[j].mass;
		}
	}
//	cout << "V: " << V << endl;
	V += rho->ModifyPotential(part);
//	cout << "V: " << V << endl;
	return V;
}

string V_UCHO::GetType() {
	string name = "UCHO";
	return name;
}
