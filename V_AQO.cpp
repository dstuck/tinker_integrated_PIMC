/*
 * VAQO.cpp
 *
 *  Created on: May 24, 2012
 *      Author: dstuck
 *	Anharmonic quartic oscillator: V(x) = 1/2*m*w^2*x^2 + lambda*x^4
 */

#include "V_AQO.h"

V_AQO::V_AQO() {
}

V_AQO::V_AQO(double w, double l, int N) {
	for(int j=0; j<N; j++) {
		omega.push_back(w);
		lambda.push_back(l);
	}
}

V_AQO::V_AQO(vector<double> w, vector<double> l) {
	if(w.size()!=l.size()){
		cout << "!!!Error: omega and lambda dimensions different in V_AQO initialization!!!" << endl;
	}
	omega = w;
	lambda = l;
}

double V_AQO::GetV(vector<Particle> part, Propagator * rho) {
	V = 0;
        V += GetV(part);
	V += rho->ModifyPotential(part);
	return V;
}
double V_AQO::GetV(vector<Particle> part) {
	V = 0;
	int dim = part[0].pos.size();
	for(int j=0; j<(int)part.size(); j++) {
		for(int k=0; k<dim; k++) {
			V += (part[j].pos[k])*(part[j].pos[k])/2.0*omega[j]*omega[j]*part[j].mass + (part[j].pos[k])*(part[j].pos[k])*(part[j].pos[k])*(part[j].pos[k])*lambda[j];
		}
	}
	return V;
}

string V_AQO::GetType() {
	string name = "AQO";
	return name;
}

CoordUtil* V_AQO::GetCoordUtil() {
        cout << "ERROR: CoordUtil not supported in UCHO" << endl;
        return NULL;
}

V_AQO::~V_AQO() {
	// TODO Auto-generated destructor stub
}

