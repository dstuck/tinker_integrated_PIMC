/*
 * V_Morse.cpp
 *
 *  Created on: May 24, 2012
 *      Author: dstuck
 */

#include "V_Morse.h"

V_Morse::V_Morse(CoordUtil* coords, double deVar, double aVar, double mVar, int N) : coordKeeper(coords) {
   vector<double> mass;
   vector<double> omega;
	for(int j=0; j<N; j++) {
                mass.push_back(mVar);
		de.push_back(deVar);
		a.push_back(aVar);
                omega.push_back(aVar*sqrt(2.0*deVar/mVar));
	}
   coords->reducedMass = mass;
   coords->omega = omega;
}

V_Morse::V_Morse(CoordUtil* coords, vector<double> deVar, vector<double> aVar, vector<double> mVar) : coordKeeper(coords) {
	if(deVar.size()!=aVar.size()){
		cout << "!!!Error: De and a dimensions different in V_Morse initialization!!!" << endl;
	}
	de = deVar;
	a = aVar;
   vector<double> omega;
   for(int i=0; i<deVar.size(); i++) {
      omega.push_back(aVar[i]*sqrt(2.0*deVar[i]/mVar[i]));
   }
   coords->reducedMass = mVar;
   coords->omega = omega;
}

V_Morse::~V_Morse() {
	// TODO Auto-generated destructor stub
}

double V_Morse::GetV(vector<Particle> part, Propagator * rho) {
	double V = 0;
        V += GetV(part);
//	cout << "preV: " << V << endl;
	V += rho->ModifyPotential(part);
//	cout << "V: " << V << endl;

	return V;
}

double V_Morse::GetV(vector<Particle> part) {
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

	return V;
}

string V_Morse::GetType() {
	string name = "Morse";
	return name;
}

CoordUtil* V_Morse::GetCoordUtil() {
        return coordKeeper;
}
