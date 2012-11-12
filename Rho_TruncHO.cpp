/*
 * Rho_TruncHO.cpp
 *
 *  Created on: Jun 13, 2012
 *      Author: dstuck
 */

#include "Rho_TruncHO.h"

Rho_TruncHO::Rho_TruncHO(vector<double> w) {
	omega = w;
}

Rho_TruncHO::Rho_TruncHO(double w, int N) {
	for(int j=0; j<N; j++){
		omega.push_back(w);
	}
}

Rho_TruncHO::~Rho_TruncHO() {
	// TODO Auto-generated destructor stub
}

double Rho_TruncHO::GetRho(vector<Particle> slice1, vector<Particle> slice2, double eps) {
// 	This function is not called if using Levy Flights
	cout << "This shouldn't be called" << endl;
	double rho = 0;
	vector<double> gamma;
	if(omega.size() != slice1.size()) {
		cout << "Error! Number of particles not equal to number of frequencies! Exiting with -99999." << endl;
		return -99999;
	}
	for(int j=0; j<(int)slice1.size(); j++) {
		for(int k=0; k<(int)slice1[j].pos.size(); k++) {
			rho += slice1[j].mass*omega[j]/2.0/eps * ((slice1[j].pos[k] - slice2[j].pos[k])*(slice1[j].pos[k] - slice2[j].pos[k])/sinh(omega[j]*eps));
//			rho += slice1[j].mass*omega[j]/2.0/eps * tanh(eps*omega[j]/2.0)/omega[j]/eps*(slice1[j].pos[k]*slice1[j].pos[k]+slice2[j].pos[k]*slice2[j].pos[k]);		//TODO: remove this
		}
	}
	return rho;
}

double Rho_TruncHO::ModifyPotential(vector<Particle> part) {
//	return 0;
	double dV = 0;
	int dim = part[0].pos.size();
	for(int j=0; j<(int)part.size(); j++) {
		for(int k=0; k<dim; k++) {
			dV -= (part[j].pos[k])*(part[j].pos[k])/2.0*omega[j]*omega[j]*part[j].mass;
		}
	}
	return dV;
}

double Rho_TruncHO::Estimate(vector<Particle> slice1, vector<Particle> slice2, double eps, int P) {
/*
 * 	Estimator explained in Whitfield and Martyna, J. Chem. Phys. 126, 074104 (2007) with missing factor of 1/2 added to term 2
 */
	if(omega.size() != slice1.size()) {
		cout << "Error! Number of particles not equal to number of frequencies! Exiting with -99999." << endl;
		return -99999;
	}

	double est = 0;
	for(int j=0; j<(int)slice1.size(); j++) {
		for(int k=0; k<(int)slice1[j].pos.size(); k++) {
			est += omega[j]/2.0/tanh(eps*omega[j]);
			est += -slice1[j].mass*omega[j]*omega[j]/2.0/tanh(eps*omega[j])/sinh(eps*omega[j])*(slice1[j].pos[k]-slice2[j].pos[k])*(slice1[j].pos[k]-slice2[j].pos[k]);
			est += slice1[j].mass*omega[j]*omega[j]/2.0/cosh(eps*omega[j]/2.0)/cosh(eps*omega[j]/2.0)*(slice1[j].pos[k]*slice1[j].pos[k]);			//Divergent part of energy
		}
	}
//	Test
//	int dim = slice1[0].pos.size();
//	for(int j=0; j<(int)slice1.size(); j++) {
//		for(int k=0; k<dim; k++) {
//			est -= (slice1[j].pos[k])*(slice1[j].pos[k])/2.0*omega[j]*omega[j]*slice1[j].mass;
//		}
//	}
	est /= (double)P;
	return est;
}


vector<double> Rho_TruncHO::GetSpringLength(vector<Particle> part, double eps) {
	vector<double> r0;
	for(int j=0; j<(int)part.size(); j++) {
		cout << "This shouldn't be called" << endl;
		r0.push_back(sqrt(sinh(eps*omega[j])/part[j].mass/omega[j]));
	}
	return r0;
}

vector<double> Rho_TruncHO::GetLevyMean(Particle partI, Particle partF, double delta, double eps, int partNum) {
	vector<double> mean;
	bool fullLevy = true;
	if(fullLevy) {
		double tDelI = tanh(eps*omega[partNum]);
		double tDelF = tanh(delta*eps*omega[partNum]);
		double sDelI = sinh(eps*omega[partNum]);
		double sDelF = sinh(delta*eps*omega[partNum]);
		for(int k = 0; k<(int)partI.pos.size(); k++) {
			mean.push_back( (partI.pos[k]/sDelI + partF.pos[k]/sDelF) * tDelI*tDelF/(tDelI+tDelF) );
		}
	}
	else {
		for(int k = 0; k<(int)partI.pos.size(); k++) {
			mean.push_back( (partI.pos[k]*delta + partF.pos[k]) / (delta+1.0) );
		}
	}
	return mean;
}

double Rho_TruncHO::GetLevySigma(double delta, double eps, int partNum) {
	double sigma;
	bool fullLevy = true;
	if(fullLevy) {
		double tDelI = tanh(eps*omega[partNum]);
		double tDelF = tanh(delta*eps*omega[partNum]);
		sigma = sqrt(tDelI*tDelF/(tDelI+tDelF) / omega[partNum]);
	}
	else {
		sigma = sqrt(sinh(eps*omega[partNum])/(1.0+1.0/delta));
	}
	return sigma;
}


string Rho_TruncHO::GetType() {
	string name = "TruncHO";
	return name;
}
