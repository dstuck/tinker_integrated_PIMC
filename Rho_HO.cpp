/*
 * RhoHO.cpp
 *
 *  Created on: May 17, 2012
 *      Author: dstuck
 */

#include "Rho_HO.h"

Rho_HO::Rho_HO(vector<double> w) {
	omega = w;
//      for(int i=0; i<(int)omega.size(); i++) {
//         cout << omega[i]/0.00000455633 << endl;
//      }
}

Rho_HO::Rho_HO(double w, int N) {
        omega.clear();
	for(int j=0; j<N; j++){
		omega.push_back(w);
	}
}

Rho_HO::~Rho_HO() {
	// TODO Auto-generated destructor stub
}

double Rho_HO::GetRho(vector<Particle> slice1, vector<Particle> slice2, double eps) {
// 	This function is not called if using Levy Flights
	double rho = 0;
	vector<double> gamma;
	if(omega.size() != slice1.size()) {
		cout << "Error! Number of particles not equal to number of frequencies! Exiting with -99999." << endl;
		return -99999;
	}
	for(int j=0; j<(int)slice1.size(); j++) {
		for(int k=0; k<(int)slice1[j].pos.size(); k++) {
//			rho += slice1[j].mass*omega[j]/2.0/eps * ((slice1[j].pos[k] - slice2[j].pos[k])*(slice1[j].pos[k] - slice2[j].pos[k])/eps/omega[j]/sinh(omega[j]*eps) + tanh(eps*omega[j]/2.0)/omega[j]/eps*(slice1[j].pos[k]*slice1[j].pos[k]+slice2[j].pos[k]*slice2[j].pos[k]));
			rho += slice1[j].mass*omega[j]/2.0/eps * ((slice1[j].pos[k] - slice2[j].pos[k])*(slice1[j].pos[k] - slice2[j].pos[k])/sinh(omega[j]*eps) + tanh(eps*omega[j]/2.0)*(slice1[j].pos[k]*slice1[j].pos[k]+slice2[j].pos[k]*slice2[j].pos[k]));
		}
	}
	return rho;
}

double Rho_HO::ModifyPotential(vector<Particle> part) {
	double dV = 0.0;
	int dim = part[0].pos.size();
	for(int j=0; j<(int)part.size(); j++) {
		for(int k=0; k<dim; k++) {
			dV -= (part[j].pos[k])*(part[j].pos[k])/2.0*omega[j]*omega[j]*part[j].mass;
		}
	}
//	cout << "\tV_HO =\t" << dV << endl;
//	cout << dV << endl;
	return dV;
}

double Rho_HO::Estimate(vector<Particle> slice1, vector<Particle> slice2, double eps, int P) {
/*
 * 	Estimator taken from Whitfield and Martyna, J. Chem. Phys. 126, 074104 (2007) with missing factor of 1/2 added to term 2
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
			est += slice1[j].mass*omega[j]*omega[j]/2.0/cosh(eps*omega[j]/2.0)/cosh(eps*omega[j]/2.0)*(slice1[j].pos[k]*slice1[j].pos[k]);
		}
	}
	est /= (double)P;
	return est;
}


vector<double> Rho_HO::GetSpringLength(vector<Particle> part, double eps) {
	cout << "Should not be calling this function!" << endl;
	vector<double> r0;
	for(int j=0; j<(int)part.size(); j++) {
		r0.push_back(sqrt(tanh(eps*omega[j])/part[j].mass/omega[j]));
//		r0.push_back(sqrt(sinh(eps*omega[j])/part[j].mass/omega[j]));
	}
	return r0;
}


vector<double> Rho_HO::GetLevyMean(Particle partI, Particle partF, double delta, double eps, int partNum) {
	vector<double> mean;
	double tDelI = tanh(eps*omega[partNum]);
	double tDelF = tanh(delta*eps*omega[partNum]);
	double sDelI = sinh(eps*omega[partNum]);
	double sDelF = sinh(delta*eps*omega[partNum]);
	for(int k = 0; k<(int)partI.pos.size(); k++) {
		mean.push_back( (partI.pos[k]/sDelI + partF.pos[k]/sDelF) * tDelI*tDelF/(tDelI+tDelF) );
	}
	return mean;
}

double Rho_HO::GetLevySigma(double delta, double eps, int partNum) {
	double sigma;
	double tDelI = tanh(eps*omega[partNum]);
	double tDelF = tanh(delta*eps*omega[partNum]);
	sigma = sqrt(tDelI*tDelF/(tDelI+tDelF) / omega[partNum]);
	return sigma;
}


string Rho_HO::GetType() {
	string name = "HO";
	return name;
}
