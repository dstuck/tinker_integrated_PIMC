/*
 * RhoFree.cpp
 *
 *  Created on: May 9, 2012
 *      Author: dstuck
 */

#include "Rho_Free.h"

Rho_Free::Rho_Free() {
	// TODO Auto-generated constructor stub

}

Rho_Free::~Rho_Free() {
	// TODO Auto-generated destructor stub
}

double Rho_Free::GetRho(vector<Particle> slice1,vector<Particle> slice2, double eps) {
	double delta;
	double rho = 0;
	for(int j=0; j<(int)slice1.size(); j++) {
		delta = 0;
		for(int k=0; k<(int)slice1[j].pos.size(); k++) {
			delta += (slice1[j].pos[k]-slice2[j].pos[k])*(slice1[j].pos[k]-slice2[j].pos[k]);
		}
		rho += slice1[j].mass/2.0/eps/eps*(delta);
	}
//	cout << "delta: " << delta << endl;
	return rho;
}

double Rho_Free::ModifyPotential(vector<Particle> part) {
	return 0;
}

double Rho_Free::Estimate(vector<Particle> slice1,vector<Particle> slice2, double eps, int P) {
	double delta;
	double kin;
	double est = 0;
	for(int j=0; j<(int)slice1.size(); j++) {
		delta = 0;
		kin = 0;
		for(int k=0; k<(int)slice1[j].pos.size(); k++) {
			delta += -(slice1[j].pos[k]-slice2[j].pos[k])*(slice1[j].pos[k]-slice2[j].pos[k]);
			kin += 1.0/2.0/eps;
		}
		est += (slice1[j].mass/2.0/eps/eps*(delta)+kin)/(double)P;
	}
	return est;
}

vector<double> Rho_Free::GetSpringLength(vector<Particle> part, double eps) {
	cout << "Shouldn't be calling GetSpringLength" << endl;
	vector<double> r0;
	for(int j=0; j<(int)part.size(); j++) {
		r0.push_back(sqrt(eps/part[j].mass));
	}
	return r0;
}

vector<double> Rho_Free::GetLevyMean(Particle partI, Particle partF, double delta, double eps, int unused) {
	vector<double> mean;
	for(int k = 0; k<(int)partI.pos.size(); k++) {
		mean.push_back( (partI.pos[k]*delta + partF.pos[k]) / (delta+1.0) );
	}
	return mean;
}

double Rho_Free::GetLevySigma(double delta, double eps, int unused) {
	double sigma;
	sigma = sqrt(eps/(1.0+1.0/delta));
	return sigma;
}


string Rho_Free::GetType() {
	string name = "Free";
	return name;
}
