/*
 * ClassicalSys.cpp
 *
 *  Created on: May 7, 2012
 *      Author: dstuck
 */

#include "TestSystem.h"

TestSystem::TestSystem(int nPart, int pSlice, int dim, vector<string> aType, double beta, vector<double> m, CoordUtil coords, PhysicsUtil physics) : System() {
//	cout << "nPart = " << nPart << "\npSlice = " << pSlice << "\ndim = " << dim << "\nbeta = " << beta << endl;
//	cout <<"m = ";
//	for(int i = 0; i < int(m.size()); i++) {
//		cout << m[i] << ", ";
//	}
//	cout << endl;
	N = nPart;
	P = pSlice;
	eps = beta/((double)P);
	coorDim = dim;
	//vector<Particle> * part_slice = new vector<Particle>(N);
	vector< vector<Particle> > part_init(P,vector<Particle>(N));
	part = part_init;
	oldPart = part_init;
	upToDate.resize(P, true);
	sliceV.resize(P, 0.0);
	oldSliceV.resize(P, 0.0);
	avgV = 0.0;
	numSteps = 0;
      cout << "Size of omega: " << (int)coords.omega.size() << endl;
      cout << "Size of coords.omega: " << (int)coords.omega.size() << endl;
//        for(int i=0; i<(int)w.size(); i++) {
//            coords.omega[i] = w[i];
//        }
//	V = new V_UCHO(w);
//	vector<double> lambda;
//	for(int j=0; j<N; j++){
//		lambda.push_back(5);
//	}
//	V = new V_AQO(w,lambda);
//	V = new V_Morse(12.5,2,N);
//	V = new V_Morse(0.033,2.0,N);		//Corresponds to reasonable molecular values (actually too stiff a spring!)
//	V = new V_Morse(0.00298,2.0,N);		//Corresponds to reasonable molecular values
//	V = new V_Morse(0.176, 1.02, N);	//Corresponds to H2 morse potential from http://www1.uprh.edu/rbaretti/MorsePotential15mar2011.htm but 'a' corresponds to 3000 cm-1
//	V = new V_Morse(0.176, 1.4886, N);	//Corresponds to H2 with omega 4522.42
//	V = new V_TinkerExecutable(coords);
//	rho = new Rho_Free();
//	rho = new Rho_HO(w);
//	rho = new Rho_TruncHO(w);

	if(!physics.vType.compare("V_TinkerExecutable")) {
		V = new V_TinkerExecutable(coords);
	}
	else if(!physics.vType.compare("V_Tinker")) {
//		cout << "Unimplemented type" << endl;
//		exit(-1);
LINE
		V = new V_Tinker(coords, m, eps, beta);
	}
	else if(!physics.vType.compare("V_Morse")) {
		if(physics.morseAlpha != 0.0 && physics.morseDE != 0.0) {
			V = new V_Morse(physics.morseDE, physics.morseAlpha, N);
		}
		else {
			V = new V_Morse(0.176, 1.4886, N);		//H2
		}
	}
	else if(!physics.vType.compare("V_UCHO")) {
		V = new V_UCHO(coords.omega);
	}

	if(!physics.rhoType.compare("Rho_HO")) {
		rho = new Rho_HO(coords.omega);
	}
	else if(!physics.rhoType.compare("Rho_Free")) {
		rho = new Rho_Free();
	}
	else if(!physics.rhoType.compare("Rho_TruncHO")) {
		rho = new Rho_TruncHO(coords.omega);
	}
	else {
		cout << "Invalid rhoType:\t" << physics.rhoType << endl;
		exit(-1);
	}


	for(int i=0; i<P; i++){
		for(int j=0; j<N; j++){
			part[i][j].mass = m[j];
//			part[i][j].atomType = aType[j];
			oldPart[i][j].mass = m[j];
//			oldPart[i][j].atomType = aType[j];
		}
	}

	vector<double> r0;

//*******************************************************
//*	Initialize Bead Polymer Positions with Random Walks *
//*******************************************************

//	Initialize all beads to center
	for(int j=0; j<N; j++) {
		for(int i=0; i<P; i++) {
			for(int k=0; k<coorDim; k++) {
				part[i][j].pos.push_back(0.0);
				oldPart[i][j].pos.push_back(0.0);
			}
		}
	}

//	bool randomWalks = true;
	bool randomWalks = false;
	static int initRanSeed = -time(0);
	int bead, beadm1;
	vector <double> levyMean;
	double levySigma;
	if(randomWalks) {
		for(int j=0; j<N; j++) {
			for(int s=0; s<P; s++) {
				bead = (s+1)%P;
				beadm1 = (s)%P;
				levyMean = rho->GetLevyMean(part[beadm1][j], part[P-1][j], (double)(P-s), eps, j);
				levySigma = rho->GetLevySigma((double)(P-s), eps, j) / sqrt(part[0][j].mass);
//				cout << "levySigma = " << levySigma << endl;
				for(int k=0; k<coorDim; k++){
					part[bead][j].pos[k] = RandomNum::rangau(levyMean[k], levySigma, &initRanSeed);
					oldPart[bead][j].pos[k] = part[bead][j].pos[k];
				}
				upToDate[s] = false;
			}
		}
	}
	else {
		int levyInit = physics.numInit;
		int snipBegin;
		int snipEnd;
		for (int j=0; j<N; j++) {
			for(int round=0; round<P/(levyInit+1); round++) {
				snipBegin = round*(levyInit+1);
//				cout << "snipBegin = " << snipBegin << endl;
				snipEnd = (snipBegin+levyInit+1)%P;
//				cout << "snipEnd = " << snipEnd << endl;
//				cout << "round = " << round << "\t j = " << j << endl;
				for (int s=0; s<levyInit; s++) {
					bead = (snipBegin+s+1)%P;
					beadm1 = (snipBegin+s)%P;
					levyMean = rho->GetLevyMean(part[beadm1][j], part[snipEnd][j], (double)(levyInit-s), eps, j);
					levySigma = rho->GetLevySigma((double)(levyInit-s), eps, j) / sqrt(part[0][j].mass);
					for(int k=0; k<coorDim; k++){
						part[bead][j].pos[k] = RandomNum::rangau(levyMean[k], levySigma, &initRanSeed);
					}
					upToDate[bead] = false;
				}
			}
			int levyLength = P%(levyInit+1) - 2;
			snipBegin = P-P%(levyInit+1);
			snipEnd = (snipBegin+levyLength+1);
			for (int s=0; s<levyLength; s++) {
				bead = (snipBegin+s+1)%P;
				beadm1 = (snipBegin+s)%P;
				levyMean = rho->GetLevyMean(part[beadm1][j], part[snipEnd][j], (double)(levyLength-s), eps, j);
				levySigma = rho->GetLevySigma((double)(levyLength-s), eps, j) / sqrt(part[0][j].mass);
				for(int k=0; k<coorDim; k++){
					part[bead][j].pos[k] = RandomNum::rangau(levyMean[k], levySigma, &initRanSeed);
				}
				upToDate[bead] = false;
			}
		}
	}


	oldEnergy = 1000000;		// TODO: clean this up
	oldPotE = 100000;
	ECheckFlag = false;
}

TestSystem::~TestSystem() {
	// TODO Auto-generated destructor stub
}

double TestSystem::GetWeight() {
	if(ECheckFlag) {
		cout << "I shouldn't be here!" << endl;
		CalcEnergy();
		return energy;
	}
	else {
		CalcPotential();
		return potE;
	}
}

double TestSystem::GetOldWeight() {
//	TODO: Fix this so it doesn't break if using different ECheckFlag on different steps
	if(ECheckFlag) {
		cout << "I shouldn't be here!" << endl;
		return oldEnergy;
	}
	else {
		return oldPotE;
	}
}

double TestSystem::EstimatorV() {
	CalcPotential();			//Might not need this
	return potE/((double)P);
//	return exp(double(-eps)*(potE));
}

double TestSystem::EstimatorE() {
//	avgV = (avgV*double(numSteps) + potE)/double(numSteps+1);
//	numSteps++;

	double est=0;
	for(int i=0; i<P; i++) {
		est += rho->Estimate(part[i], part[(i+1)%P], eps, P);
	}
	est += potE/(double)P;
//	cout << 1.0 - 1.0/(1.0+double(eps)*(avgV - potE)) << endl;
//	est *= (1.0 - 1.0/(1.0+double(eps)*(avgV - potE)));
//	est *= (1.0 - (1.0+double(eps)*potE)/(1.0+double(eps)*avgV));
//	est += potE/double(P);
//	cout << est << endl;
	return est;
}

vector < vector<Particle> > TestSystem::GetParticle() {
	return part;
}

void TestSystem::CalcEnergy(){
//	TODO: Consider storing V and only updating those that moved!!
	cout << "I shouldn't be here!" << endl;
	energy = 0;
	for(int i=0; i<P; i++){
		if(!upToDate[i]) {
			sliceV[i] = V->GetV(part[i], rho);
		}
		upToDate[i] = true;
		energy += sliceV[i];
	}
	potE = energy;
	for(int i=0; i<P; i++){
		energy += rho->GetRho(part[i], part[(i+1)%P], eps);
	}
//	cout << "Energy = " << energy << endl;
}

void TestSystem::CalcPotential(){
	potE = 0;
	for(int i=0; i<P; i++){
		if(!upToDate[i]) {
			sliceV[i] = V->GetV(part[i], rho);
		}
		upToDate[i] = true;
		potE += sliceV[i];
	}
}


void TestSystem::Move(vector<double> prob, int levyNum){
//	cout << "Take a step" << endl;
/*
//	Random Walks
	ECheckFlag = true;
	double stepSize = 1.0;
	int iPart = static_cast<int>(N*P*prob[1]);
	int iCoor = static_cast<int>(coorDim*prob[2]);
	part[(iPart-iPart%N)/N][(iPart%N)].pos[iCoor] += (stepSize*(prob[0]-0.5));		// For number (0,1)
//	part[(iPart-iPart%N)/N][(iPart%N)].pos[iCoor] += (stepSize*(prob[0]/2));		// For symm random number
*/
///*
	bool debugging = false;
	if(debugging) {
//		Debug movement
		ECheckFlag = false;
		for(int i=0; i<P; i++) {
			part[i][1].pos[0] -= 1.0;
			upToDate[i] = true;
		}
	}
	else {
//	Levy-Flight move
		if(P<3) {
			cout << "Warning! You should not be using Levy-Flight with P < 3!" << endl;
		}
		ECheckFlag = false;		// Check only against potential energy!
		static int seed = -time(0);
		vector <double> levyMean;
		double levySigma;
		int bead, beadm1;
//	Snip length must be at least 3 and shouldn't be more than P/2+1 (or 20) or else you'll just collapse it
//	snipLength-2 is the number of beads moved.
//		int maxLevy = 30;
//		int levyLength = min(static_cast<int>(prob[0]*(P-2)/2+1), maxLevy);
//	*******************************************************
		levyLength = levyNum;
//	*******************************************************
		int iPart = static_cast<int>(N*P*prob[1]);
//		int iCoor = static_cast<int>(coorDim*prob[2]);
		int pickedPart = int(iPart/P);
		int snipBegin = iPart%P;
		int snipEnd = (snipBegin+levyLength+1)%P;

//		cout << "SnipBegin: " << snipBegin << "\tSnipEnd" << snipEnd << "\tSnipLength:" << snipLength << endl;

		for (int s=0; s<levyLength; s++) {
			bead = (snipBegin+s+1)%P;
			beadm1 = (snipBegin+s)%P;
			levyMean = rho->GetLevyMean(part[beadm1][pickedPart], part[snipEnd][pickedPart], (double)(levyLength-s), eps, pickedPart);
			levySigma = rho->GetLevySigma((double)(levyLength-s), eps, pickedPart) / sqrt(part[0][pickedPart].mass);
			for(int k=0; k<coorDim; k++){
				part[bead][pickedPart].pos[k] = RandomNum::rangau(levyMean[k], levySigma, &seed);
			}
			upToDate[bead] = false;
		}
	}
	//*/
}

void TestSystem::Forget() {
	oldEnergy = energy;
	oldPotE = potE;
	for(int i=0; i<P; i++) {
		for(int j=0; j<N; j++) {
			oldPart[i][j].pos = part[i][j].pos;
		}
		oldSliceV[i] = sliceV[i];
	}
}

void TestSystem::Undo() {
	energy = oldEnergy;
	potE = oldPotE;
	for(int i=0; i<P; i++) {
		for(int j=0; j<N; j++) {
			part[i][j].pos = oldPart[i][j].pos;
		}
		sliceV[i] = oldSliceV[i];
	}
}


double TestSystem::Debug() {
	return potE;
}

string TestSystem::GetVType() {
	return V->GetType();
}

string TestSystem::GetRhoType() {
	return rho->GetType();
}
