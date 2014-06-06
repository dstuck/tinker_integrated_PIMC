/*
 * Simulation.h
 *
 *  Created on: May 4, 2012
 *      Author: dstuck
 */

#ifndef SIMULATION_H_
#define SIMULATION_H_
#include <iostream>
#include <fstream>
#include <math.h>
#include "debug.h"
#include <time.h>
#include <string>
#include "System.h"
#include "Stats.h"
#include "TestSystem.h"
#include "Random.h"
#include "CoordUtil.h"
#include "PhysicsUtil.h"
using namespace std;

class Simulation {
	public:
		Simulation(string, string, string);
		virtual ~Simulation();
		void TakeStep();
		bool Check();
		void Update();
		void Revert();
		void Sample();
		void Run();
		void FinalPrint();     //TODO: Remove
		void WritePosToFile();
		void Log();
		void FinalLog();
		void TILog(double,double);
		void Store();
		void Tokenize(const string&, vector<string>&, const string& = " ");
        void GetGaussianQuad(int, vector<double>&, vector<double>&);
        void GetLobattoQuad(int, vector<double>&, vector<double>&);
        void GetLinearQuad(int, vector<double>&, vector<double>&);

        int maxSim;
		int stepNum;      //TODO: Remove
		int maxStep;
		int sampleStart;
		int sampleFreq;
		int convFreq;
		int storeFreq;
		int levyNum;
        int numTI;
		int * idum;
//		int * idum2;
		double beta;			//TODO: Remove
		double epsTemp;
		double stepSize;
		std::string posFileName;
		std::string outFileName;
		ofstream posFile;
		ofstream logFile;
		Stats * simStats;
		Stats * simPotStats;
//		Stats * simComboStats;      //TODO: Remove
		Stats * energyStats;
		Stats * potentialStats;
//		Stats * comboStats;
		Stats * convergenceStats;
		Stats * acceptanceStats;
//		Stats * xStats;		//TODO: Remove
		System * sys;
};

#endif /* SIMULATION_H_ */
