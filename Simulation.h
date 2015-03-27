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
#include "Est_AutoCorr.h"
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
      void WritePosToFile();
      void Log();
      void FinalLog();
      void TILog(double,double);
      void Store();
      void Tokenize(const string&, vector<string>&, const string& = " ");
      void GetGaussianQuad(int, vector<double>&, vector<double>&);
      void GetLobattoQuad(int, vector<double>&, vector<double>&);
      void GetLinearQuad(int, vector<double>&, vector<double>&);

      int initLen;
      int autoCorrLen;     //if set, overwrites sampleStart, numSamples, maxStep
      int maxSim;        //Old Var
      int stepNum;
      int maxStep;
      int sampleStart;
      int numSamples;      //temp variable to overwrite maxStep
      int sampleFreq;
      int convFreq;        //TODO: Remove
      int storeFreq;
      int storeNum;        //temp variable to overwrite storeFreq
      int levyNum;
      int levyModes;
      int numTI;
      int * idum;
//		int * idum2;
      double tau;
      double errorThresh;
      double beta;			//TODO: Remove
      double epsTemp;
      double stepSize;
      std::string posFileName;
      std::string outFileName;
      ofstream posFile;
      ofstream logFile;
      Est_AutoCorr autoCorr;
      Stats * simStats;
      Stats * simPotStats;
      Stats * energyStats;
      Stats * potentialStats;
      Stats * convergenceStats;
      Stats * acceptanceStats;
      Stats * xStats;		//TODO: Remove
      System * sys;

//DES Temp:
//                ofstream vFile;      //TODO: Remove
//                int tempNum;      //TODO: Remove
};

#endif /* SIMULATION_H_ */
