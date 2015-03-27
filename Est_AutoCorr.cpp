/*
 * Est_AutoCorr.cpp
 *
 *  Created on: March 3, 2015
 *      Author: dstuck
 */

#include "Est_AutoCorr.h"
#include "armadillo"

Est_AutoCorr::Est_AutoCorr() {
}

Est_AutoCorr::Est_AutoCorr(int maxH) : maxHist(maxH) {
   //hist = new deque<double>();
   //totalStats = Stats();
   Stats tempStats = Stats();
   for(int i=0; i<maxH; i++) {
      corrStats.push_back(tempStats);
   }
   nSamples=0;
}

Est_AutoCorr::~Est_AutoCorr() {
   //delete hist;
}

void Est_AutoCorr::Reset() {
   for(int i=0;i<corrStats.size();i++) {
      corrStats[i].Reset();
   }
   totalStats.Reset();
   //delete hist;
   //hist = new deque<double>();
   hist.clear();
   nSamples = 0;
   
}

void Est_AutoCorr::AddVal(double newVal) {
// corrStats is maintaining summ x_0*x_t for t=0::1000
//Add val to history
   hist.push_front(newVal);
   if(hist.size()>maxHist) {
      hist.pop_back();
   }
   nSamples++;

   totalStats.AddVal(newVal);
   //double newMean = totalStats.GetMean();
   //newVal -= newMean;
   int count=0;
   for(deque<double>::iterator iter=hist.begin(); iter!=hist.end(); iter++) {
      //corrStats[count].AddVal((*iter-newMean) * newVal);
      corrStats[count].AddVal((*iter) * newVal);
      count++;
   }
}

vector<double> Est_AutoCorr::GetCorr() {
   if(nSamples<maxHist) {
      cout << "DES Warning: number of samples less than max history in Est_AutoCorr" << endl;
      cout << nSamples << " < " << maxHist << endl;
   }
   vector<double> corrFunc(maxHist, 0.0);
   //double var = totalStats.GetVariance();
   double meanSquared = pow(totalStats.GetMean(),2);
   double var = (corrStats[0].GetMean() - meanSquared)/double(nSamples);
//   cout << "total stats mean = " << totalStats.GetMean() << endl;
//   cout << "total stats stdev = " << totalStats.GetStDev() << endl;
   //var /= min(maxHist,nSamples);
//   cout << "totalStats.GetVariance() =" << totalStats.GetVariance() << endl;
//   cout << "corrStats[0].GetVariance() = " << corrStats[0].GetVariance() << endl;
//   cout << "var = " << var << endl;
//   cout << "corrStats[0] mean = " << corrStats[0].GetMean() << endl;;
   for(int i=0; i<min(maxHist,nSamples); i++) {
      corrFunc[i] = (corrStats[i].GetMean() - meanSquared)/totalStats.GetVariance()/double(nSamples-i)*double(nSamples);
   }
   return corrFunc;
}

double Est_AutoCorr::GetTau() {
//   if(nSamples<maxHist) {
//      cout << "DES Warning: number of samples less than max history in Est_AutoCorr" << endl;
//      cout << nSamples << " < " << maxHist << endl;
//   }

// Compute corrFunction
   arma::vec corrFunc(GetCorr());
   ofstream outLog;
   outLog.open("autoCorr.txt");
   corrFunc.print(outLog);
   //corrFunc.print("autoCorr");
   //cout << "nSamples = " << nSamples << endl;
   //find(corrFunc < 0.1, 1).print("first < 0.1")

/* Old linear regression approach
   arma::uvec firstIncrease = find(corrFunc < 0.1, 1);
   int stopIndex;
   if(firstIncrease.size() == 1) {
      stopIndex = arma::as_scalar(find(corrFunc < 0.1, 1));
   }
   else {
      cout << "Warning: AutoCorr never increases in Est_AutoCorr!" << endl;
      stopIndex = min(nSamples,maxHist);
   }
   cout << "Below 0.1 at " << stopIndex << endl;
   bool found=false;
   int i = 0;
   while(!found && i<corrStats.size()-1) {
      i++;
      found = corrFunc(i-1) < corrFunc(i);
   }
   cout << "first turn up at " << i << endl;
   if(i < stopIndex) {
      if(corrFunc(i) < 0.3) {
         stopIndex = i;
      }
      else {
         cout << "Warning: Autocorrelation function turns up while above 0.3 in Est_AutoCorr.cpp" << endl;
         stopIndex = (stopIndex+i)/2;
      }
   }
 
   //stopIndex += i;
   //stopIndex /= 2;
   cout << "Final stop index = " << stopIndex << endl;
   corrFunc = corrFunc.subvec(0,stopIndex);
   int nHist = corrFunc.size();
   corrFunc = log(abs(corrFunc));
   arma::vec tStep = corrFunc;
   for(int i=0; i<nHist; i++) {
      tStep(i) = double(i);
   }
   //tStep.print("index");
   //corrFunc.print("logGorr");
   double tau = -1.0/arma::as_scalar(tStep.t()*corrFunc)*double(nHist*(nHist-1)*(2*nHist-1))/6.0;
*/
   int maxCheck=500;
   arma::uvec firstZero = find(corrFunc < 0.0, 1);
   if(firstZero.size() == 1) {
      maxCheck = min(maxCheck,int(arma::as_scalar(firstZero)));
      //cout << "First zero at " << arma::as_scalar(firstZero) << endl;
   }
   double tau = arma::sum(corrFunc.subvec(0,maxCheck));
   //cout << "first tau = " <<  tau << endl;
   if(maxCheck > 5*tau) {
      //cout << "Recalculating" << endl;
      tau = arma::sum(corrFunc.subvec(0,int(tau*5)));
   }
   
   //cout << "tau = " <<  tau << endl;
   
   //exit(-1);

   return tau + 0.5;
}

double Est_AutoCorr::GetTotalMean() {
   return totalStats.GetMean();
}

double Est_AutoCorr::GetTotalVariance() {
   return totalStats.GetVariance();
}

double Est_AutoCorr::GetTotalStDev() {
   return totalStats.GetStDev();
}
