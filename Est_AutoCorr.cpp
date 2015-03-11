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
   //delete hist;
   //hist = new deque<double>();
   hist.clear();
   nSamples = 0;
   
}

void Est_AutoCorr::AddVal(double newVal) {
//Add val to history
   hist.push_front(newVal);
   if(hist.size()>maxHist) {
      hist.pop_back();
   }
   nSamples++;

   totalStats.AddVal(newVal);
   double newMean = totalStats.GetMean();
   newVal -= newMean;
   int count=0;
   for(deque<double>::iterator iter=hist.begin(); iter!=hist.end(); iter++) {
      corrStats[count].AddVal((*iter-newMean) * newVal);
      count++;
   }
}

vector<double> Est_AutoCorr::GetCorr() {
   if(nSamples<maxHist) {
      cout << "DES Warning: number of samples less than max history in Est_AutoCorr" << endl;
      cout << nSamples << " < " << maxHist << endl;
   }
   vector<double> corrFunc (maxHist, 0.0);
   //double var = totalStats.GetVariance();
   double var = corrStats[0].GetMean()/double(nSamples);
   //var /= min(maxHist,nSamples);
//   cout << "totalStats.GetVariance() =" << totalStats.GetVariance() << endl;
//   cout << "corrStats[0].GetVariance() = " << corrStats[0].GetVariance() << endl;
//   cout << "var = " << var << endl;
//   cout << "corrStats[0] mean = " << corrStats[0].GetMean() << endl;;
   for(int i=0; i<min(maxHist,nSamples); i++) {
      corrFunc[i] = corrStats[i].GetMean()/var/double(nSamples-i);
   }
   return corrFunc;
}

double Est_AutoCorr::GetTau() {
//   if(nSamples<maxHist) {
//      cout << "DES Warning: number of samples less than max history in Est_AutoCorr" << endl;
//      cout << nSamples << " < " << maxHist << endl;
//   }

// Compute corrFunction
   arma::vec logCorr(GetCorr());
   logCorr.print("logCorr");
   cout << "nSamples = " << nSamples << endl;
   //find(logCorr < 0.1, 1).print("first < 0.1")
   arma::uvec firstIncrease = find(logCorr < 0.1, 1);
   int stopIndex;
   if(firstIncrease.size() == 1) {
      stopIndex = arma::as_scalar(find(logCorr < 0.1, 1));
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
      found = logCorr(i-1) < logCorr(i);
   }
   cout << "first turn up at " << i << endl;
   if(i < stopIndex) {
      if(logCorr(i) < 0.3) {
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
   logCorr = logCorr.subvec(0,stopIndex);
   int nHist = logCorr.size();
   logCorr = log(abs(logCorr));
   arma::vec tStep = logCorr;
   for(int i=0; i<nHist; i++) {
      tStep(i) = double(i);
   }
   tStep.print("index");
   logCorr.print("logCorr");
   double tau = -1.0/arma::as_scalar(tStep.t()*logCorr)*double(nHist*(nHist-1)*(2*nHist-1))/6.0;
   cout << "tau = " <<  tau << endl;
   
   exit(-1);

   return tau;
}
