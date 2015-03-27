/*
 * Est_AutoCorr.h
 *
 *  Created on: March 3, 2015
 *      Author: dstuck
 */

#ifndef EST_AUTOCORR_H_
#define EST_AUTOCORR_H_

#include <cmath>
#include <deque>
#include <vector>
#include <fstream>
#include "armadillo"
#include "Stats.h"

using namespace std;

class Est_AutoCorr {
   public:
      Est_AutoCorr();
      Est_AutoCorr(int);
      virtual ~Est_AutoCorr();
      void Reset();
      void AddVal(double);
      vector<double> GetCorr();
      double GetTau();
      double GetTotalMean();
      double GetTotalVariance();
      double GetTotalStDev();

      int maxHist;
      int nSamples;
      deque<double> hist;
      Stats totalStats;
      vector<Stats> corrStats;
};

#endif /* EST_AUTOCORR_H_ */
