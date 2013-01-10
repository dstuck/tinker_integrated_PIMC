/*
 * Stats.h
 *
 *  Created on: May 9, 2012
 *      Author: dstuck
 */

#ifndef STATS_H_
#define STATS_H_

#include<cmath>
#include<vector>
using namespace std;

class Stats {
public:
	Stats();
//	Stats(bool);		// Could implement if I'd like this class to store data in vals as well...
	virtual ~Stats();
	void AddVal(double);
	double GetMean();
	double GetMeanAbs();
	double GetRMS();
	double GetVariance();
	double GetStDev();
	double GetMax();
	double GetMin();
        void Reset();

//	bool store;
	int numVals;
	double sum;
	double normSum;
	double squareSum;
	double maxVal;
	double minVal;
//	double maxSigned;
//	double minSigned;
//	vector<double> vals;
};

#endif /* STATS_H_ */
