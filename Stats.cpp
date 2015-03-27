/*
 * Stats.cpp
 *
 *  Created on: May 9, 2012
 *      Author: dstuck
 */

#include "Stats.h"

Stats::Stats() {
	numVals = 0;
	sum = 0;
	normSum = 0;
	squareSum = 0;
	maxVal = 0;
	minVal = 0;
}

Stats::~Stats() {
	// TODO Auto-generated destructor stub
}

void Stats::AddVal(double newVal) {
	numVals++;
	sum += newVal;
	normSum += abs(newVal);
	squareSum += newVal*newVal;
	maxVal = max(maxVal, abs(newVal));
	minVal = min(minVal, abs(newVal));
}

double Stats::GetMean() {
	if(numVals>0) {
		return sum/(double)numVals;
	}
	else {
		return 0.0;
	}
}

double Stats::GetMeanAbs() {
	if(numVals>0) {
		return normSum/(double)numVals;
	}
	else {
		return 0.0;
	}
}

double Stats::GetRMS() {
	if(numVals>0) {
		return sqrt(squareSum/(double)numVals);
	}
	else {
		return 0.0;
	}
}

double Stats::GetVariance() {
	if(numVals>0) {
		return squareSum/(double)numVals - (sum/(double)numVals)*(sum/(double)numVals);
		//return (squareSum/(double)numVals - (sum/(double)numVals)*(sum/(double)numVals))/(double)numVals;
	}
	else {
		return 0.0;
	}
}

double Stats::GetStDev() {
	if(numVals>0) {
		return sqrt((squareSum/(double)numVals - (sum/(double)numVals)*(sum/(double)numVals)));
//		return sqrt((squareSum/(double)numVals - (sum/(double)numVals)*(sum/(double)numVals))/(double)numVals);
	}
	else {
		return 0.0;
	}
}

double Stats::GetMax() {
	return maxVal;
}

double Stats::GetMin() {
	return minVal;
}

void Stats::Reset() {
	numVals = 0;
	sum = 0;
	normSum = 0;
	squareSum = 0;
	maxVal = 0;
	minVal = 0;
}
