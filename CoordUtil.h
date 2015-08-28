/*
 * CoordUtil.h
 *
 *  Created on: Sep 16, 2012
 *      Author: dstuck
 */

#ifndef COORDUTIL_H_
#define COORDUTIL_H_

#include "debug.h"
#include "Particle.h"
#include <math.h>
#include <iostream>
#include <vector>
#include <map>
#include "armadillo"
#include "optking/molecule.h"
using namespace std;

class CoordUtil {
private:
	vector< vector<double> > LinearStep(arma::vec,vector< vector<double> >, bool equib = false);

public:
	CoordUtil();
	CoordUtil(int, int, bool, vector< vector< vector<double> > >, vector<double>, vector<double>, vector< vector<double> >, vector<string>, vector<int>, vector< vector<int> >, string, string,bool,bool);
	virtual ~CoordUtil();

        void initInternals();
	vector< vector<double> > normalModeToCart(vector<Particle>);
	vector< vector<double> > normalModeToCart(vector<double>);      //Inefficient, for debug
	vector<double> cartToNormalMode(vector< vector<double> >);
        double GetHarmV(vector<Particle>);      //TODO: Make this vector position...
        void MatchModes(CoordUtil*);
        void MakeScalingVec(CoordUtil*);
        void initMass();
        int atomStrToInt(string);
        CoordUtil* Clone();

      bool useGuess;
        bool internals;
        bool readOmega;
        bool readGeom;
        bool scaleGeom;
	int numModes;
	int numPart;
	int dim;
	std::string outFileName;
	std::string prmName;
        map<string,double> atomToMass;
	vector<string> atomType;
	vector<int> paramType;
//        vector<double> mass;
        vector<double> reducedMass;
	vector<double> omega;
	vector< vector<int> > connectivity;
	vector< vector<double> > initCart;
	vector< vector<double> > guessCarts;
	vector<Particle> guessPart;
	vector< vector< vector<double> > > normModes;    //TODO: remove this in favor of armaNormModes
      arma::vec scalingVec;
      arma::vec totalCOM;
      arma::vec atomMass;
      arma::vec armaInitCart;
      arma::mat fragExtraB;
      arma::mat armaNormModes;
      arma::mat armaInvModes;
      arma::mat sqrtInvMass;  //TODO: Don't really need this as class data
      //arma::mat equibBInv;
      arma::mat equibBMat;
      arma::mat delocIntCos;
      //arma::mat wMat;         // Weight matrix to normalize equibBInv, not sure if this is right...
      opt::MOLECULE * workingMol;
};

#endif /* COORDUTIL_H_ */
