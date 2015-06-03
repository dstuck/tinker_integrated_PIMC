/*
 * Simulation.cpp
 *
 *  Created on: May 4, 2012
 *      Author: dstuck
 */

#include "Simulation.h"

Simulation::Simulation(string inFileName, string logFileName, string prmFile) : autoCorr(1000) {

//*******************************************************
//*		  Default Parameters			*
//*******************************************************
/*
	stepNum = 0;
	maxStep = 1000000;
	sampleStart = 1000;
	sampleFreq = 5;
	convFreq = 10000;
	storeFreq = 100000;
//	beta = 21000;	//21000 corresponds to 15 K
	beta = 10500;	//30 K
//	beta = 5250;	//60 K
 */
// DES Temp:
//        tempNum = 0;
        maxSim = 1;
        epsInit = -1.0;
        bool readOmega = false;
        bool readGeom = false;
        bool internalCoords = false;
        bool scaleGeom = false;
	int nPart = 0;
	int pSlice = 0;
	int nMode;
	int coorDim;
        bool deltaAI;
//	I had been running this previously with the wrong intended mass (3500 instead of 35000)
//	double m = 3500.0;				// 35000 me corresponds to ~20 nuclei
//	The line below is false, w = 0.00868 corresponds to 11976 wavenumber
//	double w = 2.0*sqrt(2.0*0.033/3500.0);		//.0026 corresponds to 3600 wavenumber
	vector<double> kFreq;
	vector<double> omega;
//	for(int j=0; j<nPart; j++){
//		omega.push_back(w);
//	}
//	int a = 1;	// Atomic number
	levyNum = 30;
        levyModes = 0;
	levyNumInit = -1;
        levyModesInit = -1;
        numTI = 0;
        numGeomPrint = 0;
        errorThresh = 0.05;
        tau = -1.0;
        storeNum = -1;
        maxStep = -1;
        sampleStart = -1;
        initLen = -1;
        autoCorrLen = 20000;
        //autoCorrLen = -1;
        sampleFreq = 1;
        noEnergy = false;
	vector<double> mass;
	vector<string> atomType;
	vector<int> paramType;
	vector< vector<int> > connectivity;
	vector< vector<double> > initPos;
	vector< vector< vector<double> > > modes;				//modes[i][j][k] will the the kth coordinate of the jth atom of the ith normal mode
	PhysicsUtil * physicsParams = new PhysicsUtil();

//*******************************************************
//*				   Reading in input					    *
//*******************************************************

//	Open the infile
	ifstream input;
	input.open(inFileName.c_str());
// DES Temp:
//        string tempstr = "potential.txt";
//        vFile.open(tempstr.c_str());
//	Read in and parse the file
	string lineRead;
	while(std::getline(input, lineRead)) {		// Reads next line until it reaches end of file
//		Read molecule section from .pin file
		if(lineRead.find("$molecule") != std::string::npos) {		//If line contains $molecule
			while(!(std::getline(input, lineRead).eof()) && (lineRead.find("$end") == std::string::npos)) {			//Reads next line always, goes through parsing if line isn't $end
//				cout << lineRead << endl;
				vector<string> lineTokens;
				Tokenize(lineRead, lineTokens);
				int tokenShift = 0;
//				Ignore first token if it is just an integer
				if(atoi(lineTokens[0].c_str())) {
					tokenShift++;
				}
				if(int(lineTokens.size()) < 4+tokenShift) {
					cout << "Error in input file: $molecule" << endl;
					cout << "Found line: " << lineRead << endl;
					exit(-1);
				}
				atomType.push_back(lineTokens[tokenShift]);
				tokenShift++;
				vector<double> p0;
				for(int k=0; k< 3; k++) {
					p0.push_back(atof(lineTokens[k+tokenShift].c_str()));
				}
				tokenShift += 3;
				initPos.push_back(p0);
				vector<int> con0;
				for(int c=0; c<int(lineTokens.size()-tokenShift); c++) {
					con0.push_back(atoi(lineTokens[c+tokenShift].c_str()));
				}
				connectivity.push_back(con0);
			}
		} // end of reading molecule section

//		Read in parameters from .pin file
		if(lineRead.find("$params") != std::string::npos) {
			while(!(std::getline(input, lineRead).eof()) && (lineRead.find("$end") == std::string::npos)) {
//				cout << lineRead << endl;
				if(lineRead.find("!") != std::string::npos) {
				}
				else if(lineRead.find_first_not_of(" \t\r\n") == std::string::npos) {
				}
				else {
					vector<string> lineTokens;
					Tokenize(lineRead, lineTokens);
					if(lineTokens.size() != 2) {
						cout << "Error in input file: $params" << endl;
						cout << "Found line: " << lineRead << endl;
						exit(-1);
					}
					if(lineTokens[0].find("maxSim") != std::string::npos) {
						maxSim = atoi(lineTokens[1].c_str());
					}
					else if(lineTokens[0].find("maxStep") != std::string::npos) {
						maxStep = atoi(lineTokens[1].c_str());
					}
					else if(lineTokens[0].find("initLen") != std::string::npos) {
						initLen = atoi(lineTokens[1].c_str());
					}
					else if(lineTokens[0].find("autoCorrLen") != std::string::npos) {
						autoCorrLen = atoi(lineTokens[1].c_str());
					}
					else if(lineTokens[0].find("sampleStart") != std::string::npos) {
						sampleStart = atoi(lineTokens[1].c_str());
					}
					else if(lineTokens[0].find("sampleFreq") != std::string::npos) {
						sampleFreq = atoi(lineTokens[1].c_str());
					}
					else if(lineTokens[0].find("convFreq") != std::string::npos) {
						convFreq = atoi(lineTokens[1].c_str());
					}
					else if(lineTokens[0].find("storeFreq") != std::string::npos) {
						storeFreq = atoi(lineTokens[1].c_str());
					}
					else if(lineTokens[0].find("storeNum") != std::string::npos) {
						storeNum = atoi(lineTokens[1].c_str());
					}
					else if(lineTokens[0].find("errorThresh") != std::string::npos) {
						errorThresh = atof(lineTokens[1].c_str())/100.0;
					}
					else if(lineTokens[0].find("beta") != std::string::npos) {
						beta = atof(lineTokens[1].c_str());
					}
					else if(lineTokens[0].find("temperature") != std::string::npos) {
// Atomic temperature units from http://en.wikipedia.org/wiki/Atomic_units 4/18/15
						beta = 315774.65/atof(lineTokens[1].c_str());
					}
					else if(lineTokens[0].find("tempInit") != std::string::npos) {
						epsInit = 315774.65/atof(lineTokens[1].c_str());
					}
					else if(lineTokens[0].find("initTemp") != std::string::npos) {
						epsInit = 315774.65/atof(lineTokens[1].c_str());
					}
	//				TODO: merge nPart read with $molecule section
					else if(lineTokens[0].find("nPart") != std::string::npos) {
						nPart = atoi(lineTokens[1].c_str());
					}
					else if(lineTokens[0].find("pSlice") != std::string::npos) {
						pSlice = atoi(lineTokens[1].c_str());
					}
					else if(lineTokens[0].find("levyNum") != std::string::npos) {
                                           if(lineTokens[0].find("levyNumInit") != std::string::npos) {
                                                levyNumInit = atoi(lineTokens[1].c_str());
                                           }
                                           else {
						levyNum = atoi(lineTokens[1].c_str());
                                           }
					}
					else if(lineTokens[0].find("levyModes") != std::string::npos) {
                                           if(lineTokens[0].find("levyModesInit") != std::string::npos) {
                                                levyModesInit = atoi(lineTokens[1].c_str());
                                           }
                                           else {
						levyModes = atoi(lineTokens[1].c_str());
                                           }
					}
					else if(lineTokens[0].find("levyInit") != std::string::npos) {
						physicsParams->numInit = atoi(lineTokens[1].c_str());
					}
					else if(lineTokens[0].find("rhoType") != std::string::npos) {
						physicsParams->rhoType = lineTokens[1];
					}
					else if(lineTokens[0].find("vType") != std::string::npos) {
						physicsParams->vType = lineTokens[1];
					}
					else if(lineTokens[0].find("morseDE") != std::string::npos) {
						physicsParams->morseDE = atof(lineTokens[1].c_str());
					}
					else if(lineTokens[0].find("morseAlpha") != std::string::npos) {
						physicsParams->morseAlpha = atof(lineTokens[1].c_str());
					}
					else if(lineTokens[0].find("morseMass") != std::string::npos) {
						physicsParams->morseMass = atof(lineTokens[1].c_str());
					}
					else if(lineTokens[0].find("freqCutoff") != std::string::npos) {
						physicsParams->freqCutoff = atof(lineTokens[1].c_str())*0.00000455633;
					}
					else if(lineTokens[0].find("numFrozModes") != std::string::npos) {
						physicsParams->lowFrozModes = atof(lineTokens[1].c_str());
					}
					else if(lineTokens[0].find("lowFrozModes") != std::string::npos) {
						physicsParams->lowFrozModes = atof(lineTokens[1].c_str());
					}
					else if(lineTokens[0].find("highFrozModes") != std::string::npos) {
						physicsParams->highFrozModes = atof(lineTokens[1].c_str());
					}
					else if(lineTokens[0].find("numTI") != std::string::npos) {
					        numTI=atoi(lineTokens[1].c_str());
					}
					else if(lineTokens[0].find("numGeomPrint") != std::string::npos) {
					        numGeomPrint=atoi(lineTokens[1].c_str());
					}
					else if(lineTokens[0].find("charge") != std::string::npos) {
						physicsParams->charge = atoi(lineTokens[1].c_str());
					}
					else if(lineTokens[0].find("multiplicity") != std::string::npos || lineTokens[0].find("spin") != std::string::npos) {
						physicsParams->multiplicity = atoi(lineTokens[1].c_str());
					}
					else if(lineTokens[0].find("deltaAI") != std::string::npos) {
                                            if(atof(lineTokens[1].c_str())>0) {
                                                physicsParams->deltaAbInit = true;
                                            }
					}
					else if(lineTokens[0].find("noEnergy") != std::string::npos) {
                                            if(atof(lineTokens[1].c_str())>0) {
                                                noEnergy = true;
                                            }
					}
					else if(lineTokens[0].find("readGeom") != std::string::npos) {
                                            if(atof(lineTokens[1].c_str())>0) {
                                                readGeom = true;
                                            }
					}
					else if(lineTokens[0].find("internalCoords") != std::string::npos) {
                                            if(atof(lineTokens[1].c_str())>0) {
                                                internalCoords = true;
                                            }
					}
					else if(lineTokens[0].find("scaleGeom") != std::string::npos) {
                                            if(atof(lineTokens[1].c_str())>0) {
                                                scaleGeom = true;
                                            }
					}
					else if(lineTokens[0].find("readOmega") != std::string::npos) {
                                            if(atof(lineTokens[1].c_str())>0) {
                                                readOmega = true;
                                            }
					}
					else {
						cout << "Unrecognized input name " << lineTokens[0] << " skipped." << endl;
					}

				}
			}
		}	// end of reading $params section

//		Read in normal modes from .pin file
		if(lineRead.find("$modes") != std::string::npos) {
			while(!(std::getline(input, lineRead).eof()) && (lineRead.find("$end") == std::string::npos)) {
				vector<string> lineTokens;
				Tokenize(lineRead, lineTokens);
				if(lineTokens.size() >= 1) {
	//				Ignore mode number
					if(lineTokens[0].find("Mode:") == std::string::npos) {
						cout << "Error in input file: $modes" << endl;
						cout << "Found line: " << lineRead << endl;
						exit(-1);
					}
	//				Set up frequencies
					std::getline(input, lineRead);
					lineTokens.clear();
					Tokenize(lineRead, lineTokens);
					if(lineTokens[0].find("Frequency:") == std::string::npos) {
						cout << "Error in input file: $modes" << endl;
						cout << "Found line: " << lineRead << endl;
						exit(-1);
					}
					for(int w = 1; w < int(lineTokens.size()); w++) {
						omega.push_back(atof(lineTokens[w].c_str())*0.00000455633);			// 0.00000455633 is the conversion factor from cm^-1 to (atomic seconds)^-1
					}
	//				Ignore force constants
					std::getline(input, lineRead);
					lineTokens.clear();
					Tokenize(lineRead, lineTokens);
					if(lineTokens[0].find("Force") == std::string::npos) {
						cout << "Error in input file: $modes" << endl;
						cout << "Found line: " << lineRead << endl;
						exit(-1);
					}


	//				Set up reduced masses of oscillators
					std::getline(input, lineRead);
					lineTokens.clear();
					Tokenize(lineRead, lineTokens);
					if(lineTokens[0].find("Red.") == std::string::npos) {
						cout << "Error in input file: $modes" << endl;
						cout << "Found line: " << lineRead << endl;
						exit(-1);
					}
					for(int w = 2; w < int(lineTokens.size()); w++) {
						mass.push_back(atof(lineTokens[w].c_str())*1822.8886);			// 1822.8886 is the conversion factor from amu to atomic mass units (i.e. electron masses)
					}

	//				Ignore IR Active
					std::getline(input, lineRead);
					lineTokens.clear();
					Tokenize(lineRead, lineTokens);
					if(lineTokens[0].find("IR") == std::string::npos) {
						cout << "Error in input file: $modes" << endl;
						cout << "Found line: " << lineRead << endl;
						exit(-1);
					}

	//				Ignore IR Intens
					std::getline(input, lineRead);
					lineTokens.clear();
					Tokenize(lineRead, lineTokens);
					if(lineTokens[0].find("IR") == std::string::npos) {
						cout << "Error in input file: $modes" << endl;
						cout << "Found line: " << lineRead << endl;
						exit(-1);
					}

	//				Ignore Raman Active
					std::getline(input, lineRead);
					lineTokens.clear();
					Tokenize(lineRead, lineTokens);
					if(lineTokens[0].find("Raman") == std::string::npos) {
						cout << "Error in input file: $modes" << endl;
						cout << "Found line: " << lineRead << endl;
						exit(-1);
					}


	//				Skip XYZ line
					std::getline(input, lineRead);
					lineTokens.clear();
					Tokenize(lineRead, lineTokens);
					if(lineTokens[0].find("X") == std::string::npos) {
						cout << "Error in input file: $modes" << endl;
						cout << "Found line: " << lineRead << endl;
						exit(-1);
					}

	//				Store harmonic modes
					int columns = 0;
					vector<double> tempModeInner;
					vector< vector<double> > tempModeMid[3];
					while(!(std::getline(input, lineRead).eof()) && (lineRead.find("TransDip") == std::string::npos)) {			// Important: Keeps reading until TransDip line
						lineTokens.clear();
						Tokenize(lineRead, lineTokens);
						if(lineTokens.size()%3 != 1) {
							cout << "Error in input file: $modes" << endl;
							cout << "Found line: " << lineRead << endl;
							exit(-1);
						}
						columns = (lineTokens.size()-1)/3;
						for(int j = 0; j < columns; j++) {
							for(int i = 0; i < 3; i++) {
								tempModeInner.push_back(atof(lineTokens[1+i+3*j].c_str()));
							}
							tempModeMid[j].push_back(tempModeInner);
							tempModeInner.clear();
						}
					}
					for(int j=0; j<columns; j++){
						modes.push_back(tempModeMid[j]);
						tempModeMid[j].clear();
					}
				} //end of not empty line check
			}
		} //end of reading $modes section
	} // end of reading .pin file
//**********************
//* Variable Cleanup
//**********************
// Set up job length parameters(
      if(autoCorrLen<=0 || maxSim > 1) {
         cout << "Warning: Not using new standard of autoCorrLen > 0 and maxSim = 1. Continue at your own risk!" << endl;
      }
      if(initLen<0 && autoCorrLen > 0) {
         initLen = autoCorrLen/4;
      }
      if(sampleStart>=0) {
         if(initLen < 0) {
            initLen = sampleStart;
         }
         else if((initLen+autoCorrLen) != sampleStart) {
            cout << "Error: sampleStart != initLen+autoCorrLen" << endl;
            exit(-1);
         }
      }
      if(autoCorrLen > 0) {
         sampleStart = initLen + autoCorrLen;
      }
      else {
         sampleStart = initLen;
         numSamples = maxStep - initLen;
      }
// For autocorrelation based jobs, we don't know maxStep yet so set it after init
      maxStep = max(maxStep,sampleStart+1);
// )

	if(int(omega.size())!=int(mass.size())){
		cout << "Error: Number of reduced masses not equal to number of frequencies" << endl;
		cout << "omega dim = " << omega.size() << " but mass dim = " << mass.size() << endl;
		exit (-1);
	}
        if(nPart == 0) {
            nPart = connectivity.size();
        }
        else if(int(connectivity.size()) != nPart) {
		cout << "Error: Number of atoms not equal to number of coordinates" << endl;
		cout << "connectivity dim = " << connectivity.size() << " but nPart = " << nPart << endl;
		exit (-1);
        }
//DES Temp
//      cout << "DES: connectivity" << endl;
//      for(int i=0; i<connectivity.size(); i++) {
//         for(int j=0; j< connectivity[i].size(); j++) {
//            cout << connectivity[i][j] << "\t";
//         }
//         cout << endl;
//      }
//      exit(-1);
	if(nPart==2) {
		nMode = 1;		//Linear case!
	}
	else {
// TODO: Doesn't work for linear molecules!
		nMode = nPart*3-6;
	}

      if(levyModes==0) {
         levyModes = nMode - (physicsParams->lowFrozModes+physicsParams->highFrozModes);
      }
      else if(levyModes>nMode - (physicsParams->lowFrozModes+physicsParams->highFrozModes)) {
         cout << "Error: Total modes moved > number of modes!" << endl;
         exit(-1);
      }
      if(levyModesInit==0) {
         levyModesInit = nMode - (physicsParams->lowFrozModes+physicsParams->highFrozModes);
      }
      else if(levyModesInit>nMode - (physicsParams->lowFrozModes+physicsParams->highFrozModes)) {
         cout << "Error: Total modes moved in initialization > number of modes!" << endl;
         exit(-1);
      }
      if(physicsParams->lowFrozModes+physicsParams->highFrozModes > nMode) {
         cout << "Error: Total frozen modes > number of modes!" << endl;
         exit(-1);
      }
	coorDim = 1;
//		nMode = nPart;
//		coorDim = 3;
//		TODO: fix the following hack
//		mass.clear();
//		for(int i=0; i<nPart; i++) {
//			if(!atomType[i].compare("H")) {
//				mass.push_back(1.0*1822.8886+1.0);
//			}
//			else if(!atomType[i].compare("O")) {
//				mass.push_back(16.0*1822.8886+8.0);
//			}
//			else if(!atomType[i].compare("S")) {
//				mass.push_back(32.0*1822.8886+16.0);
//			}
//			else {
//				cout << "Unimplemented atom type. See Simulation.cpp" << endl;
//				exit(-1);
//			}
//		}
	if(physicsParams->numInit<0) {
		physicsParams->numInit = pSlice-1;
	}
	if(physicsParams->isDeltaAI() && numTI<1) {
		numTI=1;
               cout << "DES: Setting numTI to 1 for deltaAI" << endl;
	}
	outName = logFileName.substr(0, logFileName.find_first_of("."));
	epsilon = beta/((double)pSlice);
        if(epsInit > 0) {
	    epsInit /= ((double)pSlice);
	}
        if(levyNumInit < 1) {
	    levyNumInit = levyNum;
	}
        if(levyModesInit < 1) {
	    levyModesInit = levyModes;
	}


//	TODO: Remove this
//	double eZPE=0;
//	double eFull=0;
//	for(int i=0; i<nMode; i++){
//		eZPE += omega[i]/2.0;
//                eFull += omega[i]/2.0/tanh(beta*omega[i]/2.0);
//	}
//	cout << "Harmonic ZPE:\t" << eZPE << endl;
//	cout << "Harmonic Energy:\t" << eFull << endl;
/*
	for(int i = 0; i < nPart; i++) {
		for (int k = 0; k < coorDim; k++) {
			p0.push_back(0.0);
		}
		initPos.push_back(p0);
		mass.push_back(m);
		atomType.push_back(a);
	}
*/

	CoordUtil * coords = new CoordUtil(nMode, nPart, internalCoords,  modes, omega, mass, initPos, atomType, paramType, connectivity, outName, prmFile, readGeom, readOmega);
        coords->scaleGeom = scaleGeom;
	sys = new TestSystem(pSlice, beta, coords, physicsParams);
	simStats = new Stats();
	simPotStats = new Stats();
//	simComboStats = new Stats();
	energyStats = new Stats();
	potentialStats = new Stats();
//	comboStats = new Stats();
	convergenceStats = new Stats();
	acceptanceStats = new Stats();
//	xStats = new Stats();

//	posFileName = "SHO.xyz";
//	posFile.open(posFileName.c_str());
	outFileName = logFileName;
	logFile.open(outFileName.c_str(), ios::out | ios::app);			//Append logFile

//	Write Parameters to logFile
	time_t t = time(0);
	struct tm * current = localtime( & t );
//	cout << outFileName << endl;
	logFile << "-------------------------------------------------------------------------------" << endl;
	logFile << (current->tm_mon+1) << "/" << (current->tm_mday) << "/" << (current->tm_year+1900) << "  " << current->tm_hour << ":" << current->tm_min << ":" << current->tm_sec << endl;
	logFile << inFileName << endl;
	logFile << "Potential: " << sys->GetVType() << "\t Propagator: " << sys->GetRhoType() << endl;
	//logFile << "P: "<< pSlice << ", N: "<< nPart << ", Beta: " << beta << ", levyModes: " << levyModes << ", levyNum: " << levyNum << "\n" << "maxSim: " << maxSim << ", maxStep: " << maxStep << ", sampleFreq: " << sampleFreq << ", convFreq: " << convFreq << ", storeFreq: " << storeFreq << endl;
	logFile << "P: "<< pSlice << ", N: "<< nPart << ", Beta: " << beta << ", Internals: " << internalCoords << ", levyNum: " << levyNum << ", levyModes: " << levyModes << "\n" << "initLen: " << initLen << ", autoCorrLen: " << autoCorrLen << ", errorThresh: " << errorThresh << ", storeNum: " << storeNum;
        if(autoCorrLen > 0) {
            logFile << endl;
        } else {
            logFile << ", maxStep: " << maxStep << endl;
        }
        if(physicsParams->lowFrozModes+physicsParams->highFrozModes>0) {
            logFile << "Frozen Modes:" << endl;
        }
        for(int i=0; i<physicsParams->lowFrozModes; i++) {
            logFile << coords->omega[i]/0.00000455633 << "\t";
            if((i)%10==9) {
               logFile << "\n";
            }
        }
        for(int i=(int)coords->omega.size()-physicsParams->highFrozModes; i<(int)coords->omega.size(); i++) {
            logFile << coords->omega[i]/0.00000455633 << "\t";
            if((i)%10==9) {
               logFile << "\n";
            }
        }
        if(physicsParams->lowFrozModes+physicsParams->highFrozModes>0) {
            logFile << endl;
        }
        logFile << "Active Modes:" << endl;
        for(int i=physicsParams->lowFrozModes; i<(int)coords->omega.size()-physicsParams->highFrozModes; i++) {
            logFile << coords->omega[i]/0.00000455633 << "\t";
            if((i-physicsParams->lowFrozModes)%10==9) {
               logFile << "\n";
            }
        }
        logFile << endl;
        if(!coords->readGeom) {
            logFile << "Coords after opt" << endl;
            for(int j=0; j<coords->numPart; j++) {
               for(int k=0; k<3; k++) {
                  logFile << coords->initCart[j][k] << "\t";
               }
               logFile << endl;
            }
        }
        logFile << "Harmonic Energy is: " << sys->GetHarmonicE() << endl;

//	srand(time(NULL));
	idum = new int;
	*idum = -time(0);

	input.close();

//	cout << "Initialization successful! Exiting for now" << endl;
//	exit(0);
}

Simulation::~Simulation() {
//	posFile.close();
	logFile.close();
// DES Temp:
//	vFile.close();
   delete simStats;
   delete simPotStats;
   delete energyStats;
   delete potentialStats;
   delete convergenceStats;
   delete acceptanceStats;
   delete idum;
   delete sys;
}

void Simulation::TakeStep(){
//	rand();
	vector<double> x;
//	x.push_back(RandomNum::rangau(0,1,idum));
//	for(int n = 0; n < 4; n++){
//	for(int n = 0; n < 3+levyModes; n++){
        if(epsInit>0 && stepNum < initLen) {
           for(int n = 0; n < levyModesInit+1; n++){
                   x.push_back(RandomNum::rand3(idum));
           }
	    sys->Move(x, epsInit, levyNumInit, levyModesInit);
	}
        else {
           for(int n = 0; n < levyModes+1; n++){
                   x.push_back(RandomNum::rand3(idum));
           }
	    sys->Move(x, epsilon, levyNum, levyModes);
	}
}

bool Simulation::Check(double epsVal){
//	double pTest = ((double)rand()/(double)RAND_MAX);
	double pTest = (RandomNum::rand3(idum));
//	cout<<"Weight: " << sys->GetWeight() << endl;
	double checkVal = exp(-epsVal*(sys->GetWeight() - sys->GetOldWeight()));
//	double checkVal = 1.0;
//	cout << "V diff\t" << sys->GetWeight() - sys->GetOldWeight() << endl;
	if (checkVal > pTest) {
		acceptanceStats->AddVal(1.0);
//		cout << "\tT" << endl;
		return true;
	} else {
		acceptanceStats->AddVal(0.0);
//		cout << "\tF" << endl;
//		cout << "False: " << checkVal << endl;
		return false;
	}
}

void Simulation::Update() {
	sys->Forget();
}

void Simulation::Revert() {
	sys->Undo();
}

void Simulation::Sample() {
        if(numTI==0) {
// If not TI use E
           if(!noEnergy){
               autoCorr.AddVal(sys->EstimatorE());
           }
           else{
               autoCorr.AddVal(0.0);
           }
        }
        else {
           autoCorr.AddVal(sys->EstimatorAnharmonicV());
        }
        if(stepNum == initLen) {
            autoCorr.Reset();
            acceptanceStats->Reset();
            if(autoCorrLen <= 0) {
               if(storeNum>0) {
                  storeFreq = int((maxStep-initLen)/storeNum);
               }
            }
        }
// Estimate remaining steps once autocorrelation estimated
        if(stepNum == sampleStart && autoCorrLen>0) {
            //sampleFreq = int(autoCorr.GetTau()*2+1);
            tau = autoCorr.GetTau();
            sampleFreq = 1;
            //cout << "DES: Working out sampling" << endl;
            //cout << "HarmoincE = " << sys->GetHarmonicE() << endl;
            if(numTI < 1 ) {
               numSamples = int(autoCorr.GetTotalVariance()/pow(((autoCorr.GetTotalMean()-sys->GetHarmonicE())*errorThresh),2)*(2*tau)*1.96*1.96);
            }
            else {
               numSamples = int(autoCorr.GetTotalVariance()/pow((autoCorr.GetTotalMean()*errorThresh),2)*(2*tau)*1.96*1.96);
            }
            if(maxStep != sampleStart+1) {
               numSamples = min(numSamples,maxStep);
               if(numSamples == maxStep) {
                  cout << "Error: required number of steps > maxStep. Continuing for " << maxStep << " steps." << endl;
               }
            }
            if(storeNum>0) {
               storeFreq = int(numSamples/storeNum);
            }
//            if(maxStep == sampleStart) {
//               maxStep = sampleStart + numSamples*sampleFreq;
//            }
//            else {
//               maxStep = sampleStart + min(numSamples*sampleFreq,maxStep);
//            }
            if(maxStep == sampleStart+1) {
               maxStep = initLen + numSamples;
            }
            else {
               maxStep = sampleStart + min(numSamples-autoCorrLen,maxStep);
            }
// Setting numSamples to the number of remaining samples here
            numSamples = max(0,numSamples-autoCorrLen);
            time_t t = time(0);
            struct tm * current = localtime( & t );
            logFile << "****Error Estimation****" << endl;
            if(numTI < 1 ) {
               logFile << "Tau: " << tau << ", Approx. Mean: " << autoCorr.GetTotalMean()-sys->GetHarmonicE() << ", Approx. St. Dev.: " << autoCorr.GetTotalStDev() << endl;
            }
            else {
               logFile << "Tau: " << tau << ", Approx. Mean: " << autoCorr.GetTotalMean() << ", Approx. St. Dev.: " << autoCorr.GetTotalStDev() << endl;
            }
               logFile << "numSamples: " << numSamples+autoCorrLen << ", maxStep: " << maxStep << ", sampleFreq: " << sampleFreq << endl;
            logFile << (current->tm_mon+1) << "/" << (current->tm_mday) << "/" << (current->tm_year+1900) << "  " << current->tm_hour << ":" << current->tm_min << ":" << current->tm_sec << endl;
            logFile << "************************" << endl;
        }
// DES TODO: Don't estimate E if doing TI
        if(autoCorrLen <=0) {
           if((stepNum >= sampleStart)&&((stepNum-sampleStart)%sampleFreq==0)) {
                   //cout << "Sampling Step = " << (stepNum-sampleStart) << endl;
                   //energyStats->AddVal(sys->EstimatorE());
                   if(!noEnergy){
                       energyStats->AddVal(sys->EstimatorE());
                   }
                   else{
                       energyStats->AddVal(0.0);
                   }
                   if(!sys->GetPhysics()->isDeltaAI()) {
                       potentialStats->AddVal(sys->EstimatorAnharmonicV());
                   }
                   else {
                       double tempV = sys->EstimatorAnharmonicV();
                       potentialStats->AddVal(tempV);
   //                    vFile << tempNum << "\t" << tempV << endl;
   //                    tempNum++;
                   }
                   vector< vector<Particle>  > part = sys->GetParticle();
   /* xstats
                   double tempX = 0.0;
                   for(int i=0; i<(int)part.size(); i++) {
                           for(int j=0; j<(int)part[0].size(); j++) {
                                   tempX += part[i][j].pos[0];
                                   //xStats->AddVal(part[i][j].pos[0]-part[(i+1)%part.size()][j].pos[0]);
                           }
                   }
                   xStats->AddVal(tempX/part.size()/part[0].size());
   */
   // To restore function uncomment posFile above
   //		if((stepNum-sampleStart)/(double)sampleFreq < 5000) {
   //			WritePathToFile(sys->GetParticle(),posFileName);
   //		}
           }
           if(numGeomPrint > 0) {
              if(stepNum >= sampleStart && (stepNum-sampleStart+1)%((maxStep-sampleStart)/numGeomPrint)==0) {
                  //cout << "DES: sampling point = " << (maxStemp-sampleStart) << endl;
/*
                  vector< vector<Particle>  > tempParts = sys->GetParticle();
                  vector<Particle> averagePart = tempParts[0];
                  for(int j=0; j<tempParts[0].size(); j++) {
                     for(int k=0; k<tempParts[0][0].pos.size(); k++) {
                        averagePart[j].pos[k] = 0.0;
                        for(int i=0; i<tempParts.size(); i++) {
                           averagePart[j].pos[k] += tempParts[i][j].pos[k];
                        }
                        averagePart[j].pos[k] /= tempParts.size();
                     }
                  }
*/
                  ostringstream convert;
                  convert << int((stepNum-sampleStart+1)/((maxStep-sampleStart)/numGeomPrint));
                  string geomNum = convert.str();
                  string partFileName = outName + "." + geomNum + ".xyz";
                  //WritePartToFile(averagePart,partFileName);
                  sys->WritePartToFile(partFileName);
              }
           }
        }
 //       if(stepNum == initLen || stepNum == sampleStart) {
 //           autoCorr.Reset();
 //       }
}

void Simulation::WritePathToFile(vector< vector<Particle>  > part, string pFileName) {
	//vector< vector<Particle>  > part = sys->GetParticle();
        ofstream pFile;
        pFile.open(pFileName.c_str());
	if((int)part[0][0].pos.size()>3){
		cout << "Warning! Turn off write to file if dim > 3!!!!" << endl;
	}
	else{
		pFile << part.size()*part[0].size() << "\n" << endl;
		for(int i=0; i<(int)part.size(); i++) {
			for(int j=0; j<(int)part[0].size(); j++) {
//				pFile << 1 << "\t";
				pFile << part[i][j].atomType << "\t";
				for(int k=0; k<(int)part[0][0].pos.size(); k++) {
					pFile << part[i][j].pos[k] << "\t";
				}
				for(int k=0; k<(3-(int)part[0][0].pos.size()); k++) {
					pFile << 0 << "\t";
				}
				pFile << endl;
			}
		}
	}
        pFile.close();
}

void Simulation::WritePartToFile(vector<Particle> part, string pFileName) {
   ofstream pFile;
   pFile.open(pFileName.c_str());
   int dim = (int)part[0].pos.size();
   if(dim>3){
      cout << "Warning! Turn off write to file if dim > 3!!!!" << endl;
      exit(-1);
   }
   else{
      pFile << part.size() << "\n" << endl;
      for(int i=0; i<(int)part.size(); i++) {
         pFile << part[i].atomType << "\t";
         for(int k=0; k<dim; k++) {
            pFile << part[i].pos[k] << "\t";
         }
// For printing carts if dim == 2
         for(int k=0; k<(3-dim); k++) {
            pFile << 0 << "\t";
         }
         pFile << endl;
      }
   }
   pFile.close();
}

void Simulation::Log() {
   time_t t = time(0);
   struct tm * current = localtime( & t );
   if(numTI <= 0) {
      if(autoCorrLen > 0) {
         logFile << "Mean Anharmonic Energy: " << autoCorr.GetTotalMean() - sys->GetHarmonicE() << endl;
         logFile << "Energy St. Dev.: " << autoCorr.GetTotalStDev() << endl;
         if(stepNum >= maxStep) {
            logFile << "95% Confidence Bounds: " << autoCorr.GetTotalStDev()/sqrt(double(numSamples+autoCorrLen))*sqrt(tau*2)*1.96 << endl;
         }
      }
      else {
         logFile << "Mean Energy: " << energyStats->GetMean() << endl;
         logFile << "Energy St. Dev.: " << energyStats->GetStDev() << endl;
         logFile << "Mean Potential: " << potentialStats->GetMean() << endl;
         logFile << "Potential St. Dev.: " << potentialStats->GetStDev() << endl;
      }
   }
   else {
      if(autoCorrLen > 0) {
         logFile << "Mean Potential: " << autoCorr.GetTotalMean() << endl;
         logFile << "Potential St. Dev.: " << autoCorr.GetTotalStDev() << endl;
         if(stepNum >= maxStep) {
            logFile << "95% Confidence Bounds: " << autoCorr.GetTotalStDev()/sqrt(double(numSamples+autoCorrLen))*sqrt(tau*2)*1.96 << endl;
         }
      }
      else {
         logFile << "Mean Potential: " << potentialStats->GetMean() << endl;
         logFile << "Potential St. Dev.: " << potentialStats->GetStDev() << endl;
      }
   }
   //	logFile << "Mean X Position: " << xStats->GetMean() << endl;
   //	logFile << "RMS X Position: " << xStats->GetRMS() << endl;
   logFile << "Acceptance probability: " << acceptanceStats->GetMean() << endl;
   logFile << (current->tm_mon+1) << "/" << (current->tm_mday) << "/" << (current->tm_year+1900) << "  " << current->tm_hour << ":" << current->tm_min << ":" << current->tm_sec << endl;
}

void Simulation::FinalLog() {
	time_t t = time(0);
	struct tm * current = localtime( & t );
        logFile << "****** " << maxSim << " Jobs Finished ******" << endl;
	logFile << "Mean Energy: " << simStats->GetMean() << endl;
	logFile << "Mean Potential: " << simPotStats->GetMean() << endl;
	logFile << "Standard Deviation: " << simStats->GetStDev() << endl;
	logFile << "Potential Standard Deviation: " << simPotStats->GetStDev() << endl;
	logFile << (current->tm_mon+1) << "/" << (current->tm_mday) << "/" << (current->tm_year+1900) << "  " << current->tm_hour << ":" << current->tm_min << ":" << current->tm_sec << endl;
}

void Simulation::TILog(double dG, double dGstdev, double taudGstdev) {
	time_t t = time(0);
	struct tm * current = localtime( & t );
        logFile << "****** " << numTI << " TI Points Finished ******" << endl;
	logFile << "Free Energy: " << dG << endl;
	logFile << "Standard Deviation: " << dGstdev << endl;
        if(maxSim == 1) {
	    logFile << "95% Confidence Bounds: " << taudGstdev/sqrt(double(numSamples+autoCorrLen))*sqrt(2.0)*1.96 << endl;
	    //logFile << "95% Confidence Bounds: " << dGstdev/sqrt(double(numSamples+autoCorrLen))*sqrt(tau*2)*1.96 << endl;
        }
	logFile << (current->tm_mon+1) << "/" << (current->tm_mday) << "/" << (current->tm_year+1900) << "  " << current->tm_hour << ":" << current->tm_min << ":" << current->tm_sec << endl;
}

void Simulation::Store() {
   if (storeFreq!=0 && (stepNum > sampleStart) && (stepNum-sampleStart+1)%storeFreq==0){
//      time_t t = time(0);
//      struct tm * current = localtime( & t );
//      logFile << "Mean Energy: " << energyStats->GetMean() << endl;
//      logFile << "Energy Convergence: " << convergenceStats->GetStDev() << endl;
//      logFile << "Mean Potential: " << potentialStats->GetMean() << endl;
//      logFile << "Acceptance probability: " << acceptanceStats->GetMean() << endl;
//      logFile << (current->tm_mon+1) << "/" << (current->tm_mday) << "/" << (current->tm_year+1900) << "  " << current->tm_hour << ":" << current->tm_min << ":" << current->tm_sec << endl;
      Log();
   }
}


void Simulation::Run() {
   double savedMax = maxStep;
   if(numTI==0) {
      double runEps = epsilon;
      for(int simNum = 0; simNum<maxSim; simNum++) {
         if(epsInit>0) {
            runEps = epsInit;
            //cout << "DES Temp: epsInit = " << epsInit << "\t epsilon = " << epsilon << endl;
         }
         for(stepNum=0; stepNum<maxStep; stepNum++){
            if(stepNum==initLen) {
               runEps = epsilon;
            }
            TakeStep();
            if (Check(runEps)) {
               Update();
            } else {
               Revert();
            }
            Sample();
//                   if (storeFreq!=0 && (stepNum+1)%storeFreq==0 && (stepNum > (sampleStart))){
//                           Store();		//Write statistics to the log.txt file every storeFreq step
//                   }
            Store();
         }
         logFile << "****Final Tau****" << endl;
         tau = autoCorr.GetTau();
         logFile << "Tau: " << tau << "\t numSamples: " << numSamples + autoCorrLen << endl;
         if(maxSim==1) {
            logFile << " ****Simulation Finished****" << endl;
         }
         else {
            logFile << "****** Finished Simulation " << simNum << " ******" << endl;
         }
         Log();
         simStats->AddVal(energyStats->GetMean());
         simPotStats->AddVal(potentialStats->GetMean());
//           simComboStats->AddVal(comboStats->GetMean());
         if(maxSim>1) {
            sys->Reset();
            maxStep = savedMax;
         }
         energyStats->Reset();
         potentialStats->Reset();
//           comboStats->Reset();
      }
      if(maxSim>1) {
         FinalLog();
      }
   }
// Running thermodynamic integration
   else {
      vector<double> gPoints;
      vector<double> gWeights;
      if(numTI==11) {
         GetLinearQuad(numTI,gPoints,gWeights);
      }
      else {
         GetGaussianQuad(numTI,gPoints,gWeights);
      }
//        double energyPoints = 0.0;          //TODO: Remove this
//        double energyStDevPoints = 0.0;     //TODO: Remove this
      double potentialPoints = 0.0;
      double potentialStDevPoints = 0.0;
      double tauScaledPoints = 0.0;
      for(int pLam=0; pLam<numTI; pLam++) {
//      Set the lambda value remembering to convert from [-1,1] range to [0,1] range
         sys->GetPhysics()->lambdaTI=0.5*(gPoints[pLam]+1);
         double runEps = epsilon;
         for(int simNum = 0; simNum<maxSim; simNum++) {
            if(epsInit>0) {
               runEps = epsInit;
               //cout << "DES Temp: epsInit = " << epsInit << "\t epsilon = " << epsilon << endl;
            }
            for(stepNum=0; stepNum<maxStep; stepNum++) {
               if(stepNum==initLen) {
                  runEps = epsilon;
               }
               TakeStep();
               if (Check(runEps)) {
                  Update();
               } else {
                  Revert();
               }
               Sample();
//                       if (storeFreq!=0 && (stepNum+1)%storeFreq==0 && (stepNum > (sampleStart))){
//                            Store();		//Write statistics to the log.txt file every storeFreq step
//                       }
               Store();		//Write statistics to the log.txt file every storeFreq step
            }
            logFile << "****Final Tau****" << endl;
            tau = autoCorr.GetTau();
            logFile << "Tau: " << tau << "\t numSamples: " << numSamples + autoCorrLen << endl;
            if(maxSim==1) {
               logFile << " ****Simulation Finished****" << endl;
            }
            else {
               logFile << "****** Finished Simulation " << simNum << " ******" << endl;
            }
            Log();
            if(maxSim==1) {
//               if(autoCorrLen > 0) {
               potentialPoints += autoCorr.GetTotalMean() * 0.5 * gWeights[pLam];
               potentialStDevPoints += autoCorr.GetTotalStDev() * 0.5 * gWeights[pLam];
               tauScaledPoints += sqrt(tau)*autoCorr.GetTotalStDev() * 0.5 * gWeights[pLam];
//               }
//               else {
//                  potentialPoints += potentialStats->GetMean() * 0.5 * gWeights[pLam];
//                  potentialStDevPoints += potentialStats->GetStDev() * 0.5 * gWeights[pLam];
//               }
            }
            else {
               simStats->AddVal(energyStats->GetMean());
               simPotStats->AddVal(potentialStats->GetMean());
            }
            if(maxSim > 1 || numTI > 1) {
               sys->Reset();
            }
            energyStats->Reset();
            potentialStats->Reset();
//               comboStats->Reset();
         }
         logFile << "Run with lambda = " << sys->GetPhysics()->lambdaTI << endl;
         if(maxSim>1) {
            FinalLog();
            potentialPoints += simPotStats->GetMean() * 0.5 * gWeights[pLam];
            potentialStDevPoints += simPotStats->GetStDev() * 0.5 * gWeights[pLam];
            simStats->Reset();
            simPotStats->Reset();
            maxStep = savedMax;
         }
//            energyPoints += simStats->GetMean() * 0.5 * gWeights[pLam];
//            energyStDevPoints += simStats->GetStDev() * 0.5 * gWeights[pLam];
//            simComboStats->Reset();
         maxStep = savedMax;
      }
      TILog(potentialPoints, potentialStDevPoints, tauScaledPoints);
   }
}




// For given number of points, nti, return points and weight for numerical integration
// using Gauss-Legendre quadrature http://en.wikipedia.org/wiki/Gaussian_quadrature.
// Uses limited table for now, but could implement Gautschi's theorem for general soln.
// Gives exact integral for polynomials of degree nti.
void Simulation::GetGaussianQuad(int nti, vector<double>& points, vector<double>& weights) {
    points.clear();
    weights.clear();
    if(nti==0) {
        cout << "ERROR: nti = 0 in Simulation::GetGaussianQuad" << endl;
        exit(-1);
    }
    else if(nti==1) {
        points.push_back(0.0);
        weights.push_back(2.0);
    }
    else if(nti==2) {
        points.push_back(-0.577350);
        weights.push_back(1.0);
        points.push_back(0.577350);
        weights.push_back(1.0);
    }
    else if(nti==3) {
        points.push_back(-0.774596);
        weights.push_back(0.5555556);
        points.push_back(0.0);
        weights.push_back(0.888889);
        points.push_back(0.774596);
        weights.push_back(0.5555556);
    }
    else if(nti==4) {
        points.push_back(-0.861136);
        weights.push_back(0.34785);
        points.push_back(-0.339981);
        weights.push_back(0.652145);
        points.push_back(0.339981);
        weights.push_back(0.652145);
        points.push_back(0.861136);
        weights.push_back(0.34785);
    }
    else if(nti==5) {
      weights.push_back(0.2369268850561891);
      points.push_back(-0.9061798459386640);
      weights.push_back(0.4786286704993665);
      points.push_back(-0.5384693101056831);
      weights.push_back(0.5688888888888889);
      points.push_back(0.0000000000000000);
      weights.push_back(0.4786286704993665);
      points.push_back(0.5384693101056831);
      weights.push_back(0.2369268850561891);
      points.push_back(0.9061798459386640);
    }
    else {
        cout << "ERROR: nti > 5 not implemented yet in Simulation::GetGaussianQuad" << endl;
        exit(-1);
    }
}

void Simulation::GetLobattoQuad(int nti, vector<double>& points, vector<double>& weights) {
    points.clear();
    weights.clear();
    if(nti<3) {
        cout << "ERROR: nti < 3 in Simulation::GetLobattoQuad" << endl;
        exit(-1);
    }
    else if(nti==3) {
        points.push_back(-1.0);
        weights.push_back(0.333333);
        points.push_back(0.0);
        weights.push_back(1.333333);
        points.push_back(1.0);
        weights.push_back(0.333333);
    }
    else if(nti==4) {
        points.push_back(-1.0);
        weights.push_back(0.166667);
        points.push_back(-0.447214);
        weights.push_back(0.833333);
        points.push_back(0.447214);
        weights.push_back(0.833333);
        points.push_back(1);
        weights.push_back(0.166667);
    }
    else if(nti==5) {
        points.push_back(-1.0);
        weights.push_back(0.1);
        points.push_back(-0.654654);
        weights.push_back(0.544444);
        points.push_back(0.0);
        weights.push_back(0.711111);
        points.push_back(0.654654);
        weights.push_back(0.544444);
        points.push_back(1);
        weights.push_back(0.1);
    }
    else if(nti==6) {
        points.push_back(-1.0);
        weights.push_back(0.066667);
        points.push_back(-0.765055);
        weights.push_back(0.378475);
        points.push_back(-0.285232);
        weights.push_back(0.554858);
        points.push_back(0.285232);
        weights.push_back(0.554858);
        points.push_back(0.765055);
        weights.push_back(0.378475);
        points.push_back(1);
        weights.push_back(0.066667);
    }
    else {
        cout << "ERROR: nti > 4 not implemented yet in Simulation::GetGaussianQuad" << endl;
        exit(-1);
    }
}

void Simulation::GetLinearQuad(int nti, vector<double>& points, vector<double>& weights) {
    points.clear();
    weights.clear();
    if(nti!=11) {
        cout << "ERROR: nti != 11 in Simulation::GetLobattoQuad" << endl;
        exit(-1);
    }
    else {
        points.push_back(-1.0);
        weights.push_back(1.0/11.0);
        points.push_back(-0.8);
        weights.push_back(2.0/11.0);
        points.push_back(-0.6);
        weights.push_back(2.0/11.0);
        points.push_back(-0.4);
        weights.push_back(2.0/11.0);
        points.push_back(-0.2);
        weights.push_back(2.0/11.0);
        points.push_back(0.0);
        weights.push_back(2.0/11.0);
        points.push_back(0.2);
        weights.push_back(2.0/11.0);
        points.push_back(0.4);
        weights.push_back(2.0/11.0);
        points.push_back(0.6);
        weights.push_back(2.0/11.0);
        points.push_back(0.8);
        weights.push_back(2.0/11.0);
        points.push_back(1.0);
        weights.push_back(1.0/11.0);
    }
}


//	Utility function stolen from http://www.oopweb.com/CPP/Documents/CPPHOWTO/Volume/C++Programming-HOWTO-7.html
void Simulation::Tokenize(const string& str, vector<string>& tokens, const string& delimiters)
{
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}
