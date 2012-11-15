/*
 * Simulation.cpp
 *
 *  Created on: May 4, 2012
 *      Author: dstuck
 */

#include "Simulation.h"

Simulation::Simulation(string inFileName, string logFileName, string prmFile) {

//*******************************************************
//*				  Default Parameters				    *
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
        maxSim = 1;
	int nPart = 0;
	int pSlice = 0;
	int nMode;
	int coorDim;
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
					else if(lineTokens[0].find("beta") != std::string::npos) {
						beta = atof(lineTokens[1].c_str());
					}
	//				TODO: merge nPart read with $molecule section
					else if(lineTokens[0].find("nPart") != std::string::npos) {
						nPart = atoi(lineTokens[1].c_str());
					}
					else if(lineTokens[0].find("pSlice") != std::string::npos) {
						pSlice = atoi(lineTokens[1].c_str());
					}
					else if(lineTokens[0].find("levyNum") != std::string::npos) {
						levyNum = atoi(lineTokens[1].c_str());
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

// A few conversions:
	if(int(omega.size())!=int(mass.size())){
		cout << "Error, number of reduced masses not equal to number of frequencies" << endl;
		cout << "omega dim = " << omega.size() << " but mass dim = " << mass.size() << endl;
		exit (-1);
	}
	if(nPart==2) {
		nMode = 1;		//Linear case!
	}
	else {
		nMode = nPart*3-6;
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
	string outName = logFileName.substr(0, logFileName.find_first_of("."));
	epsTemp = beta/((double)pSlice);

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

	CoordUtil * coords = new CoordUtil(nMode, nPart, modes, omega, mass, initPos, atomType, paramType, connectivity, outName, prmFile);
	sys = new TestSystem(pSlice, beta, coords, physicsParams);
	simStats = new Stats();
	energyStats = new Stats();
	potentialStats = new Stats();
	convergenceStats = new Stats();
	acceptanceStats = new Stats();
	xStats = new Stats();

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
	logFile << "P: "<< pSlice << ", N: "<< nPart << ", Beta: " << beta << ", levyNum: " << levyNum << "\n" << "maxSim: " << maxSim << ", maxStep: " << maxStep << ", sampleFreq: " << sampleFreq << ", convFreq: " << convFreq << ", storeFreq: " << storeFreq << endl;


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
}

void Simulation::TakeStep(){
//	rand();
	vector<double> x;
//	x.push_back(RandomNum::rangau(0,1,idum));
//	for(int n = 0; n < 2; n++){
	for(int n = 0; n < 4; n++){
		x.push_back(RandomNum::rand3(idum));
	}
	sys->Move(x, levyNum);
}

bool Simulation::Check(){
//	double pTest = ((double)rand()/(double)RAND_MAX);
	double pTest = (RandomNum::rand3(idum));
//	cout<<"Weight: " << sys->GetWeight() << endl;
	double checkVal = exp(-epsTemp*(sys->GetWeight() - sys->GetOldWeight()));
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

void Simulation::Update(){
	sys->Forget();
}

void Simulation::Revert(){
	sys->Undo();
}

void Simulation::Sample(){
	if((stepNum > sampleStart)&&(stepNum%convFreq==0)){
		convergenceStats->AddVal(sys->EstimatorE());
//		cout << "Energy: " << sys->EstimatorE() << endl;
	}
	if((stepNum > sampleStart)&&(stepNum%sampleFreq==0)) {
		energyStats->AddVal(sys->EstimatorE());
		potentialStats->AddVal(sys->EstimatorV());
		vector< vector<Particle>  > part = sys->GetParticle();
		for(int i=0; i<(int)part.size(); i++) {
			for(int j=0; j<(int)part[0].size(); j++) {
				xStats->AddVal(part[i][j].pos[0]-part[(i+1)%part.size()][j].pos[0]);
			}
		}
// To restore function uncomment posFile above
//		if((stepNum-sampleStart)/(double)sampleFreq < 5000) {
//			WritePosToFile();
//		}
	}
}

void Simulation::FinalPrint() {
	cout << " ****Simulation Finished****" << endl;
	cout << "Mean Energy: " << energyStats->GetMean() << endl;
//	cout << "Energy RMS: " << energyStats->GetRMS() << endl;
//	cout << "Energy Variance: " << energyStats->GetVariance() << endl;
//	cout << "Energy Convergence: " << convergenceStats->GetVariance() << endl;
	cout << "Energy Convergence: " << convergenceStats->GetStDev() << endl;
	cout << "Percent Convergence: " << convergenceStats->GetStDev()/convergenceStats->GetMeanAbs()*100.0 << endl;
	cout << "Mean Potential: " << potentialStats->GetMean() << endl;
	cout << "Mean X Position: " << xStats->GetMean() << endl;
	cout << "RMS X Position: " << xStats->GetRMS() << endl;
	cout << "Acceptance probability: " << acceptanceStats->GetMean() << endl;
}

void Simulation::WritePosToFile() {
	vector< vector<Particle>  > part = sys->GetParticle();
	if((int)part[0][0].pos.size()>3){
		cout << "Warning! Turn off write to file if dim > 3!!!!" << endl;
	}
	else{
		posFile << part.size()*part[0].size() << "\n" << endl;
		for(int i=0; i<(int)part.size(); i++) {
			for(int j=0; j<(int)part[0].size(); j++) {
//				posFile << 1 << "\t";
				posFile << part[i][j].atomType << "\t";
				for(int k=0; k<(int)part[0][0].pos.size(); k++) {
					posFile << part[i][j].pos[k] << "\t";
	//				cout << part[i][j].pos[k] << endl;
				}
				for(int k=0; k<(3-(int)part[0][0].pos.size()); k++) {
					posFile << 0 << "\t";
				}
				posFile << endl;
			}
		}
	}
}

void Simulation::Log() {
	time_t t = time(0);
	struct tm * current = localtime( & t );
	logFile << "Mean Energy: " << energyStats->GetMean() << endl;
	logFile << "Energy Convergence: " << convergenceStats->GetStDev() << endl;
	logFile << "Mean Potential: " << potentialStats->GetMean() << endl;
	logFile << "Mean X Position: " << xStats->GetMean() << endl;
	logFile << "RMS X Position: " << xStats->GetRMS() << endl;
	logFile << "Acceptance probability: " << acceptanceStats->GetMean() << endl;
	logFile << (current->tm_mon+1) << "/" << (current->tm_mday) << "/" << (current->tm_year+1900) << "  " << current->tm_hour << ":" << current->tm_min << ":" << current->tm_sec << endl;
}

void Simulation::FinalLog() {
	time_t t = time(0);
	struct tm * current = localtime( & t );
        logFile << "****** " << maxSim << " Jobs Finished ******" << endl;
	logFile << "Mean Energy: " << simStats->GetMean() << endl;
	logFile << "Standard Deviation: " << simStats->GetStDev() << endl;
	logFile << (current->tm_mon+1) << "/" << (current->tm_mday) << "/" << (current->tm_year+1900) << "  " << current->tm_hour << ":" << current->tm_min << ":" << current->tm_sec << endl;
}

void Simulation::Store() {
	time_t t = time(0);
	struct tm * current = localtime( & t );
	logFile << "Mean Energy: " << energyStats->GetMean() << endl;
	logFile << "Energy Convergence: " << convergenceStats->GetStDev() << endl;
	logFile << "Mean Potential: " << potentialStats->GetMean() << endl;
	logFile << "Acceptance probability: " << acceptanceStats->GetMean() << endl;
	logFile << (current->tm_mon+1) << "/" << (current->tm_mday) << "/" << (current->tm_year+1900) << "  " << current->tm_hour << ":" << current->tm_min << ":" << current->tm_sec << endl;
}


void Simulation::Run(){
      for(int simNum = 0; simNum<maxSim; simNum++) {
           for(stepNum=0; stepNum<maxStep; stepNum++){
                   TakeStep();
                   if (Check()) {
                           Update();
                   } else {
                           Revert();
                   }
                   Sample();
                   if (storeFreq!=0 && (stepNum+1)%storeFreq==0){					//Write statistics to the log.txt file every storeFreq step
                           Store();
                   }
           }
           if(maxSim==1) {
               FinalPrint();
           }
           else {
               logFile << "****** Finished Simulation " << simNum << " ******" << endl;
           }
           Log();
           simStats->AddVal(energyStats->GetMean());
           sys->Reset();
      }
      if(maxSim>1) {
         FinalLog();
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
