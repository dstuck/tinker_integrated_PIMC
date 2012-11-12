#include <iostream>
#include <fstream>
#include "debug.h"
#include "Simulation.h"
using namespace std;

int main(int argc, char *argv[]){
//	Make sure only one input
	if(argc<2) {
		cout << "Error: Requires an input file!" << endl;
		exit(-1);
	}
//	Check if input exists
	std::string inFile(argv[1]);
	ifstream inp(inFile.c_str());
	if(!inp) {
		cout << "Error: File does not exist!" << endl;
		exit(-1);
	}
	std::string outFile;
	if(argc==3) {
		outFile = std::string(argv[2]);
	}
	else {
		outFile = std::string("log.txt");
	}
	std::string prmFile;
	if(argc==4) {
		prmFile = std::string(argv[2]);
	}
	else {
		prmFile = std::string("~/projects/PIMC/sulfate_water.prm");
	}
	Simulation sim(inFile, outFile, prmFile);
	sim.Run();
	return 0;
}
