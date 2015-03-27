/*
 * ClassicalSys.cpp
 *
 *  Created on: May 7, 2012
 *      Author: dstuck
 */

#include "TestSystem.h"

TestSystem::TestSystem(int pSlice, double beta, CoordUtil* coords, PhysicsUtil * phys) : System(), physics(phys) {
// DES Temp:
   tempNum = 0;
   string tempstr = "potential.txt";
   if(phys->isDeltaAI()) {
      vFile.open(tempstr.c_str());
   }

   N = coords->numModes;
   P = pSlice;
   eps = beta/((double)P);
   coorDim = 1;
   vector< vector<Particle> > part_init(P,vector<Particle>(N));
   part = part_init;
   oldPart = part_init;

   upToDate.resize(P, true);
   sliceV.resize(P, 0.0);
   oldSliceV.resize(P, 0.0);
   avgV = 0.0;
   numSteps = 0;
//        for(int i=0; i<(int)w.size(); i++) {
//            coords.omega[i] = w[i];
//        }
//	V = new V_UCHO(w);
//	vector<double> lambda;
//	for(int j=0; j<N; j++){
//		lambda.push_back(5);
//	}
//	V = new V_AQO(w,lambda);
//	V = new V_Morse(12.5,2,N);
//	V = new V_Morse(0.033,2.0,N);		//Corresponds to reasonable molecular values (actually too stiff a spring!)
//	V = new V_Morse(0.00298,2.0,N);		//Corresponds to reasonable molecular values
//	V = new V_Morse(0.176, 1.02, N);	//Corresponds to H2 morse potential from http://www1.uprh.edu/rbaretti/MorsePotential15mar2011.htm but 'a' corresponds to 3000 cm-1
//	V = new V_Morse(0.176, 1.4886, N);	//Corresponds to H2 with omega 4522.42
//	V = new V_TinkerExecutable(coords);
//	rho = new Rho_Free();
//	rho = new Rho_HO(w);
//	rho = new Rho_TruncHO(w);
   CoordUtil * qchemCoords;
   if(phys->isDeltaAI()) {
// We clone coords here since unless readOmega and readGeom are set, they will be overwritten by V_Tinker constructor
      qchemCoords = coords->Clone();
   }

   if(!physics->vType.compare("V_TinkerExecutable")) {
      V = new V_TinkerExecutable(coords);
   }
   else if(!physics->vType.compare("V_Tinker")) {
      V = new V_Tinker(coords, eps, beta);
   }
   else if(!physics->vType.compare("V_Morse")) {
      if(physics->morseAlpha != 0.0 && physics->morseDE != 0.0) {
         V = new V_Morse(coords, physics->morseDE, physics->morseAlpha, N);
      }
      else {
         V = new V_Morse(coords, 0.176, 1.4886, N);		//H2
      }
   }
   else if(!physics->vType.compare("V_UCHO")) {
      V = new V_UCHO(coords->omega);
   }
   else {
      cout << "Error in selecing V" << endl;
      exit(-1);
   }
   harmonicE = 0.0;
   for(int i=0; i<N; i++) {
      harmonicE += coords->omega[i]/2.0;
   }
   //cout << "Harmonic Energy is: " << harmonicE << endl;

   coords->initInternals();

   if(phys->isDeltaAI()) {
      V2 = new V_QChem(qchemCoords,coords,eps,beta,physics->charge,physics->multiplicity);
      qchemCoords->initInternals();
      if(coords->scaleGeom) {
         qchemCoords->MakeScalingVec(coords);
      }
   }

   if(!physics->rhoType.compare("Rho_HO")) {
      rho = new Rho_HO(coords->omega,physics->lowFrozModes,physics->highFrozModes);
   }
   else if(!physics->rhoType.compare("Rho_Free")) {
      rho = new Rho_Free();
   }
   else if(!physics->rhoType.compare("Rho_TruncHO")) {
      rho = new Rho_TruncHO(coords->omega);
   }
   else {
      cout << "Invalid rhoType:\t" << physics->rhoType << endl;
      LINE
      exit(-1);
   }

   for(int i=0; i<P; i++){
      for(int j=0; j<N; j++){
         part[i][j].mass = coords->reducedMass[j];
         oldPart[i][j].mass = coords->reducedMass[j];
      }
   }
   vector< vector<double> > cart0;
   cart0 = V->GetCoordUtil()->initCart;
   for(int i=0; i<P; i++) {
      oldCarts.push_back(cart0);
      //newCarts.push_back(cart0);
   }
   newCarts = oldCarts;


//    Initialize bead polymers with random walks
   Reset();

//      for(int i=0; i<(int)coords->omega.size(); i++) {
//         cout << i << "\t" << coords->omega[i]/0.00000455633 << endl;
//      }
}

TestSystem::~TestSystem() {
// DES Temp:
   if(physics->isDeltaAI()) {
      vFile.close();
   }
   delete V;
   delete rho;
   delete physics;
}

void TestSystem::Reset() {
//*******************************************************
//*	Initialize Bead Polymer Positions with Random Walks *
//*******************************************************

//	Initialize all beads to center
   for(int j=0; j<N; j++) {
      for(int i=0; i<P; i++) {
         part[i][j].pos.clear();
         oldPart[i][j].pos.clear();
         for(int k=0; k<coorDim; k++) {
            part[i][j].pos.push_back(0.0);
            oldPart[i][j].pos.push_back(0.0);
         }
      }
   }
   for(int i=0; i<P; i++) {
      oldCarts[i] = V->GetCoordUtil()->initCart;
      newCarts[i] = oldCarts[i];
   }

//// DES Test( PES Scan
//   //double scanBegin = -0.2;
//   //double scanEnd = 0.2;
//   double scanBegin = -0.2;   //Bohr, so twice that in Angrstrom
//   double scanEnd = 0.2;
//   int scanSteps = 400;
//   double curX = 0.0;
//   int pickedMode = 0;
//   //int pickedMode2 = 3;
//// Print out frequencies
//   cout << "   \tw_Tinker\t\tw_QChem" << endl;
//   for(int i=0;i<N; i++) {
//      cout << i << "\t" << V->GetCoordUtil()->omega[i] << "\t" << V2->GetCoordUtil()->omega[i] << endl;
//   }
//// Print out masses
//   cout << "   \tm_Tinker\tm_QChem" << endl;
//   for(int i=0;i<N; i++) {
//      cout << i << "\t" << V->GetCoordUtil()->reducedMass[i] << "\t" << V2->GetCoordUtil()->reducedMass[i] << endl;
//   }
//   //V2->GetCoordUtil()->MakeScalingVec(V->GetCoordUtil());
//
//   cout << "x\tV_full\tV_HO" << endl;
/////*
//   for(int n=0; n<scanSteps; n++) {
//      curX = scanBegin + (scanEnd-scanBegin)/(double)(scanSteps-1)*(double)(n);
//      part[0][pickedMode].pos[0]=curX;
////      part[0][pickedMode2].pos[0]=curX/2;
//      cout << curX << "\t" << V->GetV(part[0]) << "\t" << -rho->ModifyPotential(part[0]) << endl;
//      V2->GetV(part[0]);
//   }
////*/
//   exit(-1);
//// DES)

// DES ( independent mode deltaAI
   if(physics->isDeltaAI()) {
      vector<double> scanEnd;
      for(int i=0; i<N; i++) {
//       Set xMax to the classical turning point for n=5 or ~3.5 times the groundstate turning point
         double minW = min(V->GetCoordUtil()->omega[i],V2->GetCoordUtil()->omega[i]);
         scanEnd.push_back(sqrt(2.0/minW/V->GetCoordUtil()->reducedMass[i]*(5.0+0.5)));
               //Bohr, so twice that in Angrstrom
      }
      int scanSteps = 20;
      double curX = 0.0;
      string filepref = V->GetCoordUtil()->outFileName;
      ofstream theFile;
      ofstream theFile2;
      ofstream theFile3;
      theFile.open((filepref+".nGrid").c_str());
      theFile << "400" << endl;
      theFile.close();
      theFile.open((filepref+".omegaMM").c_str());
      theFile2.open((filepref+".omegaQM").c_str());
      theFile3.open((filepref+".xMax").c_str());
      for(int i=0; i<N; i++) {
         theFile << V->GetCoordUtil()->omega[i] << "\t";
         theFile2 << V2->GetCoordUtil()->omega[i] << "\t";
         theFile3 << scanEnd[i] << "\t";
      }
      theFile  << endl;
      theFile2 << endl;
      theFile3 << endl;
      theFile.close();
      theFile2.close();
      theFile3.close();
      theFile.open((filepref+".massMM").c_str());
      theFile2.open((filepref+".massQM").c_str());
      for(int i=0; i<N; i++) {
         theFile << V->GetCoordUtil()->reducedMass[i] << "\t";
         theFile2<< V2->GetCoordUtil()->reducedMass[i] << "\t";
      }
      theFile  << endl;
      theFile2 << endl;
      theFile.close();
      theFile2.close();
      theFile.open((filepref+".beta").c_str());
      theFile << eps*P << endl;
      theFile.close();

      theFile.open((filepref+".fullMM").c_str());
      theFile2.open((filepref+".hoMM").c_str());
      theFile3.open((filepref+".hoQM").c_str());
      V2->GetV(part[0]);
      if(rename("0_pimc.in","equib.in") != 0) {
         cout << "Error writing qchem file in TestSystem!" << endl;
      }
      
      for(int n=0; n<scanSteps; n++) {
         for(int i=0; i<N; i++) {
            curX = -scanEnd[i] + (2*scanEnd[i])/(double)(scanSteps-1)*(double)(n);
            part[0][i].pos[0]=curX;
            theFile << V->GetV(part[0]) << "\t";
            theFile2 << -rho->ModifyPotential(part[0]) << "\t";
            theFile3 << static_cast<V_QChem*>(V2)->GetVHO(part[0]) << "\t";
            V2->GetV(part[0]);
            //cout << "mv: " +toString(n)+"_pimc.in" << endl;
            //cout << "to: "+ toString(n)+"_"+toString(i)+"_"+filepref+".in" << endl;
            if(rename((toString(i+n*N+1)+"_pimc.in").c_str(),(toString(n)+"_"+toString(i)+"_"+filepref+".in").c_str()) != 0) {
               cout << "Error writing qchem file in TestSystem!" << endl;
            }
            part[0][i].pos[0]=0.0;
         }
         theFile << endl;
         theFile2 << endl;
         theFile3 << endl;
      }
      theFile.close();
      theFile2.close();
      theFile3.close();
      exit(-1);
   }
// DES)


   static int initRanSeed = -time(0);
   int bead, beadm1;
   vector <double> levyMean;
   double levySigma;
   int levyInit = physics->numInit;
   int snipBegin;
   int snipEnd;
//        for (int j=0; j<N; j++) {}
   for (int j=physics->lowFrozModes; j<N-physics->highFrozModes; j++) {
      for(int round=0; round<P/(levyInit+1); round++) {
         snipBegin = round*(levyInit+1);
//				cout << "snipBegin = " << snipBegin << endl;
         snipEnd = (snipBegin+levyInit+1)%P;
//				cout << "snipEnd = " << snipEnd << endl;
//				cout << "round = " << round << "\t j = " << j << endl;
         for (int s=0; s<levyInit; s++) {
            bead = (snipBegin+s+1)%P;
            beadm1 = (snipBegin+s)%P;
            levyMean = rho->GetLevyMean(part[beadm1][j], part[snipEnd][j], (double)(levyInit-s), eps, j);
            levySigma = rho->GetLevySigma((double)(levyInit-s), eps, j) / sqrt(part[0][j].mass);
            for(int k=0; k<coorDim; k++){
               part[bead][j].pos[k] = RandomNum::rangau(levyMean[k], levySigma, &initRanSeed);
            }
            upToDate[bead] = false;
         }
      }
      int levyLength = P%(levyInit+1) - 2;
      snipBegin = P-P%(levyInit+1);
      snipEnd = (snipBegin+levyLength+1);
      for (int s=0; s<levyLength; s++) {
         bead = (snipBegin+s+1)%P;
         beadm1 = (snipBegin+s)%P;
         levyMean = rho->GetLevyMean(part[beadm1][j], part[snipEnd][j], (double)(levyLength-s), eps, j);
         levySigma = rho->GetLevySigma((double)(levyLength-s), eps, j) / sqrt(part[0][j].mass);
         for(int k=0; k<coorDim; k++){
            part[bead][j].pos[k] = RandomNum::rangau(levyMean[k], levySigma, &initRanSeed);
         }
         upToDate[bead] = false;
      }
   }
   oldEnergy = 1000000;		// TODO: clean this up
   oldPotE = 100000;
   ECheckFlag = false;
}


double TestSystem::GetWeight() {
   if(ECheckFlag) {
      cout << "I shouldn't be here! Non-Levy flights have been deprecated" << endl;
      exit(-1);
      CalcEnergy();
      return energy;
   }
   else {
      CalcPotential();
      return potE*physics->lambdaTI;
//    return potE;
   }
}

double TestSystem::GetOldWeight() {
//	TODO: Fix this so it doesn't break if using different ECheckFlag on different steps
   if(ECheckFlag) {
      cout << "I shouldn't be here!" << endl;
      return oldEnergy;
   }
   else {
      return oldPotE*physics->lambdaTI;
   }
}

double TestSystem::EstimatorV() {
   if(!physics->isDeltaAI()) {
      //CalcPotential();			//Might not need this
      return potE/((double)P);
   }
   else {
// DES: running ab initio thermodynamic integration correction to molecular mechanics
      static int seed = -time(0);
      int pickedSlice = P*(RandomNum::rand3(&seed));
      double pickedV = V->GetV(part[pickedSlice],rho);
      double harmV = -rho->ModifyPotential(part[pickedSlice]);
      V2->GetV(part[pickedSlice]);

//DES TODO: Make this work for numTI > 1
      vFile << tempNum << "\t" << harmV << "\t" << pickedV << endl;
// DES Temp!!!!
      //vFile << tempNum << "\t" << part[pickedSlice][0].pos[0] << "\t" << harmV << "\t" << pickedV << "\t" << endl;
      tempNum++;
      
      return pickedV;
   }
//	return exp(double(-eps)*(potE));
}

double TestSystem::EstimatorE() {
//	avgV = (avgV*double(numSteps) + potE)/double(numSteps+1);
//	numSteps++;

   double est=0;
   for(int i=0; i<P; i++) {
      est += rho->Estimate(part[i], part[(i+1)%P], eps, P);
   }
   est += potE/(double)P;
//	cout << 1.0 - 1.0/(1.0+double(eps)*(avgV - potE)) << endl;
//	est *= (1.0 - 1.0/(1.0+double(eps)*(avgV - potE)));
//	est *= (1.0 - (1.0+double(eps)*potE)/(1.0+double(eps)*avgV));
//	est += potE/double(P);
//	cout << est << endl;
   return est;
}

vector < vector<Particle> > TestSystem::GetParticle() {
   return part;
}

void TestSystem::CalcEnergy(){
//	TODO: Consider storing V and only updating those that moved!!
   cout << "I shouldn't be here!" << endl;
   energy = 0;
   for(int i=0; i<P; i++){
      if(!upToDate[i]) {
         sliceV[i] = V->GetV(part[i], rho);
      }
      upToDate[i] = true;
      energy += sliceV[i];
   }
   potE = energy;
   for(int i=0; i<P; i++){
      energy += rho->GetRho(part[i], part[(i+1)%P], eps);
   }
//	cout << "Energy = " << energy << endl;
}

void TestSystem::CalcPotential(){
   potE = 0;
   for(int i=0; i<P; i++){
      if(!upToDate[i]) {
//DES: Set guessCarts and guessModes
         V->GetCoordUtil()->guessCarts = oldCarts[i];
         V->GetCoordUtil()->guessPart = oldPart[i];
         sliceV[i] = V->GetV(part[i], rho);
         //V->GetCoordUtil()->useGuess = true;
         V->GetCoordUtil()->useGuess = false;
//DES: Read out new carts from guessCarts
         newCarts[i] = V->GetCoordUtil()->guessCarts;
      }
      upToDate[i] = true;
      potE += sliceV[i];
   }
}


void TestSystem::Move(vector<double> prob, int levyNum, int levyModes){
//	cout << "Take a step" << endl;
   /*
//	Random Walks
   ECheckFlag = true;
   double stepSize = 1.0;
   int iPart = static_cast<int>(N*P*prob[1]);
   int iCoor = static_cast<int>(coorDim*prob[2]);
   part[(iPart-iPart%N)/N][(iPart%N)].pos[iCoor] += (stepSize*(prob[0]-0.5));		// For number (0,1)
//	part[(iPart-iPart%N)/N][(iPart%N)].pos[iCoor] += (stepSize*(prob[0]/2));		// For symm random number
    */
///*
   bool debugging = false;
   if(debugging) {
//		Debug movement
      ECheckFlag = false;
      for(int i=0; i<P; i++) {
         part[i][1].pos[0] -= 1.0;
         upToDate[i] = true;
      }
   }
   else {
//	Levy-Flight move
      if(P<3) {
         cout << "Warning! You should not be using Levy-Flight with P < 3!" << endl;
      }
      ECheckFlag = false;		// Check only against potential energy!
      static int seed = -time(0);
      vector <double> levyMean;
      double levySigma;
      int bead, beadm1;
      vector <int> remainingPart;
      for(int i=0; i<N; i++) {
         remainingPart.push_back(i);
      }
//	Snip length must be at least 3 and shouldn't be more than P/2+1 or else you'll just collapse it
//	snipLength-2 is the number of beads moved.
//		int maxLevy = 30;
//		int levyLength = min(static_cast<int>(prob[0]*(P-2)/2+1), maxLevy);
//	*******************************************************
      levyLength = levyNum;
//	*******************************************************
      int snipBegin = static_cast<int>(prob[levyModes]*P);
      int snipEnd = (snipBegin+levyLength+1)%P;
      int iPart, pickedPart;
      for(int i=0; i<levyModes; i++) {
         iPart = static_cast<int>((N-i-physics->lowFrozModes-physics->highFrozModes)*prob[i]);
         //int pickedPart = int(iPart/P)+physics->lowFrozModes;
         pickedPart = remainingPart[iPart+physics->lowFrozModes];
         remainingPart.erase(remainingPart.begin()+iPart+physics->lowFrozModes);
         //cout << "DES: PickedPart = " <<  pickedPart << endl;
   /*
         for(int i=0; i<N; i++) {
            double tempAvg=0.0;
            for(int j=0; j<P; j++) {
               tempAvg+=part[j][i].pos[0];
            }
            tempAvg/=P;
            cout << tempAvg << "\t";
         }
         cout << endl;
   */

   //		cout << "SnipBegin: " << snipBegin << "\tSnipEnd" << snipEnd << "\tSnipLength:" << snipLength << endl;
         if(N > physics->lowFrozModes+physics->highFrozModes) {
            for (int s=0; s<levyLength; s++) {
               bead = (snipBegin+s+1)%P;
               beadm1 = (snipBegin+s)%P;
               levyMean = rho->GetLevyMean(part[beadm1][pickedPart], part[snipEnd][pickedPart], (double)(levyLength-s), eps, pickedPart);
               levySigma = rho->GetLevySigma((double)(levyLength-s), eps, pickedPart) / sqrt(part[0][pickedPart].mass);
               for(int k=0; k<coorDim; k++){
                  part[bead][pickedPart].pos[k] = RandomNum::rangau(levyMean[k], levySigma, &seed);
               }
               upToDate[bead] = false;
            }
         }
      }
   }
//*/
}

void TestSystem::Forget() {
   oldEnergy = energy;
   oldPotE = potE;
   oldCarts = newCarts;
   for(int i=0; i<P; i++) {
      for(int j=0; j<N; j++) {
         oldPart[i][j].pos = part[i][j].pos;
      }
      oldSliceV[i] = sliceV[i];
   }
}

void TestSystem::Undo() {
   energy = oldEnergy;
   potE = oldPotE;
   newCarts = oldCarts;
   for(int i=0; i<P; i++) {
      for(int j=0; j<N; j++) {
         part[i][j].pos = oldPart[i][j].pos;
      }
      sliceV[i] = oldSliceV[i];
   }
}


double TestSystem::Debug() {
   return potE;
}

double TestSystem::GetHarmonicE() {
//   cout << "Harmonic Energy is: " << harmonicE << endl;
   return harmonicE;
}

string TestSystem::GetVType() {
   return V->GetType();
}

string TestSystem::GetRhoType() {
   return rho->GetType();
}

PhysicsUtil* TestSystem::GetPhysics() {
   return physics;
}

string TestSystem::toString(int i) {
   ostringstream convert;
   convert << i;
   return convert.str();
}

