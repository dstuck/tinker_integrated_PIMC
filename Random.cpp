//	DES: Random functions from Numerical Recipes	5/11/12

#include "Random.h"

//using namespace std;

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)


// Long period (> 2 x 1018) random number generator of L’Ecuyer with Bays-Durham shuffle
// and added safeguards. Returns a uniform random deviate between 0.0 and 1.0
// (exclusive of the endpoint values). Call with idum a negative integer to initialize;
// thereafter, do not alter idum between successive deviates in a sequence. RNMX should approximate the
// largest floating value that is less than 1.
double RandomNum::ran2(int *idum) {
  int j;
  long k;
  static long idum2=123456789;
  static long iy=0;
  static long iv[NTAB];
  double temp;
  if (*idum <= 0) { // Initialize.
    if (-(*idum) < 1) *idum=1; // Be sure to prevent idum = 0.
    else *idum = -(*idum);
    idum2=(*idum);
    for (j=NTAB+7;j>=0;j--) { // Load the shuffle table (after 8 warm-ups).
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1; // Start here when not initializing.
  *idum=IA1*(*idum-k*IQ1)-k*IR1; // Compute idum=(IA1*idum) % IM1 without
  if (*idum < 0) *idum += IM1; // overflows by Schrage’s method.
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2; // Compute idum2=(IA2*idum) % IM2 likewise.
  if (idum2 < 0) idum2 += IM2;
  j=iy/NDIV; // Will be in the range 0..NTAB-1.
  iy=iv[j]-idum2; // Here idum is shuffled, idum and idum2 are
  iv[j] = *idum; // combined to generate output.
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX; // Because users don’t expect endpoint values.
  else return temp;
}


//Translated from random.F in qchem
//C     ******************************************************************
//C     *                                                                *
//C     *  Random number generators                                      *
//C     *                                                                *
//C     *  rand3  generates random numbers on [0,1]                      *
//C     *  rangau generates gaussian-distributed random numbers with a   *
//C     *         specified mean and standard deviation                  *
//C     *                                                                *
//C     *  To use either function, initialize with one call to rand3     *
//C     *  with iseed < 0.  Thereafter do not modify iseed between       *
//C     *  calls.  Calling from C, a useful choice is iseed = -time(0)   *
//C     *  to initialize.                                                *
//C     *                                                                *
//C     *  JMH (3/05)                                                    *
//C     *                                                                *
//C     ******************************************************************

double RandomNum::rand3(int * iseed) {
	static int iseed2=123456789;
	static int iy=0;
	static int iv[NTAB];
	int j,k;
	double ret;
	if(*iseed<0) {
		*iseed=std::max(-*iseed,1);
		iseed2 = *iseed;
		for(j=NTAB+8; j>0; j--) {
			k = *iseed/IQ1;
			*iseed=IA1*(*iseed-k*IQ1)-k*IR1;
			if(*iseed<0) {
				*iseed = *iseed+IM1;
			}
			if(j<NTAB){
				iv[j-1] = *iseed;
			}
		}
		iy = iv[0];
	}
	k = *iseed/IQ1;
	*iseed=IA1*(*iseed-k*IQ1)-k*IR1;
	if(*iseed<0) {
		*iseed=*iseed+IM1;
	}
	k=iseed2/IQ2;
	iseed2=IA2*(iseed2-k*IQ2)-k*IR2;
	if(iseed2<0) {
		iseed2=iseed2+IM2;
	}
	j=1+iy/NDIV;
	iy=iv[j-1]-iseed2;
	iv[j-1] = *iseed;
	if(iy<1) {
		iy=iy+IMM1;
	}
	ret = std::min(AM*iy,RNMX);
	return ret;
}

double RandomNum::rangau(double avg, double stddev, int * iseed) {
	double rangau, fac, rsq, v1, v2, gasdev;
	static int iset;
	static double gset;
	if(*iseed<0){
		iset=0;
	}
	if(iset==0) {
		rsq=10;
		while((rsq>1)||rsq==0){
			v1=2*rand3(iseed)-1;
			v2=2*rand3(iseed)-1;
			rsq=v1*v1+v2*v2;
		}
		fac=sqrt(-2*log(rsq)/rsq);
		gset=v1*fac;
		gasdev=v2*fac;
		iset=1;
	} else {
		gasdev=gset;
		iset=0;
	}
	rangau = avg + stddev*gasdev;
	return rangau;
}
