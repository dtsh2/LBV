/* pomp model file: seir no mda/seas lbv */

#include <R.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>
#include <C:/Users/David Hayman/Documents/R/win-library/3.0/pomp/include/pomp.h>

// define parameters

#define BETA        (p[parindex[0]]) // transmission rate
#define MU          (p[parindex[1]]) // adult death rate
#define DELTA       (p[parindex[2]]) // juvenile mortality rate
#define ALPHA		(p[parindex[3]]) // dis induced mortality
#define RHO			(p[parindex[4]]) // prob infectious
#define SIGMA		(p[parindex[5]]) // incubation period
#define K			(p[parindex[6]]) // carrying capacity
#define EPSILON		(p[parindex[7]]) // rate of juvenile aging
#define TAU			(p[parindex[8]]) // seroconversion
#define	KAPPA		(p[parindex[9]]) // birth rate (peak size)
#define S			(p[parindex[10]]) // synchrony
#define OMEGA		(p[parindex[11]]) // pulse per yr
#define PHI			(p[parindex[12]]) // timing during year

// define states

#define SUSJ       (x[stateindex[0]]) // number of susceptible juveniles
#define EIJ			(x[stateindex[1]]) // number of exposed to infectious juveniles
#define ERJ			(x[stateindex[2]]) // number of exposed to recovered juveniles
#define INFJ       (x[stateindex[3]]) // number of infected juveniles
#define RECJ       (x[stateindex[4]]) // number of recovered juveniles
#define SUSA       (x[stateindex[5]]) // number of susceptibles adults
#define EIA			(x[stateindex[6]]) // number of exposed to infectious adults
#define ERA			(x[stateindex[7]]) // number of exposed to recovered adults
#define INFA       (x[stateindex[8]]) // number of infected adults
#define RECA       (x[stateindex[9]]) // number of recovered adults

// the process model:
// an SIR model with Euler-multinomial step,

void sir_euler_simulator (double *x, const double *p, 
			  const int *stateindex, const int *parindex, const int *covindex,
			  int covdim, const double *covar, 
			  double t, double dt)
{
  int nrate = 22; 			// number of rates
  double rate[nrate];		// transition rates
  double trans[nrate];		// transition numbers
  double N = x[0]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7]+x[8]+x[9];		// population size
  void (*reulmult)(int,double,double*,double,double*);

  // to evaluate the basis functions and compute the transmission rate, use some of 
  // pomp's C-level eulermultinomial simulator
  reulmult = (void (*)(int,double,double*,double,double*)) R_GetCCallable("pomp","reulermultinom");

  // define pi

const double pi = 3.14159265359;

// in C --- pow(a,b) to do a^b 

  // compute the transition rates
  rate[0] = (KAPPA*(1/sqrt((1/S)*pi)*exp(-pow((cos(pi*OMEGA*t-PHI)),2)/(1/S))))*(SUSA+ERA+RECA);		// approx delta function birth into susceptible class
  rate[1] = DELTA*(N/K);			// density dept mortality - Juvenile 
  rate[2] = EPSILON;		// aging from juv to ad
  rate[3] = RHO*BETA*(INFJ+INFA)/N;			// sus to exposed- infectious route
  rate[4] = (1-RHO)*BETA*(INFJ+INFA)/N; 		// sus to exposed- recovered route
  rate[5] = DELTA*(N/K);			// density dept mortality - Juvenile   
  rate[6] = SIGMA;		// incubation rate
  rate[7] = DELTA*(N/K);			// density dept mortality - Juvenile 
  rate[8] = TAU; // seroconversion rate 
  rate[9] = EPSILON;		// aging from juv to ad
  rate[10] = ALPHA;		// disease induced mortality
  rate[11] = DELTA*(N/K);			// density dept mortality - Juvenile
  rate[12] = EPSILON;		// aging from juv to ad
  rate[13] = MU*(N/K);			// density dept mortality - Adult
  rate[14] = RHO*BETA*(INFJ+INFA)/N;			// sus to exposed- infectious route
  rate[15] = (1-RHO)*BETA*(INFJ+INFA)/N; 		// sus to exposed- recovered route
  rate[16] = MU*(N/K);			// density dept mortality - Adult
  rate[17] = SIGMA;		// incubation rate
  rate[18] = MU*(N/K);			// density dept mortality - adult
  rate[19] = TAU; // seroconversion rate 
  rate[20] = ALPHA;		// disease induced mortality
  rate[21] = MU*(N/K);			// density dept mortality - adult
    
  // compute the transition numbers // in reulmult, first # is transitions, state name, rate # trans start, trans name at start
  trans[0] = rpois(rate[0]*dt);	// births are Poisson // is this correct still with ADF birth pulse?
  (*reulmult)(4,SUSJ, &rate[1],dt,&trans[1]); // euler-multinomial exits from SUSJ class - 4, death, aging, infection (2 routes, rho, 1-rho)
  (*reulmult)(2,EIJ,&rate[5],dt,&trans[5]); // euler-multinomial exits from EIJ class
  (*reulmult)(3,ERJ,&rate[7],dt,&trans[7]); // euler-multinomial exits from ERJ class
  (*reulmult)(1,INFJ,&rate[10],dt,&trans[10]); // euler-multinomial exits from INFJ class
  (*reulmult)(2,RECJ,&rate[11],dt,&trans[11]); // euler-multinomial exits from RECJ class
  (*reulmult)(3,SUSA,&rate[13],dt,&trans[13]); // euler-multinomial exits from SUSA class
  (*reulmult)(2,EIA,&rate[16],dt,&trans[16]); // euler-multinomial exits from EIA class
  (*reulmult)(2,ERA,&rate[18],dt,&trans[18]); // euler-multinomial exits from ERA class
  (*reulmult)(1,INFA,&rate[20],dt,&trans[20]); // euler-multinomial exits from INFA class
  (*reulmult)(1,RECA,&rate[21],dt,&trans[21]); // euler-multinomial exits from RECA class

  // balance the equations
  SUSJ += trans[0]-trans[1]-trans[2]-trans[3]-trans[4]; // IN births; OUT juv mort, aging, inf * 2 routes
  EIJ += trans[3]-trans[5]-trans[6]; // IN inf from SUSJ class exposed; OUT incubation and juv mort
  ERJ += trans[4]-trans[7]-trans[8]-trans[9]; // IN inf from SUSJ class; OUT serocon, juv mort, aging
  INFJ += trans[6]-trans[10]; // IN incubation; OUT dis induced mortality
  RECJ += trans[8]-trans[11]-trans[12]; // IN seroconversion; OUT aging and juv mortality
  SUSA += trans[2]-trans[13]-trans[14]-trans[15]; // IN aging from SUSJ class; OUT adult mort, inf * 2 classes
  EIA += trans[14]-trans[16]-trans[17]; // IN from adult SUS exposed; OUT incubation and adult mortality
  ERA += trans[15]+trans[9]-trans[18]-trans[19]; // IN from adult SUS exposed & aging; OUT serconverstion and adult mortality
  INFA += trans[17]-trans[20]; // IN from incubation; OUT dis induced mortality
  RECA += trans[19]+trans[12]-trans[21]; // IN from serconverstion & aging; OUT from death

}