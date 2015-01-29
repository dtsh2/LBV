/* pomp model file: seir no seas lbv */

#include <R.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>
#include <C:/Users/David Hayman/Documents/R/win-library/3.0/pomp/include/pomp.h>

// define parameters

#define BETA        (p[parindex[0]]) // transmission rate
#define RHO        (p[parindex[1]]) // transmission rate

// define states

#define SUSJ       (x[stateindex[0]]) // number of susceptible juveniles
#define MDAJ		(x[stateindex[1]]) // number of mda juveniles
#define SUSJM		(x[stateindex[2]]) // number of susceptible juveniles from mda category
#define EIJ			(x[stateindex[3]]) // number of exposed to infectious juveniles
#define ERJ			(x[stateindex[4]]) // number of exposed to recovered juveniles
#define INFJ       (x[stateindex[5]]) // number of infected juveniles
#define RECJ       (x[stateindex[6]]) // number of recovered juveniles
#define SUSA       (x[stateindex[7]]) // number of susceptibles adults
#define EIA			(x[stateindex[8]]) // number of exposed to infectious adults
#define ERA			(x[stateindex[9]]) // number of exposed to recovered adults
#define INFA       (x[stateindex[10]]) // number of infected adults
#define RECA       (x[stateindex[11]]) // number of recovered adults
#define SPA			(x[stateindex[12]]) // adult seroprevalence
#define SPJ			(x[stateindex[13]]) // juvenile seroprevalence

// define observations

#define DRECA	(y[obsindex[0]]) // data adult sero pos
#define DRECJ	(y[obsindex[1]]) // data juvenile sero pos
#define DPOPA	(y[obsindex[2]]) // data adult sampled
#define DPOPJ	(y[obsindex[3]]) // data juvenile sampled
// #define DSPA	(y[obsindex[4]]) // adult seroprevalence
// #define DSPJ	(y[obsindex[5]]) // juvenile seroprevalence

// binomial measurement error density

void binomial_dmeasure (double *lik, double *y, double *x, double *p, int give_log,
			  int *obsindex, int *stateindex, int *parindex, int *covindex,
			  int ncovars, double *covars, double t)
{
  double ppa = (SPA);
  double ppj = (SPJ);
  double f = 0.0;
   f += dbinom(DRECA,DPOPA,ppa,TRUE); // we want this to be on the log scale
   f += dbinom(DRECJ,DPOPJ,ppj,TRUE);
 *lik = (give_log) ? f : exp(f);
}


// binomial measurement error simulator

void binomial_rmeasure (double *y, double *x, double *p,
			  int *obsindex, int *stateindex, int *parindex, int *covindex,
			  int ncovars, double *covars, double t)
{
  double ppa = (SPA);
  double ppj = (SPJ); 
  DRECA = rbinom(DPOPA,ppa);
  DRECJ = rbinom(DPOPJ,ppj);
}

// the process model:
// an SIR model with Euler-multinomial step,

void sir_euler_simulator (double *x, const double *p, 
			  const int *stateindex, const int *parindex, const int *covindex,
			  int covdim, const double *covar, 
			  double t, double dt)
{
  int nrate = 29; 			// number of rates
  double rate[nrate];		// transition rates
  double trans[nrate];		// transition numbers
  double N = x[0]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7]+x[8]+x[9]+x[10]+x[11];		// population size
  void (*reulmult)(int,double,double*,double,double*);

  // to evaluate the basis functions and compute the transmission rate, use some of 
  // pomp's C-level eulermultinomial simulator
  reulmult = (void (*)(int,double,double*,double,double*)) R_GetCCallable("pomp","reulermultinom");

  // define constant values

const double pi = 3.14159265359;
const double MU=0.000510492;
const double DELTA=0.002312247;
const double ALPHA=0.2;
const double SIGMA=0.02083333;
const double K=1000000;
const double EPSILON=0.002739726;
const double TAU=0.04166667;
const double PSI=0.01;
const double KAPPA=0.004109589;
const double S=14.35;
const double OMEGA=0.002739726;
const double PHI=0.0;
const double GAMMA=0.003773585;
  
// in C --- pow(a,b) to do a^b 

  // compute the transition rates
  rate[0] = (KAPPA*(1/sqrt((1/S)*pi)*exp(-pow((cos(pi*OMEGA*t-PHI)),2)/(1/S))))*(RECA);		// approx delta function birth into maternally-derived antibody class 
  rate[1] = PSI;			// loss of maternally derived antibody
  rate[2] = DELTA*(N/K);			// density dept mortality - Juvenile 
  rate[3] = GAMMA; // aging from SUSJM to adult SUS
  rate[4] = DELTA*(N/K); // density dept mortality - juvenile
  rate[5] = RHO*BETA*(INFJ+INFA)/N;			// sus to exposed- infectious route
  rate[6] = (1-RHO)*BETA*(INFJ+INFA)/N; 		// sus to exposed- recovered route
  rate[7] = (KAPPA*(1/sqrt((1/S)*pi)*exp(-pow((cos(pi*OMEGA*t-PHI)),2)/(1/S))))*(SUSA+ERA);		// approx delta function birth into susceptible class
  rate[8] = DELTA*(N/K);			// density dept mortality - Juvenile 
  rate[9] = EPSILON;		// aging from juv to ad
  rate[10] = RHO*BETA*(INFJ+INFA)/N;			// sus to exposed- infectious route
  rate[11] = (1-RHO)*BETA*(INFJ+INFA)/N; 		// sus to exposed- recovered route
  rate[12] = DELTA*(N/K);			// density dept mortality - Juvenile   
  rate[13] = SIGMA;		// incubation rate
  rate[14] = DELTA*(N/K);			// density dept mortality - Juvenile 
  rate[15] = TAU; // seroconversion rate 
  rate[16] = EPSILON;		// aging from juv to ad
  rate[17] = ALPHA;		// disease induced mortality
  rate[18] = DELTA*(N/K);			// density dept mortality - Juvenile
  rate[19] = EPSILON;		// aging from juv to ad
  rate[20] = MU*(N/K);			// density dept mortality - Adult
  rate[21] = RHO*BETA*(INFJ+INFA)/N;			// sus to exposed- infectious route
  rate[22] = (1-RHO)*BETA*(INFJ+INFA)/N; 		// sus to exposed- recovered route
  rate[23] = MU*(N/K);			// density dept mortality - Adult
  rate[24] = SIGMA;		// incubation rate
  rate[25] = MU*(N/K);			// density dept mortality - adult
  rate[26] = TAU; // seroconversion rate 
  rate[27] = ALPHA;		// disease induced mortality
  rate[28] = MU*(N/K);			// density dept mortality - adult
    
  // compute the transition numbers // in reulmult, first # is transitions, state name, rate # trans start, trans name at start
  trans[0] = rpois(rate[0]*dt);	// births are Poisson // is this correct still with ADF birth pulse?
  (*reulmult)(2,MDAJ,&rate[1],dt,&trans[1]); // euler-multinomial exits from MDAJ class - 2, loss MDA, death
  (*reulmult)(4,SUSJM, &rate[3],dt,&trans[3]); // euler-multinomial exits from SUSJM class - 4, death, aging, infection (2 routes, rho, 1-rho)
  trans[7] = rpois(rate[7]*dt);	// births are Poisson // is this correct still with ADF birth pulse?
  (*reulmult)(4,SUSJ, &rate[8],dt,&trans[8]); // euler-multinomial exits from SUSJ class - 4, death, aging, infection (2 routes, rho, 1-rho)
  (*reulmult)(2,EIJ,&rate[12],dt,&trans[12]); // euler-multinomial exits from EIJ class
  (*reulmult)(3,ERJ,&rate[14],dt,&trans[14]); // euler-multinomial exits from ERJ class
  (*reulmult)(1,INFJ,&rate[17],dt,&trans[17]); // euler-multinomial exits from INFJ class
  (*reulmult)(2,RECJ,&rate[18],dt,&trans[18]); // euler-multinomial exits from RECJ class
  (*reulmult)(3,SUSA,&rate[20],dt,&trans[20]); // euler-multinomial exits from SUSA class
  (*reulmult)(2,EIA,&rate[23],dt,&trans[23]); // euler-multinomial exits from EIA class
  (*reulmult)(2,ERA,&rate[25],dt,&trans[25]); // euler-multinomial exits from ERA class
  (*reulmult)(1,INFA,&rate[27],dt,&trans[27]); // euler-multinomial exits from INFA class
  (*reulmult)(1,RECA,&rate[28],dt,&trans[28]); // euler-multinomial exits from RECA class

  // balance the equations
  MDAJ += trans[0]-trans[1]-trans[2]; // IN births; OUT loss mda, juv mortality
  SUSJM += trans[1]-trans[3]-trans[4]-trans[5]-trans[6]; // IN loss mda; OUT aging (own rate so age to adults same as SUSJ), juv mort, inf * 2 routes
  SUSJ += trans[7]-trans[8]-trans[9]-trans[10]-trans[11]; // IN births; OUT juv mort, aging, inf * 2 routes
  EIJ += trans[5]+trans[10]-trans[12]-trans[13]; // IN inf from both SUSJ classes exposed; OUT incubation and juv mort
  ERJ += trans[6]+trans[11]-trans[14]-trans[15]-trans[16]; // IN inf from both SUSJ classes; OUT serocon, juv mort, aging
  INFJ += trans[13]-trans[17]; // IN incubation; OUT dis induced mortality
  RECJ += trans[15]-trans[18]-trans[19]; // IN seroconversion; OUT aging and juv mortality
  SUSA += trans[3]+trans[9]-trans[20]-trans[21]-trans[22]; // IN aging from both SUSJ classes; OUT adult mort, inf * 2 classes
  EIA += trans[21]-trans[23]-trans[24]; // IN from adult SUS exposed; OUT incubation and adult mortality
  ERA += trans[22]+trans[16]-trans[25]-trans[26]; // IN from adult SUS exposed & aging; OUT serconverstion and adult mortality
  INFA += trans[24]-trans[27]; // IN from incubation; OUT dis induced mortality
  RECA += trans[26]+trans[18]-trans[28]; // IN from serconverstion & aging; OUT from death
  SPA = RECA/(SUSA+EIA+ERA+RECA); // adult seroprevalence
  SPJ = RECJ/(SUSJ+EIJ+ERJ+RECJ); // juvenile seroprevalence
}