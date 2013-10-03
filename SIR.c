/* pomp model file: SIR */

#include <pomp.h>
#include <R_ext/Rdynload.h>

#define gamma	(__p[__parindex[0]])
#define mu	(__p[__parindex[1]])
#define iota	(__p[__parindex[2]])
#define beta1	(__p[__parindex[3]])
#define beta2	(__p[__parindex[4]])
#define beta3	(__p[__parindex[5]])
#define beta_sd	(__p[__parindex[6]])
#define popsize	(__p[__parindex[7]])
#define rho	(__p[__parindex[8]])
#define theta	(__p[__parindex[9]])
#define S_0	(__p[__parindex[10]])
#define I_0	(__p[__parindex[11]])
#define R_0	(__p[__parindex[12]])
#define seas1	(__covars[__covindex[0]])
#define seas2	(__covars[__covindex[1]])
#define seas3	(__covars[__covindex[2]])
#define S	(__x[__stateindex[0]])
#define I	(__x[__stateindex[1]])
#define R	(__x[__stateindex[2]])
#define incid	(__x[__stateindex[3]])
#define W	(__x[__stateindex[4]])
#define cases	(__y[__obsindex[0]])
#define DS	(__f[__stateindex[0]])
#define DI	(__f[__stateindex[1]])
#define DR	(__f[__stateindex[2]])
#define Dincid	(__f[__stateindex[3]])
#define DW	(__f[__stateindex[4]])
#define Tgamma	(__pt[__parindex[0]])
#define Tmu	(__pt[__parindex[1]])
#define Tiota	(__pt[__parindex[2]])
#define Tbeta1	(__pt[__parindex[3]])
#define Tbeta2	(__pt[__parindex[4]])
#define Tbeta3	(__pt[__parindex[5]])
#define Tbeta_sd	(__pt[__parindex[6]])
#define Tpopsize	(__pt[__parindex[7]])
#define Trho	(__pt[__parindex[8]])
#define Ttheta	(__pt[__parindex[9]])
#define TS_0	(__pt[__parindex[10]])
#define TI_0	(__pt[__parindex[11]])
#define TR_0	(__pt[__parindex[12]])
#define lik	(__lik[0])

void SIR_par_trans (double *__pt, double *__p, int *__parindex)
{

  double sum;
  Tgamma = exp(gamma);
  Tmu = exp(mu);
  Tiota = exp(iota);
  Tbeta1 = exp(beta1);
  Tbeta2 = exp(beta2);
  Tbeta3 = exp(beta3);
  Tbeta_sd = exp(beta_sd);
  Trho = expit(rho);
  Ttheta = exp(theta);
  TS_0 = exp(S_0);
  TI_0 = exp(I_0);
  TR_0 = exp(R_0);
  sum = TS_0+TI_0+TR_0;
  TS_0 /= sum;
  TI_0 /= sum;
  TR_0 /= sum;
 
}


void SIR_par_untrans (double *__pt, double *__p, int *__parindex)
{

  double sum;
  Tgamma = log(gamma);
  Tmu = log(mu);
  Tiota = log(iota);
  Tbeta1 = log(beta1);
  Tbeta2 = log(beta2);
  Tbeta3 = log(beta3);
  Tbeta_sd = log(beta_sd);
  Trho = logit(rho);
  Ttheta = log(theta);
  sum = S_0+I_0+R_0;
  TS_0 = log(S_0/sum);
  TI_0 = log(I_0/sum);
  TR_0 = log(R_0/sum);
 
}


void SIR_rmeasure (double *__y, double *__x, double *__p, int *__obsindex, int *__stateindex, int *__parindex, int *__covindex, int __ncovars, double *__covars, double t)
{
 
  double prob = theta/(theta+rho*incid);
  cases = rnbinom(theta,prob);
 
}


void SIR_dmeasure (double *__lik, double *__y, double *__x, double *__p, int give_log, int *__obsindex, int *__stateindex, int *__parindex, int *__covindex, int __ncovars, double *__covars, double t)
{
 
  double prob = theta/(theta+rho*incid);
  lik = dnbinom(cases,theta,prob,give_log);
 
}


void SIR_stepfn (double *__x, const double *__p, const int *__stateindex, const int *__parindex, const int *__covindex, int __covdim, const double *__covars, double t, double dt)
{
	void (*reulermultinom)(int,double,double*,double,double*);
reulermultinom = (void (*)(int,double,double*,double,double*)) R_GetCCallable("pomp","reulermultinom");

  int nrate = 6;
  double rate[nrate];		// transition rates
  double trans[nrate];		// transition numbers
  double beta;			// transmission rate
  double dW;			// white noise increment
  int k;

  beta = beta1*seas1+beta2*seas2+beta3*seas3;

  // compute the environmental stochasticity
  dW = rgammawn(beta_sd,dt);

  // compute the transition rates
  rate[0] = mu*popsize;		// birth into susceptible class
  rate[1] = (iota+beta*I*dW/dt)/popsize; // force of infection
  rate[2] = mu;			// death from susceptible class
  rate[3] = gamma;		// recovery
  rate[4] = mu;			// death from infectious class
  rate[5] = mu; 		// death from recovered class

  // compute the transition numbers
  trans[0] = rpois(rate[0]*dt);	// births are Poisson
  reulermultinom(2,S,&rate[1],dt,&trans[1]);
  reulermultinom(2,I,&rate[3],dt,&trans[3]);
  reulermultinom(1,R,&rate[5],dt,&trans[5]);

  // balance the equations
  S += trans[0]-trans[1]-trans[2];
  I += trans[1]-trans[3]-trans[4];
  R += trans[3]-trans[5];
  incid += trans[3];		// cases are cumulative recoveries
  if (beta_sd > 0.0) W += (dW-dt)/beta_sd; // increment has mean = 0, variance = dt
 
}


void SIR_skelfn (double *__f, double *__x, double *__p, int *__stateindex, int *__parindex, int *__covindex, int __ncovars, double *__covars, double t)
{

  int nrate = 6;
  double rate[nrate];		// transition rates
  double term[nrate];		// transition numbers
  double beta;			// transmission rate
  double dW;			// white noise increment
  int k;
  
  beta = beta1*seas1+beta2*seas2+beta3*seas3;

  // compute the transition rates
  rate[0] = mu*popsize;		// birth into susceptible class
  rate[1] = (iota+beta*I)/popsize; // force of infection
  rate[2] = mu;			// death from susceptible class
  rate[3] = gamma;		// recovery
  rate[4] = mu;			// death from infectious class
  rate[5] = mu; 		// death from recovered class

  // compute the several terms
  term[0] = rate[0];
  term[1] = rate[1]*S;
  term[2] = rate[2]*S;
  term[3] = rate[3]*I;
  term[4] = rate[4]*I;
  term[5] = rate[5]*R;

  // assemble the differential equations
  DS = term[0]-term[1]-term[2];
  DI = term[1]-term[3]-term[4];
  DR = term[3]-term[5];
  Dincid = term[3];		// accumulate the new I->R transitions
  DW = 0;
 
}

#undef gamma
#undef mu
#undef iota
#undef beta1
#undef beta2
#undef beta3
#undef beta_sd
#undef popsize
#undef rho
#undef theta
#undef S_0
#undef I_0
#undef R_0
#undef seas1
#undef seas2
#undef seas3
#undef S
#undef I
#undef R
#undef incid
#undef W
#undef cases
#undef DS
#undef DI
#undef DR
#undef Dincid
#undef DW
#undef Tgamma
#undef Tmu
#undef Tiota
#undef Tbeta1
#undef Tbeta2
#undef Tbeta3
#undef Tbeta_sd
#undef Tpopsize
#undef Trho
#undef Ttheta
#undef TS_0
#undef TI_0
#undef TR_0
