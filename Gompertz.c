/* pomp model file: Gompertz */

#include <pomp.h>
#include <R_ext/Rdynload.h>

#define r	(__p[__parindex[0]])
#define K	(__p[__parindex[1]])
#define sigma	(__p[__parindex[2]])
#define X_0	(__p[__parindex[3]])
#define tau	(__p[__parindex[4]])
#define X	(__x[__stateindex[0]])
#define Y	(__y[__obsindex[0]])
#define DX	(__f[__stateindex[0]])
#define Tr	(__pt[__parindex[0]])
#define TK	(__pt[__parindex[1]])
#define Tsigma	(__pt[__parindex[2]])
#define TX_0	(__pt[__parindex[3]])
#define Ttau	(__pt[__parindex[4]])
#define lik	(__lik[0])

void Gompertz_par_trans (double *__pt, double *__p, int *__parindex)
{

  Tr = log(r);
  TK = log(K);
  Tsigma = log(sigma);
  TX_0 = log(X_0);
  Ttau = log(tau);
 
}


void Gompertz_par_untrans (double *__pt, double *__p, int *__parindex)
{

  Tr = exp(r);
  TK = exp(K);
  Tsigma = exp(sigma);
  TX_0 = exp(X_0);
  Ttau = exp(tau);
 
}


void Gompertz_rmeasure (double *__y, double *__x, double *__p, int *__obsindex, int *__stateindex, int *__parindex, int *__covindex, int __ncovars, double *__covars, double t)
{
 
    Y = rlnorm(log(X),tau);
 
}


void Gompertz_dmeasure (double *__lik, double *__y, double *__x, double *__p, int give_log, int *__obsindex, int *__stateindex, int *__parindex, int *__covindex, int __ncovars, double *__covars, double t)
{
 
    lik = dlnorm(Y,log(X),tau,give_log);
 
}


void Gompertz_stepfn (double *__x, const double *__p, const int *__stateindex, const int *__parindex, const int *__covindex, int __covdim, const double *__covars, double t, double dt)
{

  double S = exp(-r*dt);
  double logeps = (sigma > 0.0) ? rnorm(0,sigma) : 0.0;
  /* note that X is over-written by the next line */
  X = pow(K,(1-S))*pow(X,S)*exp(logeps); 
 
}


void Gompertz_skelfn (double *__f, double *__x, double *__p, int *__stateindex, int *__parindex, int *__covindex, int __ncovars, double *__covars, double t)
{

  double dt = 1.0;
  double S = exp(-r*dt);
  /* note that X is not over-written in the skeleton function */
  DX = pow(K,1-S)*pow(X,S); 
 
}

#undef r
#undef K
#undef sigma
#undef X_0
#undef tau
#undef X
#undef Y
#undef DX
#undef Tr
#undef TK
#undef Tsigma
#undef TX_0
#undef Ttau
