## LBV SEIRS model with Frequency-dept transmission
## rho/beta estimation
## clean up

rm(list=ls())

# instal r developer toolbox first (Rtools from Cran R)

## pomp test run lbv
getwd()
setwd("~/GitHub/LBV")

library(pomp)

#Compiling C code and loading the dll
system("R CMD SHLIB lbvseirsFreqRMeasure.c")

dyn.load("lbvseirsFreqRMeasure.dll")

dat <- read.csv("lbv_data.csv",header=T)
head(dat)
class(dat)

pomp(
  data = dat,
  times="cumulative_time",
  t0=0,
  ## native routine for the process simulator:
  rprocess=euler.sim(
    step.fun="sir_euler_simulator",
    delta.t=1,
    #PACKAGE="pomp"  ## may be include if does not work - this is where to look for the c file 
  ## name of the shared-object library containing the PACKAGE="pomp",
  ),
  dmeasure="binomial_dmeasure",
  rmeasure="binomial_rmeasure",
  ## the order of the state variables assumed in the native routines:
  statenames=c("SUSJ","MDAJ", "SUSJM","EIJ","ERJ","INFJ", "RECJ", "SUSA", "EIA","ERA","INFA", "RECA"),
  ## the order of the parameters assumed in the native routines:
  paramnames=c("BETA","MU","DELTA","ALPHA","RHO","SIGMA","K","EPSILON","TAU","PSI","KAPPA","S","OMEGA","PHI","GAMMA","LAMBDA","ZETA",
               "SUSJ.0","MDAJ.0", "SUSJM.0","EIJ.0","ERJ.0","INFJ.0", "RECJ.0", "SUSA.0", "EIA.0","ERA.0","INFA.0", "RECA.0"),
  initializer=function(params,t0,statenames,...){
    x0<-params[paste(statenames,".0",sep="")]
    names(x0)<-statenames
    return(x0)
  }
) -> sir

params <- c(
  BETA=10, # 0
  MU=0.000510492, #1
  DELTA=0.002312247, #2
  ALPHA=0.2, #3
  RHO=0.05, #4 
  SIGMA=1/48, #5
  K=1000000, #6
  EPSILON=1/365, #7
  TAU=1/24, #8
  PSI=0.01, #9
  KAPPA=1.5/365, #10
  S=77.82, #11
  OMEGA=1/365, #12
  PHI=0.5, #13
  GAMMA=0.003225806, # check data #14
  LAMBDA=0.001, # check data #15
  ZETA=0.5, #16
  OBSER=0.1, #17
  SUSJ.0=4000,MDAJ.0=4000, SUSJM.0=1000,EIJ.0=100,ERJ.0=1000,INFJ.0=100,
  RECJ.0=10000,SUSA.0=50000, EIA.0=100,
  ERA.0=1000,INFA.0=500, RECA.0=50000)

sim <- simulate(sir,params=c(params),seed=3493885L,nsim=1,states=T,obs=F,as.data.frame=T) # double check seed in this

plot(sim$time,sim$SUSA,type="l",ylim=c(0,max(sim$RECA)))
points(sim$time,sim$RECA,col="green",type="l")
points(sim$time,sim$INFA,col="red",type="l")
points(sim$time,sim$SUSJ,lty=2,type="l")
points(sim$time,sim$RECJ,col="green",lty=2,type="l")
points(sim$time,sim$MDA,col="brown",lty=2,type="l")
points(sim$time,sim$INFJ,col="red",lty=2,type="l")

## to check can run for >1 simulation
# sims <- simulate(sir,params=c(params),seed=3493885L,nsim=10,states=T,obs=F,as.data.frame=T) # 100 simulations of 1 parameter set.

class(sim)
class(sir)

# try pfilter
#
pf<-pfilter(sir,params=params,Np=1000)
# loglik<-logLik(pf)

# now try dprocess
# pars wish to estimate
estpars1<-c("BETA","RHO","BIAS") # what I wish to estimate
par.start1<-c(BETA=10,RHO=0.05,BIAS=0.1)
# but try with all pars
estpars2=c("BIAS","BETA","MU","DELTA","ALPHA","RHO","SIGMA","K","EPSILON","TAU","PSI","KAPPA","S","OMEGA","PHI","GAMMA","LAMBDA","ZETA",
             "SUSJ.0","MDAJ.0", "SUSJM.0","EIJ.0","ERJ.0","INFJ.0", "RECJ.0", "SUSA.0", "EIA.0","ERA.0","INFA.0", "RECA.0")
par.start2=c(
  BIAS=0.1,
  BETA=10,
  MU=0.000510492,
  DELTA=0.002312247,
  ALPHA=0.2,
  RHO=0.05,
  SIGMA=1/48,
  K=1000000,
  EPSILON=1/365,
  TAU=1/24,
  PSI=0.01,
  KAPPA=1.5/365,
  S=77.82,
  OMEGA=1/365,
  PHI=0.5,
  GAMMA=0.003225806, # check data
  LAMBDA=0.001, # check data
  ZETA=0.5,
  SUSJ.0=4000,MDAJ.0=4000, SUSJM.0=1000,EIJ.0=100,ERJ.0=1000,INFJ.0=100,
  RECJ.0=10000,SUSA.0=50000, EIA.0=100,
  ERA.0=1000,INFA.0=500, RECA.0=50000)
  
replicate(
  n=2,      
      {
        mif(
          sir,
          Nmif=10,
          start=par.start2,
          pars<-estpars2,
          rw.sd=c(
            BIAS=0.01,
            BETA=0.1,
            MU=0.0001,
            DELTA=0.001,
            ALPHA=0.1,
            RHO=0.01,
            SIGMA=0.001,
            K=10,
            EPSILON=0.0001,
            TAU=0.001,
            PSI=0.011,
            KAPPA=0.001,
            S=1,
            OMEGA=0.001,
            PHI=0.01,
            GAMMA=0.001, # check data
            LAMBDA=0.001, # check data
            ZETA=0.1,
            SUSJ.0=10,MDAJ.0=10, SUSJM.0=10,EIJ.0=10,ERJ.0=10,INFJ.0=10,
            RECJ.0=10,SUSA.0=10, EIA.0=10,
            ERA.0=10,INFA.0=10, RECA.0=10),
          Np=100,
          var.factor=4,
          ic.lag=10,
          cooling.factor=0.999,
          max.fail=10) # 
})->mf
