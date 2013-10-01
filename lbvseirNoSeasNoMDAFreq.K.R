## LBV SEIR model with Frequency-dept transmission
## no MDA
## K
## clean up

rm(list=ls())

# instal r developer toolbox first (Rtools from Cran R)

## pomp test run lbv
getwd()
setwd("~/Cambridge/CSU 2013/LBV Model/lbvmodels/new_models_sept")

library(pomp)

#Compiling C code and loading the dll
system("R CMD SHLIB lbvseirNoSeasNoMDAFreq.c")

dyn.load("lbvseirNoSeasNoMDAFreq.dll")

pomp(
  data = data.frame(
    time=seq(from=0,to=365*30,by=1),  # time for simulations to run
    X = NA # dummy variables
    ),
  times="time",
  t0=0,
  ## native routine for the process simulator:
  rprocess=euler.sim(
    step.fun="sir_euler_simulator",
    delta.t=1,
    #PACKAGE="pomp"  ## may be include if does not work - this is where to look for the c file 
  ## name of the shared-object library containing the PACKAGE="pomp",
  ),
  
  ## the order of the state variables assumed in the native routines:
  statenames=c("SUSJ","EIJ","ERJ","INFJ", "RECJ", "SUSA", "EIA","ERA","INFA", "RECA"),
  ## the order of the parameters assumed in the native routines:
  paramnames=c("BETA","MU","DELTA","ALPHA","RHO","SIGMA","K","EPSILON","TAU","KAPPA","S","OMEGA","PHI",
               "SUSJ.0","EIJ.0","ERJ.0","INFJ.0", "RECJ.0", "SUSA.0", "EIA.0","ERA.0","INFA.0", "RECA.0"),
  initializer=function(params,t0,statenames,...){
    x0<-params[paste(statenames,".0",sep="")]
    names(x0)<-statenames
    return(x0)
  }
) -> sir

params <- c(
  BETA=8,
  MU=0.000510492,
  DELTA=0.002312247,
  ALPHA=0.2,
  RHO=0.02,
  SIGMA=1/48,
  K=1000000,
  EPSILON=1/365,
  TAU=1/24,
  KAPPA=1.5/365,
  S=77.82,
  OMEGA=1/365,
  PHI=0.5,
  SUSJ.0=4000,EIJ.0=1000,ERJ.0=1000,INFJ.0=1000,
  RECJ.0=10000,SUSA.0=50000, EIA.0=100,
  ERA.0=1000,INFA.0=5000, RECA.0=50000)

sim <- simulate(sir,params=c(params),seed=3493885L,nsim=1,states=T,obs=F,as.data.frame=T) # double check seed in this

plot(sim$time,sim$SUSJ,type="l",ylim=c(0,max(sim$SUSJ)))
points(sim$time,sim$RECJ,col="green",type="l")
points(sim$time,sim$INFJ,col="red",type="l")

plot(sim$time,sim$SUSA,type="l",ylim=c(0,max(sim$SUSA)))
points(sim$time,sim$RECA,col="green",type="l")
points(sim$time,sim$INFA,col="red",type="l")

## to check can run for >1 simulation
# sims <- simulate(sir,params=c(params),seed=3493885L,nsim=10,states=T,obs=F,as.data.frame=T) # 100 simulations of 1 parameter set.

## for LHS parameter set
## Calling requisite libraries for parallel computing

library(foreach)
library(doSNOW)

#Setting up "parallel backend"

w<-makeCluster(3,type="SOCK") # makes the cluster, i.e. no of cores ABC = 8 cores, DEF = 12 see performance to see # of cores
registerDoSNOW(w) # 

#Checks that the number of workers is set up correctly.

getDoParWorkers()

## ~~~~~~~~~~ LHS SAMPLING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#############

# sims <- simulate(sir,params=c(params),seed=3493885L,nsim=10,states=T,obs=F,as.data.frame=T) # 100 simulations of 1 parameter set.

# NB not testing K - carrying cap, or rate of aging..
# no LHS needed


nonVarying = matrix(c(
  BETA=10,
  MU=0.000510492,
  DELTA=0.002312247,
  ALPHA=0.2,
  RHO=0.02,
  SIGMA=1/48,
  #K=1000000,
  EPSILON=1/365,
  TAU=1/24,
  KAPPA=1.5/365,
  S=77.82,
  OMEGA=1/365,
  PHI=0.5,
  SUSJ.0=4000,EIJ.0=1000,ERJ.0=1000,INFJ.0=1000,
  RECJ.0=10000,SUSA.0=50000, EIA.0=100,
  ERA.0=1000,INFA.0=5000, RECA.0=50000),
                    ncol=22,
                    nrow=40,
                    byrow=T) #binded with non-varying parameters


dimnames(nonVarying)[[2]]=c("BETA",
                            "MU",
                            "DELTA",
                            "ALPHA",
                            "RHO",
                            "SIGMA",
                            "EPSILON",
                            "TAU",
                            "KAPPA",
                            "S",
                            "OMEGA",
                            "PHI",
                            "SUSJ.0",
                            "EIJ.0","ERJ.0","INFJ.0",
                            "RECJ.0","SUSA.0", "EIA.0",
                            "ERA.0","INFA.0","RECA.0") # naming non-varying columns

## from other code
Kset=seq(from = 100, to=200000, by =5000)

fullParamSets = cbind(nonVarying,Kset) # full parameter set
#fullParamSets[,1] <- fullParamSets[,1]*fullParamSets[,29]
head(fullParamSets)
dim(fullParamSets)

dimnames(fullParamSets)[[2]]=c("BETA","MU","DELTA","ALPHA","RHO",
                               "SIGMA","EPSILON","TAU",
                               "KAPPA","S","OMEGA","PHI",
                               "SUSJ.0","EIJ.0",
                               "ERJ.0","INFJ.0","RECJ.0","SUSA.0",
                               "EIA.0","ERA.0","INFA.0","RECA.0",
                               "K")

# order for pomp/C model:  
BETA = fullParamSets[,1]
MU = fullParamSets[,2]
DELTA = fullParamSets[,3]
ALPHA = fullParamSets[,4]
RHO = fullParamSets[,5]
SIGMA = fullParamSets[,6]
K = fullParamSets[,23]
EPSILON = fullParamSets[,7]
TAU = fullParamSets[,8]
KAPPA = fullParamSets[,9]
S = fullParamSets[,10]
OMEGA = fullParamSets[,11]
PHI = fullParamSets[,12]
SUSJ.0 = fullParamSets[,13]
EIJ.0 = fullParamSets[,14]
ERJ.0 = fullParamSets[,15]
INFJ.0 = fullParamSets[,16]
RECJ.0 = fullParamSets[,17]
SUSA.0 = fullParamSets[,18]
EIA.0 = fullParamSets[,19]
ERA.0 = fullParamSets[,20]
INFA.0 = fullParamSets[,21]
RECA.0 = fullParamSets[,22]

paramset<-cbind(BETA,MU,DELTA,ALPHA,RHO,SIGMA,K,EPSILON,TAU,KAPPA,S,OMEGA,PHI,
                SUSJ.0,EIJ.0,ERJ.0,INFJ.0, RECJ.0, SUSA.0, EIA.0,ERA.0,INFA.0, RECA.0)
######################################################################################

## Calling requisite libraries for parallel computing

library(foreach)
library(doSNOW)

#Setting up "parallel backend"

w<-makeCluster(3,type="SOCK") # makes the cluster, i.e. no of cores ABC = 8 cores, DEF = 12 see performance to see # of cores
registerDoSNOW(w) # 

#Checks that the number of workers is set up correctly.

getDoParWorkers()

#######################################################3
results<-array(NA,dim=c(40,1,5))

# for one parameter set....
out1 <-simulate(sir,params=c(paramset[1,]),
                seed=1493885L,nsim=100,states=T,obs=F,as.data.frame=T) #

outres1 <- out1[seq(from=10951,to=1095100,by=10951),] # select last #s
N1 = array(0,c(100,5)) # same dimensions as No. runs * outputs I want
for (i in 1:100){ # each stochastic run
  N1[i,1]<-sum(outres1[i,1:10])
  N1[i,2]<-(sum(outres1[i,4],outres1[i,9])/sum(outres1[i,1:10]))*100 # prevalence; total
  N1[i,3]<-((outres1[i,10])/(sum(outres1[i,6:10])))*100 # adult seroprevalence; total
  N1[i,4]<-ifelse(sum(outres1[i,1:10])>0,1,0) # population extinct for each run
  N1[i,5]<-ifelse(sum(outres1[i,4],outres1[i,9])>0,1,0) # pathogen extinction for each run
}
N1[is.na(N1)]<- 0
## now average
M1 = array(0,c(1,5))
M1[1] = mean(N1[1:100,1]) # population size
M1[2] = mean(N1[1:100,2]) # prevalence
M1[3] = mean(N1[1:100,3]) # adult seroprevalence
M1[4] = mean(N1[1:100,4]) # adult seroprevalence
M1[5] = mean(N1[1:100,5]) # adult seroprevalence
rm(out1)
M1
results[1,,]<-M1
results[1,,]

##########################################################33
# for all parameter sets....
results<-array(NA,dim=c(40,1,5))

for (j in 1:length(paramset[,1])){
  out <-simulate(sir,params=c(paramset[j,]),
                 seed=1493885L,nsim=100,states=T,obs=F,as.data.frame=T) #
  outres <- out[seq(from=10951,to=1095100,by=10951),] # select last #s
  N = array(0,c(100,5)) # same dimensions as No. runs * outputs I want
  for (i in 1:100){ # each stochastic run
    N[i,1]<-sum(outres[i,1:10])
    N[i,2]<-(sum(outres[i,4],outres[i,9])/sum(outres[i,1:10]))*100 # prevalence; total
    N[i,3]<-((outres[i,10])/(sum(outres[i,6:10])))*100 # adult seroprevalence; total
    N[i,4]<-ifelse(sum(outres[i,1:10])>0,1,0) # population extinct for each run
    N[i,5]<-ifelse(sum(outres[i,4],outres[i,9])>0,1,0) # pathogen extinction for each run
  }
  N[is.na(N)]<- 0
  ## now average
  M = array(0,c(1,5))
  M[1] = mean(N[1:100,1]) # population size
  M[2] = mean(N[1:100,2]) # prevalence
  M[3] = mean(N[1:100,3]) # adult seroprevalence
  M[4] = mean(N[1:100,4]) # mean pop extinction
  M[5] = mean(N[1:100,5]) # mean path extinction
  rm(out)
  results[j,,]<-M
}
#
#########################################################33

w<-makeCluster(1,type="SOCK") # return to one core

############################################################3
## need matrix of results...
X<-aperm(results,c(1,2,3))
dim(X)<-c(40,5)
head(X)
tail(X)

################################################################################
# below for K

plot(X[,1],X[,2])
# plot

par(mar=c(5, 6, 4, 4) + 0.1)

plot(X[,1],X[,2],pch=16,
     ylab="Mean prevalence",xlab="Population size",
     col="grey25", cex.lab=1.2)

plot(X[,1],X[,3],pch=16,
     ylab="Mean seroprevalence",xlab="Population size",
     col="grey25", cex.lab=1.2)

plot(X[,1],X[,5],pch=16,
     ylab="P[persist]",xlab="Population size",
     col="grey25", cex.lab=1.2)
