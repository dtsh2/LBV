## LBV SEIR model with Frequency-dept transmission
## K
## clean up

rm(list=ls())

# instal r developer toolbox first (Rtools from Cran R)

## pomp test run lbv
getwd()
setwd("~/Cambridge/CSU 2013/LBV Model/lbvmodels/new_models_sept")

library(pomp)

#Compiling C code and loading the dll
system("R CMD SHLIB lbvseirNoSeasFreq.c")

dyn.load("lbvseirNoSeasFreq.dll")

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
  statenames=c("SUSJ","MDAJ", "SUSJM","EIJ","ERJ","INFJ", "RECJ", "SUSA", "EIA","ERA","INFA", "RECA"),
  ## the order of the parameters assumed in the native routines:
  paramnames=c("BETA","MU","DELTA","ALPHA","RHO","SIGMA","K","EPSILON","TAU","PSI","KAPPA","S","OMEGA","PHI","GAMMA",
               "SUSJ.0","MDAJ.0", "SUSJM.0","EIJ.0","ERJ.0","INFJ.0", "RECJ.0", "SUSA.0", "EIA.0","ERA.0","INFA.0", "RECA.0"),
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
  PSI=0.01,
  KAPPA=1.5/365,
  S=77.82,
  OMEGA=1/365,
  PHI=0.5,
  GAMMA=0.003225806, # check data
  SUSJ.0=4000,MDAJ.0=4000, SUSJM.0=1000,EIJ.0=1000,ERJ.0=1000,INFJ.0=1000,
  RECJ.0=10000,SUSA.0=50000, EIA.0=100,
  ERA.0=1000,INFA.0=5000, RECA.0=50000)

sim <- simulate(sir,params=c(params),seed=3493885L,nsim=1,states=T,obs=F,as.data.frame=T) # double check seed in this

plot(sim$time,sim$SUSJ,type="l")
points(sim$time,sim$RECJ,col="green",type="l")
points(sim$time,sim$MDA,col="brown",type="l")
points(sim$time,sim$INFJ,col="red",type="l")

plot(sim$time,sim$SUSA,type="l")
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
# not testing rate of aging or # birth pulses
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
  PSI=55,
  KAPPA=1.5/365,
  S=77.82,
  OMEGA=1/365,
  PHI=0.5,
  SUSJ.0=4000,MDAJ.0=4000, SUSJM.0=1000,EIJ.0=1000,ERJ.0=1000,INFJ.0=1000,
  RECJ.0=10000,SUSA.0=50000, EIA.0=100,
  ERA.0=1000,INFA.0=5000, RECA.0=50000),
                    ncol=25,
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
                            "PSI",
                            "KAPPA",
                            "S",
                            "OMEGA",
                            "PHI",
                            #"GAMMA", # check data
                            "SUSJ.0",
                            "MDAJ.0",
                            "SUSJM.0","EIJ.0","ERJ.0","INFJ.0",
                            "RECJ.0","SUSA.0", "EIA.0",
                            "ERA.0","INFA.0","RECA.0") # naming non-varying columns

gamma.day <- (365-nonVarying[,9])
gamma <-1/gamma.day
nonVarying[,9] <-1/nonVarying[,9]
## from other code
Kset=seq(from = 100, to=200000, by =5000)

fullParamSets = cbind(nonVarying,gamma,Kset) # full parameter set
#fullParamSets[,1] <- fullParamSets[,1]*fullParamSets[,29]
head(fullParamSets)
dim(fullParamSets)

dimnames(fullParamSets)[[2]]=c("BETA","MU","DELTA","ALPHA","RHO",
                               "SIGMA","EPSILON","TAU","PSI",
                               "KAPPA","S","OMEGA","PHI",
                               "SUSJ.0","MDAJ.0","SUSJM.0","EIJ.0",
                               "ERJ.0","INFJ.0","RECJ.0","SUSA.0",
                               "EIA.0","ERA.0","INFA.0","RECA.0",
                               "GAMMA","K")

# order for pomp/C model:  
BETA = fullParamSets[,1]
MU = fullParamSets[,2]
DELTA = fullParamSets[,3]
ALPHA = fullParamSets[,4]
RHO = fullParamSets[,5]
SIGMA = fullParamSets[,6]
K = fullParamSets[,27]
EPSILON = fullParamSets[,7]
TAU = fullParamSets[,8]
PSI = fullParamSets[,9]
KAPPA = fullParamSets[,10]
S = fullParamSets[,11]
OMEGA = fullParamSets[,12]
PHI = fullParamSets[,13]
GAMMA = fullParamSets[,26]
SUSJ.0 = fullParamSets[,14]
MDAJ.0 = fullParamSets[,15]
SUSJM.0 = fullParamSets[,16]
EIJ.0 = fullParamSets[,17]
ERJ.0 = fullParamSets[,18]
INFJ.0 = fullParamSets[,19]
RECJ.0 = fullParamSets[,20]
SUSA.0 = fullParamSets[,21]
EIA.0 = fullParamSets[,22]
ERA.0 = fullParamSets[,23]
INFA.0 = fullParamSets[,24]
RECA.0 = fullParamSets[,25]

paramset<-cbind(BETA,MU,DELTA,ALPHA,RHO,SIGMA,K,EPSILON,TAU,PSI,KAPPA,S,OMEGA,PHI,GAMMA,
                SUSJ.0, MDAJ.0, SUSJM.0,EIJ.0,ERJ.0,INFJ.0, RECJ.0, SUSA.0, EIA.0,ERA.0,INFA.0, RECA.0)
######################################################################################
results<-array(NA,dim=c(40,1,5))

# for one parameter set....
out1 <-simulate(sir,params=c(paramset[1,]),
                seed=1493885L,nsim=100,states=T,obs=F,as.data.frame=T) #

outres1 <- out1[seq(from=10951,to=1095100,by=10951),] # select last #s
N1 = array(0,c(100,5)) # same dimensions as No. runs * outputs I want
for (i in 1:100){ # each stochastic run
  N1[i,1]<-sum(outres1[i,1:12])
  N1[i,2]<-(sum(outres1[i,6],outres1[i,11])/sum(outres1[i,1:12]))*100 # prevalence; total
  N1[i,3]<-((outres1[i,12])/(sum(outres1[i,8:12])))*100 # adult seroprevalence; total
  N1[i,4]<-ifelse(sum(outres1[i,1:12])>0,1,0) # population extinct for each run
  N1[i,5]<-ifelse(sum(outres1[i,6],outres1[i,11])>0,1,0) # pathogen extinction for each run
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
  N[i,1]<-sum(outres[i,1:12])
  N[i,2]<-(sum(outres[i,6],outres[i,11])/sum(outres[i,1:12]))*100 # prevalence; total
  N[i,3]<-((outres[i,12])/(sum(outres[i,8:12])))*100 # adult seroprevalence; total
  N[i,4]<-ifelse(sum(outres[i,1:12])>0,1,0) # population extinct for each run
  N[i,5]<-ifelse(sum(outres[i,6],outres[i,11])>0,1,0) # pathogen extinction for each run
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
