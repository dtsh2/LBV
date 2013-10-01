## LBV SEIRS model with Frequency-dept transmission
## rho/beta simulations
## clean up

rm(list=ls())

# instal r developer toolbox first (Rtools from Cran R)

## pomp test run lbv
getwd()
setwd("~/Cambridge/CSU 2013/LBV Model/lbvmodels/new_models_sept")

library(pomp)

#Compiling C code and loading the dll
system("R CMD SHLIB lbvseirsFreq.c")

dyn.load("lbvseirsFreq.dll")

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
  paramnames=c("BETA","MU","DELTA","ALPHA","RHO","SIGMA","K","EPSILON","TAU","PSI","KAPPA","S","OMEGA","PHI","GAMMA","LAMBDA","ZETA",
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
  LAMBDA=0.001, # check data
  ZETA=0.5,
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

nonVarying = matrix(c(
  #  BETA=0.0001,
  MU=0.000510492,
  DELTA=0.002312247,
  ALPHA=0.2,
  #  RHO=0.05,
  SIGMA=1/48,
  K=1000000,
  EPSILON=1/365,
  TAU=1/24,
  PSI=55,
  KAPPA=1.5/365,
  S=77.82,
  OMEGA=1/365,
  PHI=0.5,
  LAMBDA=0.001, # check data
  ZETA=0.5,
  SUSJ.0=4000,MDAJ.0=4000, SUSJM.0=1000,EIJ.0=1000,ERJ.0=1000,INFJ.0=1000,
  RECJ.0=10000,SUSA.0=50000, EIA.0=100,
  ERA.0=1000,INFA.0=5000, RECA.0=50000),
                    ncol=26,
                    nrow=100,
                    byrow=T) #binded with non-varying parameters


dimnames(nonVarying)[[2]]=c(#"BETA",
  "MU",
  "DELTA",
  "ALPHA",
  #"RHO",
  "SIGMA",
  "K",
  "EPSILON",
  "TAU",
  "PSI",
  "KAPPA",
  "S",
  "OMEGA",
  "PHI",
  #"GAMMA", # check data
  "LAMBDA", # check data
  "ZETA",
  "SUSJ.0",
  "MDAJ.0",
  "SUSJM.0","EIJ.0","ERJ.0","INFJ.0",
  "RECJ.0","SUSA.0", "EIA.0",
  "ERA.0","INFA.0","RECA.0") # naming non-varying columns

gamma.day <- (365-nonVarying[,8])
gamma <-1/gamma.day
nonVarying[,8] <-1/nonVarying[,8]
## from other code

library(lhs)

nspaces=100 ## how many bins/intervals, 100 needed for PRCC

hypercube=optimumLHS(n=nspaces, k=2) ## function with N columns
dimnames(hypercube)[[2]]=c("beta","rho")  # named columns with parameter names

mins = c(       	    ## set mins for each parameters-exclude those not to be varied if any
  beta= 1,		     # transmission 
  rho= 0.0001				# probability that exposure/infection will lead to infection & infectiousness (and dead)
)       			

maxs = c( 				    ## set mins for each parameters-exclude those not to be varied if any
  beta= 10,           # transmission
  rho= 0.1)                      

diffs=maxs-mins ## range of each variable

hypercubeadj = hypercube # create copy of hypercube samples to modify, hypercube adjusted; i.e. new matrix
for (i in 1:ncol(hypercube)){
  hypercubeadj[,i]=hypercube[,i]*diffs[i]+mins[i] # scale samples to difference and add minimum
}

head(hypercubeadj)
dim(hypercubeadj)

dimnames(hypercubeadj)[[2]]=c("beta","rho")

#beta=seq(from = 0, to=1/190, by =5.263158e-05)
#rho=seq(from = 0.0001,to=0.4, by=0.00399999)

#fullParamSets = cbind(nonVarying,gamma,beta,rho) # full parameter set

fullParamSets = cbind(nonVarying,gamma,hypercubeadj)
head(fullParamSets)
dim(fullParamSets)

dimnames(fullParamSets)[[2]]=c("MU","DELTA","ALPHA",
                               "SIGMA","K","EPSILON","TAU","PSI",
                               "KAPPA","S","OMEGA","PHI","LAMBDA","ZETA",
                               "SUSJ.0","MDAJ.0","SUSJM.0","EIJ.0",
                               "ERJ.0","INFJ.0","RECJ.0","SUSA.0",
                               "EIA.0","ERA.0","INFA.0","RECA.0",
                               "GAMMA","BETA","RHO")

# order for pomp/C model:  
BETA = fullParamSets[,28]
MU = fullParamSets[,1]
DELTA = fullParamSets[,2]
ALPHA = fullParamSets[,3]
RHO = fullParamSets[,29]
SIGMA = fullParamSets[,4]
K = fullParamSets[,5]
EPSILON = fullParamSets[,6]
TAU = fullParamSets[,7]
PSI = fullParamSets[,8]
KAPPA = fullParamSets[,9]
S = fullParamSets[,10]
OMEGA = fullParamSets[,11]
PHI = fullParamSets[,12]
GAMMA = fullParamSets[,27]
LAMBDA = fullParamSets[,13]
ZETA = fullParamSets[,14]
SUSJ.0 = fullParamSets[,15]
MDAJ.0 = fullParamSets[,16]
SUSJM.0 = fullParamSets[,17]
EIJ.0 = fullParamSets[,18]
ERJ.0 = fullParamSets[,19]
INFJ.0 = fullParamSets[,20]
RECJ.0 = fullParamSets[,21]
SUSA.0 = fullParamSets[,22]
EIA.0 = fullParamSets[,23]
ERA.0 = fullParamSets[,24]
INFA.0 = fullParamSets[,25]
RECA.0 = fullParamSets[,26]

paramset<-cbind(BETA,MU,DELTA,ALPHA,RHO,SIGMA,K,EPSILON,TAU,PSI,KAPPA,S,OMEGA,PHI,GAMMA,LAMBDA,ZETA,
                SUSJ.0, MDAJ.0, SUSJM.0,EIJ.0,ERJ.0,INFJ.0, RECJ.0, SUSA.0, EIA.0,ERA.0,INFA.0, RECA.0)
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
results<-array(NA,dim=c(100,1,5))

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
results<-array(NA,dim=c(100,1,5))

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
dim(X)<-c(100,5)
head(X)
tail(X)

################################################################################
library(lattice)
library(akima)
#install.packages("scatterplot3d", dependencies = TRUE)
library(scatterplot3d)

par(mfrow=c(1,1))
scatterplot3d(x = paramset[,1], y = paramset[,5], z = X[,2],
              angle=35,pch=16, highlight.3d=T,main="",
              xlab=expression(beta),ylab=(expression(rho)),
              zlab=expression("Prevalence"))

scatterplot3d(x = paramset[,1], y = paramset[,5], z = X[,3],
              angle=45,pch=16, highlight.3d=T,main="",
              xlab=expression(beta),ylab=(expression(rho)),
              zlab=expression("Adult seroprevalence"))

scatterplot3d(x = paramset[,1], y = paramset[,5], z = X[,5],
              angle=35,pch=16, highlight.3d=T,main="",
              xlab=expression(beta),ylab=(expression(rho)),
              zlab=expression("P[LBV extinct]"))

#install.packages("rgl", dependencies = TRUE)
library(rgl)
#plot3d(x = paramset[,1], y = paramset[,5], z = a[,2])

par(mai=c(1.5,1.5,1,1.5))
zzp<- interp(paramset[,1],paramset[,5],X[,2], duplicate=T)
#image(zz,col=topo.colors(12),main="prevalence")
#contour(zz,add=T,cex=2)

filled.contour(zzp, col=topo.colors(30),
               ylab=expression(rho),xlab=expression(beta),
               cex.lab=1.2)
mtext("Prevalence", side=4,line=3,cex=1.2)

zzs<- interp(paramset[,1],paramset[,5],X[,3], duplicate=T)
#image(zz,col=topo.colors(12),main="seroprevalence")
#contour(zz,add=T,cex=2)

filled.contour(zzs, col=topo.colors(30),
               ylab=expression(rho),xlab=expression(beta),
               cex.lab=1.2)
mtext("Adult seroprevalence", side=4,line=3,cex=1.2)

## 

## reconvert to beta
betaorig <-fullParamSets[,28]/fullParamSets[,5]
zzp<- interp(log(betaorig),paramset[,5],X[,2], duplicate=T)
filled.contour(zzp, col=topo.colors(24),
               ylab=expression(rho),xlab= #expression(beta),
                 expression(paste(plain(Log)," (", beta,")")),
               cex.lab=1.2)
mtext("Prevalence", side=4,line=1,cex=1.2)

zzs<- interp(betaorig,paramset[,5],X[,3], duplicate=T)
#image(zz,col=topo.colors(12),main="seroprevalence")
#contour(zz,add=T,cex=2)

filled.contour(zzs, col=topo.colors(30),
               ylab=expression(rho),xlab= expression(beta),#expression(paste(plain(Log)," (", beta,")")),
               cex.lab=1.2)
mtext("Adult seroprevalence", side=4,line=1,cex=1.2)

plot(X[,2],X[,3])
plot(X[,1],X[,2])
plot(X[,1],X[,3])

## with default params

params <- c(
  BETA= 10,
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
  LAMBDA=0.001, # check data
  ZETA=0.5,
  SUSJ.0=4000,MDAJ.0=4000, SUSJM.0=1000,EIJ.0=1000,ERJ.0=1000,INFJ.0=1000,
  RECJ.0=10000,SUSA.0=50000, EIA.0=100,
  ERA.0=1000,INFA.0=5000, RECA.0=50000)

sim <- simulate(sir,params=c(params),seed=3493885L,nsim=1,states=T,obs=F,as.data.frame=T) # double check seed in this
#
#
#par(mfrow=c(2,1))
plot(sim$SUSJM[6800:7300], main ="Juvenile",ylim=c(0,200000),col="blue",
     xlab="Time (Days)",ylab="Numbers",type="l")
lines(sim$SUSJ[6800:7300],col="blue")
lines(sim$MDAJ[6800:7300],col="brown")
lines(sim$EIJ[6800:7300],col="yellow")
lines(sim$ERJ[6800:7300],col="grey")
lines(sim$INFJ[6800:7300],col="red")
lines(sim$RECJ[6800:7300],col="green")

plot(sim$SUSA[6800:7300], main ="Adults",ylim=c(0,500000),col="blue",
     xlab="Time (Days)",ylab="Numbers",type="l")
lines(sim$EIA[6800:7300],col="yellow")
lines(sim$ERA[6800:7300],col="grey")
lines(sim$INFA[6800:7300],col="red")
lines(sim$RECA[6800:7300],col="green")

## add susj
par(mfrow=c(1,1))
par(mai=c(2,2,1,1))

susjt<-rowSums(sim[,c(1,3:5)])

susat<-rowSums(sim[,c(8:10)])

new.sim <-cbind(sim,susjt,susat)
# plot
plot(new.sim$susjt[6000:7300], main ="",ylim=c(0,300000),col="blue",
     xlab="Time (Days)",ylab="",type="l",lwd=2,cex.lab=1.2)
mtext("Numbers", side=2,las=0,at=1.5e+5,line=5,cex=1.2)
lines(new.sim$MDAJ[6000:7300],col="brown",lwd=2)
#lines(new.sim$EIJ[6000:7300],col="yellow",lwd=2)
#lines(new.sim$ERJ[6000:7300],col="grey",lwd=2)
lines(new.sim$INFJ[6000:7300],col="red",lwd=2)
lines(new.sim$RECJ[6000:7300],col="green",lwd=2)

plot(new.sim$susat[6000:7300], main ="",ylim=c(0,600000),col="blue",
     xlab="Time (Days)",ylab="",type="l",lwd=2,cex.lab=1.2)
mtext("Numbers", side=2,las=0,at=250000,line=5,cex=1.2)
#lines(new.sim$EIA[6000:7300],col="yellow",lwd=2)
#lines(new.sim$ERA[6000:7300],col="grey",lwd=2)
lines(new.sim$INFA[6000:7300],col="red",lwd=2)
lines(new.sim$RECA[6000:7300],col="green",lwd=2)

#
dim(new.sim)

par(mfrow=c(1,1))
plot(new.sim$INFA[6000:7300], main ="",ylim=c(0,100),col="red",
     xlab="Time (Days)",ylab="",type="l",lwd=2,cex.lab=1.2)
mtext("Numbers", side=2,las=0,at=50,line=4,cex=1.2)
lines(new.sim$INFJ[6000:7300],col="orange",lwd=2)
##############3
