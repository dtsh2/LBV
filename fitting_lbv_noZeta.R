# libraries
rm(list=ls())

require(ggplot2)
require(plyr)
require(deSolve)
require(reshape)
require(bbmle)

setwd("~/Cambridge/CSU 2013/LBV model/lbvmodels/new_models_sept/deterministic")
## data

lbv<-read.csv("lbv_data.csv")


juv_sp<-lbv$DRECJ/(lbv$DRECJ+lbv$DSUSJ)
ad_sp<-lbv$DRECA/(lbv$DRECA+lbv$DSUSA)
lbv.new<-cbind(lbv,juv_sp,ad_sp)
dim(lbv.new)
lbv.sp<-cbind(lbv.new[,c(1,6,7)])
lbv.sp
###################3

## sir-model-defn
## freq dept

lbv.model <- function (t, x, params) {
  ## first extract the state variables
  SUSJ <- x[1]
  MDA  <- x[2]
  SUSJM<- x[3]
  SUSA <- x[4]
  EIJ  <- x[5]
  EIA  <- x[6]
  ERJ  <- x[7]
  ERA  <- x[8]
  INFJ <- x[9]
  INFA <- x[10]
  RECJ <- x[11]
  RECA <- x[12]
  N    <- x[13]
  SPJ  <- x[14]
  SPA  <- x[15]
  N    <- unname((SUSJ+MDA+SUSJM+SUSA+EIJ+EIA+ERJ+ERA+INFJ+INFA+RECJ+RECA))
  SPJ  <- unname((RECJ)/(SUSJ+SUSJM+EIJ+ERJ))
  SPA  <- unname((RECA)/(SUSA+EIA+ERA+RECA))
  ## now extract the parameters
  beta <- params["beta"]
  gamma <- params["gamma"]
  delta <- params["delta"]
  mu    <- params["mu"]
  alpha <- params["alpha"]
  rho   <- params["rho"]
  epsilon <- params["epsilon"]
  sigma <- params["sigma"]
  tau   <- params["tau"]
  psi   <- params["psi"]
  phi   <- params["phi"]
  omega <- params["omega"]
  kappa <- params["kappa"]
  s     <- params["s"]
  K     <- params["K"]
  lambda<- params["lambda"]
  zeta  <- params["zeta"]
  ## now code the model equations
  dSUSJ.dt <-(kappa*(1/sqrt((1/s)*pi)*exp(-((cos(pi*omega*t-phi))^2)/(1/s))))*(SUSA+ERA)-beta*SUSJ*((INFJ+INFA)/N)-delta*SUSJ*(N/K)-epsilon*SUSJ
  dMDA.dt  <- (kappa*(1/sqrt((1/s)*pi)*exp(-((cos(pi*omega*t-phi))^2)/(1/s))))*(RECA)- delta*MDA*(N/K)- psi*MDA
  dSUSJM.dt<- psi*MDA - delta*SUSJM*(N/K) - beta*SUSJM*((INFJ+INFA)/N) - gamma*SUSJM
  dSUSA.dt <- epsilon*SUSJ + gamma*SUSJM- mu*SUSA*(N/K)+ lambda*RECA- beta*SUSA*((INFJ+INFA)/N)
  dEIJ.dt  <- rho*beta*SUSJ*((INFJ+INFA)/N) + rho*beta*SUSJM*((INFJ+INFA)/N)- delta*EIJ*(N/K)- sigma*EIJ
  dEIA.dt  <- rho*beta*SUSA*((INFJ+INFA)/N) - mu*EIA*(N/K)- sigma*EIA
  dERJ.dt  <- (1-rho)*beta*SUSJ*((INFJ+INFA)/N) +(1-rho)*beta*SUSJM*((INFJ+INFA)/N)- epsilon*ERJ- delta*ERJ*(N/K)- tau*ERJ
  dERA.dt  <- (1-rho)*beta*SUSA*((INFJ+INFA)/N)- tau*ERA - mu*ERA*(N/K)+epsilon*ERJ
  dINFJ.dt <- sigma*EIJ - alpha*INFJ
  dINFA.dt <- sigma*EIA - alpha*INFA
  dRECJ.dt <- tau*ERJ - delta*RECJ*(N/K)- epsilon*RECJ
  dRECA.dt <- tau*ERA - mu*RECA*(N/K)+ epsilon*RECJ - lambda*RECA
  ## combine results into a single vector

  dxdt<-c(dSUSJ.dt,
          dMDA.dt,
          dSUSJM.dt,
          dSUSA.dt,
          dEIJ.dt,
          dEIA.dt,
          dERJ.dt,
          dERA.dt,
          dINFJ.dt,
          dINFA.dt,
          dRECJ.dt,
          dRECA.dt) 
    ## return result as a list
  list(dxdt,N,SPJ,SPA)
  }


times <- seq(from=0,to=15000,by=1) ## returns a sequence
#
# init cond
xstart <- c(SUSJ=10000,
            MDA=10,
            SUSJM=1000,
            SUSA=800000,
            EIJ=10,
            EIA=10,
            ERJ=10,
            ERA=10,
            INFJ=10,
            INFA=10,
            RECJ=1000,
            RECA=200000) # initial cond & order passing to y below

## sir-model-likelihood

lbv.nll <- function (beta,
                     mu,
                     delta,
                     alpha,
                     rho,
                     sigma,
                     K,
                     epsilon,
                     tau,
                     psi,
                     kappa,
                     s,
                     omega,
                     phi,
                     gamma, # check data
                     lambda, # check data
                     zeta,
                     p, k) {
  times <- c(lbv.sp$cumulative_time[1]-1,lbv.sp$cumulative_time)
  ode.params <- c(
    beta=beta,#5,
    mu=mu,#0.000510492,
    delta=delta,#0.002312247,
    alpha=alpha,#0.2,
    rho=rho,#0.05,
    sigma=sigma,#0.02083333,
    K=K,#1000000,
    epsilon=epsilon,#0.002739726,
    tau=tau,#0.04166667,
    psi=psi,#0.01,
    kappa=kappa,#0.004109589,
    s=s,#77.82,
    omega=omega,#0.002739726,
    phi=phi,#0.5,
    gamma=gamma,#0.02092154, # check data
    lambda=lambda,#0.001, # check data
    zeta=zeta#0.5
    )
  xstart <- xstart
  out <- as.data.frame(ode(
             func=lbv.model,
             y=xstart,
             times=times,
             parms=ode.params
             ))
  colnames(out)[14:16]<-c("N","SpJ","SpA") # 
  ## 'out' is a matrix
  ll1 <- dnorm(x=lbv.sp$ad_sp,as.numeric(as.matrix(out["SpA"])),log=TRUE)
  ll2 <- dnorm(x=lbv.sp$juv_sp,as.numeric(as.matrix(out["SpJ"])),log=TRUE)
  ll<-ll1+ll2
  -sum(ll)
}


#################################################
# beta, rho -loglik

## but major identifiability issues

nll <- function (par) {
  lbv.nll(beta=par[1],
          mu=0.000510492,
          delta=0.002312247,
          alpha=0.2,
          rho=par[2],
          sigma=0.02083333,
          K=1000000,
          epsilon=0.002739726,
          tau=0.04166667,
          psi=0.01,
          kappa=1.5/365,
          s=14.35,
          omega=0.002739726,
          phi=0.0,
          gamma=0.02092154, # check data
          lambda=0.001,
          zeta=0.0#,
  ) # 
}

# using optim
fit <- optim(fn=nll,par=c(18,0.012),method="Nelder-Mead")
fit

## values:
#beta/rho$par
#[1] 1.546018e+02 1.945056e-03
#
#LL$value
#[1] 24.0809

###########################################
## mle2-beta, rho
## Nelder-Mead

fit <- mle2(lbv.nll,
            start=list(beta=18,rho=0.01),
            method="Nelder-Mead",
            fixed=list(mu=0.000510492,
                       delta=0.002312247,
                       alpha=0.2,
                       sigma=0.02083333,
                       K=1000000,
                       epsilon=0.002739726,
                       tau=0.04166667,
                       psi=0.01,
                       kappa=1.5/365,
                       s=14.35,
                       omega=0.002739726,
                       phi=0.0,
                       gamma=0.02092154, # check data
                       lambda=0.001,
                       zeta=0.0))
coef(fit)

# pfit <- profile(fit,c("beta","rho"))
# plot(pfit)
# confint(pfit)

## mle2-beta, rho
## L-BFGS-B

fit <- mle2(lbv.nll,
            start=list(beta=18,rho=0.015),
            method="L-BFGS-B",
            lower=c(0,0),
            upper=c(30,1),
            fixed=list(mu=0.000510492,
                       delta=0.002312247,
                       alpha=0.2,
                       sigma=0.02083333,
                       K=1000000,
                       epsilon=0.002739726,
                       tau=0.04166667,
                       psi=0.01,
                       kappa=1.5/365,
                       s=14.35,
                       omega=0.002739726,
                       phi=0.0,
                       gamma=0.02092154, # check data
                       lambda=0.001,
                       zeta=0.0)  
            )
coef(fit)

# pfit <- profile(fit,c("beta","rho"))
# plot(pfit)
# confint(pfit)

##########################################
## contour plots
######################################################################
## grid of {beta, rho} combinations. # time!!

beta = seq(from=0,to=30,by=0.05)  # range of beta
rho = seq(from=0.001,to=0.04, by=0.005) # range of rho
pars.grid <- expand.grid(beta,rho)
dim(pars.grid)

# for all parameter sets....
results<-array(NA,dim=c(length(pars.grid[,1]),1))

for (j in 1:length(pars.grid[,1])){
  results[j,]<-nll(par=c(pars.grid[j,1],pars.grid[j,2]))
}

library(akima) # sometimes problems with interp when data was regular
library(tgp) # for the interp.loess function for regular data
#library(lattice)

zz <- interp.loess( pars.grid[,1],
                    pars.grid[,2],
                    results[,1],
                    duplicate=T)
zz <- interp( pars.grid[,1],
              pars.grid[,2],
              results[,1],
              duplicate=T)
min(results)
image(zz,xlab="beta",ylab="rho")
contour(zz,add=T,levels=c(#24.1,
                          24.5,
                          25,
                          30,
                          #40,
                          50,
                          #60,70,
                          #80,90,
                          100,
                          #110,
                          #120,130,
                          150))#seq(24.2,40,1))

points(x=18,y=0.0172,pch=16,col="pink")
text(x=15,y=0.015,"MLE",col="pink")
###################################################
##
## for each par value alone 
##
# beta-loglik

nll <- function (par) {
  lbv.nll(beta=par[1]
          ,
          mu=0.000510492,
          delta=0.002312247,
          alpha=0.2,
          rho=0.0172,
          sigma=0.02083333,
          K=1000000,
          epsilon=0.002739726,
          tau=0.04166667,
          psi=0.01,
          kappa=1.5/365,
          s=14.35,
          omega=0.002739726,
          phi=0.0,
          gamma=0.02092154, # check data
          lambda=0.001,
          zeta=0.0#,
  ) # 
}

betacurve <- data.frame(beta=seq(1,30,length=100))
within(betacurve,nll <- sapply(c(beta),nll)) -> betacurve

head(betacurve)

ggplot(data=betacurve,mapping=aes(x=beta,y=nll))+geom_line()

### optim-beta
#par0<-c(beta=18.9)
#
#fit <- optim(fn=nll,par=par0,method="Brent",lower=0.01,upper=30)
#fit
#
##############################3## 
#
#crit.lr <- pchisq(q=0.05,df=1,lower.tail=FALSE)
#betacurve <- data.frame(beta=seq(1,30,length=100))
#betacurve <- within(betacurve,nll <- sapply(beta,nll))
#ggplot(data=betacurve,mapping=aes(x=beta,y=nll))+geom_line()+
#  ylim(fit$value+c(0,30))+
#  geom_vline(xintercept=fit$par,color='red')+
#  geom_hline(yintercept=fit$value+crit.lr,color='blue')

## mle for rho

rhonll <- function (par) {
  lbv.nll(beta=18.2,
          mu=0.000510492,
          delta=0.002312247,
          alpha=0.2,
          rho=par[1],
          sigma=0.02083333,
          K=1000000,
          epsilon=0.002739726,
          tau=0.04166667,
          psi=0.01,
          kappa=1.5/365,
          s=14.35,
          omega=0.002739726,
          phi=0.0,
          gamma=0.02092154, # check data
          lambda=0.001,
          zeta=0.0#,
  ) # 
}

rhocurve <- data.frame(rho=seq(0.001,0.2,length=100))
within(rhocurve,rhonll <- sapply(c(rho),rhonll)) -> rhocurve

head(rhocurve)

ggplot(data=rhocurve,mapping=aes(x=rho,y=rhonll))+geom_line()+ylab("nll")

## optim-rho
#par0<-c(rho=0.01)
#
#fit <- optim(fn=nll,par=par0,method="Brent",lower=0.001,upper=0.1)
#fit
#
#crit.lr <- pchisq(q=0.05,df=1,lower.tail=FALSE)
#rhocurve <- data.frame(rho=seq(0.001,0.2,length=100))
#within(rhocurve,rhonll <- sapply(c(rho),rhonll)) -> rhocurve
#ggplot(data=rhocurve,mapping=aes(x=rho,y=rhonll))+geom_line()+
#  ylim(fit$value+c(0.001,0.2))+
#  geom_vline(xintercept=fit$par,color='red')+
#  geom_hline(yintercept=fit$value+crit.lr,color='blue')

##### plot with new values and compare results and data
#ggplot(lbv, aes(cumulative_time)) + 
#  geom_line(aes(y = DSUSJ, colour = "dsusj")) + 
#  geom_line(aes(y = DSUSA, colour = "dsusa")) +
#  geom_line(aes(y = DRECJ, colour = "drecj")) + 
#  geom_line(aes(y = DRECA, colour = "dreca"))
#
#lbv_long <- melt(lbv, id="cumulative_time")  # convert to long format
#
#ggplot(data=lbv_long,
#       aes(x=cumulative_time, y=value, colour=variable)) +
#  geom_line()

juv_sp<-lbv$DRECJ/(lbv$DRECJ+lbv$DSUSJ)
ad_sp<-lbv$DRECA/(lbv$DRECA+lbv$DSUSA)
lbv.new<-cbind(lbv,juv_sp,ad_sp)
dim(lbv.new)
lbv.sp<-cbind(lbv.new[,c(1,6,7)])
lbv.sp

lbv_sp <- melt(lbv.sp, id="cumulative_time")  # convert to long format

ggplot(data=lbv_sp,
       aes(x=cumulative_time, y=value, colour=variable)) +
  geom_line()

## binomial CI for SP
## each sample get CIs
for(ii in 1:12)
{
  lbv.new$SpA.ci.l[ii] <- binom.test(lbv.new$DRECA[ii], sum(lbv.new$DRECA[ii],lbv.new$DSUSA[ii]))$conf.int[1]
  lbv.new$SpA.ci.u[ii] <- binom.test(lbv.new$DRECA[ii], sum(lbv.new$DRECA[ii],lbv.new$DSUSA[ii]))$conf.int[2]    
}

for(ii in 1:12)
{
  lbv.new$SpJ.ci.l[ii] <- binom.test(lbv.new$DRECJ[ii], sum(lbv.new$DRECJ[ii],lbv.new$DSUSJ[ii]))$conf.int[1]
  lbv.new$SpJ.ci.u[ii] <- binom.test(lbv.new$DRECJ[ii], sum(lbv.new$DRECJ[ii],lbv.new$DSUSJ[ii]))$conf.int[2]    
}

## Plot data vs the true underlying epidemic.
plot(lbv.new$cumulative_time, lbv.new$ad_sp, type="l", col="red", bty = "n",ylim=c(0,1),
     #ylim = c(0, max(lbv.new$SpA.ci.u)), 
     xlab = "days", ylab = "seroprevalence")
## Add data
points(lbv.new$cumulative_time, lbv.new$ad_sp, col = "black", pch = 19)
## Add CIs
arrows(lbv.new$cumulative_time, lbv.new$SpA.ci.l, lbv.new$cumulative_time, lbv.new$SpA.ci.u, length = .01,
       angle = 90, code = 3)
lines(lbv.new$cumulative_time, lbv.new$juv_sp, type="l", col="blue")
## Add data
points(lbv.new$cumulative_time, lbv.new$juv_sp, col = "black", pch = 19)
## Add CIs
arrows(lbv.new$cumulative_time, lbv.new$SpJ.ci.l, lbv.new$cumulative_time, lbv.new$SpJ.ci.u, length = .01,
       angle = 90, code = 3)
## Add legends
legend("topright", c("Adult", "Juvenile"),
       col=c("red","blue"),lty=1,bty="n")

##########################
## must add other data here
## closed-params

params <- c(
  beta=18,#, mle
  mu=0.000510492,
  delta=0.002312247,
  alpha=0.2,
  rho=0.0172,#mle
  sigma=0.02083333,
  K=1000000,
  epsilon=0.002739726,
  tau=0.04166667,
  psi=0.01,
  kappa=1.5/365,
  s=14.35,
  omega=0.002739726,
  phi=0.0,
  gamma=0.02092154, # check data
  lambda=0.001, # check data
  zeta=0.0)

# set-times-ics

times <- seq(from=0,to=15000,by=1) ## returns a sequence

##  to check works
out <- as.data.frame(
  ode(
    func=lbv.model,
    y=xstart,
    times=times,
    parms=params
  )
)
head(out)
tail(out)
colnames(out)[14:16]<-c("N","SpJ","SpA") # 
head(out[,14:16])
head(out)
## epi-curve-ggplot
res <- subset(out,select=c(SUSJ,
                           MDA,
                           SUSJM,
                           SUSA,
                           EIJ,
                           EIA,
                           ERJ,
                           ERA,
                           INFJ,
                           INFA,
                           RECJ,
                           RECA, # initial cond
                           time))
ggplot(data=melt(res,id.var="time"),
       mapping=aes(x=time,y=value,group=variable,color=variable))+
  geom_line()

res.sp <- subset(out,select=c(SpJ,
                              SpA,
                              time))
ggplot(data=melt(res.sp,id.var="time"),
       mapping=aes(x=time,y=value,group=variable,color=variable))+
  geom_line()
## using old plots

## Plot data vs the true underlying epidemic.
plot(lbv.new$cumulative_time, lbv.new$ad_sp, type="l", col="red", bty = "n",
     ylim = c(0, max(lbv.new$SpA.ci.u)), xlab = "days", ylab = "seroprevalence")
## Add data
points(lbv.new$cumulative_time, lbv.new$ad_sp, col = "black", pch = 19)
## Add CIs
arrows(lbv.new$cumulative_time, lbv.new$SpA.ci.l, lbv.new$cumulative_time, lbv.new$SpA.ci.u, length = .01,
       angle = 90, code = 3)
lines(lbv.new$cumulative_time, lbv.new$juv_sp, type="l", col="blue")
## Add data
points(lbv.new$cumulative_time, lbv.new$juv_sp, col = "black", pch = 19)
## Add CIs
arrows(lbv.new$cumulative_time, lbv.new$SpJ.ci.l, lbv.new$cumulative_time, lbv.new$SpJ.ci.u, length = .01,
       angle = 90, code = 3)

## simulated data

plot(out$time[13846:15001], out$SpJ[13846:15001], type="l",lty=2, col="blue",ylim=c(0,1),
     xlab="time",xaxt="n",ylab="seroprevalence")
lines(out$time[13846:15001], out$SpA[13846:15001], type="l",lty=2, col="red")
# data
#lines(lbv.new$cumulative_time+13846, lbv.new$ad_sp, col="red",pch=16)
points(lbv.new$cumulative_time+13846, lbv.new$ad_sp, col = "red", pch = 19)
## Add CIs
arrows(lbv.new$cumulative_time+13846, lbv.new$SpA.ci.l, lbv.new$cumulative_time+13846, lbv.new$SpA.ci.u, length = .01,
       angle = 90, code = 3,col="red")
#lines(lbv.new$cumulative_time+13846, lbv.new$juv_sp, type="l", col="blue")
points(lbv.new$cumulative_time+13840, lbv.new$juv_sp, col = "blue", pch = 19)
## Add CIs
arrows(lbv.new$cumulative_time+13840, lbv.new$SpJ.ci.l, lbv.new$cumulative_time+13846, lbv.new$SpJ.ci.u, length = .01,
       angle = 90, code = 3,col="blue")

## Add legends
legend("topright", c("Adult", "Juvenile"),
       col=c("red","blue"),lty=1,bty="n")
legend("top", c("data", NA),
       pch=c(16,NA),bty="n")
legend("top", c(NA, "simulation"),
       lty=c(NA,2),bty="n")

## cor coef.
tTest<-c(lbv.new$cumulative_time+13846)
tTest
tTest<-tTest[-1] # remove first point
cor(out$SpA[tTest],lbv.new$ad_sp[-1])
cor(out$SpJ[tTest],lbv.new$juv_sp[-1])
cor.test(out$SpA[tTest],lbv.new$ad_sp[-1])
cor.test(out$SpJ[tTest],lbv.new$juv_sp[-1])

plot(out$SpJ[tTest],lbv.new$juv_sp[-1],xlim=c(0,0.5),
     ylim=c(0,0.5),col="blue",pch=16,
     ylab="data",xlab="predicted",main="Seroprevalence")
points(out$SpA[tTest],lbv.new$ad_sp[-1],col="red",pch=16)
abline(a=0,b=1,lty=2)
legend("topleft",c("Juvenile","Adult"),pch=16,col=c("blue","red"),bty="n")

###############################################
## to here for plotting using fit data ##########################3

layout(1:2)
plot(log(out$N),out$INFA,type="l")
plot(log(out$N),out$INFJ,type="l")
#plot(out$INFA,out$INFJ,type="l")
layout(1:3)
plot(out$INFJ[10000:15000],log(out$SUSJ[10000:15000]),type="l",ylab="log(SUSJ)",xlab="INFJ")
plot(out$INFA[10000:15000],log(out$SUSA[10000:15000]),type="l",ylab="log(SUSA)",xlab="INFA")
plot(out$INFA[10000:15000],log(out$N[10000:15000]),type="l",ylab="log(N)",xlab="INFA")

layout(1:3)
plot(out$INFJ,log(out$SUSJ),type="l",ylab="log(SUSJ)",xlab="INFJ")
plot(out$INFA,log(out$SUSA),type="l",ylab="log(SUSA)",xlab="INFA")
plot(out$INFA,log(out$N),type="l",ylab="log(N)",xlab="INFA")

layout(1:3)
plot((out$INFA[10000:15000]+out$INFJ[10000:15000]),log(out$N[10000:15000]),type="l",ylab="log(N)",xlab="INF")
plot((out$INFA[10000:15000]+out$INFJ[10000:15000]),log(out$SUSJ[10000:15000]),type="l",ylab="log(SUSJ)",xlab="INF")
plot((out$INFA[10000:15000]+out$INFJ[10000:15000]),log(out$SUSA[10000:15000]),type="l",ylab="log(SUSA)",xlab="INF")

layout(1:3)
plot((out$INFA+out$INFJ),log(out$N),type="l",ylab="log(N)",xlab="INF")
plot((out$INFA+out$INFJ),log(out$SUSJ),type="l",ylab="log(SUSJ)",xlab="INF")
plot((out$INFA+out$INFJ),log(out$SUSA),type="l",ylab="log(SUSA)",xlab="INF")
