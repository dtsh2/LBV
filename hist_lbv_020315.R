################################################################################
##### Writing the pomp objects to estimate the parameters using particle filtering
##
## LBV SEIR model with Frequency-dept transmission
## 25 years
## Seasonal birthing
## maternally-derived antibody
##
## Sensitivity analysis below
##
## clean up
rm(list=ls())
# instal r developer toolbox first (Rtools from Cran R)
## pomp test run lbv
getwd()
setwd("~/GitHub/LBV") # revise as necessary
library(pomp)
#Compiling C code and loading the dll
# system("R CMD SHLIB lbvseirNoSeasFreqTwoParsMeasureBinom.c")
dyn.load("lbvseirNoSeasFreqTwoParsMeasureBinom.dll")
# dyn.unload("lbvseirNoSeasFreqTwoParsMeasureBinom.dll")
pomp(
data = data.frame(
time=seq(from=0,to=365*25,by=1),  # time for simulations to run
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
statenames=c("SUSJ","MDAJ", "SUSJM","EIJ","ERJ","INFJ", "RECJ", "SUSA", "EIA","ERA","INFA", "RECA","SPA","SPJ"),
## the order of the parameters assumed in the native routines:
paramnames=c("BETA",
"RHO",
# "ETA",
"SUSJ.0","MDAJ.0", "SUSJM.0","EIJ.0","ERJ.0","INFJ.0", "RECJ.0", "SUSA.0", "EIA.0","ERA.0","INFA.0", "RECA.0","SPA.0","SPJ.0"),
initializer=function(params,t0,statenames,...){
x0<-params[paste(statenames,".0",sep="")]
names(x0)<-statenames
return(x0)
}
) -> sir
params <- c(
BETA=6.25,
RHO=0.0635,
SUSJ.0=4000,MDAJ.0=4000, SUSJM.0=1000,EIJ.0=1000,ERJ.0=1000,INFJ.0=1000,
RECJ.0=10000,SUSA.0=50000, EIA.0=100,
ERA.0=1000,INFA.0=5000, RECA.0=50000,
SPA.0=0.4994506,SPJ.0=0.5882353)
lbvd<-read.csv("lbv_data_plustime.csv")
DSPJ<-lbvd$DRECJ/(lbvd$DRECJ+lbvd$DSUSJ)
DSPA<-lbvd$DRECA/(lbvd$DRECA+lbvd$DSUSA)
DPOPA<-lbvd$DRECA+lbvd$DSUSA
DPOPJ<-lbvd$DRECJ+lbvd$DSUSJ
DRECA<-lbvd$DRECA
DRECJ<-lbvd$DRECJ
times<-lbvd$cumulative_time
#
lbv.new<-cbind(times,DSPJ,DSPA,DRECJ,DRECA,DPOPJ,DPOPA)
lbv.new<-as.data.frame(lbv.new)
for(ii in 1:12)
{
lbv.new$SpA.ci.l[ii] <- binom.test(lbv.new$DRECA[ii], lbv.new$DPOPA[ii])$conf.int[1]
lbv.new$SpA.ci.u[ii] <- binom.test(lbv.new$DRECA[ii], lbv.new$DPOPA[ii])$conf.int[2]
}
for(ii in 1:12)
{
lbv.new$SpJ.ci.l[ii] <- binom.test(lbv.new$DRECJ[ii], lbv.new$DPOPJ[ii])$conf.int[1]
lbv.new$SpJ.ci.u[ii] <- binom.test(lbv.new$DRECJ[ii], lbv.new$DPOPJ[ii])$conf.int[2]
}
pomp(
data = data.frame(
time=lbv.new[,1],  # time for simulations to run
#  DSPJ = lbv.new[,2],
#  DSPA = lbv.new[,3],
DRECJ = lbv.new[,4],
DRECA = lbv.new[,5],
DPOPJ = lbv.new[,6],
DPOPA = lbv.new[,7]
),
times='time',
t0=0,
## native routine for the process simulator:
rprocess=euler.sim(
step.fun="sir_euler_simulator",
delta.t=1,
#PACKAGE="pomp"  ## may be include if does not work - this is where to look for the c file
## name of the shared-object library containing the PACKAGE="pomp",
),
rmeasure="binomial_rmeasure",
dmeasure="binomial_dmeasure",
## the order of the state variables assumed in the native routines:
statenames=c("SUSJ","MDAJ", "SUSJM","EIJ","ERJ","INFJ", "RECJ", "SUSA", "EIA","ERA","INFA", "RECA","SPA","SPJ"),
obsnames=c(#"DSPJ","DSPA",
"DRECJ","DRECA","DPOPJ","DPOPA"),
## the order of the parameters assumed in the native routines:
paramnames=c("BETA","RHO",#"ETA",
"SUSJ.0","MDAJ.0", "SUSJM.0","EIJ.0","ERJ.0","INFJ.0", "RECJ.0", "SUSA.0", "EIA.0","ERA.0","INFA.0", "RECA.0","SPA.0","SPJ.0"),
initializer=function(params,t0,statenames,...){
x0<-params[paste(statenames,".0",sep="")]
names(x0)<-statenames
return(x0)
}
) -> lbvdat
BetaV = seq(from=0.001,to=40,by=2.5)  # range of beta
RhoV = seq(from=0.001,to=1, by=0.0625) # range of rho
#
BetaV = seq(from=0.01,to=10,by=0.125)  # range of beta
length(BetaV)
RhoV = seq(from=0.001,to=0.5, by=0.03125) # range of rho
length(RhoV)
RhoV
0.03125
0.03125/2
RhoV = seq(from=0.001,to=0.5, by=0.015625) # range of rho
length(RhoV)
0.03125/4
RhoV = seq(from=0.001,to=0.5, by=0.0078125) # range of rho
parametset<- expand.grid(BetaV,RhoV)
dim(parametset)
#EtaV<-rep(0.1,length(parametset[,1]))
nonV = matrix(c(
SUSJ.0=4000,MDAJ.0=4000, SUSJM.0=1000,EIJ.0=1000,ERJ.0=1000,INFJ.0=1000,
RECJ.0=10000,SUSA.0=50000, EIA.0=100,
ERA.0=1000,INFA.0=5000, RECA.0=50000,
SPA.0=0.4994506,SPJ.0=0.5882353),
ncol=14,
nrow=length(parametset[,1]),
byrow=T) #binded with non-varying parameters
dimnames(nonV)[[2]]=c("SUSJ.0","MDAJ.0","SUSJM.0","EIJ.0",
"ERJ.0","INFJ.0","RECJ.0","SUSA.0",
"EIA.0","ERA.0","INFA.0","RECA.0",
"SPA.0","SPJ.0") # naming non-varying columns
parsV<-cbind(parametset,nonV)
BETA = as.numeric(parsV[,1])
RHO = as.numeric(parsV[,2])
#ETA = as.numeric(parsV[,3])
SUSJ.0 = as.numeric(parsV[,3])
MDAJ.0 = as.numeric(parsV[,4])
SUSJM.0 = as.numeric(parsV[,5])
EIJ.0 = as.numeric(parsV[,6])
ERJ.0 = as.numeric(parsV[,7])
INFJ.0 = as.numeric(parsV[,8])
RECJ.0 = as.numeric(parsV[,9])
SUSA.0 = as.numeric(parsV[,10])
EIA.0 = as.numeric(parsV[,11])
ERA.0 = as.numeric(parsV[,12])
INFA.0 = as.numeric(parsV[,13])
RECA.0 = as.numeric(parsV[,14])
SPA.0 = as.numeric(parsV[,15])
SPJ.0 = as.numeric(parsV[,16])
params<-cbind(
BETA,
RHO,
#ETA,
SUSJ.0,
MDAJ.0,
SUSJM.0,
EIJ.0,
ERJ.0,
INFJ.0,
RECJ.0,
SUSA.0,
EIA.0,
ERA.0,
INFA.0,
RECA.0,
SPA.0,
SPJ.0
)
results<-array(NA,dim=c(length(parametset[,1]),3))
## nb check # particles - reduced for training
for (j in 1:length(params[,1])){
results[j,1]<-logLik(pfilter(lbvdat,params=c(params[j,]),Np=100,max.fail=100,tol=1e-20))
#pf<-pfilter(lbvdat,params=c(params[j,]),Np=6000,max.fail=1000,tol=1e-20)
results[j,2:3]<-#c(logLik(pf))}
#
c(params[j,1],params[j,2])
}
library(akima)
library(lattice)
library(tgp)
library(rgl)
library(fields)
write.csv(results,file="mll_surface_fine.csv")
rholab<-expression(symbol(rho))
betalab<-expression(symbol(beta))
zzg <- interp(x=results[,2], #
y=results[,3], #
z=results[,1],
duplicate=T)#,grid.len=c(50,50))#,span=0.1)
max(results[,1])
results[results[,1]==max(results[,1]),]
maxLL<-as.data.frame(t(results[results[,1]==max(results[,1]),]))
names(maxLL)<-c("negll","Beta","Rho")
image(zzg,ann=T,xlim=c(0,max(results[,2])),ylim=c(0,max(results[,3])),
ylab=rholab,xlab=betalab)
#contour(zzg,add=T,labcex=1,drawlabels=T,nlevels=10)
#contour(zzg,add=F,labcex=1,drawlabels=T,nlevels=100)
points(x=maxLL$Beta,y=maxLL$Rho,pch=16,col="black")
points(x=maxLL[2,],y=maxLL[3,],pch=16,col="black")
par(omi=c(1,1,0.5,1))
par(mai=c(1,1,0.8,0.8))
surface(zzg,#col ="#FFFFFF",
ylab=rholab,xlab=betalab,
#zlim=c(0,10),
labcex=1)
points(x=maxLL$Beta,y=maxLL$Rho,pch=16,col="black")
points(x=maxLL$Beta,y=maxLL$Rho,pch=16,col="white")
split.screen(rbind(c(0.1,0.45,0.55, 0.9),
c(0.55, 0.9, 0.55, 0.9),
c(0.55, 0.9, 0.1, 0.45)))
#plot(1:100, rnorm(100), xaxs = "i", ylim = c(-3, 3), xaxt = "n")
#axis(1, at = seq(0, 100, 20))
#par(omi=c(1,1,0.5,1))
#par(mai=c(1,1,0.8,0.8))
screen(2)
par(mar = c(0, 0, 0, 0))
#surface(zzg,#col ="#FFFFFF",
#        ylab=rholab,xlab=betalab,
#        #zlim=c(0,10),
#        labcex=1)
image(zzg,ann=T,xlim=c(0,max(results[,2])),ylim=c(0,max(results[,3])),
ylab=rholab,xlab=betalab,col=tim.colors(10))
abline(h=maxLL[3,],col="white",lty=1)
abline(v=maxLL[2,],col="white",lty=1)
maxLL[3,]
maxLL
abline(h=maxLL$Rho,col="white",lty=1)
abline(v=maxLL$Beta,col="white",lty=1)
screen(1)
par(mar = c(0, 0, 0, 0))
rhop<-which(results[,3]==maxLL[3,])
#results[rhop,]
plot(results[betap,1],results[betap,3],type="l",ylab="Negative log-likelihood",xlab=rholab,
xlim = rev(range(results[betap,1])))
rhop<-which(results[,3]==maxLL[3,])
betap<-which(results[,2]==maxLL[2,])
rhop
betap
rhop<-which(results[,3]==maxLL$Rho)
betap<-which(results[,2]==maxLL$Beta)
betap
rhop
plot(results[betap,1],results[betap,3],type="l",ylab="Negative log-likelihood",xlab=rholab,
xlim = rev(range(results[betap,1])))
screen(3)
par(mar = c(0, 0, 0, 0))
#results[betap,]
plot(results[rhop,2],results[rhop,1],type="l",ylab="Negative log-likelihood",xlab=betalab,
ylim = rev(range(results[rhop,1])))
close.screen(all.screens = TRUE)
tiff("ll_beta_rho.tiff",width=8,height=8,units='in',res=300, compression = "lzw")
image(zzg,ann=T,xlim=c(0,max(results[,2])),ylim=c(0,max(results[,3])),
ylab=rholab,xlab=betalab)
#contour(zzg,add=T,labcex=1,drawlabels=T,nlevels=10)
#contour(zzg,add=F,labcex=1,drawlabels=T,nlevels=100)
points(x=maxLL$Beta,y=maxLL$Rho,pch=16,col="black")
points(x=maxLL[2,],y=maxLL[3,],pch=16,col="black")
dev.off()
tiff("ll_beta_rho_surf.tiff",width=8,height=8,units='in',res=300, compression = "lzw")
par(omi=c(1,1,0.5,1))
par(mai=c(1,1,0.8,0.8))
surface(zzg,#col ="#FFFFFF",
ylab=rholab,xlab=betalab,
#zlim=c(0,10),
labcex=1)
# contour(zzg,add=T,labcex=1,drawlabels=F,nlevels=10)
#surface(zzg,#col ="#FFFFFF",
#        ylab=rholab,xlab=betalab,
#        xlim=c(0,10),ylim=c(0,0.5),
#        labcex=1)
points(x=maxLL$Beta,y=maxLL$Rho,pch=16,col="white")
#points(x=maxLL[2,],y=maxLL[3,],pch=16,col="white")
dev.off()
#########################
# resdf<-as.data.frame(results)#,value=cut(results[,1],breaks=seq(min(results[,1]),max(results[,1]),25)))
# names(resdf)<-c("nll","beta","rho")#,"NegLL")
# library(ggplot2)
# p <- ggplot(resdf) +
#  geom_tile(aes(beta,rho,fill=nll)) +
#  geom_contour(aes(x=beta,y=rho,z=nll), colour="white",bins=25)
# p
# p = p + geom_point(aes(y=0.33, x=32.4),colour="red")
# p
## to here..
#write.csv(resdf,file="pfilterres.csv")
#test<-read.csv("pfilterres.csv",header=T)
#results<-test[2:4]
#head(results)
########
tiff("ll_beta_rho_surf_v1.tiff",width=8,height=8,units='in',res=300, compression = "lzw")
split.screen(rbind(c(0.1,0.45,0.55, 0.9),
c(0.55, 0.9, 0.55, 0.9),
c(0.55, 0.9, 0.1, 0.45)))
#plot(1:100, rnorm(100), xaxs = "i", ylim = c(-3, 3), xaxt = "n")
#axis(1, at = seq(0, 100, 20))
#par(omi=c(1,1,0.5,1))
#par(mai=c(1,1,0.8,0.8))
screen(2)
par(mar = c(0, 0, 0, 0))
#surface(zzg,#col ="#FFFFFF",
#        ylab=rholab,xlab=betalab,
#        #zlim=c(0,10),
#        labcex=1)
image(zzg,ann=T,xlim=c(0,max(results[,2])),ylim=c(0,max(results[,3])),
ylab=rholab,xlab=betalab,col=tim.colors(10))
abline(h=maxLL$Rho,col="white",lty=1)
abline(v=maxLL$Beta,col="white",lty=1)
screen(1)
par(mar = c(0, 0, 0, 0))
rhop<-which(results[,3]==maxLL$Rho)
betap<-which(results[,2]==maxLL$Beta)
#results[rhop,]
plot(results[betap,1],results[betap,3],type="l",ylab="Negative log-likelihood",xlab=rholab,
xlim = rev(range(results[betap,1])))
screen(3)
par(mar = c(0, 0, 0, 0))
#results[betap,]
plot(results[rhop,2],results[rhop,1],type="l",ylab="Negative log-likelihood",xlab=betalab,
ylim = rev(range(results[rhop,1])))
close.screen(all.screens = TRUE)
dev.off()
savehistory("~/GitHub/LBV/hist_lbv_020315.R")
