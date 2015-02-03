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
  BETA=7.635,
  RHO=0.0556875,
  SUSJ.0=4000,MDAJ.0=4000, SUSJM.0=1000,EIJ.0=1000,ERJ.0=1000,INFJ.0=1000,
  RECJ.0=10000,SUSA.0=50000, EIA.0=100,
 ERA.0=1000,INFA.0=5000, RECA.0=50000,
  SPA.0=0.4994506,SPJ.0=0.5882353)
#
sim <- simulate(sir,params=c(params),seed=3593885L,
                nsim=1,states=T,obs=F,as.data.frame=T) # 
class(sir) # pomp object
class(sim) # data frame - change states, obs and data.frame if want pomp obj
#
tiff("juv_dynamics.tiff",width=8,height=8,units='in',res=300, compression = "lzw")

plot(sim$time,sim$SUSJ,type="l",xlab="days",ylab="numbers")
points(sim$time,sim$RECJ,col="green",type="l")
points(sim$time,sim$MDA,col="brown",type="l")
points(sim$time,sim$INFJ,col="red",type="l")
legend("topleft",c("Susceptible","Recovered","Maternal-antibody","Infected"),
       col=c("black","green","brown","red"),lty=1,bty="n")
dev.off()

tiff("juv_2yrs_dynamics.tiff",width=8,height=8,units='in',res=300, compression = "lzw")

plot(sim$time[8396:9126],sim$SUSJ[8396:9126],type="l",
     ylim=c(0,max(c(sim$SUSJ,sim$RECJ,sim$MDA))),xlab="days",ylab="numbers")
points(sim$time[8396:9126],sim$RECJ[8396:9126],col="green",type="l")
points(sim$time[8396:9126],sim$MDA[8396:9126],col="brown",type="l")
points(sim$time[8396:9126],sim$INFJ[8396:9126],col="red",type="l")
legend("topleft",c("Susceptible","Recovered","Maternal-antibody","Infected"),
       col=c("black","green","brown","red"),lty=1,bty="n")
dev.off()
#

tiff("adult_dynamics.tiff",width=8,height=8,units='in',res=300, compression = "lzw")

plot(sim$time,sim$RECA,col="green",type="l",ylim=c(0,max(c(sim$RECA,sim$SUSA)))
     ,xlab="days",ylab="numbers")
points(sim$time,sim$SUSA,type="l")
points(sim$time,sim$INFA,col="red",type="l")
legend("topleft",c("Susceptible","Recovered","Infected"),
       col=c("black","green","red"),lty=1,bty="n")
dev.off()

tiff("adult_2yrs_dynamics.tiff",width=8,height=8,units='in',res=300, compression = "lzw")

plot(sim$time[8396:9126],sim$RECA[8396:9126],col="green",type="l",ylim=c(0,max(c(sim$RECA,sim$SUSA)))
     ,xlab="days",ylab="numbers")
points(sim$time[8396:9126],sim$SUSA[8396:9126],type="l")
points(sim$time[8396:9126],sim$INFA[8396:9126],col="red",type="l")
legend("left",c("Susceptible","Recovered","Infected"),
       col=c("black","green","red"),lty=1,bty="n")
dev.off()
#
tiff("sp_dynamics.tiff",width=8,height=8,units='in',res=300, compression = "lzw")

plot(sim$time,sim$SPA,type="l",col="darkgreen",ylim=c(0,1)
     ,xlab="days",ylab="seroprevalence (%)")
points(sim$time,sim$SPJ,type="l",col="darkblue")
legend("topright",c("Adult","Juvenile"),
       col=c("darkgreen","darkblue"),lty=1,bty="n")
dev.off()

tiff("sp_2yrs_dynamics.tiff",width=8,height=8,units='in',res=300, compression = "lzw")

plot(sim$time[8396:9126],sim$SPA[8396:9126],type="l",col="darkgreen",ylim=c(0,1)
     ,xlab="days",ylab="seroprevalence (%)")
points(sim$time[8396:9126],sim$SPJ[8396:9126],type="l",col="darkblue")
legend("topright",c("Adult","Juvenile"),
       col=c("darkgreen","darkblue"),lty=1,bty="n")
dev.off()
#
tiff("inf_dynamics.tiff",width=8,height=8,units='in',res=300, compression = "lzw")

plot(sim$time,sim$INFA,type="l",col="orange",ylim=c(0,500)
     ,xlab="days",ylab="numbers")
points(sim$time,sim$INFJ,type="l",col="red")
legend("topright",c("adults","juveniles"),
       col=c("orange","red"),lty=1,bty="n")
dev.off()

tiff("inf_2yrs_dynamics.tiff",width=8,height=8,units='in',res=300, compression = "lzw")

plot(sim$time[8396:9126],sim$INFA[8396:9126],type="l",col="orange",ylim=c(0,200)
     ,xlab="days",ylab="numbers")
points(sim$time[8396:9126],sim$INFJ[8396:9126],type="l",col="red")
legend("topright",c("adults","juveniles"),
       col=c("orange","red"),lty=1,bty="n")
dev.off()
#########################################################
## multiple sims...
#
sim <- simulate(sir,params=c(params),#seed=3593885L,
                nsim=100,states=T,obs=F,as.data.frame=T) # 
class(sir) # pomp object
class(sim) # data frame - change states, obs and data.frame if want pomp obj

with(sim,plot(time,SPA,type="n"))
for(i in unique(sim$sim))
  with(sim[ sim$sim==i,],lines(time,SPA,col=adjustcolor("grey", alpha=0.5)))
test<-tapply(sim$SPA, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="red")

simspa<-test[c(7301,7666,7725,7756,7970,8031,8090,8121,8212,8335,8396,8455)]

for(i in unique(sim$sim))
  with(sim[ sim$sim==i,],lines(time,SPJ,col=adjustcolor("grey", alpha=0.5)))
test<-tapply(sim$SPJ, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="black")

simspj<-test[c(7301,7666,7725,7756,7970,8031,8090,8121,8212,8335,8396,8455)]

## plots Ns

tiff("mean_numbers.tiff",width=8,height=8,units='in',res=300, compression = "lzw")

with(sim,plot(time,SUSA,type="n",xlim=c(7000,9000),ylab="numbers"))
#for(i in unique(sim$sim))
  sim$RECA[sim$RECA==0] <- NA
#with(sim[ sim$sim==i,],points(time,RECA,col=adjustcolor("green", alpha=0.5)))
test<-tapply(sim$RECA, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="black",lwd=2)

#for(i in unique(sim$sim))
  sim$SUSA[sim$SUSA==0] <- NA
#with(sim[ sim$sim==i,],points(time,SUSA,col=adjustcolor("blue", alpha=0.5)))
test<-tapply(sim$SUSA, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="blue",lwd=2)

#dev.off()

##
#tiff("juv_points.tiff",width=8,height=8,units='in',res=300, compression = "lzw")

#with(sim,plot(time,SUSJ,type="n",xlim=c(6000,9000),ylab="numbers"))
#for(i in unique(sim$sim))
#  with(sim[ sim$sim==i,],points(time,RECJ,col=adjustcolor("darkgreen", alpha=0.5)))
sim$RECJ[sim$RECJ==0] <- NA
test<-tapply(sim$RECJ, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="black",lty=3,lwd=2)

#for(i in unique(sim$sim))
  sim$MDA[sim$MDA==0] <- NA
#with(sim[ sim$sim==i,],points(time,MDA,col=adjustcolor("grey", alpha=0.5)))
test<-tapply(sim$MDA, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="green",lty=3,lwd=2)

#for(i in unique(sim$sim))
  sim$SUSJ[sim$SUSJ==0] <- NA
#with(sim[ sim$sim==i,],points(time,SUSJ,col=adjustcolor("blue", alpha=0.5)))
test<-tapply(sim$SUSJ, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="blue",lty=3,lwd=2)

legend("topleft",c("recovered adult","susceptible adult","recovered juvenile","maternal immunity","susceptibe juvenile"),
       col=c("black","blue","black","green","blue"),lty=c(1,1,3,3,3),lwd=2,
       cex=0.9,bg="white",box.col="white")
box(col="black")
dev.off()

# N infected
tiff("inf_points.tiff",width=8,height=8,units='in',res=300, compression = "lzw")

with(sim,plot(time,INFA,type="n",xlim=c(8000,9000),ylim=c(0,300),ylab="numbers"))
#sim$INFA[sim$INFA==0] <- NA
for(i in unique(sim$sim))
  with(sim[ sim$sim==i,],points(time, INFA,col=adjustcolor("orange", alpha=0.5)))
#sim$INFA[sim$INFA==0] <- NA
test<-tapply(sim$INFA, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="black")

for(i in unique(sim$sim))
#sim$INFJ[sim$INFJ==0] <- NA
  with(sim[ sim$sim==i,],points(time, INFJ,col=adjustcolor("red", alpha=0.5)))
#sim$INFJ[sim$INFJ==0] <- NA
test<-tapply(sim$INFJ, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="white")

legend("topright",c("juvenile","adult"),col=c("orange","red"),pch=16,bty="n")
dev.off()
##
## library(ggplot2)
#plot
## ggplot(data=data,aes(times,norm,colour=as.factor(sim))) +
##  geom_line()

#
#######

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

cr<-cor(c(simspa,simspj),c(DSPA,DSPJ))
cr
crt<-cor.test(c(simspa,simspj),c(DSPA,DSPJ))
crt
lmp<-lm(c(DSPA,DSPJ)~(c(simspa,simspj)))

tiff("cor_pred_vs_data.tiff",width=8,height=8,units='in',res=300, compression = "lzw")
plot(c(simspa,simspj),c(DSPA,DSPJ),ylim=c(0,max(c(DSPA,DSPJ))),xlim=c(0,max(c(simspa,simspj))),col=c(rep("red",12),rep("blue",12)),
     ylab="Mean seroprevalence (%)",xlab="Predicted seroprevalence (%)",pch=19)
legend("topleft",c("juvenile","adult"),col=c("blue","red"),pch=19,bty="n")
abline(lmp)
corval<-format(round(cr[1],3))
corlab<-bquote(plain(R==.(corval)))
text(c(0.25),(0.45),corlab)
dev.off()

#sim <- simulate(sir,params=c(params),seed=3593885L,
#                nsim=1,states=T,obs=F,as.data.frame=T) # 
## correlation
# plot(c(sim$SPJ[c(7301,7666,7725,7756,7970,8031,8090,8121,8212,8335,8396,8455)],
#       sim$SPA[c(7301,7666,7725,7756,7970,8031,8090,8121,8212,8335,8396,8455)]),
#     c(DSPJ,DSPA),ylim=c(0,0.6),xlim=c(0,0.7),col=c(rep("red",12),rep("blue",12)),
#     ylab="Mean seroprevalence (%)",xlab="Predicted seroprevalence (%)",pch=19)
#legend("topleft",c("juvenile","adult"),col=c("red","blue"),pch=19,bty="n")

##
#cr<-cor(c(sim$SPJ[c(7301,7666,7725,7756,7970,8031,8090,8121,8212,8335,8396,8455)],
#       sim$SPA[c(7301,7666,7725,7756,7970,8031,8090,8121,8212,8335,8396,8455)]),
#     c(DSPJ,DSPA))
#
#cor.test(c(sim$SPJ[c(7301,7666,7725,7756,7970,8031,8090,8121,8212,8335,8396,8455)],
#      sim$SPA[c(7301,7666,7725,7756,7970,8031,8090,8121,8212,8335,8396,8455)]),
#    c(DSPJ,DSPA))
#lmp<-lm(c(DSPJ,DSPA)~(c(sim$SPJ[c(7301,7666,7725,7756,7970,8031,8090,8121,8212,8335,8396,8455)],
#      sim$SPA[c(7301,7666,7725,7756,7970,8031,8090,8121,8212,8335,8396,8455)])))
#abline(lmp)
#text(0.45,0.45,round(cr[1],3))

##
## binomial CI for SP
## each sample get CIs

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

## Plot data vs the true underlying epidemic.
tiff("sp_data_sim.tiff",width=8,height=8,units='in',res=300, compression = "lzw")
#plot(lbv.new$times, lbv.new$DSPA, type="l", col="red", bty = "n",ylim=c(0,1),
#     #ylim = c(0, max(lbv.new$SpA.ci.u)), 
#     xlab = "days", ylab = "seroprevalence")
## Add data
plot(lbv.new$times, lbv.new$DSPA, col = "red", pch = 19,ylim=c(0,1.5),
     xlab = "days", ylab = "seroprevalence (%)",yaxt="n")
axis(2,at=c(0,0.5,1))
## Add CIs
arrows(lbv.new$times, lbv.new$SpA.ci.l, lbv.new$times, lbv.new$SpA.ci.u, length = .01,
       angle = 90, code = 3,lty=3,col="red")
# lines(lbv.new$times, lbv.new$DSPJ, type="l", col="blue")
## Add data
points(lbv.new$times, lbv.new$DSPJ, col = "blue", pch = 19)
## Add CIs
arrows(lbv.new$times, lbv.new$SpJ.ci.l, lbv.new$times, lbv.new$SpJ.ci.u, length = .01,
       angle = 90, code = 3,lty=3,col="blue")
## Add legends
legend("topleft", c("data", "simulation mean",
                "adult", "juvenile","95% CI data"),
       col=c("black","black","red","blue","black"),pch=c(16,1,16,16,NA),
       lty=c(0,0,0,0,3),lwd=c(0,0,0,0,1),bty="n")

points(lbv.new$times, simspa, col = "red", pch = 1)
points(lbv.new$times, simspj, col = "blue", pch = 1)

par(fig = c(0.5, 1, 0.4, 1), #mar=c(0,0,0,0),
    new=TRUE)
#par(mgp = c(0, 1, 0))
plot(c(simspa,simspj),c(DSPA,DSPJ),ylim=c(0,max(c(DSPA,DSPJ))),
     xlim=c(0,max(c(simspa,simspj))),col=c(rep("red",12),rep("blue",12)),
     ylab="",xlab="",pch=19)
#legend("topleft",c("juvenile","adult"),col=c("blue","red"),pch=19,bty="n")
mtext("predicted (%)",side=1,line=2)
mtext("data (%)",side=2,line=2)
abline(lmp)
corval<-format(round(cr[1],3))
corlab<-bquote(plain(R==.(corval)))
text(c(0.25),(0.45),corlab)

dev.off()

## for all parameter sets....
# res1<-array(NA,dim=c(100,12,2))
# for (j in 1:100){
## change # sims
#  sim <- simulate(sir,params=c(params),#seed=3493885L,
#                  nsim=1,states=T,obs=F,as.data.frame=T) # 
#  outres <- sim[c(7301,7666,7725,7756,7970,8031,8090,8121,8212,8335,8396,8455),] # select last #s
#  N = array(NA,c(12,2)) # same dimensions as No. runs * outputs I want
#  for (i in 1:12){  
#  N[i,1]<-outres[i,13]
#  N[i,2]<-outres[i,14]
#  }
#  res1[j,,]<-N
#}
#res2<-array(NA,dim=c(12,2))
#res2[,1]<-colMeans(res1[,,1])
#res2[,2]<-colMeans(res1[,,2])
#plot(res2[,2])
#
#
#cr<-cor(c(res2[,1],res2[,2]),c(DSPA,DSPJ))
#cr
#crt<-cor.test(c(res2[,1],res2[,2]),c(DSPA,DSPJ))
#crt
#lmp<-lm(c(DSPA,DSPJ)~(c(res2[,1],res2[,2])))
#
#tiff("cor_pred_vs_data.tiff",width=8,height=8,units='in',res=300, compression = "lzw")
#plot(c(res2[,1],res2[,2]),c(DSPA,DSPJ),ylim=c(0,max(c(DSPA,DSPJ))),xlim=c(0,max(c(res2[,1],res2[,2]))),col=c(rep("red",12),rep("blue",12)),
#     ylab="Mean seroprevalence (%)",xlab="Predicted seroprevalence (%)",pch=19)
#legend("topleft",c("juvenile","adult"),col=c("blue","red"),pch=19,bty="n")
#abline(lmp)
#corval<-format(round(cr[1],3))
#corlab<-bquote(plain(R==.(corval)))
#text(c(0.15),(0.05),corlab)
#dev.off()

#sim <- simulate(sir,params=c(params),seed=3593885L,
#                nsim=1,states=T,obs=F,as.data.frame=T) # 
## correlation
# plot(c(sim$SPJ[c(7301,7666,7725,7756,7970,8031,8090,8121,8212,8335,8396,8455)],
#       sim$SPA[c(7301,7666,7725,7756,7970,8031,8090,8121,8212,8335,8396,8455)]),
#     c(DSPJ,DSPA),ylim=c(0,0.6),xlim=c(0,0.7),col=c(rep("red",12),rep("blue",12)),
#     ylab="Mean seroprevalence (%)",xlab="Predicted seroprevalence (%)",pch=19)
#legend("topleft",c("juvenile","adult"),col=c("red","blue"),pch=19,bty="n")

##
#cr<-cor(c(sim$SPJ[c(7301,7666,7725,7756,7970,8031,8090,8121,8212,8335,8396,8455)],
#       sim$SPA[c(7301,7666,7725,7756,7970,8031,8090,8121,8212,8335,8396,8455)]),
#     c(DSPJ,DSPA))
#
#cor.test(c(sim$SPJ[c(7301,7666,7725,7756,7970,8031,8090,8121,8212,8335,8396,8455)],
#      sim$SPA[c(7301,7666,7725,7756,7970,8031,8090,8121,8212,8335,8396,8455)]),
#    c(DSPJ,DSPA))
#lmp<-lm(c(DSPJ,DSPA)~(c(sim$SPJ[c(7301,7666,7725,7756,7970,8031,8090,8121,8212,8335,8396,8455)],
#      sim$SPA[c(7301,7666,7725,7756,7970,8031,8090,8121,8212,8335,8396,8455)])))
#abline(lmp)
#text(0.45,0.45,round(cr[1],3))

##
## binomial CI for SP
## each sample get CIs

#lbv.new<-as.data.frame(lbv.new)
#for(ii in 1:12)
#{
#  lbv.new$SpA.ci.l[ii] <- binom.test(lbv.new$DRECA[ii], lbv.new$DPOPA[ii])$conf.int[1]
#  lbv.new$SpA.ci.u[ii] <- binom.test(lbv.new$DRECA[ii], lbv.new$DPOPA[ii])$conf.int[2]    
#}
#
#for(ii in 1:12)
#{
#  lbv.new$SpJ.ci.l[ii] <- binom.test(lbv.new$DRECJ[ii], lbv.new$DPOPJ[ii])$conf.int[1]
#  lbv.new$SpJ.ci.u[ii] <- binom.test(lbv.new$DRECJ[ii], lbv.new$DPOPJ[ii])$conf.int[2]    
#}

## Plot data vs the true underlying epidemic.
#tiff("sp_data.tiff",width=8,height=8,units='in',res=300, compression = "lzw")
#
#plot(lbv.new$times, lbv.new$DSPA, type="l", col="red", bty = "n",ylim=c(0,1),
     #ylim = c(0, max(lbv.new$SpA.ci.u)), 
#     xlab = "days", ylab = "seroprevalence")
## Add data
#points(lbv.new$times, lbv.new$DSPA, col = "red", pch = 19)
## Add CIs
#arrows(lbv.new$times, lbv.new$SpA.ci.l, lbv.new$times, lbv.new$SpA.ci.u, length = .01,
#       angle = 90, code = 3,col="red")
#lines(lbv.new$times, lbv.new$DSPJ, type="l", col="blue")
## Add data
#points(lbv.new$times, lbv.new$DSPJ, col = "blue", pch = 19)
## Add CIs
#arrows(lbv.new$times, lbv.new$SpJ.ci.l, lbv.new$times, lbv.new$SpJ.ci.u, length = .01,
#       angle = 90, code = 3,col="blue")
## Add legends
#legend("topright", c("Adult", "Juvenile"),
#       col=c("red","blue"),lty=1,bty="n")
#dev.off()

#points(lbv.new$times,
#     sim$SPJ[c(7301,7666,7725,7756,7970,8031,8090,8121,8212,8335,8396,8455)], 
#     pch=19)
#points(lbv.new$times,
#       sim$SPA[c(7301,7666,7725,7756,7970,8031,8090,8121,8212,8335,8396,8455)], 
#       pch=19)

################
##
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

plot(lbvdat)
#########

pf<-pfilter(lbvdat,params=c(params),Np=100,max.fail=100,tol=1e-20)
logLik(pf)
coef(pf)

BetaV = seq(from=0.001,to=40,by=2.5)  # range of beta
RhoV = seq(from=0.001,to=1, by=0.0625) # range of rho
#
BetaV = seq(from=0.001,to=10,by=1.25)  # range of beta
RhoV = seq(from=0.001,to=0.5, by=0.03125) # range of rho
#
BetaV = seq(from=0.01,to=10,by=0.125)  # range of beta
RhoV = seq(from=0.001,to=0.5, by=0.0078125) # range of rho
#
parametset<- expand.grid(BetaV,RhoV)
dim(parametset)
#EtaV<-rep(0.1,length(parametset[,1]))
#paramsV<-cbind(parametset,EtaV)
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

##
#
#
write.csv(results,file="mll_surface_fine.csv")

##

library(akima)
library(lattice)
library(tgp)
library(rgl)
library(fields)


rholab<-expression(symbol(rho))
betalab<-expression(symbol(beta))

zzg <- interp(x=results[,2], #
              y=results[,3], # 
              z=results[,1],
              duplicate=T)#,grid.len=c(50,50))#,span=0.1)
## narrow figure

#plot(results[,1])
max(results[,1])
results[results[,1]==max(results[,1]),]
maxLL<-as.data.frame(t(results[results[,1]==max(results[,1]),]))
names(maxLL)<-c("negll","Beta","Rho")

tiff("ll_beta_rho.tiff",width=8,height=8,units='in',res=300, compression = "lzw")

image(zzg,ann=T,xlim=c(0,max(results[,2])),ylim=c(0,max(results[,3])),
      ylab=rholab,xlab=betalab)
#contour(zzg,add=T,labcex=1,drawlabels=T,nlevels=10)
#contour(zzg,add=F,labcex=1,drawlabels=T,nlevels=100)
points(x=maxLL$Beta,y=maxLL$Rho,pch=16,col="black")
points(x=maxLL[2,],y=maxLL[3,],pch=16,col="black")
dev.off()

########################
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
##
####################################
## for LHS parameter set
## Calling requisite libraries for parallel computing
#
#library(foreach)
#library(doSNOW)
#
#Setting up "parallel backend"
#
#w<-makeCluster(3,type="SOCK") # makes the cluster, i.e. no of cores ABC = 8 cores, DEF = 12 see performance to see # of cores
#registerDoSNOW(w) # 
#
#Checks that the number of workers is set up correctly.
#
#getDoParWorkers()
#
## ~~~~~~~~~~ LHS SAMPLING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#############
# NB not testing K - carrying cap, or rate of aging..
dyn.unload("lbvseirNoSeasFreqTwoParsMeasure.dll")
dyn.unload("lbvseirNoSeasFreqTwoParsMeasureBinom.dll")

dyn.unload("lbvseirNoSeasFreqTwoParsMeasure.dll")
dyn.unload("lbvseirNoSeasFreqMeasure.dll")

# system("R CMD SHLIB lbvseirNoSeasFreqMeasure.c")

dyn.load("lbvseirNoSeasFreqMeasure.dll")

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
  paramnames=c("BETA","MU","DELTA","ALPHA","RHO","SIGMA","K","EPSILON",
               "TAU","PSI","KAPPA","S","OMEGA","PHI","GAMMA",'ETA',
               "SUSJ.0", 'MDAJ.0', "SUSJM.0","EIJ.0","ERJ.0","INFJ.0",
               "RECJ.0", "SUSA.0", "EIA.0","ERA.0","INFA.0", "RECA.0","SPA.0","SPJ.0"),
  initializer=function(params,t0,statenames,...){
    x0<-params[paste(statenames,".0",sep="")]
    names(x0)<-statenames
    return(x0)
  }
) -> sir


library(lhs)

nspaces=100 ## how many bins/intervals, 100 needed for PRCC

hypercube=randomLHS(n=nspaces, k=10) ## function with N columns
dimnames(hypercube)[[2]]=c("beta","mu","delta","alpha","rho",
                           "sigma","K",#"epsilon",
                           "tau","psi",#"k",
                           "s")#,
#"phi")  # named columns with parameter names

mins = c(   			    ## set mins for each parameters-exclude those not to be varied if any
  beta= 0.6,		     # transmission 
  mu= 0.0000510492,  			# natural mortality from my CMR study
  delta= 0.0002312247, 			# juvenile mortality rate
  alpha= 0.02,			# dis induced mortality
  rho= 0.03,				# probability that exposure/infection will lead to infection & infectiousness (and dead)
  sigma= 0.002083333,			# incubation period 
  K= 100000,			      # K, however, population size fluctuates...up to >1*10^6
  #  epsilon= 1/10,				# rate of aging for those juveniles, should be ~ annual
  tau= 0.004166667, 	      # rate of seroconversion
  psi= 10,			######### this will need to be 1/psi for analysis 
  #k=,        			# nb this is birth rate which halved
  s=14.83/10)#,   # very synchronous
#  phi=0.01)

maxs = c( 				    ## set mins for each parameters-exclude those not to be varied if any
  beta= 60,           # transmission
  mu= 0.00510492,  	          # natural mortality from my CMR study
  delta= 0.02312247,            # juvenile mortality rate
  alpha= 0.9,		    # dis induced mortality
  rho= 0.99,			    # probability 
  sigma= 0.2,          # incubation period 
  K= 100000*10,			    # K, however, population size fluctuates...up to >1*10^6
  #  epsilon= 1*10,	                # rate of aging for those juveniles, should be ~ annual
  tau= 0.4167,         # rate of seroconversion
  psi= 200,      ###### will need to be 1/psi for the model
  #k=1.5,                       # nb this is peak
  s=14.83*10)#,
#  phi=0.99)                      

diffs=maxs-mins ## range of each variable

hypercubeadj = hypercube # create copy of hypercube samples to modify, hypercube adjusted; i.e. new matrix
for (i in 1:ncol(hypercube)){
  hypercubeadj[,i]=hypercube[,i]*diffs[i]+mins[i] # scale samples to difference and add minimum
}

head(hypercubeadj)
dim(hypercubeadj)

dimnames(hypercubeadj)[[2]]=c("beta","mu","delta","alpha","rho",
                              "sigma","K",#"epsilon",
                              "tau","psi",#"k",
                              "s")#,"phi")

nonVarying = matrix(c(
  #K= 1000000,			# K, however, population size fluctuates...up to >1*10^6
  epsilon= 1/365,	# rate of aging for those juveniles, should be ~ annual - per cap per year
  omega= 1/365, # for pb per year - 1 per year
  phi= 0.0, # for timing - ...
  k=1.5/365,
  eta=0.1,
  SUSJ.0=4000,MDAJ.0=4000, SUSJM.0=1000,EIJ.0=1000,ERJ.0=1000,INFJ.0=1000,
  RECJ.0=10000,SUSA.0=50000, EIA.0=100,
  ERA.0=1000,INFA.0=5000, RECA.0=50000,
  SPA.0=0.4994506,SPJ.0=0.5882353),
  ncol=19,
  nrow=100,
  byrow=T) #binded with non-varying parameters


dimnames(nonVarying)[[2]]=c("omega","epsilon","phi","k","eta","SUSJ.0","MDAJ.0", "SUSJM.0","EIJ.0","ERJ.0","INFJ.0", "RECJ.0", "SUSA.0", "EIA.0","ERA.0","INFA.0", "RECA.0","SPA.0","SPJ.0") # naming non-varying columns

gamma.day <- (365-hypercubeadj[,9])
gamma <-1/gamma.day
hypercubeadj[,9] <-1/hypercubeadj[,9]
#hypercubeadj[,1] <-hypercubeadj[,1]*hypercubeadj[,7]
fullParamSets = cbind(nonVarying,hypercubeadj,gamma) # full parameter set

head(fullParamSets)
dim(fullParamSets)

dimnames(fullParamSets)[[2]]=c("OMEGA", # 1
                               "EPSILON",
                               "PHI",
                               "KAPPA",
                               "ETA", # 5
                               "SUSJ.0",
                               "MDAJ.0",
                               "SUSJM.0",
                               "EIJ.0",
                               "ERJ.0", # 10
                               "INFJ.0",
                               "RECJ.0",
                               "SUSA.0",
                               "EIA.0",
                               "ERA.0", # 15
                               "INFA.0",
                               "RECA.0",
                               "SPA.0",
                               "SPJ.0",
                               "BETA", # 20
                               "MU",
                               "DELTA",
                               "ALPHA",
                               "RHO",
                               "SIGMA", # 25
                               "K",
                               "TAU",
                               "PSI",
                               # "KAPPA",
                               "S",
                               "GAMMA") # 30

# order for pomp/C model:  
BETA = fullParamSets[,20]
MU = fullParamSets[,21]
DELTA = fullParamSets[,22]
ALPHA = fullParamSets[,23]
RHO = fullParamSets[,24]
SIGMA = fullParamSets[,25]
K = fullParamSets[,26]
EPSILON = fullParamSets[,2]
TAU = fullParamSets[,27]
PSI = fullParamSets[,28]
KAPPA = fullParamSets[,4]
S = fullParamSets[,29]
OMEGA = fullParamSets[,1]
PHI = fullParamSets[,3]
GAMMA = fullParamSets[,30]
ETA = fullParamSets[,5]
SUSJ.0 = fullParamSets[,6]
MDAJ.0 = fullParamSets[,7]
SUSJM.0 = fullParamSets[,8]
EIJ.0 = fullParamSets[,9]
ERJ.0 = fullParamSets[,10]
INFJ.0 = fullParamSets[,11]
RECJ.0 = fullParamSets[,12]
SUSA.0 = fullParamSets[,13]
EIA.0 = fullParamSets[,14]
ERA.0 = fullParamSets[,15]
INFA.0 = fullParamSets[,16]
RECA.0 = fullParamSets[,17]
SPA.0 = fullParamSets[,18]
SPJ.0 = fullParamSets[,19]


paramset<-cbind(BETA,MU,DELTA,ALPHA,RHO,SIGMA,K,EPSILON,TAU,PSI,KAPPA,S,OMEGA,PHI,GAMMA,ETA,
                SUSJ.0, MDAJ.0, SUSJM.0,EIJ.0,ERJ.0,INFJ.0, RECJ.0, SUSA.0, EIA.0,ERA.0,INFA.0, RECA.0,SPA.0,SPJ.0)

######################################################################################
results<-array(NA,dim=c(100,1,5))

# for one parameter set....
out1 <-simulate(sir,params=c(paramset[1,]),
                seed=1493885L,nsim=100,states=T,obs=F,as.data.frame=T) #

outres1 <- out1[seq(from=9126,to=912600,by=9126),] # select last #s
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
M1[4] = mean(N1[1:100,4]) # pop ext
M1[5] = mean(N1[1:100,5]) # lbv ext
rm(out1)
M1
results[1,,]<-M1
results[1,,]

##########################################################33
# for all parameter sets....
results<-array(NA,dim=c(100,1,5))
## change # sims
for (j in 1:length(paramset[,1])){
  out <-simulate(sir,params=c(paramset[j,]),
                 seed=1493885L,nsim=100,states=T,obs=F,as.data.frame=T) #
  outres <- out[seq(from=9126,to=912600,by=9126),] # select last #s
  N = array(0,c(100,5)) # same dimensions as No. runs * outputs I want
  for (i in 1:100){ # each stochastic run
    N[i,1]<-sum(outres[i,1:12])
    N[i,2]<-(sum(outres[i,6],outres[i,11])/sum(outres[i,1:12]))*100 # prevalence; total
    N[i,3]<-((outres[i,12])/(sum(outres[i,8:12])))*100 # adult seroprevalence; total
    N[i,4]<-ifelse(sum(outres[i,1:12])>0,1,0) # population persistence for each run
    N[i,5]<-ifelse(sum(outres[i,6],outres[i,11])>0,1,0) # pathogen persistence for each run
  }
  N[is.na(N)]<- 0
  ## now average
  M = array(0,c(1,5))
  M[1] = mean(N[1:100,1]) # population size
  M[2] = mean(N[1:100,2]) # prevalence
  M[3] = mean(N[1:100,3]) # adult seroprevalence
  M[4] = mean(N[1:100,4]) # mean pop persistence
  M[5] = mean(N[1:100,5]) # mean path persistence
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
#####Functions to calculate and plot partial-rank correlation coefficients
#####(PRCCs) between parameter values and model output.
#####
#####Written by: Michael Buhnerkempe
#####      Date: Oct. 7, 2011
#####
##### Functions:
#####          prcc - calculate the partial rank correlation coefficient
#####     plot.prcc - plot a prcc object
#####
##### A brief example is presented at the end
################################################################################

################################################################################
## prcc - function to calculate partial rank correlation coefficients between
##           each of p parameters and k model outputs using n different
##           observations (number of parameter sets)
##
##    Arguments:
##          par.mat = n x p matrix containing the parameter values obtained
##                        from Latin Hypercube Sampling
##      model.output = n x k matrix containing the model outputs
##           routine = how should the PRCCs be calculated? One of:
##                        "blower" - calculated according to Appendix A in
##                                   Blower and Dowlatabadi (1994). DEFAULT.
##                        "regression" - calculated using a regression approach.
##                                   Here, the partial correlation coefficient
##                                   is defined as the correlation between the
##                                   residuals after regressing a model output
##                                   on all of the parameters except for the
##                                   parameter of interest and the residuals
##                                   after regressing the parameter of interest
##                                   on all of the other parameters. This can
##                                   be interpreted as the correlation between
##                                   of a parameter and the model output when
##                                   the effects of all other parameters have
##                                   been removed.
##         par.names = names of parameters
##      output.names = names of model outputs
##               ... = additional arguments to be passed to functions called
##                     within this function
##
##
##    Output attributes:
##      $par.matrix = original matrix of parameter set
##      $model.output = original model output
##      $(model output name) = prcc results for the named model output

prcc = function( par.mat, model.output, routine = "blower",
                 par.names = NA,output.names = NA, ...){
  
  #Make sure par.mat and model.output are matrices
  par.mat = as.matrix(par.mat)
  model.output = as.matrix(model.output)
  
  #How many parameter sets are there?
  n = length(par.mat[,1])
  
  #How many parameters are there?
  p = length(par.mat[1,])
  
  #How many model outputs are we calculating PRCCs for?
  k = length(model.output[1,])
  
  #Find the ranks of the parameter values and the model output
  par.rank = apply(par.mat,2,rank,...)
  output.rank = apply(model.output,2,rank,...)
  
  #What is the average rank?
  ave.rank = (1 + n)/2
  
  #Create a list object to store the PRCCs
  results = list()
  
  results$num.sets = n
  
  #Try to automatically get parameter and output names if they are not
  #given
  if( sum(is.na(par.names)) > 0){par.names=dimnames(par.mat)[[2]]}
  if( sum(is.na(output.names)) > 0){output.names=dimnames(model.output)[[2]]}
  
  ########################################################################
  #Calculate the PRCCs using Appendix A from Blower and Dowlatabadi (1994)
  ########################################################################
  if( routine == "blower" ){
    
    #Do the calculation for each model output
    for( i in 1:k ){
      
      work.mat = cbind(par.rank,output.rank[,i])
      
      C.temp = matrix(0,nrow=p+1,ncol=p+1)
      
      #Calculate the C matrix
      for( j in 1:(p+1) ){
        for( l in 1:(p+1) ){
          
          C.temp[j,l]=sum((work.mat[,j]-ave.rank)*(work.mat[,l]-ave.rank))/
            sqrt(sum((work.mat[,j]-ave.rank)^2)*
                   sum((work.mat[,l]-ave.rank)^2))
        }
      }
      
      #Calculate the B matrix (qr helps with inversion)
      B.temp = solve(qr(C.temp))
      
      coeff.val = rep(0,p)
      
      #Calculate the PRCC
      for( j in 1:p ){
        
        coeff.val[j] = -B.temp[j,p+1]/sqrt(B.temp[j,j]*B.temp[p+1,p+1])
        
      }
      
      #Calculate the t-test statistics and p-values
      t.val = coeff.val*sqrt((n-2)/1-coeff.val)
      p.val = 2*pt(abs(t.val),df=(n-2),lower.tail=F)
      
      #Output the results
      results[[output.names[i]]] = data.frame(
        prcc = coeff.val,
        t.value = t.val,
        p.value = p.val,
        row.names = par.names)
    }
    
    return(results)
  }
  
  ########################################################################
  #Calculate the PRCCs using regression methods
  ########################################################################
  else if( routine == "regression" ){
    
    #Do the calculation for each model output
    for( i in 1:k ){
      
      coeff.val = rep(0,p)
      
      #Calculate the PRCC
      for( j in 1:p ){
        
        #Variation in output that can not be explained by all other predictors
        #(except the predictor of interest)
        fit.y = lm(output.rank[,i] ~ par.rank[,-j])
        
        #Variation in the predictor of interest that can not be explained
        #by the other predictors
        fit.x = lm(par.rank[,j] ~ par.rank[,-j])
        
        #PRCC is the correlation between the residuals of the two
        #regressions above
        coeff.val[j] = cor(fit.y$residuals,fit.x$residuals)
        
      }
      
      #Calculate the t-test statistics and p-values
      t.val = coeff.val*sqrt((n-2)/1-coeff.val)
      p.val = 2*pt(abs(t.val),df=(n-2),lower.tail=F)
      
      #Output the results
      results[[output.names[i]]] = data.frame(
        prcc = coeff.val,
        t.value = t.val,
        p.value = p.val,
        row.names = par.names)
    }
    
    return(results)
  }
  
  else{ return("Error: Calculation type is invalid. Must be either 'blower' or 'regression'") }
  
}


################################################################################
## plot.prcc - function to plot a prcc object
##
##    Arguments:
##          prcc.obj = a prcc object from the 'prcc' function
##             alpha = level of significance desired for cut-off lines
##               ... = additional arguments to be passed to functions called
##                     within this function
##
##
##    Output:
##      A plot that has a bar graph for each output variable giving the PRCCs.
##      Dashed red lines give the cutoff values for significance. It parameter
##      names are specified correctly, axis labels will be smart.

plot.prcc = function(prcc.obj,alpha=0.05,...){
  
 # x11()
  par(mfrow=c(ceiling((length(prcc.obj)-1)/3),min(c(length(prcc.obj)-1,3))))
  
  for( i in 2:length(prcc.obj) ){
    
    names.list=dimnames(results[[i]])[[1]]
    
    #Bar graph with the prcc values. The function in names.arg converts character
    #strings containing the names of greek letters to actual greek letters
    barplot(prcc.obj[[i]][,"prcc"],
            names.arg=sapply(names.list,
                             function(x) as.expression(substitute(list(a),list(a=as.symbol(x))))),
            main=names(results)[i],ylab="PRCC",cex.lab=1.1,...)
    
    #Plot lines to show the cutoff values for the alpha level of significance
    t.cutoff=qt(alpha/2,df=prcc.obj$num.sets-2)
    sig.cutoffs=c( (-(t.cutoff^2)-sqrt((t.cutoff^4) + 4*(prcc.obj$num.sets-2)*(t.cutoff^2)))/(2*(prcc.obj$num.sets-2)),
                   (-(t.cutoff^2)+sqrt((t.cutoff^4) + 4*(prcc.obj$num.sets-2)*(t.cutoff^2)))/(2*(prcc.obj$num.sets-2)))
    abline(h=sig.cutoffs,lty=2,col="red")
    abline(h=0,lty=1,col="black")
  }
}

# prcc with stoch simulation results...
res<-cbind(X[,1],X[,2],
           X[,3],
           X[,4], # pop extinction
           X[,5])
dimnames(res)[[2]]<-c("Population","Prevalence",
                      "Adult seroprevalence",
                      "Population persistence",
                      "LBV persistence")

write.csv(res, "results_25yr1000Sens_PRCC.csv", row.names=F, na="")

res<-cbind(#X[,1],#X[,2],
  #X[,3],
  #            X[,4], # pop extinction
  X[,5])

results=prcc(par.mat=hypercube,model.output=res ## results matrix here...
             ,routine="blower" # NB removed par names so uses symbols, add [par.names="",]
             ,output.names=c(#"Population size",
               #"Prevalence",
               #"Adult seroprevalence",
               #"Population persistence",
               "LBV persistence"))
tiff("prcc_lbv_pers.tiff",width=6,height=8,units='in',res=300, compression = "lzw")
plot.prcc(results,ylim=c(-1,1),cex.sub=1.5,cex.axis=1.25,cex.names=1.5)
dev.off()

## all works if enough paramets sets.... [otherwise get singular matrix warning..]

# NB testing K - carrying cap.
# no LHS needed

nonVarying = matrix(c(
  BETA=7.635,
  MU=0.000510492,
  DELTA=0.002312247,
  ALPHA=0.2,
  RHO=0.0556875,
  SIGMA=1/48,
  #K=1000000,
  EPSILON=1/365,
  TAU=1/24,
  KAPPA=1.5/365,
  PSI = 0.1,
  S=14.83,
  OMEGA=1/365,
  PHI=0.5,
  GAMMA=0.0037758,
  ETA=0.1,
  SUSJ.0=4000,MDAJ.0=4000, SUSJM.0=1000,EIJ.0=1000,ERJ.0=1000,INFJ.0=1000,
  RECJ.0=10000,SUSA.0=50000, EIA.0=100,
  ERA.0=1000,INFA.0=5000, RECA.0=50000,
  SPA.0=0.4994506,SPJ.0=0.5882353),
  ncol=29,
  nrow=100,
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
                            "PSI",
                            "S",
                            "OMEGA",
                            "PHI",
                            "GAMMA",
                            "ETA",
                            "SUSJ.0","MDAJ.0", "SUSJM.0","EIJ.0","ERJ.0","INFJ.0", "RECJ.0", "SUSA.0", "EIA.0","ERA.0","INFA.0", "RECA.0","SPA.0","SPJ.0") # naming non-varying columns

## from other code
Kset=seq(from = 100, to=1000000, by =10000)

fullParamSets = cbind(nonVarying,Kset) # full parameter set
#fullParamSets[,1] <- fullParamSets[,1]*fullParamSets[,29]
head(fullParamSets)
dim(fullParamSets)

dimnames(fullParamSets)[[2]]=c("BETA","MU","DELTA","ALPHA","RHO",
                               "SIGMA","EPSILON","TAU",
                               "KAPPA","PSI","S","OMEGA","PHI","GAMMA","ETA",
                               "SUSJ.0","MDAJ.0", "SUSJM.0","EIJ.0","ERJ.0","INFJ.0", "RECJ.0", "SUSA.0", "EIA.0","ERA.0","INFA.0", "RECA.0","SPA.0","SPJ.0",
                               "K")

# order for pomp/C model:  
BETA = fullParamSets[,1]
MU = fullParamSets[,2]
DELTA = fullParamSets[,3]
ALPHA = fullParamSets[,4]
RHO = fullParamSets[,5]
SIGMA = fullParamSets[,6]
K = fullParamSets[,30]
EPSILON = fullParamSets[,7]
TAU = fullParamSets[,8]
KAPPA = fullParamSets[,9]
PSI = fullParamSets[,10]
S = fullParamSets[,11]
OMEGA = fullParamSets[,12]
PHI = fullParamSets[,13]
GAMMA = fullParamSets[,14]
ETA = fullParamSets[,15]
SUSJ.0 = fullParamSets[,16]
MDAJ.0 = fullParamSets[,17]
SUSJM.0 = fullParamSets[,18]
EIJ.0 = fullParamSets[,19]
ERJ.0 = fullParamSets[,20]
INFJ.0 = fullParamSets[,21]
RECJ.0 = fullParamSets[,22]
SUSA.0 = fullParamSets[,23]
EIA.0 = fullParamSets[,24]
ERA.0 = fullParamSets[,25]
INFA.0 = fullParamSets[,26]
RECA.0 = fullParamSets[,27]
SPA.0 = fullParamSets[,28]
SPJ.0 = fullParamSets[,29]

paramset<-cbind(BETA,
                MU,
                DELTA,
                ALPHA,
                RHO,
                SIGMA,
                K,
                EPSILON,
                TAU,
                KAPPA,
                PSI,
                S,
                OMEGA,
                PHI,
                GAMMA,
                ETA,
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
                SPJ.0)
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
## change # sims
out1 <-simulate(sir,params=c(paramset[1,]),
                seed=1493885L,nsim=100,states=T,obs=F,as.data.frame=T) #

outres1 <- out1[seq(from=9126,to=912600,by=9126),] # select last #s
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
  outres <- out[seq(from=9126,to=912600,by=9126),] # select last #s
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
# below for K
# plot(X[,1],X[,2])
# plot
par(mfrow=c(1,1))
par(mar=c(5, 6, 4, 4) + 0.1)

tiff("k_lbv_prev.tiff",width=8,height=8,units='in',res=300, compression = "lzw")

plot(X[,1],X[,2],pch=16,
     ylab="Mean prevalence (%)",xlab="Population size",
     col="grey25", cex.lab=1.2)
dev.off()

tiff("k_lbv_serop.tiff",width=8,height=8,units='in',res=300, compression = "lzw")

plot(X[,1],X[,3],pch=16,
     ylab="Mean seroprevalence (%)",xlab="Population size",
     col="grey25", cex.lab=1.2)
dev.off()

tiff("k_lbv_pers.tiff",width=8,height=8,units='in',res=300, compression = "lzw")

plot(X[,1],X[,5],pch=16,
     ylab="P[persist]",xlab="Population size",
     col="grey25", cex.lab=1.2)
dev.off()
##
