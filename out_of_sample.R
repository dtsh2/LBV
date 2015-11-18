# system("R CMD SHLIB lbvseirNoSeasFreqTwoParsMeasureBinom.c")

dyn.load("lbvseirNoSeasFreqTwoParsMeasureBinom.dll")
# dyn.unload("lbvseirNoSeasFreqTwoParsMeasureBinom.dll")

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
#

for (i in 1:length(lbv.new[,1])){
  lbv.n <- lbv.new[-i,]
  
pomp(
  data = data.frame(
    time=lbv.n[,1],  # time for simulations to run
    #  DSPJ = lbv.n[,2],
    #  DSPA = lbv.n[,3],
    DRECJ = lbv.n[,4],
    DRECA = lbv.n[,5],
    DPOPJ = lbv.n[,6],
    DPOPA = lbv.n[,7]
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

#plot(lbvdat)
#########
# 
# pf<-pfilter(lbvdat,params=c(params),Np=100,max.fail=100,tol=1e-20)
# logLik(pf)
# coef(pf)

BetaV = seq(from=0.001,to=15,by=2.5)  # range of beta
RhoV = seq(from=0.001,to=0.3, by=0.0625) # range of rho
#
BetaV = seq(from=0.001,to=20,by=1.25)  # range of beta
RhoV = seq(from=0.001,to=0.5, by=0.03125) # range of rho

#
parametset<- expand.grid(BetaV,RhoV)
dim(parametset)

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
#
write.csv(results,file=paste("mll_surface_sample_",i,".csv",sep="")) # file=paste(x, ".mean", sep="")
#
}

r<-list()

for (i in 1:length(lbv.new[,1])){
  file_n=paste("mll_surface_sample_",i,".csv",sep="")
  r[[i]]<-read.csv(file=file_n,header=T)
  print(r)
 }


library(akima)
library(lattice)
library(tgp)
library(rgl)
library(fields)

rholab<-expression(symbol(rho))
betalab<-expression(symbol(beta))

zzg<-list()
for (i in 1:length(lbv.new[,1])){

zzg[[i]] <- interp(x=r[[i]][,3], #
              y=r[[i]][,4], # 
              z=r[[i]][,2],
              duplicate=T)#,grid.len=c(50,50))#,span=0.1)
}
## narrow figure
#plot(results[,1])
# max(results[,1])
# results[results[,1]==max(results[,1]),]
# maxLL<-as.data.frame(results[results[,1]==max(results[,1]),])
# names(maxLL)<-c("negll","Beta","Rho")
# 
# tiff("ll_beta_rho.tiff",width=8,height=8,units='in',res=300, compression = "lzw")

for (i in 1:length(lbv.new[,1])){
  tiff(paste("mll_surface_sample_",i,".tiff",sep=""),width=8,height=8,units='in',res=300, compression = "lzw")  
image(zzg[[i]],ann=T,#xlim=c(0,max(results[,2])),ylim=c(0,max(results[,3])),
      ylab=rholab,xlab=betalab)
dev.off()
}
#contour(zzg,add=T,labcex=1,drawlabels=T,nlevels=10)
#contour(zzg,add=F,labcex=1,drawlabels=T,nlevels=100)
# points(x=maxLL$Beta,y=maxLL$Rho,pch=16,col="black")
# points(x=maxLL[2,],y=maxLL[3,],pch=16,col="black")
# dev.off()

## now estimate the variance
zzg_v <-matrix(NA,nrow=nrow(zzg[[3]]$z),ncol=ncol(zzg[[3]]$z))
dat_v<-vector()
for (j in 1:nrow(zzg[[3]]$z)){
  for (k in 1:ncol(zzg[[3]]$z)){
    for (i in 1:length(r)){
      dat_v[i]<-c(zzg[[i]]$z[j,k])}
    zzg_v[j,k] <-sd(dat_v)
  }}

tiff(paste("mll_sd.tiff",sep=""),width=8,height=8,units='in',res=300, compression = "lzw")  
image(zzg_v,ann=T,ylab=rholab,xlab=betalab, xaxt='n',yaxt='n')
axis(side=1, at=BetaV/max(BetaV),
     labels=BetaV)
axis(side=2, at=RhoV/max(RhoV),
     labels=RhoV)
dev.off()

tiff(paste("mll_sd_contour.tiff",sep=""),width=8,height=6,units='in',res=300, compression = "lzw")  
par(mar=c(5,7,4,2)+0.1)
filled.contour(zzg_v,color = terrain.colors,
              # ylab=rholab,xlab=betalab, 
               plot.axes={c(axis(side=1, at=BetaV/max(BetaV),
                                labels=BetaV),
                           axis(side=2, at=RhoV/max(RhoV),
                                              labels=RhoV))},
               plot.title={
                 title(xlab=betalab,cex.lab=1.5)
                 mtext(rholab,2,cex=1.5,line=5,las=1)})
dev.off()

## find best values and the fit and plot them together

# dyn.load("lbvseirNoSeasFreqTwoParsMeasureBinom.dll")
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

# select best parameter pairs
dat_ll<-matrix(NA,nrow=length(r),ncol=4)
sim_sp_t<-matrix(NA,nrow=length(r),ncol=7)

for (i in 1:length(r)){
  dat_ll[i,]<-as.numeric(r[[i]][r[[i]][,2]==max(r[[i]][,2]),])
}
#

for (i in 1:nrow(dat_ll)){
  params <- c(
    BETA=dat_ll[i,3],
    RHO=dat_ll[i,4],
    SUSJ.0=4000,MDAJ.0=4000, SUSJM.0=1000,EIJ.0=1000,ERJ.0=1000,INFJ.0=1000,
    RECJ.0=10000,SUSA.0=50000, EIA.0=100,
    ERA.0=1000,INFA.0=5000, RECA.0=50000,
    SPA.0=0.4994506,SPJ.0=0.5882353)
  #    print(params)
  #
  sim <- simulate(sir,params=c(params),seed=3593885L,
                  nsim=100,states=T,obs=F,as.data.frame=T) # 
  #########################################################
  ## multiple sims...
  #
  # sim <- simulate(sir,params=c(params),#seed=3593885L,
  #                 nsim=100,states=T,obs=F,as.data.frame=T) # 
  test<-tapply(sim$SPA, sim$time, mean, na.rm=T)
  simspa<-test[c(7301,7666,7725,7756,7970,8031,8090,8121,8212,8335,8396,8455)]
  test<-tapply(sim$SPJ, sim$time, mean, na.rm=T)
  simspj<-test[c(7301,7666,7725,7756,7970,8031,8090,8121,8212,8335,8396,8455)]
  
  # times<-lbvd$cumulative_time
  #
  # lbv.new<-cbind(times,DSPJ,DSPA,DRECJ,DRECA,DPOPJ,DPOPA)
  
#   crt<-cor.test(c(simspa,simspj),c(DSPA,DSPJ))
#   cr<-cor(c(simspa,simspj),c(DSPA,DSPJ))
  sim_sp_t[i,1]<-simspa[i]
  sim_sp_t[i,2]<-simspj[i]
  sim_sp_t[i,3]<-DSPA[i]
  sim_sp_t[i,4]<-DSPJ[i]
  sim_sp_t[i,5]=dat_ll[i,3] # beta
  sim_sp_t[i,6]=dat_ll[i,4] # rho
  sim_sp_t[i,7]=dat_ll[i,2] # ll
}

colnames(sim_sp_t)<-c("Sim-SPA","Sim-SPJ","Data-SPA","Data-SPJ",'Beta',"Rho","MLL")
write.csv(sim_sp_t,file=paste("out_of_sample.csv",sep="")) #

tiff(paste("pred_vs_data.tiff",sep=""),width=8,height=8,units='in',res=300, compression = "lzw")  
plot(c(sim_sp_t[,1],sim_sp_t[,2]),c(sim_sp_t[,3],sim_sp_t[,4]),
     ylim=c(0,0.6),xlim=c(0,0.4),col=c(rep("red",12),rep("blue",12)),
     pch=16,ylab="Data",xlab="Simulation")
crt<-cor.test(c(sim_sp_t[,1],sim_sp_t[,2]),c(sim_sp_t[,3],sim_sp_t[,4]))
cr<-cor(c(sim_sp_t[,1],sim_sp_t[,2]),c(sim_sp_t[,3],sim_sp_t[,4]))
vall <- format(cr,digits=3)
vallp <- format(crt$p.value,digits=2)
cr_res<-bquote(r == .(vall))
text(x=0.3,y=0.2,cr_res,
     cex = 1)
cr_res_p<-bquote(p-value == .(vallp))
text(x=0.3,y=0.1,cr_res_p,
     cex = 1)
lmp<-lm(c(sim_sp_t[,3],sim_sp_t[,4])~c(sim_sp_t[,1],sim_sp_t[,2]))
abline(lmp)
legend('topright',c("Adult","Juvenile"),col=c("red","blue"),pch=16)
dev.off()