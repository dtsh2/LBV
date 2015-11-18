library(akima)

data(akima)

r<-list()
r[[1]]<-cbind(akima$x,akima$y,akima$z)
random.stuff <- (runif(prod(length(r[[1]][,3])), min = -0.00001, max = 0.0001))
r[[2]]<-r[[1]]#$z+random.stuff
r[[3]]<-r[[1]]#$z+2*random.stuff
r[[2]][,3]<-r[[2]][,3]+random.stuff
r[[3]][,3]<-r[[3]][,3]-random.stuff

zzg<-list()
for (i in 1:length(r)){
  
  zzg[[i]] <- interp(x=r[[i]][,1], #
                     y=r[[i]][,2], # 
                     z=r[[i]][,3],
                     duplicate=T)#,grid.len=c(50,50))#,span=0.1)
}

## now estimate the variance
zzg_v <-matrix(NA,nrow=nrow(zzg[[3]]$z),ncol=ncol(zzg[[3]]$z))

for (j in 1:nrow(zzg[[3]]$z)){
    for (k in 1:ncol(zzg[[3]]$z)){
      for (i in 1:length(r)){
        dat_v[i]<-c(zzg[[i]]$z[j,k])}
            zzg_v[j,k] <-var(dat_v)
}}


zzg_vp<-list(zzg[[3]]$x,zzg[[3]]$y,zzg_v)

image(zzg_v,ann=T)#,#xlim=c(0,max(results[,2])),ylim=c(0,max(results[,3])),
      ylab=rholab,xlab=betalab)
filled.contour(zzg_v,color = terrain.colors)

## find best values and the fit and plot them together

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

# select best parameter pairs
dat_ll<-matrix(NA,nrow=length(r),ncol=3)
sim_sp_t<-matrix(NA,nrow=length(r),ncol=6)

for (i in 1:length(r)){
dat_ll[i,]<-r[[i]][r[[i]][,3]==max(r[[i]][,3]),]
}
#
dat_ll[2,2]<-5

for (i in 1:nrow(dat_ll)){
  params <- c(
    BETA=dat_ll[i,2],
    RHO=dat_ll[i,1],
    SUSJ.0=4000,MDAJ.0=4000, SUSJM.0=1000,EIJ.0=1000,ERJ.0=1000,INFJ.0=1000,
    RECJ.0=10000,SUSA.0=50000, EIA.0=100,
    ERA.0=1000,INFA.0=5000, RECA.0=50000,
    SPA.0=0.4994506,SPJ.0=0.5882353)
#    print(params)
#
sim <- simulate(sir,params=c(params),seed=3593885L,
                nsim=1,states=T,obs=F,as.data.frame=T) # 
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

cr<-cor(c(simspa,simspj),c(DSPA,DSPJ))
sim_sp_t[i,1]<-cr
crt<-cor.test(c(simspa,simspj),c(DSPA,DSPJ))
sim_sp_t[i,2]<-crt$p.value
sim_sp_t[i,3]=dat_ll[i,2]
sim_sp_t[i,4]=dat_ll[i,1]
sim_sp_t[i,2]=dat_ll[i,3]
}
