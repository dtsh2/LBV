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

lbv.new<-cbind(times,DSPJ,DSPA,DRECJ,DRECA,DPOPJ,DPOPA)

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
    ## PACKAGE="pomp"  
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
