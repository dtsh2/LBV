## fitting birth pulse data
## David Hayman 9/24/2013
## email davidtshayman@gmail.com
#######################################################
##

# set working directory - chose the folder that the results were saved to

rm(list=ls())
setwd("C:/Users/dtshayma/Documents/GitHub/LBV")

# getwd()
library(nlme)

dat <- read.csv("eidolonbirthdata.csv",header=T)
head(dat)
attach(dat)
# plot
# plot(Month,RateChange)

x=Month/12 #since function is set up so 1=year
y=RateChange
omega=1 # no pulses/yr
####start loop to do all 6 columns of dataframe
lfn=function(params){
    k=params["k"]
    s=params["s"]
    #omega=params["omega"]
    phi=params["phi"]
    ybar=k*(1/sqrt((1/s)*pi)*exp(-(cos(pi*omega*x-phi))^2/(1/s))) #this is the approximate delta function
    
    n=length(ybar)
    ss=sum((y-ybar)^2) # sum of squares- we want to minimize
    log(ss/n)
  }
  
  params=c(k=0.1,s=70,phi=0.5)

  lfn(params)
  
  out=optim(params, lfn,control=list(trace=2), method="L-BFGS-B",
            lower=c(0,0,0),upper=c(1,100,1),hessian=TRUE)
  #trace allows you to see what the computer is doing
  #hessian makes it save the slopes of the tries - this is important for CIs
  
  out
  
  k=out$par["k"]
  s=out$par["s"]
  phi=out$par["phi"]
  #omega=out$par["omega"]  # 
  
  #the prediction using the estimates
  pred=k*(1/sqrt((1/s)*pi)*exp(-(cos(pi*omega*x-phi))^2/(1/s)))
  
par(mfrow=c(1,2))
  plot(x,y,ylab="Proportion change in pregnancy",xlab="1 year",bty="n")
  lines(x,pred,col="red")
  
  text(x=0.7,y=max(y)-0.1,labels=paste(
" k=",signif(k,digits=4),
",\n s=",signif(s,digits=4),
",\n omega=",signif(omega,digits=4),
",\n phi=",signif(phi,digits=4),sep=""),cex=0.9,pos=1)

# x<-seq(0,1,by=0.01)
# pred=k*(1/sqrt((1/s)*pi)*exp(-(cos(pi*omega*x-phi))^2/(1/s)))
# plot(x,pred,type="l",xlab="time",ylab="")
# x<-Month/12
# points(x,y,pch=14) # appalling fit

# 
# b <-function(t,s=14.3,omega=1,phi=0,k=0.21){ ## adding birth pulse
#   kappa<-1/s
#   x<-cos(pi*omega*t-phi)
#   k*(1/sqrt(kappa*pi)*exp(-x^2/kappa))
# }
# 
# plot(b) ## no "pulsed" over the 0-1
# plot(b,xlab="",ylab="",xaxt="n",yaxt="n",ylim=c(0,0.5))
# # points(x,y,pch=16) # appalling fit
# 
# legend("topright",c("fit: k 1.5, s 14.3"),
#        col=c("black"),lty=1,bty="n")
# legend("topleft",c("data"),
#        col=c("black"),pch=16,bty="n")
# 
# ## NB s == controls synchrony - 
# ##    omega == # pulses
# ##    phi == position
# ##    k == heat of pulse
# ## t == time obviously

intBt=integrate((adf<-function(t,s=14.3,omega=1,phi=0,
                               k=1.51){
  epsilon<-1/s
  x<-cos(pi*omega*t-phi)
  k*(1/sqrt(epsilon*pi)*exp(-x^2/epsilon))
}),0,1)

intBt # to get value for integral of B(t) this should be annual birth rate

##
## plotting birth pulse

library(deSolve)

## this is the BP function to add....
# 
# b <-function(t,s=50,omega=2,phi=0,k=1.2){ ## adding birth pulse
#   kappa<-1/s
#   x<-cos(pi*omega*t-phi)
#   k*(1/sqrt(kappa*pi)*exp(-x^2/kappa))
# }
# 
# plot(b)#,add=T) ## no "pulsed" over the 0-1

# ss<-seq(20,70,by=5)
# test<-for (i in 1:length(ss)){
#   res<-integrate((adf<-function(t,s=ss[i],omega=2,phi=0,k=1.2){
#     epsilon<-1/s
#     x<-cos(pi*omega*t-phi)
#     k*(1/sqrt(epsilon*pi)*exp(-x^2/epsilon))
#   }),0,1)
#   res
# }
# 
# ## NB s == controls synchrony - 
# ##    omega == # pulses
# ##    phi == position
# ##    k == heat of pulse
# ## t == time obviously
# 
# intBt=integrate((adf<-function(t,s=30,omega=2,phi=0,k=1.2){
#   epsilon<-1/s
#   x<-cos(pi*omega*t-phi)
#   k*(1/sqrt(epsilon*pi)*exp(-x^2/epsilon))
# }),0,1)
# 
# intBt # to get value for integral of B(t) this should be annual birth rate
# 
# B.t = function(t,s,k,omega=2,phi=0) {
#   k*sqrt(s/pi)*exp(-s*cos(pi*t*omega-phi)^2)
# }
# 
# # with the birth rate set so that it matches expected values
# setK = function(t,s,omega=2,phi=0,br){
#   k=br/integrate(B.t,0,1,s,1,omega=2,phi=0)$value  
#   # then calculate B(t) with this calculated value of k
#   B.t(t,s,k,omega,phi)
# }
# 
# res<-setK(t=1,s=30,omega=2,phi=0,br=0.48)
# 
# b <-function(t=1,s=30,omega=2,phi=0,k=res){ ## adding birth pulse
#   kappa<-1/s
#   x<-cos(pi*omega*t-phi)
#   k*(1/sqrt(kappa*pi)*exp(-x^2/kappa))
# }
# plot(b)
# # Birth pulse function
 B.t = function(t,s,k,omega=2,phi=0) {
   k*sqrt(s/pi)*exp(-s*cos(pi*t*omega-phi)^2)
 }
# 
 # Function to calculate the birth pulse function, with the birth rate set so that it balances the death rate, m
 setB.t = function(t,s,m,omega=2,phi=0){
   k=m/integrate(B.t,0,1,s,1,omega,phi)$value
   k/365
 }
 
 ss<-seq(1.435,143.5,by=10)
 res<-matrix(NA,ncol=1,nrow=length(ss))
 for (i in 1:nrow(res)){
   res[i,]<-setB.t(t=1,s=ss[i],m=0.48)  
   res
 }
 res
# 
 par.plot<-cbind(res,ss,rep(2,length(ss)),rep(0,length(ss)))
 dimnames(par.plot)[[2]]<-c("k","s","omega","phi")
 par.plot
# 
# par(omi=c(1,1,0.5,0.5))
# par(mai=c(0.8,0.8,0.8,0.8))
# par(mar=c(1, 4, 1, 1) + 0.1)
# par(mfrow=c(1,1))
# 
# b <-function(t=1,s=par.plot[15,2],omega=2,phi=0,k=par.plot[15,1]){ ## adding birth pulse
#   kappa<-1/s
#   x<-cos(pi*omega*t-phi)
#   k*(1/sqrt(kappa*pi)*exp(-x^2/kappa))
# }
# plot(b,xlab="",ylab="b(t)",#xaxt="n",yaxt="n",
#      col="black",ylim=c(0,0.05),lwd=1.2,lty=2,bty="n")
# 
# b <-function(t=1,s=par.plot[1,2],omega=2,phi=0,k=par.plot[1,1]){ ## adding birth pulse
#   kappa<-1/s
#   x<-cos(pi*omega*t-phi)
#   k*(1/sqrt(kappa*pi)*exp(-x^2/kappa))
# }
# plot(b,add=T,col="black",lty=3,lwd=1.2)
# 
# b <-function(t=1,s=par.plot[15,2],omega=1,phi=0,k=par.plot[1,1]){ ## adding birth pulse
#   kappa<-1/s
#   x<-cos(pi*omega*t-phi)
#   k*(1/sqrt(kappa*pi)*exp(-x^2/kappa))
# }
# plot(b,add=T,col="darkgrey",lty=2,lwd=1.2)
# 
# 
# b <-function(t=1,s=par.plot[1,2],omega=1,phi=0,k=par.plot[1,1]){ ## adding birth pulse
#   kappa<-1/s
#   x<-cos(pi*omega*t-phi)
#   k*(1/sqrt(kappa*pi)*exp(-x^2/kappa))
# }
# 
# plot(b,add=T,col="darkgrey",lty=3,lwd=1.2)
# 
# ## values used
# 
# b <-function(t=1,s=14.35,omega=2,phi=0,k=4.1*10^-3){ ## adding birth pulse
#   kappa<-1/s
#   x<-cos(pi*omega*t-phi)
#   k*(1/sqrt(kappa*pi)*exp(-x^2/kappa))
# }
# plot(b,add=T,col="black",lty=1,lwd=1.2)
# 
# b <-function(t=1,s=14.35,omega=1,phi=0,k=4.1*10^-3){ ## adding birth pulse
#   kappa<-1/s
#   x<-cos(pi*omega*t-phi)
#   k*(1/sqrt(kappa*pi)*exp(-x^2/kappa))
# }
# plot(b,add=T,col="darkgrey",lty=1,lwd=1.2)
# 
# ##
# 
# legend(x=0.25,y=0.05,lty=rep(c(1:3),2),lwd=rep(1.2,6),
#        legend=c(expression(paste(omega==1, ", k = 0.0041, s = 14.35"),paste(omega==1, ", k = 0.0041, s = 143.5"),
#                            paste(omega==1, ", k = 0.0035, s = 1.435"),paste(omega==2, ", k = 0.0041, s = 14.35"),
#                            paste(omega==2, ", k = 0.0041, s = 143.5"),paste(omega==2, ", k = 0.0035, s = 1.43"))),
#        bty="n",col=c(rep(1,3),rep("darkgrey",3)),cex=0.8)
# 
# mtext("1 year",side=1,outer=F,line=3)
# 
###

## for lattice plot

ss<-seq(1.43,14.3,by=(14.3-1.43)/9)
res<-matrix(NA,ncol=1,nrow=length(ss))
for (i in 1:nrow(res)){
  res[i,]<-setB.t(t=1,s=ss[i],m=0.48)  
  res
}
res

par.plot<-cbind(res,ss,rep(2,length(ss)),rep(0,length(ss)))
dimnames(par.plot)[[2]]<-c("k","s","omega","phi")
par.plot

#par(omi=c(1,1,0.5,0.5))
# par(mai=c(0.8,0.8,0.8,0.8))
# par(mar=c(1, 4, 1, 1) + 0.1)
#par(mfrow=c(1,1))

b <-function(t=1,s=par.plot[1,2],omega=1,phi=0,k=par.plot[1,1]){ ## adding birth pulse
  kappa<-1/s
  x<-cos(pi*omega*t-phi)
  k*(1/sqrt(kappa*pi)*exp(-x^2/kappa))
}

plot(b,xlab="",ylab="b(t)",#xaxt="n",yaxt="n",
     col="black",ylim=c(0,0.01),lwd=1.2,lty=2,bty="n")

for (i in 2:length(par.plot[,2])){
  b<-function(t=1,s=par.plot[i,2],omega=1,phi=0,k=par.plot[i,1]){ ## adding birth pulse
    kappa<-1/s
    x<-cos(pi*omega*t-phi)
    k*(1/sqrt(kappa*pi)*exp(-x^2/kappa))
  }
  plot(b,add=T,lwd=1.2,lty=2)}

mtext("1 year",side=1,outer=F,line=3)

b <-function(t=1,s=par.plot[10,2],omega=1,phi=0,k=par.plot[10,1]){ ## adding birth pulse
  kappa<-1/s
  x<-cos(pi*omega*t-phi)
  k*(1/sqrt(kappa*pi)*exp(-x^2/kappa))
}

plot(b,add=T,lwd=2,col="black")

# 
# b <-function(t=1,s=par.plot[3,2],omega=1,phi=0,k=par.plot[3,1]){ ## adding birth pulse
#   kappa<-1/s
#   x<-cos(pi*omega*t-phi)
#   k*(1/sqrt(kappa*pi)*exp(-x^2/kappa))
# }
# 
# plot(b,add=T,lwd=2,col="orange")

###

