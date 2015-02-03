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
dev.off()
getwd()
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
sim <- simulate(sir,params=c(params),#seed=3593885L,
nsim=100,states=T,obs=F,as.data.frame=T) #
class(sir) # pomp object
class(sim) # data frame - change states, obs and data.frame if want pomp obj
with(sim,plot(time,SPA,type="n"))
for(i in unique(sim$sim))
with(sim[ sim$sim==i,],lines(time,SPA,col=adjustcolor("grey", alpha=0.5)))
test<-tapply(sim$SPA, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="red")
with(sim,plot(time,RECA,type="n"))
with(sim,plot(time,RECA,type="n"),xlim=c(7000,9000))
with(sim,plot(time,RECA,type="n"),xlim=c(7000,9000))
for(i in unique(sim$sim))
with(sim[ sim$sim==i,],points(time,RECA,col=adjustcolor("green", alpha=0.5)))
for(i in unique(sim$sim))
with(sim[ sim$sim==i,],points(time,SUSA,col=adjustcolor("blue", alpha=0.5)))
with(sim,plot(time,SUSA,type="n",xlim=c(8000,9000)))
with(sim,plot(time,SUSA,type="n",xlim=c(8000,9000),ylab="numbers"))
for(i in unique(sim$sim))
with(sim[ sim$sim==i,],points(time,RECA,col=adjustcolor("green", alpha=0.5)))
with(sim,plot(time,SUSA,type="n",xlim=c(8000,9000),ylab="numbers"))
for(i in unique(sim$sim))
with(sim[ sim$sim==i,],points(time,INFA,col=adjustcolor("grey", alpha=0.5)))
test<-tapply(sim$INFA, sim$time, mean, na.rm=T)
head(test)
tail(test)
nzmean <- function(x) {
zvals <- x==0
if (all(zvals)) 0 else mean(x[!zvals])
}
with(sim,plot(time,INFA,type="n",xlim=c(8000,9000),ylab="numbers"))
for(i in unique(sim$sim))
with(sim[ sim$sim==i,],points(time, INFA,col=adjustcolor("orange", alpha=0.5)))
test<-tapply(sim$INFA, sim$time, nzmean, na.rm=T)
nzmean <- function(x) {
if (all(x==0)) 0 else mean(x[x!=0])
}
nzmean <- function(x) {
if (all(x==0)) 0 else mean(x[x!=0])
}
test<-tapply(sim$INFA, sim$time, nzmean, na.rm=T)
with(sim,plot(time,INFA,type="n",xlim=c(8000,9000),ylim=c(500),ylab="numbers"))
with(sim,plot(time,INFA,type="n",xlim=c(8000,9000),ylim=c(0,500),ylab="numbers"))
for(i in unique(sim$sim))
with(sim[ sim$sim==i,],points(time, INFA,col=adjustcolor("orange", alpha=0.5)))
with(sim,plot(time,INFA,type="n",xlim=c(8200,9000),ylim=c(0,400),ylab="numbers"))
for(i in unique(sim$sim))
with(sim[ sim$sim==i,],points(time, INFA,col=adjustcolor("orange", alpha=0.5)))
test<-tapply(sim$INFA, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="yellow")
test<-tapply(sim$INFJ, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="darkred")
nzmean <- function(x) {
if (x==0) 0 else mean(x[x!=0])
}
test<-tapply(sim$INFJ, sim$time, nzmean, na.rm=T)
class(sim)
sim$INFA[sim$INFA==0] <- NA
test<-tapply(sim$INFJ, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="pink")
test<-tapply(sim$INFA, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="pink")
sim$INFA[sim$INFA==NA] <- 0
test<-tapply(sim$INFA, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="red")
with(sim,plot(time,INFA,type="n",xlim=c(8200,9000),ylim=c(0,400),ylab="numbers"))
for(i in unique(sim$sim))
with(sim[ sim$sim==i,],points(time, INFA,col=adjustcolor("orange", alpha=0.5)))
sim$INFJ[sim$INFJ==0] <- NA
with(sim,plot(time,INFA,type="n",xlim=c(8200,9000),ylim=c(0,400),ylab="numbers"))
for(i in unique(sim$sim))
with(sim[ sim$sim==i,],points(time, INFA,col=adjustcolor("orange", alpha=0.5)))
#sim$INFA[sim$INFA==0] <- NA
test<-tapply(sim$INFA, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="black")
for(i in unique(sim$sim))
with(sim[ sim$sim==i,],points(time, INFJ,col=adjustcolor("red", alpha=0.5)))
test<-tapply(sim$INFJ, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="red")
lines(1:length(test),test,col="white")
legend("topright",c("juvenile","adult"),col=c("orange","red"))
legend("topright",c("juvenile","adult"),col=c("orange","red"),pch=1)
legend("topright",c("juvenile","adult"),col=c("orange","red"),pch=16)
with(sim,plot(time,SUSA,type="n",xlim=c(8900,9000),ylab="numbers"))
for(i in unique(sim$sim))
sim$RECA[sim$RECA==0] <- NA
with(sim[ sim$sim==i,],points(time,RECA,col=adjustcolor("green", alpha=0.5)))
test<-tapply(sim$RECA, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="black")
for(i in unique(sim$sim))
sim$SUSA[sim$SUSA==0] <- NA
with(sim[ sim$sim==i,],points(time,SUSA,col=adjustcolor("blue", alpha=0.5)))
test<-tapply(sim$SUSA, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="white")
with(sim,plot(time,SUSA,type="n",xlim=c(8900,9000),ylab="numbers"))
for(i in unique(sim$sim))
sim$RECA[sim$RECA==0] <- NA
with(sim[ sim$sim==i,],points(time,RECA,col=adjustcolor("green", alpha=0.5)))
test<-tapply(sim$RECA, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="black")
for(i in unique(sim$sim))
sim$SUSA[sim$SUSA==0] <- NA
with(sim[ sim$sim==i,],points(time,SUSA,col=adjustcolor("blue", alpha=0.5)))
test<-tapply(sim$SUSA, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="white")
with(sim,plot(time,SUSJ,type="n",xlim=c(8000,9000),ylab="numbers"))
for(i in unique(sim$sim))
with(sim[ sim$sim==i,],points(time,RECJ,col=adjustcolor("green", alpha=0.5)))
with(sim,plot(time,SUSJ,type="n",xlim=c(8990,9000),ylab="numbers"))
for(i in unique(sim$sim))
with(sim[ sim$sim==i,],points(time,RECJ,col=adjustcolor("darkgreen", alpha=0.5)))
lines(1:length(test),test,col="black")
test<-tapply(sim$RECJ, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="black")
with(sim[ sim$sim==i,],points(time,MDA,col=adjustcolor("lightbrown", alpha=0.5)))
sim$MDA[sim$MDA==0] <- NA
with(sim[ sim$sim==i,],points(time,MDA,col=adjustcolor("lightbrown", alpha=0.5)))
with(sim[ sim$sim==i,],points(time,MDA,col=adjustcolor("lightbrown", alpha=0.5)))
with(sim[ sim$sim==i,],points(time,MDA,col=adjustcolor("grey", alpha=0.5)))
test<-tapply(sim$MDA, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="black")
with(sim[ sim$sim==i,],points(time,SUSJ,col=adjustcolor("blue", alpha=0.5)))
test<-tapply(sim$SUSJ, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="white")
test<-tapply(sim$SUSJ, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="white")
plot(1:length(test),test,col="white")
plot(1:length(test),test,col="blue")
test<-tapply(sim$RECJ, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="black")
test<-tapply(sim$MDA, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="black")
test<-tapply(sim$RECJ, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="black")
test<-tapply(sim$MDA, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="black")
test<-tapply(sim$SUSJ, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="white")
with(sim,plot(time,SUSJ,type="n",xlim=c(8990,9000),ylab="numbers"))
test<-tapply(sim$RECJ, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="black")
test<-tapply(sim$MDA, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="grey")
test<-tapply(sim$SUSJ, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="white")
lines(1:length(test),test,col="blue")
with(sim,plot(time,SUSJ,type="n",xlim=c(0,9000),ylab="numbers"))
#for(i in unique(sim$sim))
#  with(sim[ sim$sim==i,],points(time,RECJ,col=adjustcolor("darkgreen", alpha=0.5)))
#sim$RECJ[sim$RECJ==0] <- NA
test<-tapply(sim$RECJ, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="black")
#for(i in unique(sim$sim))
#  sim$MDA[sim$MDA==0] <- NA
#with(sim[ sim$sim==i,],points(time,MDA,col=adjustcolor("grey", alpha=0.5)))
test<-tapply(sim$MDA, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="grey")
#for(i in unique(sim$sim))
#  sim$SUSJ[sim$SUSJ==0] <- NA
#with(sim[ sim$sim==i,],points(time,SUSJ,col=adjustcolor("blue", alpha=0.5)))
test<-tapply(sim$SUSJ, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="blue")
with(sim,plot(time,SUSJ,type="n",xlim=c(8000,9000),ylab="numbers"))
#for(i in unique(sim$sim))
#  with(sim[ sim$sim==i,],points(time,RECJ,col=adjustcolor("darkgreen", alpha=0.5)))
#sim$RECJ[sim$RECJ==0] <- NA
test<-tapply(sim$RECJ, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="black")
#for(i in unique(sim$sim))
#  sim$MDA[sim$MDA==0] <- NA
#with(sim[ sim$sim==i,],points(time,MDA,col=adjustcolor("grey", alpha=0.5)))
test<-tapply(sim$MDA, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="grey")
#for(i in unique(sim$sim))
#  sim$SUSJ[sim$SUSJ==0] <- NA
#with(sim[ sim$sim==i,],points(time,SUSJ,col=adjustcolor("blue", alpha=0.5)))
test<-tapply(sim$SUSJ, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="blue")
with(sim,plot(time,SUSA,type="n",xlim=c(8900,9000),ylab="numbers"))
#for(i in unique(sim$sim))
sim$RECA[sim$RECA==0] <- NA
#with(sim[ sim$sim==i,],points(time,RECA,col=adjustcolor("green", alpha=0.5)))
test<-tapply(sim$RECA, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="black")
#for(i in unique(sim$sim))
sim$SUSA[sim$SUSA==0] <- NA
#with(sim[ sim$sim==i,],points(time,SUSA,col=adjustcolor("blue", alpha=0.5)))
test<-tapply(sim$SUSA, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="blue")
with(sim,plot(time,SUSJ,type="n",xlim=c(8000,9000),ylab="numbers"))
#for(i in unique(sim$sim))
#  with(sim[ sim$sim==i,],points(time,RECJ,col=adjustcolor("darkgreen", alpha=0.5)))
sim$RECJ[sim$RECJ==0] <- NA
test<-tapply(sim$RECJ, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="black")
#for(i in unique(sim$sim))
sim$MDA[sim$MDA==0] <- NA
#with(sim[ sim$sim==i,],points(time,MDA,col=adjustcolor("grey", alpha=0.5)))
test<-tapply(sim$MDA, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="grey")
#for(i in unique(sim$sim))
sim$SUSJ[sim$SUSJ==0] <- NA
#with(sim[ sim$sim==i,],points(time,SUSJ,col=adjustcolor("blue", alpha=0.5)))
test<-tapply(sim$SUSJ, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="blue")
with(sim,plot(time,SUSJ,type="n",xlim=c(6000,9000),ylab="numbers"))
#for(i in unique(sim$sim))
#  with(sim[ sim$sim==i,],points(time,RECJ,col=adjustcolor("darkgreen", alpha=0.5)))
sim$RECJ[sim$RECJ==0] <- NA
test<-tapply(sim$RECJ, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="black")
#for(i in unique(sim$sim))
sim$MDA[sim$MDA==0] <- NA
#with(sim[ sim$sim==i,],points(time,MDA,col=adjustcolor("grey", alpha=0.5)))
test<-tapply(sim$MDA, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="grey")
#for(i in unique(sim$sim))
sim$SUSJ[sim$SUSJ==0] <- NA
#with(sim[ sim$sim==i,],points(time,SUSJ,col=adjustcolor("blue", alpha=0.5)))
test<-tapply(sim$SUSJ, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="blue")
sim$RECA[sim$RECA==0] <- NA
#with(sim[ sim$sim==i,],points(time,RECA,col=adjustcolor("green", alpha=0.5)))
test<-tapply(sim$RECA, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="black")
#for(i in unique(sim$sim))
sim$SUSA[sim$SUSA==0] <- NA
#with(sim[ sim$sim==i,],points(time,SUSA,col=adjustcolor("blue", alpha=0.5)))
test<-tapply(sim$SUSA, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="blue")
with(sim,plot(time,SUSA,type="n",xlim=c(8900,9000),ylab="numbers"))
#for(i in unique(sim$sim))
sim$RECA[sim$RECA==0] <- NA
#with(sim[ sim$sim==i,],points(time,RECA,col=adjustcolor("green", alpha=0.5)))
test<-tapply(sim$RECA, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="black")
#for(i in unique(sim$sim))
sim$SUSA[sim$SUSA==0] <- NA
#with(sim[ sim$sim==i,],points(time,SUSA,col=adjustcolor("blue", alpha=0.5)))
test<-tapply(sim$SUSA, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="blue")
#for(i in unique(sim$sim))
sim$SUSA[sim$SUSA==0] <- NA
#with(sim[ sim$sim==i,],points(time,SUSA,col=adjustcolor("blue", alpha=0.5)))
lines(1:length(test),test,col="blue")
test<-tapply(sim$SUSA, sim$time, mean, na.rm=T)
#dev.off()
##
#tiff("juv_points.tiff",width=8,height=8,units='in',res=300, compression = "lzw")
#with(sim,plot(time,SUSJ,type="n",xlim=c(6000,9000),ylab="numbers"))
#for(i in unique(sim$sim))
#  with(sim[ sim$sim==i,],points(time,RECJ,col=adjustcolor("darkgreen", alpha=0.5)))
sim$RECJ[sim$RECJ==0] <- NA
test<-tapply(sim$RECJ, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="black")
#for(i in unique(sim$sim))
sim$MDA[sim$MDA==0] <- NA
#with(sim[ sim$sim==i,],points(time,MDA,col=adjustcolor("grey", alpha=0.5)))
test<-tapply(sim$MDA, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="grey")
#for(i in unique(sim$sim))
sim$SUSJ[sim$SUSJ==0] <- NA
#with(sim[ sim$sim==i,],points(time,SUSJ,col=adjustcolor("blue", alpha=0.5)))
test<-tapply(sim$SUSJ, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="blue")
with(sim,plot(time,SUSA,type="n",xlim=c(6000,9000),ylab="numbers"))
#for(i in unique(sim$sim))
sim$RECA[sim$RECA==0] <- NA
#with(sim[ sim$sim==i,],points(time,RECA,col=adjustcolor("green", alpha=0.5)))
test<-tapply(sim$RECA, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="black")
#for(i in unique(sim$sim))
sim$SUSA[sim$SUSA==0] <- NA
#with(sim[ sim$sim==i,],points(time,SUSA,col=adjustcolor("blue", alpha=0.5)))
test<-tapply(sim$SUSA, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="blue")
#dev.off()
##
#tiff("juv_points.tiff",width=8,height=8,units='in',res=300, compression = "lzw")
#with(sim,plot(time,SUSJ,type="n",xlim=c(6000,9000),ylab="numbers"))
#for(i in unique(sim$sim))
#  with(sim[ sim$sim==i,],points(time,RECJ,col=adjustcolor("darkgreen", alpha=0.5)))
sim$RECJ[sim$RECJ==0] <- NA
test<-tapply(sim$RECJ, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="black")
#for(i in unique(sim$sim))
sim$MDA[sim$MDA==0] <- NA
#with(sim[ sim$sim==i,],points(time,MDA,col=adjustcolor("grey", alpha=0.5)))
test<-tapply(sim$MDA, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="grey")
#for(i in unique(sim$sim))
sim$SUSJ[sim$SUSJ==0] <- NA
#with(sim[ sim$sim==i,],points(time,SUSJ,col=adjustcolor("blue", alpha=0.5)))
test<-tapply(sim$SUSJ, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="blue")
legend("topright",c("recovered adult","susceptible adult","recovered juvenile","maternal immunity"," susceptibe juvenile"),
col=c("black","blue","black","grey","blue"),lty=c(1,1,1,3,3))
legend("topright",c("recovered adult","susceptible adult","recovered juvenile","maternal immunity"," susceptibe juvenile"),
col=c("black","blue","black","grey","blue"),lty=c(1,1,1,3,3),bty="n")
with(sim,plot(time,SUSA,type="n",xlim=c(7000,9000),ylab="numbers"))
#for(i in unique(sim$sim))
sim$RECA[sim$RECA==0] <- NA
#with(sim[ sim$sim==i,],points(time,RECA,col=adjustcolor("green", alpha=0.5)))
test<-tapply(sim$RECA, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="black")
#for(i in unique(sim$sim))
sim$SUSA[sim$SUSA==0] <- NA
#with(sim[ sim$sim==i,],points(time,SUSA,col=adjustcolor("blue", alpha=0.5)))
test<-tapply(sim$SUSA, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="blue")
#dev.off()
##
#tiff("juv_points.tiff",width=8,height=8,units='in',res=300, compression = "lzw")
#with(sim,plot(time,SUSJ,type="n",xlim=c(6000,9000),ylab="numbers"))
#for(i in unique(sim$sim))
#  with(sim[ sim$sim==i,],points(time,RECJ,col=adjustcolor("darkgreen", alpha=0.5)))
sim$RECJ[sim$RECJ==0] <- NA
test<-tapply(sim$RECJ, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="black",lty=3)
#for(i in unique(sim$sim))
sim$MDA[sim$MDA==0] <- NA
#with(sim[ sim$sim==i,],points(time,MDA,col=adjustcolor("grey", alpha=0.5)))
test<-tapply(sim$MDA, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="grey",lty=3)
#for(i in unique(sim$sim))
sim$SUSJ[sim$SUSJ==0] <- NA
#with(sim[ sim$sim==i,],points(time,SUSJ,col=adjustcolor("blue", alpha=0.5)))
test<-tapply(sim$SUSJ, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="blue",lty=3)
legend("topright",c("recovered adult","susceptible adult","recovered juvenile","maternal immunity"," susceptibe juvenile"),
col=c("black","blue","black","grey","blue"),lty=c(1,1,1,3,3),bty="n")
test<-tapply(sim$MDA, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="green",lty=3)
legend("topright",c("recovered adult","susceptible adult","recovered juvenile","maternal immunity"," susceptibe juvenile"),
col=c("black","blue","black","green","blue"),lty=c(1,1,1,3,3),bty="n")
legend("top",c("recovered adult","susceptible adult","recovered juvenile","maternal immunity"," susceptibe juvenile"),
col=c("black","blue","black","green","blue"),lty=c(1,1,1,3,3),bty="n")
legend("top",c("recovered adult","susceptible adult","recovered juvenile","maternal immunity"," susceptibe juvenile"),
col=c("black","blue","black","green","blue"),lty=c(1,1,1,3,3),bty="n",cex=0.4)
with(sim,plot(time,SUSA,type="n",xlim=c(7000,9000),ylab="numbers"))
#for(i in unique(sim$sim))
sim$RECA[sim$RECA==0] <- NA
#with(sim[ sim$sim==i,],points(time,RECA,col=adjustcolor("green", alpha=0.5)))
test<-tapply(sim$RECA, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="black")
#for(i in unique(sim$sim))
sim$SUSA[sim$SUSA==0] <- NA
#with(sim[ sim$sim==i,],points(time,SUSA,col=adjustcolor("blue", alpha=0.5)))
test<-tapply(sim$SUSA, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="blue")
#dev.off()
##
#tiff("juv_points.tiff",width=8,height=8,units='in',res=300, compression = "lzw")
#with(sim,plot(time,SUSJ,type="n",xlim=c(6000,9000),ylab="numbers"))
#for(i in unique(sim$sim))
#  with(sim[ sim$sim==i,],points(time,RECJ,col=adjustcolor("darkgreen", alpha=0.5)))
sim$RECJ[sim$RECJ==0] <- NA
test<-tapply(sim$RECJ, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="black",lty=3)
#for(i in unique(sim$sim))
sim$MDA[sim$MDA==0] <- NA
#with(sim[ sim$sim==i,],points(time,MDA,col=adjustcolor("grey", alpha=0.5)))
test<-tapply(sim$MDA, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="green",lty=3)
#for(i in unique(sim$sim))
sim$SUSJ[sim$SUSJ==0] <- NA
#with(sim[ sim$sim==i,],points(time,SUSJ,col=adjustcolor("blue", alpha=0.5)))
test<-tapply(sim$SUSJ, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="blue",lty=3)
legend("top",c("recovered adult","susceptible adult","recovered juvenile","maternal immunity"," susceptibe juvenile"),
col=c("black","blue","black","green","blue"),lty=c(1,1,1,3,3),bty="n",cex=0.9)
?legend
legend("topleft",c("recovered adult","susceptible adult","recovered juvenile","maternal immunity","susceptibe juvenile"),
col=c("black","blue","black","green","blue"),lty=c(1,1,3,3,3),bty="n",cex=0.9,fill="white")
legend("topleft",c("recovered adult","susceptible adult","recovered juvenile","maternal immunity","susceptibe juvenile"),
col=c("black","blue","black","green","blue"),lty=c(1,1,3,3,3),
bty="n",cex=0.9,bg="white")
dev.off()
with(sim,plot(time,SUSA,type="n",xlim=c(7000,9000),ylab="numbers"))
sim$SUSA[sim$SUSA==0] <- NA
#with(sim[ sim$sim==i,],points(time,SUSA,col=adjustcolor("blue", alpha=0.5)))
test<-tapply(sim$SUSA, sim$time, mean, na.rm=T)
lines(1:length(test),test,col="blue")
legend("topleft",c("recovered adult","susceptible adult","recovered juvenile","maternal immunity","susceptibe juvenile"),
col=c("black","blue","black","green","blue"),lty=c(1,1,3,3,3),
bty="n",cex=0.9,bg="white")
legend("topleft",c("recovered adult","susceptible adult","recovered juvenile","maternal immunity","susceptibe juvenile"),
col=c("black","blue","black","green","blue"),lty=c(1,1,3,3,3),
cex=0.9,bg="white")
legend("topleft",c("recovered adult","susceptible adult","recovered juvenile","maternal immunity","susceptibe juvenile"),
col=c("black","blue","black","green","blue"),lty=c(1,1,3,3,3),
cex=0.9,bg="white",box.col="white")
box(col="black")
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
savehistory("~/GitHub/LBV/hist_lbv_020315.R")
