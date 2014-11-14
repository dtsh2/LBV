result1<-rep(1:3,3)
result2<-rep(1:3,3)
result3<-rep(1:3,3)
pars<-expand.grid(result1,result2,result3)
dim(pars)
res<-rnorm(length(pars[,1]))
results<-cbind(pars,res)

results$Var1[which(results[,4]==min(results[,4]))]
results$Var2[which(results[,4]==min(results[,4]))]
results$Var3[which(results[,4]==min(results[,4]))]

parsplot <- expand.grid(result1,result2)
parsplotres<-cbind(parsplot,result1)*NA


library(data.table)

d <- data.table(results)
rest<-d[, min(res, na.rm=TRUE), by=c("Var1","Var2")]

library(akima)
library(lattice)
library(tgp)
library(rgl)
library(fields)

rholab<-expression(symbol(rho))
betalab<-expression(symbol(beta))

zzg <- interp(rest$Var1,rest$Var2,rest$V1)

image(zzg,ann=T,ylab=rholab,xlab=betalab)
contour(zzg,add=T,labcex=1,drawlabels=T,nlevels=10)
