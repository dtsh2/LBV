######################################################
# aged teeth sera anti-LBV analysis               ####
######################################################
	rm(list=ls())
## get data

setwd("~/GitHub/LBV") 

# setwd("~/Cambridge/CSU 2013/LBV model/lbvmodels/new_models_sept/deterministic")
## data

age<-read.table("Ageofinfection.txt",header=T)

#setwd("~/GitHub/LBV") 

head(age)
class(age)
attach(age)
Number<-rep(1,80)
new.data<-cbind(age,Number)
head(new.data)
data<-aggregate(new.data[,c(2,7,8)],by=list(Age=Age),sum)
data
class(data)
colnames(data)<-c("Age","sumd","Pos","Tot")
head(data)

tiff("dist_titres_aged_bats.tiff",width=8,height=8,units='in',res=300, compression = "lzw")
hist(Log.Titre.,breaks=100,#xlim=c(0.5,4),
     main="",xlab="Log(Reciprocal titre)",col="grey")#,add=T)
par(fig = c(0.6, 1, 0.4, 1), #mar=c(0,0,0,0),
    new=TRUE)
hist(RecTitre,breaks=1000,main="",xlab="Reciprocal titre")
dev.off()
## fit Muench's model, p90 Hens et al
model1<- glm(cbind(Tot-Pos,Pos)~-1+Age,data=data,family=binomial(link="log"))
summary(model1)

model2<- glm(cbind(Pos,Tot-Pos)~1,data=data, offset=log(Age), family=binomial(link="cloglog"))
summary(model2)
exp(coef(model2))
## FOI est 0.1255391

## fit Griffith's model with linear FOI and quadratic
model3<- glm(cbind(Tot-Pos,Pos)~-1+Age+I(Age^2),data=data, 
             family=binomial(link="log"))
summary(model3)

## fit Grenfell and Anderson's quadratic

model4<- glm(cbind(Tot-Pos,Pos)~-1+Age+I(Age^2)+I(Age^3),data=data, 
             family=binomial(link="log"))
summary(model4)

### Farrington's model
#farr<-function(alpha,beta,gamma){
#  p=1-exp((alpha/beta)*data$Age*exp(-beta*data$Age)
#          +(1/beta)*((alpha/beta)-gamma)*(exp(-beta*data$Age)-1)-gamma*data$Age)
#ll=data$Pos*log(p)+(data$Tot-data$Pos)*log(1-p)
## see alternative ll
#return(-sum(ll))
#}

#library(stats4)
#model5 = mle(farr,start=list(alpha=0.07,beta=0.1,gamma=0.03))
#summary(model5)
#
#foiFar=1-exp((coef(model5)[1]/coef(model5)[2])*data$Age*exp(-coef(model5)[2]*data$Age)
#        +(1/coef(model5)[2])*((coef(model5)[1]/coef(model5)[2])-coef(model5)[3])*(exp(-coef(model5)[2]*data$Age)-1)-coef(model5)[3]*data$Age)
#points(data$Age,foiFar,col="orange",lty=2,type="l")
#
AIC(model2,model3,model4)

## 
tiff("foi.tiff",width=8,height=12,units='in',res=300, compression = "lzw")
par(mfrow=c(2,1))
par(cex.axis=1.5)
par(cex.lab=1.5)
symbols(data$Age,(data$Pos/data$Tot),(data$Tot),
        bg="#0000FF0A",
        #fg="white",bg="red",
         xlim=c(-1,15),ylim=c(-0.2,1.2),
        ylab="Seroprevalence",xlab="Age")
#text(data$Age, (data$Pos/data$Tot), data$Tot, cex=1)
text(data$Age, (data$Pos/data$Tot), ".", cex=1)
#for(i in 1:13) points(i,1.5,cex=i)
symbols(-0.5,0.9,0.2,inches=F,
        bg="#0000FF0A",
        #fg="white",bg="red",
        add=T)
text(1,0.9, "size = 1", cex=1.5)

foi<-exp(coef(model2))
plot(data$Age,rep(foi,14),lty=2,type="l",xlim=c(-1,15),
     ylab="Force of infection",xlab="Age",ylim=c(0,0.3),lwd=1.5)
#legend(8,0,"constant force of infection",lty=2, cex=1,bty="n")

## to plot others
## but don't, less clear figure
(coef(model3))
foi2<--((coef(model3)[1])+2*(coef(model3)[2])*data$Age)
points(data$Age,foi2,col="red",lty=2,type="l",lwd=1.5)

(coef(model4))
foi3<--(((coef(model4))[1])+2*((coef(model4))[2])*data$Age+3*((coef(model4))[3])*data$Age^2)
points(data$Age,foi3,col="blue",lty=2,type="l",lwd=1.5)

abline(h=0,col="grey")

legend("topleft",c("Constant","Quadratic","Polynomial"),lty=2,col=c("black","red","blue"),bty="n",lwd=1.5,cex=1.5)

dev.off()