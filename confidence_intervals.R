######## R code from lecture 11 #########

### profile likelihood CI for a two-parameter log-likelihod ###

#use aphid data to fit NB model 
num.stems<-c(6,8,9,6,6,2,5,3,1,4)
#generate raw data from tabulated values
aphid.data<-rep(0:9,num.stems)
NB.LL<-function(mu,theta) sum(log(dnbinom(aphid.data,mu=mu,size=theta)))
#for nlm we need negative LL and a function of a vector
negNB.LL<-function(p) -NB.LL(p[1],p[2])
#obtain MLEs
nlm(negNB.LL,c(4,4))->out.NB

#generate grid of values
mu<-seq(2,5,.1)
theta<-seq(1,9,.1)
g<-expand.grid(mu,theta)

#log-likelihood functions
NB.LL<-function(mu,theta) sum(log(dnbinom(aphid.data,mu=mu,size=theta)))
NB.LL.p<-function(p) sum(log(dnbinom(aphid.data,mu=p[1],size=p[2])))

#generate z-coordinate on grid
g$z<-apply(g,1,NB.LL.p)

#reorganize z as a matrix same shape as grid
zmat<-matrix(g$z,nrow=length(seq(2,5,.1)))
contour(mu,theta,zmat,xlab=expression(mu),ylab=expression(theta))

#add MLE to graph
points(out.NB$estimate[1],out.NB$estimate[2],pch=16,col=2,cex=1.2)

#lower bound for joint LR confidence interval
lowerlimit2<- -out.NB$minimum - .5*qchisq(.95,2)
lowerlimit2

#lower bound for profile likelihood confidence intervals
lowerlimit1<- -out.NB$minimum - .5*qchisq(.95,1)

#add the lower limits as contours to contour plot
contour(mu,theta,zmat,xlab=expression(mu),ylab=expression(theta),levels=c(lowerlimit2,lowerlimit1),lwd=c(2,1))
points(out.NB$estimate[1],out.NB$estimate[2],pch=16,col=2,cex=1.2)

#add bounds for profile likelhood confidence intervals
abline(v=out.NB$estimate[1],lty=2,col=2)
abline(h=out.NB$estimate[2],lty=2,col=2)

#calculate profile likelihood confidence intervals using Bhat package
library(Bhat)
x.in<-list(label=c('mu','theta'),est=out.NB$estimate,low=c(2,1),high=c(5,9))
plkhci(x.in,negNB.LL,'mu')->profile.mu
profile.mu
plkhci(x.in,negNB.LL,'theta')->profile.theta
profile.theta

#add these intervals to the graph
abline(h=profile.theta,col=4,lty=2)
abline(v=profile.mu,col=4,lty=2)

#calculate Wald confidence intervals for comparison
out.NB$estimate
#extract Hessian
nlm(negNB.LL,c(4,4),hessian=T)->out.NB
#calculate confidence limits
lower.bound<-out.NB$estimate+qnorm(.025)*sqrt(diag(solve(out.NB$hessian)))
upper.bound<-out.NB$estimate+qnorm(.975)*sqrt(diag(solve(out.NB$hessian)))
rbind(lower.bound,upper.bound)
#compare to profile likelihood CI
profile.mu
profile.theta

### Poisson and negative binomial regression

slugs<-read.table('http://www.bio.ic.ac.uk/research/mjcraw/statcomp/data/slugsurvey.txt',header=T)
slugtable<-data.frame(table(slugs$slugs,slugs$field))

#fit NB and Poisson models using nlm function
negNB.LL<-function(p){
  mu<-p[1]
  theta<-p[2]
  LL<-sum(log(dnbinom(slugs$slugs,mu=mu,size=theta)))
  -LL
}

negNB.LL1<-function(p){
  z<-as.numeric(slugs$field)-1
  mu<-p[1]+p[2]*z
  theta<-p[3]
  LL<-sum(log(dnbinom(slugs$slugs,mu=mu,size=theta)))
  -LL
}

#Estimate NB models
nlm(negNB.LL,c(2,1))->out.NB
nlm(negNB.LL1,c(2,1,1))->out.NB1

#Poisson neg LL functions from last time
negpois2.LL<-function(p){
  z<-as.numeric(slugs$field)-1
  mu<-p[1]+p[2]*z
  LL<-sum(log(dpois(slugs$slugs,lambda=mu)))
  -LL
}

negpois1.LL<-function(p){
  z<-as.numeric(slugs$field)-1
  mu<-p[1]
  LL<-sum(log(dpois(slugs$slugs,lambda=mu)))
  -LL
}

#Estimate Poisson models
nlm(negpois1.LL,2)->out.pois1
nlm(negpois2.LL,c(2,1))->out.pois2
out.pois1
out.pois2
names(slugs)

### fit Poisson regression models using glm function ###

#common mean model
glm(slugs~1,data=slugs,family=poisson(link=identity))->out0

#separate means model
glm(slugs~field,data=slugs,family=poisson(link=identity))->out1

#compare log-likelihoods with nlm
-out.pois1$minimum
logLik(out0)
-out.pois2$minimum
logLik(out1)

#carry out likelihood ratio test nlm output
2*(out.pois1$minimum-out.pois2$minimum)->LRstat
LRstat
1-pchisq(LRstat,1)

#LR test using anova function
anova(out0,out1,test='Chisq')

#get Hessian for Wald tests
nlm(negpois2.LL,c(2,1),hessian=T)->out.pois2

#compare standard errors
summary(out1)
sqrt(diag(solve(out.pois2$hessian)))
sqrt(diag(solve(out.pois2$hessian)))->std.err

#carry out Wald tests
out.pois2$estimate/std.err
2*(1-pnorm(abs(out.pois2$estimate/std.err)))

#compare identity link to default log link for mean
glm(slugs~field,data=slugs,family=poisson(link=identity))->out1
glm(slugs~field,data=slugs,family=poisson)->out2
coef(out1)
coef(out2)
exp(coef(out2))

#for such a simple model there is no difference in log-likelihood with different links
logLik(out1)
logLik(out2)

### fit negative binomial regression model using glm.nb function from MASS ###

out.NB
out.NB1
library(MASS)

#fit with both identity and log links
glm.nb(slugs~1,data=slugs)->out0a
glm.nb(slugs~1,data=slugs,link=identity)->out0b
glm.nb(slugs~field,data=slugs)->out1a
glm.nb(slugs~field,data=slugs,link=identity)->out1b
coef(out0b)
coef(out1b)

#theta is stored separately
out1b$theta
out.NB1$estimate

#likelihood ratio test
anova(out0b,out1b,test='Chisq')
2*(out.NB$minimum-out.NB1$minimum)
summary(out1b)
logLik(out1b)

#examine what is stored in model object
names(out1b)

#summary table can also be stored
summary(out1b)->out.summary
names(out.summary)
out.summary$coefficients

#standard errors
out.summary$coefficients[,2]
out1b$coefficients
coef(out1b)

##### negative binomial model with a continuous predictor #####

gala<-read.table('http://www.unc.edu/courses/2010fall/ecol/563/001/data/lectures/galapagos.txt', header=T)
dim(gala)
gala[1:10,]

#model with log link
glm.nb(Species~log(Area),data=gala)->out0
summary(out0)

#intercept-only model
glm.nb(Species~1,data=gala)->out0a

#likelihood ratio test for log(area) predictor
anova(out0a,out0,test='Chisq')

#try to fit model with identity link
glm.nb(Species~log(Area),data=gala,link=identity)->out0b
#need to supply an initial guess for theta
glm.nb(Species~log(Area),data=gala,link=identity,init.theta=1)->out0b

#model with log link yields larger log-likelihood
logLik(out0)
logLik(out0b)

#predictions on level of response
fitted(out0)
#predictions on level of link
predict(out0)
gala$Species

#predicted counts
gala$mu <- fitted(out0)
#predicted probability
gala$z <- dnbinom(gala$Species, mu=gala$mu, size=out0$theta)

#largest observed richness value
max(gala$Species)

library(lattice)

#plot NB distributions for each island separately. Include observed richness value
xyplot(z~Species|Island, data=gala, xlab='# of species', ylab='Probability',
       ylim=c(0,.03), xlim=c(-25, max(gala$Species)+25), panel=function(x,y,subscripts) {
         panel.xyplot(0:450,dnbinom(0:450,mu=gala$mu[subscripts],size=out0$theta), type='h', col='grey')
         panel.abline(v=gala$Species[subscripts],col=2,lty=2)
       })

#segregate islands by size
grp1<-gala[gala$Area<5,]
grp2<-gala[gala$Area>=5,]

#plot small islands only
xyplot(z~Species|Island, data=grp1, xlab='# of species', ylab='Probability',
       ylim=c(0,.15), xlim=c(-25, max(grp1$Species)+25), panel=function(x,y,subscripts) {
         panel.xyplot(0:450,dnbinom(0:450,mu=grp1$mu[subscripts],size=out0$theta), type='h', col='grey')
         panel.abline(v=grp1$Species[subscripts],col=2,lty=2)
       })

#plot big islands only
xyplot(z~Species|Island, data=grp2, xlab='# of species', ylab='Probability',
       ylim=c(0,.02), xlim=c(-25,525), panel=function(x,y,subscripts) {
         panel.xyplot(0:500,dnbinom(0:500,mu=grp2$mu[subscripts],size=out0$theta), type='h', col='grey')
         panel.abline(v=grp2$Species[subscripts],col=2,lty=2)
       })

#prob(observed species richness value given negative binomial model for this island)
1-pnbinom(gala$Species,mu=gala$mu,size=out0$theta)

#upper-tailed p-value
dotplot(gala$Island~1-pnbinom(gala$Species,mu=gala$mu,size=out0$theta),xlab='p-value')

#add vertical line at alpha=.05
dotplot(gala$Island~1-pnbinom(gala$Species,mu=gala$mu,size=out0$theta),xlab='upper p-value',panel=function(x,y) {
  panel.dotplot(x,y)
  panel.abline(v=.05,col=2,lty=2) })

#is it unusual to obtain two significant results?
dim(gala)
dim(gala)[1]*.05

#lower-tailed p-value plot
dotplot(gala$Island~pnbinom(gala$Species,mu=gala$mu,size=out0$theta),xlab='lower p-value',panel=function(x,y) {
  panel.dotplot(x,y)
  panel.abline(v=.05,col=2,lty=2) })