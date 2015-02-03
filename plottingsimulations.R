times<-c(1:100000)
norm<-rnorm(100000)
sim<-seq(from=1, to = 1000, by = 100)
hist(norm)
data<-as.data.frame(cbind(times,norm,sim))
#plot blank
with(data,plot(times,norm,type="n"))
#plot lines
for(i in unique(data$sim))
  with(data[ data$sim==i,],lines(times,norm,col=i))

library(ggplot2)
#plot
ggplot(data=data,aes(times,norm,colour=as.factor(sim))) +
  geom_line()
