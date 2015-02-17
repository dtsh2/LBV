
rm(list=ls())
## pomp test run lbv
getwd()
setwd("~/GitHub/LBV") # revise as necessary
tiff("hist_sero.tiff",width=8,height=8,units='in',res=300, compression = "lzw")
res<-c(31.7,
60,
33,
43,
80,
26.7,
40.4,
43.8,
5.8,
32.7)
hist(res,breaks=6,col="grey",xlab="Seroprevalence (%)",
     main="")
summary(res)
dev.off()

