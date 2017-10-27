rm(list=ls())
########################
### Inference via ML ###
########################
library(gmm)
library(ReIns)
library(parallel)
library(methods)
source("Raux.R")
source("Distributions.R");source("auxillary.R");source("function_inference.R")

data<-read.csv("CYA_newdata.txt")$Ultimate
x<-tail(data,200)

#x<-ReIns::rburr(100,alpha=2,rho=-1)
#x<-rbeta(100,2,4)

#ms<-seq(5,20, by = 5)
#system.time(res<-mclapply(ms,egpd_s,data=x,omega=10,mc.cores=detectCores()))
system.time(EGP_penm5<-egpd_s(m=25,x,model ="Pa",omega=10,plot=T,lwd=2,ylim=c(-1,2)))
ReIns::Hill(x,add=T,col=2)
ReIns::EPD(x,add=T,col=4,direct=T)
ReIns::GPDmle(x,add=T,col=14)
abline(h=0.5)
# plot(res[[1]]$k,res[[1]]$gamma,type="l",ylim=c(0,2),lwd=2)
# I<-1:length(ms)
# for(i in I) lines(res[[i]]$k,res[[i]]$gamma,col=i,lwd=2)
# legend("topright",legend=paste("m =",ms),col=I,lwd=c(2,2,2,2,2))

pdf(file="CY2_MeanExcess.pdf",width = 5, height = 5)
par(mar=c(6,5,5,2))
ReIns::MeanExcess(data)
dev.off()
