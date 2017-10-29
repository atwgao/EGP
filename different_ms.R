library(gmm)
library(ReIns)
library(parallel)
library(methods)
source("./headers/Raux.R")
source("./headers/Distributions.R")
source("./headers/auxillary.R")
source("./headers/function_inference_pareto.R")
source("./headers/parametric_functions_Inference.R")

n<-100
xi<-0.5
x <- rburr(n,1/xi,-0.5)#abs(rt(n,2))#rfrechet(n,1/xi) 
ms<-3

out<-function(i,o){
  return(egpd_s(m=i,x,model ="Pa",omega=o,badu=F)$gamma)
}

estimators<-data.frame(matrix(NA,ncol=1+ms,nrow=(n-1)))
#plot(estimators[,1]<-RevEPD(x)[1:(n-1)],col=1,lwd=3,type="l",ylim=c(0,2),ylab="EVI",main="Estimates of the Extreme value index",xlab="K")

#EPD(x,add=T,direct=T,col=1,lwd=5)
Hill(x,plot=T,ylim=c(0,1))
lines(RevEPD(x)[1:(n-1)],col=1,lwd=3)
abline(h=xi)
result<-mclapply(c(1:ms),out,o=0.01,mc.cores=64)
for(i in 1:ms){
K<-i:(n-1)
Sys.sleep(1)
lines(K,(estimators[K,i]<-result[[i]]),col=i+1,lwd=2)
print(i)
}


# plot(estimators[,2],col=1,lwd=3,type="l",ylim=c(0,2),ylab="EVI",main="Estimates of the Extreme value index",xlab="K")

bads<-egpd_s(m=1,x,model ="Pa",omega=0.01,badu=T)$gamma
lines(bads,lwd=3,col="steelblue")
# abline(h=xi)
