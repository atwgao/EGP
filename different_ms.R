library(gmm)
library(ReIns)
library(parallel)
library(methods)
source("./headers/Raux.R")
source("./headers/Distributions.R")
source("./headers/auxillary.R")
source("./headers/function_inference_pareto.R")
source("./headers/parametric_functions_Inference.R")

n<-500
xi<-0.5
x <- rburr(n,1/xi,-1)
ms<-16

out<-function(i){
  return(egpd_s(m=i,x,model ="Pa",omega=0,badu=F)$gamma)
}

estimators<-data.frame(matrix(NA,ncol=ms,nrow=(n-1)))
plot(estimators[,1]<-RevEPD(x)[1:(n-1)],col=1,lwd=3,type="l",ylim=c(0,2),ylab="EVI",main="Estimates of the Extreme value index",xlab="K")
#EPD(x,add=T,direct=T,col=1,lwd=5)
abline(h=xi)
result<-mclapply(c(1:ms),out,mc.cores=detectCores())
for(i in 2:ms){
K<-i:(n-1)
lines(K,(estimators[K,i]<-result[[i]]),col=i+1,lwd=2)
print(i)
}


#plot(estimators[,1],col=1,lwd=3,type="l",ylim=c(0,2),ylab="EVI",main="Estimates of the Extreme value index",xlab="K")
