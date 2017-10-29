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
ms<-8

out<-function(i,o){
  return(egpd_s(m=i,x,model ="Pa",omega=o,badu=F,mfun=1)$gamma)
}

estimators<-data.frame(matrix(NA,ncol=ms,nrow=(n-1)))
#EPD(x,add=T,direct=T,col=1,lwd=5)
Hill(x,plot=T,ylim=c(0,1))
lines(RevEPD(x,-0.5)[1:(n-1)],col=1,lwd=3)
abline(h=xi)

result<-mclapply(c(1:ms),out,o=0.01,mc.cores=detectCores())
for(i in 1:ms){
K<-i:(n-1)
Sys.sleep(1)
lines(K,(estimators[K,i]<-result[[i]]),col=i+1,lwd=2)
print(i)
}

K<-1:(n-1)
M<-n/log(n)
mfun1=ceiling(K/log(K))
mfun2=ceiling(M*(K/n)^3)
bads<-egpd_s(m=1,x,model ="Pa",omega=0.01,badu=T,mfun=mfun1)$gamma
lines(bads,lwd=3,col="steelblue")
# abline(h=xi)
