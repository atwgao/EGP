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
x <- rburr(n,1/xi,-1)


plot(RevEPD(x),col=1,lwd=3,type="l",ylim=c(0,2))
Hill(x,add=T,lwd="orange")
#EPD(x,add=T,direct=T,col=1,lwd=5)
abline(h=xi)
for(i in 1:1){
K<-i:(n-1)
lines(K,egpd_s(m=i,x,model ="Pa",omega=0,badu=F)$gamma,col=i+1,lwd=2)
}
