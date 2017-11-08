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

x <- rburr(n,1/xi,-0.5)
ms<-ceiling(n/log(n))

.cores<-2#detectCores()

K<-1:(n-1)
M<-n/log(n)
.mfun <- matrix(c(ceiling(M*(K/n)^0.25),ceiling(M*(K/n)^0.5),ceiling(M*(K/n)^0.75),ceiling(K/log(K)),
                  ceiling(M*(K/n)^1),ceiling(M*(K/n)^2),ceiling(M*(K/n)^3)), ncol=7,nrow=(n-1))

out.1<-function(i,o,model,beta){
  return(egpd_s(m=i,x,model = model,omega=o,badu=F,mfun=1,beta)$gamma)
}

out.2<-function(i,o,model,beta,.mfun){
  return(egpd_s(m=1,x,model = model,omega=o,badu=T,mfun=.mfun[,i],beta)$gamma)
}

.dims <- c((n-1),(ms+8),6)
estimators_Pa<-data.frame(array(dim=.dims[1:2]))
estimators_EP<-array(dim=.dims)

t1<-system.time({
result_Pa1<-mclapply(c(1:ms),out.1,o=0.01,model="Pa",mc.cores=.cores)
result_Pa2<-mclapply(c(1:7),out.2,o=0.01,model="Pa",.mfun=.mfun,mc.cores=.cores)
for(i in 1:ms){K<-i:(n-1);estimators_Pa[K,i]<-result_Pa1[[i]]}
for(i in (ms+1):(ms+7)){K<-1:(n-1);estimators_Pa[K,i]<-result_Pa2[[i-ms]]}
estimators_Pa[,ms+8]<-Hill(x)$gamma

.betas<-c(-0.5,-1,-2,-4,-8,-16)

for(j in 1:6){
  result_EP1<-mclapply(c(1:ms),out.1,o=0.01,model="EP",beta=.betas[j],mc.cores=.cores)
  result_EP2<-mclapply(c(1:7),out.2,o=0.01,model="EP",beta=.betas[j],.mfun=.mfun,mc.cores=.cores)
  for(i in 1:ms){K<-i:(n-1);estimators_EP[K,i,j]<-result_EP1[[i]]}
  for(i in (ms+1):(ms+7)){K<-1:(n-1);estimators_EP[K,i,j]<-result_EP2[[i-ms]]}
  estimators_EP[K,(ms+8),j]<-RevEPD(x,.betas[j])[1:(n-1)]
  print(j)
}
})

save(ms,xi,.betas,.dims,estimators_Pa,estimators_EP,file="data.RData")
print(t1)

