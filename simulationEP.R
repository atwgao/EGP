#simulation
##########
rm(list=ls())
########################
### Inference via ML ###
########################
library(gmm)
library(ReIns)
library(parallel)
library(methods)
source("./headers/Raux.R")
source("./headers/Distributions.R")
source("./headers/auxillary.R")
source("./headers/function_inference.R")
source("./headers/parametric_functions_Inference.R")

simul.r <- function(r,n,m,babu,sigma=NA,xi=NA,empiric=FALSE,seed=1){
  set.seed(seed+r)
  #x <- rgamma(n,1,1)
  #x <- rlnorm(n,0,1)
  #x<-runif(n)
  x <- rbeta(n,2,2)
  #x <- ReIns::rburr(n,alpha=1/xi,rho=-1)
  #x <- rweibull(n,10,0.5)
  #x <- ReIns::rfrechet(n,1/xi)
  #x<- abs(rt(n,1/xi))
  #x <- exp(rgamma(n,4,2))
  x <- x[x>0]
  bb<-m
  return(list(EGP1=egpd_s(m=m,x,model ="GP",omega=0,badu=babu)$gamma,EGP2=egpd_s(m=m,x,model ="GP",omega=100,badu=babu)$gamma,
              H=ReIns:::GPDmle(x)$gamma[bb:(n-1)],EPD=Moment(x)$gamma[bb:(n-1)]))
}

nsim<-40
m=1
babu <- TRUE
n<-100
xi<--0.5
t1<-system.time(res <- mclapply(c(1:nsim),simul.r,n=n,m=ifelse(babu==T,1,m),babu=babu,sigma=sigma,xi=xi,seed=1,mc.cores=detectCores()))

return.arg <- function(l,i){
  return(l[[i]])
}

return.arg2 <- function(l){
  if(!is.null(ncol(l[[1]]))) {return(array(unlist(l), dim = c(nrow(l[[1]]), ncol(l[[1]]), length(l))))
  }else{ return(matrix(unlist(l),ncol=length(l)))}
}

EGP1<-rowMeans(return.arg2(lapply(res,return.arg,i=1)),na.rm=TRUE)
EGP2<-rowMeans(return.arg2(lapply(res,return.arg,i=2)),na.rm=TRUE)
H   <-rowMeans(return.arg2(lapply(res,return.arg,i=3)),na.rm=TRUE)
EPD <-rowMeans(return.arg2(lapply(res,return.arg,i=4)),na.rm=TRUE)

rmse_EGP1<-rowMeans(sqrt((return.arg2(lapply(res,return.arg,i=1))-xi)^2),na.rm=TRUE)
rmse_EGP2<-rowMeans(sqrt((return.arg2(lapply(res,return.arg,i=2))-xi)^2),na.rm=TRUE)
rmse_H<-rowMeans(sqrt((return.arg2(lapply(res,return.arg,i=3))-xi)^2),na.rm=TRUE)
rmse_EPD<-rowMeans(sqrt((return.arg2(lapply(res,return.arg,i=4))-xi)^2),na.rm=TRUE)

var_EGP1<-apply(return.arg2(lapply(res,return.arg,i=1)),1,var,na.rm=TRUE)
var_EGP2<-apply(return.arg2(lapply(res,return.arg,i=2)),1,var,na.rm=TRUE)
var_H<-apply(return.arg2(lapply(res,return.arg,i=3)),1,var,na.rm=TRUE)
var_EPD<-apply(return.arg2(lapply(res,return.arg,i=4)),1,var,na.rm=TRUE)

K<-m:(n-1)
#pdf(file="EP_frechet_evi_m32.pdf",width = 6, height = 6)
#par(mar=c(5,5,5,2))
plot(K,EGP1,type="l",xlab="K",ylab=bquote(EVI),ylim=c(-1.5,xi+0.5),lwd=2.5,lty=1)
lines(K,EGP2,col=4,lwd=2.5,lty=1)
lines(K,H,col=8,lwd=2,lty=2)
lines(K,EPD,col="darkred",lwd=3,lty=1)
abline(h=xi)
#legend("bottomright",col=c(1,4,8,"darkred"),lwd=c(2.5,2.5,2,2),lty=c(1,1,2,1),
#       legend=c(expression(paste(H[xi]^Pareto,(w==0)),paste(H[xi]^Pareto,(w==0.01)),paste(H[k]),paste(EPD))),ncol=1,cex=1)
#dev.off()

#pdf(file="EP_frechet_rmse_m32.pdf",width = 6, height = 6)
#par(mar=c(5,5,5,2))
plot(K,rmse_EGP1,type="l",xlab="K",ylab=bquote(RMSE),ylim=c(-0.3,xi+0.3),lwd=2.5,lty=1)
lines(K,rmse_EGP2,col=4,lwd=2.5,lty=1)
lines(K,rmse_H,col=8,lwd=2,lty=2)
lines(K,rmse_EPD,col="darkred",lwd=2,lty=1)
abline(h=0)
#legend("bottom",col=c(1,4,8,"darkred"),lwd=c(2.5,2.5,2,2),lty=c(1,1,2,1),
#       legend=c(expression(paste(H[xi]^Pareto,(w==0)),paeste(H[xi]^Pareto,(w==0.01)),paste(H[k]),paste(EPD))),ncol=2,cex=1)
#dev.off()
# 
# pdf(file="EP_burr_var_m32.pdf",width = 6, height = 6)
# par(mar=c(5,5,5,2))
# plot(K,var_EGP1,type="l",xlab="K",ylab=bquote(Variance),ylim=c(-0.01,0.01),lwd=2.5,lty=1)
# lines(K,var_EGP2,col=4,lwd=2.5,lty=1)
# lines(K,var_H,col=8,lwd=2,lty=1)
# lines(K,var_EPD,col="darkred",lwd=2,lty=2)
# abline(h=0)
# legend("bottom",col=c(1,4,8,"darkred"),lwd=c(2.5,2.5,2,2),lty=c(1,1,2,1),
#        legend=c(expression(paste(H[xi]^Pareto,(w==0)),paste(H[xi]^Pareto,(w==0.01)),paste(H[k]),paste(EPD))),ncol=2,cex=1)
# dev.off()
# 
# save.image(file="../../results_rdata/burr_imagefileEP_m32.RData")