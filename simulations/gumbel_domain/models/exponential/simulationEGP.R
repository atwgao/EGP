#simulation

rm(list=ls())
########################
### Inference via ML ###
########################
library(gmm)
library(ReIns)
library(parallel)
library(methods)
source("../../../../headers/Raux.R")
source("../../../../headers/Distributions.R")
source("../../../../headers/auxillary.R")
source("../../../../headers/function_inference.R")
source("../../../../headers/parametric_functions_Inference.R")

simul.r <- function(r,n,m,sigma=NA,xi=NA,empiric=FALSE,seed=1){
  set.seed(seed+r)
  x <- rexp(n)
  x <- x[x>0]
  return(list(EGP1=egpd_s(m=m,x,model ="GP",omega=0,badu=babu)$gamma,EGP2=egpd_s(m=m,x,model ="GP",omega=100,badu=babu)$gamma,
              GP=ReIns:::GPDmle(x)$gamma,M=Moment(x)$gamma))
}

nsim<-500
m=1
babu <- TRUE
n<-200
xi<-0
sigma<-1
t1<-system.time(res <- mclapply(c(1:nsim),simul.r,n=n,m=m,sigma=sigma,xi=xi,seed=1,mc.cores=detectCores()))

return.arg <- function(l,i){
  return(l[[i]])
}

return.arg2 <- function(l){
  if(!is.null(ncol(l[[1]]))) {return(array(unlist(l), dim = c(nrow(l[[1]]), ncol(l[[1]]), length(l))))
  }else{ return(matrix(unlist(l),ncol=length(l)))}
}

EGP1<-rowMeans(return.arg2(lapply(res,return.arg,i=1)))
EGP2<-rowMeans(return.arg2(lapply(res,return.arg,i=2)))
GP<-rowMeans(return.arg2(lapply(res,return.arg,i=3)))
Mom<-rowMeans(return.arg2(lapply(res,return.arg,i=4)))


rmse_EGP1<-rowMeans(sqrt((return.arg2(lapply(res,return.arg,i=1))-xi)^2))
rmse_EGP2<-rowMeans(sqrt((return.arg2(lapply(res,return.arg,i=2))-xi)^2))
rmse_GP<-rowMeans(sqrt((return.arg2(lapply(res,return.arg,i=3))-xi)^2))
rmse_Mom<-rowMeans(sqrt((return.arg2(lapply(res,return.arg,i=4))-xi)^2))


var_EGP1<-apply(return.arg2(lapply(res,return.arg,i=1)),1,var)
var_EGP2<-apply(return.arg2(lapply(res,return.arg,i=2)),1,var)
var_GP<-apply(return.arg2(lapply(res,return.arg,i=3)),1,var)
var_Mom<-apply(return.arg2(lapply(res,return.arg,i=4)),1,var)

K<-m:(n-1)
pdf(file="../../results_plots/EGP_exponential_xi.pdf",width = 6, height = 6)
par(mar=c(3,5,1,5))
plot(K,EGP1,type="l",xlab="K",ylab=bquote(EVI),ylim=c(xi-1,xi+0.5),lwd=2.5,lty=1)
lines(K,EGP2,col=4,lwd=2.5,lty=1)
lines(GP,col="darkgrey",lwd=2,lty=2)
lines(K,Mom,col="red4",lwd=3,lty=1)
abline(h=xi)
legend("topleft",col=c(1,4,"darkgrey","red4"),lwd=c(2.5,2.5,2,3),lty=c(1,1,2,1),
       legend=c(expression(paste(H[xi]^GP,(w==0)),paste(H[xi]^GP,(w==1)),paste(GPD),paste(Moment))),ncol=2,cex=1)
dev.off()


pdf(file="../../results_plots/EGP_exponential_rmse.pdf",width = 6, height = 6)
par(mar=c(5,5,5,2))
plot(K,rmse_EGP1,type="l",xlab="K",ylab=bquote(RMSE),ylim=c(-0.3,xi+0.5),lwd=2.5,lty=1)
lines(K,rmse_EGP2,col=4,lwd=2.5,lty=1)
lines(K,rmse_GP,col=8,lwd=2,lty=2)
lines(K,rmse_Mom,col="red4",lwd=3,lty=1)
abline(h=0)
legend("topleft",col=c(1,4,"darkgrey","red4"),lwd=c(2.5,2.5,2,3),lty=c(1,1,2,1),
       legend=c(expression(paste(H[xi]^GP,(w==0)),paste(H[xi]^GP,(w==1)),paste(GPD),paste(Moment))),ncol=2,cex=1)
dev.off()
 

pdf(file="../../results_plots/EGP_exponential_var.pdf",width = 6, height = 6)
par(mar=c(5,5,5,2))
plot(K,var_EGP1,type="l",xlab="K",ylab=bquote(Variance),ylim=c(-0.05,0.1),lwd=2.5,lty=1)
lines(K,var_EGP2,col=4,lwd=2.5,lty=1)
lines(K,var_GP,col=8,lwd=2,lty=2)
lines(K,var_Mom,col="red4",lwd=3,lty=1)
abline(h=0)
legend("topleft",col=c(1,4,"darkgrey","red4"),lwd=c(2.5,2.5,2,3),lty=c(1,1,2,1),
     legend=c(expression(paste(H[xi]^GP,(w==0)),paste(H[xi]^GP,(w==1)),paste(GPD),paste(Moment))),ncol=2,cex=1)
dev.off()

save.image(file="../../results_rdata/exponential_imagefileEG.RData")
