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
ms<-8

out<-function(i,o){
  return(egpd_s(m=i,x,model ="Pa",omega=o,badu=F,mfun=1)$gamma)
}

estimators<-data.frame(matrix(NA,ncol=ms+5,nrow=(n-1)))

result<-mclapply(c(1:ms),out,o=0.01,mc.cores=detectCores())

K<-1:(n-1)
M<-n/log(n)
estimators[,ms+1]<-egpd_s(m=1,x,model ="Pa",omega=0.01,badu=T,mfun=ceiling(K/log(K)))$gamma
estimators[,ms+2]<-egpd_s(m=1,x,model ="Pa",omega=0.01,badu=T,mfun=ceiling(M*(K/n)^1))$gamma
estimators[,ms+3]<-egpd_s(m=1,x,model ="Pa",omega=0.01,badu=T,mfun=ceiling(M*(K/n)^2))$gamma
estimators[,ms+4]<-egpd_s(m=1,x,model ="Pa",omega=0.01,badu=T,mfun=ceiling(M*(K/n)^3))$gamma
estimators[,ms+5]<-egpd_s(m=1,x,model ="Pa",omega=0.01,badu=T,mfun=ceiling(M*(K/n)^4))$gamma
estimators[,ms+6]<-RevEPD(x,-1)[1:(n-1)]
estimators[,ms+7]<-RevEPD(x,-2)[1:(n-1)]
estimators[,ms+8]<-RevEPD(x,-3)[1:(n-1)]
estimators[,ms+9]<-RevEPD(x,-4)[1:(n-1)]


for(i in 1:ms){
K<-i:(n-1)
#Sys.sleep(1)
estimators[K,i]<-result[[i]]
#lines(K,result[[i]],col=i+1,lwd=2)
print(i)
}

# #EPD(x,add=T,col=1,lwd=5)
# Hill(x,plot=T,ylim=c(0,1),main="Estimatates of the EVI")
# lines(estimators[,ms+9],col=2,lwd=3)
# lines(i:(n-1),result[[i]],col=i+1,lwd=2)
# abline(h=xi)
# lines(estimators[,ms+3],lwd=4,col="red4")
# # abline(h=xi)
# 
# plot(mfun1,type="l",lwd=4)
# lines(mfun2,col="red4",lwd=4)
# 
# a<- data.frame(estimators)
# names(a)<-c(paste("m_",c(1:ms),sep=""),paste("badu_",c(1:5),sep=""),paste("EPD_",c(1:4)))
# b<-melt(a,id="K")
# names(b)<-c("K","Estimator","Gamma")
# p<-ggplot(b,aes(x=K,y=Gamma,col=Estimator))+geom_line()
# library(plotly)
# ggplotly(p)
# 
# save(b,file="data.RData")


