library(gmm)
library(ReIns)
library(parallel)
library(methods)
source("./headers/Raux.R")
source("./headers/Distributions.R")
source("./headers/auxillary.R")
source("./headers/function_inference_pareto.R")
source("./headers/parametric_functions_Inference.R")

n<-200
xi<--0.75

x<- model3(n,0.5,-3/2)#model4(n,-0.25)#rbeta(n,2,2)#rlnorm(n,0,1)#rgamma(n,1,1)#rburr(n,1/xi,-1) #exp(rgamma(n,4,2))#rfrechet(n,1/xi)#
ms<-38
t1<-system.time({
  out<-function(i,o){
    return(egpd_s(m=i,x,model ="GP",omega=o,badu=F,mfun=1)$gamma)
  }
  
  estimators<-data.frame(matrix(NA,ncol=ms+6,nrow=(n-1)))
  
  result<-mclapply(c(1:ms),out,o=0,mc.cores=detectCores())
  
  K<-1:(n-1)
  M<-n/log(n)
  estimators[,ms+1]<-egpd_s(m=1,x,model ="GP",omega=0,badu=T,mfun=ceiling(K/log(K)))$gamma
  estimators[,ms+2]<-egpd_s(m=1,x,model ="GP",omega=0,badu=T,mfun=ceiling(M*(K/n)^1))$gamma
  estimators[,ms+3]<-egpd_s(m=1,x,model ="GP",omega=0,badu=T,mfun=ceiling(M*(K/n)^2))$gamma
  estimators[,ms+4]<-egpd_s(m=1,x,model ="GP",omega=0,badu=T,mfun=ceiling(M*(K/n)^3))$gamma
  estimators[,ms+5]<-egpd_s(m=1,x,model ="GP",omega=0,badu=T,mfun=ceiling(M*(K/n)^4))$gamma
  estimators[,ms+6]<-Moment(x)$gamma[1:(n-1)]

  
  for(i in 1:ms){
    K<-i:(n-1)
    #Sys.sleep(1)
    estimators[K,i]<-result[[i]]
    #lines(K,result[[i]],col=i+1,lwd=2)
    print(i)
  }
  
  # #EPD(x,add=T,col=1,lwd=5)
   GPDmle(x,plot=T,ylim=c(-1,2),main="Estimatates of the EVI")
   #  lines(estimators[,1],col=2,lwd=3)
   # lines(i:(n-1),result[[i]],col=i+1,lwd=2)
   # abline(h=xi)
   lines(estimators[,ms+3],lwd=4,col="red4")
  abline(h=xi)
  # 
  # plot(mfun1,type="l",lwd=4)
  # lines(mfun2,col="red4",lwd=4)
  # 
  df<- data.frame(estimators)
  names(df)<-c(paste("m_",c(1:ms),sep=""),paste("badu_",c(1:5),sep=""),paste("Moment"))
  # b<-melt(a,id="K")
  # names(b)<-c("K","Estimator","Gamma")
  # p<-ggplot(b,aes(x=K,y=Gamma,col=Estimator))+geom_line()
  # library(plotly)
  # ggplotly(p)
  # 
})
#save(df,file="./results/frechet_data_m.RData")
print(t1)


#nothing happening here
