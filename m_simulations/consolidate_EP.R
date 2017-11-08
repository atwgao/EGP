#consolidate data
nsimul<-500
distr<-"burr"
results_Pa<-array(dim=c(.dims[1:2],nsimul))
results_EP<-array(dim=c(.dims,nsimul))

for(i in 1:500){
  load(paste(distr,"/","GP_",distr,"_data_m_",i,".RData",sep=""))
  results_Pa[,,i]<-as.matrix(estimators_Pa)
  results_EP[,,,i]<-estimators_EP
  rm(ms,.betas,.dims,estimators_Pa,estimators_EP)
  print(i)
}

bias_Pa<-apply(results_Pa,c(1,2),mean,na.rm=T)
rmse_Pa<-(apply((results_Pa-0.5),c(1,2),mean,na.rm=T)^2)

bias_EP<-apply(results_EP,c(1,2,3),mean,na.rm=T)
rmse_EP<-(apply((results_EP-xi),c(1,2,3),mean,na.rm=T)^2)

save(bias,rmse,file=paste("EP_",distr,"_simulation500.RData",sep=""))
