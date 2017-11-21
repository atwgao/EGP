#consolidate data
distr<-"frechet"
dat.files <- dir(paste("results/Pa/",distr,sep=''))
load(paste("results/Pa/",distr,"/",dat.files[1],sep=''))
nsimul<-length(dat.files)
results_Pa<-array(dim=c(.dims[1:2],nsimul))
results_EP<-array(dim=c(.dims,nsimul))
rm(ms,.betas,.dims,estimators_Pa,estimators_EP)
for(i in 1:nsimul){
  load(paste("results/Pa/",distr,"/",dat.files[i],sep=''))
  results_Pa[,,i]<-as.matrix(estimators_Pa)
  results_EP[,,,i]<-estimators_EP
  rm(ms,.betas,.dims,estimators_Pa,estimators_EP)
  print(i)
}

bias_Pa<-apply(results_Pa,c(1,2),mean,na.rm=T)
rmse_Pa<-(apply((results_Pa-0.5)^2,c(1,2),mean,na.rm=T))

bias_EP<-apply(results_EP,c(1,2,3),mean,na.rm=T)
rmse_EP<-(apply((results_EP-0.5)^2,c(1,2,3),mean,na.rm=T))

save(bias_Pa,rmse_Pa,bias_EP,rmse_EP,file=paste("./consolidated_data/EP_",distr,"_simulation500.RData",sep=""))
