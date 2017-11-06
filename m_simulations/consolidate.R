#consolidate data

distr<-"burr"
results<-array(dim=c(199,44,500))

for(i in 1:500){
load(paste(distr,"/","GP_",distr,"_data_m_",i,".RData",sep=""))
  results[,,i]<-as.matrix(df)
  rm(df)
  print(i)
}

bias<-apply(results,c(1,2),mean,na.rm=T)
rmse<-(apply((results-0.5),c(1,2),mean,na.rm=T)^2)

save(bias,rmse,file=paste("GP_",distr,"_simulation500.RData",sep=""))
