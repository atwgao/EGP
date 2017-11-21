text <- readLines("../different_ms_Pa.R")
distr<-"burr15"
simFun <- switch(distr,
  loggamma      = paste("xi = 0.5;  x <- exp(rgamma(n,4,2))"),
  burr          = paste("xi = 0.5;  x <- rburr(n,2,-1)"),
  burr05        = paste("xi = 0.5;  x <- rburr(n,2,-0.5)"),
  burr15        = paste("xi = 0.5;  x <- rburr(n,2,-1.5)"),
  frechet       = paste("xi = 0.5;  x <- rfrechet(n,2)"),
  stop("unkown distribution")
  )
text[13]<-simFun
  text[5:9]<-c("source(\"../../../../headers/Raux.R\")",                        
               "source(\"../../../../headers/Distributions.R\")",                 
               "source(\"../../../../headers/auxillary.R\")",                     
               "source(\"../../../../headers/function_inference_pareto.R\")",     
               "source(\"../../../../headers/parametric_functions_Inference.R\")")
for(i in 1:250){
  dir.create(paste("Pa/",distr,sprintf("/case%03d",i,sep=""),sep=""))
  text[55]<-paste("save(ms,xi,.betas,.dims,estimators_Pa,estimators_EP,file=\"../../../results/Pa/",distr,"/EP_",distr,"_data_m_",i,".RData\")",sep="")
  file.create(paste("./Pa/",distr,"/",sprintf("case%03d",i),"/different_ms.R",sep=""))
  write(file=paste("./Pa/",distr,"/",sprintf("case%03d",i),"/different_ms.R",sep=""),text)
  print(i)
}
