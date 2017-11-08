text <- readLines("../different_ms.R")
distr<-"burr"
for(i in 1:500){
  dir.create(paste(distr,sprintf("/case%03d",i,sep=""),sep=""))
  text[58]<-paste("save(ms,xi,.betas,.dims,estimators_Pa,estimators_EP,file=\"../../../results/Pa/",distr,"/EP_",distr,"_data_m_",i,".RData\")",sep="")
  text[5:9]<-c("source(\"../../../../headers/Raux.R\")",                        
               "source(\"../../../../headers/Distributions.R\")",                 
               "source(\"../../../../headers/auxillary.R\")",                     
               "source(\"../../../../headers/function_inference_pareto.R\")",     
               "source(\"../../../../headers/parametric_functions_Inference.R\")")
  
  file.create(paste("./Pa/",distr,"/",sprintf("case%03d",i),"/different_ms.R",sep=""))
  write(file=paste("./Pa/",distr,"/",sprintf("case%03d",i),"/different_ms.R",sep=""),text)
  print(i)
}