text <- readLines("../different_ms_GP.R")
distr<-"model3"
for(i in 1:500){
  dir.create(paste(distr,sprintf("/case%03d",i,sep=""),sep=""))
  text[65]<-paste("save(df,file=\"../../results/GP/",distr,"/GP_",distr,"_data_m_",i,".RData\")",sep="")
  text[5:9]<-c("source(\"../../../headers/Raux.R\")",                        
               "source(\"../../../headers/Distributions.R\")",                 
               "source(\"../../../headers/auxillary.R\")",                     
               "source(\"../../../headers/function_inference_pareto.R\")",     
               "source(\"../../../headers/parametric_functions_Inference.R\")")
  
  file.create(paste("./",distr,"/",sprintf("case%03d",i),"/different_ms_GP.R",sep=""))
  write(file=paste("./",distr,"/",sprintf("case%03d",i),"/different_ms_GP.R",sep=""),text)
print(i)
  }