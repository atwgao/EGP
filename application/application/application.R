#application
library(ReIns)
library(gmm)
library(ReIns)
library(parallel)
library(methods)
library(dplyr)
source("../../headers/Raux.R")
source("../../headers/Distributions.R")
source("../../headers/auxillary.R")
source("../../headers/function_inference.R")
source("../../headers/parametric_functions_Inference.R")
source("QQplots.R")

alburq <-read.csv("Albuquerque.txt",header=F)$V1
desmoin<-read.csv("DesMoines.txt",header=F)$V1
grand  <-read.csv("GrandRapids.txt",header=F)$V1
zaven  <-read.csv("zaventem.txt",header=F)$V1
massmax<-read.csv("maasmax.txt",header=F,sep=" ")$V1
ecad   <-read.csv("ecad00045TX.txt",header=F,sep=" ")
secura <-as.data.frame(read.csv("secura1.txt",header=F,sep=""))
secura <-secura %>%  mutate_all( funs(as.character(.)), names( .[,sapply(., is.factor)] )) %>% select(V2) %>% .$V2
secura <- as.numeric(gsub(",","",secura))
soa <- read.csv("soa.txt",header=F,sep="")$V1
aon <- read.csv("aon.txt",header=F,sep="")$V2
sp500<-read.csv("sp500.txt",header=F,sep=",")$V2
cpt <-read.csv("cape_town_wind.txt",header=F,sep="")$V1
load("rain30339001_Seasons.RData")
rain <- read.csv("rain.csv",header=T)$pr
temp <- read.csv("temparature.csv",header=T)$tas

x<-alburq
x<-x[x>=quantile(x)[4]]

x<-x[x>30]
n<-length(x);print(n)
ReIns::ExpQQ(x)
ReIns::ParetoQQ(x)
WeibullQQ(x)
MeanExcess(x)
Moment(x,plot=T,ylim=c(-2,2))
GPDmle(x,add=T,col=2,lwd=2)
K<-1:(n-1)
lines(100:(length(x)-1),(res<-egpd_s(x,model="Pa",omega=0,m=100,badu=F,mfun=K^0.5,beta=-2)$gamma),col="red4",lwd=2)
theta_GP = suppressWarnings(EGP.fitPWM(x=x,type=1,kappa0=2,sigma0=3,xi0=-0.5))
abline(h=theta_GP[3])