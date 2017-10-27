dEGP.BB_squared = function(xg,param,model){
  return(dEGP.BB(xg,param,model)^2)
}

par_optm<-function(x,mGrid,model){
  kappa0 = 3; sigma0 = 2; xi0 = 0.15;
  n<-length(x)
  xi0<-mean(ReIns::Hill(x)$gamma[4:(0.1*n)])
  if(model=="Pa") startv<-xi0 else startv<-c(2,xi0) 
  shrink.coef<-1
  theta0 <- startv
  xx<-sort(x)
  if(model=="Pa") x<-(xx/min(xx))[-1]
  n = length(x);
  t = proc.time()
  LSCV = numeric(length(mGrid)); results=list();
  bounds = list(inits=theta0,lower_gam=1,upper_gam=10)
  for (mi in 1:length(mGrid)){
    cat(paste("m=",mGrid[mi],"/",sep=""))
    results[[mi]] = EGPBBnon.fitMLE.dq(x,mGrid[mi],xgridMax = max(x),model,shrink.coef=0,bounds=bounds)
    LSCV[mi] = integrate(dEGP.BB_squared,lower = 0, upper = Inf,param=results[[mi]]$paramBnon,model, subdivisions = 500)$value
    f <- dEGP.BB(x[1],EGPBBnon.fitMLE.dq(x=x[-1], m=mGrid[mi],xgridMax = max(x),model,shrink.coef,bounds=bounds)$paramBnon,model)
    for(i in 2:n){
      print(i)
      f <- f + dEGP.BB(x[i],EGPBBnon.fitMLE.dq(x=x[-i], m=mGrid[mi],xgridMax = max(x),model,shrink.coef,bounds=bounds)$paramBnon,model)
    }
    LSCV[mi] = LSCV[mi] - 2*f/n
  }
  minLSCV<-which.min(LSCV)
  return(list(m=mGrid[minLSCV],results=results))
}

#Inference functions
neg.log.lik.G <- function(x,theta,type,model){
  if(theta[1]<0 | theta[1]>1 | theta[2]<=0 | theta[3]<theta[2] | theta[3]<=0 | theta[4]<=0 |theta[5]<=10^(-6)){
    obj<-10^7
  }else{
    obj<--sum(dExt(x,prob=theta[1],kappa=theta[2], delta=theta[3],sigma=theta[4],xi=theta[5],log=TRUE,model=model,type=type),na.rm=TRUE) 
  }
  return(obj)
}

neg.log.lik.G.delta <- function(x,theta,type,model){
  if(theta[1]<0 | theta[1]>1 | theta[2]<=0 | theta[3]<theta[2] | theta[3]<=0 | theta[4]<=0 |theta[5]<=10^(-6) ){
    obj<-10^7
  }else{
    obj<--sum(dExt.delta(x,prob=theta[1],kappa=theta[2], delta=theta[3],sigma=theta[4],xi=theta[5],delta.k=theta[6],log=TRUE,model=model,type=type),na.rm=TRUE) 
  }
  return(obj)
}

fitW_non = function(paramGPD,m,xg,kap,model){
  n = length(xg)
  w = numeric(m)
  if(model!="GP" && model!="Pa") stop("Extensions limited to Pa and GP")
  if(model=="GP"){
    H = evd:::pgpd(q=xg, scale=paramGPD[1], shape=paramGPD[2])
  }else{
    H = ReIns::ppareto(xg,shape=paramGPD,scale=1)#ppareto(xg, gamma=paramGPD[1])
  }
  for (k in 1:m){
    Fn = sum(H<=(k/m))/n;
    Fn_1 = sum(H<=((k-1)/m))/n;
    w[k] = Fn-Fn_1;
    if ((k==m) && w[k]==0){
      w[k] = max(Fn-Fn_1,kap/m);
      w = w/sum(w)
    }
  }
  return(w)
}

func_MLEnon = function(theta,xg,m,kap,model,shrink.coef){
  if(model!="GP" && model!="Pa") stop("Extensions limited to Pa and GP")
  if(model=="GP"){
    if(theta[1]<=10^-6|theta[2] < (0)|theta[2] > (0.99)) return(10^7)
    H = evd:::pgpd(q=xg, scale=theta[1], shape=theta[2])
    h = evd:::dgpd(x=xg, scale=theta[1], shape=theta[2])
  }else{
    if(theta[1]<=10^-6) return(10^7)
    H = ReIns::ppareto(xg,shape=theta[1],scale=1)#ppareto(xg, gamma=theta[1])
    h = ReIns::dpareto(xg,shape=theta[1],scale=1)#dpareto(xg, gamma=theta[1])
  }
  w.hat <- fitW_non(theta,m,xg,kap,model)
  r = -sum(log(dBB(w=w.hat,u=H)) + log(h)) + shrink.coef*(1-m*w.hat[m])^2 
  #r <- ifelse(is.finite(r),r,10^6)
  return(r)
}

func_repar_non = function(theta){
  xi = theta[2]
  sigma = theta[1];
  
  #### xi
  xi = xi;
  
  #### sigma
  sigma = sigma;
  
  return(param=c(sigma,xi))
}


func_repar_non2 = function(theta){
  xi = theta[1]
  #### xi
  xi = xi;
  
  return(param=xi)
}

EGPBBnon.fitMLE.dq = function(x,m,xgridMax=max(x),model,kap0,shrink.coef,params=NA,bounds) {
  n = length(x);
  xgrid = seq(0.001,xgridMax, length.out=1000)
  p = c(1:length(x))/(length(x)+1)
  if(model!="GP" && model!="Pa") stop("Extensions limited to Pa and GP")

    theta0 =  bounds$inits#mean(Hill(x)$gamma[1:(0.25*n)])
  # init algorithm
  if(model=="GP"){
    par.GPD = numeric(2);
    #res <- tryCatch({optim(par=theta0, fn=function(theta,x,m,kap,model,shrink.coef) func_MLEnon(func_repar_non(theta),x,m,kap,model,shrink.coef), x=x, m=m, model=model, kap=fit.parEGP.PWM[1],shrink.coef=shrink.coef)},error=function(e) return(list(par=theta0,value=1)))
    res <- tryCatch({optim(par=theta0, fn=function(theta,x,m,kap,model,shrink.coef) func_MLEnon(func_repar_non(theta),x,m,kap,model,shrink.coef), 
                           x=x, m=m, model=model, kap=kap0,shrink.coef=shrink.coef,
                           method="Nelder-Mead")},error=function(e) return(list(par=theta0,value=1)))
    go = TRUE; tol=1e-6;
    if(length(res)==2) go = FALSE
    t1<-proc.time()
    while (go){
      if(matrix(proc.time()-t1)[3,]>2) go=FALSE
      cat(".")
      #res1 <- tryCatch({optim(par=res$par, fn=function(theta,x,m,kap,model,shrink.coef) func_MLEnon(func_repar_non(theta),x,m,kap,model,shrink.coef), x=x, m=m, model=model,kap=fit.parEGP.PWM[1],shrink.coef=shrink.coef)},error=function(e) return(list(par=theta0,value=1)))
      res1 <- tryCatch({optim(par=res$par, fn=function(theta,x,m,kap,model,shrink.coef) func_MLEnon(func_repar_non(theta),x,m,kap,model,shrink.coef), 
                              x=x, m=m, model=model, kap=kap0,shrink.coef=shrink.coef,
                              method="Nelder-Mead")},error=function(e) return(list(par=theta0,value=1)))
      try(if(abs(res1$value-res$value)<tol)  {
        go=FALSE
      })
      res = res1
    } 
    par = res$par;
    par.GPD = func_repar_non(par);
  }else{
    par.GPD = numeric(1);
    res <- tryCatch({optim(par=theta0, fn=function(theta,x,m,kap,model,shrink.coef) func_MLEnon(func_repar_non2(theta),x,m,kap,model,shrink.coef), x=x, m=m, shrink.coef=shrink.coef, model=model, kap=kap0,method="Brent",lower=bounds$lower_gam,upper=bounds$upper_gam)},error=function(e) return(list(par=theta0,value=1)))
    go = TRUE; tol=1e-6;
    if(length(res)==2) go = FALSE
    t1<-proc.time()
    while (go){
      if(matrix(proc.time()-t1)[3,]>2) go=FALSE
      cat(".")
      res1 <- tryCatch({optim(par=res$par, fn=function(theta,x,m,kap,model,shrink.coef) func_MLEnon(func_repar_non2(theta),x,m,kap,model,shrink.coef), x=x, m=m, shrink.coef=shrink.coef, model=model, kap=kap0,method="Brent",lower=bounds$lower_gam,upper=bounds$upper_gam)},error=function(e) return(list(par=theta0,value=1)))
      try(if(abs(res1$value-res$value)<tol)  {
        go=FALSE
      })
      res = res1
    }
    par<-res$par
    par.GPD = func_repar_non2(par)
  }
  par.W = fitW_non(par.GPD,m,x,kap=kap0,model)
  fit.nonparEGP = c(par.W,par.GPD)
  if(model=="Pa") fit.nonparEGP = c(par.W,1/par.GPD)
  #q.nonparEGP = qEGP.BB(p,fit.nonparEGP,model);
  
  return(list(paramBnon = fit.nonparEGP))
}

# EGP.fitMLE.dq <- function(x,type=1,model=NA){
#   p<-c(1:length(x))/(length(x)+1)
#   try1 <- try( opt <- optim(par=theta0,fn=neg.log.lik.G,x=x,model=model,type=type,method="Nelder-Mead"))
#   if(!is(try1,"try-error")){
#     theta.hat <- opt$par
#     theta_quant<-qExt(p,prob=theta.hat[1],kappa=theta.hat[2],delta=theta.hat[3],sigma=theta.hat[4],xi=theta.hat[5],type=type,model=model)
#   } else{
#     theta.hat <- c(NA,NA,NA,NA,NA)
#   }
#   if(model=="Pa")  theta.hat<-c(theta.hat[-length(theta.hat)],1/theta.hat[length(theta.hat)])
#   return(list(paramB=theta.hat,q.parEGP=theta_quant))
# }

egpd.fit <- function(data,omega,model,badu,m){
  x<-data
  k1=n#floor(n^(0.995)) # high k value if not chosen adaptively
  n<-length(x)
  p<-c(1:length(x))/(length(x)+1)
  M<-ifelse(badu==T,ceiling(n/log(n)),m)
  bb<-ifelse(badu==T,1,M)
  K<-bb:(n-1)
  xx<-sort(x)
  colsize<-M+2
  theta0_Pa =  1/mean(Hill(x)$gamma[1:(0.4*n)])
  theta0_GP = tryCatch({suppressWarnings(EGP.fitPWM(x=x,type=1,kappa0=2,sigma0=3,xi0=0.5))},error = function(e) {return(c(1,2,1/theta0_Pa))})
  if(model=="Pa") {bounds<-list(inits=theta0_Pa,lower_gam=1,upper_gam=5)
  }else bounds<-list(inits=theta0_GP[2:3],lower_gam=-1.5,upper_gam=1)
  FixWeights<-FALSE #estimate 
  if(model=="Pa") colsize <- M+1
  fitnonparEGP<-fitnonparEGP_pen<-matrix(NA,ncol=colsize,nrow=(n-1))
  for(k in 1:(n-1)){mk<-ceiling(M*k/n);mk=ifelse(badu==T,ifelse(mk==Inf,1,mk),M)#floor(k/log(k))
  if(model=="GP"){
    y1<-xx[(n-k+1):n]-xx[n-k]
  }else{
    y1<-xx[(n-k+1):n]/xx[n-k]
  }
  shrink.coef<-0.5*omega*(k/n)^(-2)
  fitnonparEGP_pen[k,c((1:mk),((M+1):colsize))]<-EGPBBnon.fitMLE.dq(y1,mk,max(y1),model,kap0=theta0_GP[1],shrink.coef=shrink.coef,bounds=bounds)$par
  print(k);
  if(model=="Pa") {
    bounds<-list(inits=fitnonparEGP_pen[k,(M+1)],lower_gam=1, upper_gam = 5)#          
    }
  }
  return(list(gamma=fitnonparEGP_pen[K,],K=K,cs=colsize,k1=k1))
}

egpd_s<-function(data,omega=0,model ="GP",badu,m,warnings = FALSE, plot = FALSE, add = FALSE, main = "EGPD-BB estimates of EVI", ...){
  gammas<-egpd.fit(data, omega,model,badu,m)
  M<-ifelse(badu==T,ceiling(length(data)/log(length(data))),m)
  K<-gammas$K
  .end<-gammas$cs
  k1<-gammas$k1
  if(plot){
    .plotfun(K,  gammas$gamma[,.end], type="l", xlab="k", ylab=bquote("gamma"), main=main, plot=plot, add=add, ...)
  }
  if(add) {
    lines(K,gammas$gamma[,.end],lwd=2,...)
  }
  .output(list(k=K, weights = gammas$gamma[,(1:M)], gamma=gammas$gamma[,.end],k1=k1), plot=plot, add=add)
}


nloglike <- function(p,rho,z1,k) {
  n <- length(z1)
  g <- p[1];   del<-p[2];   tau<-rho/p[1]
  if (g > 0 && tau < 0 && del > max(-0.999,1/tau) && all((del*(1-z1^tau)+1) > 0)&&all(1+del*(1-(1+tau)*(z1^tau))>0)) {
    objf2 <- n*log(g)+(1/g+1)*sum(log(z1))+(1/g+1)*sum(log(1+del*(1-z1^tau)))-sum(log(1+del*(1-(1+tau)*(z1^tau))))
  } else { objf2 <- 1000000 }
  return((objf2))
}

RevEPD<-function(data,rho=-1,mif=1){
  Z<-sort(data)
  Hill.x<-rep(mean(ReIns::Hill(Z)$gamma[1:(0.5*length(Z))]),length(Z))
  startv<-cbind(Hill.x,-0.01)
  p<-startv
  n <- length(Z)
  Gamma <- matrix(nrow=n, ncol=2)
  
  for (k in (n-1):1){
    #print(k)
    z1<-Z[(n-k+1):n]/Z[n-k]
    optEPD <-  tryCatch({optim(p=p[k,],nloglike,z1=z1,k=k,rho=rho,
                               control=list(trace=0))$par},error=function(e) {return(c(NA,NA))})
    parsEPD<-as.vector(optEPD)
    Gamma[k,]<-parsEPD
  }
  return(Gamma[,1])
}
########################################################################################################################################################
EGPBBnon.fitMLE.boot = function(data,i,mopt,model){
  return( EGPBBnon.fitMLE.dq(data[i],mopt,model))
}

EGP.fitML.boot <- function(data,i,type=1,prob0=NA,kappa0=NA,delta0=NA,sigma0=NA,xi0=NA,censoring=c(0,Inf),rounded=0.1,print=FALSE){
  return( EGP.fitML(data[i],type=type,prob0=prob0,kappa0=kappa0,delta0=delta0,sigma0=sigma0,xi0=xi0,censoring=censoring,rounded=rounded,print=print) )
}