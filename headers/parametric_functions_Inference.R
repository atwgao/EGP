#########################
###                   ###
### Inference via PMW ###
###                   ###
#########################

PWM <- function(orders,prob=NA,kappa=NA,delta=NA,sigma=NA,xi=NA,type=1,empiric=FALSE,Unif=NULL,NbSamples=10^4,N=200){
  if(!empiric){
    if(type==0){
      return( sigma/((1+orders)*(1+orders-xi)) )
    } else if(type==1){
      be <- beta(c(1:(max(orders)+1))*kappa,1-xi)
      expect <- c()
      for(i in 1:length(orders)){
        expect[i] <- sum((-1)^c(0:orders[i])*choose(orders[i],c(0:orders[i]))*be[1:(orders[i]+1)])
      }
      return( (sigma/xi)*(kappa*expect - 1/(1+orders)) )
    } else if(type==2){
      res <- c()
      for(i in 1:length(orders)){
        l.vals <- c(0:orders[i])
        res[i] <- (sigma/xi)*( sum( (-1)^l.vals*choose(orders[i],l.vals)*(1+delta)^(orders[i]-l.vals+1)/((1+delta*l.vals+orders[i]-xi)*(1+delta*(l.vals+1)+orders[i]-xi)) )/delta^orders[i] - 1/(1+orders[i]) )
      }
      return( res )
    } else if(type==3){
      Expect <- function(kappa,delta,xi,s,N){
        l.vals <- matrix(c(0:s),nrow=N+1,ncol=s+1,byrow=TRUE)
        n.vals <- matrix(c(0:N),nrow=N+1,ncol=s+1)
        
        P1 <- kappa/(2*delta)
        P2 <- (-1)^(n.vals+l.vals)
        P3 <- choose(s,l.vals)
        P4 <- choose((kappa/2)*(l.vals+1)-1,n.vals)
        P5 <- (1+1/delta)^(n.vals+1)
        P6 <- (1+delta)^((n.vals-xi+1)/delta)
        P7 <- P6*beta((n.vals-xi+1)/delta,n.vals+1)*pbeta(1/(1+delta),(n.vals-xi+1)/delta,n.vals+1)
        P8 <- (1+delta)*P6*beta((n.vals-xi+1)/delta+1,n.vals+1)*pbeta(1/(1+delta),(n.vals-xi+1)/delta+1,n.vals+1)
        
        ToSum <- rowSums(P2*P3*P4*P5*(P7-P8))
        
        if(!is.nan(ToSum[N+1]/ToSum[N])){
          power <- log(ToSum[N+1]/ToSum[N])/log(1+1/N)
          add <- as.numeric(ToSum[N+1]*(zeta(-power)*N^(-power)-sum((c(1:N)/N)^power)))
          return( P1*(sum(P2*P3*P4*P5*(P7-P8))+add) )
        } else{
          return( P1*sum(P2*P3*P4*P5*(P7-P8)) )	
        }
      }
      
      res <- c()
      for(i in 1:length(orders)){
        res[i] <- (sigma/xi)*(Expect(kappa,delta,xi,orders[i],N) - 1/(1+orders[i]))
      }
      
      return( res )
    } else if(type==4){
      res <- c()
      for(i in 1:length(orders)){
        l.vals <- matrix(c(0:orders[i]),nrow=orders[i]+1,ncol=orders[i]+1,byrow=TRUE)
        m.vals <- matrix(c(0:orders[i]),nrow=orders[i]+1,ncol=orders[i]+1)
        inds <- l.vals<m.vals
        l.vals[inds] <- m.vals[inds] <- NA	
        const <- choose(orders[i],l.vals)*choose(l.vals,m.vals)*(-1)^l.vals*prob^m.vals*(1-prob)^(l.vals-m.vals)
        be1 <- beta(kappa*(m.vals+1)+delta*(l.vals-m.vals),1-xi)
        be2 <- beta(kappa*m.vals+delta*(l.vals-m.vals+1),1-xi)
        expect <- sum(const*(prob*kappa*be1 + (1-prob)*delta*be2),na.rm=TRUE)
        res[i] <- (sigma/xi)*(expect - 1/(1+orders[i]))
      }
      return(res)
    }
  } else{
    if(is.null(Unif)){
      Unif <- runif(NbSamples)
    }
    XG <- rG(length(Unif),prob,kappa,delta,type,Unif,direct=TRUE)
    X <- evd:::qgpd(XG,scale=sigma,shape=xi)
    res <- c()
    for(i in 1:length(orders)){
      res[i] <- mean(X*(1-Unif)^orders[i],na.rm=TRUE)
    }
    return( res ) 
  }
}

EGP.fitPWM <- function(x,type=1,prob0=NA,kappa0=NA,delta0=NA,sigma0=NA,xi0=NA,empiric=FALSE,Unif=NULL,NbSamples=10^4,N=200){
  Fn = ecdf(x)
  mu0hat = mean(x)
  mu1hat = mean(x*(1-Fn(x)))
  mu2hat = mean(x*(1-Fn(x))^2)
  mu3hat = mean(x*(1-Fn(x))^3)
  mu4hat = mean(x*(1-Fn(x))^4)
  
  if(type==0){
    xihat <- (mu0hat-4*mu1hat)/(mu0hat-2*mu1hat)
    sigmahat <- mu0hat*(1-xihat)
    
    thetahat <- c(sigmahat,xihat)
    names(thetahat) <- c("sigma","xi")
    return(thetahat)
  } else if(type==1){
    fct <- function(theta,x){
      pwm.theor <- PWM(orders=c(0:2),kappa=theta[1],sigma=theta[2],xi=theta[3],type=1,empiric=empiric,Unif=Unif,NbSamples=NbSamples)
      pwm.empir <- c(mu0hat,mu1hat,mu2hat)
      return(matrix(pwm.theor - pwm.empir,ncol=3))
    }
    theta0 <- c(kappa0,sigma0,xi0)
    res <- gmm(fct, x, theta0, optfct = "nlminb", lower = c(0.0001, 0.0001, 10^(-6)),  upper = c(Inf,Inf, .99),vcov="iid")
    thetahat <- res$coefficients
    names(thetahat) <- c("kappa","sigma","xi")
    return(thetahat)
  } else if(type==2){
    fct <- function(theta,x){
      pwm.theor <- PWM(orders=c(0:2),delta=theta[1],sigma=theta[2],xi=theta[3],type=2,empiric=empiric,Unif=Unif,NbSamples=NbSamples)
      pwm.empir <- c(mu0hat,mu1hat,mu2hat)
      return(matrix(pwm.theor - pwm.empir,ncol=3))
    }
    theta0 <- c(delta0,sigma0,xi0)
    res <- gmm(fct, x, theta0, optfct = "nlminb", lower = c(0.0001, 0.0001, 10^(-6)),  upper = c(100,Inf, .99),vcov="iid")
    thetahat <- res$coefficients
    names(thetahat) <- c("delta","sigma","xi")
    return(thetahat)
  } else if(type==3){
    fct <- function(theta,x){
      pwm.theor <- PWM(orders=c(0:3),kappa=theta[1],delta=theta[2],sigma=theta[3],xi=theta[4],type=3,empiric=empiric,Unif=Unif,NbSamples=NbSamples,N=N)
      pwm.empir <- c(mu0hat,mu1hat,mu2hat,mu3hat)
      return(matrix(pwm.theor - pwm.empir,ncol=4))
    }
    theta0 <- c(kappa0,delta0,sigma0,xi0)
    if(empiric){
      res <- gmm(fct, x, theta0, optfct = "nlminb", lower = c(0.0001, 0.0001, 0.0001, 10^(-6)),  upper = c(Inf,100,Inf, .99),vcov="iid")
    } else{
      res <- gmm(fct, x, theta0, optfct = "nlminb", lower = c(1.5, 0.5, 0.01, 0.05),  upper = c(10,10,100, 0.5),vcov="iid")
    }
    thetahat <- res$coefficients
    names(thetahat) <- c("kappa","delta","sigma","xi")
    return(thetahat)
  } else if(type==4){
    fct <- function(theta,x){
      pwm.theor <- PWM(orders=c(0:4),prob=theta[1],kappa=theta[2],delta=theta[2]+theta[3],sigma=theta[4],xi=theta[5],type=4,empiric=empiric,Unif=Unif,NbSamples=NbSamples)
      pwm.empir <- c(mu0hat,mu1hat,mu2hat,mu3hat,mu4hat)
      return(matrix(pwm.theor - pwm.empir,ncol=5))
    }
    res0 <- EGP.fitPWM(x[x>quantile(x,0.95)]-quantile(x,0.95),type=0)
    if(is.na(prob0)){ prob0 <- 0.5 }
    if(is.na(kappa0)){ kappa0 <- mean(x[x<quantile(x,0.05)])/(quantile(x,0.05)-mean(x[x<quantile(x,0.05)])) }
    if(is.na(delta0)){ delta0 <- kappa0+0.01}
    Ddelta0 <- delta0-kappa0 
    if(is.na(sigma0)){ sigma0 <- res0[1] }
    if(is.na(xi0)){ xi0 <- res0[2] }
    theta0 <- c(prob0,kappa0,Ddelta0,sigma0,xi0)
    names(theta0) <- c("prob","kappa","Ddelta","sigma","xi")
    #print(theta0)
    res <- gmm(fct, x, theta0, optfct = "nlminb", lower = c(0, 0.0001, 0, 0.0001, 10^(-6)),  upper = c(1, Inf, Inf,Inf, .99),vcov="iid")
    thetahat <- res$coefficients; thetahat[3] <- thetahat[2] + thetahat[3]
    names(thetahat) <- c("prob","kappa","delta","sigma","xi")
    return(thetahat)
  }
}


EGP.fitPWM.boot <- function(data,i,type=1,prob0=NA,kappa0=NA,delta0=NA,sigma0=NA,xi0=NA,empiric=FALSE,Unif=NULL,NbSamples=10^4,N=200){
  return( EGP.fitPWM(data[i],type=type,prob0=prob0,kappa0=kappa0,delta0=delta0,sigma0=sigma0,xi0=xi0,empiric=empiric,Unif=Unif,NbSamples=NbSamples,N=N) )
}


########################
###                  ###
### Inference via ML ###
###                  ###
########################

neg.log.lik <- function(theta,x,censoring,rounded,type){
  if(type==0){
    if(theta[1]<=0 | theta[2]<=10^(-6) | theta[2]>0.99 ){
      return(Inf)
    } else{
      if(rounded==0){
        censor<-0
        if(censoring[1]>0){
          if(censoring[2]<Inf){
            censor <- pEGP(censoring[2],sigma=theta[1],xi=theta[2],type=0)-pEGP(censoring[1],sigma=theta[1],xi=theta[2],type=0)
          } else{
            censor <- 1-pEGP(censoring[1],sigma=theta[1],xi=theta[2],type=0)
          }
        } else{
          if(censoring[2]<Inf){
            censor <- pEGP(censoring[2],sigma=theta[1],xi=theta[2],type=0)
          } else{
            censor<- 1
          }
        }
        return( -sum(dEGP(x,sigma=theta[1],xi=theta[2],type=0,log=TRUE)-log(censor),na.rm=TRUE) )
      } else if(rounded>0){
        return( -sum(log(pEGP(x+rounded,sigma=theta[1],xi=theta[2],type=0)-pEGP(x,sigma=theta[1],xi=theta[2],type=0))-log(1-pEGP(rounded,sigma=theta[1],xi=theta[2],type=0))) )
      }
    }
  } else if(type==1){
    if(theta[1]<=0 | theta[2]<=0  | theta[3]<=10^(-6) | theta[3]>0.99){
      return(Inf)
    } else{
      if(rounded==0){
        censor<-0
        if(censoring[1]>0){
          if(censoring[2]<Inf){
            censor <- pEGP(censoring[2],kappa=theta[1],sigma=theta[2],xi=theta[3],type=1)-pEGP(censoring[1],kappa=theta[1],sigma=theta[2],xi=theta[3],type=1)
          } else{
            censor <- 1-pEGP(censoring[1],kappa=theta[1],sigma=theta[2],xi=theta[3],type=1)
          }
        } else{
          if(censoring[2]<Inf){
            censor <- pEGP(censoring[2],kappa=theta[1],sigma=theta[2],xi=theta[3],type=1)
          } else{
            censor<- 1
          }
        }
        return( -sum(dEGP(x,kappa=theta[1],sigma=theta[2],xi=theta[3],type=1,log=TRUE)-log(censor),na.rm=TRUE) )
      } else if(rounded>0){
        return( -sum(log(pEGP(x+rounded,kappa=theta[1],sigma=theta[2],xi=theta[3],type=1)-pEGP(x,kappa=theta[1],sigma=theta[2],xi=theta[3],type=1))-log(1-pEGP(rounded,kappa=theta[1],sigma=theta[2],xi=theta[3],type=1))) )
      }
    }
  } else if(type==2){
    if(theta[1]<=0 | theta[1]>100 | theta[2]<=0 | theta[3]<=10^(-6) | theta[3]>0.99){
      return(Inf)
    } else{
      if(rounded==0){
        censor<-0
        if(censoring[1]>0){
          if(censoring[2]<Inf){
            censor <- pEGP(censoring[2],delta=theta[1],sigma=theta[2],xi=theta[3],type=2)-pEGP(censoring[1],delta=theta[1],sigma=theta[2],xi=theta[3],type=2)
          } else{
            censor <- 1-pEGP(censoring[1],delta=theta[1],sigma=theta[2],xi=theta[3],type=2)
          }
        } else{
          if(censoring[2]<Inf){
            censor <- pEGP(censoring[2],delta=theta[1],sigma=theta[2],xi=theta[3],type=2)
          } else{
            censor<- 1
          }
        }
        return( -sum(dEGP(x,delta=theta[1],sigma=theta[2],xi=theta[3],type=2,log=TRUE)-log(censor),na.rm=TRUE) )
      } else if(rounded>0){
        return( -sum(log(pEGP(x+rounded,delta=theta[1],sigma=theta[2],xi=theta[3],type=2)-pEGP(x,delta=theta[1],sigma=theta[2],xi=theta[3],type=2))-log(1-pEGP(rounded,delta=theta[1],sigma=theta[2],xi=theta[3],type=2))) )
      }
    }
  } else if(type==3){
    if(theta[1]<=0 | theta[2]<=0 | theta[2]>100 | theta[3]<=0 | theta[4]<=10^(-6) | theta[4]>0.99){
      return(Inf)
    } else{
      if(rounded==0){
        censor<-0
        if(censoring[1]>0){
          if(censoring[2]<Inf){
            censor <- pEGP(censoring[2],kappa=theta[1],delta=theta[2],sigma=theta[3],xi=theta[4],type=3)-pEGP(censoring[1],kappa=theta[1],delta=theta[2],sigma=theta[3],xi=theta[4],type=3)
          } else{
            censor <- 1-pEGP(censoring[1],kappa=theta[1],delta=theta[2],sigma=theta[3],xi=theta[4],type=3)
          }
        } else{
          if(censoring[2]<Inf){
            censor <- pEGP(censoring[2],kappa=theta[1],delta=theta[2],sigma=theta[3],xi=theta[4],type=3)
          } else{
            censor<- 1
          }
        }
        return( -sum(dEGP(x,kappa=theta[1],delta=theta[2],sigma=theta[3],xi=theta[4],type=3,log=TRUE)-log(censor),na.rm=TRUE) )
      } else if(rounded>0){
        return( -sum(log(pEGP(x+rounded,kappa=theta[1],delta=theta[2],sigma=theta[3],xi=theta[4],type=3)-pEGP(x,kappa=theta[1],delta=theta[2],sigma=theta[3],xi=theta[4],type=3))-log(1-pEGP(rounded,kappa=theta[1],delta=theta[2],sigma=theta[3],xi=theta[4],type=3))) )
      }
    }
  } else if(type==4){
    if(theta[1]<0 | theta[1]>1 | theta[2]<=0 | theta[3]<=0 | theta[4]<=0 | theta[5]<=10^(-6) | theta[5]>0.99 | theta[3]<theta[2]){
      return(Inf)
    } else{
      if(rounded==0){
        censor<-0
        if(censoring[1]>0){
          if(censoring[2]<Inf){
            censor <- pEGP(censoring[2],prob=theta[1],kappa=theta[2],delta=theta[3],sigma=theta[4],xi=theta[5],type=4)-pEGP(censoring[1],prob=theta[1],kappa=theta[2],delta=theta[3],sigma=theta[4],xi=theta[5],type=4)
          } else{
            censor <- 1-pEGP(censoring[1],prob=theta[1],kappa=theta[2],delta=theta[3],sigma=theta[4],xi=theta[5],type=4)
          }
        } else{
          if(censoring[2]<Inf){
            censor <- pEGP(censoring[2],prob=theta[1],kappa=theta[2],delta=theta[3],sigma=theta[4],xi=theta[5],type=4)
          } else{
            censor<- 1
          }
        }
        return( -sum(dEGP(x,prob=theta[1],kappa=theta[2],delta=theta[3],sigma=theta[4],xi=theta[5],type=4,log=TRUE)-log(censor),na.rm=TRUE) )
      } else if(rounded>0){
        return( -sum(log(pEGP(x+rounded,prob=theta[1],kappa=theta[2],delta=theta[3],sigma=theta[4],xi=theta[5],type=4)-pEGP(x,prob=theta[1],kappa=theta[2],delta=theta[3],sigma=theta[4],xi=theta[5],type=4))-log(1-pEGP(rounded,prob=theta[1],kappa=theta[2],delta=theta[3],sigma=theta[4],xi=theta[5],type=4))) )
      }
    }
  }
}

EGP.fitML <- function(x,type=1,prob0=NA,kappa0=NA,delta0=NA,sigma0=NA,xi0=NA,censoring=c(0,Inf),rounded=0.1,print=FALSE){
  x <- x[x>=censoring[1] & x<=censoring[2]]
  if(type==0){
    theta0 <- c(sigma0,xi0)
    opt <- optim(par=theta0,fn=neg.log.lik,x=x,censoring=censoring,rounded=rounded,type=type,method="Nelder-Mead",control=list(maxit=1000),hessian=FALSE)
    if(print){
      print(opt)
    }
    thetahat <- opt$par
    names(thetahat) <- c("sigma","xi")
    return(thetahat)
  } else if(type==1){
    theta0 <- c(kappa0,sigma0,xi0)
    opt <- optim(par=theta0,fn=neg.log.lik,x=x,censoring=censoring,rounded=rounded,type=type,method="Nelder-Mead",control=list(maxit=1000),hessian=FALSE)
    if(print){
      print(opt)
    }
    thetahat <- opt$par
    names(thetahat) <- c("kappa","sigma","xi")
    return(thetahat)
  } else if(type==2){
    theta0 <- c(delta0,sigma0,xi0)
    opt <- optim(par=theta0,fn=neg.log.lik,x=x,censoring=censoring,rounded=rounded,type=type,method="Nelder-Mead",control=list(maxit=1000),hessian=FALSE)
    if(print){
      print(opt)
    }
    thetahat <- opt$par
    names(thetahat) <- c("delta","sigma","xi")
    return(thetahat)
  } else if(type==3){
    theta0 <- c(kappa0,delta0,sigma0,xi0)
    opt <- optim(par=theta0,fn=neg.log.lik,x=x,censoring=censoring,rounded=rounded,type=type,method="Nelder-Mead",control=list(maxit=1000),hessian=FALSE)
    if(print){
      print(opt)
    }
    thetahat <- opt$par
    names(thetahat) <- c("kappa","delta","sigma","xi")
    return(thetahat)
  } else if(type==4){
    theta0 <- c(prob0,kappa0,delta0,sigma0,xi0)
    opt <- optim(par=theta0,fn=neg.log.lik,x=x,censoring=censoring,rounded=rounded,type=type,method="Nelder-Mead",control=list(maxit=1000),hessian=FALSE)
    if(print){
      print(opt)
    }
    thetahat <- opt$par
    names(thetahat) <- c("prob","kappa","delta","sigma","xi")
    return(thetahat)
  }
}

EGP.fitML.boot <- function(data,i,type=1,prob0=NA,kappa0=NA,delta0=NA,sigma0=NA,xi0=NA,censoring=c(0,Inf),rounded=0.1,print=FALSE){
  return( EGP.fitML(data[i],type=type,prob0=prob0,kappa0=kappa0,delta0=delta0,sigma0=sigma0,xi0=xi0,censoring=censoring,rounded=rounded,print=print) )
}

















