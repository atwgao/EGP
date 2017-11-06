###########Parametric Distributions########
###########################################


model1<-function(n,t1,t2){
  res<-((1-runif(n))^(-1/t2)-1)^(-1/t1)
  return(max(res)+200-res)}

model2<- function(n,l){return(1-exp(-rgamma(n,2,l)))}
model3 <- function(n,t1,t2){return((((1-runif(n))^(-1/t2)-1)^(-1/t1)+1)^(-1))}
model4 <- function(n,g){return(1-(1-runif(n))^(-g))}


# Pareto

dpareto <- function(x, gamma, log = FALSE) {
  if (gamma<=0) {
    stop("gamma should be strictly positive.")
  }

  d <- 1/gamma*(1/x)^(1/gamma+1)
  
  if (log) d <- log(d)
  
  return(d)
  
}
###
ppareto <- function(x, gamma, lower.tail = TRUE, log.p = FALSE) {
  if (gamma<=0) {
    stop("shape should be strictly positive.")
  }


  p <-  1-(1/x)^(1/gamma)
  
  if (!lower.tail) p <- 1-p
  
  if (log.p) p <- log(p)
  
  
  return(p)
  
}

### GPD
dgpd <- function(x, loc = 0, scale = 1, shape = 0, log.d = FALSE) {
  if(min(scale) <= 0)
    stop("Invalid scale")
  cond1 <- (length(x) > 1) &
    (((length(loc) != length(x)) & (length(loc) != 1)) |
       ((length(scale) != length(x)) & (length(scale) != 1)) |
       ((length(shape) != length(x)) & (length(shape) != 1)))
  cond2 <- (length(x) == 1) &
    (length(unique(c(length(x), length(loc), length(scale), length(shape)))) > 2)
  if(cond1 | cond2)
    stop("Invalid parameter length!")
  if(length(shape) == 1)
    shape <- rep(shape, max(length(x), length(loc), length(scale)))
  below.support <- x < loc
  x <- pmax(x, loc)
  x <- ifelse(shape >= 0, x, pmin(x, (loc - scale/shape)))
  w <- (x - loc) / scale
  log.density <- -log(scale) - ifelse(shape == 0, w, ((1/shape) + 1) * log1p(w * shape))
  log.density[is.nan(log.density) | is.infinite(log.density) | below.support] <- -Inf
  if(!log.d)
    log.density <- exp(log.density)
  log.density
}

rgpd <- function(n, loc = 0, scale = 1, shape = 0) {
  if(min(scale) <= 0)
    stop("Invalid scale")
  cond1 <- (n > 1) &
    (((length(loc) != n) & (length(loc) != 1)) |
       ((length(scale) != n) & (length(scale) != 1)) |
       ((length(shape) != n) & (length(shape) != 1)))
  cond2 <- (n == 1) &
    (length(unique(c(n, length(loc), length(scale), length(shape)))) > 2)
  if(cond1 | cond2)
    stop("Invalid parameter length!")
  qgpd(runif(n), loc, scale, shape)
}


qgpd <- function(p, loc = 0, scale = 1, shape = 0, lower.tail = TRUE, log.p = FALSE) {
  if(log.p)
    p <- exp(p)
  if((min(p, na.rm = TRUE) <= 0) || (max(p, na.rm = TRUE) >= 1))
    stop("`p' must contain probabilities in (0,1)")
  if(min(scale) <= 0)
    stop("Invalid scale")
  cond1 <- (length(p) > 1) &
    (((length(loc) != length(p)) & (length(loc) != 1)) |
       ((length(scale) != length(p)) & (length(scale) != 1)) |
       ((length(shape) != length(p)) & (length(shape) != 1)))
  cond2 <- (length(p) == 1) &
    (length(unique(c(length(p), length(loc), length(scale), length(shape)))) > 2)
  if(cond1 | cond2)
    stop("Invalid parameter length!")
  if(lower.tail)
    p <- 1 - p
  if(length(shape) == 1)
    shape <- rep(shape, max(length(p), length(loc), length(scale)))
  ifelse(shape == 0, loc - scale * log(p), loc + scale * expm1(-shape * log(p)) / shape)
}


pgpd <- function(q, loc = 0, scale = 1, shape = 0, lower.tail = TRUE, log.p = FALSE) {
  if(min(scale) <= 0)
    stop("Invalid scale")
  cond1 <- (length(q) > 1) &
    (((length(loc) != length(q)) & (length(loc) != 1)) |
       ((length(scale) != length(q)) & (length(scale) != 1)) |
       ((length(shape) != length(q)) & (length(shape) != 1)))
  cond2 <- (length(q) == 1) &
    (length(unique(c(length(q), length(loc), length(scale), length(shape)))) > 2)
  if(cond1 | cond2)
    stop("Invalid parameter length!")
  if(length(shape) == 1)
    shape <- rep(shape, max(length(q), length(loc), length(scale)))
  q <- pmax(q, loc)
  q <- ifelse(shape >= 0, q, pmin(q, (loc - scale/shape)))
  w <- (q - loc) / scale
  p <- ifelse(shape == 0, 1 - exp(-w), 1 - exp((-1/shape)*log1p(w*shape)))
  if(!lower.tail)
    p <- 1 - p
  if(log.p)
    p <- log(p)
  p
}

qgpd <- function (p, loc = 0, scale = 1, shape = 0, lower.tail = TRUE) 
{
  inds <- p <=0 | p>=1
  res <- rep(NA,length(p))
  if (min(scale) < 0) 
    stop("invalid scale")
  if (length(shape) != 1) 
    stop("invalid shape")
  if (lower.tail) 
    p <- 1 - p
  if (shape == 0)
    res[!inds] <- loc - scale * log(p[!inds])
  else res[!inds] <- loc + scale * (p[!inds]^(-shape) - 1)/shape
  
  return(res)
}

qgpd.fullrange <- function (p, loc = 0, scale = 1, shape = 0, prob.loc=0.95) 
{
  res <- loc + qgpd((p-prob.loc)/(1-prob.loc),scale=scale,shape=shape)
  return(res)
}


pG <- function(H,prob=NA,kappa=NA,delta=NA,type=1){
  if(type==1){
    G<-H^kappa
  } else if(type==4){
    G <- prob*H^kappa + (1-prob)*H^delta
  }
  return(G)
}


dG <- function(H,prob=NA,kappa=NA,delta=NA,type=1,log=FALSE){
  if(log==FALSE){
    if(type==1){
      g<-kappa*H^(kappa-1)
    } else if(type==4){
      g <- prob*kappa*H^(kappa-1) + (1-prob)*delta*H^(delta-1)
    }
  } else{
    if(type==1){
      g<-log(kappa) + (kappa-1)*log(H)
    } else if(type==4){
      g<-log(prob*kappa*H^(kappa-1) + (1-prob)*delta*H^(delta-1))
    }
  }
  return(g)
}

qG <- function(H,prob=NA,kappa=NA,delta=NA,type=1){
  if(type==1){
    q<-H^(1/kappa)
  } else if(type==4){
    dummy.func <- function(H,p,prob=NA,kappa=NA,delta=NA){
      return(pG(H=H,prob=prob,kappa=kappa,delta=delta,type=4)-p )
    }
    find.root <- function(H,prob=NA,kappa=NA,delta=NA){
      return( uniroot(dummy.func,interval=c(0,1),p=H,prob=prob,kappa=kappa,delta=delta)$root )
    }
    q<-sapply(H,FUN=find.root,prob=prob,kappa=kappa,delta=delta)
  }
  return(q)
}

pExt <- function(x,prob=NA,kappa=NA,delta=NA,sigma=NA,xi=NA,type=1,model="GP"){
  if(model!="GP" && model!="Pa") stop("Extensions limited to Pa and GP")
  if(model=="GP"){
    H<-pgpd(x,scale=sigma,shape=xi)
  }else{
    H<-actuar::ppareto(x,shape=xi,scale=1) 
  } 
  return( pG(H,prob,kappa,delta,type) )
}

dExt <- function(x,prob=NA,kappa=NA,delta=NA,sigma=NA,xi=NA,type=1,log=FALSE,model="GP"){
  if(model!="GP" && model!="Pa") stop("Extensions limited to Pa and GP")
  if(model=="GP"){
    H<-pgpd(x,scale=sigma,shape=xi)
    h<-dgpd(x,scale=sigma,shape=xi)
  }else{
    H<-actuar::ppareto(x,shape=xi,scale=1) 
    h<-actuar::dpareto(x,shape=xi,scale=1)
  } 
  
  if(log==FALSE){
    return(dG(H,prob,kappa,delta,type)*h )
  } else{
    return(dG(H,prob,kappa,delta,type,log=TRUE) + log(h) )
  }
}

qExt <- function(p,prob=NA,kappa=NA,delta=NA,sigma=NA,xi=NA,type=1,model="GP"){
  qG <- qG(p,prob,kappa,delta,type)
  if(model!="GP" && model!="Pa") stop("Extensions limited to Pa and GP")
  if(model=="GP"){
    qH <- qgpd(qG,scale=sigma,shape=xi)
  }else{
    qH <- actuar::qpareto(qG,shape=xi,scale=1) 
  } 
  return(qH)
}

rExt <- function(n,prob=NA,kappa=NA,sigma=NA,delta=NA,xi=NA,type=1,Unif=NULL,model="GP"){
  rGx<-rG(n,prob,kappa,delta,type,Unif)
  if(model!="GP" && model!="Pa") stop("Extensions limited to Pa and GP")
  if(model=="GP"){
    x<-qgpd(rGx,scale=sigma,shape=xi) 
  }else{
    x<-actuar::qpareto(rGx,shape=1/xi,scale=1) 
  } 
  return(x)
}

rG <- function(n,prob=NA,kappa=NA,delta=NA,type=1,Unif=NULL,direct=FALSE){
  if(is.null(Unif)){
    Unif <- runif(n)
  }
  if(type!=4 | (type==4 & direct)){
    return( qG(Unif,prob,kappa,delta,type) )		
  } else if(type==4 & !direct){
    components <- sample(x=c(1,2),size=n,replace=TRUE,prob=c(prob,1-prob))
    res <- c()
    res[components==1] <- qG(Unif[components==1],prob=NA,kappa=kappa,delta=NA,type=1)
    res[components==2] <- qG(Unif[components==2],prob=NA,kappa=delta,delta=NA,type=1)
    return(res)
  }
}

####################################WITH DELTA.K ########################################
pExt.delta <- function(x,prob=NA,kappa=NA,delta=NA,sigma=NA,xi=NA,delta.k=NA,type=1,model="GP"){
  if(model!="GP" && model!="Pa") stop("Extensions limited to Pa and GP")
  if(model=="GP"){
    H<-pgpd(x,scale=sigma,shape=xi)
  }else{
    H<-1-actuar::ppareto(x,shape=xi,scale=1) 
  } 
  return(1-delta.k+delta.k*H*pG(H,prob,kappa,delta,type) )
}

dExt.delta <- function(x,prob=NA,kappa=NA,delta=NA,sigma=NA,xi=NA,delta.k=NA,type=1,log=FALSE,model="GP"){
  if(model!="GP" && model!="Pa") stop("Extensions limited to Pa and GP")
  if(model=="GP"){
    H<-pgpd(x,scale=sigma,shape=xi)
    h<-dgpd(x,scale=sigma,shape=xi)
  }else{
    H<-1-actuar::ppareto(x,shape=xi,scale=1) 
    h<-actuar::dpareto(x,shape=xi,scale=1)
  } 
  
  if(log==FALSE){
    return(1-delta.k+delta.k*dG(H,prob,kappa,delta,type)*h )
  } else{
    if(all((1-delta.k+delta.k*dG(H,prob,kappa,delta,type,log=FALSE))>0)){
    return(log(1-delta.k+delta.k*dG(H,prob,kappa,delta,type,log=FALSE)) + log(h) )
    }else{
      return(dG(H,prob,kappa,delta,type,log=TRUE)+log(h))
    }
  }
}



#######################################################
###############Non Parametric distributions############
########################################################


dBB = function(w,u){
  m_w = length(w);  
  fct = numeric(length(u));
  for (i in 1:m_w){
    fct = fct + w[i] * dbeta(u,shape1 = i, shape2 = m_w-i+1)
  }
  return(fct)
}

dBBnon= function(m,u,kap){
  n = length(u)
  fct = numeric(n);
  w = numeric(m);
  for (i in 1:m){
    Fn = sum(u<=(i/m))/n;
    Fn_1 = sum(u<=((i-1)/m))/n;
    w[i] = Fn-Fn_1;
    if ((i==m) && w[i]==0){
      w[i] = max(Fn-Fn_1,kap/m);
      w = w/sum(w)  
    }
  }
  for (i in 1:m)
    fct = fct + w[i] * dbeta(u,shape1 = i, shape2 = m-i+1)
  
  return(fct)
}

pBB = function(w,u){
  m_w = length(w)
  fct = numeric(length(u))
  for (i in 1:m_w){
    fct = fct + w[i] * pbeta(u,shape1 = i, shape2 = m_w-i+1)
  }
  return(fct)
}
qBB = function(w,u){
  dummy.func = function(u,p,w=NA){
    return( pBB(w=w,u=u)-p)
  }
  find.root = function(u,w=NA){
    return(uniroot(dummy.func,interval=c(0,1),p=u,w=w,tol=1e-08)$root)
  }
  s = sapply(u,FUN=find.root,w=w)
  s[s==0] = min(s[s!=0])
  return(s)
}

dEGP.BB = function(x,param,model){
  pn = length(param);
  if(model!="GP" && model!="Pa") stop("Extensions limited to Pa and GP")
  if(model=="GP"){
  H = evd:::pgpd(q=x, scale=param[pn-1], shape=param[pn])
  h = evd:::dgpd(x=x, scale=param[pn-1], shape=param[pn])
  g = dBB(w=param[1:(pn-2)],u=H)
  }else{
    H = actuar:::ppareto(q=x, shape=param[pn], scale=1)
    h = actuar:::dpareto(x=x, shape=param[pn], scale=1)
    g = dBB(w=param[1:(pn-1)],u=(H))
  }
  fx = g*h
  return(fx=fx)
}

pEGP.BB = function(x,param,model){
  pn = length(param);
  if(model!="GP" && model!="Pa") stop("Extensions limited to Pa and GP")
  if(model=="GP"){
  H = evd:::pgpd(q=x, scale=param[pn-1], shape=param[pn])
  G = pBB(w=param[1:(pn-2)],u=H)
  }else{
    H = actuar:::ppareto(x, shape=param[pn], scale=1)
    G = pBB(w=param[1:(pn-1)],u=(1-H))
  }
  return(Fx=G)
}
qEGP.BB = function(x,param,model){
  pn = length(param)
  if(model!="GP" && model!="Pa") stop("Extensions limited to Pa and GP")
  if(model=="GP"){
  qx = evd:::qgpd(qBB(w=param[1:(pn-2)],u = x),scale=param[pn-1],shape=param[pn])
  }else{
  qx = actuar:::qpareto(qBB(w=param[1:(pn-1)],u = x),scale=1,shape=param[pn])
  }
  return(qx=qx)
}
