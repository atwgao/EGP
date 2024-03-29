
# This file contains implementation of QQ-plots and their derivative plots.

##############################################################################################

# Computes the quantile scores to the Exponential quantile plot
# (Section 1.2.1 in Beirlant et al. (2004)) for a numeric vector of observations (data)
#
# If plot=TRUE then the scores are plotted on a Exponential
# quantile plot

ExpQQ <- function(data, plot = TRUE, main = "Exponential QQ-plot", ...) {
  
  # Check input arguments
  .checkInput(data, pos=FALSE)
  
  X <- as.numeric(sort(data))
  n <- length(X)
  
  # calculate theoretical and empirical quantiles
  
  eqq.the <- -log(1 - (1:n)/(n+1))
  eqq.emp <- X
  
  # plots if TRUE
  
  .plotfun(eqq.the, eqq.emp, type="p", xlab="Quantiles of Standard Exponential", ylab="X", 
          main=main, plot=plot, add=FALSE, ...)
  abline(lm(eqq.emp~eqq.the))
  # output list with theoretical quantiles eqq.the and empirical quantiles eqq.emp
  .output(list(eqq.the=eqq.the, eqq.emp=eqq.emp), plot=plot, add=FALSE)
  
}

# Computes the mean excess scores e_k,n (Section 1.2.2 in Beirlant et al. (2004)) 
# for a numeric vector of observations (data)
#
# If plot=TRUE then the mean excess scores are plotted as
# a function of k
#
# If plot=TRUE and k=TRUE then the mean excess scores are
# plotted as a function of the order statistics X_n-k,n


MeanExcess <- function(data, plot = TRUE, k = FALSE, main = "Mean excess plot", ...) {
  
  # Check input arguments
  .checkInput(data, pos=FALSE)
  
  
  X <- as.numeric(sort(data))
  n <- length(X)
  K <- 1:(n-1)
  e <- numeric(n)
  
  # estimation of mean excess values
  
  e[K] <- cumsum(X[n-K+1])/K - X[n-K]
  
  # plots if TRUE
  
  if (plot) {
    if (k) {  	# as function of k
      .plotfun(K, e[K], type="p", xlab="k", ylab=bquote(e["k,n"]), main=main, plot=TRUE, add=FALSE, ...)
    } else { 	   	# as function of order statistics X_n-k,n
      .plotfun(X[n-K], e[K], type="b", xlab=bquote(X["n-k,n"]), ylab=bquote(e["k,n"]), main=main, plot=TRUE, add=FALSE, ...)
    }
  }
  
  # output list with values of k, order statistics X_n-k,n 
  # and mean excess scores e_k,n
  .output(list(k=K, X=X[n-K], e=e[K]), plot=plot, add=FALSE)
  
}


##############################################################################################

# Computes the quantile scores to the Pareto quantile plot
# (Section 1.2.1 in Beirlant et al. (2004)) for a numeric vector of observations (data)
#
# If plot=TRUE then the scores are plotted on a Pareto
# quantile plot

ParetoQQ <- function(data, plot = TRUE, main = "Pareto QQ-plot", ...) {
	
  # Check input arguments
  .checkInput(data)
  
	X <- as.numeric(sort(data))
	n <- length(X)

  # calculate theoretical and empirical quantiles
  
  pqq.the <- -log(1 - (1:n)/(n+1))
  pqq.emp <- log(X)
  
  # plots if TRUE
	.plotfun(pqq.the, pqq.emp, type="p", xlab="Quantiles of Standard Exponential", ylab="log(X)", 
          main=main, plot=plot, add=FALSE, ...)
	abline(lm(pqq.emp~pqq.the))
  # output list with theoretical quantiles pqq.the and empirical quantiles pqq.emp
	.output(list(pqq.the=pqq.the, pqq.emp=pqq.emp), plot=plot, add=FALSE)

}

# Derivative plot of Pareto QQ-plot
ParetoQQ_der <- function(data, k = FALSE, plot = TRUE, main = "Derivative plot of Pareto QQ-plot", ...) {
  
  # Check input arguments
  .checkInput(data)
  
  X <- as.numeric(sort(data))
  n <- length(X)
  K <- 1:(n-1)
  

  if (k) {
    xval <- K
    xlab <- "k"
  } else {
    xval <- log(X[n-K])
    xlab <- bquote(log(X["n-k,n"]))
  }

  # Derivative
  yval <- Hill(X, plot=FALSE)$gamma
  
  # plots if TRUE
  .plotfun(xval, yval, type="p", xlab=xlab, ylab="Derivative", 
          main=main, plot=plot, add=FALSE, ...)
  
  # output list with theoretical quantiles pqq.the and empirical quantiles pqq.emp
  .output(list(xval=xval, yval=yval), plot=plot, add=FALSE)
  
}

#######################################################################################


# Log-normal QQ-plot
LognormalQQ <- function(data, plot = TRUE, main = "Log-normal QQ-plot", ...) {
  
  # Check input arguments
  .checkInput(data)
  
  X <- as.numeric(sort(data))
  n <- length(X)
  
  # calculate theoretical and empirical quantiles
  
  lnqq.the <- qnorm((1:n)/(n+1))
  lnqq.emp <- log(X)
  
  # plots if TRUE
  .plotfun(lnqq.the, lnqq.emp, type="p", xlab="Quantiles of Standard Normal", ylab="log(X)", 
          main=main, plot=plot, add=FALSE, ...)
  abline(lm(lnqq.emp~lnqq.the))
  # output list with theoretical quantiles lnqq.the and empirical quantiles lnqq.emp
  .output(list(lnqq.the=lnqq.the, lnqq.emp=lnqq.emp), plot=plot, add=FALSE)
  
}

# Derivative plot of log-normal QQ-plot
LognormalQQ_der <- function(data, k = FALSE, plot = TRUE, main = "Derivative plot of log-normal QQ-plot", ...) {
  
  # Check input arguments
  .checkInput(data)
  
  X <- as.numeric(sort(data))
  n <- length(X)
  K <- 1:(n-1)
  
  if (k) {
    xval <- K
    xlab <- "k"
  } else {
    xval <- log(X[n-K])
    xlab <- bquote(log(X["n-k,n"]))
  }
  
  
  H <- Hill(X, plot=FALSE)$gamma
  N <- (n+1)/(K+1) * dnorm(qnorm(1-(K+1)/(n+1))) - qnorm(1-(K+1)/(n+1))
  yval <- H/N
  
  # plots if TRUE
  .plotfun(xval, yval, type="p", xlab=xlab, ylab="Derivative", 
          main=main, plot=plot, add=FALSE, ...)
  
  # output list with theoretical quantiles pqq.the and empirical quantiles pqq.emp
  .output(list(xval=xval, yval=yval), plot=plot, add=FALSE)
  
}

##############################################################################################

# Weibull QQ-plot
WeibullQQ <- function(data, plot = TRUE, main = "Weibull QQ-plot", ...) {
  
  # Check input arguments
  .checkInput(data)
  
  X <- as.numeric(sort(data))
  n <- length(X)
  
  # calculate theoretical and empirical quantiles
  wqq.the <- log(-log(1-(1:n)/(n+1)))
  wqq.emp <- log(X)
  
  # plots if TRUE
  .plotfun(wqq.the, wqq.emp, type="p", xlab="Quantiles of Standard Weibull", ylab="log(X)", 
           main=main, plot=plot, add=FALSE, ...)
  abline(lm(wqq.emp~wqq.the))
  # output list with theoretical quantiles lnqq.the and empirical quantiles lnqq.emp
  .output(list(wqq.the=wqq.the, wqq.emp=wqq.emp), plot=plot, add=FALSE)
  
}

# Derivative plot of Weibull QQ-plot
WeibullQQ_der <- function(data, k = FALSE, plot = TRUE, main = "Derivative plot of Weibull QQ-plot", ...) {
  
  # Check input arguments
  .checkInput(data)
  
  X <- as.numeric(sort(data))
  n <- length(X)
  K <- 1:(n-1)
  
  if (k) {
    xval <- K
    xlab <- "k"
  } else {
    xval <- log(X[n-K])
    xlab <- bquote(log(X["n-k,n"]))
  }
  
  
  H <- Hill(X, plot=FALSE)
  j <- 1:(n-1)
  yval <- H$gamma / (1/H$k * cumsum(log(log((n+1)/j))) - log(log((n+1)/(H$k+1))))

  # plots if TRUE
  .plotfun(xval, yval, type="p", xlab=xlab, ylab="Derivative", 
           main=main, plot=plot, add=FALSE, ...)
  
  # output list with theoretical quantiles pqq.the and empirical quantiles pqq.emp
  .output(list(xval=xval, yval=yval), plot=plot, add=FALSE)
  
}

##############################################################################################

# Computes the quantile scores to the generalised quantile plot
# (Section 5.2.3 in Beirlant et al. (2004)) for a numeric vector of observations (data)
#
# If plot=TRUE then the scores are plotted on a generalised
# quantile plot

genQQ <- function(data, gamma, plot = TRUE, main = "Generalised QQ-plot", ...) {
  
  # Check input arguments
  .checkInput(data, gamma)
  
  
  X <- as.numeric(sort(data))
  n <- length(X)
  gqq.the <- numeric(n)
  gqq.emp <- numeric(n)
  
  # calculate theoretical and empirical quantiles
  K <- 1:(n-1)
  gqq.the[K] <- log((n+1) / (K+1))
  gqq.emp[K] <- log(X[n-K] * gamma[K])
  
  
  # plots if TRUE
  .plotfun(gqq.the[K], gqq.emp[K], type="p", xlab="Quantiles of Standard Exponential", ylab="log(UH)", 
          main=main, plot=plot, add=FALSE, ...)
  abline(lm(gqq.emp~gqq.the))
  # output list with theoretical quantiles gqq.the
  # and empirical quantiles gqq.emp
  .output(list(gqq.the=gqq.the, gqq.emp=gqq.emp), plot=plot, add=FALSE)
  
}

# Add extra name for compatibility with S-Plus code
generalizedQQ <- genQQ

