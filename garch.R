#' Simulate Data from a GARCH(1,1) Model
#'
#' @param N total number of observations.
#' @param nsim number of time series to generate.
#' @param omega,alpha,beta GARCH(1,1) parameters.   
#' @param eps0,sig0 Initial values, included as first observation.
#' @param mu vector of means of the stocks
#' @return Return
garch.sim <- function(N, omega, alpha, beta, mu, R, sig20, eps0, y0) {
  library(MASS)
  #set.seed(100)
  if(N==1){
    z <- as.matrix(mvrnorm(N, rep(0, length(mu)), R))
  } else{
    z <- t(mvrnorm(N, rep(0, length(mu)), R))
  }
  # memory allocation
  eps <- matrix(NA, length(mu), N+1)
  sig2 <- matrix(NA, length(mu), N+1)
  eps[,1] <- eps0
  sig2[,1] <- sig20
  for(ii in (1:N)+1) {
    sig2[,ii] <- omega + alpha * eps[,ii-1]^2 + beta * sig2[,ii-1]
    eps[,ii] <- sqrt(sig2[,ii]) *z[,ii-1]
  }
  # remove starting values
  sig2 <- sig2[,-1,drop=FALSE]
  eps <- eps[,-1,drop=FALSE]
  xt <- mu + eps
  logp<- log(y0)+apply(xt, 1, cumsum)
  p <- exp(logp)
  p
}

#' Fit a Garch model for daily stock prices
#'
#' @param yt prices of stocks
#' @return MLE and zt (innovations).
garch.fit <- function(yt) {
  xt <- diff(log(yt))[-1] # log returns
  n <- length(xt)
  # convert to numeric format
  xt <- as.numeric(xt)
  muh <- mean(xt)
  
  #intial parameters
  eps <- xt - muh
  theta0 <- garch.clo(eps)
  eta0 <- theta0["alpha"]/theta0["omega"]
  beta0 <- theta0["beta"]
  eps20 <- mean(eps*eps)
  sig20 <- 1
  param0 <- c(eta0, beta0)
  gpfit <- optim(par = log(param0),
                 fn = function(param) {
                   -garch.profll(eta = exp(param[1]),
                                 beta = exp(param[2]),
                                 eps = eps, eps20 = eps20,
                                 sig20 = sig20)
                 })
  if(gpfit$convergence != 0) {
    warning("optim did not converge.")
  }
  
  # convert back to original GARCH parameters
  theta.mle <- garch.profmle(eta = exp(gpfit$par[1]),
                             beta = exp(gpfit$par[2]),
                             eps = eps, eps20 = eps20, sig20 = sig20)
  
  sig2 <- garch.sig2(theta.mle['omega'], theta.mle['alpha'], theta.mle['beta'], eps, eps20, sig20)
  zt <- eps/sig2
  
  c(list(theta = theta.mle,zt = zt, mu = muh, epsN = eps[n], sig2N =sig2[n]))
}

#' Fit a Garch model for weekly stock prices
#'
#' @param yt prices of stocks
#' @return MLE and zt (innovations).
garch.fitw <- function(yt) {
  xt <- diff(log(yt))[-1] # log returns
  n <- length(xt)
  # convert to numeric format
  xt <- as.numeric(xt)
  muh <- mean(xt)
  
  #intial parameters
  eps <- xt - muh
  eps20 <- mean(eps*eps)
  sig20 <- 1
  beta0 <- 0.7
  opfit <- optimize(f=garch.profll, c(0,1000),
                    beta= beta0, eps=eps, eps20 =eps20, 
                    sig20=sig20, maximum = TRUE)
  
  eta0 <- opfit$maximum
  param0 <- c(eta0, beta0)
  
  gpfit <- optim(par = log(param0),
                 fn = function(param) {
                   -garch.profll(eta = exp(param[1]),
                                 beta = exp(param[2]),
                                 eps = eps, eps20 = eps20,
                                 sig20 = sig20)
                 })
  if(gpfit$convergence != 0) {
    warning("optim did not converge.")
  }
  
  # convert back to original GARCH parameters
  theta.mle <- garch.profmle(eta = exp(gpfit$par[1]),
                             beta = exp(gpfit$par[2]),
                             eps = eps, eps20 = eps20, sig20 = sig20)
  
  sig2 <- garch.sig2(theta.mle['omega'], theta.mle['alpha'], theta.mle['beta'], eps, eps20, sig20)
  zt <- eps/sig2
  
  c(list(theta = theta.mle,zt = zt, mu = muh, epsN = eps[n], sig2N =sig2[n]))
}


#' Recover full parameter estimates from profile likelihood estimates
#' Credit to Prof Lysy
#' @param eta profile parameter equal to 
#' @param beta GARCH parameter.
#' @param eps20,sig20 initial standardized volatility
#' @return The named vector 
garch.profmle <- function(eta, beta, eps, eps20, sig20) {
  sig2 <- garch.sig2(1, eta, beta, eps, eps20, sig20) # tsig^2
  omegah <- mean(eps^2/sig2) # omega.hat(eta, beta)
  # transform back to original scale
  c(omega = as.numeric(omegah), alpha = as.numeric(eta*omegah),
    beta = as.numeric(beta))
}

#' Garch Profile Loglikelihood 
#' Credit to Prof Lysy
#' eps_t = omega * sig_t * z_t
#' z_t ~iid N(0,R)
#' sig_t^2 = 1 + eta * eps^2_t-1 + beta * sig^2_t-1,
#'
#' @param eta corresponds to alpha/beta from the original Garch model
#' @param beta parameter of original Garch model
#' @param eps vector of Garch observations 
#' @param eps2, sig20 initial value 
#' @export
garch.profll <- function(eta, beta, eps, eps20, sig20, debug = FALSE) {
  n <- length(eps)
  sig2 <- garch.sig2(1, eta, beta, eps, eps20, sig20) #Profiled out alpha, omega
  if(any(sig2 < 0)) return(-Inf) #Make sure all sigt2 is positive
  omegah <- mean(eps^2/sig2) 
  -.5 * (n + n*log(omegah) + sum(log(sig2)))
}


#' Computes Sigma^2 (volatilities) at time t from GARCH observations
#' Credit to Prof Lysy
#' @param omega, alpha, beta Garch parameters.
#' @param eps vector of Garch observations 
#' @param eps20, sig20 initial values.  
#' @return sig^2
garch.sig2 <- function(omega, alpha, beta, eps, eps20, sig20) {
  n <- length(eps)
  eps2 <- eps^2
  mu <- mean(eps2)
  if(missing(eps20)) eps20 <- mu
  if(missing(sig20)) sig20 <- mu
  fsig2 <- omega + alpha * c(eps20, eps2[-n]) #First part of the sigma^2 equation
  sig2 <-as.numeric(filter(fsig2, beta, "recursive", init = sig20))
  sig2
}

#' GARCH Loglikelihood function
#' Credit to Prof Lysy
#' @param omega,alpha,beta GARCH parameters (scalars).
#' @param eps vector of GARCH observations 
#' @param eps20,sig20 initial values.  For default values see
#' @return The GARCH log-likelihood.
garch.loglik <- function(omega, alpha, beta, eps, eps20, sig20) {
  sig2 <- garch.sig2(omega, alpha, beta, eps, eps20, sig20) # sigma^2
  if(any(sig2 < 0)) return(-Inf) # prevent parameters giving negative variances
  sum(dnorm(eps, mean = 0, sd = sqrt(sig2), log = TRUE))
}


#' Closed-form estimator by Kristensen & Linton (using notation from the paper)
#'
#' @param eps Garch observation parameter
#' @return vector of estimated Garch parameters (alpha, beta, omega).
garch.clo <- function(eps) {
  n <- length(eps)
  xt <- eps*eps
  sigma2 <- mean(xt)
  rho <- cov(xt[1:(n-1)],xt[2:n])
  phi <- cov(xt[1:(n-2)], xt[3:n])/rho
  b <-  (phi^2 + 1 - 2*rho*phi)/(phi - rho)
  theta <- (-b + sqrt(b^2 -4))/2
  alpha <- theta + phi
  if (alpha <0 ){
    alpha <- abs(alpha)
  }
  beta <- -theta
  omega <- sigma2*(alpha+beta)
  c(omega = omega, alpha = alpha, beta = beta)
}
