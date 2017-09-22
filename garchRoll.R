#Transform sample to normal looking histogram
ztrans <- function(x) {
  n <- length(x)
  z <- rank(x)/(n+1) # approximately have z[i] = Pr(x <= x[i])
  qnorm(z) # corresponding z-scores
}

#Get rank of sample
zrank <- function(x){
  n <- length(x)
  rank <- rank(x)/(n)
}

#' Simulate Data from a GARCH(1,1) Model
#'
#' @param N total number of observations.
#' @param nsim number of time series to generate.
#' @param omega,alpha,beta GARCH(1,1) parameters.   
#' @param eps0,sig0 Initial values, included as first observation.
#' @param mu vector of means of the stocks
#' @return Return
garch.rsim <- function(N, omega, alpha, beta, mu, R, sig20, eps0, y0, emp, Zt) {
  library(MASS)
  #set.seed(100)
  if(N==1){
    z <- as.matrix(mvrnorm(N, rep(0, length(mu)*5), R))
  } else{
    z <- t(mvrnorm(N, rep(0, length(mu)*5), R)) #5-day innovations
  }
  if(emp){
    pz <- pnorm(z)
    empzrank <- apply(Zt, 2, zrank)
    stockrank <- matrix(NA,length(mu)*5,N)
    empz <- matrix(NA, length(mu)*5,N)
    for(i in 1:(length(mu)*5)){
      stockrank[i,] <-vapply(pz[i,],function(x) which.min(abs(empzrank[,i]-x)),1)
      empz[i,] <- Zt[stockrank[i,],i]
      z <- empz
    }
  }
  # memory allocation for each the 5-day innovations
  eps1 <- matrix(NA, length(mu), N+1)
  eps2 <- matrix(NA, length(mu), N+1)
  eps3 <- matrix(NA, length(mu), N+1)
  eps4 <- matrix(NA, length(mu), N+1)
  eps5 <- matrix(NA, length(mu), N+1)
  sig21 <- matrix(NA, length(mu), N+1)
  sig22 <- matrix(NA, length(mu), N+1)
  sig23 <- matrix(NA, length(mu), N+1)
  sig24 <- matrix(NA, length(mu), N+1)
  sig25 <- matrix(NA, length(mu), N+1)
  #z <- cbind(eps0/sig0, z) # first innovation determined by initial values
  eps1[,1] <- eps0[,1] # 1st of 5-day epsilon (not epsilon of stock 1)
  eps2[,1] <- eps0[,2]
  eps3[,1] <- eps0[,3]
  eps4[,1] <- eps0[,4]
  eps5[,1] <- eps0[,5]
  sig21[,1] <- sig20[,1]
  sig22[,1] <- sig20[,2]
  sig23[,1] <- sig20[,3]
  sig24[,1] <- sig20[,4]
  sig25[,1] <- sig20[,5]
  
  for(ii in (1:N)+1) {
    #Simulate 5-day eps and sig2, order is z^(W,1)_1, z^(W,2)_1, ... 
    sig21[,ii] <- omega + alpha * eps1[,ii-1]^2 + beta * sig21[,ii-1]
    eps1[,ii] <- sqrt(sig21[,ii]) *z[seq(1, length(mu), 5),ii-1] 
    sig22[,ii] <- omega + alpha * eps2[,ii-1]^2 + beta * sig22[,ii-1]
    eps2[,ii] <- sqrt(sig22[,ii]) *z[seq(2, length(mu), 5),ii-1]
    sig23[,ii] <- omega + alpha * eps3[,ii-1]^2 + beta * sig23[,ii-1]
    eps3[,ii] <- sqrt(sig23[,ii]) *z[seq(3, length(mu), 5),ii-1]
    sig24[,ii] <- omega + alpha * eps4[,ii-1]^2 + beta * sig24[,ii-1]
    eps4[,ii] <- sqrt(sig24[,ii]) *z[seq(4, length(mu), 5),ii-1]
    sig25[,ii] <- omega + alpha * eps5[,ii-1]^2 + beta * sig25[,ii-1]
    eps5[,ii] <- sqrt(sig25[,ii]) *z[seq(5, length(mu), 5),ii-1]
  }
  # remove starting values
  sig21 <- sig21[,-1,drop=FALSE]
  eps1 <- eps1[,-1,drop=FALSE]
  sig22 <- sig22[,-1,drop=FALSE]
  eps2 <- eps2[,-1,drop=FALSE]
  sig23 <- sig23[,-1,drop=FALSE]
  eps3 <- eps3[,-1,drop=FALSE]
  sig24 <- sig24[,-1,drop=FALSE]
  eps4 <- eps4[,-1,drop=FALSE]
  sig25 <- sig25[,-1,drop=FALSE]
  eps5 <- eps5[,-1,drop=FALSE]
  mu <- c(mu)
  xt1 <- mu + eps1
  xt2 <- mu + eps2
  xt3 <- mu + eps3
  xt4 <- mu + eps4
  xt5 <- mu + eps5
  xt <- (xt1 + xt2 + xt3 + xt4 + xt5)/5
  logp<- log(y0)+t(apply(xt, 1, cumsum))
  p <- exp(logp)
  p
}

#' Fit a Garch model for rolling stock prices
#'
#' @param yt prices of stocks
#' @return MLE and zt (innovations).
garch.fitr <- function(yt1, yt2, yt3, yt4, yt5) {
  xt1 <- diff(log(yt1))[-1]
  xt2 <- diff(log(yt2))[-1]
  xt3 <- diff(log(yt3))[-1]
  xt4 <- diff(log(yt4))[-1]
  xt5 <- diff(log(yt5))[-1]# log returns
  n <- length(xt1)
  # convert to numeric format
  xt1 <- as.numeric(xt1)
  xt2 <- as.numeric(xt2)
  xt3 <- as.numeric(xt3)
  xt4 <- as.numeric(xt4)
  xt5 <- as.numeric(xt5)
  muh1 <- mean(xt1)
  muh2 <- mean(xt2)
  muh3 <- mean(xt3)
  muh4 <- mean(xt4)
  muh5 <- mean(xt5)
  muh <- mean(muh1, muh2, muh3, muh4, muh5)
  
  #intial parameters
  eps1 <- xt1 - muh
  eps2 <- xt2- muh
  eps3 <- xt3 - muh
  eps4 <- xt4- muh
  eps5 <- xt5 - muh
  
  theta0 <- garch.clo(eps1)+ garch.clo(eps2)+ garch.clo(eps3)+ garch.clo(eps4)+ garch.clo(eps5)
  theta0 <- theta0/5
  #eta0 <- abs(theta0["alpha"]/theta0["omega"])
  eta0 <- 300
  beta0 <- 0.65
  
  eps20 <- rowMeans(rbind(eps1*eps1, eps2*eps2, eps3*eps3, eps4*eps4, eps5*eps5), na.rm=TRUE)
  sig20 <- 1
  param0 <- c(eta0, beta0)
  
  
  gpfit <- optim(par = log(param0),
                 fn = function(param) {
                   -garch.rprofll(eta = exp(param[1]),
                                 beta = exp(param[2]),
                                 eps1 = eps1, eps2=eps2, eps3=eps3, eps4=eps4,
                                 eps5 = eps5,eps20 = eps20,
                                 sig20 = sig20)
                 })
  if(gpfit$convergence != 0) {
    warning("optim did not converge.")
  }
  
  # convert back to original GARCH parameters
  theta.mle <- garch.rprofmle(eta = exp(gpfit$par[1]),
                             beta = exp(gpfit$par[2]),
                             eps = eps1,eps2=eps2, eps3=eps3, eps4=eps4,
                             eps5=eps5, eps20 = eps20, sig20 = sig20)
  
  sig21 <- garch.sig2(theta.mle['omega'], theta.mle['alpha'], 
                     theta.mle['beta'], eps1, eps20[1], sig20)
  sig22 <- garch.sig2(theta.mle['omega'], theta.mle['alpha'], 
                      theta.mle['beta'], eps2, eps20[2], sig20)
  sig23 <- garch.sig2(theta.mle['omega'], theta.mle['alpha'], 
                      theta.mle['beta'], eps3, eps20[3], sig20)
  sig24 <- garch.sig2(theta.mle['omega'], theta.mle['alpha'], 
                      theta.mle['beta'], eps4, eps20[4], sig20)
  sig25 <- garch.sig2(theta.mle['omega'], theta.mle['alpha'], 
                      theta.mle['beta'], eps5, eps20[5], sig20)
  
  zt1 <- eps1/sig21
  zt2 <- eps2/sig22
  zt3 <- eps3/sig23
  zt4 <- eps4/sig24
  zt5 <- eps5/sig25
  
  eps <- c(eps1[n], eps2[n], eps3[n], eps4[n], eps5[n])
  sig2 <- c(sig21[n], sig22[n], sig23[n], sig24[n], sig25[n])
  c(list(theta = theta.mle,zt1 = zt1, zt2 =zt2, zt3=zt3, zt4=zt4, zt5=zt5,
         mu = muh, epsN = eps, sig2N =sig2))
}

#' Recover full parameter estimates from profile likelihood estimates
#' Credit to Prof Lysy
#' @param eta profile parameter equal to 
#' @param beta GARCH parameter.
#' @param eps20,sig20 initial standardized volatility
#' @return The named vector 
#' 
garch.rprofmle <- function(eta, beta, eps1, eps2, eps3, eps4, eps5, eps20, sig20) {
  sig21 <- garch.sig2(1, eta, beta, eps1, eps20[1], sig20) # tsig^2
  sig22 <- garch.sig2(1, eta, beta, eps2, eps20[2], sig20)
  sig23 <- garch.sig2(1, eta, beta, eps3, eps20[3], sig20)
  sig24 <- garch.sig2(1, eta, beta, eps4, eps20[4], sig20)
  sig25 <- garch.sig2(1, eta, beta, eps5, eps20[5], sig20)
  
  omegah1 <- mean(eps1^2/sig21)
  omegah2 <- mean(eps2^2/sig22)
  omegah3 <- mean(eps3^2/sig23)
  omegah4 <- mean(eps4^2/sig24)
  omegah5 <- mean(eps5^2/sig25)
  omegah <- mean(c(omegah1,omegah2,omegah3, omegah4, omegah5))
  
  # omega.hat(eta, beta)
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
garch.rprofll <- function(eta, beta, eps1, eps2, eps3, eps4, eps5, eps20, sig20, debug = FALSE) {
  n <- length(eps1)
  sig21 <- garch.sig2(1, eta, beta, eps1, eps20[1], sig20) #Profiled out alpha, omega
  sig22 <- garch.sig2(1, eta, beta, eps2, eps20[2], sig20) 
  sig23 <- garch.sig2(1, eta, beta, eps3, eps20[3], sig20) 
  sig24 <- garch.sig2(1, eta, beta, eps4, eps20[4], sig20) 
  sig25 <- garch.sig2(1, eta, beta, eps5, eps20[5], sig20) 
  
  if(any(sig21 < 0)) return(-Inf) #Make sure all sigt2 is positive
  if(any(sig22 < 0)) return(-Inf)
  if(any(sig23 < 0)) return(-Inf)
  if(any(sig24 < 0)) return(-Inf)
  if(any(sig25 < 0)) return(-Inf)
  
  omegah1 <- mean(eps1^2/sig21)
  omegah2 <- mean(eps2^2/sig22)
  omegah3 <- mean(eps3^2/sig23)
  omegah4 <- mean(eps4^2/sig24)
  omegah5 <- mean(eps5^2/sig25)
  d1 <- -.5 * (n + n*log(omegah1) + sum(log(sig21)))
  d2 <- -.5 * (n + n*log(omegah2) + sum(log(sig22)))
  d3 <- -.5 * (n + n*log(omegah3) + sum(log(sig23)))
  d4 <- -.5 * (n + n*log(omegah4) + sum(log(sig24)))
  d5 <- -.5 * (n + n*log(omegah5) + sum(log(sig25)))
  d1+d2+d3+d4+d5
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

#' GARCH Loglikelihood function for rolling window
#'
#' @param omega,alpha,beta GARCH parameters (scalars).
#' @param eps1-5 vector of GARCH observations
#' @param eps20,sig20 initial values.
#' @return The GARCH log-likelihood.
garch.rloglik <- function(omega, alpha, beta, eps1, eps2, eps3, eps4, eps5, eps20, sig20) {
  sig21 <- garch.sig2(omega, alpha, beta, eps1, eps20[1], sig20) # sigma^2
  sig22 <- garch.sig2(omega, alpha, beta, eps2, eps20[2], sig20)
  sig23 <- garch.sig2(omega, alpha, beta, eps3, eps20[3], sig20)
  sig24 <- garch.sig2(omega, alpha, beta, eps4, eps20[4], sig20)
  sig25 <- garch.sig2(omega, alpha, beta, eps5, eps20[5], sig20)
  
  if(any(sig21 < 0)) return(-Inf) # prevent parameters giving negative variances
  if(any(sig22 < 0)) return(-Inf)
  if(any(sig23 < 0)) return(-Inf)
  if(any(sig24 < 0)) return(-Inf)
  if(any(sig25 < 0)) return(-Inf)
  
  d1 <- sum(dnorm(eps1, mean = 0, sd = sqrt(sig21), log = TRUE))
  d2 <- sum(dnorm(eps2, mean = 0, sd = sqrt(sig22), log = TRUE))
  d3 <- sum(dnorm(eps3, mean = 0, sd = sqrt(sig23), log = TRUE))
  d4 <- sum(dnorm(eps4, mean = 0, sd = sqrt(sig24), log = TRUE))
  d5 <- sum(dnorm(eps5, mean = 0, sd = sqrt(sig25), log = TRUE))
  d1+d2+d3+d4+d5
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
  beta <- -theta
  omega <- sigma2*(alpha+beta)
  c(omega = omega, alpha = alpha, beta = beta)
}
