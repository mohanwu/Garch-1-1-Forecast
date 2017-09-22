#MLE Checks

source("mle-check.R")
source("garch.R")
source("garchRoll.R")

require(quantmod)

getSymbols(Symbols = c("LEA"),
           src = "yahoo", from = "1990-01-01")
dates <- c("2012-01-04", "2015-12-31")
yt <- LEA[paste0(dates, collapse = "/")]$LEA.Close


#Get eps from price
geteps <- function(yt){
  xt <- diff(log(yt))[-1] # log returns
  n <- length(xt)
  # convert to numeric format
  xt <- as.numeric(xt)
  muh <- mean(xt)
  #intial parameters
  eps <- xt - muh
  eps20 <- mean(eps^2)
  list(eps20=eps20, eps = eps)
}

#Rolling Geteps
rgeteps <-function(yt1, yt2, yt3, yt4, yt5){
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
  eps20 <- rowMeans(rbind(eps1*eps1, eps2*eps2, 
                          eps3*eps3, eps4*eps4, eps5*eps5), na.rm=TRUE)
  list(eps1=eps1, eps2=eps2, eps3=eps3, eps4 =eps4, eps5=eps5, eps20=eps20)
}

#Garch Loglikihood function for weekly/daily data
loglik <- function(theta) {
  garch.loglik(omega = theta[1], alpha = theta[2], beta = theta[3],
               eps = eps, eps20 = eps20,
               sig20 = theta[1])
}

#Garch Loglikihood function for rolling data
rloglik <- function(theta) {
  garch.rloglik(omega = theta[1], alpha = theta[2], beta = theta[3],
                eps1 = eps1, eps2=eps2, eps3=eps3, eps4=eps4, eps5=eps5,
                eps20,sig20 = theta[1])
}

#Daily data
eps <- geteps(yt)$eps
eps20 <- geteps(yt)$eps20
theta.mle <- garch.fit(yt)$theta

tnames <- expression(omega, alpha, beta)
#plot check
mle.check(loglik = loglik, theta.mle = theta.mle,
          refit = TRUE, theta.names = tnames) # pass

#Weekly Data
wt <- to.weekly(yt)$yt.Close
eps <- geteps(wt)$eps
eps20 <- geteps(wt)$eps20
theta.mle <- garch.fitw(wt)$theta

tnames <- expression(omega, alpha, beta)
#plot check
mle.check(loglik = loglik, theta.mle = theta.mle,
          refit = TRUE, theta.names = tnames) # pass

#Rolling data 

n <- length(yt)
rt1 <- yt[seq(1, n, 5)]
rt2 <- yt[seq(2, n, 5)]
rt3 <- yt[seq(3, n, 5)]
rt4 <- yt[seq(4, n, 5)]
rt5 <- yt[seq(5, n, 5)]

eps <- rgeteps(rt1, rt2, rt3, rt4, rt5)
eps1 <- eps$eps1
eps2 <- eps$eps2
eps3 <- eps$eps3
eps4 <- eps$eps4
eps5 <- eps$eps5
eps20 <- eps$eps20

theta.mle <- garch.fitr(rt1, rt2, rt3, rt4, rt5)$theta

tnames <- expression(omega, alpha, beta)
#plot check
mle.check(loglik = rloglik, theta.mle = theta.mle,
          refit = TRUE, theta.names = tnames) # pass

