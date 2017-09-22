garch.dforecast <- function(yt1, yt2, yt3, yt4, yt5, R, N){
  daily1 <- garch.fit(yt1)
  daily2 <- garch.fit(yt2)
  daily3 <- garch.fit(yt3)
  daily4 <- garch.fit(yt4)
  daily5 <- garch.fit(yt5)
  
  theta.mle <- daily1$theta
  omega <- c(daily1$theta['omega'],daily2$theta['omega'], daily3$theta['omega'],
             daily4$theta['omega'], daily5$theta['omega'])
  alpha <- c(daily1$theta['alpha'], daily2$theta['alpha'], daily3$theta['alpha'],
             daily4$theta['alpha'], daily5$theta['alpha'])
  beta <- c(daily1$theta['beta'], daily2$theta['beta'], daily3$theta['beta'],
            daily4$theta['beta'], daily5$theta['beta'])
  
  mu <- c(daily1$mu, daily2$mu, daily3$mu, daily4$mu, daily5$mu)
  if (R == 'I'){
    R <- diag(length(mu)) #Assuming Z ~ MVN(0, I)
  } else{
    Zt <- cbind(daily1$zt, daily2$zt, daily3$zt, daily4$zt, daily5$zt)
    R <- cor(Zt)
  }
  sig20 <- c(daily1$sig2N, daily2$sig2N, daily3$sig2N, daily4$sig2N, daily5$sig2N)
  eps0 <- c(daily1$epsN, daily2$epsN, daily3$epsN, daily4$epsN, daily5$epsN)
  y <-  cbind(yt1, yt2, yt3, yt4, yt5)
  y0 <- as.numeric(tail(y,n=1))
  
  price <- garch.sim(N, omega, alpha, beta, mu, R, sig20, eps0, y0)
  price[,N]
}


garch.wforecast <- function(wt1, wt2, wt3, wt4, wt5, R, N){
  weekly1 <- garch.fitw(wt1)
  weekly2 <- garch.fitw(wt2)
  weekly3 <- garch.fitw(wt3)
  weekly4 <- garch.fitw(wt4)
  weekly5 <- garch.fitw(wt5)
 
  omega <- c(weekly1$theta['omega'],weekly2$theta['omega'], weekly3$theta['omega'],
             weekly4$theta['omega'], weekly5$theta['omega'])
  alpha <- c(weekly1$theta['alpha'], weekly2$theta['alpha'], weekly3$theta['alpha'],
             weekly4$theta['alpha'], weekly5$theta['alpha'])
  beta <- c(weekly1$theta['beta'], weekly2$theta['beta'], weekly3$theta['beta'],
            weekly4$theta['beta'], weekly5$theta['beta'])
  
  mu <- c(weekly1$mu, weekly2$mu, weekly3$mu, weekly4$mu, weekly5$mu)
  if (R == 'I'){
    R <- diag(length(mu)) #Assuming Z ~ MVN(0, I)
  } else{
    Zt <- cbind(weekly1$zt, weekly2$zt, weekly3$zt, weekly4$zt, weekly5$zt)
    R <- cor(Zt)
  }
  sig20 <- c(weekly1$sig2N, weekly2$sig2N, weekly3$sig2N, weekly4$sig2N, weekly5$sig2N)
  eps0 <- c(weekly1$epsN, weekly2$epsN, weekly3$epsN, weekly4$epsN, weekly5$epsN)
  y <-  cbind(yt1, yt2, yt3, yt4, yt5)
  y0 <- as.numeric(tail(y,n=1))
  
  price <- garch.sim(N, omega, alpha, beta, mu, R, sig20, eps0, y0)
  price
}

garch.rforecast <- function(yt1, yt2, yt3, yt4, yt5, R, N){
  n <- length(yt1)
  rt11 <- yt1[seq(1, n, 5)]
  rt12 <- yt1[seq(2, n, 5)]
  rt13 <- yt1[seq(3, n, 5)]
  rt14 <- yt1[seq(4, n, 5)]
  rt15 <- yt1[seq(5, n, 5)]
  
  rt21 <- yt2[seq(1, n, 5)]
  rt22 <- yt2[seq(2, n, 5)]
  rt23 <- yt2[seq(3, n, 5)]
  rt24 <- yt2[seq(4, n, 5)]
  rt25 <- yt2[seq(5, n, 5)]
  
  rt31 <- yt3[seq(1, n, 5)]
  rt32 <- yt3[seq(2, n, 5)]
  rt33 <- yt3[seq(3, n, 5)]
  rt34 <- yt3[seq(4, n, 5)]
  rt35 <- yt3[seq(5, n, 5)]
  
  rt41 <- yt4[seq(1, n, 5)]
  rt42 <- yt4[seq(2, n, 5)]
  rt43 <- yt4[seq(3, n, 5)]
  rt44 <- yt4[seq(4, n, 5)]
  rt45 <- yt4[seq(5, n, 5)]
  
  rt51 <- yt5[seq(1, n, 5)]
  rt52 <- yt5[seq(2, n, 5)]
  rt53 <- yt5[seq(3, n, 5)]
  rt54 <- yt5[seq(4, n, 5)]
  rt55 <- yt5[seq(5, n, 5)]
  
  roll1 <- garch.fitr(rt11, rt12, rt13, rt14, rt15)
  roll2 <- garch.fitr(rt21, rt22, rt23, rt24, rt25)
  roll3 <- garch.fitr(rt31, rt32, rt33, rt34, rt35)
  roll4 <- garch.fitr(rt41, rt42, rt43, rt44, rt45)
  roll5 <- garch.fitr(rt51, rt52, rt53, rt54, rt55)
  
  omega <- c(roll1$theta['omega'],roll2$theta['omega'], roll3$theta['omega'],
             roll4$theta['omega'], roll5$theta['omega'])
  alpha <- c(roll1$theta['alpha'], roll2$theta['alpha'], roll3$theta['alpha'],
             roll4$theta['alpha'], roll5$theta['alpha'])
  beta <- c(roll1$theta['beta'], roll2$theta['beta'], roll3$theta['beta'],
            roll4$theta['beta'], roll5$theta['beta'])
  
  mu <- cbind(roll1$mu, roll2$mu, roll3$mu, roll4$mu, roll5$mu)
  if (R == 'I'){
    R <- diag(length(mu)*5) #Assuming Z ~ MVN(0, I)
    emp <- FALSE
  } else if(R=='F'){
    Zt <- cbind(roll1$zt1, roll1$zt2, roll1$zt3, roll1$zt4, roll1$zt5,
                roll2$zt1, roll2$zt2, roll2$zt3, roll2$zt4, roll2$zt5,
                roll3$zt1, roll3$zt2, roll3$zt3, roll3$zt4, roll3$zt5,
                roll4$zt1, roll4$zt2, roll4$zt3, roll4$zt4, roll4$zt5,
                roll5$zt1, roll5$zt2, roll5$zt3, roll5$zt4, roll5$zt5)
    R <- cor(Zt)
    emp <-FALSE
  }else{
    Zt <- cbind(roll1$zt1, roll1$zt2, roll1$zt3, roll1$zt4, roll1$zt5,
                roll2$zt1, roll2$zt2, roll2$zt3, roll2$zt4, roll2$zt5,
                roll3$zt1, roll3$zt2, roll3$zt3, roll3$zt4, roll3$zt5,
                roll4$zt1, roll4$zt2, roll4$zt3, roll4$zt4, roll4$zt5,
                roll5$zt1, roll5$zt2, roll5$zt3, roll5$zt4, roll5$zt5)
    Ztrans <- apply(Zt, 2, ztrans)
    R <- cor(Ztrans)
    emp <- TRUE
  }
  sig20 <- rbind(roll1$sig2N, roll2$sig2N, roll3$sig2N, roll4$sig2N, roll5$sig2N)
  eps0 <- rbind(roll1$epsN, roll2$epsN, roll3$epsN, roll4$epsN, roll5$epsN)
  y <-  cbind(yt1, yt2, yt3, yt4, yt5)
  y0 <- as.numeric(tail(y,n=1))
  
  price <- garch.rsim(N, omega, alpha, beta, mu, R, sig20, 
                      eps0, y0, emp, Zt)
  if (N==1){
    price
  }else{
    price[,N]
  }
}