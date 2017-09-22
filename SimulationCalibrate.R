#Daily 5-day Simulation

require(quantmod)
source("garchtest.R")
source("garch.R")
source("garchRoll.R")
getSymbols(Symbols = c("BIIB","MAN","LEA","DHI","NHTC"),
           src = "yahoo", from = "1990-01-01")
dates <- c("2012-01-01", "2016-02-29")
yt1 <- BIIB[paste0(dates, collapse = "/")]$BIIB.Close
yt2 <- MAN[paste0(dates, collapse = "/")]$MAN.Close
yt3 <- LEA[paste0(dates, collapse = "/")]$LEA.Close
yt4 <- DHI[paste0(dates, collapse = "/")]$DHI.Close
yt5 <- NHTC[paste0(dates, collapse = "/")]$NHTC.Close

#Real Prices
dates2 <- c("2016-03-01", "2016-03-30")
pt1 <- BIIB[paste0(dates2, collapse = "/")]$BIIB.Close
pt2 <- MAN[paste0(dates2, collapse = "/")]$MAN.Close
pt3 <- LEA[paste0(dates2, collapse = "/")]$LEA.Close
pt4 <- DHI[paste0(dates2, collapse = "/")]$DHI.Close
pt5 <- NHTC[paste0(dates2, collapse = "/")]$NHTC.Close

set.seed(100)
rprice <- pt3[5]
terror <- 0
sim <- 50
for (i in 1:sim){
  dayforecast <- garch.dforecast(yt1, yt2, yt3, yt4, yt5, 'F', 5)[3]
  terror <- terror + sum(abs(dayforecast - rprice))
}
merror <- terror/sim

#Weekly Forecast
wt1 <- to.weekly(yt1)$yt1.Close
wt2 <- to.weekly(yt2)$yt2.Close
wt3 <- to.weekly(yt3)$yt3.Close
wt4 <- to.weekly(yt4)$yt4.Close
wt5 <- to.weekly(yt5)$yt5.Close

set.seed(100)
terrorw <-0
for (i in 1:sim){
  weekforecast <- garch.wforecast(wt1, wt2, wt3, wt4, wt5, 'F', 1)[3]
  terrorw <- terrorw + sum(abs(weekforecast - rprice))
}
merrorw <- terrorw/sim


#Rolling

set.seed(100)
terrorr <- 0

for (i in 1:sim){
  rollforecast <- garch.rforecast(yt1, yt2, yt3, yt4, yt5, 'F', 5)[3]
  terrorr <- terrorr + sum(abs(rollforecast - rprice))
}
merrorr <- terrorr/sim
