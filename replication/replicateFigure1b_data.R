## ####################################################################
## ####################################################################
## Supplement for
## "Backward CUSUM for Testing and Monitoring Structural Change"
## by Sven Otto and Jörg Breitung.
## This R-script allows to reproduce the simulation for Figure 4.
## ####################################################################
## ####################################################################
rm(list=ls())
start<-Sys.time()
library("parallel")
## ##################################
## Reproducible random number
## generator for parallelization
## ##################################
RNGkind("L'Ecuyer-CMRG")
set.seed(42)
CORE <- detectCores(all.tests = FALSE, logical = TRUE)
## ##################################
## Simulation:
## ##################################
## Monte Carlo repetitions:
MC <- 10
## ##################################
## For replication please set the number of repetitions to
## MC <- 100000
## Warning! This process takes a very long time
## ##################################
T <- 1000
m <- 4
## Simulated critical values:
CRIT <- c(0.9440, 0.7957, 1.3190)
## Simulated critical values can be replicated using replicateFigure4_data_init.R
##
BrownianMotion <- function(T)  ( cumsum(rnorm(T,0,sqrt(1/T))) )
##
getSBQ <- function(rT, W, c, rstar, h){
  sT <- 1:(rT-1)
  bound <- 1+2*(rT-sT)/T
  max(abs(W[rT] - W[sT] + h[rT] - h[sT])/bound)
}
##
SIM <- function(c = 0, rstar = 1){
  r <- (m-1)*(1:T)/T
  h <- numeric(T)
  J <- r >= rstar - 1
  h[J] <- c*(rstar-1)*(log(r[J])-log(rstar-1))
  boundary <- 1+2*r
  boundaryChu <- sqrt((r+1)*(log((r+1)/0.05^2)))
  sim.dist <- function(x){
    W <- sqrt(m-1)*BrownianMotion(T)
    allQ <- abs(W + h)/boundary
    allQChu <- abs(W + h)/boundaryChu
    allSBQ <- c(0,sapply(2:T,getSBQ, W=W, c=c, rstar=rstar, h=h))
    Q <- which(allQ > CRIT[1])[1]
    QChu <- which(allQChu > CRIT[2])[1]
    SBQ <- which(allSBQ > CRIT[3])[1]
    c(SBQ,Q,QChu)
  }
  DETECTORS<-sapply(1:MC,sim.dist)
  meandelay <- numeric(dim(DETECTORS)[1])
  for(i in 1:dim(DETECTORS)[1]){
    meandelay[i] <- (mean(DETECTORS[i,which(DETECTORS[i,] >= (rstar-1)/(m-1)*T)]) - ceiling((rstar-1)/(m-1)*T))/T
  }
  results <- c(meandelay)
  results
}
##
C <- 0:30
Qnames <- c('SBQ-delay', 'Q-delay', 'QChu-delay')
## First panel data
rstar <- 1.5
delcur.mon.r15 <- data.frame(mcmapply(SIM, C, mc.cores = CORE, rstar = rstar))
colnames(delcur.mon.r15) <- C
rownames(delcur.mon.r15) <- Qnames
delcur.mon.r15
## Second panel data
rstar <- 3
delcur.mon.r3 <- data.frame(mcmapply(SIM, C, mc.cores = CORE, rstar = rstar))
colnames(delcur.mon.r3) <- C
rownames(delcur.mon.r3) <- Qnames
delcur.mon.r3
##
R <- 3*(1:39/40) + 1
## Third panel data
c <- 20
delcur.mon.c20 <- data.frame(mcmapply(SIM, R, mc.cores = CORE, c = c))
colnames(delcur.mon.c20) <- R
rownames(delcur.mon.c20) <- Qnames
delcur.mon.c20
## ##################################
## ##################################
end<-Sys.time()
end-start
