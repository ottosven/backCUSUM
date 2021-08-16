## ####################################################################
## ####################################################################
## Supplement for
## "Backward CUSUM for Testing and Monitoring Structural Change"
## by Sven Otto and Jörg Breitung.
## This R-script allows to reproduce the simulated data for the upper panels of Figure 3.
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
## Simulated critical values:
CRIT <- c(0.9376, 0.9376, 1.1875)
## Simulated critical values can be replicated using replicateFigure3_data_init_upper.R
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
  r <- (1:T)/T
  h <- numeric(T)
  J <- r>=rstar
  h[J] <- c*rstar*(log(r[J])-log(rstar))
  boundary <- 1+2*r
  sim.dist <- function(x){
    W <- BrownianMotion(T)
    Q <- max(abs(W + h)/boundary)
    BQ <- max(abs(W + h[T]-rev(h))/boundary)
    SBQ <- max(sapply(2:T,getSBQ, W=W, c=c, rstar=rstar, h=h))
    c(Q,BQ,SBQ)
  }
  DETECTORS<-sapply(1:MC,sim.dist)
  results <- numeric(dim(DETECTORS)[1])
  for(i in 1:dim(DETECTORS)[1]) ( results[i] <- length(which(DETECTORS[i,] > CRIT[i]))/MC )
  results
}
##
C <- 0:30
Qnames <- c('Q', 'BQ', 'SBQ')
## First panel data
rstar <- 0.1
powcur.ret.r01 <- data.frame(mcmapply(SIM, C, mc.cores = CORE, rstar = rstar))
colnames(powcur.ret.r01) <- C
rownames(powcur.ret.r01) <- Qnames
powcur.ret.r01
## Second panel data
rstar <- 0.3
powcur.ret.r03 <- data.frame(mcmapply(SIM, C, mc.cores = CORE, rstar = rstar))
colnames(powcur.ret.r03) <- C
rownames(powcur.ret.r03) <- Qnames
powcur.ret.r03
## Third panel data
rstar <- 0.5
powcur.ret.r05 <- data.frame(mcmapply(SIM, C, mc.cores = CORE, rstar = rstar))
colnames(powcur.ret.r05) <- C
rownames(powcur.ret.r05) <- Qnames
powcur.ret.r05
## Fourth panel data
rstar <- 0.7
powcur.ret.r07 <- data.frame(mcmapply(SIM, C, mc.cores = CORE, rstar = rstar))
colnames(powcur.ret.r07) <- C
rownames(powcur.ret.r07) <- Qnames
powcur.ret.r07
## Fifth panel data
rstar <- 0.9
powcur.ret.r09 <- data.frame(mcmapply(SIM, C, mc.cores = CORE, rstar = rstar))
colnames(powcur.ret.r09) <- C
rownames(powcur.ret.r09) <- Qnames
powcur.ret.r09
##
R <- 1:39/40
## Sixth panel data
c <- 10
powcur.ret.c10 <- data.frame(mcmapply(SIM, R, mc.cores = CORE, c = c))
colnames(powcur.ret.c10) <- R
rownames(powcur.ret.c10) <- Qnames
powcur.ret.c10
## ##################################
## ##################################
end<-Sys.time()
end-start
