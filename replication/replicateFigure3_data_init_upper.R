## ####################################################################
## ####################################################################
## Supplement for
## "Backward CUSUM for Testing and Monitoring Structural Change"
## by Sven Otto and Jörg Breitung.
## This R-script allows to reproduce the critical values for the upper panels
## CRIT <- c(0.9376, 0.9376, 1.1875) for the simulation for Figure 3.
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
## For replication of the data
## please set the number of
## Monte Carlo repetitions to
## MC <- 10000000
## ##################################
## Monte Carlo repetitions:
MC <- 1000
##
T <- 1000
##
BrownianMotion <- function(T)  ( cumsum(rnorm(T,0,sqrt(1/T))) )
##
getSBQ <- function(rT, W, c, rstar){
  sT <- 1:(rT-1)
  bound <- 1+2*(rT-sT)/T
  max(abs(W[rT] - W[sT])/bound)
}
##
SIM.initial <- function(c = 0, rstar = 1){
  r <- (1:T)/T
  boundary <- 1+2*r
  sim.dist <- function(x){
    W <- BrownianMotion(T)
    Q <- max(abs(W)/boundary)
    BQ <- max(abs(W)/boundary)
    SBQ <- max(sapply(2:T,getSBQ, W=W, c=c, rstar=rstar))
    c(Q,BQ,SBQ)
  }
  DETECTORS<-mcmapply(sim.dist, 1:MC, mc.cores = CORE)
  results <- numeric(dim(DETECTORS)[1])
  for(i in 1:dim(DETECTORS)[1]) ( results[i] <- quantile(DETECTORS[i,], 0.95) )
  results
}
##
CRIT <- SIM.initial(0)
CRIT
## ##################################
## ##################################
end<-Sys.time()
end-start
