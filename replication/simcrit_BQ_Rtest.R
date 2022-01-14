## ####################################################################
## ####################################################################
## Supplement for
## "Backward CUSUM for Testing and Monitoring Structural Change"
## by Sven Otto and Jörg Breitung.
## This R-script allows to simulate critical values
## for the Q and BQ R-tests
## ####################################################################
## ####################################################################
rm(list=ls())
start<-Sys.time()
library(backCUSUM)
library(parallel)
## ####################################################################
## Parallelization setup for both SLURM and local machine
## ####################################################################
if(is.na(strtoi(Sys.getenv(c("SLURM_NTASKS"))))){
  cl = makeCluster(detectCores()-1)
} else {
  ntasks = strtoi(Sys.getenv(c("SLURM_NTASKS")))
  nslaves = ntasks-1
  cl = makeCluster(nslaves, type="MPI")
}
## ####################################################################
## Reproducible random number generator
## ####################################################################
RNGkind("L'Ecuyer-CMRG")
set.seed(42)
snow::clusterSetupRNG(cl)
## ####################################################################
## Setting for the critical values
## ####################################################################
k = 10 # maximum number of regressors/restrictions
alphalevels = c(0.1, 0.05, 0.01, 0.001) # significance levels
## ####################################################################
## Simulation setting (set MC=100000 and T=50000 to replicate paper values)
## ####################################################################
MC = 100 # number of Monte Carlo repetitions
T <- 50000  # grid on which the Wiener process is approximated
## ##################################
sim.BQ = function(T,k){
  BrownianMotion = function(T)  ( cumsum(rnorm(T,0,sqrt(1/T))) )
  W = sapply(rep(T,k), BrownianMotion)
  Wmax.multi = apply(abs(W), 1, cummax)
  bound = 1 + 2*(1:T)/T
  bound.multi = matrix(rep(bound, k), ncol=T, byrow = TRUE)
  BQ.multi = apply(Wmax.multi/bound.multi, 1, max)
  return(BQ.multi)
}
##
out = parSapply(cl,rep(T,MC),sim.BQ, k=k)
knames = paste("k=",1:k,sep="")
anames = paste("alpha=",alphalevels,"%",sep="")
rownames(out) = knames
critval.2sided = apply(out, 1, quantile, 1-alphalevels)
critval.rsided = apply(out, 1, quantile, (1-2*alphalevels))
critval.lsided = -apply(out, 1, quantile, (1-2*alphalevels))
rownames(critval.2sided) = anames
rownames(critval.rsided) = anames
rownames(critval.lsided) = anames
## ####################################################################
## Critical values for two-sided, right-sided, and left-sided tests:
## ####################################################################
round(critval.2sided,3)
round(critval.rsided,3)
round(critval.lsided,3)
##
Sys.time()-start
stopCluster(cl)
