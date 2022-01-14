## ####################################################################
## ####################################################################
## Supplement for
## "Backward CUSUM for Testing and Monitoring Structural Change"
## by Sven Otto and Jörg Breitung.
## This R-script allows to simulate critical values
## for the Q M-test with infinite monitoring horizon
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
## Simulation setting (set MC=100000 and G=50000 to replicate paper values)
## ####################################################################
MC = 100 # number of Monte Carlo repetitions
T = 50000  # grid on which the Wiener process is approximated
## ##################################
sim.Q = function(T,k){
  BrownianMotion = function(T)  ( cumsum(rnorm(T,0,sqrt(1/T))) )
  detector.dim = function(i){
    W = BrownianMotion(T)
    B = W-((1:T)/T)*W[T]
    r=(1:(T-1))/T
    r.tr = r/(1-r)
    max(abs(B[1:(T-1)])/( (1-r)*(1+2*r.tr)))
  }
  cummax(sapply(1:k,detector.dim))
}
##
out = parSapply(cl,rep(T,MC),sim.Q, k=k)
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
