## ####################################################################
## ####################################################################
## Supplement for
## "Backward CUSUM for Testing and Monitoring Structural Change"
## by Sven Otto and Jörg Breitung.
## This R-script allows to simulate critical values
## for the SBQ M-test with finite monitoring horizon
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
MC = 10 # number of Monte Carlo repetitions
G = 5000  # grid on which the Wiener process is approximated
m = c(1.2,1.4,1.6,1.8,2,4,10) # monitoring horizons
T = floor(G/(tail(m,1)-1))
## ##################################
sim.SBQ = function(T,k,m){
  M = (tail(m,1)-1)*T
  BrownianMotion = function(T)  ( cumsum(c(rnorm(M,0,sqrt(1/T)))) )
  get.innermax = function(t,W) max(abs(W[t]-c(0,W[1:(t-1)]))/(1+2*(t/T-(0:(t-1))/T)))
  detector.dim = function(j){
    W = BrownianMotion(T)
    cummax(sapply(1:M, get.innermax, W=W))
  }
  SBQcummax = apply(sapply(1:k,detector.dim),1,cummax)
  return(SBQcummax[,ceiling((m-1)*T)])
}
##
out = parSapply(cl,rep(T,MC),sim.SBQ, k=k, m=m,simplify="array")
knames = paste("k=",1:k,sep="")
anames = paste("alpha=",alphalevels,"%",sep="")
mnames = paste("m=",m,sep="")
rownames(out) = knames
critval.2sided = apply(out, c(1,2), quantile, 1-alphalevels)
critval.rsided = apply(out, c(1,2), quantile, (1-2*alphalevels))
critval.lsided = -apply(out, c(1,2), quantile, (1-2*alphalevels))
rownames(critval.2sided) = anames
rownames(critval.rsided) = anames
rownames(critval.lsided) = anames
## ####################################################################
## Critical values for two-sided, right-sided, and left-sided tests:
## ####################################################################
converttolist = function(array){
  l = lapply(seq(dim(array)[3]), function(x) round(array[ , , x],3))
  names(l) = mnames
  l
}
critval.2sided = converttolist(critval.2sided)
critval.rsided = converttolist(critval.rsided)
critval.lsided = converttolist(critval.lsided)
critval.2sided
critval.rsided
critval.lsided
##
Sys.time()-start
stopCluster(cl)
