## ####################################################################
## ####################################################################
## Supplement for
## "Backward CUSUM for Testing and Monitoring Structural Change"
## by Sven Otto and Jörg Breitung.
## This R-script allows to reproduce Table 1.
## ####################################################################
## ####################################################################
rm(list=ls())
start<-Sys.time()
library(backCUSUM)  # install package with remotes::install_github("ottosven/backCUSUM")
library(parallel)
## ##################################
## Cluster setup
## ##################################
if(is.na(strtoi(Sys.getenv(c("SLURM_NTASKS"))))){
  cl = makeCluster(detectCores()-1)
} else {
  ntasks = strtoi(Sys.getenv(c("SLURM_NTASKS")))
  nslaves = ntasks-1
  cl = makeCluster(nslaves, type="MPI")
}
## ##################################
## Reproducible random number generator
## ##################################
RNGkind("L'Ecuyer-CMRG")
set.seed(42)
snow::clusterSetupRNG(cl)
## ##################################
## Simulation setting
## ##################################
MC <- 100000 # number of Monte Carlo repetitions
T <- 50000  # grid on which the Wiener process is approximated
k = 10 # number of dimensions
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
result = apply(out, 1, quantile, c(0.9, 0.95, 0.99))
colnames(result) = paste("k=",1:k,sep="")
rownames(result) = paste("alpha=",c(10,5,1),"%",sep="")
result
##
write.csv(round(result,3), file="./results/Table1.csv")
##
Sys.time()-start
stopCluster(cl)
