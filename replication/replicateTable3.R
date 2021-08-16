## ####################################################################
## ####################################################################
## Supplement for
## "Backward CUSUM for Testing and Monitoring Structural Change"
## by Sven Otto and Jörg Breitung.
## This R-script allows to reproduce Table 3.
## ####################################################################
## ####################################################################
rm(list=ls())
start = Sys.time()
library(backCUSUM)
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
k = 5 # number of dimensions
## ###################################
sim.SBQinf = function(T,k){
  BrownianMotion = function(T)  ( cumsum(rnorm(T,0,sqrt(1/T))) )
  get.innermax = function(t,B) max(abs((1-(0:(t-1))/T)*B[t]-(1-t/T)*c(0,B[1:(t-1)]))/((1-t/T)*(1-(0:(t-1))/T)*sqrt(1/(1-t/T))*(1+2*((1/(1-t/T))-1/(1-(0:(t-1))/T)))))
  detector.dim = function(j){
    W = BrownianMotion(T)
    B = W-((1:T)/T)*W[T]
    max(sapply(1:(T-1), get.innermax, B=B))
  }
  return(cummax(sapply(1:k,detector.dim)))
}
##
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
out.SBQ = parSapply(cl,rep(T,MC),sim.SBQinf,k=k)
result.SBQ = apply(out.SBQ,1,quantile,c(0.9, 0.95, 0.99))
out.Q = parSapply(cl,rep(T,MC),sim.Q, k=k)
result.Q = apply(out.Q, 1, quantile, c(0.9, 0.95, 0.99))
##
result = cbind(result.SBQ, result.Q)
colnames(result) = paste(rep(c("SBQ","Q"),each=k),".k=",rep(1:k,2),sep="")
rownames(result) = paste("alpha=",c(10,5,1),"%",sep="")
result
##
write.csv(round(result,3), file="./results/Table3.csv")
##
Sys.time()-start
stopCluster(cl)
