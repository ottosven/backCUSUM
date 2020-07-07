## ####################################################################
## ####################################################################
## Supplement for
## "Backward CUSUM for Testing and Monitoring Structural Change"
## by Sven Otto and Jörg Breitung.
## This R-script allows to reproduce Table 2.
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
  ntasks <- strtoi(Sys.getenv(c("SLURM_NTASKS")))
  nslaves <- ntasks-1
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
input <- commandArgs(trailingOnly = TRUE)
input <- as.numeric(input) ##convert string to number
input
k <- input[1]
m <- input[2]
## ##################################
## Please insert the specification for k and m for
## which the critical values should be simulated:
## k <- 1
## m <- Inf
## ##################################
k
m
MC <- 100000
T <- 10000
## ##################################
sim.dist <- function(T, k, B.end){
  get.innermax <- function(t, B){
    R <- (1:T)/T
    S <- (0:(T-1))/T
    max(apply(abs(B[t,] - (1-R[t])/(1-S[1:t])*t(B[1:t,])),2,max)/(1-R[t] + 2*(R[t] - S[1:t])/(1-S[1:t])))
  }
  B <- sapply(1:k, function(...) e1071::rbridge(end=1, frequency=T)[1:B.end])
  max(sapply(1:B.end, get.innermax, B=B))
}
##
sim.crit <- function(k, T, m){
  if(m == Inf){
    B.end = T
  } else {
    B.end = floor(T*(m-1)/m)
  }
  levels <- c(0.2, 0.1, 0.05, 0.02, 0.01, 0.002, 0.001)
  realizations <- parSapply(cl,rep(T,MC),sim.dist, k=k, B.end=B.end)
  quantiles <- quantile(realizations, 1-levels)
  names(quantiles) <- levels
  quantiles
}
##
crits <- sim.crit(k,T,m)
k
m
result <- matrix(c(k,m,MC,T,crits), ncol = length(crits)+4)
result
write.table(result,file=paste("./critval/SBQ-crit.csv", sep=''),append=T,col.names=F,row.names=F)
Sys.time()-start
stopCluster(cl)