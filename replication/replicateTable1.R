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
# k <- 1
# m <- Inf
## ##################################
k
m
MC <- 100000
T <- 10000
## ##################################
sim.crit.Q <- function(k, m = 2,  MC = 1000, T = 10000, levels = c(0.1, 0.05, 0.01, 0.001), alternative = "two.sided", CORE = 1){
  r <- (1:T)/T
  boundary.transf <- 1+r
  if(m == Inf){
    MRange <- 1:T
  } else {
    MRange <- 1:floor(T*(m-1)/m)
  }
  boundary.mrange <- boundary.transf[MRange]
  BrownianMotion <- function(T)  ( cumsum(rnorm(T,0,sqrt(1/T))) )
  sim.dist.Q <- function(...){
    W <- sapply(rep(T,k), BrownianMotion)
    B <- W - t(W[T,]*matrix(rep(r,k), nrow = k, byrow = TRUE))
    B.mrange <- matrix(B[MRange,], ncol = k)
    Q <- max((apply(abs(B.mrange), 1, max)/boundary.mrange))
    return(Q)
  }
  realizations <- parallel::mcmapply(sim.dist.Q,1:MC,mc.cores = getOption("mc.cores", CORE))
  if(alternative == "one.sided"){
    quantiles <- quantile(realizations, 1-levels*2)
  } else {
    quantiles <- quantile(realizations, 1-levels)
  }
  names(quantiles) <- levels
  return(quantiles)
}


sim.dist <- function(T, k, B.end){
  B <- sapply(1:k, function(...) e1071::rbridge(end=1, frequency=T)[1:B.end])
  max((apply(abs(B[1:B.end,,drop=F]), 1, max)/(1+(1:B.end)/T)))
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
write.table(result,file=paste("./critval/Q-crit.csv", sep=''),append=T,col.names=F,row.names=F)
Sys.time()-start
stopCluster(cl)