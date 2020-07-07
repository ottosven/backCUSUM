## ####################################################################
## ####################################################################
## Supplement for
## "Backward CUSUM for Testing and Monitoring Structural Change"
## by Sven Otto and Jörg Breitung.
## This R-script allows to reproduce Table 3.
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
MC <- 100000
## ## ##################################
simM1 <- function(T){
  u <- rnorm(T,0,1)
  y <- 2 + u
  model <- y ~ 1
  return(c(
    backCUSUM::Q.test(model)$statistic,
    backCUSUM::BQ.test(model)$statistic,
    backCUSUM::SBQ.test(model)$statistic
  ))
}
##
simM2 <- function(T){
  u <- rnorm(T,0,1)
  x <- rnorm(T,0,1)
  y <- 2 + x + u
  model <- y ~ 1 + x
  return(c(
    backCUSUM::Q.test(model)$statistic,
    backCUSUM::BQ.test(model)$statistic,
    backCUSUM::SBQ.test(model)$statistic
  ))
}
##
sim.size <- function(T,k){
  crit <- c(get.crit.Q(k)[2], get.crit.BQ(k)[2], get.crit.SBQ(k)[2])
  if(k == 1){
    Statistics <- parSapply(cl,rep(T,MC), simM1)
  } else {
    Statistics <- parSapply(cl,rep(T,MC), simM2)
  }
  size <- c(
    length(which(Statistics[1,] > crit[1]))/MC,
    length(which(Statistics[2,] > crit[2]))/MC,
    length(which(Statistics[3,] > crit[3]))/MC
  )
  return(size)
}
##
size.100.M1 <- sim.size(100,1)
size.200.M1 <- sim.size(200,1)
size.500.M1 <- sim.size(500,1)
size.100.M2 <- sim.size(100,2)
size.200.M2 <- sim.size(200,2)
size.500.M2 <- sim.size(500,2)
##
resultsM1 <- c(size.100.M1, size.200.M1, size.500.M1)
resultsM2 <- c(size.100.M2, size.200.M2, size.500.M2)
table3 <- matrix(c(resultsM1, resultsM2), nrow = 2, byrow = TRUE, dimnames = list(c("Model I", "Model II"), rep(c("Q", "BQ", "SBQ"),3)))
table3
write.table(table3,file="./results/Table3.csv", col.names=NA)
Sys.time()-start
stopCluster(cl)
