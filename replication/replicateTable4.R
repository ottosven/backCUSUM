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
taustars <- 1:9/10
MC <- 100000
## ## ##################################
simM1 <- function(T, tstar = 1){
  u <- rnorm(T,0,1)
  y <- 2 + u + c(rep(0,floor(tstar*T)),rep(0.8,floor((1-tstar)*T)))
  model <- y ~ 1
  return(c(
    backCUSUM::Q.test(model)$statistic,
    backCUSUM::BQ.test(model)$statistic,
    backCUSUM::SBQ.test(model)$statistic,
    backCUSUM::sup.wald(model)
  ))
}
##
simM2 <- function(T, tstar = 1){
  u <- rnorm(T,0,1)
  x <- rnorm(T,0,1)
  y <- 2 + x + u + x*c(rep(0,floor(tstar*T)),rep(0.8,floor((1-tstar)*T)))
  model <- y ~ 1 + x
  return(c(
    backCUSUM::Q.test(model)$statistic,
    backCUSUM::BQ.test(model)$statistic,
    backCUSUM::SBQ.test(model)$statistic,
    backCUSUM::sup.wald(model)
  ))
}
##
sim.crit <- function(T,k){
  if(k == 1){
    Statistics <- parSapply(cl,rep(T,MC), simM1)
    CRIT <- c(
      quantile(Statistics[1,], 0.95),
      quantile(Statistics[2,], 0.95),
      quantile(Statistics[3,], 0.95),
      quantile(Statistics[4,], 0.95)
    )
  } else {
    Statistics <- parSapply(cl,rep(T,MC), simM2)
    CRIT <- c(
      quantile(Statistics[1,], 0.95),
      quantile(Statistics[2,], 0.95),
      quantile(Statistics[3,], 0.95),
      quantile(Statistics[4,], 0.95)
    )
  }
  return(CRIT)
}
##
sim.pow <- function(taustar,T,k,crit){
  if(k == 1){
    Statistics <- parSapply(cl,rep(T,MC), simM1, tstar = taustar)
  } else {
    Statistics <- parSapply(cl,rep(T,MC), simM2, tstar = taustar)
  }
  power <- c(
    length(which(Statistics[1,] > crit[1]))/MC,
    length(which(Statistics[2,] > crit[2]))/MC,
    length(which(Statistics[3,] > crit[3]))/MC,
    length(which(Statistics[4,] > crit[4]))/MC
  )
  return(power)
}
##
crit1 <- sim.crit(100,1)
crit2 <- sim.crit(100,2)
pow1 <- lapply(taustars, T = 100, k = 1, crit = crit1, sim.pow)
pow2 <- lapply(taustars, T = 100, k = 2, crit = crit2, sim.pow)
table4 <- matrix(ncol= 8, nrow = length(taustars), dimnames = list(taustars, rep(c("Q", "BQ", "SBQ", "supW"),2)))
for(i in 1:length(taustars)){
  table4[i, 1:4] <- pow1[[i]]
  table4[i, 5:8] <- pow2[[i]]
}
table4
write.table(table4,file="./results/Table4.csv", col.names=NA)
Sys.time()-start
stopCluster(cl)
