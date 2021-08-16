## ####################################################################
## ####################################################################
## Supplement for
## "Backward CUSUM for Testing and Monitoring Structural Change"
## by Sven Otto and Jörg Breitung.
## This R-script allows to reproduce Table 7.
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
breakpoints <- c(0.5, 0.65, 0.8, 0.85, 0.9, 0.95, 0.97, 0.99)
MC <- 100000
## ## ##################################
sim <- function(T, taustar){
  u <- rnorm(T,0,1)
  mu <- rep(2,T)
  mu[(floor(taustar*T)+1):T] <- 2.8
  y <- mu + u
  c(backCUSUM::breakpoint.est(y~1, type = "BQ"), backCUSUM::breakpoint.est(y~1, type = "ML"))
}
##
sim.mse <- function(taustar, T){
  Statistics<-parSapply(cl,rep(T,MC),sim, taustar=taustar)  # FOR MPI CLUSTER ONLY
  statnames <- c('BQ', 'ML')
  Bias <- rowMeans(Statistics/T - taustar)
  names(Bias) <- statnames
  MSE <- sqrt(rowMeans((Statistics/T-taustar)^2))
  names(MSE) <- statnames
  rbind(Bias, MSE)
}
##
result <- lapply(breakpoints, T = 100, sim.mse)
BiasMSE <- matrix(ncol = 4, nrow = length(breakpoints))
rownames(BiasMSE) <- breakpoints
colnames(BiasMSE) <- c('BQ (Bias)', 'ML (Bias)', 'BQ (RMSE)', 'ML (RMSE)')
for(i in 1:length(breakpoints)){
  BiasMSE[i,c(1,2)] <- result[[i]][1,]
  BiasMSE[i,c(3,4)] <- result[[i]][2,]
}
results100 <- round(BiasMSE,2)
##
result <- lapply(breakpoints, T = 200, sim.mse)
BiasMSE <- matrix(ncol = 4, nrow = length(breakpoints))
rownames(BiasMSE) <- breakpoints
colnames(BiasMSE) <- c('BQ (Bias)', 'ML (Bias)', 'BQ (RMSE)', 'ML (RMSE)')
for(i in 1:length(breakpoints)){
  BiasMSE[i,c(1,2)] <- result[[i]][1,]
  BiasMSE[i,c(3,4)] <- result[[i]][2,]
}
results200 <- round(BiasMSE,2)
##
table7 <- cbind(results100, results200)
table7
write.table(table7,file="./results/Table7.csv", col.names=NA)
Sys.time()-start
stopCluster(cl)
