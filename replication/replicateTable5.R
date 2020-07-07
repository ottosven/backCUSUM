## ####################################################################
## ####################################################################
## Supplement for
## "Backward CUSUM for Testing and Monitoring Structural Change"
## by Sven Otto and Jörg Breitung.
## This R-script allows to reproduce Table 5.
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
## ##################################
sim1 <- function(T){
  m <- 10
  mT <- floor(m*T)
  u <- rnorm(mT,0,1)
  y <- 2 + u
  model <- y ~ 1
  SBQ <- backCUSUM::SBQ.mon(model, T)$detector.scaled
  Q.det <- backCUSUM::Q.mon.lin(model, T)$detector
  r <- (1:mT)/T
  boundary.lin <- 1+2*(r-1)[(T+1):mT]
  Q <- Q.det/boundary.lin
  boundary.CSW <- sqrt(r[(T+1):mT]*(log(r[(T+1):mT]/0.05^2)))
  CSW <- Q.det/boundary.CSW
  return(c(
    which(SBQ > backCUSUM::get.crit.SBQ.mon(1)$crit[2])[1] + T,
    which(Q > backCUSUM::get.crit.Q.mon(1)$crit[2])[1] + T,
    which(CSW > 1)[1] + T
  ))
}
##
sim2 <- function(T){
  m <- 10
  mT <- floor(m*T)
  u <- rnorm(mT,0,1)
  x <- rnorm(mT,0,1)
  y <- 2 + x + u
  model <- y ~ 1 + x
  SBQ <- backCUSUM::SBQ.mon(model, T)$detector.scaled
  r <- (1:mT)/T
  boundary.lin <- 1+2*(r-1)[(T+1):mT]
  Q <- backCUSUM::Q.mon.lin(model, T)$detector/boundary.lin
  return(c(
    which(SBQ > backCUSUM::get.crit.SBQ.mon(2)$crit[2])[1] + T,
    which(Q > backCUSUM::get.crit.Q.mon(2)$crit[2])[1] + T
  ))
}
##
sim.size <- function(T,k){
  m <- c(1.5, 2, 4, 10)
  if(k == 1){
    Statistics <- parSapply(cl,rep(T,MC), sim1)
    result <- matrix(nrow = 4, ncol = 3)
    rownames(result) <- m
    colnames(result) <- c(paste("SBQ",T,k), paste("Q",T,k), paste("CSW",T))
    for(i in 1:4){
      for(j in 1:3){
        result[i,j] <- length(which(Statistics[j,] < m[i]*T))/MC
      }
    }
  } else {
    Statistics <- parSapply(cl,rep(T,MC), sim2)
    result <- matrix(nrow = 4, ncol = 2)
    rownames(result) <- m
    colnames(result) <- c(paste("SBQ",T,k), paste("Q",T,k))
    for(i in 1:4){
      for(j in 1:2){
        result[i,j] <- length(which(Statistics[j,] < m[i]*T))/MC
      }
    }
  }
  return(result)
}
##
sizes1 <- sim.size(100,1)
sizes2 <- sim.size(300,1)
sizes3 <- sim.size(100,2)
sizes4 <- sim.size(200,2)
sizes5 <- sim.size(300,2)
##
table5 <- cbind(sizes1, sizes2, sizes3, sizes4, sizes5)
table5
write.table(table5,file="./results/Table5.csv", col.names=NA)
Sys.time()-start
stopCluster(cl)