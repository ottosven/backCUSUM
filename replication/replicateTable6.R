## ####################################################################
## ####################################################################
## Supplement for
## "Backward CUSUM for Testing and Monitoring Structural Change"
## by Sven Otto and Jörg Breitung.
## This R-script allows to reproduce Table 6.
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
##
sim1.H0 <- function(T){
  m <- 20
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
  result <- c(max(SBQ), max(Q), max(CSW))
}
##
sim1.H1 <- function(T, taustar, crit){
  m <- 20
  mT <- floor(m*T)
  u <- rnorm(mT,0,1)
  brea <- floor(taustar*T)
  y <- c(rep(2,brea),rep(2+0.8,mT-brea)) + u
  model <- y ~ 1
  SBQ <- backCUSUM::SBQ.mon(model, T)$detector.scaled
  Q.det <- backCUSUM::Q.mon.lin(model, T)$detector
  r <- (1:mT)/T
  boundary.lin <- 1+2*(r-1)[(T+1):mT]
  Q <- Q.det/boundary.lin
  boundary.CSW <- sqrt(r[(T+1):mT]*(log(r[(T+1):mT]/0.05^2)))
  CSW <- Q.det/boundary.CSW
  result <- unname(c(
    which(SBQ > crit[1])[1] + T,
    which(Q > crit[2])[1] + T,
    which(CSW > crit[3])[1] + T
  ))
  return(result)
}
##
sim2.H0<- function(T){
  m <- 20
  mT <- floor(m*T)
  u <- rnorm(mT,0,1)
  x <- rnorm(mT,0,1)
  y <- 2 + x + u
  model <- y ~ 1 + x
  SBQ <- backCUSUM::SBQ.mon(model, T)$detector.scaled
  r <- (1:mT)/T
  boundary.lin <- 1+2*(r-1)[(T+1):mT]
  Q <- backCUSUM::Q.mon.lin(model, T)$detector/boundary.lin
  return(c(max(SBQ), max(Q)))
}
##
sim2.H1<- function(T, taustar, crit){
  m <- 20
  mT <- floor(m*T)
  u <- rnorm(mT,0,1)
  x <- rnorm(mT,0,1)
  brea <- floor(taustar*T)
  y <- 2 + c(rep(1,brea),rep(1+0.8,mT-brea))*x + u
  model <- y ~ 1 + x
  SBQ <- backCUSUM::SBQ.mon(model, T)$detector.scaled
  r <- (1:mT)/T
  boundary.lin <- 1+2*(r-1)[(T+1):mT]
  Q <- backCUSUM::Q.mon.lin(model, T)$detector/boundary.lin
  return(c(
    which(SBQ > crit[1])[1] + T,
    which(Q > crit[2])[1] + T)
  )
}
##
sim.crit <- function(T,k){
  if(k == 1){
    Statistics <- parSapply(cl,rep(T,MC), sim1.H0)
    CRIT <- c(
      quantile(Statistics[1,], 0.95),
      quantile(Statistics[2,], 0.95),
      quantile(Statistics[3,], 0.95)
    )
  } else {
    Statistics <- parSapply(cl,rep(T,MC), sim2.H0)
    CRIT <- c(
      quantile(Statistics[1,], 0.95),
      quantile(Statistics[2,], 0.95)
    )
  }
  return(CRIT)
}
##
sim.pow <- function(T,k,crit){
  taustar <- c(1.5, 2, 2.5, 3, 5, 10)
  if(k == 1){
    delay <- matrix(nrow = 6, ncol = 3)
    rownames(delay) <- taustar
    colnames(delay) <- c(paste('SBQ',k) , paste('Q',k), 'CSW')
    for(i in 1:6){
      Statistics <- parSapply(cl,rep(T,MC), sim1.H1, taustar = taustar[i], crit = crit)
      for(j in 1:3){
        CDI <- which(Statistics[j,] > (taustar*T)[i]) # correct detections indices
        delay[i,j] <- mean(Statistics[j,CDI] - (taustar*T)[i])
      }
    }
  } else {
    delay <- matrix(nrow = 6, ncol = 2)
    rownames(delay) <- taustar
    colnames(delay) <- c(paste('SBQ',k) , paste('Q',k))
    for(i in 1:6){
      Statistics <- parSapply(cl,rep(T,MC), sim2.H1, taustar = taustar[i], crit = crit)
      for(j in 1:2){
        CDI <- which(Statistics[j,] > (taustar*T)[i]) # correct detections indices
        delay[i,j] <- mean(Statistics[j,CDI] - (taustar*T)[i])
      }
    }
  }
  return(delay)
}
##
crit1 <- sim.crit(100,1)
crit2 <- sim.crit(100,2)
pow1 <- sim.pow(100,1,crit1)
pow2 <- sim.pow(100,2,crit2)
table6 <- cbind(pow1[1:3,], pow2[1:3,], pow1[4:6,], pow2[4:6,])
table6
write.table(table6,file="./results/Table6.csv", col.names=NA)
Sys.time()-start
stopCluster(cl)
