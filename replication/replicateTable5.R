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
MC = 100000
## ##################################
##
simM1 <- function(T, m=10, tstar = m){
  mT <- floor(m*T)
  e = rnorm(mT,0,1)
  y = c(rep(0,floor(tstar*T)),rep(0.8,floor((m-tstar)*T))) + e
  model <- y ~ 1
  SBQ = backCUSUM::SBQ.mon(model, T, m=m)$detector.scaled
  Q.det <- backCUSUM::Q.mon.lin(model, T, m=m)$detector
  r <- (1:mT)/T
  boundary.lin <- 1+2*(r-1)[(T+1):mT]
  Q <- Q.det/boundary.lin
  return(
    c(
      which(SBQ > 1.202)[1] + T,
      which(Q > 0.947)[1] + T
    )
  )
}
##
simM2 <- function(T, m=10, tstar = m){
  mT <- floor(m*T)
  e = rnorm(mT,0,1)
  x=filter(rnorm(mT,0,1),0.5,"recursive")
  y=1 + x*c(rep(0,floor(tstar*T)),rep(0.8,floor((m-tstar)*T)))+e
  model <- y ~ 1 + x
  SBQ = backCUSUM::SBQ.mon(model, T, m=m)$detector.scaled
  Q.det <- backCUSUM::Q.mon.lin(model, T, m=m)$detector
  r <- (1:mT)/T
  boundary.lin <- 1+2*(r-1)[(T+1):mT]
  Q <- Q.det/boundary.lin
  return(
    c(
      which(SBQ > 1.274)[1] + T,
      which(Q > 1.034)[1] + T
    )
  )
}
##
simM3 <- function(T, m=10, tstar = m){
  mT <- floor(m*T)
  e = rnorm(mT,0,1)
  emu = e + c(rep(0,floor(tstar*T)),rep(0.8,floor((m-tstar)*T)))
  y=filter(emu,0.5,"recursive")
  model = y[2:mT] ~ 1 + y[1:(mT-1)]
  H = matrix(c(1,0), ncol = 1)
  SBQ = backCUSUM::SBQ.mon(model, T, m=m, H=H)$detector.scaled
  Q.det <- backCUSUM::Q.mon.lin(model, T, m=m, H=H)$detector
  r <- (1:(mT-1))/T
  boundary.lin <- 1+2*(r-1)[(T+1):(mT-1)]
  Q <- Q.det/boundary.lin
  return(
    c(
      which(SBQ > 1.202)[1] + T,
      which(Q > 0.947)[1] + T
    )
  )
}
##
sim.sizedel <- function(T,j,m,tstar){
  if(j == 1){
    Statistics = parSapply(cl,rep(T,MC), simM1, m=m, tstar = tstar)
  } else if(j==2){
    Statistics <- parSapply(cl,rep(T,MC), simM2, m=m, tstar = tstar)
  } else {
    Statistics <- parSapply(cl,rep(T,MC), simM3, m=m, tstar = tstar)
  }
  rejectionrate.SBQ = sum(!is.na(Statistics[1,]))/MC
  rejectionrate.Q = sum(!is.na(Statistics[2,]))/MC
  CDI.SBQ = which(Statistics[1,] > (tstar*T)) # correct detections indices
  CDI.Q = which(Statistics[2,] > (tstar*T)) # correct detections indices
  m.Delay.SBQ =  mean(Statistics[1,CDI.SBQ] - (tstar*T))
  m.Delay.Q =  mean(Statistics[2,CDI.Q] - (tstar*T))
  return(
    c(
      rejectionrate.SBQ,
      rejectionrate.Q,
      m.Delay.SBQ,
      m.Delay.Q
    )
  )
}
##
m=2
taustars = c(m, 1.1, 1.2, 1.3, 1.4, 1.5)
M1.l = lapply(taustars, T = 200, j=1, m=m, sim.sizedel)
M2.l = lapply(taustars, T = 200, j=2, m=m, sim.sizedel)
M3.l = lapply(taustars, T = 200, j=3, m=m, sim.sizedel)
##
rnames = c("size", taustars[-1])
cnames = c("SBQ","Q", "SBQ", "Q")
##
M1 = matrix(unlist(M1.l), nrow = length(taustars), byrow=TRUE, dimnames = list(rnames,cnames))
M2 = matrix(unlist(M2.l), nrow = length(taustars), byrow=TRUE, dimnames = list(rnames,cnames))
M3 = matrix(unlist(M3.l), nrow = length(taustars), byrow=TRUE, dimnames = list(rnames,cnames))
M1
M2
M3
##
tab1 = rbind(M1[1,1:2]*100, M1[-1,3:4])
tab2 = rbind(M2[1,1:2]*100, M2[-1,3:4])
tab3 = rbind(M3[1,1:2]*100, M3[-1,3:4])
##
table=round(cbind(tab1, tab2, tab3),1)
rownames(table) = rnames
table
write.csv(table,file="./results/Table5.csv")
##
Sys.time()-start
stopCluster(cl)
