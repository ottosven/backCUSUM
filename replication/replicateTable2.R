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
MC = 100000 # number of Monte Carlo repetitions
G = 50000  # grid on which the Wiener process is approximated
k = 8 # number of dimensions
m = c(1.2,1.4,1.6,1.8,2,4,10) # monitoring horizons
T = floor(G/(tail(m,1)-1))
## ##################################
sim.SBQ = function(T,k,m){
  M = (tail(m,1)-1)*T
  BrownianMotion = function(T)  ( cumsum(c(rnorm(M,0,sqrt(1/T)))) )
  get.innermax = function(t,W) max(abs(W[t]-c(0,W[1:(t-1)]))/(1+2*(t/T-(0:(t-1))/T)))
  detector.dim = function(j){
    W = BrownianMotion(T)
    cummax(sapply(1:M, get.innermax, W=W))
  }
  SBQcummax = apply(sapply(1:k,detector.dim),1,cummax)
  return(SBQcummax[,floor((m-1)*T)])
}
##
out = parSapply(cl,rep(T,MC),sim.SBQ, k=k, m=m,simplify="array")
result= apply(out,c(1,2),quantile,c(0.9, 0.95, 0.99))
##
uppertable = matrix(nrow = length(m), ncol = 12)
rownames(uppertable) = paste("m=",m,sep="")
colnames(uppertable) = paste("k=",rep(1:4,each=3),";",rep(c(10,5,1),4),"%",sep="")
for(i in 1:length(m)) uppertable[i,] = c(result[,1:4,i])
uppertable
lowertable = matrix(nrow = length(m), ncol = 12)
rownames(lowertable) = paste("m=",m,sep="")
colnames(lowertable) = paste("k=",rep(5:8,each=3),";",rep(c(10,5,1),4),"%",sep="")
for(i in 1:length(m)) lowertable[i,] = c(result[,5:8,i])
lowertable
##
write.csv(round(uppertable,3), file="./results/Table2upper.csv")
write.csv(round(lowertable,3), file="./results/Table2lower.csv")
##
Sys.time()-start
stopCluster(cl)
