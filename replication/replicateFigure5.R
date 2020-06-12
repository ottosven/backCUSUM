## ####################################################################
## ####################################################################
## Supplement for
## "Backward CUSUM for Testing and Monitoring Structural Change"
## by Sven Otto and Jörg Breitung.
## This R-script allows to reproduce Figure 5.
## ####################################################################
## ####################################################################
rm(list=ls())
start<-Sys.time()
library("parallel")
## ##################################
## Reproducible random number
## generator for parallelization
## ##################################
RNGkind("L'Ecuyer-CMRG")
set.seed(42)
CORE <- detectCores(all.tests = FALSE, logical = TRUE)
## ##################################
## Simulation:
## ##################################
## Monte Carlo repetitions:
MC <- 100
## ##################################
## For replication please set the number of repetitions to
## MC <- 10000000
## Warning! This process takes a very long time
## ##################################
T <- 1000
## Simulated critical values:
CRIT <- c(0.9376, 0.9376, 1.1875)
## Simulated critical values can be replicated using replicateFigure3_data_init.R
##
BrownianMotion <- function(T)  ( cumsum(rnorm(T,0,sqrt(1/T))) )
##
getSBQ <- function(rT, W){
  sT <- 1:(rT-1)
  bound <- 1+2*(rT-sT)/T
  max(abs(W[rT] - W[sT])/bound)
}
##
## Upper row
##
r <- (1:T)/T
boundary <- 1+2*r
sim.dist1 <- function(x){
  W <- BrownianMotion(T)
  allQ <- abs(W)/boundary
  which(allQ > CRIT[1])
  allBQ <- rev(allQ)
  which(allBQ > CRIT[2])
  allSBQ <- c(0,sapply(2:T,getSBQ, W=W))
  Q <- which(allQ > CRIT[1])[1]
  BQ <- rev(which(allBQ > CRIT[2]))[1]
  SBQ <- which(allSBQ > CRIT[3])[1]
  c(Q,BQ,SBQ)
}
##
REJECTIONS1<-mcmapply(sim.dist1, 1:MC, mc.cores = CORE)
##
## Botttom row
##
m <- 10
r <- (m-1)*(1:T)/T
boundary <- 1+2*r
boundaryChu <- sqrt((r+1)*(log((r+1)/0.05^2)))
CRIT2 <- c(0.9328, 0.8854, 1.4116)
sim.dist2 <- function(...){
  W <- sqrt(m-1)*BrownianMotion(T)
  allQ <- abs(W)/boundary
  allQChu <- abs(W)/boundaryChu
  allSBQ <- c(0,sapply(2:T,getSBQ, W=W))
  Q <- which(allQ > CRIT2[1])[1]
  QChu <- which(allQChu > CRIT2[2])[1]
  SBQ <- which(allSBQ > CRIT2[3])[1]
  c(SBQ,Q,QChu)
}
REJECTIONS2<-mcmapply(sim.dist2, 1:MC, mc.cores = CORE)
## ##################################
## Figure 5
## ##################################
pdf('figure5.pdf')
namesHist1 <- c('Forward CUSUM', 'Backward CUSUM', 'Stacked backward CUSUM')
for(i in 1:dim(REJECTIONS1)[1]) ( hist(REJECTIONS1[i,]/T, breaks = 0:20/20, freq = FALSE, main= namesHist1[i], xlab ='time point of rejection', ylab = 'density', cex.lab = 1.5, cex.axis=1.7, cex.main = 2) )
namesHist2 <- c('Stacked backward CUSUM', 'Forward CUSUM', 'Forward CUSUM (radical boundary)')
for(i in 1:dim(REJECTIONS2)[1]) ( hist(1+(m-1)*REJECTIONS2[i,]/T, breaks = 1 + (m-1)*(0:20)/20, freq = FALSE, main= namesHist2[i], xlab ='time point of rejection', ylab = 'density', cex.lab = 1.5, cex.axis=1.7, cex.main = 2) )
dev.off()
## ##################################
## ##################################
end<-Sys.time()
end-start
