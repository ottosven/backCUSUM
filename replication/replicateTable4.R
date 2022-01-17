## ####################################################################
## ####################################################################
## Supplement for
## "Backward CUSUM for Testing and Monitoring Structural Change"
## by Sven Otto and Jörg Breitung.
## This R-script allows to reproduce Table 4.
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
## ## ##################################
simM1 <- function(T, tstar = 1){
  sup.wald <- function(formula, eps = 0.15){
    X <- model.matrix(formula)
    y <- model.frame(formula)[,1]
    T <- dim(X)[1]
    RSS0 <- deviance(lm(formula))
    wald <- function(t){
      RSS1 <- deviance(lm(y[1:t] ~ X[1:t,]))
      RSS2 <- deviance(lm(y[(t+1):T] ~ X[(t+1):T,]))
      return(T*(RSS0 - (RSS1 + RSS2))/(RSS1 + RSS2))
    }
    return(max(sapply(floor(eps*T):(T-floor(eps*T)), wald)))
  }
  e = rnorm(T,0,1)
  y = c(rep(0,floor(tstar*T)),rep(0.8,floor((1-tstar)*T))) + e
  model <- y ~ 1
  return(c(
    backCUSUM::BQ.test(model)$statistic,
    backCUSUM::SBQ.test(model)$statistic,
    backCUSUM::Q.test(model)$statistic,
    sup.wald(model)
  ))
}
##
simM2 <- function(T, tstar = 1){
  sup.wald <- function(formula, eps = 0.15){
    X <- model.matrix(formula)
    y <- model.frame(formula)[,1]
    T <- dim(X)[1]
    RSS0 <- deviance(lm(formula))
    wald <- function(t){
      RSS1 <- deviance(lm(y[1:t] ~ X[1:t,]))
      RSS2 <- deviance(lm(y[(t+1):T] ~ X[(t+1):T,]))
      return(T*(RSS0 - (RSS1 + RSS2))/(RSS1 + RSS2))
    }
    return(max(sapply(floor(eps*T):(T-floor(eps*T)), wald)))
  }
  e = rnorm(T,0,1)
  x=filter(rnorm(T,0,1),0.5,"recursive")
  y=1 + x*c(rep(0,floor(tstar*T)),rep(0.8,floor((1-tstar)*T)))+e
  model <- y ~ 1 + x
  return(c(
    backCUSUM::BQ.test(model)$statistic,
    backCUSUM::SBQ.test(model)$statistic,
    backCUSUM::Q.test(model)$statistic,
    sup.wald(model)
  ))
}
##
simM3 <- function(T, tstar = 1){
  sup.wald <- function(formula, eps = 0.15){
    X <- model.matrix(formula)
    y <- model.frame(formula)[,1]
    T <- dim(X)[1]
    RSS0 <- deviance(lm(formula))
    wald <- function(t){
      RSS1 <- deviance(lm(y[1:t] ~ X[1:t,]))
      RSS2 <- deviance(lm(y[(t+1):T] ~ X[(t+1):T,]))
      return(T*(RSS0 - (RSS1 + RSS2))/(RSS1 + RSS2))
    }
    return(max(sapply(floor(eps*T):(T-floor(eps*T)), wald)))
  }
  e = rnorm(T,0,1)
  emu = e + c(rep(0,floor(tstar*T)),rep(0.8,floor((1-tstar)*T)))
  y=filter(emu,0.5,"recursive")
  model1 = y[2:T] ~ 1 + y[1:(T-1)]
  H = matrix(c(1,0), ncol = 1)
  return(c(
    backCUSUM::BQ.test(model1,H=H)$statistic,
    backCUSUM::SBQ.test(model1,H=H)$statistic,
    backCUSUM::Q.test(model1,H=H)$statistic,
    sup.wald(model1)
  ))
}
##
sim.sizepow <- function(T,j,tstar){
  if(j == 1){
    Statistics <- parSapply(cl,rep(T,MC), simM1, tstar = tstar)
    crit = c(0.947, 1.202, 0.947, 8.85)
  } else if(j==2){
    Statistics <- parSapply(cl,rep(T,MC), simM2, tstar = tstar)
    crit = c(1.034, 1.274, 1.034, 11.79)
  } else {
    Statistics <- parSapply(cl,rep(T,MC), simM3, tstar = tstar)
    crit = c(0.947, 1.202, 0.947, 11.79)
  }
  rejections = function(j) (length(which(Statistics[j,] > crit[j]))/MC)
  size = sapply(1:length(crit),rejections)
  return(size)
}
##
taustars <- c(0, 0.1, 0.4, 0.6, 0.9)
M1.l = lapply(taustars, T = 200, j=1,  sim.sizepow)
M2.l = lapply(taustars, T = 200, j=2,  sim.sizepow)
M3.l = lapply(taustars, T = 200, j=3,  sim.sizepow)
##
rnames = c("size", taustars[-1])
cnames = c("BQ","SBQ","Q","supW")
##
M1 = matrix(unlist(M1.l), nrow = length(taustars), byrow=TRUE, dimnames = list(rnames, cnames))
M2 = matrix(unlist(M2.l), nrow = length(taustars), byrow=TRUE, dimnames = list(rnames, cnames))
M3 = matrix(unlist(M3.l), nrow = length(taustars), byrow=TRUE, dimnames = list(rnames, cnames))
##
table = round(cbind(M1,M2,M3)*100,1)
table
write.csv(table,file="./results/Table4.csv")
##
Sys.time()-start
stopCluster(cl)
