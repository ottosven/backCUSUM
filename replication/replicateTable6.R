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
library(backCUSUM)
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
MC = 10
###############################
simM1 <- function(T, m=10, tstar=m){
  mT <- floor(m*T)
  fremdt.detector <- function(formula, T){
    X = model.matrix(formula)
    y.mod = eval(formula[[2]])
    n <- length(y.mod)  #end of monitoring
    trainingfit = lm(y.mod[1:T] ~ X[1:T] - 1)
    coef = trainingfit$coefficients
    olsresid = y.mod - X %*% coef
    olsQ = cumsum(olsresid)
    pageQ = matrix(NA, nrow = n, ncol = n)
    for(t in (T+1):n) (pageQ[(T+1):t,t] <- abs(olsQ[t] - olsQ[T:(t-1)]))
    detector = apply(pageQ[-(1:T),-(1:T)],2,max,na.rm=TRUE)
    detector/sd(trainingfit$residuals)
  }
  ##
  csw <- function(formula, T, alpha = 0.05, alternative = "two.sided"){
    n <- dim(model.matrix(formula))[1] #current time point
    k <- dim(model.matrix(formula))[2]
    H <- matrix(c(1,numeric(k-1)), ncol = 1)
    detector <- backCUSUM::Q.mon(formula, T, alternative = alternative, H = H)$detector
    # boundary function
    r <- (1:n)/T
    if(alternative == "two.sided"){
      boundary.CSW <- sqrt(r[(T+1):n]*(log(r[(T+1):n]/alpha^2)))
    } else {
      boundary.CSW <- sqrt(r[(T+1):n]*(log(r[(T+1):n]/(2*alpha)^2)))
    }
    CSW <- detector/boundary.CSW
    # maximum statistic
    statistic <- max(CSW)
    # critical values and test decision
    rejection <- statistic > 1
    # detection time point
    detectiontime <- unname(T + which(CSW > 1)[1])
    return(list(detector = round(unname(detector),6), boundary = round(boundary.CSW,6), rejection = rejection, detectiontime = detectiontime, statistic = round(statistic,6)))
  }
  ##
  sbq.detector <- function(model, T){
    boundary = matrix(NA, ncol = mT, nrow = mT)
    for(j in (T+1):mT){
      r=j/T
      S=(T:(j-1))/T
      boundary[(T+1):j,j] = (1+2*(r-S))*sqrt(r)
    }
    Q <- backCUSUM::get.cusumprocess(model, T)
    SBQ = matrix(NA, nrow = mT, ncol = mT)
    for(t in (T+1):mT) (SBQ[(T+1):t,t] <- abs(Q[t] - Q[T:(t-1)]))
    detector.array = SBQ/boundary
    SBQ = apply(detector.array[-(1:T),-(1:T)],2,max,na.rm=TRUE)
    return(SBQ)
  }
  ##
  e = rnorm(mT,0,1)
  y = c(rep(0,floor(tstar*T)),rep(0.8,floor((m-tstar)*T))) + e
  model <- y ~ 1
  ##
  SBQ = sbq.detector(model, T)
  Q.det <- backCUSUM::Q.mon(model, T)$detector
  r <- (1:mT)/T
  boundary.lin <- 1+2*(r-1)[(T+1):mT]
  Q <- Q.det/boundary.lin
  CSW.out = csw(model,T)
  CSW=CSW.out$detector/CSW.out$boundary
  fremdt = fremdt.detector(model, T)
  k.fr = 1:(length(y)-T)
  fremdt2 = fremdt/(sqrt(T)*(1+k.fr/T)*(k.fr/(k.fr+T))^0.25)
  result <- unname(c(
    which(SBQ > 0.976)[1] + T,
    which(Q > 0.958)[1] + T,
    which(CSW > 1)[1] + T,
    which(fremdt2 > 2.4296)[1] + T
  ))
  result
}
##
simM2 <- function(T, m=10, tstar=m){
  mT <- floor(m*T)
  fremdt.detector <- function(formula, T){
    X = model.matrix(formula)
    y.mod = eval(formula[[2]])
    n <- length(y.mod)  #end of monitoring
    trainingfit = lm(y.mod[1:T] ~ X[1:T,] - 1)
    coef = trainingfit$coefficients
    olsresid = y.mod - X %*% coef
    olsQ = cumsum(olsresid)
    pageQ = matrix(NA, nrow = n, ncol = n)
    for(t in (T+1):n) (pageQ[(T+1):t,t] <- abs(olsQ[t] - olsQ[T:(t-1)]))
    detector = apply(pageQ[-(1:T),-(1:T)],2,max,na.rm=TRUE)
    detector/sd(trainingfit$residuals)
  }
  ##
  csw <- function(formula, T, alpha = 0.05, alternative = "two.sided"){
    n <- dim(model.matrix(formula))[1] #current time point
    k <- dim(model.matrix(formula))[2]
    H <- matrix(c(1,numeric(k-1)), ncol = 1)
    detector <- backCUSUM::Q.mon(formula, T, alternative = alternative, H = H)$detector
    # boundary function
    r <- (1:n)/T
    if(alternative == "two.sided"){
      boundary.CSW <- sqrt(r[(T+1):n]*(log(r[(T+1):n]/alpha^2)))
    } else {
      boundary.CSW <- sqrt(r[(T+1):n]*(log(r[(T+1):n]/(2*alpha)^2)))
    }
    CSW <- detector/boundary.CSW
    # maximum statistic
    statistic <- max(CSW)
    # critical values and test decision
    rejection <- statistic > 1
    # detection time point
    detectiontime <- unname(T + which(CSW > 1)[1])
    return(list(detector = round(unname(detector),6), boundary = round(boundary.CSW,6), rejection = rejection, detectiontime = detectiontime, statistic = round(statistic,6)))
  }
  ##
  sbq.detector <- function(model, T){
    boundary = matrix(NA, ncol = mT, nrow = mT)
    for(j in (T+1):mT){
      r=j/T
      S=(T:(j-1))/T
      boundary[(T+1):j,j] = (1+2*(r-S))*sqrt(r)
    }
    Q <- backCUSUM::get.cusumprocess(model, T)
    SBQ = matrix(NA, nrow = mT, ncol = mT)
    for(t in (T+1):mT) (SBQ[(T+1):t,t] <- apply(abs(Q[,t] - Q[,T:(t-1),drop=F]),2,max))
    detector.array = SBQ/boundary
    SBQ = apply(detector.array[-(1:T),-(1:T)],2,max,na.rm=TRUE)
    return(SBQ)
  }
  ##
  e = rnorm(mT,0,1)
  x=filter(rnorm(mT,0,1),0.5,"recursive")
  y=1 + x*c(rep(0,floor(tstar*T)),rep(0.8,floor((m-tstar)*T)))+e
  model <- y ~ 1 + x
  ##
  SBQ = sbq.detector(model, T)
  Q.det <- backCUSUM::Q.mon(model, T)$detector
  r <- (1:mT)/T
  boundary.lin <- 1+2*(r-1)[(T+1):mT]
  Q <- Q.det/boundary.lin
  CSW.out = csw(model,T)
  CSW=CSW.out$detector/CSW.out$boundary
  fremdt = fremdt.detector(model, T)
  k.fr = 1:(length(y)-T)
  fremdt2 = fremdt/(sqrt(T)*(1+k.fr/T)*(k.fr/(k.fr+T))^0.25)
  result <- unname(c(
    which(SBQ > 1.036)[1] + T,
    which(Q > 1.044)[1] + T,
    which(CSW > 1)[1] + T,
    which(fremdt2 > 2.4296)[1] + T
  ))
  result
}
##
simM3 <- function(T, m=10, tstar=m){
  mT <- floor(m*T)
  fremdt.detector <- function(formula, T){
    X = model.matrix(formula)
    y.mod = eval(formula[[2]])
    n <- length(y.mod)  #end of monitoring
    trainingfit = lm(y.mod[1:T] ~ X[1:T,] - 1)
    coef = trainingfit$coefficients
    olsresid = y.mod - X %*% coef
    olsQ = cumsum(olsresid)
    pageQ = matrix(NA, nrow = n, ncol = n)
    for(t in (T+1):n) (pageQ[(T+1):t,t] <- abs(olsQ[t] - olsQ[T:(t-1)]))
    detector = apply(pageQ[-(1:T),-(1:T)],2,max,na.rm=TRUE)
    detector/sd(trainingfit$residuals)
  }
  ##
  csw <- function(formula, T, alpha = 0.05, alternative = "two.sided"){
    n <- dim(model.matrix(formula))[1] #current time point
    k <- dim(model.matrix(formula))[2]
    H <- matrix(c(1,numeric(k-1)), ncol = 1)
    detector <- backCUSUM::Q.mon(formula, T, alternative = alternative, H = H)$detector
    # boundary function
    r <- (1:n)/T
    if(alternative == "two.sided"){
      boundary.CSW <- sqrt(r[(T+1):n]*(log(r[(T+1):n]/alpha^2)))
    } else {
      boundary.CSW <- sqrt(r[(T+1):n]*(log(r[(T+1):n]/(2*alpha)^2)))
    }
    CSW <- detector/boundary.CSW
    # maximum statistic
    statistic <- max(CSW)
    # critical values and test decision
    rejection <- statistic > 1
    # detection time point
    detectiontime <- unname(T + which(CSW > 1)[1])
    return(list(detector = round(unname(detector),6), boundary = round(boundary.CSW,6), rejection = rejection, detectiontime = detectiontime, statistic = round(statistic,6)))
  }
  ##
  sbq.detector <- function(model, T){
    boundary = matrix(NA, ncol = mT-1, nrow = mT-1)
    for(j in (T+1):(mT-1)){
      r=j/T
      S=(T:(j-1))/T
      boundary[(T+1):j,j] = (1+2*(r-S))*sqrt(r)
    }
    Q <- backCUSUM::get.cusumprocess(model, T)
    Q=Q[1,]
    SBQ = matrix(NA, nrow = mT-1, ncol = mT-1)
    for(t in (T+1):(mT-1)) (SBQ[(T+1):t,t] <- abs(Q[t] - Q[T:(t-1)]))
    detector.array = SBQ/boundary
    SBQ = apply(detector.array[-(1:T),-(1:T)],2,max,na.rm=TRUE)
    return(SBQ)
  }
  ##
  e = rnorm(mT,0,1)
  emu = e + c(rep(0,floor(tstar*T)),rep(0.8,floor((m-tstar)*T)))
  y=filter(emu,0.5,"recursive")
  model = y[2:mT] ~ 1 + y[1:(mT-1)]
  H = matrix(c(1,0), ncol = 1)
  ##
  SBQ = sbq.detector(model, T)
  Q.det <- backCUSUM::Q.mon(model, T, H=H)$detector
  r <- (1:(mT-1))/T
  boundary.lin <- 1+2*(r-1)[(T+1):(mT-1)]
  Q <- Q.det/boundary.lin
  CSW.out = csw(model,T)
  CSW=CSW.out$detector/CSW.out$boundary
  fremdt = fremdt.detector(model, T)
  k.fr = 1:(length(y)-1-T)
  fremdt2 = fremdt/(sqrt(T)*(1+k.fr/T)*(k.fr/(k.fr+T))^0.25)
  result <- unname(c(
    which(SBQ > 0.976)[1] + T,
    which(Q > 0.958)[1] + T,
    which(CSW > 1)[1] + T,
    which(fremdt2 > 2.4296)[1] + T
  ))
  result
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
  rejectionrates = list()
  CDI = list() # correct detections indices
  delays = list()
  for(i in 1:4){
    rejectionrates[[i]] = sum(!is.na(Statistics[i,]))/MC
    CDI[[i]] = which(Statistics[i,] > (tstar*T))
    delays[[i]] = mean(Statistics[i,CDI[[i]]] - (tstar*T))
  }
  return(
    c(
      unlist(rejectionrates),
      unlist(delays)
    )
  )
}
##
m=10
taustars = c(m, 1.5, 2, 4, 6)
M1.l = lapply(taustars, T = 200, j=1, m=m, sim.sizedel)
M2.l = lapply(taustars, T = 200, j=2, m=m, sim.sizedel)
M3.l = lapply(taustars, T = 200, j=3, m=m, sim.sizedel)
##
rnames = c("size", taustars[-1])
cnames = rep(c("SBQ","Q","CSW","FR"),2)
##
M1 = matrix(unlist(M1.l), nrow = length(taustars), byrow=TRUE, dimnames = list(rnames,cnames))
M2 = matrix(unlist(M2.l), nrow = length(taustars), byrow=TRUE, dimnames = list(rnames,cnames))
M3 = matrix(unlist(M3.l), nrow = length(taustars), byrow=TRUE, dimnames = list(rnames,cnames))
M1
M2
M3
##
tab1 = rbind(M1[1,1:4]*100, M1[-1,5:8])
tab2 = rbind(M2[1,1:4]*100, M2[-1,5:8])
tab3 = rbind(M3[1,1:4]*100, M3[-1,5:8])
##
table=round(cbind(tab1, tab2, tab3),1)
rownames(table) = rnames
table
write.csv(table,file="./results/Table6.csv")
##
Sys.time()-start
stopCluster(cl)
