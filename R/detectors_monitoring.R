#' CUSUM monitoring detector
#'
#' @param formula Specification of the linear regression model by an object of the class "formula"
#' @param T Length of the training sample. Monitoring starts at T+1.
#' @param bound Type of the boundary. Either 'linear' or 'radical'
#'
#' @return A vector containing the Forward CUSUM monitoring detector series
#' @export
#'
#' @examples
#' T <- 100
#' u <- rnorm(T,0,1)
#' y <- c(rep(0,T/2), rep(1,T/2)) + u
#' Q.detector(y~1)
Qmon.detector <- function(formula, T, bound = 'linear'){
  wt.tail <- strucchange::recresid(formula)
  k <- lm(formula)$rank
  sigmahat <- sd(wt.tail[1:(T-k)])
  wt <- c(rep(0,k),wt.tail)
  r <- (1:length(wt))/T
  if(bound == 'linear'){
    boundary <- 1+2*(r-1)
    crit <- 0.957
  }
  if(bound == 'radical'){
    boundary <- sqrt(r*(log(r/0.05^2)))
    crit <- 1
  }
  Qt <- abs(cumsum(wt[-(1:T)]))/sqrt(T)/sigmahat/boundary[-(1:T)]
  detection <- which(Qt > crit)[1]
  return(list(Qt, detection, crit))
}

#' Stacked backward CUSUM detector
#'
#' @param formula Specification of the linear regression model by an object of the class "formula"
#' #' @param T Length of the training sample. Monitoring starts at T+1.
#'
#' @return A vector containing the Stacked backward CUSUM monitoring detector series
#' @export
#'
#' @examples
#' T <- 100
#' u <- rnorm(T,0,1)
#' y <- c(rep(0,T/2), rep(1,T/2)) + u
#' BQ.detector(y~1)
SBQmon.detector <- function(formula, T){
  wt.tail <- strucchange::recresid(formula)
  k <- lm(formula)$rank
  sigmahat <- sd(wt.tail[1:(T-k)])
  wt <- c(rep(0,k),wt.tail)
  stackedBackCUSUM <- function(t, T, sigmahat, wt){
    J <- ((t-T):1)/T
    boundSBQ <- 1+2*J
    max( abs(rev(cumsum(rev(wt[(T+1):t]))))/sqrt(T)/sigmahat/boundSBQ )
  }
  MaxSBQ <- sapply((T+1):length(wt), stackedBackCUSUM, T=T, sigmahat=sigmahat, wt = wt)
  crit <- 1.514
  detection <- which(MaxSBQ > crit)[1]
  return(list(MaxSBQ, detection, crit))
}

