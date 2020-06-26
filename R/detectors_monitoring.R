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
Qmon.detector <- function(formula, T, bound = 'linear', alternative = 'two.sided'){
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
    if(alternative == 'one.sided'){
      boundary <- sqrt(r*(log(r/0.1^2)))
    }
    crit <- 1
  }
  Qt <- cumsum(wt[-(1:T)])/sqrt(T)/sigmahat/boundary[-(1:T)]
  if(alternative != 'one.sided'){
    Qt <- abs(Qt)
  }
  detection <- which(Qt > crit)[1]
  return(list(Qt = Qt, recres = wt[-(1:T)], detectiontime = detection, criticalvalue = crit))
}

#' Stacked backward CUSUM detector
#'
#' @param formula Specification of the linear regression model by an object of the class "formula"
#' @param T Length of the training sample. Monitoring starts at T+1.
#'
#' @return A vector containing the Stacked backward CUSUM monitoring detector series
#' @export
#'
#' @examples
#' T <- 100
#' u <- rnorm(T,0,1)
#' y <- c(rep(0,T/2), rep(1,T/2)) + u
#' BQ.detector(y~1)
#'
SBQmon.detector <- function(formula, T, alternative = 'two.sided'){
  wt.tail <- strucchange::recresid(formula)
  k <- lm(formula)$rank
  sigmahat <- sd(wt.tail[1:(T-k)])
  wt <- c(rep(0,k),wt.tail)
  stackedBackCUSUM <- function(t, T, sigmahat, wt){
    J <- ((t-T):1)/T
    boundSBQ <- 1+2*J
    SBQ <- rev(cumsum(rev(wt[(T+1):t])))/sqrt(T)/sigmahat/boundSBQ
    if(alternative == 'two.sided'){
      SBQ <- abs(SBQ)
    }
    max( SBQ )
  }
  MaxSBQ <- sapply((T+1):length(wt), stackedBackCUSUM, T=T, sigmahat=sigmahat, wt = wt)
  crit <- 1.514
  if(alternative == 'one.sided'){
    crit <- 1.450
  }
  detection <- which(MaxSBQ > crit)[1]
  return(list(MaxSBQ, detection, crit))
}


#' Forward CUSUM monitoring procedure
#'
#' Performs multivariate forward CUSUM monitoring using the linear boundary of Brown, Durbin, and Evans (1975)
#'
#' @param formula Specification of the linear regression model by an object of the class "formula"
#' @param T Length of the training sample. Monitoring starts at T+1.
#' @param m The length of the relative monitoring period; default is m = Inf (infinite horizon monitoring).
#' @param alternative A character string specifying the alternative hypothesis; must be one of "two.sided" (default), "greater" or "less".
#' The detector statistic is given by the maximum norm of \eqn{Q_t} ("two.sided"), maximum entry of \eqn{Q_t} ("greater"), or maximum entry of \eqn{-Q_t} ("less"), respectively.
#' @param H An optional matrix for the partial hypothesis \eqn{H'\beta_t = H'\beta_0}, where \eqn{H'Q_t} is considered instead of \eqn{Q_t}.
#' \eqn{H} must have orthonormal columns. The full structural break test is considered as the default setting (NULL).
#'
#' @return A list containung the following components:
#' \item{statistic}{The test statistic; maximum of the detector scaled by its boundary \eqn{d(r)}}
#' \item{critical.value}{A vector containing critical values for different significance levels; NA if critical value for this specification is not implemented}
#' \item{rejection}{A logical vector containing the test decision for different significance levels; TRUE for rejection; NA if critical value is not implemented}
#' \item{detector}{The vector containing the path of the forward cusum detector}
#' \item{detectontime}{The vector containing the detection time points for different significance levels, which are the time indices of the first boundary crossing; NA if the null hypothesis is not rejected}
#' \item{boundary}{The vector containing the values of the boundary function}
#' @export
#'
#' @examples
#' T <- 100
#' t <- 5*T
#' u <- rnorm(t,0,1)
#' x <- rnorm(t,1,2)
#' y <- c(rep(0,480), rep(10,20)) + x + I(x^2) + u
#' Q.monitoring(y~1+x+I(x^2), T)
#' Q.monitoring(y~1+x+I(x^2), T, alternative = "greater")
#' H <- matrix(c(1,0,0), ncol = 1)
#' Q.monitoring(y~1+x+I(x^2), T, m=6, H = H)
Q.monitoring <- function(formula, T, m=Inf, alternative = "two.sided", H = NULL){
  t <- dim(model.matrix(formula))[1]
  k <- dim(model.matrix(formula))[2]
  Q <- get.cusumprocess(formula, T)
  # in case of a partial structural break test modify the process
  if (!is.null(H)){
    Q <- t(H) %*% Q
    k <- dim(H)[2]
  }
  # detector statistic
  if(alternative == "two.sided")( detector <- apply(abs(Q[,(T+1):t,drop=F]-Q[,T]), 2, max) )
  if(alternative == "less") ( detector <- apply(-Q[,(T+1):t,drop=F]-Q[,T], 2, max) )
  if(alternative == "greater") ( detector <- apply(Q[,(T+1):t,drop=F]-Q[,T], 2, max) )
  # boundary function
  boundary <- 1+2*((T+1):t)/T
  # maximum statistic
  statistic <- max(detector/boundary)
  # critical values and test decision
  if(alternative == "two.sided"){
    crit.val <- get.crit.Q.mon(k, m)$critical.values
  } else {
    crit.val <- get.crit.Q.mon(k, m, "one.sided")$critical.values
  }
  rejection <- statistic > crit.val
  detection <- function(crit) ( T + which(detector/boundary > crit)[1] )
  detectiontime <- apply(matrix(crit.val), 1, detection)
  names(detectiontime) <- names(rejection)
  detectiontime
  return(list(detector = round(detector,6), boundary = round(boundary,6), critical.value = crit.val, rejection = rejection, detectiontime = detectiontime, statistic = round(statistic,6)))
}
