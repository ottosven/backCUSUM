#' Forward CUSUM monitoring
#'
#' Performs the multivariate forward CUSUM monitoring procedure using the linear boundary of Brown, Durbin, and Evans (1975)
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
#' y <- c(rep(0,480), rep(5,20)) + x + I(x^2) + u
#' Q.monitoring(y~1+x+I(x^2), T)
#' Q.monitoring(y~1+x+I(x^2), T, alternative = "greater")
#' H <- matrix(c(1,0,0), ncol = 1)
#' Q.monitoring(y~1+x+I(x^2), T, m=6, H = H)
Q.monitoring <- function(formula, T, m=Inf, alternative = "two.sided", H = NULL){
  n <- dim(model.matrix(formula))[1] #current time point
  k <- dim(model.matrix(formula))[2]
  Q <- get.cusumprocess(formula, T)
  # in case of a partial structural break test modify the process
  if (!is.null(H)){
    Q <- t(H) %*% Q
    k <- dim(H)[2]
  }
  # detector statistic
  if(alternative == "two.sided")( detector <- apply(abs(Q[,(T+1):n,drop=F]-Q[,T]), 2, max) )
  if(alternative == "less") ( detector <- apply(-Q[,(T+1):n,drop=F]-Q[,T], 2, max) )
  if(alternative == "greater") ( detector <- apply(Q[,(T+1):n,drop=F]-Q[,T], 2, max) )
  # boundary function
  boundary <- 1+2*(1:(n-T))/T
  # maximum statistic
  statistic <- max(detector/boundary)
  # critical values and test decision
  if(alternative == "two.sided"){
    crit.val <- get.crit.Q.mon(k, m)$critical.values
  } else {
    crit.val <- get.crit.Q.mon(k, m, "one.sided")$critical.values
  }
  rejection <- statistic > crit.val
  # detection time point
  detection <- function(crit) ( T + which(detector/boundary > crit)[1] )
  detectiontime <- apply(matrix(crit.val), 1, detection)
  names(detectiontime) <- names(rejection)
  return(list(detector = round(detector,6), boundary = round(boundary,6), critical.value = crit.val, rejection = rejection, detectiontime = detectiontime, statistic = round(statistic,6)))
}



#' Stacked backward CUSUM monitoring
#'
#' Performs the multivariate stacked backward CUSUM monitoring procedure using the linear boundary of Brown, Durbin, and Evans (1975)
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
#' \item{detectontime}{The vector containing the detection time points for different significance levels, which are the time indices of the first boundary crossing; NA if the null hypothesis is not rejected}
#' \item{detector.array}{The matrix containing the triangular array of the stacked backward cusum detector}
#' \item{detector.scaled}{The vector containing the path of the sequential scaled stacked backward cusum detector}
#' \item{boundary}{The matrix containing the values of the triangular boundary surface}
#' @export
#'
#' @examples
#' T <- 100
#' t <- 5*T
#' u <- rnorm(t,0,1)
#' x <- rnorm(t,1,2)
#' y <- c(rep(0,480), rep(5,20)) + x + I(x^2) + u
#' SBQ.monitoring(y~1+x+I(x^2), T)
#' SBQ.monitoring(y~1+x+I(x^2), T, alternative = "greater")
#' H <- matrix(c(1,0,0), ncol = 1)
#' SBQ.monitoring(y~1+x+I(x^2), T, m=6, H = H)
SBQ.monitoring <- function(formula, T, m=Inf, alternative = "two.sided", H = NULL){
  n <- dim(model.matrix(formula))[1]  #current time point
  k <- dim(model.matrix(formula))[2]
  Q <- get.cusumprocess(formula, T)
  # in case of a partial structural break test modify the process
  if (!is.null(H)){
    Q <- t(H) %*% Q
    k <- dim(H)[2]
  }
  SBQ <- array(NA, dim=c(n,n,k), dimnames = list(colnames(Q), colnames(Q), rownames(Q)))
  for(t in (T+1):n) (SBQ[(T+1):t,t,] <- t(Q[,t] - Q[,T:(t-1)]))
  # detector statistic
  if(alternative == "two.sided")( detector.array <- apply(abs(SBQ), c(1,2), max) )
  if(alternative == "less") ( detector.array <- apply(-SBQ, c(1,2), max) )
  if(alternative == "greater") ( detector.array <- apply(SBQ, c(1,2), max) )
  # boundary function
  boundary <- matrix(NA, ncol = n, nrow = n, dimnames = list(colnames(Q), colnames(Q)))
  for(j in (T+1):n) ( boundary[(T+1):j,j] <- 1+2*(j-((T+1):j)+1)/T )
  # maximum statistic
  detector.scaled <- apply(detector.array[-(1:T),-(1:T)]/boundary[-(1:T),-(1:T)], 2, max, na.rm = TRUE)
  statistic <- max(detector.scaled)
  # critical values and test decision
  if(alternative == "two.sided"){
    crit.val <- get.crit.SBQ.mon(k,m)$critical.values
  } else {
    crit.val <- get.crit.SBQ.mon(k, m, "one.sided")$critical.values
  }
  rejection <- statistic > crit.val
  # detection time point
  detection <- function(crit) ( T + which(detector.scaled > crit)[1] )
  detectiontime <- apply(matrix(crit.val), 1, detection)
  names(detectiontime) <- names(rejection)
  return(list(
    detector.scaled = round(detector.scaled,6),
    detector.array = round(detector.array[-(1:T),-(1:T)],6),
    boundary = round(boundary[-(1:T),-(1:T)],6),
    critical.value = crit.val,
    rejection = rejection,
    detectiontime = detectiontime,
    statistic = round(statistic,6)
  ))
}
