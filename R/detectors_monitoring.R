#' Stacked backward CUSUM monitoring
#'
#' Performs the multivariate stacked backward CUSUM monitoring procedure using the linear boundary of Brown, Durbin, and Evans (1975)
#'
#' @param formula Specification of the linear regression model by an object of the class "formula"
#' @param T Length of the training sample. Monitoring starts at T+1.
#' @param m The length of the relative monitoring period m > 1; default is m = 10. Horizons larger than m=10 are not implemented.
#' @param alternative A character string specifying the alternative hypothesis; must be one of "two.sided" (default), "greater" or "less".
#' The output detector is the maximum norm of \eqn{Q_t} ("two.sided"), the maximum entry of \eqn{Q_t} ("greater"), or maximum entry of \eqn{-Q_t} ("less"), respectively.
#' @param H An optional matrix for the partial hypothesis \eqn{H'\beta_t = H'\beta_0}, where \eqn{H'Q_t} is considered instead of \eqn{Q_t}.
#' \eqn{H} must have orthonormal columns. For a test for a break in the intercept, H can also set to the string "intercept".
#' The full structural break test is considered as the default setting (NULL).
#'
#' @return A list containing the following components:
#' \item{detector.scaled}{A vector containing the path of the detector statistic scaled by the boundary function from T+1 onwards depending on the specificaton for the alternative hypothesis}
#' \item{detector.array}{A matrix containing the triangular array of the stacked backward cusum detector from T+1 onwards}
#' \item{boundary}{A matrix containing the values of the triangular boundary surface from T+1 onwards}
#' \item{statistic}{The test statistic; maximum of detector.scaled}
#' \item{detectontime}{The vector containing the detection time points for different significance levels, which are the time indices of the first boundary crossing; NA if the null hypothesis is not rejected}
#' \item{alternative}{The specification for the alternative hypothesis}
#' \item{critical.value}{A vector containing critical values for different significance levels; NA if critical value for this specification is not implemented}
#' \item{rejection}{A logical vector containing the test decision for different significance levels; TRUE for rejection; NA if critical value is not implemented}
#' @export
#'
#' @examples
#' T <- 100
#' t <- 5*T
#' u <- rnorm(t,0,1)
#' x <- rnorm(t,1,2)
#' y <- c(rep(0,480), rep(5,20)) + x + I(x^2) + u
#' SBQ.mon(y~1+x+I(x^2), T)
#' SBQ.mon(y~1+x+I(x^2), T, alternative = "greater")
#' H <- matrix(c(1,0,0), ncol = 1)
#' SBQ.mon(y~1+x+I(x^2), T, m=6, H = H)
SBQ.mon <- function(formula, T, m=10, alternative = c("two.sided", "greater", "less"), H = NULL){
  alternative=match.arg(alternative)
  if(m>10) stop("Critical values for m>10 are not implemented. Please specify a smaller monitoring horizon or use Q.mon.inf")
  n <- dim(model.matrix(formula))[1] #current time point
  if(n>m*T) stop(paste0("The maximum number of monitoring times has been reached. M=",m*T,", t=",n))
  k <- dim(model.matrix(formula))[2]
  if (is.null(H)){
    Q <- get.cusumprocess(formula, T)
  } else {
    Q <- get.partialcusum(formula, T, H)
    k <- dim(Q)[1]
  }
  SBQ <- array(NA, dim=c(n,n,k), dimnames = list(colnames(Q), colnames(Q), rownames(Q)))
  for (i in 1:k) (SBQ[T+1,(T+1):n, i] <- Q[i,(T+1):n] - Q[i,T])
  for (t in (T+2):n) (SBQ[(T+2):t, t, ] <- t(Q[, t] - Q[, (T+1):(t - 1)]))
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
    crit.val <- crit.SBQ.Mtest(k,m)$crit
  } else {
    crit.val <- crit.SBQ.Mtest(k, m, "one.sided")$crit
  }
  rejection <- statistic > crit.val
  # detection time point
  detection <- function(crit) ( T + which(detector.scaled > crit)[1] )
  detectiontime <- apply(matrix(crit.val), 1, detection)
  names(detectiontime) <- names(rejection)
  output = list(
    detector.scaled = round(unname(detector.scaled),6),
    detector.array = round(detector.array[-(1:T),-(1:T)],6),
    boundary = round(boundary[-(1:T),-(1:T)],6),
    statistic = round(statistic,6),
    detectiontime = detectiontime,
    alternative = alternative,
    critical.value = crit.val,
    rejection = rejection
  )
  return(output)
}


#' Multivariate forward CUSUM monitoring
#'
#' Performs the multivariate forward CUSUM monitoring procedure using the linear boundary of Brown, Durbin, and Evans (1975)
#'
#' @param formula Specification of the linear regression model by an object of the class "formula"
#' @param T Length of the training sample. Monitoring starts at T+1.
#' @param m The length of the relative monitoring period m > 1; default is m = 10. Horizons larger than m=10 are not implemented.
#' @param alternative A character string specifying the alternative hypothesis; must be one of "two.sided" (default), "greater" or "less".
#' The output detector is the maximum norm of \eqn{Q_t} ("two.sided"), the maximum entry of \eqn{Q_t} ("greater"), or maximum entry of \eqn{-Q_t} ("less"), respectively.
#' @param H An optional matrix for the partial hypothesis \eqn{H'\beta_t = H'\beta_0}, where \eqn{H'Q_t} is considered instead of \eqn{Q_t}.
#' \eqn{H} must have orthonormal columns. For a test for a break in the intercept, H can also set to the string "intercept".
#' The full structural break test is considered as the default setting (NULL).
#'
#' @return A list containing the following components:
#' \item{detector}{A vector containing the path of the detector statistic from T+1 onwards depending on the specificaton for the alternative hypothesis}
#' \item{boundary}{A vector containing the values of the linear boundary function from T+1 onwards}
#' \item{detector.scaled}{A vector containing the path of the detector divided by the boundary}
#' \item{statistic}{The test statistic; maximum of detector.scaled}
#' \item{detectontime}{The vector containing the detection time points for different significance levels, which are the time indices of the first boundary crossing; NA if the null hypothesis is not rejected}
#' \item{alternative}{The specification for the alternative hypothesis}
#' \item{critical.value}{A vector containing critical values for different significance levels; NA if critical value for this specification is not implemented}
#' \item{rejection}{A logical vector containing the test decision for different significance levels; TRUE for rejection; NA if critical value is not implemented}
#' @export
#'
#' @examples
#' T <- 100
#' t <- 5*T
#' u <- rnorm(t,0,1)
#' x <- rnorm(t,1,2)
#' y <- c(rep(0,480), rep(5,20)) + x + I(x^2) + u
#' Q.mon(y~1+x+I(x^2), T)
#' Q.mon(y~1+x+I(x^2), T, alternative = "greater")
#' H <- matrix(c(1,0,0), ncol = 1)
#' Q.mon(y~1+x+I(x^2), T, m=6, H = H)
Q.mon <- function(formula, T, m=10, alternative = c("two.sided", "greater", "less"), H = NULL){
  alternative=match.arg(alternative)
  if(m>10) stop("Critical values for m>10 are not implemented. Please specify a smaller monitoring horizon or use Q.mon.inf")
  n <- dim(model.matrix(formula))[1] #current time point
  if(n>m*T) stop(paste0("The maximum number of monitoring times has been reached. M=",m*T,", t=",n))
  k <- dim(model.matrix(formula))[2]
  if (is.null(H)){
    Q <- get.cusumprocess(formula, T)
  } else {
    Q <- get.partialcusum(formula, T, H)
    k <- dim(Q)[1]
  }
  # detector statistic
  m.detector <- Q[,(T+1):n,drop=F]-Q[,T]
  if(alternative == "two.sided")( detector <- apply(abs(m.detector), 2, max) )
  if(alternative == "less") ( detector <- apply(-m.detector, 2, max) )
  if(alternative == "greater") ( detector <- apply(m.detector, 2, max) )
  # boundary function
  boundary <- 1+2*(1:(n-T))/T
  # maximum statistic
  statistic <- max(detector/boundary)
  # critical values and test decision
  if(alternative == "two.sided"){
    crit.val <- crit.Q.Mtest(k, m)$crit
  } else {
    crit.val <- crit.Q.Mtest(k, m, "one.sided")$crit
  }
  rejection <- statistic > crit.val
  # detection time point
  detection <- function(crit) ( T + which(detector/boundary > crit)[1] )
  detectiontime <- apply(matrix(crit.val), 1, detection)
  names(detectiontime) <- names(rejection)
  output = list(
    detector = round(unname(detector),6),
    boundary = round(boundary,6),
    detector.scaled = round(unname(detector/boundary),6),
    statistic = round(statistic,6),
    detectiontime = detectiontime,
    alternative = alternative,
    critical.value = crit.val,
    rejection = rejection
  )
  return(output)
}
