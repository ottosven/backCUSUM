#' CUSUM detector
#'
#' @param formula Specification of the linear regression model by an object of the class "formula"
#'
#' @return A vector containing the Forward CUSUM detector series
#' @export
#'
#' @examples
#' T <- 100
#' u <- rnorm(T,0,1)
#' y <- c(rep(0,T/2), rep(1,T/2)) + u
#' Q.detector(y~1)
Q.detector <- function(formula){
  wt.tail <- strucchange::recresid(formula)
  sigmahat <- sd(wt.tail)
  k <- lm(formula)$rank
  wt <- c(rep(0,k),wt.tail)
  Qt <- abs(cumsum(wt))/sqrt(length(wt))/sigmahat
  return(Qt)
}

#' Backward CUSUM detector
#'
#' @param formula Specification of the linear regression model by an object of the class "formula"
#'
#' @return A vector containing the Backward CUSUM detector series
#' @export
#'
#' @examples
#' T <- 100
#' u <- rnorm(T,0,1)
#' y <- c(rep(0,T/2), rep(1,T/2)) + u
#' BQ.detector(y~1)
BQ.detector <- function(formula){
  wt.tail <- strucchange::recresid(formula)
  sigmahat <- sd(wt.tail)
  k <- lm(formula)$rank
  wt <- c(rep(0,k),wt.tail)
  BQt <- abs(sum(wt) - c(0,cumsum(wt)[-length(wt)]))/sqrt(length(wt))/sigmahat
  return(BQt)
}


#' Retrospective forward CUSUM test
#'
#' Performs the multivariate forward CUSUM test in the retrospective contect using the linear boundary of Brown, Durbin, and Evans (1975)
#'
#' @param formula Specification of the linear regression model by an object of the class "formula"
#' @param alternative A character string specifying the alternative hypothesis; must be one of "two.sided" (default), "greater" or "less".
#' The detector statistic is given by the maximum norm of \eqn{Q_t} ("two.sided"), maximum entry of \eqn{Q_t} ("greater"), or maximum entry of \eqn{-Q_t} ("less"), respectively.
#' @param H An optional matrix for the partial hypothesis \eqn{H'\beta_t = H'\beta_0}, where \eqn{H'Q_t} is considered instead of \eqn{Q_t}.
#' \eqn{H} must have orthonormal columns. The full structural break test is considered as the default setting (NULL).
#' @param bound An optional vector that contains the values of a user specified boundary function. The linear boundary \eqn{d(r) = 1+2r} is considered as the default setting (NULL).
#'
#' @return A list containung the following components:
#' \item{statistic}{The test statistic; maximum of the detector scaled by its boundary \eqn{d(r)}}
#' \item{critical.value}{A vector containing critical values for different significance levels; NA if critical value for this specification is not implemented}
#' \item{rejection}{A logical vector containing the test decision for different significance levels; TRUE for rejection; NA if critical value is not implemented}
#' \item{detector}{The vector containing the path of the forward cusum detector}
#' \item{boundary}{The vector containing the values of the boundary function}
#' @export
#'
#' @examples
#' T <- 100
#' u <- rnorm(T,0,1)
#' x <- rnorm(T,1,2)
#' y <- c(rep(0,T/2), rep(0.7,T/2)) + x + I(x^2) + u
#' Q.test(y~1+x+I(x^2))
#' Q.test(y~1+x+I(x^2), alternative = "greater")
#' H <- matrix(c(1,0,0), ncol = 1)
#' Q.test(y~1+x+I(x^2), H = H)
Q.test <- function(formula, alternative = "two.sided", H = NULL, bound = NULL){
  T <- dim(model.matrix(formula))[1]
  k <- dim(model.matrix(formula))[2]
  Q <- get.cusumprocess(formula, T)
  # in case of a partial structural break test modify the process
  if (!is.null(H)){
    Q <- t(H) %*% Q
    k <- dim(H)[2]
  }
  # detector statistic
  if(alternative == "two.sided")( detector <- apply(abs(Q), 2, max) )
  if(alternative == "less") ( detector <- apply(-Q, 2, max) )
  if(alternative == "greater") ( detector <- apply(Q, 2, max) )
  # boundary function
  boundary <- 1+2*(1:T)/T
  # maximum statistic
  statistic <- max(detector/boundary)
  # critical values and test decision
  if(alternative == "two.sided"){
    crit.val <- get.crit.Q(k)
  } else {
    crit.val <- get.crit.Q(k, "one.sided")
  }
  rejection <- statistic > crit.val
  # user specified boundary function
  if(!is.null(bound)){
    boundary <- bound
    crit.val <- NA
    rejection <- NA
    statistic <- max(detector/boundary)
  }
  return(list(detector = round(detector,6), boundary = round(boundary,6), critical.value = crit.val, rejection = rejection, statistic = round(statistic,6)))
}




#' Retrospective backward CUSUM test
#'
#' Performs the multivariate backward CUSUM test in the retrospective contect using the linear boundary \eqn{d(r) = 1+2r}.
#'
#' @param formula Specification of the linear regression model by an object of the class "formula"
#' @param alternative A character string specifying the alternative hypothesis; must be one of "two.sided" (default), "greater" or "less".
#' The detector statistic is given by the maximum norm of \eqn{Q_t} ("two.sided"), maximum entry of \eqn{Q_t} ("greater"), or maximum entry of \eqn{-Q_t} ("less"), respectively.
#' @param H An optional matrix for the partial hypothesis \eqn{H'\beta_t = H'\beta_0}, where \eqn{H'Q_t} is considered instead of \eqn{Q_t}.
#' \eqn{H} must have orthonormal columns. The full structural break test is considered as the default setting (NULL).
#' @param bound An optional vector that contains the values of a user specified boundary function. The linear boundary \eqn{d(r) = 1+2r} is considered as the default setting (NULL).
#'
#' @return A list containung the following components:
#' \item{statistic}{The test statistic; maximum of the detector scaled by its boundary \eqn{d(r)}}
#' \item{critical.value}{A vector containing critical values for different significance levels; NA if critical value for this specification is not implemented}
#' \item{rejection}{A logical vector containing the test decision for different significance levels; TRUE for rejection; NA if critical value is not implemented}
#' \item{detector}{The vector containing the path of the backward cusum detector}
#' \item{boundary}{The vector containing the values of the boundary function}
#' @export
#'
#' @examples
#' T <- 100
#' u <- rnorm(T,0,1)
#' x <- rnorm(T,1,2)
#' y <- c(rep(0,T/2), rep(0.7,T/2)) + x + I(x^2) + u
#' BQ.test(y~1+x+I(x^2))
#' BQ.test(y~1+x+I(x^2), alternative = "greater")
#' H <- matrix(c(1,0,0), ncol = 1)
#' BQ.test(y~1+x+I(x^2), H = H)
BQ.test <- function(formula, alternative = "two.sided", H = NULL, bound = NULL){
  T <- dim(model.matrix(formula))[1]
  k <- dim(model.matrix(formula))[2]
  Q <- get.cusumprocess(formula, T)
  # in case of a partial structural break test modify the process
  if (!is.null(H)){
    Q <- t(H) %*% Q
    k <- dim(H)[2]
  }
  BQ <- cbind(Q[,T],Q[,T] - Q[,1:(T-1)])
  # detector statistic
  if(alternative == "two.sided")( detector <- apply(abs(BQ), 2, max) )
  if(alternative == "less") ( detector <- apply(-BQ, 2, max) )
  if(alternative == "greater") ( detector <- apply(BQ, 2, max) )
  # boundary function
  boundary <- rev(1+2*(1:T)/T)
  # maximum statistic
  statistic <- max(detector/boundary)
  # critical values and test decision
  if(alternative == "two.sided"){
    crit.val <- get.crit.BQ(k)
  } else {
    crit.val <- get.crit.BQ(k, "one.sided")
  }
  rejection <- statistic > crit.val
  # user specified boundary function
  if(!is.null(bound)){
    boundary <- rev(bound)
    crit.val <- NA
    rejection <- NA
    statistic <- max(detector/boundary)
  }
  return(list(detector = round(detector,6), boundary = round(boundary,6), critical.value = crit.val, rejection = rejection, statistic = round(statistic,6)))
}



#' Retrospective stacked backward CUSUM test
#'
#' Performs the multivariate stacked backward CUSUM test in the retrospective contect using the linear boundary \eqn{d(r) = 1+2r}.
#'
#' @param formula Specification of the linear regression model by an object of the class "formula"
#' @param alternative A character string specifying the alternative hypothesis; must be one of "two.sided" (default), "greater" or "less".
#' The detector statistic is given by the maximum norm of \eqn{Q_t} ("two.sided"), maximum entry of \eqn{Q_t} ("greater"), or maximum entry of \eqn{-Q_t} ("less"), respectively.
#' @param H An optional matrix for the partial hypothesis \eqn{H'\beta_t = H'\beta_0}, where \eqn{H'Q_t} is considered instead of \eqn{Q_t}.
#' \eqn{H} must have orthonormal columns. The full structural break test is considered as the default setting (NULL).
#' @param bound An optional vector that contains the values of a user specified boundary function. The linear boundary \eqn{d(r) = 1+2r} is considered as the default setting (NULL).
#'
#' @return A list containung the following components:
#' \item{statistic}{The test statistic; maximum of the detector scaled by its boundary \eqn{d(r)}}
#' \item{critical.value}{A vector containing critical values for different significance levels; NA if critical value for this specification is not implemented}
#' \item{rejection}{A logical vector containing the test decision for different significance levels; TRUE for rejection; NA if critical value is not implemented}
#' \item{detector.array}{The matrix containing the triangular array of the stacked backward cusum detector}
#' \item{detector.scaled}{The vector containing the path of the sequential scaled stacked backward cusum detector}
#' \item{boundary}{The matrix containing the values of the triangular boundary surface}
#' @export
#'
#' @examples
#' T <- 100
#' u <- rnorm(T,0,1)
#' x <- rnorm(T,1,2)
#' y <- c(rep(0,T/2), rep(0.7,T/2)) + x + I(x^2) + u
#' SBQ.test(y~1+x+I(x^2))
#' SBQ.test(y~1+x+I(x^2), alternative = "greater")
#' H <- matrix(c(1,0,0), ncol = 1)
#' SBQ.test(y~1+x+I(x^2), H = H)
SBQ.test <- function(formula, alternative = "two.sided", H = NULL, bound = NULL){
  T <- dim(model.matrix(formula))[1]
  k <- dim(model.matrix(formula))[2]
  Q <- get.cusumprocess(formula, T)
  # in case of a partial structural break test modify the process
  if (!is.null(H)){
    Q <- t(H) %*% Q
    k <- dim(H)[2]
  }
  SBQ <- array(NA, dim=c(T,T,k), dimnames = list(colnames(Q), colnames(Q), rownames(Q)))
  for(i in 1:k) (SBQ[1,,i] <- Q[i,])
  for(t in 2:T) (SBQ[2:t,t,] <- t(Q[,t] - Q[,1:(t-1)]))
  # detector statistic
  if(alternative == "two.sided")( detector.array <- apply(abs(SBQ), c(1,2), max) )
  if(alternative == "less") ( detector.array <- apply(-SBQ, c(1,2), max) )
  if(alternative == "greater") ( detector.array <- apply(SBQ, c(1,2), max) )
  # boundary function
  boundary <- matrix(NA, ncol = T, nrow = T, dimnames = list(colnames(Q), colnames(Q)))
  for(t in 1:T) ( boundary[1:t,t] <- 1+2*(t-(1:t)+1)/T )
  # maximum statistic
  detector.scaled <- apply(detector.array/boundary, 2, max, na.rm = TRUE)
  statistic <- max(detector.scaled)
  # critical values and test decision
  if(alternative == "two.sided"){
    crit.val <- get.crit.SBQ(k)
  } else {
    crit.val <- get.crit.SBQ(k, "one.sided")
  }
  rejection <- statistic > crit.val
  # user specified boundary function
  if(!is.null(bound)){
    boundary <- matrix(NA, ncol = T, nrow = T, dimnames = list(colnames(Q), colnames(Q)))
    for(t in 1:T) ( boundary[1:t,t] <- bound[t-(1:t)+1] )
    crit.val <- NA
    rejection <- NA
    detector.scaled <- apply(detector.array/boundary, 2, max, na.rm = TRUE)
    statistic <- max(detector.scaled)
  }
  return(list(detector.scaled = round(detector.scaled,6), detector.array = round(detector.array,6), boundary = round(boundary,6), critical.value = crit.val, rejection = rejection, statistic = round(statistic,6)))
}



#' Sup-Wald Test by Andrews (1993)
#'
#' Performs the Sup-Wald Test by Andrews (1993).
#'
#' @param formula Specification of the linear regression model by an object of the class "formula"
#' @param eps Trimming parameter, eps = 0.15 is default.
#'
#' @return The sup-Wald test statistic
#' @export
#'
#' @examples
#' T <- 100
#' u <- rnorm(T,0,1)
#' x <- rnorm(T,1,2)
#' y <- c(rep(0,T/2), rep(0.7,T/2)) + x + I(x^2) + u
#' sup.wald(y~1+x+I(x^2))
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
  sapply(floor(eps*T):(T-floor(eps*T)), wald)
  return(max(sapply(floor(eps*T):(T-floor(eps*T)), wald)))
}
