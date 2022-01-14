#' Forward CUSUM R-test (retrospective)
#'
#' Performs the multivariate forward CUSUM R-test (retrospective) with linear boundary
#'
#' @param formula Specification of the linear regression model by an object of the class "formula"
#' @param alternative A character string specifying the alternative hypothesis; must be one of "two.sided" (default), "greater" or "less".
#' The output detector is the maximum norm of \eqn{Q_t} ("two.sided"), the maximum entry of \eqn{Q_t} ("greater"), or maximum entry of \eqn{-Q_t} ("less"), respectively.
#' @param H An optional matrix for the partial hypothesis \eqn{H'\beta_t = H'\beta_0}, where \eqn{H'Q_t} is considered instead of \eqn{Q_t}.
#' \eqn{H} must have orthonormal columns. For a test for a break in the intercept, H can also set to the string "intercept".
#' The full structural break test is considered as the default setting (NULL).
#'
#' @return A list containing the following components:
#' \item{detector}{A vector containing the path of the detector statistic depending on the specificaton for the alternative hypothesis}
#' \item{boundary}{A vector containing the values of the linear boundary function}
#' \item{detector.scaled}{A vector containing the path of the detector divided by the boundary}
#' \item{statistic}{The test statistic; maximum of detector.scaled}
#' \item{alternative}{The specification for the alternative hypothesis}
#' \item{critical.value}{A vector containing critical values for different significance levels; NA if critical value for this specification is not implemented}
#' \item{rejection}{A logical vector containing the test decision for different significance levels; TRUE for rejection; NA if critical value is not implemented}
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
Q.test = function(formula, alternative = c("two.sided", "greater", "less"), H = NULL){
  alternative=match.arg(alternative)
  T <- dim(model.matrix(formula))[1]
  k <- dim(model.matrix(formula))[2]
  if (is.null(H)){
    Q <- get.cusumprocess(formula, T)
  } else {
    Q <- get.partialcusum(formula, T, H)
    k <- dim(Q)[1]
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
  output = list(
    detector = round(unname(detector),6),
    boundary = round(boundary,6),
    detector.scaled = round(unname(detector/boundary),6),
    statistic = round(statistic,6),
    alternative = alternative,
    critical.value = crit.val,
    rejection = rejection
  )
  return(output)
}




#' Backward CUSUM R-test (retrospective)
#'
#' Performs the multivariate backward CUSUM R-test (retrospective) with linear boundary
#'
#' @param formula Specification of the linear regression model by an object of the class "formula"
#' @param alternative A character string specifying the alternative hypothesis; must be one of "two.sided" (default), "greater" or "less".
#' The output detector is the maximum norm of \eqn{BQ_t} ("two.sided"), the maximum entry of \eqn{BQ_t} ("greater"), or maximum entry of \eqn{-BQ_t} ("less"), respectively.
#' @param H An optional matrix for the partial hypothesis \eqn{H'\beta_t = H'\beta_0}, where \eqn{H'BQ_t} is considered instead of \eqn{BQ_t}.
#' \eqn{H} must have orthonormal columns. For a test for a break in the intercept, H can also set to the string "intercept".
#' The full structural break test is considered as the default setting (NULL).
#'
#' @return A list containing the following components:
#' \item{detector}{A vector containing the path of the detector statistic depending on the specificaton for the alternative hypothesis}
#' \item{boundary}{A vector containing the values of the linear boundary function}
#' \item{detector.scaled}{A vector containing the path of the detector divided by the boundary}
#' \item{statistic}{The test statistic; maximum of detector.scaled}
#' \item{alternative}{The specification for the alternative hypothesis}
#' \item{critical.value}{A vector containing critical values for different significance levels; NA if critical value for this specification is not implemented}
#' \item{rejection}{A logical vector containing the test decision for different significance levels; TRUE for rejection; NA if critical value is not implemented}
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
BQ.test <- function(formula, alternative = c("two.sided", "greater", "less"), H = NULL){
  alternative=match.arg(alternative)
  T <- dim(model.matrix(formula))[1]
  k <- dim(model.matrix(formula))[2]
  if (is.null(H)){
    Q <- get.cusumprocess(formula, T)
  } else {
    Q <- get.partialcusum(formula, T, H)
    k <- dim(Q)[1]
  }
  BQ <- cbind(Q[,T],matrix(Q[,T] - Q[,1:(T-1)], nrow = k))
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
  output = list(
    detector = round(unname(detector),6),
    boundary = round(boundary,6),
    detector.scaled = round(unname(detector/boundary),6),
    statistic = round(statistic,6),
    alternative = alternative,
    critical.value = crit.val,
    rejection = rejection
  )
  return(output)
}



#' Stacked backward CUSUM R-test (retrospective)
#'
#' Performs the multivariate stacked backward CUSUM R-test (retrospective) with linear boundary
#'
#' @param formula Specification of the linear regression model by an object of the class "formula"
#' @param alternative A character string specifying the alternative hypothesis; must be one of "two.sided" (default), "greater" or "less".
#' The output detector matrix are the maximum norms of \eqn{SBQ_{s,t}} ("two.sided"), the maximum entries of \eqn{SBQ_{s,t}}("greater"), or maximum entries of \eqn{-SBQ_{s,t}} ("less"), respectively.
#' @param H An optional matrix for the partial hypothesis \eqn{H'\beta_t = H'\beta_0}, where \eqn{H'BQ_t} is considered instead of \eqn{SBQ_{s,t}}.
#' \eqn{H} must have orthonormal columns. For a test for a break in the intercept, H can also set to the string "intercept".
#' The full structural break test is considered as the default setting (NULL).
#'
#' @return A list containing the following components:
#' \item{detector.array}{A matrix containing the triangular array of the stacked backward cusum detector depending on the specificaton for the alternative hypothesis}
#' \item{boundary}{The matrix containing the values of the triangular boundary surface}
#' \item{detector.scaled}{A vector containing the path of the sequential stacked backward cusum detector divided by the boundary}
#' \item{statistic}{The test statistic; maximum of detector.scaled}
#' \item{alternative}{The specification for the alternative hypothesis}
#' \item{critical.value}{A vector containing critical values for different significance levels; NA if critical value for this specification is not implemented}
#' \item{rejection}{A logical vector containing the test decision for different significance levels; TRUE for rejection; NA if critical value is not implemented}
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
SBQ.test <- function(formula, alternative = c("two.sided", "greater", "less"), H = NULL){
  alternative=match.arg(alternative)
  T <- dim(model.matrix(formula))[1]
  k <- dim(model.matrix(formula))[2]
  if (is.null(H)){
    Q <- get.cusumprocess(formula, T)
  } else {
    Q <- get.partialcusum(formula, T, H)
    k <- dim(Q)[1]
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
  output = list(
    detector.array = round(detector.array,6),
    boundary = round(boundary,6),
    detector.scaled = round(unname(detector.scaled),6),
    statistic = round(statistic,6),
    alternative = alternative,
    critical.value = crit.val,
    rejection = rejection
  )
  return(output)
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
  return(max(sapply(floor(eps*T):(T-floor(eps*T)), wald)))
}
