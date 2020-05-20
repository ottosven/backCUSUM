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

#' Recursive residuals
#'
#' @param formula Specification of the linear regression model by an object of the class "formula"
#'
#' @return A vector containing the recursive residuals
#' @export
#'
#' @examples
#' T <- 100
#' u <- rnorm(T,0,1)
#' y <- c(rep(0,T/2), rep(1,T/2)) + u
#' get.recresid(y~1)
get.recresid <- function(formula){
  k <- lm(formula)$rank
  wt <- c(rep(0,k),strucchange::recresid(formula))
  return(wt)
}
