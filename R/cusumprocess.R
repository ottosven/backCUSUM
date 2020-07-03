#' @import stats
NULL


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


#' Multivariate CUSUM process
#'
#' @param formula Specification of the linear regression model by an object of the class "formula"
#' @param T In the retropspective context: length of the sample. In the monitoring context: length of the training sample, where monitoring starts at T+1.
#'
#' @return A matrix containing the multivariate forward CUSUM process Q_T(r)
#' @export
#'
#' @examples
#' T <- 100
#' u <- rnorm(T,0,1)
#' x <- rnorm(T,0,1)
#' y <- c(rep(0,T/2), rep(1,T/2)) + x + u
#' get.cusumprocess(y~x, T)
get.cusumprocess <- function(formula, T){
  wt <- get.recresid(formula)
  sig.hat <- sd(wt[1:T])
  X <- model.matrix(formula)
  sqCTinv <- expm::sqrtm(solve((t(X[1:T,]) %*% X[1:T,])/T))
  scores <- X * wt
  Q <- (sqCTinv %*% t(apply(scores, 2, cumsum)))/sig.hat/sqrt(T)
  rownames(Q) <- colnames(X)
  return(Q)
}



#' Critical values for the forward CUSUM monitoring procedure
#'
#' Provides critical values for some selected combinations of k and m for the forward CUSUM monitoring procedure with the linear boundary \eqn{d(r) = 1+2r} for different significance levels.
#'
#' @param k Number of regressors and dimension of the CUSUM process
#' @param m The length of the relative monitoring period; m = Inf for infinite horizon monitoring (default).
#' @param alternative A character string specifying the alternative hypothesis; must be one of "two.sided" (default) or "one.sided".
#'
#' @return A list containung the following components:
#' \item{critical.values}{A vector of critical values for different significance levels}
#' \item{m}{The value for the length of the relative monitoring period m. If the critical values for the input value of m is not available, the next higher available value is considered.}
#' @export
#'
#' @examples
#' get.crit.Q.mon(1,Inf)
#' get.crit.Q.mon(3,0.5)
get.crit.Q.mon <- function(k, m = Inf, alternative = "two.sided"){
  alphas <- as.numeric(colnames(Q.crit[[1]]))
  if(alternative == "one.sided") ( alphas <- alphas/2 )
  horizons <- as.numeric(rownames(Q.crit[[1]]))
  next.m <- which(horizons >= m)[1]
  if(k <= 8){
    critical <- Q.crit[[k]][next.m,]
  } else {
    critical <- rep(NA,length(alphas))
  }
  names(critical) <- alphas
  return(list(critical.values = critical, m = horizons[next.m]))
}



#' Critical values for the retrospective forward CUSUM test
#'
#' Provides critical values for some selected values for k for the forward CUSUM test with the linear boundary \eqn{d(r) = 1+2r} for different significance levels.
#'
#' @param k Number of regressors and dimension of the CUSUM process
#' @param alternative A character string specifying the alternative hypothesis; must be one of "two.sided" (default) or "one.sided".
#'
#' @return A vector of critical values for different significance levels
#' @export
#'
#' @examples
#' get.crit.Q(1)
get.crit.Q <- function(k, alternative = "two.sided"){
  return(get.crit.Q.mon(k, m=2,alternative)[[1]])
}



#' Critical values for the retrospective backward CUSUM test
#'
#' Provides critical values for some selected values for k for the backward CUSUM test with the linear boundary \eqn{d(r) = 1+2r} for different significance levels.
#'
#' @param k Number of regressors and dimension of the CUSUM process
#' @param alternative A character string specifying the alternative hypothesis; must be one of "two.sided" (default) or "one.sided".
#'
#' @return A vector of critical values for different significance levels
#' @export
#'
#' @examples
#' get.crit.BQ(1)
get.crit.BQ <- function(k, alternative = "two.sided"){
  return(get.crit.Q.mon(k, m=2, alternative)[[1]])
}



#' Critical values for the stacked backward CUSUM monitorng procedure
#'
#' Provides critical values for some selected combinations of k and m for the stacked backward CUSUM monitoring procedure with the linear boundary \eqn{d(r) = 1+2r} for different significance levels.
#'
#' @param k Number of regressors and dimension of the CUSUM process
#' @param m The length of the relative monitoring period; m = Inf for infinite horizon monitoring (default).
#' @param alternative A character string specifying the alternative hypothesis; must be one of "two.sided" (default) or "one.sided".
#'
#' @return A list containung the following components:
#' \item{critical.values}{A vector of critical values for different significance levels}
#' \item{m}{The value for the length of the relative monitoring period m. If the critical values for the input value of m is not available, the next higher available value is considered.}
#' @export
#'
#' @examples
#' get.crit.SBQ.mon(1,Inf)
#' get.crit.SBQ.mon(3,0.5)
get.crit.SBQ.mon <- function(k, m = Inf, alternative = "two.sided"){
  alphas <- as.numeric(colnames(SBQ.crit[[1]]))
  if(alternative == "one.sided") ( alphas <- alphas/2 )
  horizons <- as.numeric(rownames(SBQ.crit[[1]]))
  next.m <- which(horizons >= m)[1]
  if(k <= 8){
    critical <- SBQ.crit[[k]][next.m,]
  } else {
    critical <- rep(NA,length(alphas))
  }
  names(critical) <- alphas
  return(list(critical.values = critical, m = horizons[next.m]))
}



#' Critical values for the retrospective stacked backward CUSUM test
#'
#' Provides critical values for some selected values for k for the stacked backward CUSUM test with the linear boundary \eqn{d(r) = 1+2r} for different significance levels.
#'
#' @param k Number of regressors and dimension of the CUSUM process
#' @param alternative A character string specifying the alternative hypothesis; must be one of "two.sided" (default) or "one.sided".
#'
#' @return A vector of critical values for different significance levels
#' @export
#'
#' @examples
#' get.crit.SBQ(1)
get.crit.SBQ <- function(k, alternative = "two.sided"){
  return(get.crit.SBQ.mon(k, m=2, alternative)[[1]])
}
