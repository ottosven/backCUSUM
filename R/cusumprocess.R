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


#' Multivariate CUSUM process \eqn{Q_T(r)}
#'
#' @param formula Specification of the linear regression model by an object of the class "formula"
#' @param T Optional length of the training sample for monitoring. For retrospective tests, T should be the sample size. If NULL, T is stet to the sample size. NULL is default.
#'
#' @return A matrix containing the multivariate forward CUSUM process \eqn{Q_T(r)}
#' @export
#'
#' @examples
#' T <- 100
#' u <- rnorm(T,0,1)
#' x <- rnorm(T,0,1)
#' y <- c(rep(0,T/2), rep(1,T/2)) + x + u
#' get.cusumprocess(y~x, T)
get.cusumprocess <- function(formula, T=NULL){
  if(is.null(T)) T = dim(model.frame(formula))[1]
  wt <- get.recresid(formula)
  sig.hat <- sd(wt[1:T])
  X <- model.matrix(formula)
  sqCTinv <- expm::sqrtm(solve((t(X[1:T,]) %*% X[1:T,])/T))
  scores <- X * wt
  Q <- (sqCTinv %*% t(apply(scores, 2, cumsum)))/sig.hat/sqrt(T)
  rownames(Q) <- colnames(X)
  return(Q)
}



#' Partial CUSUM process \eqn{Q_T^*(r)}
#'
#' @param formula Specification of the linear regression model by an object of the class "formula"
#' @param T Optional length of the training sample for monitoring. For retrospective tests, T should be the sample size. If NULL, T is stet to the sample size. NULL is default.
#' @param H An optional matrix for the partial hypothesis \eqn{H'\beta_t = H'\beta_0}. \eqn{H} must have orthonormal columns.
#' If "intercept", H is set to the setting of testing for a break only in the first regressor variable (typically the intercept). "intercept" is default.
#'
#' @return A matrix containing the partial forward CUSUM process \eqn{Q_T^*(r)}
#' @export
#'
#' @examples
#' T <- 100
#' u <- rnorm(T,0,1)
#' x <- rnorm(T,0,1)
#' y <- c(rep(0,T/2), rep(1,T/2)) + x + u
#' H <- matrix(c(1,0), ncol = 1)
#' get.partialcusum(y~x, T, H=H)
get.partialcusum <- function(formula, T = NULL, H = "intercept"){
  if(is.null(T)) T = dim(model.frame(formula))[1]
  if(any(H == "intercept")) H = matrix(c(1,rep(0,dim(model.frame(formula))[2]-1)), ncol=1)
  wt <- get.recresid(formula)
  sig.hat <- sd(wt[1:T])
  X <- model.matrix(formula)
  partialregressors <- t(H) %*% t(X)
  sqCTinv <- expm::sqrtm(solve((partialregressors[,1:T, drop=F] %*% t(partialregressors[,1:T, drop=F]))/T))
  scores <- partialregressors * wt
  Q <- (sqCTinv %*% t(apply(scores, 1, cumsum)))/sig.hat/sqrt(T)
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
#' @return A list containing the following components:
#' \item{critical.values}{A vector of critical values for different significance levels}
#' \item{m}{The value for the length of the relative monitoring period m. If the critical values for the input value of m is not available, the next higher available value is considered.}
#' @export
#'
#' @examples
#' crit.Q.Mtest(1,Inf)
#' crit.Q.Mtest(3,0.5)
crit.Q.Mtest <- function(k, m = Inf, alternative = "two.sided"){
  if(alternative == "one.sided"){
    index <- c(1,2,4,6)
  } else {
    index <- c(2,3,5,7)
  }
  alphas <- c(0.1, 0.05, 0.01, 0.001)
  horizons <- as.numeric(rownames(Q.crit[[1]]))
  next.m <- which(horizons >= m)[1]
  if(k <= 30){
    critical <- Q.crit[[k]][next.m,index]
  } else {
    critical <- rep(NA,length(alphas))
  }
  names(critical) <- alphas
  return(list(crit = critical, m = horizons[next.m]))
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
#' crit.Q.Rtest(1)
crit.Q.Rtest <- function(k, alternative = "two.sided"){
  return(crit.Q.Mtest(k, m=2,alternative)[[1]])
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
#' crit.BQ(1)
crit.BQ <- function(k, alternative = "two.sided"){
  return(crit.Q.Rtest(k, alternative)[[1]])
}



#' Critical values for the stacked backward CUSUM monitorng procedure
#'
#' Provides critical values for some selected combinations of k and m for the stacked backward CUSUM monitoring procedure with the linear boundary \eqn{d(r) = 1+2r} for different significance levels.
#'
#' @param k Number of regressors and dimension of the CUSUM process
#' @param m The length of the relative monitoring period; m = Inf for infinite horizon monitoring (default).
#' @param alternative A character string specifying the alternative hypothesis; must be one of "two.sided" (default) or "one.sided".
#'
#' @return A list containing the following components:
#' \item{critical.values}{A vector of critical values for different significance levels}
#' \item{m}{The value for the length of the relative monitoring period m. If the critical values for the input value of m is not available, the next higher available value is considered.}
#' @export
#'
#' @examples
#' crit.SBQ.Mtest(1,Inf)
#' crit.SBQ.Mtest(3,0.5)
crit.SBQ.Mtest <- function(k, m = Inf, alternative = "two.sided"){
  if(alternative == "one.sided"){
    index <- c(1,2,4,6)
  } else {
    index <- c(2,3,5,7)
  }
  alphas <- c(0.1, 0.05, 0.01, 0.001)
  horizons <- as.numeric(rownames(SBQ.crit[[1]]))
  next.m <- which(horizons >= m)[1]
  if(k <= 30){
    critical <- SBQ.crit[[k]][next.m,index]
  } else {
    critical <- rep(NA,length(alphas))
  }
  names(critical) <- alphas
  return(list(crit = critical, m = horizons[next.m]))
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
#' crit.SBQ.Rtest(1)
crit.SBQ.Rtest <- function(k, alternative = "two.sided"){
  return(crit.SBQ.Mtest(k, m=2, alternative)[[1]])
}
