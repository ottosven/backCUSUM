#' Breakpoint estimation
#'
#' @param formula Specification of the linear regression model by an object of the class "formula"
#' @param type A character string specifying the breakpoint estimator type; must be one of "BQ" (backward CUSUM estimator, default) or "ML" (maximum likelihood estimator).
#' @param H An optional matrix for the partial hypothesis \eqn{H'\beta_t = H'\beta_0}, where \eqn{H'Q_t} is considered instead of \eqn{Q_t}.
#' \eqn{H} must have orthonormal columns. For a test for a break in the intercept, H can also set to the string "intercept".
#' The full structural break test is considered as the default setting (NULL).

#'
#' @return The estimated location of the breakpoint
#' @export
#'
#' @examples
#' T <- 100
#' u <- rnorm(T,0,1)
#' y <- 2 + c(rep(0,95), rep(0.8,5)) + u
#' breakpoint.est(y~1)
#' breakpoint.est(y~1, type = "ML")
breakpoint.est <- function(formula, type = c("BQ", "ML"), H = NULL){
  type=match.arg(type)
  T <- dim(model.matrix(formula))[1]
  k <- dim(model.matrix(formula))[2]
  if(type == "BQ"){
    if (is.null(H)){
      Q <- get.cusumprocess(formula, T)
    } else {
      Q <- get.partialcusum(formula, T, H)
      k <- dim(H)[2]
    }
    BQ <- cbind(Q[,T],matrix(Q[,T] - Q[,1:(T-1)], nrow = k))
    break.est <- which.max(BQ/sqrt(T-(1:T)+1))
  }
  if(type == "ML"){
    y <- model.frame(formula)[,1]
    X <- model.matrix(formula)
    RSS <- function(t) ( deviance(lm(y[1:t] ~ X[1:t,] - 1)) + deviance(lm(y[(t+1):T] ~ X[(t+1):T,] - 1)) )
    break.est <- which.min(c(rep(Inf,k),sapply(k:(T-k), RSS),rep(Inf,k)))
  }
  return(break.est)
}
