#' Breakpoint estimation
#'
#' @param formula Specification of the linear regression model by an object of the class "formula"
#' @param type A character string specifying the breakpoint estimator type; must be one of "BQ" (backward CUSUM estimator, default) or "ML" (maximum likelihood estimator, see Bai (1997).
#' @param H An optional matrix for the partial breaks. \eqn{H} must have orthonormal columns. The full structural break model is considered as the default setting (NULL).

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
breakpoint.est <- function(formula, type = "BQ", H = NULL){
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
  } else {
    y <- model.frame(formula)[,1]
    X <- model.matrix(formula)
    RSS <- function(t) ( deviance(lm(y[1:t] ~ X[1:t,] - 1)) + deviance(lm(y[(t+1):T] ~ X[(t+1):T,] - 1)) )
    break.est <- which.min(c(rep(Inf,k),sapply(k:(T-k), RSS),rep(Inf,k)))
  }
  return(break.est)
}
