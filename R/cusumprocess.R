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
  sqCTinv <- expm::sqrtm(solve((t(X) %*% X)/T))
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





#' Simulates critical values for the multivariate forward CUSUM test
#'
#' The asymptotic critical values for the linear boundary \eqn{d(r) = 1+2r} are simulated based on \eqn{MC} Monte Carlo replications, where the limiting process under the null hypothesis is approximated on a grid of \eqn{T} equidistant points.
#' The critical values of the retrospective backward CUSUM coincide with those of the forward CUSUM for m=2.
#'
#' @param k Number of regressors and dimension of the CUSUM process.
#' @param m The length of the relative monitoring period; m = Inf for infinite horizon monitoring; for retrospective testing set m = 2 (default).
#' @param MC Number of Monte Carlo replications.
#' @param T gridsize.
#' @param levels A vector if significance levels.
#' @param alternative A character string specifying the alternative hypothesis; must be one of "two.sided" (default) or "one.sided".
#' @param CORE Number of cores used for parralell simulation
#'
#' @return A vector of critical values for different significance levels.
#' @export
#'
#' @examples
#' sim.crit.Q(1, m=Inf)
sim.crit.Q <- function(k, m = 2,  MC = 1000, T = 10000, levels = c(0.1, 0.05, 0.01, 0.001), alternative = "two.sided", CORE = 1){
  r <- (1:T)/T
  boundary.transf <- 1+r
  if(m == Inf){
    MRange <- 1:T
  } else {
    MRange <- 1:floor(T*(m-1)/m)
  }
  boundary.mrange <- boundary.transf[MRange]
  BrownianMotion <- function(T)  ( cumsum(rnorm(T,0,sqrt(1/T))) )
  sim.dist.Q <- function(...){
    W <- sapply(rep(T,k), BrownianMotion)
    B <- W - t(W[T,]*matrix(rep(r,k), nrow = k, byrow = TRUE))
    B.mrange <- matrix(B[MRange,], ncol = k)
    Q <- max((apply(abs(B.mrange), 1, max)/boundary.mrange))
    return(Q)
  }
  realizations <- parallel::mcmapply(sim.dist.Q,1:MC,mc.cores = getOption("mc.cores", CORE))
  if(alternative == "one.sided"){
    quantiles <- quantile(realizations, 1-levels*2)
  } else {
    quantiles <- quantile(realizations, 1-levels)
  }
  names(quantiles) <- levels
  return(quantiles)
}




#' Simulates critical values for the multivariate stacked backward CUSUM test
#'
#' The asymptotic critical values for the linear boundary \eqn{d(r) = 1+2r} are simulated based on \eqn{MC} Monte Carlo replications, where the limiting process under the null hypothesis is approximated on a grid of \eqn{T} equidistant points.
#'
#' @param k Number of regressors and dimension of the CUSUM process.
#' @param m The length of the relative monitoring period; m = Inf for infinite horizon monitoring; for retrospective testing set m = 2 (default).
#' @param MC Number of Monte Carlo replications.
#' @param T gridsize.
#' @param levels A vector if significance levels.
#' @param alternative A character string specifying the alternative hypothesis; must be one of "two.sided" (default) or "one.sided".
#' @param CORE Number of cores used for parralell simulation
#'
#' @return A vector of critical values for different significance levels.
#' @export
#'
#' @examples
#' sim.crit.SBQ(1, m=Inf)
sim.crit.SBQ <- function(k, m = 2, MC = 100, T = 1000, levels = c(0.1, 0.05, 0.01, 0.001), alternative = "two.sided", CORE = 1){
  r <- (1:T)/T
  # boundary <- 1+2*(m-1)*r
  # boundary.transf <- 1+r
  if(m == Inf){
    MRange <- 2:T
  } else {
    MRange <- 2:floor(T*(m-1)/m)
  }
  # boundary.mrange <- boundary.transf[MRange]
  BrownianMotion <- function(T)  ( cumsum(rnorm(T,0,sqrt(1/T))) )
  getSBQMax <- function(rT, B){
    sT <- 1:(rT-1)
    bound <- (1-rT/T)*(1-sT/T)+2*(rT-sT)/T
    Diff <- (1-sT/T)*matrix(rep(B[rT,],rT-1), nrow = rT-1, byrow = TRUE) - (1-rT/T)*B[sT,]
    Inner <- matrix(abs(Diff), nrow = rT-1)
    InnerMax <- apply(Inner, 1, max)
    SBQ <- max(InnerMax/bound)
  }
  sim.dist.SBQ <- function(...){
    W <- sapply(rep(T,k), BrownianMotion)
    B <- W - t(W[T,]*matrix(rep(r,k), nrow = k, byrow = TRUE))
    B.mrange <- matrix(B[c(1,MRange),], ncol = k)
    SBQ <- max(sapply(MRange, getSBQMax, B=B.mrange))
    return(SBQ)
  }
  realizations <- parallel::mcmapply(sim.dist.SBQ,1:MC,mc.cores = getOption("mc.cores", CORE))
  if(alternative == "one.sided"){
    quantiles <- quantile(realizations, 1-levels*2)
  } else {
    quantiles <- quantile(realizations, 1-levels)
  }
  names(quantiles) <- levels
  return(quantiles)
}
