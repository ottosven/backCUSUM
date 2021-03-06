% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/breakpoint.R
\name{breakpoint.est}
\alias{breakpoint.est}
\title{Breakpoint estimation}
\usage{
breakpoint.est(formula, type = "BQ", H = NULL)
}
\arguments{
\item{formula}{Specification of the linear regression model by an object of the class "formula"}

\item{type}{A character string specifying the breakpoint estimator type; must be one of "BQ" (backward CUSUM estimator, default) or "ML" (maximum likelihood estimator, see Bai (1997).}

\item{H}{An optional matrix for the partial breaks. \eqn{H} must have orthonormal columns. The full structural break model is considered as the default setting (NULL).}
}
\value{
The estimated location of the breakpoint
}
\description{
Breakpoint estimation
}
\examples{
T <- 100
u <- rnorm(T,0,1)
y <- 2 + c(rep(0,95), rep(0.8,5)) + u
breakpoint.est(y~1)
breakpoint.est(y~1, type = "ML")
}
