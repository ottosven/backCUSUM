% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cusumprocess.R
\name{get.cusumprocess}
\alias{get.cusumprocess}
\title{Multivariate CUSUM process}
\usage{
get.cusumprocess(formula, T)
}
\arguments{
\item{formula}{Specification of the linear regression model by an object of the class "formula"}

\item{T}{In the retropspective context: length of the sample. In the monitoring context: length of the training sample, where monitoring starts at T+1.}
}
\value{
A matrix containing the multivariate forward CUSUM process Q_T(r)
}
\description{
Multivariate CUSUM process
}
\examples{
T <- 100
u <- rnorm(T,0,1)
x <- rnorm(T,0,1)
y <- c(rep(0,T/2), rep(1,T/2)) + x + u
get.cusumprocess(y~x, T)
}
