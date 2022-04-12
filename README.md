
<!-- README.md is generated from README.Rmd. Please edit that file -->

# backCUSUM

<!-- badges: start -->

<!-- badges: end -->

The goal of backCUSUM is to provide functionality to apply the methods
developed in the paper “[Backward CUSUM for Testing and Monitoring
Structural Change with an Application to COVID-19 Pandemic Data](https://arxiv.org/abs/2003.02682)” by [Sven
Otto](https://www.svenotto.com) and [Jörg
Breitung](https://wisostat.uni-koeln.de/en/institute/professors/breitung).
The repository also provides code to replicate all simulation results in
this paper.

Paper: https://doi.org/10.1017/S0266466622000159
Preprint: https://arxiv.org/abs/2003.02682

## Installation

You can install the package using the following command:

``` r
library(remotes)
install_github("ottosven/backCUSUM")
```

## Example

This is a basic example:

``` r
library(backCUSUM)
T <- 100
breakpoint <- floor(3*T/4)
beta1 = 0
beta2 = 1
u <- rnorm(T,0,1)
y <- c(rep(beta1,breakpoint), rep(beta2,T-breakpoint)) + u
Q.test(y~1)
BQ.test(y~1)
SBQ.test(y~1)
get.recresid(y~1)
```
