
<!-- README.md is generated from README.Rmd. Please edit that file -->

# backCUSUM

<!-- badges: start -->

<!-- badges: end -->

The goal of backCUSUM is to provide functionality to apply the methods
developed in the paper “[Backward CUSUM for Testing and Monitoring
Structural Change](https://arxiv.org/abs/2003.02682)” by [Sven
Otto](https://www.svenotto.com) and [Jörg
Breitung](https://wisostat.uni-koeln.de/en/institute/professors/breitung).
The repository also provides code to replicate all simulation results in
this paper.

## Installation

You can install the package using the following command:

``` r
library('devtools')
install_github("ottosven/backCUSUM")
```

## Example

This is a basic example to reproduce Figure 2 of the paper:

``` r
library(backCUSUM)
RNGkind("L'Ecuyer-CMRG")
set.seed(42)
T <- 100
breakpoint <- floor(3*T/4)
beta1 = 0
beta2 = 1
u <- rnorm(T,0,1)
y <- c(rep(beta1,breakpoint), rep(beta2,T-breakpoint)) + u
Qt <- Q.detector(y~1)
BQt <- BQ.detector(y~1)
wt <- get.recresid(y~1)
```

This is the plot:

<img src="man/figures/README-pressure-1.png" width="100%" />
