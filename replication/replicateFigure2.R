## ####################################################################
## ####################################################################
## Supplement for
## "Backward CUSUM for Testing and Monitoring Structural Change"
## by Sven Otto and Jörg Breitung.
## This R-script allows to reproduce Figure 2.
## ####################################################################
## ####################################################################
rm(list=ls())
start<-Sys.time()
library("backCUSUM")
## ##################################
## Reproducible random number
## generator for parallelization
## ##################################
RNGkind("L'Ecuyer-CMRG")
set.seed(42)
## ##################################
## Simulation
## ##################################
T <- 100
breakpoint <- floor(3*T/4)
beta1 = 0
beta2 = 1
u <- rnorm(T,0,1)
y <- c(rep(beta1,breakpoint), rep(beta2,T-breakpoint)) + u
## Recursive residuals
Qt <- Q.detector(y~1)
BQt <- BQ.detector(y~1)
wt <- get.recresid(y~1)
# ## Forward and backward CUSUM
# Qt1 <- abs(cumsum(wt))/sqrt(T)/sd(recresid(y ~ 1))
# BQt1 <- abs(sum(wt) - c(0,cumsum(wt)[-length(wt)]))/sqrt(T)/sd(recresid(y ~ 1))
## ##################################
## Figure 2
## ##################################
pdf('figure2.pdf', width = 17, height = 7.5, pointsize = 18)
par(mfrow=c(1,2))
plot(Qt, lwd = 4, type='l', col='blue', xlab = "", ylab = expression(), ylim=c(-2, 3.2), main='Forward CUSUM')
mtext(text = "time",side = 1, line = 2)
lines((0.948*(1+2*(1:T)/T)), lty=2, col='red', lwd=3)     # critical boundaries for alpha = 0.0001
lines(wt, col='grey', lwd = 1)
lines( x=c(1,100), y=c(0, 0), col='black', lwd=1, lty = 2)
plot(BQt, lwd = 4, type='l', col='blue', xlab = "", ylab = expression(), ylim=c(-2, 3.2), main='Backward CUSUM')
mtext(text = "time",side = 1, line = 2)
lines((0.948*(1+2*(T:1)/T)), lty=2, col='red', lwd=3)
lines(wt, col='grey', lwd = 1)
lines( x=c(1,100), y=c(0, 0), col='black', lwd=1, lty = 2)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", c("detector statistic", "linear boundary (5%)", "recursive residuals"), horiz = TRUE, bty = "n", cex=1 , col = c('blue', 'red', 'grey'), lty=c(1, 2, 1), lwd = c(4, 3, 1) )
dev.off()
## ##################################
## ##################################
end<-Sys.time()
end-start
