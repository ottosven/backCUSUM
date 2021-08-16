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
library(backCUSUM)
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
beta1 = 0
beta2 = 1
u <- rnorm(T,0,1)
## ##################################
## Figure 2
## ##################################
pdf('figure2.pdf', width = 17, height = 10, pointsize = 18)
par(mfrow=c(2,2))
breakpoint <- floor(3*T/4)
y <- c(rep(beta1,breakpoint), rep(beta2,T-breakpoint)) + u
Qt <- Q.test(y~1)$detector
BQt <- BQ.test(y~1)$detector
wt <- get.recresid(y~1)
plot(Qt, lwd = 4, type='l', col='blue', xlab = "", ylab = expression(), ylim=c(-2, 4.5), main='Forward CUSUM with a break at 75')
mtext(text = "time",side = 1, line = 2)
lines((0.947*(1+2*(1:T)/T)), lty=2, col='red', lwd=3)
lines((1.55*(1+2*(1:T)/T)), lty=4, col='coral4', lwd=3)
lines(wt, col='grey', lwd = 1)
lines( x=c(1,100), y=c(0, 0), col='black', lwd=1, lty = 2)
plot(BQt, lwd = 4, type='l', col='blue', xlab = "", ylab = expression(), ylim=c(-2, 4.5), main='Backward CUSUM with a break at 75')
mtext(text = "time",side = 1, line = 2)
lines((0.947*(1+2*(T:1)/T)), lty=2, col='red', lwd=3)
lines((1.55*(1+2*(T:1)/T)), lty=4, col='coral4', lwd=3)
lines(wt, col='grey', lwd = 1)
lines( x=c(1,100), y=c(0, 0), col='black', lwd=1, lty = 2)
#
breakpoint <- floor(T/4)
y <- c(rep(beta1,breakpoint), rep(beta2,T-breakpoint)) + u
Qt <- Q.test(y~1)$detector
BQt <- BQ.test(y~1)$detector
wt <- get.recresid(y~1)
plot(Qt, lwd = 4, type='l', col='blue', xlab = "", ylab = expression(), ylim=c(-2, 4.5), main='Forward CUSUM with a break at 25')
mtext(text = "time",side = 1, line = 2)
lines((0.947*(1+2*(1:T)/T)), lty=2, col='red', lwd=3)
lines((1.55*(1+2*(1:T)/T)), lty=4, col='coral4', lwd=3)
lines(wt, col='grey', lwd = 1)
lines( x=c(1,100), y=c(0, 0), col='black', lwd=1, lty = 2)
plot(BQt, lwd = 4, type='l', col='blue', xlab = "", ylab = expression(), ylim=c(-2, 4.5), main='Backward CUSUM with a break at 25')
mtext(text = "time",side = 1, line = 2)
lines((0.947*(1+2*(T:1)/T)), lty=2, col='red', lwd=3)
lines((1.55*(1+2*(T:1)/T)), lty=4, col='coral4', lwd=3)
lines(wt, col='grey', lwd = 1)
lines( x=c(1,100), y=c(0, 0), col='black', lwd=1, lty = 2)
#
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", c("detector statistic", "recursive residuals", "linear boundary (5%)", "linear boundary (0.1%)"), horiz = TRUE, bty = "n", cex=1, seg.len=4 , col = c('blue', 'grey', 'red', 'coral4'), lty=c(1, 1, 2, 4), lwd = c(4, 1, 3, 3) )
dev.off()
## ##################################
## ##################################
end<-Sys.time()
end-start
