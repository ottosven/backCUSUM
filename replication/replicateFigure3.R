## ####################################################################
## ####################################################################
## Supplement for 
## "Backward CUSUM for Testing and Monitoring Structural Change"
## by Sven Otto and Jörg Breitung.
## This R-script allows to reproduce Figure 3.
## ####################################################################
## ####################################################################
rm(list=ls())
start<-Sys.time()
## ##################################
## Load simulated data
## ##################################
load(file = "figure3_data.RData")
## ##################################
## Figure 3
## ##################################
pdf('figure3.pdf', width=20, height=15, pointsize = 30 )
par(oma = c(1, 0, 0, 0))
par(mfrow=c(2,3))
c.range <- 0:30
##
plot(c.range, powcur.ret.r01[1,], type='l', col=1, main=expression(paste(tau, "* = 0.1")), ylab = '', xlab = '', lty=1, lwd=4, ylim = c(0.05,1), cex.lab = 1.2, cex.axis=1.2, cex.main = 1.5)
mtext(text = "rejection frequency",side = 2, line = 3, cex=0.7)
mtext(text = expression(paste("c/",sigma)),side = 1, line = 2, cex=0.7)
lines(c.range, powcur.ret.r01[2,], col="red", lty=2, lwd=4)
lines(c.range, powcur.ret.r01[3,], col="blue", lty=4, lwd=4)
##
plot(c.range, powcur.ret.r03[1,], type='l', col=1, main=expression(paste(tau, "* = 0.3")), ylab = '', xlab = '', lty=1, lwd=4, ylim = c(0.05,1), cex.lab = 1.2, cex.axis=1.2, cex.main = 1.5)
mtext(text = "rejection frequency",side = 2, line = 3, cex=0.7)
mtext(text = expression(paste("c/",sigma)),side = 1, line = 2, cex=0.7)
lines(c.range, powcur.ret.r03[2,], col="red", lty=2, lwd=4)
lines(c.range, powcur.ret.r03[3,], col="blue", lty=4, lwd=4)
##
plot(c.range, powcur.ret.r05[1,], type='l', col=1, main=expression(paste(tau, "* = 0.5")), ylab = '', xlab = '', lty=1, lwd=4, ylim = c(0.05,1), cex.lab = 1.2, cex.axis=1.2, cex.main = 1.5)
mtext(text = "rejection frequency",side = 2, line = 3, cex=0.7)
mtext(text = expression(paste("c/",sigma)),side = 1, line = 2, cex=0.7)
lines(c.range, powcur.ret.r05[2,], col="red", lty=2, lwd=4)
lines(c.range, powcur.ret.r05[3,], col="blue", lty=4, lwd=4)
##
plot(c.range, powcur.ret.r07[1,], type='l', col=1, main=expression(paste(tau, "* = 0.7")), ylab = '', xlab = '', lty=1, lwd=4, ylim = c(0.05,1), cex.lab = 1.2, cex.axis=1.2, cex.main = 1.5)
mtext(text = "rejection frequency",side = 2, line = 3, cex=0.7)
mtext(text = expression(paste("c/",sigma)),side = 1, line = 2, cex=0.7)
lines(c.range, powcur.ret.r07[2,], col="red", lty=2, lwd=4)
lines(c.range, powcur.ret.r07[3,], col="blue", lty=4, lwd=4)
##
plot(c.range, powcur.ret.r09[1,], type='l', col=1, main=expression(paste(tau, "* = 0.9")), ylab = '', xlab = '', lty=1, lwd=4, ylim = c(0.05,1), cex.lab = 1.2, cex.axis=1.2, cex.main = 1.5)
mtext(text = "rejection frequency",side = 2, line = 3, cex=0.7)
mtext(text = expression(paste("c/",sigma)),side = 1, line = 2, cex=0.7)
lines(c.range, powcur.ret.r09[2,], col="red", lty=2, lwd=4)
lines(c.range, powcur.ret.r09[3,], col="blue", lty=4, lwd=4)
##
rstar.range <- 1:39/40
plot(rstar.range, powcur.ret.c10[1,], type='l', col=1, main=expression(paste("c/",sigma," = 10")), ylab = '', xlab = '', lty=1, lwd=4, ylim = c(0.05,1), cex.lab = 1.2, cex.axis=1.2, cex.main = 1.5)
mtext(text = "rejection frequency",side = 2, line = 3, cex=0.7)
mtext(text = 'r*',side = 1, line = 2, cex=0.7)
lines(rstar.range, powcur.ret.c10[2,], col="red", lty=2, lwd=4)
lines(rstar.range, powcur.ret.c10[3,], col="blue", lty=4, lwd=4)
##
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", c("forward CUSUM", "backward CUSUM", "stacked backward CUSUM"), horiz = TRUE, bty = "n", cex=1.3 , col = c(1, "red", "blue"), lty = c(1,2,4), lwd = c(4, 4, 4), seg.len=3 )
dev.off()
## ##################################
## ##################################
end<-Sys.time()
end-start