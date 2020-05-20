## ####################################################################
## ####################################################################
## Supplement for
## "Backward CUSUM for Testing and Monitoring Structural Change"
## by Sven Otto and Jörg Breitung.
## This R-script allows to reproduce Figure 4.
## ####################################################################
## ####################################################################
rm(list=ls())
start<-Sys.time()
## ##################################
## Load simulated data
## ##################################
load(file = "figure4_data.RData")
## ##################################
## Figure 3
## ##################################
pdf('figure4.pdf', width=20, height=7.5, pointsize = 30)
par(oma = c(1, 0, 0, 0))
par(mfrow=c(1,3))
c.range <- 5:30
##
plot(c.range, delcur.mon.r15[1,-(1:5)], type='l', col=1, ylim = c(0,0.3), main=expression(paste(tau, "* = 1.5")), ylab = '', xlab = '', lty=1, lwd=4, cex.lab = 1.2, cex.axis=1.2, cex.main = 1.5)
mtext(text = expression(paste("c/",sigma)), side = 1, line = 2, cex=0.7)
mtext(text = 'relative mean delay', side = 2, line = 3, cex=0.7)
lines(c.range, delcur.mon.r15[2,-(1:5)], col="red", lty=2, lwd=4)
lines(c.range, delcur.mon.r15[3,-(1:5)], col="blue", lty=4, lwd=4)
##
plot(c.range, delcur.mon.r3[1,-(1:5)], type='l', col=1, ylim = c(0,0.3), main=expression(paste(tau, "* = 3")), ylab = '', xlab = '', lty=1, lwd=4, cex.lab = 1.2, cex.axis=1.2, cex.main = 1.5)
mtext(text = expression(paste("c/",sigma)), side = 1, line = 2, cex=0.7)
mtext(text = 'relative mean delay', side = 2, line = 3, cex=0.7)
lines(c.range, delcur.mon.r3[2,-(1:5)], col="red", lty=2, lwd=4)
lines(c.range, delcur.mon.r3[3,-(1:5)], col="blue", lty=4, lwd=4)
##
rstar.range <- 3*(1:39/40) + 1
plot(rstar.range, delcur.mon.c20[1,], type='l', col=1, main=expression(paste("c/",sigma," = 20")), ylab = '', xlab = '', lty=1, lwd=4, ylim = c(0,0.15), cex.lab = 1.2, cex.axis=1.2, cex.main = 1.5)
mtext(text = 'r*', side = 1, line = 2, cex=0.7)
mtext(text = 'relative mean delay', side = 2, line = 3, cex=0.7)
lines(rstar.range, delcur.mon.c20[2,], col="red", lty=2, lwd=4)
lines(rstar.range, delcur.mon.c20[3,], col="blue", lty=4, lwd=4)
##
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", c("stacked backward CUSUM", "forward CUSUM (linear boundary)", "forward CUSUM (radical boundary)"), horiz = TRUE, bty = "n", cex=1.2 , col = c(1, "red", "blue"), lty = c(1,2,4), lwd = c(4, 4, 4) )
dev.off()
## ##################################
## ##################################
end<-Sys.time()
end-start
