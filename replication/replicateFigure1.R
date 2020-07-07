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
load(file = "figure1a_data.RData")
## ##################################
## Figure 3
## ##################################
pdf('figure1.pdf', width=30, height=22.5, pointsize = 30 )
# par(oma = c(1, 0, 0, 0))
par(mfrow=c(3,3))
c.range <- 0:30
##
plot(c.range, powcur.ret.r01[1,], type='l', col=1, main=expression(paste(tau, "* = 0.1")), ylab = '', xlab = '', lty=1, lwd=4, ylim = c(0.05,1), cex.lab = 1.2, cex.axis=1.2, cex.main = 2)
mtext(text = "rejection frequency",side = 2, line = 2.4, cex=0.8)
mtext(text = expression(paste("c/",sigma)),side = 1, line = 2.2, cex=0.9)
lines(c.range, powcur.ret.r01[2,], lty=2, lwd=4)
lines(c.range, powcur.ret.r01[3,], lty=3, lwd=4)
##
plot(c.range, powcur.ret.r03[1,], type='l', col=1, main=expression(paste(tau, "* = 0.3")), ylab = '', xlab = '', lty=1, lwd=4, ylim = c(0.05,1), cex.lab = 1.2, cex.axis=1.2, cex.main = 2)
mtext(text = "rejection frequency",side = 2, line = 2.4, cex=0.8)
mtext(text = expression(paste("c/",sigma)),side = 1, line = 2.2, cex=0.9)
lines(c.range, powcur.ret.r03[2,], lty=2, lwd=4)
lines(c.range, powcur.ret.r03[3,], lty=3, lwd=4)
##
plot(c.range, powcur.ret.r05[1,], type='l', col=1, main=expression(paste(tau, "* = 0.5")), ylab = '', xlab = '', lty=1, lwd=4, ylim = c(0.05,1), cex.lab = 1.2, cex.axis=1.2, cex.main = 2)
mtext(text = "rejection frequency",side = 2, line = 2.4, cex=0.8)
mtext(text = expression(paste("c/",sigma)),side = 1, line = 2.2, cex=0.9)
lines(c.range, powcur.ret.r05[2,], lty=2, lwd=4)
lines(c.range, powcur.ret.r05[3,], lty=3, lwd=4)
##
plot(c.range, powcur.ret.r07[1,], type='l', col=1, main=expression(paste(tau, "* = 0.7")), ylab = '', xlab = '', lty=1, lwd=4, ylim = c(0.05,1), cex.lab = 1.2, cex.axis=1.2, cex.main = 2)
mtext(text = "rejection frequency",side = 2, line = 2.4, cex=0.8)
mtext(text = expression(paste("c/",sigma)),side = 1, line = 2.2, cex=0.9)
lines(c.range, powcur.ret.r07[2,], lty=2, lwd=4)
lines(c.range, powcur.ret.r07[3,], lty=3, lwd=4)
##
plot(c.range, powcur.ret.r09[1,], type='l', col=1, main=expression(paste(tau, "* = 0.9")), ylab = '', xlab = '', lty=1, lwd=4, ylim = c(0.05,1), cex.lab = 1.2, cex.axis=1.2, cex.main = 2)
mtext(text = "rejection frequency",side = 2, line = 2.4, cex=0.8)
mtext(text = expression(paste("c/",sigma)),side = 1, line = 2.2, cex=0.9)
lines(c.range, powcur.ret.r09[2,], lty=2, lwd=4)
lines(c.range, powcur.ret.r09[3,], lty=3, lwd=4)
##
rstar.range <- 1:39/40
plot(rstar.range, powcur.ret.c10[1,], type='l', col=1, main=expression(paste("c/",sigma," = 10")), ylab = '', xlab = '', lty=1, lwd=4, ylim = c(0.05,1), cex.lab = 1.2, cex.axis=1.2, cex.main = 2)
mtext(text = "rejection frequency",side = 2, line = 2.4, cex=0.8)
mtext(text = expression(paste(tau, "*")),side = 1, line = 2.2, cex=0.9)
lines(rstar.range, powcur.ret.c10[2,], lty=2, lwd=4)
lines(rstar.range, powcur.ret.c10[3,], lty=3, lwd=4)
##
load(file = "figure1b_data.RData")
c.range <- 5:30
##
plot(c.range, delcur.mon.r15[1,-(1:5)], type='l', col=1, ylim = c(0,0.3), main=expression(paste(tau, "* = 1.5")), ylab = '', xlab = '', lty=3, lwd=4, cex.lab = 1.2, cex.axis=1.2, cex.main = 2)
mtext(text = expression(paste("c/",sigma)), side = 1, line = 2.2, cex=0.9)
mtext(text = 'relative mean delay', side = 2, line = 2.4, cex=0.8)
lines(c.range, delcur.mon.r15[2,-(1:5)], lty=1, lwd=4)
lines(c.range, delcur.mon.r15[3,-(1:5)], lty=4, lwd=4)
##
plot(c.range, delcur.mon.r3[1,-(1:5)], type='l', col=1, ylim = c(0,0.3), main=expression(paste(tau, "* = 3")), ylab = '', xlab = '', lty=3, lwd=4, cex.lab = 1.2, cex.axis=1.2, cex.main = 2)
mtext(text = expression(paste("c/",sigma)), side = 1, line = 2.2, cex=0.9)
mtext(text = 'relative mean delay', side = 2, line = 2.4, cex=0.8)
lines(c.range, delcur.mon.r3[2,-(1:5)], lty=1, lwd=4)
lines(c.range, delcur.mon.r3[3,-(1:5)], lty=4, lwd=4)
##
rstar.range <- 3*(1:39/40) + 1
plot(rstar.range, delcur.mon.c20[1,], type='l', col=1, main=expression(paste("c/",sigma," = 20")), ylab = '', xlab = '', lty=3, lwd=4, ylim = c(0,0.15), cex.lab = 1.2, cex.axis=1.2, cex.main = 2)
mtext(text = expression(paste(tau, "*")),side = 1, line = 2.2, cex=0.9)
mtext(text = 'relative mean delay', side = 2, line = 2.4, cex=0.8)
lines(rstar.range, delcur.mon.c20[2,], lty=1, lwd=4)
lines(rstar.range, delcur.mon.c20[3,], lty=4, lwd=4)
##
dev.off()
## ##################################
## ##################################
end<-Sys.time()
end-start
