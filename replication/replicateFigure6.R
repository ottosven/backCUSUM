## ####################################################################
## ####################################################################
## Supplement for
## "Backward CUSUM for Testing and Monitoring Structural Change"
## by Sven Otto and Jörg Breitung.
## This R-script allows to reproduce Figure 6.
## ####################################################################
## ####################################################################
rm(list=ls())
start<-Sys.time()
library(backCUSUM)
##
T <- 100
breaktime <- floor(T/2)
beta1 = 0
beta2 = 1
y <- c(rep(beta1,breaktime), rep(beta2,T-breaktime))
wt <- get.recresid(y~1)
## ##################################
## Figure 6
## ##################################
pdf('figure6.pdf', width = 12, height = 7.5, pointsize = 16)
hstar <- abs(sum(wt) - c(0,cumsum(wt)[-length(wt)]))/sqrt(T)
scaledhstar <- abs(sum(wt) - c(0,cumsum(wt)[-length(wt)]))/sqrt(T-(1:T)+1)
plot((1:T)/T, scaledhstar, lwd = 4, type='l', col='blue', xlab = "", ylab = expression(), ylim=c(0, 5), main='')
lines((1:T)/T, hstar, lwd = 4, type='l', col='black', lty = 2)
mtext(text = "r",side = 1, line = 2)
dev.off()
## ##################################
## ##################################
end<-Sys.time()
end-start
