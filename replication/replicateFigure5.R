## ####################################################################
## ####################################################################
## Supplement for
## "Backward CUSUM for Testing and Monitoring Structural Change"
## by Sven Otto and Jörg Breitung.
## This R-script allows to reproduce Figure 5.
## ####################################################################
## ####################################################################
rm(list=ls())
## ##################################
## Load simulated data
## ##################################
load(file = "figure5_data.RData")
## ##################################
library(backCUSUM)
library(tsutils)
library(xts)
library(dyn)
library(cointReg)
## ##################################
total = xts(US.data$total, order.by = as.Date(US.data$dates))
newinfections = diff(total)
sdiffdata = diff(newinfections,7)
lagdata = xts(embed(na.omit(sdiffdata), 8), order.by = time(na.omit(sdiffdata)[-(1:7)]))
colnames(lagdata) = paste("lag",0:7,sep="")
newcases.pt = window(newinfections, start = "2020-03-10", end = "2020-10-31")
sdiff.pt = window(sdiffdata, start = "2020-03-10", end = "2020-10-31")
##########
## Monitoring specifications
##########
trainingstart1 = "2020-04-10"
trainingstart2 = "2020-07-20"
traininglength = 42 #6 weeks
##########
## CSW monitoring
csw <- function(formula, T, alpha = 0.05, alternative = "two.sided"){
  n <- dim(model.matrix(formula))[1] #current time point
  k <- dim(model.matrix(formula))[2]
  H <- matrix(c(1,numeric(k-1)), ncol = 1)
  detector <- backCUSUM::Q.mon(formula, T, alternative = alternative, H = H)$detector
  # boundary function
  r <- (1:n)/T
  if(alternative == "two.sided"){
    boundary.CSW <- sqrt(r[(T+1):n]*(log(r[(T+1):n]/alpha^2)))
  } else {
    boundary.CSW <- sqrt(r[(T+1):n]*(log(r[(T+1):n]/(2*alpha)^2)))
  }
  CSW <- detector/boundary.CSW
  # maximum statistic
  statistic <- max(CSW)
  # critical values and test decision
  rejection <- statistic > 1
  # detection time point
  detectiontime <- unname(T + which(CSW > 1)[1])
  return(list(detector = round(unname(detector),6), boundary = round(boundary.CSW,6), rejection = rejection, detectiontime = detectiontime, statistic = round(statistic,6)))
}
##########
## FIRST MONITORING PERIOD: DYNAMIC
##########
lagdata1 = window(lagdata, start = trainingstart1, end = "2020-12-31")
model=lagdata1[,1]~lagdata1[,3]+lagdata1[,8]
H <- matrix(c(1,0,0), ncol = 1)
recres <- xts(get.recresid(model), order.by = time(lagdata1))
QChu <- csw(model, T=traininglength, alternative = "greater", alpha = 0.05)
Q.detector <- QChu$detector/QChu$boundary
Q.scaled = xts(Q.detector, order.by = time(lagdata1)[-(1:traininglength)])
Q.detection = which(Q.detector > 1)[1] + traininglength
SBQ <- SBQ.mon(model, T=traininglength, alternative = "greater", H=H)
detect.max = SBQ$detector.array
bound = matrix(nrow = dim(detect.max)[1], ncol = dim(detect.max)[1])
for(t in 1:dim(detect.max)[1]){
  for(s in 1:t){
    bound[s,t] = sqrt((t+traininglength)/traininglength)*(1+2*(t-s+1)/traininglength)
  }
}
SBQ.detector = apply(detect.max/bound, 2, max, na.rm=TRUE)
SBQ.scaled = xts(SBQ.detector/0.911, time(lagdata1)[-(1:traininglength)])
SBQ.detection = which(SBQ.scaled > 1)[1] + traininglength
testingperiod = (traininglength+1):SBQ.detection
model2=lagdata[testingperiod,1]~lagdata[testingperiod,3]+lagdata[testingperiod,8]
BQest <- breakpoint.est(model2, type="BQ", H=H) + traininglength
modeldata1 = list(
  newcases=window(newinfections, end = "2020-11-01"),
  Q.det = Q.scaled,
  SBQ.det = SBQ.scaled,
  recresid = recres,
  train.start = time(lagdata1)[1],
  mon.start = time(lagdata1)[traininglength + 1] ,
  Q.dettime = time(Q.scaled[which(Q.scaled > 1)[1]]),
  SBQ.dettime = time(SBQ.scaled[which(SBQ.scaled > 1)[1]]),
  break.est = time(lagdata1)[BQest],
  alldates = time(lagdata1)
)
##########
##########
## SECOND MONITORING PERIOD: DYNAMIC
##########
lagdata1 = window(lagdata, start = trainingstart2, end = "2020-12-31")
model=lagdata1[,1]~lagdata1[,3]+lagdata1[,8]
H <- matrix(c(1,0,0), ncol = 1)
recres <- xts(get.recresid(model), order.by = time(lagdata1))
QChu <- csw(model, T=traininglength, alternative = "greater", alpha = 0.05)
Q.detector <- QChu$detector/QChu$boundary
Q.scaled = xts(Q.detector, order.by = time(lagdata1)[-(1:traininglength)])
Q.detection = which(Q.detector > 1)[1] + traininglength
SBQ <- SBQ.mon(model, T=traininglength, alternative = "greater", H=H)
detect.max = SBQ$detector.array
bound = matrix(nrow = dim(detect.max)[1], ncol = dim(detect.max)[1])
for(t in 1:dim(detect.max)[1]){
  for(s in 1:t){
    bound[s,t] = sqrt((t+traininglength)/traininglength)*(1+2*(t-s+1)/traininglength)
  }
}
SBQ.detector = apply(detect.max/bound, 2, max, na.rm=TRUE)
SBQ.scaled = xts(SBQ.detector/0.911, time(lagdata1)[-(1:traininglength)])
SBQ.detection = which(SBQ.scaled > 1)[1] + traininglength
testingperiod = (traininglength+1):SBQ.detection
model2=lagdata[testingperiod,1]~lagdata[testingperiod,3]+lagdata[testingperiod,8]
BQest <- breakpoint.est(model2, type="BQ", H=H) + traininglength
##
modeldata2 = list(
  newcases=window(newinfections, end = "2020-11-01"),
  Q.det = Q.scaled,
  SBQ.det = SBQ.scaled,
  recresid = recres,
  train.start = time(lagdata1)[1],
  mon.start = time(lagdata1)[traininglength + 1] ,
  Q.dettime = time(Q.scaled[which(Q.scaled > 1)[1]]),
  SBQ.dettime = time(SBQ.scaled[which(SBQ.scaled > 1)[1]]),
  break.est = time(lagdata1)[BQest],
  alldates = time(lagdata1)
)
##########
##########
## FIRST MONITORING PERIOD: ROBUST
##########
lagdata1 = window(lagdata, start = trainingstart1, end = "2020-12-31")
model = lagdata1[,1] ~ 1
recres <- xts(get.recresid(model), order.by = time(lagdata1))
QChu <- csw(model, T=traininglength, alternative = "greater", alpha = 0.05)
rec = recres[1:traininglength]
omegasigma = c(sqrt(getLongRunVar(rec)$Omega)/sd(rec))
Q.detector <- QChu$detector/QChu$boundary/omegasigma
Q.scaled = xts(Q.detector, order.by = time(lagdata1)[-(1:traininglength)])
Q.detection = which(Q.detector > 1)[1] + traininglength
SBQ <- SBQ.mon(model, T=traininglength, alternative = "greater")
detect.max = SBQ$detector.array
bound = matrix(nrow = dim(detect.max)[1], ncol = dim(detect.max)[1])
for(t in 1:dim(detect.max)[1]){
  for(s in 1:t){
    bound[s,t] = sqrt((t+traininglength)/traininglength)*(1+2*(t-s+1)/traininglength)
  }
}
SBQ.detector = apply(detect.max/bound, 2, max, na.rm=TRUE)/omegasigma
SBQ.scaled = xts(SBQ.detector/0.911, order.by = time(lagdata1)[-(1:traininglength)])
SBQ.detection = which(SBQ.scaled > 1)[1] + traininglength
testingperiod = (traininglength+1):SBQ.detection
model2=lagdata[testingperiod,1]~lagdata[testingperiod,3]+lagdata[testingperiod,8]
BQest <- breakpoint.est(model2, type="BQ", H=H) + traininglength
modeldata3 = list(
  newcases=window(newinfections, end = "2020-11-01"),
  Q.det = Q.scaled,
  SBQ.det = SBQ.scaled,
  recresid = recres,
  train.start = time(lagdata1)[1],
  mon.start = time(lagdata1)[traininglength + 1] ,
  Q.dettime = time(Q.scaled[which(Q.scaled > 1)[1]]),
  SBQ.dettime = time(SBQ.scaled[which(SBQ.scaled > 1)[1]]),
  break.est = time(lagdata1)[BQest],
  alldates = time(lagdata1)
)
##########
##########
## SECOND MONITORING PERIOD: ROBUST
##########
lagdata1 = window(lagdata, start = trainingstart2, end = "2020-12-31")
model = lagdata1[,1] ~ 1
recres <- xts(get.recresid(model), order.by = time(lagdata1))
QChu <- csw(model, T=traininglength, alternative = "greater", alpha = 0.05)
rec = recres[1:traininglength]
library(cointReg)
omegasigma = c(sqrt(getLongRunVar(rec)$Omega)/sd(rec))
Q.detector <- QChu$detector/QChu$boundary/omegasigma
Q.scaled = xts(Q.detector, order.by = time(lagdata1)[-(1:traininglength)])
Q.detection = which(Q.detector > 1)[1] + traininglength
SBQ <- SBQ.mon(model, T=traininglength, alternative = "greater")
detect.max = SBQ$detector.array
bound = matrix(nrow = dim(detect.max)[1], ncol = dim(detect.max)[1])
for(t in 1:dim(detect.max)[1]){
  for(s in 1:t){
    bound[s,t] = sqrt((t+traininglength)/traininglength)*(1+2*(t-s+1)/traininglength)
  }
}
SBQ.detector = apply(detect.max/bound, 2, max, na.rm=TRUE)/omegasigma
SBQ.scaled = xts(SBQ.detector/0.911, order.by = time(lagdata1)[-(1:traininglength)])
SBQ.detection = which(SBQ.scaled > 1)[1] + traininglength
testingperiod = (traininglength+1):SBQ.detection
model2=lagdata[testingperiod,1]~lagdata[testingperiod,3]+lagdata[testingperiod,8]
BQest <- breakpoint.est(model2, type="BQ", H=H) + traininglength
modeldata4 = list(
  newcases=window(newinfections, end = "2020-11-01"),
  Q.det = Q.scaled,
  SBQ.det = SBQ.scaled,
  recresid = recres,
  train.start = time(lagdata1)[1],
  mon.start = time(lagdata1)[traininglength + 1] ,
  Q.dettime = time(Q.scaled[which(Q.scaled > 1)[1]]),
  SBQ.dettime = time(SBQ.scaled[which(SBQ.scaled > 1)[1]]),
  break.est = time(lagdata1)[BQest],
  alldates = time(lagdata1)
)
#######################################################################
options(scipen=5)
pdf(paste("figure5.pdf", sep=""), width=30, height=35, pointsize = 30)
par(mfrow=c(4,2))
##########
## PLOT NEW INFECTIONS + DETECTIONS (DYNAMIC)
##########
plot(time(newcases.pt), newcases.pt,
     main="New infections and detection times (dynamic modeling)",
     type="l", xlab = "", ylab="",
     xaxt = "n",
     cex.lab = 1, cex.axis=1, cex.main = 1.7, lwd=3)
dates = c("2020-03-01", "2020-05-01", "2020-07-01", "2020-09-01", "2020-11-01")
dates.lab = c("Mar", "May", "Jul", "Sep", "Nov")
axis(1, at = c(as.Date(dates)), labels = dates.lab, cex.axis=1)
x <- as.Date(c(as.Date(trainingstart1), as.Date(trainingstart1)+traininglength-1, as.Date(trainingstart1)+traininglength-1, as.Date(trainingstart1)))
y <- c(-100000, -100000, 800000, 800000)
polygon(x,y,  col='grey', border = 'black', density = 15, lty = 2)
x <- as.Date(c(as.Date(trainingstart2), as.Date(trainingstart2)+traininglength-1, as.Date(trainingstart2)+traininglength-1, as.Date(trainingstart2)))
polygon(x,y,  col='grey', border = 'black', density = 15, lty = 2)
abline(v=modeldata1$SBQ.dettime, lwd=5, lty = 2)
abline(v=modeldata1$Q.dettime, lwd=5, lty = 3)
abline(v=modeldata2$SBQ.dettime, lwd=5, lty = 2)
abline(v=modeldata2$Q.dettime, lwd=5, lty = 3)
legend("topleft", c("SBQ detection time", "Q detection time"), lwd = c(5,5), bty = "n", lty=c(2,3), seg.len=2, cex=1.2)
##########
## PLOT SEASONAL DIFFERENCES + DETECTIONS (DYNAMIC)
##########
plot(time(sdiff.pt), sdiff.pt,
     main="Seasonal differences and detection times (dynamic modeling)",
     type="l", xlab = "", ylab="",
     xaxt = "n",
     cex.lab = 1, cex.axis=1, cex.main = 1.7, lwd=3)
dates = c("2020-03-01", "2020-05-01", "2020-07-01", "2020-09-01", "2020-11-01")
dates.lab = c("Mar", "May", "Jul", "Sep", "Nov")
axis(1, at = c(as.Date(dates)), labels = dates.lab, cex.axis=1)
x <- as.Date(c(as.Date(trainingstart1), as.Date(trainingstart1)+traininglength-1, as.Date(trainingstart1)+traininglength-1, as.Date(trainingstart1)))
y <- c(-100000, -100000, 1000000, 1000000)
polygon(x,y,  col='grey', border = 'black', density = 15, lty = 2)
x <- as.Date(c(as.Date(trainingstart2), as.Date(trainingstart2)+traininglength-1, as.Date(trainingstart2)+traininglength-1, as.Date(trainingstart2)))
polygon(x,y,  col='grey', border = 'black', density = 15, lty = 2)
abline(v=modeldata3$SBQ.dettime, lwd=5, lty = 2)
abline(v=modeldata3$Q.dettime, lwd=5, lty = 3)
abline(v=modeldata4$SBQ.dettime, lwd=5, lty = 2)
abline(v=modeldata4$Q.dettime, lwd=5, lty = 3)
##########
## PLOT DETECTORS
##########
SBQ.pt1 = window(modeldata1$SBQ.det, start = as.Date(trainingstart1)+traininglength, end="2020-06-26")
Q.pt1 = window(modeldata1$Q.det, start = as.Date(trainingstart1)+traininglength, end="2020-06-26")
plot(time(SBQ.pt1), SBQ.pt1,
     main="Monitoring second wave (dynamic modeling)",
     type="b", xlab = "", ylab="",
     xaxt="n",
     cex.lab = 1, cex.axis=1, cex.main = 1.7, lwd = 3, lty = 1, pch=20, cex=0.9)
lines(time(Q.pt1), Q.pt1, lwd = 1, lty=1, pch=1, type="b", cex=0.7)
abline(h=1, lty=5)
dates = c("2020-05-22", "2020-05-29", "2020-06-05", "2020-06-12", "2020-06-19", "2020-06-26")
dates.lab = c("May 22", "May 29", "Jun 05", "Jun 12", "Jun 19", "Jun 26")
axis(1, at = c(as.Date(dates)), labels = dates.lab, cex.axis=1)
abline(v=modeldata1$break.est, lwd=5, lty=4)
legend("topleft", c("scaled SBQ detector", "scaled Q detector", "breakpoint estimator", "critical value"),
       lwd = c(3,1,4,1), bty = "n", lty=c(1,1,4,5), seg.len=2, cex=1.2, pch=c(20,1, NA, NA), pt.cex=c(0.9,0.8,1,1))
#
SBQ.pt3 = window(modeldata3$SBQ.det, start = as.Date(trainingstart1)+traininglength, end="2020-06-26")
Q.pt3 = window(modeldata3$Q.det, start = as.Date(trainingstart1)+traininglength, end="2020-06-26")
plot(time(SBQ.pt3), SBQ.pt3,
     main="Monitoring second wave (autocorrelation robust)",
     type="b", xlab = "", ylab="",
     xaxt="n",
     cex.lab = 1, cex.axis=1, cex.main = 1.7, lwd = 3, lty = 1, pch=20, cex=0.9)
lines(time(Q.pt3), Q.pt3, lwd = 1, lty=1, pch=1, type="b", cex=0.7)
abline(h=1, lty=5)
dates = c("2020-05-22", "2020-05-29", "2020-06-05", "2020-06-12", "2020-06-19", "2020-06-26")
dates.lab = c("May 22", "May 29", "Jun 05", "Jun 12", "Jun 19", "Jun 26")
axis(1, at = c(as.Date(dates)), labels = dates.lab, cex.axis=1)
abline(v=modeldata3$break.est, lwd=5, lty=4)
#
SBQ.pt1 = window(modeldata2$SBQ.det, start = as.Date(trainingstart2)+traininglength, end="2020-10-12")
Q.pt1 = window(modeldata2$Q.det, start = as.Date(trainingstart2)+traininglength, end="2020-10-12")
plot(time(SBQ.pt1), SBQ.pt1,
     main="Monitoring third wave (dynamic modeling)",
     type="b", xlab = "", ylab="",
     xaxt="n",
     cex.lab = 1, cex.axis=1, cex.main = 1.7, lwd = 3, pch=20, cex=0.9)
lines(time(Q.pt1), Q.pt1, lwd = 1, pch=1, type="b", cex=0.7)
abline(h=1, lty=5)
dates = c("2020-08-31", "2020-09-07", "2020-09-14", "2020-09-21", "2020-09-28", "2020-10-05", "2020-10-12")
dates.lab = c("Aug 31", "Sep 07", "Sep 14", "Sep 21", "Sep 28", "Oct 05", "Oct 12")
axis(1, at = c(as.Date(dates)), labels = dates.lab, cex.axis=1)
abline(v=modeldata2$break.est, lwd=5, lty=4)
###################
## MODEL4
###################
SBQ.pt1 = window(modeldata4$SBQ.det, start = as.Date(trainingstart2)+traininglength, end="2020-10-12")
Q.pt1 = window(modeldata4$Q.det, start = as.Date(trainingstart2)+traininglength, end="2020-10-12")
plot(time(SBQ.pt1), SBQ.pt1,
     main="Monitoring third wave (autocorrelation robust)",
     type="b", xlab = "", ylab="",
     xaxt="n",
     cex.lab = 1, cex.axis=1, cex.main = 1.7, lwd = 3, pch=20, cex=0.9)
lines(time(Q.pt1), Q.pt1, lwd = 1, pch=1, type="b", cex=0.7)
abline(h=1, lty=5)
dates = c("2020-08-31", "2020-09-07", "2020-09-14", "2020-09-21", "2020-09-28", "2020-10-05", "2020-10-12")
dates.lab = c("Aug 31", "Sep 07", "Sep 14", "Sep 21", "Sep 28", "Oct 05", "Oct 12")
axis(1, at = c(as.Date(dates)), labels = dates.lab, cex.axis=1)
abline(v=modeldata4$break.est, lwd=5, lty=4)
###################
## RECURSIVE RESIDUALS AND DETETORS
###################
###################
## MODEL1
###################
recresid.pt1 = window(modeldata1$recresid, start = as.Date(trainingstart1)+traininglength, end="2020-06-26")
recresid.pt2 = window(modeldata3$recresid, start = as.Date(trainingstart1)+traininglength, end="2020-06-26")
plot(time(recresid.pt2), recresid.pt2,
     main="Recursive residuals (second wave monitoring)",
     type="l", xlab = "", ylab="",
     xaxt="n",
     lwd=4,
     cex.lab = 1, cex.axis=1, cex.main = 1.7, col ="grey")
lines(time(recresid.pt1), recresid.pt1, lwd = 1, lty=1)
abline(h=0)
dates = c("2020-05-22", "2020-05-29", "2020-06-05", "2020-06-12", "2020-06-19", "2020-06-26")
dates.lab = c("May 22", "May 29", "Jun 05", "Jun 12", "Jun 19", "Jun 26")
axis(1, at = c(as.Date(dates)), labels = dates.lab, cex.axis=1)
legend("topleft", c("recursive residuals (dynamic modeling)", "recursive residuals (autocorrelation robust)"),
       lwd = c(1,4), bty = "n", col=c("black", "grey"), seg.len=2, cex=1.2)
#
recresid.pt1 = window(modeldata2$recresid, start = as.Date(trainingstart2)+traininglength, end="2020-10-12")
recresid.pt2 = window(modeldata4$recresid, start = as.Date(trainingstart2)+traininglength, end="2020-10-12")
plot(time(recresid.pt2), recresid.pt2,
     main="Recursive residuals (third wave monitoring)",
     type="l", xlab = "", ylab="",
     xaxt="n",
     lwd=4,
     cex.lab = 1, cex.axis=1, cex.main = 1.7, col ="grey")
lines(time(recresid.pt1), recresid.pt1, lwd = 1, lty=1)
abline(h=0)
abline(v=modeldata1$SBQ.dettime, lwd=3, lty = 2)
abline(v=modeldata1$Q.dettime, lwd=3, lty = 3)
dates = c("2020-08-31", "2020-09-07", "2020-09-14", "2020-09-21", "2020-09-28", "2020-10-05", "2020-10-12")
dates.lab = c("Aug 31", "Sep 07", "Sep 14", "Sep 21", "Sep 28", "Oct 05", "Oct 12")
axis(1, at = c(as.Date(dates)), labels = dates.lab, cex.axis=1)
#
dev.off()
