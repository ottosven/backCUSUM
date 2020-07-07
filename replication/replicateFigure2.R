## ####################################################################
## ####################################################################
## Supplement for
## "Backward CUSUM for Testing and Monitoring Structural Change"
## by Sven Otto and Jörg Breitung.
## This R-script allows to reproduce Figure 2.
## ####################################################################
## ####################################################################
rm(list=ls())
## ##################################
## Load simulated data
## ##################################
load(file = "figure2_data.RData")
## ##################################
library(backCUSUM)
library(tsutils)
## ##################################
## Load country data
countrydata.raw <- read.csv(url("https://data.humdata.org/hxlproxy/api/data-preview.csv?url=https%3A%2F%2Fraw.githubusercontent.com%2FCSSEGISandData%2FCOVID-19%2Fmaster%2Fcsse_covid_19_data%2Fcsse_covid_19_time_series%2Ftime_series_covid19_confirmed_global.csv&filename=time_series_covid19_confirmed_global.csv"))
countrydata <- data.frame(t(countrydata.raw[,-c(1:4)]))
colnames(countrydata) <- countrydata.raw[,2]
##
monitoring.dynamicmodel <- function(data, trainingstart, traininglength){
  startindex <- which(rownames(data)==trainingstart)
  todayindex <- dim(data)[1]
  data.mon <- data[startindex:todayindex,]
  maxlag <- 7
  T <- traininglength - maxlag
  lags <- lagmatrix(data.mon$newcases, c(0,1,7))
  rownames(lags) <- rownames(data.mon)
  lagsclean <- na.omit(lags)
  dates <- as.Date(rownames(lagsclean))
  model <- lagsclean[,1] ~ lagsclean[,-1]
  recres <- get.recresid(model)
  QChu <- Q.mon.csw(model, T=T, alternative = "greater", alpha = 0.05)
  H <- matrix(c(1,0,0), ncol = 1)
  SBQ <- SBQ.mon(model, T=T, alternative = "greater", H=H)
  Q.detector <- QChu$detector/QChu$boundary
  SBQ.detector <- SBQ$detector.scaled/SBQ$critical.value["0.05"]
  Q.detection <- QChu$detectiontime
  SBQ.detection<- SBQ$detectiontime["0.05"]
  #
  if(!is.na(SBQ.detection)){
    model2 <- lagsclean[1:SBQ.detection,1] ~ lagsclean[1:SBQ.detection,-1]
    BQest <- breakpoint.est(model2, type="BQ", H=H)
    break.est <- as.Date(rownames(lagsclean)[BQest])
  } else {
    break.est <- NA
  }
  newcases <- data$newcases[-1]
  alldates <- as.Date(rownames(data)[-1])
  train.start <- as.Date(trainingstart)
  mon.start <- as.Date(train.start) + traininglength
  Q.dettime <- dates[Q.detection]
  SBQ.dettime <- dates[SBQ.detection]
  list(
    newcases=newcases,
    Q.det = Q.detector,
    SBQ.det = SBQ.detector,
    recresid = recres,
    train.start = train.start,
    mon.start = mon.start,
    Q.dettime = Q.dettime,
    SBQ.dettime = SBQ.dettime,
    break.est = break.est,
    alldates = alldates,
    maxlag = maxlag
  )
}
##
getplot.newcases <- function(modeldata, name, startdate = "2020-03-10", enddate = NULL){
  if(is.null(enddate)) ( enddate = modeldata$alldates[length(modeldata$alldates)] )
  startindex <- which(modeldata$alldates == startdate)
  endindex <- which(modeldata$alldates == enddate)
  plotdata <- modeldata$newcases[startindex:endindex]
  plot(modeldata$alldates[startindex:endindex], modeldata$newcases[startindex:endindex], type='l', xlab = "", ylab="", main=paste("New infections", name), cex.lab = 1.2, cex.axis=1.2, cex.main = 1.5)
  abline(v=modeldata$break.est, lwd=2, lty=5)
  abline(v=modeldata$SBQ.dettime, lwd=2, lty = 1)
  abline(v=modeldata$Q.dettime, lwd=2, lty = 4)
  x <- as.Date(c(modeldata$train.start - modeldata$maxlag, modeldata$mon.start-1, modeldata$mon.start-1, modeldata$train.start - modeldata$maxlag))
  y <- c(-10000, -10000, 1000000, 1000000)
  polygon(x,y,  col='grey', border = 'black', density = 15, lty = 2)
}
##
getplot.recresid <- function(modeldata, name, enddate = NULL){
  if(is.null(enddate)) ( enddate = modeldata$alldates[length(modeldata$alldates)] )
  mon.dates <- modeldata$alldates[which(modeldata$alldates == modeldata$mon.start):length(modeldata$alldates)]
  T <- as.numeric(modeldata$mon.start - modeldata$train.start - modeldata$maxlag)
  plotdata <- modeldata$recresid[-(1:T)]
  endindex <- which(mon.dates == enddate)
  plot(mon.dates[1:endindex], plotdata[1:endindex], type='l', main=paste("Recursive residuals", name), xlab="", ylab="",cex.lab = 1.2, cex.axis=1.2, cex.main = 1.5)
  abline(h=0, lty = 5, lwd = 1, col="darkgrey")
}
##
getplot.detectors <- function(modeldata, name, enddate = NULL){
  if(is.null(enddate)) ( enddate = modeldata$alldates[length(modeldata$alldates)] )
  mon.dates <- modeldata$alldates[which(modeldata$alldates == modeldata$mon.start):length(modeldata$alldates)]
  T <- as.numeric(modeldata$mon.start - modeldata$train.start - modeldata$maxlag)
  plotdata1 <- modeldata$SBQ.det
  plotdata2 <- modeldata$Q.det
  endindex <- which(mon.dates == enddate)
  plot(mon.dates[1:endindex], plotdata1[1:endindex], type='l', ylim = c(min(modeldata$SBQ.det,modeldata$Q.det)-0.1, max(modeldata$SBQ.det,modeldata$Q.det,1)+0.3), main=paste("Scaled detectors", name), xlab="", ylab="",cex.lab = 1.2, cex.axis=1.2, cex.main = 1.5)
  lines(mon.dates[1:endindex], plotdata2[1:endindex], lty=3)
  abline(h=1, lty = 2, lwd = 2)
}
US.modeldata <- monitoring.dynamicmodel(US.data, "2020-04-10", 28)
getplot.newcases(US.modeldata, "US total", enddate = "2020-06-30")
getplot.recresid(US.modeldata, "US total")
getplot.detectors(US.modeldata, "US total")
##
AZ.modeldata <- monitoring.dynamicmodel(AZ.data, "2020-04-10", 28)
getplot.newcases(AZ.modeldata, "Arizona")
getplot.recresid(AZ.modeldata, "Arizona")
getplot.detectors(AZ.modeldata, "Arizona")
##
FL.modeldata <- monitoring.dynamicmodel(FL.data, "2020-04-10", 28)
getplot.newcases(FL.modeldata, "Florida")
getplot.recresid(FL.modeldata, "Florida")
getplot.detectors(FL.modeldata, "Florida")
##
NV.modeldata <- monitoring.dynamicmodel(NV.data, "2020-04-10", 28)
getplot.newcases(NV.modeldata, "Nevada")
getplot.recresid(NV.modeldata, "Nevada")
getplot.detectors(NV.modeldata, "Nevada")
##
TX.modeldata <- monitoring.dynamicmodel(TX.data, "2020-04-10", 28)
getplot.newcases(TX.modeldata, "Texas")
getplot.recresid(TX.modeldata, "Texas")
getplot.detectors(TX.modeldata, "Texas")
##
enddate = "2020-06-25"
pdf(paste("figure2.pdf", sep=""), width=30, height=37.5, pointsize = 30)
par(mfrow=c(5,3))
getplot.newcases(US.modeldata, "US total", enddate = enddate)
getplot.recresid(US.modeldata, "US total", enddate = enddate)
getplot.detectors(US.modeldata, "US total", enddate = enddate)
getplot.newcases(AZ.modeldata, "Arizona", enddate = enddate)
getplot.recresid(AZ.modeldata, "Arizona", enddate = enddate)
getplot.detectors(AZ.modeldata, "Arizona", enddate = enddate)
getplot.newcases(FL.modeldata, "Florida", enddate = enddate)
getplot.recresid(FL.modeldata, "Florida", enddate = enddate)
getplot.detectors(FL.modeldata, "Florida", enddate = enddate)
getplot.newcases(NV.modeldata, "Nevada", enddate = enddate)
getplot.recresid(NV.modeldata, "Nevada", enddate = enddate)
getplot.detectors(NV.modeldata, "Nevada", enddate = enddate)
getplot.newcases(TX.modeldata, "Texas", enddate = enddate)
getplot.recresid(TX.modeldata, "Texas", enddate = enddate)
getplot.detectors(TX.modeldata, "Texas", enddate = enddate)
dev.off()
