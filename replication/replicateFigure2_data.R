## ####################################################################
## ####################################################################
## Supplement for
## "Backward CUSUM for Testing and Monitoring Structural Change"
## by Sven Otto and Jörg Breitung.
## This R-script allows to reproduce the data for Figure 2.
## ####################################################################
## ####################################################################
rm(list=ls())
## ##################################
## Load country data
countrydata.raw <- read.csv(url("https://data.humdata.org/hxlproxy/api/data-preview.csv?url=https%3A%2F%2Fraw.githubusercontent.com%2FCSSEGISandData%2FCOVID-19%2Fmaster%2Fcsse_covid_19_data%2Fcsse_covid_19_time_series%2Ftime_series_covid19_confirmed_global.csv&filename=time_series_covid19_confirmed_global.csv"))
countrydata <- data.frame(t(countrydata.raw[,-c(1:4)]))
colnames(countrydata) <- countrydata.raw[,2]
##
getdata.country <- function(country){
  tot <- countrydata[[country]]
  startdate <- as.Date("2020-01-22")
  alldates <- seq(startdate, startdate + length(tot) - 1, by="day")
  negativevalues <- which(diff(tot) < 0)
  for(i in negativevalues){ 
    tot[i] <- tot[i+1]
    for(j in 1:i) ( if(tot[j] > tot[i]) ( tot[j] <- tot[i] ) )
  }
  new <- c(NA,diff(tot))
  data <- data.frame(tot, new)
  rownames(data) <- alldates
  colnames(data) <- c("total", "newcases")
  return(data)
}
##
getdata.usstate <- function(acronym){
  url <- paste("https://covidtracking.com/api/v1/states/",acronym,"/daily.csv", sep="")
  df <- read.csv(url)
  negativevalues <- which(diff(df$positive) > 0)
  for(i in negativevalues){
    for(j in i:dim(df)[1]) ( if(df$positive[j] > df$positive[i]) ( df$positive[j] <- df$positive[i] ) )
  }
  Rdate <- as.Date(paste(substring(df$date,1,4),
                         substring(df$date,5,6),
                         substring(df$date,7,8)), tryFormats="%Y %m %d")
  startdate <- as.Date("2020-01-22")
  alldates <- seq(startdate, max(Rdate), by="day")
  tot <- numeric(length(alldates))
  for(i in 1:length(alldates)){
    tot[i] <- sum(df$positive[which(Rdate == alldates[i])])
  }
  new <- c(NA,diff(tot))
  data <- data.frame(tot, new)
  rownames(data) <- alldates
  colnames(data) <- c("total", "newcases")
  return(data)
}
##
US.data <- getdata.country("US")
AZ.data <- getdata.usstate("az")
FL.data <- getdata.usstate("fl")
NV.data <- getdata.usstate("nv")
TX.data <- getdata.usstate("tx")
