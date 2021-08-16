## ####################################################################
## ####################################################################
## Supplement for
## "Backward CUSUM for Testing and Monitoring Structural Change"
## by Sven Otto and Jörg Breitung.
## This R-script allows to reproduce the data for Figure 5.
## ####################################################################
## ####################################################################
rm(list=ls())
## ##################################
## Load country data
countrydata.raw <- read.csv(url("https://data.humdata.org/hxlproxy/api/data-preview.csv?url=https%3A%2F%2Fraw.githubusercontent.com%2FCSSEGISandData%2FCOVID-19%2Fmaster%2Fcsse_covid_19_data%2Fcsse_covid_19_time_series%2Ftime_series_covid19_confirmed_global.csv&filename=time_series_covid19_confirmed_global.csv"))
countrydata <- data.frame(t(countrydata.raw[,-c(1:4)]))
colnames(countrydata) <- countrydata.raw[,2]
##
getdata.country <- function(country, from = "2020-01-22", to = "2021-06-09"){
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
  return(data[which(rownames(data) == from):which(rownames(data) == to),])
}
USdata = getdata.country("US")
USdata
