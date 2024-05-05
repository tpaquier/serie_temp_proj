#==============================================================================#
##### Packages ####
#==============================================================================#

library(zoo)
library(tseries)
library(fUnitRoots)
library(purrr)


#==============================================================================#
##### Import + Traitement ####
#==============================================================================#

data <- read.csv("data/serie_010767848_01052024/tmpZipSerieCsv3541213619101186198/valeurs_mensuelles.csv", header=FALSE, sep=";")
str(data)

data <- data[5:nrow(data), 0:2]
colnames(data) <- c('Date', 'Valeur')
data$Valeur <- as.numeric(data$Valeur)
data$Date <- as.Date(paste(data$Date, "-01", sep = ""), format = "%Y-%m-%d")
data$Valeur <- rev(data$Valeur)
dates <- as.yearmon(seq(from=1990,to=2024+1/12,by=1/12))
dates_diff <- as.yearmon(seq(from=1990,to=2024,by=1/12))

ts <- zoo(data[[2]])
diff_ts <- ts - lag(ts, -1)

#Check for trends using LS
summary(lm(ts ~ dates))#Strong evidence supporting trend but need to check autocorr first
summary(lm(diff_ts ~ dates_2))

#Unit root tests
#ADF
#Hypothesis to check : Residuals are uncorrelated
adf <- adfTest(ts, lag=0, type="ct")

Qtests <- function(series, k, fitdf=0) {
  pvals <- apply(matrix(1:k), 1, FUN=function(l) {
    pval <- if (l<=fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
    return(c("lag"=l,"pval"=pval))
  })
  return(t(pvals))
}

Qtests(adf@test$lm$residuals[0:24],length(adf@test$lm$coefficients))


adfTest_valid <-function(series,kmax,type){ #ADF tests until no more autocorrelated residuals
  k <- 0
  noautocorr <- 0
  while (noautocorr==0){
    cat(paste0("ADF with ",k, " lags: residuals OK? "))
    adf <- adfTest(series,lags=k,type=type)
    pvals <- Qtests(adf@test$lm$residuals,24,fitdf=length(adf@test$lm$coefficients))[,2]
    if (sum(pvals<0.05,na.rm=T) == 0) {
      noautocorr <- 1; cat("OK \n")}
    else cat("nope \n")
    k <- k + 1
  }
  return(adf)
}

adf <- adfTest_valid(ts,24,"ct")
adf <- adfTest(ts, lag=2, type="ct")
