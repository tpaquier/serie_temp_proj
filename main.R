#==============================================================================#
##### Packages ####
#==============================================================================#

library(zoo)
library(tseries)
library(fUnitRoots)
library(purrr)
library(astsa)
library(forecast)
library(portes)

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


adf <- adfTest(ts, lag=2, type="ct")

###Graphiques--------------------------------------------------------

png(file="graphs/plot_initial_de_la_serie.png",width=600, height=350)
plot(x=dates, y=ts, type='l',
     xlab = 'date', 
     ylab = "valeur de l'indice (base 100 en 2021)",
     main = 'plot de la série de 1990 à 2024')
ticks <- seq(min(dates), max(dates), by = 2)
axis(1, at = ticks)
dev.off()


png(file="graphs/plot_de_la_serie_differenciee.png",width=600, height=350)
plot(x=dates, y=diff_ts, xlab = 'date', 
     ylab = "valeur de X_t - X_t-1",
     main = 'plot de la série différenciée (ordre1) de 1990 à 2024')
ticks <- seq(min(dates), max(dates), by = 2)
axis(1, at = ticks)
dev.off()

#Tests------------------------------------------------------------

stats::acf(ts)
stats::pacf(ts)

stats::acf(diff_ts) ###MA(2)
stats::pacf(diff_ts) ###AR(2)

res22 <- astsa::sarima(xdata = ts, p = 2, d = 1, q = 2) #on a pas besoin de prendre 
#la diff pcq qu'on colle le lag dedans
res22$fit

#seems a little bit too complicated as the p_values are high
#as the p_values of the MA part are the worst, we start by checking if there is a model with 
#better p_values going down on the MA part

res21 <- astsa::sarima(xdata=ts,p=2,d=1,q=1)
res21$fit

#as the coefficients are not that good, we might check if 
#it is not the AR part that was too complicated

res12 <- astsa::sarima(xdata = ts,p=1,d=1,q=2)
res12$fit

#well, it's a bit better but nothing from the other world
#we might also check the results for an arima(1,1,1)

res11 <- astsa::sarima(xdata = ts, p=1,d=1,q=1)
res11$fit

#the model is, at first glance, better
#we now can do checks to see how it is 'really'

res11$ttable
res12$ttable
res21$ttable
res22$ttable


portes::LjungBox(obj = res11$fit)
#p_values are quite high, we thus cannot reject the fact that this model might 
#be good
















