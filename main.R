#==============================================================================#
##### Packages ####
#==============================================================================#

library(zoo)
library(tseries)
library(fUnitRoots)
library(purrr)
library(urca)
library(astsa)
library(forecast)
library(portes)

#==============================================================================#
##### Import ####
#==============================================================================#

data <- read.csv("data/serie_010767848_01052024/tmpZipSerieCsv3541213619101186198/valeurs_mensuelles.csv", header=FALSE, sep=";")
str(data)

data <- data[5:nrow(data), 0:2]
colnames(data) <- c('Date', 'Valeur')
data$Valeur <- as.numeric(data$Valeur)
data$Date <- as.Date(paste(data$Date, "-01", sep = ""), format = "%Y-%m-%d")
data$Valeur <- rev(data$Valeur)

dates <- as.yearmon(seq(from=1990,to=2024+1/12,by=1/12))#Set of dates to plot the series
dates_diff <- as.yearmon(seq(from=1990,to=2024,by=1/12))#To plot the differenciated series

time_serie <- zoo(data[[2]])


#==============================================================================#
##### Question 2 ####
#==============================================================================#
#We consider the differenciated serie of order 1.
diff_ts <- time_serie - lag(time_serie, -1)

#Check for trends using LS
model <- lm(time_serie ~ dates)

summary(lm(time_serie ~ dates))#Strong evidence supporting trend
summary(lm(diff_ts ~ dates_diff))


#Unit root tests
#Augmented Dickey-Fuller test
#Necessary hypothesis to check : Residuals are uncorrelated
adf <- adfTest(time_serie, lag=0, type="ct")

#Qtests -> function to perform L.-B. tests on the residuals
Qtests <- function(series, k, fitdf=0) {
  pvals <- apply(matrix(0:k), 1, FUN=function(l) {
    pval <- if (l<=fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
    return(c("lag"=l,"pval"=pval))
  })
  return(t(pvals))
}
Qtests(adf@test$lm$residuals[0:24],length(adf@test$lm$coefficients))

#adfTest_valid -> function to check the lag needed for the ADF test to be performed
adfTest_valid <-function(series,kmax,type){
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

adfTest_valid(time_serie, 24, 'ct')#2 lags are needed for the ADF test
#P-value = 0.61 ; we cannot reject the trend stationnarity

adfTest_valid(diff_ts, 24, 'c')#1 lag is needed for the ADF test on the differenciated series
#P-value = 0.01 ; we now have evidence supporting stationnarity


#Using KPSS TESTS
ur.kpss(time_serie, type=c('tau'))#tau means we test the constant trend scenario
ur.kpss(time_serie, type=c('tau'))@cval
#We strongly reject trend_stationnarity (1.267 >> 0.216)

ur.kpss(diff_ts, type=c('mu'))#mu means we test the constant scenario
ur.kpss(diff_ts, type=c('mu'))@cval
#We cannot reject stationnarity (0.196 << 0.347)


#At this point, the differenciated series looks stationnary


#==============================================================================#
##### Question 3 ####
#==============================================================================#
#Here, we plot the series before and after differenciating it

png(file="graphs/initial_ts_plot.png",width=1200, height=700)
plot(x=dates, y=time_serie, type='l',
     xlab = 'Date', 
     ylab = "Index value (set as 100 in 2021)",
     main = 'Time serie plot from 1990 to 2024')
ticks <- seq(min(dates), max(dates), by = 2)
axis(1, at = ticks)
dev.off()


png(file="graphs/plot_de_la_serie_differenciee.png",width=1200, height=700)
plot(x=dates_diff, y=diff_ts, xlab = 'date', type='l',
     ylab = "valeur de X_t - X_t-1",
     main = 'plot de la série différenciée (ordre1) de 1990 à 2024'
     )
ticks <- seq(min(dates), max(dates), by = 2)
axis(1, at = ticks)
dev.off()


#Both on the same graph
png(file="graphs/comparison_plots.png", width=2400, height=700) # Increased width for two plots
par(mfrow=c(1, 2))

plot(x=dates, y=time_serie, type='l',
     xlab = 'Date', 
     ylab = "Index value (set as 100 in 2021)",
     main = 'Time serie plot from 1990 to 2024')
ticks <- seq(min(dates), max(dates), by = 2)
axis(1, at = ticks)

plot(x=dates_diff, y=diff_ts, type='l',
     xlab = 'Date', 
     ylab = "Value",
     main = 'Differentiated time serie plot from 1990 to 2024')
ticks <- seq(min(dates), max(dates), by = 2)
axis(1, at = ticks)

dev.off()
#==============================================================================#
##### Question 4 ####
#==============================================================================#
#First, let's pick an ARMA model for the differenciated series

stats::acf(diff_ts) ###MA(2)
stats::pacf(diff_ts) ###AR(2)

png(file="graphs/acf_pacf_diff_ts.png", width=1200, height=700)
par(mfrow=c(1, 2))
acf(diff_ts, main='ACF of differentiated time series')
pacf(diff_ts, main='PACF of differentiated time series')
dev.off()
#ACF and PACF steer us at testing MA(2), AR(2), and mixed ARMA models.

#ARMA(2,2)
res22 <- astsa::sarima(xdata = time_serie, p = 2, d = 1, q = 2, details=FALSE)
#seems a little bit too complicated as the p_values are high
#as the p_values of the MA part are the worst, we start by checking if there is a model with 
#better p_values going down on the MA part

#ARMA(2,1)
res21 <- astsa::sarima(xdata=time_serie,p=2,d=1,q=1, details=FALSE)
#as the coefficients are not that good, we might check if 
#it is not the AR part that was too complicated

#ARMA(1,2)
res12 <- astsa::sarima(xdata = time_serie,p=1,d=1,q=2, details=FALSE)
#it's a bit better but nothing from the other world
#we might also check the results for an arima(1,1,1)

#ARMA(1,1)
res11 <- astsa::sarima(xdata = time_serie, p=1,d=1,q=1, details=FALSE)
#Maybe the best fit of all for now, but still the AR1 coeff is not significant


#Now lets check for 'extreme cases' of an AR(2) or 1 and of a MA(2) or 1
res20 <- astsa::sarima(xdata = time_serie, p=2,d=1,q=0, details=FALSE)
#The fit is good

res02 <- astsa::sarima(xdata = time_serie, p=0,d=1,q=2, details=FALSE)
#these two are really not better than the ones before


#At this point, we still consider AR(2) and MA(2)


portes::LjungBox(obj=res20$fit)
portes::LjungBox(obj=res02$fit)

#p_values are also quite high, which indicate that we cannot reject the null hyp
#of absence of correlation

#To decide which model to keep, let's compare the ICs of each model

res20$ICs
res02$ICs

#Both AIC and BIC choose the AR(2) => this will be our model of interest

#==============================================================================#
##### Question 5 ####
#==============================================================================#
#Now that we know that the best fit for the differenciated series is an AR(2)
#We can express our model as an ARIMA(2,1,0)

res20$fit$coef
res20$ttable
-1.96*sqrt(var(res20$fit$residuals))/length(res20$fit$residuals)


#==============================================================================#
##### Question 6-9 ####
#==============================================================================#
#Let's represent the confidence interval for the future values X_T+1 and X_T+2

ts_data <- ts(data = data$Valeur)

astsa::sarima.for(xdata = ts_data, n.ahead = 2, p = 2,d = 1,q = 0,
                  main='Prediction of future values at T+1 and T+2')

#The sarima.for function plots confidence intervals of +- 1*se and +- 2*se 
#Let's plot our own confidence intervals
#Refer to the report to have the general idea behind the construction of these intervals

prediction <- astsa::sarima.for(xdata = ts_data, n.ahead = 2, p = 2,d = 1,q = 0,
                         main='Prediction of future values at T+1 and T+2')
prediction$pred
est_se <- sqrt(res20$fit$sigma2)


#For T+1
first_upper_bound <- prediction$pred[1] + 1.96*est_se
first_lower_bound <- prediction$pred[1] - 1.96*est_se

#For T+2
sec_upper_bound <- prediction$pred[2] + 1.96*est_se*sqrt(1+res20$fit$coef[1]**2)
sec_lower_bound <- prediction$pred[2] - 1.96*est_se*sqrt(1+res20$fit$coef[1]**2)


#We keep the last 50 observations
df_pred <- data.frame(ts_data, dates)[311:410 ,]


new_dates <- as.yearmon(seq(from=2024+2/12,to=2024+3/12,by=1/12))
prediction_df <- data.frame(
  ts_data = c(prediction$pred[1], prediction$pred[2]),
  dates = new_dates
  )

#We assign new columns to help us plot
df_pred <- rbind(df_pred, prediction_df)
df_pred$is_last_two <- FALSE
df_pred$is_last_two[(nrow(df_pred)-1):nrow(df_pred)] <- TRUE

df_pred$upper_bound <- NA
df_pred$lower_bound <- NA

df_pred$upper_bound[(nrow(df_pred)-1)] <- first_upper_bound
df_pred$lower_bound[(nrow(df_pred)-1)] <- first_lower_bound
df_pred$upper_bound[nrow(df_pred)] <- sec_upper_bound
df_pred$lower_bound[nrow(df_pred)] <- sec_lower_bound

#The final plot !!
png(file="graphs/predictons_plot.png",width=1200, height=700)
ggplot(data=df_pred, aes(x=dates, y=ts_data)) +
  geom_line() +
  geom_point(aes(color=is_last_two), size=0.9) +
  geom_ribbon(data=df_pred[df_pred$is_last_two, ], aes(ymin=lower_bound, ymax=upper_bound), fill='red', alpha=0.2) +
  scale_color_manual(values = c('FALSE' = 'black', 'TRUE' = 'red')) +
  theme_minimal() +
  labs(title="T+1 and T+2 predictions and associated confidence bounds") +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  guides(color = "none")
dev.off()









