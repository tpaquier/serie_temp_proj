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
dates <- as.yearmon(seq(from=1990,to=2024+1/12,by=1/12))
dates_diff <- as.yearmon(seq(from=1990,to=2024,by=1/12))

ts <- zoo(data[[2]])
diff_ts <- ts - lag(ts, -1)

#Check for trends using LS
summary(lm(ts ~ dates))


res_adf <- adfTest(diff_ts, lag=0, type="nc")

png(file="graphs/plot_initial_de_la_serie.png",width=600, height=350)
plot(x=dates, y=ts, type='l',
     xlab = 'date', 
     ylab = "valeur de l'indice (base 100 en 2021)",
     main = 'plot de la série de 1990 à 2024')
ticks <- seq(min(dates), max(dates), by = 2)
axis(1, at = ticks)
dev.off()


png(file="graphs/plot_de_la_serie_differenciee.png",width=600, height=350)
plot(ts, xlab = 'date', 
     ylab = "valeur de X_t - X_t-1",
     main = 'plot de la série différenciée (ordre1) de 1990 à 2024')
ticks <- seq(min(dates), max(dates), by = 2)
axis(1, at = ticks)
dev.off()

plot(x = dates, y = ts, type='l')

