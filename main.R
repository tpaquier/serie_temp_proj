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
summary(lm(diff_ts ~ dates_2))

res_adf <- adfTest(diff_ts, lag=0, type="nc")


