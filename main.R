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

data <- read.csv("C:/Users/grnma/OneDrive/Bureau/ENSAE/2A/S2/Séries temp/Projet Michelin/serie_010767848_01052024/tmpZipSerieCsv3541213619101186198/valeurs_mensuelles.csv", header=FALSE, sep=";")
str(data)

data <- data[5:nrow(data), 0:2]
colnames(data) <- c('Date', 'Valeur')
data$Date <- as.Date(paste(data$Date, "-01", sep = ""), format = "%Y-%m-%d")
