# Panel Regression

The code and description below outlines the steps for conducting a regression using panel data (for leveraging repeated observations of a given unit over time as a way to address bias due to confounders). The basic steps will involve: 

1. Loading the panel data and libraries needed for analysis
2. Setting up the regression model
3. Visualizing output

The dataset used here ('MosquitoAbundance_PanelData') pertains to observations of various mosquito species across California counties from 2012 - 2020.

#### 1. Load packages and data ####
library(data.table)
library(mltools)

DataA = read.csv("~/Documents/Current_Projects/DogHeartworm/Q1analysis_AllPredictors_MosqAbundance.csv", header =T)[,-1]
DataADT = data.table(DataA)
DataA2 = one_hot(DataADT)

# set NAs in mosquito abundance data to 0
DataA2$Lag_Ae.sierrensis[is.na(DataA2$Lag_Ae.sierrensis)] <- 0
DataA2$Lag_Ae.albopictus[is.na(DataA2$Lag_Ae.albopictus)] <- 0
DataA2$Lag_Ae.aegypti[is.na(DataA2$Lag_Ae.aegypti)] <- 0
DataA2$Lag_Ae.vexans[is.na(DataA2$Lag_Ae.vexans)] <- 0
DataA2$Lag_An.freeborni[is.na(DataA2$Lag_An.freeborni)] <- 0
DataA2$Lag_Cs.incidencs[is.na(DataA2$Lag_Cs.incidencs)] <- 0
DataA2$Lag_Cs.inornata[is.na(DataA2$Lag_Cs.inornata)] <- 0
DataA2$Lag_Cx.tarsalis[is.na(DataA2$Lag_Cx.tarsalis)] <- 0
DataA2$Lag_Cx.quinquefasciatus[is.na(DataA2$Lag_Cx.quinquefasciatus)] <- 0

#### 2. Separate Northern, Central, & Southern CA bioregions #####

NorthBior = c("BayDelta", "Klamath", "Sierra", "SacramentoValley")
DataCentral =  DataA2[DataA2$Bioregion %in% CentralBior, ]
SouthBior = c("SouthCoast", "ColoradoDesert")

DataNorth = DataA2[DataA2$Bioregion %in% NorthBior, ]
DataCentral =  DataA2[DataA2$Bioregion %in% CentralBior, ]
DataSouth = DataA2[DataA2$Bioregion %in% SouthBior, ]


#### 3. Run panel mdoels for Northern and Southern CA ####
# scale predictors so they are directly comparable

# Full NorCal Panel Model
NorCalpm = lm(TotalPositive ~  scale(Lag_Ae.aegypti) + scale(Lag_Ae.albopictus) + scale(Lag_Ae.sierrensis) +
                scale(Lag_Ae.vexans) + scale(Lag_An.freeborni) + scale(Lag_Cs.incidencs) + scale(Lag_Cs.inornata) +
                scale(Lag_Cx.quinquefasciatus) + scale(Lag_Cx.tarsalis) + 
                scale(Lagged_DogDensity) + scale(Lagged_Income) 
                factor(Year) + factor(Bioregion) - 1, data = DataNorth)
summary(NorCalpm)

# Full SoCal Panel model #
SoCalpm = lm(TotalPositive ~  scale(Lag_Ae.aegypti) + scale(Lag_Ae.albopictus) + scale(Lag_Ae.sierrensis) +
               scale(Lag_Ae.vexans) + scale(Lag_An.freeborni) + scale(Lag_Cs.incidencs) + scale(Lag_Cs.inornata) +
               scale(Lag_Cx.quinquefasciatus) + scale(Lag_Cx.tarsalis) + 
               scale(Lagged_DogDensity) + scale(Lagged_Income) 
               factor(Year) + factor(Bioregion) - 1, data = DataSouth)
summary(SoCalpm)

# Full CentralCal Panel Model
CentralCalpm = lm(TotalPositive ~  scale(Lag_Ae.aegypti) + scale(Lag_Ae.albopictus) + scale(Lag_Ae.sierrensis) +
                scale(Lag_Ae.vexans) + scale(Lag_An.freeborni) + scale(Lag_Cs.incidencs) + scale(Lag_Cs.inornata) +
                scale(Lag_Cx.quinquefasciatus) + scale(Lag_Cx.tarsalis) + 
                scale(Lagged_DogDensity) + scale(Lagged_Income) 
                factor(Year) + factor(Bioregion) - 1, data = DataCentral)
summary(CentralCalpm)

# NorCal fixed effects only
NorCalpmFEonly = lm(TotalPositive ~   
                      factor(Year) + factor(Bioregion) - 1, data = DataNorth)
summary(NorCalpmFEonly)

# SoCal fixed effects only
SoCalpmFEonly = lm(TotalPositive ~   
                      factor(Year) + factor(Bioregion) - 1, data = DataSouth)
summary(SoCalpmFEonly)

# CentralCal fixed effects only
CentralCalpmFEonly = lm(TotalPositive ~   
                      factor(Year) + factor(Bioregion) - 1, data = DataCentral)
summary(CentralCalpmFEonly)

# NorCal year fixed effect only
NorCalpmYearonly = lm(TotalPositive ~   
                        factor(Year) - 1, data = DataNorth)
summary(NorCalpmYearonly)

# SoCal year fixed effect only
SoCalpmYearonly = lm(TotalPositive ~   
                        factor(Year) - 1, data = DataSouth)
summary(SoCalpmYearonly)

# CentralCal year fixed effect only
CentralCalpmYearonly = lm(TotalPositive ~   
                        factor(Year) - 1, data = DataCentral)
summary(CentralCalpmYearonly)


# NorCal bioregion fixed effect only
NorCalpmBioronly = lm(TotalPositive ~   
                        factor(Bioregion) - 1, data = DataNorth)
summary(NorCalpmBioronly)

# SoCal bioregion fixed effect only
SoCalpmBioronly = lm(TotalPositive ~   
                        factor(Bioregion) - 1, data = DataSouth)
summary(SoCalpmBioronly)

# CentralCal bioregion fixed effect only
CentralCalpmBioronly = lm(TotalPositive ~   
                        factor(Bioregion) - 1, data = DataCentral)
summary(CentralCalpmBioronly)
