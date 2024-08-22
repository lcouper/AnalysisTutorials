#### Script to accompany panel regression tutorial #####

##### 1. Load packages and data #####


setwd("~/Downloads/PanelRegression")

library(data.table)
library(mltools)
library(regclass)

DataA = read.csv("MosquitoAbundance_PanelData.csv", header =T)[,-1]

##### 2. Test for multicollinearity in predictors #####

#Calculate the variance inflation factor (VIF) to assess multicollinearity of mosquito abundacne predictors:
  
testformulti <- lm(TotalPositive ~  Lag_Ae.aegypti + Lag_Ae.albopictus + Lag_Ae.sierrensis +
                       Lag_Ae.vexans + Lag_An.freeborni + Lag_Cs.incidencs + Lag_Cs.inornata +
                       Lag_Cx.quinquefasciatus + Lag_Cx.tarsalis, data = DataA)
VIF(testformulti)


##### 3. Separate out data based on bioregion #####

NorthBior = c("BayDelta", "Klamath", "Sierra", "SacramentoValley")
CentralBior = c("SanJoaquinValley", "CentralCoast")
SouthBior = c("SouthCoast", "ColoradoDesert")

DataNorth = DataA[DataA$Bioregion %in% NorthBior, ]
DataCentral =  DataA[DataA$Bioregion %in% CentralBior, ]
DataSouth = DataA[DataA$Bioregion %in% SouthBior, ]

##### 4. Run panel models for each region #####

SoCalpm = lm(TotalPositive ~  scale(Lag_Ae.aegypti) + scale(Lag_Ae.albopictus) + scale(Lag_Ae.sierrensis) +
               scale(Lag_Ae.vexans) + scale(Lag_An.freeborni) + scale(Lag_Cs.incidencs) + scale(Lag_Cs.inornata) +
               scale(Lag_Cx.quinquefasciatus) + scale(Lag_Cx.tarsalis) + 
               scale(Lagged_DogDensity) + scale(Lagged_Income) +
               factor(Year) + factor(Bioregion) - 1, data = DataSouth)
summary(SoCalpm)

CentralCalpm = lm(TotalPositive ~  scale(Lag_Ae.aegypti)  + scale(Lag_Ae.sierrensis) +
                    scale(Lag_Ae.vexans) + scale(Lag_An.freeborni) + scale(Lag_Cs.incidencs) + scale(Lag_Cs.inornata) +
                    scale(Lag_Cx.quinquefasciatus) + scale(Lag_Cx.tarsalis) + 
                    scale(Lagged_DogDensity) + scale(Lagged_Income) +
                    factor(Year) + factor(Bioregion) - 1, data = DataCentral)
summary(CentralCalpm)

NorCalpm = lm(TotalPositive ~  scale(Lag_Ae.aegypti) + scale(Lag_Ae.albopictus) + scale(Lag_Ae.sierrensis) +
                scale(Lag_Ae.vexans) + scale(Lag_An.freeborni) + scale(Lag_Cs.incidencs) + scale(Lag_Cs.inornata) +
                scale(Lag_Cx.quinquefasciatus) + scale(Lag_Cx.tarsalis) + 
                scale(Lagged_DogDensity) + scale(Lagged_Income) +
                factor(Year) + factor(Bioregion) - 1, data = DataNorth)
summary(NorCalpm)

#### 5. Repeat with fixed effects only ####

#To understand how much of the variation in dog heartworm cases may be explained by the predictors of interest (i.e., mosquito abundance, dog density, income) we repeat the model but including only the year and bioregion fixed effets:
  
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

