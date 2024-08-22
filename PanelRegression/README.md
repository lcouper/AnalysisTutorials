# Panel Regression

The code and description below outlines the steps for conducting a regression using panel data (for leveraging repeated observations of a given unit over time as a way to address bias due to confounders). The basic steps will involve: 

1. Loading the panel data and libraries needed for analysis
2. Setting up the regression model
3. Visualizing output

The dataset used here ('MosquitoAbundance_PanelData') pertains to observations of various mosquito species across California counties from 2012 - 2020, as well as cases of dog heartworm for those counties and years. We are using a panel regression approach to investigate whether variation in the abundance of any particular mosquito species is associated with variation in dog heartworm cases.

#### 1. Load packages and data ####
```
library(data.table)
library(mltools)

setwd("~/PanelRegression")
data = read.csv("MosquitoAbundance_PanelData.csv", header =T)[,-1]
```
#### 2. Separate data into bioregions #####

Given broad climatic and ecological differences between California bioregions, we may want to run separate regressions for each. Here, we separate out data into Northern, Central, and Southern bioregions:

```
NorthBior = c("BayDelta", "Klamath", "Sierra", "SacramentoValley")
DataCentral =  DataA2[DataA2$Bioregion %in% CentralBior, ]
SouthBior = c("SouthCoast", "ColoradoDesert")

DataNorth = DataA2[DataA2$Bioregion %in% NorthBior, ]
DataCentral =  DataA2[DataA2$Bioregion %in% CentralBior, ]
DataSouth = DataA2[DataA2$Bioregion %in% SouthBior, ]
```

#### 3. Run panel mdoels for each bioregion ####

Notes: 
- In the model call,  we scale each predictor so that their effects are more directly comparable.
- We are using 1-year lagged version of the mosquito abundances to account for the lag between transmission and potential case detection
- To define the panel data structure, we set year and bioregion as fixed effects
- We also include dog density and county-level income as predictors in attempt to control for year-to-year and within-bioregion variation in host abundance and reporting tendencies (the year and bioregion fixed effects should account for yearly variation that affects all bioregions, and spatial variation that is constant over time within a given bioregion)

**Full Southern California panel model**
```
SoCalpm = lm(TotalPositive ~  scale(Lag_Ae.aegypti) + scale(Lag_Ae.albopictus) + scale(Lag_Ae.sierrensis) +
               scale(Lag_Ae.vexans) + scale(Lag_An.freeborni) + scale(Lag_Cs.incidencs) + scale(Lag_Cs.inornata) +
               scale(Lag_Cx.quinquefasciatus) + scale(Lag_Cx.tarsalis) + 
               scale(Lagged_DogDensity) + scale(Lagged_Income) +
               factor(Year) + factor(Bioregion) - 1, data = DataSouth)
summary(SoCalpm)
```

**Repeat for Central California**
CentralCalpm = lm(TotalPositive ~  scale(Lag_Ae.aegypti) + scale(Lag_Ae.albopictus) + scale(Lag_Ae.sierrensis) +
                scale(Lag_Ae.vexans) + scale(Lag_An.freeborni) + scale(Lag_Cs.incidencs) + scale(Lag_Cs.inornata) +
                scale(Lag_Cx.quinquefasciatus) + scale(Lag_Cx.tarsalis) + 
                scale(Lagged_DogDensity) + scale(Lagged_Income) +
                factor(Year) + factor(Bioregion) - 1, data = DataCentral)
summary(CentralCalpm)

**Repeat for Northern California**
```
NorCalpm = lm(TotalPositive ~  scale(Lag_Ae.aegypti) + scale(Lag_Ae.albopictus) + scale(Lag_Ae.sierrensis) +
                scale(Lag_Ae.vexans) + scale(Lag_An.freeborni) + scale(Lag_Cs.incidencs) + scale(Lag_Cs.inornata) +
                scale(Lag_Cx.quinquefasciatus) + scale(Lag_Cx.tarsalis) + 
                scale(Lagged_DogDensity) + scale(Lagged_Income) +
                factor(Year) + factor(Bioregion) - 1, data = DataNorth)
summary(NorCalpm)
```


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
