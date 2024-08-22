# Panel Regression

The code and description below outlines the steps for conducting a regression using panel data (for leveraging repeated observations of a given unit over time as a way to address bias due to confounders). The basic steps will involve: 

1. Loading the panel data and libraries needed for analysis
2. Setting up the regression model
3. Visualizing output

The dataset used here ('MosquitoAbundance_PanelData') pertains to observations of various mosquito species across California counties from 2012 - 2020, as well as cases of dog heartworm for those counties and years. We are using a panel regression approach to investigate whether variation in the abundance of any particular mosquito species is associated with variation in dog heartworm cases.

This tutorial accompanies the manuscript "Ecological drivers of dog heartworm transmission in California" available [here.](https://link.springer.com/article/10.1186/s13071-022-05526-x#Sec13)

#### 1. Load packages and data ####
```
library(data.table)
library(mltools)
library(regclass)

setwd("~/PanelRegression")
data = read.csv("MosquitoAbundance_PanelData.csv", header =T)[,-1]
```
#### 2. Test for multicollinearity in predictors #####

Calculate the variance inflation factor (VIF) to assess multicollinearity of mosquito abundacne predictors:
```
testformulti <- lm(TotalPositive ~  Lag_Ae.aegypti + Lag_Ae.albopictus + Lag_Ae.sierrensis +
                   Lag_Ae.vexans + Lag_An.freeborni + Lag_Cs.incidencs + Lag_Cs.inornata +
                   Lag_Cx.quinquefasciatus + Lag_Cx.tarsalis, data = DataA)
VIF(testformulti)
```
None of the mosquito predictors have a VIF value > 3. Ok to proceed without excluding any species from the model 

#### 3. Separate out data based on bioregion #####

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
- To define the panel data structure, we set year and bioregion as fixed effects (or 'dummy variables'). This should help control for any unobserved heterogeneity that might influence dog heartworm cases in a particular bioregion across all years (e.g. geographic features, number of veterinary clinics) or influence cases in all bioregions in a given year (e.g. an influx of shelter dogs to the state due to natural disaster, higher case reporting). 
- We also include dog density and county-level income as predictors in attempt to control for variation across time within a given bioregion in host abundance and reporting tendencies.

**Full Southern California panel model**
```
SoCalpm = lm(TotalPositive ~  scale(Lag_Ae.aegypti) + scale(Lag_Ae.albopictus) + scale(Lag_Ae.sierrensis) +
               scale(Lag_Ae.vexans) + scale(Lag_An.freeborni) + scale(Lag_Cs.incidencs) + scale(Lag_Cs.inornata) +
               scale(Lag_Cx.quinquefasciatus) + scale(Lag_Cx.tarsalis) + 
               scale(Lagged_DogDensity) + scale(Lagged_Income) +
               factor(Year) + factor(Bioregion) - 1, data = DataSouth)
summary(SoCalpm)
```

<img width="650" alt="image" src="https://github.com/user-attachments/assets/b92454eb-9636-4b95-86a6-13e9a89b0f65">

The scaled coefficient estimates shown here denote the change in dog heartworm cases from a one standard deviation change in mosquito abundances.

**Repeat for Central California**
Note: Aedes albopictus was not included as a predictor in the model for Central California as it was not found in either Central California bioregion in any year.

```
CentralCalpm = lm(TotalPositive ~  scale(Lag_Ae.aegypti)  + scale(Lag_Ae.sierrensis) +
                scale(Lag_Ae.vexans) + scale(Lag_An.freeborni) + scale(Lag_Cs.incidencs) + scale(Lag_Cs.inornata) +
                scale(Lag_Cx.quinquefasciatus) + scale(Lag_Cx.tarsalis) + 
                scale(Lagged_DogDensity) + scale(Lagged_Income) +
                factor(Year) + factor(Bioregion) - 1, data = DataCentral)
summary(CentralCalpm)
```

**Repeat for Northern California**
```
NorCalpm = lm(TotalPositive ~  scale(Lag_Ae.aegypti) + scale(Lag_Ae.albopictus) + scale(Lag_Ae.sierrensis) +
                scale(Lag_Ae.vexans) + scale(Lag_An.freeborni) + scale(Lag_Cs.incidencs) + scale(Lag_Cs.inornata) +
                scale(Lag_Cx.quinquefasciatus) + scale(Lag_Cx.tarsalis) + 
                scale(Lagged_DogDensity) + scale(Lagged_Income) +
                factor(Year) + factor(Bioregion) - 1, data = DataNorth)
summary(NorCalpm)
```

To understand how much of the variation in dog heartworm cases may be explained by the predictors of interest (i.e., mosquito abundance, dog density, income) we repeat the model but including only the year and bioregion fixed effets:

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
```

Coefficient estimates across models:

![image](https://github.com/user-attachments/assets/4185386f-ed50-4a1e-b924-904beea3ba80)

#### 4. Interpretations ####

Given the results generated above, we may infer: 



