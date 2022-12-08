
# Tutorial on Bayesian Curve Fitting 

```{r, setup, include=FALSE}
knitr::opts_knit$set(root.dir = '~/Documents/Current_Projects/LifeHistoryTraitExp/Analysis_TraitFits')
```

## Overview

In this tutorial, we will go through the steps involved in fitting Bayesian models to mosquito thermal performance curves. Specifically, we will fit temperature response functions to larval mosquito development rate ("LDR") for 10 populations of Aedes sierrensis, the western tree hole mosquito. The basic steps will involve:

1) Loading the data and libraries needed for analysis
2) Plotting the raw data to determine the appropriate functional form to fit
3) Setting up the JAGS model
4) Setting up the MCMC simulation
5) Fitting using low information priors (uniform priors bounded by biologically realistic constraints)
6) Generating informative priors using a leave-one-out approach
7) Fitting using these informative priors
8) Visually checking the model fits
9) Calculating additional thermal performance characteristics
10) Plotting thermal response curves and parameter estimates 

The data used in this tutorial can be found [here.](https://github.com/lcouper/AnalysisTutorials/blob/main/BayesianCurveFitting/LifeHistoryTraitExp_TutorialData.csv) 

Note this tutorial builds on code and approaches used by [Shocket et al. 2020](https://elifesciences.org/articles/58511) and [Mordecai et al. 2013](https://onlinelibrary.wiley.com/doi/abs/10.1111/ele.12015)[& 2017](https://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0005568). For more background on mosquito thermal biology and Bayesian curve fitting, see: [Mordecai et al. 2019. Ecology Letters.](https://onlinelibrary.wiley.com/doi/full/10.1111/ele.13335)

### 1. Loading data and libraries 

```{r, warning = FALSE, message = FALSE}
# Set working directory (Note: this will differ for you. Just make sure it matches where you downloaded the data)
setwd("~/Documents/Current_Projects/LifeHistoryTraitExp/Analysis_TraitFits")

# Load libraries for fitting traits
library('R2jags')
library('mcmcplots')
library('MASS')
library('HDInterval')

# Load the data and remove any missing values
data.LDR = read.csv("LifeHistoryTraitExp_Data.csv")[,-1]
data.LDR = data.LDR[!is.na(data.LDR$LarvalDevRate),]

# Subset the data by mosquito population. Here, 'population' refers to a single tree hole from which Aedes sierrensis larvae were collected. Populations used here ranged from Southern California to Eugene, Orgeon

data.LDR.HOP <- subset(data.LDR, Population == "HOP")
data.LDR.MAR1 <- subset(data.LDR, Population == "MARIN35")
data.LDR.MAR2 <- subset(data.LDR, Population == "MARIN29")
data.LDR.WAW <- subset(data.LDR, Population == "WAW")
data.LDR.EUG <- subset(data.LDR, Population == "EUG")
data.LDR.PLA <- subset(data.LDR, Population == "PLA")
data.LDR.SB <- subset(data.LDR, Population == "SB")
data.LDR.JRA <- subset(data.LDR, Population == "JRA")
data.LDR.PAR <- subset(data.LDR, Population == "PAR")
data.LDR.POW <- subset(data.LDR, Population == "POW")
```
### 2. Plot raw trait data

```{r}
plot(LarvalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.15), pch = 16,
     data = data.LDR, ylab = "Larval dev rate", xlab = "Temperature")
```

The raw experimental data appears left-skewed, which is typical for rate traits of ectotherms such as development, fecundity, and biting. These are well described by Briere functions: 
$$cT(T-T0)\sqrt{Tm−T}$$ 
where T0 and Tm are the critical thermal minimum and maximum, respectively and c is a positive rate constant.

The other commonly used function is the quadratic, which is symmetrical and better describes probability traits such as survival at a given life stage. The quadratic function is given as:
$$-c(T–T0)(T–Tm))$$
