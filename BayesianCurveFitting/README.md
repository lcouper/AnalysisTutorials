
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

![LDR fit](./Figures/LDR_Hopland_Uniform_Fit.jpeg)

The raw experimental data appears left-skewed, which is typical for rate traits of ectotherms such as development, fecundity, and biting. These are well described by Briere functions: 
$$cT(T-T0)\sqrt{Tm−T}$$ 
where T0 and Tm are the critical thermal minimum and maximum, respectively and c is a positive rate constant.

The other commonly used function is the quadratic, which is symmetrical and better describes probability traits such as survival at a given life stage. The quadratic function is given as:
$$-c(T–T0)(T–Tm))$$

### 3. Setting up the JAGS model

Next we set-up the JAGS ('Just Another Gibbs Sampler') model, which we will use to run Bayesian models using simulation.
We are using this approach to estimate T0, Tm, and c (the critical thermal minimum, maximum, and rate constant shown in the equation above, and predicted trait values, needed to generate continuous thermal response curves.

In the set-up, we: 

1) Specify the 'priors' -- a probability distribution capturing our beliefs about what the parameters may be, using our knowledge of the system. In the first round of model fitting, we are using 'low information' priors. These are uniform priors, bounded by biologically realistic temperature constraints. For T0, we set these bounds to be 0-20, meaning the probability distribution will only take values in this range, and all values within these bounds are equally likely. For Tm, we set the bounds as 28-35. These bounds were based on prior experiments and information from our mosquito species. We also specify sigma and tau, which expresses the standard deviation used when making trait predictions.
2) Specify the likelihood function, which specifies how the raw data link to the parameters you are trying to estimate. Here, we use the Briere function from above, and set trait performance to 0 if  T0 > T > Tm.
3) Specify any additional quantities we want the model to estimate, which for us are the predicted trait values at other temperatures. This will enable us to generate a continuous thermal performance curve, from experimental trait data measured at 6 temperatures.

```{r, results = 'hide'}
{
  setwd("~/Documents/Current_Projects/LifeHistoryTraitExp/Analysis_TraitFits")
  sink("briere.txt") # This drives the text below to a text file, which will be saved in your working directory
cat("
    model{  
    
    # 1. Specify the Priors
    cf.q ~ dunif(0, 1) # c = a positive rate constant
    cf.T0 ~ dunif(0, 20) 
    cf.Tm ~ dunif(28, 35) # 35 set based on highest CTmax for larval survival
    cf.sigma ~ dunif(0, 1000) # standard deviation 
    cf.tau <- 1 / (cf.sigma * cf.sigma)
    
    # 2. Specify the Likelihood function (here, Briere)
    for(i in 1:N.obs){
    trait.mu[i] <- cf.q * temp[i] * (temp[i] - cf.T0) * sqrt((cf.Tm - temp[i]) * (cf.Tm > temp[i])) * (cf.T0 < temp[i])
    trait[i] ~ dnorm(trait.mu[i], cf.tau) # trait values drawn from normal dist with variance from ct tau
    }
    
    # 3. Specify the Derived Quantities and Predictions
    for(i in 1:N.Temp.xs){
    z.trait.mu.pred[i] <- cf.q * Temp.xs[i] * (Temp.xs[i] - cf.T0) * sqrt((cf.Tm - Temp.xs[i]) * (cf.Tm > Temp.xs[i])) * (cf.T0 < Temp.xs[i])
    }
    
    } # close model
    ",fill=T)
sink()
}
```

### 4. Setting up the MCMC simulation

Next we set-up the Markov Chain Monte Carlo sampler. This is a way to *approximate* the posterior distribution when you can't solve for it analytically (i.e., because the shape of the prior and likelihood distributions are funky and can't be easily combined). The way the simulation works is it start with some initial parameter value, for example '32$^\circ$C' for Tm. Then it picks another random parameter value, compares it with the previous one, and keeps whichever has the *higher likelihood* of being in the posterior distribution. This value gets added to the 'chain', and the process is repeated for as many iterations as you specify. Eventually your chain will converge on a single parameter value. We 'burn-in', or drop from the chain, some number of initial entries in the chain, as these are more heavily influenced by our choice of initial value, which we don't want. We can also 'thin' our chain by removing every xth value. This can help reduce the autocorrelation of points in the chain, since each point depends on the prior point. You can also repeat this process several times (i.e., build multiple chains), which will help with convergence and enable more precise parameter estimates.

```{r, results = 'hide'}
# For each parameter, we specify the initial values to use in the MCMC simulation
inits<-function(){list(
  cf.q = 0.01,
  cf.Tm = 35,
  cf.T0 = 5,
  cf.sigma = rlnorm(1))}

# Parameters to Estimate 
parameters <- c("cf.q", "cf.T0", "cf.Tm","cf.sigma", "z.trait.mu.pred")

# Specify the MCMC settings 
ni <- 25000 # number of iterations in each chain
nb <- 5000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nth iterations in each chain
nc <- 3 # number of chains

# Temperature sequence for derived quantity calculations
# This enables us to build the continuous thermal performance curve
Temp.xs <- seq(1, 45, 0.1)
N.Temp.xs <-length(Temp.xs)
# For priors - fewer temperatures for derived calculations makes it go faster
Temp.xs <- seq(5, 45, 0.5)
N.Temp.xs <-length(Temp.xs)
```
## 5. For each population, fit the model using low information priors
```{r, warning = FALSE, results = 'hide'}
# Below, we fit the model for a single population -- Hopland, or 'HOP'

# We first organize the data for the JAGS model
data <- data.LDR.HOP
trait <- data$LarvalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Then we run the JAGS model
LDR.HOP.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())
```
#### Now we can examine the output...

The output summary shows the mean, standard deviation, and quartiles for each of the specified parameters.
You can also view values from each iteration of the chain and model diagnostics
```{r, warning = FALSE}
LDR.HOP.out$BUGSoutput$summary[1:5,]
# View (LDR.HOP.out$BUGSoutput) # To look at values from each iteration
# mcmcplot(LDR.HOP.prior.out) # To view model diagnostic plots (e.g., convergence, autocorrelation)
```

#### ...and plot the raw data with the model fits. Here we're plotting the mean, 2.5%, and 97.5% quantiles for the predicted trait values across temperature
```{r, fig = TRUE}
plot(LarvalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.12), data = data.LDR.HOP, ylab = "LDR for Hop", xlab = "Temperature")
lines(LDR.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(LDR.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
```


The model fit appears to describe the raw data well for this trait and population.

The above code (starting from "data <- data.LDR.HOP") was just for the 'HOP' population. Try modifying the code to repeat the process for the nine remaining populations: MAR1, MAR2, WAW, EUG, PLA, SB, JRA, PAR, and POW. (Below is the code to do so, if needed)

```{r, warning = FALSE, results = 'hide'}
# MAR1
data <- data.LDR.MAR1 
trait <- data$LarvalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)
LDR.MAR1.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                     n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# MAR2
data <- data.LDR.MAR2 
trait <- data$LarvalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)
LDR.MAR2.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                     n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# WAW
data <- data.LDR.WAW 
trait <- data$LarvalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)
LDR.WAW.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# EUG
data <- data.LDR.EUG 
trait <- data$LarvalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)
LDR.EUG.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# PLA
data <- data.LDR.PLA 
trait <- data$LarvalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)
LDR.PLA.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# SB
trait <- data$LarvalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)
LDR.SB.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                   n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# JRA
data <- data.LDR.JRA 
trait <- data$LarvalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)
LDR.JRA.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# PAR 
trait <- data$LarvalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)
LDR.PAR.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# POW 
data <- data.LDR.POW 
trait <- data$LarvalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)
LDR.POW.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())
```

Lets visually check the fits for each population by plotting the raw data and model fits :
```{r, echo = FALSE, warning = FALSE}

par(mfrow = c(5,2))
par(mar = c(2.5, 2.5, 1, 1))

plot(LarvalDevRate ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.12), data = data.LDR.HOP, ylab = "LDR for Hop", xlab = "Temperature")
lines(LDR.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(LDR.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
text("HOP", x = 36.5, y = 0.105, cex = 1.5)

plot(LarvalDevRate ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.12), data = data.LDR.MAR1, ylab = "LDR for MAR1", xlab = "Temperature")
lines(LDR.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(LDR.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
text("MAR1", x = 36.5, y = 0.105, cex = 1.5)

plot(LarvalDevRate ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.12), data = data.LDR.MAR2 , ylab = "LDR for MAR2 prior", xlab = "Temperature")
lines(LDR.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(LDR.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
text("MAR2", x = 36.5, y = 0.105, cex = 1.5)

plot(LarvalDevRate ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.12), data = data.LDR.WAW , ylab = "LDR for WAW prior", xlab = "Temperature")
lines(LDR.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(LDR.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
text("WAW", x = 36.5, y = 0.105, cex = 1.5)

plot(LarvalDevRate ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.12), data = data.LDR.EUG , ylab = "LDR for EUG prior", xlab = "Temperature")
lines(LDR.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(LDR.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
text("EUG", x = 36.5, y = 0.105, cex = 1.5)

plot(LarvalDevRate ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.12), data = data.LDR.PLA , ylab = "LDR for PLA prior", xlab = "Temperature")
lines(LDR.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(LDR.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
text("PLA", x = 36.5, y = 0.105, cex = 1.5)

plot(LarvalDevRate ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.12), data = data.LDR.SB , ylab = "LDR for SB prior", xlab = "Temperature")
lines(LDR.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(LDR.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
text("SB", x = 36.5, y = 0.105, cex = 1.5)

plot(LarvalDevRate ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.12), data = data.LDR.JRA , ylab = "LDR for JRA prior", xlab = "Temperature")
lines(LDR.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(LDR.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
text("JRA", x = 36.5, y = 0.105, cex = 1.5)

plot(LarvalDevRate ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.12), data = data.LDR.PAR , ylab = "LDR for PAR prior", xlab = "Temperature")
lines(LDR.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(LDR.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
text("PAR", x = 36.5, y = 0.105, cex = 1.5)

plot(LarvalDevRate ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.12), data = data.LDR.POW , ylab = "LDR for POW prior", xlab = "Temperature")
lines(LDR.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(LDR.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
text("POW", x = 36.5, y = 0.105, cex = 1.5)
```
