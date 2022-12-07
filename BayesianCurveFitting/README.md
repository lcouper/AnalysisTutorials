
## Tutorial on Bayesian Curve Fitting 
This tutorial go through the steps involved in fitting Bayesian models to mosquito thermal performance curves. Specifically, we will fit temperature response functions to larval mosquito development rate ("LDR") for 10 populations of Aedes sierrensis, the western tree hole mosquito. The basic steps will involve:

1) Loading the data and libraries needed for analysis
2) Plotting the raw data to determine the appropriate functional form to fit
2) Setting up the JAGS model
3) Fitting using low information priors (uniform priors bounded by biologically realistic constraints)
4) Generating informative priors using a leave-one-out approach
5) Fitting using these informative priors
6) Checking the fits from the low information and informative priors
6) Plotting thermal response curves and parameter estimates 

The data can be found in this folder, titled "LifeHistoryTraitExp_TutorialData.csv"
