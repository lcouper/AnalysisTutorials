# R Script Accompanying the "Bayesian Curve Fitting" analysis Tutorial

setwd("~/Documents/Current_Projects/LifeHistoryTraitExp/Analysis_TraitFits")

#### Load libraries and data ####
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

plot(LarvalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.15), pch = 16,
     data = data.LDR, ylab = "Larval dev rate", xlab = "Temperature")

#### Set up JAGS model ####
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

#### Set up MCMC simulation ####

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


#### Run JAGS model for each population ####

data <- data.LDR.HOP
trait <- data$LarvalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Then we run the JAGS model
LDR.HOP.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())
# Examine output
LDR.HOP.out$BUGSoutput$summary[1:5,]

plot(LarvalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.12), data = data.LDR.HOP, ylab = "LDR for Hop", xlab = "Temperature")
lines(LDR.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(LDR.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

# Repeat for other 9 populations 

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

#### Visually check fits from low informative priors ####

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

##### Generate informative priors using leave-one-out #####

data.LDR.HOP.prior <- subset(data.LDR, Population != "HOP")
data.LDR.MAR1.prior <- subset(data.LDR, Population != "MARIN35")
data.LDR.MAR2.prior <- subset(data.LDR, Population != "MARIN29")
data.LDR.WAW.prior <- subset(data.LDR, Population != "WAW")
data.LDR.EUG.prior <- subset(data.LDR, Population != "EUG")
data.LDR.PLA.prior <- subset(data.LDR, Population != "PLA")
data.LDR.SB.prior <- subset(data.LDR, Population != "SB")
data.LDR.JRA.prior <- subset(data.LDR, Population != "JRA")
data.LDR.PAR.prior <- subset(data.LDR, Population != "PAR")
data.LDR.POW.prior <- subset(data.LDR, Population != "POW")

# HOP
data <- data.LDR.HOP.prior
trait <- data$LarvalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)
LDR.HOP.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# MAR1
data <- data.LDR.MAR1.prior
trait <- data$LarvalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)
LDR.MAR1.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                           n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# MAR2
data <- data.LDR.MAR2.prior
trait <- data$LarvalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)
LDR.MAR2.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                           n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# WAW
data <- data.LDR.WAW.prior
trait <- data$LarvalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)
LDR.WAW.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# EUG
data <- data.LDR.EUG.prior
trait <- data$LarvalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)
LDR.EUG.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# PLA
data <- data.LDR.PLA.prior
trait <- data$LarvalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)
LDR.PLA.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# SB
data <- data.LDR.SB.prior
trait <- data$LarvalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)
LDR.SB.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# JRA
data <- data.LDR.JRA.prior
trait <- data$LarvalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)
LDR.JRA.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# PAR
data <- data.LDR.PAR.prior
trait <- data$LarvalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)
LDR.PAR.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# POW
data <- data.LDR.POW.prior
trait <- data$LarvalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)
LDR.POW.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Save output 

LDR.HOP.prior.cf.dists <- data.frame(q = as.vector(LDR.HOP.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(LDR.HOP.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(LDR.HOP.prior.out$BUGSoutput$sims.list$cf.Tm))
LDR.HOP.prior.gamma.fits = apply(LDR.HOP.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

LDR.MAR1.prior.cf.dists <- data.frame(q = as.vector(LDR.MAR1.prior.out$BUGSoutput$sims.list$cf.q),
                                      T0 = as.vector(LDR.MAR1.prior.out$BUGSoutput$sims.list$cf.T0), 
                                      Tm = as.vector(LDR.MAR1.prior.out$BUGSoutput$sims.list$cf.Tm))
LDR.MAR1.prior.gamma.fits = apply(LDR.MAR1.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

LDR.MAR2.prior.cf.dists <- data.frame(q = as.vector(LDR.MAR2.prior.out$BUGSoutput$sims.list$cf.q),
                                      T0 = as.vector(LDR.MAR2.prior.out$BUGSoutput$sims.list$cf.T0), 
                                      Tm = as.vector(LDR.MAR2.prior.out$BUGSoutput$sims.list$cf.Tm))
LDR.MAR2.prior.gamma.fits = apply(LDR.MAR2.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

LDR.WAW.prior.cf.dists <- data.frame(q = as.vector(LDR.WAW.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(LDR.WAW.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(LDR.WAW.prior.out$BUGSoutput$sims.list$cf.Tm))
LDR.WAW.prior.gamma.fits = apply(LDR.WAW.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

LDR.EUG.prior.cf.dists <- data.frame(q = as.vector(LDR.EUG.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(LDR.EUG.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(LDR.EUG.prior.out$BUGSoutput$sims.list$cf.Tm))
LDR.EUG.prior.gamma.fits = apply(LDR.EUG.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

LDR.PLA.prior.cf.dists <- data.frame(q = as.vector(LDR.PLA.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(LDR.PLA.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(LDR.PLA.prior.out$BUGSoutput$sims.list$cf.Tm))
LDR.PLA.prior.gamma.fits = apply(LDR.PLA.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

LDR.SB.prior.cf.dists <- data.frame(q = as.vector(LDR.SB.prior.out$BUGSoutput$sims.list$cf.q),
                                    T0 = as.vector(LDR.SB.prior.out$BUGSoutput$sims.list$cf.T0), 
                                    Tm = as.vector(LDR.SB.prior.out$BUGSoutput$sims.list$cf.Tm))
LDR.SB.prior.gamma.fits = apply(LDR.SB.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

LDR.JRA.prior.cf.dists <- data.frame(q = as.vector(LDR.JRA.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(LDR.JRA.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(LDR.JRA.prior.out$BUGSoutput$sims.list$cf.Tm))
LDR.JRA.prior.gamma.fits = apply(LDR.JRA.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

LDR.PAR.prior.cf.dists <- data.frame(q = as.vector(LDR.PAR.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(LDR.PAR.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(LDR.PAR.prior.out$BUGSoutput$sims.list$cf.Tm))
LDR.PAR.prior.gamma.fits = apply(LDR.PAR.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

LDR.POW.prior.cf.dists <- data.frame(q = as.vector(LDR.POW.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(LDR.POW.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(LDR.POW.prior.out$BUGSoutput$sims.list$cf.Tm))
LDR.POW.prior.gamma.fits = apply(LDR.POW.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

LDR.hypers <- list(LDR.HOP.prior.gamma.fits, LDR.MAR1.prior.gamma.fits, LDR.MAR2.prior.gamma.fits,
                   LDR.WAW.prior.gamma.fits, LDR.EUG.prior.gamma.fits, LDR.PLA.prior.gamma.fits,
                   LDR.SB.prior.gamma.fits, LDR.JRA.prior.gamma.fits, LDR.PAR.prior.gamma.fits, LDR.POW.prior.gamma.fits)
save(LDR.hypers, file = "~/Documents/Current_Projects/LifeHistoryTraitExp/Analysis_TraitFits/LDRhypers.Rsave")


#### Specify updated JAGS model #####
sink("briere_inf.txt")
cat("
    model{
    
    # 1. Specify the informative priors as gamma distributed
    cf.q ~ dgamma(hypers[1,1], hypers[2,1]) 
    cf.T0 ~ dgamma(hypers[1,2], hypers[2,2])
    cf.Tm ~ dgamma(hypers[1,3], hypers[2,3])
    cf.sigma ~ dunif(0, 1000)
    cf.tau <- 1 / (cf.sigma * cf.sigma)
    
    # 2. Specify the Likelihood function (here, Briere)
    for(i in 1:N.obs){
    trait.mu[i] <- cf.q * temp[i] * (temp[i] - cf.T0) * sqrt((cf.Tm - temp[i]) * (cf.Tm > temp[i])) * (cf.T0 < temp[i])
    trait[i] ~ dnorm(trait.mu[i], cf.tau)
    }
    
    # 3. Specify the Derived Quantities and Predictions
    for(i in 1:N.Temp.xs){
    z.trait.mu.pred[i] <- cf.q * Temp.xs[i] * (Temp.xs[i] - cf.T0) * sqrt((cf.Tm - Temp.xs[i]) * (cf.Tm > Temp.xs[i])) * (cf.T0 < Temp.xs[i])
    }
    
    } # close model
    ",fill=T)
sink()

##### Fit models for each population using informative priors #####
# First load saved hyperparameters
load("LDRhypers.Rsave")
LDR.HOP.prior.gamma.fits <- LDR.hypers[[1]]
LDR.MAR1.prior.gamma.fits <- LDR.hypers[[2]]
LDR.MAR2.prior.gamma.fits <- LDR.hypers[[3]]
LDR.WAW.prior.gamma.fits <- LDR.hypers[[4]]
LDR.EUG.prior.gamma.fits <- LDR.hypers[[5]]
LDR.PLA.prior.gamma.fits <- LDR.hypers[[6]]
LDR.SB.prior.gamma.fits <- LDR.hypers[[7]]
LDR.JRA.prior.gamma.fits <- LDR.hypers[[8]]
LDR.PAR.prior.gamma.fits <- LDR.hypers[[9]]
LDR.POW.prior.gamma.fits <- LDR.hypers[[10]]

# HOP
data <- data.LDR.HOP
hypers <- LDR.HOP.prior.gamma.fits * 0.1 # Note this is where you can adjust the weight of the priors. Decreasing this value (e.g., to '0.01' will increase the variance of the priors and reduce their effect on the posterior distribution)
trait <- data$LarvalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)
LDR.HOP.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# MAR1
data <- data.LDR.MAR1
hypers <- LDR.MAR1.prior.gamma.fits * 0.1
trait <- data$LarvalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)
LDR.MAR1.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# MAR2
data <- data.LDR.MAR2
hypers <- LDR.MAR2.prior.gamma.fits * 0.1
trait <- data$LarvalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)
LDR.MAR2.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# WAW
data <- data.LDR.WAW
hypers <- LDR.WAW.prior.gamma.fits * 0.1
trait <- data$LarvalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)
LDR.WAW.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# EUG
data <- data.LDR.EUG
hypers <- LDR.EUG.prior.gamma.fits * 0.1
trait <- data$LarvalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)
LDR.EUG.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# PLA
data <- data.LDR.PLA
hypers <- LDR.PLA.prior.gamma.fits * 0.1
trait <- data$LarvalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)
LDR.PLA.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# SB
data <- data.LDR.SB
hypers <- LDR.SB.prior.gamma.fits * 0.1
trait <- data$LarvalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)
LDR.SB.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                       n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# JRA
data <- data.LDR.JRA
hypers <- LDR.JRA.prior.gamma.fits * 0.1
trait <- data$LarvalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)
LDR.JRA.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# PAR
data <- data.LDR.PAR
hypers <- LDR.PAR.prior.gamma.fits * 0.1
trait <- data$LarvalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)
LDR.PAR.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# POW
data <- data.LDR.POW
hypers <- LDR.POW.prior.gamma.fits * 0.1
trait <- data$LarvalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)
LDR.POW.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())


#### Visually check model fits from informative priors #####

par(mfrow = c(5,2))
par(mar = c(2.5, 2.5, 1, 1))

plot(LarvalDevRate ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.12), data = data.LDR.HOP, ylab = "LDR for Hop", xlab = "Temperature")
lines(LDR.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(LDR.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col ="red")
lines(LDR.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("HOP", x = 36.5, y = 0.105, cex = 1.5)

plot(LarvalDevRate ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.12), data = data.LDR.MAR1, ylab = "LDR for MAR1", xlab = "Temperature")
lines(LDR.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(LDR.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(LDR.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("MAR1", x = 36.5, y = 0.105, cex = 1.5)

plot(LarvalDevRate ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.12), data = data.LDR.MAR2 , ylab = "LDR for MAR2 prior", xlab = "Temperature")
lines(LDR.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(LDR.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(LDR.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
text("MAR2", x = 36.5, y = 0.105, cex = 1.5)

plot(LarvalDevRate ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.12), data = data.LDR.WAW , ylab = "LDR for WAW prior", xlab = "Temperature")
lines(LDR.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(LDR.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(LDR.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
text("WAW", x = 36.5, y = 0.105, cex = 1.5)

plot(LarvalDevRate ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.12), data = data.LDR.EUG , ylab = "LDR for EUG prior", xlab = "Temperature")
lines(LDR.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(LDR.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(LDR.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
text("EUG", x = 36.5, y = 0.105, cex = 1.5)

plot(LarvalDevRate ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.12), data = data.LDR.PLA , ylab = "LDR for PLA prior", xlab = "Temperature")
lines(LDR.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(LDR.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(LDR.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
text("PLA", x = 36.5, y = 0.105, cex = 1.5)

plot(LarvalDevRate ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.12), data = data.LDR.SB , ylab = "LDR for SB prior", xlab = "Temperature")
lines(LDR.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(LDR.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(LDR.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
text("SB", x = 36.5, y = 0.105, cex = 1.5)

plot(LarvalDevRate ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.12), data = data.LDR.JRA , ylab = "LDR for JRA prior", xlab = "Temperature")
lines(LDR.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(LDR.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(LDR.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
text("JRA", x = 36.5, y = 0.105, cex = 1.5)

plot(LarvalDevRate ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.12), data = data.LDR.PAR , ylab = "LDR for PAR prior", xlab = "Temperature")
lines(LDR.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(LDR.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(LDR.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
text("PAR", x = 36.5, y = 0.105, cex = 1.5)

plot(LarvalDevRate ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.12), data = data.LDR.POW , ylab = "LDR for POW prior", xlab = "Temperature")
lines(LDR.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(LDR.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(LDR.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
text("POW", x = 36.5, y = 0.105, cex = 1.5)


#### Calculate additional thermal response parameters #####
Topt = function(x) {
  Matrix = x[["BUGSoutput"]][["sims.matrix"]]
  ToptVec = rep(NA, nrow(Matrix))
  for (i in 1:nrow(Matrix))
  {ToptVec[i] = Temp.xs[which.max(Matrix[i,6:86])]}
  return(c(mean(ToptVec), quantile(ToptVec, c(0.025, 0.975))))
}

# Function to calculate Pmax of LDR for each population ######

Pmax = function(x) {
  Matrix = x[["BUGSoutput"]][["sims.matrix"]]
  PmaxVec = rep(NA, nrow(Matrix))
  for (i in 1:nrow(Matrix))
  {PmaxVec[i] = max(Matrix[i,6:86])}
  return(c(mean(PmaxVec), quantile(PmaxVec, c(0.025, 0.975))))
}

# Function to calculate Tbreadth for each population #####

# Tbreadth defined as the temperature range where 
# performance remains >= 50% peak performance

Tbreadth = function(x) {
  Matrix = x[["BUGSoutput"]][["sims.matrix"]]
  TbreadthVec = rep(NA, nrow(Matrix))
  for (i in 1:nrow(Matrix))
  {fiftypeak = (max(Matrix[i,6:86])/2)
  fiftypeakindexes = which(Matrix[i,6:86] >= fiftypeak)
  mintempindex = fiftypeakindexes[1]
  maxtempindex = fiftypeakindexes[length(fiftypeakindexes)]
  TbreadthVec[i] = Temp.xs[maxtempindex] - Temp.xs[mintempindex]
  }
  return(c(mean(TbreadthVec), quantile(TbreadthVec, c(0.025, 0.975))))
}

#### Compile parameter estimates and save #####
# Function to compile mean & 95% credible intervals for all parameters for each population 
ParamCompile = function(x){
  DF = as.data.frame(matrix(,nrow = 5, ncol =4))
  colnames(DF) = c("mean", "lower95", "upper95", "param")
  DF[,4] = c("Tmin", "Tmax", "Topt", "Tbreadth", "Pmax")
  DF[1,1] = x$BUGSoutput$summary[1,1]
  DF[1,2:3] = hdi(x$BUGSoutput$sims.list$cf.T0, 0.95)[c(1,2)]
  DF[2,1] = x$BUGSoutput$summary[2,1]
  DF[2,2:3] = hdi(x$BUGSoutput$sims.list$cf.Tm, 0.95)[c(1,2)]
  DF[3,1:3] = as.numeric(Topt(x))
  DF[4,1:3] = as.numeric(Tbreadth(x))
  DF[5,1:3] = as.numeric(Pmax(x))
  return(DF)
}

# Store in a single dataframe
LDRdf_inf <- rbind.data.frame(ParamCompile(LDR.HOP.out.inf), ParamCompile(LDR.MAR1.out.inf),
                              ParamCompile(LDR.MAR2.out.inf), ParamCompile(LDR.WAW.out.inf),
                              ParamCompile(LDR.EUG.out.inf), ParamCompile(LDR.PLA.out.inf),
                              ParamCompile(LDR.SB.out.inf), ParamCompile(LDR.JRA.out.inf),
                              ParamCompile(LDR.PAR.out.inf), ParamCompile(LDR.POW.out.inf))
# Append the population names
LDRdf_inf$Population = c(rep(c("HOP", "MAR1", "MAR2", "WAW", "EUG", "PLA", "SB", "JRA", "PAR", "POW"), each = 5))
# Save for later use
save(LDRdf_inf, file = "LDR_meansd_inf.Rsave")

#### Plot thermal response curves and parameter estimates ####

library(ggplot2)

# The first plot will be the thermal response curves for each population:
load("LDR_meansd_inf.Rsave") # to bypass above code of generating parameter estimates & credible intervals

plot(LarvalDevRate ~ Temp.Treatment, 
     xlim = c(5, 45), ylim = c(0,0.11), data = data.LDR.HOP, type = "n", bty = "n",
     ylab = "rate (1/day)", xlab = "Temperature (\u00B0C)", pch = 1,
     main = "Larval development rates", cex.main = 2, cex.lab = 1.4, cex.axis = 1.2)
lines(LDR.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "#313695", lwd = 1.5)
lines(LDR.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "#4575b4", lwd = 1.5)
lines(LDR.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "#74add1", lwd = 1.5)
lines(LDR.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "#abd9e9", lwd = 1.5)
lines(LDR.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "#e0f3f8", lwd = 1.5)
lines(LDR.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "#fee090", lwd = 1.5)
lines(LDR.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "#fdae61", lwd = 1.5)
lines(LDR.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "#f46d43", lwd = 1.5)
lines(LDR.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "#ec3c30", lwd = 1.5)
lines(LDR.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "#ab041b", lwd = 1.5)

# The second plot will be show the T0, Topt, and Tm for each population with the points and error bars as the mean and 95% credible intervals, respectively
# First I'm ordering populations based on the latitude of their collection site
LatOrder = c("EUG", "HOP", "PLA", "MAR2", "MAR1", "JRA", "WAW", "PAR", "SB", "POW")
# I spent way too much time manually picking out these colors for each population, but use whatever color palette suits you! 
LatColors = rep(rev(c("#ab041b", "#ec3c30", "#f46d43", "#fdae61", "#fee090", "#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695")), each = 3)

LDRdf_inf$Population = factor(LDRdf_inf$Population, levels = LatOrder)
LDRdf_inf = LDRdf_inf[order(LDRdf_inf$Population),]
# I'm removing the Tbreadth and Pmax parameters for now
# But note that you could follow the below steps to create a similar plot for any individual or subset of the parameters
LDRdf_inf = LDRdf_inf[LDRdf_inf$param != "Tbreadth",] 
LDRdf_inf = LDRdf_inf[LDRdf_inf$param != "Pmax",]

ggplot(LDRdf_inf, aes(x=Population, y=mean, col = LatColors)) + 
  scale_y_continuous(limits = c(0, 45)) +
  geom_point(stat="identity", col=LatColors, 
             position=position_dodge(width = 1)) +
  geom_errorbar(aes(ymin=lower95, ymax=upper95), col = LatColors,
                position=position_dodge(.9), width = 0.3) + 
  ggtitle("Larval dev rates (1/day)") + 
  theme_minimal() +  coord_flip()  + 
  labs(x = "Population", y = "Temperature (C)") + 
  theme(axis.text=element_text(size=15), 
        axis.title = element_text(size = 16),
        legend.text=element_text(size=15),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 24),
        legend.position= "none") 