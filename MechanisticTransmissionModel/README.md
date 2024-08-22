# Mechanistic Transmission Model 

The code and description below outlines the steps for constructing a temperature-dependent, trait-based R0 model. In this example, we will specifically be estimating tempearure-based transmission suitability for dengue transmitted by *Aedes aegypti* 🦟

An overview of the steps below:
1. Load packages and climate data
2. Define functions for each temperature-dependent mosquito life history trait
3. Calculate R0 for the given month/year
4. Rescale and plot the data

Note this approach is outlined in [Mordecai et al. 2017](https://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0005568), and Dr. Erin Mordecai helped guide the construction of the code below

### 1. Load packages and climate data 

```
library(ncdf4)
library(raster)
library(scales)
library(tidyr)
library(sf)
library(maptools)
library(RNetCDF)
library(tigris)
library(RColorBrewer)
library(tigris)

setwd("TransmissionModel")
source("Temperature_R0_Functions.R")
```

### 2. Define functions for life history trait 
