## Case-Crossover 😷 

The code and description below outlines the steps for conducting a case-crossover -- a commnon study design in epidemiology.

**Overview**   
The case-crossover approach is useful for examining associations between exposures to some risk factor that varies over time (e.g., air pollution, alcohol consumption) and an acute outcome (e.g., heart attack, car crash).  
The study population consists entirely of individuals that have experienced the outcome, and inference is based on comparing that individuals exposure to the risk factor before the outcome and during some control period.  
Because comparisons are made at the level of a given individual, confounders that don't change over time such as sex or genetics are inherently controlled for.  

In the code and explanation below, we will walk through an example of a case-crossover analysis to investigate the impact of oil and gas well development on coccidioidomycosis risk in Kern County, California, following the [manuscript available here](https://www.medrxiv.org/content/10.1101/2025.09.19.25336198v1). The tutorial uses the mock dataets 'mock_cases.csv' and 'mock_well_data.csv' uploaded here, as the real data are protected health information.  

This code and tutorial was made in collaboration with Dr. Jennifer Head and Philip Collender, MPH  


**Step 1. Load libraries**
Here, we will first define the exposure window (fill in what this is), and cutoff distance (fill in what this is)
```
source("YourPathToMockDataFilesAndScript/functions_plotting.R")

library(sf); library(lubridate); library(plotrix); library(dplyr); library(ggplot2); library(tidyr)
library(stringr); library(survival); library(parallel); library(doParallel); library(foreach); library(splines); library(paletteer); library(cowplot); library(data.table); library(usethis)
```

**Step 2. Bring in (mock) datasets**  
The representative (mock) datasets provided here contain randomly generated values and rows but maintain the same data frame format and necessary columns as used in the analysis. 
```
mock_wells <- fread("mock_well_data.csv")
mock_cases <- fread("mock_cases.csv")
```

**Step 3. Set up variables**  
Here we will define the exposure window (fill in what this is) and cutoff distance (fill in what this is).    
```
# exposure window
days <- 30*3

# cutoff distance (in meters)
limit <- 5000 # 5KM
```

**Step 4. Data preparation**
(desription of how we are preparing the data for this analysis)

```
mock_wells <- mock_wells %>% mutate(
    INTPTLAT = as.numeric(gsub("^\\+", "", INTPTLAT)),
    INTPTLON = as.numeric(gsub("^\\+", "", INTPTLON)),
    GEOID = as.character(GEOID),  # keep as character to preserve leading zeros
    STATEFP = as.character(STATEFP),
    COUNTYFP = as.character(COUNTYFP),
    TRACTCE = as.character(TRACTCE))

# convert date columns
date_cols <- c("spud_date","completion_date","first_prod_date",
               "last_prod_date","date_earliest","date_latest")
mock_wells[ , (date_cols) := lapply(.SD, as.Date), .SDcols = date_cols]

# convert to sf object
mock_wells_sf <- st_as_sf(mock_wells, coords = c("INTPTLON","INTPTLAT"),
                    crs = 4326, remove = FALSE)

# establish buffers for the sake of cropping cases later (part 2)
mock_wells_sf <- st_transform(mock_wells_sf, crs = 3310)
mock_wells_buffers <- st_union(st_buffer(mock_wells_sf, dist = limit)) %>% st_as_sf()

mock_cases <- mock_cases %>% mutate(
    INTPTLAT20 = as.numeric(gsub("^\\+", "", INTPTLAT20)),
    INTPTLON20 = as.numeric(gsub("^\\+", "", INTPTLON20)))

# Make points in WGS84 (EPSG:4326), keep the INTPT* columns
mock_cases_sf <- st_as_sf(mock_cases, coords = c("INTPTLON20","INTPTLAT20"),
                          crs = 4326, remove = FALSE)

# Transform back to California Albers to match originals
mock_cases_sf <- st_transform(mock_cases_sf, 3310)

mock_cases_buffer <- st_union(st_buffer(mock_cases_sf, dist = limit))

# crop to the well buffers area
mock_cases_sf <- st_join(mock_cases_sf, mock_wells_buffers, left = F)

# correct the diagnosis dates
mock_cases_sf$diagDate <- as.Date(mock_cases_sf$Est_DtOnset, origin = '1960-01-01')
mock_cases_sf <- mock_cases_sf %>% subset(diagDate >= ymd('2007-01-01') & Est_DtOnset <= ymd('2022-12-31') ) # subset to 2007-2022
nobs <- nrow(mock_cases_sf) # 1770 cases in mock dataset

# Sample whether you will be sampling ahead of behind
mock_cases_sf$dir <- 2*rbinom(n = nobs, size = 1, prob = 0.5) - 1

# subset to wells spud within exposure period of any case
earliest_date = min(mock_cases_sf$diagDate) - days - (365 + 49)
latest_date = max(mock_cases_sf$diagDate) + 365 - 49

mock_wells_sf_sub = subset(mock_wells_sf, spud_date >= earliest_date & 
                     spud_date <= latest_date)

# subset to wells within 5 km of any case
mock_wells_sf_sub = mock_wells_sf_sub[unlist(st_contains(mock_cases_buffer, mock_wells_sf_sub)),] 
```

**Step 5. Calculate distances between all wells and cases**
(say something about why we're doing this)
```
dists <- units::drop_units(st_distance(mock_wells_sf_sub, mock_cases_sf,))
saveRDS(dists, file = "Calculated_distances.RDS")
dists <- readRDS("Calculated_distances.RDS") # each row: a well, each column: a case

# construct logical matrix for spudding within exposure period for each case (each row: a well, each column: a case)
spud_dates = matrix(mock_wells_sf_sub$spud_date, nrow = nrow(dists), ncol = ncol(dists)) 

# construct logical matrix for completion within exposure period for each case (each row: a well, each column: a case)
# add 7 to completion date to account for one week post-completion dust
mock_wells_sf_sub$completion_date <- mock_wells_sf_sub$completion_date + 7
comp_dates = matrix(mock_wells_sf_sub$completion_date, nrow = nrow(dists), ncol = ncol(dists)) 

# calculate difference between case and spud date, and case and completion date
diagDatenum <- as.numeric(mock_cases_sf$diagDate)
spud_dates_reg = t(diagDatenum - t(spud_dates)) # reg for regular
comp_dates_reg = t(diagDatenum - t(comp_dates)) # reg for regular

spud_dates_reg[1:4,1:4] # Here, each row represents a well, each column represents a case
comp_dates_reg[1:4,1:4] 
```

**Step 6. Calculate exposures during the hazard period**
```
# For cases: keep distances if they fall in the right time period
distsCase <- dists 
distsCase[spud_dates_reg > (days + 49) & comp_dates_reg > (days + 49)] = NA
distsCase[spud_dates_reg < (0 + 49) & comp_dates_reg < (0 + 49)] = NA 

# Tally wells within distances on each row
mock_cases_sf$case_di0_1 = colSums(distsCase < 1000,na.rm = T)
mock_cases_sf$case_di1_2 = colSums(distsCase >= 1000 & distsCase < 2000,na.rm = T)
mock_cases_sf$case_di2_3 = colSums(distsCase >= 2000 & distsCase < 3000,na.rm = T)
mock_cases_sf$case_di3_4 = colSums(distsCase >= 3000 & distsCase < 4000,na.rm = T)
mock_cases_sf$case_di4_5 = colSums(distsCase >= 4000 & distsCase <= 5000,na.rm = T)

mock_cases_sf <- mock_cases_sf %>% mutate(case_di0_5 = case_di0_1 + case_di1_2 +
                                case_di2_3 + case_di3_4 + case_di4_5)
summary(case_sf$case_di0_5)

# how many exposed during hazard period
sum(mock_cases_sf$case_di0_5 >= 1) # 148
```

**Step 7. Calculate exposures during the control period**

```
# For controls: keep distances if they fall in the right time period
distsControlpre <- distsControlpost <- dists

distsControlpost[spud_dates_reg > (days + 49) + 365 & comp_dates_reg > (days + 49) + 365] = NA
distsControlpost[spud_dates_reg < (0 + 49) + 365 & comp_dates_reg < (0 + 49) + 365] = NA

distsControlpre[spud_dates_reg > (days + 49) - 365 & comp_dates_reg > (days + 49) - 365] = NA
distsControlpre[spud_dates_reg < (0 + 49) - 365 & comp_dates_reg < (0 + 49) - 365] = NA

# Tally wells within distances on each row
mock_cases_sf$contrpost_di0_1 = colSums(distsControlpost < 1000,na.rm = T)
mock_cases_sf$contrpost_di1_2 = colSums(distsControlpost >= 1000 & distsControlpost < 2000,na.rm = T)
mock_cases_sf$contrpost_di2_3 = colSums(distsControlpost >= 2000 & distsControlpost < 3000,na.rm = T)
mock_cases_sf$contrpost_di3_4 = colSums(distsControlpost >= 3000 & distsControlpost < 4000,na.rm = T)
mock_cases_sf$contrpost_di4_5 = colSums(distsControlpost >= 4000 & distsControlpost <= 5000,na.rm = T)

mock_cases_sf$contrpre_di0_1 = colSums(distsControlpre < 1000,na.rm = T)
mock_cases_sf$contrpre_di1_2 = colSums(distsControlpre >= 1000 & distsControlpre < 2000,na.rm = T)
mock_cases_sf$contrpre_di2_3 = colSums(distsControlpre >= 2000 & distsControlpre < 3000,na.rm = T)
mock_cases_sf$contrpre_di3_4 = colSums(distsControlpre >= 3000 & distsControlpre < 4000,na.rm = T)
mock_cases_sf$contrpre_di4_5 = colSums(distsControlpre >= 4000 & distsControlpre <= 5000,na.rm = T)

mock_cases_sf <- mock_cases_sf %>% mutate(contr_di0_1 = ifelse(dir == 1, contrpost_di0_1,
                                                   contrpre_di0_1),
                              contr_di1_2 = ifelse(dir == 1, contrpost_di1_2,
                                                   contrpre_di1_2),
                              contr_di2_3 = ifelse(dir == 1, contrpost_di2_3,
                                                   contrpre_di2_3),
                              contr_di3_4 = ifelse(dir == 1, contrpost_di3_4,
                                                   contrpre_di3_4),
                              contr_di4_5 = ifelse(dir == 1, contrpost_di4_5,
                                                   contrpre_di4_5))

mock_cases_sf <- mock_cases_sf %>% mutate(contr_di0_5 = contr_di0_1 + contr_di1_2 +
                                contr_di2_3 + contr_di3_4 + contr_di4_5)

# how many exposed during control hazard period
sum(mock_cases_sf$contr_di0_5 >= 1) # 132

# Convert wide to long.
# Note this uses a custom function written in the accompanying 'functions_plotting.R' script provided here  
cclong <- make_cclong(mock_cases_sf)
```
