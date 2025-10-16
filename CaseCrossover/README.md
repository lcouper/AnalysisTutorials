## Case-Crossover 😷 

The R code and description below outline the steps for conducting a case-crossover analysis -- a commnon study design in epidemiology. If you'd like to run the code on your own, please download the accompanying mock datasets. 


**Overview**   
The case-crossover approach is useful for examining associations between exposures to some risk factor that varies over time (e.g., air pollution, alcohol consumption) and an acute outcome (e.g., heart attack, car crash).  
The study population consists entirely of individuals that have experienced the outcome, and inference is based on comparing that individuals' exposure to the risk factor before the outcome and during some control period.  
Because comparisons are made at the level of a given individual, confounders that don't change over time such as sex or genetics are inherently controlled for.  

In the code and explanation below, we will walk through an example of a case-crossover analysis to investigate the impact of oil and gas well development on coccidioidomycosis risk in Kern County, California, following the [manuscript available here](https://www.medrxiv.org/content/10.1101/2025.09.19.25336198v1). The tutorial uses the mock dataets 'mock_cases.csv' and 'mock_well_data.csv' uploaded here, as the real data are protected health information.  

This code and tutorial was made in collaboration with Dr. Jennifer Head and Philip Collender, MPH  


**Step 1. Load libraries**   
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
Here we will define the exposure window (here, the full number of days across which exposures to wells are considered relevant) and cutoff distance (here, the maximum radius around each case residence within which wells are considered relevant for exposure. In this analysis, we defined a 90 day exposure window to account for variability in the time between exposure and symptom onset, and time between symptom onset and diagnosis. We set the cutoff distance at 5km given prior work finding that air pollution from oil and gas well development was detectable up to this distance [Gonzalez et al. 2022](https://www.sciencedirect.com/science/article/pii/S0048969721053754) 
```
# exposure window
days <- 30*3

# cutoff distance (in meters)
limit <- 5000 # 5KM
```

**Step 4. Prepare the data**  
Here we will parse dates, convert well data to sf objects, and unionize buffers around all cases to trim wells to the relevant spatial extent. We will also restrict cases to the study period (2007-2022) and randomly assign a direction for pre/post control periods (i.e. applying a semi-symmetric bi-directional design).   
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
In order to count exposures by distance bands, here we calculate the well-to-case distances and build matrices encoding whether well activity fell inside each case’s hazard/control windows. 
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
Here, we will tabulate the number of wells with activity dates overlapping the hazard window. In our days, the hazard window is the 49-139 days prior to the case date (again to account for delays in symptom presentation, health-care seeking, diagnosis, and reporting). 
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

# Tip: case_di0_5 is the total wells within 0–5 km in the hazard window.

# Calculate how many cases had exposures during their hazard period
sum(mock_cases_sf$case_di0_5 >= 1) # 148
```

**Step 7. Calculate exposures during the control period**   
The logic here is to mirror the hazard window, but exactly one year before or one year after the case (bidirectional controls). This helps control for seasonality in exposures. Here, each case gets one matched control (but this could be increased). We then reshape the data to long format for conditional logistic regression (clogit, which estimates within-person odds ratios comparing hazard vs control exposures.   
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

# Convert wide to long. Note this uses a custom function written in the accompanying 'functions_plotting.R' script provided here 
cclong <- make_cclong(mock_cases_sf)
```

**Step 8. Generate negative control**   
While the case-crossover study design, and the approach to control period selection implemented here, controls for confounding from time-invariant factors, seasonality, and time trends in the exposure, we can further reduce residual confounding from any unmeasured, time-varying factors by including a negative control exposure in our model. Here, we have defined the negative control exposure as a “hazard” period that occurred in the 7-97 days *following* a case report. This serves as an ideal negative control as the construction of new wells occurring *after* the estimated disease onset is expected to be associated with the true exposure (as certain locations may experience more frequent oil and gas well construction), yet conditionally independent of disease onset in the absence of confounding, model misspecification, and measurement error.   
```
# Here, we exposure occurs in the 3.5 month period AFTER the event 
# For cases: keep distances if they fall in the right time period
distsCaseNC <- dists 
distsCaseNC[spud_dates_reg > (0-7)] = NA
distsCaseNC[spud_dates_reg < (-1*days)-7] = NA

# Tally wells within distances on each row
mock_cases_sfNC2 <- mock_cases_sf
mock_cases_sfNC2$case_di0_1 = colSums(distsCaseNC < 1000,na.rm = T)
mock_cases_sfNC2$case_di1_2 = colSums(distsCaseNC >= 1000 & distsCaseNC < 2000,na.rm = T)
mock_cases_sfNC2$case_di2_3 = colSums(distsCaseNC >= 2000 & distsCaseNC < 3000,na.rm = T)
mock_cases_sfNC2$case_di3_4 = colSums(distsCaseNC >= 3000 & distsCaseNC < 4000,na.rm = T)
mock_cases_sfNC2$case_di4_5 = colSums(distsCaseNC >= 4000 & distsCaseNC <= 5000,na.rm = T)

mock_cases_sfNC2 <- mock_cases_sfNC2 %>% mutate(case_di0_5 = case_di0_1 + case_di1_2 +
                                      case_di2_3 + case_di3_4 + case_di4_5)
# how many exposed during "hazard" period
sum(mock_cases_sfNC2$case_di0_5 >= 1) # 121

# Do it for the control period : keep distances if they fall in the right time period
distsControlpreNC <- distsControlpostNC <- dists 

distsControlpostNC[spud_dates_reg > (365 + 7 + days)] = NA 
distsControlpostNC[spud_dates_reg < 365 + 7] = NA

distsControlpreNC[spud_dates_reg > (0 - 365 + 7 + days)] = NA
distsControlpreNC[spud_dates_reg < 0 - 365 + 7] = NA

# Tally wells within distances on each row
mock_cases_sfNC2$contrpost_di0_1 = colSums(distsControlpostNC < 1000,na.rm = T)
mock_cases_sfNC2$contrpost_di1_2 = colSums(distsControlpostNC >= 1000 & distsControlpostNC < 2000,na.rm = T)
mock_cases_sfNC2$contrpost_di2_3 = colSums(distsControlpostNC >= 2000 & distsControlpostNC < 3000,na.rm = T)
mock_cases_sfNC2$contrpost_di3_4 = colSums(distsControlpostNC >= 3000 & distsControlpostNC < 4000,na.rm = T)
mock_cases_sfNC2$contrpost_di4_5 = colSums(distsControlpostNC >= 4000 & distsControlpostNC <= 5000,na.rm = T)

mock_cases_sfNC2$contrpre_di0_1 = colSums(distsControlpreNC < 1000,na.rm = T)
mock_cases_sfNC2$contrpre_di1_2 = colSums(distsControlpreNC >= 1000 & distsControlpreNC < 2000,na.rm = T)
mock_cases_sfNC2$contrpre_di2_3 = colSums(distsControlpreNC >= 2000 & distsControlpreNC < 3000,na.rm = T)
mock_cases_sfNC2$contrpre_di3_4 = colSums(distsControlpreNC >= 3000 & distsControlpreNC < 4000,na.rm = T)
mock_cases_sfNC2$contrpre_di4_5 = colSums(distsControlpreNC >= 4000 & distsControlpreNC <= 5000,na.rm = T)

mock_cases_sfNC2 <- mock_cases_sfNC2 %>% mutate(contr_di0_1 = ifelse(dir == 1, contrpost_di0_1,
                                                         contrpre_di0_1),
                                    contr_di1_2 = ifelse(dir == 1, contrpost_di1_2,
                                                         contrpre_di1_2),
                                    contr_di2_3 = ifelse(dir == 1, contrpost_di2_3,
                                                         contrpre_di2_3),
                                    contr_di3_4 = ifelse(dir == 1, contrpost_di3_4,
                                                         contrpre_di3_4),
                                    contr_di4_5 = ifelse(dir == 1, contrpost_di4_5,
                                                         contrpre_di4_5))

mock_cases_sfNC2 <- mock_cases_sfNC2 %>% mutate(contr_di0_5 = contr_di0_1 + contr_di1_2 +
                                      contr_di2_3 + contr_di3_4 + contr_di4_5)

# how many exposed during control, 'hazard' period
sum(mock_cases_sfNC2$contr_di0_5 >= 1) # 112

# run the conditional logistic regression model
cclongNC2 <- make_cclong(mock_cases_sfNC2)

# Combine observed and negative control data
cclongNC2sel <- cclongNC2 %>% select(names(cclongNC2)[c(10, 15, 18, 21, 24, 27:39)])
colnames(cclongNC2sel) <- paste0(colnames(cclongNC2sel), "NC")
cclongF <- cbind(cclong, cclongNC2sel)

coefsNC <- run_models_with_neg_control(cclongF, model5 = T)
colnames(coefsNC)[1:4] <- c("OR", "logOR", "lower", "upper")
```

**Step 9. Plot coefficient estimates**   
Belo, we plot odds ratios (OR) and confidence intervals for the association between well activity within distance rings and case onset, comparing hazard vs. control windows. OR > 1 suggests higher odds of coccidioidomycosis given exposure to well activity.
```
# First, plot without negative control
coefsNC = read.csv("Mock_CoefEstimates.csv", row.names = 1)

make_plots_obsvOnly(coefsNC, model_num = 1, xlab = "Distance from case", ystart = 0.2, h = 0.8)
make_plots_obsvOnly(coefsNC, model_num = 3, xlab = "Distance from case", ystart = 0.2, h = 0.8)
make_plots_obsvOnly(coefsNC, model_num = 5, xlab = "Quartile", ystart = 0.2, h = 0.8)

make_plots_obsvOnly(coefsNC, model_num = 2, xlab = "Distance from case", ystart = 0.2, h = 0.8)
make_plots_obsvOnly(coefsNC, model_num = 4, xlab = "Distance from case", ystart = 0.2, h = 0.8)

# Now, plot with negative controls
# main model
coefsNC = read.csv("Mock_CoefEstimates.csv", row.names = 1)

# continuous
make_plots(coefsNC, model_num = 1, xlab = "Distance from case", ystart = 0.2, h = 0.8)
# binary 
make_plots(coefsNC, model_num = 3, xlab = "Distance from case", ystart = 0.2, h = 0.8)
# quartile
make_plots(coefsNC, model_num = 5, xlab = "Quartile", ystart = 0.2, h = 0.8)
```

**Step 10. To stratify analyses by season of exposure**   
To investigate whether the association between well activity and coccidioidomycosis infection varies based on season of exposure, we can stratify our analyses by season. The logic here is that soil moisture, wind, disturbances, and other outdoor activity may vary by season, which may modify the association between wells and coccidioidomycosis risk. 
```
cclongF <- fread("MockExposures.csv")

##### SEASON #####

cclongF$Month = month(as.POSIXlt(cclongF$diagDate, format = "%Y-%m-%d"))
cclongF$Season = recode(cclongF$Month, '1' = "Winter", '2' = "Winter", '12' = "Winter",
                        '3' = "Spring", '4' = "Spring", '5' = "Spring",
                        '6' = "Summer", '7' = "Summer", '8' = "Summer",
                        '9' = "Fall", '10' = "Fall", '11' = "Fall")

# Defining season based on when case occurred, not on when well was spud
table(cclongF$Season)/2

coefsWinter <- run_models_with_neg_control(cclongF %>% subset(Season == "Winter"), model5 = F)
coefsSpring <- run_models_with_neg_control(cclongF %>% subset(Season == "Spring"), model5 = F) # not enough unique breaks to run models here
coefsSummer <- run_models_with_neg_control(cclongF %>% subset(Season == "Summer"), model5 = F)
coefsFall <- run_models_with_neg_control(cclongF %>% subset(Season == "Fall"), model5 = F)

coefsWinter$Estimate <- "Winter"
coefsSpring$Estimate <- "Spring"
coefsSummer$Estimate <- "Summer"
coefsFall$Estimate <- "Fall"

coefsSeason <- rbind(coefsWinter, coefsSpring, coefsSummer, coefsFall)
colnames(coefsSeason)[1:4] <- c("OR", "logOR", "lower", "upper")

# Remove estimates for negative controls for plotting:
coefsSeason = coefsSeason[!str_detect(rownames(coefsSeason), "NC"),]

# continuous, binary, quantile
make_plots_strata(coefsSeason, model_num = 1, xlab = "Distance from case", ystart = 0.2, h = 0.8)
make_plots_strata(coefsSeason, model_num = 3, xlab = "Distance from case", ystart = 0.2, h = 0.8)
# Note: not enough unique breaks within each season to run quartile model
make_plots_strata(coefsSeason, model_num = 5, xlab = "Quartile", ystart = 0.2, h = 0.8)

# Wald tests:
cclongWinter <- cclongF %>% subset(Season == "Winter")
cclongFall <- cclongF %>% subset(Season == "Fall")
cclongSpring <- cclongF %>% subset(Season == "Spring")
cclongSummer <- cclongF %>% subset(Season == "Summer")

waldtest(cclongFall, cclongSpring, model = "C") # sig diff at 0-5 , p : 0.006117326
waldtest(cclongWinter, cclongSummer, model = "B") # not sig
```


