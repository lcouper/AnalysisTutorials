## Case-Crossover 😷 

The code and description below outlines the steps for conducting a case-crossover -- a commnon study design in epidemiology.

**Overview**   
The case-crossover approach is useful for examining associations between exposures to some risk factor that varies over time (e.g., air pollution, alcohol consumption) and an acute outcome (e.g., heart attack, car crash).  
The study population consists entirely of individuals that have experienced the outcome, and inference is based on comparing that individuals exposure to the risk factor before the outcome and during some control period.  
Because comparisons are made at the level of a given individual, confounders that don't change over time such as sex or genetics are inherently controlled for.  

In the code and explanation below, we will walk through an example of a case-crossover analysis to investigate the impact of oil and gas well development on coccidioidomycosis risk in Kern County, California, following the [manuscript available here](https://www.medrxiv.org/content/10.1101/2025.09.19.25336198v1). The tutorial uses the mock dataets 'mock_cases.csv' and 'mock_well_data.csv' uploaded here, as the real data are protected health information.  


Step 1. 
