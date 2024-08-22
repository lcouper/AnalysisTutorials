# Gradient Boosted Machines 💻

The code and description below outlines the steps for conducting a gradient boosted machine approach using the [XGboost package](https://xgboost.readthedocs.io/en/stable/python/) in R.
Briefly: 
- Gradient boosting is a supervised machine learning approach in which regression or classification trees are sequentially built from the prediction errors of the prior tree. 
- GBM algorithms allow for complex, nonlinear relationships among predictor and outcome variables and collinearity between predictors, making them well suited for ecological analyses.
- Extreme gradient boosting is a scalable and efficient GBM implementation that minimizes overfitting and has been shown to achieve high predictive accuracy.

In this example, we develop an XGboost classification model to predict the presence or absence of various mosquito species based on ecological predictors (climate, land cover) and surveillance characteristics, while controlling for spatial and temporal confounders (lat/long, surveillance year)

This example accompanies the manuscript "Ecological drivers of dog heartworm transmission in California" available [here.](https://link.springer.com/article/10.1186/s13071-022-05526-x#Sec2)

**The basic steps are:**

1. Load packages and data
2. One-hot encode data
3. Remove highly correlated ecological predictors
4. Split data into testing (20%) and training (80%)
5. Define predictors & response variables
6. Specify xgboost model
7. Tune hyperparamters using (a) Bayesian optimization and (b) grid search
8. Run model with optimal hyperparameters
9. Evaluate model performance on test data
10. Examine feature importance
11. Bootstrap xgboost models

### 1. Load packages and data ###
```
library(xgboost)
library(mltools)
library(caret)
library(vip)
library(pdp)
library(rBayesianOptimization)
library(data.table)
library(pROC) 

setwd("~/Downloads/GBT")
Species = read.csv("Q2_AllSpeciesPresAbs_AllPredictors_1000.csv", header = T)[,-1]
```

### 2. One-hot encode the data ###

Categorical data can not be handled by many ML models and must be converted to numerical. One-hot encoding is a way of doing this conversion, in which each value of the categorical variable is separated into its own, binary variable (e.g., red, yellow, blue --> red (0/1), blue (0/1), yellow (0/1))

```
# first = convert the county name & trap type predictors (confounds) to factors
# then convert df to data table for one_hot function to work

Species$name = as.factor(Species$name)
Species$trap = as.factor(Species$trap)
SpeciesDT = data.table(Species)
Species2 = one_hot(SpeciesDT)
```

### 3. Remove highly correlated ecological predictors ###

ML models can handle correlated predictors, but including many collinear predictors can minimize interpretability.
Here, we remove any predictor with a > 0.90 pairwise correlation with another predictor (previously determined)
```
removeInds = c(81,82,84,85,87,88,93:97,99,100,102,103,122,123,124,128:139,152)
Species3 = subset(Species2, select = -c(removeInds))
```

### 4. Split data into training (80%) and testing set (20%) ###
This gives a way to evaluate model accuracy. Here, each observation was a surveillance record from a single collection date and location
```
set.seed(1234)
parts = sort(sample(nrow(Species3), nrow(Species3)*.80))
train = data.matrix(Species3[parts, ])
test = data.matrix(Species3[-parts, ])
```

### 5: Define predictors & response variables in testing & training data sets ###

Here, the predictors are:
- all climate & land cover variables, 
- trap features (type, # of nights, vector control agency), and
- potential spatial & temporal confounders (lat / long, surveillance year)
The response is presence in trap (0/1)

```
# doing first just Aedes sierrensis (column 160)
train_x = train[, -c(1,58,61,69:72,75:79,158:166)]
train_y = train[,160]
test_x = test[, -c(1,58,61,69:72, 75:79,158, 158:166)]
test_y = test[, 160]

# Create matrix for use in XGBoost 
xgb_train = xgb.DMatrix(data = train_x, label = train_y)
xgb_test = xgb.DMatrix(data = test_x, label = test_y)
```

### 6. Specify the model  #### 
A note on some of the key arguments:
- early_stopping_rounds : trains the model until train_rmse hasn't improved in x rounds
- max_dpeth: the maximum tree depth (or number of nodes)
- eta: controls learning rate. Smaller = less overfitting
  
```
minitial_xgb = xgboost(data = xgb_train, nrounds = 40, objective = "binary:logistic", eval_metric = "logloss",
                       early_stopping_rounds = 10
                       max.depth = 4
                       eta = 0.23) 
```
