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

Note: ML models can handle correlated predictors, but including many collinear predictors can minimize interpretability.
In this dataset, we have already removed any predictor with a > 0.90 pairwise correlation with another predictor. 
But noting here that you may need to perform this step for your own data.


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
The response is presence in trap (0/1) (i.e. a binary variable)

```
# doing first just Aedes sierrensis (column 146)
train_x = train[, -c(144:152)]
train_y = train[,146]
test_x = test[, -c(144:152)]
test_y = test[, 146]

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

### 7. Tune hyperparameters  ###

There are 2 main options for hyperparmater tuning: grid search (uninformed, slower) and Bayesian optimization (faster, 'informed' based on prior results)

### 7a. Tune using Bayesian Optimization 

```
ntrees.max = 200
xgb_cv_bayes <- function(eta, max.depth,gamma) {
  cv <- xgb.cv(params = list(booster = "gbtree", eta = eta,
                             max_depth = max.depth, 
                             gamma = gamma,
                             objective = "binary:logistic",
                             eval_metric = "logloss",
                             seed = 25),
               data = xgb_train,
               nrounds = ntrees.max,
               nfold = 5, 
               early_stopping_rounds = 10, 
               scale_pos_weight = 4,
               verbose = T)
  list(Score = -unlist(cv$evaluation_log[cv$best_iteration, "test_logloss_mean"]), # Ensure score is negative, since optimization maximizes
       Pred = cv$pred,
       cb.print.evaluation(period = 1))}

best_params <- BayesianOptimization(xgb_cv_bayes,
                                    bounds = list(eta = c(0.1, 0.3),
                                                  max.depth = c(2L, 10L),
                                                  gamma = c(10, 20)),
                                    init_grid_dt = NULL, init_points = 10,
                                    n_iter = 3, acq = "ucb", kappa = 3,
                                    eps = 1.5, verbose = T)

```

### 7b. Tune using grid search
Note, this approach was not used in the manuscript, but followed the example [here.](https://www.r-bloggers.com/2020/11/r-xgboost-regression/)

```
# Set up grid with range of hyperparameters to try
hyper_grid <- expand.grid(max_depth = seq(3, 6, 1), eta = seq(.2, .35, .01))  

xgb_train_rmse = vector()
xgb_test_rmse = vector()

for (j in 1:nrow(hyper_grid)) {
  set.seed(1234) 
  m_xgb_untuned <- xgb.cv(
    data = train_x,
    label = train_y,
    nrounds = 10, # INCREASE LATER
    metrics = "rmse",
    early_stopping_rounds = 10,
    nfold = 5, # 5-fold cross validation
    max_depth = hyper_grid$max_depth[j],
    eta = hyper_grid$eta[j])
  # if you want to store all the RMSE somewhere (not necessary for rest of loop to work)
  xgb_train_rmse[j] <- m_xgb_untuned$evaluation_log$train_rmse_mean[m_xgb_untuned$best_iteration]
  xgb_test_rmse[j] <- m_xgb_untuned$evaluation_log$test_rmse_mean[m_xgb_untuned$best_iteration]
  cat(j, "\n")}    

# Find the best parameters from the grid seach above

df = cbind(xgb_train_rmse, xgb_test_rmse, 1:64) # 64 is based on grid size (4 depth values x 16 eta values)
which.min(df[,2]) 
hyper_grid[which.min(df[,2]),]
```

### 8. Run model with optimal hyper parameters ###

```
watchlist <- list(train = xgb_train, test = xgb_test)
mfit_xgb = xgboost(data = xgb_train, 
                   nrounds = 100, 
                   eta = best_params$Best_Par[1],
                   max_depth = best_params$Best_Par[2],
                   gamma = best_params$Best_Par[3],
                   objective = "binary:logistic",
                   eval_metric = "logloss",
                   scale_pos_weight = 4)
```

### 9. Evaluate model performance on test data ###

Here we use AUC for evaluation since the outcome is binary 

```
xgbpred = predict(mfit_xgb, xgb_test)
true_vals = as.data.frame(test)[,146] # This number will vary depending on species. 146 is for Aedes sierrensis
PredsTrue = cbind.data.frame(xgbpred, true_vals)
# calculate AUC and plot ROC curve
aucoutput = auc(PredsTrue$true_vals, PredsTrue$xgbpred) 
plot(roc(PredsTrue[,2], PredsTrue[,1]), main = "ROC curve")
```

### 10. Examine feature importance 

We can calculate a few metrics of "importance" for each predictor:
1) gain : how much the tree improves by adding a split based on a given variable
2) cover : how many observations related to this feature? 
3) frequency : relative % of times a features was used in trees 

```
importance_matrixFit <- xgb.importance(model = mfit_xgb)

# Plot "Importance" based on Gain
xgb.plot.importance(importance_matrixFit, xlab = "Feature Importance",
                    top_n = 12) # how many predictors to plot
```

<img width="631" alt="image" src="https://github.com/user-attachments/assets/2cfd58ea-57ed-48d3-9057-299fb9ce7828">


### 11. Bootstrap xgboost models ###

Now we can run the optimal model 100 times, each time on different 80% subset of data.  
This allows us to calculate confidence intervals for:  feature importance for each of the predictors and AUC stats

```
# 1. Importance matrix with rows for each predictor 
importanceMatrixSubSamps = data.frame(matrix(nrow = 143, ncol = 101))
importanceMatrixSubSamps[,1] = colnames(Species2)[-c(144:152)]
colnames(importanceMatrixSubSamps)[1] = "Feature"
# re-order alphabetically based on predictor name
importanceMatrixSubSamps = importanceMatrixSubSamps[order(importanceMatrixSubSamps[,'Feature']),]

# 2. AUC value from each run
AUC = vector(length = 100)

## Bootstrap run ##

for (i in 1:100){ 
  # select random 80% subsample of original data & set-up xgboost data
  # still splitting into testing and training so I can calculate AUC each time
  SubSampRows = sort(sample(nrow(Species2), nrow(Species2)*.80))
  SubSampTrain = data.matrix(Species2[SubSampRows, ])
  SubSampTest = data.matrix(Species2[-SubSampRows, ])
  SubSampTrain_x = SubSampTrain[, -c(144:152)]
  SubSampTrain_y = SubSampTrain[,146]
  SubSampTest_x = SubSampTest[, -c(144:152)]
  SubSampTest_y = SubSampTest[,146]
  
  xgb_SubSampTrain = xgb.DMatrix(data = SubSampTrain_x, label = SubSampTrain_y)
  xgb_SubSampTest = xgb.DMatrix(data = SubSampTest_x, label = SubSampTest_y)
  
  # run xgboost model
  mfit_xgb_SubSamp = xgboost(data = xgb_SubSampTrain, 
                             nrounds = 2,
                             eta = best_params$Best_Par[1],
                             max_depth = best_params$Best_Par[2],
                             gamma = best_params$Best_Par[3],
                             early_stopping_rounds = 10,
                             objective = "binary:logistic",
                             eval_metric = "logloss",
                             scale_pos_weight = 4)
  
  # Pull out feature importance
  temp = xgb.importance(model = mfit_xgb_SubSamp)
  # Reorder alphabetically using merge
  temp2 = merge(importanceMatrixSubSamps, temp[,1:2], by = "Feature", all.x = TRUE)
  # Add this re-ordered column (in column 102 of temp2, the merged df) to matrix
  importanceMatrixSubSamps[,(i+1)] = temp2[,102]
  
  # calculate AUC
  xgbpredSubSamp = predict(mfit_xgb_SubSamp, xgb_SubSampTest)
  true_valsSubSamp = as.data.frame(SubSampTest)[,146]
  PredsTrueSubSamp = cbind.data.frame(xgbpredSubSamp, true_valsSubSamp)
  # calculate AUC and plot ROC curve
  AUC[i] = auc(PredsTrueSubSamp[,2], PredsTrueSubSamp[,1]) 
}
```

Plot showing results of variable importance in predicting vector species presence/absence. For each species, predictors are ranked based on their mean gain across the 100 model iterations

![Image 1](https://github.com/user-attachments/assets/41a38c74-4200-4ead-8828-216594870d17)

