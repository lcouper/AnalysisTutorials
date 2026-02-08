# Gradient Boosted Machines 💻

The code and description below outlines the steps for conducting a gradient boosted machine approach using the [XGboost package](https://xgboost.readthedocs.io/en/stable/python/) in R.

Briefly: 
- Gradient boosting is a supervised machine learning approach in which regression or classification trees are sequentially built from the prediction errors of the prior tree.
- GBM algorithms allow for complex, nonlinear relationships among predictor and outcome variables and can tolerate collinearity, making them well suited for ecological analyses.
- Extreme gradient boosting (XGBoost) is a scalable and efficient GBM implementation that minimizes overfitting and has been shown to achieve high predictive accuracy.
- In this example, we develop an XGBoost classification model to predict the presence or absence of a mosquito species based on ecological predictors (climate, land cover) and surveillance characteristics, while controlling for spatial confounders (e.g. latitude/longitude).

This example accompanies the manuscript "Ecological drivers of dog heartworm transmission in California" available [here.](https://link.springer.com/article/10.1186/s13071-022-05526-x#Sec2)

Note that help constructing this code was provided by Dr. Caroline Glidden

**Overview of analysis steps**
1. Load packages and data
2. Define outcome, predictors, and categorical variables
3. Remove highly correlated ecological predictors
4. One-hot encode categorical data
5. Split data into training (80%) and testing (20%) sets
6. Define predictor and response matrices
7. Tune hyperparameters using Bayesian optimization
8. Train the final XGBoost model
9. Evaluate model performance on test data
10. Examine feature importance
11. Bootstrap model performance and feature importance
12. Summarize bootstrap uncertainty
    
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

# Set working directory and seed for randomization 
setwd("~/Downloads")
set.seed(1234)

Species <- read.csv("Q2_AllSpeciesPresAbs_AllPredictors_1000.csv", header = TRUE)
Species <- Species[, -1]  # drop row index column if present
```

### 2. Define outcome, predictors, and categorical variables ###

Remove any columns in the dataset that are free-text fields (e.g., raw dates) or metadata fields that behave like identifiers (e.g., county), or generally non-meaningful for the prediction task.   
This is to avoid leakage and to keep interpretation focused on ecological predictors

```
response_var <- "AedesSierrensis"
id_vars <- c("name", "trap", "latitude", "longitude")
categorical_vars <- c("name", "trap")

drop_extra <- c("county", "collection_date", "trap_problem_bit", "sex", "collection_type", "closest_longitude", "closest_latitude", "trap_problem_bit", "collection_longitude", "collection_latitude")
Species <- Species %>% select(-any_of(drop_extra))
```

### 3. Remove highly correlated ecological predictors ###

Note: ML models can handle correlated predictors, but including many collinear predictors can minimize interpretability.
Here, we remove any predictor with a > 0.90 pairwise correlation with another predictor. 

```
# identify numerical predictors
candidate_vars <- setdiff(colnames(Species), c(response_var, id_vars, categorical_vars))
num_vars <- candidate_vars[sapply(Species[, candidate_vars, drop = FALSE], is.numeric)]

# remove any with no variation 
sds <- sapply(Species[, num_vars, drop = FALSE], sd, na.rm = TRUE)
num_vars_filt <- num_vars[is.finite(sds) & sds > 0]

# compute pairwise correlations
cor_mat <- cor(Species[, num_vars_filt, drop = FALSE], use = "pairwise.complete.obs")
cor_mat[is.na(cor_mat)] <- 0
diag(cor_mat) <- 1
high_corr <- findCorrelation(cor_mat, cutoff = 0.90, names = TRUE, exact = TRUE)

# How many highly correlated?
length(high_corr) # should be 36

# Remove these
Species <- Species[, !colnames(Species) %in% high_corr]
```


### 4. One-hot encode categorical data for XGBoost ###

Because XGBoost requires all numeric inputs, categorical predictors must be converted using one-hot encoding.

```
Species[categorical_vars] <- lapply(Species[categorical_vars], as.factor)
SpeciesDT <- as.data.table(Species)
SpeciesOH <- one_hot(SpeciesDT)
SpeciesOH <- as.data.frame(SpeciesOH)
```


### 5: Split data into training (80%) and testing set (20%) ###

```
n <- nrow(SpeciesOH)
train_idx <- sort(sample(seq_len(n), size = floor(0.8 * n)))

train_df <- SpeciesOH[train_idx, ]
test_df  <- SpeciesOH[-train_idx, ]
```

### 6. Create testing and training matrices ####

Here, the predictors are:
- all climate & land cover variables, 
- some trap features 
- potential spatial & temporal confounders (lat / long, surveillance year)
The response is presence in trap (0/1) (i.e. a binary variable)

```
train_y <- as.numeric(train_df[[response_var]])
test_y  <- as.numeric(test_df[[response_var]])

drop_vars <- c(response_var, id_vars)

pred_cols <- setdiff(colnames(train_df), drop_vars)
pred_cols <- pred_cols[sapply(train_df[, pred_cols, drop = FALSE], is.numeric)]

train_x <- as.matrix(train_df[, pred_cols, drop = FALSE])
test_x  <- as.matrix(test_df[, pred_cols, drop = FALSE])

xgb_train <- xgb.DMatrix(data = train_x, label = train_y)
xgb_test  <- xgb.DMatrix(data = test_x,  label = test_y)
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

