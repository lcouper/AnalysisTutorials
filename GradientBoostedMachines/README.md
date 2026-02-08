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


Necessary data can be downloaded [here.](https://github.com/lcouper/AnalysisTutorials/blob/main/GradientBoostedMachines/Q2_AllSpeciesPresAbs_AllPredictors_1000.csv)

To run this more easily on your own machine, the R script can be found [here.](https://github.com/lcouper/AnalysisTutorials/blob/main/GradientBoostedMachines/GradientBoostedMachine_RScript.R)
    
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

There are 2 main options for hyperparmater tuning: grid search (uninformed, slower) and Bayesian optimization (faster, 'informed' based on prior results). Here, we will do just Bayesian optimization.

```
# Maximum number of boosting iterations evaluated during cross-validation
ntrees_max <- 200 

# Function evaluated by the Bayesian optimizer
# Each call fits an XGBoost model using k-fold CV on the training data
# and returns the (negative) test log loss at the optimal iteration

xgb_cv_bayes <- function(eta, max_depth, gamma) {
  cv <- xgb.cv(params = list(
    booster = "gbtree",                    # tree-based boosting
    eta = eta,                             # learning rate (smaller = slower, less overfitting)
    max_depth = as.integer(max_depth),     # maximum depth of each tree
    gamma = gamma,                         # minimum loss reduction required for a split
    objective = "binary:logistic",          # binary classification (0/1 outcome)
    eval_metric = "logloss"),                 # evaluation metric minimized during CV
  data = xgb_train,                         
  nrounds = ntrees_max,                     # upper bound on number of trees
  nfold = 5,                                # 5-fold cross-validation
  early_stopping_rounds = 10,               # stop if log loss does not improve
  verbose = FALSE)
  
  list(Score = -cv$evaluation_log[cv$best_iteration, test_logloss_mean])}

# Run Bayesian optimization to efficiently explore hyperparameter space
opt <- BayesianOptimization(
  FUN = xgb_cv_bayes,
  bounds = list(
    eta = c(0.05, 0.3),                       # plausible learning rate range
    max_depth = c(2L, 8L),                    # shallow to moderately deep trees
    gamma = c(1, 10)),                        # from weak to strong split penalty
  init_points = 10,                           # random initial evaluations
  n_iter = 5,                                 # number of optimization steps
  acq = "ucb",                                # upper confidence bound acquisition
  kappa = 2.576,                              # exploration vs exploitation tradeoff
  verbose = TRUE)

# Best-performing hyperparameter combination
best_params <- opt$Best_Par
```



### 8. Run model with optimal hyper parameters ###

```
params <- list(booster = "gbtree",
          eta = unname(best_params["eta"]),
          max_depth = as.integer(best_params["max_depth"]),
          gamma = unname(best_params["gamma"]),
          objective = "binary:logistic",  
          eval_metric = "logloss")

mfit_xgb <- xgb.train(params = params,
          data = xgb_train,
          nrounds = ntrees_max,
          watchlist = list(train = xgb_train, test = xgb_test),
          early_stopping_rounds = 10,
          verbose = 1)
```

### 9. Evaluate model performance on test data ###

Here we use AUC for evaluation since the outcome is binary 

```
xgb_pred <- predict(mfit_xgb, xgb_test)
roc_obj <- roc(test_y, xgb_pred)
auc_val <- auc(roc_obj)
plot(roc_obj, main = sprintf("ROC curve (AUC = %.3f)", auc_val))

# Confusion matrix; sensitivity & specificity
pred_class <- ifelse(xgb_pred >= 0.5, 1, 0)
confusionMatrix(factor(pred_class, levels = c(0, 1)),
  factor(test_y, levels = c(0, 1)))
```

### 10. Examine feature importance 

We can calculate a few metrics of "importance" for each predictor:
1) gain : how much the tree improves by adding a split based on a given variable
2) cover : how many observations related to this feature? 
3) frequency : relative % of times a features was used in trees 

```
imp <- xgb.importance(model = mfit_xgb)

xgb.plot.importance(imp, top_n = 15, xlab = "Feature importance (Gain)")
```

### 11. Bootstrap xgboost models ###

Now we can run the optimal model 100 times, each time on different 80% subset of data.  
This allows us to calculate confidence intervals for:  feature importance for each of the predictors and AUC stats

```
# Repeat model fitting on 100 random 80% subsamples to quantify uncertainty
# Note this can take a while to run

B <- 100
auc_boot <- numeric(B)

imp_boot <- data.frame(Feature = imp$Feature, matrix(NA, nrow = nrow(imp), ncol = B))

for (b in seq_len(B)) {
  idx <- sort(sample(seq_len(n), size = floor(0.8 * n)))
  
  train_b <- SpeciesOH[idx, ]
  test_b  <- SpeciesOH[-idx, ]
  
  train_x_b <- as.matrix(train_b[, !colnames(train_b) %in% drop_vars])
  train_y_b <- as.numeric(train_b[[response_var]])
  
  test_x_b  <- as.matrix(test_b[, !colnames(test_b) %in% drop_vars])
  test_y_b  <- as.numeric(test_b[[response_var]])
  
  xgb_train_b <- xgb.DMatrix(train_x_b, label = train_y_b)
  xgb_test_b  <- xgb.DMatrix(test_x_b,  label = test_y_b)
  
  m_b <- xgboost(
    data = xgb_train_b,
    nrounds = mfit_xgb$best_iteration,
    eta = best_params["eta"],
    max_depth = as.integer(best_params["max_depth"]),
    gamma = best_params["gamma"],
    objective = "binary:logistic",
    eval_metric = "logloss",
    verbose = 0)
  
  pred_b <- predict(m_b, xgb_test_b)
  auc_boot[b] <- auc(test_y_b, pred_b)
  
  imp_b <- xgb.importance(model = m_b)
  imp_boot[, b + 1] <- imp_b$Gain[match(imp_boot$Feature, imp_b$Feature)]}

```

### 12. Create bootstrap summaries for model performance and feature importance ####

```
# AUC mean + 95% CI
auc_summary <- c(
  mean = mean(auc_boot, na.rm = TRUE),
  lwr  = quantile(auc_boot, 0.025, na.rm = TRUE),
  upr  = quantile(auc_boot, 0.975, na.rm = TRUE))

auc_summary

# show distribution of AUC values 
hist(auc_boot, breaks = 20, main = "Bootstrap AUC distribution", xlab = "AUC")
abline(v = auc_summary, col = c("black","red","red"), lwd = 2, lty = c(1,2,2))

# Feature importance mean + 95% CI
imp_mat <- imp_boot[, -1, drop = FALSE]
imp_boot$selection_rate <- rowMeans(!is.na(imp_mat))
imp_boot$mean_gain <- rowMeans(imp_mat, na.rm = TRUE)
imp_boot$lwr_gain <- apply(imp_mat, 1, quantile, probs = 0.025, na.rm = TRUE)
imp_boot$upr_gain <- apply(imp_mat, 1, quantile, probs = 0.975, na.rm = TRUE)

# Keep features with finite importance and reasonable stability
imp_boot <- imp_boot[is.finite(imp_boot$mean_gain) & imp_boot$mean_gain > 0 & imp_boot$selection_rate >= 0.2, ]
imp_boot <- imp_boot[order(-imp_boot$mean_gain), ]

# Plot top features with uncertainty (most important at top)
top_imp <- head(imp_boot, 15)
ypos <- rev(seq_len(nrow(top_imp)))
xlims <- range(c(top_imp$lwr_gain, top_imp$upr_gain), finite = TRUE, na.rm = TRUE)

par(mar = c(5, 10, 2, 2))
plot(top_imp$mean_gain, ypos, xlim = xlims, yaxt = "n",
     xlab = "Feature importance (Gain)", ylab = "", pch = 16)
segments(top_imp$lwr_gain, ypos, top_imp$upr_gain, ypos)
axis(2, at = ypos, labels = top_imp$Feature, las = 2)
```
