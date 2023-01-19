##### 
##### title: BRT for identifying wild hosts likely to be exposed to parasites in the leishmania genus - dont use imputed data
##### author: Caroline Glidden
##### version: 09/10/2021

##this script uses Bayesian optimization and nested-CV to select final parameter values 
##as well as evaluate model performance via calculating in/out sample AUC and target shuffling

#------------------------------------------------------
# Set up
#------------------------------------------------------

library(xgboost) 
library(rBayesianOptimization)
library(caret)
library(rsample) #to split stratified data
library(dplyr)

#------------------------------------------------------
#read in final trait data -- that was cleaned and reduced in the variable selection script
#------------------------------------------------------
#analysis_data <- readRDS('../cleaned data/analysis data/trait_data_small.rds')
analysis_data <- readRDS('trait_data_small.rds')
#table(trait_data_small$leish.infection) #102 positive, 1368 negative

#------------------------------------------------------#
#------------------------------------------------------#
#------------------------------------------------------#
# Nested cross-validation model -- test model performance
#------------------------------------------------------#
#------------------------------------------------------#
#------------------------------------------------------#

#initialize outputs
cv_preds_all <- c()
auc_list_all <- c()
sensitivity_list_all <- c()
parameter_table_all <- c()

 for(i in 1:25){
   set.seed(i*8120)
   
   ##remove species name
   loop_data <- analysis_data[,-2]
   
   #shuffle before creating kfolds
   rows <- sample(nrow(loop_data))
   loop_data <- loop_data[rows,]
   
   ##create three stratified folds for out nested cv loop
   train.index <- createFolds(loop_data[,1], k=3, list = FALSE)
   
   #inititalize outputs
   cv_preds <- as.data.frame(matrix(ncol=3, nrow=0)); names(cv_preds) <- c('MSW05_Binomial','xgbpred','true_vals')
   auc_list <- c()
   sensitivity_list <- c()
   parameter_table <- as.data.frame(matrix(nrow=9, ncol=3))
   names(parameter_table) <- c('outerfold_1','outerfold_2','outerfold_3')
   parameter_table$parameter <- c('test mean logloss','nrows','eta','max_depth','min_child',
                                  'subsample','colsample_bytree','gamma','ntrees')
   #nested cv
   for(j in 1:3){
     tryCatch({
     set.seed(i*j*982)
     
     test <- loop_data[which(train.index==j),]
     test <- as.matrix(test)
     d_test <- xgb.DMatrix(test[,-1], label=test[,1])
     
     train <- loop_data[which(train.index!=j),]
     train <- as.matrix(train)
     d_train <- xgb.DMatrix(train[,-1], label=train[,1])
     
     #------------------------------------------------------
     # New cross-validation function as the engine of Bayesian optimization
     #------------------------------------------------------
     ntrees.max = 200
     xgb_cv_bayes <- function(eta, max.depth, min.child.weight, subsample, colsample_bytree, gamma) {
       cv <- xgb.cv(params = list(booster = "gbtree",
                                  eta = eta,
                                  max_depth = max.depth,
                                  min_child_weight = min.child.weight,
                                  subsample = subsample,
                                  colsample_bytree = colsample_bytree,
                                  gamma = gamma,
                                  objective = "binary:logistic",
                                  eval_metric = "logloss",
                                  seed = 25),
                    data = d_train,
                    nrounds = ntrees.max,
                    nfold = 5, #this then uses 5 fold CV within this function
                    early_stopping_rounds = 10,
                    scale_pos_weight = 4, #sqrt(1368/102)
                    verbose = T)
       list(Score = -unlist(cv$evaluation_log[cv$best_iteration, "test_logloss_mean"]), # Ensure score is negative, since optimization maximizes
            Pred = cv$pred,
            cb.print.evaluation(period = 1))
     }
     
     #------------------------------------------------------
     # Acquire optimal parameters with Bayesian optimization (maximization function) via the R package "rBayesianOptimization"
     #------------------------------------------------------
     best_params <- BayesianOptimization(xgb_cv_bayes,
                                         bounds = list(eta = c(0.1, 0.3),
                                                       max.depth = c(2L, 10L),
                                                       min.child.weight = c(10L, 25L),
                                                       subsample = c(0.5, 0.8),
                                                       colsample_bytree = c(0.5, 0.8),
                                                       gamma = c(10, 20)),
                                         init_grid_dt = NULL,
                                         init_points = 10,
                                         n_iter = 40,
                                         acq = "ucb",
                                         kappa = 3,
                                         eps = 1.5,
                                         verbose = T)
     
     #------------------------------------------------------
     # Using the tuned hyperparameters, run a second cross-validation to acquire nrounds
     #------------------------------------------------------
     xgb_cv <- xgb.cv(params = best_params,
                      data = d_train,
                      nrounds = ntrees.max,
                      nfold = 5,
                      scale_pos_weight = 4,
                      early_stopping_rounds = 10,
                      objective = "binary:logistic",
                      eval_metric = "logloss",
                      verbose = T)
     
     best_params$nrounds <- xgb_cv$best_ntreelimit
     
     #------------------------------------------------------
     # Check evaluation log, to see that testing and training errors are declining -- need to make a data frame to save this
     # Ensure that optimized hyper parameters are within the pre specified bounds
     #------------------------------------------------------
     
     parameter_table[1,j] <- xgb_cv$evaluation_log$test_logloss_mean[1]
     parameter_table[2,1] <- xgb_cv$evaluation_log$test_logloss_mean[nrow(xgb_cv$evaluation_log)]
     
     parameter_table[3,j] <- xgb_cv$params$Best_Par[1] %>% round(4)
     parameter_table[4,j] <- xgb_cv$params$Best_Par[2] %>% round(4)
     parameter_table[5,j] <- xgb_cv$params$Best_Par[3] %>% round(4)
     parameter_table[6,j] <- xgb_cv$params$Best_Par[4] %>% round(4)
     parameter_table[7,j] <- xgb_cv$params$Best_Par[5] %>% round(4)
     parameter_table[8,j] <- xgb_cv$params$Best_Par[6] %>% round(4)
     parameter_table[9,j] <- best_params$nrounds
     
     #------------------------------------------------------
     # Run the full xgb model with the suite of optimal parameters
     #------------------------------------------------------
     watchlist <- list(train = d_train, test = d_test)
     xgb.fit <- xgboost(data = d_train,
                        eta = best_params$Best_Par[1],
                        max_depth = best_params$Best_Par[2],
                        min_child_weight = best_params$Best_Par[3],
                        subsample = best_params$Best_Par[4],
                        colsample_bytree = best_params$Best_Par[5],
                        gamma = best_params$Best_Par[6],
                        nrounds = best_params$nrounds,
                        scale_pos_weight = 4,
                        objective = "binary:logistic",
                        eval_metric = "logloss")
     
     ###prediction test
     xgbpred <- predict(xgb.fit, d_test); true_vals <- as.data.frame(test); true_vals <- true_vals$leish.infection
     xgbpred_df <- as.data.frame(cbind(loop_data[which(train.index==j),]$MSW05_Binomial, xgbpred, true_vals))
     cv_preds <- rbind(cv_preds, xgbpred_df)
     
     auc <- pROC::roc(response=true_vals, predictor=xgbpred, levels=c(0,1), auc = TRUE, plot = FALSE)
     auc_list <- c(auc_list, auc$auc)
     
     threshold <- pROC::coords(auc, "best", ret = "threshold")
     CM <- confusionMatrix(as.factor(ifelse(xgbpred >= threshold$threshold, 1, 0)), as.factor(true_vals), positive = "1")
     CM_table <- as.data.frame(CM[[4]])
     sensitivity <- CM_table[1,]
     sensitivity_list <- c(sensitivity_list, sensitivity)
     
     }, error=function(e){})
   }
   
   cv_preds_all <- rbind(cv_preds_all, cv_preds)
   auc_list_all <- c(auc_list_all, auc_list) 
   sensitivity_list_all <- c(sensitivity_list_all, sensitivity_list)
   parameter_table_all <- rbind(parameter_table_all, parameter_table)
   
 }
 
write.csv(cv_preds_all, "nested cv predictions leishmania july 8 2022.csv")
write.csv(auc_list_all, "nested cv auc leishmania july 8 2022.csv") 
write.csv(sensitivity_list_all, "nested cv sensitivity leishmania july 8 2022.csv") 
write.csv(parameter_table_all, "parameter table leishmania july 8 2022.csv")
 
#use mean of parameters for final model
parameter_table_rename <- parameter_table_all; names(parameter_table_rename) <- c('value','value','value','parameter')
reshape_params <- rbind(parameter_table_rename[,c(1,4)], parameter_table_rename[,c(2,4)], parameter_table_rename[,c(3,4)])
mean_params <-aggregate(reshape_params, by=list(reshape_params$parameter), 
                        FUN=mean, na.rm=TRUE)
write.csv(mean_params, "mean parameter table leishmania july 8 2022.csv")


#------------------------------------------------------#
#------------------------------------------------------#
#------------------------------------------------------#
# Nested cross-validation model w/target shuffling (shuffle response to see if model is fitting spurious correlations in the data -- I did this bc my data is super unbalanced)
#------------------------------------------------------#
#------------------------------------------------------#
#------------------------------------------------------#

#initialize outputs
auc_list_null <- c()

for(i in 1:25){
  set.seed(i*8120)
  
  ##remove species name
  loop_data <- analysis_data[,-2]
  null_data <- transform(loop_data, leish.infection = sample(leish.infection))
  
  #randomize rows
  rows <- sample(nrow(null_data))
  null_data <- null_data[rows,]
  
  ##create three stratified folds for out nested cv loop
  train.index <- createFolds(null_data[,1], k=3, list = FALSE)
  
  #inititalize outputs
  auc_list <- c()
  
  #nested cv
  for(j in 1:3){
    tryCatch({
    set.seed(i*j*982)
    
    train <- null_data[which(train.index!=j),]
    train <- as.matrix(train)
    d_train <- xgb.DMatrix(train[,-1], label=train[,1])
    
    test <- null_data[which(train.index==j),]
    true_vals <- test$leish.infection
    test <- as.matrix(test)
    d_test <- xgb.DMatrix(test[,-1], label=test[,1])
    
    #------------------------------------------------------
    # New cross-validation function as the engine of Bayesian optimization
    #------------------------------------------------------
    ntrees.max = 200
    xgb_cv_bayes <- function(eta, max.depth, min.child.weight, subsample, colsample_bytree, gamma) {
      cv <- xgb.cv(params = list(booster = "gbtree",
                                 eta = eta,
                                 max_depth = max.depth,
                                 min_child_weight = min.child.weight,
                                 subsample = subsample,
                                 colsample_bytree = colsample_bytree,
                                 gamma = gamma,
                                 objective = "binary:logistic",
                                 eval_metric = "logloss",
                                 seed = 25),
                   data = d_train,
                   nrounds = ntrees.max,
                   nfold = 5,
                   early_stopping_rounds = 10,
                   scale_pos_weight = 4,
                   verbose = T)
      list(Score = -unlist(cv$evaluation_log[cv$best_iteration, "test_logloss_mean"]), # Ensure score is negative, since optimization maximizes
           Pred = cv$pred,
           cb.print.evaluation(period = 1))
    }
    
    #------------------------------------------------------
    # Acquire optimal parameters with Bayesian optimization (maximization function) via the R package "rBayesianOptimization"
    #------------------------------------------------------
    best_params <- BayesianOptimization(xgb_cv_bayes,
                                        bounds = list(eta = c(0.1, 0.3),
                                                      max.depth = c(2L, 10L),
                                                      min.child.weight = c(10L, 25L),
                                                      subsample = c(0.5, 0.8),
                                                      colsample_bytree = c(0.5, 0.8),
                                                      gamma = c(10, 20)),
                                        init_grid_dt = NULL,
                                        init_points = 10,
                                        n_iter = 40,
                                        acq = "ucb",
                                        kappa = 3,
                                        eps = 1.5,
                                        verbose = T)
    
    #------------------------------------------------------
    # Using the tuned hyperparameters, run a second cross-validation to acquire nrounds
    #------------------------------------------------------
    xgb_cv <- xgb.cv(params = best_params,
                     data = d_train,
                     nrounds = ntrees.max,
                     nfold = 5,
                     scale_pos_weight = 4,
                     early_stopping_rounds = 10,
                     objective = "binary:logistic",
                     eval_metric = "logloss",
                     verbose = T)
    
    best_params$nrounds <- xgb_cv$best_ntreelimit
    
    #------------------------------------------------------
    # Run the full xgb model with the suite of optimal parameters
    #------------------------------------------------------
    watchlist <- list(train = d_train, test = d_test)
    xgb.fit <- xgboost(data = d_train,
                       eta = best_params$Best_Par[1],
                       max_depth = best_params$Best_Par[2],
                       min_child_weight = best_params$Best_Par[3],
                       subsample = best_params$Best_Par[4],
                       colsample_bytree = best_params$Best_Par[5],
                       gamma = best_params$Best_Par[6],
                       nrounds = best_params$nrounds,
                       scale_pos_weight = 4,
                       objective = "binary:logistic",
                       eval_metric = "logloss")
    
    ###prediction test
    xgbpred <- predict(xgb.fit, d_test)

    auc <- pROC::roc(response=true_vals, predictor=xgbpred, levels=c(0,1), auc = TRUE, plot = FALSE)
    auc_list <- c(auc_list, auc$auc)
    }, error=function(e){})
  }
  
  auc_list_null <- c(auc_list_null, auc_list) 
  
}

write.csv(auc_list_null, "nested cv target shuffled auc leishmania july 8 2022.csv")


#------------------
#get in sample AUC  (i.e., fit model to all of data) to check for over-fitting
#------------------
full_data <- analysis_data[,-2]
full_data <- as.matrix(full_data)
d_all <- xgb.DMatrix(full_data[,-1], label=full_data[,1])

full_auc <- c()
for(i in 1:100){ #iterate 100 times to get variation around model algorithm
  set.seed(i*9)
  
  xgb.full <- xgboost(data = d_all,
                      eta = mean_params[mean_params$Group.1=='eta','value'],
                      max_depth = round(mean_params[mean_params$Group.1=='max_depth','value']),
                      min_child_weight = mean_params[mean_params$Group.1=='min_child','value'],
                      subsample = mean_params[mean_params$Group.1=='subsample','value'],
                      colsample_bytree = mean_params[mean_params$Group.1=='colsample_bytree','value'],
                      gamma = mean_params[mean_params$Group.1=='gamma','value'],
                      nrounds = 15,
                      scale_pos_weight = 4,
                      objective = "binary:logistic",
                      eval_metric = "logloss")
  
  ###prediction test
  xgbpred <- predict(xgb.full, d_all); true_vals <- as.data.frame(analysis_data); true_vals <- true_vals$leish.infection
  auc <- pROC::roc(response=true_vals, predictor=xgbpred, levels=c(0,1), auc = TRUE, plot = TRUE)
  
  full_auc <- c(full_auc, auc$auc)  
  
}

write.csv(full_auc, "in sample auc leishmania july 8 2022.csv")

