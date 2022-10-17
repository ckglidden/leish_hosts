##### 
##### title: BRT for testing trait profile & performance of study effort model
##### author: Caroline Glidden
##### version: 03/16/2022

##this script uses Bayesian optimization and nested-CV to test model performance of study effort model 
##as well as determine trait profile

#------------------------------------------------------
# Set up
#------------------------------------------------------
#rm(list=ls())

library(xgboost) 
library(rBayesianOptimization)
library(caret)
library(rsample) #to split stratified data
library(dplyr)

#setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #set working directory to current directory

#------------------------------------------------------
# Final steps of data clean up
#------------------------------------------------------

#analysis_data <- readRDS("../cleaned data/analysis data/pubmed_analysis_data.rds")
analysis_data <- readRDS("pubmed_analysis_data.rds")

#------------------------------------------------------#
#------------------------------------------------------#
#------------------------------------------------------#
# Nested cross-validation model -- test model performance
#------------------------------------------------------#
#------------------------------------------------------#
#------------------------------------------------------#
rsq <- function (x, y) cor(x, y) ^ 2

#initialize outputs
#cv_preds_all <- c()
r2_list_all <- c()
parameter_table_all <- c()

for(i in 1:2){
  tryCatch({
    set.seed(i*8120)
    
    ##remove species name
    loop_data <- analysis_data[,-2]
    
    ##create three stratified folds for out nested cv loop
    train.index <- createFolds(loop_data[,1], k=3, list = FALSE)
    
    #inititalize outputs
    #cv_preds <- as.data.frame(matrix(ncol=3, nrow=0)); names(cv_preds) <- c('MSW05_Binomial','xgbpred','true_vals')
    r2_list <- c()
    parameter_table <- as.data.frame(matrix(nrow=9, ncol=3))
    names(parameter_table) <- c('outerfold_1','outerfold_2','outerfold_3')
    parameter_table$parameter <- c('test mae','nrows','eta','max_depth','min_child',
                                   'subsample','colsample_bytree','gamma','ntrees')
    #nested cv
    for(j in 1:3){
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
                                   objective = "reg:squarederror",
                                   eval_metric = "mae",
                                   seed = 25),
                     data = d_train,
                     nrounds = ntrees.max,
                     nfold = 5,
                     early_stopping_rounds = 10,
                     #scale_pos_weight = 4,
                     verbose = T)
        list(Score = -unlist(cv$evaluation_log[cv$best_iteration, "test_mae_mean"]), # Ensure score is negative, since optimization maximizes
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
                       #scale_pos_weight = 4,
                       early_stopping_rounds = 10,
                       objective = "reg:squarederror",
                       eval_metric = "mae",
                       verbose = T)
      
      best_params$nrounds <- xgb_cv$best_ntreelimit
      
      #------------------------------------------------------
      # Check evaluation log, to see that testing and training errors are declining -- need to make a data frame to save this
      # Ensure that optimized hyper parameters are within the pre specified bounds
      #------------------------------------------------------
      
      parameter_table[1,j] <- xgb_cv$evaluation_log$test_mae_mean[1]
      parameter_table[2,1] <- xgb_cv$evaluation_log$test_mae_mean[nrow(xgb_cv$evaluation_log)]
      
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
                         #scale_pos_weight = 4,
                         objective = "reg:squarederror",
                         eval_metric = "mae")
      
      ###prediction test
      xgbpred <- predict(xgb.fit, d_test); true_vals <- as.data.frame(test); true_vals <- true_vals$pubmed.count
      
      r2 <- rsq(xgbpred, true_vals)
      r2_list <- c(r2_list, r2)
    }
    
    r2_list_all <- c(r2_list_all, r2_list) 
    parameter_table_all <- rbind(parameter_table_all, parameter_table)
  }, error=function(e){})
}

write.csv(r2_list_all, "study effort nested cv r2 july 6 2022.csv")
write.csv(parameter_table_all, "study effort parameter table july2 2022.csv")

#use mean of parameters for final model
parameter_table_rename <- parameter_table_all; names(parameter_table_rename) <- c('value','value','value','parameter')
reshape_params <- rbind(parameter_table_rename[,c(1,4)], parameter_table_rename[,c(2,4)], parameter_table_rename[,c(3,4)])
mean_params <-aggregate(reshape_params, by=list(reshape_params$parameter), 
                        FUN=mean, na.rm=TRUE)
write.csv(mean_params, "study mean parameter table july 6 2022.csv")
#mean_params <- read.csv("study mean parameter table march 16 2022.csv")
# 

