##### 
##### title: BRT for identifying hosts likely to be exposed to parasites in the leishmania genus - dont use imputed data
##### author: Caroline Glidden
##### version: 09/10/2021

########calculates feature importance
########calculates pdps
########uses bootstrap values to get host predictions

#------------------------------------------------------
# Set up
#------------------------------------------------------
rm(list=ls())

library(xgboost) 
library(rBayesianOptimization)
library(caret)
library(rsample) #to split stratified data
library(dplyr)
library(SHAPforxgboost)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #set working directory to current directory

#------------------------------------------------------
#read in final trait data -- that was cleaned and reduced in the variable selection script
#------------------------------------------------------

analysis_data <- readRDS("../cleaned data/analysis data/trait_data_small.rds")
#analysis_data <- readRDS("trait_data_small.rds")
mean_params <- read.csv("../output/mean parameter table leishmania july 8 2022.csv")
#mean_params <- read.csv("mean parameter table leishmania july 8 2022.csv")

#------------------------------------------------------#
#Feature importance - shapley                          #
#------------------------------------------------------#
####var importance --- use this to decide pdfs to run
means_df <- c()
dependence_df <- c()

for(i in 1:100){
  set.seed(i*1226)
  
  #shuffle before creating kfolds
  rows <- sample(nrow(analysis_data))
  analysis_data_v2 <- analysis_data[rows,]
  feature_data <- analysis_data_v2[,-2] #remove species names
  
  ##create test & train split
  train_index <- createDataPartition(feature_data$leish.infection, p = .7, list = FALSE)
  train <- feature_data[train_index,]; train_matrix <- as.matrix(train); d_train_matrix <- xgb.DMatrix(train_matrix[,-1], label= train_matrix[,1])
  test <- feature_data[-train_index,]; test_matrix <- as.matrix(test); d_test_matrix <- xgb.DMatrix(test_matrix[,-1], label= test_matrix[,1])
  
  #get AUC with all the features
  xgb.base <- xgboost(data = d_train_matrix,
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
  
  
  # To return the SHAP values and ranked features by mean|SHAP|
  shap_values <- shap.values(xgb_model = xgb.base, X_train = test_matrix[,-1])
  
  # The ranked features by mean |SHAP| for feature importance
  means <- as.data.frame(shap_values$mean_shap_score) #this gets global feature importance - this is the mean SHAP value per feature
  means$Feature <- rownames(means) #this just makes sure that the SHAP values is associated with each feature 
  means$iter <- i #saving iteration to potentially aggregate and plot better
  means_df <- rbind(means_df, means)
  
  #shap scores per observation for partial dependence (shap value X feature value)
  dependence_plots <- as.data.frame(cbind(analysis_data_v2[-train_index, ], shap_values$shap_score)) #this gets a dataframe of SHAP values per obs per feature (y axis), and binds it to the original dataframe so you retain the associated feature value (x axis)
  dependence_plots$iter <- i

  dependence_df <- rbind(dependence_df, dependence_plots)
  
}

names(means_df) <- c("mean_shap_score", "Feature", "iter")
write.csv(means_df, "../output/shapley variable importance leishmania july 11 2022.csv")
write.csv(dependence_df, "../output/shapley dependence leishmania july 11 2022.csv")

#######
#aggregate per group of features for each iteration
#######

#make a list of features to iterate through
features <- names(feature_data)[2:ncol(feature_data)]
feature_names <- features

pc <- c("mean_crop_cover", "mean_urban_cover","mean_tree_cover")
hab <- c("forest","savanna","shrubland","grassland","wetlands","rocky_areas","desert","artificial")
diet <- c("MammalEater","Insectivore","Frugivore","Granivore","Folivore")
activity <- c("ActivityCycle.1", "ActivityCycle.2", "ActivityCycle.3")
strata <- c("ForStrat.ValueA", "ForStrat.ValueAr", "ForStrat.ValueG", "ForStrat.ValueS")
#order <- c("MSW05_OrderPrimates","MSW05_OrderRodentia","MSW05_OrderCarnivora","MSW05_OrderChiroptera","MSW05_OrderDidelphimorphia","MSW05_OrderSoricomorpha","MSW05_OrderLagomorpha","MSW05_OrderCingulata","MSW05_OrderPilosa","MSW05_OrderArtiodactyla",
# #           "MSW05_OrderMicrobiotheria", "MSW05_OrderPaucituberculata", "MSW05_OrderPerissodactyla")
phylo <- c("V1","V2","V3","V4","V5","V6")
trophic <- c("TrophicLevel.yOmnivore","TrophicLevel.yCarnivore","TrophicLevel.yHerbivore")
drop_traits <- c(pc, hab, diet, activity, strata, phylo, trophic)

features <- features[!features %in% drop_traits]
list_features <- as.list(features)
list_features[[26]] <- pc
list_features[[27]] <- hab
list_features[[28]] <- diet
list_features[[29]] <- activity
list_features[[30]] <- strata
list_features[[31]] <- phylo
list_features[[32]] <- trophic
# #list_features[[33]] <- order

feature_names <- c(features, "pc", "hab", "diet", "activity", "strata", "phylo", "trophic") #took out order

sum_mean_shap_scores <- c()
for(j in 1:100) {
  
  sub_data <- subset(means_df, iter == j)
  
  new_df <- c()
  
  for(i in 1:length(list_features)) { #32 features
    
    new_dat <- sub_data[sub_data$Feature %in% list_features[[i]], ]
    shap_sum <-  sum(new_dat$mean_shap_score)
    
    new_df0 <- data.frame(Feature = feature_names[i], mean_shap_score = shap_sum, iter = j)
    new_df <- rbind(new_df, new_df0)
    
  }
  
  sum_mean_shap_scores <- rbind(sum_mean_shap_scores, new_df)
    
}


# #get top variables for further analysis
aggregate_mean <- aggregate(mean_shap_score ~ Feature, data = sum_mean_shap_scores, 
                              FUN = function(x) c(mean = mean(x), 
                                                  median = quantile(x, probs = 0.5),
                                                  lower.q = quantile(x, probs = 0.05),
                                                  upper.q = quantile(x, probs = 0.95)))

aggregate_mean <- do.call(data.frame, aggregate_mean)

write.csv(aggregate_mean, "mean sum shap scores leishmania july 18.csv")

imp_df <- subset(aggregate_mean, mean_shap_score.lower.q.5. > 0)

importance_plot_all <- ggplot(imp_df, aes(x = reorder(Feature, mean_shap_score.lower.q.5.), y = mean_shap_score.mean)) +
  geom_point(size = 3) + xlab('feature') + ylab('mean |shap value|') +
  geom_errorbar(aes(ymin = mean_shap_score.lower.q.5., ymax = mean_shap_score.upper.q.95.), position = "dodge", width = 0.4, size = 1.5) +
  #scale_color_manual(values=c('grey',colors[c(3,5,8)])) + 
  coord_flip() + theme_bw(base_size = 14) ##convert y to sd from the mean

 
# ##plot var importance
# #ggsave('plots/all wildlife importance plotapril 7.pdf',importance_plot_all, units='in', dpi=600)

#------------------------------------------------------#
#SHAPLEY hosts                                         #
#------------------------------------------------------#
mean_dependence_df <- data.frame(MSW05_Binomial = dependence_df$species, 
                                 iter = dependence_df$iter, 
                                 shap_prediction = as.numeric(rowSums(dependence_df[,60:116]))) #make sure tg

summary_shap_predictions <- aggregate(shap_prediction ~ MSW05_Binomial, data = mean_dependence_df, 
                                      FUN = function(x) c(mean = mean(x), 
                                                          median = quantile(x, probs = 0.5),
                                                          lower.q = quantile(x, probs = 0.05),
                                                          upper.q = quantile(x, probs = 0.95)))

summary_shap_predictions <- do.call(data.frame, summary_shap_predictions)
summary_shap_predictions <- left_join(summary_shap_predictions, 
                                      analysis_data[,c("leish.infection", "MSW05_Binomial")], by = "MSW05_Binomial")

write.csv(summary_shap_predictions, "../output/leishmania summary shap predictions july 19.csv")
