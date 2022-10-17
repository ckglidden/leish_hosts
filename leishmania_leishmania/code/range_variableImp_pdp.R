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

new_hosts <- subset(summary_shap_predictions, shap_prediction.lower.q.5. > 0) #161 observations
new_hosts <- subset(new_hosts, leish.infection == 0) #98 unrecognized hosts

get_orders <- readRDS('/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_leishmania/cleaned data/analysis data/total data.rds')
orders <- unique(get_orders[,c("MSW05_Binomial", "MSW05_Order")])

new_hosts <- left_join(new_hosts, orders, by = "MSW05_Binomial")
table(new_hosts$MSW05_Order)



































#------------------------------------------------------#
#PDP data & host predictions                           #
#------------------------------------------------------#
##create list for top variables named imp_df
cv_preds <- c()
# pd_df_all <- c() #store partial dependence data
# 
# names_imp_df <- feature_names[1:2]
# 
# for(i in 1:100){
#   set.seed(i*8120)
#   
#   pred_data <- analysis_data
#   
#   #shuffle before creating kfolds
#   rows <- sample(nrow(pred_data))
#   pred_data <- pred_data[rows,]
#   
#   ###remove species name
#   loop_data <- pred_data[,-2]
#   
#   ##create test & train split
#   train_index <- createDataPartition(loop_data$leish.infection, p = .8, list = FALSE)
#   train <- loop_data[train_index,]; train_matrix <- as.matrix(train); d_train_matrix <- xgb.DMatrix(train_matrix[,-1], label= train_matrix[,1])
# 
#   xgb.fit <- xgboost(data = d_train_matrix,
#                        eta = mean_params[mean_params$Group.1=='eta','value'],
#                        max_depth = round(mean_params[mean_params$Group.1=='max_depth','value']),
#                        min_child_weight = mean_params[mean_params$Group.1=='min_child','value'],
#                        subsample = mean_params[mean_params$Group.1=='subsample','value'],
#                        colsample_bytree = mean_params[mean_params$Group.1=='colsample_bytree','value'],
#                        gamma = mean_params[mean_params$Group.1=='gamma','value'],
#                        nrounds = 15,
#                        scale_pos_weight = 4,
#                        objective = "binary:logistic",
#                        eval_metric = "logloss")
#     
#     ###prediction test
#     xgbpred <- predict(xgb.fit, d_train_matrix); true_vals <- train$leish.infection
#     xgbpred_df <- as.data.frame(cbind(pred_data[train_index,]$MSW05_Binomial, xgbpred, true_vals))
#     names(xgbpred_df) <- c("MSW05_Binomial", "predicted_probability", "leish.infection")
#     cv_preds <- rbind(cv_preds, xgbpred_df)
#     
#     #pdps
#     pd_df = data.frame(matrix(vector(), 0, 3, dimnames=list(c(), c('variable', 'value','yhat'))),
#                        row.names = NULL, stringsAsFactors=F)
#     
#     for (k in 1:length(names_imp_df)) { #loop through each variable
#       
#       output <- as.data.frame(pdp::partial(xgb.fit, pred.var = names_imp_df[k], prob = TRUE, train = bts_data[,-1]))
#       
#       loop_pdp_df <- data.frame(matrix(vector(), nrow(output), 4,
#                                    dimnames=list(c(), c('variable', 'value','yhat','iter'))), stringsAsFactors=F,
#                             row.names=NULL)
#       
#       loop_pdp_df$variable <- names_imp_df[k]
#       loop_pdp_df$value <- output[[1]]
#       loop_pdp_df$yhat <- output[[2]]
#       loop_pdp_df$iter <- i
#       
#       pd_df <- rbind(pd_df, loop_pdp_df)
#       
#     }
#     
#     pd_df_all <- rbind(pd_df_all, pd_df)
#     
#   }
#   
# 
# #names(cv_predictions_5x) <- c('MSW05_Binomial', 'predicted_probability')
# #cv_predictions_5x <- merge(cv_predictions_5x, analysis_data[,c('leish.infection','MSW05_Binomial')])
# write.csv(cv_predictions_5x, "../output/leishmania july 7 bootstrap predictions.csv")
# write.csv(pd_df_iter, "../output/leishmania july 7 bootstrap pdp prob to use.csv")
# 
# 
# ##plot pdp
# rug_data <- analysis_data[,names_imp_df[c(2:8, 15:17)]]
# rug_data <- rug_data %>%
#   tidyr::pivot_longer(cols= 1:10,names_to = "variable", values_to = "value")
# rug_data <- as.data.frame(rug_data)
# 
# pdp_all_wildlife <- ggplot(pd_df_iter, aes(x=value, y=yhat)) + 
#   stat_smooth(aes(group=iter), color='lightgrey', method='loess', size=0.5, se=FALSE) + 
#   stat_smooth(aes(), method='loess', size=2, se=FALSE) +
#   geom_rug(data=rug_data, aes(x=value), alpha=0.3, length=unit(0.05, "npc"), inherit.aes = FALSE, sides='b') +
#   facet_wrap(~variable, scales='free', ncol=3) +
#   #scale_color_manual(values=c('darkgrey',colors[c(3,5,8)])) +
#   theme_bw(base_size = 12)
# 
# #ggsave('plots/all wildlife pdp plot april 7.pdf', pdp_all_wildlife, dpi=600)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ###############################
# ###############extra code    #
# ###############################
# 
# #####get figure for phylogenetic traits
# phylo_traits <- c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10")
# pd_df_iter_phylo <- c() #store partial dependence data
# 
# for(i in 1:100){
#   set.seed(i*8120)
#   
#   ###remove species name
#   loop_data <- analysis_data[,-2]
#   
#   ##create 5 stratified folds
#   bootstrap_folds <- createFolds(loop_data[,1], k=5, list = FALSE)
#   
#   #inititalize outputs
#   cv_preds <- as.data.frame(matrix(ncol=3, nrow=0)); names(cv_preds) <- c('MSW05_Binomial','xgbpred','true_vals')
#   
#   #initialize df to save pdp?
#   pd_df_all <- c()
#   
#   #bootstrap (80% of data)
#   for(j in 1:5){
#     set.seed(i*j*982)
#     
#     bts_data <- loop_data[which(bootstrap_folds!=j),]
#     bts_matrix <- as.matrix(bts_data)
#     d_bts_data <- xgb.DMatrix(bts_matrix[,-1], label=bts_matrix[,1])
#     
#     xgb.fit <- xgboost(data = d_bts_data,
#                        eta = mean_params[mean_params$Group.1=='eta','value'],
#                        max_depth = round(mean_params[mean_params$Group.1=='max_depth','value']),
#                        min_child_weight = mean_params[mean_params$Group.1=='min_child','value'],
#                        subsample = mean_params[mean_params$Group.1=='subsample','value'],
#                        colsample_bytree = mean_params[mean_params$Group.1=='colsample_bytree','value'],
#                        gamma = mean_params[mean_params$Group.1=='gamma','value'],
#                        nrounds = 15,
#                        scale_pos_weight = 4,
#                        objective = "binary:logistic",
#                        eval_metric = "logloss")
#     
#     
#     #pdps
#     pd_df = data.frame(matrix(vector(), 0, 3, dimnames=list(c(), c('variable', 'value','yhat'))),
#                        row.names = NULL, stringsAsFactors=F)
#     
#     for (k in 1:10) { #loop through each variable
#       
#       output <- as.data.frame(pdp::partial(xgb.fit, pred.var = phylo_traits[k], center = TRUE, train = bts_data[,-1]))
#       
#       loop_pdp_df <- data.frame(matrix(vector(), nrow(output), 4,
#                                        dimnames=list(c(), c('variable', 'value','yhat','iter'))), stringsAsFactors=F,
#                                 row.names=NULL)
#       
#       loop_pdp_df$variable <- phylo_traits[k]
#       loop_pdp_df$value <- output[[1]]
#       loop_pdp_df$yhat <- output[[2]]
#       loop_pdp_df$iter <- i
#       
#       pd_df <- rbind(pd_df, loop_pdp_df)
#       
#     }
#     
#     pd_df_all <- rbind(pd_df_all, pd_df)
#     
#   }
#   
#   pd_df_iter_phylo <- rbind(pd_df_iter_phylo, pd_df_all)
# }
# 
# pdp_phylo <- ggplot(pd_df_iter_phylo, aes(x=value, y=yhat)) + 
#   stat_smooth(aes(group=iter), color='lightgrey', method='loess', size=0.5, se=FALSE) + 
#   stat_smooth(aes(), method='loess', size=2, se=FALSE) +
#   #geom_rug(data=rug_data, aes(x=value), alpha=0.3, length=unit(0.05, "npc"), inherit.aes = FALSE, sides='b') +
#   facet_wrap(~variable, scales='free', ncol=3) +
#   #scale_color_manual(values=c('darkgrey',colors[c(3,5,8)])) +
#   theme_bw(base_size = 12)
# ggsave('../output/plots/pdp phylo plot april 29.png', pdp_phylo, dpi=600)
# 
# 
# #------------------------------------------------------#
# #Feature importance - gain                             #
# #------------------------------------------------------#
# feature_data <- analysis_data[,-2] #remove species names
# 
# ####var importance --- use this to decide pdfs to run
# importance_all <- c()
# 
# for(i in 1:20){
#   set.seed(i*8120)
#   
#   #shuffle before creating kfolds
#   rows <- sample(nrow(feature_data))
#   feature_data <- feature_data[rows,]
#   
#   ##create 5 stratified folds
#   bootstrap_folds <- createFolds(feature_data[,1], k=5, list = FALSE)
#   
#   #initialize outputs
#   importance_df <- c()
#   
#   #bootstrap (80% of data)
#   for(j in 1:5){
#     set.seed(i*j*982)
#     
#     bts_data <- feature_data[which(bootstrap_folds!=j),]
#     bts_data <- as.matrix(bts_data)
#     d_bts_data <- xgb.DMatrix(bts_data[,-1], label=bts_data[,1])
#     
#     xgb.fit <- xgboost(data = d_bts_data,
#                        eta = mean_params[mean_params$Group.1=='eta','value'],
#                        max_depth = round(mean_params[mean_params$Group.1=='max_depth','value']),
#                        min_child_weight = mean_params[mean_params$Group.1=='min_child','value'],
#                        subsample = mean_params[mean_params$Group.1=='subsample','value'],
#                        colsample_bytree = mean_params[mean_params$Group.1=='colsample_bytree','value'],
#                        gamma = mean_params[mean_params$Group.1=='gamma','value'],
#                        nrounds = 15,
#                        scale_pos_weight = 4,
#                        objective = "binary:logistic",
#                        eval_metric = "logloss")
#     
#     
#     ###var importance
#     #importance <- EIX::importance(xgb.fit, d_bts_data, option = "both")
#     #importance <- as.data.frame(matrix(unlist(importance), ncol = 7, nrow = nrow(importance)))
#     importance <- xgb.importance(model=xgb.fit, feature_names = names(bts_data[,-1]))
#     
#     a <- colnames(bts_data[,-1]); b <- importance$Feature; add <- a[!a %in% b]
#     add_df <- data.frame(Feature = add, Gain=rep(0, length(add)), Cover=rep(0, length(add)), Frequency=rep(0, length(add)))
#     importance <- rbind(importance, add_df)
#     
#     pc <- c("mean_crop_cover", "mean_urban_cover","mean_tree_cover"); pc_df <- importance %>% filter(importance$Feature %in% pc)
#     add_lulc <- data.frame(Feature = "LULC", Gain=sum(pc_df$Gain), Cover=sum(pc_df$cover), Frequency=sum(pc_df$Frequency))
#     
#     hab <- c("forest","savanna","shrubland","grassland","wetlands","rocky_areas","desert","artificial"); hab_df <- importance %>% filter(importance$Feature %in% hab)
#     add_hab <- data.frame(Feature = "habitat", Gain=sum(hab_df$Gain), Cover=sum(hab_df$cover), Frequency=sum(hab_df$Frequency))
#     
#     diet <- c("MammalEater","Insectivore","Frugivore","Granivore","Folivore"); diet_df <- importance %>% filter(importance$Feature %in% diet)
#     add_diet <- data.frame(Feature = "diet", Gain=sum(diet_df$Gain), Cover=sum(diet_df$cover), Frequency=sum(diet_df$Frequency))
#     
#     activity <- c("ActivityCycle.1", "ActivityCycle.2", "ActivityCycle.3"); activity_df <- importance %>% filter(importance$Feature %in% activity)
#     add_activity <- data.frame(Feature = "activity", Gain=sum(activity_df$Gain), Cover=sum(activity_df$cover), Frequency=sum(activity_df$Frequency))
#     
#     strata <- c("ForStrat.ValueA", "ForStrat.ValueAr", "ForStrat.ValueG", "ForStrat.ValueS"); strata_df <- importance %>% filter(importance$Feature %in% strata)
#     add_strata <- data.frame(Feature = "strata", Gain=sum(strata_df$Gain), Cover=sum(strata_df$cover), Frequency=sum(strata_df$Frequency))
#     
#     order <- c("MSW05_OrderPrimates","MSW05_OrderRodentia","MSW05_OrderCarnivora","MSW05_OrderChiroptera","MSW05_OrderDidelphimorphia","MSW05_OrderSoricomorpha","MSW05_OrderLagomorpha","MSW05_OrderCingulata","MSW05_OrderPilosa","MSW05_OrderArtiodactyla",
#                "MSW05_OrderMicrobiotheria", "MSW05_OrderPaucituberculata", "MSW05_OrderPerissodactyla"); order_df <- importance %>% filter(importance$Feature %in% order)
#     add_order <- data.frame(Feature = "order", Gain=sum(order_df$Gain), Cover=sum(order_df$cover), Frequency=sum(order_df$Frequency))
#     
#     phylo <- c("V1","V2","V3","V4","V5","V6"); phylo_df <- importance %>% filter(importance$Feature %in% phylo)
#     add_phylo <- data.frame(Feature = "phylo", Gain=sum(phylo_df$Gain), Cover=sum(phylo_df$cover), Frequency=sum(phylo_df$Frequency))
#     
#     trophic <- c("TrophicLevel.yOmnivore","TrophicLevel.yCarnivore","TrophicLevel.yHerbivore"); trophic_df <- importance %>% filter(importance$Feature %in% trophic)
#     add_trophic <- data.frame(Feature = "trophic", Gain=sum(trophic_df$Gain), Cover=sum(trophic_df$cover), Frequency=sum(trophic_df$Frequency))
#     
#     importance <- rbind(importance, add_lulc, add_hab, add_diet, add_activity, add_strata, add_order, add_phylo, add_trophic)
#     
#     importance_df <- rbind(importance_df, importance)
#     
#   }
#   
#   importance_df$iter <- i
#   importance_all <- rbind(importance_all, importance_df)
#   
# }
# 
# write.csv(importance_all, "../output/bootstrap variable importance leishmania july 7 2022.csv")
# 
# 
# #get top variables for pdp plots
# mean_importance <- aggregate(Gain ~ Feature, data = importance_all, 
#                              FUN = function(x) c(mean = Rmisc::CI(x, 0.95)[2], 
#                                                  lowerCI = Rmisc::CI(x, 0.95)[3], 
#                                                  upperCI = Rmisc::CI(x, 0.95)[1],
#                                                  median = quantile(x, probs = 0.5),
#                                                  lower.q = quantile(x, probs = 0.025),
#                                                  upper.q = quantile(x, probs = 0.975)))
# 
# mean_importance <- do.call(data.frame,mean_importance)
# mean_importance <- mean_importance %>% arrange(desc(Gain.median.50.))
# 
# write.csv(mean_importance, "../output/bootstrap mean variable importance leishmania july 7 2022.csv")
# 
# # ##get mean of categorical variables: pc_land, habitat, diet, activity, strata, order, phylogenetic
# # 
# # pc <- c("mean_crop_cover", "mean_urban_cover","mean_tree_cover"); pc_df <- mean_importance %>% filter(mean_importance$Feature %in% pc); LULC <- sum(pc_df$Gain.mean.mean)
# # pc_df$sq.ci <- I(pc_df$Gain.upperCI.upper - pc_df$Gain.lowerCI.lower)^2; sqrt_ci <- sqrt(sum(pc_df$sq.ci)); LULC.lowerCI <- LULC - sqrt_ci;  LULC.upperCI <- LULC + sqrt_ci
# # 
# # hab <- c("forest","savanna","shrubland","grassland","wetlands","rocky_areas","desert","artificial"); hab_df <- mean_importance %>% filter(mean_importance$Feature %in% hab); HAB <- sum(hab_df$Gain.mean.mean)
# # diet <- c("MammalEater","Insectivore","Frugivore","Granivore","Folivore"); diet_df <- mean_importance %>% filter(mean_importance$Feature %in% diet); DIET <- sum(diet_df$Gain.mean.mean)
# # activity <- c("ActivityCycle.1", "ActivityCycle.2", "ActivityCycle.3"); activity_df <- mean_importance %>% filter(mean_importance$Feature %in% activity); ACT <- sum(activity_df$Gain.mean.mean)
# # strata <- c("ForStrat.ValueA", "ForStrat.ValueAr", "ForStrat.ValueG", "ForStrat.ValueS"); strata_df <- mean_importance %>% filter(mean_importance$Feature %in% strata); SRT <- sum(strata_df$Gain.mean.mean)
# # order <- c("MSW05_OrderPrimates","MSW05_OrderRodentia","MSW05_OrderCarnivora","MSW05_OrderChiroptera","MSW05_OrderDidelphimorphia","MSW05_OrderSoricomorpha","MSW05_OrderLagomorpha","MSW05_OrderCingulata","MSW05_OrderPilosa","MSW05_OrderArtiodactyla",
# #            "MSW05_OrderMicrobiotheria", "MSW05_OrderPaucituberculata", "MSW05_OrderPerissodactyla"); order_df <- mean_importance %>% filter(mean_importance$Feature %in% order); ORD <- sum(order_df$Gain.mean.mean)
# # phylo <- c("V1","V2","V3","V4","V5","V6"); phylo_df <- mean_importance %>% filter(mean_importance$Feature %in% phylo); PHY <- sum(phylo_df$Gain.mean.mean)
# # phylo_df$sq.ci <- I(phylo_df$Gain.upperCI.upper - phylo_df$Gain.lowerCI.lower)^2; sqrt_ci <- sqrt(sum(phylo_df$sq.ci)); PHY.lowerCI <- PHY - sqrt_ci;  PHY.upperCI <- PHY + sqrt_ci
# # 
# # trophic <- c("TrophicLevel.yOmnivore","TrophicLevel.yCarnivore","TrophicLevel.yHerbivore"); troph_df <- mean_importance %>% filter(mean_importance$Feature %in% trophic); TRP <- sum(troph_df$Gain.mean.mean)
# # 
# # sum_gains <- c(LULC, PHY, HAB, DIET, ACT, SRT, ORD, TRP)
# # sum_gains_names <- c("range_lulc","phylo_distance","habitat","diet","activity_cycle","strata","order","trophic_level")
# # categorical <- data.frame(Feature = sum_gains_names, Gain.mean.mean = sum_gains, 
# #                           Gain.lowerCI.lower = c(LULC.lowerCI, PHY.lowerCI, rep(NA, 6)), Gain.upperCI.upper = c(LULC.upperCI, PHY.upperCI, rep(NA, 6)), Gain.n = rep(500, length(sum_gains)), Gain.se = rep(NA, length(sum_gains)))
# # 
# # mean_importance <- rbind(mean_importance, categorical)
# 
# #drop single traits
# pc <- c("mean_crop_cover", "mean_urban_cover","mean_tree_cover")
# hab <- c("forest","savanna","shrubland","grassland","wetlands","rocky_areas","desert","artificial")
# diet <- c("MammalEater","Insectivore","Frugivore","Granivore","Folivore")
# activity <- c("ActivityCycle.1", "ActivityCycle.2", "ActivityCycle.3")
# strata <- c("ForStrat.ValueA", "ForStrat.ValueAr", "ForStrat.ValueG", "ForStrat.ValueS")
# order <- c("MSW05_OrderPrimates","MSW05_OrderRodentia","MSW05_OrderCarnivora","MSW05_OrderChiroptera","MSW05_OrderDidelphimorphia","MSW05_OrderSoricomorpha","MSW05_OrderLagomorpha","MSW05_OrderCingulata","MSW05_OrderPilosa","MSW05_OrderArtiodactyla",
#            "MSW05_OrderMicrobiotheria", "MSW05_OrderPaucituberculata", "MSW05_OrderPerissodactyla")
# phylo <- c("V1","V2","V3","V4","V5","V6")
# trophic <- c("TrophicLevel.yOmnivore","TrophicLevel.yCarnivore","TrophicLevel.yHerbivore")
# drop_traits <- c(pc, hab, diet, activity, strata, order, phylo, trophic)
# 
# mean_importance <- mean_importance[!mean_importance$Feature %in% drop_traits, ]
# 
# write.csv(mean_importance, "../output/bootstrap sum mean variable importance leishmania july 7 2022.csv")
# 
# imp_df <- subset(mean_importance, Gain.lower.q.2.5. > 0.005)
# importance_plot_all <- ggplot(imp_df, aes(x = reorder(Feature, Gain.lower.q.2.5.), y = Gain.median.50.)) +
#   geom_point(size = 3) + xlab('feature') + ylab('importance') +
#   geom_errorbar(aes(ymin = Gain.lower.q.2.5., ymax = Gain.upper.q.97.5.), position = "dodge", width = 0.4, size = 1.5) +
#   #scale_color_manual(values=c('grey',colors[c(3,5,8)])) + 
#   coord_flip() +
#   theme_bw(base_size = 14) ##convert y to sd from the mean
# 
# ##plot var importance
# #ggsave('plots/all wildlife importance plotapril 7.pdf',importance_plot_all, units='in', dpi=600)
# 
# # pc <- c("mean_crop_cover", "mean_urban_cover","mean_tree_cover"); pc_df <- importance_df %>% filter(importance_df$Feature %in% pc)
# # add_lulc <- data.frame(Feature = "LULC", delta_auc=sum(pc_df$delta_auc), delta_sensitivity=sum(pc_df$delta_sensitivity),
# #                        perm_auc = NA, base_auc = NA)
# # 
# # hab <- c("forest","savanna","shrubland","grassland","wetlands","rocky_areas","desert","artificial"); hab_df <- importance_df %>% filter(importance_df$Feature %in% hab)
# # add_hab <- data.frame(Feature = "habitat", delta_auc=sum(hab_df$delta_auc), delta_sensitivity=sum(hab_df$delta_sensitivity),
# #                       perm_auc = NA, base_auc = NA)
# # 
# # diet <- c("MammalEater","Insectivore","Frugivore","Granivore","Folivore"); diet_df <- importance_df %>% filter(importance_df$Feature %in% diet)
# # add_diet <- data.frame(Feature = "diet", delta_auc=sum(diet_df$delta_auc), delta_sensitivity=sum(diet_df$delta_sensitivity),
# #                        perm_auc = NA, base_auc = NA)
# # 
# # activity <- c("ActivityCycle.1", "ActivityCycle.2", "ActivityCycle.3"); activity_df <- importance_df %>% filter(importance_df$Feature %in% activity)
# # add_activity <- data.frame(Feature = "activity", delta_auc=sum(activity_df$delta_auc), delta_sensitivity=sum(activity_df$delta_sensitivity),
# #                            perm_auc = NA, base_auc = NA)
# # 
# # strata <- c("ForStrat.ValueA", "ForStrat.ValueAr", "ForStrat.ValueG", "ForStrat.ValueS"); strata_df <- importance_df %>% filter(importance_df$Feature %in% strata)
# # add_strata <- data.frame(Feature = "strata", delta_auc=sum(strata_df$delta_auc), delta_sensitivity=sum(strata_df$delta_sensitivity),
# #                          perm_auc = NA, base_auc = NA)
# # 
# # order <- c("MSW05_OrderPrimates","MSW05_OrderRodentia","MSW05_OrderCarnivora","MSW05_OrderChiroptera","MSW05_OrderDidelphimorphia","MSW05_OrderSoricomorpha","MSW05_OrderLagomorpha","MSW05_OrderCingulata","MSW05_OrderPilosa","MSW05_OrderArtiodactyla",
# #            "MSW05_OrderMicrobiotheria", "MSW05_OrderPaucituberculata", "MSW05_OrderPerissodactyla"); order_df <- importance_df %>% filter(importance_df$Feature %in% order)
# # add_order <- data.frame(Feature = "order", delta_auc=sum(order_df$delta_auc), delta_sensitivity=sum(order_df$delta_sensitivity),
# #                         perm_auc = NA, base_auc = NA)
# # 
# # phylo <- c("V1","V2","V3","V4","V5","V6"); phylo_df <- importance_df %>% filter(importance_df$Feature %in% phylo)
# # add_phylo <- data.frame(Feature = "phylo", delta_auc=sum(phylo_df$delta_sensitivity), delta_sensitivity=sum(phylo_df$delta_sensitivity),
# #                         perm_auc = NA, base_auc = NA)
# # 
# # trophic <- c("TrophicLevel.yOmnivore","TrophicLevel.yCarnivore","TrophicLevel.yHerbivore"); trophic_df <- importance_df %>% filter(importance_df$Feature %in% trophic)
# # add_trophic <- data.frame(Feature = "trophic", delta_auc=sum(trophic_df$delta_auc), delta_sensitivity=sum(trophic_df$delta_sensitivity),
# #                           perm_auc = NA, base_auc = NA)
# # 
# # importance_df <- rbind(importance_df, add_lulc, add_hab, add_diet, add_activity, add_strata, add_order, add_phylo, add_trophic)
# 

