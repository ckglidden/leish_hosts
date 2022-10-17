######Figures for leish host paper - leish study effort
library(ggplot2)
library(tidyr)
library(dplyr)
library(sf)
library(raster)
library(foreach)
library(speciesgeocodeR)

###########
##figure 2: variable importance
###########
colors0 <- nationalparkcolors::park_palette("Badlands", 3)
colors <- c(colors0[c(1,3,2)], "grey")

mean_importance_leish <- read.csv("/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_leishmania/study effort supp/mean sum shap scores study effort leishmania july 20.csv")
#imp_df_leish <- subset(mean_importance_leish, mean_shap_score.lower.q.5. > 0)
imp_df_leish <- mean_importance_leish

importance_plot_leish <- ggplot(imp_df_leish, aes(x = reorder(Feature, mean_shap_score.upper.q.95.), y = mean_shap_score.mean)) +
  geom_point(size = 3) + xlab('feature') + ylab('mean |SHAP value|') + 
  ggtitle(expression(~italic(L. (Leishmania))~" study effort")) +
  geom_errorbar(aes(ymin = mean_shap_score.lower.q.5., ymax = mean_shap_score.upper.q.95.), position = "dodge", width = 0.4, size = 1.5) +
  coord_flip() +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        legend.key=element_blank()) 

ggsave("/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_leishmania/study effort supp/leish_studyEff_importance.png", importance_plot_leish, dpi=600)

###########
##figure 3: shap pdps
###########

pd_df_iter_leish <- read.csv('/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_leishmania/study effort supp/shapley dependence study leishmania sept 20 2022.csv'); names(pd_df_iter_leish)[2] <- 'Feature'
pd_df_iter_leish$row_code <- seq(1, nrow(pd_df_iter_leish), by = 1)

#split df into value & shap value; make each of them in long format; merge by feature and iteration
feature_value <- pd_df_iter_leish[,c(3:58,115, 116)] %>%
  tidyr::pivot_longer(cols = 1:55, #check this
               names_to = "feature",
               values_to = "value",
               values_drop_na = FALSE)
  
feature_shap <- pd_df_iter_leish[,c(59:114, 115, 116)] %>%
  tidyr::pivot_longer(cols = 1:55, #check this
                      names_to = "feature",
                      values_to = "shap_val",
                      values_drop_na = FALSE)
feature_shap$feature <- substr(feature_shap$feature,1,nchar(feature_shap$feature)-2)

shap_pdp_df <- left_join(feature_value, feature_shap, by = c("feature","iter", "row_code")) #merge dfs again

vars_to_plot <- c("LitterSize","GR_Area_km2", "GestationLen_d", "NeonateBodyMass_g", "WeaningAge_d", "GR_MaxLong_dd") #subset to significant variables

shap_pdp_plot_df <- shap_pdp_df[shap_pdp_df$feature %in% vars_to_plot, ]

trait_names <- list(
  "mean_tree_cover" = "% forest cvr",
  "mean_crop_cover" = "% crop cvr", 
  "mean_urban_cover" = "log(% urban cvr)", 
  "GestationLen_d" = "gest length (d)", 
  "bio10_temp_warmest_qt" = "tmp warmest qt (C)", 
  "LitterSize" = "litter size (n)",
  "GR_Area_km2" = "range area (km)",
  "PopulationDensity_n.km2" = "log(pop density) (n/km2)",
  "GR_MinLong_dd" = "min longitude (dd)",
  "GR_MaxLong_dd" = "max longitude (dd)"
)

trait_labeller <- function(variable,value){
  return(trait_names[value])
}

pdp_leish_se <- ggplot(shap_pdp_plot_df, aes(x=value, y=shap_val)) +
  #geom_point() +
  ggtitle(expression(~italic(L. (Leishmania))~" study effort")) +
  stat_smooth(aes(group=iter), color='lightgrey', method='loess', size=0.5, se=FALSE) + 
  stat_smooth(aes(), method='loess', size=2, se=FALSE) +
  geom_rug(data=shap_pdp_plot_df, aes(x=value), alpha=0.3, length=unit(0.05, "npc"), inherit.aes = FALSE, sides='b') +
  facet_wrap(~feature, scales='free', ncol=4) +
  theme_bw(base_size = 12) + 
  ylab("Shapley score") +
  theme(legend.position = "none")

ggsave("/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_leishmania/study effort supp/leish_studyEff_pdps.png", pdp_leish_se, dpi=600)

