######Figures for leish host paper
library(ggplot2)
library(tidyr)
library(dplyr)
library(sf)
library(raster)
library(foreach)
library(speciesgeocodeR)
###########
##figure 1 or 2?: model performance - shapley predictions
###########

nested_auc_leish <- read.csv('/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_leishmania/output/nested cv auc leishmania july 8 2022.csv')
Rmisc::CI(nested_auc_leish$x, 0.95)
#   upper      mean     lower 
#0.8912433 0.8856315 0.8800197 

target_shuffle_leish <- read.csv('/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_leishmania/output/nested cv target shuffled auc leishmania july 8 2022.csv')
Rmisc::CI(target_shuffle_leish$x, 0.95)
#upper      mean     lower 
#0.5464865 0.5372712 0.5280559 

in_sample_leish <- read.csv("/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_leishmania/output/in sample auc leishmania july 8 2022.csv")
Rmisc::CI(in_sample_leish$x, 0.95)
#upper      mean     lower 
#0.9752981 0.9748173 0.9743364

predictions_leish <- read.csv('/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_leishmania/output/leishmania summary shap predictions july 19.csv')

performance_plot_leish <- ggplot(predictions_leish, aes(x = shap_prediction.mean, y = stat(width*density), fill = as.factor(leish.infection))) +
  ylab('proportion of hosts') + xlab("mean sum shap value") + ggtitle("B. Leishmania (Leishmania)") +
  geom_histogram(binwidth= 0.3, position = "dodge", alpha=0.75) + 
  theme_bw(base_size = 14) + 
  theme(legend.position = c(0.5, 0.85)) + 
  labs(fill = "documented host") + 
  scale_fill_manual(values=RColorBrewer::brewer.pal(10, "RdYlBu")[c(2,8)])

ggsave("/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_leishmania/output/plots/predicted_prob_model_performance.pdf", 
       performance_plot_leish, dpi=600)


###########
##figure 2: variable importance
###########
colors0 <- nationalparkcolors::park_palette("Badlands", 3)
colors <- c(colors0[c(1,3,2)], "grey")

mean_importance_leish <- read.csv("/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_leishmania/output/mean sum shap scores leishmania july 18.csv")
imp_df_leish <- subset(mean_importance_leish, mean_shap_score.lower.q.5. > 0)

imp_df_leish$trait_type <- c("biogeography", "life-history", "biogeography", "biogeography", "biogeography",
                             "life-history", "biogeography", "phylogenetic distance", "study effort")

importance_plot_leish <- ggplot(imp_df_leish, aes(x = reorder(Feature, mean_shap_score.lower.q.5.), y = mean_shap_score.mean, color=trait_type)) +
  geom_point(size = 3) + xlab('feature') + ylab('mean |SHAP value|') + 
  ggtitle(expression(~bold(b.)~" "~italic(L. (Leishmania))~" ")) +
  geom_errorbar(aes(ymin = mean_shap_score.lower.q.5., ymax = mean_shap_score.upper.q.95.), position = "dodge", width = 0.4, size = 1.5) +
  scale_color_manual(values=colors, 'trait type') + 
  coord_flip() +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        legend.key=element_blank()) + 
  scale_x_discrete(breaks=imp_df_leish$Feature, 
                   labels=c("tmp warmest qt", "gestation len", "range area",  "max longitude", "min longitude", 
                   "litter size", "range % land cvr", "phylo distance", "study effort"))

importance_both <- cowplot::plot_grid(importance_plot_vianna, importance_plot_leish, ncol = 2)

ggsave("/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_leishmania/output/plots/trait_importance_performance.pdf", importance_plot_leish, dpi=600)
ggsave("/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_leishmania/output/plots/trait_importance_performance.png", importance_plot_leish, dpi=600)

ggsave("/Users/carolineglidden/Desktop/reservoir hosts - FINAL/all_trait_importance_performance.pdf", importance_both, dpi=600)

###########
##figure 3: shap pdps
###########

pd_df_iter_leish <- read.csv('/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_leishmania/output/shapley dependence leishmania july 11 2022.csv'); names(pd_df_iter_leish)[2] <- 'Feature'
pd_df_iter_leish$row_code <- seq(1, nrow(pd_df_iter_leish), by = 1)

#split df into value & shap value; make each of them in long format; merge by feature and iteration
feature_value <- pd_df_iter_leish[,c(4:60,118, 119)] %>%
  tidyr::pivot_longer(cols = 1:57, #check this
               names_to = "feature",
               values_to = "value",
               values_drop_na = FALSE)
  
feature_shap <- pd_df_iter_leish[,c(61:118, 119)] %>%
  tidyr::pivot_longer(cols = 1:57, #check this
                      names_to = "feature",
                      values_to = "shap_val",
                      values_drop_na = FALSE)
feature_shap$feature <- substr(feature_shap$feature,1,nchar(feature_shap$feature)-2)

shap_pdp_df <- left_join(feature_value, feature_shap, by = c("feature","iter", "row_code")) #merge dfs again

vars_to_plot <- c("mean_tree_cover","mean_crop_cover", "mean_urban_cover", "GestationLen_d", "LitterSize", "bio10_temp_warmest_qt", "GR_Area_km2", "GR_MaxLong_dd", "GR_MinLong_dd") #subset to significant variables

shap_pdp_plot_df <- shap_pdp_df[shap_pdp_df$feature %in% vars_to_plot, ]
names(imp_df_leish)[2] <- "feature"
trait_type_df <- imp_df_leish[,c("feature", "trait_type")]
trait_type_df <- rbind(trait_type_df, 
                       data.frame(feature = c("mean_tree_cover","mean_crop_cover", "mean_urban_cover"),
                                  trait_type = rep("biogeography", 3)))

shap_pdp_plot_df <- left_join(shap_pdp_plot_df, trait_type_df, by = "feature")

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

pdp_leish <- ggplot(shap_pdp_plot_df, aes(x=value, y=shap_val, color=trait_type)) +
  #geom_point() +
  ggtitle(expression(~bold(b.)~" "~italic(L. (Leishmania))~" ")) +
  stat_smooth(aes(group=iter), color='lightgrey', method='loess', size=0.5, se=FALSE) + 
  stat_smooth(aes(), method='loess', size=2, se=FALSE) +
  scale_color_manual(values=colors[1:2], 'trait type') + #fix colors
  geom_rug(data=shap_pdp_plot_df, aes(x=value), alpha=0.3, length=unit(0.05, "npc"), inherit.aes = FALSE, sides='b') +
  facet_wrap(~feature, scales='free', ncol=4, labeller = trait_labeller) +
  theme_bw(base_size = 12) + 
  ylab("Shapley score") +
  theme(legend.position = "none")

ggsave("/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_leishmania/output/plots/prob_pdps.pdf", pdp_leish, dpi=600)
ggsave("/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_leishmania/output/plots/prob_pdps.png", pdp_leish, dpi=600)

pdp_both <- cowplot::plot_grid(pdp_vianna, pdp_leish, ncol = 1)
ggsave("/Users/carolineglidden/Desktop/reservoir hosts - FINAL/shapley_pdps.pdf", pdp_both, dpi=600)
ggsave("/Users/carolineglidden/Desktop/reservoir hosts - FINAL/shapley_pdps.png", pdp_both, dpi=600)

###########
##figure 4: host per order
###########

summary_shap_predictions_leish <- read.csv("/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_leishmania/output/leishmania summary shap predictions july 19.csv")

abv_avg_leish <- subset(summary_shap_predictions_leish, shap_prediction.lower.q.5. > 0) #161 observations
new_leish <- subset(abv_avg_leish, leish.infection == 0) #98 unrecognized hosts

get_orders <- readRDS('/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_leishmania/cleaned data/analysis data/total data.rds')
orders <- unique(get_orders[,c("MSW05_Binomial", "MSW05_Order")])

new_leish <- left_join(new_leish, orders, by = "MSW05_Binomial")
tbl_leish <- as.data.frame(table(new_leish$MSW05_Order)); names(tbl_leish) <- c("order", "freq")
tbl_leish$subgenera <- rep("leishmania", 11)
tbl_leish$status <- rep("new", 11)

#old leish
old_leish <- subset(summary_shap_predictions_leish, leish.infection == 1)
old_leish <- left_join(old_leish, orders, by = "MSW05_Binomial")
tbl_old_leish <- as.data.frame(table(old_leish$MSW05_Order)); names(tbl_old_leish) <- c("order", "freq")
tbl_old_leish$subgenera <- rep("leishmania", 7)
tbl_old_leish$status <- rep("old", 7)

tbl_all_leish <- rbind(tbl_leish, tbl_old_leish)

# Barplot
order_colors <- c(nationalparkcolors::park_palette("GeneralGrant", 8), nationalparkcolors::park_palette("Saguaro", 3))

leish_ab_bar <- ggplot(tbl_all_leish, aes(x=order, y=freq, fill = order, color = as.factor(status))) +
  scale_fill_manual(values = order_colors) +  
  scale_color_manual(values = c("white", "white")) +
  ylim(0, 75) +
  ggtitle("b. L. (Leishmania)") +
  ylab("number of species") +
  geom_bar(stat = "identity") + 
  theme_classic(base_size = 12) +
  theme(legend.position="none",
        axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) 

leish_rel_bar <- ggplot(tbl_leish, aes(x=subgenera, y=freq, fill = order)) +
  scale_fill_manual(values = order_colors) +
  ggtitle("d.") +
  ylab("proportion of species") +
  geom_bar(position="fill", stat="identity") + 
  theme_classic(base_size = 12)

tbl_comb <- rbind(tbl_vianna, tbl_leish)
tbl_comb$order <- factor(tbl_comb$order, levels=c("Artiodactyla", "Carnivora", "Chiroptera", "Cingulata", "Didelphimorphia",
                                                      "Lagomorpha", "Perissodactyla", "Pilosa",          
                                                      "Primates", "Rodentia" , "Soricomorpha"))
tbl_comb$subgenera <- factor(tbl_comb$subgenera, levels=c("vianna", "leishmania"))

combined_rel_bar <- ggplot(tbl_comb, aes(x=subgenera, y=freq, fill = order)) +
  scale_fill_manual(values = order_colors) +
  ggtitle("c.") +
  ylab("proportion of species") +
  geom_bar(position="fill", stat="identity") + 
  theme_classic(base_size = 12) +
  theme(legend.position="none",
        axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) 


new_host_summary <- cowplot::plot_grid(vianna_ab_bar, combined_rel_bar, leish_ab_bar,
                                       ncol = 2, rel_widths = c(2,1))

ggsave("/Users/carolineglidden/Desktop/reservoir hosts - FINAL/new_host_summary.pdf", new_host_summary, dpi=600)
ggsave("/Users/carolineglidden/Desktop/reservoir hosts - FINAL/new_prop_host_summary.pdf", combined_rel_bar, dpi=600)







###########
##figure 4: host ranges
###########

###################
##figure 5 host locations
###################
library(sf)
library(RWmisc)
library(raster)
library(dplyr)

countries <- read_sf("/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_leishmania/raw data/leishmania_ranges/ne_10m_admin_0_countries/ne_10m_admin_0_countries.shp", layer="ne_10m_admin_0_countries")
latin_america <- c("Argentina", "Belize","Bolivia", "Brazil", "Chile", "Colombia", "Costa Rica", "Cuba", "Dominican Republic", 
                   "Ecuador", "El Salvador", "Guatemala", "Honduras", "Mexico", "Nicaragua", "Panama", "Paraguay", "Peru",
                   "Puerto Rico", "Uruguay", "Venezuela", "Guyana","Suriname","United States of America")
m_ca_sa <- countries %>% filter(countries$GEOUNIT %in% latin_america)
shape_v0 <- st_union(m_ca_sa)
shape <- st_crop(shape_v0, c(xmin= -127, ymin = -57, xmax = -31, ymax = 50))

ranges <- st_read("/Users/carolineglidden/Desktop/reservoir host - DRAFTS/habitattypes/cleaning feature collection/leishmania_animals_all")
for(i in 1:length(ntl)){
  ranges$BINOMIAL[ranges$BINOMIAL==ntl[i]] <- pantheria[i]
}

predicted_host_ranges_leish_v0 <- ranges[ranges$BINOMIAL %in% new_leish$MSW05_Binomial, ]
known_host_ranges_leish_v0 <- ranges[ranges$BINOMIAL %in% old_leish$MSW05_Binomial, ]

#r <- raster(as(shape, "Spatial"), ncols = 1000, nrows = 1000)
#rr <- rasterize(as(shape, "Spatial"), r, res=1000000, getCover = TRUE, progress = "text"); box <- extent(rr)

##crop rasters to Latin America, reduce shapefiles to one shapefile per species, count number of shapefile per grid cell, save for future plotting

#predicted hosts
predicted_host_ranges_leish_v1 <- st_crop(predicted_host_ranges_leish_v0, box)
predicted_host_ranges_leish_v2 <- foreach(i = new_leish$MSW05_Binomial, #make sure there is only one polygon per species
                                           .packages = c("sf", "tidyverse", "magrittr", "raster", "readxl"),
                                           .combine = "rbind") %do% {
                                             
                                             SHP <- predicted_host_ranges_leish_v1 %>%
                                               filter(BINOMIAL == i) %>%
                                               group_by(BINOMIAL) %>%
                                               summarise()
                                             
                                           }

predicted_host_ranges_leish_v2 = as(predicted_host_ranges_leish_v2, "Spatial")
predicted_leish_raster <- RangeRichness(predicted_host_ranges_leish_v2, res = 0.08)

writeRaster(predicted_leish_raster, '/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_leishmania/output/leish_predicted_host_richness.tif', format="GTiff", overwrite=TRUE)

##repeat for known hosts
#sf::sf_use_s2(FALSE)
known_host_ranges_leish_v1 <- st_crop(known_host_ranges_leish_v0, box)
known_host_ranges_leish_v2 <- foreach(i = old_leish$MSW05_Binomial, #make sure there is only one polygon per species
                                       .packages = c("sf", "tidyverse", "magrittr", "raster", "readxl"),
                                       .combine = "rbind") %do% {
                                         
                                         SHP <- known_host_ranges_leish_v1 %>%
                                           filter(BINOMIAL == i) %>%
                                           group_by(BINOMIAL) %>%
                                           summarise()
                                         
                                       }
known_host_ranges_leish_v2 = as(known_host_ranges_leish_v2, "Spatial")
known_leish_raster <- RangeRichness(known_host_ranges_leish_v2, res = 0.08)
writeRaster(known_leish_raster, '/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_leishmania/output/plots/leish_known_host_richness.tif', format="GTiff", overwrite=TRUE)

######################################
#############joint host ranges figure

vianna_known <- raster('/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_vianna/output/plots/vianna_known_host_richness.tif')
vianna_new <- raster('/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_vianna/output/plots/vianna_predicted_host_richness.tif')

leish_known <- raster('/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_leishmania/output/plots/leish_known_host_richness.tif')
leish_new <- raster('/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_leishmania/output/leish_predicted_host_richness.tif')

countries <- read_sf("/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_leishmania/raw data/leishmania_ranges/ne_10m_admin_0_countries/ne_10m_admin_0_countries.shp", layer="ne_10m_admin_0_countries")
latin_america <- c("Argentina", "Belize","Bolivia", "Brazil", "Chile", "Colombia", "Costa Rica", "Cuba", "Dominican Republic", 
                   "Ecuador", "El Salvador", "Guatemala", "Honduras", "Mexico", "Nicaragua", "Panama", "Paraguay", "Peru",
                   "Puerto Rico", "Uruguay", "Venezuela", "Guyana","Suriname","United States of America")
m_ca_sa <- countries %>% filter(countries$GEOUNIT %in% latin_america)
shape_v0 <- st_union(m_ca_sa)
america_shape <- st_crop(shape_v0, c(xmin= -127, ymin = -57, xmax = -31, ymax = 50))

leish_shp <- st_read('/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_leishmania/raw data/leishmania_distributions/leishmania_genus_range')
leish_shp_clipped <-st_intersection(st_make_valid(leish_shp), america_shape)
vianna_shp <- st_read('/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_vianna/raw data/leishmania_distributions/vianna_subgenus_range')
vianna_shp_clipped <-st_intersection(st_make_valid(vianna_shp), america_shape)

#####get host shape files
top_hosts <- c("Hylaeamys megacephalus", "Calomys callosus", "Dasypus novemcinctus", "Leopardus wiedii")
ranges <- st_read("/Users/carolineglidden/Desktop/reservoir host - DRAFTS/habitattypes/cleaning feature collection/leishmania_animals_all")
#^^make sure names are converted
top_host_shp <- ranges[ranges$BINOMIAL %in% top_hosts, ]

top_host_shp_v2 <- foreach(i = unique(top_host_shp$BINOMIAL), #make sure there is only one polygon per species
                           .packages = c("sf", "tidyverse", "magrittr", "raster", "readxl"),
                           .combine = "rbind") %do% {
                             
                             SHP <- top_host_shp %>%
                               filter(BINOMIAL == i) %>%
                               group_by(BINOMIAL) %>%
                               summarise()
                             
                           }

####vianna figures
library(cartomisc)
gplot_vianna_known<- gplot_data(vianna_known); gplot_vianna_known$value[gplot_vianna_known$value == 0] <- NA
gplot_vianna_predicted <- gplot_data(vianna_new); gplot_vianna_predicted$value[gplot_vianna_predicted$value == 0] <- NA

map_a <- ggplot() +
  ggtitle(expression(~bold(a.)~~italic(L. (Viannia))~'known')) +
  geom_tile(data = gplot_vianna_known, 
            aes(x = x, y = y, fill = value)) +
  scale_fill_viridis(name = "no. of viannia hosts", option = "C", na.value = "white") +
  geom_sf(data = america_shape, fill = "NA",
          colour = "black", size = 0.2) +
  geom_sf(data = vianna_shp_clipped, fill = NA,
          colour = "grey", size = 1.5) +
  theme_void() +
  theme(legend.position = c(0.2,0.2))

map_b <- ggplot() +
  ggtitle(expression(~bold(b.)~~italic(L. (Viannia))~'predicted')) +
  geom_tile(data = gplot_vianna_predicted, 
            aes(x = x, y = y, fill = value)) +
  scale_fill_viridis(option = "C", na.value = "white") +
  geom_sf(data = america_shape, fill = "NA",
          colour = "black", size = 0.2) +
  geom_sf(data = top_host_shp_v2[3,], fill = NA,
          colour = "darkseagreen1", size = 1) +
  geom_sf(data = top_host_shp_v2[4,], fill = NA,
          colour = "darkgreen", size = 1) +
  theme_void() +
  theme(legend.position = "none")


####leishmania figures
#library(cartomisc)
gplot_leish_known<- gplot_data(leish_known); gplot_leish_known$value[gplot_leish_known$value == 0] <- NA
gplot_leish_predicted <- gplot_data(leish_new); gplot_leish_predicted$value[gplot_leish_predicted$value == 0] <- NA

map_c <- ggplot() +
  ggtitle(expression(~bold(c.)~~italic(L. (Leishmania))~'known')) +
  geom_tile(data = gplot_leish_known, 
            aes(x = x, y = y, fill = value)) +
  scale_fill_viridis(name = "no. leish hosts", option = "C", na.value = "white") +
  geom_sf(data = america_shape, fill = "NA",
          colour = "black", size = 0.2) +
  geom_sf(data = leish_shp_clipped, fill = NA,
          colour = "grey", size = 1.5) +
  theme_void() +
  theme(legend.position = c(0.2,0.3))

map_d <- ggplot() +
  ggtitle(expression(~bold(d.)~~italic(L. (Leishmania))~'predicted')) +
  geom_tile(data = gplot_leish_predicted, 
            aes(x = x, y = y, fill = value)) +
  scale_fill_viridis(option = "C", na.value = "white") +
  geom_sf(data = america_shape, fill = "NA",
          colour = "black", size = 0.2) +
  geom_sf(data = top_host_shp_v2[1,], fill = NA,
          colour = "lightblue", size = 1) +
  geom_sf(data = top_host_shp_v2[2,], fill = NA,
          colour = "blue", size = 1) +
  theme_void() +
  theme(legend.position = "none")


gridExtra::grid.arrange(map_a, map_b, map_c, map_d, ncol = 2)


final_maps <- gridExtra::arrangeGrob(map_a, map_b, map_c, map_d, ncol = 2)
ggsave("/Users/carolineglidden/Desktop/reservoir hosts - FINAL/final_maps.pdf", final_maps, dpi = 300)
ggsave("/Users/carolineglidden/Desktop/reservoir hosts - FINAL/final_maps.png", final_maps, dpi = 300)


#############PCoA supplemental figure
phylo_importance <- read.csv("/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_leishmania/output/shapley variable importance leishmania july 11 2022.csv")
phylo_importance <- phylo_importance[phylo_importance$Feature %in% c("V1", "V2", "V3", "V4", "V5", "V6"), ]
mean_phylo_importance <- aggregate(mean_shap_score ~ Feature, data = phylo_importance, 
                                   FUN = function(x) c(mean = mean(x),
                                                       median = quantile(x, probs = 0.5),
                                                       lower.q = quantile(x, probs = 0.05),
                                                       upper.q = quantile(x, probs = 0.95)))

pcoa_data <- readRDS("/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_leishmania/cleaned data/analysis data/total data.rds")

pcoa_data <- pcoa_data[, c("MSW05_Order", "V1", "V2", "V3", "V4", "V5", "V6")]

#V1

v1 <- ggplot(pcoa_data, aes(x=V1)) + 
  geom_histogram(aes(fill = MSW05_Order)) # +
 # xlim(min(pcoa_data$V1), max(pcoa_data$V1)) +
 # facet_wrap(vars(MSW05_Order), scales = "free")

v2 <- ggplot(pcoa_data, aes(x=V2)) + 
  geom_histogram(aes(fill = MSW05_Order)) #+
  #xlim(min(pcoa_data$V2), max(pcoa_data$V2)) +
  #facet_wrap(vars(MSW05_Order), scales = "free")

v3 <- ggplot(pcoa_data, aes(x=V3)) + 
  geom_histogram(aes(fill = MSW05_Order)) #+
 # xlim(min(pcoa_data$V3), max(pcoa_data$V3)) +
#  facet_wrap(vars(MSW05_Order), scales = "free")

v4 <- ggplot(pcoa_data, aes(x=V4)) + 
  geom_histogram(aes(fill = MSW05_Order)) #+
  #xlim(min(pcoa_data$V4), max(pcoa_data$V4)) +
  #facet_wrap(vars(MSW05_Order), scales = "free")

v5 <- ggplot(pcoa_data, aes(x=V5)) + 
  geom_histogram(aes(fill = MSW05_Order)) #+
  #xlim(min(pcoa_data$V5), max(pcoa_data$V5)) +
  #facet_wrap(vars(MSW05_Order), scales = "free")

v6 <- ggplot(pcoa_data, aes(x=V6)) + 
  geom_histogram(aes(fill = MSW05_Order)) #+
  #xlim(min(pcoa_data$V6), max(pcoa_data$V6)) +
  #facet_wrap(vars(MSW05_Order), scales = "free")


pcoa_pdp <- read.csv("/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_leishmania/output/shapley dependence leishmania july 11 2022.csv")

pcoa_pdp$row_code <- seq(1, nrow(pcoa_pdp), by = 1)

#split df into value & shap value; make each of them in long format; merge by feature and iteration
feature_value <- pcoa_pdp[,c(4:60,118, 119)] %>%
  tidyr::pivot_longer(cols = 1:57, #check this
                      names_to = "feature",
                      values_to = "value",
                      values_drop_na = FALSE)

feature_shap <- pcoa_pdp[,c(61:118, 119)] %>%
  tidyr::pivot_longer(cols = 1:57, #check this
                      names_to = "feature",
                      values_to = "shap_val",
                      values_drop_na = FALSE)
feature_shap$feature <- substr(feature_shap$feature,1,nchar(feature_shap$feature)-2)

pcoa_pdp_df <- left_join(feature_value, feature_shap, by = c("feature","iter", "row_code")) #merge dfs again
pcoa_pdp_df <- pcoa_pdp_df[pcoa_pdp_df$feature %in% c("V1", "V2", "V3", "V4", "V5", "V6"), ]

pcoa_pdp_plot <- ggplot(pcoa_pdp_df, aes(x=value, y=shap_val)) +
  #geom_point() +
  #ggtitle(expression(~bold(b.)~" "~italic(L. (Leishmania))~" ")) +
  stat_smooth(aes(group=iter), color='lightgrey', method='loess', size=0.5, se=FALSE) + 
  stat_smooth(aes(), method='loess', size=2, se=FALSE) +
  #scale_color_manual(values=colors[1:2], 'trait type') + #fix colors
  #geom_rug(data=shap_pdp_plot_df, aes(x=value), alpha=0.3, length=unit(0.05, "npc"), inherit.aes = FALSE, sides='b') +
  facet_wrap(~feature, scales='free', ncol=4) +
  theme_bw(base_size = 12) + 
  ylab("Shapley score") +
  theme(legend.position = "none")



















#####extra code
hosts_to_plot <- subset(mean_predictions_leish, lower_ci > 0.5)
prediction_plot_leish <- predictions_leish[predictions_leish$MSW05_Binomial %in% hosts_to_plot$MSW05_Binomial,]
prediction_plot_leish <- subset(prediction_plot_leish, leish.infection < 1)

get_orders <- readRDS('/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_leishmania/cleaned data/analysis data/total data.rds')
orders <- unique(get_orders[,c("MSW05_Binomial", "MSW05_Order")])

prediction_plot_leish <- merge(prediction_plot_leish, orders, by = "MSW05_Binomial")

#pred_probs_leish <- ggplot(prediction_plot_leish, aes(as.numeric(predicted_probability), fill=MSW05_Order)) + 
#  geom_histogram(alpha=0.75) + facet_wrap(~MSW05_Order + MSW05_Binomial, ncol=4) + xlab('probability') +
#  scale_fill_brewer(palette = 'Spectral') +
#  theme_bw(base_size = 12) +
#  theme(legend.position = "none")
#ggsave("../output/plots/pred_probs.pdf", pred_probs, dpi=600)
#ggsave("../output/plots/pred_probs.png", pred_probs, dpi=600)

##pred plots w/ci intervals

hosts_to_plot_new_leish <- subset(hosts_to_plot, leish.infection < 1)
hosts_to_plot_old_leish <- subset(mean_predictions_leish, leish.infection > 0)

hosts_to_plot_new_leish <- merge(hosts_to_plot_new_leish, orders, by = "MSW05_Binomial")

pred_plot_new_leish <- ggplot(hosts_to_plot_new_leish, aes(x = fct_reorder(fct_reorder(MSW05_Binomial, mean_prob), MSW05_Order), y = mean_prob, color=MSW05_Order)) +
  geom_point(size = 3) + xlab('species') + ylab('probability') + ggtitle("A. Leishmania (Leishmania)") +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), position = "dodge", width = 0.4, size = 1.5) +
  scale_color_brewer(palette = 'Spectral') + 
  coord_flip() +
  theme_bw(base_size = 10) +
  theme(legend.position = c(0.85, 0.1))
ggsave("/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_leishmania/pred_probs_line_leish.pdf", pred_plot_new_leish, dpi=600)
ggsave("/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_leishmania/pred_probs_line_leish.png", pred_plot_new_leish, dpi=600)

combined_preds <- gridExtra::arrangeGrob(pred_plot_new_leish, pred_plot_vianna, nrow = 1)
                                         
ggsave("/Users/carolineglidden/Desktop/reservoir hosts - FINAL/predictions_new_hosts.png", combined_preds, dpi = 300)



#####new combined plots
all_hosts_to_plot <- unique(c(unique(hosts_to_plot_new_leish$MSW05_Binomial), unique(hosts_to_plot_new_vianna$MSW05_Binomial)))
###48 total hosts

final_leish_plot <- mean_predictions_leish[mean_predictions_leish$MSW05_Binomial %in% all_hosts_to_plot, ]; final_leish_plot$parasite <- "leishmania"

final_vianna_plot <- mean_predictions_vianna[mean_predictions_vianna$MSW05_Binomial %in% all_hosts_to_plot, ]; final_vianna_plot$parasite <- "vianna"
joint_plot_data <- rbind(final_leish_plot, final_vianna_plot)
joint_plot_data$sig <- ifelse(joint_plot_data$lower_ci > 0.5, "sig", "non")
joint_plot_data <- merge(joint_plot_data, orders, by = "MSW05_Binomial")
joint_plot_data$numeric_order <- ifelse(joint_plot_data$MSW05_Order == "Rodentia", 3,
                                        ifelse(joint_plot_data$MSW05_Order == "Chiroptera", 2,
                                               ifelse(joint_plot_data$MSW05_Order == "Pilosa", 3,
                                                      ifelse(joint_plot_data$MSW05_Order == "Didelphimorphia", 4,
                                                             ifelse(joint_plot_data$MSW05_Order == "Carnivora", 1,
                                                                    ifelse(joint_plot_data$MSW05_Order == "Cingulata", 6,
                                                                           ifelse(joint_plot_data$MSW05_Order == "Artiodactyla", 7,8)))))))


#write.csv(joint_plot_data, "/Users/carolineglidden/Desktop/reservoir hosts - FINAL/prediction_plot_data.csv")

pointplot_joined <- ggplot(joint_plot_data, aes(x = fct_reorder(fct_reorder(MSW05_Binomial, mean_prob), numeric_order, .desc=TRUE), 
                            y = mean_prob, color = interaction(parasite, leish.infection), shape = parasite)) +
  geom_point() + coord_flip() + xlab("species") +
  ylab("probability") +
  ylim(0.5, 0.9) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci, group=interaction(parasite, leish.infection)))+
  scale_color_brewer(palette="Paired", values = c(1,2,5,6)) +
    theme_minimal() +
  theme(legend.position = "none")
 # facet_wrap(~MSW05_Order, scales = 'free')
ggsave("/Users/carolineglidden/Desktop/reservoir hosts - FINAL/predictions_new_hosts_pointplot.pdf", pointplot_joined, dpi = 300)



####try ridgeline plot
# library
library(ggridges)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
my_palette = c(brewer.pal(6, "Paired")[c(1,2,5,6)])

ridge_leish_plot <- predictions_leish[predictions_leish$MSW05_Binomial %in% all_hosts_to_plot, ]; ridge_leish_plot$parasite <- "leishmania"
ridge_vianna_plot <- predictions_vianna[predictions_vianna$MSW05_Binomial %in% all_hosts_to_plot, ]; ridge_vianna_plot$parasite <- "vianna"
ridge_plot <- rbind(ridge_leish_plot, ridge_vianna_plot)
ridge_plot <- left_join(ridge_plot, orders, by = "MSW05_Binomial"); ridge_plot <- left_join(ridge_plot, mean_predictions_leish[1:2], by = "MSW05_Binomial")
ridge_plot$numeric_order <- ifelse(ridge_plot$MSW05_Order == "Rodentia", 3,
                                                    ifelse(ridge_plot$MSW05_Order == "Chiroptera", 2,
                                                           ifelse(ridge_plot$MSW05_Order == "Pilosa", 3,
                                                                  ifelse(ridge_plot$MSW05_Order == "Didelphimorphia", 4,
                                                                         ifelse(ridge_plot$MSW05_Order == "Carnivora", 1,
                                                                                ifelse(ridge_plot$MSW05_Order == "Cingulata", 6,
                                                                                       ifelse(ridge_plot$MSW05_Order == "Artiodactyla", 7,8)))))))

#write.csv(ridge_plot, "/Users/carolineglidden/Desktop/reservoir hosts - FINAL/prediction_ridge_plot_data.csv")

ridge_figure <- ggplot(ridge_plot, aes(y = fct_reorder(fct_reorder(MSW05_Binomial, mean_prob), numeric_order, .desc = TRUE), 
                       x = predicted_probability, 
                       fill=interaction(parasite, leish.infection),
                       color=interaction(parasite, leish.infection))) +
  geom_density_ridges(rel_min_height = 0.05, alpha=0.5, scale = 0.9) +
  theme_ridges(font_size = 9) + 
  geom_vline(xintercept = 0.5) +
  ylab("species") + 
  xlab("predicted probability") +
  xlim(0.04, 0.93) + 
  scale_fill_manual(values = my_palette) +
  scale_color_manual(values = my_palette) +
  theme(legend.position = "none",
        axis.title=element_text(size=9),
        axis.text = element_text(size=9))

ggsave("/Users/carolineglidden/Desktop/reservoir hosts - FINAL/predictions_new_hosts_ridgeplot.pdf", 
       device = "pdf",
       ridge_figure, dpi = 300)


max(ridge_plot$predicted_probability)






###################
##figure 5 host locations
###################
library(sf)
library(RWmisc)
library(raster)

countries <- read_sf("/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_leishmania/raw data/leishmania_ranges/ne_10m_admin_0_countries/ne_10m_admin_0_countries.shp", layer="ne_10m_admin_0_countries")
latin_america <- c("Argentina", "Belize","Bolivia", "Brazil", "Chile", "Colombia", "Costa Rica", "Cuba", "Dominican Republic", 
                   "Ecuador", "El Salvador", "Guatemala", "Honduras", "Mexico", "Nicaragua", "Panama", "Paraguay", "Peru",
                   "Puerto Rico", "Uruguay", "Venezuela", "Guyana","Suriname","", "United States of America")
m_ca_sa <- countries %>% filter(countries$GEOUNIT %in% latin_america)
shape <- st_union(m_ca_sa)

ranges <- st_read("/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_leishmania/raw data/leishmania_animals_all_clipped")
predicted_host_ranges <- ranges[ranges$BINOMIAL %in% hosts_to_plot_new_leish$MSW05_Binomial, ]
known_host_ranges <- ranges[ranges$BINOMIAL %in% hosts_to_plot_old_leish$MSW05_Binomial, ]

r <- raster(as(shape, "Spatial"), ncols = 1000, nrows = 1000)
rr <- rasterize(as(shape, "Spatial"), r, res=10000, getCover = TRUE, progress = "text"); extent(rr) <- extent(predicted_host_ranges)

predicted_host_raster <- overlap.weight(rr, predicted_host_ranges, count = TRUE, warn = TRUE)
writeRaster(predicted_host_raster, '/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_leishmania/output/plots/predicted_host_richness.tif', format="GTiff", overwrite=TRUE)

known_host_raster <- overlap.weight(rr, known_host_ranges, count = TRUE, warn = TRUE)
know_host_raster <- crop(known_host_raster, extent(predicted_host_raster))
writeRaster(known_host_raster, '/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_leishmania/output/plots/known_host_richness.tif', format="GTiff", overwrite=TRUE)

plot(predicted_host_raster)  
plot(known_host_raster)  

library(RColorBrewer)

tiff('range_map_draft_plot.tiff', width = 5, height = 4, units = 'in',                       # fixed
     compression = c( "lzw"), res=300)
par(mfrow=c(1,2))
plot(known_host_raster,
     col=brewer.pal(n=5, name='BuPu'),
     main="Known Hosts",
     axes = FALSE, box = FALSE,
     xlim = c(-124.587, -32),
     ylim = c(-55.02719, 50.65755),
     cex.main = 1, cex.lab=1)
plot(shape, bg="transparent", 
     xlim = c(-124.587, -32),
     ylim = c(-55.02719, 50.65755), add=TRUE)

plot(predicted_host_raster,
     col=brewer.pal(n=5, name='BuPu'),
     main="Predicted Hosts",
     axes = FALSE, box = FALSE,
     xlim = c(-124.587, -32),
     ylim = c(-55.02719, 50.65755),
     cex.main = 1, cex.lab=1)
plot(shape, bg="transparent", 
     xlim = c(-124.587, -32),
     ylim = c(-55.02719, 50.65755),
     add=TRUE)
dev.off()

library(cartomisc)
gplot_known<- gplot_data(known_host_raster)
gplot_predicted <- gplot_data(predicted_host_raster)

ggplot() +
  geom_sf(data = shape, fill = "white",
          colour = "black", size = 0.2) +
  #xlim = c(-124.587, -32) +
 # ylim = c(-55.02719, 50.65755) +
  geom_tile(data = gplot_known, 
            aes(x = x, y = y, fill = value)) +
  scale_fill_brewer(palette = "BuPu")

ggplot() +
  geom_sf(data = shape, fill = "white",
          colour = "black", size = 0.2) +
  geom_tile(data = gplot_predicted, 
            aes(x = x, y = y, fill = value)) +
  #xlim = c(-124.587, -32) +
  #ylim = c(-55.02719, 50.65755) +
  scale_fill_brewer(palette = "BuPu")














####old code 
predictions_leish$predicted_probability <- as.numeric(predictions_leish$predicted_probability)
predictions_leish$MSW05_Binomial <- as.factor(predictions_leish$MSW05_Binomial)

mean_predictions_leish <- aggregate(as.numeric(predicted_probability) ~ as.factor(MSW05_Binomial), data = cv_preds, 
                                    FUN = function(x) c(mean = mean(x), 
                                                        lowerCI = Rmisc::CI(x, 0.95)[3], 
                                                        upperCI = Rmisc::CI(x, 0.95)[1],
                                                        lower.q = quantile(x, probs = 0.025),
                                                        upper.q = quantile(x, probs = 0.975)))
mean_predictions_leish <- do.call(data.frame, mean_predictions_leish); names(mean_predictions_leish) <- c('MSW05_Binomial','mean_prob','lower_ci','upper_ci', "lower_q", "upper_q")

mean_predictions_leish <- merge(mean_predictions_leish, analysis_data[,c('MSW05_Binomial','leish.infection')], by = 'MSW05_Binomial')
mean_predictions_leish <- unique(mean_predictions_leish)


predicted <- ifelse(mean_predictions_leish$predicted_probability > 0.7390197, 1, 0)
known <- mean_predictions_leish$leish.infection
caret::confusionMatrix(as.factor(predicted), as.factor()) #90th percentile: sensitivity = 0.94, specificity = 0.91


#pd_df_iter <- merge(pd_df_iter, imp_df[c(2,3,8)], by = 'Feature')
new_labels_leish <- data.frame(Feature = feature_names_leish, 
                               new_feat = c("gestation length (d)", "study effort", "mean temp (C)", 
                                            "litter size", "range area (km2)", "mean vpd",
                                            "home range (km2)",  "minimum long (dd)",
                                            "phylo dim 1", "phylo dim 2", "phylo dim 3", "phylo dim 4", "phylo dim 5", "phlyo dim 6",
                                            "% crop cover ", "% grass cover", "log(% urban cover)", "% tree cover (0.006)"),
                               trait_type = c('life-history', 'study effort', 
                                              'biogeography', 'life-history', 'biogeography', 'biogeography', 'life-history',
                                              'biogeography', rep("phylogeny",6), rep("biogeography",4)))
pd_df_iter_leish <- merge(pd_df_iter_leish, new_labels_leish, by = 'Feature')

#trait_data <- readRDS("/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_leishmania/cleaned data/analysis data/trait_data_small.rds")

#pd_df_iter <- pd_df_iter %>%
#  mutate(across(Feature, factor, levels=imp_df$Feature))

rug_data <- trait_data[,feature_names_leish]
rug_data <- rug_data %>%
  tidyr::pivot_longer(cols= 1:18,names_to = "Feature", values_to = "value")
rug_data <- as.data.frame(rug_data); rug_data <- merge(rug_data, new_labels_leish, by = 'Feature')


range_pdp <- pd_df_iter_leish[pd_df_iter_leish$Feature %in% c("mean_crop_cover","mean_grass_cover","mean_urban_cover", "mean_tree_cover"), ]
range_rug <- rug_data[rug_data$Feature %in% c("mean_crop_cover","mean_grass_cover","mean_urban_cover", "mean_tree_cover"),]

pdp_leish_range <- ggplot(range_pdp, aes(x=value, y=yhat)) + 
  stat_smooth(aes(group=iter), color='lightgrey', method='loess', size=0.5, se=FALSE) + 
  stat_smooth(aes(), method='loess', size=2, se=FALSE, color=colors) +
  #scale_color_manual(values=colors, 'trait type') + 
  geom_rug(data=range_rug, aes(x=value), alpha=0.3, length=unit(0.05, "npc"), inherit.aes = FALSE, sides='b') +
  facet_wrap(~new_feat, scales='free', ncol=4) +
  #scale_color_manual(values=c('darkgrey',colors[c(3,5,8)])) +
  theme_bw(base_size = 11) + 
  theme(legend.position = "none")

ggsave("/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_leishmania/output/plots/pdps.pdf", pdp_leish, dpi=600)
ggsave("/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_leishmania/output/plots/pdps.png", pdp_leish, dpi=600)

ggsave("/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_leishmania/output/plots/pdps_range.pdf", pdp_leish_range, dpi=600)
ggsave("/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_leishmania/output/plots/pdps_range.png", pdp_leish_range, dpi=600)

