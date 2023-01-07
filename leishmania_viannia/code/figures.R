######Figures for leish host paper
library(ggplot2)
library(forcats)
library(dplyr)
library(tidyr)
library(sf)

###########
##summary 1: model performance
###########

nested_auc <- read.csv('/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_viannia/output/nested cv auc viannia july 8 2022.csv')
Rmisc::CI(nested_auc$x, 0.95)
#upper      mean     lower 
#0.8561613 0.8468799 0.8375984

target_shuffle <- read.csv('/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_viannia/output/nested cv target shuffled auc viannia july 8 2022.csv')
#t.test(target_shuffle$x, mu=0.5)
Rmisc::CI(target_shuffle$x, 0.95)
#upper      mean     lower 
#0.5506404 0.5406717 0.5307030  

in_sample <- read.csv("/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_viannia/output/in sample auc viannia july 8 2022.csv")
Rmisc::CI(in_sample$x, 0.95)
#upper      mean     lower 
#0.9707527 0.9701160 0.9694793 


###########
##figure 2: variable importance
###########
colors0 <- nationalparkcolors::park_palette("Badlands", 3)
colors <- c(colors0[c(1,3,2)], "grey")

mean_importance_viannia <- read.csv("/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_viannia/output/mean sum shap scores leishmania july 21.csv")
imp_df_viannia <- subset(mean_importance_viannia, mean_shap_score.lower.q.5. > 0)

imp_df_viannia$trait_type <- c("biogeography", "life-history", "biogeography", "biogeography", "phylogenetic distance", "life-history", "study effort")

importance_plot_viannia <- ggplot(imp_df_viannia, aes(x = reorder(Feature, mean_shap_score.lower.q.5.), y = mean_shap_score.mean, color=trait_type)) +
  geom_point(size = 3) + xlab('feature') + ylab('mean |SHAP value|') + 
  ggtitle(expression(~bold(a.)~" "~italic(L. (viannia))~" ")) +
  geom_errorbar(aes(ymin = mean_shap_score.lower.q.5., ymax = mean_shap_score.upper.q.95.), position = "dodge", width = 0.4, size = 1.5) +
  scale_color_manual(values=colors, 'trait type') + 
  coord_flip() +
  theme_bw(base_size = 12) +
  theme(legend.position = c(0.9, 0.1),
        legend.key=element_blank()) + scale_x_discrete(breaks=imp_df_viannia$Feature, 
                                                       labels=c("tmp warmest qt", "gestation len", "min longitude", 
                                                                "range % land cvr", "phylo distance", "pop density", "study effort"))
                                                              

#ggsave("/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_viannia/output/plots/trait_importance_performance.pdf", importance_plot_viannia, dpi=600)
#ggsave("/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_viannia/output/plots/trait_importance_performance.png", importance_plot_viannia, dpi=600)


###########
##figure 3: shap pdps
###########

pd_df_iter_viannia <- read.csv('/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_viannia/output/shapley dependence viannia july 21 2022.csv'); names(pd_df_iter_viannia)[2] <- 'Feature'
pd_df_iter_viannia$row_code <- seq(1, nrow(pd_df_iter_viannia), by = 1)

#split df into value & shap value; make each of them in long format; merge by feature and iteration
feature_value <- pd_df_iter_viannia[,c(4:58,114, 115)] %>%
  tidyr::pivot_longer(cols = 1:55, #check this
                      names_to = "feature",
                      values_to = "value",
                      values_drop_na = FALSE)

feature_shap <- pd_df_iter_viannia[,c(59:115)] %>%
  tidyr::pivot_longer(cols = 1:55, #check this
                      names_to = "feature",
                      values_to = "shap_val",
                      values_drop_na = FALSE)
feature_shap$feature <- substr(feature_shap$feature,1,nchar(feature_shap$feature)-2)

shap_pdp_df <- left_join(feature_value, feature_shap, by = c("feature","iter", "row_code")) #merge dfs again

vars_to_plot <- c("GR_MinLong_dd", "mean_tree_cover","mean_crop_cover", "mean_urban_cover", "GestationLen_d", "bio10_temp_warmest_qt", "PopulationDensity_n.km2") #subset to significant variables

shap_pdp_plot_df <- shap_pdp_df[shap_pdp_df$feature %in% vars_to_plot, ]
names(imp_df_viannia)[2] <- "feature"
trait_type_df <- imp_df_viannia[,c("feature", "trait_type")]
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
  "GR_MinLong_dd" = "min longitude (dd)"
)

trait_labeller <- function(variable,value){
  return(trait_names[value])
}

pdp_viannia <- ggplot(shap_pdp_plot_df, aes(x=value, y=shap_val, color=trait_type)) +
  #geom_point() +
  ggtitle(expression(~bold(a.)~" "~italic(L. (viannia))~" ")) +
  stat_smooth(aes(group=iter), color='lightgrey', method='loess', size=0.5, se=FALSE) + 
  stat_smooth(aes(), method='loess', size=2, se=FALSE) +
  scale_color_manual(values=colors[1:2], 'trait type') + #fix colors
  geom_rug(data=shap_pdp_plot_df, aes(x=value), alpha=0.3, length=unit(0.05, "npc"), inherit.aes = FALSE, sides='b') +
  facet_wrap(~feature, scales='free', ncol=4, labeller = trait_labeller) +
  theme_bw(base_size = 12) + 
  ylab("Shapley score") +
  theme(legend.position = "none")

ggsave("/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_viannia/output/plots/prob_pdps.pdf", pdp_viannia, dpi=600)
ggsave("/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_viannia/output/plots/prob_pdps.png", pdp_viannia, dpi=600)


###########
##figure 4: host per order
###########

summary_shap_predictions_viannia <- read.csv("/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_viannia/output/leishmania summary shap predictions july 21.csv")

abv_avg_viannia <- subset(summary_shap_predictions_viannia, shap_prediction.lower.q.5. > 0) #126 observations
new_viannia <- subset(abv_avg_viannia, leish.infection == 0) #83 unrecognized hosts

get_orders <- readRDS('/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_leishmania/cleaned data/analysis data/total data.rds')
orders <- unique(get_orders[,c("MSW05_Binomial", "MSW05_Order")])

new_viannia <- left_join(new_viannia, orders, by = "MSW05_Binomial")
tbl_viannia <- as.data.frame(table(new_viannia$MSW05_Order)); names(tbl_viannia) <- c("order", "freq")
tbl_viannia <- rbind(tbl_viannia, data.frame(order = "Perissodactyla", freq = 0))
tbl_viannia$subgenera <- rep("viannia", 11)
tbl_viannia$order <- factor(tbl_viannia$order, levels=c("Artiodactyla", "Carnivora", "Chiroptera", "Cingulata", "Didelphimorphia",
                                                      "Lagomorpha", "Perissodactyla", "Pilosa",          
                                                      "Primates", "Rodentia" , "Soricomorpha"))
tbl_viannia$status <- rep("new", 11)

#old viannia
old_viannia <- subset(summary_shap_predictions_viannia, leish.infection == 1)
old_viannia <- left_join(old_viannia, orders, by = "MSW05_Binomial")
tbl_old_viannia <- as.data.frame(table(old_viannia$MSW05_Order)); names(tbl_old_viannia) <- c("order", "freq")
tbl_old_viannia$subgenera <- rep("viannia", 9)
tbl_old_viannia$status <- rep("old", 9)

tbl_all_viannia <- rbind(tbl_viannia, tbl_old_viannia)

# Barplot
order_colors <- c(nationalparkcolors::park_palette("GeneralGrant", 8), nationalparkcolors::park_palette("Saguaro", 3))

viannia_ab_bar <- ggplot(tbl_all_viannia, aes(x=order, y=freq, fill = order, color = as.factor(status))) +
  scale_fill_manual(values = order_colors) +
  scale_color_manual(values = c("white", "white")) +
  ylim(0, 75) +
  ggtitle("a. L. (viannia)") +
  ylab("number of species") +
  geom_bar(stat = "identity") + 
  theme_classic(base_size = 12) +
  theme(legend.position="none",
        axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) 

viannia_rel_bar <- ggplot(tbl_viannia, aes(x=subgenera, y=freq, fill = order)) +
  scale_fill_manual(values = order_colors) +
  ggtitle("b.") +
  ylab("proportion of species") +
  geom_bar(position="fill", stat="identity") + 
  theme_classic(base_size = 12)

#gridExtra::grid.arrange(viannia_ab_bar, viannia_rel_bar, ncol = 2)



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
shape <- st_crop(shape, c(xmin= -127, ymin = -57, xmax = -31, ymax = 50))

ranges <- st_read("/Users/carolineglidden/Desktop/reservoir host - DRAFTS/habitattypes/cleaning feature collection/leishmania_animals_all")
for(i in 1:length(ntl)){
  ranges$BINOMIAL[ranges$BINOMIAL==ntl[i]] <- pantheria[i]
}

predicted_host_ranges_viannia_v0 <- ranges[ranges$BINOMIAL %in% new_viannia$MSW05_Binomial, ]
known_host_ranges_viannia_v0 <- ranges[ranges$BINOMIAL %in% old_viannia$MSW05_Binomial, ]

#r <- raster(as(shape, "Spatial"), ncols = 1000, nrows = 1000)
#rr <- rasterize(as(shape, "Spatial"), r, res=1000000, getCover = TRUE, progress = "text"); box <- extent(rr)

##crop rasters to Latin America, reduce shapefiles to one shapefile per species, count number of shapefile per grid cell, save for future plotting

#predicted hosts
predicted_host_ranges_viannia_v1 <- st_crop(predicted_host_ranges_viannia_v0, box)
predicted_host_ranges_viannia_v2 <- foreach(i = new_viannia$MSW05_Binomial, #make sure there is only one polygon per species
                                                 .packages = c("sf", "tidyverse", "magrittr", "raster", "readxl"),
                                                 .combine = "rbind") %do% {
  
  SHP <- predicted_host_ranges_viannia_v1 %>%
    filter(BINOMIAL == i) %>%
    group_by(BINOMIAL) %>%
    summarise()
                                                 }

predicted_host_ranges_viannia_v2 = as(predicted_host_ranges_viannia_v2, "Spatial")
predicted_viannia_raster <- RangeRichness(predicted_host_ranges_viannia_v2, res = 0.08)
writeRaster(predicted_viannia_raster, '/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_viannia/output/plots/viannia_predicted_host_richness.tif', format="GTiff", overwrite=TRUE)

##repeat for known hosts
#sf::sf_use_s2(FALSE)
known_host_ranges_viannia_v1 <- st_crop(known_host_ranges_viannia_v0, box)
known_host_ranges_viannia_v2 <- foreach(i = old_viannia$MSW05_Binomial, #make sure there is only one polygon per species
                                           .packages = c("sf", "tidyverse", "magrittr", "raster", "readxl"),
                                           .combine = "rbind") %do% {
                                             
                                             SHP <- known_host_ranges_viannia_v1 %>%
                                               filter(BINOMIAL == i) %>%
                                               group_by(BINOMIAL) %>%
                                               summarise()
                                             
                                           }

known_host_ranges_viannia_v2 = as(known_host_ranges_viannia_v2, "Spatial")
known_viannia_raster <- RangeRichness(known_host_ranges_viannia_v2, res = 0.08)
writeRaster(known_viannia_raster, '/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_viannia/output/plots/viannia_known_host_richness.tif', format="GTiff", overwrite=TRUE)




#############PCoA supplemental figure
phylo_importance <- read.csv("/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_viannia/output/shapley variable importance viannia july 21 2022.csv")
phylo_importance <- phylo_importance[phylo_importance$Feature %in% c("V1", "V2", "V3", "V4", "V5", "V6"), ]
mean_phylo_importance <- aggregate(mean_shap_score ~ Feature, data = phylo_importance, 
                                   FUN = function(x) c(mean = mean(x),
                                                       median = quantile(x, probs = 0.5),
                                                       lower.q = quantile(x, probs = 0.05),
                                                       upper.q = quantile(x, probs = 0.95))) ##plot V2 & V6
mean_phylo_importance <- do.call(data.frame, mean_phylo_importance)


phylo_dimensions_imp_plot <- ggplot(mean_phylo_importance, aes(x = reorder(Feature, mean_shap_score.lower.q.5.), y = mean_shap_score.mean)) +
  geom_point(size = 3) + xlab('PCoA dimension') + ylab('mean |SHAP value|') + 
  #ggtitle(expression(~bold(a.)~" "~italic(L. (viannia))~" ")) +
  geom_errorbar(aes(ymin = mean_shap_score.lower.q.5., ymax = mean_shap_score.upper.q.95.), position = "dodge", width = 0.4, size = 1.5) +
  #scale_color_manual(values=colors, 'trait type') + 
  coord_flip() +
  theme_bw(base_size = 12)
ggsave("viannia_phylo_trait_imp.png", phylo_dimensions_imp_plot)

pcoa_data <- readRDS("/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_leishmania/cleaned data/analysis data/total data.rds")

pcoa_data <- pcoa_data[, c("MSW05_Order", "V2", "V6")]

#V1

# v1 <- ggplot(pcoa_data, aes(x=V1)) + 
#   geom_histogram(aes(fill = MSW05_Order)) # +
# xlim(min(pcoa_data$V1), max(pcoa_data$V1)) +
# facet_wrap(vars(MSW05_Order), scales = "free")

v2 <- ggplot(pcoa_data, aes(x=V2)) + 
  geom_histogram(aes(fill = MSW05_Order)) #+
#xlim(min(pcoa_data$V2), max(pcoa_data$V2)) +
#facet_wrap(vars(MSW05_Order), scales = "free")

# v3 <- ggplot(pcoa_data, aes(x=V3)) + 
#   geom_histogram(aes(fill = MSW05_Order)) #+
# xlim(min(pcoa_data$V3), max(pcoa_data$V3)) +
#  facet_wrap(vars(MSW05_Order), scales = "free")

# v4 <- ggplot(pcoa_data, aes(x=V4)) + 
#   geom_histogram(aes(fill = MSW05_Order)) #+
#xlim(min(pcoa_data$V4), max(pcoa_data$V4)) +
#facet_wrap(vars(MSW05_Order), scales = "free")

# v5 <- ggplot(pcoa_data, aes(x=V5)) + 
#   geom_histogram(aes(fill = MSW05_Order)) #+
#xlim(min(pcoa_data$V5), max(pcoa_data$V5)) +
#facet_wrap(vars(MSW05_Order), scales = "free")

v6 <- ggplot(pcoa_data, aes(x=V6)) + 
  geom_histogram(aes(fill = MSW05_Order)) #+
#xlim(min(pcoa_data$V6), max(pcoa_data$V6)) +
#facet_wrap(vars(MSW05_Order), scales = "free")


pcoa_pdp <- read.csv("/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_viannia/output/shapley dependence viannia july 21 2022.csv")

pcoa_pdp$row_code <- seq(1, nrow(pcoa_pdp), by = 1)

#split df into value & shap value; make each of them in long format; merge by feature and iteration
feature_value <- pcoa_pdp[,c(4:58,114, 115)] %>%
  tidyr::pivot_longer(cols = 1:55, #check this
                      names_to = "feature",
                      values_to = "value",
                      values_drop_na = FALSE)

feature_shap <- pcoa_pdp[,c(59:115)] %>%
  tidyr::pivot_longer(cols = 1:55, #check this
                      names_to = "feature",
                      values_to = "shap_val",
                      values_drop_na = FALSE)
feature_shap$feature <- substr(feature_shap$feature,1,nchar(feature_shap$feature)-2)

pcoa_pdp_df <- left_join(feature_value, feature_shap, by = c("feature","iter", "row_code")) #merge dfs again

pcoa_pdp_df <- pcoa_pdp_df[pcoa_pdp_df$feature %in% c("V2", "V6"), ]

pcoa_pdp_plot <- ggplot(pcoa_pdp_df, aes(x=value, y=shap_val)) +
  #geom_point() +
  #ggtitle(expression(~bold(b.)~" "~italic(L. (Leishmania))~" ")) +
  stat_smooth(aes(group=iter), color='lightgrey', method='loess', size=0.5, se=FALSE) + 
  stat_smooth(aes(), method='loess', size=2, se=FALSE) +
  #scale_color_manual(values=colors[1:2], 'trait type') + #fix colors
  #geom_rug(data=shap_pdp_plot_df, aes(x=value), alpha=0.3, length=unit(0.05, "npc"), inherit.aes = FALSE, sides='b') +
  facet_wrap(~feature, scales='free', ncol=1) +
  theme_bw(base_size = 12) + 
  ylab("Shapley score") +
  theme(legend.position = "none")

#V2 = -25
#V6 = 0

pcoa_data$v2_cutoff <- ifelse(pcoa_data$V2 > -25, 1, 0)
pcoa_data$v6_cutoff <- ifelse(pcoa_data$V6 > 0, 1, 0)

v2_table_below <- as.data.frame(table(pcoa_data[pcoa_data$v2_cutoff == 0,]$MSW05_Order)); v2_table_below$threshold <- "below"
v2_table_above <- as.data.frame(table(pcoa_data[pcoa_data$v2_cutoff == 1,]$MSW05_Order)); v2_table_above$threshold <- "above"
v2_table <- rbind(v2_table_below, v2_table_above)
v2_table$dimension <- "V2"

v6_table_below <- as.data.frame(table(pcoa_data[pcoa_data$v6_cutoff == 0,]$MSW05_Order)); v6_table_below$threshold <- "below"
v6_table_above <- as.data.frame(table(pcoa_data[pcoa_data$v6_cutoff == 1,]$MSW05_Order)); v6_table_above$threshold <- "above"
v6_table <- rbind(v6_table_below, v6_table_above)
v6_table$dimension <- "V6"

phlyo_table <- rbind(v2_table, v6_table)

phlyo_table <- phlyo_table[-c(7, 18),]

colorz <- viridis::turbo(13)

phylo_plot <- ggplot(phlyo_table, aes(x=threshold, y=Freq, fill = Var1)) +
  scale_fill_manual(values = colorz) +
  #ggtitle("b.") +
  ylab("proportion of species") +
  geom_bar(position="fill", stat="identity") + 
  facet_wrap(~dimension, scales='free', ncol=1) +
  theme_classic(base_size = 12)


combined_pcoa <- gridExtra::arrangeGrob(pcoa_pdp_plot, phylo_plot, ncol = 2)

ggsave("pcoa_full_supplementary_figure.pdf", combined_pcoa)




























#####try plotting rasters again
library(raster); library(rasterVis); library(ggplot2); library(viridis)

viannia_predicted <- raster('/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_viannia/output/plots/predicted_host_richness.tif')
viannia_predicted[viannia_predicted$predicted_host_richness == 0, ] <- NA
viannia_known <- raster('/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_viannia/output/plots/known_host_richness.tif')
viannia_known[viannia_known$known_host_richness == 0, ] <- NA

leish_predicted <- raster('/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_leishmania/output/plots/predicted_host_richness.tif')
leish_predicted[leish_predicted$predicted_host_richness == 0, ] <- NA
leish_known <- raster('/Users/carolineglidden/Desktop/reservoir hosts - FINAL/leishmania_leishmania/output/plots/known_host_richness.tif')
leish_known[leish_known$known_host_richness == 0, ] <- NA


d <- gplot(viannia_predicted) + 
  geom_tile(aes(fill = value))+
  ylim(-60, 50) +
  xlim(-130, -10) +
  ggtitle("D. viannia predicted") +
  scale_fill_viridis(na.value="white") + #, limits = c(1, 24)) +
  theme_void() +
  theme(#legend.key.height = unit(0.5, 'cm'),
        #legend.key.width = unit(0.5, 'cm'),
        legend.position = c(0.25, 0.40),
        legend.title = element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        plot.margin=unit(c(0.5,-2,0.5,0),"mm"))

c <- gplot(viannia_known) + 
  geom_tile(aes(fill = value)) +
  ylim(-60, 50) +
  xlim(-130, -10) + 
  ggtitle("C. viannia known") +
  scale_fill_viridis(na.value="white") + #, limits = c(1, 46)) +
  theme_void() +
  theme(#legend.key.height = unit(0.5, 'cm'),
    #legend.key.width = unit(0.5, 'cm'),
    legend.position = c(0.25, 0.40),
    legend.title = element_blank(),
    axis.text.x=element_blank(),
    axis.text.y=element_blank(),
    plot.margin=unit(c(0.5,-2,0.5,0),"mm"))

a <- gplot(leish_known) + 
  geom_tile(aes(fill = value)) +
  ylim(-60, 50) +
  xlim(-130, -10) + 
  ggtitle("A. Leishmania known") +
  scale_fill_viridis(na.value="white") + #, limits = c(1, 46)) +
  theme_void() +
  theme(#legend.key.height = unit(0.5, 'cm'),
    #legend.key.width = unit(0.5, 'cm'),
    legend.position = c(0.25, 0.40),
    legend.title = element_blank(),
    axis.text.x=element_blank(),
    axis.text.y=element_blank(),
    plot.margin=unit(c(0.5,-2,0.5,0),"mm"))

b <- gplot(leish_predicted) + 
  geom_tile(aes(fill = value)) +
  scale_fill_viridis(na.value="white") + #, limits = c(1, 24)) +
  ylim(-60, 50) +
  xlim(-130, -10) + 
  ggtitle("B. Leishmania predicted") +
  theme_void() +
  theme(#legend.key.height = unit(0.5, 'cm'),
    #legend.key.width = unit(0.5, 'cm'),
    legend.position = c(0.25, 0.40),
    legend.title = element_blank(),
    axis.text.x=element_blank(),
    axis.text.y=element_blank(),
    plot.margin=unit(c(0.5,-2,0.5,0),"mm"))


gridExtra::grid.arrange(a, b, c, d, ncol = 2)
ggsave('/Users/carolineglidden/Desktop/reservoir hosts - FINAL/joint_maps_diff_scales.png', gridExtra::arrangeGrob(a, b, c, d, ncol = 2))





