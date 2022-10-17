##### Cleaning and merging reservoir host trait data
##### title: Cleaning reservoir host data
##### author: Aisling Murran, Caroline Glidden
##### version: 06/01/2021


######### READ_ME: 
##### DESCRIPTION: 
  ### This program accesses reservoir host trait data from several database files, cleans 
  ### this data, and combines it into one complete data file. This large set of data is stored as 
  ### a dataframe named total_data.

##### SECTIONS:

### Formatting:
  # Formatting Pantheria Data
  # Formatting Neotropic and Invasive Species Data
  # Formatting Leish Reservoir and Infection Status Data (infected by any Leish: 0/1)
  # Formatting Biogeography Data and Forest Integrity Data
  # Formatting MammalDiet Data
  # Formatting Latitude and Longitude Data
  # Formatting Pubmed Data (details of data collection in pubmed_data_collection.R)
  # Formatting Extra Land-use Traits Data

### Merging to total_data:
  # Merging Pantheria Data
  # Merging Neotropic and Invasive Species Data
  # Merging Leish Reservoir and Infection Status Data (infected by any Leish: 0/1)
  # Merging Biogeography data and Forest Integrity Data
  # Merging MammalDiet Data
  # Merging Latitude and Longitude Data
  # Merging Pubmed Data
  # Merging Extra Land-use Traits Data

### Reformatting total_data to create total_trait_data:
  # Round 1 - removing non-life history traits, repeated traits, or generally irrelevant traits
  # Round 2 - converting categorical variables to binomial variables
  # Round 3 - removing extinct species
  # Round 4 - removing colinearly related traits


#######################
####################### SET UP
#######################

library(caret) #for finding correlations
library(dplyr) #for data clean up
library(tidyverse)
library(purrr) #for merging multiple dataframes at once
library(sf) #this is for reading in shapefiles
library(raster) #for extent function
library(geosphere)
library(foreach)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #set working directory to current directory

###Important variables
#total_data is the compete collection of all our database information saved as a dataframe organized by species. (Includes everything)
total_data <- data.frame()

#total_trait_data is the compete collection of all our traits saved as a dataframe organized by species. (Includes correlated traits, but not repeated traits or inapplicable traits)
total_trait_data <- data.frame()

#trait_data is the refined collection of all our traits saved as a dataframe organized by species. (Traits to train the model on)


#######################
####################### FORMATTING DATA
#######################


####################### Formatting Pantheria Data
#######################

###Import all the pantheria data
Pantheria_1_df <- read.delim("../raw data/Pantheria & Species Data/PanTHERIA_1-0_WR05_Aug2008.txt") #read in txt 1 file separating by tabs
Pantheria_1_df$taxonomic.year <- 'MSW05' #keep track of what binomial we are using
Pantheria_2_df <- read.delim("../raw data/Pantheria & Species Data/PanTHERIA_1-0_WR93_Aug2008.txt") #read in txt 2 file separating by tabs
Pantheria_2_df$taxonomic.year <- 'MSW93' #keep track of what binomial we are using

#read in the species that weren't in pantheria until the update
#cat("\n", file = file.choose(), append = TRUE) #if you get the annoying warning
Pantheria_3_df <- read.delim('../raw data/Pantheria & Species Data/pantheria_updated.txt')
Pantheria_3_df$taxonomic.year <- NA #to keep cols consistent with the other Pantheria data


###Combine the pantheria data into one dataframe

#1: Remove species already included in pantheria WR05 from pantheria WR93
duplicate.names <- Pantheria_1_df$MSW05_Binomial
Pantheria_2_df_reduced <- Pantheria_2_df[!(Pantheria_2_df$MSW93_Binomial %in% duplicate.names), ]

#2: Coercing the column names to be the same in order to merge. (The names of the cols are not identical because diff taxonomic book editions were used, but they do contain the same info. For rbind they must be the same so we are changing them.)
names(Pantheria_2_df_reduced) <- names(Pantheria_1_df) 
names(Pantheria_3_df) <- names(Pantheria_1_df)

#3: Merge the three pantheria data sets
Pantheria_df <- rbind(Pantheria_1_df, Pantheria_2_df_reduced, Pantheria_3_df) #Combining the dataframes
Pantheria_df[Pantheria_df == -999] <- NA #replace -999 w/NA (-999 was used to signify that data wasn't collected)


####################### Formatting Neotropical and Invasive Species Data
#######################

###Import the neotropical and invasive (including domestic animals) species list (also includes zoonotic host status & endemic vs invasive)
endemic_species <- read.csv("../raw data/Pantheria & Species Data/NeotropicalMammalHostStatus_cleaned.csv"); endemic_species <- endemic_species[,-8]
invasive_species <- read.csv("../raw data/Pantheria & Species Data/list of domestic and invasive species GISD.csv")


###Combine the species data into one species dataframe
species <- rbind(endemic_species, invasive_species)


###Renaming columns to match total_data (pantheria)
names(species) <- c("MSW05_Order","MSW05_Family","MSW05_Genus","MSW05_Species","MSW05_Binomial","zoonotic_host_status", "native") #change names to match pantheria

###Reformatting taxonomy to match total_data (pantheria)
fix.tax <- read.csv("../raw data/Pantheria & Species Data/fixing taxonomy.csv")
ntl <- fix.tax$MSW05_Binomial.ntl
pantheria <- fix.tax$MSW05_Binomial.pan
for(i in 1:length(ntl)){
  species$MSW05_Binomial[species$MSW05_Binomial==ntl[i]] <- pantheria[i]
}

#compare species names between data-sets, identify which ones don't match up
#list <- setdiff(species$MSW05_Binomial, Pantheria_df$MSW05_Binomial) #282 species that are in IUCN list but not in Pantheria
#write.csv(list, '../raw data/Extra Data/list of species not in pantheria.csv') #saved in folder

####################### Formatting Leish Reservoir Data to convert to "infected by any leish"
#######################

###Import leish reservoir data
reservoir.data <- read.csv("../raw data/Pantheria & Species Data/reservoir host status_using rank_final.csv")
reservoir.data <- reservoir.data[,c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29)]
reservoir.data[is.na(reservoir.data)] <- 0 #change NAs to 0s

reservoir.data$amazonensis.host <- ifelse(reservoir.data$amazonensis.reservoir > 0, 1,0)
#reservoir.data$braziliensis.host <- ifelse(reservoir.data$braziliensis.reservoir > 0, 1,0)
#reservoir.data$enriettii.host <- ifelse(reservoir.data$enriettii.reservoir > 0, 1,0)
#reservoir.data$guyanensis.host <- ifelse(reservoir.data$guyanensis.reservoir > 0, 1,0)
#reservoir.data$lainsoni.host <- ifelse(reservoir.data$lainsoni.reservoir > 0, 1,0)
#reservoir.data$lindenbergi.host <- ifelse(reservoir.data$lindenbergi.reservoir > 0, 1,0)
reservoir.data$mexicana.host <- ifelse(reservoir.data$mexicana.reservoir > 0, 1,0)
#reservoir.data$naiffi.host <- ifelse(reservoir.data$naiffi.reservoir > 0, 1,0)
#reservoir.data$panamensis.host <- ifelse(reservoir.data$panamensis.reservoir > 0, 1,0)
#reservoir.data$peruviana.host <- ifelse(reservoir.data$peruviana.reservoir > 0, 1,0)
#reservoir.data$shawi.host <- ifelse(reservoir.data$shawi.reservoir > 0, 1,0)
reservoir.data$venezuelensis.host <- ifelse(reservoir.data$venezuelensis.reservoir > 0, 1,0)
#reservoir.data$colombiensis.host <- ifelse(reservoir.data$colombiensis.reservoir > 0, 1,0)
reservoir.data$infantum.host <- ifelse(reservoir.data$infantum.reservoir > 0, 1,0)

infection.data <- reservoir.data[,c(1,16:19)]
infection.data$sum <- rowSums(infection.data[,2:5])
infection.data$leish.infection <- ifelse(infection.data$sum > 0, 1, 0)

###reduce just to infectioun status ("has species every been infected with Leish")
leish.status <- infection.data[,c(1,7)] #my rewritten code added extra rows, prob bc of the NA replace part so I just delete them here
names(leish.status)[1] <- 'MSW05_Binomial'

###Reformatting taxonomy to match total_data 
for(i in 1:length(ntl)){
  leish.status$MSW05_Binomial[leish.status$MSW05_Binomial==ntl[i]] <- pantheria[i]
}


####################### Formatting Biogeography Data
#######################

chelsa_climate <- read.csv("../raw data/ACL_GEE_FINAL/leishmania_hosts_climate.csv")

chelsa_climate_agg <- aggregate(cbind(bio1_mean_annual_temp, bio4_temp_seasonality,
                                      bio10_temp_warmest_qt, 
                                      bio11_temp_coldest_qt, bio12_annual_precip, bio15_precip_seasonality,
                                      bio16_precip_wettest_qt, bio17_precip_driest_qt, 
                                      cmi_range, cmi_mean,
                                      hurs_range, hurs_mean) ~ BINOMIAL, data = chelsa_climate, FUN = mean, na.rm = TRUE)
names(chelsa_climate_agg)[1] <- "MSW05_Binomial"

human.popsize <- read.csv("../raw data/ACL_GEE_FINAL/mean_pop_1km2_2020.csv")
pop.agg <- aggregate(by=list(human.popsize$BINOMIAL), x=human.popsize$mean, FUN=mean) #take mean for each species
names(pop.agg) <- c("MSW05_Binomial",'mean_humanpop_km2')

#elevation <- read.csv("../raw data/Biogeography Data/mean_range_elevation.csv")
#elev.agg <- aggregate(by=list(elevation$BINOMIAL), x=elevation$mean, FUN=mean) #take mean for each species
#names(elev.agg) <- c("MSW05_Binomial",'mean_elevation_m')

flii <- read.csv("../raw data/ACL_GEE_FINAL/mean_range_integrity_updated.csv")
flii.agg <- aggregate(by=list(flii$BINOMIAL), x=flii$mean, FUN=mean) 
names(flii.agg) <- c("MSW05_Binomial",'flii')

#habitat_breadth <- read.csv("../raw data/Biogeography Data/neotropical species habitat breadth.csv")
#habitat_breadth <- habitat_breadth[,c(2:4)] #only include columns we need
#habitat_breadth.agg <- aggregate(.~MSW05_Binomial, habitat_breadth, mean) #reduce to one line per species with mean

tree_cover <- read.csv("../raw data/ACL_GEE_FINAL/mean_tree_cover.csv")
tree_cover.agg <- aggregate(by=list(tree_cover$BINOMIAL), x=tree_cover$mean, FUN=mean) #reduce to one line per species with mean
names(tree_cover.agg) <- c("MSW05_Binomial",'mean_tree_cover')

grass_cover <- read.csv("../raw data/ACL_GEE_FINAL/mean_grass_cover.csv")
grass_cover.agg <- aggregate(by=list(grass_cover$BINOMIAL), x=grass_cover$mean, FUN=mean) #reduce to one line per species with mean
names(grass_cover.agg) <- c("MSW05_Binomial",'mean_grass_cover')

crop_cover <- read.csv("../raw data/ACL_GEE_FINAL/mean_range_crop_cover.csv")
crop_cover.agg <- aggregate(by=list(crop_cover$BINOMIAL), x=crop_cover$mean, FUN=mean) #take mean for each species
names(crop_cover.agg) <- c("MSW05_Binomial",'mean_crop_cover')

urban_cover <- read.csv("../raw data/ACL_GEE_FINAL/mean_range_urban_cover.csv")
urban_cover.agg <- aggregate(by=list(urban_cover$BINOMIAL), x=urban_cover$mean, FUN=mean) #take mean for each species
names(urban_cover.agg) <- c("MSW05_Binomial",'mean_urban_cover')

canopyHeight <- read.csv("../raw data/ACL_GEE_FINAL/mean_canopy_cover.csv")
canopyHeight.agg <- aggregate(by=list(canopyHeight$BINOMIAL), x=canopyHeight$mean, FUN=mean) #take mean for each species
names(canopyHeight.agg) <- c("MSW05_Binomial",'mean_canopy_height')

gHM <- read.csv("../raw data/ACL_GEE_FINAL/mean_gHM.csv") 
gHM.agg <- aggregate(by=list(gHM$BINOMIAL), x=gHM$mean, FUN=mean) #take mean for each species
names(gHM.agg) <- c("MSW05_Binomial",'mean_gHM')

habitat <- read.csv('../raw data/aoh habitat data/iucn_habitat_classifications_april42022.csv'); names(habitat)[2] <- "MSW05_Binomial"; habitat <- habitat[,-1]
habitat$habitat_breadth <- rowSums(habitat[, 2:ncol(habitat)])

biogeography <- list(chelsa_climate_agg, flii.agg, pop.agg, tree_cover.agg, grass_cover.agg, crop_cover.agg, urban_cover.agg, gHM.agg, canopyHeight.agg, habitat) %>% reduce(left_join, by = "MSW05_Binomial")


###Reformatting taxonomy to match total_data 
for(i in 1:length(ntl)){
  biogeography$MSW05_Binomial[biogeography$MSW05_Binomial==ntl[i]] <- pantheria[i] #make sure names match pantheria
}


####################### Formatting MammalDiet Data
#######################

#import mammalDiet data
MammalDiet_df <- read.delim("../raw data/MammalDiet Data/MammalDIET_v1.0 2.txt") #read in txt 1 file separating by tabs

MammalDiet_df$TaxonID <- c(paste(MammalDiet_df$Genus, MammalDiet_df$Species))
names(MammalDiet_df)[1] <- 'MSW05_Binomial'

###Reformatting taxonomy to match total_data 
for(i in 1:length(ntl)){
  MammalDiet_df$MSW05_Binomial[MammalDiet_df$MSW05_Binomial==ntl[i]] <- pantheria[i] #make sure names match pantheria
}


####################### Formatting Latitude and Longitude Data & Getting area of range
#######################

#import latitude and longitude data
#Lat_Long_shape_file <- st_read("../raw data/Neotropical Species Range Data/data_0.shx") #read in shape file
#range_file <- as_Spatial(Lat_Long_shape_file)
#area_sqkm <- areaPolygon(range_file)/1000
#range_size <- cbind(Lat_Long_shape_file$BINOMIAL, area_sqkm)

#Converting from a matrix to a dataframe 
#range_size <- data.frame(range_size)
#names(range_size)[1] <- 'MSW05_Binomial'

#reduce to just one range per species by aggregating by mean
#range_size$area_sqkm <- as.numeric(range_size$area_sqkm)
#range_size.agg <- aggregate(.~MSW05_Binomial, range_size, mean)

###Reformatting taxonomy to match total_data 
#for(i in 1:length(ntl)){
#  range_size.agg$MSW05_Binomial[range_size.agg$MSW05_Binomial==ntl[i]] <- pantheria[i] #make sure names match pantheria
#}

#extent_lapply <- lapply(st_geometry(Lat_Long_shape_file), st_bbox)
#species_num <- length(extent_lapply)
#Lat_Long_matrix <- matrix(ncol = 5)

#Converting from a shape file to a matrix
#for (i in 1:species_num) #for each species in file
#{
#  name <- Lat_Long_shape_file [i,]$BINOMIAL
#  xmin <- extent_lapply[[i]][1]
#  ymin <- extent_lapply[[i]][2]
#  xmax <- extent_lapply[[i]][3]
#  ymax <- extent_lapply[[i]][4]
#  individual_vector <- c(name, xmin, xmax, ymin, ymax)
#  Lat_Long_matrix <- rbind(Lat_Long_matrix, individual_vector)
#}

#Lat_Long_matrix <- Lat_Long_matrix[2:species_num,]

#Converting from a matrix to a dataframe - this messes up repeated species
#Lat_Long_df <- data.frame(Lat_Long_matrix)
#names(Lat_Long_df)[1] <- 'MSW05_Binomial'

#reduce to just one range per species by aggregating by mean
#Lat_Long_df <- Lat_Long_df %>% mutate_at(c(2:5), as.numeric)
#Lat_Long_df.agg <- aggregate(.~MSW05_Binomial, Lat_Long_df, mean)

###Reformatting taxonomy to match total_data 
#for(i in 1:length(ntl)){
#  Lat_Long_df$MSW05_Binomial[Lat_Long_df$MSW05_Binomial==ntl[i]] <- pantheria[i] #make sure names match pantheria
#}

#names(Lat_Long_df) <- c("MSW05_Binomial", 'xmin','ymin','xmax','ymax')

####################### Formatting Elton Trait Data
#######################

#Convert to dataframe
EltonTraits_df <- read.delim("../raw data/Elton Traits Data/MamFuncDat.txt") #read in shape file

names(EltonTraits_df)[2] <- 'MSW05_Binomial'

###Reformatting taxonomy to match total_data 
for(i in 1:length(ntl)){
  EltonTraits_df$MSW05_Binomial[EltonTraits_df$MSW05_Binomial==ntl[i]] <- pantheria[i] #make sure names match pantheria
}


####################### Formatting Pubmed Trait
#######################

pubMed <- read.csv("../raw data/Extra Data/neotropical species pubmed citations.csv")
#####ADD BINOMIAL


#######################Formatting phylogeny data
phylogenetic <- read.csv('../cleaned data/extra/pcoa_phylogenetic_distances.csv'); phylogenetic <- phylogenetic[,-1] #read in pcoa results
phylogenetic <- phylogenetic[,1:7]

####Identifying wild versus domestic hosts
domestic <- read.csv("../cleaned data/extra/MSW_order_family_domestic_status.csv")

#######################
####################### MERGING DATA
#######################

##### Pantheria Data Merging
total_data <- Pantheria_df


##### Neotropical and Invasive Species Merging - make sure to keep all species listed in the species data file
total_data <- left_join(species[5:7],total_data, by='MSW05_Binomial') #only use binomial from species list, high tax groups all from pantheria
colnames(total_data)[8:57]<-sub("^[^_]*_","",colnames(total_data)[8:57]) #remove unnecessary characters in column names


##### Leish Reservoir and Infection Data Merging
total_data <- left_join(total_data, leish.status, by = 'MSW05_Binomial', all=TRUE)
total_data$leish.infection[is.na(total_data$leish.infection)] <- 0

##### Biogeography Data Merging
total_data <- left_join(total_data, biogeography, by="MSW05_Binomial")


##### MammalDiet Data Merging
total_data <- left_join(total_data, MammalDiet_df, by="MSW05_Binomial")

##### updated area
#total_data <- left_join(total_data, range_size.agg, by="MSW05_Binomial")

##### Latitude and Longitude Data Merging
#total_data <- left_join(total_data, Lat_Long_df, by="MSW05_Binomial")


##### Elton Traits Data Merging
total_data <- left_join(total_data, EltonTraits_df, by="MSW05_Binomial")

##### Extra Land-Use Data Merging
#total_data <- left_join(total_data, biogeography_2, by="MSW05_Binomial")

##### PubMed Merging
total_data <- left_join(total_data, pubMed, by="MSW05_Binomial")

##### Phylogeny Merging
total_data <- left_join(total_data, phylogenetic, by = "MSW05_Binomial")

##### Domestic Merging
total_data <- left_join(total_data, domestic[,c(1,4)], by = "MSW05_Binomial")


####Formating/Saving total_data after merging
total_data <- unique(total_data) #remove duplicate rows, there are still some species repeated but we can remove those in the next step (cutting down traits)

saveRDS(total_data, '../cleaned data/analysis data/total data.rds')

#######################
#######################  REFORMATTING TOTAL_DATA TO CREATE TOTAL_TRAIT_DATA
#######################

####################### Round One: creating total_trait_data without unhelpful traits
####################### 

##### Pantheria Data Cleaning

#Create combined Pantheria traits
total_data$RelAgetoSexualMaturity <- total_data$SexualMaturityAge_d/(total_data$MaxLongevity_m*30.42) #30.42 = average # of days in a month in a non-leap year
total_data$RelAgetoDisperal <- total_data$DispersalAge_d/(total_data$MaxLongevity_m*30.42)
total_data$RelAgetoWeaning <- total_data$WeaningAge_d/(total_data$MaxLongevity_m*30.42)
total_data$metabloic_standardized <- total_data$BasalMetRate_mLO2hr/total_data$BasalMetRateMass_g

#Designate what Pantheria traits are not helpful and will be removed using the trait description data table on GDrive
pantheria_unhelpful <- c('MSW05_Genus',	'MSW05_Species',	
                         'HabitatBreadth',	'HomeRange_Indiv_km2',	'InterbirthInterval_d',	
                         'Terrestriality',	'TrophicLevel.x',	'References',	
                         'AdultBodyMass_g_EXT',	'LittersPerYear_EXT',	'NeonateBodyMass_g_EXT',	'WeaningBodyMass_g_EXT',	
                         'HuPopDen_Min_n.km2',	'HuPopDen_Mean_n.km2',	'HuPopDen_5p_n.km2',	'HuPopDen_Change',
                         'Precip_Mean_mm',	'Temp_Mean_01degC',	'AET_Mean_mm',	'PET_Mean_mm',	
                         'taxonomic.year', "BasalMetRate_mLO2hr", "BasalMetRateMass_g")


#### EltonTrait Data Cleaning

#Designate what EltonTrait data are not helpful and will be removed (traits to keep on description table on GDrive)
elton_unhelpful <- c('MSW3_ID',	'MSWFamilyLatin',	
                     'Diet.Inv',	'Diet.Vend',	'Diet.Vect',	'Diet.Vfish',	'Diet.Vunk',	'Diet.Scav',	'Diet.Fruit',	'Diet.Nect',
                     'Diet.Seed',	'Diet.PlantO',	'Diet.Source',	'Diet.Certainty',	
                     'ForStrat.Certainty',	'ForStrat.Comment',	"Activity.Nocturnal", "Activity.Crepuscular", "Activity.Diurnal", 'Activity.Source',	'Activity.Certainty',	
                     'BodyMass.Value',	'BodyMass.Source',	'BodyMass.SpecLevel')


#### MammalDiet Data Cleaning

#Designate what MammalDiet traits are not helpful and will be removed
mammalDiet_unhelpful <- c('Order',	'Family',	'Genus',	'Species',	'Animal',	'Vertebrate',	'Mammal',	'Bird',	
                         'Herptile',	'Fish',	'Invertebrate',	'Plant',	'Seed',	'Fruit',	'Nectar',	'Root',	'Leaf',	'Woody',	
                         'Herbaceous',	'Other',	'TaxonomicNote',	'FillCode',	'DataSource')

#### Other Data Cleaning
#Designate what other traits are not helpful and will be removed
extra_unhelpful <- c('X', "X.x", "X.x.x","X.y","X.y.y", "habitatB.biomes")

unhelpful_traits <- append(pantheria_unhelpful, elton_unhelpful)
unhelpful_traits <- append(unhelpful_traits, mammalDiet_unhelpful)
unhelpful_traits <- append(unhelpful_traits, extra_unhelpful)
  
total_trait_data <- total_data[,!(names(total_data) %in% unhelpful_traits)] #we should have 88 traits (variables) left!
names(total_trait_data) #manually double check that traits are correct  
total_trait_data <- unique(total_trait_data) #try to get rid of multiple rows for a species again
total_trait_data <- total_trait_data %>% distinct(MSW05_Binomial, .keep_all = TRUE) #removes final rows for multiple species, keeps first row of duplicated species


#MAKE SURE TO RUN CODE BELOW!!!
####################### Round Two: changing categorical traits into binomial traits (each category becomes its own trait)
####################### 

###First we will import a list of all the trait names and their data type
trait_types <- read.csv("../raw data/Extra Data/trait_types.csv")
categorical_traits <- trait_types %>% filter(trait_types$type == 'categorical')
#binomial_traits <- trait_types %>% filter(trait_types$type == 'binary')
#numeric_traits <- trait_types %>% filter(trait_types$type == 'numeric')

#####one hot encode variables
categorical <- total_trait_data[categorical_traits$Trait_name]
categorical <- categorical %>% mutate_if(is.integer,as.factor)
dmy <- dummyVars(" ~ .", data = categorical)
dmy.rewrite <- data.frame(predict(dmy, newdata = categorical))
dmy.rewrite$MSW05_Binomial <- c(total_trait_data$MSW05_Binomial)
rewritten_categorical <- dmy.rewrite

#remove the categorical variables. Replace them with the new ones
total_trait_data <- total_trait_data[, !(names(total_trait_data) %in% categorical_traits$Trait_name)]
total_trait_data <- left_join(total_trait_data, rewritten_categorical, by = 'MSW05_Binomial')

#remove marine orders and foraging strata, trophic level not assigned
total_trait_data <- total_trait_data[,-c(82, 159, 164)]

###We also reordered our traits to make dealing with data type easier in the future -> categorical, binomial, and then numeric
#First we have to remove MSW05_Binomial from the names in rewritten_categorical because we want this trait first
#drop <- c('MSW05_Binomial')
#rewritten_categorical = rewritten_categorical[,!(names(rewritten_categorical) %in% drop)]

#Now we can reorder by col names
#col_order <- c('MSW05_Binomial', binomial_traits$Trait_name, numeric_traits$Trait_name, colnames(rewritten_categorical))
#total_trait_data <- total_trait_data[, col_order]


####################### Round Three: removing extinct species
####################### 

###First we will import a list of all the extinct mammals from IUCN

extinct_species <- read.delim("../raw data/Extra Data/list of extinct hosts.txt")

#Remove the rows of the species in the extinct_species dataframe
total_trait_data <- total_trait_data[!(total_trait_data$MSW05_Binomial %in% extinct_species$names),]
total_trait_data <- subset(total_trait_data, MSW05_Binomial != 'Gracilinanus ignitus')
total_trait_data <- subset(total_trait_data, MSW05_Binomial != 'Brotomys voratus')
total_trait_data <- subset(total_trait_data, MSW05_Binomial != 'Inia geoffrensis')

########################
#Remove species with < 10% overlap with Leishmania ranges
positive_hosts <- subset(total_trait_data, leish.infection == 1)
negative_hosts <- subset(total_trait_data, leish.infection == 0)

#read in animal ranges
species_ranges <- st_read("../raw data/leishmania_animals_all_clipped/leishmania_animals_all_clipped.shp")
for(i in 1:length(ntl)){
  species_ranges$BINOMIAL[species_ranges$BINOMIAL==ntl[i]] <- pantheria[i]
}

names(species_ranges)[1] <- "MSW05_Binomial"

#####read in leish ranges
leishmania_ranges <- st_read("../raw data/leishmania_distributions/leishmania_distributions.shp")
leishmania_ranges <- leishmania_ranges[c(1, 7, 10),]
leish_genus_range <- leishmania_ranges %>% summarise()
#st_write(leish_genus_range, "../raw data/leishmania_distributions/leishmania_genus_range.shp")

#percent ovelap per animal

#area overlap
area_data <- foreach(i = unique(negative_hosts$MSW05_Binomial), 
                     .packages = c("sf", "tidyverse", "R.utils","units","doParallel"),
                     .combine = rbind,
                     .errorhandling = 'remove') %dopar% {
                       
                       SHP1 <- species_ranges %>%
                         filter(MSW05_Binomial == i) %>%
                         group_by(MSW05_Binomial) %>%
                         summarise()
                       
                       SHP2 <- leish_genus_range
                       
                       AREA1 <- set_units(st_area(SHP1),km^2)
                       AREA2 <- set_units(st_area(SHP2),km^2)
                       
                       OVR <- st_intersection(SHP1, SHP2)
                       
                       # add in areas in km2
                       
                       AREA <- ifelse(nrow(OVR) > 0 , set_units(st_area(OVR),km^2), 0)
                       
                       mammal_overlap <- as.numeric(AREA/AREA1)
                       
                       area_data = data.frame(MSW05_Binomial = unique(SHP1$MSW05_Binomial), 
                                              mammal_area = AREA1,
                                              leish_area = AREA2,
                                              area_km2 = AREA,
                                              percent_overlap = mammal_overlap
                                              )
                       
                     }

###only retain animals w/ > 10% of range in Leish range
no_overlap_hosts <- subset(area_data, percent_overlap < 0.1); id_no_overlap_hosts <- unique(no_overlap_hosts$MSW05_Binomial) 

"%ni%" <- Negate("%in%")
negative_hosts <- negative_hosts[negative_hosts$MSW05_Binomial %ni% id_no_overlap_hosts,]

total_trait_data <- rbind(positive_hosts, negative_hosts) #now down to 1480 hosts b/c removed 205 hosts


####################### SAVING THE FINAL FILE WITH ALL TRAIT DATA!

#This is the data file imputations were preformed on.
saveRDS(total_trait_data, '../cleaned data/analysis data/total trait data.rds')

####START HERE!
####################### Round Four: creating trait_data without low coverage & correlated traits
####################### 
trait_data <- total_trait_data[,-c(93:152)] #remove family

# how many NA?
na.amount = apply(trait_data, 2, FUN = function(x) length(which(is.na(x))))

na.df = data.frame(Variable=names(na.amount),
                   NAlength=as.numeric(na.amount)) #Lists num of NA by trait

na.df = na.df[order(na.df[,2], decreasing=T),] #Reorders the data so that the traits with the most NA are at the front
row.names(na.df) = NULL #resets row names to match order

na.df$percent_coverage <- 1-na.df$NAlength/1480

## take columns with more than 1332 non-NA entries ((1332  = 90% of the data (so 90% can be missing)))
which.in = with(na.df, Variable[which(NAlength < 1332)])
trait_data_small = trait_data[,which(names(trait_data) %in% which.in)]

trait_data_small <- trait_data_small[,c(25, 1:24, 26:91)]

#Look at correlation among variables by creating a correlation matrix
#num_cols <- unlist(lapply(total_trait_data[,c(4:32, 34:44, 60:65,67:70)], is.numeric)) #extract col numbers that are numeric (continuous variables)
cor.matrix <- cor(trait_data_small[,c(5:45, 54, 60:66)], use="pairwise.complete.obs") #make full correlation matrix
cor.matrix[is.na(cor.matrix)] <- 0 # change NA to 0 to apply cutoff
write.csv(cor.matrix, "../output/corelation_matrix.csv")
high.cor <- findCorrelation(cor.matrix, cutoff=0.7, names=TRUE)

#[1] "WeaningHeadBodyLen_mm"  "DispersalAge_d"         "NeonateHeadBodyLen_mm"  "AgeatFirstBirth_d"      "SexualMaturityAge_d"   
#[6] "RelAgetoDisperal"       "AdultHeadBodyLen_mm"    "mean_tree_cover"        "NeonateBodyMass_g"      "SocialGrpSize"         
#[11] "RelAgetoSexualMaturity" "WeaningBodyMass_g"      "GR_MaxLat_dd"           "GR_MidRangeLong_dd"     "GR_MinLat_dd"          
#[16] "mean_humanpop_km2" 

#high.cor.matrix <-cor.matrix[,high.cor] #subset to traits that are highly correlative
#high.cor.matrix[high.cor.matrix < 0.7] <- NA 

#correlations (without double counting)
#adult body mass correlated with: weaning body length, neonate body length, weaning body mass
#adult forearm length correlated with: weaning body length, neonate body length, adult body length, sexual maturity age
#age at first birth correlated with: neonate body length, max longevity, sexual maturity age
#dispersal age correlated with: gestation length
#litter per year correlated with: weaning body length
#litter size correlated with: weaning body length
#max longevity correlated with: age at first birth, sexual maturity age
#neonate body mass correlated with: neonate body length, weaning body mass
#Population group size correlate with: relative age to dispersal
#weaning age correlated with: weaning body legnth, age at first birth, sexual maturity age
#weaning body mass correlated with: weaning body length, neonate body mass
#low coverage traits to exclude: age at first birth, dispersal age, neonate body length, weaning body mass, weaning body length, age to dispersal, relative age to dispersal
# also exclude: adult forearm length, adult body length, sexual maturity age

#Choose variables to minimize correlation - consider excluding one if correlation > 0.8, can base which ones to remove off of low coverage traits
#correlated_traits <- c('AgeatFirstBirth_d', "DispersalAge_d", "NeonateHeadBodyLen_mm", "NeonateBodyMass_g", "WeaningBodyMass_g", "WeaningHeadBodyLen_mm", "DispersalAge_d",
#                       "RelAgetoDisperal", "AdultForearmLen_mm", "AdultHeadBodyLen_mm", "GR_Area_km2", "GR_MaxLat_dd", "GR_MinLat_dd", "GR_MaxLong_dd",
#                       "GR_MinLong_dd", "mean_tree_cover")

#correlated_traits <- names(as.data.frame(high.cor.matrix)); correlated_traits <- correlated_traits[c(1:9, 11:16, 26, 28)] #make sure just numeric traits are thrown out

#then removed correlated traits
high.cor.v1 <- high.cor[-1] #remove tree cover, swap with canopy height
high.cor.v2 <- c(high.cor.v1, "mean_canopy_height", "flii", "mean_grass_cover")

trait_data <- trait_data_small[,!(names(trait_data_small) %in% high.cor.v2)]
saveRDS(trait_data, '../cleaned data/analysis data/trait data.rds')
#32 traits
####################### Round Five: transform variables, remove wild hosts
####################### 

#subset down to wild hosts and remove domestic column & trophic level y
trait_data_small <- subset(trait_data, domestic == 0); trait_data_small <- trait_data_small[,-49]

#remove order b/c prob colinear w/ phylo
trait_data_small <- trait_data_small[,-c(50:62)]

#log transform extreme cases
#trait_data_small$AdultBodyMass_g <- log(trait_data_small$AdultBodyMass_g + 1)
trait_data_small$PopulationDensity_n.km2 <- log(trait_data_small$PopulationDensity_n.km2 + 1)
trait_data_small$GR_Area_km2  <- log(trait_data_small$GR_Area_km2 + 1)
trait_data_small$mean_urban_cover  <- log(trait_data_small$mean_urban_cover + 0.01)
trait_data_small$pubmed.count  <- 1/ (trait_data_small$pubmed.count + 1) # more extreme transformation

saveRDS(trait_data_small, "../cleaned data/analysis data/trait_data_small.rds")


#####################
###now save pubmed data

pubmed_data <- trait_data_small[,c(42,3:41, 43:59)]
saveRDS(pubmed_data, "../cleaned data/analysis data/pubmed_analysis_data.rds")

table(trait_data_small$leish.infection)
#0 = 1368, 1 = 102










################################OLD CODE - ANALYSES FROM USING RANGES INSTEAD OF AOH
####################### Formatting Extra Land-use Traits
#######################

#import Extra Land-Use data

#canopyHeight <- read.csv("../raw data/Extra Land-use Data/mean_canopy_height_v2.csv")
#canopyHeight.agg <- aggregate(by=list(canopyHeight$BINOMIAL), x=canopyHeight$mean, FUN=mean) #take mean for each species
#names(canopyHeight.agg) <- c("MSW05_Binomial",'mean_canopy_height')

#gHM <- read.csv("../raw data/Extra Land-use Data/mean_gHM.csv") 
#gHM.agg <- aggregate(by=list(gHM$BINOMIAL), x=gHM$mean, FUN=mean) #take mean for each species
#names(gHM.agg) <- c("MSW05_Binomial",'mean_gHM')

#ndvi <- read.csv("../raw data/Extra Land-use Data/mean_ndvi_2009to2019.csv") #read in txt 1 file separating by tabs
#ndvi.agg <- aggregate(by=list(ndvi$BINOMIAL), x=ndvi$mean, FUN=mean) #take mean for each species
#names(ndvi.agg) <- c("MSW05_Binomial",'mean_ndvi_2009to2019')

#crop_cover <- read.csv("../raw data/Extra Land-use Data/mean_range_crop_cover.csv") #read in txt 1 file separating by tabs
#crop_cover.agg <- aggregate(by=list(crop_cover$BINOMIAL), x=crop_cover$mean, FUN=mean) #take mean for each species
#names(crop_cover.agg) <- c("MSW05_Binomial",'mean_range_crop_cover')

#urban_cover <- read.csv("../raw data/Extra Land-use Data/mean_range_urban_cover.csv") #read in txt 1 file separating by tabs
#urban_cover.agg <- aggregate(by=list(urban_cover$BINOMIAL), x=urban_cover$mean, FUN=mean) #take mean for each species
#names(urban_cover.agg) <- c("MSW05_Binomial",'mean_range_urban_cover')

#biogeography_2 <- list(gHM.agg, ndvi.agg, crop_cover.agg, urban_cover.agg, canopyHeight.agg) %>% reduce(left_join, by = "MSW05_Binomial")

###Reformatting taxonomy to match total_data 
#for(i in 1:length(ntl)){
#  biogeography_2$MSW05_Binomial[biogeography_2$MSW05_Binomial==ntl[i]] <- pantheria[i] #make sure names match pantheria
#}