##### Cleaning and merging reservoir host trait data
##### title: Cleaning reservoir host data
##### author: Caroline Glidden
##### version: 08/30/2021


######### READ_ME: 
##### DESCRIPTION: 
### This program reads in a newick file from TimeTree, converts to phylogenetic distances, 
###and reduces dimensionality in a PCoA

library(ape)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #set working directory to current directory

tree <- read.tree('../raw data/phylogenetic distance time tree/time tree species list.nwk') #from TimeTree website
#tree2 <- read.nexus('../raw data/phylogenetic distance time tree/MamPhy_fullPosterior_BDvr_DNAonly_4098sp_topoFree_NDexp_MCC_v2_target.tre') #from Upham et al. 2019
#tree3 <- read.nexus('../raw data/phylogenetic distance time tree/Bininda-emonds 2007 mammals.nex') #from Bininda-emonds et al. 2007

dist.mat <- cophenetic.phylo(tree) #get branch length between tips

isSymmetric(dist.mat) #check to make sure matrix is symmetric
any(is.complex(dist.mat)) #check to make sure there are no complex numbers

pcoa <- cmdscale (dist.mat, k=10, eig = TRUE) #run pcoa, save ten axes
round(pcoa$eig*100/sum(pcoa$eig),1) #get % of variance explained by each axes to get number of axes to include; first 5 axes explain 78.5% of variation
barplot(pcoa$eig, names = paste ('PCoA', 1:21), las = 3, ylab = 'eigenvalues')

newdata <- as.data.frame(pcoa$points) #make into new dataframe
species <- row.names(newdata) #make sure species names are attached
newdata$species <- species #make sure species names are attached
newdata$species <- gsub("_", " ", newdata$species) #remove underscore and replace with space for species names

vegan::ordiplot(pcoa$points, choices=c(1,3)) #data viz

#read in binomial names to make dataframe easy to merge with total trait data
#species_names <- read.csv('../raw data/phylogenetic distance time tree/time tree species list.csv', header = FALSE)
#names(species_names) <- 'MSW05_Binomial'
#species_names$species <- species_names$MSW05_Binomial

#test <- merge(species_names, newdata, by='species', all=TRUE) #see which names don't merge
#write.csv(test,'../raw data/phylogenetic distance time tree/update timetree species names.csv')

#upload data frame to fix names
new_names <- read.csv('../raw data/phylogenetic distance time tree/update timetree species names.csv')

newdataMSW <- merge(new_names, newdata, by='species') #merge so that pantheria nomenclature is maintained

write.csv(newdataMSW, '../cleaned data/extra/pcoa_phylogenetic_distances.csv')
