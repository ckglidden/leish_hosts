# leish_hosts
> Code for: *Phylogenetic and biogeographical traits predict unrecognized hosts of zoonotic leishmaniasis* 2023. Caroline K Glidden, Aisling R. Murran, Rafaella Albuquerque Silva, Adrian A. Castellanos,  Barbara A. Han, Erin A. Mordecai. *PLoS NTD*: in Reiviion. (preprint : https://www.biorxiv.org/content/10.1101/2022.10.11.511693v2)
> Code written by Caroline K. Glidden & Aisling R. Murran

This github contains the following R files (for each *Leishmania* genus:

1. timetree_phylo_distance.R: code to get phylogenetic distance PCoA axes
2. pubmed_data_collection.R: code to build variable importance plots and tables for each model in the manuscript
3. clean_merge_data_4.1.R : code to clean and merge host status data and trait data included in the raw_data folder
4. species_range_xgboost_model_performance.R : codel to evaluate model performance using nested cross-validation
5. range_variableImp_pdp.R: code to bootstrap variable importance and partial dependence plot analysis
6. figures.R : code to generate figures for final model
 *The folder study_effort_supp includes the same analyses but for the supplemental study effort analysis*

The final data included in the analysis (cleaned_data/analysis_data/trait_data_small.rds) includes data from the following sources: 

* Host status data: data sources can be found in raw_data/reservoir_host_status_using_rank_final.csv
* Species life-history traits: Pantheria https://doi.org/10.1890/08-1494.1
* Foraging data: EltonTraits https://doi.org/10.1890/13-1917.1
* Climate data: CHELSA https://chelsa-climate.org/
* Main habitat (categorical variables): IUCN https://www.iucnredlist.org/
* Phylogenetic distances: TimeTree http://timetree.org/
* Landcover data: Copernicus Global Land Cover Layers: CGLS-LC100 Collection 3 https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_Landcover_100m_Proba-V-C3_Global
* Global human modification index: CSP gHM: Global Human Modification https://developers.google.com/earth-engine/datasets/catalog/CSP_HM_GlobalHumanModification
* Zoonotic host status: GIDEON https://www.gideononline.com/
 
 Species range shapefiles for cleaning data, building figures, and aggregating habitat traits for species ranges can be found here: https://www.iucnredlist.org/
 
 \n 
 
 Occurrence points to generate polygons for invasive species in the Americas were downloaded from da Rosa et al. 2020:  https://doi.org/10.1002/ecy.3115

\n 

Occurrence points to generate polygons for zoonotic *Leishmania* occurrence were downloaded from Herrera et al. 2020: https://doi.org/10.1038/s41597-020-0451-5
