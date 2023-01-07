#!/bin/bash
#SBATCH --job-name=range_study_effort
#SBATCH --time=96:00:00
#SBATCH --qos long
#SBATCH -p normal
#SBATCH -c 16
#SBATCH --mem=128GB

ml R/4.0

Rscript species_range_xgboost_study_effort.R
