#!/bin/bash
#SBATCH --job-name=range_reservoirs_hosts
#SBATCH --time=96:00:00
#SBATCH --qos long
#SBATCH -p normal
#SBATCH -c 16
#SBATCH --mem=128GB

ml R/4.0

Rscript range_variableImp_pdp.R
