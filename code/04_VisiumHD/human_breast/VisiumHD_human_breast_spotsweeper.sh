#!/bin/bash
#SBATCH --job-name=SW_Human_Breast    # Job name
#SBATCH --output=./logs/SW_Human_Breast.out  # Output file
#SBATCH --error=./logs/SW_Human_Breast.err   # Error file
#SBATCH --ntasks=1                                # Run on a single CPU
#SBATCH --mem=128G                                


# Load R module (if necessary, adjust this to match your system)
module load conda_R/4.4

# Run the R script
Rscript VisiumHD_human_breast_spotsweeper.R
