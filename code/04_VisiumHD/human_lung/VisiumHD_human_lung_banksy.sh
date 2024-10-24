#!/bin/bash
#SBATCH --job-name=Banksy_Human_lung    # Job name
#SBATCH --output=./logs/Banksy_Human_lung.out  # Output file
#SBATCH --error=./logs/Banksy_Human_lung.err   # Error file
#SBATCH --ntasks=1                                # Run on a single CPU
#SBATCH --mem=128G                                


# Load R module (if necessary, adjust this to match your system)
module load conda_R/4.4

# Run the R script
Rscript VisiumHD_human_lung_banksy.R
