#!/bin/bash
#SBATCH --job-name=BreastCancer_sw    # Job name
#SBATCH --output=BreastCancer_sw.out  # Output file
#SBATCH --error=BreastCancer_sw.err   # Error file
#SBATCH --ntasks=1                                # Run on a single CPU
#SBATCH --mem=64G                                


# Load R module (if necessary, adjust this to match your system)
module load conda_R/4.4

# Run the R script
Rscript breast_cancer_spotsweeper.R
