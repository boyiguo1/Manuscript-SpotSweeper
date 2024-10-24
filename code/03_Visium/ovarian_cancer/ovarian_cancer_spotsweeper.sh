#!/bin/bash
#SBATCH --job-name=OvarianCancer_sw    # Job name
#SBATCH --output=OvarianCancer_sw.out  # Output file
#SBATCH --error=OvarianCancer_sw.err   # Error file
#SBATCH --ntasks=1                                # Run on a single CPU
#SBATCH --mem=64G                                


# Load R module (if necessary, adjust this to match your system)
module load conda_R/4.4

# Run the R script
Rscript ovarian_cancer_spotsweeper.R
