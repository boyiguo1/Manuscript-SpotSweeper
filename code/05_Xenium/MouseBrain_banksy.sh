#!/bin/bash
#SBATCH --job-name=Xenium_Banksy    # Job name
#SBATCH --output=Xenium_Banksy.out  # Output file
#SBATCH --error=Xenium_Banksy.err   # Error file
#SBATCH --ntasks=1                                # Run on a single CPU
#SBATCH --mem=128G                                


# Load R module (if necessary, adjust this to match your system)
module load conda_R/4.4

# Run the R script
Rscript MouseBrain_banksy.R
