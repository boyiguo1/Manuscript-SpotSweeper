#!/bin/bash
#SBATCH --job-name=Xenium_spotsweeper    # Job name
#SBATCH --output=./logs/Xenium_spotsweeper.out  # Output file
#SBATCH --error=./logs/Xenium_spotsweeper.err   # Error file
#SBATCH --ntasks=1                                # Run on a single CPU
#SBATCH --mem=128G                                


# Load R module (if necessary, adjust this to match your system)
module load conda_R/4.4

# Run the R script
Rscript sfe_xenium.R