#!/bin/bash
#SBATCH --job-name=stomics_spotsweeper    # Job name
#SBATCH --output=stomics_spotsweeper.out  # Output file
#SBATCH --error=stomics_spotsweeper.err   # Error file
#SBATCH --ntasks=1                                # Run on a single CPU
#SBATCH --mem=64G                                


# Load R module (if necessary, adjust this to match your system)
module load conda_R/4.4

# Run the R script
Rscript stomics_spotsweeper.R