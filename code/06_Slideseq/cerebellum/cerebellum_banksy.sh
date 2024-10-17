#!/bin/bash
#SBATCH --job-name=Cerebellum_Banksy    # Job name
#SBATCH --output=Cerebellum_Banksy.out  # Output file
#SBATCH --error=Cerebellum_Banksy.err   # Error file
#SBATCH --ntasks=1                                # Run on a single CPU
#SBATCH --mem=64G                                


# Load R module (if necessary, adjust this to match your system)
module load conda_R/4.4

# Run the R script
Rscript cerebellum_banksy.R
