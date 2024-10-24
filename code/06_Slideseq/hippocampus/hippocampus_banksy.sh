#!/bin/bash
#SBATCH --job-name=HPC_Banksy    # Job name
#SBATCH --output=HPC_Banksy.out  # Output file
#SBATCH --error=HPC_Banksy.err   # Error file
#SBATCH --ntasks=1                                # Run on a single CPU
#SBATCH --mem=90G                                


# Load R module (if necessary, adjust this to match your system)
module load conda_R/4.4

# Run the R script
Rscript hippocampus_banksy.R
