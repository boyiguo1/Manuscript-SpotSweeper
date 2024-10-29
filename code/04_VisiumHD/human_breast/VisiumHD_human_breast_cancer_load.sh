#!/bin/bash
#SBATCH --job-name=HD_SPE_HumanBreast    # Job name
#SBATCH --output=./logs/HD_SPE_HumanBreast.out  # Output file
#SBATCH --error=./logs/HD_SPE_HumanBreast.err   # Error file
#SBATCH --ntasks=1                                # Run on a single CPU
#SBATCH --mem=128G                                


# Load R module (if necessary, adjust this to match your system)
module load conda_R/4.4

# Run the R script
Rscript VisiumHD_human_breast_cancer_load.R
