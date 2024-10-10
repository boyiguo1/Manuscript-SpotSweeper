#!/bin/bash
#SBATCH --job-name=BayesSpace_mouse_brain   # Job name
#SBATCH --output=./logs/BayesSpace_mouse_brain.out  # Output file
#SBATCH --error=./logs/BayesSpace_mouse_brain.err   # Error file
#SBATCH --ntasks=1                                # Run on a single CPU
#SBATCH --mem=128G                                


# Load R module (if necessary, adjust this to match your system)
module load conda_R/4.4

# Run the R script
Rscript VisiumHD_mouse_brain_bayesspace.R
