#!/bin/bash
#SBATCH --job-name=VisiumHD_Banksy_008    # Job name
#SBATCH --output=VisiumHD_Banksy_008.out  # Output file
#SBATCH --error=VisiumHD_Banksy_008.err   # Error file
#SBATCH --ntasks=1                                # Run on a single CPU
#SBATCH --mem=128G                                


# Load R module (if necessary, adjust this to match your system)
module load conda_R/4.4

# Run the R script
Rscript VisiumHD_mouse_brain_banksy_008.R

