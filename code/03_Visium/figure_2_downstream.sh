#!/bin/bash
#SBATCH --job-name=Visium_downstream    # Job name
#SBATCH --output=Visium_downstreamout  # Output file
#SBATCH --error=VisiumH_downstream.err   # Error file
#SBATCH --ntasks=1                                # Run on a single CPU
#SBATCH --mem=128G                                


# Load R module (if necessary, adjust this to match your system)
module load conda_R/4.4

# Run the R script
Rscript figure_2_downstream.R

