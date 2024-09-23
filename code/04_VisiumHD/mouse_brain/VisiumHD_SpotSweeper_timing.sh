#!/bin/bash
#SBATCH --job-name=VisiumHD_SpotSweeper_timing    # Job name
#SBATCH --output=VisiumHD_SpotSweeper_timing.out  # Output file
#SBATCH --error=VisiumHD_SpotSweeper_timing.err   # Error file
#SBATCH --ntasks=1                                # Run on a single CPU
#SBATCH --mem=32G                                 # Memory request (32 GB)
#SBATCH --time=12:00:00                           # Time limit (12 hours)


# Load R module (if necessary, adjust this to match your system)
module load conda_R/4.4

# Run the R script
Rscript VisiumHD_SpotSweeper_timing.R
