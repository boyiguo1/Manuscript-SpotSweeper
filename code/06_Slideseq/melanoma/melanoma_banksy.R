library(SpotSweeper)
library(RANN)
library(here)
library(SpatialExperiment)
library(spatialLIBD)
library(escheR)
library(patchwork)
library(scran)
library(scater)
library(tidyr)
library(dplyr)


plot_dir <- here("plots","Slideseq")

# open .csv
location <- read.csv(here("processed-data", "Slideseq","Melanoma","melanoma_spatial.csv"))
counts <- read.csv(here("processed-data", "Slideseq","Melanoma","melanoma_counts.csv"))
# "sp_count" "location"

spe <- SpatialExperiment(assay=counts,
                           colData=location,
                           spatialCoordsNames=c("xcoord", "ycoord"))
spe

counts(spe) <- assay(spe) 