library(SpotSweeper)
library(RANN)
library(here)
library(SpatialExperiment)
library(spatialLIBD)
library(escheR)
library(patchwork)
library(scran)
library(scater)


plot_dir <- here("plots","VisiumHD","mouse_brain")


spe <- readRDS(here("processed-data", "VisiumHD","VisiumHD_MouseBrain_002.rds"))
spe


library(data.table)
barcodes = fread(here("raw-data","VisiumHD", "Visium_HD_Mouse_Brain", "output","binned_outputs", "square_002um", "raw_feature_bc_matrix", "barcodes.tsv.gz"))
head(barcodes)