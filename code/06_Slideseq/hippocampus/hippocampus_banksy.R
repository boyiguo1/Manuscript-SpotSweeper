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


load(here("processed-data", "Slideseq","SlideseqV2Hippocampus", "slideseq_hippocampus.RData"))
# "countmat" "location"


spe <- SpatialExperiment(assay=countmat,
                           colData=location,
                           spatialCoordsNames=c("xcoord", "ycoord"))


counts(spe) <- assay(spe) 
                     

# ========== Banksy clustering ===========
library(Banksy)

lambda <- 0.8
k_geom <- 30
npcs <- 50
aname <- "counts"
spe <- Banksy::computeBanksy(spe, assay_name = aname, k_geom = k_geom)

set.seed(1000)
spe <- Banksy::runBanksyPCA(spe, lambda = lambda, npcs = npcs)

set.seed(1000)
spe <- Banksy::clusterBanksy(spe, lambda = lambda, npcs = npcs, resolution = 0.8)


saveRDS(spe, here("processed-data",  "Slideseq", "SlideseqV2Hippocampus", "slideseq_hippocampus_banksy.rds"))