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


load(here("processed-data", "Slideseq","SlideseqCerebellum", "slideseq_cerebellum.Rdata"))
# "sp_count" "location"

spe <- SpatialExperiment(assay=sp_count,
                           colData=location,
                           spatialCoordsNames=c("xcoord", "ycoord"))

counts(spe) <- assay(spe) 

# ========== Banksy clustering ===========
library(Banksy)

lambda <- 0.3
k_geom <- 30
npcs <- 50
aname <- "counts"
spe <- Banksy::computeBanksy(spe, assay_name = aname, k_geom = k_geom)

set.seed(1000)
spe <- Banksy::runBanksyPCA(spe, lambda = lambda, npcs = npcs)

set.seed(1000)
spe <- Banksy::clusterBanksy(spe, lambda = lambda, npcs = npcs, resolution = 0.7)


saveRDS(spe, here("processed-data",  "Slideseq","SlideseqCerebellum", "slideseq_cerebellum_banksy.rds"))

