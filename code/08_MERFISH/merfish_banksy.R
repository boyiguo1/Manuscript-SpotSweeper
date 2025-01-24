library(SpotSweeper)
library(RANN)
library(here)
library(SpatialExperiment)
library(spatialLIBD)
library(escheR)
library(patchwork)
library(scran)
library(scater)


spe <- readRDS(here("processed-data",  "MERFISH", "merfish_spotsweeper.rds"))
spe


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


saveRDS(spe, here("processed-data",  "MERFISH",  "merfish_spotsweeper_banksy.rds"))
