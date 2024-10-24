library(SpotSweeper)
library(RANN)
library(here)
library(SpatialExperiment)
library(spatialLIBD)
library(escheR)
library(patchwork)
library(scran)
library(scater)


plot_dir <- here("plots","ovarian_cancer")


spe <- readRDS(here("processed-data", "spe_humanOvarian_preprocessed.rds"))
spe
# class: SpatialExperiment 
# dim: 33538 20123 
# metadata(0):
# assays(2): counts logcounts
# rownames(33538): ENSG00000243485 ENSG00000237613 ... ENSG00000277475
#   ENSG00000268674
# rowData names(2): gene_id gene_name
# colnames(20123): AAACAAGTATCTCCCA-1 AAACATTTCCCGGATT-1 ...
#   TTGTTTCCATACAACT-1 TTGTTTGTATTACACG-1
# colData names(5): in_tissue array_row array_col sample_id sizeFactor
# reducedDimNames(0):
# mainExpName: NULL
# altExpNames(0):
# spatialCoords names(2) : pxl_col_in_fullres pxl_row_in_fullres
# imgData names(4): sample_id image_id data scaleFactor

# ===== Add QC Metrics =====

# drop out of tissue
spe <- spe[, spe$in_tissue]



# ========== Banksy clustering ===========
library(Banksy)

lambda <- 0.2
k_geom <- 18
npcs <- 50
aname <- "counts"
spe <- Banksy::computeBanksy(spe, assay_name = aname, k_geom = k_geom)

set.seed(1000)
spe <- Banksy::runBanksyPCA(spe, lambda = lambda, npcs = npcs)

set.seed(1000)
spe <- Banksy::clusterBanksy(spe, lambda = lambda, npcs = npcs, resolution = 0.8)


saveRDS(spe, here("processed-data", "Visium","ovarian_cancer","OvarianCancer_banksy.rds"))
