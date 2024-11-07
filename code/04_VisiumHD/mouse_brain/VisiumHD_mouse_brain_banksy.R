library(SpotSweeper)
library(RANN)
library(here)
library(SpatialExperiment)
library(spatialLIBD)
library(escheR)
library(patchwork)
library(scran)
library(scater)
library(Seurat)


plot_dir <- here("plots","VisiumHD","mouse_brain")


spe <- readRDS(here("processed-data", "VisiumHD","VisiumHD_MouseBrain_008.rds"))
spe


# ===== Add QC Metrics =====

# drop out of tissue
#spe <- spe[, spe$in_tissue]

# get mito genes
rownames(spe) <- rowData(spe)$symbol
is.mito <- grepl("^mt-", rownames(spe))

# subset to SPE to is mito
spe.mito <- spe[rownames(spe) %in% is.mito,]

# get qc metrics
df <- scuttle::perCellQCMetrics(spe, subsets=list(Mito=is.mito))

# add to colData
spe$sum <- df$sum
spe$detected <- df$detected
spe$subsets_mito_percent <- df$subsets_Mito_percent
spe$subsets_mito_sum <- df$subsets_Mito_sum

spe$sum_discard <- isOutlier(spe$sum, nmads=3, type="lower", log=TRUE)
spe$detected_discard <- isOutlier(spe$detected, nmads=3, type="lower", log=TRUE)
spe$subsets_mito_percent_discard <- isOutlier(spe$subsets_mito_percent, nmads=3, type="higher")


# ========== Banksy clustering ===========
library(Banksy)

set.seed(1000)
spe <- logNormCounts(spe)

# get only top 1000 HVGs
dec <- modelGeneVar(spe)
chosen <- getTopHVGs(dec, n=1000)

spe.hvg <- spe[chosen,]

lambda <- 0.8
k_geom <- 30
npcs <- 50
aname <- "counts"

spe.hvg <- Banksy::computeBanksy(spe.hvg, assay_name = aname, k_geom = k_geom)
spe.hvg <- Banksy::runBanksyPCA(spe.hvg, lambda = lambda, npcs = npcs)
spe.hvg <- Banksy::clusterBanksy(spe.hvg, lambda = lambda, npcs = npcs, resolution = 0.8)


saveRDS(spe.hvg, here("processed-data", "VisiumHD","VisiumHD_MouseBrain_008_banksy_hvgs.rds"))
