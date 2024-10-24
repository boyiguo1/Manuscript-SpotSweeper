library(SpotSweeper)
library(RANN)
library(here)
library(SpatialExperiment)
library(spatialLIBD)
library(escheR)
library(patchwork)
library(scran)
library(scater)


plot_dir <- here("plots","VisiumHD","human_lung")


spe <- readRDS(here("processed-data", "VisiumHD", "VisiumHD_HumanLung_016.rds"))
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

lambda <- 0.8
k_geom <- 30
npcs <- 50
aname <- "counts"
spe <- Banksy::computeBanksy(spe, assay_name = aname, k_geom = k_geom)

set.seed(1000)
spe <- Banksy::runBanksyPCA(spe, lambda = lambda, npcs = npcs)

set.seed(1000)
spe <- Banksy::clusterBanksy(spe, lambda = lambda, npcs = npcs, resolution = 0.8)


saveRDS(spe, here("processed-data", "VisiumHD","VisiumHD_MouseBrain_016_banksy.rds"))
