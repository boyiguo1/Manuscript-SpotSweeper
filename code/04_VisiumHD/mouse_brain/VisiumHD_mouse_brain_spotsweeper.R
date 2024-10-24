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

spe <- readRDS(here("processed-data", "VisiumHD", "mouse_brain", "VisiumHD_MouseBrain_016_banksy.rds"))
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

# ====== Add discards by fixed thresholds ========
spe$detected_threshold <- spe$detected < 500
spe$sum_threshold <- spe$sum < 500
spe$subsets_mito_percent_threshold <- spe$subsets_mito_percent > 30

# ======= SpotSweeper =========
spe <- localOutliers(spe, metric = "sum", direction = "lower", log = TRUE)
spe <- localOutliers(spe, metric = "detected", direction = "lower", log = TRUE)
spe <- localOutliers(spe, metric = "subsets_mito_percent", direction = "higher", log = TRUE)

saveRDS(spe, here("processed-data", "VisiumHD", "mouse_brain", "VisiumHD_MouseBrain_016_banksy_spotsweeper.rds"))
