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


spe <- readRDS(here("processed-data", "VisiumHD", "VisiumHD_MouseBrain_016.rds"))
spe

# ===== Add QC Metrics =====

# drop out of tissue
#spe <- spe[, spe$in_tissue]

# get mito genes
rownames(spe) <- rowData(spe)$symbol
is.mito <- grepl("^MT-", rownames(spe))

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
library(BayesSpace)

## Add BayesSpace metadata
spe <- spatialPreprocess(spe, platform="Visium", skip.PCA=FALSE)

colData(spe)$row <- spe$array_row
colData(spe)$col <- spe$array_col

metadata(spe)$BayesSpace.data <- list(platform = "Visium", is.enhanced = FALSE)

d <- 15  # Number of PCs

## Run BayesSpace clustering
set.seed(104)
dlpfc <- spatialCluster(spe, q=10, d=d, platform='ST',
                        nrep=10000, gamma=3, save.chain=TRUE)
spe$BS_k10 <- as.factor(dlpfc$spatial.cluster)

# save 
saveRDS(spe, here("processed-data", "VisiumHD","mouse_brain", "VisiumHD_MouseBrain_016_bayesspace.rds"))