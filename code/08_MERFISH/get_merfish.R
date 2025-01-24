library(SpotSweeper)
library(RANN)
library(here)
library(SpatialExperiment)
library(spatialLIBD)
library(escheR)
library(patchwork)
library(scran)
library(scater)
library(SpatialFeatureExperiment)
library(SFEData)
library(RColorBrewer)


plot_dir <- here("plots","Xenium")

spe <- STexampleData::STARmapPLUS_mouseBrain()
spe
# class: SpatialExperiment 
# dim: 1022 46184 
# metadata(0):
# assays(1): counts
# rownames(1022): A2M ABCC9 ... ZIC1 ZMYM1
# rowData names(0):
# colnames(46184): well05_0 well05_1 ... well05_46182 well05_46183
# colData names(7): NAME Main_molecular_cell_type ...
#   Molecular_spatial_cell_type sample_id
# reducedDimNames(0):
# mainExpName: NULL
# altExpNames(0):
# spatialCoords names(3) : X Y Z
# imgData names(0):


# ========= Basic QC ========
outs <- perCellQCMetrics(spe)

spe$nCounts <- outs$sum
spe$nGenes <- outs$detected

spe$nCounts_discard_threshold <- spe$nCounts < 30
spe$nGenes_discard_threshold  <- spe$nGenes < 50

spe$nCounts_discard_3mad <- isOutlier(spe$nCounts, nmads=3, type="lower", log=TRUE)
spe$nGenes_discard_3mad <- isOutlier(spe$nGenes, nmads=3, type="lower", log=TRUE)

# local outlier: physical distance
spe <- localOutliers(spe, metric="nCounts", direction="lower", log=TRUE)
spe <- localOutliers(spe, metric="nGenes", direction="lower", log=TRUE)


# save
saveRDS(spe, here("processed-data", "MERFISH", "merfish_spotsweeper.rds"))
