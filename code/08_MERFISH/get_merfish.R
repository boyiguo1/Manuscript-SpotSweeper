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
library(MerfishData)


plot_dir <- here("plots","Xenium")

spe <- MouseColonIbdCadinu2024()
spe
# class: SpatialExperiment 
# dim: 943 1348582 
# metadata(3): colors_tier1 colors_tier2 colors_tier3
# assays(2): counts logcounts
# rownames(943): Ackr1 Ackr2 ... Zc3h12a Znrf3
# rowData names(0):
# colnames: NULL
# colData names(12): cell_id sample_id ... tier3 leiden_neigh
# reducedDimNames(5): PCA UMAP UMAP_Tier1 UMAP_Tier2 UMAP_Tier3
# mainExpName: NULL
# altExpNames(1): Blank
# spatialCoords names(2) : x y
# imgData names(0):

colnames(colData(spe))
#  [1] "cell_id"                 "sample_id"              
#  [3] "sample_type"             "mouse_id"               
#  [5] "technical_repeat_number" "slice_id"               
#  [7] "fov"                     "n_genes"                
#  [9] "tier1"                   "tier2"                  
# [11] "tier3"                   "leiden_neigh"  

# subset to mouse_id = "082421_D0_m6" and spe$sample_type == "Healthy" & spe$technical_repeat_number == "1" & spe$slice_id == "2"
spe.subset <- spe[, spe$mouse_id == "082421_D0_m6" & spe$sample_type == "Healthy" & spe$technical_repeat_number == "1" & spe$slice_id == "2"]
spe.subset
# class: SpatialExperiment 
# dim: 943 27785 
# metadata(3): colors_tier1 colors_tier2 colors_tier3
# assays(2): counts logcounts
# rownames(943): Ackr1 Ackr2 ... Zc3h12a Znrf3
# rowData names(0):
# colnames: NULL
# colData names(12): cell_id sample_id ... tier3 leiden_neigh
# reducedDimNames(5): PCA UMAP UMAP_Tier1 UMAP_Tier2 UMAP_Tier3
# mainExpName: NULL
# altExpNames(1): Blank
# spatialCoords names(2) : x y
# imgData names(0):

spe <- spe.subset

# ========= Basic QC ========
outs <- perCellQCMetrics(spe)

spe$nCounts <- outs$sum
spe$nGenes <- outs$detected

spe$nCounts_discard_3mad <- isOutlier(spe$nCounts, nmads=2, type="lower", log=TRUE)
spe$nGenes_discard_3mad <- isOutlier(spe$nGenes, nmads=2, type="lower", log=TRUE)

# local outlier: physical distance
spe <- localOutliers(spe, metric="nCounts", direction="lower", log=TRUE, n_neighbors=50)
spe <- localOutliers(spe, metric="nGenes", direction="lower", log=TRUE, n_neighbors=50)


# save
saveRDS(spe, here("processed-data", "MERFISH", "merfish_spotsweeper.rds"))
