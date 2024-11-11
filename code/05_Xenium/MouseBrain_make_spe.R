library(SpotSweeper)
library(here)
library(SpatialExperiment)
library(spatialLIBD)
library(escheR)
library(patchwork)
library(scran)
library(scater)
library(RColorBrewer)
library(SpatialExperimentIO)

plot_dir <- here("plots","Xenium", "mouse_brain")

spe <- readXeniumSXE(here("raw-data", "Xenium", "MouseBrain_5k", "output"))
spe
# class: SpatialExperiment 
# dim: 13780 63173 
# metadata(0):
# assays(1): counts
# rownames(13780): ENSMUSG00000052595 ENSMUSG00000030111 ...
#   DeprecatedCodeword_16059 DeprecatedCodeword_18533
# rowData names(3): ID Symbol Type
# colnames(63173): aaaalabd-1 aaabcnhh-1 ... oilgmbbn-1 oilhbhom-1
# colData names(13): cell_id transcript_counts ... segmentation_method
#   sample_id
# reducedDimNames(0):
# mainExpName: NULL
# altExpNames(0):
# spatialCoords names(2) : x_centroid y_centroid
# imgData names(0):

colnames(colData(spe)) 
#  [1] "cell_id"                    "transcript_counts"         
#  [3] "control_probe_counts"       "genomic_control_counts"    
#  [5] "control_codeword_counts"    "unassigned_codeword_counts"
#  [7] "deprecated_codeword_counts" "total_counts"              
#  [9] "cell_area"                  "nucleus_area"              
# [11] "nucleus_count"              "segmentation_method"       
# [13] "sample_id"  

unique(rowData(spe)$Type)
# [1] "Gene Expression"           "Negative Control Probe"   
# [3] "Genomic Control"           "Negative Control Codeword"
# [5] "Unassigned Codeword"       "Deprecated Codeword"

# subset to only "Gene Expression" genes
spe <- spe[rowData(spe)$Type == "Gene Expression",]
spe
# class: SpatialExperiment 
# dim: 5006 63173 
# metadata(0):
# assays(1): counts
# rownames(5006): ENSMUSG00000052595 ENSMUSG00000030111 ...
#   ENSMUSG00000055670 ENSMUSG00000027596
# rowData names(3): ID Symbol Type
# colnames(63173): aaaalabd-1 aaabcnhh-1 ... oilgmbbn-1 oilhbhom-1
# colData names(13): cell_id transcript_counts ... segmentation_method
#   sample_id
# reducedDimNames(0):
# mainExpName: NULL
# altExpNames(0):
# spatialCoords names(2) : x_centroid y_centroid
# imgData names(0):


# ============== Adding QC metrics ===============


# get qc metrics
df <- scuttle::perCellQCMetrics(spe)

# add to colData
spe$sum <- df$sum
spe$detected <- df$detected

# drop any sum or detected = 0 in spe
spe <- spe[,spe$sum > 0]
spe <- spe[,spe$detected > 0]

spe$sum_discard <- isOutlier(spe$sum, nmads=3, type="lower", log=TRUE)
spe$detected_discard <- isOutlier(spe$detected, nmads=3, type="lower", log=TRUE)


# ====== Add discards by fixed thresholds ========
spe$detected_threshold <- spe$detected < 100
spe$sum_threshold <- spe$sum < 100


# ======= SpotSweeper =========
spe <- localOutliers(spe,  n_neighbors = 36, metric = "sum", direction = "lower", log = TRUE)
spe <- localOutliers(spe,  n_neighbors = 36, metric = "detected", direction = "lower", log = TRUE)

# save
saveRDS(spe, here("processed-data",  "Xenium", "MouseBrain_5k", "Xenium_MouseBrain_spotsweeper.rds"))