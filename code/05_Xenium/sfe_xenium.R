remotes::install_github("MicTott/SpotSweeper")

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

(sfe <- JanesickBreastData(dataset = "rep2"))
sfe
# class: SpatialFeatureExperiment
# dim: 541 118708
# metadata(1): Samples
# assays(1): counts
# rownames(541): ABCC11 ACTA2 ... BLANK_0497 BLANK_0499
# rowData names(6): ID Symbol ... vars cv2
# colnames: NULL
# colData names(10): Sample Barcode ... nCounts nGenes
# reducedDimNames(0):
#   mainExpName: NULL
# altExpNames(0):
#   spatialCoords names(2) : x_centroid y_centroid
# imgData names(1): sample_id
#
# unit:
#   Geometries:
#   colGeometries: centroids (POINT), cellSeg (POLYGON), nucSeg (GEOMETRY)
#
# Graphs:
#   sample01:

# ======== Convert to SpatialExperiment ========
sfe$xcoord <- spatialCoords(sfe)[,1]
sfe$ycoord <- spatialCoords(sfe)[,2]

spe <- SpatialExperiment(assay=counts(sfe),
                         colData=colData(sfe),
                         spatialCoordsNames=c("xcoord", "ycoord"))

rowData(spe) <- rowData(sfe)
counts(spe) <- assay(sfe)

# ========= Read in cell type annotations from 10x ========
# read in excel
library(readxl)

# read in the excel file - sheet "Xenium R1 Fig1-5 (supervised)"
file <- here("processed-data","Xenium", "Cell_Barcode_Type_Matrices.xlsx")
cell_barcode_type <- readxl::read_excel(file, sheet="Xenium R2 Fig1-5 (supervised)")

# add Cluster to colData based on matching Barcode
spe$Cluster <- cell_barcode_type$Cluster[match(spe$Barcode, cell_barcode_type$Barcode)]

# drop single NA cell
spe <- spe[,!is.na(spe$Cluster)]
colnames(colData(spe))
#  [1] "Sample"                  "Barcode"                
#  [3] "transcript_counts"       "control_probe_counts"   
#  [5] "control_codeword_counts" "cell_area"              
#  [7] "nucleus_area"            "sample_id"              
#  [9] "nCounts"                 "nGenes"                 
# [11] "Cluster" 


# ========= Basic QC ========
spe$nCounts_discard_threshold <- spe$nCounts < 25
spe$nGenes_discard_threshold  <- spe$nGenes < 10

spe$nCounts_discard_3mad <- isOutlier(spe$nCounts, nmads=3, type="lower", log=TRUE)
spe$nGenes_discard_3mad <- isOutlier(spe$nGenes, nmads=3, type="lower", log=TRUE)

# local outlier: physical distance
spe <- localOutliers(spe, metric="nCounts", direction="lower", log=TRUE)
spe <- localOutliers(spe, metric="nGenes", direction="lower", log=TRUE)


# save
saveRDS(spe, here("processed-data", "Xenium", "Xenium_spotsweeper.rds"))

