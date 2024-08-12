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


plot_dir <- here("plots","xenium")

(sfe <- BiermannMelaMetasData(dataset = "MBM05_rep1"))
sfe
# class: SpatialFeatureExperiment
# dim: 27566 29536
# metadata(0):
#   assays(1): counts
# rownames(27566): A1BG A1BG-AS1 ... ZZZ3 snoZ196
# rowData names(3): means vars cv2
# colnames(29536): ACCACTCATTTCTC-1 GTTCANTCCACGTA-1 ... ACGCGCAATCGTAG-1 TTGTTCCGTTCATA-1
# colData names(4): sample_id nCounts nGenes prop_mito
# reducedDimNames(0):
#   mainExpName: NULL
# altExpNames(0):
#   spatialCoords names(2) : xcoord ycoord
# imgData names(1): sample_id
#
# unit: full_res_image_pixels
# Geometries:
#   colGeometries: centroids (POINT)
#
# Graphs:
#   sample01:

sfe$xcoord <- spatialCoords(sfe)[,1]
sfe$ycoord <- spatialCoords(sfe)[,2]

spe <- SpatialExperiment(assay=counts(sfe),
                         colData=colData(sfe),
                         spatialCoordsNames=c("xcoord", "ycoord"))

rowData(spe) <- rowData(sfe)
counts(spe) <- assay(sfe)



p1 <- plotQC(spe, metric="nCounts_log", point_size=1)
p2 <- plotQC(spe, metric="nGenes_log", point_size=1)
p1+p2


p1 <- plotColData(spe, x="sample_id", y="nCounts") +
  scale_y_log10()
p2 <- plotColData(spe, x="sample_id", y="nGenes") +
  scale_y_log10()
p1+p2

# ======== SpotSweeper ===========

spe <- localOutliers(spe, metric="nCounts", direction="lower", log=TRUE, n_neighbors=18)
spe <- localOutliers(spe, metric="nGenes", direction="lower", log=TRUE, n_neighbors=18)
colnames(colData(spe))

png(file.path(plot_dir, "spotplot_spotsweeper_nCounts_discard.png"))
p1 <- plotQC(spe, metric="nCounts_log", outliers="nCounts_outliers", point_size=1)
p1
dev.off()



plotColData(spe, x="sample_id", y="nCounts_z", color="nCounts_outliers")
plotColData(spe, x="sample_id", y="nGenes_z", color="nGenes_outliers")
