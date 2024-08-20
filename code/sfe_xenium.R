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
library(ggVennDiagram)
library(RColorBrewer)


plot_dir <- here("plots","xenium")

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
file <- here("processed-data","Cell_Barcode_Type_Matrices.xlsx")
cell_barcode_type <- readxl::read_excel(file, sheet="Xenium R2 Fig1-5 (supervised)")

# add Cluster to colData based on matching Barcode
spe$Cluster <- cell_barcode_type$Cluster[match(spe$Barcode, cell_barcode_type$Barcode)]

# drop single NA cell
spe <- spe[,!is.na(spe$Cluster)]


# ========= Plotting clusters =========
# unique palette for clusters
cluster_pal <- colorRampPalette(brewer.pal(9, "Set1"))(length(unique(spe$Cluster)))
png(here(plot_dir, "breastcancer_spotplot_clusters.png"), width=12, height=10, units="in", res=300)
make_escheR(spe) |>
  add_ground(var="Cluster", point_size=.3) +
  scale_color_manual(values=cluster_pal)
dev.off()

# ========= QC metrics =========
png(here(plot_dir, "breastcancer_spotplot_qc_metrics.png"), width=15, height=7.5, units="in", res=300)
p1 <- plotQC(spe, metric="nCounts", point_size=.7)
p2 <- plotQC(spe, metric="nGenes", point_size=.7)
p1+p2
dev.off()

spe$nCounts_log <- log1p(spe$nCounts)
png(here(plot_dir, "breastcancer_violin_qc_metrics.png"), width=6, height=9, units="in", res=300)
p1 <- plotColData(spe, x="Cluster", y="nCounts_log", color_by="Cluster") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position="none") +
  labs(title="Library Size",
       x="Cell Type") +
  coord_flip()
p2 <- plotColData(spe, x="Cluster", y="nGenes", color_by="Cluster") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position="none") +
  labs(title="Num Genes",
       x=NULL) +
  coord_flip()
p1+p2
dev.off()


# ======== SpotSweeper ========
spe <- localOutliers(spe, metric="nCounts", direction="lower", log=TRUE)
spe <- localOutliers(spe, metric="nGenes", direction="lower", log=TRUE)

png(here(plot_dir, "breastcancer_spotplot_discard_spotsweeper_newSamples.png"), width=15, height=7.5, units="in", res=300)
p1 <- plotQC(spe, metric="nCounts_log", outliers="nCounts_outliers", point_size=1)
p2 <- plotQC(spe, metric="nGenes_log", outliers="nGenes_outliers", point_size=1)
p1+p2
dev.off()

png(here(plot_dir, "breastcancer_violin_discard_spotsweeper.png"), width=8, height=9, units="in", res=300)
p1 <- plotColData(spe, x="Cluster", y="nCounts_log", color_by="nCounts_outliers") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values=c("TRUE"="red", "FALSE"="grey")) +
  theme(legend.position="none") +
  labs(title="Library Size",
       x="Cell Type") +
  coord_flip()
p2 <- plotColData(spe, x="Cluster", y="nGenes", color_by="nGenes_outliers") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values=c("TRUE"="red", "FALSE"="grey")) +
  #theme(legend.position="none") +
  labs(title="Num Genes",
       x=NULL) +
  coord_flip()
p1+p2
dev.off()



# ======== Scater qc =======

# ====== arbitrary thresholding ======
spe$nCounts_discard_threshold <- spe$nCounts < 25
spe$nGenes_discard_threshold  <- spe$nGenes < 10

png(here(plot_dir, "breastcancer_spotplot_discard_threshold.png"), width=15, height=7.5, res=300, units="in")
p1 <- plotQC(spe, metric="nCounts_log", outliers="nCounts_discard_threshold", point_size=.7) +
  ggtitle("Sum Counts")

p2 <- plotQC(spe, metric="nGenes_log", outliers="nGenes_discard_threshold", point_size=.7) +
  ggtitle("Genes Detected")

p1+p2
dev.off()

png(here(plot_dir, "breastcancer_violin_discard_threshold.png"), width=8, height=9, units="in", res=300)
p1 <- plotColData(spe, x="Cluster", y="nCounts_log", color_by="nCounts_discard_threshold") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values=c("TRUE"="red", "FALSE"="grey")) +
  theme(legend.position="none") +
  labs(title="Library Size",
       x="Cell Type") +
  coord_flip()
p2 <- plotColData(spe, x="Cluster", y="nGenes", color_by="nGenes_discard_threshold") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values=c("TRUE"="red", "FALSE"="grey")) +
  #theme(legend.position="none") +
  labs(title="Num Genes",
       x=NULL) +
  coord_flip()
p1+p2
dev.off()



# ====== 3 MAD ======
spe$nCounts_discard_3mad <- isOutlier(spe$nCounts, nmads=3, type="lower", log=TRUE)
spe$nGenes_discard_3mad <- isOutlier(spe$nGenes, nmads=3, type="lower", log=TRUE)

png(here(plot_dir, "breastcancer_spotplot_discard_3mad.png"), width=15, height=7.5, res=300, units="in")
p1 <- plotQC(spe, metric="nCounts_log", outliers="nCounts_discard_3mad", point_size=.7) +
  ggtitle("Sum Counts")

p2 <- plotQC(spe, metric="nGenes_log", outliers="nGenes_discard_3mad", point_size=.7) +
  ggtitle("Genes Detected")

p1+p2
dev.off()


png(here(plot_dir, "breastcancer_violin_discard_3mad.png"), width=8, height=9, units="in", res=300)
p1 <- plotColData(spe, x="Cluster", y="nCounts_log", color_by="nCounts_discard_3mad") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values=c("TRUE"="red", "FALSE"="grey")) +
  theme(legend.position="none") +
  labs(title="Library Size",
       x="Cell Type") +
  coord_flip()
p2 <- plotColData(spe, x="Cluster", y="nGenes", color_by="nGenes_discard_3mad") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values=c("TRUE"="red", "FALSE"="grey")) +
  #theme(legend.position="none") +
  labs(title="Num Genes",
       x=NULL) +
  coord_flip()
p1+p2
dev.off()




# ========= Violin plots ==========


png(here(plot_dir, "breastcancer_violins_discard.png"), width=10, height=10, units="in", res=300)
p1 <- plotColData(spe, x="sample_id", y="nCounts_log", color="nCounts_discard_threshold")
p2 <- plotColData(spe, x="sample_id", y="nGenes_log", color="nGenes_discard_threshold")

p3 <- plotColData(spe, x="sample_id", y="nCounts_log", color="nCounts_discard_3mad")
p4 <- plotColData(spe, x="sample_id", y="nGenes_log", color="nGenes_discard_3mad")

p5 <- plotColData(spe, x="sample_id", y="nCounts_z", color="nCounts_outliers")
p6 <- plotColData(spe, x="sample_id", y="nGenes_z", color="nGenes_outliers")
(p1+p2) / (p3+p4) / (p5+p6)
dev.off()



# ========= Venn diagram ========

# nCounts
discard_list <- list(
  Threshold = which(spe$nCounts_discard_threshold),
  ThreeMAD = which(spe$nCounts_discard_3mad),
  SpotSweeper = which(spe$nCounts_outliers)
)

png(here(plot_dir, "breastcancer_vennDiagram_nCounts_discard.png"), width=5, height=5, units="in", res=300)
ggVennDiagram(discard_list) +
  scale_fill_gradient(low="white",high = "red") +
  ggtitle("nCounts")
dev.off()

# nCounts
discard_list <- list(
  Threshold = which(spe$nGenes_discard_threshold),
  ThreeMAD = which(spe$nGenes_discard_3mad),
  SpotSweeper = which(spe$nGenes_outliers)
)

png(here(plot_dir, "breastcancer_vennDiagram_nGenes_discard.png"), width=5, height=5, units="in", res=300)
ggVennDiagram(discard_list) +
  scale_fill_gradient(low="white",high = "red") +
  ggtitle("nGenes")
dev.off()




# ========= Norm counts ==========

spe$nCounts_norm <- spe$nCounts / spe$cell_area

plotQC(spe, metric="nCounts_norm")
plotQC(spe, metric="nCounts")




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

png(here(plot_dir, "breastcancer_banksy_spotplots.png"), width=10, height=10, units="in", res=300)
make_escheR(spe) |>
  add_ground(var="clust_M0_lam0.8_k50_res0.8", point_size=.3)
dev.off()

spe$banksy_clusters <- spe$clust_M0_lam0.8_k50_res0.8

# violin plots
png(here(plot_dir, "breastcancer_banksy_violins_discard_threshold.png"), width=8, height=9, units="in", res=300)
p1 <- plotColData(spe, x="clust_M0_lam0.8_k50_res0.8", y="nCounts_log", color_by="nCounts_discard_threshold") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values=c("TRUE"="red", "FALSE"="grey")) +
  theme(legend.position="none") +
  labs(title="Library Size",
       x="Cell Type") +
  coord_flip()
p2 <- plotColData(spe, x="clust_M0_lam0.8_k50_res0.8", y="nGenes", color_by="nGenes_discard_threshold") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values=c("TRUE"="red", "FALSE"="grey")) +
  #theme(legend.position="none") +
  labs(title="Num Genes",
       x=NULL) +
  coord_flip()
p1+p2
dev.off()


png(here(plot_dir, "breastcancer_banksy_violins_discard_3mad.png"), width=8, height=9, units="in", res=300)
p1 <- plotColData(spe, x="clust_M0_lam0.8_k50_res0.8", y="nCounts_log", color_by="nCounts_discard_3mad") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values=c("TRUE"="red", "FALSE"="grey")) +
  theme(legend.position="none") +
  labs(title="Library Size",
       x="Cell Type") +
  coord_flip()
p2 <- plotColData(spe, x="clust_M0_lam0.8_k50_res0.8", y="nGenes", color_by="nGenes_discard_3mad") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values=c("TRUE"="red", "FALSE"="grey")) +
  #theme(legend.position="none") +
  labs(title="Num Genes",
       x=NULL) +
  coord_flip()
p1+p2
dev.off()


png(here(plot_dir, "breastcancer_banksy_violins_discard_spotsweeper.png"), width=8, height=9, units="in", res=300)
p1 <- plotColData(spe, x="clust_M0_lam0.8_k50_res0.8", y="nCounts_log", color_by="nCounts_outliers") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values=c("TRUE"="red", "FALSE"="grey")) +
  theme(legend.position="none") +
  labs(title="Library Size",
       x="Cell Type") +
  coord_flip()
p2 <- plotColData(spe, x="clust_M0_lam0.8_k50_res0.8", y="nGenes", color_by="nGenes_outliers") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values=c("TRUE"="red", "FALSE"="grey")) +
  #theme(legend.position="none") +
  labs(title="Num Genes",
       x=NULL) +
  coord_flip()
p1+p2
dev.off()



# ======= UMAPs =======

plotReducedDim(spe, dimred="PCA_M0_lam0.8", color_by="clust_M0_lam0.8_k50_res0.8")

spe <- runUMAP(spe, assay.type="counts", ncomponents = 30, ntop = 500)

png(here(plot_dir, "breastcancer_umap_banksy_clusters.png"), width=5, height=5, units="in", res=300)
plotReducedDim(spe, dimred="UMAP", color_by="clust_M0_lam0.8_k50_res0.8")
dev.off()

png(here(plot_dir, "breastcancer_umap_banksy_nCounts.png"),  width=5, height=5, units="in", res=300)
plotReducedDim(spe, dimred="UMAP", color_by="nCounts_log")
dev.off()

png(here(plot_dir, "breastcancer_umap_banksy_nCounts_threshold.png"),  width=5, height=5, units="in", res=300)
plotReducedDim(spe, dimred="UMAP", color_by="nCounts_discard_threshold") +
  scale_color_manual("color", values=c("TRUE"="red", "FALSE"="grey"))
dev.off()

png(here(plot_dir, "breastcancer_umap_banksy_nCounts_3mad.png"),  width=5, height=5, units="in", res=300)
plotReducedDim(spe, dimred="UMAP", color_by="nCounts_discard_3mad")+
  scale_color_manual("color", values=c("TRUE"="red", "FALSE"="grey"))
dev.off()

png(here(plot_dir, "breastcancer_umap_banksy_nCounts_spotsweeper.png"),  width=5, height=5, units="in", res=300)
plotReducedDim(spe, dimred="UMAP", color_by="nCounts_outliers")+
  scale_color_manual("color", values=c("TRUE"="red", "FALSE"="grey"))
dev.off()
