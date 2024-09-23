library(SpotSweeper)
library(RANN)
library(here)
library(SpatialExperiment)
library(spatialLIBD)
library(escheR)
library(patchwork)
library(scran)
library(scater)


plot_dir <- here("plots","ovarian_cancer")


spe <- readRDS(here("processed-data", "spe_humanOvarian_preprocessed.rds"))
spe
# class: SpatialExperiment 
# dim: 33538 20123 
# metadata(0):
# assays(2): counts logcounts
# rownames(33538): ENSG00000243485 ENSG00000237613 ... ENSG00000277475
#   ENSG00000268674
# rowData names(2): gene_id gene_name
# colnames(20123): AAACAAGTATCTCCCA-1 AAACATTTCCCGGATT-1 ...
#   TTGTTTCCATACAACT-1 TTGTTTGTATTACACG-1
# colData names(5): in_tissue array_row array_col sample_id sizeFactor
# reducedDimNames(0):
# mainExpName: NULL
# altExpNames(0):
# spatialCoords names(2) : pxl_col_in_fullres pxl_row_in_fullres
# imgData names(4): sample_id image_id data scaleFactor

# ===== Add QC Metrics =====

# drop out of tissue
spe <- spe[, spe$in_tissue]

# get mito genes
rownames(spe) <- rowData(spe)$gene_name
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

colnames(colData(spe))
# [1] "in_tissue"            "array_row"            "array_col"            "sample_id"            "sizeFactor"
# [6] "sum"                  "detected"             "subsets_mito_percent" "subsets_mito_sum"

# ===== Spatially naive QC ======
# plot metrics
png(here(plot_dir, "cancer_spotplot_metrics.png"), width=12, height=5, res=300, units="in")
p1 <- plotQC(spe, metric="detected")
p2 <- plotQC(spe, metric="sum")
p3 <- plotQC(spe, metric="subsets_mito_percent")

p1+p2+p3
dev.off()

png(here(plot_dir, "cancer_spotplot_metrics_spatialLIBD.png"), width=5, height=5, res=300, units="in")
vis_gene(spe, geneid="detected")
dev.off()

png(here(plot_dir, "cancer_spotplot_detected_grid.png"), width=15, height=15, res=300, units="in")
p_list <- vis_grid_gene(spe,
        geneid="detected",
        spatial = FALSE,
        return_plots = TRUE
    )

cowplot::plot_grid(plotlist = p_list, ncol = 3)
dev.off()

png(here(plot_dir, "cancer_spotplot_sum_grid.png"), width=15, height=15, res=300, units="in")
p_list <- vis_grid_gene(spe,
        geneid="sum",
        spatial = FALSE,
        return_plots = TRUE
    )

cowplot::plot_grid(plotlist = p_list, ncol = 3)
dev.off()


png(here(plot_dir, "cancer_spotplot_percentMito_grid.png"), width=15, height=15, res=300, units="in")
p_list <- vis_grid_gene(spe,
        geneid="subsets_mito_percent",
        spatial = FALSE,
        return_plots = TRUE
    )

cowplot::plot_grid(plotlist = p_list, ncol = 3)
dev.off()
# run QC using 3 MAD

spe$sum_discard <- isOutlier(spe$sum, nmads=3, type="lower", log=TRUE)
spe$detected_discard <- isOutlier(spe$detected, nmads=3, type="lower", log=TRUE)
spe$subsets_mito_percent_discard <- isOutlier(spe$subsets_mito_percent, nmads=3, type="higher")


png(here(plot_dir, "cancer_violins_discarded_3mad.png"), width=15, height=3, res=300, units="in")
p1 <- plotColData(spe, x="sample_id", y="sum", color="sum_discard") +
   scale_y_log10() +
  ggtitle("Sum UMI") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p2 <- plotColData(spe, x="sample_id", y="detected", color="detected_discard") +
  scale_y_log10() +
  ggtitle("Gene Detected") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p3 <- plotColData(spe, x="sample_id", y="subsets_mito_percent", color="subsets_mito_percent_discard") +
  ggtitle("Mito Percent") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p1+p2+p3
dev.off()


png(here(plot_dir, "cancer_spotplot_discarded_3mad.png"), width=15, height=5, res=300, units="in")
p1 <- plotQC(spe, metric="sum", outliers="sum_discard") +
  ggtitle("Sum UMI")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p2 <- plotQC(spe, metric="detected", outliers="detected_discard") +
  ggtitle("Genes Detected")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p3 <- plotQC(spe, metric="subsets_mito_percent", outliers="subsets_mito_percent_discard") +
  ggtitle("Mito Percent")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p1+p2+p3
dev.off()


# make TRUE FALSE into DISCARD KEEP
spe$subsets_mito_percent_discard <- ifelse(spe$subsets_mito_percent_discard, "DISCARD", "KEEP")
spe$detected_discard <- ifelse(spe$detected_discard, "DISCARD", "KEEP")
spe$sum_discard <- ifelse(spe$sum_discard, "DISCARD", "KEEP")


png(here(plot_dir, "cancer_spotplot_discarded_3mad_grid_mito.png"), width=15, height=15, res=300, units="in")
p_list <- vis_grid_clus(spe,
        clustervar="subsets_mito_percent_discard",
        spatial = FALSE,
        return_plots = TRUE
    )

cowplot::plot_grid(plotlist = p_list, ncol = 3)
dev.off()

png(here(plot_dir, "cancer_spotplot_discarded_3mad_grid_detected.png"), width=15, height=15, res=300, units="in")
p_list <- vis_grid_clus(spe,
        clustervar="detected_discard",
        spatial = FALSE,
        return_plots = TRUE
    )

cowplot::plot_grid(plotlist = p_list, ncol = 3)
dev.off()

png(here(plot_dir, "cancer_spotplot_discarded_3mad_grid_sum.png"), width=15, height=15, res=300, units="in")
p_list <- vis_grid_clus(spe,
        clustervar="sum_discard",
        spatial = FALSE,
        return_plots = TRUE
    )

cowplot::plot_grid(plotlist = p_list, ncol = 3)
dev.off()




# run QC using arbitrary thresholds

spe$sum_discard_threshold <- spe$sum < 2000
spe$detected_discard_threshold  <- spe$detected < 1000
spe$subsets_mito_percent_discard_threshold  <- spe$subsets_mito_percent > 10

png(here(plot_dir, "cancer_spotplot_discarded_threshold.png"), width=12, height=5, res=300, units="in")
p1 <- plotQC(spe, metric="sum", outliers="sum_discard_threshold") +
  ggtitle("Sum UMI")

p2 <- plotQC(spe, metric="detected", outliers="detected_discard_threshold") +
  ggtitle("Genes Detected")

p3 <- plotQC(spe, metric="subsets_mito_percent", outliers="subsets_mito_percent_discard_threshold") +
  ggtitle("Mito Percent")

p1+p2+p3
dev.off()


# ====== SpotSweeper =====

# run SpotSweeper

spe <- localOutliers(spe, metric = "sum", direction = "lower", log = TRUE)
spe <- localOutliers(spe, metric = "detected", direction = "lower", log = TRUE)
spe <- localOutliers(spe, metric = "subsets_mito_percent", direction = "higher", log = TRUE)

# violins
png(here(plot_dir, "cancer_violins_discarded_spotsweeper.png"), width=9, height=3, res=300, units="in")
p1 <- plotColData(spe, x="sample_id", y="sum_z", color="sum_outliers") +
  ggtitle("Sum UMI")

p2 <- plotColData(spe, x="sample_id", y="detected_z", color="detected_outliers") +
  ggtitle("Genes Detected")

p3 <- plotColData(spe, x="sample_id", y="subsets_mito_percent_z", color="subsets_mito_percent_outliers") +
  ggtitle("Mito Percent")

p1+p2+p3
dev.off()

# spot plots
png(here(plot_dir, "cancer_spotplot_discarded_spotsweeper.png"), width=12, height=5, res=300, units="in")
p1 <- plotQC(spe, metric="sum", outliers="sum_outliers") +
  ggtitle("Sum UMI")

p2 <- plotQC(spe, metric="detected", outliers="detected_outliers") +
  ggtitle("Genes Detected")

p3 <- plotQC(spe, metric="subsets_mito_percent", outliers="subsets_mito_percent_outliers") +
  ggtitle("Mito Percent")

p1+p2+p3
dev.off()



# ======== Venn Diagram ==========

# counts
discard_list <- list(
  Threshold = which(spe$sum_discard_threshold),
  ThreeMAD = which(spe$sum_discard),
  SpotSweeper = which(spe$sum_outliers)
)

png(here(plot_dir, "ovariancancer_vennDiagram_counts_discard.png"), width=5, height=5, units="in", res=300)
ggVennDiagram(discard_list) +
  scale_fill_gradient(low="white",high = "red") +
  ggtitle("Sum Counts")
dev.off()


# genes
discard_list <- list(
  Threshold = which(spe$detected_discard_threshold),
  ThreeMAD = which(spe$detected_discard),
  SpotSweeper = which(spe$detected_outliers)
)

png(here(plot_dir, "ovariancancer_vennDiagram_detected_discard.png"), width=5, height=5, units="in", res=300)
ggVennDiagram(discard_list) +
  scale_fill_gradient(low="white",high = "red") +
  ggtitle("Detected Genes")
dev.off()


# mito ratio
discard_list <- list(
  Threshold = which(spe$subsets_mito_percent_discard_threshold),
  ThreeMAD = which(spe$subsets_mito_percent_discard),
  SpotSweeper = which(spe$subsets_mito_percent_outliers)
)

png(here(plot_dir, "ovariancancer_vennDiagram_mito_discard.png"), width=5, height=5, units="in", res=300)
ggVennDiagram(discard_list) +
  scale_fill_gradient(low="white",high = "red") +
  ggtitle("Mito Ratio")
dev.off()






