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


spe <- readRDS(here("processed-data", "VisiumHD","VisiumHD_HumanLung_016.rds"))
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


colnames(colData(spe))
#  [1] "barcode"                      "in_tissue"                   
#  [3] "array_row"                    "array_col"                   
#  [5] "sample_id"                    "sum"                         
#  [7] "detected"                     "subsets_mito_percent"        
#  [9] "subsets_mito_sum"             "sum_discard"                 
# [11] "detected_discard"             "subsets_mito_percent_discard"




# ===== Spatially naive QC ======
# plot metrics
png(here(plot_dir, "HD_HumanLung_016_spotplot_metrics.png"), width=12, height=5, res=300, units="in")
p1 <- plotQC(spe, metric="detected", point_size=0.4)
p2 <- plotQC(spe, metric="sum", point_size=0.4)
p3 <- plotQC(spe, metric="subsets_mito_percent", point_size=0.4)

p1+p2+p3
dev.off()

png(here(plot_dir, "HD_HumanLung_016_spotplot_metrics_discard.png"), width=12, height=5, res=300, units="in")
p4 <- plotQC(spe,metric="detected", outliers="detected_discard", point_size=0.4, stroke=.5)
p5 <- plotQC(spe, metric="sum", outliers="sum_discard", point_size=0.4, stroke=.5)
p6 <- plotQC(spe,  metric="subsets_mito_percent", outliers="subsets_mito_percent_discard", point_size=0.4, stroke=.5)

p4+p5+p6
dev.off()

# scatter plot of mito and sum
png(here(plot_dir, "HD_HumanLung_016_scatter_sum_mito.png"), width=5, height=5, res=300, units="in")
plotColData(spe, x="sum", y="subsets_mito_percent")
dev.off()


# ======= SpotSweeper =========
spe <- localOutliers(spe, metric = "sum", direction = "lower", log = TRUE)
spe <- localOutliers(spe, metric = "detected", direction = "lower", log = TRUE)
spe <- localOutliers(spe, metric = "subsets_mito_percent", direction = "higher", log = TRUE)

# violins
png(here(plot_dir, "HD_HumanLung_016_violins_discarded_spotsweeper.png"), width=9, height=3, res=300, units="in")
p1 <- plotColData(spe, x="sample_id", y="sum_z", color="sum_outliers") +
  ggtitle("Sum UMI")

p2 <- plotColData(spe, x="sample_id", y="detected_z", color="detected_outliers") +
  ggtitle("Genes Detected")

p3 <- plotColData(spe, x="sample_id", y="subsets_mito_percent_z", color="subsets_mito_percent_outliers") +
  ggtitle("Mito Percent")

p1+p2+p3
dev.off()


# spot plots
png(here(plot_dir, "HD_HumanLung_016_spotplot_discarded_spotsweeper.png"), width=12, height=5, res=1200, units="in")
p1 <- plotQC(spe, metric="sum", outliers="sum_outliers", point_size=.3) +
  ggtitle("Sum UMI")

p2 <- plotQC(spe, metric="detected", outliers="detected_outliers", point_size=.3) +
  ggtitle("Genes Detected")

p3 <- plotQC(spe, metric="subsets_mito_percent", outliers="subsets_mito_percent_outliers", point_size=.3) +
  ggtitle("Mito Percent")

p1+p2+p3
dev.off()



# ====== Spot plot for main figure =======
png(here(plot_dir, "HD_HumanLung_016_spotplot_main_figure.png"), width=12, height=5, res=300, units="in")
p1 <- plotQC(spe, metric="detected", point_size=0.4)
p2 <- plotQC(spe, metric="detected", outliers="sum_discard", point_size=0.4, stroke=.5)
p3 <- plotQC(spe, metric="detected", outliers="sum_outliers", point_size=.4, stroke=.5)
p1+p2+p3
dev.off()
