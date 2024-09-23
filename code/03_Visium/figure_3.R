library(SpotSweeper)
library(RANN)
library(here)
library(SpatialExperiment)
library(spatialLIBD)
library(ggspavis)
library(BiocNeighbors)
library(escheR)
library(patchwork)
library(tictoc)
library(dplyr)
library(tidyr)
library(scran)
library(scater)
library(ggridges)
library(escheR)
library(STexampleData)

#library(ExperimentHub)

#hub <- ExperimentHub()

## Select rows in the table
#spe.cancer <- hub[['EH6695']]



plot_dir = here('plots',"figure_3")
processed_dir = here('processed-data')

# example mouse data from 10x
#load(here("raw-data", "spe_raw_mouse.Rdata"))
#spe.mouse <- spe
#colnames(colData(spe.mouse))
# [1] "sample_id"              "in_tissue"              "array_row"              "array_col"
# [5] "10x_graphclust"         "10x_kmeans_10_clusters" "10x_kmeans_2_clusters"  "10x_kmeans_3_clusters"
# [9] "10x_kmeans_4_clusters"  "10x_kmeans_5_clusters"  "10x_kmeans_6_clusters"  "10x_kmeans_7_clusters"
# [13] "10x_kmeans_8_clusters"  "10x_kmeans_9_clusters"  "key"                    "sum_umi"
# [17] "sum_gene"               "expr_chrM"              "expr_chrM_ratio"        "ManualAnnotation"

spe.mouse <- Visium_mouseCoronal()
spe.mouse
# class: SpatialExperiment
# dim: 32285 4992
# metadata(0):
#   assays(1): counts
# rownames(32285): ENSMUSG00000051951 ENSMUSG00000089699 ... ENSMUSG00000095019 ENSMUSG00000095041
# rowData names(3): gene_id gene_name feature_type
# colnames(4992): AAACAACGAATAGTTC-1 AAACAAGTATCTCCCA-1 ... TTGTTTGTATTACACG-1 TTGTTTGTGTAAATTC-1
# colData names(5): barcode_id sample_id in_tissue array_row array_col
# reducedDimNames(0):
#   mainExpName: NULL
# altExpNames(0):
#   spatialCoords names(2) : pxl_col_in_fullres pxl_row_in_fullres
# imgData names(4): sample_id image_id data scaleFactor


# n=30 DLFPC
load(here("raw-data","dlPFC_raw.RData"))
spe.pfc <- spe_raw
colnames(colData(spe.pfc))
# [1] "sample_id"              "in_tissue"              "array_row"              "array_col"
# [5] "10x_graphclust"         "10x_kmeans_10_clusters" "10x_kmeans_2_clusters"  "10x_kmeans_3_clusters"
# [9] "10x_kmeans_4_clusters"  "10x_kmeans_5_clusters"  "10x_kmeans_6_clusters"  "10x_kmeans_7_clusters"
# [13] "10x_kmeans_8_clusters"  "10x_kmeans_9_clusters"  "key"                    "sum_umi"
# [17] "sum_gene"               "expr_chrM"              "expr_chrM_ratio"        "ManualAnnotation"
# [21] "subject"                "region"                 "sex"                    "age"
# [25] "diagnosis"              "sample_id_complete"     "count"

# Maynard et al
spe <- fetch_data(type = "spe")
spe.maynard <- spe
colnames(colData(spe.maynard))
# [1] "sample_id"                   "Cluster"                     "sum_umi"
# [4] "sum_gene"                    "subject"                     "position"
# [7] "replicate"                   "subject_position"            "discard"
# [10] "key"                         "cell_count"                  "SNN_k50_k4"
# [13] "SNN_k50_k5"                  "SNN_k50_k6"                  "SNN_k50_k7"
# [16] "SNN_k50_k8"                  "SNN_k50_k9"                  "SNN_k50_k10"
# [19] "SNN_k50_k11"                 "SNN_k50_k12"                 "SNN_k50_k13"
# [22] "SNN_k50_k14"                 "SNN_k50_k15"                 "SNN_k50_k16"
# [25] "SNN_k50_k17"                 "SNN_k50_k18"                 "SNN_k50_k19"
# [28] "SNN_k50_k20"                 "SNN_k50_k21"                 "SNN_k50_k22"
# [31] "SNN_k50_k23"                 "SNN_k50_k24"                 "SNN_k50_k25"
# [34] "SNN_k50_k26"                 "SNN_k50_k27"                 "SNN_k50_k28"
# [37] "GraphBased"                  "Maynard"                     "Martinowich"
# [40] "layer_guess"                 "layer_guess_reordered"       "layer_guess_reordered_short"
# [43] "expr_chrM"                   "expr_chrM_ratio"             "SpatialDE_PCA"
# [46] "SpatialDE_pool_PCA"          "HVG_PCA"                     "pseudobulk_PCA"
# [49] "markers_PCA"                 "SpatialDE_UMAP"              "SpatialDE_pool_UMAP"
# [52] "HVG_UMAP"                    "pseudobulk_UMAP"             "markers_UMAP"
# [55] "SpatialDE_PCA_spatial"       "SpatialDE_pool_PCA_spatial"  "HVG_PCA_spatial"
# [58] "pseudobulk_PCA_spatial"      "markers_PCA_spatial"         "SpatialDE_UMAP_spatial"
# [61] "SpatialDE_pool_UMAP_spatial" "HVG_UMAP_spatial"            "pseudobulk_UMAP_spatial"
# [64] "markers_UMAP_spatial"        "spatialLIBD"                 "ManualAnnotation"
# [67] "in_tissue"                   "array_row"                   "array_col"


# drop out of tissue spots
spe.pfc <- spe.pfc[,colData(spe.pfc)$in_tissue == TRUE]
spe.mouse <- spe.mouse[,colData(spe.mouse)$in_tissue == TRUE]
spe.maynard <- spe.maynard[,colData(spe.maynard)$in_tissue == TRUE]


# ======== Local Outlier detection =========

# ===== Find

# === PFC ===
spe.pfc <- localOutliers(spe.pfc,
                         metric="sum_umi",
                         direction="both"
)

spe.pfc <- localOutliers(spe.pfc,
                         metric="sum_gene",
                         direction="both"
)

spe.pfc <- localOutliers(spe.pfc,
                         metric="expr_chrM_ratio",
                         direction="higher",
                         log=FALSE
)

spe.pfc$local_outliers <- as.logical(spe.pfc$sum_umi_outliers) |
  as.logical(spe.pfc$sum_gene_outliers) |
  as.logical(spe.pfc$expr_chrM_ratio_outliers)

# === Maynard ===

spe.maynard <- localOutliers(spe.maynard,
                             metric="sum_umi",
                             direction="both"
)

spe.maynard <- localOutliers(spe.maynard,
                             metric="sum_gene",
                             direction="both"
)

spe.maynard <- localOutliers(spe.maynard,
                             metric="expr_chrM_ratio",
                             direction="higher",
                             log=FALSE
)

spe.maynard$local_outliers <- as.logical(spe.maynard$sum_umi_outliers) |
  as.logical(spe.maynard$sum_gene_outliers) |
  as.logical(spe.maynard$expr_chrM_ratio_outliers)

# === Mouse ===
rownames(spe.mouse) <- rowData(spe.mouse)$gene_name
is.mito <- grepl("^mt-", rownames(counts(spe.mouse)))

spe.mouse <- addPerCellQCMetrics(spe.mouse, subsets=list(Mito=is.mito))
spe.mouse$sum_umi <- spe.mouse$sum
spe.mouse$sum_gene <- spe.mouse$detected
spe.mouse$expr_chrM_ratio <- spe.mouse$subsets_Mito_percent/100

spe.mouse <- localOutliers(spe.mouse,
                           metric="sum_umi",
                           direction="both"
)

spe.mouse <- localOutliers(spe.mouse,
                           metric="sum_gene",
                           direction="both"
)

spe.mouse <- localOutliers(spe.mouse,
                           metric="expr_chrM_ratio",
                           direction="higher",
                           log=FALSE
)

spe.mouse$local_outliers <- as.logical(spe.mouse$sum_umi_outliers) |
  as.logical(spe.mouse$sum_gene_outliers) |
  as.logical(spe.mouse$expr_chrM_ratio_outliers)


# ======== Check for reoccurring outliers in DLFPC datasets ========

Z <- 3

spe.pfc$sum_umi_high <- spe.pfc$sum_umi_z > Z
spe.pfc$sum_umi_low <- spe.pfc$sum_umi_z < -Z
spe.mouse$sum_umi_low <- spe.mouse$sum_umi_z < -Z

spe.maynard$sum_umi_high <- spe.maynard$sum_umi_z > Z
spe.maynard$sum_umi_low <- spe.maynard$sum_umi_z < -Z
spe.mouse$sum_umi_low <- spe.mouse$sum_umi_z < -Z

# drop any NAs
spe.pfc$sum_umi_high[is.na(spe.pfc$sum_umi_high)] <- FALSE
spe.pfc$sum_umi_low[is.na(spe.pfc$sum_umi_low)] <- FALSE
spe.mouse$sum_umi_low[is.na(spe.mouse$sum_umi_low)] <- FALSE

spe.maynard$sum_umi_high[is.na(spe.maynard$sum_umi_high)] <- FALSE
spe.maynard$sum_umi_low[is.na(spe.maynard$sum_umi_low)] <- FALSE
spe.mouse$sum_umi_low[is.na(spe.mouse$sum_umi_low)] <- FALSE

# subset to only outliers
pfc.outliers_high <- spe.pfc[,spe.pfc$sum_umi_high]
pfc.outliers_low <- spe.pfc[,spe.pfc$sum_umi_low]

maynard.outliers_high <- spe.maynard[,spe.maynard$sum_umi_high]
maynard.outliers_low <- spe.maynard[,spe.maynard$sum_umi_low]

# create dataframes
pfc_outliers_high <- data.frame(sample_id = pfc.outliers_high$sample_id,
                           row = pfc.outliers_high$array_row,
                           col = pfc.outliers_high$array_col,
                           group="high outliers",
                           expt="DLFPC")
pfc_outliers_low <- data.frame(sample_id = pfc.outliers_low$sample_id,
                           row = pfc.outliers_low$array_row,
                           col = pfc.outliers_low$array_col,
                           group="low outliers",
                           expt="DLFPC")

maynard_outliers_high <- data.frame(sample_id = maynard.outliers_high$sample_id,
                           row = maynard.outliers_high$array_row,
                           col = maynard.outliers_high$array_col,
                           group="high outliers",
                           expt="Maynard")
maynard_outliers_low <- data.frame(sample_id = maynard.outliers_low$sample_id,
                           row = maynard.outliers_low$array_row,
                           col = maynard.outliers_low$array_col,
                           group="low outliers",
                           expt="Maynard")

# join the all dataframes
all_outliers <- rbind(pfc_outliers_high, pfc_outliers_low,
                      maynard_outliers_high, maynard_outliers_low)


# get the number of samples each spot is an outlier in, per group
spots_in_multiple_samples <- all_outliers %>%
  group_by(group, row, col) %>%  # Group by 'group', 'row', and 'col'
  summarise(samples = n_distinct(sample_id)) %>%
  arrange(desc(samples))

# plot boxplot with jitter by group
p <- ggplot(spots_in_multiple_samples, aes(x = group, y = samples)) +
  geom_jitter(width = 0.2, size=3) +
  geom_boxplot(outlier.shape = NA) +
  labs(y = "Total number of samples",
       x=NULL) +
  ggtitle("Outliers in multiple samples") +
  theme_bw() +
  theme(text=element_text(size=18),
        axis.text=element_text(size=17))

png(here(plot_dir, "Repeat_outliers_boxplot_3z.png"), width=5, height=5.5, units="in", res=300)
p
dev.off()


# get row and col for spots with counts > 20 samples
biased_spots <- spots_in_multiple_samples[spots_in_multiple_samples$samples > 20,]
biased_spots
# # A tibble: 6 × 4
# # Groups:   row, col [6]
# row   col group        count
# <int> <int> <chr>        <int>
# 1 low (z < -3)    57    99      42
# 2 low (z < -3)    38    88      39
# 3 low (z < -3)    47    77      39
# 4 low (z < -3)    49   103      35
# 5 low (z < -3)    28    32      29
# 6 low (z < -3)    53   117      25


biased_spots <- biased_spots %>%
  mutate(combined = paste(row, col, sep = "-")) %>%
  group_by(combined) %>%
  filter(row_number() == 1) %>%
  ungroup()

spe.pfc$array_combined <- paste(spe.pfc$array_row, spe.pfc$array_col, sep = "-")
spe.maynard$array_combined <- paste(spe.maynard$array_row, spe.maynard$array_col, sep = "-")
spe.mouse$array_combined <- paste(spe.mouse$array_row, spe.mouse$array_col, sep = "-")

spe.pfc$biased_spots <- spe.pfc$array_combined %in% biased_spots$combined
spe.maynard$biased_spots <- spe.maynard$array_combined %in% biased_spots$combined
spe.mouse$biased_spots <- spe.mouse$array_combined %in% biased_spots$combined

# ================================================
# Getting barcodes of bad spots
# ================================================

# get barcodes by splitting keys
spe.pfc$barcode <- sapply(strsplit(spe.pfc$key, "-"), "[[", 1)
spe.mouse$barcode <- sapply(strsplit(spe.mouse$barcode_id, "-"), "[[", 1)
spe.maynard$barcode <- sapply(strsplit(spe.maynard$key, "-"), "[[", 1)
spe.maynard$barcode <- sapply(strsplit(spe.maynard$barcode, "_"), "[[", 2)

# spot plots of bad spots
png(here(plot_dir, "Spotplot_badSpots_Maynard.png"), width=4.5, height=5, units="in", res=300)
p1 <- SpotSweeper::plotQC(spe.maynard, metric="sum_umi_log",outliers="biased_spots",
                          point_size=1.65) +
  ggtitle("DLFPC (Maynard et al)") +
  guides(color = guide_legend(title = "Low quality")) +
  theme(text=element_text(size=15))
p1
dev.off()

png(here(plot_dir, "Spotplot_badSpots_PFC.png"), width=5, height=5, units="in", res=300)
p2 <- SpotSweeper::plotQC(spe.pfc, metric="sum_umi_log",outliers="biased_spots",
                          point_size=1.8) +
  ggtitle("DLPFC (Huuki-Myers et al)") +
  guides(color = guide_legend(title = "Low quality")) +
  theme(text=element_text(size=15))
p2
dev.off()



# drop duplicate colData from spe.mouse
colData(spe.mouse) <- colData(spe.mouse)[, !duplicated(colnames(colData(spe.mouse)))]


png(here(plot_dir, "Spotplot_badSpots_Mouse.png"), width=5, height=5, units="in", res=300)
p3 <- SpotSweeper::plotQC(spe.mouse, metric="sum_umi_log",outliers="biased_spots",
                          point_size=2.4) +
  ggtitle("Mouse Coronal (10x Genomics)") +
  guides(color = guide_legend(title = "Low quality")) +
  theme(text=element_text(size=18))
p3
dev.off()





# === export csv of bad spot barcodes and col and row ===
# subset to only bad spots
spe.temp <- spe.maynard[,spe.maynard$biased_spots]

biased_spots_df <- data.frame(row = spe.temp$array_row,
                           col = spe.temp$array_col,
                           barcode = spe.temp$barcode)

# export to csv
write.csv(unique(biased_spots_df), here(processed_dir, "barcodes_for_biased_spots.csv"), row.names = FALSE)

# ===============================
# Characterizing bad spots
# ===============================

# sum umi of bad spots
p1 <- plotColData(spe.maynard, x="biased_spots",y="sum_umi", color_by="sample_id") +
  ggtitle("Maynard et al.")+
  labs(y="Sum UMI",
       x="Consistent outliers")+
  theme(legend.position = "none",
        text=element_text(size=15),
        axis.text=element_text(size=15))

p2 <- plotColData(spe.pfc, x="biased_spots",y="sum_umi", color_by="sample_id") +
  ggtitle("Huuki-Myers et al.")+
  labs(y="Sum UMI",
       x="Consistent outliers")+
  theme(legend.position = "none",
        text=element_text(size=15),
        axis.text=element_text(size=15))

p3 <- plotColData(spe.mouse, x="biased_spots",y="sum_umi", color_by="sample_id") +
  ggtitle("Mouse Coronal (10x Genomics)")+
  labs(y="Sum UMI",
       x="Consistent outliers")+
  theme(legend.position = "none",
        text=element_text(size=15),
        axis.text=element_text(size=15)) +
  # rename


png(here(plot_dir, "bad_spots_sum_umi.png"), width=15, height=5, units="in", res=300)
(p1+p2+p3)
dev.off()


# sum gene of bad spots
p1 <- plotColData(spe.maynard, x="biased_spots",y="sum_gene", color_by="sample_id") +
  labs(y="Unique genes",
       x="Consistent outliers")+
  theme(legend.position = "none",
        text=element_text(size=15),
        axis.text=element_text(size=15))

p2 <- plotColData(spe.pfc, x="biased_spots",y="sum_gene", color_by="sample_id") +
  labs(y="Unique genes",
       x="Consistent outliers")+
  theme(legend.position = "none",
        text=element_text(size=15),
        axis.text=element_text(size=15))

p3 <- plotColData(spe.mouse, x="biased_spots",y="sum_gene", color_by="sample_id") +
  labs(y="Unique genes",
       x="Consistent outliers")+
  theme(legend.position = "none",
        text=element_text(size=15),
        axis.text=element_text(size=15))

png(here(plot_dir, "bad_spots_sum_gene.png"), width=15, height=5, units="in", res=300)
(p1+p2+p3)
dev.off()


# mito ratio of bad spots
p1 <- plotColData(spe.maynard, x="biased_spots",y="expr_chrM_ratio", color_by="sample_id") +
  labs(y="Mito ratio",
       x="Consistent outliers")+
  theme(legend.position = "none",
        text=element_text(size=15),
        axis.text=element_text(size=15))

p2 <- plotColData(spe.pfc, x="biased_spots",y="expr_chrM_ratio", color_by="sample_id") +
  labs(y="Mito ratio",
       x="Consistent outliers")+
  theme(legend.position = "none",
        text=element_text(size=15),
        axis.text=element_text(size=15))

p3 <- plotColData(spe.mouse, x="biased_spots",y="expr_chrM_ratio", color_by="sample_id") +
  labs(y="Mito ratio",
       x="Consistent outliers")+
  theme(legend.position = "none",
        text=element_text(size=15),
        axis.text=element_text(size=15))

png(here(plot_dir, "bad_spots_mito_ratio.png"), width=15, height=5, units="in", res=300)
(p1+p2+p3)
dev.off()


# z-score of bad spots
p1 <- plotColData(spe.maynard, x="biased_spots",y="sum_umi_z", color_by="sum_umi_low") +
  ggtitle("Maynard et al.")+
  labs(y="Local z-score",
       x="Consistent outliers") +
  theme(legend.position = "none",
        text=element_text(size=23),
        axis.text=element_text(size=20),
        plot.title = element_text(face="plain")) +
  scale_color_manual(values = c("grey", "red"))

p2 <- plotColData(spe.pfc, x="biased_spots",y="sum_umi_z", color_by="sum_umi_low") +
  ggtitle("Huuki-Myers et al.") +
  labs(y="Local z-score",
       x="Consistent outliers") +
  theme(legend.position = "none",
        text=element_text(size=23),
        axis.text=element_text(size=20),
        plot.title = element_text(face="plain"))+
  scale_color_manual(values = c("grey", "red"))


p3 <- plotColData(spe.mouse, x="biased_spots",y="sum_umi_z", color_by="sum_umi_low") +
  ggtitle("Mouse Coronal") +
  labs(y="Local z-score",
       x="Consistent outliers") +
  scale_color_manual(values = c("grey", "red")) +
  guides(color = guide_legend(title = "Outlier")) +
  theme(text=element_text(size=23),
      axis.text=element_text(size=20),
      plot.title = element_text(face="plain"))


png(here(plot_dir, "bad_spots_z.png"), width=15, height=5, units="in", res=300)
(p1+p2+p3)
dev.off()

# =======================================
# K-mer analysis of bad barcodes
# =======================================

# make dataframe with barcodes and bad spots
barcodes_df <- data.frame(biased_spots = spe.maynard$biased_spots, barcode = spe.maynard$barcode)

# replace TRUE with Bad, FALSE with Good
barcodes_df$biased_spots <- ifelse(barcodes_df$biased_spots, "Bad", "Good")


# ===== prepping data frames =====
# get dataframe with barcode sequence and read count
barcode_counts_mouse <- data.frame(read_count = spe.mouse$sum_umi, barcode = spe.mouse$barcode)
barcode_counts_pfc <- data.frame(read_count = spe.pfc$sum_umi, barcode = spe.pfc$barcode)
barcode_counts_maynard <- data.frame(read_count = spe.maynard$sum_umi, barcode = spe.maynard$barcode)

# concatenate dfs
barcode_counts <- rbind(#barcode_counts_mouse,
                        barcode_counts_pfc,
                        barcode_counts_maynard
                        )

# get barcodes for bad spots
bad_barcodes <- barcodes_df[barcodes_df$biased_spots == "Bad", "barcode"]
good_barcodes <- barcodes_df[barcodes_df$biased_spots == "Good", "barcode"]

# average read counts with standard deviation for each unique barcode with dplyr
barcode_counts_avg <- barcode_counts %>%
  group_by(barcode) %>%
  summarise(mean = mean(read_count))


# sort by highest to lowest
barcode_counts_avg <- barcode_counts_avg[order(barcode_counts_avg$mean, decreasing = TRUE),]

# add an integer for row number
barcode_counts_avg$row <- 1:nrow(barcode_counts_avg)


# ggplot lineplot with SEM shaded error bars
p1 <- ggplot(barcode_counts_avg, aes(x = row, y = mean)) +
  geom_line()  +
  labs(x = "Barcodes", y = "UMI count") +
  ggtitle("Average library size") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        text=element_text(size=20))


png(here(plot_dir, "barcode_read_count_ranked.png"), width=5, height=5, units="in", res=300)
p1
dev.off()

# get the top 10 and bottom 10 barcodes
top_barcodes <- barcode_counts_avg[1:6, "barcode"]
bottom_barcodes <- barcode_counts_avg[(nrow(barcode_counts_avg)-9):nrow(barcode_counts_avg), "barcode"]
random_barcodes <- barcode_counts_avg[sample(1:nrow(barcode_counts_avg), 1000), "barcode"]
bad_barcodes <- unique(biased_spots_df)$barcode

library(Biostrings)
library(msa)
library(ggmsa)


# bad spots
dss.bad <- DNAStringSet(unique(bad_barcodes))
test <- msa(dss.bad)
msaPrettyPrint(test, output="pdf", file=here("alignment_biased_spots.pdf"),askForOverwrite=FALSE)

# bottom spots
dss.bottom <- DNAStringSet(bottom_barcodes$barcode)
test <- msa(dss.bottom)
msaPrettyPrint(test, output="pdf", file="alignment_bottom.pdf",askForOverwrite=FALSE)

# top spots
dss.top <- DNAStringSet(top_barcodes$barcode)
test <- msa(dss.top)
msaPrettyPrint(test, output="pdf", file="alignment_top.pdf",askForOverwrite=FALSE)


# random spots
#dss.random <- DNAStringSet(random_barcodes$barcode)
#test <- msa(dss.random)
#msaPrettyPrint(test, output="pdf", file="alignment_random.pdf",askForOverwrite=FALSE)



library(stringr)

calculate_gc_content <- function(sequence) {
  num_g <- str_count(sequence, "G")
  num_c <- str_count(sequence, "C")

  gc_content <- (num_g + num_c) / str_length(sequence) * 100
  return(gc_content)
}

bad_barcodes_gc <- calculate_gc_content(bad_barcodes)
top_barcodes_gc <- calculate_gc_content(pull(top_barcodes))
random_barcodes_gc <- calculate_gc_content(pull(random_barcodes))

bad_barcodes_gc
#[1] 62.50 43.75 68.75 62.50 43.75 50.00

top_barcodes_gc
# [1] 31.25 56.25 25.00 43.75 50.00 37.50

random_barcodes_gc
# [1] 43.75 62.50 50.00 37.50 43.75 50.00


# Combine data into a single dataframe
gc_content_data <- data.frame(
  GC_Content = c(bad_barcodes_gc, top_barcodes_gc),
  Barcode_Type = factor(c(rep("Bad", length(bad_barcodes_gc)),
                          rep("Top", length(top_barcodes_gc))
                          )
                        )
)

# Create the box plot with geom_points
ggplot(gc_content_data, aes(x = Barcode_Type, y = GC_Content)) +
  geom_boxplot(aes(fill=Barcode_Type)) +
  geom_point(position = position_jitter(width = 0.2, height = 0)) +
  labs(title = "GC Content Comparison of Barcodes", x = "Barcode Type", y = "GC Content (%)") +
  theme_bw() +
  theme(text=element_text(size=28),
        axis.text=element_text(size=28),
        plot.title = element_text(size=28)) +
  scale_color_manual(values=c("firebrick1", "dodgerblue"))






# make a function to get_GC_content
get_GC_content <- function(sequence,string){
  num_g <- str_count(sequence, string[1])
  gc_content <- (num_g / str_length(sequence)) * 100
  return(gc_content)
}

# apply over all barcodes
barcode_counts_avg$a_content <- sapply(barcode_counts_avg$barcode, get_GC_content, string = c("A"))
barcode_counts_avg$t_content <- sapply(barcode_counts_avg$barcode, get_GC_content, string = c("T"))
barcode_counts_avg$c_content <- sapply(barcode_counts_avg$barcode, get_GC_content, string = c("C"))
barcode_counts_avg$g_content <- sapply(barcode_counts_avg$barcode, get_GC_content, string = c("G"))

p1 <- ggplot(barcode_counts_avg, aes(x = a_content, y = mean)) +
  geom_point() +
  ggtitle("A content vs read count") +
  labs(x = "A content (%)", y = "Read count") +
  #trend line
  geom_smooth(method = "lm")

p2 <- ggplot(barcode_counts_avg, aes(x = t_content, y = mean)) +
  geom_point() +
  ggtitle("T content vs read count") +
  labs(x = "T content (%)", y = "Read count") +
  #trend line
  geom_smooth(method = "lm")

p3 <- ggplot(barcode_counts_avg, aes(x = c_content, y = mean)) +
  geom_point() +
  ggtitle("C content vs read count") +
  labs(x = "C content (%)", y = "Read count") +
  #trend line
  geom_smooth(method = "lm")

p4 <- ggplot(barcode_counts_avg, aes(x = g_content, y = mean)) +
  geom_point() +
  ggtitle("G content vs read count") +
  labs(x = "G content (%)", y = "Read count") +
  #trend line
  geom_smooth(method = "lm")

png(here(plot_dir, "atcg_correlations.png"), width=10, height=10, units="in", res=300)
(p1+p2)/(p3+p4)
dev.off()

# get A content in top 10 vs bottom 10
top_barcodes_a <- sapply(top_barcodes$barcode, get_GC_content, string = c("A"))
bottom_barcodes_a <- sapply(bottom_barcodes$barcode, get_GC_content, string = c("A"))
a_df <- data.frame(a_content = c(top_barcodes_a, bottom_barcodes_a),
                   group = c(rep("Top", length(top_barcodes_a)),
                             rep("Bottom", length(bottom_barcodes_a))))

# get T content in top 10 vs bottom 10
top_barcodes_t <- sapply(top_barcodes$barcode, get_GC_content, string = c("T"))
bottom_barcodes_t <- sapply(bottom_barcodes$barcode, get_GC_content, string = c("T"))
t_df <- data.frame(t_content = c(top_barcodes_t, bottom_barcodes_t),
                   group = c(rep("Top", length(top_barcodes_t)),
                             rep("Bottom", length(bottom_barcodes_t))))

# get C content in top 10 vs bottom 10
top_barcodes_c <- sapply(top_barcodes$barcode, get_GC_content, string = c("C"))
bottom_barcodes_c <- sapply(bottom_barcodes$barcode, get_GC_content, string = c("C"))
c_df <- data.frame(c_content = c(top_barcodes_c, bottom_barcodes_c),
                   group = c(rep("Top", length(top_barcodes_c)),
                             rep("Bottom", length(bottom_barcodes_c))))

# get G content in top 10 vs bottom 10
top_barcodes_g <- sapply(top_barcodes$barcode, get_GC_content, string = c("G"))
bottom_barcodes_g <- sapply(bottom_barcodes$barcode, get_GC_content, string = c("G"))
g_df <- data.frame(g_content = c(top_barcodes_g, bottom_barcodes_g),
                   group = c(rep("Top", length(top_barcodes_g)),
                             rep("Bottom", length(bottom_barcodes_g))))


# plot boxplot of top vs bottom
p1 <- ggplot(a_df, aes(x = group, y = a_content, fill=group)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge()) +
  ggtitle("A content") +
  labs(x = "Group", y = "A content (%)")


p2 <- ggplot(t_df, aes(x = group, y = t_content, fill=group)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge()) +
  ggtitle("T content") +
  labs(x = "Group", y = "T content (%)")


p3 <- ggplot(c_df, aes(x = group, y = c_content, fill=group)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge()) +
  ggtitle("C content") +
  labs(x = "Group", y = "C content (%)")

p4 <- ggplot(g_df, aes(x = group, y = g_content, fill=group)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge()) +
  ggtitle("G content") +
  labs(x = "Group", y = "G content (%)")

png(here(plot_dir, "atcg_boxplot.png"), width=10, height=10, units="in", res=300)
(p1+p2)/(p3+p4)
dev.off()



get_GC_sum <- function(sequence,string){
  num_g <- str_count(sequence, string[1])
  gc_content <- num_g
  return(gc_content)
}

# GT sequences
top_barcodes_gt <- sapply(top_barcodes$barcode, get_GC_sum, string = c("GT"))
bottom_barcodes_gt <- sapply(bottom_barcodes$barcode, get_GC_sum, string = c("GT"))
gt_df <- data.frame(gt_content = c(top_barcodes_gt, bottom_barcodes_gt),
                    group = c(rep("Top", length(top_barcodes_gt)),
                              rep("Bottom", length(bottom_barcodes_gt))))

ggplot(gt_df, aes(x = group, y = gt_content, fill=group)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge()) +
  ggtitle("GT count") +
  labs(x = "Group", y = "GT count")


top_barcodes_aa <- sapply(top_barcodes$barcode, get_GC_sum, string = c("AAA"))
bottom_barcodes_aa <- sapply(bottom_barcodes$barcode, get_GC_sum, string = c("AAA"))
aa_df <- data.frame(aa_content = c(top_barcodes_aa, bottom_barcodes_aa),
                    group = c(rep("Top", length(top_barcodes_aa)),
                              rep("Bottom", length(bottom_barcodes_aa))))

ggplot(aa_df, aes(x = group, y = aa_content, fill=group)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge()) +
  ggtitle("AAA count") +
  labs(x = "Group", y = "AAA count")





barcode_counts_avg$gt_content <- sapply(barcode_counts_avg$barcode, get_GC_sum, string = c("AAA"))


# ========= Multiple linear regression ===========
# let's see how much of this we can regress out. I'm betting.
ggplot(barcode_counts_avg, aes(x = gt_content, y = mean)) +
  geom_point() +
  ggtitle("A content vs read count") +
  labs(x = "A content (%)", y = "Read count") +
  #trend line
  geom_smooth(method = "lm")


# perform multiple linear regression of a, g, t, c
model <- lm(mean ~t_content + a_content + g_content, data = barcode_counts_avg)
summary(model)

plot(resid(model))


# Fit the logistic regression model
model <- glm(mean ~t_content + a_content + g_content + c_content, family = poisson(link="identity"), data = barcode_counts_avg)
summary(model)

plot(resid(model))


# add bad barcodes to barcode_counts_avg
barcode_counts_avg$biased_spots <- barcode_counts_avg$barcode %in% bad_barcodes

# ================  differential kmer analysis ==============
library(Biostrings)

compute_kmers <- function(dss, k) {
  oligonucleotideFrequency(dss, k)
}

max_k <- 4

dss.all <- DNAStringSet(unique(barcode_counts_avg$barcode))
kmer.all_list <- lapply(1:max_k, compute_kmers, dss = dss.all)
kmer.all_df <- do.call(cbind, kmer.all_list)
# drop 0s
kmer.all_df <- kmer.all_df[rowSums(kmer.all_df) != 0,]

sce <- SingleCellExperiment(assays = list(counts = t(kmer.all_df)))

sce$barcode <- barcode_counts_avg$barcode
sce$biased_spots <- barcode_counts_avg$biased_spots

# add top and bottom barcode coldata
top_barcodes <- barcode_counts_avg[1:6, "barcode"]
bottom_barcodes <- barcode_counts_avg[(nrow(barcode_counts_avg)-9):nrow(barcode_counts_avg), "barcode"]
bad_barcodes <- unique(biased_spots_df)$barcode

# make colData of Top, Bottom, NA barcodes
sce$top_barcodes <- NA

sce$top_barcodes[sce$barcode %in% top_barcodes$barcode] <- "Top"
sce$top_barcodes[sce$barcode %in% unique(biased_spots_df)$barcode] <- "Bad"
unique(sce$top_barcodes)
#[1] "Top"    NA       "Bottom"


# subset to just top and bottom
sce.subset <- sce[,!is.na(sce$top_barcodes)]


# differential expression of kmers
# library("DESeq2")
# library("EnhancedVolcano")
#
# dds <- DESeqDataSet(sce.subset, design = ~ top_barcodes)
# keep <- rowSums(counts(dds) > 0)
#
# dds <- DESeq(dds, test="LRT", reduced=~1)
# res <- results(dds)
# res
#
# png(here(plot_dir, "kmer_volcano_DEseq2.png"), width=5, height=7, units="in", res=300)
# EnhancedVolcano(res,
#                 lab = rownames(res),
#                 x = 'log2FoldChange',
#                 y = 'pvalue',
#                 pCutoff = 0.05)
# dev.off()

# DE using scran findMarkers - binomial
logcounts(sce.subset) <- counts(sce.subset)
markers <- scran::findMarkers(sce.subset, groups = sce.subset$top_barcodes, test.type="binom")
Bottom_markers <- markers[[1]]

png(here(plot_dir, "kmer_volcano_k4_binom.png"), width=5, height=7, units="in", res=300)
EnhancedVolcano(Bottom_markers, lab=rownames(Bottom_markers), x='logFC.Top', y='p.value',
                pCutoff = 0.05,
                FCcutoff = 1) +
  xlim(-3, 3) +
  ylim(0, 7)
dev.off()


# DE using scran findMarkers - t-test
markers <- scran::findMarkers(sce.subset, groups = sce.subset$top_barcodes, test.type="t")
Bottom_markers <- markers[[1]]
png(here(plot_dir, "kmer_volcano_k4_t.png"), width=5, height=7, units="in", res=300)
EnhancedVolcano(Bottom_markers, lab=rownames(Bottom_markers), x='logFC.Top', y='p.value',
                pCutoff = 0.05,
                FCcutoff = .5,
                title = 'K-mer differential expression',
                subtitle=NULL,
                pointSize = 4.0,
                labSize = 6.0,
                colAlpha = 1,
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75) +
  xlim(-3, 3) +
  ylim(0, 4) +
  theme(text=element_text(size=20),
        plot.title=element_text(face="plain"))
dev.off()

# get markers < pval .05
Bottom_markers_sig <- Bottom_markers[Bottom_markers$p.value < .05,]

# ======== DE when subsetting for k-mer size =======
# subset to kmers length 1
sce.kmer1 <- sce[nchar(rownames(sce)) == 1,]
sce.kmer2 <- sce[nchar(rownames(sce)) == 2,]
sce.kmer3 <- sce[nchar(rownames(sce)) == 3,]
sce.kmer4 <- sce[nchar(rownames(sce)) == 4,]

logcounts(sce.kmer1) <- counts(sce.kmer1)
logcounts(sce.kmer2) <- counts(sce.kmer2)
logcounts(sce.kmer3) <- counts(sce.kmer3)
logcounts(sce.kmer4) <- counts(sce.kmer4)

# find marker
markers.kmer1 <- scran::findMarkers(sce.kmer1, groups = sce.kmer1$top_barcodes)
markers.kmer1[[1]]

markers.kmer2 <- scran::findMarkers(sce.kmer2, groups = sce.kmer2$top_barcodes)
markers.kmer2[[1]]

markers.kmer3 <- scran::findMarkers(sce.kmer3, groups = sce.kmer3$top_barcodes)
markers.kmer3[[1]]

markers.kmer4 <- scran::findMarkers(sce.kmer4, groups = sce.kmer4$top_barcodes)
markers.kmer4[[1]]

# volcano plots
png(here(plot_dir, "kmer1_volcano.png"), width=5, height=7, units="in", res=300)
EnhancedVolcano(markers.kmer1[[1]], lab=rownames(markers.kmer1[[1]]), x='logFC.Top', y='FDR',
                pCutoff=.01,
                FCcutoff = .5)+
  ylab("-Log10P (FDR)")
dev.off()

png(here(plot_dir, "kmer2_volcano.png"), width=5, height=7, units="in", res=300)
EnhancedVolcano(markers.kmer2[[1]], lab=rownames(markers.kmer2[[1]]), x='logFC.Top', y='FDR',
                pCutoff=.01,
                FCcutoff = .5) +
  ylab("-Log10P (FDR)")
dev.off()

png(here(plot_dir, "kmer3_volcano.png"), width=5, height=7, units="in", res=300)
EnhancedVolcano(markers.kmer3[[1]], lab=rownames(markers.kmer3[[1]]), x='logFC.Top', y='FDR',
                pCutoff=.01,
                FCcutoff = .5)+
  ylab("-Log10P (FDR)")
dev.off()

png(here(plot_dir, "kmer4_volcano.png"), width=5, height=7, units="in", res=300)
EnhancedVolcano(markers.kmer4[[1]], lab=rownames(markers.kmer4[[1]]), x='logFC.Top', y='p.value',
                pCutoff=.01,
                FCcutoff = .5)+
  ylab("-Log10P (FDR)")
dev.off()


# ======= Heatmaps ======
# ComplexHeatmap
library("ComplexHeatmap")
library(RColorBrewer)

# subset to only Topbarcodes
sce.subset.sig <- sce.subset[rownames(Bottom_markers_sig),]
sce.subset.sig$barcodes <- sce.subset.sig$top_barcodes

png(here(plot_dir, "kmer_sig_heatmap.png"), width=7, height=4.5, units="in", res=300)
p <- plotHeatmap(sce.subset.sig, features=unique(rownames(sce.subset.sig)),
            #group="top_barcodes",
            assay="counts",
            center=TRUE,
            scale=TRUE,
            colour_columns_by="barcodes",
            color=colorRampPalette(c("blue", "white", "red"))(101)
)

# Draw the plot
grid.draw(p$gtable)

# Add color scale label
grid.text("Color Scale Label", x=0.92, y=0.5, rot=90, gp=gpar(fontsize=10))

dev.off()



library(ComplexHeatmap)
# Subset the SCE object to only include Topbarcodes
sce.subset.sig <- sce.subset[rownames(Bottom_markers_sig),]
sce.subset.sig$barcodes <- sce.subset.sig$top_barcodes
sce.subset.sig$barcodes
sce.subset.sig$barcodes[sce.subset.sig$barcodes == "Bad"] <- "Consistent_outlier"

# Extract counts matrix
counts_matrix <- assay(sce.subset.sig, "counts")

# Center and scale the data
centered_scaled_matrix <- t(scale(t(counts_matrix), center=TRUE, scale=TRUE))

# Create a color mapping for barcodes
unique_barcodes <- unique(sce.subset.sig$barcodes)
barcode_colors <- setNames(c("darkgreen","orange"), unique_barcodes)

# Create column annotations
col_annotation <- HeatmapAnnotation(
  barcodes = sce.subset.sig$barcodes,
  col = list(barcodes = barcode_colors)
)

# Generate the heatmap with ComplexHeatmap
png(here(plot_dir, "kmer_sig_heatmap.png"), width=7, height=4.5, units="in", res=300)

ht <- Heatmap(centered_scaled_matrix,
        name = "K-mer counts\n(centered and scaled)",
        col = colorRampPalette(c("blue", "white", "red"))(101),
        show_row_names = TRUE,
        show_column_names = FALSE,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        rect_gp = gpar(col = "darkgrey", lwd = 1),
        top_annotation = col_annotation)

draw(ht, merge_legend = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right")

dev.off()




# ========= Box plots comparing k-mer expression ============
# get all barcodes with 'AAA' or 'CGT'
AAA_barcodes <- barcode_counts_avg[grep("AAC", barcode_counts_avg$barcode),]
CGT_barcodes <- barcode_counts_avg[grep("GTGT", barcode_counts_avg$barcode),]

# combine with a new column "sequence"
AAA_barcodes$sequence <- "AAC"
CGT_barcodes$sequence <- "GTGT"
barcode_counts_avg$sequence <- "other"

# combine
combined_barcodes <- rbind(AAA_barcodes, CGT_barcodes, barcode_counts_avg)

# box plots
png(here(plot_dir, "AAC_GTGT_boxplot.png"), width=8, height=8, units="in", res=300)
ggplot(combined_barcodes, aes(x = sequence, y = mean, fill=sequence)) +
  geom_point(position = position_jitterdodge()) +
  geom_boxplot() +
  ggtitle("AAC vs GTGT library size") +
  labs(x = "Sequence", y = "UMI count") +
  theme_bw() +
  theme(text=element_text(size=28),
        axis.text=element_text(size=28),
        plot.title = element_text(size=28)) +
  scale_fill_manual(values=c("firebrick1",  "dodgerblue", "grey"))
dev.off()

png(here(plot_dir, "AAC_GTGT_violin.png"), width=8, height=8, units="in", res=300)
ggplot(combined_barcodes, aes(x = sequence, y = mean, fill=sequence)) +
  geom_violin(aes(fill = sequence)) +
  scale_fill_manual(values=c("firebrick1",  "dodgerblue", "grey")) +
  geom_boxplot(width = 0.2, fill = "white") +
  ggtitle("AAC vs GTGT library size") +
  labs(x = "Sequence", y = "UMI count") +
  theme_bw() +
  theme(text=element_text(size=28),
        axis.text=element_text(size=28),
        plot.title = element_text(size=28))
dev.off()




# figuring out which barcode goes to which spot
biased_spots
# # A tibble: 6 × 5
# group          row   col samples combined
# <chr>        <int> <int>   <int> <chr>
#   1 low (z < -3)    57    99      42 57-99
# 2 low (z < -3)    38    88      39 38-88
# 3 low (z < -3)    47    77      39 47-77
# 4 low (z < -3)    49   103      35 49-103
# 5 low (z < -3)    28    32      29 28-32
# 6 low (z < -3)    53   117      25 53-117

biased_spots$barcodes <- bad_barcodes

# plot row and col as x and y coordinates. scatter plot. color by combined

p <- ggplot(biased_spots, aes(x = col, y = row, color = barcodes)) +
  geom_point(size=3) +
  ggtitle("Consistent outliers") +
  labs(x = "Column", y = "Row") +
  theme_bw() +
  theme(text=element_text(size=28),
        axis.text=element_text(size=28),
        plot.title = element_text(size=28)) +
  scale_color_manual(values = c("red", "blue", "green", "purple", "orange", "black"))

png(here(plot_dir, "BadSpots_scatter_color_by_barcode.png"), width=10, height=6, units="in", res=300)
p
dev.off()



