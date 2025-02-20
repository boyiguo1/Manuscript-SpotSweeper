library(SpotSweeper)
library(here)
library(ggspavis)
library(escheR)
library(patchwork)
library(dplyr)
library(tidyr)
library(scran)
library(scater)

processed_dir = here('processed-data')

spe <- readRDS(here("processed-data", "Xenium", "Xenium_spotsweeper.rds"))
spe
# class: SpatialExperiment 
# dim: 541 118650 
# metadata(0):
# assays(3): '' counts logcounts
# rownames(541): ABCC11 ACTA2 ... BLANK_0497 BLANK_0499
# rowData names(6): ID Symbol ... vars cv2
# colnames: NULL
# colData names(26): Sample Barcode ... nCounts_outliers_pca
#   nGenes_outliers_pca
# reducedDimNames(1): PCA
# mainExpName: NULL
# altExpNames(0):
# spatialCoords names(2) : xcoord ycoord
# imgData names(1): sample_id


# ========= Spatial outlier detection for Sum UMI =======
spe <- localOutliers(spe, metric="nCounts", direction="lower", log=TRUE ,n_neighbors=10)
spe$nCounts_outliers_k10 <- spe$nCounts_outliers
spe$nCounts_outliers_z_k10 <- spe$nCounts_z

spe <- localOutliers(spe, metric="nCounts", direction="lower", log=TRUE ,n_neighbors=36)
spe$nCounts_outliers_k36 <- spe$nCounts_outliers
spe$nCounts_outliers_z_k36 <- spe$nCounts_z

spe <- localOutliers(spe, metric="nCounts", direction="lower", log=TRUE ,n_neighbors=50)
spe$nCounts_outliers_k50 <- spe$nCounts_outliers
spe$nCounts_outliers_z_k50 <- spe$nCounts_z

spe <- localOutliers(spe, metric="nCounts", direction="lower", log=TRUE ,n_neighbors=75)
spe$nCounts_outliers_k75 <- spe$nCounts_outliers
spe$nCounts_outliers_z_k75 <- spe$nCounts_z

spe <- localOutliers(spe, metric="nCounts", direction="lower", log=TRUE ,n_neighbors=100)
spe$nCounts_outliers_k100 <- spe$nCounts_outliers
spe$nCounts_outliers_z_k100 <- spe$nCounts_z




# ========= Ggplot box plots of average outliers per sample over different K =========

# Color pallete
color_palette <- c("k10" = "#E41A1C", "k36" = "#377EB8", "k50" = "#4DAF4A", "k75" = "#984EA3", "k100" = "#FF7F00")

# Calculate the number of outliers per sample_id for each k
outlier_counts <- colData(spe) %>%
  as.data.frame() %>%
  group_by(sample_id) %>%
  summarise(across(starts_with("nCounts_outliers_k"), ~ sum(.x, na.rm = TRUE), .names = "total_{.col}"))


# Convert to long format and simplify K names
outlier_counts_long <- outlier_counts %>%
  pivot_longer(
    cols = starts_with("total_nCounts_outliers_k"),
    names_to = "K",
    values_to = "total_outliers"
  ) %>%
  mutate(K = sub("total_nCounts_outliers_", "", K))  # Simplify K names to just the values

# Reorder K as a factor
outlier_counts_long <- outlier_counts_long %>%
  mutate(K = factor(K, levels = c("k10", "k36", "k50", "k75", "k100")))

png(here("plots", "parameter_selection", "knn_selection_xenium", "barplot.png"), unit="in", width=3.5, height=3, res=300)
ggplot(outlier_counts_long, aes(x = K, y = total_outliers, fill=K)) +
  geom_col() + # Boxplot without outlier shapes
  geom_text(aes(label = total_outliers), vjust = -0.5, size = 2.5) +
  labs(
    title = "Number of Outliers Per Sample",
    x = "K-nearest neighbors",
    y = "Total Outliers"
  ) +
  ggpubr::theme_pubr() +
  #ylim 0
  #  ylim(0, 40) +
  theme(legend.position = "none") +
  scale_fill_manual(values = color_palette) 
dev.off()





# ============== Quantifying change in outliers ==============
# Extract outliers per sample
outliers_per_k <- colData(spe) %>%
  as.data.frame() %>%
  select(starts_with("nCounts_outliers_k")) %>%
  mutate(spot = rownames(.), sample_id = colData(spe)$sample_id)  # Add spot and sample info

# Get reference outliers for K = 36
reference_outliers <- outliers_per_k %>%
  filter(nCounts_outliers_k50 == 1) %>%
  pull(spot)  # Reference spots for K = 36

# Calculate shared and additional spots per sample
comparison_per_sample <- outliers_per_k %>%
  pivot_longer(
    cols = starts_with("nCounts_outliers_k"),
    names_to = "K",
    values_to = "is_outlier"
  ) %>%
  mutate(K = sub("nCounts_outliers_", "", K)) %>%  # Simplify K values
  group_by(sample_id, K) %>%
  summarise(
    total_spots = sum(is_outlier == 1, na.rm = TRUE),
    shared_spots = sum(spot %in% reference_outliers & is_outlier == 1, na.rm = TRUE),
    additional_spots = sum(!spot %in% reference_outliers & is_outlier == 1, na.rm = TRUE)
  ) %>%
  mutate(
    percent_shared = 100 * shared_spots / total_spots  # Percentage shared
  ) %>%
  filter(K != "k50")  # Exclude K = 36 from the comparison

# Reorder K as a factor
comparison_per_sample <- comparison_per_sample %>%
  mutate(K = factor(K, levels = c("k10", "k36", "k75", "k100")))

# box plot of shared outliers
png(here("plots", "parameter_selection", "knn_selection_xenium", "overlap.png"), unit="in", width=3, height=3, res=300)
ggplot(comparison_per_sample, aes(x = K, y = percent_shared, fill=K)) +
  geom_col() + # Boxplot without outlier shapes
  geom_text(aes(label = sprintf("%.1f%%", percent_shared)), vjust = -0.5, size = 2.5) +
  labs(
    title = "Shared Outliers w/ K=50",
    x = "K-nearest neighbors",
    y = "Percentage of Shared Outliers"
  ) +
  ggpubr::theme_pubr() +
  ylim(0, 100) +
  scale_fill_manual(values = color_palette) +
  theme(legend.position = "none") 
dev.off()

# box plot of additional outliers
png(here("plots", "parameter_selection", "knn_selection_xenium", "overlap_additional.png"), unit="in", width=3, height=3, res=300)
ggplot(comparison_per_sample, aes(x = K, y = additional_spots, fill=K)) +
  geom_col() + # Boxplot without outlier shapes
  geom_text(aes(label =  additional_spots), vjust = -0.5, size = 2.5) +
  labs(
    title = "Additional Outliers",
    x = "K-nearest neighbors",
    y = "Number of Additional Spots"
  ) +
  ggpubr::theme_pubr() +
  #ylim(0, 31) +
  scale_fill_manual(values = color_palette)  +
  theme(legend.position = "none") 
dev.off()





# ==== Cumulative distribution of Sum UMI z-thresholds ====

# Extract z-scores for each K
z_scores <- colData(spe) %>%
  as.data.frame() %>%
  select(starts_with("nCounts_outliers_z")) %>%
  mutate(spot = rownames(.), sample_id = colData(spe)$sample_id)  # Add spot and sample info

# subset z-scores
z_scores <- z_scores %>%
  filter(sample_id == unique(sample_id)[1])

# Reshape to long format
z_scores_long <- z_scores %>%
  pivot_longer(
    cols = starts_with("nCounts_outliers_z_k"),
    names_to = "K",
    values_to = "z_score"
  ) %>%
  mutate(K = sub("nCounts_outliers_z_", "", K))  # Simplify K names

# Reorder K as a factor
z_scores_long <- z_scores_long %>%
  mutate(K = factor(K, levels = c("k10", "k36", "k50", "k75", "k100")))

# Plot the cumulative distribution
png(here("plots", "parameter_selection", "knn_selection_xenium", "cumulative_distribution.png"), unit="in", width=4, height=3, res=300)
ggplot(z_scores_long, aes(x = z_score, color = K)) +
  stat_ecdf(geom = "step") +
  labs(
    title = "UMI Z-scores",
    x = "Z-score",
    y = "Cumulative Probability"
  ) +
  ggpubr::theme_pubr() +
  scale_color_brewer(palette = "Set1") +
  xlim(-4, 4) +
  theme(legend.position = "right") +
  guides(color = guide_legend(title = NULL)) +
  geom_vline(xintercept = c(-3), linetype = "dashed", color = "black")
dev.off()


# zoom in to -3 to highlight potential differences in the distributions
png(here("plots", "parameter_selection", "knn_selection_xenium", "cumulative_distribution_zoom.png"), unit="in", width=2.5, height=2.5, res=300)
ggplot(z_scores_long, aes(x = z_score, color = K)) +
  stat_ecdf(geom = "step") +
  labs(
    title = NULL,
    x = NULL,
    y = NULL
  ) +
  ggpubr::theme_pubr() +
  scale_color_brewer(palette = "Set1") +
  coord_cartesian(xlim = c(-4, -2), ylim =c(0,0.03)) +
  theme(legend.position = "none") +
  guides(color = guide_legend(title = NULL)) +
  geom_vline(xintercept = c(-3), linetype = "dashed", color = "black")
dev.off()





# ===============================================================
#  Cumulative distributions for detected genes and mito ratio 
# ===============================================================

# ======== Detected Genes ========
spe <- localOutliers(spe, metric="nGenes", direction="lower", log=TRUE ,n_neighbors=10)
spe$nGenes_outliers_k10 <- spe$nGenes_outliers
spe$nGenes_z_k10 <- spe$nGenes_z

spe <- localOutliers(spe, metric="nGenes", direction="lower", log=TRUE ,n_neighbors=36)
spe$nGenes_outliers_k36 <- spe$nGenes_outliers
spe$nGenes_z_k36 <- spe$nGenes_z

spe <- localOutliers(spe, metric="nGenes", direction="lower", log=TRUE ,n_neighbors=50)
spe$nGenes_outliers_k50 <- spe$nGenes_outliers
spe$nGenes_z_k50 <- spe$nGenes_z

spe <- localOutliers(spe, metric="nGenes", direction="lower", log=TRUE ,n_neighbors=75)
spe$nGenes_outliers_k75 <- spe$nGenes_outliers
spe$nGenes_z_k75 <- spe$nGenes_z

spe <- localOutliers(spe, metric="nGenes", direction="lower", log=TRUE ,n_neighbors=100)
spe$nGenes_outliers_k100 <- spe$nGenes_outliers
spe$nGenes_z_k100 <- spe$nGenes_z

# Extract z-scores for each K
z_scores <- colData(spe) %>%
  as.data.frame() %>%
  select(starts_with("nGenes_z")) %>%
  mutate(spot = rownames(.), sample_id = colData(spe)$sample_id)  # Add spot and sample info

  # subset z-scores
z_scores <- z_scores %>%
  filter(sample_id == unique(sample_id)[1])

# Reshape to long format
z_scores_long <- z_scores %>%
  pivot_longer(
    cols = starts_with("nGenes_z_k"),
    names_to = "K",
    values_to = "z_score"
  ) %>%
  mutate(K = sub("nGenes_z_", "", K))  # Simplify K values

# Reorder K as a factor
z_scores_long <- z_scores_long %>%
  mutate(K = factor(K, levels = c("k10", "k36", "k50", "k75", "k100")))

# Plot the cumulative distribution
png(here("plots", "parameter_selection", "knn_selection_xenium", "cumulative_distribution_genes.png"), unit="in", width=4, height=3, res=300)
ggplot(z_scores_long, aes(x = z_score, color = K)) +
  stat_ecdf(geom = "step") +
  labs(
    title = "Detected Genes Z-scores",
    x = "Z-score",
    y = "Cumulative Probability"
  ) +
  ggpubr::theme_pubr() +
  scale_color_brewer(palette = "Set1") +
  xlim(-4, 4) +
  theme(legend.position = "right") + 
  guides(color = guide_legend(title = NULL)) +
  geom_vline(xintercept = c(-3), linetype = "dashed", color = "black")
dev.off()

# zoom
png(here("plots", "parameter_selection", "knn_selection_xenium", "cumulative_distribution_genes_zoom.png"), unit="in", width=2.5, height=2.5, res=300)
ggplot(z_scores_long, aes(x = z_score, color = K)) +
  stat_ecdf(geom = "step") +
  labs(
    title = NULL,
    x = NULL,
    y = NULL
  ) +
  ggpubr::theme_pubr() +
  scale_color_brewer(palette = "Set1") +
  coord_cartesian(xlim = c(-4, -2), ylim =c(0,0.03)) +
  theme(legend.position = "none") +
  guides(color = guide_legend(title = NULL)) +
  geom_vline(xintercept = c(-3), linetype = "dashed", color = "black")
dev.off()




# ==== Sum UMI ====

# extract nCounts for all outliers across Ks
nCounts_outliers <- colData(spe) %>%
  as.data.frame() %>%
  select(starts_with("nCounts_outliers_k")) %>%
  mutate(spot = rownames(.), sample_id = colData(spe)$sample_id)  # Add spot and sample info

# add nCounts
nCounts_outliers$nCounts <- colData(spe)$nCounts

# Reshape to long format
nCounts_outliers_long <- nCounts_outliers %>%
  pivot_longer(
    cols = starts_with("nCounts_outliers_k"),
    names_to = "K",
    values_to = "is_outlier"
  ) %>%
  mutate(K = sub("nCounts_outliers_", "", K))  # Simplify K names

# keep only outliers
nCounts_outliers_long <- nCounts_outliers_long %>%
  filter(is_outlier == 1)

# Reorder K as a factor
nCounts_outliers_long <- nCounts_outliers_long %>%
  mutate(K = factor(K, levels = c("k10", "k36", "k50", "k75", "k100")))

# Plot the violin plot
png(here("plots", "parameter_selection", "knn_selection_xenium", "violin_nCounts.png"), unit="in", width=3, height=3, res=300)
ggplot(nCounts_outliers_long, aes(x = K, y = nCounts, fill=K)) +
  geom_violin(scale = "width", alpha = 0.7) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.5, color="black") +
  labs(
    title = "Sum UMI Outliers",
    x = "K-nearest neighbors",
    y = "Sum UMI"
  ) +
  ggpubr::theme_pubr() +
  scale_fill_manual(values = color_palette) +
  theme(legend.position = "none")
dev.off()


# ==== Detected Genes ====
# extract outliers for all outliers across Ks
nGenes_outliers <- colData(spe) %>%
  as.data.frame() %>%
  select(starts_with("nGenes_outliers_k")) %>%
  mutate(spot = rownames(.), sample_id = colData(spe)$sample_id)  # Add spot and sample info

# add nGenes
nGenes_outliers$nGenes <- colData(spe)$nGenes

# Reshape to long format
nGenes_outliers_long <- nGenes_outliers %>%
  pivot_longer(
    cols = starts_with("nGenes_outliers_k"),
    names_to = "K",
    values_to = "is_outlier"
  ) %>%
  mutate(K = sub("nGenes_outliers_", "", K))  # Simplify K names

# keep only outliers
nGenes_outliers_long <- nGenes_outliers_long %>%
  filter(is_outlier == 1)

# Reorder K as a factor
nGenes_outliers_long <- nGenes_outliers_long %>%
  mutate(K = factor(K, levels = c("k10", "k36", "k50", "k75", "k100")))

# Plot the violin plot
png(here("plots", "parameter_selection", "knn_selection_xenium", "violin_nGenes.png"), unit="in", width=3, height=3, res=300)
ggplot(nGenes_outliers_long, aes(x = K, y = nGenes, fill=K)) +
  geom_violin(scale = "width", alpha = 0.7) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.5, color="black") +
  labs(
    title = "Gene Outliers",
    x = "K-nearest neighbors",
    y = "Detected Genes"
  ) +
  ggpubr::theme_pubr() +
  scale_fill_manual(values = color_palette) +
  theme(legend.position = "none")
dev.off()
