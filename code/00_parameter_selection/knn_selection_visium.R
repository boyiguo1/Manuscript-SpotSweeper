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

## Download the spot-level data
spe <- spatialLIBD::fetch_data(type = "spe")
spe


# ========= Spatial outlier detection for Sum UMI =======
spe <- localOutliers(spe, metric="sum_umi", direction="lower", log=TRUE ,n_neighbors=6)
spe$sum_umi_outliers_k6 <- spe$sum_umi_outliers
spe$sum_umi_outliers_z_k6 <- spe$sum_umi_z

spe <- localOutliers(spe, metric="sum_umi", direction="lower", log=TRUE ,n_neighbors=18)
spe$sum_umi_outliers_k18 <- spe$sum_umi_outliers
spe$sum_umi_outliers_z_k18 <- spe$sum_umi_z

spe <- localOutliers(spe, metric="sum_umi", direction="lower", log=TRUE ,n_neighbors=36)
spe$sum_umi_outliers_k36 <- spe$sum_umi_outliers
spe$sum_umi_outliers_z_k36 <- spe$sum_umi_z

spe <- localOutliers(spe, metric="sum_umi", direction="lower", log=TRUE ,n_neighbors=60)
spe$sum_umi_outliers_k60 <- spe$sum_umi_outliers
spe$sum_umi_outliers_z_k60 <- spe$sum_umi_z

spe <- localOutliers(spe, metric="sum_umi", direction="lower", log=TRUE ,n_neighbors=90)
spe$sum_umi_outliers_k90 <- spe$sum_umi_outliers
spe$sum_umi_outliers_z_k90 <- spe$sum_umi_z

colnames(colData(spe))

table(spe$sum_umi_outliers_k6)
# FALSE  TRUE 
# 47292   389 

table(spe$lsum_umi_outliers_k18)
# FALSE  TRUE 
# 47539   142 

table(spe$sum_umi_outliers_k36)
# FALSE  TRUE 
# 47567   114 

table(spe$sum_umi_outliers_k60)
# FALSE  TRUE 
# 47571   110 

table(spe$sum_umi_outliers_k90)
# FALSE  TRUE 
# 47564   117 


# ========= Ggplot box plots of average outliers per sample over different K =========

# Color pallete
color_palette <- c("k6" = "#E41A1C", "k18" = "#377EB8", "k36" = "#4DAF4A", "k60" = "#984EA3", "k90" = "#FF7F00")

# Calculate the number of outliers per sample_id for each k
outlier_counts <- colData(spe) %>%
  as.data.frame() %>%
  group_by(sample_id) %>%
  summarise(across(starts_with("sum_umi_outliers_k"), ~ sum(.x, na.rm = TRUE), .names = "total_{.col}"))


# Convert to long format and simplify K names
outlier_counts_long <- outlier_counts %>%
  pivot_longer(
    cols = starts_with("total_sum_umi_outliers_k"),
    names_to = "K",
    values_to = "total_outliers"
  ) %>%
  mutate(K = sub("total_sum_umi_outliers_", "", K))  # Simplify K names to just the values

# Reorder K as a factor
outlier_counts_long <- outlier_counts_long %>%
  mutate(K = factor(K, levels = c("k6", "k18", "k36", "k60", "k90")))

png(here("plots", "parameter_selection", "knn_selection", "boxplot.png"), unit="in", width=3.5, height=3, res=300)
ggplot(outlier_counts_long, aes(x = K, y = total_outliers, fill=K)) +
  geom_boxplot(outlier.shape = NA,  alpha = 0.7) + # Boxplot without outlier shapes
  geom_jitter(width = 0.2, size = 2, alpha = 0.8, color = "black") + # Individual data points
  labs(
    title = "Number of Outliers Per Sample",
    x = "K-nearest neighbors",
    y = "Total Outliers"
  ) +
  ggpubr::theme_pubr() +
  #ylim 0
    ylim(0, 40) +
  theme(legend.position = "none") +
  scale_fill_manual(values = color_palette) 
dev.off()


# ============== Quantifying change in outliers ==============
# Extract outliers per sample
outliers_per_k <- colData(spe) %>%
  as.data.frame() %>%
  select(starts_with("sum_umi_outliers_k")) %>%
  mutate(spot = rownames(.), sample_id = colData(spe)$sample_id)  # Add spot and sample info

# Get reference outliers for K = 36
reference_outliers <- outliers_per_k %>%
  filter(sum_umi_outliers_k36 == 1) %>%
  pull(spot)  # Reference spots for K = 36

# Calculate shared and additional spots per sample
comparison_per_sample <- outliers_per_k %>%
  pivot_longer(
    cols = starts_with("sum_umi_outliers_k"),
    names_to = "K",
    values_to = "is_outlier"
  ) %>%
  mutate(K = sub("sum_umi_outliers_", "", K)) %>%  # Simplify K values
  group_by(sample_id, K) %>%
  summarise(
    total_spots = sum(is_outlier == 1, na.rm = TRUE),
    shared_spots = sum(spot %in% reference_outliers & is_outlier == 1, na.rm = TRUE),
    additional_spots = sum(!spot %in% reference_outliers & is_outlier == 1, na.rm = TRUE)
  ) %>%
  mutate(
    percent_shared = 100 * shared_spots / total_spots  # Percentage shared
  ) %>%
  filter(K != "k36")  # Exclude K = 36 from the comparison

# Reorder K as a factor
comparison_per_sample <- comparison_per_sample %>%
  mutate(K = factor(K, levels = c("k6", "k18", "k60", "k90")))

# box plot of shared outliers
png(here("plots", "parameter_selection", "knn_selection", "overlap.png"), unit="in", width=3, height=3, res=300)
ggplot(comparison_per_sample, aes(x = K, y = percent_shared, fill=K)) +
  geom_boxplot( alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7, color = "black") +
  labs(
    title = "Shared Outliers w/ K=36",
    x = "K-nearest neighbors",
    y = "Percentage of Shared Outliers"
  ) +
  ggpubr::theme_pubr() +
  ylim(0, 100) +
  scale_fill_manual(values = color_palette) +
  theme(legend.position = "none") 
dev.off()

# box plot of additional outliers
png(here("plots", "parameter_selection", "knn_selection", "overlap_additional.png"), unit="in", width=3, height=3, res=300)
ggplot(comparison_per_sample, aes(x = K, y = additional_spots, fill=K)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7, color="black") +
  labs(
    title = "Additional Outliers to K=36",
    x = "K-nearest neighbors",
    y = "Number of Additional Spots"
  ) +
  ggpubr::theme_pubr() +
  ylim(0, 31) +
  scale_fill_manual(values = color_palette)  +
  theme(legend.position = "none") 
dev.off()



# =========== Spot plots ==========

# subset to first sample, unique(spe$sample_id)[1]
spe_sub <- spe[, spe$sample_id == unique(spe$sample_id)[1]]
png(here("plots", "parameter_selection", "knn_selection", "spotplots.png"), unit="in", width=10, height=3, res=300)
p1 <- plotSpotQC(spe_sub, plot_type = "spot", in_tissue = "in_tissue", 
                 annotate = "sum_umi_outliers_k6", point_size = 0.2) +
                guides(color = "none") +
                ggtitle("K = 6 (1st order)")

p2 <- plotSpotQC(spe_sub, plot_type = "spot", in_tissue = "in_tissue", 
                 annotate = "sum_umi_outliers_k18", point_size = 0.2) +
                guides(color = "none") +
                ggtitle("K = 18 (2nd order)")

p3 <- plotSpotQC(spe_sub, plot_type = "spot", in_tissue = "in_tissue", 
                 annotate = "sum_umi_outliers_k36", point_size = 0.2) +
                guides(color = "none") +
                ggtitle("K = 36 (3rd order)")

p4 <- plotSpotQC(spe_sub, plot_type = "spot", in_tissue = "in_tissue", 
                 annotate = "sum_umi_outliers_k60", point_size = 0.2) +
                guides(color = "none") +
                ggtitle("K = 60 (4th order)")

p5 <- plotSpotQC(spe_sub, plot_type = "spot", in_tissue = "in_tissue", 
                 annotate = "sum_umi_outliers_k90", point_size = 0.2) +
                guides(color = "none") +
                ggtitle("K = 90 (5th order)")

p1 | p2 | p3 | p4 | p5
dev.off()




# ==== Cumulative distribution of Sum UMI z-thresholds ====

# Extract z-scores for each K
z_scores <- colData(spe) %>%
  as.data.frame() %>%
  select(starts_with("sum_umi_outliers_z")) %>%
  mutate(spot = rownames(.), sample_id = colData(spe)$sample_id)  # Add spot and sample info

# subset z-scores
z_scores <- z_scores %>%
  filter(sample_id == unique(sample_id)[1])

# Reshape to long format
z_scores_long <- z_scores %>%
  pivot_longer(
    cols = starts_with("sum_umi_outliers_z_k"),
    names_to = "K",
    values_to = "z_score"
  ) %>%
  mutate(K = sub("sum_umi_outliers_z_", "", K))  # Simplify K names

# Reorder K as a factor
z_scores_long <- z_scores_long %>%
  mutate(K = factor(K, levels = c("k6", "k18", "k36", "k60", "k90")))

# Plot the cumulative distribution
png(here("plots", "parameter_selection", "knn_selection", "cumulative_distribution.png"), unit="in", width=4, height=3, res=300)
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
png(here("plots", "parameter_selection", "knn_selection", "cumulative_distribution_zoom.png"), unit="in", width=2.5, height=2.5, res=300)
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
spe <- localOutliers(spe, metric="sum_gene", direction="lower", log=TRUE ,n_neighbors=6)
spe$sum_gene_outliers_k6 <- spe$sum_gene_outliers
spe$sum_gene_z_k6 <- spe$sum_gene_z

spe <- localOutliers(spe, metric="sum_gene", direction="lower", log=TRUE ,n_neighbors=18)
spe$sum_gene_outliers_k18 <- spe$sum_gene_outliers
spe$sum_gene_z_k18 <- spe$sum_gene_z

spe <- localOutliers(spe, metric="sum_gene", direction="lower", log=TRUE ,n_neighbors=36)
spe$sum_gene_outliers_k36 <- spe$sum_gene_outliers
spe$sum_gene_z_k36 <- spe$sum_gene_z

spe <- localOutliers(spe, metric="sum_gene", direction="lower", log=TRUE ,n_neighbors=60)
spe$sum_gene_outliers_k60 <- spe$sum_gene_outliers
spe$sum_gene_z_k60 <- spe$sum_gene_z

spe <- localOutliers(spe, metric="sum_gene", direction="lower", log=TRUE ,n_neighbors=90)
spe$sum_gene_outliers_k90 <- spe$sum_gene_outliers
spe$sum_gene_z_k90 <- spe$sum_gene_z

# Extract z-scores for each K
z_scores <- colData(spe) %>%
  as.data.frame() %>%
  select(starts_with("sum_gene_z")) %>%
  mutate(spot = rownames(.), sample_id = colData(spe)$sample_id)  # Add spot and sample info

  # subset z-scores
z_scores <- z_scores %>%
  filter(sample_id == unique(sample_id)[1])

# Reshape to long format
z_scores_long <- z_scores %>%
  pivot_longer(
    cols = starts_with("sum_gene_z_k"),
    names_to = "K",
    values_to = "z_score"
  ) %>%
  mutate(K = sub("sum_gene_z_", "", K))  # Simplify K values

# Reorder K as a factor
z_scores_long <- z_scores_long %>%
  mutate(K = factor(K, levels = c("k6", "k18", "k36", "k60", "k90")))

# Plot the cumulative distribution
png(here("plots", "parameter_selection", "knn_selection", "cumulative_distribution_genes.png"), unit="in", width=4, height=3, res=300)
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
png(here("plots", "parameter_selection", "knn_selection", "cumulative_distribution_genes_zoom.png"), unit="in", width=2.5, height=2.5, res=300)
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




# ========= mito ratio ============
spe <- localOutliers(spe, metric="expr_chrM_ratio", direction="higher", log=TRUE ,n_neighbors=6)
spe$expr_chrM_ratio_outliers_k6 <- spe$expr_chrM_ratio_outliers
spe$expr_chrM_ratio_z_k6 <- spe$expr_chrM_ratio_z

spe <- localOutliers(spe, metric="expr_chrM_ratio", direction="higher", log=TRUE ,n_neighbors=18)
spe$expr_chrM_ratio_outliers_k18 <- spe$expr_chrM_ratio_outliers
spe$expr_chrM_ratio_z_k18 <- spe$expr_chrM_ratio_z

spe <- localOutliers(spe, metric="expr_chrM_ratio", direction="higher", log=TRUE ,n_neighbors=36)
spe$expr_chrM_ratio_outliers_k36 <- spe$expr_chrM_ratio_outliers
spe$expr_chrM_ratio_z_k36 <- spe$expr_chrM_ratio_z

spe <- localOutliers(spe, metric="expr_chrM_ratio", direction="higher", log=TRUE ,n_neighbors=60)
spe$expr_chrM_ratio_outliers_k60 <- spe$expr_chrM_ratio_outliers
spe$expr_chrM_ratio_z_k60 <- spe$expr_chrM_ratio_z

spe <- localOutliers(spe, metric="expr_chrM_ratio", direction="higher", log=TRUE ,n_neighbors=90)
spe$expr_chrM_ratio_outliers_k90 <- spe$expr_chrM_ratio_outliers
spe$expr_chrM_ratio_z_k90 <- spe$expr_chrM_ratio_z

# Extract z-scores for each K
z_scores <- colData(spe) %>%
  as.data.frame() %>%
  select(starts_with("expr_chrM_ratio_z")) %>%
  mutate(spot = rownames(.), sample_id = colData(spe)$sample_id)  # Add spot and sample info

# subset z-scores
z_scores <- z_scores %>%
  filter(sample_id == unique(sample_id)[1])

# Reshape to long format
z_scores_long <- z_scores %>%
  pivot_longer(
    cols = starts_with("expr_chrM_ratio_z_k"),
    names_to = "K",
    values_to = "z_score"
  ) %>%
  mutate(K = sub("expr_chrM_ratio_z_", "", K))  # Simplify K values

# Reorder K as a factor
z_scores_long <- z_scores_long %>%
  mutate(K = factor(K, levels = c("k6", "k18", "k36", "k60", "k90")))

# Plot the cumulative distribution
png(here("plots", "parameter_selection", "knn_selection", "cumulative_distribution_mito.png"), unit="in", width=4, height=3, res=300)
ggplot(z_scores_long, aes(x = z_score, color = K)) +
  stat_ecdf(geom = "step") +
  labs(
    title = "Mito Ratio Z-scores",
    x = "Z-score",
    y = "Cumulative Probability"
  ) +
  ggpubr::theme_pubr() +
  scale_color_brewer(palette = "Set1") +
  xlim(-4, 4) +
  theme(legend.position = "right") +
  guides(color = guide_legend(title = NULL)) +
  geom_vline(xintercept = c(3), linetype = "dashed", color = "black")
dev.off()

# zoom in to -3 to highlight potential differences in the distributions
png(here("plots", "parameter_selection", "knn_selection", "cumulative_distribution_mito_zoom.png"), unit="in", width=2.5, height=2.5, res=300)
ggplot(z_scores_long, aes(x = z_score, color = K)) +
  stat_ecdf(geom = "step") +
  labs(
    title = NULL,
    x = NULL,
    y = NULL
  ) +
  ggpubr::theme_pubr() +
  scale_color_brewer(palette = "Set1") +
  coord_cartesian(xlim = c(2, 4), ylim =c(0.97,1)) +
  theme(legend.position = "none") +
  guides(color = guide_legend(title = NULL)) +
  geom_vline(xintercept = c(3), linetype = "dashed", color = "black")
dev.off()





# ==================================================================
# Violins of QC metrics across K outliers
# ==================================================================

# ==== Sum UMI ====

# extract sum_umi for all outliers across Ks
sum_umi_outliers <- colData(spe) %>%
  as.data.frame() %>%
  select(starts_with("sum_umi_outliers_k")) %>%
  mutate(spot = rownames(.), sample_id = colData(spe)$sample_id)  # Add spot and sample info

# add sum_umi
sum_umi_outliers$sum_umi <- colData(spe)$sum_umi

# Reshape to long format
sum_umi_outliers_long <- sum_umi_outliers %>%
  pivot_longer(
    cols = starts_with("sum_umi_outliers_k"),
    names_to = "K",
    values_to = "is_outlier"
  ) %>%
  mutate(K = sub("sum_umi_outliers_", "", K))  # Simplify K names

# keep only outliers
sum_umi_outliers_long <- sum_umi_outliers_long %>%
  filter(is_outlier == 1)

# Reorder K as a factor
sum_umi_outliers_long <- sum_umi_outliers_long %>%
  mutate(K = factor(K, levels = c("k6", "k18", "k36", "k60", "k90")))

# Plot the violin plot
png(here("plots", "parameter_selection", "knn_selection", "violin_sum_umi.png"), unit="in", width=3, height=3, res=300)
ggplot(sum_umi_outliers_long, aes(x = K, y = sum_umi, fill=K)) +
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
sum_gene_outliers <- colData(spe) %>%
  as.data.frame() %>%
  select(starts_with("sum_gene_outliers_k")) %>%
  mutate(spot = rownames(.), sample_id = colData(spe)$sample_id)  # Add spot and sample info

# add sum_gene
sum_gene_outliers$sum_gene <- colData(spe)$sum_gene

# Reshape to long format
sum_gene_outliers_long <- sum_gene_outliers %>%
  pivot_longer(
    cols = starts_with("sum_gene_outliers_k"),
    names_to = "K",
    values_to = "is_outlier"
  ) %>%
  mutate(K = sub("sum_gene_outliers_", "", K))  # Simplify K names

# keep only outliers
sum_gene_outliers_long <- sum_gene_outliers_long %>%
  filter(is_outlier == 1)

# Reorder K as a factor
sum_gene_outliers_long <- sum_gene_outliers_long %>%
  mutate(K = factor(K, levels = c("k6", "k18", "k36", "k60", "k90")))

# Plot the violin plot
png(here("plots", "parameter_selection", "knn_selection", "violin_sum_gene.png"), unit="in", width=3, height=3, res=300)
ggplot(sum_gene_outliers_long, aes(x = K, y = sum_gene, fill=K)) +
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

# ==== Mito ratio ====
# extract outliers for all outliers across Ks
expr_chrM_ratio_outliers <- colData(spe) %>%
  as.data.frame() %>%
  select(starts_with("expr_chrM_ratio_outliers_k")) %>%
  mutate(spot = rownames(.), sample_id = colData(spe)$sample_id)  # Add spot and sample info

# add expr_chrM_ratio
expr_chrM_ratio_outliers$expr_chrM_ratio <- colData(spe)$expr_chrM_ratio

# Reshape to long format
expr_chrM_ratio_outliers_long <- expr_chrM_ratio_outliers %>%
  pivot_longer(
    cols = starts_with("expr_chrM_ratio_outliers_k"),
    names_to = "K",
    values_to = "is_outlier"
  ) %>%
  mutate(K = sub("expr_chrM_ratio_outliers_", "", K))  # Simplify K names

# keep only outliers
expr_chrM_ratio_outliers_long <- expr_chrM_ratio_outliers_long %>%
  filter(is_outlier == 1)

# Reorder K as a factor
expr_chrM_ratio_outliers_long <- expr_chrM_ratio_outliers_long %>%
  mutate(K = factor(K, levels = c("k6", "k18", "k36", "k60", "k90")))

# Plot the violin plot
png(here("plots", "parameter_selection", "knn_selection", "violin_mito.png"), unit="in", width=3, height=3, res=300)
ggplot(expr_chrM_ratio_outliers_long, aes(x = K, y = expr_chrM_ratio, fill=K)) +
  geom_violin(scale = "width", alpha = 0.7) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.5, color="black") +
  labs(
    title = "Mito Ratio Outliers",
    x = "K-nearest neighbors",
    y = "Mito Ratio"
  ) +
  ggpubr::theme_pubr() +
  scale_fill_manual(values = color_palette) +
  theme(legend.position = "none")
dev.off()
