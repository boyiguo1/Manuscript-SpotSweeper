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

plot_dir = here('plots',"figure_2")
processed_dir = here('processed-data')

## Download the spot-level data
spe <- fetch_data(type = "spe")
spe


# ========= Spatial outlier detection =======
spe <- localOutliers(spe,
                     metric="sum_umi",
                     direction="lower"
)

spe <- localOutliers(spe,
                     metric="sum_gene",
                     direction="lower"
)

spe <- localOutliers(spe,
                     metric="expr_chrM_ratio",
                     direction="higher",
                     log=FALSE
)

spe$local_outliers <- as.logical(spe$sum_umi_outliers) |
  as.logical(spe$sum_gene_outliers) |
  as.logical(spe$expr_chrM_ratio_outliers)



# ========= Standard outlier detection (MAD) =======
library(scuttle)

qc.umi_mad <- isOutlier(spe$sum_umi, log=TRUE, type="lower")
table(qc.umi_mad)
# FALSE  TRUE
# 47230   451

qc.gene_mad <- isOutlier(spe$sum_gene, log=TRUE, type="lower")
table(qc.gene_mad)
# FALSE  TRUE
# 47171   510

qc.mito_mad <- isOutlier(spe$expr_chrM_ratio, type="higher")
table(qc.mito_mad)
# qc.mito_mad
# FALSE  TRUE
# 47598    83

colData(spe)$qc_mito_mad <- qc.mito_mad
colData(spe)$qc_umi_mad <- qc.umi_mad
colData(spe)$qc_gene_mad <- qc.gene_mad

mad_outliers <- qc.umi_mad | qc.mito_mad | qc.gene_mad
spe$discard_mad <- mad_outliers
colnames(colData(spe))

# ======== Arbitrary threshold QC =======

qc_mito <- colData(spe)$expr_chrM_ratio > 0.275 # this is arbitrary
table(qc_mito)
# FALSE  TRUE
# 47401   280

qc_umi <- colData(spe)$sum_umi < 500 # this is arbitrary
table(qc_umi)
# FALSE  TRUE
# 47372   309

qc_gene <- colData(spe)$sum_gene < 500 # this is arbitrary
table(qc_gene)

colData(spe)$qc_mito_threshold <- qc_mito
colData(spe)$qc_umi_threshold <- qc_umi
colData(spe)$qc_gene_threshold <- qc_gene

threshold_outliers <- qc_umi | qc_mito | qc_gene
spe$discard_threshold <- threshold_outliers


# # ===== Average percentage discarded across samples =====
# Extract the required columns
discard_mad_col <- spe$discard_mad
discard_local_col <- spe$local_outliers
discard_threshold_col <- spe$discard_threshold
layer_guess_ordered_col <- spe$layer_guess_reordered
sample_id_col <- spe$sample_id
subject_col <- spe$subject

# Create a data frame
df <- data.frame(discard_mad = discard_mad_col,
                 discard_threshold = discard_threshold_col,
                 discard_local = as.logical(discard_local_col),
                 layer_guess_ordered = layer_guess_ordered_col,
                 sample_id = sample_id_col,
                 subject = subject_col)

df <- df[!is.na(df$layer_guess_ordered), ]

discarded_df <- df %>%
  group_by(sample_id, layer_guess_ordered, subject) %>%
  summarise(
    total_spots = n(),  # Total number of spots in each domain for each sample
    discarded_mad = sum(discard_mad, na.rm = TRUE),  # Counting TRUE as 1 and FALSE as 0
    discarded_local = sum(discard_local, na.rm = TRUE),  # Counting TRUE as 1 and FALSE as 0
    discarded_threshold = sum(discard_threshold, na.rm = TRUE),  # Counting TRUE as 1 and FALSE as 0
    .groups = 'keep'  # Dropping the grouping after summarise
  ) %>%
  mutate(
    percentage_mad = (discarded_mad / total_spots) * 100  # Calculating the percentage
  ) %>%
  mutate(
    percentage_local = (discarded_local/ total_spots) * 100  # Calculating the percentage
  )%>%
  mutate(
    percentage_threshold = (discarded_threshold/ total_spots) * 100  # Calculating the percentage
  )



head(discarded_df)
# # A tibble: 6 × 10
# # Groups:   sample_id, layer_guess_ordered, subject [6]
# sample_id layer_guess_ordered subject total_spots discarded_mad discarded_local
# <chr>     <fct>               <chr>         <int>         <int>           <int>
#   1 151507    Layer1              Br5292          817            14               1
# 2 151507    Layer2              Br5292          305             0               0
# 3 151507    Layer3              Br5292         1215            10               4
# 4 151507    Layer4              Br5292          369             9               1
# 5 151507    Layer5              Br5292          675             3               1
# 6 151507    Layer6              Br5292          486             0               0
# # ℹ 4 more variables: discarded_threshold <int>, percentage_mad <dbl>, percentage_local <dbl>,
# #   percentage_threshold <dbl>


# ===== Make df of QC metrics across layers =====

# Extract data
expr_data <- colData(spe)$expr_chrM
expr_ratio <- colData(spe)$expr_chrM_ratio
expr_ratio_z <- colData(spe)$expr_chrM_ratio_z
layer_data <- colData(spe)$layer_guess_reordered
cell_count <- colData(spe)$cell_count
sum_umi <- colData(spe)$sum_umi
sum_umi_z <- colData(spe)$sum_umi_z
sum_gene <- colData(spe)$sum_gene
sum_gene_z <- colData(spe)$sum_gene_z

# Combine the expression data with the layer data
qc_df <- data.frame(
  expr_chrM = expr_data,
  expr_ratio = expr_ratio,
  expr_ratio_z = expr_ratio_z,
  layer = layer_data,
  sum_umi = sum_umi,
  sum_umi_z = sum_umi_z,
  cell_count = cell_count,
  sum_gene = sum_gene,
  sum_gene_z = sum_gene_z,
  cell_count = cell_count,
  discard_local = as.logical(discard_local_col)
)


# drop "NA" layer
qc_df <- qc_df[!is.na(qc_df$layer),]

# reverse the y order so that layer 1 is at the top
qc_df$layer <- factor(qc_df$layer, levels = rev(levels(qc_df$layer)))




# ================== From here we can start to generate plots ==================

# subset sample
spe.subset <- subset(spe, ,sample_id == unique(spe$sample_id)[1])

# ===== Panel A =====
# visualizng the manual layer annotations and common QC metrics

p1 <- make_escheR(spe.subset) |>
  add_fill(var = "layer_guess_reordered", stroke = 1) +
  scale_fill_manual(
    name = "", # turn off legend name for ground_truth
    values = libd_layer_colors) +
  ggtitle("Manual Annotations") +
  theme(plot.title = element_text(size = 20),
        legend.text = element_text(size=12))

pdf(height = 5, width=5, here(plot_dir, 'Figure2_spotplot_layers.pdf'))
p1
dev.off()


# ===== Panels B + C =====


p2 <- make_escheR(spe.subset) |>
  add_fill(var = "expr_chrM_ratio") |>
  add_ground(var = "qc_mito_threshold", stroke = 1) +
  scale_color_manual(
    name = "", # turn off legend name for ground_truth
    values = c(
      "TRUE" = "red2",
      "FALSE" = "transparent")
  ) +
  scale_fill_gradient(low ="white",high =  "black") +
  labs(fill = "Mito ratio") +
  guides(color = guide_legend(title = "Discarded")) +
  ggtitle("Mitochondrial Ratio") +
  theme(plot.title = element_text(size = 20),
        legend.text = element_text(size=11))

pdf(height = 5, width=5, here(plot_dir, 'Figure2_spotplot_mito_global.pdf'))
p2
dev.off()

p3 <- make_escheR(spe.subset) |>
  add_fill(var = "sum_gene") |>
  add_ground(var = "qc_gene_threshold", stroke = 1) +
  scale_color_manual(
    name = "", # turn off legend name for ground_truth
    values = c(
      "TRUE" = "red2",
      "FALSE" = "transparent")
  ) +
  scale_fill_gradient(low ="white",high =  "black") +
  labs(fill = "Sum Genes") +
  guides(color = guide_legend(title = "Discarded")) +
  ggtitle("Global outliers") +
  theme(plot.title = element_text(size = 20),
        legend.text = element_text(size=11))

pdf(height = 5, width=5, here(plot_dir, 'Figure2_spotplot_gene_global.pdf'))
p3
dev.off()


p4 <- make_escheR(spe.subset) |>
  add_fill(var = "sum_umi") |>
  add_ground(var = "qc_umi_threshold", stroke = 1) +
  scale_color_manual(
    name = "", # turn off legend name for ground_truth
    values = c(
      "TRUE" = "red2",
      "FALSE" = "transparent")
  ) +
  scale_fill_gradient(low ="white",high =  "black") +
  labs(fill = "Sum UMI") +
  guides(color = guide_legend(title = "Discarded")) +
  ggtitle("Total UMI") +
  theme(plot.title = element_text(size = 20),
        legend.text = element_text(size=11))

pdf(height = 5, width=5, here(plot_dir, 'Figure2_spotplot_umi_global.pdf'))
p4
dev.off()

# ===== Panel D ======

p <- ggplot(discarded_df, aes(x = layer_guess_ordered, y = percentage_threshold, fill = layer_guess_ordered)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, aes(color = layer_guess_ordered), alpha = 0.7) +  # Display all data points with matching layer colors
  scale_fill_manual(values = spatialLIBD::libd_layer_colors) +
  scale_color_manual(values = spatialLIBD::libd_layer_colors) +
  scale_y_continuous(limits = c(0, 60)) +
  labs(title = "Global outliers",
       x = "Spatial Domain",
       y = "Percent discarded (%)") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size=12),
        plot.title = element_text(size = 20),
        text = element_text(size=15))

pdf(height = 5, width=5, here(plot_dir, 'Figure2_boxplot_global.pdf'))
p
dev.off()



# ===== Panel E-G =====

# ridge plot of sum_umi with a threshold of 500
p4 <- ggplot(qc_df, aes(x = sum_umi, y = layer, fill = layer)) +
  geom_density_ridges(alpha = 0.8, scale = 1.5) +
  geom_vline(xintercept = 500, linetype = "dashed", color = "red", size=1) +
  scale_x_continuous(trans='log10') +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(size = 20),
        text = element_text(size = 15)) +
  scale_fill_manual(values = libd_layer_colors) +
  labs(title = "Library size per layer",
       x = "Library Size",
       y = "Spatial Domain")

# ridge plot of expr_ratio with a threshold of 0.3
p5 <- ggplot(qc_df, aes(x = expr_ratio, y = layer, fill = layer)) +
  geom_density_ridges(alpha = 0.8, scale = 1.5) +
  geom_vline(xintercept = 0.275, linetype = "dashed", color = "red", size=1) +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(size = 20),
        text = element_text(size = 15)) +
  scale_fill_manual(values = libd_layer_colors) +
  labs(title = "Mito percent per layer",
       x = "Mitochondrial Ratio",
       y = "Spatial Domain")

# ridge plot of sum_gene with a threshold of 500
p6 <- ggplot(qc_df, aes(x = sum_gene, y = layer, fill = layer)) +
  geom_density_ridges(alpha = 0.8, scale = 1.5) +
  geom_vline(xintercept = 500, linetype = "dashed", color = "red", size=1) +
  theme_bw() +
  scale_x_continuous(trans='log10') +
  theme(legend.position = "none",
        plot.title = element_text(size = 20),
        text = element_text(size = 15)) +
  scale_fill_manual(values = libd_layer_colors) +
  labs(title = "Unique genes per layer",
       x = "Number of genes",
       y = "Spatial Domain")

pdf(height = 5, width=5, here(plot_dir, 'Figure2_ridge_umi.pdf'))
p4
dev.off()

pdf(height = 5, width=5, here(plot_dir, 'Figure2_ridge_mito.pdf'))
p5
dev.off()

pdf(height = 5, width=5, here(plot_dir, 'Figure2_ridge_gene.pdf'))
p6
dev.off()





# ===== Panel E ======

# p <- ggplot(discarded_df, aes(x = layer_guess_ordered, y = percentage_mad, fill = layer_guess_ordered)) +
#   geom_boxplot(alpha = 0.7, outlier.shape = NA) +
#   geom_jitter(width = 0.2, size = 2, aes(color = layer_guess_ordered), alpha = 0.7) +  # Display all data points with matching layer colors
#   scale_fill_manual(values = spatialLIBD::libd_layer_colors) +
#   scale_color_manual(values = spatialLIBD::libd_layer_colors) +
#   labs(title = "3 MAD",
#        x = "Spatial Domain",
#        y = "Percent discarded (%)") +
#   theme_bw(base_size=20) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12),
#         legend.position = "none")
#
# pdf(height = 5, width=5, here(plot_dir, 'Figure2_Panel_E.pdf'))
# p
# dev.off()




# ======== Visualizing outliers via UMAP ========
# set.seed(555)
#
# # change to gene names
# rownames(spe) <- rowData(spe)$gene_name
#
# # get top 5000 HVGs
# top.hvg <- getTopHVGs(spe, n=5000)
#
# #run PCA
# spe <- runPCA(spe, ncomponents=50, subset_row=top.hvg)
#
# # batch correction using fastMNN
# library(batchelor)
# mnn.output <- fastMNN(spe, batch = spe$sample_id)
#
# # copy corrected PCA
# reducedDims(spe)$PCA_mnn <-reducedDims(mnn.output)$corrected
#
# # run UMAP
# spe <- runUMAP(spe, dimred = "PCA_mnn")
#
# # plot UMAP
# p1 <- plotReducedDim(spe, "UMAP", colour_by = "layer_guess_reordered", point_size=0.2) +
#   scale_color_manual(values = libd_layer_colors) +
#   labs(title = "Manual Annotations") +
#   theme(plot.title = element_text(size = 20))
#
# p2 <- plotReducedDim(spe, "UMAP", colour_by = "discard_mad", point_size=0.2)+
#   scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
#   labs(title = "Global Outliers (3 MADs)") +
#   theme(plot.title = element_text(size = 20))
#
# p3 <- plotReducedDim(spe, "UMAP", colour_by = "local_outliers", point_size=0.2) +
#   scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
#   labs(title = "Local Outliers (3 Zs)") +
#   theme(plot.title = element_text(size = 20))
#
# p4 <- plotReducedDim(spe, "UMAP", colour_by = "expr_chrM_ratio", point_size=0.2) +
#   labs(title = "Mitochondrial Ratio") +
#   theme(plot.title = element_text(size = 20))
#
# p5 <- plotReducedDim(spe, "UMAP", colour_by = "sum_umi", point_size=0.2) +
#   labs(title = "Library Size") +
#   theme(plot.title = element_text(size = 20))
#
# p6 <- plotReducedDim(spe, "UMAP", colour_by = "sum_gene", point_size=0.2) +
#   labs(title = "Unique Genes") +
#   theme(plot.title = element_text(size = 20))
#
# pdf(height = 7.5, width=15, here(plot_dir, 'Figure2_UMAPs.pdf'))
# (p1+p2+p3)/(p4+p5+p6)
# dev.off()



# Get the percent of local outliers that are also considered global outliers
local_outliers <- as.logical(as.character(spe$local_outliers))
discard_threshold <-  as.logical(as.character(spe$discard_threshold))
discard_mad <-  as.logical(as.character(spe$discard_mad))

# First, find the observations that are local outliers
local_outliers_indices <- which(local_outliers)

# Then, check which of these are also mad outliers
mad_outliers_among_local <- discard_mad[local_outliers_indices]
threshold_outliers_among_local <- discard_threshold[local_outliers_indices]

# Calculate the percentage
percentage_mad_among_local <- sum(mad_outliers_among_local) / length(mad_outliers_among_local) * 100
percentage_threshold_among_local <- sum(threshold_outliers_among_local) / length(threshold_outliers_among_local) * 100

# visualize using bar plot. Use nice colors
p <- ggplot(data.frame(x = c("MAD", "Threshold"), y = c(percentage_mad_among_local, percentage_threshold_among_local)), aes(x, y)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(title = "Local Outliers that are also Global Outliers",
       x = "Outlier Type",
       y = "Percent") +
  theme_bw(base_size=20) +
  theme(plot.title = element_text(size = 20))
p

# NOTES: Roughly 60% of current local outliers are also global outliers.

# ========= Visualing outliers like miQC =========
#library(miQC)


p1 <- plotColData(spe, x="sum_gene", y="expr_chrM_ratio", color_by="discard_threshold", point_size=1) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
  theme_bw() +
  labs(title = "Global Outliers",
       x = "Number of genes",
       y = "Mitochondrial Ratio") +
  theme(plot.title = element_text(size = 20),
        text = element_text(size = 15)) +
  guides(color = guide_legend(title = "Discarded"))

png(height = 5, width=7, units="in",res=300, here(plot_dir, 'Figure2_gene_vs_mito.png'))
p1
dev.off()




# ============= Figure 1 - Part 2 =============


p1 <- make_escheR(spe.subset) |>
  add_fill(var = "sum_gene") |>
  add_ground(var = "sum_gene_outliers", stroke = 1) +
  scale_color_manual(
    name = "", # turn off legend name for ground_truth
    values = c(
      "TRUE" = "red2",
      "FALSE" = "transparent")
  ) +
  scale_fill_gradient(low ="white",high =  "black") +
  labs(fill = "Sum Genes") +
  guides(color = guide_legend(title = "Discarded")) +
  ggtitle("Local outliers") +
  theme(plot.title = element_text(size = 20),
        legend.text = element_text(size=11))

pdf(height = 5, width=5, here(plot_dir, 'Figure2_spotplot_gene_outliers.pdf'))
p1
dev.off()

p2 <- make_escheR(spe.subset) |>
  add_fill(var = "expr_chrM_ratio") |>
  add_ground(var = "expr_chrM_ratio_outliers", stroke = 1) +
  scale_color_manual(
    name = "", # turn off legend name for ground_truth
    values = c(
      "TRUE" = "red2",
      "FALSE" = "transparent")
  ) +
  scale_fill_gradient(low ="white",high =  "black") +
  labs(fill = "Mito Ratio)") +
  guides(color = guide_legend(title = "Discarded")) +
  ggtitle("Mitochondrial Ratio") +
  theme(plot.title = element_text(size = 20),
        legend.text = element_text(size=11))

pdf(height = 5, width=5, here(plot_dir, 'Figure2_spotplot_mito_outliers.pdf'))
p2
dev.off()

p3 <- make_escheR(spe.subset) |>
  add_fill(var = "sum_gene_log") |>
  add_ground(var = "sum_gene_outliers", stroke = 1) +
  scale_color_manual(
    name = "", # turn off legend name for ground_truth
    values = c(
      "TRUE" = "red2",
      "FALSE" = "transparent")
  ) +
  scale_fill_gradient(low ="white",high =  "black") +
  labs(fill = "log1p(Sum Genes)") +
  guides(color = guide_legend(title = "Discarded")) +
  ggtitle("Unique Genes") +
  theme(plot.title = element_text(size = 20),
        legend.text = element_text(size=11))

pdf(height = 5, width=5, here(plot_dir, 'Figure2_spotplot_geneLog_outliers.pdf'))
p3
dev.off()



# ============ Plots for schematic =============

p1 <- make_escheR(spe.subset) |>
  add_fill(var = "sum_gene_outliers") +
  scale_fill_manual(
    name = "", # turn off legend name for ground_truth
    values = c(
      "TRUE" = "red2",
      "FALSE" = "grey")
  ) +
  labs(fill = "Outliers") +
  ggtitle("Outliers") +
  theme(plot.title = element_text(size = 20),
        legend.text = element_text(size=11))

pdf(height = 5, width=5, here(plot_dir, 'Figure2_spotplot_schematic.pdf'))
p1
dev.off()



p1 <- plotColData(spe.subset, x="sample_id", y="sum_gene_z", color_by="sum_gene_outliers", point_size=4) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
  #horizontal line at -3
  geom_hline(yintercept = -3, linetype = "dashed", color = "red", size=1) +
  #theme(text=element_text(size=21)) +
  theme_classic(base_size = 21) +
  xlab("Sample") +
  ylab("z-score") +
  # remove x tick
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")

pdf(height = 3, width=3, here(plot_dir, 'Figure2_schematic_violin.pdf'))
p1
dev.off()




# ====== Panel L =====
p <- ggplot(discarded_df, aes(x = layer_guess_ordered, y = percentage_local, fill = layer_guess_ordered)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, aes(color = layer_guess_ordered), alpha = 0.7) +  # Display all data points with matching layer colors
  scale_fill_manual(values = spatialLIBD::libd_layer_colors) +
  scale_color_manual(values = spatialLIBD::libd_layer_colors) +
  scale_y_continuous(limits = c(0, 60)) +
  labs(title = "Local outliers",
       x = "Spatial Domain",
       y = "Percent discarded (%)") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size=12),
        plot.title = element_text(size = 20),
        text = element_text(size=15))

pdf(height = 5, width=5, here(plot_dir, 'Figure2_boxplot_local.pdf'))
p
dev.off()



discarded_df
# # Groups:   sample_id, layer_guess_ordered, subject [76]
# sample_id layer_guess_ordered subject total_spots discarded_mad discarded_local discarded_threshold percentage_mad percentage_local
# <chr>     <fct>               <chr>         <int>         <int>           <int>               <int>          <dbl>            <dbl>
#   1 151507    Layer1              Br5292          817            14               1                  61          1.71             0.122
# 2 151507    Layer2              Br5292          305             0               0                   0          0                0
# 3 151507    Layer3              Br5292         1215            10               4                  14          0.823            0.329
# 4 151507    Layer4              Br5292          369             9               1                  10          2.44             0.271


discarded_grouped <- pivot_longer(discarded_df,
                                  cols = c("percentage_local", "percentage_threshold", "percentage_mad"),
                                  names_to = "method",
                                  values_to = "percentage_discarded")

# sample_id layer_guess_ordered subject total_spots discarded_mad discarded_local discarded_threshold method               percentage_discarded
# <chr>     <fct>               <chr>         <int>         <int>           <int>               <int> <chr>                               <dbl>
#   1 151507    Layer1              Br5292          817            14               1                  61 percentage_local                    0.122
# 2 151507    Layer1              Br5292          817            14               1                  61 percentage_threshold                7.47
# 3 151507    Layer1              Br5292          817            14               1                  61 percentage_mad                      1.71

# change the method names to Local, Manual, and MAD
discarded_grouped$method <- factor(discarded_grouped$method, levels = c("percentage_local", "percentage_threshold", "percentage_mad"),
                                   labels = c("Local", "Manual", "MAD"))



p1 <- ggplot(discarded_grouped, aes(x = layer_guess_ordered, y = percentage_discarded, fill = layer_guess_ordered)) +
  geom_boxplot() +
  facet_wrap(~method) +  # Optional: Add facets if there are multiple subjects to compare within
  labs(title = "Comparison of Discard Percentage by method",
       x = "Layer",
       y = "Percentage Discarded",
       fill = "Layers") +
  scale_fill_manual(values = spatialLIBD::libd_layer_colors) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size = 15))

png(height = 5, width=10, units="in",res=300, here(plot_dir, 'Figure2_boxplot_by_method.png'))
p1
dev.off()


# average discarded for Layer1, Layer3, and WM per group
test<- discarded_grouped %>%
  filter(layer_guess_ordered %in% c("Layer1", "Layer3", "WM")) %>%
  group_by(layer_guess_ordered, method) %>%
  summarise(mean_discarded = mean(percentage_discarded)) %>%
  pivot_wider(names_from = method, values_from = mean_discarded)
test
# layer_guess_ordered Local Manual   MAD
# 1 Layer1              0.207   9.34 2.19
# 2 Layer3              0.274   4.70 0.553
# 3 WM                  0.401   9.74 5.30


filtered_df <- discarded_df %>%
  filter(layer_guess_ordered %in% c("Layer1", "Layer3", "WM"))

# Summarize the total spots and discarded spots by method and layer
summary_df <- filtered_df %>%
  group_by(sample_id, layer_guess_ordered, subject) %>%
  summarize(
    total_spots = sum(total_spots),
    discarded_local = sum(discarded_local),
    discarded_threshold = sum(discarded_threshold)
  ) %>%
  ungroup()

# Calculate the retained spots by method
summary_df <- summary_df %>%
  mutate(
    retained = discarded_threshold - discarded_local
  )

# Summarize across all samples
final_summary <- summary_df %>%
  group_by(layer_guess_ordered) %>%
  summarize(
    total_retained = sum(retained),
  )

# Display the final summary
print(final_summary)


# ===== Panel M-O =====

# ridge plot of sum_umi with a threshold of 500
p4 <- ggplot(qc_df, aes(x = sum_umi_z, y = layer, fill = layer)) +
  geom_density_ridges(alpha = 0.8, scale = 1.5) +
  geom_vline(xintercept = -3, linetype = "dashed", color = "red", size=1) +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(size = 20),
        text = element_text(size = 15)) +
  scale_fill_manual(values = libd_layer_colors) +
  labs(title = "Library size (local z-score)",
       x = "Z-score",
       y = "Spatial Domain")

# ridge plot of expr_ratio with a threshold of 0.3
p5 <- ggplot(qc_df, aes(x = expr_ratio_z, y = layer, fill = layer)) +
  geom_density_ridges(alpha = 0.8, scale = 1.5) +
  geom_vline(xintercept = 3, linetype = "dashed", color = "red", size=1) +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(size = 20),
        text = element_text(size = 15)) +
  scale_fill_manual(values = libd_layer_colors) +
  labs(title = "Mito ratio (local z-score)",
       x = "Z-score",
       y = "Spatial Domain")

# ridge plot of sum_gene with a threshold of 500
p6 <- ggplot(qc_df, aes(x = sum_gene_z, y = layer, fill = layer)) +
  geom_density_ridges(alpha = 0.8, scale = 1.5) +
  geom_vline(xintercept = -3, linetype = "dashed", color = "red", size=1) +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(size = 20),
        text = element_text(size = 15)) +
  scale_fill_manual(values = libd_layer_colors) +
  labs(title = "Unique genes (local z-score)",
       x = "Z-score",
       y = "Spatial Domain")

pdf(height = 5, width=5, here(plot_dir, 'Figure2_ridge_umi_z.pdf'))
p4
dev.off()

pdf(height = 5, width=5, here(plot_dir, 'Figure2_ridge_mito_z.pdf'))
p5
dev.off()

pdf(height = 5, width=5, here(plot_dir, 'Figure2_ridge_gene_z.pdf'))
p6
dev.off()




# ==== Mi QC ======


p1 <- plotColData(spe, x="sum_gene", y="expr_chrM_ratio", color_by="local_outliers", point_size=1) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
  theme_bw() +
  labs(title = "Local Outliers",
       x = "Number of genes",
       y = "Mitochondrial Ratio") +
  theme(plot.title = element_text(size = 20),
        text = element_text(size = 15)) +
  guides(color = guide_legend(title = "Discarded"))

png(height = 5, width=5, units="in",res=300, here(plot_dir, 'Figure2_geneVSmito_local.png'))
p1
dev.off()


