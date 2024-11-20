library(SpotSweeper)
library(RANN)
library(here)
library(SpatialExperiment)
library(spatialLIBD)
library(ggspavis)
library(BiocNeighbors)
library(escheR)
library(patchwork)
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
spe <- localOutliers(spe, metric="sum_umi", direction="lower", log=TRUE)
spe <- localOutliers(spe, metric="sum_gene", direction="lower", log=TRUE)
spe <- localOutliers(spe, metric="expr_chrM_ratio", direction="higher", log=FALSE)

spe$local_outliers <- as.logical(spe$sum_umi_outliers) |
  as.logical(spe$sum_gene_outliers) |
  as.logical(spe$expr_chrM_ratio_outliers)


colnames(colData(spe))
# [69] "array_col"                   "sum_umi_log"                
# [71] "sum_umi_outliers"            "sum_umi_z"    


# # consider outliers if theabs(z-score) by dividing 0.6745 is > 3
# spe$sum_umi_outliers_corrected <- abs(spe$sum_umi_z)*0.6745 > 3
# spe$sum_gene_outliers_corrected <- abs(spe$sum_gene_z)*0.6745 > 3
# spe$expr_chrM_ratio_outliers_corrected <- abs(spe$expr_chrM_ratio_z)*0.6745 > 3



# #subset to first sample
# spe.subset <- subset(spe, ,sample_id == unique(spe$sample_id)[1])
# png(here(plot_dir, "Figure2_spotplot_local_umi.png"), width=10, height=5, units="in", res=300)
# p1 <- make_escheR(spe.subset ) |> 
#   add_fill(var = "sum_umi_outliers") +
#   scale_fill_manual(values = c("TRUE" = "red2", "FALSE" = "grey")) 

# p2 <- make_escheR(spe.subset ) |> 
#   add_fill(var = "sum_umi_outliers_corrected") +
#   scale_fill_manual(values = c("TRUE" = "red2", "FALSE" = "grey"))
#   p1+p2
# dev.off()


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





# ============= miQC =============

library(miQC)

rownames(spe) <- rowData(spe)$gene_name

mt_genes <- grepl("^MT-",  rownames(spe))
feature_ctrls <- list(mito = rownames(spe)[mt_genes])

spe <- addPerCellQC(spe, subsets = feature_ctrls)

model <- mixtureModel(spe, model_type = "spline")

png(here(plot_dir, "miQC_mixtureModel.png"), width=4, height=4, units="in", res=300)
plotModel(spe, model) + ggtitle("Visium DLPFC")
dev.off()

png(here(plot_dir, "miQC_mixtureModelFiltering.png"), width=3.5, height=3.5, units="in", res=300)
plotFiltering(spe, model) + ggtitle("miQC outliers") +
     scale_color_manual(values = c("TRUE" = "grey", "FALSE" = "red"))
dev.off()

# NOTE: miQC does not offer a way to simply flag spots, so I'm importing part of the filter function here to do so manually
library(flexmix)

posterior_cutoff = 0.75
keep_all_below_boundary = TRUE
enforce_left_cutoff = TRUE

metrics <- as.data.frame(colData(spe))

    if (is.null(model)) {
        warning("call 'mixtureModel' explicitly to get stable model features")
        model <- mixtureModel(spe)
    }

    intercept1 <- parameters(model, component = 1)[1]
    intercept2 <- parameters(model, component = 2)[1]
    if (intercept1 > intercept2) {
        compromised_dist <- 1
        intact_dist <- 2
    } else {
        intact_dist <- 1
        compromised_dist <- 2
    }

    post <- posterior(model)
    metrics$prob_compromised <- post[, compromised_dist]
    spe$prob_compromised <- metrics$prob_compromised
    metrics$keep <- metrics$prob_compromised <= posterior_cutoff

    if (sum(metrics$keep)==nrow(metrics)){
        stop("all cells passed posterior probability filtering. One 
              cause of this is the model selecting two near-identical
              distributions. Try rerunning mixtureModel() and/or 
              setting a different random seed.")
    }
    
    if (keep_all_below_boundary == TRUE) {
        predictions <- fitted(model)[, intact_dist]
        metrics$intact_prediction <- predictions
        metrics[metrics$subsets_mito_percent <
                    metrics$intact_prediction, ]$keep <- TRUE
    }

    if (enforce_left_cutoff == TRUE) {
        min_discard <- min(metrics[!metrics$keep, ]$subsets_mito_percent)
        min_index <- which(metrics$subsets_mito_percent == min_discard)[1]
        lib_complexity <- metrics[min_index, ]$detected
        metrics[metrics$detected <= lib_complexity &
                    metrics$subsets_mito_percent >= min_discard, ]$keep <- FALSE
    }

# add "outliers" to SPE
spe$miQC_keep <- !metrics$keep


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
                 discard_miQC = spe$miQC_keep,
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
    discarded_miQC = sum(discard_miQC, na.rm = TRUE),  # Counting TRUE as 1 and FALSE as 0
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
  ) %>%
  mutate(
    percentage_miQC = (discarded_miQC/ total_spots) * 100  # Calculating the percentage
  )


discarded_long <- discarded_df %>%
  pivot_longer(
    cols = starts_with("percentage_"),
    names_to = "QC_method",
    values_to = "percentage_discarded"
  ) %>%
  mutate(QC_method = recode(QC_method,
                            "percentage_miQC" = "miQC",
                            "percentage_local" = "SpotSweeper",
                            "percentage_mad" = "MAD",
                            "percentage_threshold" = "Threshold"
))

# Reorder QC_method factor
discarded_long <- discarded_long %>%
  mutate(QC_method = factor(QC_method, levels = c("miQC",  "Threshold", "MAD","SpotSweeper")))



# Reverse the order of the layers
discarded_long <- discarded_long %>%
  mutate(layer_guess_ordered = factor(layer_guess_ordered, levels = rev(levels(layer_guess_ordered))))

# Plot 
p <- ggplot(discarded_long, aes(x = layer_guess_ordered, y = percentage_discarded, fill = QC_method)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, position = position_dodge(width = 0.75)) +  # Dodge to separate methods
  scale_y_continuous(limits = c(0, 60)) +
  labs(
    title = NULL,
    x = "Spatial Domain",
    y = "Percent outliers (%)",
    fill = "Method" 
  ) +
  theme_classic() +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    axis.text.x = element_text( hjust = 1, size = 12),
    plot.title = element_text(size = 18),
    text = element_text(size = 16)) +
  scale_fill_manual(values = c("#459395", "#eb7c69", "#fda638", "#b82ac9"), guide = guide_legend(nrow = 2)) +
  coord_flip()

# Save the plot to PNG
png(file = here(plot_dir, 'Figure2_boxplot_all_methods_colored.png'), height = 8, width = 4.5, units = "in", res = 300)
print(p)
dev.off()

# write csv
write.csv(discarded_df, here(processed_dir,"outputs_for_paper", 'Figure2_boxplots.csv'), row.names = FALSE)


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
  guides(color = guide_legend(title = "Outlier")) +
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
  guides(color = guide_legend(title = "Outlier")) +
  ggtitle("Threshold") +
  theme(plot.title = element_text(size = 20),
        legend.text = element_text(size=11))

pdf(height = 5, width=5, here(plot_dir, 'Figure2_spotplot_gene_threshold.pdf'))
p3
dev.off()

p <- make_escheR(spe.subset) |>
  add_fill(var = "sum_gene") |>
  add_ground(var = "discard_mad", stroke = 1) +
  scale_color_manual(
    name = "", # turn off legend name for ground_truth
    values = c(
      "TRUE" = "red2",
      "FALSE" = "transparent")
  ) +
  scale_fill_gradient(low ="white",high =  "black") +
  labs(fill = "Sum Genes") +
  guides(color = guide_legend(title = "Outlier")) +
  ggtitle("MAD") +
  theme(plot.title = element_text(size = 20),
        legend.text = element_text(size=11))

pdf(height = 5, width=5, here(plot_dir, 'Figure2_spotplot_gene_MAD.pdf'))
p
dev.off()


p <- make_escheR(spe.subset) |>
  add_fill(var = "sum_gene") |>
  add_ground(var = "local_outliers", stroke = 1) +
  scale_color_manual(
    name = "", # turn off legend name for ground_truth
    values = c(
      "TRUE" = "red2",
      "FALSE" = "transparent")
  ) +
  scale_fill_gradient(low ="white",high =  "black") +
  labs(fill = "Sum Genes") +
  guides(color = guide_legend(title = "Outlier")) +
  ggtitle("Local outliers") +
  theme(plot.title = element_text(size = 20),
        legend.text = element_text(size=11))

pdf(height = 5, width=5, here(plot_dir, 'Figure2_spotplot_gene_local.pdf'))
p
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
  guides(color = guide_legend(title = "Outlier")) +
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

# Extract the sum values that are marked as outliers
outliers <- spe$sum_umi[spe$qc_umi_mad]
sum_3MAD <- max(outliers)
sum_3MAD

outliers <- spe$expr_chrM_ratio[spe$qc_mito_mad]
mito_3MAD <- min(outliers)
mito_3MAD

outliers <- spe$sum_gene[spe$qc_gene_mad]
detected_3MAD <- max(outliers)
detected_3MAD



# ridge plot of sum_umi with a threshold of 500
p4 <- ggplot(qc_df, aes(x = sum_umi, y = layer, fill = layer)) +
  geom_density_ridges(alpha = 0.8, scale = 1.5) +
  geom_vline(xintercept = 500, linetype = "dashed", color = "#eb7c69", size=1.75) +
    geom_vline(xintercept = sum_3MAD, linetype = "dashed", color = "#fda638", size=1.75) +
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
  geom_vline(xintercept = 0.275, linetype = "dashed", color = "#eb7c69", size=1.75) +
    geom_vline(xintercept = mito_3MAD, linetype = "dashed", color = "#fda638", size=1.75) +
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
  geom_vline(xintercept = 500, linetype = "dashed", color = "#eb7c69", size=1.75) +
    geom_vline(xintercept = detected_3MAD, linetype = "dashed", color = "#fda638", size=1.75) +
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




# ======= Local outlier ridge plots =======

# ridge plot of sum_umi with a threshold of 500
p4 <- ggplot(qc_df, aes(x = sum_umi_z, y = layer, fill = layer)) +
  geom_density_ridges(alpha = 0.8, scale = 1.5) +
  geom_vline(xintercept = -3, linetype = "dashed", color = "#b82ac9", size=1.75) +
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
  geom_vline(xintercept = 3, linetype = "dashed", color = "#b82ac9", size=1.75) +
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
  geom_vline(xintercept = -3, linetype = "dashed", color = "#b82ac9", size=1.75) +
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




# ====== spotplot grid of all samples ======
# loop through all samples

samples <- unique(spe$sample_id)

for (sample in samples) {
  spe.subset <- subset(spe, ,sample_id == sample)
  p1 <- make_escheR(spe.subset) |>
  add_fill(var = "sum_umi") |>
  add_ground(var = "local_outliers", stroke = 1) +
  scale_color_manual(
    name = "", # turn off legend name for ground_truth
    values = c(
      "TRUE" = "red2",
      "FALSE" = "transparent")
  ) +
  scale_fill_gradient(low ="white",high =  "black") +
  labs(fill = "Sum UMI") +
  guides(color = guide_legend(title = "Outlier")) +
  ggtitle(paste0("Sample: ", sample)) +
  theme(plot.title = element_text(size = 20),
        legend.text = element_text(size=11))

  png(here(plot_dir,"spotplots", paste0('Figure2_spotplot_', sample, '.png')), width = 5, height = 5, units = "in", res = 300)
  print(p1)
  dev.off()
}


# save spe
saveRDS(spe, here(processed_dir,"figure_2", 'spe_spotsweeper.rds'))












# ===== Other stuff below ======


# == Plots for schematic ==

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
