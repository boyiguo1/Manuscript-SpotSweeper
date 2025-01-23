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

# large DLPFC dataset
load(here("processed-data","Visium","dlPFC_raw.Rdata"))
spe.pfc <- spe_raw
spe.pfc

# drop out of tissue spots
spe.pfc <- spe.pfc[,colData(spe.pfc)$in_tissue == TRUE]
spe.hangnail <- spe.pfc[,colData(spe.pfc)$sample_id == "Br8325_ant"]



# ======= Artifact detection using various multiscale resolutions =======

spe.hangnail <- findArtifacts(spe.hangnail,
                            mito_percent="expr_chrM_ratio",
                            mito_sum="expr_chrM",
                            n_rings=2,
                            name="artifact_n2"
                            )
spe.hangnail$artifact_n2 <- spe.hangnail$artifact


spe.hangnail <- findArtifacts(spe.hangnail,
                            mito_percent="expr_chrM_ratio",
                            mito_sum="expr_chrM",
                            n_rings=3,
                            name="artifact_n3"
                            )
spe.hangnail$artifact_n3 <- spe.hangnail$artifact


spe.hangnail <- findArtifacts(spe.hangnail,
                            mito_percent="expr_chrM_ratio",
                            mito_sum="expr_chrM",
                            n_rings=4,
                            name="artifact_n4"
                            )
spe.hangnail$artifact_n4 <- spe.hangnail$artifact


spe.hangnail <- findArtifacts(spe.hangnail,
                            mito_percent="expr_chrM_ratio",
                            mito_sum="expr_chrM",
                            n_rings=5,
                            name="artifact_n5"
                            )
spe.hangnail$artifact_n5 <- spe.hangnail$artifact

spe.hangnail <- findArtifacts(spe.hangnail,
                            mito_percent="expr_chrM_ratio",
                            mito_sum="expr_chrM",
                            n_rings=6,
                            name="artifact_n6"
                            )
spe.hangnail$artifact_n6 <- spe.hangnail$artifact


spe.hangnail <- findArtifacts(spe.hangnail,
                            mito_percent="expr_chrM_ratio",
                            mito_sum="expr_chrM",
                            n_rings=7,
                            name="artifact_n7"
                            )  
spe.hangnail$artifact_n7 <- spe.hangnail$artifact


# ======= Compare the results of different resolutions =======


png(here("plots", "parameter_selection", "artifact_detection", "spotplots.png"), unit="in", width=6, height=4, res=300)
p1 <- plotSpotQC(spe.hangnail, plot_type = "spot", in_tissue = "in_tissue", 
                 annotate = "artifact_n2", point_size = 0.2) +
                guides(color = "none") +
                ggtitle("1st-2nd order")

p2 <- plotSpotQC(spe.hangnail, plot_type = "spot", in_tissue = "in_tissue", 
                 annotate = "artifact_n3", point_size = 0.2) +
                guides(color = "none") +
                ggtitle("1st-3rd order")


p3 <- plotSpotQC(spe.hangnail, plot_type = "spot", in_tissue = "in_tissue", 
                 annotate = "artifact_n4", point_size = 0.2) +
                guides(color = "none") +
                ggtitle("1st-4th order")

p4 <- plotSpotQC(spe.hangnail, plot_type = "spot", in_tissue = "in_tissue", 
                 annotate = "artifact_n5", point_size = 0.2) +
                guides(color = "none") +
                ggtitle("1st-5th order")

p5 <- plotSpotQC(spe.hangnail, plot_type = "spot", in_tissue = "in_tissue", 
                 annotate = "artifact_n6", point_size = 0.2) +
                guides(color = "none") +
                ggtitle("1st-6th order")

p6 <- plotSpotQC(spe.hangnail, plot_type = "spot", in_tissue = "in_tissue", 
                 annotate = "artifact_n7", point_size = 0.2) +
                guides(color = "none") +
                ggtitle("1st-7th order")

(p1 | p2 | p3) / (p4 | p5 | p6)
dev.off()


# ======= Quantifying change in artifacts ========

# Extract artifacts per sample
artifacts_per_n <- colData(spe.hangnail) %>%
  as.data.frame() %>%
  select(starts_with("artifact_n")) %>%
  mutate(spot = rownames(.), sample_id = colData(spe.hangnail)$sample_id)  # Add spot and sample info

# Get reference artifacts for n_order = 5
reference_artifacts <- artifacts_per_n %>%
  filter(artifact_n5 == 1) %>%
  pull(spot)  

# Calculate shared and additional spots per sample
comparison_per_sample_artifacts <- artifacts_per_n %>%
  pivot_longer(
    cols = starts_with("artifact_n"),
    names_to = "n_order",
    values_to = "is_artifact"
  ) %>%
  mutate(n_order = sub("artifact_", "", n_order)) %>%  
  group_by(sample_id, n_order) %>%
  summarise(
    total_spots = sum(is_artifact == 1, na.rm = TRUE),
    shared_spots = sum(spot %in% reference_artifacts & is_artifact == 1, na.rm = TRUE),
    additional_spots = sum(!spot %in% reference_artifacts & is_artifact == 1, na.rm = TRUE)
  ) %>%
  mutate(
    percent_shared = 100 * shared_spots / total_spots  
  ) %>%
  filter(n_order != "n5")

# Reorder n_order as a factor
comparison_per_sample_artifacts <- comparison_per_sample_artifacts %>%
  mutate(n_order = factor(n_order, levels = c("n2", "n3", "n4", "n6", "n7")))

# BAR plot of shared artifacts
png(here("plots", "parameter_selection", "artifact_detection", "artifact_overlap.png"), unit="in", width=3.5, height=3.5, res=300)
ggplot(comparison_per_sample_artifacts, aes(x = n_order, y = percent_shared, fill=n_order)) +
  geom_bar(stat = "identity", alpha = 0.7) +
  geom_text(aes(label = paste0(round(percent_shared, 1), "%")), vjust = -0.3, size = 3) +
  labs(
    title = "Shared Spots w/ n=5",
    x = "Highest N-order neighbors",
    y = "Percentage of Shared Spots"
  ) +
  ggpubr::theme_pubr() +
  ylim(0, 100) +
  scale_fill_manual(values = color_palette) +
  theme(legend.position = "none")
dev.off()

# BAR plot of additional artifacts
png(here("plots", "parameter_selection", "artifact_detection", "artifact_overlap_additional.png"), unit="in", width=3.5, height=3.5, res=300)
ggplot(comparison_per_sample_artifacts, aes(x = n_order, y = additional_spots, fill=n_order)) +
  geom_bar(stat = "identity", alpha = 0.7) +
  geom_text(aes(label = additional_spots), vjust = -0.3, size = 3) +
  labs(
    title = "Additional Spots to n=5",
    x = "Highest N-order neighbors",
    y = "Number of Additional Spots"
  ) +
  ggpubr::theme_pubr() +
  scale_fill_manual(values = color_palette) +
  theme(legend.position = "none")
dev.off()