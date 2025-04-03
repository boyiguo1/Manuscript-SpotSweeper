library(SpotSweeper)
library(RANN)
library(here)
library(SpatialExperiment)
library(spatialLIBD)
library(escheR)
library(patchwork)
library(scran)
library(scater)


plot_dir <- here("plots","STOmics")


spe <- readRDS(here("processed-data", "STOmics","STOmics_prenatal.rds"))
spe

# ===== Plot ground truth spatial domains =====

# make color palette for number of domains
n_domains <- length(unique(colData(spe)$annotation))
domain_pal <- colorRampPalette(RColorBrewer::brewer.pal(n_domains, "Set1"))

pdf(here(plot_dir, "stomics_annotations.pdf"), width=5, height=5)
make_escheR(spe) |>
    add_fill(var="annotation", point_size=.7) +
    scale_fill_manual(values=domain_pal(n_domains)) +
    theme(legend.position="none")
dev.off()


# make df and export to csv
banksy_df <- data.frame(
    spatialCoords = spatialCoords(spe),
    sample_id = colData(spe)$sample_id,
    annotation = colData(spe)$annotation
)

write.csv(banksy_df, 
          file = here("processed-data", "outputs_for_paper", "figure_4", "Figure4_D.csv"),
          row.names = FALSE, quote = FALSE)
