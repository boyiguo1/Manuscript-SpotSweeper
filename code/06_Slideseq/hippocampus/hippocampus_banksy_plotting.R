library(SpotSweeper)
library(RANN)
library(here)
library(SpatialExperiment)
library(spatialLIBD)
library(escheR)
library(patchwork)
library(scran)
library(scater)


plot_dir <- here("plots","Slideseq", "hippocampus")


spe <- readRDS(here("processed-data", "Slideseq","SlideseqV2Hippocampus","slideseq_hippocampus_banksy.rds"))
spe

# ===== Plot Banksy spatial domains =====

# make color palette for number of domains
n_domains <- length(unique(colData(spe)$clust_M0_lam0.8_k50_res0.8))
domain_pal <- colorRampPalette(RColorBrewer::brewer.pal(n_domains, "Set1"))

pdf(here(plot_dir, "hippocampus_banksy_domains.pdf"), width=5, height=5)
make_escheR(spe) |>
    add_fill(var="clust_M0_lam0.8_k50_res0.8", point_size=.7) +
    scale_fill_manual(values=domain_pal(n_domains)) +
    theme(legend.position="none")
dev.off()

# make df and export to csv
banksy_df <- data.frame(
    spatialCoords = spatialCoords(spe),
    sample_id = colData(spe)$sample_id,
    banksy_domain = colData(spe)$clust_M0_lam0.8_k50_res0.8
)

write.csv(banksy_df, 
          file = here("processed-data", "outputs_for_paper", "figure_4", "Figure4_C.csv"),
          row.names = FALSE, quote = FALSE)
