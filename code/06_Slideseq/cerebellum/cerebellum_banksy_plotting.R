library(SpotSweeper)
library(RANN)
library(here)
library(SpatialExperiment)
library(spatialLIBD)
library(escheR)
library(patchwork)
library(scran)
library(scater)


plot_dir <- here("plots","Slideseq", "cerebellum")


spe <- readRDS(here("processed-data",  "Slideseq","SlideseqCerebellum", "slideseq_cerebellum_banksy.rds"))
spe

# ===== Plot Banksy spatial domains =====

# make color palette for number of domains
n_domains <- length(unique(colData(spe)$clust_M0_lam0.3_k50_res0.7))
domain_pal <- colorRampPalette(RColorBrewer::brewer.pal(n_domains, "Set1"))

png(here(plot_dir, "banksy_domains_new.png"), width=5, height=5, units="in", res=300)
make_escheR(spe) |>
    add_fill(var="clust_M0_lam0.3_k50_res0.7", point_size=1) +
    scale_fill_manual(values=domain_pal(n_domains)) +
    theme(legend.position="none")
dev.off()