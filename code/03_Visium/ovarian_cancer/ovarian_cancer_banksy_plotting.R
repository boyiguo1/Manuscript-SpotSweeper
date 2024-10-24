library(SpotSweeper)
library(RANN)
library(here)
library(SpatialExperiment)
library(spatialLIBD)
library(escheR)
library(patchwork)
library(scran)
library(scater)


plot_dir <- here("plots","Visium", "ovarian_cancer")


spe <- readRDS(here("processed-data", "Visium","ovarian_cancer","OvarianCancer_banksy.rds"))
spe


# ===== Plot Banksy spatial domains =====
unique(spe$sample_id)
# [1] "GSM6506110" "GSM6506111" "GSM6506112" "GSM6506113" "GSM6506114"
# [6] "GSM6506115" "GSM6506116" "GSM6506117"

# subset to first sample
spe <- spe[,spe$sample_id == "GSM6506110"]

# make color palette for number of domains
n_domains <- length(unique(colData(spe)$clust_M0_lam0.2_k50_res0.8))
domain_pal <- colorRampPalette(RColorBrewer::brewer.pal(n_domains, "Set1"))

png(here(plot_dir, "ovarian_banksy_domains.png"), width=5, height=5, units="in", res=300)
make_escheR(spe) |>
    add_fill(var="clust_M0_lam0.2_k50_res0.8", point_size=2.5) +
    scale_fill_manual(values=domain_pal(n_domains)) +
    theme(legend.position="none")
dev.off()