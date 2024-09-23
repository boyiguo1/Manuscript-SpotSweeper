
library(SpotSweeper)
library(RANN)
library(here)
library(SpatialExperiment)
library(spatialLIBD)
library(escheR)
library(patchwork)
library(scran)
library(scater)


plot_dir <- here("plots","VisiumHD","mouse_brain")


spe <- readRDS(here("processed-data", "VisiumHD","VisiumHD_MouseBrain_016_banksy.rds"))
spe




# ======== spot plots ========
png(here(plot_dir, "VisiumHD_MB_banksy_spotplots.png"), width=10, height=10, units="in", res=300)
make_escheR(spe) |>
  add_ground(var="clust_M0_lam0.8_k50_res0.8", point_size=.3)
dev.off()



# ========= violin plots =========
spe$sum_log <- log2(spe$sum+1)
spe$detected_log <- log2(spe$detected+1)

png(here(plot_dir, "VisiumHD_MB_banksy_violins_nCount.png"), width=8, height=9, units="in", res=300)
p1 <- plotColData(spe, x="clust_M0_lam0.8_k50_res0.8", y="sum_log", color_by="clust_M0_lam0.8_k50_res0.8") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #scale_color_manual(values=c("TRUE"="red", "FALSE"="grey")) +
  theme(legend.position="none") +
  labs(title="Library Size",
       x="Cell Type") +
  coord_flip()
p2 <- plotColData(spe, x="clust_M0_lam0.8_k50_res0.8", y="detected_log", color_by="clust_M0_lam0.8_k50_res0.8") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #scale_color_manual(values=c("TRUE"="red", "FALSE"="grey")) +
  #theme(legend.position="none") +
  labs(title="Num Genes",
       x=NULL) +
  coord_flip()
p1+p2
dev.off()

# ======= UMAPs =======

plotReducedDim(spe, dimred="PCA_M0_lam0.8", color_by="clust_M0_lam0.8_k50_res0.8")

spe <- runUMAP(spe, assay.type="counts", ncomponents = 30, ntop = 500)

png(here(plot_dir, "VisiumHD_MB_umap_banksy_clusters.png"), width=5, height=5, units="in", res=300)
plotReducedDim(spe, dimred="UMAP", color_by="clust_M0_lam0.8_k50_res0.8")
dev.off()

png(here(plot_dir, "VisiumHD_MB_umap_banksy_nCounts.png"),  width=5, height=5, units="in", res=300)
plotReducedDim(spe, dimred="UMAP", color_by="nCounts_log")
dev.off()

png(here(plot_dir, "VisiumHD_MB_umap_banksy_nCounts_threshold.png"),  width=5, height=5, units="in", res=300)
plotReducedDim(spe, dimred="UMAP", color_by="nCounts_discard_threshold") +
  scale_color_manual("color", values=c("TRUE"="red", "FALSE"="grey"))
dev.off()

png(here(plot_dir, "VisiumHD_MB_umap_banksy_nCounts_3mad.png"),  width=5, height=5, units="in", res=300)
plotReducedDim(spe, dimred="UMAP", color_by="nCounts_discard_3mad")+
  scale_color_manual("color", values=c("TRUE"="red", "FALSE"="grey"))
dev.off()

png(here(plot_dir, "VisiumHD_MB_umap_banksy_nCounts_spotsweeper.png"),  width=5, height=5, units="in", res=300)
plotReducedDim(spe, dimred="UMAP", color_by="nCounts_outliers")+
  scale_color_manual("color", values=c("TRUE"="red", "FALSE"="grey"))
dev.off()
