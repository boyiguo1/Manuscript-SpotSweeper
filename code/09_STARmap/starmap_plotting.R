
library(SpotSweeper)
library(RANN)
library(here)
library(SpatialExperiment)
library(spatialLIBD)
library(escheR)
library(patchwork)
library(scran)
library(scater)
library(dplyr)

plot_dir <- here("plots","STARmap")

spe <-  readRDS(here("processed-data",  "STARmap",  "starmap_spotsweeper_banksy.rds"))
spe

spe$nCounts_discard_threshold <- spe$nCounts < 35
spe$nGenes_discard_threshold  <- spe$nGenes < 25

# total local outliers
spe$local_outliers <- as.logical(spe$nCounts_outliers) |
  as.logical(spe$nGenes_outliers)

# drop z in spatialCoords
spatialCoords(spe) <- spatialCoords(spe)[,1:2]

colnames(colData(spe))
#  [1] "NAME"                         "Main_molecular_cell_type"    
#  [3] "Sub_molecular_cell_type"      "Main_molecular_tissue_region"
#  [5] "Sub_molecular_tissue_region"  "Molecular_spatial_cell_type" 
#  [7] "sample_id"                    "nCounts"                     
#  [9] "nGenes"                       "nCounts_discard_threshold"   
# [11] "nGenes_discard_threshold"     "nCounts_discard_3mad"        
# [13] "nGenes_discard_3mad"          "nCounts_log"                 
# [15] "nCounts_outliers"             "nCounts_z"                   
# [17] "nGenes_log"                   "nGenes_outliers"             
# [19] "nGenes_z"                     "clust_M0_lam0.8_k50_res0.8"  
# [21] "local_outliers"    

# 3MAD outliers discard
spe$MAD_outliers <- as.logical(spe$nCounts_discard_3mad) |
  as.logical(spe$nGenes_discard_3mad)

# total threshold outliers
spe$threshold_outliers <- as.logical(spe$nCounts_discard_threshold) |
    as.logical(spe$nGenes_discard_threshold)


# palette for length spe$clust_M0_lam0.8_k50_res0.8
pal <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(length(unique(spe$clust_M0_lam0.8_k50_res0.8)))


# ========= Create dataframe of outliers by domain ========
outlier_df <- data.frame(MAD=spe$MAD_outliers,
                         Threshold=spe$threshold_outliers, 
                         SpotSweeper=spe$local_outliers,
                         domain=factor(spe$clust_M0_lam0.8_k50_res0.8),
                         nCounts = spe$nCounts,
                         nGenes = spe$nGenes,
                         #subsets_mito_percent = spe$subsets_mito_percent,
                         nCounts_z = spe$nCounts_z,
                         nGenes_z = spe$nGenes_z
                         #subsets_mito_percent_z = spe$subsets_mito_percent_z
                        )

# get the number of outliers per domain
outlier_df <- outlier_df %>%
    group_by(domain) %>%
    summarise(MAD=sum(MAD),
              SpotSweeper=sum(SpotSweeper),
              Threshold=sum(Threshold)) %>% 
    tidyr::pivot_longer(cols=c(Threshold, MAD, SpotSweeper), names_to="method", values_to="count")

# change plotting order to Threhsold, Mad, SpotSweeper
outlier_df$method <- factor(outlier_df$method, levels=c("Threshold", "MAD", "SpotSweeper"))

# bar plot sw vs mad outliers per domain
png(here(plot_dir, "STARmap_outliers_per_domain.png"), width=15, height=5, res=300, units="in")
ggplot(outlier_df, aes(x=domain, y=count, fill=method)) +
    geom_bar(stat="identity", position="dodge") +
    labs(x="Spatial Domain", y="Number of Outliers", fill="Method") +
    theme_classic() +
    scale_fill_manual(values=c("#459395", "#eb7c69","#fda638")) +
    theme(legend.position=c(.7, .9)) +
    theme(text = element_text(size=16)) 
dev.off()




# ============ Plotting here =============
png(here(plot_dir, "STARmap_Banksy_domains.png"), width=4, height=5, units="in", res=300)
make_escheR(spe) |>
    add_fill(var="clust_M0_lam0.8_k50_res0.8", point_size=.5) +
    scale_fill_manual(values=pal,
    name="") +
      theme(plot.title = element_text(size = 24, hjust = 0.1),
        legend.text = element_text(size=14)) +
    theme(legend.position="none") +
    ggtitle("Spatial Domains") +
        coord_flip()
dev.off()


png(here(plot_dir, "STARmap_nCounts_3MAD.png"), width=5, height=5, res=300, units="in")
plotQC(spe, metric="nGenes", outliers="nCounts_discard_3mad", point_size=0.6, stroke=.8) +
    ggtitle("3 MAD") +
  theme(plot.title = element_text(size = 20, hjust = 0.1),
        legend.text = element_text(size=11))+
        coord_flip()
dev.off()

png(here(plot_dir, "STARmap_nCounts_threshold.png"), width=5, height=5, res=300, units="in")
plotQC(spe, metric="nGenes", outliers="nCounts_discard_threshold", point_size=0.6, stroke=.8) +
    ggtitle("Fixed Threshold") +
  theme(plot.title = element_text(size = 20, hjust = 0.1),
        legend.text = element_text(size=11))+
        coord_flip()
dev.off()

png(here(plot_dir, "STARmap_nCounts_spotsweeper.png"), width=5, height=5, res=300, units="in")
plotQC(spe, metric="nGenes", outliers="nCounts_outliers", point_size=0.6, stroke=.8) +
    ggtitle("SpotSweeper") +
  theme(plot.title = element_text(size = 20, hjust = 0.1),
        legend.text = element_text(size=11)) +
        coord_flip()
dev.off()



# ===== Threshold Ridge plots =====
outlier_df <- data.frame(MAD=spe$MAD_outliers,
                         Threshold=spe$threshold_outliers, 
                         SpotSweeper=spe$local_outliers,
                         domain=factor(spe$clust_M0_lam0.8_k50_res0.8),
                         nCounts = spe$nCounts,
                         nGenes = spe$nGenes,
                         #subsets_mito_percent = spe$subsets_mito_percent,
                         nCounts_z = spe$nCounts_z,
                         nGenes_z = spe$nGenes_z
                         #subsets_mito_percent_z = spe$subsets_mito_percent_z
                        )



library(ggridges)
# ridge plot of nCounts_umi with a threshold of 500
png(here(plot_dir, "STARmap_ridge_nCounts_umi.png"), width=4, height=6, res=300, units="in")
ggplot(outlier_df, aes(x = nCounts, y = domain, fill = domain)) +
  geom_density_ridges(alpha = 0.8, scale = 1.5) +
  geom_vline(xintercept = 35, linetype = "dashed", color = "red", size=1) +
  scale_x_continuous(trans='log10') +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(size = 18),
        text = element_text(size = 15)) +
  scale_fill_manual(values = pal) +
  labs(title = "Library size per domain",
       x = "Library Size",
       y = "Spatial Domain") +
       coord_cartesian(xlim = c(20, NA))
dev.off()


png(here(plot_dir, "STARmap_ridge_unique_genes.png"), width=4, height=6, res=300, units="in")
ggplot(outlier_df, aes(x = nGenes, y = domain, fill = domain)) +
  geom_density_ridges(alpha = 0.8, scale = 1.5) +
  geom_vline(xintercept = 35, linetype = "dashed", color = "red", size=1) +
  theme_bw() +
  scale_x_continuous(trans='log10') +
  theme(legend.position = "none",
        plot.title = element_text(size = 18),
        text = element_text(size = 15)) +
  scale_fill_manual(values = pal) +
  labs(title = "Unique genes per domain",
       x = "Number of genes",
       y = "Spatial Domain")+
       coord_cartesian(xlim = c(20, NA))
dev.off()




# ===== SpotSweeper ridge plots =====
# ridge plot of nCounts_umi with a threshold of 500
png(here(plot_dir, "STARmap_ridge_nCounts_umi_spotsweeper.png"), width=4, height=6, res=300, units="in")
ggplot(outlier_df, aes(x = nCounts_z, y = domain, fill = domain)) +
  geom_density_ridges(alpha = 0.8, scale = 1.5) +
  geom_vline(xintercept = -3, linetype = "dashed", color = "red", size=1) +
  #scale_x_continuous(trans='log10') +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(size = 18),
        text = element_text(size = 15)) +
  scale_fill_manual(values = pal) +
  labs(title = "Library size (local z-score)",
       x = "Z-score",
       y = "Spatial Domain") 
dev.off()


png(here(plot_dir, "STARmap_ridge_unique_genes_spotsweeper.png"), width=4, height=6, res=300, units="in")
ggplot(outlier_df, aes(x = nGenes_z, y = domain, fill = domain)) +
  geom_density_ridges(alpha = 0.8, scale = 1.5) +
  geom_vline(xintercept = -3, linetype = "dashed", color = "red", size=1) +
  theme_bw() +
  #scale_x_continuous(trans='log10') +
  theme(legend.position = "none",
        plot.title = element_text(size = 18),
        text = element_text(size = 15)) +
  scale_fill_manual(values = pal) +
  labs(title = "Unique genes (local z-score)",
       x = "Z-score",
       y = "Spatial Domain")
dev.off()


