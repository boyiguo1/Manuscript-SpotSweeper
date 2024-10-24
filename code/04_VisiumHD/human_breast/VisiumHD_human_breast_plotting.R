library(SpotSweeper)
library(RANN)
library(here)
library(SpatialExperiment)
library(spatialLIBD)
library(escheR)
library(patchwork)
library(scran)
library(scater)
library(tidyr)
library(dplyr)


plot_dir <- here("plots","VisiumHD","human_breast","016")


spe <- readRDS(here("processed-data", "VisiumHD","human_breast", "VisiumHD_HumanBreast_016_bspace_spotsweeper.rds"))
spe

# total local outliers
spe$local_outliers <- as.logical(spe$sum_outliers) |
  as.logical(spe$detected_outliers) |
  as.logical(spe$subsets_mito_percent_outliers)


colnames(colData(spe))
#  [1] "barcode"                        "in_tissue"                     
#  [3] "array_row"                      "array_col"                     
#  [5] "sample_id"                      "sum"                           
#  [7] "detected"                       "subsets_mito_percent"          
#  [9] "subsets_mito_sum"               "sum_discard"                   
# [11] "detected_discard"               "subsets_mito_percent_discard"  
# [13] "sizeFactor"                     "row"                           
# [15] "col"                            "BS_k10"                        
# [17] "detected_threshold"             "sum_threshold"                 
# [19] "subsets_mito_percent_threshold" "sum_log"                       
# [21] "sum_outliers"                   "sum_z"                         
# [23] "detected_log"                   "detected_outliers"             
# [25] "detected_z"                     "subsets_mito_percent_log"      
# [27] "subsets_mito_percent_outliers"  "subsets_mito_percent_z"        
# [29] "local_outliers"  

# 3MAD outliers discard
spe$MAD_outliers <- as.logical(spe$sum_discard) |
  as.logical(spe$detected_discard) |
  as.logical(spe$subsets_mito_percent_discard)

# total threshold outliers
spe$threshold_outliers <- as.logical(spe$sum_threshold) |
    as.logical(spe$detected_threshold) |
    as.logical(spe$subsets_mito_percent_threshold)



# named color palette of 1-10 for BS clusters
pal <- c("1" = "#FF0000", "2" = "#FFA500", "3" = "#FFFF00", "4" = "#008000", "5" = "#0000FF", "6" = "#4B0082", "7" = "#EE82EE", "8" = "#A52A2A", "9" = "#808080", "10" = "#000000")



# ========= Create dataframe of outliers by domain ========
outlier_df <- data.frame(MAD=spe$MAD_outliers,
                         Threshold=spe$threshold_outliers, 
                         SpotSweeper=spe$local_outliers,
                         domain=factor(spe.sw$BS_k10),
                         sum = spe$sum,
                         detected = spe$detected,
                         subsets_mito_percent = spe$subsets_mito_percent,
                         sum_z = spe$sum_z,
                         detected_z = spe$detected_z,
                         subsets_mito_percent_z = spe$subsets_mito_percent_z
                        )

# get the number of outliers per domain
outlier_df <- outlier_df %>%
    group_by(domain) %>%
    summarise(MAD=sum(MAD),
              SpotSweeper=sum(SpotSweeper),
              Threshold=sum(Threshold)) %>% 
    pivot_longer(cols=c(Threshold, MAD, SpotSweeper), names_to="method", values_to="count")

# change plotting order to Threhsold, Mad, SpotSweeper
outlier_df$method <- factor(outlier_df$method, levels=c("Threshold", "MAD", "SpotSweeper"))

# bar plot sw vs mad outliers per domain
png(here(plot_dir, "HD_HumanBreast_016_outliers_per_domain.png"), width=5, height=5, res=300, units="in")
ggplot(outlier_df, aes(x=domain, y=count, fill=method)) +
    geom_bar(stat="identity", position="dodge") +
    labs(x="Spatial Domain", y="Number of Outliers", fill="Method") +
    theme_classic() +
    scale_fill_manual(values=c("#459395", "#eb7c69","#fda638")) +
    theme(legend.position="top",
    legend.direction="horizontal") +
    theme(text = element_text(size=16))
dev.off()




# ============ Plotting here =============
png(here(plot_dir, "HD_HumanBreast_016_BayesSpace_domains.png"), width=5, height=5, units="in", res=300)
make_escheR(spe.sw) |>
    add_fill(var="BS_k10", point_size=.6) +
    scale_fill_manual(values=pal,
    name="") +
      theme(plot.title = element_text(size = 24, hjust = 0.1),
        legend.text = element_text(size=14)) +
    theme(legend.position="none") +
    ggtitle("Spatial Domains") 
dev.off()


png(here(plot_dir, "HD_HumanBreast_016_sum_3MAD.png"), width=5, height=5, res=300, units="in")
plotQC(spe.sw, metric="detected", outliers="sum_discard", point_size=0.6, stroke=.8) +
    ggtitle("3 MAD") +
  theme(plot.title = element_text(size = 20, hjust = 0.1),
        legend.text = element_text(size=11))
dev.off()

png(here(plot_dir, "HD_HumanBreast_016_sum_threshold.png"), width=5, height=5, res=300, units="in")
plotQC(spe.sw, metric="detected", outliers="sum_threshold", point_size=0.6, stroke=.8) +
    ggtitle("Fixed Threshold") +
  theme(plot.title = element_text(size = 20, hjust = 0.1),
        legend.text = element_text(size=11))
dev.off()

png(here(plot_dir, "HD_HumanBreast_016_sum_spotsweeper.png"), width=5, height=5, res=300, units="in")
plotQC(spe.sw, metric="detected", outliers="sum_outliers", point_size=0.6, stroke=.8) +
    ggtitle("SpotSweeper") +
  theme(plot.title = element_text(size = 20, hjust = 0.1),
        legend.text = element_text(size=11))
dev.off()



# ===== Threshold Ridge plots =====
library(ggridges)
# ridge plot of sum_umi with a threshold of 500
png(here(plot_dir, "HD_HumanBreast_016_ridge_sum_umi.png"), width=5, height=5, res=300, units="in")
ggplot(outlier_df, aes(x = sum, y = domain, fill = domain)) +
  geom_density_ridges(alpha = 0.8, scale = 1.5) +
  geom_vline(xintercept = 500, linetype = "dashed", color = "red", size=1) +
  scale_x_continuous(trans='log10') +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(size = 20),
        text = element_text(size = 15)) +
  scale_fill_manual(values = pal) +
  labs(title = "Library size per domain",
       x = "Library Size",
       y = "Spatial Domain") +
       coord_cartesian(xlim = c(50, NA))
dev.off()

png(here(plot_dir, "HD_HumanBreast_016_ridge_mito_percent.png"), width=5, height=5, res=300, units="in")
ggplot(outlier_df, aes(x = subsets_mito_percent, y = domain, fill = domain)) +
  geom_density_ridges(alpha = 0.8, scale = 1.5) +
  geom_vline(xintercept = 10, linetype = "dashed", color = "red", size=1) +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(size = 20),
        text = element_text(size = 15)) +
  scale_fill_manual(values = pal) +
  labs(title = "Mito percent per domain",
       x = "Mitochondrial Ratio",
       y = "Spatial Domain")+
  coord_cartesian(xlim = c(NA, 20))
dev.off()

png(here(plot_dir, "HD_HumanBreast_016_ridge_unique_genes.png"), width=5, height=5, res=300, units="in")
ggplot(outlier_df, aes(x = detected, y = domain, fill = domain)) +
  geom_density_ridges(alpha = 0.8, scale = 1.5) +
  geom_vline(xintercept = 500, linetype = "dashed", color = "red", size=1) +
  theme_bw() +
  scale_x_continuous(trans='log10') +
  theme(legend.position = "none",
        plot.title = element_text(size = 20),
        text = element_text(size = 15)) +
  scale_fill_manual(values = pal) +
  labs(title = "Unique genes per domain",
       x = "Number of genes",
       y = "Spatial Domain")+
       coord_cartesian(xlim = c(50, NA))
dev.off()




# ===== SpotSweeper ridge plots =====
# ridge plot of sum_umi with a threshold of 500
png(here(plot_dir, "HD_HumanBreast_016_ridge_sum_umi_spotsweeper.png"), width=5, height=5, res=300, units="in")
ggplot(outlier_df, aes(x = sum_z, y = domain, fill = domain)) +
  geom_density_ridges(alpha = 0.8, scale = 1.5) +
  geom_vline(xintercept = -3, linetype = "dashed", color = "red", size=1) +
  #scale_x_continuous(trans='log10') +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(size = 20),
        text = element_text(size = 15)) +
  scale_fill_manual(values = pal) +
  labs(title = "Library size (local z-score)",
       x = "Z-score",
       y = "Spatial Domain") 
dev.off()

png(here(plot_dir, "HD_HumanBreast_016_ridge_mito_percent_spotsweeper.png"), width=5, height=5, res=300, units="in")
ggplot(outlier_df, aes(x = subsets_mito_percent_z, y = domain, fill = domain)) +
  geom_density_ridges(alpha = 0.8, scale = 1.5) +
  geom_vline(xintercept = 3, linetype = "dashed", color = "red", size=1) +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(size = 20),
        text = element_text(size = 15)) +
  scale_fill_manual(values = pal) +
  labs(title = "Mito percent  (local z-score)",
       x = "Z-score",
       y = "Spatial Domain")
dev.off()

png(here(plot_dir, "HD_HumanBreast_016_ridge_unique_genes_spotsweeper.png"), width=5, height=5, res=300, units="in")
ggplot(outlier_df, aes(x = detected_z, y = domain, fill = domain)) +
  geom_density_ridges(alpha = 0.8, scale = 1.5) +
  geom_vline(xintercept = -3, linetype = "dashed", color = "red", size=1) +
  theme_bw() +
  #scale_x_continuous(trans='log10') +
  theme(legend.position = "none",
        plot.title = element_text(size = 20),
        text = element_text(size = 15)) +
  scale_fill_manual(values = pal) +
  labs(title = "Unique genes (local z-score)",
       x = "Z-score",
       y = "Spatial Domain")
dev.off()





# ====== UMAP ======
set.seed(555)

# run UMAP
spe <- runUMAP(spe, dimred = "PCA")

# plot UMAP
p1 <- plotReducedDim(spe, "UMAP", colour_by = "BS_k10", point_size=0.2) +
  scale_color_manual(values = pal) +
  labs(title = "Clusters") +
  theme(plot.title = element_text(size = 20))

p2 <- plotReducedDim(spe, "UMAP", colour_by = "MAD_outliers", point_size=0.2)+
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
  labs(title = "Global Outliers (3 MADs)") +
  theme(plot.title = element_text(size = 20))

p3 <- plotReducedDim(spe, "UMAP", colour_by = "local_outliers", point_size=0.2) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
  labs(title = "Local Outliers (3 Zs)") +
  theme(plot.title = element_text(size = 20))

p4 <- plotReducedDim(spe, "UMAP", colour_by = "subsets_mito_percent", point_size=0.2) +
  labs(title = "Mitochondrial Ratio") +
  theme(plot.title = element_text(size = 20))

p5 <- plotReducedDim(spe, "UMAP", colour_by = "sum", point_size=0.2) +
  labs(title = "Library Size") +
  theme(plot.title = element_text(size = 20))

p6 <- plotReducedDim(spe, "UMAP", colour_by = "detected", point_size=0.2) +
  labs(title = "Unique Genes") +
  theme(plot.title = element_text(size = 20))

png(here(plot_dir, 'HD_HumanBreast_016_umap.png'), width=15, height=10, units="in", res=300)
(p1+p2+p3)/(p4+p5+p6)
dev.off()





