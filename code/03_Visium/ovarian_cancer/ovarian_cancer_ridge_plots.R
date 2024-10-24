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
library(tidyr)


plot_dir <- here("plots","Visium", "ovarian_cancer")


spe <- readRDS(here("processed-data", "Visium","ovarian_cancer","OvarianCancer_banksy_spotsweeper.rds"))
spe


colnames(colData(spe))
#  [1] "in_tissue"                      "array_row"                     
#  [3] "array_col"                      "sample_id"                     
#  [5] "sizeFactor"                     "clust_M0_lam0.8_k50_res0.8"    
#  [7] "sum"                            "detected"                      
#  [9] "subsets_mito_sum"               "subsets_mito_detected"         
# [11] "subsets_mito_percent"           "total"                         
# [13] "sum_discard"                    "detected_discard"              
# [15] "subsets_mito_percent_discard"   "detected_threshold"            
# [17] "sum_threshold"                  "subsets_mito_percent_threshold"
# [19] "prob_compromised"               "miQC_keep"                     
# [21] "sum_log"                        "sum_outliers"                  
# [23] "sum_z"                          "detected_log"                  
# [25] "detected_outliers"              "detected_z"                    
# [27] "subsets_mito_percent_log"       "subsets_mito_percent_outliers" 
# [29] "subsets_mito_percent_z"        

# ===== combining outliers form the different methods ======

# total local outliers
spe$local_outliers <- as.logical(spe$sum_outliers) |
  as.logical(spe$detected_outliers) |
  as.logical(spe$subsets_mito_percent_outliers)


# 3MAD outliers discard
spe$MAD_outliers <- as.logical(spe$sum_discard) |
  as.logical(spe$detected_discard) |
  as.logical(spe$subsets_mito_percent_discard)

# total threshold outliers
spe$threshold_outliers <- as.logical(spe$sum_threshold) |
    as.logical(spe$detected_threshold) |
    as.logical(spe$subsets_mito_percent_threshold)




# ========= Create dataframe of outliers by domain ========
outlier_df <- data.frame(MAD=spe$MAD_outliers,
                         miQC=spe$miQC_keep,
                         Threshold=spe$threshold_outliers, 
                         SpotSweeper=spe$local_outliers,
                         domain=factor(spe$clust_M0_lam0.2_k50_res0.8),
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
              Threshold=sum(Threshold),
              miQC=sum(miQC)) %>% 
    pivot_longer(cols=c(Threshold, MAD, SpotSweeper, miQC), names_to="method", values_to="count")

# change plotting order to Threhsold, Mad, SpotSweeper
outlier_df$method <- factor(outlier_df$method, levels=c( "miQC", "Threshold", "MAD", "SpotSweeper"))

# bar plot sw vs mad outliers per domain
png(here(plot_dir, "overian_cancer_outliers_per_domain.png"), width=4.5, height=5, res=300, units="in")
ggplot(outlier_df, aes(x=domain, y=count, fill=method)) +
    geom_bar(stat="identity", position="dodge") +
    labs(x="Spatial Domain", y="Number of Outliers", fill="Method") +
    theme_classic() +
    scale_fill_manual(values=c("#459395", "#eb7c69","#fda638", "#b82ac9"),
                      guide = guide_legend(nrow=2)) +
    theme(legend.position="top",
    legend.direction="horizontal") +
    theme(text = element_text(size=16, face="bold")) 
dev.off()




# ===== Threshold Ridge plots =====

# Extract the sum values that are marked as outliers
outliers <- spe$sum[spe$sum_discard]
sum_3MAD <- min(outliers)
sum_3MAD

outliers <- spe$subsets_mito_percent[spe$subsets_mito_percent_discard]
mito_3MAD <- min(outliers)
mito_3MAD

outliers <- spe$detected[spe$detected_discard]
detected_3MAD <- min(outliers)
detected_3MAD

# make color palette for number of domains
n_domains <- length(unique(colData(spe)$clust_M0_lam0.2_k50_res0.8))
domain_pal <- colorRampPalette(RColorBrewer::brewer.pal(n_domains, "Set1"))
pal <- domain_pal(n_domains)

library(ggridges)
# ridge plot of sum_umi with a threshold of 500
png(here(plot_dir, "ovarian_cancer_ridge_sum_umi.png"), width=4, height=4, res=300, units="in")
ggplot(outlier_df, aes(x = sum, y = domain, fill = domain)) +
  geom_density_ridges(alpha = 0.8, scale = 1.5) +
  geom_vline(xintercept = 500, linetype = "dashed", color = "#eb7c69", size=1.75) +
  geom_vline(xintercept = sum_3MAD, linetype = "dashed", color = "#fda638", size=1.75) +
  scale_x_continuous(trans='log10') +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(size = 20),
        text = element_text(size = 16, face = "bold")) +
  scale_fill_manual(values = pal) +
  labs(title = "Library size per domain",
       x = "Library Size",
       y = "Spatial Domain")
dev.off()

png(here(plot_dir, "ovarian_cancer_ridge_mito_percent.png"), width=4, height=4, res=300, units="in")
ggplot(outlier_df, aes(x = subsets_mito_percent, y = domain, fill = domain)) +
  geom_density_ridges(alpha = 0.8, scale = 1.5) +
  geom_vline(xintercept = 10, linetype = "dashed", color = "#eb7c69", size=1.75) +
    geom_vline(xintercept = mito_3MAD, linetype = "dashed", color = "#fda638", size=1.75) +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(size = 18),
        text = element_text(size = 16, face = "bold")) +
  scale_fill_manual(values = pal) +
  labs(title = "Mito percent per domain",
       x = "Mitochondrial Ratio",
       y = "Spatial Domain")
  #coord_cartesian(xlim = c(NA, 2))
dev.off()

png(here(plot_dir, "ovarian_cancer_ridge_unique_genes.png"), width=4, height=4, res=300, units="in")
ggplot(outlier_df, aes(x = detected, y = domain, fill = domain)) +
  geom_density_ridges(alpha = 0.8, scale = 1.5) +
  geom_vline(xintercept = 500, linetype = "dashed", color = "#eb7c69", size=1.75) +
    geom_vline(xintercept = detected_3MAD, linetype = "dashed", color = "#fda638", size=1.75) +
  theme_bw() +
  scale_x_continuous(trans='log10') +
  theme(legend.position = "none",
        plot.title = element_text(size = 18),
        text = element_text(size = 16, face = "bold")) +
  scale_fill_manual(values = pal) +
  labs(title = "Unique genes per domain",
       x = "Number of genes",
       y = "Spatial Domain")
       #coord_cartesian(xlim = c(50, NA))
dev.off()




# ===== SpotSweeper ridge plots =====
# ridge plot of sum_umi with a threshold of 500
png(here(plot_dir, "ovarian_cancer_ridge_sum_umi_spotsweeper.png"), width=4, height=4, res=300, units="in")
ggplot(outlier_df, aes(x = sum_z, y = domain, fill = domain)) +
  geom_density_ridges(alpha = 0.8, scale = 1.5) +
  geom_vline(xintercept = -3, linetype = "dashed", color = "#b82ac9", size=1.75) +
  #scale_x_continuous(trans='log10') +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(size = 19),
        text = element_text(size = 16, face = "bold")) +
  scale_fill_manual(values = pal) +
  labs(title = "Library size (local z-score)",
       x = "Z-score",
       y = "Spatial Domain") 
dev.off()

png(here(plot_dir, "ovarian_cancer_ridge_mito_percent_spotsweeper.png"), width=4, height=4, res=300, units="in")
ggplot(outlier_df, aes(x = subsets_mito_percent_z, y = domain, fill = domain)) +
  geom_density_ridges(alpha = 0.8, scale = 1.5) +
  geom_vline(xintercept = 3, linetype = "dashed", color = "#b82ac9", size=1.75) +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(size = 18),
        text = element_text(size = 16, face = "bold")) +
  scale_fill_manual(values = pal) +
  labs(title = "Mito percent  (local z-score)",
       x = "Z-score",
       y = "Spatial Domain")
dev.off()

png(here(plot_dir, "ovarian_cancer_ridge_unique_genes_spotsweeper.png"), width=4, height=4, res=300, units="in")
ggplot(outlier_df, aes(x = detected_z, y = domain, fill = domain)) +
  geom_density_ridges(alpha = 0.8, scale = 1.5) +
  geom_vline(xintercept = -3, linetype = "dashed", color = "#b82ac9", size=1.75) +
  theme_bw() +
  #scale_x_continuous(trans='log10') +
  theme(legend.position = "none",
        plot.title = element_text(size = 18),
        text = element_text(size = 16, face = "bold")) +
  scale_fill_manual(values = pal) +
  labs(title = "Unique genes (local z-score)",
       x = "Z-score",
       y = "Spatial Domain")
dev.off()

