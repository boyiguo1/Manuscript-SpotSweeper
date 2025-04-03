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
library(RColorBrewer)
library(BayesSpace)

plot_dir = here('plots',"Visium", "figure_4_5")
processed_dir = here('raw-data')


# large DLPFC dataset
load(here("processed-data","Visium","dlPFC_raw.Rdata"))
spe.pfc <- spe_raw
spe.pfc




# drop out of tissue spots
spe.pfc <- spe.pfc[,colData(spe.pfc)$in_tissue == TRUE]


# ========= Artifact examples =========

# Define the function for plotting hitology
plotHistology <- function(spe, sample_id) {
  spe_subset <- spe[,colData(spe)$sample_id == sample_id]

  # Make dataframe combining column data and spatial coordinates
  df <- as.data.frame(cbind(colData(spe_subset), spatialCoords(spe_subset)), optional = TRUE)

  # Plotting
  p <- vis_gene_p(
    spe = spe_subset,
    d = df,
    title = sample_id,
    spatial = TRUE
  )

  # clean up and return
  rm(spe_subset)
  return(p)
}

# ==== histology ====

pdf(here(plot_dir,"dryspot_histology.pdf"), width=5, height=5)
p1 <- plotHistology(spe.pfc, "Br3942_mid")
p1
dev.off()

pdf(here(plot_dir,"hangnail_histology.pdf"), width=5, height=5)
p2 <- plotHistology(spe.pfc, "Br8325_ant")
p2
dev.off()

# ==== library size ====

spe.pfc$sum_umi_log <- log2(spe.pfc$sum_umi)
spe.pfc$sum_gene_log <- log2(spe.pfc$sum_gene)

spe.dryspot <- spe.pfc[,colData(spe.pfc)$sample_id == "Br3942_mid"]
spe.hangnail <- spe.pfc[,colData(spe.pfc)$sample_id == "Br8325_ant"]

pdf(here(plot_dir, "dryspot_library_size.pdf"), width=5, height=5)
p3 <- SpotSweeper::plotQC(spe.dryspot, metric="sum_umi_log",point_size=2.1) +
  ggtitle("Library size") +
  theme(text=element_text(size=18),
        plot.title = element_text(size=24),
        legend.position = "bottom",
        legend.title = element_text(hjust=0.5),
        legend.text=element_text(size=14),
        legend.key.width = unit(1, "cm"))
p3
dev.off()

pdf(here(plot_dir,"hangnail_library_size.pdf"), width=5, height=5)
p4 <- SpotSweeper::plotQC(spe.hangnail, metric="sum_umi_log",point_size=2.1) +
  ggtitle("Library size") +
  theme(text=element_text(size=18),
        plot.title = element_text(size=24),
        legend.position = "bottom",
        legend.title = element_text(hjust=0.5),
        legend.text=element_text(size=14),
        legend.key.width = unit(1, "cm"))
p4
dev.off()

# make dataframe of library size with spatalCoords an export to csv
hangnail_df <- data.frame(
  sample_id = colData(spe.hangnail)$sample_id,
  sum_umi = colData(spe.hangnail)$sum_umi,
  row = colData(spe.hangnail)$array_row,
  col = colData(spe.hangnail)$array_col
)

dryspot_df <- data.frame(
  sample_id = colData(spe.dryspot)$sample_id,
  sum_umi = colData(spe.dryspot)$sum_umi,
  row = colData(spe.dryspot)$array_row,
  col = colData(spe.dryspot)$array_col
)

# Export to CSV
write.csv(dryspot_df, file=here("processed-data", "outputs_for_paper", "figure_5", "Figure5_D.csv"), row.names = FALSE)
write.csv(hangnail_df, file=here("processed-data", "outputs_for_paper", "figure_5", "Figure5_I.csv"), row.names = FALSE)

# ==== unique genes ====

pdf(here(plot_dir, "dryspot_unique_genes.pdf"), width=5, height=5)
p3 <- SpotSweeper::plotQC(spe.dryspot, metric="sum_gene_log",point_size=2.1) +
  ggtitle("Unique Genes") +
  theme(text=element_text(size=18),
        plot.title = element_text(size=24),
        legend.position = "bottom",
        legend.title = element_text(hjust=0.5),
        legend.text=element_text(size=14),
        legend.key.width = unit(1, "cm"))
p3
dev.off()

pdf(here(plot_dir,"hangnail_unique_genes.pdf"), width=5, height=5)
p4 <- SpotSweeper::plotQC(spe.hangnail, metric="sum_gene_log",point_size=2.1) +
  ggtitle("Unique Genes") +
  theme(text=element_text(size=18),
        plot.title = element_text(size=24),
        legend.position = "bottom",
        legend.title = element_text(hjust=0.5),
        legend.text=element_text(size=14),
        legend.key.width = unit(1, "cm"))
p4
dev.off()


# ==== mito percent ====

pdf(here(plot_dir, "dryspot_mito_ratio.pdf"), width=5, height=5)
p3 <- SpotSweeper::plotQC(spe.dryspot, metric="expr_chrM_ratio",point_size=2.1) +
  ggtitle("Mito Ratio") +
  theme(text=element_text(size=18),
        plot.title = element_text(size=24),
        legend.position = "bottom",
        legend.title = element_text(hjust=0.5),
        legend.text=element_text(size=14),
        legend.key.width = unit(1, "cm"))

p3
dev.off()

pdf(here(plot_dir,"hangnail_mito_ratio.pdf"), width=5, height=5)
p4 <- SpotSweeper::plotQC(spe.hangnail, metric="expr_chrM_ratio", point_size=2.1) +
  ggtitle("Mito Ratio") +
  theme(text=element_text(size=18),
        plot.title = element_text(size=24),
        legend.position = "bottom",
        legend.title = element_text(hjust=0.5),
        legend.text=element_text(size=14),
        legend.key.width = unit(1, "cm"))
p4
dev.off()


# make dataframe of library size with spatalCoords an export to csv
hangnail_df <- data.frame(
  sample_id = colData(spe.hangnail)$sample_id,
  mito_ratio = colData(spe.hangnail)$expr_chrM_ratio,
  row = colData(spe.hangnail)$array_row,
  col = colData(spe.hangnail)$array_col
)

dryspot_df <- data.frame(
  sample_id = colData(spe.dryspot)$sample_id,
  mito_ratio = colData(spe.dryspot)$expr_chrM_ratio,
  row = colData(spe.dryspot)$array_row,
  col = colData(spe.dryspot)$array_col
)

# Export to CSV
write.csv(dryspot_df, file=here("processed-data", "outputs_for_paper", "figure_5", "Figure5_E.csv"), row.names = FALSE)
write.csv(hangnail_df, file=here("processed-data", "outputs_for_paper", "figure_5", "Figure5_J.csv"), row.names = FALSE)

# ========================================
#             Clustering
# ========================================


# ==== improper clustering ====

library(ggthemes)
library(scales)
colors<- colorblind_pal()(7)


# ===== BS clustering of the dryspot samples ======
spe.dryspot <- logNormCounts(spe.dryspot)

set.seed(101)
dec <- scran::modelGeneVar(spe.dryspot)
top <- scran::getTopHVGs(dec, n = 3000)

set.seed(102)
spe.dryspot <- scater::runPCA(spe.dryspot, subset_row=top)

## Add BayesSpace metadata
spe.dryspot <- spatialPreprocess(spe.dryspot, platform="Visium", skip.PCA=TRUE)

colData(spe.dryspot)$row <- spe.dryspot$array_row
colData(spe.dryspot)$col <- spe.dryspot$array_col

metadata(spe.dryspot)$BayesSpace.data <- list(platform = "Visium", is.enhanced = FALSE)

d <- 15  # Number of PCs
[
## Run BayesSpace clustering
set.seed(104)
dlpfc <- spatialCluster(spe.dryspot, q=2, d=d, platform='Visium',
                        nrep=10000, gamma=3, save.chain=TRUE)
spe.dryspot$BS_k2 <- as.factor(dlpfc$spatial.cluster)

dlpfc <- spatialCluster(spe.dryspot, q=7, d=d, platform='Visium',
                        nrep=10000, gamma=3, save.chain=TRUE)
spe.dryspot$BS_k7 <- as.factor(dlpfc$spatial.cluster)

library(ggsci)]


pdf(here(plot_dir,"dryspot_clustering.pdf"), width=5, height=5)
p <- make_escheR(spe.dryspot) |>
  add_fill("BS_k7") +
  ggtitle("BayesSpace k=7") +
  theme(text=element_text(size=18),
        plot.title = element_text(size=24),
        legend.position = "bottom",
        legend.text=element_text(size=14)) +
  scale_fill_manual(values=colors)
p
dev.off()

# to csv
dryspot_df <- data.frame(
  sample_id = colData(spe.dryspot)$sample_id,
  BS_k7 = as.character(colData(spe.dryspot)$BS_k7),
  row = colData(spe.dryspot)$array_row,
  col = colData(spe.dryspot)$array_col
)

# Export to CSV
write.csv(dryspot_df, file=here("processed-data", "outputs_for_paper", "figure_5", "Figure5_F.csv"), row.names = FALSE)


# ===== BS clustering of the hangnail samples ======
spe.hangnail <- logNormCounts(spe.hangnail)

set.seed(101)
dec <- scran::modelGeneVar(spe.hangnail)
top <- scran::getTopHVGs(dec, n = 3000)

set.seed(102)
spe.hangnail <- scater::runPCA(spe.hangnail, subset_row=top)

## Add BayesSpace metadata
spe.hangnail <- spatialPreprocess(spe.hangnail, platform="Visium", skip.PCA=TRUE)

colData(spe.hangnail)$row <- spe.hangnail$array_row
colData(spe.hangnail)$col <- spe.hangnail$array_col

metadata(spe.hangnail)$BayesSpace.data <- list(platform = "Visium", is.enhanced = FALSE)

d <- 15  # Number of PCs

## Run BayesSpace clustering
set.seed(104)
dlpfc <- spatialCluster(spe.hangnail, q=2, d=d, platform='Visium',
                        nrep=10000, gamma=3, save.chain=TRUE)

spe.hangnail$BS_k2 <- as.factor(dlpfc$spatial.cluster)

dlpfc <- spatialCluster(spe.hangnail, q=7, d=d, platform='Visium',
                        nrep=10000, gamma=3, save.chain=TRUE)

spe.hangnail$BS_k7 <- as.factor(dlpfc$spatial.cluster)

pdf(here(plot_dir,"hangnail_BayesSpace_k7.pdf"), width=5, height=5)
p1 <- make_escheR(spe.hangnail) |>
  add_fill("BS_k7") +
  ggtitle("BayesSpace k=7") +
  theme(text=element_text(size=18),
        plot.title = element_text(size=24),
        legend.position = "bottom",
        legend.text=element_text(size=14)) +
  #change colors
  scale_fill_manual(values=colors)

p1
dev.off()

# to csv
hangnail_df <- data.frame(
  sample_id = colData(spe.hangnail)$sample_id,
  BS_k7 = as.character(colData(spe.hangnail)$BS_k7),
  row = colData(spe.hangnail)$array_row,
  col = colData(spe.hangnail)$array_col
)

# Export to CSV
write.csv(hangnail_df, file=here("processed-data", "outputs_for_paper", "figure_5", "Figure5_K.csv"), row.names = FALSE)


# ========= Dry-spot artifact ==========

library(dbscan)

x <- data.frame(umi=log2(spe.dryspot$sum_umi),
               gene=log2(spe.dryspot$sum_gene))

plot(x)


outs <- dbscan(x, eps = .5, minPts = 20)
spe.dryspot$DBSCAN_clusters <- as.factor(outs$cluster)


clus1_mean <- mean(colData(spe.dryspot)[["sum_umi_log"]][spe.dryspot$DBSCAN_clusters == 1])
clus2_mean <- mean(colData(spe.dryspot)[["sum_umi_log"]][spe.dryspot$DBSCAN_clusters == 2])

artifact_clus <- which.min(c(clus1_mean, clus2_mean))
spe.dryspot$DBSCAN_artifact <- FALSE
spe.dryspot$DBSCAN_artifact[spe.dryspot$DBSCAN_clusters == artifact_clus] <- TRUE


p1 <- make_escheR(spe.dryspot) |>
  add_fill(var="DBSCAN_artifact") +
  ggtitle("Detected artifact") +
  theme(text=element_text(size=18),
        plot.title = element_text(size=24),
        legend.position = "bottom",
        legend.text=element_text(size=14))  +
  scale_fill_manual(values=c("grey", "red"),
                    name="Discard")

pdf(here(plot_dir,"dryspot_dbscan_spotplots.pdf"), width=5, height=5)
p1
dev.off()


df <- data.frame(umi=log2(spe.dryspot$sum_umi),
                gene=log2(spe.dryspot$sum_gene),
                cluster=spe.dryspot$DBSCAN_artifact)

#ggplot of x colored by clusters
pdf(here(plot_dir,"dryspot_dbscan_scatter.pdf"), width=6, height=5)
p <- ggplot(df, aes(x=umi, y=gene, color=cluster)) +
  geom_point(alpha=.2) +
  scale_color_manual(values=c("grey", "red"),
                     name="Discard")  +
  theme_bw() +
  theme(text=element_text(size=18),
        plot.title = element_text(size=24),
        #legend.position = "bottom",
        legend.text=element_text(size=14)) +
  ggtitle("DBSCAN clustering") +
  ylab("log2(Total UMI)") +
  xlab("log2(Unique genes)")
p
dev.off()

# to csv
write.csv(df, file=here("processed-data", "outputs_for_paper", "figure_6", "Figure6_B.csv"), row.names = FALSE)

# add spatialCoords and export
dryspot_df <- data.frame(
  sample_id = colData(spe.dryspot)$sample_id,
  DBSCAN_artifact = as.character(colData(spe.dryspot)$DBSCAN_artifact),
  row = colData(spe.dryspot)$array_row,
  col = colData(spe.dryspot)$array_col
)

# Export to CSV
write.csv(dryspot_df, file=here("processed-data", "outputs_for_paper", "figure_6", "Figure6_C.csv"), row.names = FALSE)


# ========== Hangnail artifact ============

spe.hangnail <- findArtifacts(spe.hangnail,
                            mito_percent="expr_chrM_ratio",
                            mito_sum="expr_chrM",
                            n_rings=5,
                            name="artifact"
)


# ======= Testing PCA clustering vs clustering on single local variance ========
library(mclust)
#
# pdf(here(plot_dir,"hangnail_localvariance.pdf"), width=5, height=5, units="in", res=300)
# p1 <- make_escheR(spe.hangnail) |>
#   add_fill("k18", point_size=2.1) +
#   ggtitle("Local variance") +
#   theme(text=element_text(size=18),
#         plot.title = element_text(size=24),
#         legend.position = "bottom",
#         legend.text=element_text(size=14)) +
#   scale_fill_gradient(low = "white", high = "black")
# p1
# dev.off()



spe.hangnail$artifact_pc1 <- reducedDim(spe.hangnail, "PCA_artifacts")[,1]
pdf(here(plot_dir,"hangnail_localvariance.pdf"), width=5, height=5)
p1 <- make_escheR(spe.hangnail) |>
  add_fill("artifact_pc1", point_size=2.1) +
  ggtitle("Multiscale variance") +
  theme(text=element_text(size=18),
        plot.title = element_text(size=24),
        legend.position = "bottom",
        legend.text=element_text(size=14)) +
  scale_fill_gradient(low = "white", high = "black", name="PC1")
p1
dev.off()

# to csv
hangnail_df <- data.frame(
  sample_id = colData(spe.hangnail)$sample_id,
  artifact_pc1 = colData(spe.hangnail)$artifact_pc1,
  row = colData(spe.hangnail)$array_row,
  col = colData(spe.hangnail)$array_col
)

# Export to CSV
write.csv(hangnail_df, file=here("processed-data", "outputs_for_paper", "figure_6", "Figure6_E.csv"), row.names = FALSE)



# local variance plus mito
y <- data.frame(variance=spe.hangnail$k18,
               mito=log2(spe.hangnail$expr_chrM)
                )


# GMM
clust <- Mclust(y, G=2)

spe.hangnail$Mclust_clusters <- as.factor(clust$classification)

clus1_mean <- mean(colData(spe.hangnail)[[paste0("k",18
)]][spe.hangnail$Mclust_clusters == 1])
clus2_mean <- mean(colData(spe.hangnail)[[paste0("k",18
)]][spe.hangnail$Mclust_clusters == 2])

artifact_clus <- which.min(c(clus1_mean, clus2_mean))
spe.hangnail$Mclust_artifact <- FALSE
spe.hangnail$Mclust_artifact[spe.hangnail$Mclust_clusters == artifact_clus] <- TRUE

p1 <- make_escheR(spe.hangnail) |>
  add_fill("Mclust_artifact", point_size=2.3) +
  #ggtitle("Mclust clusters (var)") +
  theme(text=element_text(size=18),
        plot.title = element_text(size=24),
        legend.position = "bottom",
        legend.text=element_text(size=14)) +
  scale_fill_manual(values=c("grey", "red"), name="Discard")


pdf(here(plot_dir,"hangnail_mclust_variance.pdf"), width=5, height=5)
p1
dev.off()

# dbscan
outs_hangnail <- dbscan(y, eps = .2, minPts = 10)

spe.hangnail$DBSCAN_clusters <- as.factor(outs_hangnail$cluster)

clus1_mean <- mean(colData(spe.hangnail)[[paste0("k",18
)]][spe.hangnail$DBSCAN_clusters == 1])
clus2_mean <- mean(colData(spe.hangnail)[[paste0("k",18
)]][spe.hangnail$DBSCAN_clusters == 2])

artifact_clus <- which.min(c(clus1_mean, clus2_mean))
spe.hangnail$DBSCAN_artifact <- FALSE
spe.hangnail$DBSCAN_artifact[spe.hangnail$DBSCAN_clusters == artifact_clus] <- TRUE

p2 <- make_escheR(spe.hangnail) |>
  add_fill(var="DBSCAN_artifact", point_size=2.3) +
  #ggtitle("DBSCAN clusters (var)") +
  theme(text=element_text(size=18),
        plot.title = element_text(size=24),
        legend.position = "bottom",
        legend.text=element_text(size=14)) +
  scale_fill_manual(values=c("grey", "red"), name="Discard")


pdf(here(plot_dir,"hangnail_dbscan_variance.pdf"), width=5, height=5)
p2
dev.off()

# kmeans
clust <- kmeans(y, centers=2)

spe.hangnail$Kmeans_clusters <- as.factor(clust$cluster)

clus1_mean <- mean(colData(spe.hangnail)[[paste0("k",18
)]][spe.hangnail$Kmeans_clusters == 1])
clus2_mean <- mean(colData(spe.hangnail)[[paste0("k",18
)]][spe.hangnail$Kmeans_clusters == 2])

artifact_clus <- which.min(c(clus1_mean, clus2_mean))
spe.hangnail$Kmeans_artifact <- FALSE
spe.hangnail$Kmeans_artifact[spe.hangnail$Kmeans_clusters == artifact_clus] <- TRUE

p3 <- make_escheR(spe.hangnail) |>
  add_fill("Kmeans_artifact", point_size=2.3) +
  #ggtitle("Kmeans clusters (var)") +
  theme(text=element_text(size=18),
        plot.title = element_text(size=24),
        legend.position = "bottom",
        legend.text=element_text(size=14)) +
  scale_fill_manual(values=c("grey", "red"), name="Discard")


pdf(here(plot_dir,"hangnail_kmeans_variance.pdf"), width=5, height=5)
p3
dev.off()



# ========== various clustering on PCA of local variances ===========

pdf(here(plot_dir,"hangnail_artifact.pdf"), width=5, height=5)
p1 <- make_escheR(spe.hangnail) |>
  add_fill("artifact", name=dis) +
  ggtitle("Detected artifact") +
  theme(text=element_text(size=18),
        plot.title = element_text(size=24),
        legend.position = "bottom",
        legend.text=element_text(size=14)) +
  scale_fill_manual(values=c("grey", "red"), name="Discard")
p1
dev.off()

# to csv
hangnail_df <- data.frame(
  sample_id = colData(spe.hangnail)$sample_id,
  artifact = as.character(colData(spe.hangnail)$artifact),
  row = colData(spe.hangnail)$array_row,
  col = colData(spe.hangnail)$array_col
)
# Export to CSV
write.csv(hangnail_df, file=here("processed-data", "outputs_for_paper", "figure_6", "Figure6_F.csv"), row.names = FALSE)




x_pca <- reducedDim(spe.hangnail,"PCA_artifacts")[,1:2]

pdf(here(plot_dir,"hangnail_pca_scatter.pdf"), width=6, height=5)
ggplot(x_pca, aes(x=PC1, y=PC2)) +
  geom_point(aes(color=as.factor(spe.hangnail$artifact)), alpha=.5) +
  scale_color_manual(values=c("grey", "red"),
                     name="Discard") +
  theme_bw() +
  theme(text=element_text(size=18),
        plot.title = element_text(size=22),
        #legend.position = "bottom",
        legend.text=element_text(size=14)) +
  ggtitle("K-means clustering") +
  ylab("PC2") +
  xlab("PC1")
dev.off()

# to csv
hangnail_pca_df <- data.frame(
  sample_id = colData(spe.hangnail)$sample_id,
  PC1 = x_pca[,1],
  PC2 = x_pca[,2],
  artifact = as.character(colData(spe.hangnail)$artifact),
  row = colData(spe.hangnail)$array_row,
  col = colData(spe.hangnail)$array_col
)

# Export to CSV
write.csv(hangnail_pca_df, file=here("processed-data", "outputs_for_paper", "figure_6", "Figure6_G.csv"), row.names = FALSE)



# GMM
clust <- Mclust(x_pca, G=2)

spe.hangnail$Mclust_clusters <- as.factor(clust$classification)

clus1_mean <- mean(colData(spe.hangnail)[[paste0("k",18
)]][spe.hangnail$Mclust_clusters == 1])
clus2_mean <- mean(colData(spe.hangnail)[[paste0("k",18
)]][spe.hangnail$Mclust_clusters == 2])

artifact_clus <- which.min(c(clus1_mean, clus2_mean))
spe.hangnail$Mclust_artifact <- FALSE
spe.hangnail$Mclust_artifact[spe.hangnail$Mclust_clusters == artifact_clus] <- TRUE

p1 <- make_escheR(spe.hangnail) |>
  add_fill("Mclust_artifact", point_size=2.3) +
  #ggtitle("Mclust clusters (PCA)") +
  theme(text=element_text(size=18),
        plot.title = element_text(size=24),
        legend.position = "bottom",
        legend.text=element_text(size=14)) +
  scale_fill_manual(values=c("grey", "red"), name="Discard")

pdf(here(plot_dir,"hangnail_mclust_pca.pdf"), width=5, height=5)
p1
dev.off()

# dbscan
outs_hangnail <- dbscan(x_pca, eps = .5, minPts = 20)

spe.hangnail$DBSCAN_clusters <- as.factor(outs_hangnail$cluster)

clus1_mean <- mean(colData(spe.hangnail)[[paste0("k",18
)]][spe.hangnail$DBSCAN_clusters == 1])
clus2_mean <- mean(colData(spe.hangnail)[[paste0("k",18
)]][spe.hangnail$DBSCAN_clusters == 2])

artifact_clus <- which.min(c(clus1_mean, clus2_mean))
spe.hangnail$DBSCAN_artifact <- FALSE
spe.hangnail$DBSCAN_artifact[spe.hangnail$DBSCAN_clusters == artifact_clus] <- TRUE

p2 <- make_escheR(spe.hangnail) |>
  add_fill(var="DBSCAN_artifact", point_size=2.3) +
  #ggtitle("DBSCAN clusters (PCA)") +
  theme(text=element_text(size=18),
        plot.title = element_text(size=24),
        legend.position = "bottom",
        legend.text=element_text(size=14)) +
  scale_fill_manual(values=c("grey", "red"), name="Discard")



pdf(here(plot_dir,"hangnail_dbscan_pca.pdf"), width=5, height=5)
p2
dev.off()


# kmeans
clust <- kmeans(x_pca, centers=2)
spe.hangnail$Kmeans_clusters <- as.factor(clust$cluster)

clus1_mean <- mean(colData(spe.hangnail)[[paste0("k",18
)]][spe.hangnail$Kmeans_clusters == 1])
clus2_mean <- mean(colData(spe.hangnail)[[paste0("k",18
)]][spe.hangnail$Kmeans_clusters == 2])

artifact_clus <- which.min(c(clus1_mean, clus2_mean))
spe.hangnail$Kmeans_artifact <- FALSE
spe.hangnail$Kmeans_artifact[spe.hangnail$Kmeans_clusters == artifact_clus] <- TRUE

p3 <- make_escheR(spe.hangnail) |>
  add_fill("Kmeans_artifact", point_size=2.3) +
  #ggtitle("Kmeans clusters (PCA)") +
  theme(text=element_text(size=18),
        plot.title = element_text(size=24),
        legend.position = "bottom",
        legend.text=element_text(size=14)) +
  scale_fill_manual(values=c("grey", "red"), name="Discard")

pdf(here(plot_dir,"hangnail_kmeans_pca.pdf"), width=5, height=5)
p3
dev.off()




# ============= Testing layer marker genes after QC  =============

library(viridis)

rownames(spe.hangnail) <- rowData(spe.hangnail)$gene_name
rownames(spe.dryspot) <- rowData(spe.dryspot)$gene_name


known_L6_markers <- c( "NR4A2", "CCK", "KRT17", "SCN3B", "SYTB", "SEMA3E")
# testing layer markers with and without discarding artifact
spe.no.hangnail <- spe.hangnail[,colData(spe.hangnail)$artifact == FALSE]

#find markers
markers <- scran::findMarkers(spe.no.hangnail, groups=spe.no.hangnail$BS_k7, test.type="wilcox", direction="up")
L6_markers <- markers[[5]]

# ranking of known markers
L6_markers <- L6_markers[rownames(L6_markers) %in% known_L6_markers,]

no.hangnail_ranks <- L6_markers["Top"]
no.hangnail_ranks
# CCK            1
# SCN3B          5
# KRT17          6
# SEMA3E        10
# NR4A2         16

#find markers
markers <- scran::findMarkers(spe.hangnail, groups=spe.hangnail$BS_k7, test.type="wilcox", direction="up")
L6_markers <- markers[[5]]

# ranking of known markers
L6_markers <- L6_markers[rownames(L6_markers) %in% known_L6_markers,]

hangnail_ranks <- L6_markers["Top"]
hangnail_ranks
# CCK            1
# KRT17          9
# SCN3B         10
# SEMA3E       122
# NR4A2        178

hangnail_ranks$artifact <- "With"
no.hangnail_ranks$artifact <- "Without"

# ggplot boxplot
df <- rbind(hangnail_ranks, no.hangnail_ranks)

# to csv
write.csv(df, file=here("processed-data", "outputs_for_paper", "figure_6", "Figure6_H.csv"), row.names = TRUE)

# boxplot with geom_points of with and without groups
library(ggpubr)
pdf(here(plot_dir,"layer_marker_ranking_boxplot.pdf"), width=5, height=5)
p <- ggplot(df, aes(x=artifact, y=Top, fill=artifact)) +
  geom_boxplot() +
  geom_point(position=position_dodge(width=0.75), size=3) +
  geom_line(aes(group=rownames(df))) +
  scale_fill_manual(values=c("red", "grey")) +
  theme_bw() +
  theme(text=element_text(size=18),
        plot.title = element_text(size=24),
        #legend.position = "bottom",
        legend.text=element_text(size=14)) +
  ggtitle("L6 marker genes") +
  ylab("Ranking") +
  xlab("Hangnail") +
  stat_compare_means(method = "wilcox.test", paired = TRUE)
p
dev.off()


# sort rownames by ranking descending
df <- df[order(df$Top, decreasing=TRUE),]

# make a bar plot instead
pdf(here(plot_dir,"layer_marker_ranking_barplot.pdf"), width=5, height=5)
p <- ggplot(df, aes(x=Top, y=reorder(rownames(df),Top), fill=reorder(artifact,Top))) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_text(aes(label=Top), position=position_dodge(width=0.9), hjust=-.15) +
  scale_fill_manual(values=c("grey", "red"), "Hangnail") +
  theme_bw() +
  theme(text=element_text(size=18),
        plot.title = element_text(size=24),
        #legend.position = "bottom",
        legend.text=element_text(size=14)) +
  ggtitle("L6 marker genes") +
  ylab("Genes") +
  xlab("Ranking") +
  xlim(0,200)
p
dev.off()



# ==== reclustering dryspot sample after artifact removal ====

# ===== BS clustering of the dryspot samples ======
spe.no.dryspot <- spe.dryspot[,colData(spe.dryspot)$DBSCAN_artifact == FALSE]
spe.no.dryspot <- logNormCounts(spe.no.dryspot)

set.seed(101)
dec <- scran::modelGeneVar(spe.no.dryspot)
top <- scran::getTopHVGs(dec, n = 3000)

set.seed(102)
spe.no.dryspot <- scater::runPCA(spe.no.dryspot, subset_row=top)

## Add BayesSpace metadata
spe.no.dryspot <- spatialPreprocess(spe.no.dryspot, platform="Visium", skip.PCA=TRUE)

colData(spe.no.dryspot)$row <- spe.no.dryspot$array_row
colData(spe.no.dryspot)$col <- spe.no.dryspot$array_col

metadata(spe.no.dryspot)$BayesSpace.data <- list(platform = "Visium", is.enhanced = FALSE)

d <- 15  # Number of PCs

## Run BayesSpace clustering
set.seed(104)
dlpfc <- spatialCluster(spe.no.dryspot, q=2, d=d, platform='Visium',
                        nrep=10000, gamma=3, save.chain=TRUE)
spe.no.dryspot$BS_k2 <- as.factor(dlpfc$spatial.cluster)

dlpfc <- spatialCluster(spe.no.dryspot, q=7, d=d, platform='Visium',
                        nrep=10000, gamma=3, save.chain=TRUE)
spe.no.dryspot$BS_k7 <- as.factor(dlpfc$spatial.cluster)


pdf(here(plot_dir,"dryspot_noArtifact_clustering.pdf"), width=5, height=5)
p <- make_escheR(spe.no.dryspot) |>
  add_fill("BS_k7") +
  ggtitle("BayesSpace k=7") +
  theme(text=element_text(size=18),
        plot.title = element_text(size=24),
        legend.position = "bottom",
        legend.text=element_text(size=14)) +
  scale_fill_manual(values=colors)
p
dev.off()

# to csv
dryspot_noartifact_df <- data.frame(
  sample_id = colData(spe.no.dryspot)$sample_id,
  BS_k7 = as.character(colData(spe.no.dryspot)$BS_k7),
  row = colData(spe.no.dryspot)$array_row,
  col = colData(spe.no.dryspot)$array_col
)

# Export to CSV
write.csv(dryspot_noartifact_df, file=here("processed-data", "outputs_for_paper", "figure_6", "Figure6_D.csv"), row.names = FALSE)



# ========== Local outliers for schematic ==========

spe.hangnail <- localOutliers(spe.hangnail,
                     metric="sum_umi",
                     direction="lower"
)

spe.hangnail <- localOutliers(spe.hangnail,
                     metric="sum_gene",
                     direction="lower"
)

spe.hangnail <- localOutliers(spe.hangnail,
                     metric="expr_chrM_ratio",
                     direction="higher",
                     log=FALSE
)

spe.hangnail$local_outliers <- as.logical(spe.hangnail$sum_umi_outliers) |
  as.logical(spe.hangnail$sum_gene_outliers) |
  as.logical(spe.hangnail$expr_chrM_ratio_outliers)


pdf(here(plot_dir,"hangnail_local_outliers.pdf"), width=5.5, height=5.5)
SpotSweeper::plotQC(spe.hangnail, metric="sum_umi_log", outliers="local_outliers", point_size=2.1) +
  ggtitle("Local outliers") +
  theme(text=element_text(size=18),
        plot.title = element_text(size=24),
        legend.text=element_text(size=14))
dev.off()


# ========================================================
# Quantifying Qc metrics between artifact and good tissue
# ========================================================

spe.dryspot <- findArtifacts(spe.dryspot,
                              mito_percent="expr_chrM_ratio",
                              mito_sum="expr_chrM",
                              n_rings=5,
                              name="artifact"
)

# hangnail
df <- data.frame(spe.dryspot$k6,
                 spe.dryspot$k18,
                 spe.dryspot$k36,
                 spe.dryspot$k60,
                 spe.dryspot$k90) %>% rowMeans()

spe.dryspot$multiscale_var <- df

# make dataframe of the three QC metrics from DBSCAN_artifact
df.dryspot <- data.frame(umi=spe.dryspot$sum_umi,
                 gene=spe.dryspot$sum_gene,
                 mito=spe.dryspot$expr_chrM_ratio,
                 var=spe.dryspot$multiscale_var,
                 artifact=spe.dryspot$DBSCAN_artifact
                 )

# library size box plot
p1 <- ggplot(df.dryspot, aes(x=artifact, y=umi)) +
  geom_boxplot(aes(fill=artifact)) +
  theme_bw() +
  theme(text=element_text(size=18),
        plot.title = element_text(size=24),
        legend.position = "none",
        legend.text=element_text(size=14)) +
  ggtitle("Library size") +
  ylab("Total UMI") +
  xlab("Dryspot") +
  scale_fill_manual(values=c("grey", "red"))

p2 <- ggplot(df.dryspot, aes(x=artifact, y=gene)) +
  geom_boxplot(aes(fill=artifact)) +
  theme_bw() +
  theme(text=element_text(size=18),
        plot.title = element_text(size=24),
        legend.position = "none",
        legend.text=element_text(size=14)) +

  ggtitle("Unique Genes") +
  ylab("Number of  unique genes") +
  xlab("Dryspot") +
  scale_fill_manual(values=c("grey", "red"))

p3 <- ggplot(df.dryspot, aes(x=artifact, y=mito)) +
  geom_boxplot(aes(fill=artifact)) +
  theme_bw() +
  theme(text=element_text(size=18),
        plot.title = element_text(size=24),
        legend.position = "none",
        legend.text=element_text(size=14)) +
  ggtitle("Mito Ratio") +
  ylab("Mitochondrial Ratio") +
  xlab("Dryspot") +
  scale_fill_manual(values=c("grey", "red"))

p4 <- ggplot(df.dryspot, aes(x=artifact, y=var)) +
  geom_boxplot(aes(fill=artifact)) +
  theme_bw() +
  theme(text=element_text(size=18),
        plot.title = element_text(size=24),
        #legend.position = "bottom",
        legend.text=element_text(size=14)) +
  ggtitle("Multiscale variance") +
  ylab("log(Variance of Mito Ratio)") +
  xlab("Dryspot") +
  scale_fill_manual(values=c("grey", "red"))

pdf(here(plot_dir,"dryspot_qc_metrics.pdf"), width=20, height=6)
(p1+p2+p3+p4) +
  plot_layout(ncol = 4)
dev.off()

# hangnail
df <- data.frame(spe.hangnail$k6,
                 spe.hangnail$k18,
                 spe.hangnail$k36,
                 spe.hangnail$k60,
                 spe.hangnail$k90) %>% rowMeans()


df.hangnail <- data.frame(umi=spe.hangnail$sum_umi,
                 gene=spe.hangnail$sum_gene,
                 mito=spe.hangnail$expr_chrM_ratio,
                 var=spe.hangnail$multiscale_var,
                 artifact=spe.hangnail$artifact
                 )

# library size box plot
p1 <- ggplot(df.hangnail, aes(x=artifact, y=umi)) +
  geom_boxplot(aes(fill=artifact)) +
  theme_bw() +
  theme(text=element_text(size=18),
        plot.title = element_text(size=24),
        legend.position = "none",
        legend.text=element_text(size=14)) +
  ggtitle("Library size") +
  ylab("Total UMI") +
  xlab("Hangnail") +
  scale_fill_manual(values=c("grey", "red"))

p2 <- ggplot(df.hangnail, aes(x=artifact, y=gene)) +
  geom_boxplot(aes(fill=artifact)) +
  theme_bw() +
  theme(text=element_text(size=18),
        plot.title = element_text(size=24),
        legend.position = "none",
        legend.text=element_text(size=14)) +

  ggtitle("Unique Genes") +
  ylab("Number of  unique genes") +
  xlab("Hangnail") +
  scale_fill_manual(values=c("grey", "red"))

p3 <- ggplot(df.hangnail, aes(x=artifact, y=mito)) +
  geom_boxplot(aes(fill=artifact)) +
  theme_bw() +
  theme(text=element_text(size=18),
        plot.title = element_text(size=24),
        legend.position = "none",
        legend.text=element_text(size=14)) +
  ggtitle("Mito Ratio") +
  ylab("Mitochondrial Ratio") +
  xlab("Hangnail") +
  scale_fill_manual(values=c("grey", "red"))

p4 <- ggplot(df.hangnail, aes(x=artifact, y=var)) +
  geom_boxplot(aes(fill=artifact)) +
  theme_bw() +
  theme(text=element_text(size=18),
        plot.title = element_text(size=24),
        #legend.position = "bottom",
        legend.text=element_text(size=14)) +
  ggtitle("Multiscale variance") +
  ylab("log(Variance of Mito Ratio)") +
  xlab("Hangnail") +
  scale_fill_manual(values=c("grey", "red"))

pdf(here(plot_dir,"hangnail_qc_metrics.pdf"), width=20, height=6)
(p1+p2+p3+p4) +
  plot_layout(ncol = 4)
dev.off()


