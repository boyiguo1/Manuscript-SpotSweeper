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

plot_dir = here('plots',"Visium", "figure_4_5", "other_clustering")
processed_dir = here('raw-data')


# large DLPFC dataset
load(here("processed-data","Visium","dlPFC_raw.Rdata"))
spe.pfc <- spe_raw
spe.pfc



spe.hangnail <- spe.pfc[,colData(spe.pfc)$sample_id == "Br8325_ant"]

# drop out of tissue spots $in_tissue
spe.hangnail <- spe.hangnail[,spe.hangnail$in_tissue == TRUE]


library(ggthemes)
library(scales)
colors<- colorblind_pal()(7)

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

png(here(plot_dir,"hangnail_BayesSpace_k7.png"), width=5, height=5, units="in", res=300)
p1 <- make_escheR(spe.hangnail) |>
  add_fill("BS_k7") +
  ggtitle("BayesSpace") +
  theme(text=element_text(size=18),
        plot.title = element_text(size=24),
        legend.position = "bottom",
        legend.text=element_text(size=14)) +
  #change colors
  scale_fill_manual(values=colors)

p1
dev.off()



# ====== Banksy ======
# remove duplicate colnames
colnames(colData(spe.hangnail)) <- make.unique(colnames(colData(spe.hangnail)))

library(Banksy)

spe.hangnail.hvg <- spe.hangnail[top,]

lambda <- 0.2
k_geom <- 18
npcs <- 50
aname <- "logcounts"
spe.hangnail.hvg <- Banksy::computeBanksy(spe.hangnail.hvg, assay_name = aname, k_geom = k_geom)

set.seed(1000)
spe.hangnail.hvg <- Banksy::runBanksyPCA(spe.hangnail.hvg, lambda = lambda, npcs = npcs)

set.seed(1000)
spe.hangnail.hvg <- Banksy::clusterBanksy(spe.hangnail.hvg, lambda = lambda, npcs = npcs, resolution = 0.8)

# drop duplicate colData
colData(spe.hangnail.hvg) <- colData(spe.hangnail.hvg)[!duplicated(colnames(colData(spe.hangnail.hvg)))]

colors<- colorblind_pal()(length(unique(colData(spe.hangnail)$clust_M0_lam0.2_k50_res0.8)))
png(here(plot_dir, "hangnail_banksy.png"), width=5, height=5, units="in", res=300)
make_escheR(spe.hangnail.hvg) |>
    add_fill(var="clust_M0_lam0.2_k50_res0.8", point_size=2.5) +
    ggtitle("BANKSY") +
    theme(text=element_text(size=18),
        plot.title = element_text(size=24),
        legend.position = "bottom",
        legend.text=element_text(size=14)) +
        guides(fill = guide_legend(title = "BANKSY")) +
        scale_fill_manual(values=colors)
dev.off()



# ============ PRECAST ============
library(PRECAST)
library(Seurat)
library(purrr)

spe <- spe.hangnail

colnames(spe) <- spe$key

seuList <- list()
sample_ids <- unique(spe$sample_id)

for (id in sample_ids) {
    tmp_spe <- spe[, spe$sample_id == id]
    
    tmp_spe$row <- tmp_spe$array_row
    tmp_spe$col <- tmp_spe$array_col
    
    # Ensure unique row names
    rownames(tmp_spe) <- make.unique(rownames(tmp_spe))
    
    seuList[[id]] <- CreateSeuratObject(
        counts = counts(tmp_spe),
        meta.data = data.frame(colData(tmp_spe))
    )
}

set.seed(1)
preobj <- CreatePRECASTObject(seuList = seuList, gene.number=2000, selectGenesMethod='HVGs',
                              premin.spots = 1, premin.features=1, postmin.spots=1, postmin.features=1)
preobj@seulist

PRECASTObj <- AddAdjList(preobj, platform = "Visium")

PRECASTObj <- AddParSetting(PRECASTObj, Sigma_equal = FALSE,  maxIter = 30, verbose = TRUE)

PRECASTObj <- PRECAST(PRECASTObj, K = 7)
PRECASTObj

PRECASTObj <- SelectModel(PRECASTObj)
seuInt <- IntegrateSpaData(PRECASTObj, species = "Human")

# Merge with spe object
cluster_df <- seuInt@meta.data |>
    mutate(cluster = factor(cluster)) |>
    rename_with(~ paste0("PRECAST_", .x)) |>
    tibble::rownames_to_column(var = "key")

col_data_df <- colData(spe) |>
    data.frame() |>
    left_join(cluster_df, by="key")

rownames(col_data_df) <- colnames(spe)
colData(spe)$PRECAST_cluster <- col_data_df$PRECAST_cluster


# subset sample
spe.subset <- subset(spe, ,sample_id == unique(spe$sample_id)[1])

colors<- colorblind_pal()(7)
png(here(plot_dir, "hangnail_precast.png"), width = 5, height = 5, units = 'in', res = 300)
make_escheR(spe.subset) |>
    add_fill(var='PRECAST_cluster') +
    theme(text=element_text(size=18),
        plot.title = element_text(size=24),
        legend.position = "bottom",
        legend.text=element_text(size=14)) +
        guides(fill = guide_legend(title = "PRECAST")) +
        scale_fill_manual(values=colors) +
        ggtitle("PRECAST")
dev.off()

#saveRDS(spe, file = here(processed_dir, 'figure_2', "spe_spotsweeper_precast.rds"))


