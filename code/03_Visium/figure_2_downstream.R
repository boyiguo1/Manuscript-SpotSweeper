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

plot_dir = here('plots',"figure_2")
processed_dir = here('processed-data')

## Download the spot-level data
spe <- readRDS(here(processed_dir,"figure_2", 'spe_spotsweeper.rds'))

# get only top 1000 HVGs
dec <- modelGeneVar(spe)
chosen <- getTopHVGs(dec, n=1000)

spe.hvg <- spe[chosen,]

# ======= Banksy clustering =======
library(Banksy)

lambda <- 0.2
k_geom <- 18
npcs <- 50
aname <- "counts"
spe.hvg <- Banksy::computeBanksy(spe.hvg, assay_name = aname, k_geom = k_geom)

set.seed(1000)
spe.hvg <- Banksy::runBanksyPCA(spe.hvg, lambda = lambda, npcs = npcs)

set.seed(1000)
spe.hvg <- Banksy::clusterBanksy(spe.hvg, lambda = lambda, npcs = npcs, resolution = 0.8)


# ======== BayesSpace clustering ========
library(BayesSpace)

## Add BayesSpace metadata
spe.hvg <- spatialPreprocess(spe.hvg, platform="Visium", skip.PCA=TRUE)

colData(spe.hvg)$row <- spe.hvg$array_row
colData(spe.hvg)$col <- spe.hvg$array_col

metadata(spe.hvg)$BayesSpace.data <- list(platform = "Visium", is.enhanced = FALSE)

d <- 15  # Number of PCs

## Run BayesSpace clustering
dlpfc <- spatialCluster(spe.hvg, q=7, d=d, platform='Visium',
                        nrep=10000, gamma=3, save.chain=TRUE)
spe.hvg$BS_k7 <- as.factor(dlpfc$spatial.cluster)

# save
saveRDS(spe.hvg, here(processed_dir, 'figure_2', 'spe_spotsweeper_clusters.rds'))

# load
spe.hvg <- readRDS(here(processed_dir, 'figure_2', 'spe_spotsweeper_clusters.rds'))

spe$neighbors <- spe.hvg$neighbors
spe$BS_k7 <- spe.hvg$BS_k7
spe$Banksy <- spe.hvg$clust_M0_lam0.2_k50_res0.8

# ============ PRECAST ============
library(PRECAST)
library(CreateSeuratObject)
library(purrr)

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
# An object of class PRECASTObj 
#  with 12 datasets and  47681 spots in total, with spots for each dataset:  4226 4384 4789 4634 3661 3498 4110 4015 3639 3673 3592 3460 
#  2000 common variable genes selected



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

png(here(plot_dir, "figure_2_downstream_3.png"), width = 5, height = 5, units = 'in', res = 300)
make_escheR(spe.subset) |>
    add_fill(var='PRECAST_cluster')
dev.off()

saveRDS(spe, file = here(processed_dir, 'figure_2', "spe_spotsweeper_precast.rds"))

# load




# ========= Get how much neighbor share the same cluster =========
# start with bayesspace
spe$BS_k7 <- factor(spe$BS_k7, levels = c("1", "2", "3", "4", "5", "6", "7"))

# Function to calculate the percentage of neighbors in the same cluster
calculate_percent_neighbors_same_cluster <- function(spe, cluster_col) {
    # Initialize a list to store the results for each sample
    all_percent_neighbors_same_cluster <- list()
    
    # Loop through each sample
    samples <- unique(colData(spe)$sample_id)
    for (sample in samples) {
        # Subset to the current sample
        spe.subset <- spe[, colData(spe)$sample_id == sample]
        
        # Get first order neighbors of all spots 
        dnn1 <- BiocNeighbors::findKNN(spatialCoords(spe.subset), k = 6, warn.ties = FALSE)$index
        
        # Calculate the percentage of neighbors in the same cluster for each spot
        percent_neighbors_same_cluster <- sapply(seq_len(nrow(dnn1)), function(i) {
            neighbors <- dnn1[i, ]
            same_cluster <- sum(spe.subset[[cluster_col]][neighbors] == spe.subset[[cluster_col]][i])
            percent_same_cluster <- (same_cluster / length(neighbors)) * 100
            return(percent_same_cluster)
        })
        
        # Store the results in the list
        all_percent_neighbors_same_cluster[[sample]] <- percent_neighbors_same_cluster
    }
    
    # Combine the results into a single vector
    combined_percent_neighbors_same_cluster <- unlist(all_percent_neighbors_same_cluster)
    
    # Add the calculated percentages to the spe object
    spe[[paste0(cluster_col, "_k6_same_neighbors")]] <- combined_percent_neighbors_same_cluster
    return(spe)
}

# Example usage
spe <- calculate_percent_neighbors_same_cluster(spe, "BS_k7")
spe <- calculate_percent_neighbors_same_cluster(spe, "Banksy")
spe <- calculate_percent_neighbors_same_cluster(spe, "GraphBased")
spe <- calculate_percent_neighbors_same_cluster(spe, "SNN_k50_k7")

# Combine the results into a single dataframe
spe_df <- as.data.frame(colData(spe)) %>%
    pivot_longer(cols = ends_with("_k6_same_neighbors"), names_to = "method", values_to = "percent_neighbors") %>%
    mutate(method = gsub("_k6_same_neighbors", "", method))

# Calculate the average percent neighbors for local_outliers per sample and method
average_percent_neighbors <- spe_df %>%
    mutate(method = recode(method, BS_k7 = "BayesSpace", SNN_k50_k7 = "PRECAST")) %>%
    group_by(sample_id, local_outliers, method) %>%
    summarize(avg_percent_neighbors = mean(percent_neighbors, na.rm = TRUE))

# Reorder the methods
average_percent_neighbors$method <- factor(average_percent_neighbors$method, levels = c("GraphBased", "PRECAST", "BayesSpace", "Banksy"))

# Generate boxplots using ggplot2
library(ggplot2)

png(here(plot_dir, "figure_2_downstream_2.png"), width = 7, height = 3, units = 'in', res = 300)
ggplot(average_percent_neighbors, aes(x = local_outliers, y = avg_percent_neighbors, fill = local_outliers)) +
    geom_boxplot() +
    geom_jitter() +
    scale_fill_manual(values = c("grey", "red"), name = "Local Outliers") +
    labs(y = "Percent Neighbors in Same Cluster", x=NULL) +
    scale_y_continuous(limits = c(0, 100)) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    facet_wrap(~ method, nrow = 1) +
    theme(strip.text = element_text(size = 14), text = element_text(size = 12)) 
dev.off()


# save results to csv
write.csv(average_percent_neighbors, here(processed_dir, 'outputs_for_paper', 'Figure_2_neighbors_boxplot.csv'), row.names = FALSE)

