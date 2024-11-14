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

# get snRNA-seq data
sce_path_zip <- fetch_data("spatialDLPFC_snRNAseq")
sce_path <- unzip(sce_path_zip, exdir = tempdir())
sce <- HDF5Array::loadHDF5SummarizedExperiment(
    file.path(tempdir(), "sce_DLPFC_annotated")
)
sce
# class: SingleCellExperiment 
# dim: 36601 77604 
# metadata(3): Samples cell_type_colors cell_type_colors_broad
# assays(2): counts logcounts
# rownames(36601): MIR1302-2HG FAM138A ... AC007325.4 AC007325.2
# rowData names(7): source type ... gene_type binomial_deviance
# colnames(77604): 1_AAACCCAAGTTCTCTT-1 1_AAACCCACAAGGTCTT-1 ...
#   19_TTTGTTGTCTCATTGT-1 19_TTTGTTGTCTTAAGGC-1
# colData names(32): Sample Barcode ... cellType_layer layer_annotation
# reducedDimNames(4): GLMPCA_approx TSNE UMAP HARMONY
# mainExpName: NULL
# altExpNames(0):

colnames(colData(sce))
#  [1] "Sample"                "Barcode"               "key"                  
#  [4] "SAMPLE_ID"             "pos"                   "BrNum"                
#  [7] "round"                 "Position"              "age"                  
# [10] "sex"                   "diagnosis"             "sum"                  
# [13] "detected"              "subsets_Mito_sum"      "subsets_Mito_detected"
# [16] "subsets_Mito_percent"  "total"                 "high_mito"            
# [19] "low_sum"               "low_detected"          "discard_auto"         
# [22] "doubletScore"          "prelimCluster"         "collapsedCluster"     
# [25] "kmeans"                "sizeFactor"            "cellType_broad_k"     
# [28] "cellType_k"            "cellType_broad_hc"     "cellType_hc"          
# [31] "cellType_layer"        "layer_annotation" 

unique(sce$cellType_broad_hc)
# [1] Inhib     Oligo     OPC       Excit     Astro     EndoMural Ambiguous
# [8] Micro    
# Levels: Astro EndoMural Micro Oligo OPC Excit Inhib Ambiguous

# drop Ambiguous and refactor
sce <- sce[,sce$cellType_broad_hc != "Ambiguous"]
sce$cellType_broad_hc <- factor(sce$cellType_broad_hc)

# ========== Deconvolution ===========
library(SPOTlight)

# get top 3000 HVGs
dec <- modelGeneVar(sce)
hvg <- getTopHVGs(dec, n = 3000)


colLabels(sce) <- colData(sce)$cellType_broad_hc

# Get vector indicating which genes are neither ribosomal or mitochondrial
genes <- !grepl(pattern = "^Rp[l|s]|MT", x = rownames(sce))

# ======== Get marker genes for seeding NMF =============
mgs <- findMarkers(sce, groups=sce$cellType_broad_hc)

mgs_fil <- lapply(names(mgs), function(i) {
    x <- mgs[[i]]
    # Filter and keep relevant marker genes, those with AUC > 0.8
    x <- x[x$summary.logFC > 2, ]
    # Sort the genes from highest to lowest weight
    x <- x[order(x$summary.logFC, decreasing = TRUE), ]
    # Add gene and cluster id to the dataframe
    x$gene <- rownames(x)
    x$cluster <- i
    # Ensure all data frames have the same column names
    x <- x[, c("summary.logFC", "gene", "cluster")]
    data.frame(x)
})

mgs_df <- do.call(rbind, mgs_fil)



# ========= Cell downsampling =======
# split cell indices by identity
idx <- split(seq(ncol(sce)), sce$cellType_broad_hc)
# downsample to at most 20 per identity & subset
n_cells <- 20
cs_keep <- lapply(idx, function(i) {
    n <- length(i)
    if (n < n_cells)
        n_cells <- n
    sample(i, n_cells)
})
sce.dsp <- sce[, unlist(cs_keep)]


# ========= Run SPOTlight ===========


res <- SPOTlight(
    x = sce.dsp,
    y = spe,
    groups = sce.dsp$cellType_broad_hc,
    mgs = mgs_df,
    hvg = hvg,
    weight_id = "summary.logFC",
    group_id = "cluster",
    gene_id = "gene")

# Scaling count matrix
# Seeding initial matrices
# Training NMF model
# Time for training: 1.97min

# Extract deconvolution matrix
head(mat <- res$mat)[, seq_len(3)]
#                         Astro EndoMural Excit_ambig
# AAACAACGAATAGTTC-1 0.03829211 0.1761871  0.05229192
# AAACAAGTATCTCCCA-1 0.00000000 0.0000000  0.00000000
# AAACAATCTACTAGCA-1 0.00000000 0.0000000  0.00000000
# AAACACCAATAACTGC-1 0.00000000 0.0000000  0.00000000
# AAACAGCTTTCAGAAG-1 0.00000000 0.0000000  0.00000000
# AAACAGGGTCTATATT-1 0.00000000 0.0000000  0.00000000

# Extract NMF model fit
mod <- res$NMF

png(here(plot_dir, "SPOTlight_topic_profiles.png"), width=4, height=4, units="in", res=300)
plotTopicProfiles(
    x = mod,
    y = sce.dsp$cellType_broad_hc,
    facet = FALSE,
    min_prop = 0.01,
    ncol = 1) +
    theme(aspect.ratio = 1)
dev.off()

png(here(plot_dir, "SPOTlight_topic_profile_celltypes.png"), width=10, height=10, units="in", res=300)
plotTopicProfiles(
    x = mod,
    y = sce.dsp$cellType_broad_hc,
    facet = TRUE,
    min_prop = 0.01,
    ncol = 6)
dev.off()





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





# ========= Scatter pie chart ===========
ct <- colnames(mat)
#mat[mat < 0.1] <- 0

# palette for length spe$clust_M0_lam0.8_k50_res0.8
pal <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(length(unique(sce$cellType_broad_hc)))
names(pal) <- ct

#subset to first sample of spe
spe.subset <- subset(spe, ,sample_id == unique(spe$sample_id)[1])

# subset mat to only the spots in spe.subset
mat.subset <- mat[colnames(spe.subset),]

png(here(plot_dir, "SPOTlight_scatter_pie.png"), width=6, height=6, units="in", res=300)
plotSpatialScatterpie(
    x = spe.subset,
    y = mat.subset,
    cell_types = colnames(mat),
    img = FALSE,
    scatterpie_alpha = 1,
    pie_scale = 0.4) +
    scale_fill_manual(
        values = pal,
        breaks = names(pal))
dev.off()






# Get the average proportions of cell types in local outlier spots vs neighboring spots
# Get first order neighbors of all spots
dnn1 <- BiocNeighbors::findKNN(spatialCoords(spe.subset), k = 6, warn.ties = FALSE)$index

# Calculate the average proportions of cell types in local outlier spots vs neighboring spots
avg_proportion <- matrix(0, nrow = nrow(dnn1), ncol = ncol(mat.subset))
avg_proportion <- t(sapply(seq_len(nrow(dnn1)), function(i) {
    neighbors <- dnn1[i, ]
    colMeans(mat.subset[neighbors,])
}))


library(BiocNeighbors)

calculate_avg_celltype_proportion_neighbors <- function(spe, mat, sample_id_col = "sample_id", k = 6) {
    # Initialize a list to store results for each sample
    all_avg_proportions <- list()
    
    # Loop through each sample
    samples <- unique(colData(spe)[[sample_id_col]])
    for (sample in samples) {
        # Subset for the current sample
        spe.subset <- spe[, colData(spe)[[sample_id_col]] == sample]
        mat.subset <- mat[colData(spe)[[sample_id_col]] == sample, ]
        
        # Get first-order neighbors of each spot
        dnn1 <- BiocNeighbors::findKNN(spatialCoords(spe.subset), k = k, warn.ties = FALSE)$index
        
        # Calculate the average cell type proportions in neighbors for each spot
        avg_proportion_neighbors <- t(sapply(seq_len(nrow(dnn1)), function(i) {
            neighbors <- dnn1[i, ]
            colMeans(mat.subset[neighbors, , drop = FALSE])
        }))
        
        # Add the results for the sample to the list
        all_avg_proportions[[sample]] <- avg_proportion_neighbors
    }
    
    # Combine the results into a single matrix
    combined_avg_proportions <- do.call(rbind, all_avg_proportions)
    
    return(combined_avg_proportions)
}

# Calculate the average cell type proportions in neighbors for each spot
avg_proportion_neighbors <- calculate_avg_celltype_proportion_neighbors(spe.subset, mat.subset, sample_id_col = "sample_id", k = 6)


library(BiocNeighbors)

calculate_avg_proportion_outliers_vs_neighbors <- function(spe, mat, outliers, k = 6) {
    # Find neighbors
    dnn1 <- BiocNeighbors::findKNN(spatialCoords(spe), k = k, warn.ties = FALSE)$index
    
    # Identify outlier spots and neighbor spots that are also outliers
    outlier_indices <- which(outliers)
    outlier_neighbors_indices <- unique(unlist(dnn1[outlier_indices, ]))
    
    # Filter to neighbors that are also outliers
    outlier_neighbors_indices <- outlier_neighbors_indices[outliers[outlier_neighbors_indices]]
    
    # Calculate average cell type proportions for outliers and their neighboring outliers
    avg_proportion_outliers <- colMeans(mat[outlier_indices, , drop = FALSE])
    avg_proportion_outlier_neighbors <- colMeans(mat[outlier_neighbors_indices, , drop = FALSE])
    
    # Combine into a data frame for plotting
    proportions_df <- data.frame(
        CellType = rep(colnames(mat), 2),
        Proportion = c(avg_proportion_outliers, avg_proportion_outlier_neighbors),
        Group = rep(c("Outliers", "Outlier Neighbors"), each = ncol(mat))
    )
    
    return(proportions_df)
}


library(ggplot2)

plot_proportions <- function(proportions_df) {
    ggplot(proportions_df, aes(x = Group, y = Proportion, fill = CellType)) +
        geom_bar(stat = "identity") +
        labs(title = "Average Cell Type Proportions in Outliers vs Outlier Neighbors",
             x = "Group", y = "Average Proportion") +
        theme_minimal() +
        scale_fill_brewer(palette = "Set3") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Assuming `proportions_df` was created by the above function
proportions_df <- calculate_avg_proportion_outliers_vs_neighbors(spe, mat, spe$local_outliers, k = 6)

png(here(plot_dir, "SPOTlight_avg_proportions_outliers_vs_neighbors.png"), width=6, height=4, units="in", res=300)
plot_proportions(proportions_df)
dev.off()





calculate_avg_proportion_outliers_vs_neighbors <- function(spe, mat, outliers, k = 6) {
    # Find neighbors
    dnn1 <- BiocNeighbors::findKNN(spatialCoords(spe), k = k, warn.ties = FALSE)$index
    
    # Identify outlier spots and neighbor spots that are also outliers
    outlier_indices <- which(outliers)
    outlier_neighbors_indices <- unique(unlist(dnn1[outlier_indices, ]))
    
    # Filter to neighbors that are also outliers
    outlier_neighbors_indices <- outlier_neighbors_indices[outliers[outlier_neighbors_indices]]
    
    # Calculate average cell type proportions for outliers and their neighboring outliers
    avg_proportion_outliers <- colMeans(mat[outlier_indices, , drop = FALSE])
    avg_proportion_outlier_neighbors <- colMeans(mat[outlier_neighbors_indices, , drop = FALSE])
    
    # Normalize to ensure they sum to 1
    avg_proportion_outliers <- avg_proportion_outliers / sum(avg_proportion_outliers)
    avg_proportion_outlier_neighbors <- avg_proportion_outlier_neighbors / sum(avg_proportion_outlier_neighbors)
    
    # Combine into a data frame for plotting
    proportions_df <- data.frame(
        CellType = rep(colnames(mat), 2),
        Proportion = c(avg_proportion_outliers, avg_proportion_outlier_neighbors),
        Group = rep(c("Outliers", "Outlier Neighbors"), each = ncol(mat))
    )
    
    return(proportions_df)
}

# Plotting the normalized proportions
proportions_df <- calculate_avg_proportion_outliers_vs_neighbors(spe, mat, spe$local_outliers, k = 6)

png(here(plot_dir, "SPOTlight_avg_proportions_outliers_vs_neighbors.png"), width=3, height=4, units="in", res=300)
plot_proportions(proportions_df)
dev.off()







calculate_avg_proportion_per_sample <- function(spe, mat, outliers, sample_id_col = "sample_id", k = 6) {
    # Initialize a list to store the results
    all_sample_proportions <- list()
    
    # Loop through each sample
    samples <- unique(colData(spe)[[sample_id_col]])
    for (sample in samples) {
        # Subset to the current sample
        spe.subset <- spe[, colData(spe)[[sample_id_col]] == sample]
        mat.subset <- mat[colData(spe)[[sample_id_col]] == sample, ]
        outliers.subset <- outliers[colData(spe)[[sample_id_col]] == sample]
        
        # Get first-order neighbors for the current sample
        dnn1 <- BiocNeighbors::findKNN(spatialCoords(spe.subset), k = k, warn.ties = FALSE)$index
        
        # Identify outlier spots and outlier neighbors within this sample
        outlier_indices <- which(outliers.subset)
        outlier_neighbors_indices <- unique(unlist(dnn1[outlier_indices, ]))
        
        # Filter to neighbors that are also outliers
        outlier_neighbors_indices <- outlier_neighbors_indices[outliers.subset[outlier_neighbors_indices]]
        
        # Calculate average proportions for outliers and outlier neighbors
        avg_proportion_outliers <- colMeans(mat.subset[outlier_indices, , drop = FALSE])
        avg_proportion_outlier_neighbors <- colMeans(mat.subset[outlier_neighbors_indices, , drop = FALSE])
        
        # Normalize proportions
        avg_proportion_outliers <- avg_proportion_outliers / sum(avg_proportion_outliers)
        avg_proportion_outlier_neighbors <- avg_proportion_outlier_neighbors / sum(avg_proportion_outlier_neighbors)
        
        # Store results in a data frame with sample information
        sample_df <- data.frame(
            CellType = rep(colnames(mat), 2),
            Proportion = c(avg_proportion_outliers, avg_proportion_outlier_neighbors),
            Group = rep(c("Outliers", "Outlier Neighbors"), each = ncol(mat)),
            Sample = sample
        )
        
        # Append to the list
        all_sample_proportions[[sample]] <- sample_df
    }
    
    # Combine all samples into one data frame
    combined_proportions_df <- do.call(rbind, all_sample_proportions)
    
    return(combined_proportions_df)
}

# Generate the data frame with sample-level proportions
proportions_df <- calculate_avg_proportion_per_sample(spe, mat, outliers, sample_id_col = "sample_id", k = 6)


library(ggplot2)

plot_boxplots_per_celltype <- function(proportions_df) {
    ggplot(proportions_df, aes(x = Group, y = Proportion, fill = Group)) +
        geom_boxplot() +
        facet_wrap(~ CellType, scales = "free_y") +
        labs(title = "Average Cell Type Proportions in Outliers vs Outlier Neighbors Across Samples",
             x = "Group", y = "Average Proportion") +
        theme_minimal() +
        scale_fill_brewer(palette = "Set3") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Plotting the box plots
png(here(plot_dir, "SPOTlight_avg_proportions_outliers_vs_neighbors_boxplots.png"), width=10, height=10, units="in", res=300)
plot_boxplots_per_celltype(proportions_df)

dev.off()












calculate_avg_proportion_per_sample <- function(spe, mat, sample_id_col = "sample_id", k = 6) {
    # Initialize a list to store the results
    all_sample_proportions <- list()
    
    # Loop through each sample
    samples <- unique(colData(spe)[[sample_id_col]])
    for (sample in samples) {
        # Subset to the current sample
        spe.subset <- spe[, colData(spe)[[sample_id_col]] == sample]
        mat.subset <- mat[colData(spe)[[sample_id_col]] == sample, ]
        
        # Get first-order neighbors for the current sample
        dnn1 <- BiocNeighbors::findKNN(spatialCoords(spe.subset), k = k, warn.ties = FALSE)$index
        
        # Identify outlier spots and outlier neighbors within this sample
        outlier_indices <- spe.subset$local_outliers
        outlier_neighbors_indices <- unique(unlist(dnn1[outlier_indices, ]))
        
        # Calculate average proportions for outliers
        avg_proportion_outliers <- colMeans(mat.subset[outlier_indices, , drop = FALSE])
        avg_proportion_outliers <- avg_proportion_outliers / sum(avg_proportion_outliers)  # Normalize
        
        # Calculate average proportions for outlier neighbors if they exist
        if (length(outlier_neighbors_indices) > 0) {
            avg_proportion_outlier_neighbors <- colMeans(mat.subset[outlier_neighbors_indices, , drop = FALSE])
            avg_proportion_outlier_neighbors <- avg_proportion_outlier_neighbors / sum(avg_proportion_outlier_neighbors)  # Normalize
        } else {
            avg_proportion_outlier_neighbors <- rep(NaN, ncol(mat.subset))
        }
        
        # Store results in a data frame with sample information
        sample_df <- data.frame(
            CellType = rep(colnames(mat), 2),
            Proportion = c(avg_proportion_outliers, avg_proportion_outlier_neighbors),
            Group = rep(c("Outliers", "Outlier Neighbors"), each = ncol(mat)),
            Sample = sample
        )
        
        # Append to the list
        all_sample_proportions[[sample]] <- sample_df
    }
    
    # Combine all samples into one data frame
    combined_proportions_df <- do.call(rbind, all_sample_proportions)
    
    return(combined_proportions_df)
}

# Generate the data frame with sample-level proportions
proportions_df <- calculate_avg_proportion_per_sample(spe, mat, sample_id_col = "sample_id", k = 6)


plot_boxplots_per_celltype <- function(proportions_df) {
    ggplot(proportions_df, aes(x = Group, y = Proportion, fill = Group)) +
        geom_boxplot() +
        facet_wrap(~ CellType, scales = "free_y") +
        labs(title = "Average Cell Type Proportions in Outliers vs Outlier Neighbors Across Samples",
             x = "Group", y = "Average Proportion") +
        theme_minimal() +
        scale_fill_brewer(palette = "Set3") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Plotting the box plots
png(here(plot_dir, "SPOTlight_avg_proportions_outliers_vs_neighbors_boxplots.png"), width=10, height=10, units="in", res=300)
plot_boxplots_per_celltype(proportions_df)
dev.off()