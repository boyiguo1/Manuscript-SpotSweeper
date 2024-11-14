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



# === Updating localVariance function to include mean-var correct variables ===

localVariance <- function(spe, n_neighbors = 36,
                          metric = c("expr_chrM_ratio"),
                          samples = "sample_id", log = FALSE, name = NULL,
                          workers = 1, correct_mean_variance = TRUE) {

  # ===== start function =====
  # log1p transform specified metric
  if (log) {
    metric_log <- paste0(metric, "_log")
    colData(spe)[metric_log] <- log1p(colData(spe)[[metric]])
    metric_to_use <- metric_log
  } else {
    metric_to_use <- metric
  }

  # Get a list of unique sample IDs
  unique_sample_ids <- unique(colData(spe)[[samples]])

  # Initialize list to store each columnData dataframe
  columnData_list <- sapply(unique_sample_ids, FUN = function(x) NULL)

  # Loop through each unique sample ID
  for (sample_id in seq_along(unique_sample_ids)) {
    # Subset the data for the current sample
    sample <- unique_sample_ids[sample_id]
    spe_subset <- subset(spe, , sample_id == sample)

    # Create a list of spatial coordinates and qc metric
    columnData <- colData(spe_subset)
    columnData$coords <- spatialCoords(spe_subset)

    # Find nearest neighbors
    dnn <- BiocNeighbors::findKNN(spatialCoords(spe_subset),
                                  k = n_neighbors,
                                  BPPARAM = BiocParallel::MulticoreParam(workers),
                                  warn.ties = FALSE
    )$index

    #  === Compute local variance ===
    # get neighborhood metric
    neighborhoods <- lapply(seq_len(nrow(dnn)), function(i) {
      indices <- dnn[i, ]
      indices <- indices[indices != 0]
      indices <- c(i, indices)

      columnData[indices, metric_to_use]
    })

    # Compute variance and mean
    stats_matrix <- t(sapply(neighborhoods, function(x) {
      c(var = var(x, na.rm = TRUE), mean = mean(x, na.rm = TRUE))
    }))

    # Handle non-finite values
    stats_matrix[!is.finite(stats_matrix)] <- 0

    # Log2 transform the variance
    stats_matrix[, "var"] <- log2(stats_matrix[, "var"])

    if (correct_mean_variance) {
      # Perform robust linear regression to regress out mean-var bias
      fit.irls <- MASS::rlm(var ~ mean,
                            data = as.data.frame(stats_matrix))

      var_resid <- resid(fit.irls) # Get residuals

      # add local variance to columnData dataframe
      if (!is.null(name)) {
        columnData[name] <- var_resid
      } else {
        metric_var <- paste0(metric, "_var")
        columnData[metric_var] <- var_resid
      }
    } else {
      # add local variance to columnData dataframe without correction
      if (!is.null(name)) {
        columnData[name] <- stats_matrix[, "var"]
      } else {
        metric_var <- paste0(metric, "_var")
        columnData[metric_var] <- stats_matrix[, "var"]
      }
    }

    # Store the modified columnData dataframe in the list
    columnData_list[[sample_id]] <- columnData
  }

  # rbind the list of dataframes
  columnData_aggregated <- do.call(rbind, columnData_list)

  # replace SPE column data with aggregated data
  colData(spe) <- columnData_aggregated

  return(spe)
}


# ==== artifact detection with no mean-var option ====
findArtifacts <- function(
    spe, mito_percent = "expr_chrM_ratio",
    mito_sum = "expr_chrM", samples = "sample_id", n_order = 5,
    shape = "hexagonal", log = TRUE, name = "artifact", var_output = TRUE,
    correct_mean_variance = TRUE) {

  # Check if spe contains more than one sample_id
  unique_sample_ids <- unique(colData(spe)[[samples]])
  if (length(unique_sample_ids) > 1) {
    stop("The input SPE must contain only one sample. Please ensure only samples containing artifacts are passed.")
  }

  # Initialize a list to store spe for each sample
  columnData_list <- sapply(unique_sample_ids, FUN = function(x) NULL)

  for (sample in unique_sample_ids) {
    # subset by sample
    spe.temp <- spe[, colData(spe)[[samples]] == sample]

    # ======= Calculate local mito variance ========
    # Use vapply to iterate over order_seq
    var_matrix <- vapply(seq_len(n_order), function(i) {
      # Calculate n_neighbors for the current order
      if (shape == "hexagonal") {
        n_neighbors <- 3 * i * (i + 1)
      } else if (shape == "square") {
        n_neighbors <- 4 * i * (i + 1)
      }
      tmp.name <- paste0("k", n_neighbors)

      # Apply local variance calculation
      spe.temp <<- localVariance(spe.temp,
                                 metric = mito_percent,
                                 n_neighbors = n_neighbors,
                                 name = tmp.name,
                                 log = log,
                                 correct_mean_variance = correct_mean_variance)

      # Extract and return the column corresponding to tmp.name from colData
      colData(spe.temp)[[tmp.name]]
    }, numeric(length(spe.temp[[1]])))


    # ========== PCA and clustering ==========
    var_matrix <- cbind(
      var_matrix,
      colData(spe.temp)[[mito_percent]],
      colData(spe.temp)[[mito_sum]]
    )

    var_df <- data.frame(var_matrix)

    # Run PCA and add to reduced dims
    pc <- prcomp(var_df, center = TRUE, scale. = TRUE)
    rownames(pc$x) <- colnames(spe.temp) # assign rownames to avoid error
    reducedDim(spe.temp, "PCA_artifacts") <- pc$x

    # Cluster using kmeans and add to temp sce
    clus <- kmeans(pc$x, centers = 2, nstart = 25)
    spe.temp$Kmeans <- clus$cluster

    # =========== Artifact annotation ===========

    # calculate average local variance of the two clusters
    clus1_mean <- mean(colData(spe.temp)[[paste0("k", if (shape == "hexagonal") 18 else 24)]][spe.temp$Kmeans == 1])
    clus2_mean <- mean(colData(spe.temp)[[paste0("k", if (shape == "hexagonal") 18 else 24)]][spe.temp$Kmeans == 2])

    artifact_clus <- which.min(c(clus1_mean, clus2_mean))

    # create a new $artifact column; if clus1 mean < clus2 - rename Kmeans
    # 1 to 'TRUE'
    spe.temp$artifact <- FALSE
    spe.temp$artifact[spe.temp$Kmeans == artifact_clus] <- TRUE
    spe.temp$Kmeans <- NULL

    # add to sample list
    columnData_list[[sample]] <- colData(spe.temp)
  }

  # rbind the list of dataframes
  columnData_aggregated <- do.call(rbind, columnData_list)

  # replace SPE column data with aggregated data
  colData(spe) <- columnData_aggregated
  reducedDims(spe) <- reducedDims(spe.temp)

  # return sce
  return(spe)
}




# ======= Subset to hangnail sample ========

spe.hangnail <- spe.pfc[,colData(spe.pfc)$sample_id == "Br8325_ant"]

# drop out of tissue spots
spe.hangnail <- spe.hangnail[,colData(spe.hangnail)$in_tissue == TRUE]


# ======= Calculating local variance ========

# with correction
spe.hangnail <- localVariance(spe.hangnail, n_neighbors = 36,
                              metric = c("expr_chrM_ratio"),
                              samples = "sample_id", log = TRUE,
                              name = "expr_chrM_ratio_var_corrected",
                              workers = 1, correct_mean_variance = TRUE)

# without correction
spe.hangnail <- localVariance(spe.hangnail,  n_neighbors = 36,
                              metric = c("expr_chrM_ratio"),
                              samples = "sample_id", log = TRUE,
                              name = "expr_chrM_ratio_var_uncorrected",
                              workers = 1, correct_mean_variance = FALSE)

colnames(colData(spe.hangnail))

# dataframe
df <- data.frame(expr_chrM_ratio_var_corrected = spe.hangnail$expr_chrM_ratio_var_corrected,
                 expr_chrM_ratio_var_uncorrected = spe.hangnail$expr_chrM_ratio_var_uncorrected,
                 expr_chrM_ratio_log = spe.hangnail$expr_chrM_ratio)

png(here(plot_dir, "mean_variance_relationship_corrected.png"), width = 5, height = 5, units = "in", res = 300)
ggplot(df, aes(x = expr_chrM_ratio_log, y = expr_chrM_ratio_var_corrected)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(title = "Corrected Relationship",
       x = "log2(expr_chrM_ratio)",
       y = "IRLS residuals") +
  theme(text=element_text(size=18),
        plot.title = element_text(size=22),
        legend.text=element_text(size=14))  +
ggtitle("Corrected mean-var")
dev.off()

png(here(plot_dir, "mean_variance_relationship_uncorrected.png"), width = 5, height = 5, units = "in", res = 300)
ggplot(df, aes(x = expr_chrM_ratio_log, y = expr_chrM_ratio_var_uncorrected)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(title = "Uncorrected Relationship",
       x = "log2(expr_chrM_ratio)",
       y = "log2(Local Variance)") +
  theme(text=element_text(size=18),
        plot.title = element_text(size=22),
        legend.text=element_text(size=14)) +
        ggtitle("Uncorrected mean-var")
dev.off()


# ====== Artifact detection with and without correction ======

# with correction
spe.hangnail <- findArtifacts(spe.hangnail, mito_percent = "expr_chrM_ratio",
                              mito_sum = "expr_chrM", samples = "sample_id", n_order = 5,
                              shape = "hexagonal", log = TRUE,
                              var_output = TRUE, correct_mean_variance = TRUE)

spe.hangnail$artifact_corrected <- spe.hangnail$artifact
reducedDim(spe.hangnail, "PCA_artifact_corrected") <- reducedDim(spe.hangnail, "PCA_artifacts")

# without correction
spe.hangnail <- findArtifacts(spe.hangnail, mito_percent = "expr_chrM_ratio",
                              mito_sum = "expr_chrM", samples = "sample_id", n_order = 5,
                              shape = "hexagonal", log = TRUE, name = "artifact",
                              var_output = TRUE, correct_mean_variance = FALSE)
spe.hangnail$artifact_uncorrected <- spe.hangnail$artifact    
reducedDim(spe.hangnail, "PCA_artifact_uncorrected") <- reducedDim(spe.hangnail, "PCA_artifacts")    

colnames(colData(spe.hangnail))

# ===== plot clusters =====
png(here(plot_dir, "artifact_detection_corrected.png"), width = 4, height = 4, units = "in", res = 300)
make_escheR(spe.hangnail) |>
    add_fill(var="artifact_corrected", point_size=1.75) +
    scale_fill_manual(values=c("grey", "red"),
    name="Artifact") +
    theme(legend.position = "bottom") +
    ggtitle("SpotPlot: Corrected")
dev.off()

png(here(plot_dir, "artifact_detection_uncorrected.png"), width = 4, height = 4, units = "in", res = 300)
make_escheR(spe.hangnail) |>
    add_fill(var="artifact_uncorrected", point_size=1.75) +
    scale_fill_manual(values=c("grey", "red"),
    name="Artifact") +
    theme(legend.position = "bottom") +
    ggtitle("SpotPlot: Uncorrected")
dev.off()

# ===== plot variance =====
png(here(plot_dir, "artifact_variance_corrected.png"), width = 4.75, height = 4.75, units = "in", res = 300)
make_escheR(spe.hangnail) |>
    add_fill(var="k36") +
    scale_fill_gradient(low = "grey", high = "black")
dev.off()

png(here(plot_dir, "artifact_variance_uncorrected.png"), width = 4.75, height = 4.75, units = "in", res = 300)
make_escheR(spe.hangnail) |>
    add_fill(var="k36") +
    scale_fill_gradient(low = "grey", high = "black")
dev.off()

# ===== PCA =====
pca_corrected <- reducedDim(spe.hangnail,"PCA_artifact_corrected")[,1:2]
pca_uncorrected <- reducedDim(spe.hangnail,"PCA_artifact_uncorrected")[,1:2]

png(here(plot_dir,"artifact_detection_corrected_pca.png"), width=6, height=5, units="in", res=300)
ggplot(pca_corrected, aes(x=PC1, y=PC2)) +
  geom_point(aes(color=as.factor(spe.hangnail$artifact_corrected)), alpha=.5) +
  scale_color_manual(values=c("grey", "red"),
                     name="Artifact") +
  theme_bw() +
  theme(text=element_text(size=18),
        plot.title = element_text(size=22),
        #legend.position = "bottom",
        legend.text=element_text(size=14)) +
  ggtitle("K-means clustering: Corrected") +
  ylab("PC2") +
  xlab("PC1")
dev.off()

png(here(plot_dir,"artifact_detection_uncorrected_pca.png"), width=6, height=5, units="in", res=300)
ggplot(pca_uncorrected, aes(x=PC1, y=PC2)) +
  geom_point(aes(color=as.factor(spe.hangnail$artifact_uncorrected)), alpha=.5) +
  scale_color_manual(values=c("grey", "red"),
                     name="Artifact") +
  theme_bw() +
  theme(text=element_text(size=18),
        plot.title = element_text(size=22),
        #legend.position = "bottom",
        legend.text=element_text(size=14)) +
  ggtitle("K-means clustering: Uncorrected") +
  ylab("PC2") +
  xlab("PC1")
dev.off()

