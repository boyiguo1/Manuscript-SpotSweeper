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
library(ggspavis)


plot_dir <- here("plots","parameter_selection")

spe <- readRDS(here("processed-data", "VisiumHD","human_breast", "VisiumHD_HumanBreast_016_bayesspace.rds"))
spe
# dim: 18085 119082 
# metadata(1): BayesSpace.data
# assays(3): '' counts logcounts
# rownames(18085): SAMD11 NOC2L ... MT-ND6 MT-CYB
# rowData names(2): symbol is.HVG
# colnames(119082): s_016um_00107_00066-1 s_016um_00126_00213-1 ...
#   s_016um_00361_00190-1 s_016um_00242_00130-1
# colData names(16): barcode in_tissue ... col BS_k10
# reducedDimNames(1): PCA
# mainExpName: NULL
# altExpNames(0):
# spatialCoords names(2) : xcoord ycoord
# imgData names(1): sample_id


# Find 1st-order neighbors (k = 8)
dnn1 <- BiocNeighbors::findKNN(spatialCoords(spe), k = 8, warn.ties = FALSE)$index
spe$neighbors <- "Non-neighbor"
spe$neighbors[1] <- "Center"  # Center
spe$neighbors[dnn1[1,]] <- "1st Order (k=8)"  # First-order neighbors

# Find 2nd-order neighbors (k = 24) based on first-order neighbors
dnn2 <- BiocNeighbors::findKNN(spatialCoords(spe), k = 24, warn.ties = FALSE)$index

second_order_indices <- dnn2[1,!(dnn2[1,] %in% dnn1[1,]) ]
spe$neighbors[second_order_indices] <- "2nd Order (k=24)"  # Second-order neighbors

# Find 3rd-order neighbors (k = 48) based on second-order neighbors
dnn3 <- BiocNeighbors::findKNN(spatialCoords(spe), k = 48, warn.ties = FALSE)$index

third_order_indices <- dnn3[1, !(dnn3[1,] %in% dnn2[1,])]
spe$neighbors[third_order_indices] <- "3rd Order (k=48)"  # Third-order neighbors

# Create data frame for plotting
neighbors_df <- data.frame(spatialCoords(spe), neighbors = spe$neighbors)

# Determine the min distance for coordinate scaling
min_distance <- min(dist(spatialCoords(spe)[dnn1[1,], ])) + 1

# Set x and y limits for the plot, centered around the first spot
xlim <- c(spatialCoords(spe)[1, 1] - 5 * min_distance, spatialCoords(spe)[1, 1] + 5 * min_distance)
ylim <- c(spatialCoords(spe)[1, 2] - 5 * min_distance, spatialCoords(spe)[1, 2] + 5 * min_distance)
neighbors2plot <- neighbors_df[
  neighbors_df$xcoord > xlim[1] & neighbors_df$xcoord < xlim[2] &
  neighbors_df$ycoord > ylim[1] & neighbors_df$ycoord < ylim[2], ]

# Scale coordinates for easier plotting
neighbors2plot$xcoord <- round((neighbors2plot$xcoord - min(neighbors2plot$xcoord)) / min_distance)
neighbors2plot$ycoord <- round((neighbors2plot$ycoord - min(neighbors2plot$ycoord)) / min_distance)

# Plot with distinct colors for each neighbor level
png(here(plot_dir, "Square_neighbor_parameters.png"), width=6, height=4, units="in", res=300)
ggplot(neighbors2plot, aes(x = xcoord, y = ycoord)) +
  geom_tile(aes(fill = neighbors, width = 0.95, height = 0.95)) +
  scale_fill_manual(values = c( "lightblue", "dodgerblue", "darkblue", "red","grey")) +
  labs(title = "Square grids (VisiumHD)") +
  ggpubr::theme_pubr() +
  theme(legend.position = "right") +
    theme(plot.title = element_text(size = 16))
dev.off()
