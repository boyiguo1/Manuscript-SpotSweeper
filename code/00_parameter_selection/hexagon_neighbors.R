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

spe <- readRDS(here("processed-data", "Visium", "breast_cancer", "BreastCancer_banksy_spotsweeper.rds"))
spe

# rename spatialCoords to "xcoord" and "ycoord"
colnames(spatialCoords(spe)) <- c("xcoord", "ycoord")

# Find 1st-order neighbors (k = 8)
dnn1 <- BiocNeighbors::findKNN(spatialCoords(spe), k = 6, warn.ties = FALSE)$index
spe$neighbors <- "Non-neighbor"
spe$neighbors[62] <- "Center"  # Center
spe$neighbors[dnn1[62,]] <- "1st Order (k=6)"  # First-order neighbors

# Find 2nd-order neighbors (k = 24) based on first-order neighbors
dnn2 <- BiocNeighbors::findKNN(spatialCoords(spe), k = 18, warn.ties = FALSE)$index

second_order_indices <- dnn2[62,!(dnn2[62,] %in% dnn1[62,]) ]
spe$neighbors[second_order_indices] <- "2nd Order (k=18)"  # Second-order neighbors

# Find 3rd-order neighbors (k = 48) based on second-order neighbors
dnn3 <- BiocNeighbors::findKNN(spatialCoords(spe), k = 36, warn.ties = FALSE)$index

third_order_indices <- dnn3[62, !(dnn3[62,] %in% dnn2[62,])]
spe$neighbors[third_order_indices] <- "3rd Order (k=36)"  # Third-order neighbors

# Create data frame for plotting
neighbors_df <- data.frame(spatialCoords(spe), neighbors = spe$neighbors)

# Determine the min distance for coordinate scaling
min_distance <- min(dist(spatialCoords(spe)[dnn1[62,], ])) + 1

# Set x and y limits for the plot, centered around the first spot
xlim <- c(spatialCoords(spe)[62, 1] - 6 * min_distance, spatialCoords(spe)[62, 1] + 6 * min_distance)
ylim <- c(spatialCoords(spe)[62, 2] - 6 * min_distance, spatialCoords(spe)[62, 2] + 6 * min_distance)
neighbors2plot <- neighbors_df[
  neighbors_df$xcoord > xlim[1] & neighbors_df$xcoord < xlim[2] &
  neighbors_df$ycoord > ylim[1] & neighbors_df$ycoord < ylim[2], ]

# Scale coordinates for easier plotting
neighbors2plot$xcoord <- (neighbors2plot$xcoord - min(neighbors2plot$xcoord)) / min_distance
neighbors2plot$ycoord <- (neighbors2plot$ycoord - min(neighbors2plot$ycoord)) / min_distance

# Plot with distinct colors for each neighbor level
png(here(plot_dir, "Hexagon_neighbor_parameters.png"), width=6, height=4, units="in", res=300)
ggplot(neighbors2plot, aes(x = xcoord, y = ycoord)) +
  geom_point(aes(color = neighbors), size=7.5) +
  scale_color_manual(values = c("lightblue", "dodgerblue", "darkblue", "red","grey")) +
  labs(title = "Hexagonal grids (Visium)") +
  ggpubr::theme_pubr() +
  theme(legend.position = "right") +
  theme(plot.title = element_text(size = 16))
dev.off()
