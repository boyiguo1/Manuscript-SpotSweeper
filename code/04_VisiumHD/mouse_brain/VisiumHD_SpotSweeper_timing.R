library(SpotSweeper)
library(RANN)
library(here)
library(SpatialExperiment)
library(spatialLIBD)
library(escheR)
library(patchwork)
library(scran)
library(scater)
library(microbenchmark)
library(dplyr)

# Set directories
plot_dir <- here("plots", "VisiumHD", "mouse_brain")
processed_dir <- here("processed-data", "VisiumHD")

# Load the data
spe <- readRDS(here("processed-data", "VisiumHD", "VisiumHD_MouseBrain_008.rds"))

# Add QC metrics as in your original code
rownames(spe) <- rowData(spe)$symbol
is.mito <- grepl("^mt-", rownames(spe))

# Get QC metrics
df <- scuttle::perCellQCMetrics(spe, subsets = list(Mito = is.mito))

# Add to colData
spe$sum <- df$sum
spe$detected <- df$detected
spe$subsets_mito_percent <- df$subsets_Mito_percent

# Function to subset SPE and time the localOutliers execution
run_spotsweeper_timing <- function(spe, n_spots, n_runs) {
  times <- numeric(n_runs)  # Store runtimes for each repetition
  
  for (i in 1:n_runs) {
    # Subset the SPE object
    if (n_spots < ncol(spe)) {
      spe_subset <- spe[, sample(ncol(spe), n_spots)]
    } else {
      spe_subset <- spe
    }
    
    # Time the localOutliers function for the "sum" metric
    time_taken <- microbenchmark(
      localOutliers(spe_subset, metric = "sum", direction = "lower", log = TRUE),
      times = n_runs  # Run once for this iteration
    )
    
    # Store the time in seconds
    times[i] <- time_taken$time / 1e9  # Convert from nanoseconds to seconds
  }
  
  # Return all runtimes
  return(times)
}


# List of dataset sizes to test
dataset_sizes <- c(1000, 10000, 50000, 100000, 400000)

# Initialize a data frame to store all timing results
timing_results <- data.frame()

# Loop through dataset sizes and collect timing information
for (size in dataset_sizes) {
  times <- run_spotsweeper_timing(spe, size, n_runs=5)
  size_results <- data.frame(DatasetSize = size, TimeTaken = times)
  timing_results <- rbind(timing_results, size_results)
}

# Save the results to a .csv file
write.csv(timing_results, file = file.path(processed_dir, "spotsweeper_timing_results_batch_008.csv"), row.names = FALSE)

# read in csv
timing_results <- read.csv(file.path(processed_dir, "spotsweeper_timing_results_batch_008.csv"))
# Output timing results
print(timing_results)


# drop the 1000 dataset size
timing_results <- timing_results[timing_results$DatasetSize != 1000, ]

# Create the plot with log scale for x-axis, vertical x-axis labels, and only major ticks
png(file.path(plot_dir, "spotsweeper_timing_results_008.png"), width = 4, height = 4, units = "in", res = 300)
ggplot(timing_results, aes(x = DatasetSize, y = TimeTaken)) +
  geom_jitter(size = 5, color = "red", width=200, alpha=.2) +  # Red points for individual samples
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "black") +  # Dashed linear trend line
  scale_x_continuous(breaks = dataset_sizes, labels = scales::comma) +  # Log scale for x-axis with custom breaks
  scale_y_continuous() +  # Continuous y-axis for runtime
  labs(title = "Scalability: SpotSweeper, single metric",
       x = "Number of spots",
       y = "Runtime (sec)") +
  ggpubr::theme_pubr() +  # Use a cleaner classic theme
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),  # Centered and bold title
    axis.title = element_text(size = 14),  # Axis title size
    axis.text = element_text(size = 12),   # Axis text size
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),  # Rotate x-axis labels
    panel.grid.major = element_line(linetype=3,size = 0.5, color = "grey"),  # Major gridlines only
    panel.grid.minor = element_blank()  # Remove minor gridlines
  )

dev.off()
