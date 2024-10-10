
library(here)
library(patchwork)
library(ggplot2)


# Set directories
plot_dir <- here("plots", "VisiumHD", "mouse_brain")
processed_dir <- here("processed-data", "VisiumHD")

# load
timing_results <- read.csv(file.path(processed_dir, "spotsweeper_timing_results_batch_008.csv"))
timing_results 
#    DatasetSize   TimeTaken
# 1        1e+03   1.1005584
# 2        1e+03   1.0701233
# 3        1e+03   0.9950546
# 4        1e+03   1.0105378
# 5        1e+03   1.0213042

# 1000 dataset size
timing_results_1000 <- timing_results[!timing_results$DatasetSize == 1e+03,]

# List of dataset sizes to test
dataset_sizes <- c(10000, 50000, 100000, 400000)

# Create the plot with log scale for x-axis, vertical x-axis labels, and only major ticks
png(file.path(plot_dir, "spotsweeper_timing_results_008.png"), width = 4.5, height = 4.5, units = "in", res = 300)
ggplot(timing_results, aes(x = DatasetSize, y = TimeTaken)) +
  geom_jitter(size = 5, color = "red", width=200, alpha=.2) +  # Red points for individual samples
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "black") +  # Dashed linear trend line
  scale_x_continuous(breaks = dataset_sizes, labels = scales::comma) +  # Log scale for x-axis with custom breaks
  scale_y_continuous() +  # Continuous y-axis for runtime
  labs(title = "Scalability:",
      subtitle =" SpotSweeper, single metric",
       x = "Number of spots",
       y = "Runtime (sec)") +
  ggpubr::theme_pubr() +  # Use a cleaner classic theme
  theme(
    plot.title = element_text( size = 18),  # Centered and bold title
    plot.subtitle = element_text( size = 14),  # Centered subtitle
    axis.title = element_text(size = 14),  # Axis title size
    axis.text = element_text(size = 12),   # Axis text size
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),  # Rotate x-axis labels
    panel.grid.major = element_line(linetype=3,size = 0.5, color = "grey"),  # Major gridlines only
    panel.grid.minor = element_blank()  # Remove minor gridlines
  )

dev.off()
