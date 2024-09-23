
library(here)
library(patchwork)
library(ggplot2)


# Set directories
plot_dir <- here("plots", "VisiumHD", "mouse_brain")
processed_dir <- here("processed-data", "VisiumHD")

# load
timing_results <- read.csv(file.path(processed_dir, "spotsweeper_timing_results_batch.csv"))

# Create the plot
png(file.path(plot_dir, "spotsweeper_timing_results_batch.png"), width = 5, height = 5, units = "in", res = 300)
ggplot(timing_results, aes(x = DatasetSize, y = TimeTaken)) +
  geom_jitter(size = 2, color = "red", width=50, alpha=.6) +  # Red points for individual samples
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "black") +  # Dashed linear trend line
  scale_x_continuous() +  # Continuous x-axis for number of spots
  scale_y_continuous() +  # Continuous y-axis for runtime
  labs(title = "Scalability: SpotSweeper, single metric (sum)",
       x = "Number of spots",
       y = "Runtime (sec)") +
  theme_minimal() +  # Minimal theme for clean look
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Centered and bold title
    axis.title = element_text(size = 12),  # Axis title size
    axis.text = element_text(size = 10),   # Axis text size
    panel.grid.major = element_line(size = 0.5, color = "grey"),  # Major gridlines
    panel.grid.minor = element_line(size = 0.25, color = "lightgrey")  # Minor gridlines
  )

dev.off()