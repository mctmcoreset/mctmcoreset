# ==========================================================================
# GGpairs Visualization for Original Covertype Dataset
# ==========================================================================

# --- Load required packages ---
library(ggplot2)
library(GGally)
library(dplyr)
library(reshape2)
library(viridis)

# --- Load and prepare data ---
covtype_data <- read.csv("covtype.data.gz", header = FALSE)

# Set column names for continuous variables
col_names <- c(
  "Elevation", "Aspect", "Slope", "Horizontal_Distance_To_Hydrology", 
  "Vertical_Distance_To_Hydrology", "Horizontal_Distance_To_Roadways", 
  "Hillshade_9am", "Hillshade_Noon", "Hillshade_3pm", 
  "Horizontal_Distance_To_Fire_Points",
  paste0("V", 11:55)
)
colnames(covtype_data) <- col_names

# Extract continuous variables
continuous_vars <- covtype_data[, 1:10]

# Set sample size and seed for reproducibility
sample_size <- 2000  # Using 2000 sample points
set.seed(123)
sample_indices <- sample(1:nrow(continuous_vars), sample_size)
sample_data <- continuous_vars[sample_indices, ]

# ==========================================================================
# Basic GGpairs Plot
# ==========================================================================

# Create basic ggpairs plot
ggpairs_basic <- ggpairs(sample_data)

# Save the image
ggsave("covtype_ggpairs_basic.png", ggpairs_basic, width = 15, height = 15, dpi = 300)

# ==========================================================================
# Customized GGpairs Plots
# ==========================================================================

# Create more customized ggpairs plots, showing only selected variables
selected_vars <- c("Elevation", "Aspect", "Slope", "Horizontal_Distance_To_Hydrology",
                   "Vertical_Distance_To_Hydrology")
sample_selected <- sample_data[, selected_vars]

# Custom ggpairs for first 5 variables
custom_ggpairs <- ggpairs(
  sample_data[, 1:5],
  upper = list(continuous = wrap("cor", method = "spearman", size = 3)),
  lower = list(continuous = wrap("points", alpha = 0.3, size = 0.5)),
  diag = list(continuous = wrap("densityDiag", fill = "lightblue")),
  title = "Pairwise Relationships Between Terrain Variables") + 
  theme_bw() + 
  theme(
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 9),
    strip.text = element_text(size = 8),
    plot.title = element_text(size = 14, hjust = 0.5)
  )

# Save customized plot
ggsave("covtype_ggpairs_custom.png", custom_ggpairs, width = 12, height = 10, dpi = 300)

# Custom ggpairs for variables 6-10
custom_ggpairs2 <- ggpairs(
  sample_data[, 6:10],
  upper = list(continuous = wrap("cor", method = "spearman", size = 3)),
  lower = list(continuous = wrap("points", alpha = 0.3, size = 0.5)),
  diag = list(continuous = wrap("densityDiag", fill = "lightblue")),
  title = "Pairwise Relationships Between Terrain Variables") + 
  theme_bw() + 
  theme(
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 9),
    strip.text = element_text(size = 8),
    plot.title = element_text(size = 14, hjust = 0.5)
  )

# Save second customized plot
ggsave("covtype_ggpairs_custom2.png", custom_ggpairs2, width = 12, height = 10, dpi = 300)

# ==========================================================================
# GGpairs Colored by Cover Type
# ==========================================================================

# Extract cover type and create a dataset with cover type
cover_type <- covtype_data[, 55]
sample_data_with_type <- continuous_vars[sample_indices, selected_vars]
sample_data_with_type$Cover_Type <- factor(cover_type[sample_indices])

# Create ggpairs with coloring by cover type
colored_ggpairs <- ggpairs(
  sample_data_with_type,
  columns = 1:5,  
  mapping = aes(color = Cover_Type, alpha = 0.5),
  upper = list(continuous = wrap("cor", size = 3)),
  lower = list(continuous = wrap("points", size = 0.5)),
  diag = list(continuous = wrap("densityDiag")),
  title = "Terrain Variables Colored by Forest Cover Type") + 
  theme_bw() + 
  theme(
    legend.position = "bottom",
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 9),
    strip.text = element_text(size = 8),
    plot.title = element_text(size = 14, hjust = 0.5)
  )

# Save plot colored by cover type
ggsave("covtype_ggpairs_by_cover_type.png", colored_ggpairs, width = 12, height = 12, dpi = 300)

# ==========================================================================
# Correlation Heatmap
# ==========================================================================

# Calculate correlation matrix
corr_matrix <- cor(sample_data, method = "spearman")
corr_data <- melt(corr_matrix)
colnames(corr_data) <- c("Var1", "Var2", "Correlation")

# Create heatmap
correlation_heatmap <- ggplot(corr_data, aes(x = Var1, y = Var2, fill = Correlation)) +
  geom_tile() +
  scale_fill_viridis() +
  geom_text(aes(label = round(Correlation, 2)),
            color = ifelse(abs(corr_data$Correlation) > 0.7, "white", "black"),
            size = 3) +
  labs(title = "Correlation Heatmap of Terrain Variables",
       x = "", y = "") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    plot.title = element_text(hjust = 0.5, size = 14)
  )

# Save heatmap
ggsave("covtype_correlation_heatmap.png", correlation_heatmap, width = 10, height = 8, dpi = 300)

# ==========================================================================
# Parallel Coordinate Plot
# ==========================================================================

# Standardize data for parallel coordinate plot
scaled_data <- as.data.frame(scale(sample_data))
scaled_data$Cover_Type <- factor(cover_type[sample_indices])

# Create parallel coordinate plot
parallel_plot <- ggparcoord(
  scaled_data,
  columns = 1:10,  # Use all 10 continuous variables
  groupColumn = "Cover_Type",
  scale = "uniminmax",  # Unified min-max scaling
  alphaLines = 0.3) +
  labs(title = "Parallel Coordinate Plot of Terrain Variables by Cover Type",
       x = "Variables", y = "Scaled Value") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    legend.position = "bottom"
  )

# Save parallel coordinate plot
ggsave("covtype_parallel_plot.png", parallel_plot, width = 12, height = 8, dpi = 300)