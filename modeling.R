# ==========================================================================
# High-Dimensional MCTM Model Fitting
# ==========================================================================

# Load required libraries
library(tram)
library(dplyr)

#' Fit a high-dimensional MCTM model
#'
#' @param data Data frame containing the variables to be modeled
#' @param formula Model formula (default: ~ 1 for unconditional model)
#' @param var_names Names of variables to include in the model (optional)
#' @param dimensions Number of dimensions (optional if var_names is provided)
#' @param order Order of the Bernstein polynomial (default: 6)
#' @return A fitted mmlt model object, or NULL if fitting fails
fit_high_dim_mctm <- function(data, formula = ~ 1, var_names = NULL, 
                              dimensions = NULL, order = 6) {
  start_time <- proc.time()  
  tryCatch({
    # Determine variable names and dimensions
    if (is.null(var_names)) {
      if (is.null(dimensions)) {
        stop("Either var_names or dimensions must be provided")
      }
      var_names <- paste0("y", 1:dimensions)
    } else {
      dimensions <- length(var_names)
    }
    
    # Check if all variables are present in the data
    missing_vars <- setdiff(var_names, names(data))
    if (length(missing_vars) > 0) {
      stop("The following variables are missing from the data: ", paste(missing_vars, collapse=", "))
    }
    
    # Create marginal BoxCox models
    margin_models <- list()
    for (i in 1:dimensions) {
      var_name <- var_names[i]
      
      # Calculate support range
      support_range <- c(min(data[[var_name]], na.rm = TRUE), 
                         max(data[[var_name]], na.rm = TRUE))
      padding <- 0.1 * (support_range[2] - support_range[1])
      support_range <- support_range + c(-padding, padding)
      
      # Create model formula
      if (deparse(formula[[2]]) == "1") {
        model_formula <- as.formula(paste(var_name, "~", "1"))
      } else {
        model_formula <- as.formula(paste(var_name, "~", deparse(formula[[2]])))
      }
      
      # Create model
      margin_models[[i]] <- BoxCox(
        model_formula, 
        data = data,
        support = support_range,
        order = order,
        add = c(-0.5, 0.5)
      )
    }
    
    # Fit multivariate model
    args <- c(margin_models, list(formula = formula, data = data))
    mctm_model <- do.call(mmlt, args)
    
    # Store variable names in the model object for future reference
    attr(mctm_model, "var_names") <- var_names
    
    end_time <- proc.time()
    elapsed_time <- end_time - start_time
    print(sprintf("Function execution time: %.2f seconds", elapsed_time["elapsed"]))
    
    return(mctm_model)
  }, error = function(e) {
    message("Error during model fitting: ", e$message)
    return(NULL)
  })
}

# ==========================================================================
# Data Preparation for Covertype Dataset
# ==========================================================================

set.seed(123) # Set seed for reproducibility

# Read Covertype dataset
cat("Loading Covertype dataset...\n")
covtype_data <- read.csv("covtype.data.gz", header = FALSE)

# Define column names according to UCI Covertype dataset documentation
col_names <- c(
  # Continuous variables
  "Elevation", "Aspect", "Slope", "Horizontal_Distance_To_Hydrology", 
  "Vertical_Distance_To_Hydrology", "Horizontal_Distance_To_Roadways", 
  "Hillshade_9am", "Hillshade_Noon", "Hillshade_3pm", 
  "Horizontal_Distance_To_Fire_Points",
  # Remaining variables
  paste0("V", 11:55))

# Set column names
colnames(covtype_data) <- col_names

# Keep only the first 10 continuous variables
covtype_continuous <- covtype_data[, 1:10]

# Standardize continuous variables
covtype_scaled <- as.data.frame(scale(covtype_continuous))

# Add Cover_Type column
covtype_scaled <- covtype_scaled %>%
  mutate(Cover_Type = covtype_data$V55)

# View data summary
cat("Covertype dataset summary:\n")
print(dim(covtype_scaled))
print(summary(covtype_scaled[, 1:5]))

# Sample from the data - using a smaller sample size for computational efficiency
sample_size <- 30000 # Reduced sample size for faster processing
covtype_sample <- covtype_scaled[sample(nrow(covtype_scaled), sample_size), ]

# ==========================================================================
# Model Fitting
# ==========================================================================

# Fit unconditional model
cat("\nFitting unconditional model...\n")
unconditional_model <- fit_high_dim_mctm(
  data = covtype_sample,
  var_names = names(covtype_sample)[1:10],
  formula = ~ 1
)

# Fit conditional model (based on Cover_Type)
cat("\nFitting conditional model based on Cover_Type...\n")
conditional_model <- fit_high_dim_mctm(
  data = covtype_sample,
  var_names = names(covtype_sample)[1:10],
  formula = ~ Cover_Type
)
