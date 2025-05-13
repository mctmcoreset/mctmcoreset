.models <- function (...) 
{
  m <- lapply(list(...), function(x) as.mlt(x))
  nm <- abbreviate(sapply(m, function(x) x$model$response), 
                   4)
  J <- length(m)
  Jp <- J * (J - 1)/2
  normal <- sapply(m, function(x) x$todistr$name == "normal")
  w <- lapply(m, weights)
  out <- lapply(w, function(x) stopifnot(isTRUE(all.equal(x, 
                                                          w[[1]]))))
  w <- w[[1L]]
  if (isTRUE(all.equal(unique(w), 1))) 
    w <- FALSE
  mm <- lapply(m, function(mod) {
    eY <- get("eY", environment(mod$parm))
    iY <- get("iY", environment(mod$parm))
    list(eY = eY, iY = iY)
  })
  cmod <- sapply(mm, function(x) !is.null(x$eY))
  dmod <- sapply(mm, function(x) !is.null(x$iY))
  stopifnot(all(xor(cmod, dmod)))
  stopifnot(all(diff(cmod) <= 0))
  stopifnot(all(diff(dmod) >= 0))
  nobs <- unique(sapply(m, nobs))
  stopifnot(length(nobs) == 1L)
  nobs <- nobs[[1L]]
  P <- sapply(m, function(x) length(coef(x)))
  fpar <- factor(rep(1:J, P))
  parm <- function(par) {
    mpar <- par[1:sum(P)]
    split(mpar, fpar)
  }
  constr <- lapply(mm, function(m) {
    if (is.null(m$eY)) 
      return(attr(m$iY$Yleft, "constraint"))
    return(attr(m$eY$Y, "constraint"))
  })
  ui <- do.call("bdiag", lapply(constr, function(x) x$ui))
  ci <- do.call("c", lapply(constr, function(x) x$ci))
  ui <- as(ui[is.finite(ci), , drop = FALSE], "matrix")
  ci <- ci[is.finite(ci)]
  mf <- lapply(1:J, function(j) {
    mf <- m[[j]]$data
    if (cmod[j]) 
      return(mf)
    yl <- m[[j]]$response$cleft
    yr <- m[[j]]$response$cright
    rp <- m[[j]]$model$response
    ml <- mr <- mf
    ml[[rp]] <- yl
    mr[[rp]] <- yr
    return(list(left = ml, right = mr))
  })
  nn <- sapply(1:J, function(j) {
    !is.null(m[[j]]$fixed) || !isTRUE(all.equal(unique(m[[j]]$offset), 
                                                0)) || m[[j]]$model$scale_shift
  })
  type <- lapply(1:J, function(j) mlt:::.type_of_response(m[[j]]$response))
  return(list(models = m, mf = mf, cont = cmod, type = type, 
              normal = normal, nobs = nobs, weights = w, nparm = P, 
              parm = parm, ui = ui, ci = ci, mm = mm, names = nm, nn = nn))
}


# ==========================================================================
# Coreset Comparison Experiment Pipeline - High-Dimensional MCTM (Covertype Dataset)
# ==========================================================================
# --- 0. Load Required Libraries ---
# Make sure these packages are installed: install.packages(c("tram", "mvtnorm", "ggplot2", "Matrix", "dplyr", "geometry", "parallel", "gridExtra"))
suppressPackageStartupMessages(library(tram))       # For fitting the MCTM model
suppressPackageStartupMessages(library(mvtnorm))    # Possibly used for computations, although handled internally by mmlt
suppressPackageStartupMessages(library(ggplot2))    # For plotting
suppressPackageStartupMessages(library(Matrix))     # For sparse matrix operations
suppressPackageStartupMessages(library(dplyr))      # For data manipulation
suppressPackageStartupMessages(library(geometry))   # For computing convex hulls
suppressPackageStartupMessages(library(parallel))   # For parallel computation to speed up experiments
suppressPackageStartupMessages(library(gridExtra))  # For arranging ggplot figures

# --- Global Settings ---
# Set working directory (if needed)
# setwd("/path/to/your/directory")

# Define number of cores for parallel computing (adjust according to your machine)
# num_cores <- detectCores() - 1
num_cores <- 2# Can be set manually

#' Precompute Convex Hull vertex indices for all dimensions
#'
#' @param A_prime_list A list containing A'_j matrices for all dimensions
#' @param n Total number of data points (used for fallback)
#' @param k_hull_fallback_size Number of points to sample if fallback is needed
#' @return A list where each element is a vector of convex hull vertex indices for one dimension
precompute_all_hull_indices <- function(A_prime_list, n, 
                                        k_hull_fallback_size = 50) {
  dimensions <- length(A_prime_list)
  print(dimensions)
  cat("Starting precomputation of Convex Hull indices for all dimensions...\n")
  
  # Use mclapply to compute in parallel
  hull_indices_list <- mclapply(1:dimensions, function(j) {
    cat(sprintf("  Computing Convex Hull for dimension %d...\n", j))
    A_prime_j <- A_prime_list[[j]]
    
    # Check if input is valid
    if (is.null(A_prime_j) || nrow(A_prime_j) < ncol(A_prime_j) + 1) {
      warning(sprintf("  Dimension %d: Matrix A' is invalid or has insufficient points, using random fallback.", j))
      return(sample(1:n, min(n, k_hull_fallback_size)))
    }
    
    # Try to compute Convex Hull with error handling
    hull_j_indices <- tryCatch({
      # Adding tiny noise may help avoid collinearity issues
      noise <- matrix(rnorm(nrow(A_prime_j) * ncol(A_prime_j), sd = 1e-9), 
                      nrow = nrow(A_prime_j))
      unique(as.vector(chull(A_prime_j + noise)))
    }, error = function(e) {
      warning(sprintf("  Dimension %d: Convex Hull computation failed: %s. Using random fallback.", j, e$message))
      return(sample(1:n, min(n, k_hull_fallback_size)))
    })
    
    cat(sprintf("  Dimension %d: Convex Hull computed successfully with %d vertices.\n", j, length(hull_j_indices)))
    return(hull_j_indices)
  }, mc.cores = num_cores) # Parallel computation
  
  cat("Precomputation of Convex Hull indices for all dimensions completed.\n")
  names(hull_indices_list) <- paste0("dim_", 1:dimensions)
  return(hull_indices_list)
}


# ==========================================================================
# Part 1: Load and Prepare Covertype Data
# ==========================================================================

#' Load and preprocess the Covertype dataset
#'
#' @param n_sample Number of samples to draw
#' @param seed Random seed
#' @param file_path Path to the data file (default is "covtype.data.gz")
#' @return A data frame containing 10 standardized continuous variables, with column names y1 to y10
load_prepare_covertype <- function(n_sample = 100000, seed = 123, file_path = "covtype.data.gz") {
  cat("Loading Covertype dataset from:", file_path, "...\n")
  if (!file.exists(file_path)) {
    stop("Dataset file not found: ", file_path,
         "\nPlease download it from: https://archive.ics.uci.edu/ml/machine-learning-databases/covtype/covtype.data.gz")
  }
  # Use data.table::fread for faster reading
  if (!requireNamespace("data.table", quietly = TRUE)) {
    message("For faster loading, it is recommended to install the 'data.table' package. Using read.csv instead.")
    covtype_data <- read.csv(gzfile(file_path), header = FALSE)
  } else {
    covtype_data <- data.table::fread(file_path, header = FALSE)
  }
  
  # Define column names (according to UCI documentation)
  col_names <- c(
    "Elevation", "Aspect", "Slope", "Horizontal_Distance_To_Hydrology",
    "Vertical_Distance_To_Hydrology", "Horizontal_Distance_To_Roadways",
    "Hillshade_9am", "Hillshade_Noon", "Hillshade_3pm",
    "Horizontal_Distance_To_Fire_Points",
    paste0("Wilderness_Area_", 1:4), # 4 Wilderness Area binary columns
    paste0("Soil_Type_", 1:40),     # 40 Soil Type binary columns
    "Cover_Type"                   # Target variable (categorical)
  )
  colnames(covtype_data) <- col_names
  
  # Select only the first 10 continuous variables
  covtype_continuous <- covtype_data[, 1:10, with = FALSE] # Using data.table syntax
  cat("Selected the first 10 continuous variables.\n")
  
  # Standardization (Scaling)
  cat("Standardizing variables...\n")
  covtype_scaled <- as.data.frame(scale(covtype_continuous))
  
  # Free memory
  rm(covtype_data, covtype_continuous)
  gc()
  
  # Random sampling
  cat(sprintf("Randomly sampling %d observations...\n", n_sample))
  set.seed(seed)
  if (n_sample >= nrow(covtype_scaled)) {
    warning("Requested sample size is greater than or equal to the dataset size. Using the full dataset.")
    sample_indices <- 1:nrow(covtype_scaled)
  } else {
    sample_indices <- sample(nrow(covtype_scaled), n_sample)
  }
  covtype_final_sample <- covtype_scaled[sample_indices, ]
  
  cat("Data preparation complete.\n")
  cat("Final sample dimensions:", dim(covtype_final_sample), "\n")
  cat("Summary of the first 5 columns:\n")
  print(summary(covtype_final_sample[, 1:5]))
  
  # Rename columns to y1, y2, ..., y10
  new_colnames <- paste0("y", 1:10)
  colnames(covtype_final_sample) <- new_colnames
  cat("Columns renamed to y1 through y10.\n")
  
  return(covtype_final_sample)
}


# ==========================================================================
# Part 2: Core MCTM Functions and Matrix Extraction
# ==========================================================================

#' Fit the full multivariate MCTM model (used as benchmark)
#'
#' @param data A data frame containing the response variables (y1, y2, ...)
#' @param formula Model formula (e.g., ~ 1 for unconditional)
#' @param order Order of the Bernstein polynomial
#' @return A fitted mmlt object, or NULL if fitting fails
fit_full_mctm <- function(data, formula = ~ 1, order = 6) {
  cat("Fitting full MCTM model (benchmark model)...\n")
  start_time <- proc.time()
  
  dimensions <- ncol(data)
  var_names <- names(data)
  
  # Create marginal models in parallel
  margin_models <- mclapply(1:dimensions, function(i) {
    var_name <- var_names[i]
    # Compute support range and add padding
    support_range <- range(data[[var_name]], na.rm = TRUE)
    padding <- 0.1 * diff(support_range)
    # Handle the case where padding is zero
    if (padding == 0) padding <- 1 # or another small positive number
    support_range_padded <- support_range + c(-padding, padding)
    
    # Handle Inf/NaN in support
    if(any(!is.finite(support_range_padded))) {
      warning(paste("Support range of variable", var_name, "contains Inf/NaN, will use mean +/- 10*SD."))
      m <- mean(data[[var_name]], na.rm=TRUE)
      s <- sd(data[[var_name]], na.rm=TRUE)
      if (is.na(s) || s == 0) s <- 1 # Fallback if sd is NA or 0
      support_range_padded <- c(m - 10*s, m + 10*s)
    }
    
    model_formula <- as.formula(paste(var_name, "~", deparse1(formula[[2]])))
    
    # Use tryCatch to handle fitting errors for individual models
    tryCatch({
      BoxCox(model_formula, data = data, support = support_range_padded,
             order = order, add = c(-0.5, 0.5)) # use default 'add'
    }, error = function(e) {
      message(paste("Error fitting marginal model for variable", var_name, ":", e$message))
      return(NULL) # Return NULL to indicate failure
    })
  }, mc.cores = num_cores) # Parallelization using mclapply
  
  # Filter out models that failed to fit (NULL)
  valid_models_idx <- !sapply(margin_models, is.null)
  if(!all(valid_models_idx)) {
    failed_vars <- var_names[!valid_models_idx]
    warning("Marginal model fitting failed for the following variables: ", paste(failed_vars, collapse=", "), ". Remaining models will be used.")
    margin_models <- margin_models[valid_models_idx]
    if(length(margin_models) < 2) {
      stop("Too many marginal model fitting failures (<2). Cannot proceed with mmlt fitting.")
    }
    # Note: dimensions of the data remain unchanged, only the number of valid models changes
  }
  
  # Fit the mmlt model
  args <- c(margin_models, list(formula = formula, data = data))
  mctm_model <- tryCatch({
    do.call(mmlt, args)
  }, error = function(e) {
    message("Error fitting mmlt multivariate model: ", e$message)
    return(NULL)
  })
  
  elapsed_time <- (proc.time() - start_time)[["elapsed"]]
  if (!is.null(mctm_model)) {
    cat(sprintf("Full MCTM model successfully fitted in %.2f seconds.\n", elapsed_time))
  } else {
    cat(sprintf("Fitting full MCTM model failed. Time elapsed: %.2f seconds.\n", elapsed_time))
  }
  return(mctm_model)
}

#' Extract basis function matrices A_j and derivative matrices A'_j from an MCTM model
#'
#' This function builds temporary marginal models only to extract
#' internal matrix structures needed for coreset sampling.
#'
#' @param data The full dataset (used for support range and temporary fitting)
#' @param formula Model formula (default ~ 1)
#' @param order Bernstein polynomial order (default 6)
#' @return A list containing:
#'         $A : list of A_j matrices (basis functions)
#'         $A_prime : list of A'_j matrices (derivatives)
#'         $dimensions : number of dimensions for which matrices were successfully extracted
extract_mctm_matrices <- function(data, formula = ~ 1, order = 6) {
  cat("Extracting basis function matrices A_j and derivative matrices A'_j...\n")
  dimensions <- ncol(data)
  var_names <- names(data)
  
  # 1. Build temporary marginal models (in parallel)
  #    The goal is to extract internal structure (eY$Y, eY$Yprime)
  cat("  Building temporary marginal models (parallel)...\n")
  margin_models_for_structure <- mclapply(1:dimensions, function(i) {
    var_name <- var_names[i]
    support_range <- range(data[[var_name]], na.rm = TRUE)
    padding <- 0.1 * diff(support_range)
    if (is.na(padding) || padding == 0) padding <- 1
    support_range_padded <- support_range + c(-padding, padding)
    
    # Handle Inf/NaN in support range
    if (any(!is.finite(support_range_padded))) {
      m <- mean(data[[var_name]], na.rm = TRUE)
      s <- sd(data[[var_name]], na.rm = TRUE)
      if (is.na(s) || s == 0) s <- 1
      support_range_padded <- c(m - 10 * s, m + 10 * s)
      message(sprintf("  Info: Non-finite support for variable %s. Using range: %s",
                      var_name, paste(round(support_range_padded, 2), collapse = " - ")))
    }
    
    model_formula <- as.formula(paste(var_name, "~", deparse1(formula[[2]])))
    
    # Fit temporary BoxCox model using tryCatch
    tmp_model <- tryCatch({
      tram::BoxCox(model_formula, data = data, support = support_range_padded,
             order = order, add = c(-0.5, 0.5))
    }, error = function(e) {
      message(sprintf("  Error while building temporary marginal model for %s: %s", var_name, e$message))
      return(NULL)
    })
    return(tmp_model)
  }, mc.cores = num_cores)
  
  # 2. Filter out failed models
  valid_models_idx <- !sapply(margin_models_for_structure, is.null)
  if (!any(valid_models_idx)) {
    stop("Failed to build ALL temporary marginal models. Cannot extract matrices.")
  }
  if (!all(valid_models_idx)) {
    failed_vars <- var_names[!valid_models_idx]
    warning("Failed to build temporary marginal models for variables: ",
            paste(failed_vars, collapse = ", "), ". Matrices will only be extracted for valid models.")
    margin_models_for_structure <- margin_models_for_structure[valid_models_idx]
  }
  
  # 3. Extract model structure via do.call
  cat("  Extracting internal model structures...\n")
  model_structures <- tryCatch({
    do.call(.models, margin_models_for_structure)
  }, error = function(e) {
    message("Error calling tram:::.models via do.call: ", e$message)
    return(NULL)
  })
  
  if (is.null(model_structures) || length(model_structures$models) == 0) {
    stop("Failed to extract matrix structure from valid marginal models.")
  }
  
  # 4. Extract A and A' matrices from structure
  A_list <- list()
  A_prime_list <- list()
  actual_dimensions <- length(model_structures$models)
  
  for (j in 1:actual_dimensions) {
    model_mm <- model_structures$mm[[j]]
    if (!is.null(model_mm) && !is.null(model_mm$eY)) {
      if (!is.null(model_mm$eY$Y)) {
        A_list[[j]] <- model_mm$eY$Y
      } else {
        warning(sprintf("Matrix A (eY$Y) not found for dimension %d (relative to valid models).", j))
        A_list[[j]] <- NULL
      }
      if (!is.null(model_mm$eY$Yprime)) {
        A_prime_list[[j]] <- model_mm$eY$Yprime
      } else {
        warning(sprintf("Matrix A' (eY$Yprime) not found for dimension %d (relative to valid models).", j))
        A_prime_list[[j]] <- NULL
      }
    } else {
      warning(sprintf("Structure 'mm' or 'eY' missing for dimension %d (relative to valid models).", j))
      A_list[[j]] <- NULL
      A_prime_list[[j]] <- NULL
    }
  }
  
  # 5. Finalize and return
  A_list_final <- Filter(Negate(is.null), A_list)
  A_prime_list_final <- Filter(Negate(is.null), A_prime_list)
  final_dimensions <- length(A_list_final)
  
  if (final_dimensions == 0) {
    stop("Failed to extract any valid A or A' matrices.")
  }
  if (final_dimensions < dimensions) {
    original_var_names_valid <- var_names[valid_models_idx]
    extracted_var_names <- names(A_list_final)
    if (is.null(extracted_var_names)) extracted_var_names <- original_var_names_valid
    warning(sprintf("Only %d matrix dimensions extracted (less than original %d). Variables: %s",
                    final_dimensions, dimensions, paste(extracted_var_names, collapse = ", ")))
  }
  
  cat("Matrix extraction complete. Valid dimensions:", final_dimensions, "\n")
  return(list(A = A_list_final, A_prime = A_prime_list_final, dimensions = final_dimensions))
}


#' # ==========================================================================
#' # Part 3: Coreset Sampling Methods (Adapted for High Dimensions)
#' # ==========================================================================

#' Combine probabilities from multiple sources
#'
#' @param prob_sources_list A list of probability vectors, each of length n
#' @param n Total number of observations
#' @return A combined probability vector of length n
combine_probabilities <- function(prob_sources_list, n) {
  combined_probs <- numeric(n)
  num_sources <- length(prob_sources_list)
  if (num_sources == 0) return(rep(1/n, n)) # If no sources, use uniform
  
  # Use log for numerical stability: log(1 - p_i) = sum_j log(1 - p_ij)
  log_one_minus_p_matrix <- matrix(0, nrow = n, ncol = num_sources)
  
  for (j in 1:num_sources) {
    probs_j <- prob_sources_list[[j]]
    # Ensure correct length and clean values
    if (length(probs_j) != n) stop("Probability vector length mismatch")
    probs_j[is.na(probs_j) | !is.finite(probs_j) | probs_j < 0] <- 0
    probs_j[probs_j >= 1] <- 1 - 1e-12 # Avoid log(0)
    log_one_minus_p_matrix[, j] <- log(1 - probs_j)
  }
  
  # Sum the logs
  sum_log_one_minus_p <- rowSums(log_one_minus_p_matrix)
  
  # Combined probability p_i = 1 - exp(sum_j log(1 - p_ij))
  combined_probs <- 1 - exp(sum_log_one_minus_p)
  
  # Remove small negative values due to precision
  combined_probs[combined_probs < 0] <- 0
  
  # Normalize to sum to 1
  total_prob <- sum(combined_probs)
  if (total_prob > 1e-12) { # Use small tolerance
    combined_probs <- combined_probs / total_prob
  } else {
    warning("Total combined probability is close to zero. Reverting to uniform.")
    return(rep(1/n, n))
  }
  
  # Ensure no NA/NaN/Inf after normalization
  combined_probs[!is.finite(combined_probs)] <- 0
  total_prob_final <- sum(combined_probs)
  if (total_prob_final > 1e-12) {
    combined_probs <- combined_probs / total_prob_final
  } else {
    return(rep(1/n, n))
  }
  
  return(combined_probs)
}

#' Uniform Sampling
#' @param n Total number of observations
#' @param k Coreset size
#' @param seed Random seed
#' @return A list containing indices, weights, and probs
sample_uniform <- function(n, k, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  indices <- sample(1:n, k, replace = FALSE)
  weights <- rep(n / k, k)
  probs <- rep(1/n, n)
  return(list(indices = sort(indices), weights = weights, probs = probs))
}

#' L2 Leverage Score Sampling Only (High-Dimensional)
#' @param A_list List of A_j matrices
#' @param k Coreset size
#' @param seed Random seed
#' @return A list containing indices, weights, and combined probs
sample_l2_only <- function(A_list, k, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (length(A_list) == 0) stop("A_list cannot be empty.")
  n <- nrow(A_list[[1]])
  dimensions <- length(A_list)
  
  cat("  Computing L2 leverage scores (parallel)...\n")
  # Compute L2 probabilities for each dimension in parallel
  prob_l2_list <- mclapply(1:dimensions, function(j) {
    Aj <- A_list[[j]]
    if (is.null(Aj) || nrow(Aj) == 0 || ncol(Aj) == 0) {
      warning(paste("Matrix A for dimension", j, "is invalid. Using uniform."))
      return(rep(1/n, n))
    }
    # QR decomposition with tryCatch
    qr_decomp <- tryCatch(qr(Aj), error = function(e) {
      warning(paste("QR decomposition failed for dimension", j, ":", e$message, ". Using uniform."))
      return(NULL)
    })
    if (is.null(qr_decomp)) return(rep(1/n, n))
    
    leverage_j <- tryCatch(rowSums(qr.Q(qr_decomp)^2), error = function(e) {
      warning(paste("Leverage computation failed for dimension", j, ":", e$message, ". Using uniform."))
      return(rep(0, n)) # Fallback to zero (will become uniform after normalization)
    })
    
    prob_j <- (leverage_j + 1/n) / (sum(leverage_j) + 1)
    prob_j[is.na(prob_j) | !is.finite(prob_j)] <- 0
    return(prob_j)
  }, mc.cores = num_cores)
  
  # Remove NULL results
  prob_l2_list <- Filter(Negate(is.null), prob_l2_list)
  if (length(prob_l2_list) == 0) stop("All L2 probability computations failed.")
  
  cat("  Combining L2 probabilities...\n")
  combined_probs <- combine_probabilities(prob_l2_list, n)
  
  cat("  Performing L2-only sampling...\n")
  indices <- sample(1:n, k, replace = FALSE, prob = combined_probs)
  
  # Compute weights: w_i = 1 / (p_i * k), normalized to sum(w_i) = n
  sampled_probs <- combined_probs[indices]
  weights <- rep(0, k)
  non_zero_prob_idx <- sampled_probs > 1e-12
  weights[non_zero_prob_idx] <- 1 / (sampled_probs[non_zero_prob_idx] * k)
  
  # Normalize weights
  current_sum_weights <- sum(weights)
  if (current_sum_weights > 1e-12) {
    weights <- weights * (n / current_sum_weights)
  } else {
    warning("Sum of L2 weights is close to zero. Using uniform weights.")
    weights <- rep(n / k, k)
  }
  
  return(list(indices = sort(indices), weights = weights, probs = combined_probs))
}

#' L2 Leverage Score + Convex Hull Sampling (High Dimension - Using Precomputed Indices)
#'
#' @param A_list List of A_j matrices
#' @param A_prime_list List of A'_j matrices (derivatives) - Actually no longer used here but kept for consistency
#' @param k Total coreset size
#' @param k_hull_prop Proportion of coreset size allocated to hull points
#' @param seed Random seed
#' @param precomputed_hull_indices List containing precomputed convex hull vertex indices for each dimension (output of precompute_all_hull_indices)
#' @return A list containing indices, weights, and combined probabilities
sample_l2_hull <- function(A_list, A_prime_list, k, k_hull_prop = 0.2, seed = NULL, precomputed_hull_indices) {
  # --- Validate Inputs ---
  if (!is.null(seed)) set.seed(seed)
  if (missing(precomputed_hull_indices)) {
    stop("Argument 'precomputed_hull_indices' must be provided.")
  }
  if (length(A_list) == 0) stop("A_list cannot be empty.")
  n <- nrow(A_list[[1]])
  dimensions <- length(A_list)
  # Ensure dimensions match
  if (length(precomputed_hull_indices) != dimensions) {
    warning("Number of dimensions in precomputed_hull_indices does not match A_list. Using the smaller one.")
    dimensions <- min(length(A_list), length(precomputed_hull_indices))
    if (dimensions == 0) stop("No matching dimensions between A_list and precomputed_hull_indices.")
  }
  
  # --- 1. Compute L2 Probabilities (same as sample_l2_only) ---
  cat("    Computing L2 leverage probabilities (parallel)...\n")
  prob_l2_list <- mclapply(1:dimensions, function(j) {
    Aj <- A_list[[j]]
    if (is.null(Aj) || nrow(Aj) == 0 || ncol(Aj) == 0) return(rep(1/n, n))
    qr_decomp <- tryCatch(qr(Aj), error = function(e) NULL)
    if (is.null(qr_decomp)) return(rep(1/n, n))
    leverage_j <- tryCatch(rowSums(qr.Q(qr_decomp)^2), error = function(e) rep(0, n))
    prob_j <- (leverage_j + 1/n) / (sum(leverage_j) + 1)
    prob_j[is.na(prob_j) | !is.finite(prob_j)] <- 0
    return(prob_j)
  }, mc.cores = num_cores)
  prob_l2_list <- Filter(Negate(is.null), prob_l2_list)
  if (length(prob_l2_list) == 0) stop("Failed to compute any L2 probabilities.")
  
  # --- 2. Use Precomputed Convex Hull Indices ---
  cat("    Using precomputed convex hull indices...\n")
  all_unique_hull_indices <- sort(unique(unlist(precomputed_hull_indices[1:dimensions])))
  n_unique_hull <- length(all_unique_hull_indices)
  cat(sprintf("      Total %d unique convex hull points from precomputation.\n", n_unique_hull))
  
  # --- 3. Build Hull Probability Vector from Precomputed Indices ---
  k_hull_target <- ceiling(k * k_hull_prop)
  k_hull_actual <- min(k_hull_target, n_unique_hull)
  
  prob_hull_vector <- numeric(n)
  if (n_unique_hull > 0 && k_hull_actual > 0) {
    individual_hull_prob <- k_hull_actual / n_unique_hull
    prob_hull_vector[all_unique_hull_indices] <- individual_hull_prob
  } else {
    cat("      No unique hull points or target hull size is zero from precomputation.\n")
  }
  
  # --- 4. Combine L2 and Hull Probabilities ---
  cat("    Combining L2 and Hull probabilities...\n")
  all_prob_sources <- c(prob_l2_list, list(prob_hull_vector))
  combined_probs <- combine_probabilities(all_prob_sources, n)
  
  # --- 5. Sample Based on Combined Probabilities ---
  cat("    Performing L2+Hull sampling...\n")
  indices <- tryCatch({
    sample(1:n, k, replace = FALSE, prob = combined_probs)
  }, error = function(e) {
    warning(paste("Sampling L2+Hull failed:", e$message, ". Falling back to uniform sampling."))
    return(sample(1:n, k, replace = FALSE))
  })
  
  # --- 6. Compute Weights ---
  sampled_combined_probs <- combined_probs[indices]
  weights <- rep(0, k)
  non_zero_prob_idx <- sampled_combined_probs > 1e-12
  
  weights[non_zero_prob_idx] <- 1 / (sampled_combined_probs[non_zero_prob_idx] * k)
  
  # Normalize weights
  current_sum_weights <- sum(weights[is.finite(weights)])
  if (current_sum_weights > 1e-12) {
    weights <- weights * (n / current_sum_weights)
  } else {
    warning("Total L2+Hull weights close to zero. Using uniform weights.")
    weights <- rep(n / k, k)
  }
  weights[!is.finite(weights)] <- mean(weights[is.finite(weights)], na.rm = TRUE)
  if (all(!is.finite(weights))) weights <- rep(n/k, k)
  
  return(list(indices = sort(indices), weights = weights, probs = combined_probs))
}

# ==========================================================================
# Part 4: Fitting the Coreset Model
# ==========================================================================

#' Fit MCTM model on coreset data with weights
#' @param data Full dataset
#' @param indices Coreset indices
#' @param weights Coreset weights
#' @param formula Model formula
#' @param order Bernstein order
#' @return mmlt model fitted on the coreset, or NULL if failed
fit_coreset_mctm <- function(data, indices, weights, formula = ~ 1, order = 6) {
  cat("  Fitting MCTM model on coreset...\n")
  start_time <- proc.time()
  
  # Ensure indices are valid
  if(any(indices > nrow(data)) || any(indices < 1)) {
    stop("Invalid coreset indices.")
  }
  coreset_data <- data[indices, , drop = FALSE]
  
  # Ensure weights are valid
  if(length(weights) != length(indices) || any(!is.finite(weights)) || any(weights < 0)) {
    warning("Invalid coreset weights. Using uniform weights instead.")
    weights <- rep(nrow(data) / length(indices), length(indices))
  }
  
  dimensions <- ncol(coreset_data)
  var_names <- names(coreset_data)
  
  # Create marginal models (must be rebuilt for the coreset data to ensure correct support)
  margin_models <- mclapply(1:dimensions, function(i) {
    var_name <- var_names[i]
    support_range <- range(coreset_data[[var_name]], na.rm = TRUE)
    padding <- 0.1 * diff(support_range)
    if (padding == 0) padding <- 1
    support_range_padded <- support_range + c(-padding, padding)
    
    if(any(!is.finite(support_range_padded))) {
      m <- mean(coreset_data[[var_name]], na.rm=TRUE)
      s <- sd(coreset_data[[var_name]], na.rm=TRUE)
      if (is.na(s) || s == 0) s <- 1
      support_range_padded <- c(m - 10*s, m + 10*s)
    }
    
    model_formula <- as.formula(paste(var_name, "~", deparse1(formula[[2]])))
    tryCatch({
      BoxCox(model_formula, data = coreset_data, 
             support = support_range_padded, order = order, add=c(-0.5, 0.5))
    }, error = function(e) {
      message(paste("Failed to create coreset marginal model for", var_name, ":", e$message))
      return(NULL)
    })
  }, mc.cores = num_cores)
  
  valid_models_idx <- !sapply(margin_models, is.null)
  if (!all(valid_models_idx)) {
    warning("Failed to create some coreset marginal models.")
    margin_models <- margin_models[valid_models_idx]
    if(length(margin_models) < 2) {
      message("  Failed to fit coreset: too few valid marginal models.")
      return(NULL)
    }
  }
  
  # Fit the mmlt model with weights
  # Ensure weights are passed correctly
  args <- c(margin_models, list(formula = formula, data = coreset_data))
                                #, 
                                #weights = weights))
  mctm_model <- tryCatch({
    do.call(mmlt, args)
  }, error = function(e) {
    message("  Failed to fit coreset mmlt model: ", e$message)
    return(NULL)
  })
  
  elapsed_time <- (proc.time() - start_time)[["elapsed"]]
  if (!is.null(mctm_model)) {
    cat(sprintf("  Coreset model successfully fitted in %.2f seconds.\n", elapsed_time))
  } else {
    cat(sprintf("  Failed to fit coreset model after %.2f seconds.\n", elapsed_time))
  }
  return(mctm_model)
}


# ==========================================================================
# Part 5: Evaluation Metrics
# ==========================================================================

#' Evaluate performance metrics of coreset model compared to full model
#' @param coreset_model mmlt model fitted on the coreset
#' @param full_model mmlt model fitted on the full dataset (baseline)
#' @return A list containing metrics: logLik_ratio, param_l2_dist, lambda_l2_dist
evaluate_metrics <- function(coreset_model, full_model) {
  # Initialize metrics to NA
  metrics <- list(
    logLik_ratio = NA_real_,
    param_l2_dist = NA_real_,
    lambda_l2_dist = NA_real_,
    coreset_ll = NA_real_ # Add raw coreset log-likelihood
  )
  
  if (is.null(coreset_model) || is.null(full_model)) {
    warning("One of the models (coreset or full) is NULL. Cannot compute metrics.")
    return(metrics)
  }
  
  # 1. Log Likelihood Ratio
  ll_coreset <- tryCatch(as.numeric(full_model$ll(coef(coreset_model))), error = function(e) NA_real_)
  
  ll_full <- tryCatch(as.numeric(full_model$ll(coef(full_model))), error = function(e) NA_real_)
  metrics$coreset_ll <- ll_coreset # Save raw coreset LL
  if (!is.na(ll_coreset) && !is.na(ll_full) && abs(ll_full) > 1e-9) {
    metrics$logLik_ratio <- ll_coreset / ll_full
  } else {
    warning("Cannot compute LogLik Ratio because full or coreset LL is NA or full LL is zero.")
  }
  
  # 2. Parameter L2 Distance (All Parameters)
  coef_coreset <- tryCatch(coef(coreset_model), error = function(e) NULL)
  coef_full <- tryCatch(coef(full_model), error = function(e) NULL)
  
  if (!is.null(coef_coreset) && !is.null(coef_full)) {
    # Find common parameters by name
    common_params <- intersect(names(coef_coreset), names(coef_full))
    if (length(common_params) > 0) {
      diff_params <- coef_coreset[common_params] - coef_full[common_params]
      metrics$param_l2_dist <- sqrt(sum(diff_params^2))
      
      # 3. Lambda L2 Distance (Only Dependency Parameters)
      # More precise regex pattern for lambda_jl (e.g., y2.y1.(Intercept), y3.y1.(Intercept), y3.y2.(Intercept))
      lambda_param_names <- grep("^y[0-9]+\\.y[0-9]+\\.\\(Intercept\\)$", common_params, value = TRUE)
      
      if (length(lambda_param_names) > 0) {
        diff_lambda <- coef_coreset[lambda_param_names] - coef_full[lambda_param_names]
        metrics$lambda_l2_dist <- sqrt(sum(diff_lambda^2))
      } else {
        # If no lambda parameters matched (e.g., 2D model with no complex dependencies)
        metrics$lambda_l2_dist <- 0 # Or NA? 0 is better if truly no lambda terms
        warning("No matching lambda parameters found between coreset and full model.")
      }
    } else {
      warning("No common parameters found between coreset and full model.")
    }
  } else {
    warning("Failed to extract coefficients from one of the models.")
  }
  
  return(metrics)
}


# ==========================================================================
# Part 6: Main Experiment Pipeline
# ==========================================================================


#' @param data Preprocessed Covertype dataset
#' @param matrices List containing A and A' matrices
#' @param full_model Baseline fitted mmlt model on full data
#' @param min_k Minimum coreset size
#' @param max_k Maximum coreset size
#' @param step_k Step size for increasing coreset size
#' @param num_trials Number of trials per size/method
#' @param k_hull_prop Proportion of hull points in l2_hull method
#' @param seed Main random seed
#' @return A list containing complete experiment results
run_covertype_comparison <- function(data, matrices, full_model,
                                     min_k, max_k, step_k, num_trials,
                                     k_hull_prop = 0.2,
                                     seed = 123) {
  
  n_full <- nrow(data)
  if (n_full == 0) stop("Input data is empty.")
  
  sizes <- seq(min_k, max_k, by = step_k)
  sizes <- sizes[sizes <= n_full]
  if (length(sizes) == 0) stop("Minimum coreset size exceeds dataset size.")
  
  methods <- c("uniform", "l2_only", "l2_hull")
  all_results <- list()
  
  A_list <- matrices$A
  A_prime_list <- matrices$A_prime
  if (is.null(A_list) || is.null(A_prime_list) || length(A_list) == 0 || length(A_prime_list) == 0) {
    stop("A or A' matrices are invalid.")
  }
  
  # --- New: Precompute all convex hull indices ---
  precomputed_hull_indices <- precompute_all_hull_indices(A_prime_list, n_full)
  # ------------------------------------------------
  
  # --- Iterate over Coreset Sizes ---
  for (k in sizes) {
    cat(sprintf("\n--- Running for Coreset Size k = %d ---\n", k))
    size_results <- list()
    
    # --- Iterate over Sampling Methods ---
    for (method in methods) {
      cat(sprintf("  Method: %s\n", method))
      method_trial_results <- list()
      
      # --- Iterate over Trials ---
      trial_seeds <- seed + k * 1000 + match(method, methods) * 100 + (1:num_trials)
      
      trial_data <- mclapply(1:num_trials, function(trial) {
        trial_seed <- trial_seeds[trial]
        cat(sprintf("    Trial %d (seed: %d)...\n", trial, trial_seed))
        
        # 1. Sampling
        sampling_start_time <- proc.time()
        sampler_output <- tryCatch({
          switch(method,
                 "uniform" = sample_uniform(n_full, k, trial_seed),
                 "l2_only" = sample_l2_only(A_list, k, trial_seed),
                 # Pass precomputed hull indices to l2_hull
                 "l2_hull" = sample_l2_hull(A_list, A_prime_list, 
                                            k, k_hull_prop, trial_seed,
                                            precomputed_hull_indices = precomputed_hull_indices)
          )
        }, error = function(e) {
          message(paste("    Sampling failed:", e$message))
          return(NULL)
        })
        sampling_time <- (proc.time() - sampling_start_time)[["elapsed"]]
        
        if (is.null(sampler_output)) {
          cat("      Sampling failed, skipping trial.\n")
          return(NULL)
        }
        cat(sprintf("      Sampling completed in %.2f seconds. Unique indices: %d\n",
                    sampling_time, length(unique(sampler_output$indices))))
        
        # 2. Fitting
        fitting_start_time <- proc.time()
        coreset_model <- fit_coreset_mctm(data, 
                                          sampler_output$indices, 
                                          sampler_output$weights,
                                          formula = ~ 1, order = 6)
        fitting_time <- (proc.time() - fitting_start_time)[["elapsed"]]
        
        if (is.null(coreset_model)) {
          cat("      Coreset fitting failed, saving partial results.\n")
          return(list(
            indices = sampler_output$indices, weights = sampler_output$weights,
            sampling_time = sampling_time, fitting_time = fitting_time, 
            total_time = sampling_time + fitting_time,
            logLik_ratio = NA, param_l2_dist = NA, 
            lambda_l2_dist = NA, coreset_ll = NA, 
            coreset_coef = NULL
          ))
        }
        cat(sprintf("      Fitting completed in %.2f seconds.\n", fitting_time))
        
        # 3. Evaluation
        metrics <- evaluate_metrics(coreset_model, full_model)
        cat(sprintf("      Evaluation: LL Ratio = %.4f, Param L2 = %.4f, Lambda L2 = %.4f\n",
                    metrics$logLik_ratio, 
                    metrics$param_l2_dist, 
                    metrics$lambda_l2_dist))
        
        return(list(
          indices = sampler_output$indices, weights = sampler_output$weights,
          sampling_time = sampling_time, 
          fitting_time = fitting_time, 
          total_time = sampling_time + fitting_time,
          logLik_ratio = metrics$logLik_ratio, 
          param_l2_dist = metrics$param_l2_dist, 
          lambda_l2_dist = metrics$lambda_l2_dist,
          coreset_ll = metrics$coreset_ll,
          coreset_coef = tryCatch(coef(coreset_model), error = function(e) NULL)
        ))
      }, mc.cores = num_cores)
      
      size_results[[method]] <- Filter(Negate(is.null), trial_data)
    }
    all_results[[paste0("k_", k)]] <- size_results
  }
  
  return(all_results)
}



# ==========================================================================
# Part 7: Analysis and Visualization of Results
# ==========================================================================

#' Analyze coreset experiment results across multiple trials
analyze_covertype_results <- function(results) {
  analysis <- list() # List to store analysis per coreset size
  
  for (size_name in names(results)) {
    if (!grepl("^k_", size_name)) next
    k <- as.numeric(gsub("k_", "", size_name))
    method_analysis <- list()
    
    for (method_name in names(results[[size_name]])) {
      trials <- results[[size_name]][[method_name]]
      valid_trials <- Filter(Negate(is.null), trials)
      num_total_trials <- length(trials)
      num_valid_trials <- length(valid_trials)
      
      if (num_valid_trials == 0) {
        warning(paste("No valid trial results for", size_name, method_name))
        method_analysis[[method_name]] <- list(k = k, method = method_name, num_valid_trials = 0)
        next
      }
      
      safe_metric_extractor <- function(trial_list, metric_name) {
        values <- sapply(trial_list, function(t) {
          val <- t[[metric_name]]
          if (is.null(val) || length(val) == 0 || !is.finite(val)) NA else val
        })
        return(as.numeric(values))
      }
      
      logLik_ratios <- safe_metric_extractor(valid_trials, "logLik_ratio")
      param_l2_dists <- safe_metric_extractor(valid_trials, "param_l2_dist")
      lambda_l2_dists <- safe_metric_extractor(valid_trials, "lambda_l2_dist")
      sampling_times <- safe_metric_extractor(valid_trials, "sampling_time")
      fitting_times <- safe_metric_extractor(valid_trials, "fitting_time")
      total_times <- safe_metric_extractor(valid_trials, "total_time")
      
      safe_summary_stats <- function(values) {
        valid_values <- values[!is.na(values)]
        n_valid <- length(valid_values)
        if (n_valid == 0) {
          return(list(mean = NA_real_, sd = NA_real_, median = NA_real_,
                      min = NA_real_, max = NA_real_, n_valid = 0))
        }
        sd_val <- if (n_valid > 1) sd(valid_values) else 0
        min_sd_threshold <- 1e-9
        computed_sd <- max(sd_val, min_sd_threshold, na.rm = TRUE)
        
        list(
          mean = mean(valid_values),
          sd = computed_sd,
          median = median(valid_values),
          min = min(valid_values),
          max = max(valid_values),
          n_valid = n_valid
        )
      }
      
      method_summary <- list(
        k = k,
        method = method_name,
        num_total_trials = num_total_trials,
        num_valid_trials = num_valid_trials,
        logLik_ratio = safe_summary_stats(logLik_ratios),
        param_l2_dist = safe_summary_stats(param_l2_dists),
        lambda_l2_dist = safe_summary_stats(lambda_l2_dists),
        sampling_time = safe_summary_stats(sampling_times),
        fitting_time = safe_summary_stats(fitting_times),
        total_time = safe_summary_stats(total_times)
      )
      method_analysis[[method_name]] <- method_summary
    }
    analysis[[size_name]] <- method_analysis
  }
  return(analysis)
}

#' Plot coreset performance comparison results
plot_covertype_results <- function(analysis, data_gen_name = "Covertype 10D",
                                   xlim = NULL, ylim_logLik = NULL, ylim_param_l2 = NULL,
                                   ylim_lambda_l2 = NULL, ylim_total_time = NULL) {
  
  metrics_df_list <- list()
  for (size_name in names(analysis)) {
    for (method_name in names(analysis[[size_name]])) {
      res <- analysis[[size_name]][[method_name]]
      if (res$num_valid_trials == 0) next
      
      metrics_df_list[[length(metrics_df_list) + 1]] <- data.frame(
        k = res$k,
        method = res$method,
        logLik_ratio_mean = res$logLik_ratio$mean,
        logLik_ratio_sd = res$logLik_ratio$sd,
        logLik_ratio_median = res$logLik_ratio$median,
        logLik_ratio_min = res$logLik_ratio$min,
        logLik_ratio_max = res$logLik_ratio$max,
        param_l2_dist_mean = res$param_l2_dist$mean,
        param_l2_dist_sd = res$param_l2_dist$sd,
        param_l2_dist_median = res$param_l2_dist$median,
        param_l2_dist_min = res$param_l2_dist$min,
        param_l2_dist_max = res$param_l2_dist$max,
        lambda_l2_dist_mean = res$lambda_l2_dist$mean,
        lambda_l2_dist_sd = res$lambda_l2_dist$sd,
        lambda_l2_dist_median = res$lambda_l2_dist$median,
        lambda_l2_dist_min = res$lambda_l2_dist$min,
        lambda_l2_dist_max = res$lambda_l2_dist$max,
        total_time_mean = res$total_time$mean,
        total_time_sd = res$total_time$sd,
        total_time_median = res$total_time$median,
        total_time_min = res$total_time$min,
        total_time_max = res$total_time$max,
        stringsAsFactors = FALSE
      )
    }
  }
  
  if (length(metrics_df_list) == 0) {
    warning("No valid summary data found for plotting.")
    return(list())
  }
  metrics_df <- do.call(rbind, metrics_df_list)
  metrics_df <- subset(metrics_df, method %in% c("uniform", "l2_hull"))
  plots <- list()
  base_theme <- theme_minimal(base_size = 11) +
    theme(
      legend.title = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      axis.title = element_text(size = 12)
    )
  
  add_zoom_text <- function(p, x_zoom = FALSE, y_zoom = FALSE, xlim_plot = NULL, ylim_plot = NULL, df_plot) {
    if (x_zoom || y_zoom) {
      zoom_text <- ifelse(x_zoom && y_zoom, "Zoom X-Y", ifelse(x_zoom, "Zoom X", "Zoom Y"))
      x_pos <- ifelse(!is.null(xlim_plot), xlim_plot[1], min(df_plot$k, na.rm = TRUE))
      y_col_name <- rlang::quo_text(p$mapping$y)
      y_min_val <- if (!is.null(ylim_plot)) ylim_plot[1] else min(df_plot[[y_col_name]], na.rm = TRUE)
      
      p <- p + annotate("text", x = x_pos, y = y_min_val, label = zoom_text,
                        hjust = 0.05, vjust = 1.5, color = "darkred", size = 3, fontface = "italic")
    }
    return(p)
  }
  
  x_zoom_active <- !is.null(xlim)
  
  # --- Plot 1: LogLik Ratio ---
  y_zoom_ll <- !is.null(ylim_logLik)
  plots$logLik <- ggplot(metrics_df, aes(x = k, y = logLik_ratio_mean, color = method, fill = method, group = method)) +
    geom_ribbon(aes(ymin = logLik_ratio_mean - logLik_ratio_sd, ymax = logLik_ratio_mean + logLik_ratio_sd), alpha = 0.25, linetype = 0) +
    geom_line(linewidth = 0.8) + geom_point(size = 1.5) +
    geom_hline(yintercept = 1.0, linetype = "dashed", color = "black", alpha = 0.7) +
    labs(title = "Log-Likelihood Ratio Comparison", subtitle = data_gen_name,
         x = "Coreset Size (k)", y = "LogLik Ratio (Coreset / Full)") +
    scale_color_brewer(palette = "Set1") + scale_fill_brewer(palette = "Set1") +
    coord_cartesian(xlim = xlim, ylim = ylim_logLik, expand = TRUE) + base_theme
  plots$logLik <- add_zoom_text(plots$logLik, x_zoom_active, y_zoom_ll, xlim, ylim_logLik, metrics_df)
  
  # --- Plot 2: Parameter L2 Distance ---
  y_zoom_pL2 <- !is.null(ylim_param_l2)
  plots$param_l2 <- ggplot(metrics_df, aes(x = k, y = param_l2_dist_mean, color = method, fill = method, group = method)) +
    geom_ribbon(aes(ymin = pmax(0, param_l2_dist_mean - param_l2_dist_sd),
                    ymax = param_l2_dist_mean + param_l2_dist_sd), alpha = 0.25, linetype = 0) +
    geom_line(linewidth = 0.8) + geom_point(size = 1.5) +
    labs(title = "Parameter L2 Distance Comparison", subtitle = data_gen_name,
         x = "Coreset Size (k)", y = "L2 Distance (Coreset vs Full)") +
    scale_y_log10() +
    scale_color_brewer(palette = "Set1") + scale_fill_brewer(palette = "Set1") +
    coord_cartesian(xlim = xlim, ylim = ylim_param_l2, expand = TRUE) + base_theme
  plots$param_l2 <- add_zoom_text(plots$param_l2, x_zoom_active, y_zoom_pL2, xlim, ylim_param_l2, metrics_df)
  
  # --- Plot 3: Lambda L2 Distance ---
  y_zoom_lL2 <- !is.null(ylim_lambda_l2)
  plots$lambda_l2 <- ggplot(metrics_df, aes(x = k, y = lambda_l2_dist_mean, color = method, fill = method, group = method)) +
    geom_ribbon(aes(ymin = pmax(0, lambda_l2_dist_mean - lambda_l2_dist_sd),
                    ymax = lambda_l2_dist_mean + lambda_l2_dist_sd), alpha = 0.25, linetype = 0) +
    geom_line(linewidth = 0.8) + geom_point(size = 1.5) +
    labs(title = "Lambda L2 Distance Comparison", subtitle = data_gen_name,
         x = "Coreset Size (k)", y = "L2 Distance (λ)") +
    scale_y_log10() +
    scale_color_brewer(palette = "Set1") + scale_fill_brewer(palette = "Set1") +
    coord_cartesian(xlim = xlim, ylim = ylim_lambda_l2, expand = TRUE) + base_theme
  plots$lambda_l2 <- add_zoom_text(plots$lambda_l2, x_zoom_active, y_zoom_lL2, xlim, ylim_lambda_l2, metrics_df)
  
  # --- Plot 4: Total Time ---
  y_zoom_time <- !is.null(ylim_total_time)
  plots$total_time <- ggplot(metrics_df, aes(x = k, y = total_time_mean, color = method, fill = method, group = method)) +
    geom_ribbon(aes(ymin = pmax(0, total_time_mean - total_time_sd),
                    ymax = total_time_mean + total_time_sd), alpha = 0.25, linetype = 0) +
    geom_line(linewidth = 0.8) + geom_point(size = 1.5) +
    labs(title = "Total Computation Time Comparison", subtitle = data_gen_name,
         x = "Coreset Size (k)", y = "Total Time (seconds)") +
    scale_color_brewer(palette = "Set1") + scale_fill_brewer(palette = "Set1") +
    coord_cartesian(xlim = xlim, ylim = ylim_total_time, expand = TRUE) + base_theme
  plots$total_time <- add_zoom_text(plots$total_time, x_zoom_active, y_zoom_time, xlim, ylim_total_time, metrics_df)
  
  return(plots)
}


#' Generate a Markdown comparison table for a selected coreset size
#'
#' @param analysis Output from analyze_covertype_results()
#' @param k_cutoff The coreset size (or closest available) for which to display the table
#' @return Invisibly returns the summary data.frame
generate_covertype_table <- function(analysis, k_cutoff) {
  available_ks <- sort(as.numeric(gsub("k_", "", names(analysis))))
  if (length(available_ks) == 0) {
    warning("No analysis results found.")
    return(invisible(NULL))
  }
  closest_k <- available_ks[which.min(abs(available_ks - k_cutoff))]
  size_key <- paste0("k_", closest_k)
  
  if (!size_key %in% names(analysis)) {
    warning(sprintf("No results for coreset size k = %d (closest to requested %d)", closest_k, k_cutoff))
    return(invisible(NULL))
  }
  
  cat(sprintf("\n--- Performance Table for Coreset Size k = %d ---\n", closest_k))
  
  table_data_list <- list()
  methods_in_size <- names(analysis[[size_key]])
  if (length(methods_in_size) == 0) {
    warning(sprintf("No methods found for coreset size k = %d", closest_k))
    return(invisible(NULL))
  }
  
  # Use 'uniform' as baseline for relative improvements
  baseline_data <- analysis[[size_key]][["uniform"]]
  if (is.null(baseline_data)) {
    warning("Method 'uniform' not found as baseline. Relative improvements won't be computed.")
    baseline_l2 <- NA; baseline_lambda <- NA; baseline_loglik_mean <- NA
  } else {
    baseline_l2 <- baseline_data$param_l2_dist$mean
    baseline_lambda <- baseline_data$lambda_l2_dist$mean
    baseline_loglik_mean <- baseline_data$logLik_ratio$mean
  }
  
  format_sd_safe <- function(mean_val, sd_val, digits = 3) {
    if (is.null(mean_val) || is.na(mean_val) || is.null(sd_val) || is.na(sd_val)) return("NA")
    sprintf(paste0("%.", digits, "f ± %.", digits, "f"), mean_val, sd_val)
  }
  
  calculate_improvement <- function(value_mean, baseline_mean, lower_is_better = TRUE) {
    if (is.na(baseline_mean) || is.na(value_mean) || abs(baseline_mean) < 1e-9) return("NA")
    if (metric_name == "logLik_ratio") {
      baseline_dist_from_1 <- abs(baseline_mean - 1)
      value_dist_from_1 <- abs(value_mean - 1)
      if (abs(baseline_dist_from_1) < 1e-9) return("0.0%")
      improvement <- (baseline_dist_from_1 - value_dist_from_1) / baseline_dist_from_1
    } else if (lower_is_better) {
      improvement <- (baseline_mean - value_mean) / abs(baseline_mean)
    } else {
      improvement <- (value_mean - baseline_mean) / abs(baseline_mean)
    }
    return(sprintf("%.1f%%", max(0, improvement * 100)))
  }
  
  for (method in methods_in_size) {
    res <- analysis[[size_key]][[method]]
    if (res$num_valid_trials == 0) next
    
    metric_name <- "param_l2_dist";   rel_impr_l2 <- calculate_improvement(res$param_l2_dist$mean, baseline_l2, TRUE)
    metric_name <- "lambda_l2_dist";  rel_impr_lambda <- calculate_improvement(res$lambda_l2_dist$mean, baseline_lambda, TRUE)
    metric_name <- "logLik_ratio";    rel_impr_loglik <- calculate_improvement(res$logLik_ratio$mean, baseline_loglik_mean, FALSE)
    
    table_data_list[[method]] <- data.frame(
      Method       = method,
      ParamL2      = format_sd_safe(res$param_l2_dist$mean, res$param_l2_dist$sd),
      LambdaL2     = format_sd_safe(res$lambda_l2_dist$mean, res$lambda_l2_dist$sd),
      LogLikRatio  = format_sd_safe(res$logLik_ratio$mean, res$logLik_ratio$sd),
      Impr_L2      = if (method == "uniform") "baseline" else rel_impr_l2,
      Impr_Lambda  = if (method == "uniform") "baseline" else rel_impr_lambda,
      Impr_LL      = if (method == "uniform") "baseline" else rel_impr_loglik,
      TotalTime    = format_sd_safe(res$total_time$mean, res$total_time$sd, 2),
      ValidTrials  = sprintf("%d/%d", res$num_valid_trials, res$num_total_trials),
      stringsAsFactors = FALSE
    )
  }
  
  if (length(table_data_list) == 0) {
    warning(sprintf("No valid method results for table at k = %d", closest_k))
    return(invisible(NULL))
  }
  
  final_table <- do.call(rbind, table_data_list)
  
  preferred_order <- c("l2_hull", "l2_only", "uniform")
  ordered_methods <- c(intersect(preferred_order, final_table$Method),
                       setdiff(final_table$Method, preferred_order))
  final_table <- final_table[match(ordered_methods, final_table$Method), ]
  
  cat("| Method   | Param L2 Dist   | Lambda L2 Dist  | LogLik Ratio    | Impr. L2 | Impr. Lambda | Impr. LL | Total Time (s) | Valid Trials     |\n")
  cat("|----------|------------------|------------------|------------------|----------|----------------|----------|------------------|------------------|\n")
  for (i in 1:nrow(final_table)) {
    cat(sprintf("| %-8s | %-16s | %-16s | %-16s | %8s | %14s | %8s | %16s | %16s |\n",
                final_table$Method[i], final_table$ParamL2[i], final_table$LambdaL2[i],
                final_table$LogLikRatio[i], final_table$Impr_L2[i], final_table$Impr_Lambda[i],
                final_table$Impr_LL[i], final_table$TotalTime[i], final_table$ValidTrials[i]))
  }
  cat("\n")
  invisible(final_table)
}

# ==========================================================================
# Part 8: Main Script to Run Covertype Coreset Experiments
# ==========================================================================

# --- Main Experiment Parameters ---
DATA_FILE_PATH <- "covtype.data.gz"  # Make sure this file exists or change the path
N_SAMPLE_DATA <- 1000              # Sample size from Covertype
MIN_K <- 10                         # Minimum coreset size
MAX_K <- 500                        # Maximum coreset size (increase if resources allow)
STEP_K <- 10                         # Step size for coreset size
NUM_TRIALS <- 5                      # Number of trials per size/method (increase for stability)
SEED <- 42                           # Global random seed
FORMULA <- ~ 1                       # Unconditional model
ORDER <- 7                          # Bernstein polynomial order
K_HULL_PROP <- 0.05                  # Proportion allocated to hull points

# --- Run the Full Experiment Pipeline ---
start_experiment_time <- Sys.time()
cat("Starting Coreset MCTM experiment pipeline for Covertype...\n")
  
  # 1. Load and prepare the data
covertype_data_sample <- load_prepare_covertype(n_sample = N_SAMPLE_DATA, 
                                                seed = SEED, 
                                                file_path = DATA_FILE_PATH)
covertype_data_sample <- covertype_data_sample

DIMENSIONS <- ncol(covertype_data_sample)  # Should be 10
dim(covertype_data_sample)
  # 2. Fit full model (baseline)
full_mctm_model <- fit_full_mctm(covertype_data_sample, formula = FORMULA,
                                   order = ORDER)

# 3. Extract A and A' matrices
length(full_mctm_model$coef)

mctm_matrices <- extract_mctm_matrices(covertype_data_sample, 
                                       formula = FORMULA, order = ORDER)
  
  
dim(mctm_matrices)
  # 4. Run the coreset com)parison experiments
comparison_results <- run_covertype_comparison(
  data = covertype_data_sample,
  matrices = mctm_matrices,
  full_model = full_mctm_model,
  min_k = MIN_K,
  max_k = MAX_K,
  step_k = STEP_K,
  num_trials = NUM_TRIALS,
  k_hull_prop = K_HULL_PROP,
  seed = SEED
)
# 5. Analyze the results
cat("\nAnalyzing experiment results...\n")
analysis_summary <- analyze_covertype_results(comparison_results)
# 6. Visualize the results
cat("Creating visualizations...\n")
plots <- plot_covertype_results(analysis_summary,
                                data_gen_name = "Covertype 10D (Unconditional)",
                                ylim_logLik = c(-1,-100),
                                ylim_total_time = c(0,50))
# Display main plots
print(plots$logLik)
print(plots$param_l2)
print(plots$lambda_l2)
print(plots$total_time)

# 7. Generate performance tables for specific k
cat("\nGenerating performance tables...\n")
performance_table_50 <- generate_covertype_table(analysis_summary, k_cutoff = 50)
performance_table_100 <- generate_covertype_table(analysis_summary, k_cutoff = 100)
performance_table_200 <- generate_covertype_table(analysis_summary, k_cutoff = 200)
performance_table_500 <- generate_covertype_table(analysis_summary, k_cutoff = 500)
  


