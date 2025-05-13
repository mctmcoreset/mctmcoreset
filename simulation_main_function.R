# Simulation_main_function.R
# Function to generate a performance comparison table at specified coreset size
generate_performance_table <- function(analysis_results, cutoff_size, data_gen_name = "Default") {
  # Find the closest size to the cutoff_size
  sizes <- sapply(names(analysis_results), function(s) as.numeric(gsub("size_", "", s)))
  closest_size_idx <- which.min(abs(sizes - cutoff_size))
  closest_size <- sizes[closest_size_idx]
  size_key <- paste0("size_", closest_size)
  
  # Check if the size exists in analysis_results
  if(!size_key %in% names(analysis_results)) {
    stop(paste("Size", cutoff_size, "not found. Closest available size is", closest_size, "but not in results."))
  }
  
  cat(paste("Using coreset size =", closest_size, "(closest to specified", cutoff_size, ")\n"))
  
  # Get methods
  methods <- names(analysis_results[[size_key]])
  
  # Ensure l2_hull, l2_only, and uniform are present
  required_methods <- c("l2_hull", "l2_only", "uniform")
  missing_methods <- required_methods[!required_methods %in% methods]
  if(length(missing_methods) > 0) {
    warning(paste("Some methods are missing:", paste(missing_methods, collapse=", ")))
  }
  
  # Prepare table data
  table_data <- data.frame(
    Method = character(),
    ParamL2Distance = character(),
    LambdaError = character(),
    LogLikelihoodRatio = character(),
    RelativeImprovement = character(),
    TotalTime = character(),
    stringsAsFactors = FALSE
  )
  
  # Get baseline values (uniform method)
  if("uniform" %in% methods) {
    uniform_data <- analysis_results[[size_key]][["uniform"]]
    baseline_l2 <- uniform_data$param_l2_mean
    baseline_lambda <- uniform_data$lambda_error_mean
    baseline_loglik <- uniform_data$logLik_mean
    baseline_time <- uniform_data$mean_total_time
  } else {
    warning("Uniform method not found. Relative improvements will not be calculated.")
    baseline_l2 <- NA
    baseline_lambda <- NA
    baseline_loglik <- NA
    baseline_time <- NA
  }
  
  # Format a value with standard deviation
  format_with_sd <- function(mean_val, sd_val, digits = 2) {
    return(sprintf("%.2f ± %.2f", round(mean_val, digits), round(sd_val, digits)))
  }
  
  # Calculate relative improvement compared to uniform
  calculate_improvement <- function(value, baseline, lower_is_better = TRUE) {
    if(is.na(baseline) || is.na(value)) return("N/A")
    
    if(lower_is_better) {
      improvement <- (baseline - value) / baseline * 100
    } else {
      # For metrics where higher is better (e.g., loglik closer to 1.0)
      # Lower absolute difference from 1.0 is better
      value_diff <- abs(value - 1.0)
      baseline_diff <- abs(baseline - 1.0)
      improvement <- (baseline_diff - value_diff) / baseline_diff * 100
      # If baseline is exactly 1.0, this would be division by zero
      if(baseline_diff == 0) return("N/A")
    }
    
    return(sprintf("%.1f%%", improvement))
  }
  
  # Process each method
  for(method in methods) {
    method_data <- analysis_results[[size_key]][[method]]
    
    # Format param L2 distance
    param_l2 <- format_with_sd(method_data$param_l2_mean, method_data$param_l2_sd)
    
    # Format lambda error
    lambda_error <- format_with_sd(method_data$lambda_error_mean, method_data$lambda_error_sd)
    
    # Format log-likelihood ratio
    loglik <- format_with_sd(method_data$logLik_mean, method_data$logLik_sd)
    
    # Calculate relative improvements
    if(method != "uniform") {
      rel_improvement_l2 <- calculate_improvement(method_data$param_l2_mean, baseline_l2)
      rel_improvement_lambda <- calculate_improvement(method_data$lambda_error_mean, baseline_lambda)
      rel_improvement_loglik <- calculate_improvement(method_data$logLik_mean, baseline_loglik, lower_is_better = FALSE)
      
      # Average improvement across metrics
      improvements <- c(
        as.numeric(gsub("%", "", rel_improvement_l2)),
        as.numeric(gsub("%", "", rel_improvement_lambda)),
        as.numeric(gsub("%", "", rel_improvement_loglik))
      )
      improvements <- improvements[!is.na(improvements)]
      
      if(length(improvements) > 0) {
        avg_improvement <- sprintf("%.1f%%", mean(improvements, na.rm = TRUE))
      } else {
        avg_improvement <- "N/A"
      }
    } else {
      # Uniform is the baseline, so relative improvement is N/A
      avg_improvement <- "baseline"
    }
    
    # Format total time
    total_time <- format_with_sd(method_data$mean_total_time, method_data$sd_total_time)
    
    # Add to table data
    table_data <- rbind(table_data, data.frame(
      Method = method,
      ParamL2Distance = param_l2,
      LambdaError = lambda_error,
      LogLikelihoodRatio = loglik,
      RelativeImprovement = avg_improvement,
      TotalTime = total_time,
      stringsAsFactors = FALSE
    ))
  }
  
  # Reorder methods to put l2_hull, l2_only, and uniform first
  preferred_order <- c("l2_hull", "l2_only", "uniform")
  other_methods <- setdiff(table_data$Method, preferred_order)
  ordered_methods <- c(preferred_order[preferred_order %in% table_data$Method], other_methods)
  table_data <- table_data[match(ordered_methods, table_data$Method),]
  
  # Print formatted table
  cat("\n")
  cat(paste("Performance Comparison at Coreset Size =", closest_size, "\n"))
  cat(paste("Data Generation Process:", data_gen_name, "\n"))
  cat("\n")
  
  # Print as markdown table
  cat("| Method | Param L2 Distance | Lambda Error | Log-Likelihood Ratio | Relative Improvement | Total Time (s) |\n")
  cat("|--------|-------------------|--------------|----------------------|----------------------|----------------|\n")
  
  for(i in 1:nrow(table_data)) {
    cat(sprintf("| %s | %s | %s | %s | %s | %s |\n",
                table_data$Method[i],
                table_data$ParamL2Distance[i],
                table_data$LambdaError[i],
                table_data$LogLikelihoodRatio[i],
                table_data$RelativeImprovement[i],
                table_data$TotalTime[i]))
  }
  cat("\n")
  
  # Return the table data invisibly
  return(invisible(table_data))
}

# Example usage:
# # For coreset size = 30
# table_30 <- generate_performance_table(analysis, 30, "Bivariate Normal")
# 
# # For coreset size = 100
# table_100 <- generate_performance_table(analysis, 100, "Bivariate Normal")
# Modified wrapper function to run multiple simulations with different datasets
run_multiple_simulations <- function(data_gen_fn, n_samples, n_sims = 5,
                                     min_size = 10, max_size = NULL, step = 10,
                                     num_trials = 5) {  # Added num_trials parameter
  # Initialize results list
  all_results <- list()
  
  # Run simulations
  for(i in 1:n_sims) {
    cat(sprintf("\nSimulation %d/%d\n", i, n_sims))
    
    # Create a list to store trial results for this simulation
    sim_trials <- list()
    
    # Run multiple trials with different random seeds
    for(j in 1:num_trials) {
      cat(sprintf("  Trial %d/%d\n", j, num_trials))
      
      # Generate absolutely unique random seed for each trial
      current_time <- as.numeric(Sys.time())
      unique_seed <- round(i*10000 + j*100 + (current_time %% 100))
      set.seed(unique_seed)
      cat(sprintf("    Using seed: %d\n", unique_seed))
      
      # Generate dataset
      data <- data_gen_fn(n = n_samples)
      
      # Print basic dataset statistics to verify variation
      cat(sprintf("    Data summary - Mean y1: %.4f, SD y1: %.4f, Mean y2: %.4f, SD y2: %.4f\n", 
                  mean(data[,1]), sd(data[,1]), mean(data[,2]), sd(data[,2])))
      
      # Fit full model for this dataset
      x <- data[,1]
      y <- data[,2]
      m_x <- tram::BoxCox(x ~ 1)
      m_y <- tram::BoxCox(y ~ 1)
      m_xy <- mmlt(m_x, m_y, formula = ~ 1)
      
      # Store the model and run comparison for this trial
      trial_result <- run_mctm_comparison(
        data = data,
        min_size = min_size,
        max_size = max_size,
        step = step,
        num_trials = 1  # Only need 1 iteration per dataset
      )
      
      # Add the full model to results
      trial_result$full_model <- m_xy
      
      # Also store the raw data for troubleshooting
      trial_result$raw_data <- data
      
      # Store trial results
      sim_trials[[j]] <- trial_result
    }
    
    # Store all trials for this simulation
    all_results[[i]] <- sim_trials
  }
  
  return(all_results)
}

# Modified analysis function to handle multiple trials per simulation
analyze_multiple_simulations <- function(all_results) {
  n_sims <- length(all_results)
  
  # Ensure at least one simulation exists
  if(n_sims == 0) {
    stop("No simulation results found")
  }
  
  # Ensure first simulation exists and has at least one trial
  if(is.null(all_results[[1]]) || length(all_results[[1]]) == 0) {
    stop("First simulation result is empty")
  }
  
  n_trials <- length(all_results[[1]])
  
  # Get size names from first simulation, first trial
  first_trial <- all_results[[1]][[1]]
  exclude_names <- c("baseline", "hull_indices", "matrices", "data", "full_model", "raw_data")
  sizes <- names(first_trial)[!names(first_trial) %in% exclude_names]
  
  if(length(sizes) == 0) {
    stop("No size information found in results")
  }
  
  analysis <- list()
  
  for(size in sizes) {
    # Extract numeric size
    size_num <- as.numeric(gsub("size_", "", size))
    if(is.na(size_num)) {
      warning(paste("Failed to convert", size, "to numeric. Skipping."))
      next
    }
    
    # Get methods list
    if(is.null(first_trial[[size]])) {
      warning(paste("Size", size, "not found in first trial. Skipping."))
      next
    }
    
    methods <- names(first_trial[[size]])
    if(length(methods) == 0) {
      warning(paste("No methods found in size", size, ". Skipping."))
      next
    }
    
    method_analysis <- list()
    for(method in methods) {
      # Lists to store all values across sims and trials
      all_logLik_values <- c()
      all_param_l2_values <- c()
      all_lambda_error_values <- c()
      all_sizes <- c()
      all_sampling_times <- c()
      all_optimization_times <- c()
      all_total_times <- c()
      
      # Process all simulations and trials
      for(sim in 1:n_sims) {
        # Check if this simulation exists
        if(is.null(all_results[[sim]])) {
          cat("Skipping missing simulation", sim, "\n")
          next
        }
        
        for(t in 1:n_trials) {
          # Check if this trial exists
          if(is.null(all_results[[sim]][[t]])) {
            cat("Skipping missing trial - sim", sim, "trial", t, "\n")
            next
          }
          
          trial_results <- all_results[[sim]][[t]]
          
          # Check if full model exists
          if(is.null(trial_results$full_model)) {
            cat("Skipping missing full model - sim", sim, "trial", t, "\n")
            next
          }
          
          # Check if this size and method exist
          if(is.null(trial_results[[size]]) || is.null(trial_results[[size]][[method]]) || 
             length(trial_results[[size]][[method]]) < 1) {
            cat("Skipping missing method results - sim", sim, "trial", t, "size", size, "method", method, "\n")
            next
          }
          
          # Extract model and trial result
          m_xy <- trial_results$full_model
          trial_result <- trial_results[[size]][[method]][[1]]
          
          # Calculate metrics
          tryCatch({
            full_loglik <- m_xy$ll(m_xy$par)
            coreset_loglik <- m_xy$ll(trial_result$coefficients)
            # Use absolute value to ensure positive values
            likelihood_ratio <- abs(coreset_loglik / full_loglik)
            all_logLik_values <- c(all_logLik_values, likelihood_ratio)
            
            # Parameter error calculation
            true_params <- m_xy$par
            coreset_params <- trial_result$coefficients
            param_diff <- coreset_params - true_params
            l2_dist <- sqrt(sum(param_diff^2))
            all_param_l2_values <- c(all_param_l2_values, l2_dist)
            
            # Lambda error calculation
            lambda_error <- abs(param_diff[length(param_diff)])
            all_lambda_error_values <- c(all_lambda_error_values, lambda_error)
            
            # Other metrics
            all_sizes <- c(all_sizes, trial_result$size)
            all_sampling_times <- c(all_sampling_times, trial_result$sampling_time)
            all_optimization_times <- c(all_optimization_times, trial_result$optimization_time)
            all_total_times <- c(all_total_times, trial_result$total_time)
          }, error = function(e) {
            cat("Error calculating metrics:", as.character(e), "- sim", sim, "trial", t, "\n")
          })
        }
      }
      
      # Check if we have enough data to analyze
      if(length(all_logLik_values) == 0) {
        warning(paste("Size", size, "method", method, "has no valid data"))
        next
      }
      
      # Print debug info to verify data variation
      cat("Size:", size_num, "Method:", method, "\n")
      cat("LogLik values range:", range(all_logLik_values), "SD:", sd(all_logLik_values), "\n")
      cat("Param L2 values range:", range(all_param_l2_values), "SD:", sd(all_param_l2_values), "\n")
      cat("Lambda error values range:", range(all_lambda_error_values), "SD:", sd(all_lambda_error_values), "\n")
      
      # Ensure non-zero standard deviation
      min_sd <- 0.001
      
      # Compute statistics with error checking
      method_analysis[[method]] <- list(
        size = size_num,
        # Log-likelihood statistics
        logLik_mean = mean(all_logLik_values),
        logLik_sd = max(min_sd, sd(all_logLik_values)), # Ensure minimum visible SD
        logLik_median = median(all_logLik_values),
        logLik_min = min(all_logLik_values),
        logLik_max = max(all_logLik_values),
        
        # Parameter L2 distance statistics
        param_l2_mean = mean(all_param_l2_values),
        param_l2_sd = max(min_sd, sd(all_param_l2_values)), # Ensure minimum visible SD
        param_l2_median = median(all_param_l2_values),
        param_l2_min = min(all_param_l2_values),
        param_l2_max = max(all_param_l2_values),
        
        # Lambda error statistics
        lambda_error_mean = mean(all_lambda_error_values),
        lambda_error_sd = max(min_sd, sd(all_lambda_error_values)), # Ensure minimum visible SD
        lambda_error_median = median(all_lambda_error_values),
        lambda_error_min = min(all_lambda_error_values),
        lambda_error_max = max(all_lambda_error_values),
        
        # Time statistics
        actual_size = mean(all_sizes),
        mean_sampling_time = mean(all_sampling_times),
        sd_sampling_time = sd(all_sampling_times),
        mean_optimization_time = mean(all_optimization_times),
        sd_optimization_time = sd(all_optimization_times),
        mean_total_time = mean(all_total_times),
        sd_total_time = sd(all_total_times)
      )
    }
    analysis[[size]] <- method_analysis
  }
  
  return(analysis)
}

# Modified plotting function to show mean and standard deviation error bars
# Modified plotting function to show mean and standard deviation error bars
# with zoom capability for both x and y axes
plot_multiple_simulation_results <- function(analysis_results, data_gen_name = "Default", 
                                             xlim = NULL,  # New parameter for x-axis zoom
                                             ylim_logLik = NULL, 
                                             ylim_param_l2 = NULL, 
                                             ylim_lambda = NULL,
                                             ylim_sampling_time = NULL,
                                             ylim_optimization_time = NULL,
                                             ylim_total_time = NULL) {
  # First, verify the structure of analysis_results
  if(length(analysis_results) == 0) {
    stop("Analysis results are empty")
  }
  
  # Extract size names
  sizes <- names(analysis_results)
  if(length(sizes) == 0) {
    stop("No size data found in analysis results")
  }
  
  # Check if first element exists and get methods
  if(is.null(analysis_results[[1]])) {
    stop("First element of analysis results is NULL")
  }
  
  methods <- names(analysis_results[[1]])
  if(length(methods) == 0) {
    stop("No methods found in analysis results")
  }
  
  # Create an empty data frame with the correct column types
  metrics_df <- data.frame(
    size = numeric(),
    method = character(),
    logLik_mean = numeric(),
    logLik_sd = numeric(),
    logLik_median = numeric(),
    logLik_min = numeric(),
    logLik_max = numeric(),
    param_l2_mean = numeric(),
    param_l2_sd = numeric(),
    param_l2_median = numeric(),
    param_l2_min = numeric(),
    param_l2_max = numeric(),
    lambda_error_mean = numeric(),
    lambda_error_sd = numeric(),
    lambda_error_median = numeric(),
    lambda_error_min = numeric(),
    lambda_error_max = numeric(),
    sampling_time = numeric(),
    sampling_time_sd = numeric(),
    optimization_time = numeric(),
    optimization_time_sd = numeric(),
    total_time = numeric(),
    total_time_sd = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Safely populate the data frame
  for(size in sizes) {
    # Convert size_X to numeric
    size_num <- as.numeric(gsub("size_", "", size))
    
    # Skip if conversion failed
    if(is.na(size_num)) {
      warning(paste("Failed to convert", size, "to numeric. Skipping."))
      next
    }
    
    for(method in methods) {
      # Check if this method exists for this size
      if(is.null(analysis_results[[size]][[method]])) {
        warning(paste("Method", method, "not found for size", size, ". Skipping."))
        next
      }
      
      res <- analysis_results[[size]][[method]]
      
      # Build a single row with all metrics
      row_data <- data.frame(
        size = size_num,
        method = method,
        stringsAsFactors = FALSE
      )
      
      # Safely add numeric metrics
      safe_add <- function(df, name, value) {
        if(!is.null(value) && !is.na(value) && is.numeric(value)) {
          df[[name]] <- value
        } else {
          df[[name]] <- NA
          warning(paste("Value for", name, "is not valid. Setting to NA."))
        }
        return(df)
      }
      
      # Safely add all metrics
      row_data <- safe_add(row_data, "logLik_mean", res$logLik_mean)
      row_data <- safe_add(row_data, "logLik_sd", res$logLik_sd)
      row_data <- safe_add(row_data, "logLik_median", res$logLik_median)
      row_data <- safe_add(row_data, "logLik_min", res$logLik_min)
      row_data <- safe_add(row_data, "logLik_max", res$logLik_max)
      
      row_data <- safe_add(row_data, "param_l2_mean", res$param_l2_mean)
      row_data <- safe_add(row_data, "param_l2_sd", res$param_l2_sd)
      row_data <- safe_add(row_data, "param_l2_median", res$param_l2_median)
      row_data <- safe_add(row_data, "param_l2_min", res$param_l2_min)
      row_data <- safe_add(row_data, "param_l2_max", res$param_l2_max)
      
      row_data <- safe_add(row_data, "lambda_error_mean", res$lambda_error_mean)
      row_data <- safe_add(row_data, "lambda_error_sd", res$lambda_error_sd)
      row_data <- safe_add(row_data, "lambda_error_median", res$lambda_error_median)
      row_data <- safe_add(row_data, "lambda_error_min", res$lambda_error_min)
      row_data <- safe_add(row_data, "lambda_error_max", res$lambda_error_max)
      
      row_data <- safe_add(row_data, "sampling_time", res$mean_sampling_time)
      row_data <- safe_add(row_data, "sampling_time_sd", res$sd_sampling_time)
      row_data <- safe_add(row_data, "optimization_time", res$mean_optimization_time)
      row_data <- safe_add(row_data, "optimization_time_sd", res$sd_optimization_time)
      row_data <- safe_add(row_data, "total_time", res$mean_total_time)
      row_data <- safe_add(row_data, "total_time_sd", res$sd_total_time)
      
      # Add the row to the data frame
      metrics_df <- rbind(metrics_df, row_data)
    }
  }
  
  # Check if we have data to plot
  if(nrow(metrics_df) == 0) {
    stop("No data to plot after processing")
  }
  
  # Create plots using ggplot2
  plots <- list()
  
  # Common theme
  base_theme <- theme_minimal() +
    theme(legend.title = element_blank(),
          legend.position = "bottom",
          panel.grid.major = element_line(color = "gray90"),
          panel.grid.minor = element_line(color = "gray95"))
  
  # Helper function to add zoom text if either x or y limits are set
  add_zoom_text <- function(p, x_zoom = FALSE, y_zoom = FALSE, x_pos = NULL, y_pos = NULL) {
    if(x_zoom || y_zoom) {
      zoom_text <- ifelse(x_zoom && y_zoom, "X-Y Zoomed view", 
                          ifelse(x_zoom, "X-Zoomed view", "Y-Zoomed view"))
      
      p <- p + annotate("text", 
                        x = if(is.null(x_pos)) min(metrics_df$size) else x_pos, 
                        y = if(is.null(y_pos)) min(metrics_df$param_l2_mean) else y_pos,
                        label = zoom_text, 
                        hjust = 0, vjust = 0, 
                        color = "darkred", fontface = "italic", size = 3)
    }
    return(p)
  }
  
  # Log-likelihood ratio plot with error bars - CONSISTENT FORMATTING
  plots$logLik <- ggplot(metrics_df) +
    geom_ribbon(aes(x = size,
                    ymin = logLik_mean - logLik_sd,
                    ymax = logLik_mean + logLik_sd,
                    fill = method,
                    group = method),
                alpha = 0.3) +
    geom_line(aes(x = size,
                  y = logLik_mean,
                  color = method,
                  group = method), 
              linewidth = 1) +
    geom_point(aes(x = size,
                   y = logLik_mean,
                   color = method),
               size = 2) +
    labs(title = "Log-Likelihood Ratio vs Coreset Size",
         subtitle = paste("Data Generation Process:", data_gen_name),
         x = "Coreset Size",
         y = "Likelihood Ratio") +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    base_theme
  
  # Parameter L2 distance plot with error bars - CONSISTENT FORMATTING
  plots$param_l2 <- ggplot(metrics_df) +
    geom_ribbon(aes(x = size,
                    ymin = pmax(0, param_l2_mean - param_l2_sd),
                    ymax = param_l2_mean + param_l2_sd,
                    fill = method,
                    group = method),
                alpha = 0.3) +
    geom_line(aes(x = size,
                  y = param_l2_mean,
                  color = method,
                  group = method),
              linewidth = 1) +
    geom_point(aes(x = size,
                   y = param_l2_mean,
                   color = method),
               size = 2) +
    labs(title = "Parameter L2 Distance vs Coreset Size",
         subtitle = paste("Data Generation Process:", data_gen_name),
         x = "Coreset Size",
         y = "L2 Distance") +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    base_theme
  
  # Lambda error plot with error bars - CONSISTENT FORMATTING
  plots$lambda <- ggplot(metrics_df) +
    geom_ribbon(aes(x = size,
                    ymin = pmax(0, lambda_error_mean - lambda_error_sd),
                    ymax = lambda_error_mean + lambda_error_sd,
                    fill = method,
                    group = method),
                alpha = 0.3) +
    geom_line(aes(x = size,
                  y = lambda_error_mean,
                  color = method,
                  group = method),
              linewidth = 1) +
    geom_point(aes(x = size,
                   y = lambda_error_mean,
                   color = method),
               size = 2) +
    labs(title = "Lambda Parameter Error vs Coreset Size",
         subtitle = paste("Data Generation Process:", data_gen_name),
         x = "Coreset Size",
         y = "Absolute Error") +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    base_theme
  
  # Time plots with error bars - FIXED THE FORMATTING TO BE CONSISTENT WITH OTHERS
  plots$sampling_time <- ggplot(metrics_df) +
    geom_ribbon(aes(x = size,
                    ymin = pmax(0, sampling_time - sampling_time_sd),
                    ymax = sampling_time + sampling_time_sd,
                    fill = method,
                    group = method),
                alpha = 0.3) +
    geom_line(aes(x = size,
                  y = sampling_time,
                  color = method,
                  group = method),
              linewidth = 1) +
    geom_point(aes(x = size,
                   y = sampling_time,
                   color = method),
               size = 2) +
    labs(title = "Sampling Time vs Coreset Size",
         subtitle = paste("Data Generation Process:", data_gen_name),
         x = "Coreset Size",
         y = "Time (seconds)") +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    base_theme
  
  plots$optimization_time <- ggplot(metrics_df) +
    geom_ribbon(aes(x = size,
                    ymin = pmax(0, optimization_time - optimization_time_sd),
                    ymax = optimization_time + optimization_time_sd,
                    fill = method,
                    group = method),
                alpha = 0.3) +
    geom_line(aes(x = size,
                  y = optimization_time,
                  color = method,
                  group = method),
              linewidth = 1) +
    geom_point(aes(x = size,
                   y = optimization_time,
                   color = method),
               size = 2) +
    labs(title = "Optimization Time vs Coreset Size",
         subtitle = paste("Data Generation Process:", data_gen_name),
         x = "Coreset Size",
         y = "Time (seconds)") +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    base_theme
  
  plots$total_time <- ggplot(metrics_df) +
    geom_ribbon(aes(x = size,
                    ymin = pmax(0, total_time - total_time_sd),
                    ymax = total_time + total_time_sd,
                    fill = method,
                    group = method),
                alpha = 0.3) +
    geom_line(aes(x = size,
                  y = total_time,
                  color = method,
                  group = method),
              linewidth = 1) +
    geom_point(aes(x = size,
                   y = total_time,
                   color = method),
               size = 2) +
    labs(title = "Total Time vs Coreset Size",
         subtitle = paste("Data Generation Process:", data_gen_name),
         x = "Coreset Size",
         y = "Time (seconds)") +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    base_theme
  
  # Apply zoom effects based on provided limits
  
  # X-axis zoom (apply to all plots if specified)
  x_zoom <- !is.null(xlim)
  
  # Define min and max x values for placing zoom indicator text
  x_min <- if(x_zoom) xlim[1] else min(metrics_df$size)
  
  # Apply Y-axis zoom for each plot
  if(!is.null(ylim_logLik)) {
    plots$logLik <- plots$logLik + 
      coord_cartesian(xlim = xlim, ylim = ylim_logLik)
    plots$logLik <- add_zoom_text(plots$logLik, x_zoom, TRUE, x_min, ylim_logLik[1])
  } else if(x_zoom) {
    plots$logLik <- plots$logLik + coord_cartesian(xlim = xlim)
    plots$logLik <- add_zoom_text(plots$logLik, TRUE, FALSE, x_min)
  }
  
  if(!is.null(ylim_param_l2)) {
    plots$param_l2 <- plots$param_l2 + 
      coord_cartesian(xlim = xlim, ylim = ylim_param_l2)
    plots$param_l2 <- add_zoom_text(plots$param_l2, x_zoom, TRUE, x_min, ylim_param_l2[1])
  } else if(x_zoom) {
    plots$param_l2 <- plots$param_l2 + coord_cartesian(xlim = xlim)
    plots$param_l2 <- add_zoom_text(plots$param_l2, TRUE, FALSE, x_min)
  }
  
  if(!is.null(ylim_lambda)) {
    plots$lambda <- plots$lambda + 
      coord_cartesian(xlim = xlim, ylim = ylim_lambda)
    plots$lambda <- add_zoom_text(plots$lambda, x_zoom, TRUE, x_min, ylim_lambda[1])
  } else if(x_zoom) {
    plots$lambda <- plots$lambda + coord_cartesian(xlim = xlim)
    plots$lambda <- add_zoom_text(plots$lambda, TRUE, FALSE, x_min)
  }
  
  if(!is.null(ylim_sampling_time)) {
    plots$sampling_time <- plots$sampling_time + 
      coord_cartesian(xlim = xlim, ylim = ylim_sampling_time)
    plots$sampling_time <- add_zoom_text(plots$sampling_time, x_zoom, TRUE, x_min, ylim_sampling_time[1])
  } else if(x_zoom) {
    plots$sampling_time <- plots$sampling_time + coord_cartesian(xlim = xlim)
    plots$sampling_time <- add_zoom_text(plots$sampling_time, TRUE, FALSE, x_min)
  }
  
  if(!is.null(ylim_optimization_time)) {
    plots$optimization_time <- plots$optimization_time + 
      coord_cartesian(xlim = xlim, ylim = ylim_optimization_time)
    plots$optimization_time <- add_zoom_text(plots$optimization_time, x_zoom, TRUE, x_min, ylim_optimization_time[1])
  } else if(x_zoom) {
    plots$optimization_time <- plots$optimization_time + coord_cartesian(xlim = xlim)
    plots$optimization_time <- add_zoom_text(plots$optimization_time, TRUE, FALSE, x_min)
  }
  
  if(!is.null(ylim_total_time)) {
    plots$total_time <- plots$total_time + 
      coord_cartesian(xlim = xlim, ylim = ylim_total_time)
    plots$total_time <- add_zoom_text(plots$total_time, x_zoom, TRUE, x_min, ylim_total_time[1])
  } else if(x_zoom) {
    plots$total_time <- plots$total_time + coord_cartesian(xlim = xlim)
    plots$total_time <- add_zoom_text(plots$total_time, TRUE, FALSE, x_min)
  }
  
  return(plots)
}

# ==========================================================================
# Implementation of Analysis for Various Data Generation Functions
# ==========================================================================

# --- Example 1: Non-linear Correlation ---
results_nonlinear <- run_multiple_simulations(
  data_gen_fn = generate_nonlinear_correlation,
  n_samples = 1000,
  n_sims = 3,        # Number of outer simulations
  min_size = 10,
  max_size = 500,    # Use smaller max size initially for faster testing
  step = 5,
  num_trials = 10     # Number of trials per simulation - increase for more variation
)

# Analyze and plot
analysis_nonlinear <- analyze_multiple_simulations(results_nonlinear)
plots_nonlinear <- plot_multiple_simulation_results(
  analysis_results = analysis_nonlinear,
  data_gen_name = "Non-linear Correlation",
  ylim_logLik = c(1.0, 1.5),          # Custom range for logLik plot
  ylim_param_l2 = c(0.5, 5.0),        # Custom range for param_l2 plot
  ylim_lambda = c(0, 1),              # Custom range for lambda plot
  ylim_total_time = c(0, 0.5),
  xlim = c(0, 200)                    # Custom range for total_time plot
)

# Display plots
plots_nonlinear$logLik
plots_nonlinear$param_l2
plots_nonlinear$lambda
plots_nonlinear$total_time

# Comparison table for coreset size = 30
table_30_nonlinear <- generate_performance_table(analysis_nonlinear, 30, "Non-linear Correlation")

# Comparison table for coreset size = 100
table_100_nonlinear <- generate_performance_table(analysis_nonlinear, 100, "Non-linear Correlation")

# --- Example 2: Bivariate Normal ---
results_normal <- run_multiple_simulations(
  data_gen_fn = generate_bivariate_normal,
  n_samples = 10000,
  n_sims = 3,        # Number of outer simulations
  min_size = 10,
  max_size = 500,    # Use smaller max size initially for faster testing
  step = 5,
  num_trials = 10     # Number of trials per simulation - increase for more variation
)

# Analyze and plot
analysis_normal <- analyze_multiple_simulations(results_normal)
plots_normal <- plot_multiple_simulation_results(
  analysis_results = analysis_normal,
  data_gen_name = "Bivariate Normal",
  ylim_logLik = c(1.0, 5),            # Custom range for logLik plot
  ylim_param_l2 = c(0.5, 10.0),       # Custom range for param_l2 plot
  ylim_lambda = c(0, 1.5),            # Custom range for lambda plot
  ylim_total_time = c(0, 5),
  xlim = c(0, 200)                    # Custom range for total_time plot
)

# Display plots
plots_normal$logLik
plots_normal$param_l2
plots_normal$lambda
plots_normal$total_time

# Comparison table for coreset size = 30
table_30_normal <- generate_performance_table(analysis_normal, 30, "Bivariate Normal")

# Comparison table for coreset size = 100
table_100_normal <- generate_performance_table(analysis_normal, 100, "Bivariate Normal")

# --- Example 3: Bivariate Normal Mixture ---
results_normal_mixture <- run_multiple_simulations(
  data_gen_fn = generate_bivariate_normal_mixture,
  n_samples = 10000,
  n_sims = 3,        # Number of outer simulations
  min_size = 10,
  max_size = 500,    # Use smaller max size initially for faster testing
  step = 5,
  num_trials = 10     # Number of trials per simulation - increase for more variation
)

# Analyze and plot
analysis_normal_mixture <- analyze_multiple_simulations(results_normal_mixture)
plots_normal_mixture <- plot_multiple_simulation_results(
  analysis_results = analysis_normal_mixture,
  data_gen_name = "Bivariate Normal Mixture",
  ylim_logLik = c(1.0, 1.5),          # Custom range for logLik plot
  ylim_param_l2 = c(0.5, 7),          # Custom range for param_l2 plot
  ylim_lambda = c(0, 1),              # Custom range for lambda plot
  ylim_total_time = c(0, 1),
  xlim = c(0, 200)                    # Custom range for total_time plot
)

# Display plots
plots_normal_mixture$logLik
plots_normal_mixture$param_l2
plots_normal_mixture$lambda
plots_normal_mixture$total_time

# Comparison table for coreset size = 30
table_30_mixture <- generate_performance_table(analysis_normal_mixture, 30, "Bivariate Normal Mixture")

# Comparison table for coreset size = 100
table_100_mixture <- generate_performance_table(analysis_normal_mixture, 100, "Bivariate Normal Mixture")

# --- Example 4: Geometric Mixed Distribution ---
results_geometric <- run_multiple_simulations(
  data_gen_fn = generate_mixture,
  n_samples = 10000,
  n_sims = 3,        # Number of outer simulations
  min_size = 10,
  max_size = 500,    # Use smaller max size initially for faster testing
  step = 5,
  num_trials = 10     # Number of trials per simulation - increase for more variation
)

analysis_geometric <- analyze_multiple_simulations(results_geometric)
plots_geometric <- plot_multiple_simulation_results(
  analysis_results = analysis_geometric,
  data_gen_name = "Geometric Mixed Distribution",
  ylim_logLik = c(1.0, 5),            # Custom range for logLik plot
  ylim_param_l2 = c(0.5, 15),         # Custom range for param_l2 plot
  ylim_lambda = c(0, 0.3),            # Custom range for lambda plot
  ylim_total_time = c(0, 1),
  xlim = c(0, 150)
)

# Display plots
plots_geometric$logLik
plots_geometric$param_l2
plots_geometric$lambda
plots_geometric$total_time

# Comparison table for coreset size = 30
table_30_geometric <- generate_performance_table(analysis_geometric, 30, "Geometric Mixed Distribution")

# Comparison table for coreset size = 100
table_100_geometric <- generate_performance_table(analysis_geometric, 100, "Geometric Mixed Distribution")

# --- Example 5: Skew-t Distribution ---
results_skewt <- run_multiple_simulations(
  data_gen_fn = generate_skewt,
  n_samples = 10000,
  n_sims = 3,        # Number of outer simulations
  min_size = 10,
  max_size = 500,    # Use smaller max size initially for faster testing
  step = 5,
  num_trials = 10
)

analysis_skewt <- analyze_multiple_simulations(results_skewt)
plots_skewt <- plot_multiple_simulation_results(
  analysis_results = analysis_skewt,
  data_gen_name = "Skew-t Distribution",
  ylim_logLik = c(1.0, 3),
  ylim_param_l2 = c(0.5, 10),
  ylim_lambda = c(0, 2),
  ylim_total_time = c(0, 1),
  xlim = c(0, 200)
)

plots_skewt$logLik
plots_skewt$param_l2
plots_skewt$lambda
plots_skewt$total_time

# Comparison table for coreset size = 30
table_30_skewt <- generate_performance_table(analysis_skewt, 30, "Skew-t Distribution")

# Comparison table for coreset size = 100
table_100_skewt <- generate_performance_table(analysis_skewt, 100, "Skew-t Distribution")

# --- Example 6: Heteroscedastic Distribution ---
results_heteroscedastic <- run_multiple_simulations(
  data_gen_fn = generate_heteroscedastic,
  n_samples = 10000,
  n_sims = 3,        # Number of outer simulations
  min_size = 10,
  max_size = 500,    # Use smaller max size initially for faster testing
  step = 5,
  num_trials = 10
)

analysis_heteroscedastic <- analyze_multiple_simulations(results_heteroscedastic)
plots_heteroscedastic <- plot_multiple_simulation_results(
  analysis_results = analysis_heteroscedastic,
  data_gen_name = "Heteroscedastic Distribution",
  xlim = c(0, 200),
  ylim_logLik = c(1.0, 1.5),
  ylim_param_l2 = c(1, 10),
  ylim_lambda = c(0, 1),
  ylim_total_time = c(0, 1)
)

plots_heteroscedastic$logLik
plots_heteroscedastic$param_l2
plots_heteroscedastic$lambda
plots_heteroscedastic$total_time

# Comparison table for coreset size = 30
table_30_heteroscedastic <- generate_performance_table(analysis_heteroscedastic, 30, "Heteroscedastic Distribution")

# Comparison table for coreset size = 100
table_100_heteroscedastic <- generate_performance_table(analysis_heteroscedastic, 100, "Heteroscedastic Distribution")

# --- Example 7: Copula Complex Distribution ---
results_copula <- run_multiple_simulations(
  data_gen_fn = generate_copula_complex,
  n_samples = 10000,
  n_sims = 3,        # Number of outer simulations
  min_size = 10,
  max_size = 500,    # Use smaller max size initially for faster testing
  step = 5,
  num_trials = 10
)

analysis_copula <- analyze_multiple_simulations(results_copula)
plots_copula <- plot_multiple_simulation_results(
  analysis_results = analysis_copula,
  xlim = c(0, 200),
  data_gen_name = "Copula Complex Distribution",
  ylim_logLik = c(1.0, 3),
  ylim_param_l2 = c(1, 10),
  ylim_lambda = c(0, 1),
  ylim_total_time = c(0, 1)
)

plots_copula$logLik
plots_copula$param_l2
plots_copula$lambda
plots_copula$total_time

# Comparison table for coreset size = 30
table_30_copula <- generate_performance_table(analysis_copula, 30, "Copula Complex Distribution")

# Comparison table for coreset size = 100
table_100_copula <- generate_performance_table(analysis_copula, 100, "Copula Complex Distribution")

# --- Example 8: Spiral Dependency ---
results_spiral <- run_multiple_simulations(
  data_gen_fn = generate_spiral_dependency,
  n_samples = 10000,
  n_sims = 3,        # Number of outer simulations
  min_size = 10,
  max_size = 500,    # Use smaller max size initially for faster testing
  step = 5,
  num_trials = 10     # Number of trials per simulation - increase for more variation
)

# Analyze and plot
analysis_spiral <- analyze_multiple_simulations(results_spiral)
plots_spiral <- plot_multiple_simulation_results(
  analysis_results = analysis_spiral,
  data_gen_name = "Spiral Dependency",
  ylim_logLik = c(1.0, 1.5),          # Custom range for logLik plot
  ylim_param_l2 = c(0.5, 7.0),        # Custom range for param_l2 plot
  ylim_lambda = c(0, 0.5),            # Custom range for lambda plot
  ylim_total_time = c(0, 1),
  xlim = c(0, 200)                    # Custom range for total_time plot
)

# Display plots
plots_spiral$logLik
plots_spiral$param_l2
plots_spiral$lambda
plots_spiral$total_time

# Comparison table for coreset size = 30
table_30_spiral <- generate_performance_table(analysis_spiral, 30, "Spiral Dependency")

# Comparison table for coreset size = 100
table_100_spiral <- generate_performance_table(analysis_spiral, 100, "Spiral Dependency")

# --- Example 9: Circular Dependency ---
results_circular <- run_multiple_simulations(
  data_gen_fn = generate_circular_dependency,
  n_samples = 10000,
  n_sims = 3,        # Number of outer simulations
  min_size = 10,
  max_size = 500,    # Use smaller max size initially for faster testing
  step = 5,
  num_trials = 10     # Number of trials per simulation - increase for more variation
)

# Analyze and plot
analysis_circular <- analyze_multiple_simulations(results_circular)
plots_circular <- plot_multiple_simulation_results(
  analysis_results = analysis_circular,
  xlim = c(0, 200),
  data_gen_name = "Circular Dependency",
  ylim_logLik = c(1.0, 1.25),         # Custom range for logLik plot
  ylim_param_l2 = c(0.5, 5.0),        # Custom range for param_l2 plot
  ylim_lambda = c(0, 0.1),            # Custom range for lambda plot
  ylim_total_time = c(0, 1)           # Custom range for total_time plot
)

# Display plots
plots_circular$logLik
plots_circular$param_l2
plots_circular$lambda
plots_circular$total_time

# Comparison table for coreset size = 30
table_30_circular <- generate_performance_table(analysis_circular, 30, "Circular Dependency")

# Comparison table for coreset size = 100
table_100_circular <- generate_performance_table(analysis_circular, 100, "Circular Dependency")

# --- Example 10: t-Copula Dependency ---
results_tcopula <- run_multiple_simulations(
  data_gen_fn = generate_t_copula_dependency,
  n_samples = 10000,
  n_sims = 3,        # Number of outer simulations
  min_size = 10,
  max_size = 500,    # Use smaller max size initially for faster testing
  step = 5,
  num_trials = 10     # Number of trials per simulation - increase for more variation
)

# Analyze and plot
analysis_tcopula <- analyze_multiple_simulations(results_tcopula)
plots_tcopula <- plot_multiple_simulation_results(
  analysis_results = analysis_tcopula,
  data_gen_name = "t-Copula Dependency",
  ylim_logLik = c(1.0, 1.5),          # Custom range for logLik plot
  ylim_param_l2 = c(0.5, 8.0),        # Custom range for param_l2 plot
  ylim_lambda = c(0, 1),              # Custom range for lambda plot
  ylim_total_time = c(0, 1),
  xlim = c(0, 200)                    # Custom range for total_time plot
)

# Display plots
plots_tcopula$logLik
plots_tcopula$param_l2
plots_tcopula$lambda
plots_tcopula$total_time

# Comparison table for coreset size = 30
table_30_tcopula <- generate_performance_table(analysis_tcopula, 30, "t-Copula Dependency")

# Comparison table for coreset size = 100
table_100_tcopula <- generate_performance_table(analysis_tcopula, 100, "t-Copula Dependency")

# --- Example 11: Piecewise Dependency ---
results_piecewise <- run_multiple_simulations(
  data_gen_fn = generate_piecewise_dependency,
  n_samples = 10000,
  n_sims = 3,        # Number of outer simulations
  min_size = 10,
  max_size = 500,    # Use smaller max size initially for faster testing
  step = 5,
  num_trials = 10     # Number of trials per simulation - increase for more variation
)

# Analyze and plot
analysis_piecewise <- analyze_multiple_simulations(results_piecewise)
plots_piecewise <- plot_multiple_simulation_results(
  analysis_results = analysis_piecewise,
  data_gen_name = "Piecewise Dependency",
  ylim_logLik = c(1.0, 1.5),          # Custom range for logLik plot
  ylim_param_l2 = c(0.5, 5.0),        # Custom range for param_l2 plot
  ylim_lambda = c(0, 0.5),            # Custom range for lambda plot
  ylim_total_time = c(0, 1),
  xlim = c(0, 200)                    # Custom range for total_time plot
)

# Display plots
plots_piecewise$logLik
plots_piecewise$param_l2
plots_piecewise$lambda
plots_piecewise$total_time

# Comparison table for coreset size = 30
table_30_piecewise <- generate_performance_table(analysis_piecewise, 30, "Piecewise Dependency")

# Comparison table for coreset size = 100
table_100_piecewise <- generate_performance_table(analysis_piecewise, 100, "Piecewise Dependency")

# --- Example 12: Hourglass Dependency ---
results_hourglass <- run_multiple_simulations(
  data_gen_fn = generate_hourglass_dependency,
  n_samples = 10000,
  n_sims = 3,        # Number of outer simulations
  min_size = 10,
  max_size = 500,    # Use smaller max size initially for faster testing
  step = 5,
  num_trials = 10     # Number of trials per simulation - increase for more variation
)

# Analyze and plot
analysis_hourglass <- analyze_multiple_simulations(results_hourglass)
plots_hourglass <- plot_multiple_simulation_results(
  analysis_results = analysis_hourglass,
  data_gen_name = "Hourglass Dependency",
  ylim_logLik = c(1.0, 1.5),          # Custom range for logLik plot
  ylim_param_l2 = c(0.5, 5.0),        # Custom range for param_l2 plot
  ylim_lambda = c(0, 1),              # Custom range for lambda plot
  ylim_total_time = c(0, 1),
  xlim = c(0, 200)                    # Custom range for total_time plot
)

# Display plots
plots_hourglass$logLik
plots_hourglass$param_l2
plots_hourglass$lambda
plots_hourglass$total_time

# Comparison table for coreset size = 30
table_30_hourglass <- generate_performance_table(analysis_hourglass, 30, "Hourglass Dependency")

# Comparison table for coreset size = 100
table_100_hourglass <- generate_performance_table(analysis_hourglass, 100, "Hourglass Dependency")

# --- Example 13: Bimodal Clusters ---
results_bimodal <- run_multiple_simulations(
  data_gen_fn = generate_bimodal_clusters,
  n_samples = 10000,
  n_sims = 3,        # Number of outer simulations
  min_size = 10,
  max_size = 500,    # Use smaller max size initially for faster testing
  step = 5,
  num_trials = 10     # Number of trials per simulation - increase for more variation
)

# Analyze and plot
analysis_bimodal <- analyze_multiple_simulations(results_bimodal)
plots_bimodal <- plot_multiple_simulation_results(
  analysis_results = analysis_bimodal,
  data_gen_name = "Bimodal Clusters",
  ylim_logLik = c(1.0, 1.5),          # Custom range for logLik plot
  ylim_param_l2 = c(0.5, 5.0),        # Custom range for param_l2 plot
  ylim_lambda = c(0, 1),              # Custom range for lambda plot
  ylim_total_time = c(0, 1),
  xlim = c(0, 200)                    # Custom range for total_time plot
)

# Display plots
plots_bimodal$logLik
plots_bimodal$param_l2
plots_bimodal$lambda
plots_bimodal$total_time

# Comparison table for coreset size = 30
table_30_bimodal <- generate_performance_table(analysis_bimodal, 30, "Bimodal Clusters")

# Comparison table for coreset size = 100
table_100_bimodal <- generate_performance_table(analysis_bimodal, 100, "Bimodal Clusters")

# --- Example 14: Sinusoidal Dependency ---
results_sinusoidal <- run_multiple_simulations(
  data_gen_fn = generate_sinusoidal_dependency,
  n_samples = 10000,
  n_sims = 3,        # Number of outer simulations
  min_size = 10,
  max_size = 500,    # Use smaller max size initially for faster testing
  step = 5,
  num_trials = 10     # Number of trials per simulation - increase for more variation
)

# Analyze and plot
analysis_sinusoidal <- analyze_multiple_simulations(results_sinusoidal)
plots_sinusoidal <- plot_multiple_simulation_results(
  analysis_results = analysis_sinusoidal,
  data_gen_name = "Sinusoidal Dependency",
  ylim_logLik = c(1.0, 1.3),          # Custom range for logLik plot
  ylim_param_l2 = c(0.5, 3.0),        # Custom range for param_l2 plot
  ylim_lambda = c(0, 1),              # Custom range for lambda plot
  ylim_total_time = c(0, 1),
  xlim = c(0, 200)                    # Custom range for total_time plot
)

# Display plots
plots_sinusoidal$logLik
plots_sinusoidal$param_l2
plots_sinusoidal$lambda
plots_sinusoidal$total_time

# Comparison table for coreset size = 30
table_30_sinusoidal <- generate_performance_table(analysis_sinusoidal, 30, "Sinusoidal Dependency")

# Comparison table for coreset size = 100
table_100_sinusoidal <- generate_performance_table(analysis_sinusoidal, 100, "Sinusoidal Dependency")

# ==========================================================================
# Summary Tables for All Analyses
# ==========================================================================

# Comprehensive analysis objects with corresponding names
all_analyses <- list(
  "Non-linear Correlation" = analysis_nonlinear,
  "Bivariate Normal" = analysis_normal,
  "Bivariate Normal Mixture" = analysis_normal_mixture,
  "Mixture Distribution" = analysis_geometric,
  "Skew-t Distribution" = analysis_skewt,
  "Heteroscedastic Distribution" = analysis_heteroscedastic,
  "Copula Complex Distribution" = analysis_copula,
  "Spiral Dependency" = analysis_spiral,
  "Circular Dependency" = analysis_circular,
  "t-Copula Dependency" = analysis_tcopula,
  "Piecewise Dependency" = analysis_piecewise,
  "Hourglass Dependency" = analysis_hourglass,
  "Bimodal Clusters" = analysis_bimodal,
  "Sinusoidal Dependency" = analysis_sinusoidal
)

# Generate and print performance tables for all analyses
for (analysis_name in names(all_analyses)) {
  cat("\n==============================\n")
  cat("Performance Table for:", analysis_name, "\n")
  cat("==============================\n\n")
  
  # Current analysis object
  current_analysis <- all_analyses[[analysis_name]]
  
  # Generate coreset size = 30 table
  table_30 <- generate_performance_table(current_analysis, 30, analysis_name)
  cat(">>> Coreset Size = 30:\n")
  print(table_30)
  cat("\n------------------------------\n")
  
  # Generate coreset size = 100 table
  table_100 <- generate_performance_table(current_analysis, 100, analysis_name)
  cat(">>> Coreset Size = 100:\n")
  print(table_100)
  cat("\n\n")
}