#sampling_function.R
# 1. Core sampling methods
library(ggplot2)

sample_uniform <- function(N, size) {
  indices <- sort(unique(sample(1:N, size, replace = TRUE)))
  weights <- rep(N/size, size)
  return(list(indices = indices, weights = weights))
}

sample_l2_only <- function(X_list, size) {
  N <- nrow(X_list[[1]])
  all_indices <- list()
  all_probs <- list()
  
  # Get L2 samples for each matrix
  for(i in seq_along(X_list)) {
    lev_result <- L2_coreset(size, X_list[[i]], 2)
    all_indices[[i]] <- lev_result$sample_index
    all_probs[[i]] <- lev_result$prob
  }
  
  indices <- unique(unlist(all_indices))
  combined <- combine_weights(all_probs, size)
  weights <- combined$weights[indices]
  
  return(list(indices = indices, weights = weights))
}

# Pre-compute hull indices for all matrices
get_separate_hull_indices <- function(X, Y, X_prime, Y_prime) {
  # Compute convex hull only for the derivative matrices
  X_prime_hull <- hull_ind(X_prime)
  Y_prime_hull <- hull_ind(Y_prime)
  
  return(list(
    X_prime_hull = X_prime_hull,
    Y_prime_hull = Y_prime_hull
  ))
}

# Modified sampling function that uses pre-computed hull indices
sample_l2_hull <- function(X, Y, X_prime, Y_prime, size, hull_indices) {
  N <- nrow(X)
  
  # 1. Perform leverage score sampling for X
  X_lev_result <- L2_coreset(size, X, 2)
  X_samples <- X_lev_result$sample_index
  X_probs <- X_lev_result$prob
  
  # 2. Perform leverage score sampling for Y
  Y_lev_result <- L2_coreset(size, Y, 2)
  Y_samples <- Y_lev_result$sample_index
  Y_probs <- Y_lev_result$prob
  
  # 3. Uniform sampling from the convex hull of X_prime
  X_hull_size <- length(hull_indices$X_prime_hull)
  X_hull_sample_size <- min(ceiling(size * 0.2), X_hull_size)
  if(X_hull_sample_size < X_hull_size) {
    X_hull_samples <- sort(sample(hull_indices$X_prime_hull, X_hull_sample_size, replace = FALSE))
    X_hull_probs <- rep(X_hull_sample_size/X_hull_size, length(X_hull_samples))
  } else {
    X_hull_samples <- hull_indices$X_prime_hull
    X_hull_probs <- rep(1, X_hull_size)
  }
  
  # 4. Uniform sampling from the convex hull of Y_prime
  Y_hull_size <- length(hull_indices$Y_prime_hull)
  Y_hull_sample_size <- min(ceiling(size * 0.2), Y_hull_size)
  if(Y_hull_sample_size < Y_hull_size) {
    Y_hull_samples <- sort(sample(hull_indices$Y_prime_hull, Y_hull_sample_size, replace = FALSE))
    Y_hull_probs <- rep(Y_hull_sample_size/Y_hull_size, length(Y_hull_samples))
  } else {
    Y_hull_samples <- hull_indices$Y_prime_hull
    Y_hull_probs <- rep(1, Y_hull_size)
  }
  
  # 5. Combine all sampled indices
  all_indices <- unique(c(X_samples, Y_samples, X_hull_samples, Y_hull_samples))
  
  # 6. Construct probability lists
  all_probs <- list()
  
  # Leverage score probability for X
  X_lev_probs <- numeric(N)
  X_lev_probs[X_samples] <- X_probs
  all_probs$X_lev <- X_lev_probs
  
  # Leverage score probability for Y
  Y_lev_probs <- numeric(N)
  Y_lev_probs[Y_samples] <- Y_probs
  all_probs$Y_lev <- Y_lev_probs
  
  # Probability for X_prime convex hull
  X_hull_full_probs <- numeric(N)
  X_hull_full_probs[X_hull_samples] <- X_hull_probs
  all_probs$X_hull <- X_hull_full_probs
  
  # Probability for Y_prime convex hull
  Y_hull_full_probs <- numeric(N)
  Y_hull_full_probs[Y_hull_samples] <- Y_hull_probs
  all_probs$Y_hull <- Y_hull_full_probs
  
  # 7. Combine probabilities and compute weights
  combined <- combine_weights(all_probs, size)
  weights <- combined$weights[all_indices]
  
  return(list(
    indices = all_indices,
    weights = weights,
    components = list(
      X_samples = X_samples,
      Y_samples = Y_samples,
      X_hull_samples = X_hull_samples,
      Y_hull_samples = Y_hull_samples
    )
  ))
}

# 2. Weight combination function
# Modified combine_weights function
combine_weights <- function(prob_list, size, N) {  # Added N parameter
  # Get all indices and compute combined probabilities
  n <- length(prob_list[[1]])
  combined_probs <- numeric(n)
  
  # Compute p_i = 1 - ∏(1 - p_i^(f))
  valid_indices <- which(rowSums(do.call(cbind, prob_list)) > 0)
  
  for(i in valid_indices) {
    prob_product <- 1
    for(probs in prob_list) {
      if(probs[i] > 0) {
        prob_product <- prob_product * (1 - probs[i])
      }
    }
    combined_probs[i] <- 1 - prob_product
  }
  
  # Compute initial weights: w_i = 1/(p_i * size)
  weights <- numeric(n)
  positive_probs <- combined_probs > 0
  weights[positive_probs] <- 1/(combined_probs[positive_probs] * size)
  
  # Normalize weights so that the sum equals N
  total_weight <- sum(weights)
  weights <- weights * (N/total_weight)
  
  return(list(probs = combined_probs, weights = weights))
}

# 3. Main sampling dispatcher

# Modified dispatcher to include hull_indices
sample_method <- function(method, X_list, size, hull_indices = NULL) {
  N <- nrow(X_list[[1]])
  
  result <- switch(method,
                   "uniform" = sample_uniform(N, size),
                   "l2_only" = sample_l2_only(X_list, size),
                   "l2_hull" = sample_l2_hull(X_list, size, hull_indices),
                   stop("Unknown sampling method")
  )
  
  return(result)
}

# 4. Main simulation function

run_mctm_comparison <- function(data, min_size = 10, max_size = NULL, step = 10, 
                                num_trials = 5) {
  N <- nrow(data)
  if(is.null(max_size)) max_size <- floor(0.5 * N)
  
  # Initialize results
  comparison_results <- list()
  
  # Get baseline results
  cat("Computing baseline...\n")
  start_time <- proc.time()
  baseline_fit <- mmlt_continuous(data, method = "BFGS",
                                  control = list(
                                    trace = 0,
                                    maxit = 1000,
                                    reltol = 1e-15
                                  ),
                                  random_init = TRUE)
  baseline_time <- (proc.time() - start_time)[3]  
  
  comparison_results$baseline <- list(
    coefficients = baseline_fit$coefficients,
    logLik = baseline_fit$logLik,
    time = baseline_time
  )
  
  # Get model matrices
  x <- data[,1]
  y <- data[,2]
  m_x <- tram::BoxCox(x ~ 1)
  m_y <- tram::BoxCox(y ~ 1)
  models <- .models(m_x, m_y)
  
  # Get the four matrices
  X <- models$mm[[1]]$eY$Y
  Y <- models$mm[[2]]$eY$Y
  X_prime <- models$mm[[1]]$eY$Yprime
  Y_prime <- models$mm[[2]]$eY$Yprime
  
  # Pre-compute hull indices
  cat("Computing convex hulls...\n")
  hull_indices <- get_separate_hull_indices(X, Y, X_prime, Y_prime)
  
  # Run simulations
  methods <- c("uniform", "l2_only", "l2_hull")
  sizes <- seq(min_size, max_size, by = step)
  
  for(size in sizes) {
    cat(sprintf("\nSize: %d\n", size))
    size_results <- list()
    
    for(method in methods) {
      cat(sprintf("Method: %s\n", method))
      method_results <- list()
      
      for(trial in 1:num_trials) {
        cat(sprintf("Trial: %d/%d\n", trial, num_trials))
        
        # Record sampling time
        sample_start_time <- proc.time()
        
        # Get samples and weights based on method
        
        if(method == "uniform") {
          indices <- sort(sample(1:N, size, replace = FALSE))
          weights <- rep(N/size, length(indices))  
          sample_result <- list(indices = indices, weights = weights)
          
        } else if(method == "l2_only") {
          # Leverage score sampling
          X_lev <- L2_coreset(size, X, 2)
          Y_lev <- L2_coreset(size, Y, 2)
          
          # Combine indices
          all_indices <- sort(unique(c(X_lev$sample_index, 
                                       Y_lev$sample_index)))
          
          # Create probability vectors
          X_probs <- numeric(N)
          Y_probs <- numeric(N)
          
          # Assign probabilities
          for(i in seq_along(X_lev$sample_index)) {
            X_probs[X_lev$sample_index[i]] <- X_lev$prob[i]
          }
          for(i in seq_along(Y_lev$sample_index)) {
            Y_probs[Y_lev$sample_index[i]] <- Y_lev$prob[i]
          }
          
          # Combine probabilities and compute normalized weights
          combined <- combine_weights(list(X = X_probs, 
                                           Y = Y_probs), size, N)
          
          sample_result <- list(
            indices = all_indices,
            weights = combined$weights[all_indices]
          )
          
        } else if(method == "l2_hull") {
          # Leverage score sampling
          X_lev <- L2_coreset(size, X, 2)
          Y_lev <- L2_coreset(size, Y, 2)
          
          # Hull sampling
          X_hull_size <- min(ceiling(size * 0.2), 
                             length(hull_indices$X_prime_hull))
          Y_hull_size <- min(ceiling(size * 0.2), 
                             length(hull_indices$Y_prime_hull))
          
          X_hull_samples <- sort(sample(hull_indices$X_prime_hull, X_hull_size, replace = FALSE))
          Y_hull_samples <- sort(sample(hull_indices$Y_prime_hull, Y_hull_size, replace = FALSE))
          
          # Combine all indices
          all_indices <- sort(unique(c(
            X_lev$sample_index, Y_lev$sample_index,
            X_hull_samples, Y_hull_samples
          )))
          
          # Create probability vectors
          X_lev_probs <- numeric(N)
          Y_lev_probs <- numeric(N)
          X_hull_probs <- numeric(N)
          Y_hull_probs <- numeric(N)
          
          # Assign leverage score probabilities
          for(i in seq_along(X_lev$sample_index)) {
            X_lev_probs[X_lev$sample_index[i]] <- X_lev$prob[i]
          }
          for(i in seq_along(Y_lev$sample_index)) {
            Y_lev_probs[Y_lev$sample_index[i]] <- Y_lev$prob[i]
          }
          
          # Assign hull probabilities
          for(i in seq_along(X_hull_samples)) {
            X_hull_probs[X_hull_samples[i]] <- X_hull_size/length(hull_indices$X_prime_hull)
          }
          for(i in seq_along(Y_hull_samples)) {
            Y_hull_probs[Y_hull_samples[i]] <- Y_hull_size/length(hull_indices$Y_prime_hull)
          }
          
          # Combine probabilities and compute normalized weights
          combined <- combine_weights(
            list(
              X_lev = X_lev_probs,
              Y_lev = Y_lev_probs,
              X_hull = X_hull_probs,
              Y_hull = Y_hull_probs
            ),
            size,
            N
          )
          
          sample_result <- list(
            indices = all_indices,
            weights = combined$weights[all_indices],
            components = list(
              X_lev = X_lev$sample_index,
              Y_lev = Y_lev$sample_index,
              X_hull = X_hull_samples,
              Y_hull = Y_hull_samples
            )
          )
        }
        
        sampling_time <- (proc.time() - sample_start_time)[3]
        
        # Record optimization time
        optim_start_time <- proc.time()
        
        # Run optimization with sampled data
        coreset_data <- data[sample_result$indices, ]
        coreset_fit <- try(mmlt_continuous(
          coreset_data, 
          method = "BFGS",
          control = list(
            trace = 0,
            maxit = 1000,
            reltol = 1e-15
          ),
          weights = sample_result$weights,
          random_init = TRUE
        ))
        
        optimization_time <- (proc.time() - optim_start_time)[3]
        total_time <- sampling_time + optimization_time
        
        # Store results
        if(!inherits(coreset_fit, "try-error")) {
          method_results[[trial]] <- list(
            coefficients = coreset_fit$coefficients,
            logLik = as.numeric(coreset_fit$logLik),
            size = length(sample_result$indices),
            indices = sample_result$indices,
            weights = sample_result$weights,
            sampling_time = sampling_time,
            optimization_time = optimization_time,
            total_time = total_time
          )
          
          if(method == "l2_hull") {
            method_results[[trial]]$components <- sample_result$components
          }
        } else {
          warning(sprintf("Optimization failed for size %d, method %s, trial %d", 
                          size, method, trial))
          method_results[[trial]] <- list(
            coefficients = rep(NA, length(baseline_fit$coefficients)),
            logLik = NA_real_,
            size = length(sample_result$indices),
            indices = sample_result$indices,
            weights = sample_result$weights,
            sampling_time = sampling_time,
            optimization_time = optimization_time,
            total_time = total_time
          )
        }
      }
      
      size_results[[method]] <- method_results
    }
    
    comparison_results[[sprintf("size_%d", size)]] <- size_results
  }
  
  comparison_results$hull_indices <- hull_indices
  comparison_results$matrices <- list(
    X = X,
    Y = Y,
    X_prime = X_prime,
    Y_prime = Y_prime
  )
  comparison_results$data <- data
  
  return(comparison_results)
}

# Wrapper function to run multiple simulations with different datasets
run_multiple_simulations <- function(data_gen_fn = NULL, dataset = NULL,n_samples, n_sims = 5, 
                                     min_size = 10, max_size = NULL, step = 10) {
  # Initialize results list
  all_results <- list()
  
  # Run simulations
  for(i in 1:n_sims) {
    cat(sprintf("\nSimulation %d/%d\n", i, n_sims))
    
    # Generate new dataset
    
    set.seed(round(runif(1,1,10000)))
    data <- data_gen_fn(n = n_samples)
    
    # Fit full model for this dataset
    x <- data[,1]
    y <- data[,2]
    m_x <- tram::BoxCox(x ~ 1)
    m_y <- tram::BoxCox(y ~ 1)
    m_xy <- mmlt(m_x, m_y, formula = ~ 1)
    
    # Store the model in the results
    sim_result <- run_mctm_comparison(
      data = data,
      min_size = min_size,
      max_size = max_size,
      step = step,
      num_trials = 1  # Only need 1 trial per dataset now
    )
    
    # Add the full model to results
    sim_result$full_model <- m_xy
    
    # Store results
    all_results[[i]] <- sim_result
  }
  
  return(all_results)
}


# Modified analysis function to include min, max, and median
analyze_multiple_simulations <- function(all_results) {
  n_sims <- length(all_results)
  
  sizes <- names(all_results[[1]])[!names(all_results[[1]]) %in% 
                                     c("baseline", "hull_indices", "matrices", "data", "full_model")]
  
  analysis <- list()
  
  for(size in sizes) {
    size_num <- as.numeric(gsub("size_", "", size))
    methods <- names(all_results[[1]][[size]])
    
    method_analysis <- list()
    for(method in methods) {
      logLik_values <- numeric(n_sims)
      param_values <- list()
      actual_sizes <- numeric(n_sims)
      sampling_times <- numeric(n_sims)
      optimization_times <- numeric(n_sims)
      total_times <- numeric(n_sims)
      
      for(sim in 1:n_sims) {
        sim_results <- all_results[[sim]]
        m_xy <- sim_results$full_model
        trial_result <- sim_results[[size]][[method]][[1]]
        
        full_loglik <- m_xy$ll(m_xy$par)
        coreset_loglik <- m_xy$ll(trial_result$coefficients)
        logLik_values[sim] <- abs(coreset_loglik / full_loglik)
        
        true_params <- m_xy$par
        coreset_params <- trial_result$coefficients
        param_values[[sim]] <- coreset_params - true_params
        
        actual_sizes[sim] <- trial_result$size
        sampling_times[sim] <- trial_result$sampling_time
        optimization_times[sim] <- trial_result$optimization_time
        total_times[sim] <- trial_result$total_time
      }
      
      param_matrix <- do.call(rbind, param_values)
      param_l2_dist <- sqrt(rowSums(param_matrix^2))
      lambda_errors <- abs(param_matrix[, ncol(param_matrix)])
      
      method_analysis[[method]] <- list(
        size = size_num,
        # Log-likelihood statistics
        logLik_median = median(logLik_values),
        logLik_min = min(logLik_values),
        logLik_max = max(logLik_values),
        # Parameter L2 distance statistics
        param_l2_median = median(param_l2_dist),
        param_l2_min = min(param_l2_dist),
        param_l2_max = max(param_l2_dist),
        # Lambda error statistics
        lambda_error_median = median(lambda_errors),
        lambda_error_min = min(lambda_errors),
        lambda_error_max = max(lambda_errors),
        # Time statistics (keeping mean/sd for these)
        actual_size = mean(actual_sizes),
        mean_sampling_time = mean(sampling_times),
        sd_sampling_time = sd(sampling_times),
        mean_optimization_time = mean(optimization_times),
        sd_optimization_time = sd(optimization_times),
        mean_total_time = mean(total_times),
        sd_total_time = sd(total_times)
      )
    }
    analysis[[size]] <- method_analysis
  }
  
  return(analysis)
}

# Modified plotting function to show min, max, and median
plot_multiple_simulation_results <- function(analysis_results) {
  sizes <- names(analysis_results)
  methods <- names(analysis_results[[1]])
  
  metrics_df <- data.frame()
  for(size in sizes) {
    size_num <- as.numeric(gsub("size_", "", size))
    for(method in methods) {
      res <- analysis_results[[size]][[method]]
      metrics_df <- rbind(metrics_df, data.frame(
        size = size_num,
        method = method,
        # Log-likelihood metrics
        logLik_median = res$logLik_median,
        logLik_min = res$logLik_min,
        logLik_max = res$logLik_max,
        # Parameter L2 metrics
        param_l2_median = res$param_l2_median,
        param_l2_min = res$param_l2_min,
        param_l2_max = res$param_l2_max,
        # Lambda error metrics
        lambda_error_median = res$lambda_error_median,
        lambda_error_min = res$lambda_error_min,
        lambda_error_max = res$lambda_error_max,
        # Time metrics (keeping mean/sd)
        sampling_time = res$mean_sampling_time,
        sampling_time_sd = res$sd_sampling_time,
        optimization_time = res$mean_optimization_time,
        optimization_time_sd = res$sd_optimization_time,
        total_time = res$mean_total_time,
        total_time_sd = res$sd_total_time
      ))
    }
  }
  
  # Create plots using ggplot2
  plots <- list()
  
  # Common theme
  base_theme <- theme_minimal() +
    theme(legend.title = element_blank(),
          legend.position = "bottom")
  
  # Log-likelihood ratio plot with min/max range
  plots$logLik <- ggplot(metrics_df, aes(x = size, color = method)) +
    geom_line(aes(y = logLik_median)) +
    geom_ribbon(aes(ymin = logLik_min,
                    ymax = logLik_max,
                    fill = method),
                alpha = 0.2) +
    labs(title = "Log-Likelihood Ratio vs Coreset Size",
         x = "Coreset Size",
         y = "Likelihood Ratio") +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    base_theme #+
  #coord_cartesian(ylim = c(0.9, 1.1))
  
  # Parameter L2 distance plot with min/max range
  plots$param_l2 <- ggplot(metrics_df, aes(x = size, color = method)) +
    geom_line(aes(y = param_l2_median)) +
    geom_ribbon(aes(ymin = param_l2_min,
                    ymax = param_l2_max,
                    fill = method),
                alpha = 0.2) +
    labs(title = "Parameter L2 Distance vs Coreset Size",
         x = "Coreset Size",
         y = "L2 Distance") +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    base_theme #+
  #coord_cartesian(ylim = c(0, 1.1))
  
  # Lambda error plot with min/max range
  plots$lambda <- ggplot(metrics_df, aes(x = size, color = method)) +
    geom_line(aes(y = lambda_error_median)) +
    geom_ribbon(aes(ymin = lambda_error_min,
                    ymax = lambda_error_max,
                    fill = method),
                alpha = 0.2) +
    labs(title = "Lambda Parameter Error vs Coreset Size",
         x = "Coreset Size",
         y = "Absolute Error") +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    base_theme #+
  #coord_cartesian(ylim = c(0, 0.5))
  
  # Time plots remain the same with mean/sd
  plots$sampling_time <- ggplot(metrics_df, aes(x = size, 
                                                y = sampling_time,
                                                color = method, 
                                                fill = method)) +
    geom_line() +
    geom_point() +
    geom_ribbon(aes(ymin = pmax(0, sampling_time - sampling_time_sd),
                    ymax = sampling_time + sampling_time_sd),
                alpha = 0.2) +
    labs(title = "Sampling Time vs Coreset Size",
         x = "Coreset Size",
         y = "Time (seconds)") +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    base_theme
  
  plots$optimization_time <- ggplot(metrics_df, aes(x = size, 
                                                    y = optimization_time,
                                                    color = method, 
                                                    fill = method)) +
    geom_line() +
    geom_point() +
    geom_ribbon(aes(ymin = pmax(0, optimization_time - optimization_time_sd),
                    ymax = optimization_time + optimization_time_sd),
                alpha = 0.2) +
    labs(title = "Optimization Time vs Coreset Size",
         x = "Coreset Size",
         y = "Time (seconds)") +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    base_theme
  
  plots$total_time <- ggplot(metrics_df, aes(x = size, 
                                             y = total_time,
                                             color = method, 
                                             fill = method)) +
    geom_line() +
    geom_point() +
    geom_ribbon(aes(ymin = pmax(0, total_time - total_time_sd),
                    ymax = total_time + total_time_sd),
                alpha = 0.2) +
    labs(title = "Total Time vs Coreset Size",
         x = "Coreset Size",
         y = "Time (seconds)") +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    base_theme
  
  return(plots)
}


results <- run_multiple_simulations(
  data_gen_fn = generate_heteroscedastic,
  n_samples = 10000,
  n_sims = 10,
  min_size = 10,
  max_size = 1000,
  step = 10
)
analysis <- analyze_multiple_simulations(results)
plots <- plot_multiple_simulation_results(analysis)
