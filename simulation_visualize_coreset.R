#simulation_visualize_coreset.R
# ==========================================================================
# Coreset Sampling Visualization Functions
# ==========================================================================

#' Visualize different coreset sampling methods
#'
#' @param data_gen_fn Function that generates bivariate data
#' @param n_samples Number of total samples to generate
#' @param coreset_size Size of the coreset to sample
#' @param plot_title Title for the overall visualization
#' @return A list containing the visualization plot and sampling details
visualize_coreset_sampling <- function(data_gen_fn, n_samples = 1000, coreset_size = 50,
                                       plot_title = "Coreset Sampling Visualization") {
  # Set seed for reproducibility
  set.seed(42)
  
  # Generate data
  data <- data_gen_fn(n = n_samples)
  
  # Setup for computing coresets
  x <- data[,1]
  y <- data[,2]
  m_x <- tram::BoxCox(x ~ 1)
  m_y <- tram::BoxCox(y ~ 1)
  models <- .models(m_x, m_y)
  
  # Get basis function matrices for X and Y
  X <- models$mm[[1]]$eY$Y
  Y <- models$mm[[2]]$eY$Y
  X_prime <- models$mm[[1]]$eY$Yprime
  Y_prime <- models$mm[[2]]$eY$Yprime
  
  # 1. Uniform sampling
  uniform_indices <- sort(sample(1:n_samples, coreset_size*2, replace = FALSE))
  
  # 2. L2-only sampling
  # Get L2 samples for X
  X_lev <- L2_coreset(coreset_size, X, 2)
  Y_lev <- L2_coreset(coreset_size, Y, 2)
  l2_indices <- sort(unique(c(X_lev$sample_index, Y_lev$sample_index)))
  
  # 3. L2-hull sampling
  # Compute convex hulls
  hull_indices <- get_separate_hull_indices(X, Y, X_prime, Y_prime)
  
  # L2 samples from above
  # Hull samples
  X_hull_size <- min(ceiling(coreset_size * 0.2), length(hull_indices$X_prime_hull))
  Y_hull_size <- min(ceiling(coreset_size * 0.2), length(hull_indices$Y_prime_hull))
  
  X_hull_samples <- sort(sample(hull_indices$X_prime_hull,
                                X_hull_size, replace = FALSE))
  Y_hull_samples <- sort(sample(hull_indices$Y_prime_hull,
                                Y_hull_size, replace = FALSE))
  
  # Combine all samples
  l2_hull_indices <- sort(unique(c(X_lev$sample_index, Y_lev$sample_index,
                                   X_hull_samples, Y_hull_samples)))
  
  # Create a data frame for plotting
  plot_data <- data.frame(
    y1 = data[,1],
    y2 = data[,2],
    uniform = FALSE,
    l2_only = FALSE,
    l2_hull = FALSE,
    x_hull = FALSE,
    y_hull = FALSE
  )
  
  # Mark sampled points
  plot_data$uniform[uniform_indices] <- TRUE
  plot_data$l2_only[l2_indices] <- TRUE
  plot_data$l2_hull[l2_hull_indices] <- TRUE
  
  # Mark hull points
  plot_data$x_hull[hull_indices$X_prime_hull] <- TRUE
  plot_data$y_hull[hull_indices$Y_prime_hull] <- TRUE
  
  # Create plots with clean white background and minimal design
  p1 <- ggplot(plot_data, aes(x = y1, y = y2)) +
    # Add contour lines without density fill
    geom_density_2d(color = "darkgray", size = 0.5) +
    # Add background points in light color
    geom_point(color = "grey70", size = 1, alpha = 0.6) +
    # Add highlighted points for uniform sampling
    geom_point(data = subset(plot_data, uniform), aes(x = y1, y = y2),
               color = "red", size = 2.5, alpha = 0.8) +
    labs(title = "Uniform Sampling") +
    theme_minimal() +
    theme(
      legend.position = "none",
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "#EEEEEE"),
      panel.grid.minor = element_blank()
    )
  
  p2 <- ggplot(plot_data, aes(x = y1, y = y2)) +
    # Add contour lines without density fill
    geom_density_2d(color = "darkgray", size = 0.5) +
    # Add background points in light color
    geom_point(color = "grey70", size = 1, alpha = 0.6) +
    # Add highlighted points for L2 sampling
    geom_point(data = subset(plot_data, l2_only), aes(x = y1, y = y2),
               color = "blue", size = 2.5, alpha = 0.8) +
    labs(title = "L2 Sensitivity Sampling") +
    theme_minimal() +
    theme(
      legend.position = "none",
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "#EEEEEE"),
      panel.grid.minor = element_blank()
    )
  
  p3 <- ggplot(plot_data, aes(x = y1, y = y2)) +
    # Add contour lines without density fill
    geom_density_2d(color = "darkgray", size = 0.5) +
    # Add background points in light color
    geom_point(color = "grey70", size = 1, alpha = 0.6) +
    # L2-hull points
    geom_point(data = subset(plot_data, l2_hull), aes(x = y1, y = y2),
               color = "purple", size = 2.5, alpha = 0.8) +
    # Highlight X hull points
    geom_point(data = subset(plot_data, x_hull & l2_hull), aes(x = y1, y = y2),
               color = "green", size = 3, alpha = 0.8, shape = 21, stroke = 1) +
    # Highlight Y hull points
    geom_point(data = subset(plot_data, y_hull & l2_hull), aes(x = y1, y = y2),
               color = "orange", size = 3, alpha = 0.8, shape = 24, stroke = 1) +
    labs(title = "L2-Hull Sampling") +
    theme_minimal() +
    theme(
      legend.position = "none",
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "#EEEEEE"),
      panel.grid.minor = element_blank()
    )
  
  # Combine plots
  combined_plot <- egg::ggarrange(p1, p2, p3, ncol = 1)
  
  # Add overall title
  final_plot <- ggpubr::annotate_figure(combined_plot,
                                        top = ggpubr::text_grob(plot_title,
                                                                face = "bold", size = 14))
  
  return(list(
    plot = final_plot,
    data = plot_data,
    uniform_indices = uniform_indices,
    l2_indices = l2_indices,
    l2_hull_indices = l2_hull_indices,
    hull_indices = hull_indices
  ))
}

# Helper function to get hull indices (copied from sampling_function.R)
get_separate_hull_indices <- function(X, Y, X_prime, Y_prime) {
  # Compute convex hull only for the derivative matrices
  X_prime_hull <- hull_ind(X_prime)
  Y_prime_hull <- hull_ind(Y_prime)
  
  return(list(
    X_prime_hull = X_prime_hull,
    Y_prime_hull = Y_prime_hull
  ))
}

# ==========================================================================
# Visualization of All Data Generation Processes
# ==========================================================================

#' Visualize all data generation processes
#'
#' @param n_samples Number of samples to generate for each process
#' @param save_plots Whether to save the plots as files
#' @param output_dir Directory to save plots if save_plots is TRUE
#' @return A list containing visualization results for all data generation processes
visualize_all_data_processes <- function(n_samples = 1000, save_plots = FALSE, output_dir = "visualization_results") {
  # Create output directory if saving plots
  if(save_plots && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # All data generation functions and their corresponding titles
  data_processes <- list(
    # Basic data generation processes
    list(fn = generate_bivariate_normal, name = "Bivariate Normal Distribution"),
    list(fn = generate_nonlinear_correlation, name = "Non-linear Correlation Structure"),
    list(fn = generate_bivariate_normal_mixture, name = "Bivariate Normal Mixture"),
    list(fn = generate_mixture, name = "Geometric Mixed Distribution"),
    
    list(fn = generate_skewt, name = "Skew-t Distribution"),
    list(fn = generate_heteroscedastic, name = "Heteroscedastic Distribution"),
    list(fn = generate_copula_complex, name = "Copula Complex Distribution"),
    
    # Additional data generation processes
    list(fn = generate_spiral_dependency, name = "Spiral Dependency"),
    list(fn = generate_circular_dependency, name = "Circular Dependency"),
    list(fn = generate_t_copula_dependency, name = "t-Copula Dependency"),
    list(fn = generate_piecewise_dependency, name = "Piecewise Dependency"),
    list(fn = generate_polynomial_dependency, name = "Polynomial Dependency"),
    list(fn = generate_hourglass_dependency, name = "Hourglass Dependency"),
    list(fn = generate_bimodal_clusters, name = "Bimodal Clusters"),
    list(fn = generate_sinusoidal_dependency, name = "Sinusoidal Dependency")
  )
  
  # Create visualization for each data generation process
  results <- list()
  
  for(i in 1:length(data_processes)) {
    process <- data_processes[[i]]
    cat("Processing", process$name, "...\n")
    
    # Create visualization using existing function
    result <- visualize_coreset_sampling(
      data_gen_fn = process$fn,
      n_samples = n_samples,
      coreset_size = 50,  # Default coreset size is 50
      plot_title = process$name
    )
    
    # Save results
    results[[process$name]] <- result
    
    # Save images if needed
    if(save_plots) {
      filename <- file.path(output_dir, paste0(gsub("[^a-zA-Z0-9]", "_", process$name), ".pdf"))
      ggplot2::ggsave(filename, result$plot, width = 10, height = 12, dpi = 300)
      cat("  Image saved to", filename, "\n")
    }
  }
  
  return(results)
}

# ==========================================================================
# Example Usage
# ==========================================================================
# 
# # Example 1: Spiral dependency
# spiral_vis <- visualize_coreset_sampling(
#   data_gen_fn = generate_spiral_dependency,
#   n_samples = 1000,
#   coreset_size = 50,
#   plot_title = "Spiral Dependency - Coreset Sampling Comparison"
# )
# print(spiral_vis$plot)
# 
# # Example 2: Bimodal clusters
# bimodal_vis <- visualize_coreset_sampling(
#   data_gen_fn = generate_bimodal_clusters,
#   n_samples = 1000,
#   coreset_size = 50,
#   plot_title = "Bimodal Clusters - Coreset Sampling Comparison"
# )
# print(bimodal_vis$plot)
# 
# # Example 3: Circular dependency
# circular_vis <- visualize_coreset_sampling(

#   data_gen_fn = generate_circular_dependency,
#   n_samples = 1000,
#   coreset_size = 50,
#   plot_title = "Circular Dependency - Coreset Sampling Comparison"
# )
# print(circular_vis$plot)

# Example 4: Visualize all data generation processes and save results
# results <- visualize_all_data_processes(n_samples = 1000, save_plots = TRUE)
# Display specific visualization
# print(results[["Spiral Dependency"]]$plot)