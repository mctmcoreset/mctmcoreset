#complicated_functions.R

# 1. Bivariate normal
generate_bivariate_normal <- function(n = 1000, mu = c(0, 0), Sigma = NULL, rho = 0.7) {
  # If Sigma is not provided, construct it using rho
  if (is.null(Sigma)) {
    # Default standard deviations
    sigma1 <- 1
    sigma2 <- 1
    
    # Construct covariance matrix from correlation
    Sigma <- matrix(c(
      sigma1^2, rho * sigma1 * sigma2,
      rho * sigma1 * sigma2, sigma2^2
    ), nrow = 2)
  }
  
  # Generate data using mvrnorm from the MASS package
  data <- MASS::mvrnorm(n = n, mu = mu, Sigma = Sigma)
  
  # Convert to data frame with column names y1 and y2 to match other functions
  df <- as.data.frame(data)
  colnames(df) <- c("y1", "y2")
  
  return(df)
}
# 2. Non-linear correlation
generate_nonlinear_correlation <- function(n = 1000) {
  # Generate base variable
  x <- seq(-3, 3, length.out = n)
  
  # Correlation coefficient changes with x
  rho <- sin(x)  # Or rho = 2 * pnorm(x) - 1
  
  # Generate correlated random variables
  y1 <- x^2 + rnorm(n, sd = 0.5)
  y2 <- numeric(n)
  
  for(i in 1:n) {
    sigma <- matrix(c(1, rho[i], rho[i], 1), 2, 2)
    temp <- MASS::mvrnorm(1, mu = c(0, 0), Sigma = sigma)
    y2[i] <- temp[2]
  }
  
  return(data.frame(y1 = y1, y2 = y2))
}

# 3. Bivariate normal mixture 
generate_mixture <- function(n = 1000) {
  # Required packages
  if(!require(MASS)) {
    stop("MASS package is required for mvrnorm function")
  }
  
  # Mixture weights
  p1 <- 0.5  # 50% from first component
  p2 <- 0.5  # 50% from second component
  
  # Determine sample sizes
  n1 <- round(p1 * n)
  n2 <- n - n1  # Ensure total is exactly n
  
  # Component 1: Circular pattern with noise
  theta1 <- runif(n1, 0, 2*pi)
  r1 <- 2 + rnorm(n1, 0, 0.2)  # Small radius with little noise
  x1 <- r1 * cos(theta1)
  y1 <- r1 * sin(theta1)
  samples1 <- cbind(x1, y1)
  
  # Component 2: Cross-shaped pattern
  # Horizontal part of cross
  n2a <- round(n2/2)
  n2b <- n2 - n2a
  
  x2a <- runif(n2a, -5, 5)
  y2a <- rnorm(n2a, 0, 0.3)
  
  # Vertical part of cross
  x2b <- rnorm(n2b, 0, 0.3)
  y2b <- runif(n2b, -5, 5)
  
  samples2 <- rbind(
    cbind(x2a, y2a),
    cbind(x2b, y2b)
  )
  
  # Combine samples
  combined_samples <- rbind(samples1, samples2)
  
  # Randomly shuffle the rows to mix the samples
  shuffled_indices <- sample(1:n, n)
  shuffled_samples <- combined_samples[shuffled_indices, ]
  
  # Return as a data frame with column names
  return(data.frame(y1 = shuffled_samples[, 1], 
                    y2 = shuffled_samples[, 2]))
}

generate_bivariate_normal_mixture <- function(n = 1000) {
  # Required packages
  if(!require(MASS)) {
    stop("MASS package is required for mvrnorm function")
  }
  
  # Parameters for the first bivariate normal distribution
  mu1 <- c(-2, 2)  # Mean vector for first component
  Sigma1 <- matrix(c(1.2, 0.6, 
                     0.6, 0.8), 2, 2)  # Covariance matrix with positive correlation
  
  # Parameters for the second bivariate normal distribution
  mu2 <- c(3, -1)  # Mean vector for second component
  Sigma2 <- matrix(c(0.9, -0.4, 
                     -0.4, 1.1), 2, 2)  # Covariance matrix with negative correlation
  
  # Mixture weights (proportion of samples from each distribution)
  p1 <- 0.4  # 40% from distribution 1
  p2 <- 0.6  # 60% from distribution 2
  
  # Determine how many samples to draw from each distribution
  n1 <- round(p1 * n)
  n2 <- n - n1  # Ensure total is exactly n
  
  # Generate samples from each distribution
  samples1 <- MASS::mvrnorm(n1, mu = mu1, Sigma = Sigma1)
  samples2 <- MASS::mvrnorm(n2, mu = mu2, Sigma = Sigma2)
  
  # Combine samples
  combined_samples <- rbind(samples1, samples2)
  
  # Randomly shuffle the rows to mix the samples
  shuffled_indices <- sample(1:n, n)
  shuffled_samples <- combined_samples[shuffled_indices, ]
  
  # Return as a data frame with column names
  return(data.frame(y1 = shuffled_samples[, 1], 
                    y2 = shuffled_samples[, 2]))
}


# 5. Skewed and heavy-tailed distribution
generate_skewt <- function(n = 1000) {
  require(sn)  # Requires sn package
  
  # Define skewed t-distribution parameters
  xi <- c(0, 0)       # Location parameter
  Omega <- matrix(c(1, 0.5, 0.5, 1), 2, 2)  # Scale matrix
  alpha <- c(5, -3)   # Skewness parameter
  nu <- 4            # Degrees of freedom
  
  # Generate skewed t-distribution samples
  y <- sn::rmst(n, xi = xi, Omega = Omega, alpha = alpha, nu = nu)
  
  return(data.frame(y1 = y[,1], y2 = y[,2]))
}

# 6. Conditional heteroscedasticity
generate_heteroscedastic <- function(n = 1000) {
  x <- seq(-3, 3, length.out = n)
  
  # Variance changes with x
  sigma1 <- exp(0.5 * x)
  sigma2 <- sqrt(abs(x))
  
  y1 <- rnorm(n, mean = x^2, sd = sigma1)
  y2 <- rnorm(n, mean = sin(x), sd = sigma2)
  
  return(data.frame(y1 = y1, y2 = y2))
}

# 7. Copula-based complex dependency structure
generate_copula_complex <- function(n = 1000) {
  require(copula)
  
  # Create a Clayton copula
  clayton_cop <- claytonCopula(param = 2)
  
  # Generate copula samples
  u <- rCopula(n, clayton_cop)
  
  # Transform to target marginal distributions
  y1 <- qgamma(u[,1], shape = 2, rate = 1)
  y2 <- qlnorm(u[,2], meanlog = 0, sdlog = 1)
  
  return(data.frame(y1 = y1, y2 = y2))
}

# 8. Spiral Pattern - creates data with a spiral dependency structure
generate_spiral_dependency <- function(n = 1000, noise_level = 0.5) {
  # Generate angle parameter
  t <- seq(0, 3*pi, length.out = n)
  
  # Create spiral pattern with noise
  r <- 0.5 * t
  x <- r * cos(t) + rnorm(n, sd = noise_level)
  y <- r * sin(t) + rnorm(n, sd = noise_level)
  
  return(data.frame(y1 = x, y2 = y))
}

# 9. Circle/Ring Pattern - creates data with a circular dependency
generate_circular_dependency <- function(n = 1000, radius = 5, thickness = 1) {
  # Generate angles
  theta <- runif(n, 0, 2*pi)
  
  # Generate distances from center (with variation for thickness)
  r <- radius + rnorm(n, 0, thickness)
  
  # Convert to cartesian coordinates
  x <- r * cos(theta)
  y <- r * sin(theta)
  
  return(data.frame(y1 = x, y2 = y))
}

# 10. Copula-Based Dependency - creates data with specific correlation structure
generate_t_copula_dependency <- function(n = 1000, rho = 0.7, df = 3) {
  # Required packages
  if(!require(copula)) {
    stop("copula package is required")
  }
  
  # Create t-copula with specified correlation and degrees of freedom
  tcop <- tCopula(param = rho, dim = 2, df = df)
  
  # Generate data from copula
  u <- rCopula(n, tcop)
  
  # Transform to desired marginal distributions
  # Using t distribution for both margins for heavier tails
  x <- qt(u[,1], df = 5) # First variable: t with 5 df
  y <- qexp(u[,2])       # Second variable: exponential
  
  return(data.frame(y1 = x, y2 = y))
}

# 11. Piecewise Linear Dependency - creates data with different correlations in different regions
generate_piecewise_dependency <- function(n = 1000) {
  # Generate base x values across the range
  x <- rnorm(n, sd = 2)
  
  # Create y with different dependency based on regions of x
  y <- numeric(n)
  
  # Implementing different correlation structures in different regions
  for (i in 1:n) {
    if (x[i] < -1) {
      # Strong positive correlation in left region
      y[i] <- 1.5 * x[i] + rnorm(1, sd = 0.5)
    } else if (x[i] < 1) {
      # Weak negative correlation in middle region
      y[i] <- -0.5 * x[i] + rnorm(1, sd = 0.8)
    } else {
      # Strong negative correlation in right region
      y[i] <- -2 * x[i] + rnorm(1, sd = 0.5)
    }
  }
  
  return(data.frame(y1 = x, y2 = y))
}

# 12. Hourglass Pattern - creates data with non-homogeneous variance structure
generate_hourglass_dependency <- function(n = 1000) {
  # Generate base x values
  x <- rnorm(n, sd = 2)
  
  # Create y with variance dependent on x (higher at extremes)
  variance <- 0.2 + 0.3 * x^2  # Variance increases quadratically away from zero
  y <- rnorm(n, mean = 0, sd = sqrt(variance))
  
  return(data.frame(y1 = x, y2 = y))
}

# 13. Bimodal Clusters with Correlation - two clusters with different correlations
generate_bimodal_clusters <- function(n = 1000) {
  # Determine sizes of each cluster
  n1 <- round(n/2)
  n2 <- n - n1
  
  # Create covariance matrices for each cluster
  Sigma1 <- matrix(c(1, 0.8, 
                     0.8, 1), 2, 2)  # Strong positive correlation
  
  Sigma2 <- matrix(c(1, -0.7, 
                     -0.7, 1), 2, 2)  # Strong negative correlation
  
  # Generate samples from each cluster
  cluster1 <- MASS::mvrnorm(n1, mu = c(-2, 2), Sigma = Sigma1)
  cluster2 <- MASS::mvrnorm(n2, mu = c(2, 2), Sigma = Sigma2)
  
  # Combine and shuffle
  data <- rbind(cluster1, cluster2)
  data <- data[sample(1:n, n),]
  
  return(data.frame(y1 = data[,1], y2 = data[,2]))
}

# 14. Sinusoidal Dependency - sine wave pattern
generate_sinusoidal_dependency <- function(n = 1000, amplitude = 2, frequency = 1, noise_sd = 0.5) {
  # Generate evenly spaced x values
  x <- seq(-3, 3, length.out = n)
  
  # Generate y with sinusoidal pattern plus noise
  y <- amplitude * sin(frequency * pi * x) + rnorm(n, sd = noise_sd)
  
  return(data.frame(y1 = x, y2 = y))
}

# Demonstration visualization function
visualize_all_dependencies <- function() {
  if(!require(ggplot2)) {
    stop("ggplot2 package is required")
  }
  
  # Set seed for reproducibility
  set.seed(123)
  
  # Generate data from all methods
  data_list <- list(
    Spiral = generate_spiral_dependency(500),
    Circular = generate_circular_dependency(500),
    TCopula = generate_t_copula_dependency(500),
    Piecewise = generate_piecewise_dependency(500),
    Hourglass = generate_hourglass_dependency(500),
    BimodalClusters = generate_bimodal_clusters(500),
    Sinusoidal = generate_sinusoidal_dependency(500)
  )
  
  # Create plots for each dataset
  plots <- list()
  for (name in names(data_list)) {
    data <- data_list[[name]]
    plots[[name]] <- ggplot(data, aes(x = y1, y = y2)) +
      geom_point(alpha = 0.5, size = 0.8) +
      geom_density_2d(color = "blue", alpha = 0.8) +
      ggtitle(paste("Dependency Pattern:", name)) +
      theme_minimal()
  }
  
  return(plots)
}

