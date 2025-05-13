#main.R
# Main function for Multivariate Conditional Transformation Models
mmlt_continuous <- function(data, random_init = TRUE, method = "BFGS",
                            control = list(
                              trace = 1,          # Print progress
                              maxit = 1000,       # Maximum iterations
                              reltol = 1e-10      # Convergence tolerance
                            ), weights = NULL) {
  # 1. Data preparation for bivariate transformation
  x <- data[, 1]
  y <- data[, 2]
  J <- 2  # Bivariate case J=2
  
  # 2. Compute marginal Box-Cox transformations
  # Get initial estimates for \tilde{h}_j
  m_x <- tram::BoxCox(x ~ 1)
  m_y <- tram::BoxCox(y ~ 1)
  
  # 3. Get model information including basis functions
  models <- .models(m_x, m_y)
  
  # 4. Extract basis function matrices
  # Y_list[[j]] contains a_j(y_j) basis functions
  # Yprime_list[[j]] contains a'_j(y_j) basis derivatives
  Y_list <- list(
    models$mm[[1]]$eY$Y,     # a_1(y_1)
    models$mm[[2]]$eY$Y      # a_2(y_2)
  )
  Yprime_list <- list(
    models$mm[[1]]$eY$Yprime,  # a'_1(y_1)
    models$mm[[2]]$eY$Yprime   # a'_2(y_2)
  )
  
  # 5. Construct parameter names
  nx <- ncol(Y_list[[1]])  # Number of basis functions for x
  ny <- ncol(Y_list[[2]])  # Number of basis functions for y
  # Names for \theta_1 coefficients
  names_x <- paste0("x.Bs", 1:nx, "(x)")
  # Names for \theta_2 coefficients  
  names_y <- paste0("y.Bs", 1:ny, "(y)")
  # Name for \lambda_{21}
  names_lambda <- "y.x.(Intercept)"
  param_names <- c(names_x, names_y, names_lambda)
  
  # 6. Parameter initialization
  if(random_init) {
    # Generate random initial parameters ensuring monotonicity
    initial_params <- generate_initial_params(Y_list, J)
  } else {
    # Use Box-Cox estimates as initial \theta_j
    theta_init <- c(m_x$par, m_y$par)
    # Initialize \lambda_{21} = 0 (independence)
    lambda_init <- 0
    initial_params <- c(theta_init, lambda_init)
  }
  names(initial_params) <- param_names
  
  # 7. Maximum likelihood optimization
  # Minimize negative log-likelihood using BFGS
  optim_result <- optim(
    par = initial_params,
    fn = log_likelihood,  # -(\log_density + \log_det_jacobian)
    #gr = score_function,
    method = method,      # Quasi-Newton method
    # method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN",
    #            "Brent"),
    control =control,
    Y_list = Y_list,
    Yprime_list = Yprime_list,
    data = data,
    J = J
  )
  
  # 8. Organize results
  # Extract estimated parameters \hat{\theta}, \hat{\lambda}
  coef <- optim_result$par
  names(coef) <- param_names
  
  # Create result object
  result <- list(
    coefficients = coef,    # \hat{\theta}, \hat{\lambda}
    logLik = -optim_result$value,  # Maximum log-likelihood
    convergence = optim_result$convergence,  # Convergence status
    hessian = optim_result$hessian,         # For standard errors
    call = match.call(),
    models = list(m_x = m_x, m_y = m_y)     # Marginal models
  )
  
  # Set class for methods dispatch
  class(result) <- c("mmlt", "mlt")
  return(result)
}

# Compute score function based on the paper's equations (2) and (3)
score_function <- function(params, Y_list, Yprime_list, data, J) {
  n <- nrow(data)
  
  # 1. Extract parameters and construct Lambda
  theta_list <- list()
  idx <- 1
  for(j in 1:J) {
    theta_length <- ncol(Y_list[[j]])
    theta_list[[j]] <- params[idx:(idx + theta_length - 1)]
    idx <- idx + theta_length
  }
  lambda <- params[idx]
  
  # Construct Lambda as ltMatrices
  Lambda_mat <- matrix(0, nrow = J, ncol = J)
  Lambda_mat[2,1] <- lambda
  Lambda_mat <- Lambda_mat + diag(J)
  Lambda <- ltMatrices(Lambda_mat[lower.tri(Lambda_mat, diag = TRUE)], 
                       diag = TRUE)
  
  # 2. Compute transformations h(y) and h'(y)
  h_tilde <- matrix(0, nrow = n, ncol = J)
  h_prime <- matrix(0, nrow = n, ncol = J)
  
  for(j in 1:J) {
    h_tilde[,j] <- Y_list[[j]] %*% theta_list[[j]]
    h_prime[,j] <- Yprime_list[[j]] %*% theta_list[[j]]
    if(any(h_prime[,j] <= 0)) return(rep(Inf, length(params)))
  }
  
  # 3. Score for theta_j
  # ∂l/∂θ_j = Σ_i [ -h_i(z_i)λ_j + (1/h'_j(y_ij)) * ∂h'_j(y_ij)/∂θ_j ]
  score_theta <- vector("list", J)
  for(j in 1:J) {
    # Transform using Lambda
    z <- Mult(Lambda, t(h_tilde))[j,]
    
    # Score matrix for each basis coefficient
    score_j <- matrix(0, nrow = ncol(Y_list[[j]]), ncol = 1)
    
    for(k in 1:ncol(Y_list[[j]])) {
      # Normal density term
      density_term <- -sum(z * Y_list[[j]][,k])
      
      # Jacobian term
      jacobian_term <- sum(Yprime_list[[j]][,k] / h_prime[,j])
      
      score_j[k] <- density_term + jacobian_term
    }
    
    score_theta[[j]] <- score_j
  }
  
  # 4. Score for lambda
  # ∂l/∂λ_{21} = -Σ_i [ h_1(y_i) * h_2(y_i) ]
  score_lambda <- -sum(h_tilde[,1] * h_tilde[,2])
  
  # 5. Combine scores
  score <- c(unlist(score_theta), score_lambda)
  names(score) <- names(params)
  
  return(score)
}

# Multivariate Conditional Transformation Model log-likelihood function
# Based on the paper notation and equations
log_likelihood <- function(params, Y_list, Yprime_list, data, J,weights = NULL) {
  # Sample size n
  
  #Implementin weights
  n <- nrow(data)
  if(is.null(weights)){
    weights <- rep(1,n)
  }
  # 1. Extract parameters theta_j and construct parameter lists
  # \theta_j: coefficients for j-th marginal transformation
  theta_list <- list()
  idx <- 1
  for(j in 1:J) {
    theta_length <- ncol(Y_list[[j]])
    theta_list[[j]] <- params[idx:(idx + theta_length - 1)]
    idx <- idx + theta_length
  }
  
  # 2. Construct Lambda matrix
  # \Lambda = I - L - L^T where L is lower triangular with \lambda_{j\ell}
  lambda <- params[idx]  # Extract \lambda_{21} for bivariate case
  Lambda_mat <- matrix(0, nrow = J, ncol = J)
  Lambda_mat[2,1] <- lambda  # Insert \lambda_{21}
  Lambda_mat <- Lambda_mat + diag(J)  # Add identity matrix
  # Convert to ltMatrices class for efficient computation
  # Using lower triangular representation
  Lambda <- try({
    ltMatrices(Lambda_mat[lower.tri(Lambda_mat, diag = TRUE)], diag = TRUE)
  }, silent = TRUE)
  
  if(inherits(Lambda, "try-error")) return(Inf)
  
  # 3. Compute transformation values
  # h_j(y) = \tilde{h}_j(y_j) + \sum_{\ell=1}^{j-1} \lambda_{j\ell} h_\ell(y)
  h_tilde <- matrix(0, nrow = n, ncol = J)  # Store h_j(y)
  h_tilde_deriv <- matrix(0, nrow = n, ncol = J)        # Store \tilde{h}'_j(y_j)
  
  # Compute basis transformations
  # \tilde{h}_j(y_j) = a_j(y_j)^T \theta_j
  for(j in 1:J) {
    h_tilde[,j] <- Y_list[[j]] %*% theta_list[[j]]  # Basic transformation
    h_tilde_deriv[,j] <- Yprime_list[[j]] %*% theta_list[[j]]   # Its derivative
    if(any(h_tilde_deriv[,j] <= 0)) return(Inf)  # Check monotonicity
  }
  
  # 4. Transform using Lambda
  # Transform using Lambda
  z_transformed <- try(Mult(Lambda, t(h_tilde)), silent = TRUE)
  if(inherits(z_transformed, "try-error")) return(Inf)
  
  # Compute log-likelihood components
  # Jacobian term with weights
  log_det_jacobian <- sum(weights * log(h_tilde_deriv))
  # Normal density term with weights
  # Note that z_transformed is a J x n matrix, need to transpose back to n x J
  z_transformed <- t(z_transformed)  #now it is n x J
  # calculate tthe weighted density.
  log_densities <- matrix(0, nrow = n, ncol = J)
  for(j in 1:J) {
    log_densities[,j] <- dnorm(z_transformed[,j], log = TRUE)
  }
  log_density <- sum(weights * rowSums(log_densities))
  # Add diagonal term if needed
  if(attr(Lambda, "diag")) {
    effective_n <- sum(weights)
    log_density <- log_density + 
      effective_n * sum(log(diagonals(Lambda)))
  }
  # Return negative weighted log-likelihood
  return(-(log_density + log_det_jacobian))
}





# test the code
set.seed(123)
N <- 1000
sigma <- 0.4
var <- 1
Sigma <- matrix(c(var, sigma, sigma, var), nrow = 2, byrow = TRUE)
data <- as.data.frame(mvtnorm::rmvt(n = N, sigma = Sigma, df = 4))
names(data) <- c("x", "y")

# 1. Using original mmlt
x <- data[,1]
y <- data[,2]
m_x <- tram::BoxCox(x ~ 1)
m_y <- tram::BoxCox(y ~ 1)
test_model <- mmlt(m_x, m_y)
test_model$par
str(test_model$par)
# 2. using our impmentation
fit <- mmlt_continuous(data, method ="BFGS",
                       control =list(
                         trace = 1,          # Print progress
                         maxit = 1000,       # Maximum iterations
                         reltol = 1e-15      # Convergence tolerance
                       ),random_init = TRUE )

# Get model information
models <- .models(m_x, m_y)

# Get basis function matrices
Y_list <- list(
  models$mm[[1]]$eY$Y,
  models$mm[[2]]$eY$Y
)

Yprime_list <- list(
  models$mm[[1]]$eY$Yprime,
  models$mm[[2]]$eY$Yprime
)
Yprime_list[[1]]

Y_list[[1]]
# 3. Verify the calculation 
verification <- verify_calculation(test_model$par, Y_list, Yprime_list, data, J=2)
verification <- verify_calculation(fit$coefficients, Y_list, Yprime_list, data, J=2)

# 4. compare the log likelihood
cat("\n=== Results Comparison ===\n")
cat("\n1. Likelihood Values:\n")
cat("Original mmlt:", test_model$ll(test_model$par), "\n")
cat("Our implementation:", fit$logLik, "\n")
#cat("Verification:", verification$total, "\n")

cat("\n2. Parameter Estimates:\n")
cat("\nOriginal parameters:\n")
print(test_model$par)
cat("\nOur parameters:\n")
print(fit$coefficients)

cat("\n3. Convergence Information:\n")
cat("Our convergence:", fit$convergence, "\n")

# 5. Lambda and correlation coefficients
cat("\n4. Lambda and Correlation Comparison:\n")
cat("Original lambda:", test_model$par["y.x.(Intercept)"], "\n")
cat("Our lambda:", fit$coefficients["y.x.(Intercept)"], "\n")


cat("\n5. Implied Correlations:\n")
cat("Original correlation:", compute_correlation(test_model$par["y.x.(Intercept)"]), "\n")
cat("Our correlation:", compute_correlation(fit$coefficients["y.x.(Intercept)"]), "\n")
cat("True correlation used in data generation:", sigma, "\n")

cat("\n=== Testing with original parameters ===\n")
original_check <- check_likelihood_calculation(test_model$par, Y_list, Yprime_list, data, J=2)

# test the optimization process
cat("\n=== Testing with our optimized parameters ===\n")
our_check <- check_likelihood_calculation(fit$coefficients, Y_list, Yprime_list, data, J=2)

# compare results
cat("\n=== Comparison ===\n")
cat("\n1. Lambda comparison:\n")
cat("Original:", test_model$par["y.x.(Intercept)"], "\n")
cat("Ours:", fit$coefficients["y.x.(Intercept)"], "\n")

cat("\n2. Jacobian term comparison:\n")
cat("Original:", original_check$log_det_jacobian, "\n")
cat("Ours:", our_check$log_det_jacobian, "\n")

cat("\n3. Normal density term comparison:\n")
cat("Original:", original_check$log_density, "\n")
cat("Ours:", our_check$log_density, "\n")

cat("\n4. Total likelihood comparison:\n")
cat("Original (verification):", test_model$ll(test_model$par), "\n")
cat("Original (our calculation):", -original_check$total, "\n")
cat("Ours (fit$logLik):", fit$logLik, "\n")
cat("Ours (our calculation):", -our_check$total, "\n")

