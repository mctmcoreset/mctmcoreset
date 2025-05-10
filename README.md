# Coresets for Multivariate Conditional Transformation Models

This repository contains the official implementation of our NeurIPS submission "Coresets for Multivariate Conditional Transformation Models". The code is organized into two main components: simulation studies and real-world applications.

## Overview

Multivariate Conditional Transformation Models (MCTMs) provide a flexible framework for modeling complex dependencies between variables while maintaining interpretable marginal distributions. This repository implements efficient coreset construction methods that dramatically reduce the computational requirements while preserving model accuracy.

## Requirements

To run the code, you'll need R with the following packages:

```r
# Core packages
install.packages(c("tram", "mvtnorm", "colorspace", "latex2exp", 
                   "Matrix", "pracma", "geometry"))

# Additional packages for visualization
install.packages(c("tidyr", "ggplot2", "egg", "bigmemory"))

# For certain data generation processes
install.packages(c("sn", "copula"))
```

## Repository Structure

The simulation part contains the following main components:

- `main.R`: Core implementation of the MCTM model
- `complicated_functions.R`: Data generation functions for various dependency structures
- `coreset.R`: Coreset construction algorithms and utility functions
- `sampling_function.R`: Methods for coreset sampling and evaluation
- `simulation_main_function.R`: Comprehensive experiments on different data distributions
- `simulation_visualize_coreset.R`: Visualization tools for coreset sampling
- `Support_functions.R`: Helper functions for model validation and verification

## Getting Started

### Running the Simulation Studies

1. Load the required packages and utility functions:

```r
# Load required packages
source("packages.R")

# Load the implementation of the MCTM model
source("main.R")

# Load data generation functions
source("complicated_functions.R")

# Load coreset construction functions
source("coreset.R")

# Load sampling and evaluation functions
source("sampling_function.R")
```

2. Run experiments with a specific data generation process:

```r
# Example: Bivariate Normal data with 10,000 samples
results_normal <- run_multiple_simulations(
  data_gen_fn = generate_bivariate_normal,
  n_samples = 10000,
  n_sims = 3,        # Number of simulations
  min_size = 10,     # Minimum coreset size
  max_size = 500,    # Maximum coreset size
  step = 5,          # Step size for coreset sizes
  num_trials = 3     # Number of trials per simulation
)

# Analyze results
analysis_normal <- analyze_multiple_simulations(results_normal)

# Visualize results
plots_normal <- plot_multiple_simulation_results(
  analysis_results = analysis_normal,
  data_gen_name = "Bivariate Normal"
)

# Display plots
plots_normal$logLik     # Log-likelihood ratio plot
plots_normal$param_l2   # Parameter L2 distance plot
plots_normal$lambda     # Lambda parameter error plot
plots_normal$total_time # Total computation time plot

# Generate performance comparison tables
table_30_normal <- generate_performance_table(analysis_normal, 30, "Bivariate Normal")
table_100_normal <- generate_performance_table(analysis_normal, 100, "Bivariate Normal")
```

3. Visualize the coreset sampling process:

```r
# Visualize coreset sampling for spiral dependency
spiral_vis <- visualize_coreset_sampling(
  data_gen_fn = generate_spiral_dependency,
  n_samples = 1000,
  coreset_size = 50,
  plot_title = "Spiral Dependency - Coreset Sampling Comparison"
)
print(spiral_vis$plot)

# Generate visualizations for all data generation processes
all_visualizations <- visualize_all_data_processes(n_samples = 1000, save_plots = TRUE)
```

## Core Model Implementation

The main implementation of the Multivariate Conditional Transformation Model is in `main.R`. The `mmlt_continuous` function is our implementation that allows for coreset-based sampling:

```r
# Fit a MCTM model using the full dataset
fit <- mmlt_continuous(data, method = "BFGS",
                       control = list(
                         trace = 1,
                         maxit = 1000,
                         reltol = 1e-15
                       ),
                       random_init = TRUE)

# Fit a MCTM model using a coreset with weights
coreset_fit <- mmlt_continuous(
  coreset_data,
  method = "BFGS",
  control = list(
    trace = 0,
    maxit = 1000,
    reltol = 1e-15
  ),
  weights = coreset_weights,
  random_init = TRUE
)
```

## Coreset Construction Methods

We implement three main coreset sampling strategies:

1. **Uniform Sampling**: Simple random sampling without replacement
   ```r
   indices <- sort(sample(1:N, size, replace = FALSE))
   weights <- rep(N/size, length(indices))
   ```

2. **L2 Sensitivity Sampling**: Leverages the sensitivity scores of data points
   ```r
   # Sample using L2 sensitivity scores
   X_lev <- L2_coreset(size, X, 2)
   Y_lev <- L2_coreset(size, Y, 2)
   l2_indices <- sort(unique(c(X_lev$sample_index, Y_lev$sample_index)))
   ```

3. **L2-Hull Sampling**: Our proposed method combining L2 sensitivity with convex hull points
   ```r
   # Get hull indices
   hull_indices <- get_separate_hull_indices(X, Y, X_prime, Y_prime)
   
   # Sample using L2 sensitivity + convex hull points
   sample_result <- sample_l2_hull(X, Y, X_prime, Y_prime, size, hull_indices)
   ```

## Data Generation Processes

We provide implementations for various data generation processes to evaluate the robustness of our coreset methods:

- Bivariate Normal (`generate_bivariate_normal`)
- Non-linear Correlation (`generate_nonlinear_correlation`)
- Bivariate Normal Mixture (`generate_bivariate_normal_mixture`)
- Geometric Mixed Distribution (`generate_mixture`)
- Skewed and Heavy-tailed Distribution (`generate_skewt`)
- Conditional Heteroscedasticity (`generate_heteroscedastic`)
- Copula-based Complex Dependency (`generate_copula_complex`)
- Spiral Dependency (`generate_spiral_dependency`)
- Circular Dependency (`generate_circular_dependency`)
- T-Copula Dependency (`generate_t_copula_dependency`)
- Piecewise Dependency (`generate_piecewise_dependency`)
- Hourglass Dependency (`generate_hourglass_dependency`)
- Bimodal Clusters (`generate_bimodal_clusters`)
- Sinusoidal Dependency (`generate_sinusoidal_dependency`)

Each data generation process creates bivariate data with different dependency structures to test the effectiveness of our coreset methods under various conditions.

## Evaluation Metrics

We evaluate the performance of coreset methods using several metrics:

1. **Log-Likelihood Ratio**: Ratio of the log-likelihood achieved by the coreset model to the full data model
2. **Parameter L2 Distance**: L2 distance between parameters estimated by the coreset model and the full data model
3. **Lambda Error**: Absolute error in the dependency parameter (λ)
4. **Computation Time**: Total time for coreset construction and model fitting

## Reproducing the Paper Results

To reproduce the comprehensive results from our paper, run the simulation experiments with all data generation processes:

```r
# Source all required files
source("packages.R")
source("main.R")
source("complicated_functions.R")
source("coreset.R")
source("sampling_function.R")
source("simulation_main_function.R")

# Run all experiments (note: this may take significant time to complete)
source("simulation_main_function.R")
```

This will run all the experiments with the 14 different data generation processes, perform analysis, and generate plots and tables for comparison.


## License

[MIT License](LICENSE)
