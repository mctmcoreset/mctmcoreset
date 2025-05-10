# Coresets for Multivariate Conditional Transformation Models

This repository contains the official implementation of our NeurIPS submission "Coresets for Multivariate Conditional Transformation Models". The code is organized into two main components: simulation studies and real-world applications with high-dimensional datasets.

## Overview

Multivariate Conditional Transformation Models (MCTMs) provide a flexible framework for modeling complex dependencies between variables while maintaining interpretable marginal distributions. This repository implements efficient coreset construction methods that dramatically reduce the computational requirements while preserving model accuracy.

The repository includes both simulation experiments with various synthetic data generation processes and real-world applications on two high-dimensional datasets:
1. **Covertype Dataset**: A forest cover type dataset with 10 continuous terrain variables
2. **Financial Data**: Stock market returns for 20 long-listed stocks

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

### Simulation Part
- `main.R`: Core implementation of the MCTM model
- `complicated_functions.R`: Data generation functions for various dependency structures
- `coreset.R`: Coreset construction algorithms and utility functions
- `sampling_function.R`: Methods for coreset sampling and evaluation
- `simulation_main_function.R`: Comprehensive experiments on different data distributions
- `simulation_visualize_coreset.R`: Visualization tools for coreset sampling
- `Support_functions.R`: Helper functions for model validation and verification

### Real-World Applications
- `coreset_covertype.R`: Implementation of coreset methods for high-dimensional forest covertype data
- `covertype_descriptive_analysis.R`: Descriptive analysis and visualization of the covertype dataset
- `financial_exp.R`: Application of coreset methods to financial stock return data
- `modeling.R`: Utilities for fitting high-dimensional MCTM models

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

### Real-World Data Analysis

#### Covertype Dataset

The Covertype dataset experiments require more computational resources due to the high dimensionality. The pipeline includes:

1. **Data Preparation**: Loading and preprocessing the 10 continuous variables from the Covertype dataset
2. **Hull Precomputation**: Computing convex hull vertices for all dimensions (a key optimization)
3. **Matrix Extraction**: Extracting basis function matrices from marginal models
4. **Coreset Construction**: Using sampling methods (uniform, L2-only, L2-hull) to build coresets
5. **Evaluation**: Comparing performance metrics across coreset sizes

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

### High-Dimensional Extensions

For high-dimensional datasets (Covertype and financial data), we extend our coreset construction methods with:

1. **Parallel Computing**: Using `mclapply` for parallel computation of sensitivity scores and convex hull vertices
   ```r
   hull_indices_list <- mclapply(1:dimensions, function(j) {
     # Compute convex hull for dimension j
     hull_j_indices <- unique(as.vector(chull(A_prime_list[[j]])))
     return(hull_j_indices)
   }, mc.cores = num_cores)
   ```

2. **Probability Combination**: Advanced probability combining method for high dimensions
   ```r
   combined_probs <- combine_probabilities(prob_sources_list, n)
   ```

3. **Precomputation Strategy**: Precomputing convex hull vertices to avoid redundant calculations
   ```r
   precomputed_hull_indices <- precompute_all_hull_indices(A_prime_list, n_full)
   ```

## Datasets

### Simulation Data Generation Processes

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

### Real-World Datasets

#### Covertype Dataset

The Covertype dataset contains cartographic variables for forest cover type prediction, from which we use 10 continuous variables:

1. **Elevation**: Elevation in meters
2. **Aspect**: Aspect in degrees azimuth
3. **Slope**: Slope in degrees
4. **Horizontal_Distance_To_Hydrology**: Horizontal distance to nearest water feature
5. **Vertical_Distance_To_Hydrology**: Vertical distance to nearest water feature
6. **Horizontal_Distance_To_Roadways**: Horizontal distance to nearest roadway
7. **Hillshade_9am**: Hillshade index at 9am
8. **Hillshade_Noon**: Hillshade index at noon
9. **Hillshade_3pm**: Hillshade index at 3pm
10. **Horizontal_Distance_To_Fire_Points**: Horizontal distance to nearest wildfire ignition point

This presents a challenging high-dimensional dataset with complex dependencies between terrain variables.

#### Financial Dataset

We use daily returns of 20 long-listed stocks from major exchanges:

- NYSE stocks: JNJ, PG, KO, XOM, WMT, IBM, GE, MMM, MCD, PFE
- NASDAQ stocks: AAPL, MSFT, INTC, CSCO, AMGN, CMCSA, COST, GILD, SBUX
- Global stocks: TOT

Financial data presents unique challenges due to its volatile nature, heavy tails, and complex dependency structures that vary over time.

## Evaluation Metrics

We evaluate the performance of coreset methods using several metrics:

1. **Log-Likelihood Ratio**: Ratio of the log-likelihood achieved by the coreset model to the full data model
2. **Parameter L2 Distance**: L2 distance between parameters estimated by the coreset model and the full data model
3. **Lambda Error**: Absolute error in the dependency parameter (λ)
4. **Computation Time**: Total time for coreset construction and model fitting

For high-dimensional models with the Covertype and financial datasets, we extend these metrics to:

1. **Lambda L2 Distance**: L2 distance specifically for the dependency parameters (λ_{jl})
2. **Sampling Time**: Time required just for the coreset construction
3. **Fitting Time**: Time required for model fitting on the coreset
4. **Relative Improvement**: Percentage improvement over the uniform sampling baseline

## Reproducing the Paper Results

### Simulation Experiments

To reproduce the simulation results from our paper, run the experiments with all data generation processes:

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

### Real-World Applications

#### Covertype Dataset Experiments

The Covertype dataset contains cartographic variables used to predict forest cover type. Our experiments focus on the 10 continuous variables, such as elevation, aspect, slope, and distances to various features.

1. Download the dataset first:
```bash
wget https://archive.ics.uci.edu/ml/machine-learning-databases/covtype/covtype.data.gz
```

2. Run the Covertype coreset experiments:
```r
# Load required packages
library(tram)
library(mvtnorm)
library(ggplot2)
library(Matrix)
library(dplyr)
library(geometry)
library(parallel)
library(gridExtra)

# Run the full experiment pipeline
source("coreset_covertype.R")
```

This script loads and preprocesses the Covertype dataset, extracts the model matrices needed for coreset construction, runs the full model as a benchmark, then compares the performance of different coreset sampling methods across various coreset sizes.

#### Financial Data Experiments

The financial data experiment analyzes daily returns of 20 long-listed stocks from major exchanges.

```r
# Load required packages
library(quantmod)
library(dplyr)
library(ggplot2)
library(lubridate)
library(tram)
library(mvtnorm)
library(Matrix)
library(parallel)

# Run the financial data experiments
source("financial_exp.R")
```

This script fetches historical stock data, computes log-returns, cleans the dataset, and then applies the same coreset methods to evaluate their performance in capturing dependencies between stock returns.

## Citation

If you find our code or paper useful for your research, please consider citing:


## License

[MIT License](LICENSE)
