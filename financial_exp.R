################################################################################
# Long-Term Listed Stock Data Acquisition and Processing
################################################################################

library(quantmod)
library(dplyr)
library(ggplot2)
library(lubridate)

# Find long-listed stocks
find_long_listed_stocks <- function(min_years = 20, exchange = NULL) {
  long_listed_stocks <- list(
    NYSE = c("JNJ", "PG", "KO", "XOM", "WMT", "IBM", "GE", "MMM", "MCD", "PFE"),
    NASDAQ = c("AAPL", "MSFT", "INTC", "CSCO", "AMGN", "ADBE", "CMCSA", "COST", "GILD", "SBUX"),
    Global = c("HSBC", "TOT", "BP", "SNY", "NVS", "TM", "HMC", "SNE", "BCS", "RDS-A")
  )
  if (!is.null(exchange) && exchange %in% names(long_listed_stocks)) {
    return(long_listed_stocks[[exchange]])
  } else {
    return(unlist(long_listed_stocks))
  }
}

# Check availability of historical stock data
check_stock_availability <- function(stocks, min_date = "2004-01-01") {
  min_date <- as.Date(min_date)
  available_stocks <- character(0)
  for (stock in stocks) {
    tryCatch({
      data <- getSymbols(stock, src = "yahoo", from = min_date, to = min_date + 30, auto.assign = FALSE)
      if (index(data)[1] <= min_date + 30) {
        available_stocks <- c(available_stocks, stock)
      }
    }, error = function(e) {})
  }
  return(available_stocks)
}

# Create full stock data with optional log-returns
create_complete_stockdat <- function(stocks, from_date = "2004-01-01", to_date = Sys.Date(), calculate_returns = TRUE) {
  logreturns <- function(prices) {
    logret <- diff(log(prices))
    return(c(NA, logret))
  }
  all_data <- list()
  date_ranges <- data.frame(stock = stocks, start_date = as.Date(NA), end_date = as.Date(NA))
  for (i in seq_along(stocks)) {
    stock <- stocks[i]
    tryCatch({
      stock_obj <- getSymbols(stock, src = "yahoo", from = from_date, to = to_date, auto.assign = FALSE)
      date_ranges[i, 2:3] <- range(index(stock_obj))
      closing_prices <- Cl(stock_obj)
      if (calculate_returns) {
        returns <- logreturns(as.numeric(closing_prices))
        data <- data.frame(date = index(closing_prices), price = as.numeric(closing_prices), return = returns)
      } else {
        data <- data.frame(date = index(closing_prices), price = as.numeric(closing_prices))
      }
      colnames(data)[2:3] <- paste0(stock, c(".price", ".return"))[1:ncol(data)-1]
      all_data[[stock]] <- data
    }, error = function(e) {})
  }
  if (length(all_data) > 0) {
    common_start <- max(date_ranges$start_date, na.rm = TRUE)
    common_end <- min(date_ranges$end_date, na.rm = TRUE)
    merged_data <- Reduce(function(x, y) merge(x, y, by = "date", all = TRUE), all_data)
    return(merged_data %>% filter(date >= common_start & date <= common_end) %>% arrange(date))
  } else {
    stop("No stock data successfully retrieved.")
  }
}

# Check completeness of price and return data
check_data_completeness <- function(data) {
  price_cols <- grep("\\.price$", colnames(data), value = TRUE)
  return_cols <- grep("\\.return$", colnames(data), value = TRUE)
  data$price_complete <- rowSums(!is.na(data[, price_cols]))
  data$return_complete <- rowSums(!is.na(data[, return_cols]))
  data$year <- year(data$date)
  return(data %>% select(-price_complete, -return_complete, -year))
}

# Keep only complete data rows
clean_stock_data <- function(data, return_only = TRUE) {
  price_cols <- grep("\\.price$", colnames(data), value = TRUE)
  return_cols <- grep("\\.return$", colnames(data), value = TRUE)
  check_cols <- if (return_only) return_cols else c(price_cols, return_cols)
  clean_data <- data %>% filter(complete.cases(select(., c("date", check_cols))))
  if (return_only) clean_data <- clean_data %>% select(date, all_of(return_cols))
  return(clean_data)
}

################################################################################
# Main Experiment Execution
################################################################################

# Step 1: Prepare data
selected_stocks <- c("JNJ", "PG", "KO", "XOM", "WMT", "IBM", "GE", "MMM", "MCD", "PFE",
                     "AAPL", "MSFT", "INTC", "CSCO", "AMGN", "CMCSA", "COST", "GILD", 
                     "SBUX", "TOT")
stock_data <- create_complete_stockdat(selected_stocks, "2015-01-01", Sys.Date())
checked_data <- check_data_completeness(stock_data)
clean_data <- clean_stock_data(checked_data)
stockdat <- clean_data %>% select(-date)
dim(stockdat)
# Step 2: Define parameters
MIN_K <- 20
MAX_K <- 300
STEP_K <- 10
NUM_TRIALS <- 5
SEED <- 42
FORMULA <- ~ 1
ORDER <- 6
K_HULL_PROP <- 0.05
DIMENSIONS <- ncol(stockdat)

# Step 3: Fit full MCTM model
colnames(stockdat) <- paste0("y", 1:ncol(stockdat))
full_mctm_model <- fit_full_mctm(stockdat, formula = FORMULA, order = ORDER)
full_mctm_model$coef

# Step 4: Extract matrices
mctm_matrices <- extract_mctm_matrices(stockdat, formula = FORMULA, order = ORDER)

# Step 5: Run coreset comparison
comparison_results <- run_covertype_comparison(
  data = stockdat,
  matrices = mctm_matrices,
  full_model = full_mctm_model,
  min_k = MIN_K,
  max_k = MAX_K,
  step_k = STEP_K,
  num_trials = NUM_TRIALS,
  k_hull_prop = K_HULL_PROP,
  seed = SEED
)

# Step 6: Analyze and visualize
analysis_summary <- analyze_covertype_results(comparison_results)

plots <- plot_covertype_results(analysis_summary, 
                                data_gen_name = "Stock returns 20D (Unconditional)",
                                ylim_logLik = c(1, -3),
                                xlim = c(40, 300),
                                ylim_total_time = c(0, 100)
)
print(plots$logLik)
print(plots$param_l2)
print(plots$lambda_l2)
print(plots$total_time)

# Step 7: Generate performance tables
performance_table_50 <- generate_covertype_table(analysis_summary, k_cutoff = 50)
performance_table_100 <- generate_covertype_table(analysis_summary, k_cutoff = 100)
performance_table_200 <- generate_covertype_table(analysis_summary, k_cutoff = 200)
performance_table_500 <- generate_covertype_table(analysis_summary, k_cutoff = 500)
