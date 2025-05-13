#coreset.R
l2_leverage<- function(X,Lp){
  result<-apply(qr.Q(qr(X)),MARGIN = 1,pracma::Norm,p = Lp)
  return(result)
}

l2_prob<- function(X,Lp){
  N <- nrow(X)
  D <- ncol(X)
  result<-apply(qr.Q(qr(X)),MARGIN = 1,
                pracma::Norm,p = Lp)
  scores <- result^2
  prob <- (scores+1/N)/(sum(scores)+1)
  return(prob)
}

sample_coreset <- function(X, prob, coreset_size){
  N <- length(prob)
  prob_adj <- prob*coreset_size
  prob_adj[prob_adj>1] <- 1
  sample_index <- sort(sample(1:N,size = coreset_size,
                              prob  = prob,replace = FALSE))
  w <-  1 / (prob*coreset_size)  
  ret <- list("sample_index"=sample_index, "weights"= w,
              "prob"=prob)
  return(ret)
}

L2_coreset <- function(coreset_size,X,Lp){
  N <- nrow(X)
  D <- ncol(X)
  leverage<- l2_leverage(X=X,
                         Lp=Lp)
  scores <- leverage^2
  prob <- (scores+1/N)/(sum(scores)+1)
  sample_index <- sort(sample(1:N,size = coreset_size,
                              prob  = prob, replace = FALSE))
  w <-  1 / (prob*coreset_size)  
  ret <- list("scores"= scores, "sample_index"=sample_index, "weights"= w,
              "prob"=prob)
  return(ret)
}


Uniform_coreset <- function(coreset_size,X){
  N <- nrow(X)
  D <- ncol(X)
  sample_index <- sort(sample(1:N,size = coreset_size,
                              replace = FALSE))
  w <-  rep(1,N)
  
  ret <- list("sample_index"=sample_index, "weights"= w )
  return(ret)
  
}

# input: x - numeric matrix containing one point per row
# output: numeric matrix containing the points forming the convex hull

hull_points <- function(x){
  stopifnot(is.matrix(x))
  hull <- chull(x)
  ind <- sort(unique(as.vector(hull)))
  return(x[ind,])
}

# input: x - numeric matrix containing one point per row
# output: numeric vector containing the indices of points forming the convex hull
hull_ind <- function(x){
  stopifnot(is.matrix(x))
  hull <- chull(x)#, options = "v")
  ind <- sort(unique(as.vector(hull)))
  return(ind)
}

heuristic_chull<- function(X, M){
  n <- nrow(X)
  d <- ncol(X)
  block_size <- M/5
  chunk_size <- M/5
  u <- matrix(rnorm(d*M),ncol=M)
  bgu <- as.big.matrix(u)
  normalized_bgu <- big.matrix(nrow = nrow(bgu), 
                               ncol = ncol(bgu), init = 0, type = 'double')
  
  # block-wise read and calculate bgu
  for (j in seq(1, ncol(bgu), by = chunk_size)) {
    chunk <- bgu[,j:(j + chunk_size - 1)]
    # apply function and store the result
    temp <- chunk/apply(chunk, 2,Norm)
    # Write the results back into the big.matrix
    normalized_bgu[, j:(j + chunk_size - 1)] <- 
      temp
  }
  #Remove the current blcok
  remove(bgu)
  #gc() is used for more memory
  gc()
  x_prime_u_bg<- big.matrix(nrow = n, ncol = M)
  
  for (start_col in seq(1, M, by = block_size)) {
    end_col = min(start_col + block_size - 1, M)
    
    # Read a submatrix from the normalized bgu
    sub_matrix <- normalized_bgu[,start_col:end_col]
    
    # multiply by X
    result_sub_matrix <- X %*% sub_matrix
    
    # Store the results in xprimeu
    x_prime_u_bg[, start_col:end_col] <- result_sub_matrix
    remove(sub_matrix)
    gc()
  }
  remove(result_sub_matrix)
  remove(normalized_bgu)
  gc()
  
  # Define the block size, e.g., process 10000 columns each time
  # Loop to perform block matrix multiplication
  
  # Initialize a vector to store the index of maximum values
  max_index_list <- numeric()
  min_index_list <- numeric()
  # Assume we process n_cols columns at a time
  
  # Calculate the total number of columns
  total_cols <- ncol(x_prime_u_bg)
  
  # Process in blocks
  for(start_col in seq(1, total_cols, by = block_size)) {
    
    end_col <- min(start_col + block_size - 1, total_cols)
    
    # Use the get_sub_matrix function to obtain a submatrix (apply your own implementation here)
    sub_matrix <- x_prime_u_bg[,start_col:end_col]
    
    # Apply which.max function to the sub-matrix
    max_index_sub <- apply(sub_matrix, 2, which.max)
    min_index_sub <- apply(sub_matrix, 2, which.min)
    
    # Save the result to max_index_list
    max_index_list <- c(max_index_list, max_index_sub)
    min_index_list <- c(min_index_list, min_index_sub)
    
    remove(sub_matrix)
    gc()
  }
  remove(x_prime_u_bg)
  gc()
  unique_index <- unique(c(max_index_list,min_index_list))
  #unique_index2 <- unique(max_index_list)
  
  return(unique_index)
}


