mmlt <- function (..., formula = ~1, data, conditional = FALSE, theta = NULL, 
                  fixed = NULL, scale = FALSE, optim = mmltoptim(), 
                  args = list(seed = 1, M = 1000), dofit = TRUE, domargins = TRUE) {
  
  # Record the function call to allow reproduction of the calling environment when necessary.
  call <- match.call() 
  
  # Use the .models function to define and initialize the model structure.
  m <- .models(...) 
  
  # Check if conditional models are allowed; conditional models are only applicable to probit-type marginal models.
  if (conditional && !all(m$normal)) 
    stop("Conditional models only available for marginal probit-type models.")
  
  # Check if conditional models can be fitted along with marginal model parameters.
  if (conditional && !domargins) 
    stop("Conditional models must fit marginal and joint parameters.")
  
  # If theta is NULL and fitting is allowed along with marginal parameter fitting.
  if (is.null(theta) && dofit && domargins) {
    cl <- match.call()
    cl$conditional <- FALSE  # Disable conditional models
    cl$domargins <- FALSE    # Disable marginal parameter fitting
    sm <- eval(cl, parent.frame())  # Re-evaluate the modified calling environment
    
    # If sm is not NULL, extract its parameters
    if (!is.null(sm)) {
      theta <- coef(sm, type = "all")  # Get all parameters
      
      # If it's a conditional model, adjust theta
      if (conditional) {
        class(sm)[1] <- "cmmlt"  # Set model class to conditional model
        d <- rowMeans(diagonals(coef(sm, newdata = data, type = "Sigma")))  # Calculate the mean of the diagonal elements of the covariance matrix
        theta[1:sum(m$nparm)] <- theta[1:sum(m$nparm)] * rep(sqrt(d), times = m$nparm)  # Adjust theta
      }
    } else {
      # If sm is NULL, use the coefficients of each model
      theta <- do.call("c", lapply(m$models, function(mod) coef(mod)))
    }
  }
  
  # Calculate the number of continuous and discrete variables
  cJ <- sum(m$cont)  # Number of continuous variables
  dJ <- sum(!m$cont)  # Number of discrete variables
  J <- cJ + dJ  # Total number of variables
  Jp <- J * (J - 1) / 2  # Number of off-diagonal elements
  
  # Initialize the log-likelihood and gradient computation functions
  llsc <- .ll(c(cJ, dJ), standardize = !conditional, args)
  
  # If there are discrete variables and weights are not provided, initialize Monte Carlo weights
  if (dJ && is.null(args$w)) 
    args$w <- .MCw(J = dJ, M = args$M, seed = args$seed)
  
  # Check if the formula only contains an intercept
  if (isTRUE(all.equal(formula, ~1))) {
    lX <- matrix(1)  # Create a design matrix with only an intercept
    colnames(lX) <- "(Intercept)"
    bx <- NULL  # No basis functions
  } else {
    bx <- formula  # Use the provided formula
    if (inherits(formula, "formula")) {
      bx <- as.basis(formula, data)  # Convert the formula to basis functions
    }
    lX <- model.matrix(bx, data = data)  # Generate the model matrix
    if (conditional) 
      warning("Conditional models with covariate-dependent correlations are order-dependent")  # Warn that covariate correlations are order-dependent
  }
  
  # Internal function to convert parameters to a parameter matrix
  .Xparm <- function(parm) {
    parm <- parm[-(1:sum(m$nparm))]  # Remove marginal model parameters
    return(matrix(parm, nrow = ncol(lX)))  # Convert remaining parameters to a matrix
  }
  
  # Initialize starting parameters
  start <- .start(m, colnames(lX), names(fixed))
  
  # Set parameter names
  parnames <- eparnames <- names(start)  # All parameter names
  lparnames <- names(start)[-(1:sum(m$nparm))]  # Non-marginal parameter names
  
  # If there are fixed parameters, update parameter names
  if (!is.null(fixed)) 
    eparnames <- eparnames[!eparnames %in% names(fixed)]
  if (!is.null(fixed)) 
    lparnames <- lparnames[!lparnames %in% names(fixed)]
  
  # If theta is provided, verify its length and set parameter names
  if (!is.null(theta)) {
    if (!is.null(fixed)) 
      theta <- theta[!names(theta) %in% names(fixed)]
    stopifnot(length(theta) == length(eparnames))  # Ensure theta has correct length
    names(theta) <- eparnames  # Set names for theta
  }
  
  # If there are continuous variables, generate model matrices
  if (cJ) {
    mm <- lapply(1:cJ, function(j) .model_matrix(m, j = j))  # Model matrices for continuous variables
    mmp <- lapply(1:cJ, function(j) .model_matrix(m, j = j, prime = TRUE))  # First derivative matrices for continuous variables
  }
  
  # If there are discrete variables, generate model matrices and handle infinite values
  if (dJ) {
    dmm <- lapply(cJ + 1:dJ, function(j) .model_matrix(m, j = j))  # Model matrices for discrete variables
    dmm <- lapply(dmm, function(x) {
      x$Yleft[!is.finite(x$Yleft[, 1]), ] <- 0  # Handle infinite values on the left boundary
      x$Yright[!is.finite(x$Yright[, 1]), ] <- 0  # Handle infinite values on the right boundary
      x
    })
  }
  
  # Initialize weights
  weights <- m$weights  # Weights provided by the model
  # if (weights) {
  #   if (cJ) {
  #     mm <- lapply(mm, function(x) x * weights)  # Apply weights to continuous variables
  #     mmp <- lapply(mmp, function(x) x * weights)  # Apply weights to continuous variable derivatives
  #   }
  #   if (dJ) {
  #     dmm <- lapply(dmm, function(x) list(Yleft = x$Yleft * weights, Yright = x$Yright * weights))  # Apply weights to discrete variables
  #   }
  # }
  
  # Create the initial Lambda matrix, initialized to zero
  LAMBDA <- ltMatrices(matrix(0, nrow = Jp, ncol = nrow(lX)), byrow = TRUE, diag = FALSE, names = names(m$models))
  
  # Define the log-likelihood function
  ll <- function(parm, newdata = NULL) {
    # If new data is provided and the formula is not constant
    if (!is.null(newdata) && !isTRUE(all.equal(formula, ~1))) 
      lX <- model.matrix(bx, data = newdata)  # Generate the model matrix for new data
    
    Lambda <- LAMBDA  # Initialize the Lambda matrix
    Lambda[] <- t(lX %*% .Xparm(parm))  # Update the Lambda matrix
    
    ret <- 0  # Initialize the return value
    
    # If there are continuous variables, calculate their log-likelihood
    if (cJ) {
      z <- .rbind(.mget(m, j = which(m$cont), 
                        parm = parm, what = "z", 
                        newdata = newdata))  # Extract z
      zp <- .rbind(.mget(m, j = which(m$cont), 
                         parm = parm, what = "zprime", 
                         newdata = newdata))  # Extract zprime
      ret <- colSums(.log(zp))  # Calculate log
      if (!dJ) 
        return(ret + llsc$logLik(obs = z, Lambda = Lambda))  # Return log-likelihood value
    }
    
    # If there are discrete variables, calculate their log-likelihood
    if (dJ) {
      lower <- .rbind(.mget(m,
                            j = which(!m$cont),
                            parm = parm, what = "zleft", 
                            newdata = newdata))  # Extract lower bounds
      upper <- .rbind(.mget(m, j = which(!m$cont), 
                            parm = parm, what = "zright", 
                            newdata = newdata))  # Extract upper bounds
      if (!cJ) 
        return(llsc$logLik(lower = lower, upper = upper, Lambda = Lambda))  # Return log-likelihood value
    }
    
    # Return the cumulative log-likelihood value
    return(ret + llsc$logLik(obs = z, lower = lower, upper = upper, Lambda = Lambda))
  }
  
  # Define the score function (gradient)
  sc <- function(parm, newdata = NULL, scores = FALSE) {
    # If individual scores are requested
    if (scores) {
      RS <- CS <- function(x) x  # Return original values
    } else {
      RS <- function(x) rowSums(x, na.rm = TRUE)  # Row sums
      CS <- function(x) colSums(x, na.rm = TRUE)  # Column sums
    }
    
    # Handle the model matrix for new data
    if (!is.null(newdata) && !isTRUE(all.equal(formula, ~1))) 
      lX <- model.matrix(bx, data = newdata)  # Generate the model matrix for new data
    
    Lambda <- LAMBDA  # Initialize the Lambda matrix
    Lambda[] <- t(lX %*% .Xparm(parm))  # Update the Lambda matrix
    
    if (cJ) {
      z <- .rbind(.mget(m, j = which(m$cont), parm = parm, what = "z", newdata = newdata))  # Extract z
      if (!dJ) 
        sc <- llsc$score(obs = z, Lambda = Lambda)  # Calculate scores for continuous variables
    }
    
    if (dJ) {
      lower <- .rbind(.mget(m, j = which(!m$cont), parm = parm, what = "zleft", newdata = newdata))  # Extract lower bounds
      upper <- .rbind(.mget(m, j = which(!m$cont), parm = parm, what = "zright", newdata = newdata))  # Extract upper bounds
      if (!cJ) 
        sc <- llsc$score(lower = lower, upper = upper, Lambda = Lambda)  # Calculate scores for discrete variables
    }
    
    if (cJ && dJ) 
      sc <- llsc$score(obs = z, lower = lower, upper = upper, Lambda = Lambda)  # Combine scores for continuous and discrete variables
    
    # Handle the lower triangular part of the Lambda matrix
    Lmat <- Lower_tri(sc$Lambda)[rep(1:Jp, each = ncol(lX)), , drop = FALSE]
    if (identical(c(lX), 1)) {
      scL <- RS(Lmat)  # Row sums
    } else {
      scL <- RS(Lmat * t(lX[, rep(1:ncol(lX), Jp), drop = FALSE]))  # Matrix multiplication followed by row sums
    }
    
    scp <- vector(mode = "list", length = cJ + dJ)  # Initialize the score vector
    
    if (cJ) {
      if (all(m$normal)) {
        zp <- .rbind(.mget(m, j = which(m$cont), parm = parm, what = "zprime", newdata = newdata))  # Extract zprime
        scp[1:cJ] <- lapply(1:cJ, function(j) {
          CS(mm[[j]] * c(sc$obs[j, ])) + CS(mmp[[j]] / c(zp[j, ]))  # Calculate scores for continuous variables
        })
      } else {
        dz <- .rbind(.mget(m, j = which(m$cont), parm = parm, what = "dtrafo", newdata = newdata))  # Extract derivatives
        ef <- lapply(which(m$cont), function(j) .mget(m, j = j, parm = parm, what = "estfun", newdata = newdata))  # Extract estimating functions
        scp[1:cJ] <- lapply(1:cJ, function(j) {
          CS(mm[[j]] * c(sc$obs[j, ] + z[j, ]) / c(dnorm(z[j, ])) * c(dz[j, ])) - CS(ef[[j]])  # Calculate scores for non-normal continuous variables
        })
      }
    }
    
    if (dJ) {
      if (all(m$normal)) {
        scp[cJ + 1:dJ] <- lapply(1:dJ, function(j) {
          CS(dmm[[j]]$Yleft * c(sc$lower[j, ])) + CS(dmm[[j]]$Yright * c(sc$upper[j, ]))  # Calculate scores for discrete variables
        })
      } else {
        dzl <- .rbind(.mget(m, j = which(!m$cont), parm = parm, what = "dzleft", newdata = newdata))  # Extract left derivatives
        dzl[!is.finite(dzl)] <- 0  # Handle infinite values
        dzr <- .rbind(.mget(m, j = which(!m$cont), parm = parm, what = "dzright", newdata = newdata))  # Extract right derivatives
        dzr[!is.finite(dzr)] <- 0  # Handle infinite values
        scp[cJ + 1:dJ] <- lapply(1:dJ, function(j) {
          return(CS(dmm[[j]]$Yleft * c(dzl[j, ]) * c(sc$lower[j, ])) + CS(dmm[[j]]$Yright * c(dzr[j, ]) * c(sc$upper[j, ])))  # Calculate scores for non-normal discrete variables
        })
      }
    }
    
    if (!scores) {
      ret <- c(do.call("c", scp), c(scL))  # Combine scores
      names(ret) <- parnames
      return(ret)
    }
    
    ret <- cbind(do.call("cbind", scp), t(scL))  # Combine score matrices
    colnames(ret) <- parnames
    return(ret)
  }
  
  # If scale is TRUE, perform scaling
  if (scale) {
    scl <- rep(apply(abs(lX), 2, max, na.rm = TRUE), times = Jp)  # Calculate the maximum absolute value for each column
    lt1 <- scl < 1.1  # Flag values less than 1.1
    gt1 <- scl >= 1.1  # Flag values greater than or equal to 1.1
    scl[gt1] <- 1/scl[gt1]  # Take the reciprocal for values greater than or equal to 1.1
    scl[lt1] <- 1  # Set values less than 1.1 to 1
    scl <- c(do.call("c", .mget(m, j = 1:J, parm = NULL, what = "scale")), scl)  # Combine scaling parameters
    names(scl) <- parnames  # Set parameter names
  } else {
    # If not scaling, set all scaling values to 1
    scl <- numeric(length(parnames))
    scl[] <- 1
    names(scl) <- parnames
  }
  
  # Initialize weights; if not defined, default to 1
  if (!weights) 
    weights <- 1
  
  # Define the objective function (negative log-likelihood)
  f <- function(par, scl, ...) {
    if (!is.null(fixed)) {
      p <- par
      names(p) <- eparnames
      p[names(fixed)] <- fixed  # Keep fixed parameters unchanged
      par <- p[parnames]
    }
    return(-sum(weights * ll(par * scl, ...)))  # Return negative log-likelihood
  }
  
  # Define the gradient function (negative scores)
  g <- function(par, scl, ...) {
    if (!is.null(fixed)) {
      p <- par
      names(p) <- eparnames
      p[names(fixed)] <- fixed  # Keep fixed parameters unchanged
      par <- p[parnames]
    }
    ret <- -sc(par * scl, ...) * scl  # Calculate negative gradient
    if (is.null(fixed)) 
      return(ret)
    if (is.matrix(ret)) 
      return(ret[, !parnames %in% names(fixed)])  # Return gradient for non-fixed parameters
    return(ret[!parnames %in% names(fixed)])
  }
  
  # If not fitting marginal parameters
  if (!domargins) {
    stopifnot(!conditional)  # Ensure it's not a conditional model
    start <- start[1:sum(m$nparm)]  # Retain only marginal parameters
    
    # Define the sub-log-likelihood function
    cll <- function(cpar) f(c(start/scl[names(start)], cpar), scl = scl)
    
    # Define the sub-gradient function
    csc <- function(cpar) {
      ret <- g(c(start/scl[names(start)], cpar), scl = scl)
      return(ret[names(lambdastart)])
    }
    
    # If theta is not defined, initialize lambda's starting values
    if (is.null(theta)) {
      lambdastart <- rep(0, length(lparnames))
      names(lambdastart) <- lparnames
    } else {
      lambdastart <- theta[lparnames]
    }
    
    # Optimize the lambda parameters
    if (length(lambdastart)) {
      for (i in 1:length(optim)) {
        op <- optim[[i]](theta = lambdastart, f = cll, g = csc)  # Optimize
        if (op$convergence == 0) 
          break  # Exit loop if converged
      }
      names(op$par) <- names(lambdastart)
      ret <- c(start, op$par * scl[names(lambdastart)])  # Combine results
      ret <- list(par = ret, value = -op$value)  # Return optimization results
    } else {
      return(NULL)  # Return NULL if lambdastart is not defined
    }
  } else {
    # If fitting marginal parameters
    ui <- m$ui  # Linear constraint matrix
    ui <- cbind(ui, matrix(0, nrow = nrow(ui), ncol = length(parnames) - ncol(ui)))
    if (!is.null(fixed)) 
      ui <- ui[, !parnames %in% names(fixed), drop = FALSE]  # Remove fixed parameters
    
    ci <- m$ci  # Right-hand side of linear constraints
    
    if (is.null(theta) && !dofit) 
      return(list(ll = function(...) f(..., scl = 1), score = function(...) g(..., scl = 1), ui = ui, ci = ci))  # Return functions and constraints when not fitting
    
    start <- theta/scl[eparnames]  # Initial parameters
    ui <- t(t(ui) * scl[eparnames])  # Scale the constraint matrix
    
    if (dofit) {
      for (i in 1:length(optim)) {
        ret <- optim[[i]](theta = start, f = function(par) f(par, scl = scl), g = function(par) g(par, scl = scl), ui = ui, ci = ci)  # Optimize
        if (ret$convergence == 0) 
          break  # Exit loop if converged
      }
      if (ret$convergence != 0) 
        warning("Optimisation did not converge")  # Warn if not converged
    } else {
      ret <- list(par = start, value = f(theta, scl = 1), convergence = NA, optim_hessian = NA)  # Return initial values
    }
    
    names(ret$par) <- eparnames
    ret$par[eparnames] <- ret$par[eparnames] * scl[eparnames]  # Scale parameters
  }
  ret$ll <- function(...) f(..., scl = 1)  # Define the returned log-likelihood function
  ret$score <- function(...) g(..., scl = 1)  # Define the returned gradient function
  ret$args <- args  # Save arguments
  ret$logLik <- -ret$value  # Save log-likelihood value
  ret$models <- m  # Save model
  ret$formula <- formula  # Save formula
  ret$bx <- bx  # Save basis functions
  ret$parm <- function(par) {
    if (!is.null(fixed)) 
      par <- c(par, fixed)[parnames]
    return(c(m$parm(par), list(.Xparm(par))))  # Return parameters and basis functions
  }
  if (!missing(data)) 
    ret$data <- data  # Save data
  ret$names <- m$names  # Save model names
  ret$call <- match.call()  # Save the call
  class(ret) <- c(ifelse(conditional, "cmmlt", "mmmlt"), "mmlt")  # Set class names
  ret$mmlt <- "Multivariate Conditional Transformation Model"  # Set model description
  ret
}

.ll <- function(dim, standardize = TRUE, args = list()) {
  # If dim has only one element, expand it to two elements (default the second to 0)
  if (length(dim) == 1L) 
    dim <- c(dim, 0L)
  
  cJ <- dim[1L]  # Number of continuous variables
  dJ <- dim[2L]  # Number of discrete variables
  
  # If there are no discrete variables
  if (!dJ) {
    # Define the log-likelihood function for continuous variables
    cll <- function(obs, Lambda) {
      # Check if Lambda is a non-diagonal matrix
      if (dim(Lambda)[2L] > 1) 
        stopifnot(!attr(Lambda, "diag"))
      if (!standardize) 
        return(ldmvnorm(obs = obs, 
                        invchol = Lambda, logLik = FALSE))  # Calculate non-standardized log-likelihood
      sLambda <- mvtnorm::standardize(invchol = Lambda)  # Standardize Lambda
      return(ldmvnorm(obs = obs, 
                      invchol = sLambda, logLik = FALSE))  # Calculate standardized log-likelihood
    }
    
    # Define the score function for continuous variables
    csc <- function(obs, Lambda) {
      if (dim(Lambda)[2L] > 1) 
        stopifnot(!attr(Lambda, "diag"))
      if (!standardize) {
        ret <- sldmvnorm(obs = obs, invchol = Lambda)  # Calculate score
        return(list(Lambda = ret$invchol, obs = ret$obs))  # Return score and transformed obs
      }
      chol <- solve(Lambda)  # Calculate the inverse matrix of Lambda
      D <- sqrt(Tcrossprod(chol, diag_only = TRUE))  # Calculate the standardization factor
      sLambda <- invcholD(Lambda, D = D)  # Standardize Lambda
      ret <- sldmvnorm(obs = obs, invchol = sLambda)  # Calculate standardized score
      ret$chol <- -vectrick(sLambda, ret$invchol)  # Calculate the correction term for the score
      dobs <- ret$obs
      ret <- destandardize(chol = chol, invchol = Lambda, score_schol = ret$chol)  # Destandardize
      return(list(Lambda = ret, obs = dobs))  # Return score and transformed obs
    }
    
    return(list(logLik = cll, score = csc))  # Return log-likelihood and score functions
  }
  
  # Define the log-likelihood function with discrete variables
  ll <- function(obs = NULL, lower, upper, Lambda) {
    if (dim(Lambda)[2L] > 1) 
      stopifnot(!attr(Lambda, "diag"))
    a <- args
    a$obs <- obs  # Set obs
    a$mean <- 0  # Set mean to 0
    a$lower <- lower  # Set lower bounds
    a$upper <- upper  # Set upper bounds
    a$logLik <- FALSE  # Disable log-likelihood flag
    if (!standardize) {
      a$invchol <- Lambda  # Use non-standardized Lambda
    } else {
      a$chol <- mvtnorm::standardize(chol = solve(Lambda))  # Use standardized Lambda
    }
    return(do.call("ldpmvnorm", a))  # Calculate log-likelihood of the multivariate normal distribution
  }
  
  # Define the score function with discrete variables
  sc <- function(obs = NULL, lower, upper, Lambda) {
    a <- args
    a$obs <- obs  # Set obs
    a$mean <- 0  # Set mean to 0
    a$lower <- lower  # Set lower bounds
    a$upper <- upper  # Set upper bounds
    a$logLik <- TRUE  # Enable log-likelihood flag
    if (!standardize) {
      a$invchol <- Lambda  # Use non-standardized Lambda
      ret <- do.call("sldpmvnorm", a)  # Calculate score
      return(list(Lambda = ret$invchol, obs = ret$obs, mean = ret$mean, lower = ret$lower, upper = ret$upper))  # Return score and related information
    }
    chol <- solve(Lambda)  # Calculate the inverse matrix of Lambda
    D <- sqrt(Tcrossprod(chol, diag_only = TRUE))  # Calculate the standardization factor
    a$invchol <- sLambda <- invcholD(Lambda, D = D)  # Standardize Lambda
    ret <- do.call("sldpmvnorm", a)  # Calculate standardized score
    ret$chol <- -vectrick(sLambda, ret$invchol)  # Calculate the correction term
    smean <- ret$mean
    sobs <- ret$obs
    slower <- ret$lower
    supper <- ret$upper
    ret <- destandardize(chol = chol, invchol = Lambda, score_schol = ret$chol)  # Destandardize
    ret <- list(Lambda = ret, mean = smean, obs = sobs, lower = slower, upper = supper)  # Return destandardized results
    return(ret)
  }
  
  return(list(logLik = ll, score = sc))  # Return log-likelihood and score functions
}

.start <- function(m, xnames, fixed = NULL) {
  J <- length(m$models)  # Get the number of models
  Jp <- J * (J - 1) / 2  # Number of off-diagonal elements
  Jnames <- m$names  # Model names
  
  # Extract marginal model parameters
  margin_par <- do.call("c", lapply(m$models, function(mod) coef(as.mlt(mod))))
  names(margin_par) <- paste(rep(Jnames, time = m$nparm), names(margin_par), sep = ".")  # Name marginal parameters
  
  # Construct Lambda parameter names
  rn <- rownames(unclass(ltMatrices(1:Jp, names = Jnames, byrow = TRUE)))
  lnames <- paste(rep(rn, each = length(xnames)), rep(xnames, length(rn)), sep = ".")
  
  # Initialize Lambda parameters to 0
  lambda_par <- rep(0, length(lnames))
  names(lambda_par) <- lnames
  
  # Combine marginal parameters and Lambda parameters
  start <- c(margin_par, lambda_par)
  
  # If there are fixed parameters, verify that all fixed parameters are in start
  if (!is.null(fixed)) 
    stopifnot(all(fixed %in% names(start)))
  
  return(start)  # Return starting parameters
}

.model_matrix <- function(models, j = 1, newdata = NULL, prime = FALSE) {
  # If new data is not provided
  if (is.null(newdata)) {
    if (models$cont[j]) {  # If it's a continuous variable
      if (prime) 
        return(models$mm[[j]]$eY$Yprime)  # Return the first derivative matrix
      return(models$mm[[j]]$eY$Y)  # Return the model matrix
    }
    stopifnot(!prime)  # Ensure discrete variables do not request derivatives
    Yleft <- models$mm[[j]]$iY$Yleft  # Get the left boundary of discrete variables
    Yright <- models$mm[[j]]$iY$Yright  # Get the right boundary of discrete variables
    return(list(Yleft = Yleft, Yright = Yright))  # Return left and right boundaries
  }
  
  # Process new data
  resp <- models$models[[j]]$model$response  # Get the name of the response variable
  y <- R(newdata[[resp]])  # Extract the response variable from new data
  
  # Determine whether it is a continuous variable or a double type
  if (models$cont[j] || mlt:::.type_of_response(y) == "double") {
    if (prime) {  # If derivative is requested
      drv <- 1L
      names(drv) <- models$models[[j]]$model$response
      return(model.matrix(models$models[[j]]$model, data = newdata, deriv = drv))  # Return the first derivative matrix
    } else {
      return(model.matrix(models$models[[j]]$model, data = newdata))  # Return the model matrix
    }
  }
  
  # If it's a discrete variable, compute the model matrix for interval boundaries
  iY <- mlt:::.mm_interval(model = models$models[[j]]$model, data = newdata, resp, y)
  Yleft <- iY$Yleft  # Left boundary
  Yright <- iY$Yright  # Right boundary
  return(list(Yleft = Yleft, Yright = Yright))  # Return left and right boundaries
}

.rbind <- function(x) {
  # If x is not a list, directly convert it to a matrix and return
  if (!is.list(x)) 
    return(matrix(x, nrow = 1))
  
  # Bind the elements of the list by rows and return the matrix
  return(do.call("rbind", x))
}

.mget <- function(models, j = 1, parm, newdata = NULL, what = c("trafo", "dtrafo", "z", "zleft", "dzleft", "zright", "dzright", "zprime", "mm", "mmprime", "estfun", "scale"), ...) {
  what <- match.arg(what)  # Ensure that 'what' is one of the allowed values
  
  # If 'j' is a vector, recursively call .mget
  if (length(j) > 1) {
    ret <- lapply(j, .mget, models = models, parm = parm, newdata = newdata, what = what)
    return(ret)
  }
  
  # If 'scale' is requested
  if (what == "scale") {
    if (models$cont[j]) {  # Continuous variable
      Y <- models$mm[[j]]$eY$Y
    } else {  # Discrete variable
      Y <- models$mm[[j]]$iY$Yleft
    }
    Ytmp <- Y
    Ytmp[!is.finite(Ytmp)] <- NA  # Replace infinite values with NA
    sc <- apply(abs(Ytmp), 2, max, na.rm = TRUE)  # Compute the maximum absolute value of each column
    lt1 <- sc < 1.1
    gt1 <- sc >= 1.1
    sc[gt1] <- 1/sc[gt1]  # For values greater than or equal to 1.1, take the reciprocal
    sc[lt1] <- 1  # For values less than 1.1, set to 1
    return(sc)  # Return the scaling parameters
  }
  
  # Extract model parameters and update the model
  prm <- models$parm(parm)[[j]]
  tmp <- models$models[[j]]
  cf <- coef(tmp)
  cf[] <- prm
  coef(tmp) <- cf
  
  # If new data is not provided
  if (is.null(newdata)) {
    if (models$nn[j]) 
      newdata <- tmp$data
  }
  
  # Handle continuous variables
  if (models$cont[j]) {
    if (is.null(newdata)) {
      tr <- c(models$mm[[j]]$eY$Y %*% prm)  # Transformed values
      trp <- c(models$mm[[j]]$eY$Yprime %*% prm)  # First derivatives
      if (!models$normal[j]) 
        trd <- tmp$todistr$d(tr) * trp  # Compute density
    } else {
      tr <- predict(tmp, newdata = newdata, type = "trafo", ...)  # Transformed values
      drv <- 1L
      names(drv) <- tmp$model$response
      trp <- predict(tmp, newdata = newdata, type = "trafo", deriv = drv, ...)  # First derivatives
      if (!models$normal[j]) 
        trd <- predict(tmp, newdata = newdata, type = "density", ...)  # Compute density
    }
  } else {
    # Handle discrete variables
    if (is.null(newdata)) {
      trl <- c(models$mm[[j]]$iY$Yleft %*% prm)  # Left boundary transformed values
      trl[!is.finite(trl)] <- -Inf  # Replace infinite values with negative infinity
      trr <- c(models$mm[[j]]$iY$Yright %*% prm)  # Right boundary transformed values
      trr[!is.finite(trr)] <- Inf  # Replace infinite values with positive infinity
    } else {
      mmj <- .model_matrix(models, j = j, newdata = newdata)
      if (is.matrix(mmj)) 
        return(c(mmj %*% prm))  # Return transformed values
      trl <- c(mmj$Yleft %*% prm)  # Left boundary transformed values
      trl[!is.finite(trl)] <- -Inf
      trr <- c(mmj$Yright %*% prm)  # Right boundary transformed values
      trr[!is.finite(trr)] <- Inf
    }
  }
  
  # Return the requested result based on 'what'
  if (what == "trafo") {
    return(tr)
  }
  if (what == "dtrafo") {
    return(tmp$todistr$d(tr))
  }
  if (what == "z") {
    if (models$normal[j]) 
      return(tr)
    return(qnorm(tmp$todistr$p(tr, log = TRUE), log.p = TRUE))
  }
  if (what == "zleft") {
    if (models$normal[[j]]) 
      return(trl)
    return(qnorm(tmp$todistr$p(trl, log = TRUE), log.p = TRUE))
  }
  if (what == "dzleft") {
    if (models$normal[[j]]) 
      return(rep(1, length(trl)))
    qn <- qnorm(tmp$todistr$p(trl, log = TRUE), log.p = TRUE)
    dn <- dnorm(qn)
    dn[!is.finite(dn)] <- 1
    return(tmp$todistr$d(trl) / dn)
  }
  if (what == "zright") {
    if (models$normal[[j]]) 
      return(trr)
    return(qnorm(tmp$todistr$p(trr, log = TRUE), log.p = TRUE))
  }
  if (what == "dzright") {
    if (models$normal[[j]]) 
      return(rep(1, length(trr)))
    qn <- qnorm(tmp$todistr$p(trr, log = TRUE), log.p = TRUE)
    dn <- dnorm(qn)
    dn[!is.finite(dn)] <- 1
    return(tmp$todistr$d(trr) / dn)
  }
  if (what == "zprime") {
    if (models$normal[[j]]) 
      return(trp)
    qn <- qnorm(tmp$todistr$p(tr, log = TRUE), log.p = TRUE)
    return(trd / dnorm(qn))
  }
  if (what == "estfun") {
    if (is.null(newdata)) 
      return(estfun(tmp))
    return(estfun(tmp, newdata = newdata))
  }
}

.log <- function(x) {
  # Compute element-wise logarithm, ensuring a minimum value to avoid log of zero
  return(log(.pmax(.Machine$double.eps, x)))
}

.pmax <- function(x, y) {
  # Compute element-wise maximum, ensuring each element is at least 'x'
  y[y < x] <- x
  return(y)
}

ltMatrices <- function(object, diag = FALSE, byrow = FALSE, names = TRUE) {
  
  # If 'object' is not a matrix, convert it to a matrix
  if (!is.matrix(object)) 
    object <- matrix(object, ncol = 1L)
  
  # If 'object' is already of type 'ltMatrices', reorder and return
  if (inherits(object, "ltMatrices")) {
    ret <- .reorder(object, byrow = byrow)  # Call internal function .reorder to reorder
    return(ret)
  }
  
  # Calculate the matrix dimension J
  J <- floor((1 + sqrt(1 + 4 * 2 * nrow(object))) / 2 - diag)
  
  # Check if the number of rows in 'object' matches the number of lower triangular elements of J
  if (nrow(object) != J * (J - 1) / 2 + diag * J) 
    stop("Dimension of object does not correspond to lower \n              triangular part of a square matrix")
  
  # Handle the 'names' parameter
  nonames <- FALSE
  if (!isTRUE(names)) {
    if (is.character(names)) 
      stopifnot(is.character(names) && length(unique(names)) == J)  # Ensure 'names' is a unique character vector
    else 
      nonames <- TRUE
  } else {
    names <- as.character(1:J)  # Default to character vector from 1 to J
  }
  
  # If 'names' are available, set row names for the matrix
  if (!nonames) {
    L1 <- matrix(names, nrow = J, ncol = J)  # Create row name matrix
    L2 <- matrix(names, nrow = J, ncol = J, byrow = TRUE)  # Create column name matrix
    L <- matrix(paste(L1, L2, sep = "."), nrow = J, ncol = J)  # Combine row and column names
    
    # Set row names based on 'byrow'
    if (byrow) 
      rownames(object) <- t(L)[upper.tri(L, diag = diag)]
    else 
      rownames(object) <- L[lower.tri(L, diag = diag)]
  }
  
  # Set attributes for 'object'
  attr(object, "J") <- J  # Set matrix dimension
  attr(object, "diag") <- diag  # Whether to include the diagonal
  attr(object, "byrow") <- byrow  # Whether stored by row
  attr(object, "rcnames") <- names  # Row and column names
  
  # Convert 'object' to 'ltMatrices' type
  class(object) <- c("ltMatrices", class(object))
  
  # Return the processed 'object'
  return(object)
}


.models <- function (...) 
{
  # Convert each input to an 'mlt' object
  m <- lapply(list(...), function(x) as.mlt(x))
  
  # Create abbreviated names based on the response variable of each model, limited to 4 characters
  nm <- abbreviate(sapply(m, function(x) x$model$response), 4)
  
  J <- length(m)  # Number of models
  Jp <- J * (J - 1)/2  # Number of off-diagonal elements for a symmetric matrix
  
  # Determine which models have a 'normal' distribution
  normal <- sapply(m, function(x) x$todistr$name == "normal")
  
  # Extract weights from each model
  w <- lapply(m, weights)
  
  # Ensure all models have identical weights
  out <- lapply(w, function(x) stopifnot(isTRUE(all.equal(x, w[[1]]))))
  
  w <- w[[1L]]  # Use the weights from the first model
  
  # If all weights are 1, set weights to FALSE (no weights)
  if (isTRUE(all.equal(unique(w), 1))) 
    w <- FALSE
  
  # Extract transformation matrices (eY for continuous, iY for discrete) from each model
  mm <- lapply(m, function(mod) {
    eY <- get("eY", environment(mod$parm))  # Extract 'eY' from the model's environment
    iY <- get("iY", environment(mod$parm))  # Extract 'iY' from the model's environment
    list(eY = eY, iY = iY)
  })
  
  # Identify which models are continuous and which are discrete
  cmod <- sapply(mm, function(x) !is.null(x$eY))  # Continuous models have 'eY'
  dmod <- sapply(mm, function(x) !is.null(x$iY))  # Discrete models have 'iY'
  
  # Ensure that each model is either continuous or discrete, but not both
  stopifnot(all(xor(cmod, dmod)))
  
  # Ensure that continuous models come before discrete models
  stopifnot(all(diff(cmod) <= 0))
  stopifnot(all(diff(dmod) >= 0))
  
  # Get the number of observations, ensuring all models have the same number
  nobs <- unique(sapply(m, nobs))
  stopifnot(length(nobs) == 1L)
  nobs <- nobs[[1L]]
  
  # Get the number of parameters for each model
  P <- sapply(m, function(x) length(coef(x)))
  
  # Create a factor to split parameters by model
  fpar <- factor(rep(1:J, P))
  
  # Define a function to split parameters based on the factor
  parm <- function(par) {
    mpar <- par[1:sum(P)]  # Extract marginal parameters
    split(mpar, fpar)  # Split parameters by model
  }
  
  # Extract constraints from each model's transformation matrices
  constr <- lapply(mm, function(m) {
    if (is.null(m$eY)) 
      return(attr(m$iY$Yleft, "constraint"))  # Constraints for discrete models
    return(attr(m$eY$Y, "constraint"))  # Constraints for continuous models
  })
  
  # Combine all constraint matrices into a block diagonal matrix
  ui <- do.call("bdiag", lapply(constr, function(x) x$ui))
  
  # Combine all constraint right-hand side values into a single vector
  ci <- do.call("c", lapply(constr, function(x) x$ci))
  
  # Remove constraints with non-finite right-hand side values
  ui <- as(ui[is.finite(ci), , drop = FALSE], "matrix")
  ci <- ci[is.finite(ci)]
  
  # Create model frames for each model
  mf <- lapply(1:J, function(j) {
    mf <- m[[j]]$data  # Extract data from the j-th model
    if (cmod[j]) 
      return(mf)  # For continuous models, return the data as is
    yl <- m[[j]]$response$cleft  # Left boundary for discrete response
    yr <- m[[j]]$response$cright  # Right boundary for discrete response
    rp <- m[[j]]$model$response  # Response variable name
    ml <- mr <- mf  # Create copies for left and right boundaries
    ml[[rp]] <- yl  # Assign left boundary to the left model frame
    mr[[rp]] <- yr  # Assign right boundary to the right model frame
    return(list(left = ml, right = mr))  # Return list of left and right model frames
  })
  
  # Identify models with fixed parameters or non-zero offsets or scaling shifts
  nn <- sapply(1:J, function(j) {
    !is.null(m[[j]]$fixed) || !isTRUE(all.equal(unique(m[[j]]$offset), 0)) || m[[j]]$model$scale_shift
  })
  
  # Determine the type of response for each model
  type <- lapply(1:J, function(j) mlt:::.type_of_response(m[[j]]$response))
  
  # Return a list containing all the extracted and processed information
  return(list(
    models = m,        # List of 'mlt' models
    mf = mf,           # List of model frames
    cont = cmod,       # Logical vector indicating continuous models
    type = type,       # List of response types
    normal = normal,   # Logical vector indicating normal distributions
    nobs = nobs,       # Number of observations
    weights = w,       # Weights (FALSE if all weights are 1)
    nparm = P,         # Number of parameters per model
    parm = parm,       # Function to split parameters by model
    ui = ui,           # Constraint matrix
    ci = ci,           # Constraint right-hand side
    mm = mm,           # List of transformation matrices
    names = nm,        # Abbreviated model names
    nn = nn            # Logical vector for models with fixed parameters or other constraints
  ))
}
