#!/usr/bin/env Rscript

# R WebAssembly Server for RJMCMC Widget
# This script creates a WebAssembly-compatible R server that can run the rjmc package

library(webr)
library(jsonlite)
library(httr)

# Load the rjmc package
if (!require("rjmc", character.only = TRUE)) {
  # If not installed, try to load from local development
  devtools::load_all("/Users/davidhodgson/Dropbox/Mac (3)/Documents/research/software/mcmc/rjmc")
}

# Define the mixture model from the vignette
mixture_model <- list(
  lowerParSupport_fitted = c(-8),
  upperParSupport_fitted = c(8),
  namesOfParameters = c("sigma"),
  
  sampleInitPrior = function(datalist) {
    rnorm(1, 0, 1)
  },
  
  sampleInitJump = function(params, datalist) {
    p <- c(0.5, 0.5)
    mu <- c(1, 4)
    sigma <- c(1, 1)
    jump_new <- matrix(c(p, mu, sigma), nrow = 3, byrow = TRUE)
    jump_new
  },
  
  evaluateLogPrior = function(params, jump, datalist) {
    p <- dnorm(params[1], 0, 1, log = TRUE)
    
    p_vec <- jump[1, ]
    for(i in 1:length(p_vec) ) {
      p <- p + dunif(p_vec[i], 0, 1, log = TRUE)
    }
    mu_vec <- jump[2, ]
    for(i in 1:length(mu_vec) ) {
      p <- p + dunif(mu_vec[i], 0, 20, log = TRUE)
    }
    sigma_vec <- jump[3, ]
    for(i in 1:length(sigma_vec) ) {
      p <- p + dunif(sigma_vec[i], 0.3, 3, log = TRUE)
    }
    N <- ncol(jump)
    p <- p + dunif(N, 1, 8)
    
    p
  },
  
  evaluateLogLikelihood = function(params, jump, datalist) {
    ll <- 0
    N <- ncol(jump)
    p_vec <- jump[1, ]
    mu_vec <- jump[2, ]
    sigma_vec <- jump[3, ]
    
    z <- datalist$obs
    N_data <- datalist$N_data
    sigma <- params[1]
    
    for (i in 1:N_data) {
      i_x = 0
      for (j in 1:N) {
        i_x <- i_x + p_vec[j] * dnorm(z[i], mu_vec[j], sigma_vec[j])
      }
      ll <- ll + log(i_x)
    }
    
    ll
  },
  
  sampleBirthProposal = function(params, jump, i_idx, datalist) {
    p_new_sample <- runif(1, 0, 1)
    p_new <- c(jump[1, ] * (1 - p_new_sample), p_new_sample)
    mu_new_sample <- runif(1, 0, 20)
    mu_new <- c(jump[2, ], mu_new_sample)
    sigma_new_sample <- runif(1, 0.3, 3)
    sigma_new <- c(jump[3, ], sigma_new_sample)
    jump_new <- matrix(c(p_new, mu_new, sigma_new), nrow = 3, byrow = TRUE)
    jump_new
  },
  
  sampleDeathProposal = function(params, jump, i_idx, datalist) {
    N <- ncol(jump)
    jump_remove <- jump[, i_idx]
    jump_new <- jump[, -i_idx]
    jump_new[1, ] <- c(jump_new[1, ] / (1 - jump_remove[1]))
    jump_new
  },
  
  evaluateBirthProposal = function(params, jump, i_idx, datalist) {
    N <- ncol(jump)
    log(1 / (N * dunif(jump[2, N], 0, 20 ) ))
  },
  
  evaluateDeathProposal = function(params, jump, i_idx, datalist) {
    N <- ncol(jump)
    log((N) * dunif(jump[2, i_idx], 0, 20 ) ) 
  },
  
  sampleJump = function(params, jump, i_idx, datalist) {
    N <- ncol(jump)
    jump_update <- jump[, i_idx]
    
    p_new <- min(max(jump_update[1] + rnorm(1, 0, 0.01), 0), 1)
    diff = (jump_update[1] - p_new) / (N - 1)
    p_new_vew <- jump[1, ] + diff
    p_new_vew[i_idx] <- p_new
    jump[1, ] <- p_new_vew
    
    jump[2, i_idx] <- jump_update[2] + rnorm(1, 0, 0.1)
    jump[3, i_idx] <- max(jump_update[3] + rnorm(1, 0, 0.1), 0.3)
    
    jump
  },
  
  sampleProposal = function(params, jump, datalist) {
    N <- ncol(jump)
    if (N == 2) {
      q <- c(0.0, 0.67, 1.0)
    } else if (N == 20) {
      q <- c(0.33, 1.0, 1.0)
    } else {
      q <- c(0.33, 0.67, 1.0)
    }
    q
  }
)

# Function to run RJMCMC with streaming results
run_rjmc_streaming <- function(data, settings, callback_url = NULL) {
  tryCatch({
    # Prepare data
    data_l <- list(
      obs = data$obs,
      N_data = length(data$obs)
    )
    
    # Run RJMCMC
    outputs <- rjmc_func(mixture_model, data_l, settings)
    
    # Process results for streaming
    results <- list(
      success = TRUE,
      outputs = outputs,
      message = "RJMCMC completed successfully"
    )
    
    # If callback URL provided, send results
    if (!is.null(callback_url)) {
      httr::POST(callback_url, 
                body = toJSON(results, auto_unbox = TRUE),
                content_type("application/json"))
    }
    
    return(results)
    
  }, error = function(e) {
    error_result <- list(
      success = FALSE,
      error = as.character(e$message),
      message = "RJMCMC failed"
    )
    
    if (!is.null(callback_url)) {
      httr::POST(callback_url, 
                body = toJSON(error_result, auto_unbox = TRUE),
                content_type("application/json"))
    }
    
    return(error_result)
  })
}

# Function to generate sample data like in vignette
generate_sample_data <- function(n_samples = 1000) {
  m <- 5
  obs <- c(sample(rnorm(10e4, 2, 1), 30 * m, replace = TRUE),
           sample(rnorm(10e4, 6, 0.5), 70 * m, replace = TRUE), 
           sample(rnorm(10e4, 12, 1), 100 * m, replace = TRUE))
  
  return(list(obs = obs))
}

# Function to get posterior summary
get_posterior_summary <- function(outputs) {
  tryCatch({
    n_chain <- 4
    tables_length <- get_lengths(outputs, n_chain)
    
    # Get most probable K
    mode_post <- names(sort(tables_length, decreasing = TRUE))[1]
    
    # Get posterior samples for the most probable K
    post_process <- get_clean_posterior(outputs, mode_post, n_chain)
    
    summary <- list(
      most_probable_K = as.numeric(mode_post),
      K_probabilities = as.list(tables_length),
      posterior_samples = post_process,
      success = TRUE
    )
    
    return(summary)
    
  }, error = function(e) {
    return(list(
      success = FALSE,
      error = as.character(e$message)
    ))
  })
}

# Export functions for WebAssembly
webr::serve({
  list(
    run_rjmc_streaming = run_rjmc_streaming,
    generate_sample_data = generate_sample_data,
    get_posterior_summary = get_posterior_summary,
    mixture_model = mixture_model
  )
}, port = 8080)
