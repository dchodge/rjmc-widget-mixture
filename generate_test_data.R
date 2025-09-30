#!/usr/bin/env Rscript

# Generate 5 different mixture distributions for testing RJMCMC
# Each mixture has random number of components (2-8) and random scale (10-100)

# Set random seed for reproducibility
set.seed(42)

# Create output directory
dir.create("sampled data for testing", showWarnings = FALSE)

# Function to generate mixture data
generate_mixture_data <- function(n_components, n_samples = 100, scale_range = c(10, 100)) {
  # Random scale factor
  scale <- runif(1, scale_range[1], scale_range[2])
  
  # Generate random mixture parameters
  # Means: spread across the scale range
  means <- sort(runif(n_components, -scale/2, scale/2))
  
  # Standard deviations: random but reasonable
  sigmas <- runif(n_components, 0.5, scale/10)
  
  # Mixing proportions: random but normalized
  p <- runif(n_components, 0.1, 0.9)
  p <- p / sum(p)  # Normalize to sum to 1
  
  # Generate samples
  samples <- numeric(n_samples)
  for (i in 1:n_samples) {
    # Choose component
    component <- sample(n_components, 1, prob = p)
    # Sample from chosen component
    samples[i] <- rnorm(1, means[component], sigmas[component])
  }
  
  return(list(
    data = samples,
    true_parameters = list(
      n_components = n_components,
      means = as.list(means),
      sigmas = as.list(sigmas),
      mixing_proportions = as.list(p),
      scale = scale
    ),
    metadata = list(
      n_samples = n_samples,
      description = paste0("Mixture with ", n_components, " components, scale=", round(scale, 1))
    )
  ))
}

# Generate 5 test datasets
cat("Generating 5 test mixture distributions...\n\n")

for (i in 1:5) {
  # Random number of components (2-8)
  n_components <- sample(2:8, 1)
  
  # Generate data
  data <- generate_mixture_data(n_components, n_samples = 100)
  
  # Save to CSV file (just the data values)
  filename <- paste0("mixture_test_", i, "_K", n_components, ".csv")
  filepath <- file.path("sampled data for testing", filename)
  
  # Write just the data values to CSV
  write.csv(data$data, filepath, row.names = FALSE, col.names = FALSE)
  
  cat("Generated", filename, ":\n")
  cat("  - Components:", n_components, "\n")
  cat("  - Scale:", round(data$true_parameters$scale, 1), "\n")
  cat("  - Means:", paste(sprintf("%.2f", data$true_parameters$means), collapse = ", "), "\n")
  cat("  - Sigmas:", paste(sprintf("%.2f", data$true_parameters$sigmas), collapse = ", "), "\n")
  cat("  - Mixing props:", paste(sprintf("%.3f", data$true_parameters$mixing_proportions), collapse = ", "), "\n")
  cat("\n")
}

cat("All test data saved to 'sampled data for testing/' directory\n")

# Generate a summary file
summary <- list(
  description = "Test mixture distributions for RJMCMC",
  generated_files = paste0("mixture_test_", 1:5, "_K", sample(2:8, 5), ".csv"),
  parameters = list(
    n_samples_per_dataset = 100,
    component_range = c(2, 8),
    scale_range = c(10, 100),
    random_seed = 42
  )
)

summary_json <- jsonlite::toJSON(summary, pretty = TRUE, auto_unbox = TRUE)
writeLines(summary_json, file.path("sampled data for testing", "README.json"))

cat("Summary saved to README.json\n")
