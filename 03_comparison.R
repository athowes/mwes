fit_centroid <- readRDS("fit_centroid.rds")
fit_integrated <- readRDS("fit_integrated.rds")

# Sampling speed

samples_per_second <- function(fit) {
  # Outside of warm-up phase
  times <- get_elapsed_time(fit)
  # Number of samples in first chain
  nsim_iter <- fit@stan_args[[1]]$iter
  nsim_iter / mean(times[, 2])
}

samples_per_second(fit_centroid)
samples_per_second(fit_integrated)

# Convergence analysis

# All Rhat < 1.1 necessary but not sufficient
rstan::summary(fit_centroid)[["summary"]][, "Rhat"]
rstan::summary(fit_integrated)[["summary"]][, "Rhat"]

# Fitted values

fitted_centroid <- summary(fit_centroid)$summary[paste0("rho[", 1:28, "]"), "mean"]
fitted_integrated <- summary(fit_integrated)$summary[paste0("rho[", 1:28, "]"), "mean"]

plot(fitted_centroid, fitted_integrated)
