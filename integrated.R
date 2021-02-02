library(sf)
library(rstan)

sf <- readRDS(file = "mw.rds")

# The number of integration points within each area
L <- 10

# The method for selecting integration points
type <- "random"

# The number of areas
n <- nrow(sf)

# Draw L samples from each area according to method type
samples <- sf::st_sample(sf, type = type, size = rep(L, n))

library(ggplot2)
plot_samples <- function(samples){
  ggplot(sf) +
    geom_sf(fill = "lightgrey") +
    geom_sf(data = samples, alpha = 0.5, shape = 4) +
    labs(x = "Longitude", y = "Latitude") +
    theme_minimal() +
    labs(fill = "") +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
}

plot_samples(samples)

# Construct an (L * n) * (L * n) matrix containing the Euclidean distance between each sample
# Note that this ^ will not be exact for some settings of type (hexagonal, regular)
# I should add code to force it to be though since Stan requires matrix dimensions
S <- sf::st_distance(samples, samples)
dim(S)

# This is required when some observations are missing (cross-validation)
ii_obs <- which(!is.na(sf$y))
ii_mis <- which(is.na(sf$y))
n_obs <- length(ii_obs)
n_mis <- length(ii_mis)

# Creating the data structure
dat <- list(n_obs = n_obs,
            n_mis = n_mis,
            ii_obs = array(ii_obs),
            ii_mis = array(ii_mis),
            n = nrow(sf),
            y_obs = sf$y[ii_obs],
            m = sf$n_obs,
            mu = rep(0, nrow(sf)),
            L = L,
            S = S)
  
# The number of warm-up iterations in the MCMC
nsim_warm <- 100

# The number of iterations in the MCMC after warm-up
nsim_iter <- 900

fit <- rstan::stan("integrated.stan",
                   data = dat,
                   warmup = nsim_warm,
                   iter = nsim_iter)

# Samples per second outside of warm-up phase
times <- get_elapsed_time(fit)
nsim_iter / mean(times[, 2])

chain <- rstan::extract(fit)

plot(chain$beta_0) # Intercept
plot(chain$l) # Length-scale

# Version which accepts unequal numbers of samples in each area

# Data structure for unequal number of points in each area
sample_index <- sf::st_intersects(sf, samples)
sample_lengths <- lengths(sample_index)
start_index <- sapply(sample_index, function(x) x[1])

dat_unequal <- list(n_obs = n_obs,
                    n_mis = n_mis,
                    ii_obs = array(ii_obs),
                    ii_mis = array(ii_mis),
                    n = nrow(sf),
                    y_obs = sf$y[ii_obs],
                    m = sf$n_obs,
                    mu = rep(0, nrow(sf)),
                    sample_lengths = sample_lengths,
                    total_samples = sum(sample_lengths),
                    start_index = start_index,
                    S = S)

fit_unequal <- rstan::stan("unequal-integrated.stan",
                           data = dat_unequal,
                           warmup = nsim_warm,
                           iter = nsim_iter)

times_unequal <- get_elapsed_time(fit_unequal)
nsim_iter / mean(times_unequal[, 2])

chain_unequal <- rstan::extract(fit_unequal)

plot(chain_unequal$beta_0) # Intercept
plot(chain_unequal$l) # Length-scale
plot(chain_unequal$sigma_phi) # Standard deviation
plot(chain_unequal$tau_phi) # Precision
