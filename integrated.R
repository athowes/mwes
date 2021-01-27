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
nsim_warm <- 10

# The number of iterations in the MCMC after warm-up
nsim_iter <- 90

fit <- rstan::stan("integrated.stan",
                   data = dat,
                   warmup = nsim_warm,
                   iter = nsim_iter)

samples <- rstan::extract(fit)

plot(samples$beta_0) # Intercept
plot(samples$l) # Length-scale
