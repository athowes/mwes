library(sf)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

sf <- readRDS(file = "mw.rds")

centroid_distance <- function(sf) {
  cent <- sf::st_centroid(sf)
  D <- sf::st_distance(cent, cent)
  return(D)
}

D <- centroid_distance(sf)

ii_obs <- which(!is.na(sf$y))
ii_mis <- which(is.na(sf$y))
n_obs <- length(ii_obs)
n_mis <- length(ii_mis)

dat <- list(n_obs = n_obs,
            n_mis = n_mis,
            ii_obs = array(ii_obs),
            ii_mis = array(ii_mis),
            n = nrow(sf),
            y_obs = sf$y[ii_obs],
            m = sf$n_obs,
            mu = rep(0, nrow(sf)),
            D = D)

nsim_warm <- 100
nsim_iter <- 400

fit <- rstan::stan("centroid.stan",
                   data = dat,
                   warmup = nsim_warm,
                   iter = nsim_iter)

saveRDS(fit, file = "fit_centroid.rds")
