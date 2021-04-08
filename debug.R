library(bsae)
library(sf)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

sf <- readRDS(file = "data/mw.rds")

l_fixed <- 1
cov <- centroid_covariance(sf, l = l_fixed)
cov <- cov / riebler_gv(cov) # Standardise so tau prior is right
C <- Matrix::solve(cov) # Precision matrix

tau_prior <- list(prec = list(prior = "logtnormal", param = c(0, 1/2.5^2), initial = 0, fixed = FALSE))
beta_prior <- list(mean.intercept = -2, prec.intercept = 1)

fit_inla <- INLA::inla(y ~ 1 + f(id, model = "generic0", Cmatrix = C, hyper = tau_prior),
                       family = "xbinomial",
                       control.family = list(control.link = list(model = "logit")),
                       control.fixed = beta_prior,
                       data = list(id = 1:nrow(sf), y = sf$y, m = sf$n_obs),
                       Ntrials = m,
                       control.predictor = list(compute = TRUE, link = 1),
                       control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE))

D <- centroid_distance(sf)

# debug.stan: fixed lengthscale

fit_stan <- rstan::stan("stan/debug.stan",
                        data = list(n = nrow(sf), y = sf$y, m = sf$n_obs, mu = rep(0, nrow(sf)), D = D, l = l_fixed),
                        warmup = 100,
                        iter = 4000)

summary_stan <- rstan::summary(fit_stan)$summary

fit_inla$summary.fixed
summary_stan["beta_0", ]

fit_inla$summary.hyperpar
summary_stan["tau_phi", ]

plot(fit_inla$summary.random$id[, "mean"], summary_stan[2:(1 + 28), "mean"])
abline(a = 0, b = 1)

# debug-01.stan: fixed lengthscale, Gram matrix passed in outright

fit_inla_01 <- fit_inla

fit_stan_01 <- rstan::stan("stan/debug-01.stan",
                           data = list(n = nrow(sf), y = sf$y, m = sf$n_obs, mu = rep(0, nrow(sf)), K = cov),
                           warmup = 100,
                           iter = 4000)

summary_stan_01 <- rstan::summary(fit_stan_01)$summary

fit_inla_01$summary.fixed
summary_stan["beta_0", ]

fit_inla_01$summary.hyperpar
summary_stan["tau_phi", ]

plot(fit_inla_01$summary.random$id[, "mean"], summary_stan[2:(1 + 28), "mean"])
abline(a = 0, b = 1)

# debug-02.stan: integers to remove xbinomial

fit_inla_02 <- INLA::inla(y ~ 1 + f(id, model = "generic0", Cmatrix = C, hyper = tau_prior),
                          family = "binomial",
                          control.family = list(control.link = list(model = "logit")),
                          control.fixed = beta_prior,
                          data = list(id = 1:nrow(sf), y = floor(sf$y), m = floor(sf$n_obs)),
                          Ntrials = m,
                          control.predictor = list(compute = TRUE, link = 1),
                          control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE))

fit_stan_02 <- rstan::stan("stan/debug-02.stan",
                           data = list(n = nrow(sf), y = floor(sf$y), m = floor(sf$n_obs), mu = rep(0, nrow(sf)), K = cov),
                           warmup = 100,
                           iter = 4000)

summary_stan <- rstan::summary(fit_stan_02)$summary

# Compare
fit_inla_02$summary.fixed
summary_stan["beta_0", ]

fit_inla_02$summary.hyperpar
summary_stan["tau_phi", ]

plot(fit_inla_02$summary.$id[, "mean"], summary_stan[2:(1 + 28), "mean"])
abline(a = 0, b = 1)
