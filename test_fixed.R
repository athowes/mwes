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

fit_inla <- INLA::inla(y ~ 1 + f(id, model = "generic0", Cmatrix = C, hyper = tau_prior),
                       family = "xbinomial",
                       control.family = list(control.link = list(model = "logit")),
                       data = list(id = 1:nrow(sf), y = sf$y, m = sf$n_obs),
                       Ntrials = m,
                       control.predictor = list(compute = TRUE, link = 1),
                       control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE))

D <- centroid_distance(sf)

fit_stan <- rstan::stan("stan/fixed-centroid.stan",
                        data = list(n = nrow(sf), y = sf$y, m = sf$n_obs, mu = rep(0, nrow(sf)), D = D, l = l_fixed),
                        warmup = 100,
                        iter = 4000)

draws <- rstan::extract(fit_stan)
summary_stan <- rstan::summary(fit_stan)$summary

# Compare
fit_inla$summary.fixed
summary_stan["beta_0", ]

fit_inla$summary.hyperpar
summary_stan["tau_phi", ]

phi_mean_inla <- fit_inla$summary.random$id[, "mean"]
phi_mean_stan <- summary_stan[2:(1 + 28), "mean"]
plot(phi_mean_inla, phi_mean_stan)
abline(a = 0, b = 1)
