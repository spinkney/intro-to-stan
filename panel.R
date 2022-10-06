library(data.table)
library(cmdstanr)
library(extrafont)
extrafont::font_import()

dat <- data.table(n_big = c(2862, 2288, 3119, NA, 4072, NA, 2638, 3981, 2773, 1827),
                  n_small     = c(55, 63, 22, 99, 16, 13, 87, 84, 49, 62),
                  hours_big = c(28620, 22440, 31190, NA, 24432, NA, 21104, 27867, 22184, 11735),
                  hours_small = c(165, 315, 66, 396, 80, 65, 435, 84, 196, 286))

total_big <- 20000
total_small <- 500

cor(dat$n_small, dat$hours_small)
dat[is.na(n_big) == F, cor(n_big, hours_big)]


data_normal <- list(
  N = nrow(dat),
  N_unknown = dat[is.na(n_big), .N],
  y_small = dat[, n_small / 10],
  
  y = dat[!is.na(n_big), n_big / 400],
  hours_big = dat[!is.na(n_big), hours_big], 
  n_big = dat[!is.na(n_big), n_big], 
  
  n_small = dat[, n_small],
  hours_small = dat[, hours_small],
  
  unknown = which(is.na(dat[, n_big])),
  known = which(!is.na(dat[, n_big]))
)

mod <- cmdstan_model("panel_fusion_student.stan")
fit <- mod$sample(
  data = data_normal,
  parallel_chains = 4,
  seed = 12312,
  iter_sampling = 500,
  iter_warmup = 500
)

fit$summary("rho_xy")

out <- as.data.table(fit$summary(c("y_out", "x_out")))
out[, variable := c(rep("hhs_impute", 2), rep("hours_impute", 2))]
out

fit$summary("beta_x")

library(bayesplot)
yrep <- fit$draws("y_rep", format = "matrix")
y <- c(dat$n_big,
           dat$n_small)

indx <- which(!is.na(y))[1:8]

ppc_dens_overlay(y = y[indx], yrep = yrep[1:400, indx])
ppc_dens_overlay(y = y[11:20], yrep = yrep[1:400, 11:20])

xrep <- fit$draws("x_rep", format = "matrix")
x <- c(dat$hours_big,
       dat$hours_small)

indx <- which(!is.na(x))[1:8]

ppc_dens_overlay(y = x[indx], yrep = xrep[1:400, indx])
ppc_dens_overlay(y = x[11:20], yrep = xrep[1:400, 11:20])

fit$summary("nu")
