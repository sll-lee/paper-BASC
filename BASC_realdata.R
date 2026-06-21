###############################################################################
# Bayesian Adaptive Synthetic Control (BASC): West Germany Real Data Example
#
# This script runs the BASC model for the West Germany application.
###############################################################################

rm(list = ls())

# -----------------------------------------------------------------------------
# 0. Packages
# -----------------------------------------------------------------------------

required_packages <- c(
  "tidysynth",
  "foreign",
  "dplyr",
  "tidyr",
  "tmvtnorm",
  "MASS",
  "MCMCpack",
  "fields",
  "coda",
  "ggplot2",
  "reshape2"
)

invisible(lapply(required_packages, require, character.only = TRUE))

# -----------------------------------------------------------------------------
# 1. User settings
# -----------------------------------------------------------------------------

# Set paths before running.
# The input data file is not included here because of data-licensing issues.
input_data_path <- "data/repgermany.dta"
output_dir <- "results/realdata"
plot_dir <- file.path(output_dir, "plot")

# MCMC settings
N <- 500000
nburn <- 500000
seeds <- c(103, 200, 305)

# Set TRUE to run the sampler. Set FALSE to load a saved result object.
run_mcmc <- TRUE
result_file <- file.path(output_dir, "BASC_realdata_results.RData")

# -----------------------------------------------------------------------------
# 2. Data preparation
# -----------------------------------------------------------------------------

prepare_realdata <- function(input_data_path) {
  data <- read.dta(input_data_path)

  vt <- data %>%
    distinct(year) %>%
    pull(year) %>%
    as.vector() - 1959

  Dt <- ifelse(vt - 31 < 0, 0, 1)

  data.y <- data %>%
    filter(country == "West Germany") %>%
    arrange(year)

  raw.y <- data.y$gdp
  year.y <- data.y$year

  pre_idx_y <- year.y < 1990
  mean_y_pre <- mean(raw.y[pre_idx_y])
  sd_y_pre <- sd(raw.y[pre_idx_y])

  base.x <- data %>%
    filter(country != "West Germany") %>%
    dplyr::select(country, year, gdp) %>%
    pivot_wider(names_from = year, values_from = gdp)

  raw.x <- base.x %>%
    dplyr::select(-country) %>%
    as.matrix()

  # Final manuscript version: original per-capita GDP scale.
  y <- raw.y
  x <- t(raw.x)

  donor_names <- base.x$country

  list(
    data = data,
    vt = vt,
    Dt = Dt,
    y = y,
    x = x,
    donor_names = donor_names,
    mean_y_pre = mean_y_pre,
    sd_y_pre = sd_y_pre
  )
}

realdata <- prepare_realdata(input_data_path)

vt <- realdata$vt
Dt <- realdata$Dt
y <- realdata$y
x <- realdata$x
donor_names <- realdata$donor_names

j <- ncol(x)
t <- length(vt)
q <- 1

# -----------------------------------------------------------------------------
# 3. Helper functions
# -----------------------------------------------------------------------------

safe_chol <- function(Sigma, base_jitter = 1e-8, max_try = 10) {
  Sigma <- (Sigma + t(Sigma)) / 2
  n <- nrow(Sigma)

  scale_S <- max(diag(Sigma), na.rm = TRUE)
  if (!is.finite(scale_S) || scale_S <= 0) scale_S <- 1

  for (k in 0:max_try) {
    jitter <- base_jitter * scale_S * 10^k
    R <- try(chol(Sigma + jitter * diag(n)), silent = TRUE)

    if (!inherits(R, "try-error")) {
      return(list(R = R, jitter = jitter))
    }
  }

  eig_min <- min(eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values)
  jitter <- abs(eig_min) + base_jitter * scale_S * 10^max_try
  R <- chol(Sigma + jitter * diag(n))

  list(R = R, jitter = jitter)
}

loglik_gamma <- function(beta, y, x, Dt, alpha.v, f, sigsq) {
  res <- as.numeric(y - x %*% beta - Dt * alpha.v - f)
  -0.5 * length(y) * log(2 * pi * sigsq) -
    0.5 * sum(res^2) / sigsq
}

log_post_l <- function(gp.l, dataList) {
  tau <- dataList$tau
  vtt <- dataList$vt
  ff <- as.numeric(dataList$f)
  alpha_l <- dataList$alpha_l
  beta_l <- dataList$beta_l

  if (!is.finite(gp.l) || gp.l <= 0) return(-Inf)
  if (!is.finite(tau) || tau <= 0) return(-Inf)
  if (any(!is.finite(ff))) return(-Inf)

  gp.Cov <- tau * exp(-0.5 * rdist(vtt)^2 / gp.l)
  gp.Cov <- (gp.Cov + t(gp.Cov)) / 2

  ch <- try(safe_chol(gp.Cov), silent = TRUE)
  if (inherits(ch, "try-error")) return(-Inf)

  R <- ch$R

  # t(R) %*% R = gp.Cov
  z <- backsolve(R, ff, transpose = TRUE)
  quad <- sum(z^2)
  logdet <- 2 * sum(log(diag(R)))

  log_prior <- -(alpha_l + 1) * log(gp.l) - beta_l / gp.l
  log_gp <- -0.5 * logdet - 0.5 * quad

  log_prior + log_gp
}

rw_mp.l <- function(nsim.l, nburn.l, delta, dataList, initsList) {
  log.curr <- log(initsList$gp.l)
  para.samples <- rep(NA_real_, nsim.l)
  acc.count <- 0

  for (iter.l in 1:(nsim.l + nburn.l)) {
    curr <- exp(log.curr)
    log.prop <- log.curr + rnorm(1, 0, delta)
    prop <- exp(log.prop)

    log_acc <- log_post_l(prop, dataList) -
      log_post_l(curr, dataList) +
      log.prop - log.curr

    accepted <- FALSE
    if (is.finite(log_acc) && log(runif(1)) < log_acc) {
      log.curr <- log.prop
      accepted <- TRUE

      if (iter.l > nburn.l) {
        acc.count <- acc.count + 1
      }
    }

    if (iter.l > nburn.l) {
      para.samples[iter.l - nburn.l] <- exp(log.curr)
    }

    if (iter.l <= 20) {
      cat(
        "iter.l =", iter.l,
        "curr =", curr,
        "prop =", prop,
        "log_acc =", log_acc,
        "accepted =", accepted,
        "\n"
      )
    }
  }

  attr(para.samples, "accept_rate") <- acc.count / nsim.l
  para.samples
}

log_post_u <- function(u, p, dataList) {
  gamma.v <- dataList$gamma.v
  alpha_u <- dataList$alpha_u
  Dt <- dataList$Dt
  alpha.v <- dataList$alpha.v
  f <- dataList$f
  x <- dataList$x
  y <- dataList$y
  sigsq <- dataList$sigsq

  u.j <- u[p]
  denom <- sum(gamma.v * u)

  if (!is.finite(u.j) || u.j <= 0) return(-Inf)
  if (!is.finite(denom) || denom <= 0) return(-Inf)
  if (!is.finite(sigsq) || sigsq <= 0) return(-Inf)

  beta.tem <- gamma.v * u / denom

  res <- as.numeric(y - x %*% beta.tem - Dt * alpha.v - f)

  log_prior <- (alpha_u - 1) * log(u.j) - alpha_u * u.j
  log_lik <- -0.5 * sum(res^2) / sigsq

  log_prior + log_lik
}

compute_beta_samples <- function(results, N, j) {
  for (chain_id in seq_along(results)) {
    results[[chain_id]]$beta.sample <- matrix(0, N, j)

    for (beta.row in 1:N) {
      for (k in 1:j) {
        results[[chain_id]]$beta.sample[beta.row, k] <-
          results[[chain_id]]$gamma.sample[beta.row, k] *
          results[[chain_id]]$u.sample[beta.row, k] /
          sum(results[[chain_id]]$gamma.sample[beta.row, ] *
                results[[chain_id]]$u.sample[beta.row, ])
      }
    }
  }

  results
}

ordinal_suffix <- function(num) {
  suffix <- ifelse(
    num %% 10 == 1 & num %% 100 != 11,
    "st",
    ifelse(
      num %% 10 == 2 & num %% 100 != 12,
      "nd",
      ifelse(num %% 10 == 3 & num %% 100 != 13, "rd", "th")
    )
  )

  paste0(num, suffix)
}

# -----------------------------------------------------------------------------
# 4. MCMC sampler
# -----------------------------------------------------------------------------

run_basc_chain <- function(seed, y, x, vt, Dt, N, nburn, q = 1) {
  j <- ncol(x)
  t <- length(vt)

  # Storage
  f.sample <- matrix(0, N, t)
  sig.sample <- rep(0, N)
  tau.sample <- rep(0, N)
  l.sample <- rep(0, N)
  w.sample <- rep(0, N)
  gamma.sample <- matrix(0, N, j)
  u.sample <- matrix(0, N, j)
  a.sample <- matrix(0, N, q)

  # Prior parameters
  alpha_sig <- 10
  beta_sig <- 5000

  alpha_tau <- 3
  beta_tau <- 2 * 10000

  alpha_l <- 3
  beta_l <- 1000

  alpha_u <- 2.5

  # Initial values
  set.seed(seed)

  tau <- rinvgamma(n = 1, alpha_tau, beta_tau)
  print(tau)

  lscale <- rinvgamma(n = 1, alpha_l, beta_l)
  print(lscale)

  Cov.l <- matrix(0, nrow = t, ncol = t)
  for (i in 1:t) {
    for (k in 1:t) {
      Cov.l[i, k] <- exp(-0.5 * (vt[i] - vt[k])^2 / lscale)
    }
  }

  mu <- rep(0, t)
  f <- mvrnorm(n = 1, mu, tau * Cov.l)
  print(f)

  w <- runif(1, 0, 1)

  gamma.v <- rep(1, j)
  if (sum(gamma.v) == 0) {
    gamma.v[sample(1:j, 1)] <- 1
  }

  u <- rgamma(j, alpha_u, scale = 1 / alpha_u)
  print(u)

  beta <- rep(0, j)
  for (k in 1:j) {
    beta[k] <- gamma.v[k] * u[k] / sum(gamma.v * u)
  }

  sigsq <- rinvgamma(1, alpha_sig, beta_sig)
  print(sigsq)

  sig.a <- rinvgamma(1, 3, 2 * 1000^2)
  mu_a <- rep(-0, q)
  alpha.v <- rtmvnorm(
    n = 1,
    mean = mu_a,
    sigma = sig.a * diag(q),
    upper = rep(0, q)
  )

  # Adaptive Metropolis-Hastings settings for u
  delta.v <- rep(3, j)
  adapt.iter <- 1000
  accept.count.u <- rep(0, j)
  proposal.count.u <- rep(0, j)
  target_acceptance_rate <- 0.23
  adapt_rate <- 0.05
  u_update_freq <- 5

  for (iter in 1:(N + nburn)) {

    # 1) Update f
    Cov.f <- tau * exp(-0.5 * rdist(vt)^2 / lscale)
    Cov.f <- (Cov.f + t(Cov.f)) / 2

    ch_Cov <- safe_chol(Cov.f)
    Cov.f_chol <- ch_Cov$R

    inv.Cov.f <- chol2inv(Cov.f_chol)
    inv.Cov.f <- (inv.Cov.f + t(inv.Cov.f)) / 2

    if (sum(gamma.v * u) <= 0) {
      stop(paste("Empty donor set in f update at iter", iter))
    }

    beta <- gamma.v * u / sum(gamma.v * u)

    Q.f <- inv.Cov.f + diag(1 / sigsq, t)
    Q.f <- (Q.f + t(Q.f)) / 2

    ch_Q <- safe_chol(Q.f)
    R_Q <- ch_Q$R

    b.f <- as.numeric(y - x %*% beta - Dt * alpha.v) / sigsq

    tmp <- forwardsolve(t(R_Q), b.f)
    mu.f <- backsolve(R_Q, tmp)

    z <- rnorm(t)
    f <- as.vector(mu.f + backsolve(R_Q, z))

    # 2) Update sigma^2
    beta <- rep(0, j)
    if (sum(gamma.v * u) != 0) {
      for (k in 1:j) {
        beta[k] <- gamma.v[k] * u[k] / sum(gamma.v * u)
      }
    } else {
      beta <- rep(0, j)
    }

    alpha_sig.post <- alpha_sig + t / 2
    beta_sig.post <- beta_sig +
      1 / 2 * sum((y - x %*% beta - Dt * alpha.v - f)^2)
    sigsq <- rinvgamma(1, alpha_sig.post, beta_sig.post)

    # 3) Update tau^2
    alpha_tau.post <- alpha_tau + t / 2
    Cov.l <- exp(-0.5 * rdist(vt)^2 / lscale)

    Cov.l <- (Cov.l + t(Cov.l)) / 2
    tryCatch(
      {
        Cov.l_chol <- chol(Cov.l)
      },
      error = function(cond) {
        Cov.l <- Cov.l +
          (-min(eigen(Cov.l)$values) + 0.001) * diag(ncol(Cov.l))
        Cov.l_chol <<- chol(Cov.l)
      }
    )
    inv.Cov.l <- chol2inv(Cov.l_chol)

    beta_tau.post <- beta_tau + 0.5 * t(f) %*% inv.Cov.l %*% f
    tau <- rinvgamma(1, alpha_tau.post, beta_tau.post)

    # 4) Update length-scale
    dataList.l <- list(
      tau = tau,
      vt = vt,
      f = f,
      alpha_l = alpha_l,
      beta_l = beta_l
    )
    delta <- 0.5
    initsList <- list(gp.l = lscale)
    lscale <- rw_mp.l(1, 0, delta, dataList.l, initsList)

    # 5) Update eta
    alpha_w.post <- sum(gamma.v) + 1
    beta_w.post <- sum(rep(1, j) - gamma.v) + 1
    w <- rbeta(1, alpha_w.post, beta_w.post)

    # 6) Update gamma_j
    for (p in 1:j) {
      gamma.tem1 <- gamma.v
      gamma.tem1[p] <- 1

      denom1 <- sum(gamma.tem1 * u)

      if (denom1 <= 0 || !is.finite(denom1)) {
        log_rj_1 <- -Inf
      } else {
        beta.tem1 <- gamma.tem1 * u / denom1
        log_rj_1 <- log(w) + loglik_gamma(
          beta = beta.tem1,
          y = y,
          x = x,
          Dt = Dt,
          alpha.v = alpha.v,
          f = f,
          sigsq = sigsq
        )
      }

      gamma.tem0 <- gamma.v
      gamma.tem0[p] <- 0

      denom0 <- sum(gamma.tem0 * u)

      if (denom0 <= 0 || !is.finite(denom0)) {
        log_rj_0 <- -Inf
      } else {
        beta.tem0 <- gamma.tem0 * u / denom0
        log_rj_0 <- log(1 - w) + loglik_gamma(
          beta = beta.tem0,
          y = y,
          x = x,
          Dt = Dt,
          alpha.v = alpha.v,
          f = f,
          sigsq = sigsq
        )
      }

      bin.p <- plogis(log_rj_1 - log_rj_0)

      if (!is.finite(bin.p) || bin.p < 0 || bin.p > 1) {
        print(list(
          p = p,
          log_rj_1 = log_rj_1,
          log_rj_0 = log_rj_0,
          sigsq = sigsq,
          sum_gamma_current = sum(gamma.v),
          denom1 = denom1,
          denom0 = denom0
        ))
        stop("Bad gamma probability")
      }

      gamma.v[p] <- rbinom(1, 1, bin.p)
    }

    # 7) Update u
    dataList.u <- list(
      gamma.v = gamma.v,
      alpha_u = alpha_u,
      Dt = Dt,
      alpha.v = alpha.v,
      f = f,
      x = x,
      y = y,
      sigsq = sigsq
    )

    initsList <- list(u = u)

    for (freq in 1:u_update_freq) {
      for (p in 1:j) {
        delta <- delta.v[p]

        u.curr <- initsList$u[p]
        u.prop <- initsList$u
        u.prop[p] <- exp(log(u.prop[p]) + rnorm(1, 0, 1) * delta)

        log_accept_prob <-
          log_post_u(u.prop, p, dataList.u) -
          log_post_u(initsList$u, p, dataList.u) +
          log(u.prop[p]) - log(u.curr)

        if (iter <= nburn) {
          proposal.count.u[p] <- proposal.count.u[p] + 1
        }

        if (is.finite(log_accept_prob) && log(runif(1)) < log_accept_prob) {
          initsList$u <- u.prop

          if (iter <= nburn) {
            accept.count.u[p] <- accept.count.u[p] + 1
          }
        }
      }
    }

    u <- initsList$u

    # 8) Update alpha
    W.inv <- sig.a^(-1) * diag(1, q) + sigsq^(-1) * t(Dt) %*% Dt
    W.a <- solve(W.inv)

    W.inv <- (W.inv + t(W.inv)) / 2
    tryCatch(
      {
        W.inv_chol <- chol(W.inv)
      },
      error = function(cond) {
        W.inv <- W.inv +
          (-min(eigen(W.inv)$values) + 0.001) * diag(ncol(W.inv))
        W.inv_chol <<- chol(W.inv)
      }
    )
    W.a <- chol2inv(W.inv_chol)

    beta <- rep(0, j)
    if (sum(gamma.v * u) != 0) {
      for (k in 1:j) {
        beta[k] <- gamma.v[k] * u[k] / sum(gamma.v * u)
      }
    } else {
      beta <- rep(0, j)
    }

    mu.a <- c(sigsq^(-1) * W.a %*% t(Dt) %*% (y - x %*% beta - f))
    alpha.v <- as.vector(rtmvnorm(1, mu.a, W.a, upper = rep(0, q)))

    # Diagnostic checks
    if (!is.finite(sigsq) || sigsq <= 0) stop(paste("Bad sigsq at iter", iter))
    if (!is.finite(tau) || tau <= 0) stop(paste("Bad tau at iter", iter))
    if (any(!is.finite(f))) stop(paste("Bad f at iter", iter))
    if (any(!is.finite(u)) || any(u <= 0)) stop(paste("Bad u at iter", iter))
    if (any(!is.finite(alpha.v))) stop(paste("Bad alpha.v at iter", iter))
    if (any(!is.finite(beta))) stop(paste("Bad beta at iter", iter))

    if (iter %% 100 == 0) {
      res <- as.numeric(y - x %*% beta - Dt * alpha.v - f)
      cat(
        "iter =", iter,
        "sigsq =", sigsq,
        "tau =", tau,
        "rss =", sum(res^2),
        "max|res| =", max(abs(res)),
        "sum_gamma =", sum(gamma.v),
        "\n"
      )
    }

    if (iter > nburn) {
      f.sample[iter - nburn, ] <- f
      sig.sample[iter - nburn] <- sigsq
      tau.sample[iter - nburn] <- tau
      l.sample[iter - nburn] <- lscale
      w.sample[iter - nburn] <- w
      gamma.sample[iter - nburn, ] <- gamma.v
      u.sample[iter - nburn, ] <- u
      a.sample[iter - nburn, ] <- alpha.v
    }

    if (iter <= nburn && iter %% adapt.iter == 0) {
      acceptance.rate <- accept.count.u / pmax(proposal.count.u, 1)

      print(acceptance.rate)
      print(delta.v)

      for (p in 1:j) {
        delta_min <- 0.05
        delta_max <- 6

        if (acceptance.rate[p] < target_acceptance_rate * 0.95) {
          delta.v[p] <- max(delta.v[p] * (1 - adapt_rate), delta_min)
        } else if (acceptance.rate[p] > target_acceptance_rate * 1.05) {
          delta.v[p] <- min(delta.v[p] * (1 + adapt_rate), delta_max)
        }
      }

      accept.count.u[] <- 0
      proposal.count.u[] <- 0
    }

    print(iter)
    print(lscale)
  }

  list(
    y = y,
    x = x,
    f.sample = f.sample,
    sig.sample = sig.sample,
    tau.sample = tau.sample,
    l.sample = l.sample,
    w.sample = w.sample,
    gamma.sample = gamma.sample,
    u.sample = u.sample,
    a.sample = a.sample
  )
}

# -----------------------------------------------------------------------------
# 5. Run or load MCMC results
# -----------------------------------------------------------------------------

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

if (run_mcmc) {
  results <- list()

  for (seed in seeds) {
    chain_id <- seed %/% 100
    results[[chain_id]] <- run_basc_chain(
      seed = seed,
      y = y,
      x = x,
      vt = vt,
      Dt = Dt,
      N = N,
      nburn = nburn,
      q = q
    )
  }

  save(results, file = result_file)
} else {
  load(result_file)
}

results <- compute_beta_samples(results, N = N, j = j)

# -----------------------------------------------------------------------------
# 6. Posterior summaries
# -----------------------------------------------------------------------------

selected_chain <- 2

beta.sample <- results[[selected_chain]]$beta.sample
gamma.sample <- results[[selected_chain]]$gamma.sample
u.sample <- results[[selected_chain]]$u.sample
sig.sample <- results[[selected_chain]]$sig.sample
a.sample <- results[[selected_chain]]$a.sample
f.sample <- results[[selected_chain]]$f.sample
y <- results[[selected_chain]]$y
x <- results[[selected_chain]]$x

beta_summary <- data.frame(
  donor = donor_names,
  beta_mean = colMeans(beta.sample),
  beta_sd = apply(beta.sample, 2, sd),
  gamma_mean = colMeans(gamma.sample)
)

print(round(beta_summary[, c("beta_mean", "beta_sd", "gamma_mean")], 2))
write.csv(
  beta_summary,
  file = file.path(output_dir, "BASC_realdata_weight_summary.csv"),
  row.names = FALSE
)

# Posterior predictive counterfactual path.
eps.sample <- matrix(0, N, t)
for (i in 1:N) {
  eps.sample[i, ] <- rnorm(t, mean = 0, sd = sig.sample[i])
}

y.sample <- matrix(0, N, t)
for (i in 1:N) {
  y.sample[i, ] <- x %*% beta.sample[i, ] + Dt * 0 + f.sample[i, ]
}

y.ci <- apply(y.sample, 2, quantile, probs = c(0.025, 0.5, 0.975))

y.predict <- data.frame(
  vt = vt,
  y_mean = y.ci[2, ],
  ymin = y.ci[1, ],
  ymax = y.ci[3, ]
)

y.true <- data.frame(vt = vt, y = y)

save(y.predict, file = file.path(output_dir, "BASC_realdata_prediction.RData"))
write.csv(
  y.predict,
  file = file.path(output_dir, "BASC_realdata_prediction_95ci.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 7. Convergence diagnostics
# -----------------------------------------------------------------------------

print(1 - rejectionRate(mcmc(results[[selected_chain]]$l.sample)))
print(1 - rejectionRate(mcmc(results[[selected_chain]]$u.sample)))

sample.names <- c(
  "f.sample",
  "sig.sample",
  "tau.sample",
  "l.sample",
  "w.sample",
  "u.sample",
  "a.sample"
)

gelman.summary <- list()

for (sample.name in sample.names) {
  mcmc.chain1 <- results[[1]][[sample.name]]
  mcmc.chain2 <- results[[2]][[sample.name]]
  mcmc.chain3 <- results[[3]][[sample.name]]

  thin <- 1

  if (is.matrix(mcmc.chain1)) {
    new.chain1 <- mcmc(mcmc.chain1[seq(1, nrow(mcmc.chain1), by = thin), ])
    new.chain2 <- mcmc(mcmc.chain2[seq(1, nrow(mcmc.chain2), by = thin), ])
    new.chain3 <- mcmc(mcmc.chain3[seq(1, nrow(mcmc.chain3), by = thin), ])
  } else {
    new.chain1 <- mcmc(mcmc.chain1[seq(1, length(mcmc.chain1), by = thin)])
    new.chain2 <- mcmc(mcmc.chain2[seq(1, length(mcmc.chain2), by = thin)])
    new.chain3 <- mcmc(mcmc.chain3[seq(1, length(mcmc.chain3), by = thin)])
  }

  mcmc.list <- mcmc.list(new.chain1, new.chain2, new.chain3)
  gelman.summary[[sample.name]] <- gelman.diag(mcmc.list)
}

print(gelman.summary)

beta.gelman <- list()
for (k in 1:j) {
  beta.gelman[[k]] <- gelman.diag(mcmc.list(
    mcmc(results[[1]]$beta.sample[, k]),
    mcmc(results[[2]]$beta.sample[, k]),
    mcmc(results[[3]]$beta.sample[, k])
  ))

  cat("\nBeta", k, "\n")
  print(beta.gelman[[k]])
}

# -----------------------------------------------------------------------------
# 8. Plotting functions
# -----------------------------------------------------------------------------

plot_counterfactual_path <- function(y.predict, y.true) {
  ggplot() +
    geom_ribbon(
      data = y.predict,
      aes(x = vt, ymin = ymin, ymax = ymax),
      fill = "skyblue",
      alpha = 0.3
    ) +
    geom_line(
      data = y.predict,
      aes(x = vt, y = y_mean),
      color = "steelblue",
      linewidth = 1
    ) +
    geom_point(
      data = y.predict,
      aes(x = vt, y = y_mean),
      color = "steelblue",
      size = 2
    ) +
    geom_line(
      data = y.true,
      aes(x = vt, y = y),
      color = "black",
      linewidth = 0.5
    ) +
    scale_y_continuous(limits = c(0, 35000)) +
    labs(x = "time", y = "y") +
    theme_minimal()
}

plot_all_traces <- function(results,
                            thin = 100,
                            selected_u = c(1, 10, 2, 7),
                            selected_f = c(1, 10, 40),
                            output_dir = "plot/") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  trace_list <- list()

  scalar_params <- c("sig.sample", "tau.sample", "l.sample", "w.sample", "a.sample")

  for (param in scalar_params) {
    chain_data <- lapply(seq_along(results), function(chain_id) {
      iter_idx <- seq(1, nrow(as.matrix(results[[chain_id]][[param]])), by = thin)
      data.frame(
        iter = iter_idx,
        value = results[[chain_id]][[param]][iter_idx],
        chain = as.factor(chain_id),
        parameter = param
      )
    })

    trace_list[[param]] <- do.call(rbind, chain_data)
  }

  if ("u.sample" %in% names(results[[1]])) {
    chain_data <- lapply(seq_along(results), function(chain_id) {
      mat <- results[[chain_id]][["u.sample"]]
      iter_idx <- seq(1, nrow(mat), by = thin)
      available_u <- ncol(mat)
      sel_u <- selected_u[selected_u <= available_u]

      df <- as.data.frame(mat[iter_idx, sel_u, drop = FALSE])
      colnames(df) <- sapply(sel_u, ordinal_suffix)
      df$iter <- iter_idx
      df$chain <- as.factor(chain_id)

      df_melt <- melt(
        df,
        id.vars = c("iter", "chain"),
        variable.name = "sub_index",
        value.name = "value"
      )

      df_melt$sub_index <- factor(
        df_melt$sub_index,
        levels = sapply(sel_u, ordinal_suffix)
      )
      df_melt$parameter <- paste0("u (", df_melt$sub_index, " control unit)")

      df_melt[, c("iter", "value", "chain", "parameter")]
    })

    trace_list[["u.sample"]] <- do.call(rbind, chain_data)
  }

  if ("f.sample" %in% names(results[[1]])) {
    chain_data <- lapply(seq_along(results), function(chain_id) {
      mat <- results[[chain_id]][["f.sample"]]
      iter_idx <- seq(1, nrow(mat), by = thin)
      available_f <- ncol(mat)
      sel_f <- selected_f[selected_f <= available_f]

      df <- as.data.frame(mat[iter_idx, sel_f, drop = FALSE])
      colnames(df) <- as.character(sel_f)
      df$iter <- iter_idx
      df$chain <- as.factor(chain_id)

      df_melt <- melt(
        df,
        id.vars = c("iter", "chain"),
        variable.name = "sub_index",
        value.name = "value"
      )

      df_melt$parameter <- paste0("f (t=", df_melt$sub_index, ")")
      df_melt[, c("iter", "value", "chain", "parameter")]
    })

    trace_list[["f.sample"]] <- do.call(rbind, chain_data)
  }

  all_trace <- do.call(rbind, trace_list)

  all_trace$parameter <- gsub("f.sample", "f", all_trace$parameter)
  all_trace$parameter <- gsub("u.sample", "u", all_trace$parameter)
  all_trace$parameter <- gsub("sig.sample", "sigma^2", all_trace$parameter)
  all_trace$parameter <- gsub("w.sample", "eta", all_trace$parameter)
  all_trace$parameter <- gsub("a.sample", "alpha", all_trace$parameter)
  all_trace$parameter <- gsub("tau.sample", "tau", all_trace$parameter)
  all_trace$parameter <- gsub("l.sample", "kappa", all_trace$parameter)

  final_levels <- c(
    "f (t=1)", "f (t=10)", "f (t=40)",
    "tau", "eta", "alpha", "sigma^2", "kappa",
    "u (1st control unit)",
    "u (10th control unit)",
    "u (2nd control unit)",
    "u (7th control unit)"
  )

  final_levels <- intersect(final_levels, unique(all_trace$parameter))
  all_trace$parameter <- factor(all_trace$parameter, levels = final_levels)

  scalar_labels <- c("sigma^2", "tau", "eta", "alpha", "kappa")
  custom_labeller <- function(x) {
    x <- as.character(x)
    lapply(x, function(lbl) {
      if (lbl %in% scalar_labels) {
        parse(text = lbl)[[1]]
      } else {
        lbl
      }
    })
  }

  n_facets <- length(unique(all_trace$parameter))
  n_rows <- ceiling(n_facets / 3)
  total_height <- n_rows * 2

  p <- ggplot(all_trace, aes(x = iter, y = value, color = chain)) +
    geom_line(alpha = 0.7) +
    scale_x_continuous(
      breaks = c(100000, 300000, 500000),
      labels = function(x) format(x, big.mark = ",", scientific = FALSE)
    ) +
    facet_wrap(
      ~ parameter,
      scales = "free_y",
      ncol = 3,
      labeller = as_labeller(custom_labeller)
    ) +
    labs(title = "", x = "Iteration", y = "Value") +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14),
      strip.text = element_text(size = 18, face = "bold")
    )

  filename <- file.path(output_dir, "trace_plots_realdata")

  ggsave(paste0(filename, ".png"), plot = p, width = 12, height = total_height, dpi = 150)
  ggsave(paste0(filename, ".pdf"), plot = p, width = 12, height = total_height)
  ggsave(
    paste0(filename, ".eps"),
    plot = p,
    width = 12,
    height = total_height,
    device = cairo_ps,
    fallback_resolution = 600
  )

  p
}

plot_all_densities <- function(results,
                               thin = 100,
                               selected_u = c(1, 10, 2, 7),
                               selected_f = c(1, 10, 25, 50),
                               output_dir = "plot/") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  density_list <- list()

  scalar_params <- c("sig.sample", "tau.sample", "l.sample", "w.sample", "a.sample")

  for (param in scalar_params) {
    chain_data <- lapply(seq_along(results), function(chain_id) {
      param_vals <- results[[chain_id]][[param]]
      thin_idx <- seq(1, length(param_vals), by = thin)

      data.frame(
        value = as.vector(param_vals[thin_idx]),
        chain = as.factor(chain_id),
        parameter = param
      )
    })

    density_list[[param]] <- do.call(rbind, chain_data)
  }

  if ("u.sample" %in% names(results[[1]])) {
    chain_data <- lapply(seq_along(results), function(chain_id) {
      mat <- results[[chain_id]][["u.sample"]]
      available_u <- ncol(mat)
      sel_u <- selected_u[selected_u <= available_u]
      thin_idx <- seq(1, nrow(mat), by = thin)

      df_u <- do.call(rbind, lapply(seq_along(sel_u), function(i) {
        data.frame(
          value = mat[thin_idx, sel_u[i]],
          chain = as.factor(chain_id),
          sub_index = ordinal_suffix(sel_u[i])
        )
      }))

      df_u$parameter <- paste0("u (", df_u$sub_index, " control unit)")
      df_u[, c("value", "chain", "parameter")]
    })

    density_list[["u.sample"]] <- do.call(rbind, chain_data)
  }

  if ("f.sample" %in% names(results[[1]])) {
    chain_data <- lapply(seq_along(results), function(chain_id) {
      mat <- results[[chain_id]][["f.sample"]]
      available_f <- ncol(mat)
      sel_f <- selected_f[selected_f <= available_f]
      thin_idx <- seq(1, nrow(mat), by = thin)

      df_f <- do.call(rbind, lapply(seq_along(sel_f), function(i) {
        data.frame(
          value = mat[thin_idx, sel_f[i]],
          chain = as.factor(chain_id),
          sub_index = as.character(sel_f[i])
        )
      }))

      df_f$parameter <- paste0("f (t=", df_f$sub_index, ")")
      df_f[, c("value", "chain", "parameter")]
    })

    density_list[["f.sample"]] <- do.call(rbind, chain_data)
  }

  all_density <- do.call(rbind, density_list)

  all_density$parameter <- gsub("sig.sample", "sigma^2", all_density$parameter)
  all_density$parameter <- gsub("tau.sample", "tau", all_density$parameter)
  all_density$parameter <- gsub("l.sample", "kappa", all_density$parameter)
  all_density$parameter <- gsub("w.sample", "eta", all_density$parameter)
  all_density$parameter <- gsub("a.sample", "alpha", all_density$parameter)

  final_levels <- c(
    "f (t=1)", "f (t=10)", "f (t=25)", "f (t=50)",
    "tau", "eta", "alpha", "sigma^2", "kappa",
    "u (1st control unit)",
    "u (10th control unit)",
    "u (2nd control unit)",
    "u (7th control unit)"
  )

  final_levels <- intersect(final_levels, unique(all_density$parameter))
  all_density$parameter <- factor(all_density$parameter, levels = final_levels)

  n_facets <- length(unique(all_density$parameter))
  n_rows <- ceiling(n_facets / 3)
  total_height <- n_rows * 2

  p <- ggplot(all_density, aes(x = value, fill = chain)) +
    geom_density(alpha = 0.5) +
    facet_wrap(~ parameter, scales = "free", ncol = 3) +
    labs(title = "", x = "Value", y = "Density") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 20, face = "bold"),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 16),
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 16),
      strip.text = element_text(size = 18, face = "bold")
    )

  filename <- file.path(output_dir, "density_plots_realdata")

  ggsave(paste0(filename, ".png"), plot = p, width = 12, height = total_height, dpi = 150)
  ggsave(paste0(filename, ".pdf"), plot = p, width = 12, height = total_height)
  ggsave(
    paste0(filename, ".eps"),
    plot = p,
    width = 12,
    height = total_height,
    device = cairo_ps,
    fallback_resolution = 600
  )

  p
}

plot_gamma_inclusion <- function(results, donor_names, selected_chain = 2, output_dir = "plot/") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  gamma_df <- data.frame(
    donor = factor(donor_names, levels = donor_names),
    gamma_mean = colMeans(results[[selected_chain]]$gamma.sample)
  )

  p <- ggplot(gamma_df, aes(x = donor, y = gamma_mean)) +
    geom_col(fill = "skyblue") +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
    labs(title = "", x = "donor pool", y = expression(hat(gamma))) +
    theme_minimal(base_size = 16) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(
    file.path(output_dir, "realdata_gamma.png"),
    plot = p,
    width = 12,
    height = 4,
    dpi = 150
  )
  ggsave(
    file.path(output_dir, "realdata_gamma.eps"),
    plot = p,
    width = 12,
    height = 4
  )
  ggsave(
    file.path(output_dir, "realdata_gamma.pdf"),
    plot = p,
    width = 12,
    height = 4
  )

  p
}

# -----------------------------------------------------------------------------
# 9. Create plots
# -----------------------------------------------------------------------------

p_path <- plot_counterfactual_path(y.predict, y.true)
ggsave(
  file.path(plot_dir, "BASC_realdata_counterfactual.png"),
  plot = p_path,
  width = 7,
  height = 5,
  dpi = 150
)

plot_all_traces(results, thin = 100, output_dir = plot_dir)
plot_all_densities(results, thin = 1, output_dir = plot_dir)
plot_gamma_inclusion(results, donor_names, selected_chain = selected_chain, output_dir = plot_dir)
