###############################################################################
# Bayesian Adaptive Synthetic Control (BASC-pre)
# West Germany application: pre-intervention fitting
#
# This script implements the BASC-pre specification for the West Germany application.
# All model parameters are estimated using 1960--1989 only.
# The post-intervention latent GP component is generated conditionally on the
# pre-intervention GP draws.
###############################################################################

rm(list = ls())

# -----------------------------------------------------------------------------
# 0. Packages
# -----------------------------------------------------------------------------

required_packages <- c(
  "foreign",
  "dplyr",
  "tidyr",
  "MASS",
  "MCMCpack",
  "fields",
  "coda"
)

invisible(lapply(required_packages, require, character.only = TRUE))

# -----------------------------------------------------------------------------
# 1. User settings
# -----------------------------------------------------------------------------

input_data_path <- "data/repgermany.dta"
output_dir <- "results/realdata/BASC_pre"

N <- 500000
nburn <- 500000
seeds <- c(103, 200, 305)

# Three chains are run for convergence assessment.
# Posterior summaries use chain 2, matching the reported analysis.
selected_chain <- 2

run_mcmc <- TRUE
result_file <- file.path(output_dir, "BASC_pre_results.RData")

# Seed used for conditional GP generation in the post-intervention period.
prediction_seed <- 20261006

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 2. Data preparation
# -----------------------------------------------------------------------------

prepare_realdata <- function(input_data_path) {
  data <- foreign::read.dta(input_data_path)

  data_y <- data %>%
    dplyr::filter(country == "West Germany") %>%
    dplyr::arrange(year)

  y <- data_y$gdp
  year_y <- data_y$year
  vt <- year_y - 1959

  base_x <- data %>%
    dplyr::filter(country != "West Germany") %>%
    dplyr::select(country, year, gdp) %>%
    tidyr::pivot_wider(names_from = year, values_from = gdp)

  x <- base_x %>%
    dplyr::select(-country) %>%
    as.matrix() %>%
    t()

  list(
    y = y,
    x = x,
    year = year_y,
    vt = vt,
    donor_names = base_x$country
  )
}

realdata <- prepare_realdata(input_data_path)

y <- realdata$y
x <- realdata$x
year_y <- realdata$year
vt <- realdata$vt
donor_names <- realdata$donor_names

pre <- which(year_y < 1990)
post <- which(year_y >= 1990)

y_fit <- y[pre]
x_fit <- x[pre, , drop = FALSE]
vt_fit <- vt[pre]

j <- ncol(x)
TT <- length(vt)
TT_fit <- length(vt_fit)

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

loglik_gamma <- function(beta, y, x, f, sigsq) {
  res <- as.numeric(y - x %*% beta - f)

  -0.5 * length(y) * log(2 * pi * sigsq) -
    0.5 * sum(res^2) / sigsq
}

log_post_l <- function(gp_l, data_list) {
  tau <- data_list$tau
  vtt <- data_list$vt
  ff <- as.numeric(data_list$f)
  alpha_l <- data_list$alpha_l
  beta_l <- data_list$beta_l

  if (!is.finite(gp_l) || gp_l <= 0) return(-Inf)
  if (!is.finite(tau) || tau <= 0) return(-Inf)
  if (any(!is.finite(ff))) return(-Inf)

  gp_cov <- tau * exp(-0.5 * fields::rdist(vtt)^2 / gp_l)
  gp_cov <- (gp_cov + t(gp_cov)) / 2

  ch <- try(safe_chol(gp_cov), silent = TRUE)
  if (inherits(ch, "try-error")) return(-Inf)

  R <- ch$R
  z <- backsolve(R, ff, transpose = TRUE)

  quad <- sum(z^2)
  logdet <- 2 * sum(log(diag(R)))

  log_prior <- -(alpha_l + 1) * log(gp_l) - beta_l / gp_l
  log_gp <- -0.5 * logdet - 0.5 * quad

  log_prior + log_gp
}

rw_mp_l <- function(nsim_l, nburn_l, delta, data_list, inits_list) {
  log_curr <- log(inits_list$gp_l)
  para_samples <- rep(NA_real_, nsim_l)
  acc_count <- 0

  for (iter_l in 1:(nsim_l + nburn_l)) {
    curr <- exp(log_curr)
    log_prop <- log_curr + rnorm(1, 0, delta)
    prop <- exp(log_prop)

    log_acc <- log_post_l(prop, data_list) -
      log_post_l(curr, data_list) +
      log_prop - log_curr

    if (is.finite(log_acc) && log(runif(1)) < log_acc) {
      log_curr <- log_prop

      if (iter_l > nburn_l) {
        acc_count <- acc_count + 1
      }
    }

    if (iter_l > nburn_l) {
      para_samples[iter_l - nburn_l] <- exp(log_curr)
    }
  }

  attr(para_samples, "accept_rate") <- acc_count / nsim_l
  para_samples
}

log_post_u <- function(u, p, data_list) {
  gamma_v <- data_list$gamma_v
  alpha_u <- data_list$alpha_u
  f <- data_list$f
  x <- data_list$x
  y <- data_list$y
  sigsq <- data_list$sigsq

  u_j <- u[p]
  denom <- sum(gamma_v * u)

  if (!is.finite(u_j) || u_j <= 0) return(-Inf)
  if (!is.finite(denom) || denom <= 0) return(-Inf)
  if (!is.finite(sigsq) || sigsq <= 0) return(-Inf)

  beta_tem <- gamma_v * u / denom
  res <- as.numeric(y - x %*% beta_tem - f)

  log_prior <- (alpha_u - 1) * log(u_j) - alpha_u * u_j
  log_lik <- -0.5 * sum(res^2) / sigsq

  log_prior + log_lik
}

compute_beta_samples <- function(results, N, j) {
  for (chain_id in seq_along(results)) {
    results[[chain_id]]$beta.sample <- matrix(0, N, j)

    for (beta_row in 1:N) {
      denom <- sum(
        results[[chain_id]]$gamma.sample[beta_row, ] *
          results[[chain_id]]$u.sample[beta_row, ]
      )

      if (!is.finite(denom) || denom <= 0) {
        stop("Empty donor set encountered while constructing posterior weights.")
      }

      results[[chain_id]]$beta.sample[beta_row, ] <-
        results[[chain_id]]$gamma.sample[beta_row, ] *
        results[[chain_id]]$u.sample[beta_row, ] / denom
    }
  }

  results
}

# -----------------------------------------------------------------------------
# 4. MCMC sampler
# -----------------------------------------------------------------------------

run_basc_pre_chain <- function(seed) {
  # Storage
  f_sample <- matrix(0, N, TT_fit)
  sig_sample <- rep(0, N)
  tau_sample <- rep(0, N)
  l_sample <- rep(0, N)
  eta_sample <- rep(0, N)
  gamma_sample <- matrix(0, N, j)
  u_sample <- matrix(0, N, j)

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

  tau <- MCMCpack::rinvgamma(n = 1, alpha_tau, beta_tau)
  lscale <- MCMCpack::rinvgamma(n = 1, alpha_l, beta_l)

  Cov_l <- exp(-0.5 * fields::rdist(vt_fit)^2 / lscale)
  f <- as.vector(
    MASS::mvrnorm(1, rep(0, TT_fit), tau * Cov_l)
  )

  eta <- runif(1, 0, 1)

  gamma_v <- rep(1, j)
  u <- rgamma(j, alpha_u, scale = 1 / alpha_u)
  beta <- gamma_v * u / sum(gamma_v * u)

  sigsq <- MCMCpack::rinvgamma(1, alpha_sig, beta_sig)

  # Adaptive MH settings for u
  delta_v <- rep(3, j)
  adapt_iter <- 1000
  accept_count_u <- rep(0, j)
  proposal_count_u <- rep(0, j)

  target_acceptance_rate <- 0.23
  adapt_rate <- 0.05
  u_update_freq <- 5

  for (iter in 1:(N + nburn)) {
    # 1) Update f
    Cov_f <- tau * exp(-0.5 * fields::rdist(vt_fit)^2 / lscale)
    Cov_f <- (Cov_f + t(Cov_f)) / 2

    Cov_f_chol <- safe_chol(Cov_f)$R
    inv_Cov_f <- chol2inv(Cov_f_chol)
    inv_Cov_f <- (inv_Cov_f + t(inv_Cov_f)) / 2

    if (sum(gamma_v * u) <= 0) {
      stop(paste("Empty donor set in f update at iteration", iter))
    }

    beta <- gamma_v * u / sum(gamma_v * u)

    Q_f <- inv_Cov_f + diag(1 / sigsq, TT_fit)
    Q_f <- (Q_f + t(Q_f)) / 2
    R_Q <- safe_chol(Q_f)$R

    b_f <- as.numeric(y_fit - x_fit %*% beta) / sigsq
    tmp <- forwardsolve(t(R_Q), b_f)
    mu_f <- backsolve(R_Q, tmp)

    z <- rnorm(TT_fit)
    f <- as.vector(mu_f + backsolve(R_Q, z))

    # 2) Update sigma^2
    beta <- gamma_v * u / sum(gamma_v * u)

    alpha_sig_post <- alpha_sig + TT_fit / 2
    beta_sig_post <- beta_sig +
      0.5 * sum((y_fit - x_fit %*% beta - f)^2)

    sigsq <- MCMCpack::rinvgamma(1, alpha_sig_post, beta_sig_post)

    # 3) Update tau
    alpha_tau_post <- alpha_tau + TT_fit / 2
    Cov_l <- exp(-0.5 * fields::rdist(vt_fit)^2 / lscale)
    Cov_l <- (Cov_l + t(Cov_l)) / 2

    tryCatch(
      {
        Cov_l_chol <- chol(Cov_l)
      },
      error = function(cond) {
        Cov_l <- Cov_l +
          (-min(eigen(Cov_l)$values) + 0.001) * diag(ncol(Cov_l))
        Cov_l_chol <<- chol(Cov_l)
      }
    )

    inv_Cov_l <- chol2inv(Cov_l_chol)

    beta_tau_post <- beta_tau + 0.5 * t(f) %*% inv_Cov_l %*% f
    tau <- MCMCpack::rinvgamma(
      1,
      alpha_tau_post,
      as.numeric(beta_tau_post)
    )

    # 4) Update GP length scale
    data_list_l <- list(
      tau = tau,
      vt = vt_fit,
      f = f,
      alpha_l = alpha_l,
      beta_l = beta_l
    )

    lscale <- rw_mp_l(
      nsim_l = 1,
      nburn_l = 0,
      delta = 0.6,
      data_list = data_list_l,
      inits_list = list(gp_l = lscale)
    )

    # 5) Update eta
    alpha_eta_post <- sum(gamma_v) + 1
    beta_eta_post <- sum(1 - gamma_v) + 1
    eta <- rbeta(1, alpha_eta_post, beta_eta_post)

    # 6) Update gamma_j
    for (p in 1:j) {
      gamma_tem1 <- gamma_v
      gamma_tem1[p] <- 1

      denom1 <- sum(gamma_tem1 * u)

      if (denom1 <= 0 || !is.finite(denom1)) {
        log_rj_1 <- -Inf
      } else {
        beta_tem1 <- gamma_tem1 * u / denom1

        log_rj_1 <- log(eta) + loglik_gamma(
          beta = beta_tem1,
          y = y_fit,
          x = x_fit,
          f = f,
          sigsq = sigsq
        )
      }

      gamma_tem0 <- gamma_v
      gamma_tem0[p] <- 0

      denom0 <- sum(gamma_tem0 * u)

      if (denom0 <= 0 || !is.finite(denom0)) {
        # The empty donor set is not allowed.
        log_rj_0 <- -Inf
      } else {
        beta_tem0 <- gamma_tem0 * u / denom0

        log_rj_0 <- log(1 - eta) + loglik_gamma(
          beta = beta_tem0,
          y = y_fit,
          x = x_fit,
          f = f,
          sigsq = sigsq
        )
      }

      bin_p <- plogis(log_rj_1 - log_rj_0)

      if (!is.finite(bin_p) || bin_p < 0 || bin_p > 1) {
        stop("Invalid donor-inclusion probability.")
      }

      gamma_v[p] <- rbinom(1, 1, bin_p)
    }

    if (sum(gamma_v) == 0) {
      stop(paste("Empty donor set after gamma update at iteration", iter))
    }

    # 7) Update u
    data_list_u <- list(
      gamma_v = gamma_v,
      alpha_u = alpha_u,
      f = f,
      x = x_fit,
      y = y_fit,
      sigsq = sigsq
    )

    inits_u <- u

    for (freq in 1:u_update_freq) {
      for (p in 1:j) {
        delta <- delta_v[p]

        u_curr <- inits_u[p]
        u_prop <- inits_u
        u_prop[p] <- exp(log(u_prop[p]) + rnorm(1, 0, 1) * delta)

        log_accept_prob <-
          log_post_u(u_prop, p, data_list_u) -
          log_post_u(inits_u, p, data_list_u) +
          log(u_prop[p]) - log(u_curr)

        if (iter <= nburn) {
          proposal_count_u[p] <- proposal_count_u[p] + 1
        }

        if (is.finite(log_accept_prob) &&
            log(runif(1)) < log_accept_prob) {
          inits_u <- u_prop

          if (iter <= nburn) {
            accept_count_u[p] <- accept_count_u[p] + 1
          }
        }
      }
    }

    u <- inits_u

    # Basic numerical checks
    if (!is.finite(sigsq) || sigsq <= 0) stop("Invalid sigma^2.")
    if (!is.finite(tau) || tau <= 0) stop("Invalid tau.")
    if (any(!is.finite(f))) stop("Invalid f.")
    if (any(!is.finite(u)) || any(u <= 0)) stop("Invalid u.")
    if (any(!is.finite(beta))) stop("Invalid donor weights.")

    # Store retained draws
    if (iter > nburn) {
      s <- iter - nburn

      f_sample[s, ] <- f
      sig_sample[s] <- sigsq
      tau_sample[s] <- tau
      l_sample[s] <- lscale
      eta_sample[s] <- eta
      gamma_sample[s, ] <- gamma_v
      u_sample[s, ] <- u
    }

    # Adapt u proposal scales during burn-in
    if (iter <= nburn && iter %% adapt_iter == 0) {
      acceptance_rate <- accept_count_u / pmax(proposal_count_u, 1)

      for (p in 1:j) {
        delta_min <- 0.05
        delta_max <- 6

        if (acceptance_rate[p] < target_acceptance_rate * 0.95) {
          delta_v[p] <- max(delta_v[p] * (1 - adapt_rate), delta_min)
        } else if (acceptance_rate[p] > target_acceptance_rate * 1.05) {
          delta_v[p] <- min(delta_v[p] * (1 + adapt_rate), delta_max)
        }
      }

      accept_count_u[] <- 0
      proposal_count_u[] <- 0
    }

    if (iter %% 10000 == 0) {
      cat("seed =", seed,
          "| iteration =", iter, "/", N + nburn, "\n")
    }
  }

  list(
    y = y,
    x = x,
    f.sample = f_sample,
    sig.sample = sig_sample,
    tau.sample = tau_sample,
    l.sample = l_sample,
    eta.sample = eta_sample,
    gamma.sample = gamma_sample,
    u.sample = u_sample
  )
}

# -----------------------------------------------------------------------------
# 5. Run or load MCMC results
# -----------------------------------------------------------------------------

if (run_mcmc) {
  results <- vector("list", length(seeds))

  for (i in seq_along(seeds)) {
    results[[i]] <- run_basc_pre_chain(seeds[i])
  }

  save(results, file = result_file)
} else {
  load(result_file)
}

results <- compute_beta_samples(results, N = N, j = j)

res <- results[[selected_chain]]

# -----------------------------------------------------------------------------
# 6. Conditional post-intervention GP draws
# -----------------------------------------------------------------------------

S <- nrow(res$f.sample)
pred_idx <- seq_len(S)

beta_pred <- res$beta.sample[pred_idx, , drop = FALSE]
f_pre_pred <- res$f.sample[pred_idx, , drop = FALSE]
tau_pred <- res$tau.sample[pred_idx]
l_pred <- res$l.sample[pred_idx]

vt_pre <- vt[pre]
vt_post <- vt[post]

set.seed(prediction_seed)

f_post_sample <- matrix(
  0,
  nrow = S,
  ncol = length(post)
)

for (s in 1:S) {
  ell <- l_pred[s]
  tau <- tau_pred[s]

  Rmm <- exp(
    -0.5 * outer(vt_pre, vt_pre, "-")^2 / ell
  )

  Rpm <- exp(
    -0.5 * outer(vt_post, vt_pre, "-")^2 / ell
  )

  Rpp <- exp(
    -0.5 * outer(vt_post, vt_post, "-")^2 / ell
  )

  Rmm <- (Rmm + t(Rmm)) / 2
  R <- safe_chol(Rmm)$R

  tmp <- forwardsolve(t(R), f_pre_pred[s, ])
  Rmm_inv_f <- backsolve(R, tmp)

  mu_post <- as.vector(
    Rpm %*% Rmm_inv_f
  )

  tmp2 <- forwardsolve(
    t(R),
    t(Rpm)
  )

  Rmm_inv_Rmp <- backsolve(
    R,
    tmp2
  )

  V_post <- tau * (
    Rpp -
      Rpm %*% Rmm_inv_Rmp
  )

  V_post <- (V_post + t(V_post)) / 2
  ch_V <- safe_chol(V_post)

  z <- rnorm(length(post))

  f_post_sample[s, ] <- as.vector(
    mu_post +
      t(ch_V$R) %*% z
  )
}

# -----------------------------------------------------------------------------
# 7. Counterfactual trajectory and summaries
# -----------------------------------------------------------------------------

x_pre <- x[pre, , drop = FALSE]
x_post <- x[post, , drop = FALSE]

y0_pre_sample <-
  beta_pred %*% t(x_pre) +
  f_pre_pred

y0_post_sample <-
  beta_pred %*% t(x_post) +
  f_post_sample

y0_sample <- matrix(
  NA_real_,
  nrow = S,
  ncol = TT
)

y0_sample[, pre] <- y0_pre_sample
y0_sample[, post] <- y0_post_sample

y0_quantile <- apply(
  y0_sample,
  2,
  quantile,
  probs = c(0.025, 0.5, 0.975)
)

counterfactual_summary <- data.frame(
  year = year_y,
  observed = y,
  posterior_mean = colMeans(y0_sample),
  posterior_median = y0_quantile[2, ],
  lower_95 = y0_quantile[1, ],
  upper_95 = y0_quantile[3, ]
)

write.csv(
  counterfactual_summary,
  file = file.path(output_dir, "BASC_pre_counterfactual.csv"),
  row.names = FALSE
)

beta_summary <- data.frame(
  donor = donor_names,
  beta_mean = colMeans(res$beta.sample),
  beta_sd = apply(res$beta.sample, 2, sd),
  gamma_mean = colMeans(res$gamma.sample)
)

write.csv(
  beta_summary,
  file = file.path(output_dir, "BASC_pre_weight_summary.csv"),
  row.names = FALSE
)

# Pre-intervention RMSE
y0_pre_mean <- colMeans(y0_pre_sample)

pre_rmse <- sqrt(
  mean((y0_pre_mean - y[pre])^2)
)

# Post-intervention effect and ATT
effect_sample <- matrix(
  y,
  nrow = S,
  ncol = length(y),
  byrow = TRUE
) - y0_sample

effect_post_sample <- effect_sample[, post, drop = FALSE]

ATT_sample <- rowMeans(effect_post_sample)
ATT_mean <- mean(ATT_sample)

ATT_ci <- quantile(
  ATT_sample,
  c(0.025, 0.5, 0.975)
)

# Donor-comparability diagnostic associated with interpolation concerns
donor_mse <- colMeans(
  (x_pre - y[pre])^2
)

beta_mean <- colMeans(res$beta.sample)

donor_discrepancy <- sqrt(
  sum(beta_mean * donor_mse)
)

summary_measures <- data.frame(
  quantity = c(
    "Pre-RMSE",
    "ATT posterior mean",
    "Weighted pre-intervention donor discrepancy"
  ),
  value = c(
    pre_rmse,
    ATT_mean,
    donor_discrepancy
  )
)

write.csv(
  summary_measures,
  file = file.path(output_dir, "BASC_pre_summary_measures.csv"),
  row.names = FALSE
)

print(beta_summary)
print(summary_measures)

# -----------------------------------------------------------------------------
# 8. Convergence diagnostics
# -----------------------------------------------------------------------------

sample_names <- c(
  "f.sample",
  "sig.sample",
  "tau.sample",
  "l.sample",
  "eta.sample",
  "u.sample"
)

gelman_summary <- list()

for (sample_name in sample_names) {
  chain1 <- coda::mcmc(results[[1]][[sample_name]])
  chain2 <- coda::mcmc(results[[2]][[sample_name]])
  chain3 <- coda::mcmc(results[[3]][[sample_name]])

  gelman_summary[[sample_name]] <- coda::gelman.diag(
    coda::mcmc.list(chain1, chain2, chain3)
  )
}

beta_gelman <- vector("list", j)

for (k in 1:j) {
  beta_gelman[[k]] <- coda::gelman.diag(
    coda::mcmc.list(
      coda::mcmc(results[[1]]$beta.sample[, k]),
      coda::mcmc(results[[2]]$beta.sample[, k]),
      coda::mcmc(results[[3]]$beta.sample[, k])
    )
  )
}

save(
  gelman_summary,
  beta_gelman,
  file = file.path(output_dir, "BASC_pre_convergence.RData")
)

print(gelman_summary)
