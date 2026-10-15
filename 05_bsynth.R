###############################################################################
# Bayesian synthetic control using the official bsynth package
# West Germany application
#
# Outcome-only specification with predictor matching OFF and Gaussian process OFF.
###############################################################################

rm(list = ls())

# -----------------------------------------------------------------------------
# 0. Packages
# -----------------------------------------------------------------------------

required_packages <- c(
  "foreign",
  "bsynth",
  "rstan"
)

invisible(lapply(required_packages, require, character.only = TRUE))

# -----------------------------------------------------------------------------
# 1. User settings
# -----------------------------------------------------------------------------

input_data_path <- "data/repgermany.dta"
output_dir <- "results/realdata/bsynth"

pre_years <- 1960:1989
post_years <- 1990:2003
all_years <- 1960:2003

stan_seed <- 20261006

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 2. Data preparation
# -----------------------------------------------------------------------------

raw <- foreign::read.dta(input_data_path)

raw$country <- as.character(raw$country)
raw$year <- as.integer(raw$year)
raw$gdp <- as.numeric(raw$gdp)

dat <- raw[
  raw$year %in% all_years,
  c("country", "year", "gdp")
]

dat <- dat[order(dat$country, dat$year), ]
row.names(dat) <- NULL

stopifnot(
  !anyNA(dat),
  !anyDuplicated(paste(dat$country, dat$year, sep = "_")),
  "West Germany" %in% dat$country,
  all(table(dat$country, dat$year) == 1L)
)

countries <- sort(unique(dat$country))
donors <- setdiff(countries, "West Germany")

# Treatment indicator: West Germany from 1990 onward.
dat$D <- as.integer(
  dat$country == "West Germany" &
    dat$year >= 1990
)

stopifnot(sum(dat$D) == length(post_years))

bmv_data <- dat[, c("year", "country", "gdp", "D")]

# -----------------------------------------------------------------------------
# 3. Fit bsynth
# -----------------------------------------------------------------------------

# This is the specification used for Appendix D.3:
#   - outcome only
#   - predictor matching OFF
#   - Gaussian process OFF
bmv <- bsynth::bayesianSynth$new(
  data = bmv_data,
  time = year,
  id = country,
  treated = D,
  outcome = gdp,
  ci_width = 0.95,
  gp = FALSE,
  predictor_match = FALSE
)

set.seed(stan_seed)

bmv$fit(
  chains = 4,
  cores = 4,
  iter = 12000,
  warmup = 6000,
  seed = stan_seed,
  control = list(
    adapt_delta = 0.95
  )
)

# -----------------------------------------------------------------------------
# 4. Counterfactual path, Pre-RMSE, and ATT
# -----------------------------------------------------------------------------

# plotData is the package's public summary object. y_synth is the posterior
# counterfactual mean and LB/UB are its uncertainty bounds.
bmv_path <- as.data.frame(bmv$plotData)

stopifnot(
  all(
    c("year", "gdp", "y_synth", "LB", "UB") %in%
      names(bmv_path)
  )
)

bmv_path <- bmv_path[order(bmv_path$year), ]

stopifnot(
  identical(
    as.integer(bmv_path$year),
    as.integer(all_years)
  )
)

bmv_path$effect <- bmv_path$gdp - bmv_path$y_synth

bmv_pre_rmse <- sqrt(
  mean(
    bmv_path$effect[
      bmv_path$year %in% pre_years
    ]^2
  )
)

bmv_att <- mean(
  bmv_path$effect[
    bmv_path$year %in% post_years
  ]
)

summary_table <- data.frame(
  method = "bsynth (outcome only; no GP)",
  pre_RMSE = bmv_pre_rmse,
  post_ATT = bmv_att
)

cat(
  "\n=== bsynth (GP OFF; predictor matching OFF) ===\n"
)
print(summary_table, row.names = FALSE)

# -----------------------------------------------------------------------------
# 5. Posterior donor weights
# -----------------------------------------------------------------------------

# In the package version used for the analysis, the underlying Stan
# fit is stored in the private fitted object. We use it only to summarize the
# posterior draws of beta; the fitted trajectory above is taken from plotData.
bmv_stan <- bmv$.__enclos_env__$private$fitted

beta_draws <- as.matrix(
  bmv_stan,
  pars = "beta"
)

stopifnot(
  ncol(beta_draws) == length(donors)
)

bmv_weight_table <- data.frame(
  country = donors,
  weight = colMeans(beta_draws),
  posterior_sd = apply(beta_draws, 2, sd)
)

cat("\nPosterior donor-weight summaries:\n")
print(
  bmv_weight_table[
    order(-bmv_weight_table$weight),
  ],
  row.names = FALSE
)

stan_beta_diagnostics <- rstan::summary(
  bmv_stan,
  pars = "beta"
)$summary[, c("mean", "n_eff", "Rhat")]

print(stan_beta_diagnostics)

# -----------------------------------------------------------------------------
# 6. Export
# -----------------------------------------------------------------------------

write.csv(
  bmv_weight_table,
  file.path(output_dir, "bsynth_weights.csv"),
  row.names = FALSE
)

write.csv(
  bmv_path,
  file.path(output_dir, "bsynth_trajectory.csv"),
  row.names = FALSE
)

write.csv(
  summary_table,
  file.path(output_dir, "bsynth_summary.csv"),
  row.names = FALSE
)

write.csv(
  stan_beta_diagnostics,
  file.path(output_dir, "bsynth_beta_diagnostics.csv"),
  row.names = TRUE
)

writeLines(
  c(
    paste(
      "bsynth:",
      as.character(utils::packageVersion("bsynth"))
    ),
    paste(
      "rstan:",
      as.character(utils::packageVersion("rstan"))
    ),
    capture.output(sessionInfo())
  ),
  file.path(output_dir, "sessionInfo.txt")
)
