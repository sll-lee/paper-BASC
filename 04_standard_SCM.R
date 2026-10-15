###############################################################################
# Standard synthetic control (s-SCM)
# West Germany application
#
# Outcome-only SCM using annual pre-intervention GDP outcomes (1960--1989).
# The weights solve the simplex-constrained least-squares problem on the
# original GDP scale.
###############################################################################

rm(list = ls())

# -----------------------------------------------------------------------------
# 0. Packages
# -----------------------------------------------------------------------------

required_packages <- c(
  "foreign",
  "Synth"
)

invisible(lapply(required_packages, require, character.only = TRUE))

# -----------------------------------------------------------------------------
# 1. User settings
# -----------------------------------------------------------------------------

input_data_path <- "data/repgermany.dta"
output_dir <- "results/realdata/standard_SCM"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

pre_years <- 1960:1989
post_years <- 1990:2003
all_years <- 1960:2003

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

dat$country_id <- match(dat$country, countries)

# -----------------------------------------------------------------------------
# 3. Standard SCM
# -----------------------------------------------------------------------------

# Each pre-intervention year enters as a separate outcome predictor.
scm_data <- Synth::dataprep(
  foo = dat,
  predictors = NULL,
  special.predictors = lapply(
    pre_years,
    function(yr) list("gdp", yr, "mean")
  ),
  dependent = "gdp",
  unit.variable = "country_id",
  unit.names.variable = "country",
  time.variable = "year",
  treatment.identifier = match("West Germany", countries),
  controls.identifier = match(donors, countries),
  time.predictors.prior = pre_years,
  time.optimize.ssr = pre_years,
  time.plot = all_years
)

# Synth internally standardizes predictor rows by their across-unit SD.
# Choosing custom.v proportional to the row variances cancels that scaling,
# so the resulting donor weights minimize the unweighted pre-intervention
# squared error on the original GDP scale.
stopifnot(nrow(scm_data$X0) == length(pre_years))

v_raw <- apply(
  cbind(scm_data$X0, scm_data$X1),
  1,
  stats::var
)

stopifnot(
  all(is.finite(v_raw)),
  all(v_raw > 0)
)

scm_fit <- Synth::synth(
  data.prep.obj = scm_data,
  custom.v = v_raw
)

# -----------------------------------------------------------------------------
# 4. Donor weights and counterfactual path
# -----------------------------------------------------------------------------

scm_w <- as.numeric(scm_fit$solution.w)
scm_ids <- as.integer(rownames(scm_fit$solution.w))

scm_weight_table <- data.frame(
  country = countries[scm_ids],
  weight = scm_w,
  row.names = NULL
)

scm_path <- data.frame(
  year = all_years,
  observed = drop(scm_data$Y1plot),
  counterfactual = drop(scm_data$Y0plot %*% scm_w)
)

scm_path$effect <- scm_path$observed - scm_path$counterfactual

scm_pre_rmse <- sqrt(
  mean(
    scm_path$effect[
      scm_path$year %in% pre_years
    ]^2
  )
)

scm_att <- mean(
  scm_path$effect[
    scm_path$year %in% post_years
  ]
)

stopifnot(
  abs(sum(scm_w) - 1) < 1e-3,
  all(scm_w >= -1e-5)
)

summary_table <- data.frame(
  method = "standard SCM (outcome only)",
  pre_RMSE = scm_pre_rmse,
  post_ATT = scm_att
)

cat("\n=== Standard SCM (outcome-only; original GDP scale) ===\n")
print(
  scm_weight_table[order(-scm_weight_table$weight), ],
  row.names = FALSE
)
print(summary_table, row.names = FALSE)

# -----------------------------------------------------------------------------
# 5. Independent numerical check
# -----------------------------------------------------------------------------

# Optional: solve the same simplex-constrained least-squares problem using
# quadprog. This is used only as a numerical cross-check and is not required
# for the reported estimator.
if (requireNamespace("quadprog", quietly = TRUE)) {
  x_pre <- scm_data$Z0
  y_pre <- drop(scm_data$Z1)

  scale_gdp <- max(abs(c(x_pre, y_pre)))
  xs <- x_pre / scale_gdp
  ys <- y_pre / scale_gdp

  J <- ncol(xs)

  qp_result <- tryCatch(
    quadprog::solve.QP(
      Dmat = crossprod(xs) + diag(1e-10, J),
      dvec = drop(crossprod(xs, ys)),
      Amat = cbind(rep(1, J), diag(J)),
      bvec = c(1, rep(0, J)),
      meq = 1
    ),
    error = function(e) e
  )

  if (inherits(qp_result, "error")) {
    warning(
      "Independent QP check unavailable: ",
      conditionMessage(qp_result)
    )
  } else {
    qp_rmse <- sqrt(
      mean(
        (y_pre - drop(x_pre %*% qp_result$solution))^2
      )
    )

    cat("Independent QP pre-RMSE:", qp_rmse, "\n")

    if (abs(qp_rmse - scm_pre_rmse) > 1) {
      warning(
        "Synth and the independent QP solution differ by ",
        "more than 1 GDP unit in pre-RMSE."
      )
    }
  }
}

# -----------------------------------------------------------------------------
# 6. Export
# -----------------------------------------------------------------------------

write.csv(
  scm_weight_table,
  file.path(output_dir, "standard_SCM_weights.csv"),
  row.names = FALSE
)

write.csv(
  scm_path,
  file.path(output_dir, "standard_SCM_trajectory.csv"),
  row.names = FALSE
)

write.csv(
  summary_table,
  file.path(output_dir, "standard_SCM_summary.csv"),
  row.names = FALSE
)

writeLines(
  c(
    paste(
      "Synth:",
      as.character(utils::packageVersion("Synth"))
    ),
    capture.output(sessionInfo())
  ),
  file.path(output_dir, "sessionInfo.txt")
)
