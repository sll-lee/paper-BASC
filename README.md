# Bayesian Adaptive Synthetic Control (BASC)

This repository provides code for the real-data application of the Bayesian Adaptive Synthetic Control (BASC) method.

The main script reproduces the West Germany application using the BASC model with donor-set selection and posterior counterfactual prediction.

## Repository structure

```text
.
├── clean_basc_realdata.R        # Main script for the West Germany application
├── data/
│   └── repgermany.dta           # Input data file; not included if restricted by license
└── results/
    └── realdata/                # Output files created by the script
        └── plot/                # Figures created by the script
```

## Requirements

The script is written in R and uses the following packages:

```r
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
```

Install missing packages before running the script, for example:

```r
install.packages(c(
  "tidysynth", "foreign", "dplyr", "tidyr", "tmvtnorm",
  "MASS", "MCMCpack", "fields", "coda", "ggplot2", "reshape2"
))
```

## Data

The real-data example uses the West Germany dataset in Stata format:

```text
data/repgermany.dta
```

If the data file is not included in this repository due to licensing restrictions, place the file manually under the `data/` directory or modify the following line in `clean_basc_realdata.R`:

```r
input_data_path <- "data/repgermany.dta"
```

## How to run

Open R and run:

```r
source("clean_basc_realdata.R")
```

By default, the script runs the MCMC sampler:

```r
run_mcmc <- TRUE
```

After the MCMC output has been saved, set

```r
run_mcmc <- FALSE
```

to load the saved result object and regenerate posterior summaries and plots without rerunning the sampler.

## Main outputs

The script creates the following outputs under `results/realdata/`:

```text
BASC_realdata_results.RData              # Saved MCMC results
BASC_realdata_weight_summary.csv         # Posterior donor-weight and inclusion summaries
BASC_realdata_prediction.RData           # Posterior counterfactual prediction object
BASC_realdata_prediction_95ci.csv        # Posterior median and 95% credible intervals
plot/BASC_realdata_counterfactual.png    # Counterfactual path plot
plot/trace_plots_realdata.*              # Trace plots
plot/density_plots_realdata.*            # Marginal posterior density plots
plot/realdata_gamma.*                    # Posterior donor-inclusion plot
```

## Notes

- The MCMC update formulas follow the original analysis code used for the paper.
- The script has been reorganized for readability and reproducibility, but the main sampling steps are intentionally preserved.
- The current setting uses three chains with seeds `103`, `200`, and `305`.
- The donor weights are computed from the posterior samples of the inclusion indicators and positive latent weights.

## Citation

If you use this code, please cite the corresponding paper.

