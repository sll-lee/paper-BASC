# Bayesian Adaptive Synthetic Control (BASC)

This repository provides replication code for the West Germany application of the Bayesian Adaptive Synthetic Control (BASC) method.

The repository contains code for BASC, standard synthetic control, and the comparison methods used in the West Germany application.

## Repository structure

```text
.
├── 01_BASC_full.R              # BASC with full-period fitting and donor-set selection
├── 02_BASC_pre.R               # BASC fitted using the pre-intervention period
├── 03_BASC_pre_NS.R            # Pre-intervention BASC without donor-set selection
├── 04_standard_SCM.R           # Standard synthetic control
├── 05_bsynth.R                 # Bayesian synthetic control using the bsynth package
├── 06_fPCA_SYNTH.ipynb         # fPCA-SYNTH / RPCA-SYNTH comparison
├── 07_ClusterSC.ipynb          # ClusterSC comparison
├── requirements_python.txt     # Python dependencies
└── data/
    └── repgermany.dta          # West Germany GDP data
```

## Data

The empirical application uses the West Germany dataset from Abadie et al. (2015), containing annual per-capita GDP and related variables for West Germany and the donor countries.

The data file used by all scripts is

```text
data/repgermany.dta
```

The analysis period is 1960--2003, with German reunification in 1990 treated as the intervention. The pre-intervention period is 1960--1989 and the post-intervention period is 1990--2003.

## BASC specifications

### `01_BASC_full.R`

Fits the main BASC specification over the full 1960--2003 period. The model performs donor-set selection and uses a two-dimensional linear post-intervention effect basis.

### `02_BASC_pre.R`

Fits BASC using only the pre-intervention observations from 1960--1989. The post-intervention latent Gaussian-process component is generated conditionally on the pre-intervention GP draws.

### `03_BASC_pre_NS.R`

Fits the no-selection counterpart of BASC-pre. All donor-inclusion indicators are fixed at one, and model parameters are estimated using 1960--1989 only.

## Comparison methods

### `04_standard_SCM.R`

Fits a standard outcome-only synthetic control using the original, unscaled GDP outcomes over the pre-intervention period.

### `05_bsynth.R`

Fits Bayesian synthetic control using the official `bsynth` package under an outcome-only specification, with predictor matching and the Gaussian-process component turned off.

### `06_fPCA_SYNTH.ipynb`

Implements the fPCA-SYNTH / RPCA-SYNTH comparison using the `mlsynth` package:

https://github.com/jgreathouse9/mlsynth

The notebook includes the executed output.

### `07_ClusterSC.ipynb`

Implements the ClusterSC comparison using the authors' `syclib` implementation:

https://github.com/srho1/ClusterSC

Donor clustering uses the donor trajectories over the full observation period, while the treated unit is assigned to a donor cluster using its pre-intervention trajectory. After donor selection, an OLS synthetic control is fitted directly to the selected donors using the original 1960--1989 GDP outcomes. The number of clusters is selected by minimizing the pre-intervention RMSE over the candidate grid.

The notebook is saved with executed output.

## R requirements

The R scripts use packages including

```r
foreign
dplyr
tidyr
tmvtnorm
MASS
MCMCpack
fields
coda
Synth
bsynth
rstan
```

For example,

```r
install.packages(c(
  "foreign", "dplyr", "tidyr", "tmvtnorm",
  "MASS", "MCMCpack", "fields", "coda", "Synth"
))
```

The `bsynth` and `rstan` dependencies should be installed according to the installation instructions for the `bsynth` package.

## Python requirements

Python dependencies for fPCA-SYNTH and ClusterSC are listed in

```text
requirements_python.txt
```

Install them from the repository root using

```bash
pip install -r requirements_python.txt
```

The Python notebooks assume that the data file is available at

```text
data/repgermany.dta
```

## Running the empirical application

Run the R scripts from the repository root:

```r
source("01_BASC_full.R")
source("02_BASC_pre.R")
source("03_BASC_pre_NS.R")
source("04_standard_SCM.R")
source("05_bsynth.R")
```

Run the following notebooks from the repository root:

```text
06_fPCA_SYNTH.ipynb
07_ClusterSC.ipynb
```

from the repository root.

The BASC MCMC scripts use long chains and may require substantial computation time.

## Outputs

The scripts save method-specific summaries, donor weights, fitted or counterfactual trajectories, and diagnostic output under

```text
results/realdata/
```

The Python notebooks also display the main numerical results directly in the saved notebook output.

## Citation

If you use this code, please cite:

> Bayesian Donor Set Selection in Synthetic Controls.

## Contact

If you encounter any issues or require additional information to reproduce the results, please contact the authors.
