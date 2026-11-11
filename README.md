# Credible Region Bayesian Maximum Entropy (CRBME)

Research code for spatial prediction and uncertainty quantification using **Credible Region Bayesian Maximum Entropy (CRBME)** in R.

## Overview

This repository contains R scripts used to investigate spatial prediction when observations are represented with uncertainty. The CRBME workflow constructs uncertainty information and incorporates it into Bayesian Maximum Entropy (BME) spatial prediction, with simulation-based comparisons and prediction maps.

**This repository is not an R package.** There is no `install_github()` step. Download or clone the repository, install its R dependencies, and load the project functions with `source()`.

The code is intended for research and reproducibility. The exact interface and required inputs of individual functions should be checked in the corresponding source files.

## Repository structure

The project is organized approximately as follows (the contents of individual directories may evolve):

```text
CRBME/
├── R/                         # Custom CRBME functions
│   ├── new_vg_results.R
│   ├── new_prob_zk.R
│   ├── new_bme_estimate.R
│   ├── new_bme_predict.R
│   └── new_bme_cv.R
├── Gumbel_Simulations/        # Simulation and plotting scripts for Gumbel data
├── Lognormal_Simulations/     # Simulation and plotting scripts for log-normal data
├── Real Data/                 # Real-data analyses (if publicly distributable)
├── figs/                      # Generated figures
├── simulation_results.csv     # Simulation results (if included)
└── CRBME.Rproj                # RStudio project file
```

Some simulation scripts may depend on additional local functions or input files. Review their `source()` and file-reading statements before running them.

## Requirements

- R (use a recent version; the project has been run with R 4.6.2)
- R packages used by the simulation scripts, including `gstat`, `sf`, `caret`, `dplyr`, `evd`, `mvtnorm`, and `ggplot2`
- `BMEmapping` for scripts that call functions from the BMEmapping package

Install the common CRBME simulation dependencies:

```r
install.packages(c(
  "gstat", "sf", "caret", "dplyr", "evd",
  "mvtnorm", "ggplot2", "viridis"
))
```

For scripts that use BMEmapping, install it separately:

```r
install.packages("remotes")
remotes::install_github("KinsprideDuah/BMEmapping")
```

Installation of `sf` and its system dependencies may require additional setup on Linux.

## Getting started

### 1. Download the repository

Clone the repository after publishing it to GitHub:

```bash
git clone https://github.com/YOUR-USERNAME/CRBME.git
cd CRBME
```

Alternatively, choose **Code > Download ZIP** on GitHub and extract the files.

### 2. Open R in the project directory

In RStudio, open `CRBME.Rproj`. In VS Code or an R terminal, set the working directory to the repository root (replace the path below):

```r
setwd("/path/to/CRBME")
getwd()
```

### 3. Load the CRBME functions

```r
source("R/new_vg_results.R")
source("R/new_prob_zk.R")
source("R/new_bme_estimate.R")
source("R/new_bme_predict.R")
source("R/new_bme_cv.R")
```

These are project scripts, **not installed package functions**. Use the actual function names defined in these files after sourcing them. If a script calls `BMEmapping::` or `BMEmapping:::`, it is explicitly using a function from that separate package; check compatibility with the version installed.

### 4. Run a simulation

For example, from the repository root:

```r
source("Lognormal_Simulations/lognormal_simulations.R")
```

For the Gumbel experiments, open the corresponding script in `Gumbel_Simulations/` and run it after loading its dependencies. Check the script for required paths and output filenames before execution.

## Methodology

The research workflow investigates BME spatial prediction with credible-region-based uncertainty representations. At a high level:

1. Generate or import spatial coordinates and observations.
2. Specify the spatial dependence model, including a variogram.
3. Construct the uncertainty representation used by CRBME.
4. Estimate model quantities and obtain spatial predictions using the project-specific functions.
5. Evaluate predictions at held-out locations and summarize errors.
6. Produce prediction maps and compare simulation configurations.

The detailed mathematical formulation, estimation algorithm, and interpretation of the credible-region parameters should be taken from the accompanying manuscript and implementation, rather than inferred from the brief overview here.

## Simulation experiments

The research includes experiments with **Gumbel** and **log-normal** simulated data. Example settings used in the Gumbel experiments include:

| Component | Example setting |
| --- | --- |
| Spatial domain | `[0, 100] × [0, 100]` |
| Number of spatial sites | 625 |
| Variogram model | Exponential |
| Partial sill | 0.5 |
| Range parameter | 60 |
| Nugget | 0.5 |
| Gumbel location | 0.5 |
| Gumbel scale | 1 |
| Hold-out sites | 100 |
| Realizations | 5 in one simulation design |

These are **example experiment settings**, not hard-coded requirements of CRBME. Consult each simulation script for its actual random seed, parameter grid, number of realizations, and data-generation procedure. For reproducibility, keep the random seed and the train/test split consistent when comparing methods.

## Validation and results

Simulation scripts may evaluate held-out predictions using metrics such as root mean squared error (RMSE) and mean absolute error (MAE). The project also produces spatial prediction visualizations and may save results in `simulation_results.csv`.

When interpreting results, distinguish between the simulated latent field, any kriged or smoothed reference surface, and predictions evaluated at the original held-out sites. These are not interchangeable definitions of ground truth.

## Reproducibility notes

- Run scripts from the repository root unless they specify otherwise.
- Inspect each script for hard-coded local paths and update them for your machine.
- Record the R version, package versions, seeds, and simulation settings used for published figures.
- Do not rerun a random realization when attempting to reproduce a figure from an earlier simulation; reuse the original simulated values and split.
- Some simulations can be computationally intensive and may be better suited to a server.
- Verify that all custom functions referenced by a script are included in the repository.

## Related software

This project uses ideas and, in some scripts, functions from [BMEmapping](https://github.com/KinsprideDuah/BMEmapping), an R package for Bayesian Maximum Entropy spatial interpolation. **CRBME is a separate research-code repository and should not be installed as though it were the BMEmapping package.**

## Citation

If you use this research code, please cite the associated CRBME manuscript when its bibliographic details become available. A formal citation will be added here.

## License

No license is specified in this README. Add a `LICENSE` file if you wish to define terms for reuse and redistribution.

## Contact

For questions, issues, or reproducibility inquiries, please use the GitHub repository's Issues section once it is available.
