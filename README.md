
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RFLC-SCP <img src="man/figures/compendium-sticker.png" align="right" style="float:right; height:120px;"/>

<!-- badges: start -->

[![License: GPL (\>=
2)](https://img.shields.io/badge/License-GPL%20%28%3E%3D%202%29-blue.svg)](https://choosealicense.com/licenses/gpl-2.0/)
<!-- badges: end -->

<p align="left">
• <a href="#overview">Overview</a><br> •
<a href="#features">Features</a><br> •
<a href="#content">Content</a><br> •
<a href="#installation">Installation</a><br> •
<a href="#usage">Usage</a><br> • <a href="#citation">Citation</a><br> •
<a href="#contributing">Contributing</a><br> •
<a href="#acknowledgments">Acknowledgments</a><br> •
<a href="#references">References</a>
</p>

## Overview

RFLC-SCP: A comprehensive framework to assess landscape connectivity for
conservation planning

This project provides all the data and codes used to produce the
methodological framework and results presented in:

Prima, M.-C., Renaud, J., Witté, I., Suarez, L., Rouveyrol, P.,
Fernando, M., Sacchi, A., Cosentino, F., Santini, L., Maiorano, L.,
Moreira, F., Dertien, J., Fernández, N., & Thuiller, W. (in prep.). A
comprehensive framework to assess multi-species landscape connectivity
for conservation planning.

## Features

The workflow is as follow,

- Step 1. Clustering species in functional groups having similar traits
  and environmental niches (analyses/Rcode01_GetFunctionalGroups.R).
- Step 2. Generating resistance and habitat suitability maps for each
  group based on an ensemble of species distribution models
  (analyses/Rcode2_GetResisSuitMaps.R).
- Step 3. Generating ecological continuities per group
  (analyses/Rcode3_GetEcologicalContinuities.R).
- Step 4. Calculating multi-scale network metrics
  (analyses/Rcode4_GetConnectivityMetrics.R).
- Step 5. (optional) Plot figures (analyses/Rcode5_GetFigures.R).

## Content

This repository is structured as follow:

- [`DESCRIPTION`](https://github.com/mcpri3/RFLC-SCP/tree/master/DESCRIPTION):
  contains project metadata (authors, date, dependencies, etc.)

- [`make.R`](https://github.com/mcpri3/RFLC-SCP/tree/master/make.R):
  main R script to load project dependencies

- [`R/`](https://github.com/mcpri3/RFLC-SCP/tree/master/R): contains R
  functions developed especially for this project

- [`man/`](https://github.com/mcpri3/RFLC-SCP/tree/master/man): contains
  help files of R functions

- [`analyses/`](https://github.com/mcpri3/RFLC-SCP/tree/master/analyses):
  contains R scripts to run each step of the workflow

- [`data/raw-data/`](https://filesender.renater.fr/?s=download&token=0b5de737-082c-478e-8f64-845251e60ab9):
  contains all raw data required to perform analyses. Due to its size,
  \[`data/raw-data/`\] folder was stored on another platform, accessible
  here:
  <https://filesender.renater.fr/?s=download&token=0b5de737-082c-478e-8f64-845251e60ab9>
  These data need to be downloaded and stored in the project folder
  before running the workflow

- [`data/derived-data/`](https://filesender.renater.fr/?s=download&token=0b5de737-082c-478e-8f64-845251e60ab9):
  contains all intermediate data created during the workflow. Due to its
  size, \[`data/derived-data/`\] folder was stored on another platform,
  accessible here:
  <https://filesender.renater.fr/?s=download&token=0b5de737-082c-478e-8f64-845251e60ab9>

- [`outputs/`](https://filesender.renater.fr/?s=download&token=0b5de737-082c-478e-8f64-845251e60ab9):
  contains final data created during the workflow. Due to its size,
  \[`outputs/`\] folder was stored on another platform, accessible here:
  <https://filesender.renater.fr/?s=download&token=0b5de737-082c-478e-8f64-845251e60ab9>

- [`figures/`](https://github.com/mcpri3/RFLC-SCP/tree/master/figures):
  contains all the figures created during the workflow

## Installation

To install this compendium:

- [Fork](https://docs.github.com/en/get-started/quickstart/contributing-to-projects)
  this repository using the GitHub interface.
- [Clone](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository)
  your fork using `git clone fork-url` (replace `fork-url` by the URL of
  your fork). Alternatively, open [RStudio
  IDE](https://posit.co/products/open-source/rstudio/) and create a New
  Project from Version Control.

## Usage

Launch the
[`make.R`](https://github.com/mcpri3/RFLC-SCP/tree/master/make.R) file
with:

``` r
source("make.R")
```

Then go to the \[`analyses/`\] folder and run successively the different
scripts.

**Notes**

- All required packages listed in the `DESCRIPTION` file will be
  installed (if necessary)
- All required packages and R functions will be loaded
- Some analyses listed in the script of the \[`analyses/`\] folder might
  take time

## Citation

Please use the following citation:

> Prima, M.-C., Renaud, J., Witté, I., Suarez, L., Rouveyrol, P.,
> Fernando, M., Sacchi, A., Cosentino, F., Santini, L., Maiorano, L.,
> Moreira, F., Dertien, J., Fernandez, N., & Thuiller, W. (in prep.). A
> comprehensive framework to assess multi-species landscape connectivity
> for conservation planning.

## Contributing

All types of contributions are encouraged and valued. For more
information, check out our [Contributor
Guidelines](https://github.com/mcpri3/RFLC-SCP/blob/main/CONTRIBUTING.md).

Please note that this project is released with a [Contributor Code of
Conduct](https://contributor-covenant.org/version/2/1/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
