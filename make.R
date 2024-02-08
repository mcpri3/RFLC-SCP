#' RFLC-SCP: A Robust Framework to assess Landscape Connectivity for Systematic Conservation Planning application
#' 
#' @description 
#' This project provides all the data and code used to produce the methodological framework and results presented in:
#' 
#' Prima, M.-C., Renaud, J., Witté, I., Suarez, L., Rouveyrol, P., Fernando, M., Sacchi, A., Cosentino, F., Santini, L., Maiorano, L., Moreira, F., Dertien, J.,
#' & Thuiller, W. (in prep.). A robust framework to assess landscape connectivity for systematic conservation planning application. 
#' 
#'
#' The workflow is as follow, 
#' 
#' Step 1. Clustering species in functional groups having similar traits and environmental niches (analyses/Rcode01_GetFunctionalGroups.R).
#' Step 2. Generating resistance and habitat suitability maps for each group based on an ensemble of species distribution models  (analyses/Rcode2_GetResisSuitMaps.R).
#' Step 3. Generating ecological continuities per group (analyses/Rcode3_GetEcologicalContinuities.R).
#' Step 4. Calculating multi-scale network metrics (analyses/Rcode4_GetConnectivityMetrics.R).
#' 
#' First run this make.R file before any run of R scripts.
#' R functions are located in the R/ folder and are automatically loaded when the make.R file is run.
#' Function help can be found by running the classic linecode help('function_name'). 
#' 
#' Initial and intermediate datasets are located is the data/ folder. 
#' Generated outputs are located in the outputs/ folder. 
#'  
#'  
#' @author Marie-Caroline Prima \email{marie-caroline.prima@univ-grenoble-alpes.fr}
#' 
#' @date 2024/01/30



## Install Dependencies (listed in DESCRIPTION) ----

devtools::install_deps(upgrade = "never")


## Load Project Addins (R Functions and Packages) ----

devtools::load_all(here::here())

# Required packages 
library(dplyr)
library(foreach)
library(igraph)

## Run Project ----

# List all R scripts in a sequential order and using the following form:
# source(here::here("analyses", "script_X.R"))

