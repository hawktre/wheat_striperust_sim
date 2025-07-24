## ---------------------------
##
## Script name: 04b_RunSim.R
##
## Purpose of script: Run the simulation
##
## Author: Trent VanHawkins
##
## Date Created: 2025-07-23
##
##
## ---------------------------

## load up the packages we will need:  (uncomment as required)
.libPaths("/home/hawkintr/R_libs/4.4.1")
options(repos = c(CRAN = "https://cloud.r-project.org"))

library(here)
library(data.table)
library(dplyr)
library(purrr)
library(stringr)
library(lubridate)
library(tidyr)
library(furrr)
library(MASS)
source(here("Code/02a_GradDescentFun.R"))
source(here("Code/03a_EMgradfun.R"))
source(here("Code/04a_SimFunc.R"))

# Read in the data --------------------------------------------------------
forward_fits <- readRDS(here("DataProcessed/results/forward_model/forward_fits.rds"))
mod_dat <- readRDS(here("DataProcessed/experimental/mod_dat_arrays.rds"))

## Set up Simulation
plan(multicore(workers = availableCores()))
nsim <- 10000

## Run the simulation
set.seed(404)
sims <- 1:nsim %>% future_map_dfr(~single_sim(.x, mod_dat, forward_fits, output_dir = here("DataProcessed/results/simulation/errors/")))


saveRDS(sims, here("DataProcessed/results/simulation/sims.rds"))