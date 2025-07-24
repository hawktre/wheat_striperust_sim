# setup_packages.R
my_lib <- "/home/hawkintr/R_libs/4.4.1"
dir.create(my_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(my_lib)

options(repos = c(CRAN = "https://cloud.r-project.org"))

pkgs <- c("here", "data.table", "dplyr", "purrr", "stringr",
          "lubridate", "tidyr", "furrr", "MASS")

install.packages(pkgs, lib = my_lib, dependencies = TRUE)
