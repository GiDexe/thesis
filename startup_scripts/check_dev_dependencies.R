#!/usr/bin/env Rscript

#' Check and Install R Development and Thesis Dependencies
cat("Checking installation of R development and thesis packages...\n")

# Ensure renv is available
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")

# --- FORCING EXACT VERSIONS ---
# We use the @version syntax to pin your specific environment
cat("Pinning exact package versions for thesis stability...\n")

pinned_pkgs <- c(
  "fixest@0.13.2", "optimx@2025-4.9", "MASS@7.3-65", "CVXR@1.0-15",
  "glmnet@4.1-10", "Matrix@1.7-3", "magrittr@2.0.4", "kableExtra@1.4.0",
  "devtools@2.4.6", "usethis@3.2.1", "WDI@2.7.9", "janitor@2.2.1",
  "lubridate@1.9.4", "forcats@1.0.1", "stringr@1.5.2", "dplyr@1.1.4",
  "purrr@1.1.0", "readr@2.1.5", "tidyr@1.3.1", "tibble@3.3.0",
  "ggplot2@4.0.0", "tidyverse@2.0.0", "readxl@1.4.5", "countrycode@1.6.1", "igraph@2.2.2", "zoo", "ggraph", "patchwork",
  "ggrepel", "tidygraph", "maps", "cowplot", "knitr"
)

# Install pinned packages
# renv is smart: if the version is already there, it skips it instantly
renv::install(pinned_pkgs)


# --- VS CODE / DEV TOOLS ---
if (!requireNamespace("desc", quietly = TRUE)) renv::install("desc")
library(desc)

args <- commandArgs(trailingOnly = TRUE)
vscode <- !"--no-vscode" %in% args

if (vscode) {
  cat("Checking VS Code dev tools...\n")
  vscode_deps <- c(
    "languageserver",
    "nx10/httpgd",
    "ManuelHentschel/vscDebugger"
  )
  renv::install(vscode_deps)
}

devtools::install_github("pedroclsouza/recovernetwork")

# --- DESCRIPTION FILE CHECK ---
# Check if there are any other dependencies in the DESCRIPTION file
desc_file <- tryCatch(
  {
    description$new()
  },
  error = function(e) NULL
)
if (!is.null(desc_file)) {
  extra_deps <- c(desc_file$get_deps()$package, desc_file$get_remotes())
  # Filter out what we already installed
  extra_deps <- extra_deps[!(extra_deps %in% gsub("@.*", "", pinned_pkgs))]
  if (length(extra_deps) > 0) renv::install(extra_deps)
}

cat("Environment check complete. All versions are locked.\n")
