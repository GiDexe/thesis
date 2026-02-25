# =============================================================================
# 02_grid_search_standardised.R

# =============================================================================

library(dplyr)
library(Matrix)
library(glmnet)
library(CVXR)
library(MASS)
library(optimx)
library(fixest)

source("/workspaces/thesis/functions.r")
set.seed(1861)
# === LOAD DATA ===============================================================

# Skipped: the script is supposed to be run after grid_search_32.r

# === LAMBDA GRID =============================================================

grid <- expand.grid(
  lambdaL1     = c(0.05, 0.1),
  lambdaL2     = c(0.01),
  lambdaL1Star = c(0.25, 0.375, 0.75)
)

# === STANDARDISE =============================================================
# Standardise y and x to unit variance for both ISO standards
# This ensures penalisation acts comparably across variables and across ISOs
# rho and W are invariant to scaling; beta and gamma become "per SD" effects

standardise <- function(data) {
  data %>%
    dplyr::mutate(
      y  = (y - mean(y)) / sd(y),
      x1 = (x1 - mean(x1)) / sd(x1),
      x2 = (x2 - mean(x2)) / sd(x2),
      x3 = (x3 - mean(x3)) / sd(x3),
      x4 = (x4 - mean(x4)) / sd(x4)
    )
}

data_37001_std <- standardise(data_37001)
data_45001_std <- standardise(data_45001)

cat("=== STANDARDISED DATA ===\n")
cat(
  "ISO 37001 - N:", length(unique(data_37001_std$id)),
  "T:", length(unique(data_37001_std$time)), "\n"
)
cat(
  "  y: mean=", round(mean(data_37001_std$y), 3),
  "sd=", round(sd(data_37001_std$y), 3), "\n"
)
cat(
  "ISO 45001 - N:", length(unique(data_45001_std$id)),
  "T:", length(unique(data_45001_std$time)), "\n"
)
cat(
  "  y: mean=", round(mean(data_45001_std$y), 3),
  "sd=", round(sd(data_45001_std$y), 3), "\n"
)

cat("\nCombinazioni nella griglia:", nrow(grid), "\n")

# === FUNCTION: RUN GRID SEARCH ===============================================

run_grid_search <- function(data, iso_label, grid) {
  n_grid <- nrow(grid)
  results <- vector("list", n_grid)
  bic_vec <- rep(NA, n_grid)
  timings <- rep(NA, n_grid)

  # ---  Wder ---
  cat("\n", iso_label, "- Calcolo Wder (run iniziale)...\n")
  t0 <- Sys.time()

  tryCatch(
    {
      first_run <- recoverNetwork(
        data = data,
        lambda = as.numeric(grid[1, ]),
        exoeffects = 1
      )

      Wder <- first_run$Wder
      results[[1]] <- first_run
      bic_vec[1] <- first_run$BIC
      timings[1] <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

      cat(
        iso_label, "- Run 1/", n_grid,
        " lambda=", as.numeric(grid[1, ]),
        " BIC=", round(bic_vec[1], 4),
        " tempo=", round(timings[1], 1), "min\n"
      )
    },
    error = function(e) {
      timings[1] <<- as.numeric(difftime(Sys.time(), t0, units = "mins"))
      cat(" ERRORE in the first run:", e$message, "\n")
      cat(" Missing Wder. Aborting ", iso_label, "\n")

      save(results, bic_vec, timings, grid,
        file = paste0("data/output_data/grid_search_", iso_label, "_FAILED.RData")
      )
      return(list(results = results, bic_vec = bic_vec, timings = timings, Wder = NULL))
    }
  )

  # Check if first run failed
  if (is.na(bic_vec[1])) {
    return(list(results = results, bic_vec = bic_vec, timings = timings, Wder = NULL))
  }

  # --- Already having Wder ---
  for (g in 2:n_grid) {
    cat(
      iso_label, "- Run", g, "/", n_grid,
      " lambda=", as.numeric(grid[g, ])
    )

    t0 <- Sys.time()

    tryCatch(
      {
        results[[g]] <- recoverNetwork(
          data = data,
          lambda = as.numeric(grid[g, ]),
          exoeffects = 1,
          Wder = Wder
        )
        bic_vec[g] <- results[[g]]$BIC
        timings[g] <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

        cat(
          " BIC=", round(bic_vec[g], 4),
          " time=", round(timings[g], 1), "min\n"
        )
      },
      error = function(e) {
        timings[g] <<- as.numeric(difftime(Sys.time(), t0, units = "mins"))
        cat(" ERROR:", e$message, "\n")
      }
    )

    # Save (safety)
    save(results, bic_vec, timings, grid, Wder,
      file = paste0("data/output_data/grid_search_", iso_label, "_PARTIAL.RData")
    )
  }

  # final save
  save(results, bic_vec, timings, grid, Wder,
    file = paste0("data/output_data/grid_search_", iso_label, ".RData")
  )

  return(list(results = results, bic_vec = bic_vec, timings = timings, Wder = Wder))
}

# === HELPER: PRINT SUMMARY ===================================================

print_summary <- function(gs, iso_label, grid) {
  valid <- which(!is.na(gs$bic_vec))
  if (length(valid) == 0) {
    cat("\n=== ", iso_label, " === All the runs failed\n")
    return(invisible(NULL))
  }
  best <- valid[which.min(gs$bic_vec[valid])]
  cat("\n===", iso_label, "===\n")
  cat("  Best lambda:", as.numeric(grid[best, ]), "\n")
  cat("  BIC:", gs$bic_vec[best], "\n")
  cat("  rho:", gs$results[[best]]$unpenalisedgmm$rho, "\n")
  cat("  N edges:", sum(abs(gs$results[[best]]$adaen$W) > 1e-5), "\n")
  cat("  Successful Runs :", length(valid), "/", nrow(grid), "\n")
}

# =============================================================================
# PART 1: STANDARDISED, N=37
# =============================================================================

cat("\n########################################\n")
cat("PART 1: STANDARDISED, N=37\n")
cat("########################################\n")
cat("Hour:", format(Sys.time(), "%H:%M"), "\n")

# --- ISO 45001 standardised ---
cat("\n>>> ISO 45001 STANDARDISED <<<\n")
t_start <- Sys.time()
gs_45001_std <- run_grid_search(data_45001_std, "ISO45001_standardised", grid)
cat(
  "\nISO 45001 std completed in",
  round(as.numeric(difftime(Sys.time(), t_start, units = "hours")), 1), "ore\n"
)
# --- ISO 37001 standardised ---
cat("\n>>> ISO 37001 STANDARDISED <<<\n")
t_start <- Sys.time()
gs_37001_std <- run_grid_search(data_37001_std, "ISO37001_standardised", grid)
cat(
  "\nISO 37001 std completed in",
  round(as.numeric(difftime(Sys.time(), t_start, units = "hours")), 1), "ore\n"
)


print_summary(gs_37001_std, "ISO 37001 STANDARDISED", grid)
print_summary(gs_45001_std, "ISO 45001 STANDARDISED", grid)

# --- Stability standardised ---
cat("\n=== STABILITY (STANDARDISED) ===\n")
for (label in c("ISO37001_standardised", "ISO45001_standardised")) {
  load(paste0("data/output_data/grid_search_", label, ".RData"))
  valid <- which(!is.na(bic_vec))
  if (length(valid) > 1) {
    all_W <- lapply(results[valid], function(r) (abs(r$adaen$W) > 1e-5) + 0)
    stability <- Reduce("+", all_W) / length(all_W)
    cat(paste0("\n", label, ":\n"))
    cat("  Edge robusti (>80%):", sum(stability > 0.8), "\n")
    cat("  Edge robusti (>60%):", sum(stability > 0.6), "\n")
    cat("  Edge robusti (>40%):", sum(stability > 0.4), "\n")
  }
}

# =============================================================================
# FINAL SUMMARY
# =============================================================================

cat("\n########################################\n")
cat("Status: Completed\n")
cat("Hour:", format(Sys.time(), "%H:%M"), "\n")
cat("########################################\n")

cat("\nFiles saved in:\n")
cat("  grid_search_ISO37001_standardised.RData\n")
cat("  grid_search_ISO45001_standardised.RData\n")
