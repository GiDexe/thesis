# =============================================================================
# ISO Diffusion in the EU: Lagged Network Recovery (1-year lag specification)
# Author: Gianluigi De Rubertis
# Method: De Paula, Rasul & Souza (2025, REStud) -- recoverNetwork
#
# LAGGED SPECIFICATION:
#   Structural model: y_{it} = rho * W * y_{it} + beta * x_{i,t-1}
#                              + gamma * W * x_{i,t-1} + alpha_i + alpha_t + eps
#   All covariates (x1..x5) are lagged one period relative to y_t.
#   y_t = Delta log(1 + n_cert_t) as in the simultaneous spec.
#   This costs one additional period (T = T_sim - 1).
#
# GRID: reduced 3-point grid per dimension (27 combinations vs 64).
# OUTPUT: output_lagged/ (separate from output/)
# =============================================================================

# ---- 0. Setup ---------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(magrittr)
  library(readxl)
  library(stringr)
  library(tidyr)
  library(janitor)
  library(purrr)
  library(readr)
  library(igraph)
  library(remotes)
  library(recoverNetwork)
  library(zoo)
})

set.seed(1861)
setwd("/workspaces/thesis")
options(readr.show_col_types = FALSE)

out_dir <- "output_lagged"
clean_dir <- "/workspaces/thesis/data/cleaned_data"
data_dir <- "/workspaces/thesis/data"

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
if (!dir.exists(clean_dir)) dir.create(clean_dir, recursive = TRUE)

# ---- 1. Country groups ------------------------------------------------------
eu <- c(
  "Austria", "Belgium", "Bulgaria", "Croatia", "Cyprus",
  "Czech Republic", "Denmark", "Estonia", "Finland", "France",
  "Germany", "Greece", "Hungary", "Ireland", "Italy", "Latvia",
  "Lithuania", "Luxembourg", "Malta", "Netherlands", "Poland",
  "Portugal", "Romania", "Slovak Republic", "Slovenia", "Spain",
  "Sweden"
)

# ---- 2. Helper functions ----------------------------------------------------
harmonise_names <- function(df) {
  if (nrow(df) == 0) {
    return(df)
  }
  df %>% dplyr::mutate(country = dplyr::case_match(country,
    "Korea (Republic of)" ~ "Korea", "Korea, Republic of" ~ "Korea",
    "Korea (the Republic of)" ~ "Korea",
    "Czechia" ~ "Czech Republic",
    "Slovakia" ~ "Slovak Republic",
    "United Kingdom of Great Britain and Northern Ireland" ~ "United Kingdom",
    "United States of America" ~ "United States",
    "Netherlands (Kingdom of the)" ~ "Netherlands",
    .default = country
  ))
}

read_iso_hist <- function(file_name, sheets_idx, my_names) {
  f_path <- file.path(data_dir, file_name)
  if (!file.exists(f_path)) stop(paste("File not found:", f_path))
  sheets <- purrr::map(sheets_idx, ~ readxl::read_excel(f_path, sheet = .x))
  expected_cols <- length(my_names)
  purrr::map(1:5, ~ sheets[[.x]] %>%
    dplyr::select(dplyr::all_of(1:expected_cols)) %>%
    dplyr::slice(-(1:2)) %>%
    purrr::set_names(my_names) %>%
    tidyr::pivot_longer(-country, names_to = "year", values_to = "n_cert")) %>%
    purrr::reduce(dplyr::bind_rows) %>%
    harmonise_names() %>%
    dplyr::mutate(n_cert = as.numeric(n_cert), year = as.integer(year)) %>%
    dplyr::group_by(country, year) %>%
    dplyr::summarise(n_cert = sum(n_cert, na.rm = TRUE), .groups = "drop")
}

extract_recent <- function(iso_code) {
  files <- list.files(data_dir,
    pattern = "ISO Survey 20(1[8-9]|2[0-9]).*\\.xlsx$",
    full.names = TRUE
  )
  purrr::map_dfr(files, function(f) {
    yr <- as.numeric(stringr::str_extract(basename(f), "20\\d{2}"))
    sheets <- readxl::excel_sheets(f)
    target <- sheets[grepl(iso_code, sheets)][1]
    if (!is.na(target)) {
      df <- readxl::read_excel(f, sheet = target)
      start_row <- which(
        apply(df, 1, function(x) any(grepl("Afghanistan|Albania", x)))
      )[1]
      df_clean <- df %>%
        dplyr::slice(start_row:dplyr::n()) %>%
        dplyr::rename(country = 1) %>%
        dplyr::mutate(dplyr::across(-1, ~ suppressWarnings(as.numeric(.)))) %>%
        dplyr::select(country, tidyselect::where(is.numeric)) %>%
        dplyr::mutate(n_cert = rowSums(dplyr::select(., -country), na.rm = TRUE)) %>%
        dplyr::select(country, n_cert) %>%
        dplyr::mutate(year = yr)
      return(df_clean)
    }
    return(NULL)
  }) %>%
    harmonise_names() %>%
    dplyr::group_by(country, year) %>%
    dplyr::summarise(n_cert = sum(n_cert, na.rm = TRUE), .groups = "drop")
}

make_covs <- function(start_yr) {
  load_w <- function(f, v) {
    readr::read_csv(file.path(data_dir, f), skip = 4, show_col_types = FALSE) %>%
      dplyr::select(country = 1, as.character(start_yr):"2024") %>%
      tidyr::pivot_longer(-country, names_to = "year", values_to = v) %>%
      dplyr::mutate(
        year  = as.integer(year),
        !!v := as.numeric(as.character(.data[[v]]))
      ) %>%
      harmonise_names()
  }
  load_w("API_NY.GDP.PCAP.KD_DS2_en_csv_v2_234.csv", "gdp_pc") %>%
    dplyr::full_join(load_w("API_NE.TRD.GNFS.ZS_DS2_en_csv_v2_334.csv", "trade"),
      by = c("country", "year")
    ) %>%
    dplyr::full_join(load_w("API_RL.EST_DS2_en_csv_v2_1905.csv", "rol"),
      by = c("country", "year")
    ) %>%
    dplyr::full_join(load_w("API_EG.FEC.RNEW.ZS_DS2_en_csv_v2_4948.csv", "renew"),
      by = c("country", "year")
    ) %>%
    dplyr::full_join(load_w("API_NY.GDP.TOTL.RT.ZS_DS2_en_csv_v2_362.csv", "rent"),
      by = c("country", "year")
    )
}

load_sbs_2017 <- function() {
  raw <- readxl::read_excel(
    file.path(data_dir, "enterprises_eu.xlsx"),
    sheet = "Sheet 1", skip = 9, col_names = FALSE
  )
  tibble::tibble(
    country    = raw[[1]],
    n_firms_17 = suppressWarnings(as.numeric(raw[[26]]))
  ) %>%
    dplyr::filter(!is.na(country), !is.na(n_firms_17)) %>%
    harmonise_names() %>%
    dplyr::distinct(country, .keep_all = TRUE)
}

# ---- 3. Panel builder (LAGGED specification) --------------------------------
#
# Difference from the simultaneous spec (estimate.R):
#   After interpolating and log-transforming covariates, each x_{i,t} is
#   replaced by its one-period lag x_{i,t-1}.  The outcome y_t = Delta log
#   n_cert_t is unchanged.  The panel is then restricted to balanced years.
#   Net cost: 1 extra period lost, so T_lag = T_sim - 1.

finalize_panel_lagged <- function(iso_full, filter_list, label) {
  res <- iso_full %>%
    dplyr::filter(country %in% filter_list) %>%
    dplyr::left_join(sbs_firms, by = "country") %>%
    dplyr::group_by(country) %>%
    dplyr::arrange(year, .by_group = TRUE) %>%
    # 1. Interpolate / fill covariates and log-transform
    dplyr::mutate(dplyr::across(c(gdp_pc, trade, rol, renew, rent), ~ {
      v <- zoo::na.approx(., na.rm = FALSE)
      v <- zoo::na.locf(zoo::na.locf(v, na.rm = FALSE),
        na.rm = FALSE,
        fromLast = TRUE
      )
      log1p(v)
    })) %>%
    # 2. Lag ALL covariates by 1 period
    #    x_{i,t-1} enters the RHS of the structural equation at time t.
    dplyr::mutate(
      x1 = dplyr::lag(gdp_pc, n = 1, order_by = year),
      x2 = dplyr::lag(trade, n = 1, order_by = year),
      x3 = dplyr::lag(rol, n = 1, order_by = year),
      x4 = dplyr::lag(renew, n = 1, order_by = year),
      x5 = dplyr::lag(rent, n = 1, order_by = year)
    ) %>%
    # 3. Outcome: Delta log n_cert (unchanged; dated t)
    dplyr::mutate(y = log1p(n_cert) - log1p(dplyr::lag(n_cert,
      n = 1,
      order_by = year
    ))) %>%
    dplyr::ungroup() %>%
    tidyr::drop_na(y, x1, x2, x3, x4, x5)

  # 4. Balance: keep only countries with the maximum T after dropping NAs
  country_T <- res %>% dplyr::count(country)
  T_max <- max(country_T$n)
  final <- res %>%
    dplyr::group_by(country) %>%
    dplyr::filter(dplyr::n() == T_max) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(id = as.integer(factor(country))) %>%
    dplyr::select(id, time = year, y, x1, x2, x3, x4, x5)

  N <- length(unique(final$id))
  cat(sprintf("Dataset %s (lagged): N=%d, T=%d\n", label, N, T_max))
  return(final)
}

# ---- 4. Build panels --------------------------------------------------------
sbs_firms <- load_sbs_2017()

h9001 <- read_iso_hist(
  "01. ISO 9001 - Number of certificates per country and industry sectors - 1993 to 2017.xlsm",
  4:8, c("country", as.character(1993:2017))
)
r9001 <- extract_recent("9001")
f9001 <- dplyr::bind_rows(h9001, r9001) %>%
  dplyr::left_join(make_covs(1993), by = c("country", "year")) %>%
  dplyr::filter(year <= 2017)

h14001 <- read_iso_hist(
  "02. ISO 14001 - Number of certificates per country and industry sectors - 1999 to 2017.xlsm",
  4:8, c("country", as.character(1999:2017))
)
r14001 <- extract_recent("14001")
f14001 <- dplyr::bind_rows(h14001, r14001) %>%
  dplyr::left_join(make_covs(1999), by = c("country", "year")) %>%
  dplyr::filter(year <= 2017)

sync_eu <- intersect(
  f9001 %>% dplyr::filter(country %in% eu) %>% dplyr::pull(country) %>% unique(),
  f14001 %>% dplyr::filter(country %in% eu) %>% dplyr::pull(country) %>% unique()
)

# Save panels (same clean_dir as main estimation for shared use downstream)
saveRDS(f9001, file.path(clean_dir, "f9001.rds"))
saveRDS(f14001, file.path(clean_dir, "f14001.rds"))
saveRDS(sync_eu, file.path(clean_dir, "sync_eu.rds"))

dp_9001_eu <- finalize_panel_lagged(f9001, sync_eu, "9001_EU")
dp_14001_eu <- finalize_panel_lagged(f14001, sync_eu, "14001_EU")

saveRDS(dp_9001_eu, file.path(clean_dir, "panel_9001_eu_lagged.rds"))
saveRDS(dp_14001_eu, file.path(clean_dir, "panel_14001_eu_lagged.rds"))
readr::write_csv(dp_9001_eu, file.path(clean_dir, "panel_9001_eu_lagged.csv"))
readr::write_csv(dp_14001_eu, file.path(clean_dir, "panel_14001_eu_lagged.csv"))

# ---- 5. Reduced lambda grid -------------------------------------------------
# Original: 4^3 = 64 combinations.
# Reduced:  3^3 = 27 combinations -- covers the same [0.05, 0.25] range
# with coarser steps.  Further reduced to the 9 most empirically relevant
# combinations if time is critical (set FAST_MODE <- TRUE below).

FAST_MODE <- FALSE # set TRUE for 9-run quick check (3 diagonal + 6 off-diag)

if (FAST_MODE) {
  # 9 targeted combinations informed by the simultaneous-spec BIC winners
  lambda_grid <- list(
    c(0.05, 0.05, 0.05), c(0.05, 0.15, 0.05), c(0.05, 0.25, 0.05),
    c(0.15, 0.05, 0.15), c(0.15, 0.15, 0.15), c(0.15, 0.25, 0.15),
    c(0.25, 0.05, 0.25), c(0.25, 0.15, 0.25), c(0.25, 0.25, 0.25)
  )
} else {
  lambda_vals <- c(0.05, 0.15, 0.25) # 3 values -> 3^3 = 27 runs
  lambda_grid <- as.list(as.data.frame(
    t(expand.grid(lambda_vals, lambda_vals, lambda_vals))
  ))
  names(lambda_grid) <- NULL
}

cat(sprintf("\nGrid size: %d combinations per panel\n", length(lambda_grid)))

# ---- 6. Grid search ---------------------------------------------------------
panels <- list(
  "9001_eu"  = dp_9001_eu,
  "14001_eu" = dp_14001_eu
)

grid_results <- list()
idx <- 0

for (pname in names(panels)) {
  df <- panels[[pname]]
  Wder <- NULL # reuse derivative across lambdas

  for (li in seq_along(lambda_grid)) {
    lam <- lambda_grid[[li]]
    idx <- idx + 1
    run_label <- sprintf("%s_lag_gridrun%02d", pname, li)

    cat(sprintf(
      "\n[%d/%d] %s | lambda = (%.3f, %.3f, %.3f) ... ",
      idx, length(lambda_grid) * length(panels),
      pname, lam[1], lam[2], lam[3]
    ))

    row <- list(
      panel = pname, lambda = lam, run_label = run_label,
      edges = NA, rho_en = NA, rho_ada = NA, rho_gmm = NA,
      bic = NA, rn = NULL
    )

    tryCatch(
      {
        set.seed(1861)
        rn <- if (is.null(Wder)) {
          recoverNetwork(df, lambda = lam)
        } else {
          recoverNetwork(df, lambda = lam, Wder = Wder)
        }
        if (is.null(Wder)) Wder <- rn$Wder

        W <- as.matrix(rn$unpenalisedgmm$W)
        diag(W) <- 0

        row$edges <- sum(abs(W) > 1e-5)
        row$rho_en <- rn$en$rho
        row$rho_ada <- rn$adaen$rho
        row$rho_gmm <- rn$unpenalisedgmm$rho
        row$bic <- rn$BIC
        row$rn <- rn

        cat(sprintf(
          "OK | edges=%d | rho_gmm=%.3f | BIC=%.2f",
          row$edges, row$rho_gmm, row$bic
        ))

        saveRDS(rn, file.path(out_dir, paste0(run_label, ".rds")))
      },
      error = function(e) {
        cat("FAILED:", conditionMessage(e))
      }
    )

    grid_results[[idx]] <- row
  }

  Wder <- NULL # reset between panels
}

# ---- 7. Save grid metadata --------------------------------------------------
grid_meta <- lapply(grid_results, function(r) {
  r$rn <- NULL
  r
})
saveRDS(grid_meta, file.path(out_dir, "grid_search_meta_lagged.rds"))
save(grid_results, file = file.path(out_dir, "grid_search_full_lagged.RData"))

# ---- 8. Summary table -------------------------------------------------------
summary_df <- do.call(rbind, lapply(grid_results, function(r) {
  data.frame(
    Panel = gsub("_", " ", r$panel),
    L1 = r$lambda[1], L2 = r$lambda[2], L1s = r$lambda[3],
    Edges = ifelse(is.na(r$edges), "FAIL", as.character(r$edges)),
    rho_EN = ifelse(is.na(r$rho_en), "FAIL", sprintf("%.3f", r$rho_en)),
    rho_Ada = ifelse(is.na(r$rho_ada), "FAIL", sprintf("%.3f", r$rho_ada)),
    rho_GMM = ifelse(is.na(r$rho_gmm), "FAIL", sprintf("%.3f", r$rho_gmm)),
    BIC = ifelse(is.na(r$bic), "FAIL", sprintf("%.2f", r$bic)),
    stringsAsFactors = FALSE
  )
}))
summary_df$Panel <- gsub("9001", "Quality (9001)", summary_df$Panel)
summary_df$Panel <- gsub("14001", "Environment (14001)", summary_df$Panel)

readr::write_csv(summary_df, file.path(out_dir, "grid_search_summary_lagged.csv"))
saveRDS(summary_df, file.path(out_dir, "grid_search_summary_lagged.rds"))

# ---- 9. Best per panel (min BIC) --------------------------------------------
best <- list()
for (pname in names(panels)) {
  panel_runs <- Filter(function(r) r$panel == pname && !is.null(r$rn), grid_results)
  if (length(panel_runs) > 0) {
    bics <- sapply(panel_runs, `[[`, "bic")
    best[[pname]] <- panel_runs[[which.min(bics)]]
    cat(sprintf(
      "\nBest [%s]: lambda=(%.3f,%.3f,%.3f) | BIC=%.2f | edges=%d | rho_gmm=%.4f\n",
      pname,
      best[[pname]]$lambda[1], best[[pname]]$lambda[2], best[[pname]]$lambda[3],
      best[[pname]]$bic, best[[pname]]$edges, best[[pname]]$rho_gmm
    ))
  }
}
saveRDS(best, file.path(out_dir, "best_per_panel_lagged.rds"))

# ---- 10. Network statistics for best networks -------------------------------
net_stats <- function(rn) {
  W <- as.matrix(rn$unpenalisedgmm$W)
  diag(W) <- 0
  B <- (abs(W) > 1e-5) * 1
  g <- igraph::graph_from_adjacency_matrix(B, mode = "directed")
  list(
    edges       = sum(B),
    clustering  = igraph::transitivity(g, "global"),
    reciprocity = sum(B * t(B)) / max(sum(B), 1),
    out_mean    = mean(igraph::degree(g, mode = "out")),
    out_sd      = sd(igraph::degree(g, mode = "out")),
    in_mean     = mean(igraph::degree(g, mode = "in")),
    in_sd       = sd(igraph::degree(g, mode = "in"))
  )
}

stats_df <- do.call(rbind, lapply(names(best), function(pname) {
  s <- net_stats(best[[pname]]$rn)
  data.frame(
    panel = pname,
    lambda1 = best[[pname]]$lambda[1],
    lambda2 = best[[pname]]$lambda[2],
    lambda1s = best[[pname]]$lambda[3],
    bic = best[[pname]]$bic,
    edges = s$edges,
    clustering = s$clustering,
    reciprocity = s$reciprocity,
    out_mean = s$out_mean, out_sd = s$out_sd,
    in_mean = s$in_mean, in_sd = s$in_sd,
    stringsAsFactors = FALSE
  )
}))
readr::write_csv(stats_df, file.path(out_dir, "network_stats_best_lagged.csv"))
saveRDS(stats_df, file.path(out_dir, "network_stats_best_lagged.rds"))

# ---- 11. Save adjacency matrices --------------------------------------------
for (pname in names(best)) {
  W <- as.matrix(best[[pname]]$rn$unpenalisedgmm$W)
  diag(W) <- 0
  saveRDS(W, file.path(out_dir, sprintf("W_best_%s_lagged.rds", pname)))
  utils::write.csv(W, file.path(out_dir, sprintf("W_best_%s_lagged.csv", pname)),
    row.names = FALSE
  )
}

cat("\n\nDone. All lagged-spec outputs written to:", normalizePath(out_dir), "\n")
