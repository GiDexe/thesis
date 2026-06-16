# =============================================================================
# ISO Diffusion in the EU: Endogenous Network Recovery for ISO 9001 / 14001
# Author: Gianluigi De Rubertis
# Method: De Paula, Rasul & Souza (2025, REStud) -- recoverNetwork
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
  library(ggplot2)
  library(readr)
  library(igraph)
  library(remotes)
  library(recoverNetwork)
  library(zoo)
})

set.seed(1861)
setwd("/workspaces/thesis")
options(readr.show_col_types = FALSE)

data_dir <- "/workspaces/thesis/data"
raw_dir <- file.path(data_dir, "data_raw")
clean_dir <- file.path(data_dir, "cleaned_data")
if (!dir.exists(clean_dir)) dir.create(clean_dir, recursive = TRUE)

out_dir <- if (Sys.getenv("RERUN") == "true") file.path(data_dir, "output_rerun") else file.path(data_dir, "output")
dir.create(out_dir, showWarnings = FALSE)

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
    "Czechia" ~ "Czech Republic", "Slovakia" ~ "Slovak Republic",
    "United Kingdom of Great Britain and Northern Ireland" ~ "United Kingdom",
    "United States of America" ~ "United States",
    "Netherlands (Kingdom of the)" ~ "Netherlands",
    .default = country
  ))
}

read_iso_hist <- function(file_name, sheets_idx, my_names) {
  f_path <- file.path(raw_dir, file_name)
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

make_covs <- function(start_yr) {
  load_w <- function(f, v) {
    readr::read_csv(file.path(raw_dir, f), skip = 4, show_col_types = FALSE) %>%
      dplyr::select(country = 1, as.character(start_yr):"2024") %>%
      tidyr::pivot_longer(-country, names_to = "year", values_to = v) %>%
      dplyr::mutate(
        year = as.integer(year),
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
  file_path <- file.path(raw_dir, "enterprises_eu.xlsx")
  raw <- readxl::read_excel(file_path, sheet = "Sheet 1", skip = 9, col_names = FALSE)
  tibble::tibble(
    country    = raw[[1]],
    n_firms_17 = suppressWarnings(as.numeric(raw[[26]]))
  ) %>%
    dplyr::filter(!is.na(country), !is.na(n_firms_17)) %>%
    harmonise_names() %>%
    dplyr::distinct(country, .keep_all = TRUE)
}

# Standardisation: log-difference of certifications
finalize_panel <- function(iso_full, filter_list, label) {
  res <- iso_full %>%
    dplyr::filter(country %in% filter_list) %>%
    dplyr::left_join(sbs_firms, by = "country") %>%
    dplyr::group_by(country) %>%
    dplyr::arrange(year, .by_group = TRUE) %>%
    dplyr::mutate(dplyr::across(c(gdp_pc, trade, rol, renew, rent), ~ {
      v <- zoo::na.approx(., na.rm = FALSE)
      v <- zoo::na.locf(zoo::na.locf(v, na.rm = FALSE), na.rm = FALSE, fromLast = TRUE)
      log1p(v)
    })) %>%
    dplyr::mutate(y = log1p(n_cert) - log1p(dplyr::lag(n_cert))) %>%
    dplyr::filter(!is.na(y)) %>%
    dplyr::ungroup()

  T_max <- max((res %>% dplyr::count(country))$n)
  final <- res %>%
    dplyr::group_by(country) %>%
    dplyr::filter(dplyr::n() == T_max) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(id = as.integer(factor(country))) %>%
    dplyr::select(id,
      time = year, y,
      x1 = gdp_pc, x2 = trade, x3 = rol, x4 = renew, x5 = rent
    ) %>%
    tidyr::drop_na()

  cat(sprintf("Dataset %s: N=%d, T=%d\n", label, length(unique(final$id)), T_max))
  return(final)
}

# ---- 3. Build panels --------------------------------------------------------
sbs_firms <- load_sbs_2017()

f9001 <- read_iso_hist(
  "01. ISO 9001 - Number of certificates per country and industry sectors - 1993 to 2017.xlsm",
  4:8, c("country", as.character(1993:2017))
)


f14001 <- read_iso_hist(
  "02. ISO 14001 - Number of certificates per country and industry sectors - 1999 to 2017.xlsm",
  4:8, c("country", as.character(1999:2017))
)


covs9001 <- make_covs(1993)
covs14001 <- make_covs(1999)
f9001 <- f9001 %>% dplyr::left_join(covs9001, by = c("country", "year"))
f14001 <- f14001 %>% dplyr::left_join(covs14001, by = c("country", "year"))

# ---- EU-27 share of global ISO certificates (role of EU in diffusion)
eu_global_share <- function(df, yr, eu_set) {
  glob <- df %>%
    dplyr::filter(year == yr) %>%
    dplyr::summarise(t = sum(n_cert, na.rm = TRUE)) %>%
    dplyr::pull(t)
  eu_tot <- df %>%
    dplyr::filter(year == yr, country %in% eu_set) %>%
    dplyr::summarise(t = sum(n_cert, na.rm = TRUE)) %>%
    dplyr::pull(t)
  data.frame(year = yr, eu = eu_tot, global = glob, share = eu_tot / glob)
}

share_9001 <- do.call(rbind, lapply(c(1999, 2010, 2017), eu_global_share, df = f9001, eu_set = eu))
share_14001 <- do.call(rbind, lapply(c(1999, 2010, 2017), eu_global_share, df = f14001, eu_set = eu))

cat("\n--- ISO 9001 EU-27 global share ---\n")
print(share_9001)
cat("\n--- ISO 14001 EU-27 global share ---\n")
print(share_14001)

saveRDS(
  list(iso9001 = share_9001, iso14001 = share_14001),
  file.path(clean_dir, "eu_global_share.rds")
)
# Truncate at 2017 (overlap of both standards' historical sheets)
f9001 <- f9001 %>% dplyr::filter(year <= 2017)
f14001 <- f14001 %>% dplyr::filter(year <= 2017)

# Sync countries across standards
sync_eu <- intersect(
  (f9001 %>% dplyr::filter(country %in% eu) %>% dplyr::pull(country) %>% unique()),
  (f14001 %>% dplyr::filter(country %in% eu) %>% dplyr::pull(country) %>% unique())
)

saveRDS(f9001, file.path(clean_dir, "f9001.rds"))
saveRDS(f14001, file.path(clean_dir, "f14001.rds"))
saveRDS(sync_eu, file.path(clean_dir, "sync_eu.rds"))

dp_9001_eu <- finalize_panel(f9001, sync_eu, "9001_EU (synced)")
dp_14001_eu <- finalize_panel(f14001, sync_eu, "14001_EU (synced)")

clean_dir <- "/workspaces/thesis/data/cleaned_data"
if (!dir.exists(clean_dir)) dir.create(clean_dir, recursive = TRUE)

saveRDS(dp_9001_eu, file.path(clean_dir, "panel_9001_eu.rds"))
saveRDS(dp_14001_eu, file.path(clean_dir, "panel_14001_eu.rds"))
readr::write_csv(dp_9001_eu, file.path(clean_dir, "panel_9001_eu.csv"))
readr::write_csv(dp_14001_eu, file.path(clean_dir, "panel_14001_eu.csv"))

# ---- 4. Grid search ---------------------------------------------------------
panels <- list(
  "9001_eu"  = dp_9001_eu,
  "14001_eu" = dp_14001_eu
)

lambda_vals <- c(0.05, 0.10, 0.15, 0.2)
lambda_grid <- as.list(as.data.frame(t(expand.grid(
  lambda_vals, lambda_vals, lambda_vals
))))
names(lambda_grid) <- NULL

grid_results <- list()
idx <- 0

for (pname in names(panels)) {
  df <- panels[[pname]]
  Wder <- NULL

  for (li in seq_along(lambda_grid)) {
    lam <- lambda_grid[[li]]
    idx <- idx + 1
    run_label <- sprintf("%s_gridrun%02d", pname, li)
    cat(sprintf(
      "\n[%d] %s | lambda = (%.3f, %.3f, %.3f) ... ",
      idx, pname, lam[1], lam[2], lam[3]
    ))

    row <- list(
      panel = pname, lambda = lam, run_label = run_label,
      edges = NA, rho_en = NA, rho_ada = NA, rho_gmm = NA, bic = NA, rn = NULL
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

  Wder <- NULL
}

# Save grid
grid_meta <- lapply(grid_results, function(r) {
  r$rn <- NULL
  r
})
saveRDS(grid_meta, file.path(out_dir, "grid_search_meta.rds"))
save(grid_results, file = file.path(out_dir, "grid_search_full.RData"))

# ---- 5. Grid summary table --------------------------------------------------
summary_df <- do.call(rbind, lapply(grid_results, function(r) {
  data.frame(
    Panel = gsub("_", " ", r$panel),
    L1 = r$lambda[1],
    L2 = r$lambda[2],
    L1s = r$lambda[3],
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

readr::write_csv(summary_df, file.path(out_dir, "grid_search_summary.csv"))
saveRDS(summary_df, file.path(out_dir, "grid_search_summary.rds"))

# ---- 6. Best run per panel (min BIC) ----------------------------------------
best <- list()
for (pname in names(panels)) {
  panel_runs <- Filter(function(r) r$panel == pname && !is.null(r$rn), grid_results)
  if (length(panel_runs) > 0) {
    bics <- sapply(panel_runs, `[[`, "bic")
    best[[pname]] <- panel_runs[[which.min(bics)]]
  }
}
saveRDS(best, file.path(out_dir, "best_per_panel.rds"))

# ---- 7. Network statistics for best networks --------------------------------
# Convention: binarize via abs(W) > 1e-5 (numerical tolerance from Caner & Zhang 2014,
# adopted by de Paula, Rasul & Souza 2025, fn. 25). Edges are defined on |W|, not W,
# so the count is sign-agnostic; in this empirical setting non-negativity is implicit
# given row-sum normalization (de Paula et al. 2025, A4').
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
    out_mean = s$out_mean,
    out_sd = s$out_sd,
    in_mean = s$in_mean,
    in_sd = s$in_sd,
    stringsAsFactors = FALSE
  )
}))

readr::write_csv(stats_df, file.path(out_dir, "network_stats_best.csv"))
saveRDS(stats_df, file.path(out_dir, "network_stats_best.rds"))

# ---- 8. Save adjacency matrices of best networks ----------------------------
for (pname in names(best)) {
  W <- as.matrix(best[[pname]]$rn$unpenalisedgmm$W)
  diag(W) <- 0
  saveRDS(W, file.path(out_dir, sprintf("W_best_%s.rds", pname)))
  utils::write.csv(W, file.path(out_dir, sprintf("W_best_%s.csv", pname)),
    row.names = FALSE
  )
}

cat("\n\nDone. All outputs written to:", normalizePath(out_dir), "\n")
