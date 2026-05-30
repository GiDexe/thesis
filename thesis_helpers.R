## =============================================================================
## thesis_helpers.R
## -----------------------------------------------------------------------------
## All helper functions used in main2.qmd, documented in one place.
##
## Sections:
##   1. Distance-matrix builders (WDI-style CSV, Hofstede culture, capitals,
##      env tax, FDI)
##   2. ABC exogenous network construction (Albuquerque, Bronnenberg & Corbett
##      2007 top-N similarity / distance networks)
##   3. BDF (2009) panel builders, outcome attachers, IV formulas
##   4. Inference helpers: pairs cluster bootstrap, robust extractors,
##      significance stars, formatters
##   5. Table-printing helpers (kableExtra wrappers)
##   6. Diagnostics: VIF, network statistics
##
## Author: Gianluigi De Rubertis
## Last updated: 2026-05-12
## =============================================================================


## =============================================================================
## 1. Distance-matrix builders
## =============================================================================

#' Convert ISO2 country code to ISO3, handling Eurostat quirks.
#'
#' Eurostat uses "EL" for Greece and "UK" for the United Kingdom rather than
#' the ISO2 codes "GR" and "GB". This helper applies those substitutions
#' before delegating to countrycode::countrycode().
#'
#' @param x character vector of ISO2 codes (possibly Eurostat-style)
#' @return character vector of ISO3 codes
iso2_to_iso3 <- function(x) {
  x <- as.character(x)
  x <- dplyr::case_match(x, "EL" ~ "GR", "UK" ~ "GB", .default = x)
  countrycode::countrycode(x,
    origin = "iso2c", destination = "iso3c", warn = FALSE
  )
}


#' Haversine distance between two points on the Earth's surface.
#'
#' @param lat1,lon1,lat2,lon2 latitudes and longitudes in degrees
#' @return great-circle distance in kilometres
haversine_km <- function(lat1, lon1, lat2, lon2) {
  R <- 6371
  to_rad <- pi / 180
  dlat <- (lat2 - lat1) * to_rad
  dlon <- (lon2 - lon1) * to_rad
  a <- sin(dlat / 2)^2 +
    cos(lat1 * to_rad) * cos(lat2 * to_rad) * sin(dlon / 2)^2
  2 * R * asin(pmin(1, sqrt(a)))
}


#' Build a dyadic distance matrix from a WDI-style country-year wide CSV.
#'
#' Reads a World Bank Indicators CSV (the standard 4-skip-row format),
#' averages the indicator over the specified year columns, and returns the
#' N x N matrix of absolute log-differences between country values. If a
#' country in country_labels is missing in the CSV, it is imputed with the
#' sample mean of available countries.
#'
#' @param path path to a WDI CSV
#' @param year_cols character vector of year column names, e.g.
#'   c("x1996", "x1997", "x1998", "x1999", "x2000")
#' @param take_log if TRUE (default), apply log(x + 1) before differencing
#' @return N x N symmetric matrix with country labels as row/col names
build_dist_wdi <- function(path, year_cols, take_log = TRUE) {
  raw <- readr::read_csv(path, skip = 4, show_col_types = FALSE) |>
    janitor::clean_names()
  available <- intersect(year_cols, names(raw))
  if (length(available) == 0) {
    available <- grep("^x[0-9]{4}$", names(raw), value = TRUE)[1:5]
  }
  d <- raw |>
    dplyr::transmute(
      iso3 = country_code,
      val = rowMeans(
        dplyr::across(
          dplyr::all_of(available),
          \(z) suppressWarnings(as.numeric(z))
        ),
        na.rm = TRUE
      )
    ) |>
    dplyr::filter(iso3 %in% cl_iso3, !is.na(val), is.finite(val)) |>
    dplyr::mutate(country = iso_to_label[iso3]) |>
    dplyr::arrange(match(country, country_labels))
  if (nrow(d) < N) {
    miss <- setdiff(country_labels, d$country)
    d <- dplyr::bind_rows(
      d,
      data.frame(
        iso3 = NA, val = mean(d$val, na.rm = TRUE),
        country = miss, stringsAsFactors = FALSE
      )
    ) |>
      dplyr::arrange(match(country, country_labels))
  }
  v <- if (take_log) log(d$val + 1) else d$val
  names(v) <- d$country
  M <- abs(outer(v, v, FUN = "-"))
  dimnames(M) <- list(country_labels, country_labels)
  M
}


## =============================================================================
## 2. ABC exogenous network construction
##    (Albuquerque, Bronnenberg & Corbett 2007, Table 3 Model 5)
## =============================================================================

#' Construct an exogenous Top-N directed network from a similarity or
#' distance matrix.
#'
#' For each country i, identifies the top_n most similar (or closest)
#' countries and assigns uniform weight 1/top_n to those edges; the diagonal
#' is set to -Inf or +Inf as appropriate so a country cannot be its own
#' neighbour. This is the ABC (2007) construction.
#'
#' @param M N x N matrix (similarity or distance)
#' @param top_n number of neighbours to keep per row (default 3)
#' @param similarity if TRUE, M is a similarity matrix and we pick the
#'   highest values; if FALSE, M is a distance and we pick the lowest values
#' @return N x N row-normalised binary-weighted network
top_n_network <- function(M, top_n = 3, similarity = FALSE) {
  W <- matrix(0, nrow(M), ncol(M), dimnames = dimnames(M))
  for (i in seq_len(nrow(M))) {
    v <- M[i, ]
    v[i] <- if (similarity) -Inf else Inf
    idx <- if (similarity) {
      order(v, decreasing = TRUE)[seq_len(top_n)]
    } else {
      order(v, decreasing = FALSE)[seq_len(top_n)]
    }
    W[i, idx] <- 1 / top_n
  }
  W
}


## =============================================================================
## 3. BDF (2009) panel builders, outcome attachers, IV formulas
## =============================================================================

#' Row-normalise a weighted matrix so each row sums to 1.
#'
#' Sets the diagonal to zero and entries with |w| < eps to zero before
#' normalisation. Rows that sum to zero are left as zero (isolated nodes).
#'
#' @param W N x N matrix
#' @param eps numeric threshold (default 1e-5)
#' @return N x N row-normalised matrix
row_normalise <- function(W, eps = 1e-5) {
  W[abs(W) < eps] <- 0
  diag(W) <- 0
  rs <- rowSums(W)
  rs[rs == 0] <- 1
  W / rs
}


#' Restrict a panel to years in which all countries are present.
#'
#' Keeps only year-rows where the number of countries equals
#' length(country_labels). Use for full-panel BDF estimation.
#'
#' @param d a data frame with columns `country`, `year`
#' @return d filtered to balanced years
keep_full <- function(d) {
  full_years <- d |>
    dplyr::count(year) |>
    dplyr::filter(n == length(country_labels)) |>
    dplyr::pull(year)
  dplyr::filter(d, year %in% full_years)
}


#' Attach intensity outcome y = n_cert / n_firms_17 * 1000.
#'
#' Intensity (certifications per thousand firms at end of sample) is the
#' Albuquerque, Bronnenberg & Corbett (2007) outcome variable.
#'
#' @param dp_raw recoverNetwork panel with id, time, x1..x5
#' @param certs ISO certification panel
#' @return tibble with country, year, x1..x5, y
attach_y_intensity <- function(dp_raw, certs) {
  dp_raw |>
    dplyr::mutate(country = id_to_country[as.character(id)], year = time) |>
    dplyr::select(country, year, x1, x2, x3, x4, x5) |>
    dplyr::inner_join(
      certs |> dplyr::filter(country %in% country_labels) |>
        dplyr::transmute(country, year, n_cert),
      by = c("country", "year")
    ) |>
    dplyr::inner_join(firm17, by = "country") |>
    dplyr::mutate(y = n_cert / n_firms_17 * 1000) |>
    dplyr::select(country, year, x1, x2, x3, x4, x5, y) |>
    tidyr::drop_na() |>
    dplyr::filter(country %in% country_labels)
}


#' Attach log-intensity outcome y = log(1 + n_cert / n_firms_17 * 1000).
attach_y_log_intensity <- function(dp_raw, certs) {
  dp_raw |>
    dplyr::mutate(country = id_to_country[as.character(id)], year = time) |>
    dplyr::select(country, year, x1, x2, x3, x4, x5) |>
    dplyr::inner_join(
      certs |> dplyr::filter(country %in% country_labels) |>
        dplyr::transmute(country, year, n_cert),
      by = c("country", "year")
    ) |>
    dplyr::inner_join(firm17, by = "country") |>
    dplyr::mutate(y = log1p(n_cert / n_firms_17 * 1000)) |>
    dplyr::select(country, year, x1, x2, x3, x4, x5, y) |>
    tidyr::drop_na() |>
    dplyr::filter(country %in% country_labels)
}


#' Attach log-difference outcome y = log(1 + n_cert_t) - log(1 + n_cert_{t-1}).
#'
#' This is the headline BDF outcome: a measure of growth in certifications
#' that absorbs country-specific level differences. Justified by
#' Albuquerque, Bronnenberg & Corbett (2007); see also Griliches (1957).
attach_y_logdiff <- function(dp_raw, certs) {
  dp_raw |>
    dplyr::mutate(country = id_to_country[as.character(id)], year = time) |>
    dplyr::select(country, year, x1, x2, x3, x4, x5) |>
    dplyr::inner_join(
      certs |> dplyr::filter(country %in% country_labels) |>
        dplyr::arrange(country, year) |>
        dplyr::group_by(country) |>
        dplyr::mutate(y = log1p(n_cert) -
          log1p(dplyr::lag(n_cert, n = 1, order_by = year))) |>
        dplyr::ungroup() |>
        dplyr::transmute(country, year, y),
      by = c("country", "year")
    ) |>
    tidyr::drop_na() |>
    dplyr::filter(country %in% country_labels)
}

#' Build the BDF (2009) regression panel from a country-year panel and a
#' network W.
#'
#' Constructs Gy = W %*% y (peer outcome), G_xk = W %*% xk (peer exogenous),
#' G2_xk = W^2 %*% xk and G3_xk = W^3 %*% xk (BDF instruments under
#' intransitivity).
#'
#' If use_lag = TRUE, applies a one-year lag uniformly to ALL regressors
#' and instruments: {Gy, x_1..x_K, G_x*, G2_x*, G3_x*}. This matches the
#' predetermined-regressor assumption E[eps_t | z_{t-1}] = 0 (Liu,
#' Patacchini & Zenou 2014, Quant. Econ. 5(1): 91-140) and avoids the
#' inconsistency of lagging only the peer outcome while leaving own and
#' contextual regressors contemporaneous. One period is lost; the effective
#' sample is T - 1.
#'
#' The contemporaneous specification (use_lag = FALSE) matches BDF (2009)
#' Proposition 1 under intransitivity.
#'
#' @param panel data frame with country, year, y, x1..x5
#' @param W N x N network matrix (will be row-normalised inside)
#' @param x_vars character vector of exogenous-variable names
#' @param use_lag logical; if TRUE applies one-year lag to ALL regressors
#'   and instruments
#' @return data frame with Gy, G_x*, G2_x*, G3_x* columns added
build_bdf <- function(panel, W, x_vars, use_lag) {
  G <- row_normalise(W)
  G2 <- G %*% G
  G3 <- G2 %*% G

  d <- panel[order(panel$year, match(panel$country, rownames(G))), ]
  yrs <- sort(unique(d$year))

  # --- peer outcome ---
  d$Gy <- unlist(lapply(yrs, function(t) {
    as.numeric(G %*% d$y[d$year == t])
  }))

  # --- contextual regressors and instruments ---
  for (v in x_vars) {
    Xt <- matrix(d[[v]], nrow = nrow(G), ncol = length(yrs))
    d[[paste0("G_", v)]] <- as.numeric(G %*% Xt)
    d[[paste0("G2_", v)]] <- as.numeric(G2 %*% Xt)
    d[[paste0("G3_", v)]] <- as.numeric(G3 %*% Xt)
  }

  d <- d[order(d$country, d$year), ]

  if (use_lag) {
    lag_cols <- c(
      "Gy",
      x_vars, # own regressors x_{t-1}
      paste0("G_", x_vars), # contextual Gx_{t-1}
      paste0("G2_", x_vars), # instruments G²x_{t-1}
      paste0("G3_", x_vars) # instruments G³x_{t-1}
    )
    d <- d |>
      dplyr::group_by(country) |>
      dplyr::mutate(
        dplyr::across(
          dplyr::all_of(lag_cols),
          \(z) dplyr::lag(z, n = 1, order_by = year)
        )
      ) |>
      dplyr::ungroup() |>
      tidyr::drop_na(dplyr::all_of(lag_cols))
  }

  d
}


## =============================================================================
## 4. Inference helpers
## =============================================================================

#' Pairs cluster bootstrap for IV models with feols().
#'
#' Resamples whole clusters (countries) with replacement and refits the IV
#' model B times. Returns the bootstrap SE and a centred two-sided p-value.
#' Recommended for samples with G < 50 clusters (Cameron, Gelbach & Miller
#' 2008, Rev. Econ. Stat. 90(3): 414-427).
#'
#' @param data data frame containing the cluster variable
#' @param fit_formula a formula compatible with fixest::feols
#' @param cluster_var name of the cluster column (default "country")
#' @param B number of bootstrap replicates (default 999)
#' @return list with `estimate`, `se`, `p_val` (named numeric vectors)
cluster_boot <- function(data, fit_formula, cluster_var = "country",
                         B = 999) {
  clusters <- unique(data[[cluster_var]])
  G <- length(clusters)
  est0 <- coef(fixest::feols(fit_formula,
    data = data, warn = FALSE, notes = FALSE
  ))
  k <- length(est0)
  boot_mat <- matrix(NA_real_,
    nrow = B, ncol = k,
    dimnames = list(NULL, names(est0))
  )
  for (b in seq_len(B)) {
    idx <- sample(clusters, G, replace = TRUE)
    boot_df <- do.call(rbind, lapply(seq_along(idx), function(j) {
      d <- data[data[[cluster_var]] == idx[j], ]
      d[[cluster_var]] <- paste0(idx[j], "_", j)
      d
    }))
    fit_b <- tryCatch(
      fixest::feols(fit_formula,
        data = boot_df, warn = FALSE, notes = FALSE
      ),
      error = function(e) NULL
    )
    if (!is.null(fit_b)) boot_mat[b, ] <- coef(fit_b)
  }
  list(
    estimate = est0,
    se = apply(boot_mat, 2, sd, na.rm = TRUE),
    p_val = sapply(seq_len(k), function(j) {
      mean(abs(boot_mat[, j] - est0[j]) >= abs(est0[j]), na.rm = TRUE)
    })
  )
}


#' Extract fit_Gy coefficient from a feols IV fit (no-FE specification).
#'
#' Returns NA values if fit_Gy is not in the coefficient table (e.g.
#' because the IV first stage collapsed under network sparsity).
extract_nofe <- function(fit) {
  ct <- summary(fit, cluster = ~country)$coeftable
  i <- which(rownames(ct) == "fit_Gy")
  if (length(i) == 0) {
    return(list(b = NA_real_, se = NA_real_, p = NA_real_))
  }
  list(
    b = ct[i, "Estimate"],
    se = ct[i, "Std. Error"],
    p = ct[i, "Pr(>|t|)"]
  )
}


#' Extract fit_Gy coefficient from a cluster_boot result (country FE
#' specification).
extract_cfe <- function(b) {
  i <- which(names(b$estimate) == "fit_Gy")
  if (length(i) == 0) {
    return(list(b = NA_real_, se = NA_real_, p = NA_real_))
  }
  list(b = b$estimate[i], se = b$se[i], p = b$p_val[i])
}


#' Convert a p-value into LaTeX significance stars.
#'
#' Returns an empty string if p is NA or zero-length.
#'
#' @param p numeric p-value
#' @return character: "$^{***}$", "$^{**}$", "$^{*}$", or ""
sig_stars <- function(p) {
  if (length(p) == 0 || is.na(p)) {
    return("")
  }
  if (p < 0.01) {
    return("$^{***}$")
  }
  if (p < 0.05) {
    return("$^{**}$")
  }
  if (p < 0.10) {
    return("$^{*}$")
  }
  ""
}


#' Format a coefficient estimate with stars; returns "---" if missing.
fmt_est <- function(b, p) {
  if (length(b) == 0 || is.na(b)) {
    return("---")
  }
  sprintf("%.3f%s", b, sig_stars(p))
}


#' Format a standard error with parentheses; returns "" if missing.
fmt_se <- function(se) {
  if (length(se) == 0 || is.na(se)) {
    return("")
  }
  sprintf("(%.3f)", se)
}


#' Build a two-row block for a regression result: estimate on top, SE on
#' bottom. Used to assemble multi-network comparison tables.
#'
#' @param panel the panel label (e.g. "Quality")
#' @param network the network label (e.g. "Recovered")
#' @param r a result list with $nofe and $cfe, each containing b, se, p
#' @return 2 x 4 character matrix
rb <- function(panel, network, r) {
  rbind(
    c(
      panel, network,
      fmt_est(r$nofe$b, r$nofe$p),
      fmt_est(r$cfe$b, r$cfe$p)
    ),
    c("", "", fmt_se(r$nofe$se), fmt_se(r$cfe$se))
  )
}


#' Single-fit BDF driver: returns no-FE and country-FE estimates.
#'
#' @param panel BDF panel built by build_bdf()
#' @param W network matrix
#' @param use_lag passes through to build_bdf
#' @param B bootstrap replicates for country FE inference
#' @return list with $nofe and $cfe
fit_one <- function(panel, W, use_lag, B = 999) {
  d <- build_bdf(panel, W, x_vars, use_lag)
  list(
    nofe = extract_nofe(feols(f_nofe, data = d, cluster = ~country)),
    cfe  = extract_cfe(cluster_boot(d, f_cfe, "country", B = B))
  )
}


#' Build a table comparing recovered vs ABC two-channel networks.
#'
#' @param panel_q,panel_e Quality and Environment outcome panels
#' @param use_lag passed to fit_one
#' @param B bootstrap replicates
#' @return data frame ready for print_kbl
table_rec_vs_abc <- function(panel_q, panel_e, use_lag, B = 999) {
  q_rec <- fit_one(panel_q, W_q, use_lag, B)
  q_abc <- fit_one(panel_q, W_geotrade, use_lag, B)
  e_rec <- fit_one(panel_e, W_e, use_lag, B)
  e_abc <- fit_one(panel_e, W_geoculture, use_lag, B)
  tab <- rbind(
    rb("Quality", "Recovered", q_rec),
    rb("Quality", "ABC: Geo + Trade", q_abc),
    rb("Environment", "Recovered", e_rec),
    rb("Environment", "ABC: Geo + Culture", e_abc)
  )
  tab <- as.data.frame(tab, stringsAsFactors = FALSE)
  colnames(tab) <- c(
    "Panel", "Network",
    "$\\hat\\beta$ (no FE)", "$\\hat\\beta$ (country FE)"
  )
  tab
}


#' Build a table for ABC single-channel networks (Geo / Trade / Culture).
table_abc_single <- function(panel_q, panel_e, use_lag, B = 999) {
  q_geo <- fit_one(panel_q, W_geo, use_lag, B)
  q_tra <- fit_one(panel_q, W_trade, use_lag, B)
  q_cul <- fit_one(panel_q, W_cul, use_lag, B)
  e_geo <- fit_one(panel_e, W_geo, use_lag, B)
  e_tra <- fit_one(panel_e, W_trade, use_lag, B)
  e_cul <- fit_one(panel_e, W_cul, use_lag, B)
  tab <- rbind(
    rb("Quality", "Geography", q_geo),
    rb("Quality", "Trade", q_tra),
    rb("Quality", "Cultural distance", q_cul),
    rb("Environment", "Geography", e_geo),
    rb("Environment", "Trade", e_tra),
    rb("Environment", "Cultural distance", e_cul)
  )
  tab <- as.data.frame(tab, stringsAsFactors = FALSE)
  colnames(tab) <- c(
    "Panel", "Network",
    "$\\hat\\beta$ (no FE)", "$\\hat\\beta$ (country FE)"
  )
  tab
}


## =============================================================================
## 5. Table-printing helpers
## =============================================================================

#' Wrapper around kableExtra::kbl with the conventions used throughout
#' the thesis: booktabs, LaTeX format, no row names, hold_position,
#' font size 10.
print_kbl <- function(tab, caption, label) {
  kbl(tab,
    booktabs = TRUE, format = "latex", escape = FALSE,
    align = "llrr", row.names = FALSE,
    caption = caption, label = label
  ) |>
    kable_styling(latex_options = "hold_position", font_size = 10)
}


## =============================================================================
## 6. Diagnostics
## =============================================================================

#' Compute basic statistics on a directed weighted network.
#'
#' Returns edge count, global transitivity (clustering), share of
#' reciprocated edges, and out- and in-degree means/SDs.
net_stats <- function(W) {
  B <- (abs(W) > 1e-5) * 1
  diag(B) <- 0
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


#' Variance Inflation Factor for one variable: regress it on the others
#' and compute 1 / (1 - R^2). VIF > 5 is the usual warning threshold,
#' VIF > 10 indicates serious collinearity (Wooldridge 2010, Econometric
#' Analysis of Cross Section and Panel Data, Sec. 4.2.4).
vif_one <- function(formula, data) {
  fit <- lm(formula, data = data)
  r2 <- summary(fit)$r.squared
  1 / (1 - r2)
}


## =============================================================================
## END OF FILE
## =============================================================================
