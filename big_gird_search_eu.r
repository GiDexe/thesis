# =============================================================================
# grid_search.R
# Symmetric grid search over lambda + De Paula Table 3 for best results
# Run AFTER master_cleaning_final_all_4.R (needs dp_* objects in environment)
# =============================================================================

library(recoverNetwork)
library(igraph)

out_dir <- "output"
if (!dir.exists(out_dir)) dir.create(out_dir)

# === GRID SETUP ===
panels <- list(
  "9001_eu"     = dp_9001_eu,
  "14001_eu"    = dp_14001_eu,
  "9001_oecd"   = dp_9001_oecd,
  "14001_oecd"  = dp_14001_oecd
)

lambda_sym <- c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30)

# === GRID SEARCH ===
grid_results <- list()
idx <- 0

for (pname in names(panels)) {
  df <- panels[[pname]]
  Wder <- NULL

  for (l in lambda_sym) {
    lam <- c(l, l, l)
    idx <- idx + 1
    run_label <- sprintf("%s_gridrun%02d", pname, idx)
    cat(sprintf("\n[%d] %s | lambda = %.2f ... ", idx, pname, l))

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

        # Stats from the FINAL W matrix (unpenalised GMM)
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

        # Save individual run
        saveRDS(rn, file.path(out_dir, paste0(run_label, ".rds")))
      },
      error = function(e) {
        cat("FAILED:", conditionMessage(e))
      }
    )

    grid_results[[idx]] <- row
  }

  Wder <- NULL # reset for next panel
}

# Save full grid results (without rn objects, to keep file small)
grid_meta <- lapply(grid_results, function(r) {
  r$rn <- NULL
  r
})
saveRDS(grid_meta, file.path(out_dir, "grid_search_meta.rds"))

# Also save full object for recomputation
save(grid_results, file = file.path(out_dir, "grid_search_full.RData"))

# === CONSOLE SUMMARY ===
cat("\n\n====== PILOT GRID SEARCH SUMMARY ======\n")
cat(sprintf(
  "%-12s %6s %6s %6s %6s %8s %8s %8s %8s\n",
  "Panel", "L1", "L2", "L1*", "Edges", "rho_EN", "rho_Ada", "rho_GMM", "BIC"
))
for (r in grid_results) {
  cat(sprintf(
    "%-12s %6.2f %6.2f %6.2f %6s %8s %8s %8s %8s\n",
    r$panel, r$lambda[1], r$lambda[2], r$lambda[3],
    ifelse(is.na(r$edges), "FAIL", as.character(r$edges)),
    ifelse(is.na(r$rho_en), "FAIL", sprintf("%.3f", r$rho_en)),
    ifelse(is.na(r$rho_ada), "FAIL", sprintf("%.3f", r$rho_ada)),
    ifelse(is.na(r$rho_gmm), "FAIL", sprintf("%.3f", r$rho_gmm)),
    ifelse(is.na(r$bic), "FAIL", sprintf("%.2f", r$bic))
  ))
}

# === SELECT BEST PER PANEL (by BIC) ===
best <- list()
for (pname in names(panels)) {
  panel_runs <- Filter(function(r) r$panel == pname && !is.null(r$rn), grid_results)
  if (length(panel_runs) > 0) {
    bics <- sapply(panel_runs, `[[`, "bic")
    b <- panel_runs[[which.min(bics)]]
    best[[pname]] <- b
    cat(sprintf(
      "\nBest %s: lambda=(%.2f,%.2f,%.2f), BIC=%.2f, edges=%d, rho_gmm=%.3f\n",
      pname, b$lambda[1], b$lambda[2], b$lambda[3],
      b$bic, b$edges, b$rho_gmm
    ))
  }
}

saveRDS(best, file.path(out_dir, "best_per_panel.rds"))

# === NETWORK STATISTICS (from final GMM W matrix) ===
net_stats <- function(rn) {
  W <- as.matrix(rn$unpenalisedgmm$W)
  diag(W) <- 0
  B <- (abs(W) > 1e-5) * 1 # binary adjacency
  g <- graph_from_adjacency_matrix(B, mode = "directed")
  list(
    edges = sum(B),
    clustering = transitivity(g, "global"),
    reciprocity = sum(B * t(B)) / max(sum(B), 1),
    out_mean = mean(degree(g, mode = "out")),
    out_sd = sd(degree(g, mode = "out")),
    in_mean = mean(degree(g, mode = "in")),
    in_sd = sd(degree(g, mode = "in"))
  )
}

# === GENERATE LATEX TABLES ===
# One table per group (EU, OECD), comparing 9001 vs 14001 — De Paula Table 3 style

write_comparison_table <- function(best_a, best_b, label_a, label_b,
                                   group_label, filename) {
  rn_a <- best_a$rn
  rn_b <- best_b$rn

  s_a <- net_stats(rn_a)
  s_b <- net_stats(rn_b)

  # Binary adjacency for overlap computation
  W_a <- (abs(as.matrix(rn_a$unpenalisedgmm$W)) > 1e-5) * 1
  diag(W_a) <- 0
  W_b <- (abs(as.matrix(rn_b$unpenalisedgmm$W)) > 1e-5) * 1
  diag(W_b) <- 0
  common <- sum(W_a & W_b)
  only_a <- sum(W_a & !W_b)
  only_b <- sum(!W_a & W_b)

  fmt <- function(x, d = 3) formatC(x, format = "f", digits = d)
  pct <- function(x) paste0(fmt(100 * x, 1), "\\%")
  degr <- function(m, s) paste0(fmt(m, 3), " (", fmt(s, 3), ")")

  tex <- c(
    "\\begin{table}[htbp]",
    "\\centering",
    sprintf("\\caption{%s: ISO 9001 versus ISO 14001 networks}", group_label),
    sprintf("\\label{tab:network_%s}", tolower(gsub(" ", "_", group_label))),
    "\\begin{tabular}{lcc}",
    "\\hline\\hline",
    sprintf("& %s & %s \\\\", label_a, label_b),
    "\\hline",
    sprintf("Number of edges & %d & %d \\\\", s_a$edges, s_b$edges),
    sprintf("Edges in both networks & %d & %d \\\\", common, common),
    sprintf("Edges in W\\textsubscript{9001} only & %d & \\\\", only_a),
    sprintf("Edges in W\\textsubscript{14001} only & & %d \\\\", only_b),
    sprintf("Clustering & %s & %s \\\\", fmt(s_a$clustering), fmt(s_b$clustering)),
    sprintf("Reciprocated edges & %s & %s \\\\", pct(s_a$reciprocity), pct(s_b$reciprocity)),
    "Degree distribution across nodes (countries) & & \\\\",
    sprintf(
      "\\quad Out-degree & %s & %s \\\\",
      degr(s_a$out_mean, s_a$out_sd), degr(s_b$out_mean, s_b$out_sd)
    ),
    sprintf(
      "\\quad In-degree & %s & %s \\\\",
      degr(s_a$in_mean, s_a$in_sd), degr(s_b$in_mean, s_b$in_sd)
    ),
    "\\hline",
    sprintf(
      "$\\hat{\\rho}$ (GMM) & %s & %s \\\\",
      fmt(rn_a$unpenalisedgmm$rho, 4), fmt(rn_b$unpenalisedgmm$rho, 4)
    ),
    sprintf(
      "$\\lambda$ & (%s) & (%s) \\\\",
      paste(fmt(best_a$lambda, 2), collapse = ", "),
      paste(fmt(best_b$lambda, 2), collapse = ", ")
    ),
    sprintf("BIC & %s & %s \\\\", fmt(best_a$bic, 2), fmt(best_b$bic, 2)),
    "\\hline\\hline",
    "\\end{tabular}",
    "\\end{table}"
  )

  writeLines(tex, file.path(out_dir, filename))
  cat("Written:", filename, "\n")
}

# EU table
if (!is.null(best[["9001_eu"]]) && !is.null(best[["14001_eu"]])) {
  write_comparison_table(
    best[["9001_eu"]], best[["14001_eu"]],
    "ISO 9001", "ISO 14001",
    "EU", "table_network_eu.tex"
  )
}

# OECD table
if (!is.null(best[["9001_oecd"]]) && !is.null(best[["14001_oecd"]])) {
  write_comparison_table(
    best[["9001_oecd"]], best[["14001_oecd"]],
    "ISO 9001", "ISO 14001",
    "OECD", "table_network_oecd.tex"
  )
}

cat("\n=== GRID SEARCH COMPLETE ===\n")
cat("Saved files:\n")
cat("  output/grid_search_full.RData     (all runs, for recomputation)\n")
cat("  output/grid_search_meta.rds       (metadata without rn objects)\n")
cat("  output/best_per_panel.rds         (best run per panel)\n")
cat("  output/<panel>_gridrunNN.rds      (individual runs)\n")
cat("  output/table_network_eu.tex       (De Paula Table 3, EU)\n")
cat("  output/table_network_oecd.tex     (De Paula Table 3, OECD)\n")
