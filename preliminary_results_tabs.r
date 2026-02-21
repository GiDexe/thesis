#to be added to container
install.packages(c("igraph"))
library(igraph)

# === LOAD ALL RESULTS ===
e_37 <- new.env(); load("data/output_data/grid_search_ISO37001.RData", envir = e_37)
e_45 <- new.env(); load("data/output_data/grid_search_ISO45001.RData", envir = e_45)
e_37s <- new.env(); load("data/output_data/grid_search_ISO37001_standardised.RData", envir = e_37s)
e_45s <- new.env(); load("data/output_data/grid_search_ISO45001_standardised.RData", envir = e_45s)

r37  <- e_37$results[[which.min(e_37$bic_vec)]]
r45  <- e_45$results[[which.min(e_45$bic_vec)]]
r37s <- e_37s$results[[which.min(e_37s$bic_vec)]]
r45s <- e_45s$results[[which.min(e_45s$bic_vec)]]

# === SAVE RDS ===
saveRDS(r37,  "data/output_data/results_37001_N32.rds")
saveRDS(r45,  "data/output_data/results_45001_N32.rds")
saveRDS(r37s, "data/output_data/results_37001_N37_std.rds")
saveRDS(r45s, "data/output_data/results_45001_N37_std.rds")

# === HELPER: format number ===
fmt <- function(x, digits = 3) formatC(x, format = "f", digits = digits)

# =============================================================================
# TABLE 1 (tex)
# =============================================================================

t1 <- file("table1.tex", "w")

writeLines(c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Estimation results}",
  "\\label{tab:estimation}",
  "\\begin{tabular}{lcccc}",
  "\\hline\\hline",
  " & \\multicolumn{2}{c}{Non-standardised ($N=32$)} & \\multicolumn{2}{c}{Standardised ($N=37$)} \\\\",
  "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5}",
  " & ISO~37001 & ISO~45001 & ISO~37001 & ISO~45001 \\\\",
  "\\hline"
), t1)

# rho
writeLines(sprintf("$\\hat{\\rho}$ & %s & %s & %s & %s \\\\",
  fmt(r37$unpenalisedgmm$rho), fmt(r45$unpenalisedgmm$rho),
  fmt(r37s$unpenalisedgmm$rho), fmt(r45s$unpenalisedgmm$rho)), t1)

# beta
covnames <- c("GDP per capita", "Rule of Law", "Trade openness", "Unemployment")
for (i in 1:4) {
  writeLines(sprintf("$\\hat{\\beta}_{\\text{%s}}$ & %s & %s & %s & %s \\\\",
    covnames[i],
    fmt(r37$unpenalisedgmm$beta[i]), fmt(r45$unpenalisedgmm$beta[i]),
    fmt(r37s$unpenalisedgmm$beta[i]), fmt(r45s$unpenalisedgmm$beta[i])), t1)
}

# gamma
for (i in 1:4) {
  writeLines(sprintf("$\\hat{\\gamma}_{\\text{%s}}$ & %s & %s & %s & %s \\\\",
    covnames[i],
    fmt(r37$unpenalisedgmm$gamma[i]), fmt(r45$unpenalisedgmm$gamma[i]),
    fmt(r37s$unpenalisedgmm$gamma[i]), fmt(r45s$unpenalisedgmm$gamma[i])), t1)
}

writeLines("\\hline", t1)

# BIC
writeLines(sprintf("BIC & %s & %s & %s & %s \\\\",
  fmt(r37$BIC, 2), fmt(r45$BIC, 2),
  fmt(r37s$BIC, 2), fmt(r45s$BIC, 2)), t1)

# Edges
writeLines(sprintf("Edges in $\\hat{W}$ & %d & %d & %d & %d \\\\",
  sum(abs(r37$adaen$W) > 1e-5), sum(abs(r45$adaen$W) > 1e-5),
  sum(abs(r37s$adaen$W) > 1e-5), sum(abs(r45s$adaen$W) > 1e-5)), t1)

# N, T
writeLines("$N$ & 32 & 32 & 37 & 37 \\\\", t1)
writeLines("$T$ & 6 & 6 & 6 & 6 \\\\", t1)
writeLines("Standardised & No & No & Yes & Yes \\\\", t1)

writeLines(c(
  "\\hline\\hline",
  "\\end{tabular}",
  "\\begin{minipage}{0.95\\textwidth}",
  "\\vspace{4pt}",
  "\\footnotesize \\textit{Notes:} Estimates from the three-step adaptive elastic net GMM procedure of De~Paula et al.\\ (2024). $\\hat{\\rho}$ is the endogenous peer effect. $\\hat{\\beta}$ and $\\hat{\\gamma}$ are the own and contextual effects of covariates. Columns~1--2 use non-standardised data on the 32 original OECD members. Columns~3--4 standardise all variables to zero mean and unit variance and include all 37 OECD members. Penalisation parameters selected by BIC over a grid search.",
  "\\end{minipage}",
  "\\end{table}"
), t1)

close(t1)
cat("Saved table1.tex\n")

# =============================================================================
# TABLE 3 (tex) — one for each specification
# =============================================================================

write_table3 <- function(r1, r2, N, std_label, filename, caption_extra) {
  
  W1 <- unname(as.matrix(r1$unpenalisedgmm$W))
  W2 <- unname(as.matrix(r2$unpenalisedgmm$W))
  W1b <- (abs(W1) > 1e-5) + 0
  W2b <- (abs(W2) > 1e-5) + 0
  
  g1 <- graph_from_adjacency_matrix(W1b, mode = "directed")
  g2 <- graph_from_adjacency_matrix(W2b, mode = "directed")
  
  e1 <- sum(W1b); e2 <- sum(W2b)
  common <- sum(W1b * W2b)
  clust1 <- transitivity(g1, "global"); clust2 <- transitivity(g2, "global")
  recip1 <- sum(W1b * t(W1b)) / max(e1, 1)
  recip2 <- sum(W2b * t(W2b)) / max(e2, 1)
  
  f <- file(filename, "w")
  
  writeLines(c(
    "\\begin{table}[htbp]",
    "\\centering",
    sprintf("\\caption{ISO~37001 versus ISO~45001 networks%s}", caption_extra),
    sprintf("\\label{tab:networks_%s}", std_label),
    "\\begin{tabular}{lcc}",
    "\\hline\\hline",
    " & ISO~37001 & ISO~45001 \\\\",
    "\\hline",
    sprintf("Number of edges & %d & %d \\\\", e1, e2),
    sprintf("Edges in both networks & %d & %d \\\\", common, common),
    sprintf("Edges in ISO~37001 only & %d & \\\\", e1 - common),
    sprintf("Edges in ISO~45001 only & & %d \\\\", e2 - common),
    sprintf("Clustering & %s & %s \\\\", fmt(clust1), fmt(clust2)),
    sprintf("Reciprocated edges & %s\\%% & %s\\%% \\\\", fmt(100*recip1, 1), fmt(100*recip2, 1)),
    "Degree distribution & & \\\\",
    sprintf("\\quad Out-degree & %s (%s) & %s (%s) \\\\",
      fmt(mean(degree(g1, mode="out"))), fmt(sd(degree(g1, mode="out"))),
      fmt(mean(degree(g2, mode="out"))), fmt(sd(degree(g2, mode="out")))),
    sprintf("\\quad In-degree & %s (%s) & %s (%s) \\\\",
      fmt(mean(degree(g1, mode="in"))), fmt(sd(degree(g1, mode="in"))),
      fmt(mean(degree(g2, mode="in"))), fmt(sd(degree(g2, mode="in")))),
    "\\hline\\hline",
    "\\end{tabular}",
    "\\begin{minipage}{0.85\\textwidth}",
    "\\vspace{4pt}",
    sprintf("\\footnotesize \\textit{Notes:} Network statistics for the BIC-optimal estimated networks. $N=%d$, $T=6$. %s The clustering coefficient is the frequency of fully connected triplets over total triplets. Reciprocated edges is the share of directed edges whose reverse also exists. Degree distribution reports mean (standard deviation).", N, ifelse(std_label=="std", "All variables standardised to zero mean and unit variance.", "Non-standardised data.")),
    "\\end{minipage}",
    "\\end{table}"
  ), f)
  
  close(f)
  cat("Saved", filename, "\n")
}

write_table3(r37, r45, 32, "nonstd", "table3_nonstd.tex", " (non-standardised, $N=32$)")
write_table3(r37s, r45s, 37, "std", "table3_std.tex", " (standardised, $N=37$)")

cat("\nDone. Files: table1.tex, table3_nonstd.tex, table3_std.tex\n")