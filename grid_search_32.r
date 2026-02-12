# =============================================================================
# 02_grid_search.R
# Grid search su lambda per ISO 37001 e ISO 45001
# Tempo stimato: ~10 ore (2 ISO × 15 combinazioni × ~20 min ciascuna)
# Lancia e vai a dormire!
# =============================================================================

library(dplyr)
library(Matrix)
library(glmnet)
library(CVXR)
library(MASS)
library(optimx)
library(fixest)

source("/workspaces/thesis/functions.r")


# === LOAD DATA ===============================================================

data_37001 <- readRDS("data_balanced_37001.rds")
data_45001 <- readRDS("data_balanced_45001.rds")

# === FILTER N=32 (exclude countries added with case_match) ===================


old_countries <- c("Australia", "Austria", "Belgium", "Canada", "Chile", "Colombia", 
  "Costa Rica", "Denmark", "Estonia", "Finland", "France", "Germany", 
  "Greece", "Hungary", "Iceland", "Ireland", "Israel", "Italy", "Japan", 
  "Latvia", "Lithuania", "Luxembourg", "Mexico", "Netherlands", "New Zealand", 
  "Norway", "Poland", "Portugal", "Slovenia", "Spain", "Sweden", "Switzerland")

cmap_37 <- readRDS("country_map_37001.rds")
old_ids_37 <- cmap_37 %>% dplyr::filter(country %in% old_countries) %>% dplyr::pull(id)
data_37001 <- data_37001 %>%
  dplyr::filter(id %in% old_ids_37) %>%
  dplyr::mutate(id = as.integer(factor(id)))

cmap_45 <- readRDS("country_map_45001.rds")
old_ids_45 <- cmap_45 %>% dplyr::filter(country %in% old_countries) %>% dplyr::pull(id)
data_45001 <- data_45001 %>%
  dplyr::filter(id %in% old_ids_45) %>%
  dplyr::mutate(id = as.integer(factor(id)))

cat("ISO 37001 - N:", length(unique(data_37001$id)), 
    "T:", length(unique(data_37001$time)), "\n")
cat("ISO 45001 - N:", length(unique(data_45001$id)), 
    "T:", length(unique(data_45001$time)), "\n")


# === LAMBDA GRID =============================================================
# Basata sulla griglia di de Paula, Rasul & Souza (2025)
# Loro usano 84 combinazioni; noi prendiamo un sottoinsieme per stare in 10h
# Con ~20-30 min per run e 2 ISO, 15 combinazioni = ~10h

# Se una run ci mette 30 min: 24 × 30 min × 2 ISO = 24h -> troppo
# Se ci mette 20 min: 15 × 20 min × 2 ISO = 10h -> ok
# Riduciamo: 3 × 2 × 4 = 24 combinazioni per ISO
# Se 30 min/run: 24 × 30 × 2 = 24h -> troppo
# Allora: 3 × 1 × 5 = 15 combinazioni per ISO  
# 15 × 30 × 2 = 15h -> borderline
grid <- expand.grid(
  lambdaL1     = c(0.05, 0.1),
  lambdaL2     = c(0.01),
  lambdaL1Star = c(0.25, 0.375, 0.75)
)


cat("\nCombinazioni nella griglia:", nrow(grid), "\n")
cat("Tempo stimato (20 min/run):", nrow(grid) * 20 * 2 / 60, "ore\n")
cat("Tempo stimato (30 min/run):", nrow(grid) * 30 * 2 / 60, "ore\n")

# === FUNCTION: RUN GRID SEARCH ===============================================

run_grid_search <- function(data, iso_label, grid) {
  
  n_grid <- nrow(grid)
  results <- vector("list", n_grid)
  bic_vec <- rep(NA, n_grid)
  timings <- rep(NA, n_grid)
  
  # --- Prima run per calcolare Wder ---
  cat("\n", iso_label, "- Calcolo Wder (run iniziale)...\n")
  t0 <- Sys.time()
  
  first_run <- recoverNetwork(
    data = data,
    lambda = as.numeric(grid[1, ]),
    exoeffects = 1
  )
  
  Wder <- first_run$Wder
  results[[1]] <- first_run
  bic_vec[1] <- first_run$BIC
  timings[1] <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  
  cat(iso_label, "- Run 1/", n_grid, 
      " lambda=", as.numeric(grid[1,]),
      " BIC=", round(bic_vec[1], 4),
      " tempo=", round(timings[1], 1), "min\n")
  
  # --- Rimanenti run con Wder riusato ---
  for (g in 2:n_grid) {
    cat(iso_label, "- Run", g, "/", n_grid, 
        " lambda=", as.numeric(grid[g,]))
    
    t0 <- Sys.time()
    
    tryCatch({
      results[[g]] <- recoverNetwork(
        data = data,
        lambda = as.numeric(grid[g, ]),
        exoeffects = 1,
        Wder = Wder
      )
      bic_vec[g] <- results[[g]]$BIC
      timings[g] <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
      
      cat(" BIC=", round(bic_vec[g], 4),
          " tempo=", round(timings[g], 1), "min\n")
      
    }, error = function(e) {
      timings[g] <<- as.numeric(difftime(Sys.time(), t0, units = "mins"))
      cat(" ERRORE:", e$message, "\n")
    })
    
    # Salva risultati intermedi dopo ogni run (safety)
    save(results, bic_vec, timings, grid, Wder,
         file = paste0("grid_search_", iso_label, "_PARTIAL.RData"))
  }
  
  # Salva risultati finali
  save(results, bic_vec, timings, grid, Wder,
       file = paste0("grid_search_", iso_label, ".RData"))
  
  return(list(results = results, bic_vec = bic_vec, timings = timings, Wder = Wder))
}

# === RUN =====================================================================

cat("\n========================================\n")
cat("INIZIO GRID SEARCH\n")
cat("Ora:", format(Sys.time(), "%H:%M"), "\n")
cat("========================================\n")

# --- ISO 37001 ---
cat("\n>>> ISO 37001 <<<\n")
t_start_37 <- Sys.time()
gs_37001 <- run_grid_search(data_37001, "ISO37001", grid)
t_end_37 <- Sys.time()
cat("\nISO 37001 completato in", 
    round(as.numeric(difftime(t_end_37, t_start_37, units = "hours")), 1), "ore\n")

# --- ISO 45001 ---
cat("\n>>> ISO 45001 <<<\n")
t_start_45 <- Sys.time()
gs_45001 <- run_grid_search(data_45001, "ISO45001", grid)
t_end_45 <- Sys.time()
cat("\nISO 45001 completato in", 
    round(as.numeric(difftime(t_end_45, t_start_45, units = "hours")), 1), "ore\n")

# === SUMMARY =================================================================

cat("\n========================================\n")
cat("GRID SEARCH COMPLETATA\n")
cat("Ora:", format(Sys.time(), "%H:%M"), "\n")
cat("Tempo totale:", round(as.numeric(difftime(Sys.time(), t_start_37, units = "hours")), 1), "ore\n")
cat("========================================\n")

# --- ISO 37001 ---
valid_37 <- which(!is.na(gs_37001$bic_vec))
best_37 <- valid_37[which.min(gs_37001$bic_vec[valid_37])]
cat("\n=== ISO 37001 ===\n")
cat("  Miglior lambda:", as.numeric(grid[best_37,]), "\n")
cat("  BIC:", gs_37001$bic_vec[best_37], "\n")
cat("  rho:", gs_37001$results[[best_37]]$unpenalisedgmm$rho, "\n")
cat("  N edges:", sum(abs(gs_37001$results[[best_37]]$adaen$W) > 1e-5), "\n")

# --- ISO 45001 ---
valid_45 <- which(!is.na(gs_45001$bic_vec))
best_45 <- valid_45[which.min(gs_45001$bic_vec[valid_45])]
cat("\n=== ISO 45001 ===\n")
cat("  Miglior lambda:", as.numeric(grid[best_45,]), "\n")
cat("  BIC:", gs_45001$bic_vec[best_45], "\n")
cat("  rho:", gs_45001$results[[best_45]]$unpenalisedgmm$rho, "\n")
cat("  N edges:", sum(abs(gs_45001$results[[best_45]]$adaen$W) > 1e-5), "\n")

# --- Stability ---
cat("\n=== STABILITY ANALYSIS ===\n")
for (iso in c("37001", "45001")) {
  gs <- get(paste0("gs_", iso))
  valid <- which(!is.na(gs$bic_vec))
  all_W <- lapply(gs$results[valid], function(r) (abs(r$adaen$W) > 1e-5) + 0)
  stability <- Reduce("+", all_W) / length(all_W)
  
  cat(paste0("\nISO ", iso, ":\n"))
  cat("  Edge robusti (>80%):", sum(stability > 0.8), "\n")
  cat("  Edge robusti (>60%):", sum(stability > 0.6), "\n")
  cat("  Edge robusti (>40%):", sum(stability > 0.4), "\n")
}

# === SALVA RISULTATI FINALI BEST =============================================

W_best_37001 <- gs_37001$results[[best_37]]$unpenalisedgmm$W
W_best_45001 <- gs_45001$results[[best_45]]$unpenalisedgmm$W
rho_best_37001 <- gs_37001$results[[best_37]]$unpenalisedgmm$rho
rho_best_45001 <- gs_45001$results[[best_45]]$unpenalisedgmm$rho

saveRDS(W_best_37001, file = "W_best_37001.rds")
saveRDS(W_best_45001, file = "W_best_45001.rds")
saveRDS(rho_best_37001, file = "rho_best_37001.rds")
saveRDS(rho_best_45001, file = "rho_best_45001.rds")

# Salva anche la tabella BIC completa
bic_table <- grid %>%
  mutate(
    BIC_37001 = gs_37001$bic_vec,
    BIC_45001 = gs_45001$bic_vec,
    time_37001 = gs_37001$timings,
    time_45001 = gs_45001$timings
  ) %>%
  arrange(BIC_37001)

saveRDS(bic_table, file = "bic_table.rds")

cat("\n=== TABELLA BIC (ordinata per ISO 37001) ===\n")
print(bic_table)

cat("\n=== FILE SALVATI ===\n")
cat("  W_best_37001.rds\n")
cat("  W_best_45001.rds\n")
cat("  rho_best_37001.rds\n")
cat("  rho_best_45001.rds\n")
cat("  bic_table.rds\n")
cat("  grid_search_ISO37001.RData  (tutti i risultati)\n")
cat("  grid_search_ISO45001.RData  (tutti i risultati)\n")

###
load("grid_search_ISO37001.RData")
rho_37_all <- sapply(results, function(r) r$unpenalisedgmm$rho)

load("grid_search_ISO45001.RData")
rho_45_all <- sapply(results, function(r) r$unpenalisedgmm$rho)

comparison <- grid %>%
  mutate(
    rho_37001 = rho_37_all,
    rho_45001 = rho_45_all,
    rho_37_greater = rho_37_all > rho_45_all
  )

print(comparison)

cat("\nrho_37001 > rho_45001 in tutte le specificazioni?", 
    all(comparison$rho_37_greater), "\n")