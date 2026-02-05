###########################################################
## SETUP GEOGRAFICO AUTOMATICO (UE-27)                   ##
###########################################################

# Carica le librerie necessarie
library(rnaturalearth) # Per i dati geografici
library(sf)            # Per manipolare dati spaziali
library(spdep)         # Per creare matrici di adiacenza
library(dplyr)

# 1. SCARICA I PAESI DEL MONDO
world <- ne_countries(scale = 50, returnclass = "sf")

# 2. DEFINISCI L'UNIONE EUROPEA (EU-27)
# Usiamo i codici ISO per essere sicuri al 100%
eu_iso <- c("AUT", "BEL", "BGR", "CZE", "DEU", "DNK", "EST", "ESP", 
           "FIN", "FRA", "GRC", "HRV", "HUN", "ITA", "LTU", "LUX", 
           "LVA", "NLD", "POL", "PRT", "ROU", "SVK", "SVN", "SWE")

eu_map <- world %>% 
  filter(iso_a3 %in% eu_iso) %>%
  dplyr::select(iso_a3, admin) %>%
  arrange(iso_a3)

# 3. CREA LA MATRICE DI ADIACENZA (CONTIGUITÀ)
# poly2nb trova chi condivide un confine
nb <- poly2nb(eu_map, queen = TRUE) 

# Nota: Isole come Irlanda, Malta e Cipro risulteranno isolate (0 vicini)


# 4. CONVERTI IN MATRICE REALE
W_geo <- nb2mat(nb, style = "W", zero.policy = TRUE)

# Visualizziamo i nomi per sicurezza
rownames(W_geo) <- eu_map$iso_a3
colnames(W_geo) <- eu_map$iso_a3

print("Matrice Geografica UE-27 generata!")
print(W_geo[1:10, 1:10]) # Anteprima dei primi 10










################################################################################
## SIMULAZIONE NETWORK EU - CODICE COMPLETO
################################################################################

# PREREQUISITI:
# - Hai caricato il package recovernetwork
# - Hai la matrice W_geo (23x23) già in memoria

library(MASS)
library(dplyr)
library(tibble)

################################################################################
## 1. FUNZIONE GENDATA MODIFICATA (copia-incolla, non modifica l'originale)
################################################################################

gendata_eu <- function(N, T, rho, beta, gamma, W_matrix, seed, FEsetting = "standard", varcov = "orthogonal") {
  
  # SET SEED
  set.seed(seed)
  
  K <- 1  # numero di covariate
  
  # Create fixed effects
  if (FEsetting == "standard") {
    FE <- rnorm(N) + 1
    TE <- rnorm(N) + 1
  } else if (FEsetting == "zero") {
    FE <- 0 * rnorm(N) + 1
    TE <- 0 * rnorm(N) + 1
  }
  
  # Define variance-covariance matrix
  if (varcov == "orthogonal") {
    Sigma <- diag(N)
  } else if (varcov == "q3") {
    q <- 0.3
    Sigma <- (1-q) * diag(N) + q
  } else if (varcov == "q5") {
    q <- 0.5
    Sigma <- (1-q) * diag(N) + q
  }
  
  # Create reduced-form matrix
  Pi0inv <- solve(diag(N) - rho * W_matrix)
  
  # Generate data for each time period
  data <- tibble(y = double(), x = double(), id = integer(), time = integer())
  
  for (t in 1:T) {
    
    idt <- 1:N
    timet <- rep(t, N)
    
    et <- MASS::mvrnorm(n = 1, mu = rep(0, N), Sigma)
    xt <- rnorm(N * K)
    dim(xt) <- c(N, K)
    
    yt <- Pi0inv %*% xt %*% rep(beta, K) + 
          Pi0inv %*% W_matrix %*% xt %*% rep(gamma, K) + 
          Pi0inv %*% FE + 
          Pi0inv %*% TE + 
          Pi0inv %*% et
    
    xt <- as_tibble(as.data.frame(xt))
    names(xt) <- paste0('x', 1:ncol(xt))
    
    datat <- bind_cols(tibble(id = idt, time = timet), 
                       tibble(y = as.vector(yt)), 
                       xt)
    data <- rbind(data, datat)
  }
  
  return(data)
}

################################################################################
## 2. PARAMETRI SIMULAZIONE
################################################################################

# La tua matrice W (già in memoria)
# W_geo <- ... (la matrice 23x23 che hai creato)

# Parametri
N <- 23          # numero paesi EU
T <- 8         # periodi temporali (usa 50, non 8 per ora!)
rho <- 0.3       # coefficiente spillover
beta <- 0.4      # coefficiente direct effect
gamma <- 0.5     # coefficiente indirect effect (WX)
seed <- 1        # seed per riproducibilità

################################################################################
## 3. UNA SINGOLA SIMULAZIONE
################################################################################

# Genera dati
cat("Generando dati simulati...\n")
data_sim <- gendata_eu(
  N = N,
  T = T,
  rho = rho,
  beta = beta,
  gamma = gamma,
  W_matrix = W_geo,
  seed = seed,
  FEsetting = "standard",  # "standard" o "zero"
  varcov = "orthogonal"    # "orthogonal", "q3", "q5"
)

# Visualizza dati
cat("Dati generati:\n")
print(head(data_sim, 20))
cat("\nDimensioni:", nrow(data_sim), "righe (", N, "paesi x", T, "periodi)\n")

# Stima il network
cat("\nStimando il network...\n")
lambda <- c(0.01, 0.01, 0.10)  # parametri regolarizzazione

risultati <- recoverNetwork(
  data = data_sim,
  lambda = lambda,
  docv = 0,      # 0 = no cross-validation
  eta = 0.01     # learning rate
)

################################################################################
## 4. RISULTATI
################################################################################

cat("\n=== RISULTATI ===\n")
cat("Parametri VERI:\n")
cat("  rho (spillover):", rho, "\n")
cat("  beta (direct):", beta, "\n")
cat("  gamma (indirect):", gamma, "\n")

cat("\nParametri STIMATI:\n")
cat("  rho_hat:", risultati$rhohat, "\n")
cat("  beta_hat:", risultati$betahat, "\n")
cat("  gamma_hat:", risultati$gammahat, "\n")

cat("\nErrori:\n")
cat("  Errore rho:", risultati$rhohat - rho, "\n")
cat("  Errore beta:", risultati$betahat - beta, "\n")
cat("  Errore gamma:", risultati$gammahat - gamma, "\n")

# Network stimato vs vero
cat("\nNetwork W (primi 5x5):\n")
cat("VERO:\n")
print(round(W_geo[1:5, 1:5], 3))
cat("\nSTIMATO:\n")
print(round(risultati$What[1:5, 1:5], 3))

################################################################################
## RISULTATI CORRETTI
################################################################################

cat("\n=== RISULTATI ===\n\n")

cat("Parametri VERI:\n")
cat("  rho (spillover):", rho, "\n")
cat("  beta (direct):", beta, "\n")
cat("  gamma (indirect):", gamma, "\n\n")

# METODO 1: Elastic Net
cat("--- Elastic Net ---\n")
cat("  rho_hat:", risultati$en$rho, "\n")
cat("  beta_hat:", risultati$en$beta, "\n")
cat("  gamma_hat:", risultati$en$gamma, "\n")
cat("  Errore rho:", risultati$en$rho - rho, "\n\n")

# METODO 2: Adaptive Elastic Net
cat("--- Adaptive Elastic Net ---\n")
cat("  rho_hat:", risultati$adaen$rho, "\n")
cat("  beta_hat:", risultati$adaen$beta, "\n")
cat("  gamma_hat:", risultati$adaen$gamma, "\n")
cat("  Errore rho:", risultati$adaen$rho - rho, "\n\n")

# METODO 3: Unpenalised GMM 
cat("--- Unpenalised GMM (migliore) ---\n")
cat("  rho_hat:", risultati$unpenalisedgmm$rho, "\n")
cat("  beta_hat:", risultati$unpenalisedgmm$beta, "\n")
cat("  gamma_hat:", risultati$unpenalisedgmm$gamma, "\n")
cat("  Errore rho:", risultati$unpenalisedgmm$rho - rho, "\n\n")

# Network stimato vs vero
cat("Network W (primi 5x5):\n\n")
cat("VERO:\n")
print(round(W_geo[1:5, 1:5], 3))
cat("\nSTIMATO (Adaptive Elastic Net):\n")
print(round(risultati$adaen$W[1:5, 1:5], 3))

# BIC
cat("\nBIC:", risultati$BIC, "\n")

# Debug: vediamo cosa contiene risultati
cat("\n=== STRUTTURA RISULTATI ===\n")
str(risultati)

cat("\n=== NOMI ELEMENTI ===\n")
names(risultati)

################################################################################
## VISUALIZZA MATRICE W COMPLETA
################################################################################

# Matrice VERA
cat("\n=== MATRICE W VERA (23x23) ===\n")
print(round(W_geo, 3))

# Matrice STIMATA (Adaptive Elastic Net)
cat("\n=== MATRICE W STIMATA (23x23) ===\n")
print(round(risultati$adaen$W, 3))

# Oppure salva in formato più leggibile
cat("\n=== CONFRONTO FIANCO A FIANCO ===\n")
cat("\nNumero di link non-zero:\n")
cat("  VERO:", sum(W_geo > 0), "link\n")
cat("  STIMATO:", sum(risultati$adaen$W > 0.001), "link (soglia 0.001)\n")

# Matrice differenza
cat("\n=== MATRICE DIFFERENZA (Stimato - Vero) ===\n")
diff_W <- risultati$adaen$W - W_geo
print(round(diff_W, 3))

# Statistiche
cat("\n=== STATISTICHE ===\n")
cat("Max errore:", max(abs(diff_W)), "\n")
cat("RMSE:", sqrt(mean(diff_W^2)), "\n")
cat("Correlazione:", cor(as.vector(W_geo), as.vector(risultati$adaen$W)), "\n")

# Visualizza solo i link più forti
cat("\n=== LINK FORTI VERI (>0.2) ===\n")
strong_links_true <- which(W_geo > 0.2, arr.ind = TRUE)
rownames(strong_links_true) <- NULL
print(data.frame(
  From = strong_links_true[,1],
  To = strong_links_true[,2],
  Weight_True = W_geo[strong_links_true],
  Weight_Estimated = risultati$adaen$W[strong_links_true]
))


###---------------------------


###########################################################
## SETUP GEOGRAFICO AUTOMATICO OECD                  ##
###########################################################

# Carica le librerie necessarie
library(rnaturalearth) # Per i dati geografici
library(sf)            # Per manipolare dati spaziali
library(spdep)         # Per creare matrici di adiacenza
library(dplyr)

# 1. SCARICA I PAESI DEL MONDO
world <- ne_countries(scale = 50, returnclass = "sf")

# 2. DEFINISCI L'UNIONE EUROPEA (EU-27)
# Usiamo i codici ISO per essere sicuri al 100%
eu_iso <- c("AUS", "AUT", "BEL", "CAN", "CHL", "COL", "CRI", "CZE", "DNK", 
  "EST", "FIN", "FRA", "DEU", "GRC", "HUN", "ISL", "IRL", "ISR", 
  "ITA", "JPN", "KOR", "LVA", "LTU", "LUX", "MEX", "NLD", "NZL", 
  "NOR", "POL", "PRT", "SVK", "SVN", "ESP", "SWE", "CHE", "TUR", 
  "GBR", "USA")

eu_map <- world %>% 
  filter(iso_a3 %in% eu_iso) %>%
  dplyr::select(iso_a3, admin) %>%
  arrange(iso_a3)

# 3. CREA LA MATRICE DI ADIACENZA (CONTIGUITÀ)
# poly2nb trova chi condivide un confine
nb <- poly2nb(eu_map, queen = TRUE) 

# Nota: Isole come Irlanda, Malta e Cipro risulteranno isolate (0 vicini)


# 4. CONVERTI IN MATRICE REALE
W_geo <- nb2mat(nb, style = "W", zero.policy = TRUE)

# Visualizziamo i nomi per sicurezza
rownames(W_geo) <- eu_map$iso_a3
colnames(W_geo) <- eu_map$iso_a3

print("Matrice Geografica UE-27 generata!")
print(W_geo[1:10, 1:10]) # Anteprima dei primi 10










################################################################################
## SIMULAZIONE NETWORK EU - CODICE COMPLETO
################################################################################

# PREREQUISITI:
# - Hai caricato il package recovernetwork
# - Hai la matrice W_geo (23x23) già in memoria

library(MASS)
library(dplyr)
library(tibble)

################################################################################
## 1. FUNZIONE GENDATA MODIFICATA (copia-incolla, non modifica l'originale)
################################################################################

gendata_eu <- function(N, T, rho, beta, gamma, W_matrix, seed, FEsetting = "standard", varcov = "orthogonal") {
  
  # SET SEED
  set.seed(seed)
  
  K <- 1  # numero di covariate
  
  # Create fixed effects
  if (FEsetting == "standard") {
    FE <- rnorm(N) + 1
    TE <- rnorm(N) + 1
  } else if (FEsetting == "zero") {
    FE <- 0 * rnorm(N) + 1
    TE <- 0 * rnorm(N) + 1
  }
  
  # Define variance-covariance matrix
  if (varcov == "orthogonal") {
    Sigma <- diag(N)
  } else if (varcov == "q3") {
    q <- 0.3
    Sigma <- (1-q) * diag(N) + q
  } else if (varcov == "q5") {
    q <- 0.5
    Sigma <- (1-q) * diag(N) + q
  }
  
  # Create reduced-form matrix
  Pi0inv <- solve(diag(N) - rho * W_matrix)
  
  # Generate data for each time period
  data <- tibble(y = double(), x = double(), id = integer(), time = integer())
  
  for (t in 1:T) {
    
    idt <- 1:N
    timet <- rep(t, N)
    
    et <- MASS::mvrnorm(n = 1, mu = rep(0, N), Sigma)
    xt <- rnorm(N * K)
    dim(xt) <- c(N, K)
    
    yt <- Pi0inv %*% xt %*% rep(beta, K) + 
          Pi0inv %*% W_matrix %*% xt %*% rep(gamma, K) + 
          Pi0inv %*% FE + 
          Pi0inv %*% TE + 
          Pi0inv %*% et
    
    xt <- as_tibble(as.data.frame(xt))
    names(xt) <- paste0('x', 1:ncol(xt))
    
    datat <- bind_cols(tibble(id = idt, time = timet), 
                       tibble(y = as.vector(yt)), 
                       xt)
    data <- rbind(data, datat)
  }
  
  return(data)
}

################################################################################
## 2. PARAMETRI SIMULAZIONE
################################################################################

# La tua matrice W (già in memoria)
# W_geo <- ... (la matrice 23x23 che hai creato)

# Parametri
N <- 36          # numero paesi EU
T <- 8         # periodi temporali (usa 50, non 8 per ora!)
rho <- 0.3       # coefficiente spillover
beta <- 0.4      # coefficiente direct effect
gamma <- 0.5     # coefficiente indirect effect (WX)
seed <- 1        # seed per riproducibilità

################################################################################
## 3. UNA SINGOLA SIMULAZIONE
################################################################################

# Genera dati
cat("Generando dati simulati...\n")
data_sim <- gendata_eu(
  N = N,
  T = T,
  rho = rho,
  beta = beta,
  gamma = gamma,
  W_matrix = W_geo,
  seed = seed,
  FEsetting = "standard",  # "standard" o "zero"
  varcov = "orthogonal"    # "orthogonal", "q3", "q5"
)

# Visualizza dati
cat("Dati generati:\n")
print(head(data_sim, 20))
cat("\nDimensioni:", nrow(data_sim), "righe (", N, "paesi x", T, "periodi)\n")

# Stima il network
cat("\nStimando il network...\n")
lambda <- c(0.01, 0.01, 0.10)  # parametri regolarizzazione

risultati <- recoverNetwork(
  data = data_sim,
  lambda = lambda,
  docv = 0,      # 0 = no cross-validation
  eta = 0.01     # learning rate
)

################################################################################
## 4. RISULTATI
################################################################################

cat("\n=== RISULTATI ===\n")
cat("Parametri VERI:\n")
cat("  rho (spillover):", rho, "\n")
cat("  beta (direct):", beta, "\n")
cat("  gamma (indirect):", gamma, "\n")

cat("\nParametri STIMATI:\n")
cat("  rho_hat:", risultati$rhohat, "\n")
cat("  beta_hat:", risultati$betahat, "\n")
cat("  gamma_hat:", risultati$gammahat, "\n")

cat("\nErrori:\n")
cat("  Errore rho:", risultati$rhohat - rho, "\n")
cat("  Errore beta:", risultati$betahat - beta, "\n")
cat("  Errore gamma:", risultati$gammahat - gamma, "\n")

# Network stimato vs vero
cat("\nNetwork W (primi 5x5):\n")
cat("VERO:\n")
print(round(W_geo[1:5, 1:5], 3))
cat("\nSTIMATO:\n")
print(round(risultati$What[1:5, 1:5], 3))

################################################################################
## RISULTATI CORRETTI
################################################################################

cat("\n=== RISULTATI ===\n\n")

cat("Parametri VERI:\n")
cat("  rho (spillover):", rho, "\n")
cat("  beta (direct):", beta, "\n")
cat("  gamma (indirect):", gamma, "\n\n")

# METODO 1: Elastic Net
cat("--- Elastic Net ---\n")
cat("  rho_hat:", risultati$en$rho, "\n")
cat("  beta_hat:", risultati$en$beta, "\n")
cat("  gamma_hat:", risultati$en$gamma, "\n")
cat("  Errore rho:", risultati$en$rho - rho, "\n\n")

# METODO 2: Adaptive Elastic Net
cat("--- Adaptive Elastic Net ---\n")
cat("  rho_hat:", risultati$adaen$rho, "\n")
cat("  beta_hat:", risultati$adaen$beta, "\n")
cat("  gamma_hat:", risultati$adaen$gamma, "\n")
cat("  Errore rho:", risultati$adaen$rho - rho, "\n\n")

# METODO 3: Unpenalised GMM 
cat("--- Unpenalised GMM (migliore) ---\n")
cat("  rho_hat:", risultati$unpenalisedgmm$rho, "\n")
cat("  beta_hat:", risultati$unpenalisedgmm$beta, "\n")
cat("  gamma_hat:", risultati$unpenalisedgmm$gamma, "\n")
cat("  Errore rho:", risultati$unpenalisedgmm$rho - rho, "\n\n")

# Network stimato vs vero
cat("Network W (primi 5x5):\n\n")
cat("VERO:\n")
print(round(W_geo[1:5, 1:5], 3))
cat("\nSTIMATO (Adaptive Elastic Net):\n")
print(round(risultati$adaen$W[1:5, 1:5], 3))

# BIC
cat("\nBIC:", risultati$BIC, "\n")

# Debug: vediamo cosa contiene risultati
cat("\n=== STRUTTURA RISULTATI ===\n")
str(risultati)

cat("\n=== NOMI ELEMENTI ===\n")
names(risultati)

################################################################################
## VISUALIZZA MATRICE W COMPLETA
################################################################################

# Matrice VERA
cat("\n=== MATRICE W VERA (23x23) ===\n")
print(round(W_geo, 3))

# Matrice STIMATA (Adaptive Elastic Net)
cat("\n=== MATRICE W STIMATA (23x23) ===\n")
print(round(risultati$adaen$W, 3))

# Oppure salva in formato più leggibile
cat("\n=== CONFRONTO FIANCO A FIANCO ===\n")
cat("\nNumero di link non-zero:\n")
cat("  VERO:", sum(W_geo > 0), "link\n")
cat("  STIMATO:", sum(risultati$adaen$W > 0.001), "link (soglia 0.001)\n")

# Matrice differenza
cat("\n=== MATRICE DIFFERENZA (Stimato - Vero) ===\n")
diff_W <- risultati$adaen$W - W_geo
print(round(diff_W, 3))

# Statistiche
cat("\n=== STATISTICHE ===\n")
cat("Max errore:", max(abs(diff_W)), "\n")
cat("RMSE:", sqrt(mean(diff_W^2)), "\n")
cat("Correlazione:", cor(as.vector(W_geo), as.vector(risultati$adaen$W)), "\n")

# Visualizza solo i link più forti
cat("\n=== LINK FORTI VERI (>0.2) ===\n")
strong_links_true <- which(W_geo > 0.2, arr.ind = TRUE)
rownames(strong_links_true) <- NULL
print(data.frame(
  From = strong_links_true[,1],
  To = strong_links_true[,2],
  Weight_True = W_geo[strong_links_true],
  Weight_Estimated = risultati$adaen$W[strong_links_true]
))