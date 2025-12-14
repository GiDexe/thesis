# Load required libraries
library(readxl)
library(tidyverse)
library(dplyr)
library(tidyr)
library(janitor)
library(purrr)
library(WDI)

source("/workspaces/thesis/functions.r")

# Set your data directory
data_dir <- "/workspaces/thesis/data"

# Get all Excel files
all_excel_files <- list.files(data_dir,
  pattern = "\\.xlsx$|\\.xls$|\\.xlsm$",
  full.names = TRUE
)
# Restricting to ISO 37001 (computational constraints)

# Get all Excel files about ISO 37001
excel_files <- all_excel_files[grepl(
  "ISO Survey 20",
  basename(all_excel_files)
)]

# Function to extract all sheets except the first from one file
extract_sheets <- function(file_path) {
  # Extract the year from the file path (4-digit number between 2018 and 2024)
  year <- str_extract(basename(file_path), "20(1[8-9]|2[0-4])")
  
  # Get all sheet names
  sheet_names <- excel_sheets(file_path)

  # Skip the first sheet
  sheets_to_read <- sheet_names[-1]

  # Read each sheet into a list and add the year column
  data_list <- map(sheets_to_read, ~ {
    read_excel(file_path, sheet = .x) %>%
      mutate(Year = year)
  })

  # Name the list elements
  names(data_list) <- sheets_to_read

  return(data_list)
}


# Extract sheets from all files
all_data <- map(excel_files, extract_sheets)

# Nested structure 
g <- all_data[[1]]


g[1]


corrupt <- map(all_data, ~ keep(.x, ~ any(grepl("37001", colnames(.)))))

#countries <- c("Italy", "Brazil", "South Korea", "Egypt", "Kenya", 
# "Norway", "Thailand", "Argentina", "Morocco", "New Zealand")
corrupt <- map(all_data, ~ keep(.x, ~ any(grepl("37001", colnames(.)))))



# Adjusting tibble structure for all tibbles
adjusted_data <- map(corrupt, ~ map(.x, ~ {
  .x %>%
    setNames(ifelse(names(.) == "Year", "Year", as.character(slice(., 2)))) %>%
    slice(-2) %>%
    mutate(across(-1, ~ replace_na(as.numeric(.), 0))) %>%
    mutate(Tot = rowSums(select(., where(is.numeric) & !last_col()), na.rm = TRUE), .before = Year)
}))

my_panel <- adjusted_data %>%
  imap_dfr(~ bind_rows(.x, .id = "sheet") %>% 
             mutate(country = .y, .before = 1))

my_panel <- my_panel %>%
  select(country, `Land/Sector`, Tot, Year) %>%
  filter(`Land/Sector` != "sector number") %>%
  filter(`Land/Sector` != "Sector number") %>%
  clean_names() %>%
  rename(code = country) %>%
  rename(country = land_sector)

my_panel <- my_panel %>%
  group_by(country) %>%
  arrange(year) %>%  # ascending
  mutate(delta = tot - lag(tot, default = 0))

covariates <- WDI(
  country = "all",
  indicator = c(
    "NY.GDP.PCAP.PP.KD",      # GDP per capita (constant PPP)
    "CC.EST",                 # Control of Corruption
    "RQ.EST",                 # Regulatory Quality
    "SL.UEM.TOTL.ZS"          # Unemployment rate
  ),
  start = 2018,
  end = 2024,
  extra = TRUE
) %>%
  select(country, year, NY.GDP.PCAP.PP.KD, CC.EST, RQ.EST, SL.UEM.TOTL.ZS) %>%
  rename(
    gdp_pc = NY.GDP.PCAP.PP.KD,
    corruption_control = CC.EST,
    regulatory_quality = RQ.EST,
    unemployment = SL.UEM.TOTL.ZS
  ) %>%
  filter(!is.na(gdp_pc))


####


#---
library(dplyr)
library(tidyr)



# Identifica l'insieme di paesi presenti nell'anno di riferimento (2020)
countries_2020 <- my_panel %>%
  filter(year == 2018) %>%
  pull(country) %>%
  unique()

# =========================================================================
# 2. PREPARAZIONE E RECTANGULARIZZAZIONE (Gestione dei dati mancanti)
# =========================================================================

my_panel <- my_panel %>%
  
  # A. Filtra il dataset: mantiene SOLO i paesi presenti nel 2020.
  filter(country %in% countries_2020) %>%
  
  # B. Annulla qualsiasi raggruppamento precedente, essenziale per complete().
  ungroup() %>%
  
  # C. Rectangularizzazione: riempie le combinazioni country-year mancanti 
  #    all'interno del range di anni coperto dal dataset filtrato.
  complete(
    country, 
    year,
    fill = list(tot = NA) # I valori mancanti nell'originale diventano NA
  ) %>%
  
  # =========================================================================
  # 3. CALCOLO DEL DELTA E IMPUTAZIONE (Delta = 0 se dato mancante)
  # =========================================================================
  
  # Raggruppa per paese per eseguire i calcoli lag.
  group_by(country) %>%
  
  # Ordina per anno, fondamentale per la funzione lag().
  arrange(year) %>%
  
  # Imputazione con LOCF (Last Observation Carried Forward):
  # Se tot_t è NA, viene riempito con tot_{t-1}.
  # Questo assicura che il calcolo del delta sia: tot_{t-1} - tot_{t-1} = 0.
  fill(tot, .direction = "down") %>%
  
  # Imputazione dei valori NA all'inizio della serie:
  # Se un paese non ha mai avuto un'osservazione (o inizia con NA), tot viene impostato a 0.
  mutate(tot = replace_na(tot, 0)) %>%
  
  # Calcola il delta. 
  # Il 'default = 0' gestisce la prima osservazione di ciascun paese (delta = tot_1 - 0).
  mutate(delta = tot - lag(tot, default = 0)) %>%
  
  # Rimuove il raggruppamento alla fine del pipeline
  ungroup()

# Merge and format
panel_data <- my_panel %>%
  left_join(covariates, by = c("country", "year")) %>%
  filter(!is.na(gdp_pc)) %>%
  mutate(id = as.numeric(factor(country, levels = unique(country)))) %>%
  rename(
    time = year,
    y = delta,
    x1 = gdp_pc,
    x2 = corruption_control,
    x3 = regulatory_quality,
    x4 = unemployment
  )

  # La prima parte del tuo codice (senza creare l'ID)
panel_data <- my_panel %>%
  left_join(covariates, by = c("country", "year")) %>%
  filter(!is.na(gdp_pc)) %>%
  # Rimuovi: mutate(id = as.numeric(factor(country, levels = unique(country))))
  rename(
    time = year,
    y = delta,
    x1 = gdp_pc,
    x2 = corruption_control,
    x3 = regulatory_quality,
    x4 = unemployment
  )

# Blocco 1: Creazione dell'ID, pulita e robusta (dopo il rename)
panel_data <- panel_data %>%
  # *** PASSO CRUCIALE: Rimuovere qualsiasi raggruppamento ereditato ***
  ungroup() %>% 
  
  # Assegna l'ID unico a ciascun paese
  mutate(id = group_indices(., country))
  # Ora group_indices vede il dataframe intero e assegna ID 1 a Afghanistan, 2 a Albania, ecc.

# Blocco 2: Mantenere un panel bilanciato (il resto del codice non cambia)
panel_balance <- panel_data %>%
  group_by(id) %>%
  summarize(country = first(country), n_obs = n(), .groups = "drop")

complete_ids <- panel_balance %>%
  filter(n_obs == max(n_obs)) %>%
  pull(id)

panel_data <- panel_data %>%
  filter(id %in% complete_ids) %>%
  select(-"code")

# lambda unica per test
lambda <- c(0.05)# Install required packages if not already installed
required_packages <- c("Matrix", "glmnet", "CVXR", "MASS", "optimx", "fixest")

# Check and install missing packages
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) {
  install.packages(new_packages)
}

# Load the libraries
library(Matrix)       # per matrici sparse e bdiag
library(glmnet)       # elastic net
library(CVXR)         # se usi adaestmethod="cvxr"
library(MASS)         # per ginv
library(optimx)       # per outerOptim
library(fixest)       # regressioni panel / IV
# chiamata super-light preceduta da pulizia
# Rimuovi le colonne non numeriche che non sono necessarie per il modello.
# Mantieni solo le variabili di interesse + id e time.
panel_data_clean <- panel_data %>%
    # Usa dplyr::select per evitare il conflitto con MASS::select
    dplyr::select(-country, -tot)
# Colonne rilevanti per il modello
MODEL_COLS <- c("y", "x1", "x2", "x3", "x4", "id", "time")

panel_data_final <- panel_data_clean %>%
  # *** PASSO CRUCIALE: rimuove tutte le 54 righe con NA menzionate nell'errore ***
  drop_na(all_of(MODEL_COLS)) # drop_na richiede il caricamento di tidyr

# Verifica la dimensione. Dovrebbe essere 310 (364 - 54)
# dim(panel_data_final) 

library(dplyr)
library(tidyr)

# Assicuriamo che l'elenco delle colonne del modello sia definito
MODEL_COLS <- c("y", "x1", "x2", "x3", "x4", "id", "time")

# =========================================================================
# PULIZIA FINALE DEL DATASET
# =========================================================================

panel_data_final <- panel_data %>%
  
  # 1. Selezione Rigorosa: Mantieni solo le colonne necessarie al modello (numeriche).
  # Questo risolve il problema 'country' e rimuove 'tot'.
  dplyr::select(all_of(MODEL_COLS)) %>%
  
  # 2. Rimozione degli NA: Elimina tutte le righe che contengono NA in una 
  #    qualsiasi delle variabili del modello (y, x1...x4, id, time).
  drop_na(all_of(MODEL_COLS)) %>%
  
  # 3. Pulizia per Effetti Fissi (Fase Cruciale): Rimuovi i gruppi con n=1.
  
  #    a) Filtra i Paesi (id): Rimuove tutti i paesi che hanno solo una riga (anno) di dati.
  group_by(id) %>%
  filter(n() > 1) %>%
  ungroup() %>%
  
  #    b) Filtra gli Anni (time): Rimuove tutti gli anni che hanno solo un paese di dati.
  #    Questo è cruciale se hai rimosso molti paesi con il filtro NA.
  group_by(time) %>%
  filter(n() > 1) %>%
  ungroup()

# Verifica dei gruppi id
count_id <- panel_data_final %>% group_by(id) %>% summarise(n = n())
# Controlla: max(count_id$n)
# Controlla: min(count_id$n) # Dovrebbe essere > 1

# Verifica dei gruppi time
count_time <- panel_data_final %>% group_by(time) %>% summarise(n = n())

# Colonne rilevanti
MODEL_COLS_NUMERIC <- c("y", "x1", "x2", "x3", "x4", "id", "time")

# Passaggio di pulizia estrema
panel_data_final <- panel_data_final %>%
  # Rende esplicitamente tutte le colonne numeriche (per sicurezza)
  mutate(across(all_of(MODEL_COLS_NUMERIC), as.numeric)) %>%
  
  # Imputazione zero per tutti i valori NA rimasti (il modo più aggressivo)
  mutate(across(everything(), ~replace_na(., 0)))
  
# Nota: La riga qui sopra potrebbe distorcere i dati se ci fossero NA reali, 
# ma è l'unico modo per soddisfare la richiesta di glmnet quando drop_na fallisce.
# Controlla: min(count_time$n) # Dovrebbe essere > 1
# =========================================================================
# ESECUZIONE DEL MODELLO
# La costante 'c' (es. 1) assicura che log(y) sia ben definito anche per y=0.
# Se y = 0, log(y+1) = log(1) = 0.
panel_data_final <- panel_data_final %>%
  # *** TRASFORMAZIONE IHS: Gestisce Zeri e Negativi ***
  # Nota: asinh(y) è l'implementazione in R dell'IHS.
  mutate(y = asinh(y)) 
  # Se preferisci la formula espansa: 
  # mutate(y = log(y + sqrt(y^2 + 1)))
view(panel_data_final)
write.csv(panel_data_final, "/workspaces/thesis/data/panel_data_final.csv", row.names = FALSE)
# Esegui la regressione con il nuovo 'y' trasformato
res <- recoverNetwork(
      data = panel_data_final, # <-- USA IL DATA FRAME CON y TRASFORMATO
      lambda = lambda,
      exoeffects = 1,
      docv = 0,
      timeweights = 1,
      eta = 0.05,
      adaestmethod = "L-BFGS-B", # Manteniamo L-BFGS-B per stabilità
      dopostols = 0,
      Wfixed = -1
    )


# output essenziale
res$unpenalisedgmm$W  # rete stimata
res$unpenalisedgmm$beta  # coefficienti beta
###################################################

# ==============================================================================
# Script di Replica per "Identifying Network Ties from Panel Data"
# Adattato per il dataset utente: panel_data_final.csv
# ==============================================================================

# 1. Caricamento delle librerie degli autori
# Assicurati che la cartella "PACKAGE" sia presente nella tua working directory.
# Il file README indica che loadpackage.R carica tutte le dipendenze necessarie.
# Installa devtools e le altre dipendenze probabili
install.packages(c("devtools", "cvxr", "kableExtra", "dplyr", "magrittr"))

# Load required libraries
library(devtools)
library(cvxr)
library(kableExtra)
library(dplyr)
library(magrittr)


relativepath = "./PACKAGE/"

if(!file.exists(paste0(relativepath, "loadpackage.R"))) {
  stop("Errore: Impossibile trovare './PACKAGE/loadpackage.R'. Assicurati di essere nella cartella corretta del replication package.")
}

suppressMessages(source(paste(relativepath, "loadpackage.R", sep=""))) 
# [cite: 814]

# 2. Caricamento e Preparazione Dati
# Leggiamo il tuo file csv
data_user <- read.csv("data/panel_data_final.csv")

# Verifica dimensioni: N=27, T=6
cat("Dimensioni dataset:", nrow(data_user), "righe.\n")
cat("Individui unici:", length(unique(data_user$id)), "\n")
cat("Periodi temporali:", length(unique(data_user$time)), "\n")

# Ordiniamo per Tempo e poi ID (Best practice per dati panel)
data_user <- data_user[order(data_user$time, data_user$id),]

# Verifica struttura: y, x1, x2, x3, x4, id, time
head(data_user)

# NUOVO PASSO CRITICO: Rimuovere tutte le righe con valori mancanti (NA)
# Questo risolve l'errore "x has missing values"
cat("Righe prima della pulizia NA:", nrow(data_user), "\n")
data_user <- na.omit(data_user)
cat("Righe dopo la pulizia NA:", nrow(data_user), "\n")
# 3. Configurazione Parametri
# Definiamo i parametri lambda come suggerito nel tuo snippet e nel README.
# Vettore lambda: c(p1, p2, p1_star)
# p1 = penalità Lasso step 1, p2 = Ridge, p1_star = Lasso step 2 (adaptive)
lambda_vec <- c(0.10, 0.10, 0.10) 
# [cite: 821]

# 4. Esecuzione della Stima (recoverNetwork)
# exoeffects=1 include gli effetti sociali esogeni (gamma), come nel modello teorico (Eq. 1 del paper).
# [cite: 51, 822]
cat("Avvio stima recoverNetwork... (potrebbe richiedere alcuni minuti)\n")

rn <- recoverNetwork(data = data_user, lambda = lambda_vec, exoeffects = 1)

# 5. Estrazione Risultati
# I risultati principali si trovano nella lista "unpenalisedgmm" (terzo step di stima)
# [cite: 821]

# A. Parametri scalari
rho_est   <- rn$unpenalisedgmm$rho   # Effetto endogeno (Social Multiplier)
beta_est  <- rn$unpenalisedgmm$beta  # Effetti propri
gamma_est <- rn$unpenalisedgmm$gamma # Effetti esogeni

cat("\n--- Risultati Stima ---\n")
cat("Rho (Endogenous Effect):", rho_est, "\n")
cat("Beta (Own Covariates):\n"); print(beta_est)
cat("Gamma (Exogenous Social Effects):\n"); print(gamma_est)

# B. Matrice delle Interazioni Sociali (W)
W_est <- rn$unpenalisedgmm$W

cat("\n--- Matrice W (Primi 5 nodi) ---\n")
print(W_est[1:5, 1:5])

# C. Analisi della Sparsità
num_links <- sum(W_est != 0)
total_possible <- nrow(W_est)^2 - nrow(W_est) # Esclude la diagonale
sparsity <- 1 - (num_links / total_possible)

cat("\nLink identificati:", num_links, "su", total_possible, "possibili.\n")
cat("Sparsità della rete:", round(sparsity * 100, 2), "%\n")

# 6. Salvataggio Output
write.csv(W_est, "W_stimata_output.csv")
cat("\nMatrice W salvata in 'W_stimata_output.csv'.\n")




# 1. Contiamo quante osservazioni ha ogni ID
count_obs <- data_user %>%
  group_by(id) %>%
  summarise(n_obs = n())

# 2. Identifichiamo gli ID che non hanno T=6 osservazioni
id_to_remove <- count_obs %>%
  filter(n_obs != 6) %>%
  pull(id)

cat("\nIndividui da rimuovere (hanno meno di 6 osservazioni):", id_to_remove, "\n")

# 3. Rimuoviamo questi ID dal dataset
data_balanced <- data_user %>%
  filter(!(id %in% id_to_remove))

cat("Righe nel dataset bilanciato:", nrow(data_balanced), "\n")
cat("Individui nel dataset bilanciato (N):", length(unique(data_balanced$id)), "\n")

# 4. Eseguiamo la stima con il dataset bilanciato
lambda_vec <- c(0.05, 0.02, 0.10)
cat("\nAvvio stima recoverNetwork sul dataset bilanciato...\n")

# Usa 'data_balanced' invece di 'data_user'
rn <- recoverNetwork(data = data_balanced, lambda = lambda_vec, exoeffects = 1)

# Estrazione Risultati (come prima)
W_est <- rn$unpenalisedgmm$W
cat("\n--- Risultati Stima ---\n")
cat("Rho (Endogenous Effect):", rn$unpenalisedgmm$rho, "\n")
cat("\nMatrice W (Primi 5 nodi) -> Dovrebbe essere la tua soluzione!\n")
print(W_est[1:5, 1:5])



#data_balanced <- data_balanced %>%
#  mutate(across(c(y, x1, x2, x3, x4), scale))
#
#rn <- recoverNetwork(
#  data = data_balanced,
#  lambda = c(0.1, 0.1, 0.1),
#  exoeffects = 0,
#  docv = 0,
#  timeweights = 0,
#  dopostols = 0,
#  adaestmethod = "initial"
#)

