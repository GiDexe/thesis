# Load libraries
library(devtools)
library(kableExtra)
library(dplyr)
library(magrittr)
library(Matrix)
library(glmnet)
library(CVXR)
library(MASS)
library(optimx)
library(fixest)

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

# Manually Solving my conflicts
select <- dplyr::select
filter <- dplyr::filter
map <- purrr::map
summarize <- dplyr::summarize
mutate <- dplyr::mutate


# Adjusting tibble structure
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
  arrange(year) %>% # ascending
  mutate(delta = tot - lag(tot, default = 0))

 #controls, old bit
#covariates <- WDI(
#  country = "all",
#  indicator = c(
#    "NY.GDP.PCAP.PP.KD", # GDP per capita (constant PPP)
#    "CC.EST", # Control of Corruption
#    "RQ.EST", # Regulatory Quality
#    "SL.UEM.TOTL.ZS" # Unemployment rate
#  ),
#  start = 2018,
#  end = 2024,
#  extra = F
#) %>%
#  select(country, year, NY.GDP.PCAP.PP.KD, CC.EST, RQ.EST, SL.UEM.TOTL.ZS) %>%
#  rename(
#    gdp_pc = NY.GDP.PCAP.PP.KD,
#    corruption_control = CC.EST,
#    regulatory_quality = RQ.EST,
#    unemployment = SL.UEM.TOTL.ZS
#  ) %>%
#  filter(!is.na(gdp_pc))
#
#saveRDS(covariates, file = "covariates.rds")
covariates <- read_rds("covariates.rds")
#---

# FROM HERE CHECK NAMES AGAIN
# Take only countries appearing in 2020 for lighter computation
#countries_2020 <- my_panel %>%
#  filter(year == 2020) %>%
#  pull(country) %>%
#  unique()


my_panel <- my_panel %>%
  #filter(country %in% countries_2020) %>%
  ungroup() %>%
  complete(
    country,
    year,
    fill = list(tot = NA)
  ) %>%
  group_by(country) %>%
  arrange(year) %>%
  # Last Observation Carried Forward):
  # Ensuring tot_{t-1} - tot_{t-1} = 0.
  fill(tot, .direction = "down") %>%
  mutate(tot = replace_na(tot, 0)) %>%
  mutate(delta = tot - lag(tot, default = 0)) %>%
  ungroup()


# this is something i just changed
panel_data <- my_panel %>%
  left_join(covariates, by = c("country", "year")) %>%
  filter(!is.na(gdp_pc)) %>%
  rename(
    time = year,
    y = tot,
    x1 = gdp_pc,
    x2 = rule_of_law,
    x3 = trade,
    x4 = unemployment
  )

# id
panel_data <- panel_data %>%
  ungroup() %>%
  mutate(id = group_indices(., country))
panel_balance <- panel_data %>%
  group_by(id) %>%
  summarize(country = first(country), n_obs = n(), .groups = "drop")

complete_ids <- panel_balance %>%
  filter(n_obs == max(n_obs)) %>%
  pull(id)

panel_data <- panel_data %>%
  filter(id %in% complete_ids) %>%
  dplyr::select(-code)

required_packages <- c("Matrix", "glmnet", "CVXR", "MASS", "optimx", "fixest")

# Check and install missing packages? Attention, conflixts are present. Need for a better structure.
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) {
  install.packages(new_packages)
}

# Create OECD countries vector
oecd <- c("Australia", "Austria", "Belgium", "Canada", "Chile", "Colombia", 
  "Costa Rica", "Czech Republic", "Denmark", "Estonia", "Finland", 
  "France", "Germany", "Greece", "Hungary", "Iceland", "Ireland", 
  "Israel", "Italy", "Japan", "Korea", "Latvia", "Lithuania", 
  "Luxembourg", "Mexico", "Netherlands", "New Zealand", "Norway", 
  "Poland", "Portugal", "Slovak Republic", "Slovenia", "Spain", 
  "Sweden", "Switzerland", "Turkey", "United Kingdom", "United States")

panel_data_clean <- panel_data %>%
  filter(country %in% oecd) %>%
  dplyr::select(-country, -delta)
MODEL_COLS <- c("y", "x1", "x2", "x3", "x4", "id", "time")

panel_data_final <- panel_data_clean %>%
  drop_na(all_of(MODEL_COLS))

MODEL_COLS <- c("y", "x1", "x2", "x3", "x4", "id", "time")

#data_user <- read.csv("data/panel_data_final.csv")
data_user <- panel_data_clean
# DImension Check
cat("Dimensioni dataset:", nrow(data_user), "righe.\n")
cat("Individui unici:", length(unique(data_user$id)), "\n")
cat("Periodi temporali:", length(unique(data_user$time)), "\n")

data_user <- data_user[order(data_user$time, data_user$id), ]

#  y, x1, x2, x3, x4, id, time
head(data_user)

data_user <- na.omit(data_user)


# Tricky: I clearly had some issues with merging some steps ago - Ideally I will remove this further section

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
rn_37 <- recoverNetwork(data = data_balanced,lambda = lambda_vec, exoeffects = 1)

# Estrazione Risultati (come prima)
W_est_37 <- rn_37$unpenalisedgmm$W
cat("\n--- Risultati Stima ---\n")
cat("Rho (Endogenous Effect):", rn_37$unpenalisedgmm$rho, "\n")
cat("\nMatrice W (Primi 5 nodi) -> Dovrebbe essere la tua soluzione!\n")
print(W_est_37[1:5, 1:5])

saveRDS(W_est_37, file = "matrice_W_stima_37001.rds")
saveRDS(rho_37, file = "rn_37$unpenalisedgmm$rho.rds")

#readRDS("matrice_W_stima.rds")


library(igraph)
library(ggraph)
library(tidygraph)
library(patchwork)

W_37001 <- readRDS("matrice_W_stima_37001.rds")
W_45001 <- readRDS("matrice_W_stima_45001.rds")

cat("\n=== CONFRONTO STATISTICO ===\n")
cat("\nISO 37001:\n")
cat("  Links (>0.01):", sum(W_37001 > 0.01), "\n")
cat("  Density:", sum(W_37001 > 0.01) / (nrow(W_37001)^2 - nrow(W_37001)), "\n")
cat("  Max weight:", max(W_37001), "\n")
cat("  Mean weight (non-zero):", mean(W_37001[W_37001 > 0.01]), "\n")

cat("\nISO 45001:\n")
cat("  Links (>0.01):", sum(W_45001 > 0.01), "\n")
cat("  Density:", sum(W_45001 > 0.01) / (nrow(W_45001)^2 - nrow(W_45001)), "\n")
cat("  Max weight:", max(W_45001), "\n")
cat("  Mean weight (non-zero):", mean(W_45001[W_45001 > 0.01]), "\n")

cat("\n=== CORRELAZIONE TRA MATRICI ===\n")
cat("Pearson:", cor(as.vector(W_37001), as.vector(W_45001)), "\n")
cat("Spearman:", cor(as.vector(W_37001), as.vector(W_45001), method = "spearman"), "\n")


cat("\n=== OVERLAP LINKS ===\n")
links_37001 <- which(W_37001 > 0.01, arr.ind = TRUE)
links_45001 <- which(W_45001 > 0.01, arr.ind = TRUE)
common_links <- nrow(merge(as.data.frame(links_37001), as.data.frame(links_45001)))
cat("Common links:", common_links, "\n")
cat("Jaccard similarity:", common_links / (nrow(links_37001) + nrow(links_14001) - common_links), "\n")

threshold <- 0.01

g_37001 <- graph_from_adjacency_matrix(W_37001 > threshold, mode = "directed")
g_14001 <- graph_from_adjacency_matrix45 > threshold, mode = "directed")

p1 <- ggraph(g_37001, layout = "fr") + 
  geom_edge_link(arrow = arrow(length = unit(2, 'mm')), 
                 end_cap = circle(2, 'mm'), alpha = 0.3) +
  geom_node_point(size = 4, color = "darkred") +
  geom_node_text(aes(label = name), size = 2.5, repel = TRUE) +
  theme_void() +
  labs(title = "ISO 37001 (Anti-Corruption)")

p2 <- ggraph(g_14001, layout = "fr") + 
  geom_edge_link(arrow = arrow(length = unit(2, 'mm')), 
                 end_cap = circle(2, 'mm'), alpha = 0.3) +
  geom_node_point(size = 4, color = "darkgreen") +
  geom_node_text(aes(label = name), size = 2.5, repel = TRUE) +
  theme_void() +
  labs(title = "ISO 14001 (Environment)")

p1 / p2

ggsave("network_comparison.png", width = 12, height = 10)