# libraries
library(readxl)
library(tidyverse)
library(janitor)
library(purrr)
library(WDI)
library(devtools)
library(kableExtra)
library(magrittr)
library(Matrix)
library(glmnet)
library(CVXR)
library(MASS)
library(optimx)
library(fixest)
library(countrycode)
# depaula's functions
source("/workspaces/thesis/functions.r")

data_dir <- "/workspaces/thesis/data"

#-----------------------------------------------------------------------

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
  data_list <- purrr::map(sheets_to_read, ~ {
    read_excel(file_path, sheet = .x) %>%
      dplyr::mutate(Year = year)
  })

  # Name the list elements
  names(data_list) <- sheets_to_read

  return(data_list)
}


# Extract sheets from all files
all_data <- purrr::map(excel_files, extract_sheets)

# Nested structure
g <- all_data[[1]]


g[1]


corrupt <- purrr::map(all_data, ~ keep(.x, ~ any(grepl("37001", colnames(.)))))


# Adjusting tibble structure
adjusted_data <- purrr::map(corrupt, ~ purrr::map(.x, ~ {
  .x %>%
    setNames(ifelse(names(.) == "Year", "Year", as.character(slice(., 2)))) %>%
    slice(-2) %>%
    dplyr::mutate(across(-1, ~ replace_na(as.numeric(.), 0))) %>%
    dplyr::mutate(Tot = rowSums(dplyr::select(., where(is.numeric) & !last_col()), na.rm = TRUE), .before = Year)
}))

my_panel <- adjusted_data %>%
  imap_dfr(~ bind_rows(.x, .id = "sheet") %>%
    dplyr::mutate(country = .y, .before = 1))

my_panel <- my_panel %>%
  dplyr::select(country, `Land/Sector`, Tot, Year) %>%
  dplyr::filter(`Land/Sector` != "sector number") %>%
  dplyr::filter(`Land/Sector` != "Sector number") %>%
  clean_names() %>%
  rename(code = country) %>%
  rename(country = land_sector)

my_panel <- my_panel %>%
  group_by(country) %>%
  arrange(year) %>% # ascending
  dplyr::mutate(delta = tot - lag(tot, default = 0))


# controls, old bit
covariates <- WDI(
  country = "all",
  indicator = c(
    "NY.GDP.PCAP.PP.KD", # GDP per capita (constant PPP)
    "CC.EST", # Control of Corruption
    "RQ.EST", # Regulatory Quality
    "SL.UEM.TOTL.ZS" # Unemployment rate
  ),
  start = 2018,
  end = 2024,
  extra = F
) %>%
  dplyr::select(country, year, NY.GDP.PCAP.PP.KD, CC.EST, RQ.EST, SL.UEM.TOTL.ZS) %>%
  rename(
    gdp_pc = NY.GDP.PCAP.PP.KD,
    corruption_control = CC.EST,
    regulatory_quality = RQ.EST,
    unemployment = SL.UEM.TOTL.ZS
  ) %>%
  dplyr::filter(!is.na(gdp_pc))

saveRDS(covariates, file = "covariates.rds")

#---

# FROM HERE CHECK NAMES AGAIN
# Take only countries appearing in 2020 for lighter computation
countries_2020 <- my_panel %>%
  dplyr::filter(year == 2020) %>%
  pull(country) %>%
  unique()


my_panel <- my_panel %>%
  dplyr::filter(country %in% countries_2020) %>%
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
  dplyr::mutate(tot = replace_na(tot, 0)) %>%
  dplyr::mutate(delta = tot - lag(tot, default = 0)) %>%
  ungroup()


# this is something i just changed
panel_data <- my_panel %>%
  left_join(covariates, by = c("country", "year")) %>%
  dplyr::filter(!is.na(gdp_pc)) %>%
  rename(
    time = year,
    y = tot,
    x1 = gdp_pc,
    x2 = corruption_control,
    x3 = regulatory_quality,
    x4 = unemployment
  )

# id
panel_data <- panel_data %>%
  ungroup() %>%
  dplyr::mutate(id = group_indices(., country))
panel_balance <- panel_data %>%
  group_by(id) %>%
  dplyr::summarize(country = first(country), n_obs = n(), .groups = "drop")

complete_ids <- panel_balance %>%
  dplyr::filter(n_obs == max(n_obs)) %>%
  pull(id)

panel_data <- panel_data %>%
  dplyr::filter(id %in% complete_ids) %>%
  dplyr::select(-code)

panel_data_clean <- panel_data %>%
  dplyr::select(-country, -delta)
MODEL_COLS <- c("y", "x1", "x2", "x3", "x4", "id", "time")

panel_data_final <- panel_data_clean %>%
  drop_na(all_of(MODEL_COLS))

MODEL_COLS <- c("y", "x1", "x2", "x3", "x4", "id", "time")

data_user <- read.csv("data/panel_data_final.csv")
data_user <- panel_data_clean
# DImension Check
cat("Dimensions of the dataset:", nrow(data_user), "righe.\n")
cat("Unique countries:", length(unique(data_user$id)), "\n")
cat("Periods:", length(unique(data_user$time)), "\n")

data_user <- data_user[order(data_user$time, data_user$id), ]

#  y, x1, x2, x3, x4, id, time
head(data_user)

data_user <- na.omit(data_user)


#  to be reviewed???

# 1. counting obs
count_obs <- data_user %>%
  group_by(id) %>%
  summarise(n_obs = n())

# 2. spotting ids not having at least 6 observations
id_to_remove <- count_obs %>%
  dplyr::filter(n_obs != 6) %>%
  pull(id)

cat("\nCountries to be removed - Less than 6 obs:", id_to_remove, "\n")

# 3. Removing those ids from the dataset
data_balanced <- data_user %>%
  dplyr::filter(!(id %in% id_to_remove))

cat("Nrows (balanced):", nrow(data_balanced), "\n")
cat("Countries in balanced dataset (N):", length(unique(data_balanced$id)), "\n")

# 4. Balanced estimate
lambda_vec <- c(0.05, 0.02, 0.10)
cat("\nEstimating recoverNetwork on balanced dataset...\n")

# Use 'data_balanced' 
rn <- recoverNetwork(data = data_balanced, lambda = lambda_vec, exoeffects = 1)

# Results
W_est <- rn$unpenalisedgmm$W
cat("\n--- Results ---\n")
cat("Rho (Endogenous Effect):", rn$unpenalisedgmm$rho, "\n")
cat("\nMatrix W (First 5 nodes)\n")
print(W_est[1:5, 1:5])

saveRDS(W_est, file = "matrice_W_stima_2.rds")

readRDS("matrice_W_stima_2.rds")
