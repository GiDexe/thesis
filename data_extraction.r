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
source("/workspaces/thesis/functions.r")

# libraries
library(readxl)
library(tidyverse)
library(janitor)
library(purrr)
library(WDI)
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

environment <- purrr::map(all_data, ~ keep(.x, ~ any(grepl("14001", colnames(.)))))


# Adjusting tibble structure
adjusted_data <- purrr::map(corrupt, ~ purrr::map(.x, ~ {
  .x %>%
    setNames(ifelse(names(.) == "Year", "Year", as.character(slice(., 2)))) %>%
    slice(-2) %>%
    dplyr::mutate(across(-1, ~ replace_na(as.numeric(.), 0))) %>%
    dplyr::mutate(Tot = rowSums(dplyr::select(., where(is.numeric) & !last_col()), na.rm = TRUE), .before = Year)
}))

eu <- c(
  "Austria", "Belgium", "Bulgaria", "Croatia", "Cyprus", "Czech Republic",
  "Denmark", "Estonia", "Finland", "France", "Germany", "Greece", "Hungary",
  "Ireland", "Italy", "Latvia", "Lithuania", "Luxembourg", "Malta",
  "Netherlands", "Poland", "Portugal", "Romania", "Slovakia", "Slovenia",
  "Spain", "Sweden"
)

my_panel <- adjusted_data %>%
  imap_dfr(~ bind_rows(.x, .id = "sheet") %>%
    dplyr::mutate(country = .y, .before = 1))%>%
  dplyr::select(country, `Land/Sector`, Tot, Year) %>%
  dplyr::filter(`Land/Sector` != "sector number") %>%
  dplyr::filter(`Land/Sector` != "Sector number") %>%
  clean_names() %>%
  rename(code = country) %>%
  rename(country = land_sector)%>%
  group_by(country) %>%
  arrange(year) %>%
  dplyr::filter(country %in% eu)%>%
  dplyr::mutate(tot = replace_na(tot, 0))%>%
  dplyr::select(-code)

# Covariates
covariates <- read_excel("/workspaces/thesis/data/cpds-1960-2023-update-2025.xlsx") %>%
  dplyr::filter(eu == 1 & year >= 2018) %>%
  dplyr::filter(country != "United Kingdom") %>%
  dplyr::select(country, year,  realgdpgr, unemp, gov_party, countryn) %>%
  rename(id = countryn)

trade_wb <- WDI(
  country = c("AT","BE","BG","HR","CY","CZ","DK","EE","FI","FR",
              "DE","GR","HU","IE","IT","LV","LT","LU","MT","NL",
              "PL","PT","RO","SK","SI","ES","SE"),
  indicator = "NE.TRD.GNFS.ZS",
  start = 2018,
  end = 2024
) %>%
dplyr::select(country, year, NE.TRD.GNFS.ZS) 

saveRDS(trade_wb, file = "tradewb.rds")


covariates <- covariates %>%
  left_join(trade_wb, by = c("country", "year")) %>%
  rename(gdp_pc = realgdpgr) %>%
  rename(unemployment = unemp) %>%
  rename(trade_gdp = NE.TRD.GNFS.ZS)

#is.na(covariates)
#covariates %>%
#      filter(is.na(trade_gdp)) %>%
#     dplyr:: select(country, year, trade_gdp)
#CZH and SLK are missing trade data. I will estimate with and without them



# this is something i just changed
data_tb <- my_panel %>%
  right_join(covariates, by = c("country", "year")) %>%
  dplyr::filter(country != "Czech Republic" & country != "Slovakia") %>%
  dplyr::mutate(tot = ifelse(is.na(tot), 0, tot)) %>%
  rename(
    time = year,
    y = tot,
    x1 = gdp_pc,
    x2 = unemployment,
    x3 = gov_party,
    x4 = trade_gdp
  ) %>%
  ungroup() %>%
  dplyr::mutate(id = group_indices(., country)) %>%
  dplyr::select(-country)

# T per id
Ti <- data_tb %>%
  group_by(id) %>%
  summarise(Ti = n_distinct(time))

Ti %>% arrange(Ti)

###

# 4. Balanced estimate
lambda_vec <- c(0.05, 0.02, 1.10)
cat("\nEstimating recoverNetwork on balanced dataset...\n")

# Use 'data_tb' 
rn <- recoverNetwork(data = data_tb, lambda = lambda_vec, exoeffects = 1)

# Results
W_est <- rn$unpenalisedgmm$W
cat("\n--- Results ---\n")
cat("Rho (Endogenous Effect):", rn$unpenalisedgmm$rho, "\n")
cat("\nMatrix W (First 5 nodes)\n")
print(W_est[1:5, 1:5])

saveRDS(data_tb, file = "data_tb.rds")

# T per id
Ti <- data_tb %>%
  group_by(id) %>%
  summarise(Ti = n_distinct(time))

Ti %>% arrange(Ti)
