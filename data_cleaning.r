# Load required libraries
library(readxl)
library(tidyverse)
library(dplyr)
library(tidyr)
library(janitor)
library(purrr)

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
  clean_names() %>%
  rename(code = country) %>%
  rename(country = land_sector)