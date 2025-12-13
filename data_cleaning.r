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


g <- all_data[[1]]


g[1]


corrupt <- map(all_data, ~ keep(.x, ~ any(grepl("37001", colnames(.)))))

colnames(corrupt[[3]][[1]])

corrupt[[3]][[1]]$Year