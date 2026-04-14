# =============================================================================
# 01_cleaning.R
# Produce i panel data bilanciati per ISO 14001 e ISO 9001, li salva come RDS
# N countries, T>=7 (2017 onwards)
# NO COVARIATES
# =============================================================================

library(dplyr)
library(readxl)
library(stringr)
library(purrr)
library(tidyr)
library(janitor)

# === SETTINGS ================================================================
set.seed(1861)
data_dir <- "/workspaces/thesis/data"
out_data_dir <- "data/output_data"

eu <- c(
  "Austria", "Belgium", "Bulgaria", "Croatia", "Cyprus", "Czech Republic",
  "Denmark", "Estonia", "Finland", "France", "Germany", "Greece",
  "Hungary", "Ireland", "Italy", "Latvia", "Lithuania", "Luxembourg",
  "Malta", "Netherlands", "Poland", "Portugal", "Romania", "Slovak Republic",
  "Slovenia", "Spain", "Sweden"
)

oecd <- c(
  "Australia", "Austria", "Belgium", "Canada", "Chile", "Colombia",
  "Costa Rica", "Czech Republic", "Denmark", "Estonia", "Finland",
  "France", "Germany", "Greece", "Hungary", "Iceland", "Ireland",
  "Israel", "Italy", "Japan", "Korea", "Latvia", "Lithuania",
  "Luxembourg", "Mexico", "Netherlands", "New Zealand", "Norway",
  "Poland", "Portugal", "Slovak Republic", "Slovenia", "Spain",
  "Sweden", "Switzerland", "Turkey", "United Kingdom", "United States"
)

# === HELPER FUNCTIONS ========================================================

extract_sheets <- function(file_path) {
  year <- str_extract(basename(file_path), "20(1[7-9]|2[0-4])")
  sheet_names <- excel_sheets(file_path)
  sheets_to_read <- sheet_names[-1]

  data_list <- purrr::map(sheets_to_read, ~ {
    tryCatch(
      {
        read_excel(file_path, sheet = .x) %>%
          dplyr::mutate(Year = year)
      },
      error = function(e) {
        NULL
      }
    )
  })

  data_list <- purrr::discard(data_list, is.null)
  names(data_list) <- sheets_to_read[seq_along(data_list)]
  return(data_list)
}

build_panel <- function(iso_code, country_filter = NULL) {
  # --- Read Excel files ---
  all_excel_files <- list.files(data_dir,
    pattern = "\\.xlsx$|\\.xls$|\\.xlsm$",
    full.names = TRUE
  )
  excel_files <- all_excel_files[grepl("ISO Survey 20", basename(all_excel_files))]
  all_data <- purrr::map(excel_files, extract_sheets)

  # --- Filter to ISO code ---
  iso_data <- purrr::map(all_data, ~ purrr::keep(.x, ~ any(grepl(iso_code, colnames(.)))))

  # --- Adjust tibble structure ---
  adjusted_data <- purrr::map(iso_data, ~ purrr::map(.x, ~ {
    if (ncol(.) < 2) {
      return(NULL)
    }

    tryCatch(
      {
        .x %>%
          setNames(ifelse(names(.) == "Year", "Year", as.character(dplyr::slice(., 2)))) %>%
          dplyr::slice(-2) %>%
          dplyr::mutate(across(-1, ~ replace_na(as.numeric(.), 0))) %>%
          dplyr::mutate(Tot = rowSums(dplyr::select(., where(is.numeric) & !last_col()), na.rm = TRUE), .before = Year)
      },
      error = function(e) {
        NULL
      }
    )
  }))

  adjusted_data <- purrr::map(adjusted_data, ~ purrr::discard(.x, is.null))

  my_panel <- adjusted_data %>%
    purrr::imap_dfr(~ bind_rows(.x, .id = "sheet") %>%
      dplyr::mutate(country = .y, .before = 1))

  my_panel <- my_panel %>%
    dplyr::select(country, `Land/Sector`, Tot, Year) %>%
    dplyr::filter(`Land/Sector` != "sector number") %>%
    dplyr::filter(`Land/Sector` != "Sector number") %>%
    clean_names() %>%
    dplyr::rename(code = country) %>%
    dplyr::rename(country = land_sector)

  my_panel <- my_panel %>%
    dplyr::group_by(country) %>%
    dplyr::arrange(year) %>%
    dplyr::mutate(delta = tot - lag(tot, default = 0))

  # --- Harmonize ISO country names to OECD standard ---
  my_panel <- my_panel %>%
    dplyr::ungroup() %>%
    dplyr::mutate(country = case_match(country,
      "Korea (Republic of)" ~ "Korea",
      "Korea (the Republic of)" ~ "Korea",
      "Czech Republic" ~ "Czech Republic",
      "Czechia" ~ "Czech Republic",
      "Slovakia" ~ "Slovak Republic",
      "Turkey" ~ "Turkey",
      "United Kingdom of Great Britain and Northern Ireland" ~ "United Kingdom",
      "United Kingdom of Great Britain and Northern Ireland (the)" ~ "United Kingdom",
      "United States of America" ~ "United States",
      "United States of America (the)" ~ "United States",
      "Netherlands (Kingdom of the)" ~ "Netherlands",
      .default = country
    ))

  # --- Complete panel ---
  my_panel <- my_panel %>%
    dplyr::ungroup() %>%
    tidyr::complete(country, year, fill = list(tot = NA)) %>%
    dplyr::group_by(country) %>%
    dplyr::arrange(year) %>%
    tidyr::fill(tot, .direction = "down") %>%
    dplyr::mutate(tot = replace_na(tot, 0)) %>%
    dplyr::mutate(delta = tot - lag(tot, default = 0)) %>%
    dplyr::ungroup()

  # --- Filter to 2017 onwards and apply country filter if specified ---
  panel_data <- my_panel %>%
    dplyr::filter(year >= 2017)

  if (!is.null(country_filter)) {
    panel_data <- panel_data %>%
      dplyr::filter(country %in% country_filter)
  }

  # --- Balance check ---
  obs_per_country <- panel_data %>%
    dplyr::group_by(country) %>%
    dplyr::summarise(n_obs = dplyr::n(), .groups = "drop")

  T_max <- max(obs_per_country$n_obs)
  countries_to_keep <- obs_per_country %>%
    dplyr::filter(n_obs == T_max) %>%
    dplyr::pull(country)

  panel_data <- panel_data %>%
    dplyr::filter(country %in% countries_to_keep)

  # --- Save country mapping, then assign consecutive IDs ---
  country_map <- panel_data %>%
    dplyr::distinct(country) %>%
    dplyr::arrange(country) %>%
    dplyr::mutate(id = dplyr::row_number())

  panel_data <- panel_data %>%
    dplyr::select(-delta) %>%
    dplyr::left_join(country_map, by = "country")

  # --- Keep only time, y, id (no covariates) ---
  data_balanced <- panel_data %>%
    dplyr::select(time = year, y = tot, id) %>%
    dplyr::arrange(time, id)

  data_balanced <- na.omit(data_balanced)

  cat(paste0("ISO ", iso_code, ":\n"))
  cat("  N:", length(unique(data_balanced$id)), "\n")
  cat("  T:", length(unique(data_balanced$time)), "\n")
  cat("  Rows:", nrow(data_balanced), "\n")
  cat("  ID range:", range(data_balanced$id), "\n")
  cat("  Years:", sort(unique(data_balanced$time)), "\n\n")

  return(list(data = data_balanced, country_map = country_map))
}

# === BUILD AND SAVE ==========================================================

cat("=== Building ISO 14001 (EU) panel ===\n")
panel_14001_eu <- build_panel("14001", country_filter = eu)

cat("=== Building ISO 14001 (OECD) panel ===\n")
panel_14001_oecd <- build_panel("14001", country_filter = oecd)

cat("=== Building ISO 9001 (EU) panel ===\n")
panel_9001_eu <- build_panel("9001", country_filter = eu)

cat("=== Building ISO 9001 (OECD) panel ===\n")
panel_9001_oecd <- build_panel("9001", country_filter = oecd)

# Save
saveRDS(panel_14001_eu$data, file = file.path(out_data_dir, "data_balanced_14001_eu.rds"))
saveRDS(panel_14001_eu$country_map, file = file.path(out_data_dir, "country_map_14001_eu.rds"))

saveRDS(panel_14001_oecd$data, file = file.path(out_data_dir, "data_balanced_14001_oecd.rds"))
saveRDS(panel_14001_oecd$country_map, file = file.path(out_data_dir, "country_map_14001_oecd.rds"))

saveRDS(panel_9001_eu$data, file = file.path(out_data_dir, "data_balanced_9001_eu.rds"))
saveRDS(panel_9001_eu$country_map, file = file.path(out_data_dir, "country_map_9001_eu.rds"))

saveRDS(panel_9001_oecd$data, file = file.path(out_data_dir, "data_balanced_9001_oecd.rds"))
saveRDS(panel_9001_oecd$country_map, file = file.path(out_data_dir, "country_map_9001_oecd.rds"))

cat("=== Files Saved ===\n")
cat("  14001 EU: data_balanced_14001_eu.rds, country_map_14001_eu.rds\n")
cat("  14001 OECD: data_balanced_14001_oecd.rds, country_map_14001_oecd.rds\n")
cat("  9001 EU: data_balanced_9001_eu.rds, country_map_9001_eu.rds\n")
cat("  9001 OECD: data_balanced_9001_oecd.rds, country_map_9001_oecd.rds\n")


### combining them

# =============================================================================
# 02_combine_data.R
# Combine 2017+ and pre-2017 data for ISO 14001 and 9001
# Output: Single long-format dataset per ISO × region
# =============================================================================

library(dplyr)
library(readr)

data_dir <- "data/output_data"

# === COMBINE DATA ============================================================

combine_datasets <- function(iso_code, region) {
  # --- Read 2017+ data (from 01_cleaning.R) ---
  data_2017plus <- readRDS(
    file.path(data_dir, paste0("data_balanced_", iso_code, "_", region, ".rds"))
  )
  country_map_2017plus <- readRDS(
    file.path(data_dir, paste0("country_map_", iso_code, "_", region, ".rds"))
  )

  # --- Read pre-2017 data (if exists) ---
  file_pre2017 <- file.path(data_dir, paste0("data_balanced_", iso_code, "_", region, "_pre2017.rds"))

  if (file.exists(file_pre2017)) {
    data_pre2017 <- readRDS(file_pre2017)
    country_map_pre2017 <- readRDS(
      file.path(data_dir, paste0("country_map_", iso_code, "_", region, "_pre2017.rds"))
    )

    # --- Merge country maps and reindex IDs ---
    all_countries <- bind_rows(
      country_map_2017plus %>% dplyr::select(country),
      country_map_pre2017 %>% dplyr::select(country)
    ) %>%
      distinct() %>%
      arrange(country) %>%
      mutate(id = row_number())

    # --- Reindex 2017+ data ---
    remap_2017plus <- country_map_2017plus %>%
      left_join(all_countries, by = "country", suffix = c(".old", ".new")) %>%
      select(id.old, id.new = id)

    data_2017plus <- data_2017plus %>%
      left_join(remap_2017plus, by = c("id" = "id.old")) %>%
      select(time, y, id = id.new)

    # --- Reindex pre-2017 data ---
    remap_pre2017 <- country_map_pre2017 %>%
      left_join(all_countries, by = "country", suffix = c(".old", ".new")) %>%
      select(id.old, id.new = id)

    data_pre2017 <- data_pre2017 %>%
      left_join(remap_pre2017, by = c("id" = "id.old")) %>%
      select(time, y, id = id.new)

    # --- Combine ---
    combined <- bind_rows(data_pre2017, data_2017plus) %>%
      arrange(time, id) %>%
      na.omit()
  } else {
    cat("Warning: No pre-2017 file found for", iso_code, "-", region, "\n")
    combined <- data_2017plus
    all_countries <- country_map_2017plus
  }

  return(list(data = combined, country_map = all_countries))
}

# === BUILD ALL COMBINATIONS ==================================================

cat("=== Combining ISO 14001 (EU) ===\n")
c14001_eu <- combine_datasets("14001", "eu")
saveRDS(c14001_eu$data, file = file.path(data_dir, "data_combined_14001_eu.rds"))
saveRDS(c14001_eu$country_map, file = file.path(data_dir, "country_map_combined_14001_eu.rds"))
cat("  N:", n_distinct(c14001_eu$data$id), "| T:", min(c14001_eu$data$time), "-", max(c14001_eu$data$time), "\n\n")

cat("=== Combining ISO 14001 (OECD) ===\n")
c14001_oecd <- combine_datasets("14001", "oecd")
saveRDS(c14001_oecd$data, file = file.path(data_dir, "data_combined_14001_oecd.rds"))
saveRDS(c14001_oecd$country_map, file = file.path(data_dir, "country_map_combined_14001_oecd.rds"))
cat("  N:", n_distinct(c14001_oecd$data$id), "| T:", min(c14001_oecd$data$time), "-", max(c14001_oecd$data$time), "\n\n")

cat("=== Combining ISO 9001 (EU) ===\n")
c9001_eu <- combine_datasets("9001", "eu")
saveRDS(c9001_eu$data, file = file.path(data_dir, "data_combined_9001_eu.rds"))
saveRDS(c9001_eu$country_map, file = file.path(data_dir, "country_map_combined_9001_eu.rds"))
cat("  N:", n_distinct(c9001_eu$data$id), "| T:", min(c9001_eu$data$time), "-", max(c9001_eu$data$time), "\n\n")

cat("=== Combining ISO 9001 (OECD) ===\n")
c9001_oecd <- combine_datasets("9001", "oecd")
saveRDS(c9001_oecd$data, file = file.path(data_dir, "data_combined_9001_oecd.rds"))
saveRDS(c9001_oecd$country_map, file = file.path(data_dir, "country_map_combined_9001_oecd.rds"))
cat("  N:", n_distinct(c9001_oecd$data$id), "| T:", min(c9001_oecd$data$time), "-", max(c9001_oecd$data$time), "\n\n")

cat("=== Files Saved ===\n")
cat("  data_combined_14001_eu.rds\n")
cat("  data_combined_14001_oecd.rds\n")
cat("  data_combined_9001_eu.rds\n")
cat("  data_combined_9001_oecd.rds\n")
cat("  (+ country_map_combined_*.rds for each)\n")
