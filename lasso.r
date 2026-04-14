library(readxl)
library(dplyr)

# File path
xlsx_path <- "data/wgidataset_with_sourcedata-2025.xlsx"

# Get all sheet names
sheet_names <- excel_sheets(xlsx_path)

# Read all sheets and add dimension column
wgi_list <- map_df(
  sheet_names,
  ~ read_excel(xlsx_path, sheet = .x) %>%
    mutate(dimension = .x, .after = `Economy (code)`)
)

# Optional: Clean up and select key columns
# Uncomment if you want a more compact version:
# wgi_list <- wgi_list %>%
#   select(
#     country_code = `Economy (code)`,
#     country_name = `Economy (name)`,
#     year = Year,
#     region = Region,
#     income = `Income classification`,
#     dimension,
#     estimate = `Governance estimate (approx. -2.5 to +2.5)`,
#     score = `Governance score (0-100)`,
#     std_error = `Standard error (estimate)`,
#     everything()
#   )

# Save to CSV
write.csv(wgi_list, "/mnt/user-data/outputs/wgi_combined_long.csv", row.names = FALSE)

# Quick check
cat("Combined dataset shape:", nrow(wgi_list), "rows x", ncol(wgi_list), "columns\n")
cat("Dimensions included:", paste(unique(wgi_list$dimension), collapse = ", "), "\n")
cat("Year range:", min(wgi_list$Year), "-", max(wgi_list$Year), "\n")
cat("Countries:", n_distinct(wgi_list$`Economy (code)`), "\n\n")
cat("First few rows:\n")
print(head(wgi_list[1:10]))

wgi_list %>%
  select(-1, -dimension, Region, `Economy (code)`) %>%
  rename(country = `Economy (name)`)
