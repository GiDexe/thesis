library(dplyr)
library(WDI)

# =============================================================================
# LABOUR FORCE DATA
# =============================================================================

oecd_iso <- c(
  "AUS", "AUT", "BEL", "CAN", "CHL", "COL", "CZE", "DNK", "EST", "FIN",
  "FRA", "DEU", "GRC", "HUN", "ISL", "IRL", "ISR", "ITA", "JPN", "KOR",
  "LVA", "LTU", "LUX", "MEX", "NLD", "NZL", "NOR", "POL", "PRT", "SVK",
  "SVN", "ESP", "SWE", "CHE", "TUR", "GBR", "USA"
)

lf <- WDI(
  country = oecd_iso,
  indicator = "SL.TLF.TOTL.IN",
  start = 2017, end = 2024
)
lf <- lf[, c("iso2c", "country", "year", "SL.TLF.TOTL.IN")]
colnames(lf) <- c("iso2c", "country", "time", "labor_force")

# Aggiungi Costa Rica
cr_lf <- WDI(
  country = "CR",
  indicator = "SL.TLF.TOTL.IN",
  start = 2017, end = 2024
)
cr_lf <- cr_lf[, c("iso2c", "country", "year", "SL.TLF.TOTL.IN")]
colnames(cr_lf) <- c("iso2c", "country", "time", "labor_force")
lf <- rbind(lf, cr_lf)

# =============================================================================
# COUNTRY MAP
# =============================================================================

country_map <- data.frame(
  country = c(
    "Australia", "Austria", "Belgium", "Canada", "Chile", "Colombia", "Costa Rica",
    "Denmark", "Estonia", "Finland", "France", "Germany", "Greece", "Hungary",
    "Iceland", "Ireland", "Israel", "Italy", "Japan", "Korea, Rep.", "Latvia",
    "Lithuania", "Luxembourg", "Mexico", "Netherlands", "New Zealand", "Norway",
    "Poland", "Portugal", "Slovak Republic", "Slovenia", "Spain", "Sweden",
    "Switzerland", "Turkiye", "United Kingdom", "United States"
  ),
  id = 1:37,
  iso2c = c(
    "AU", "AT", "BE", "CA", "CL", "CO", "CR", "DK", "EE", "FI", "FR", "DE", "GR", "HU",
    "IS", "IE", "IL", "IT", "JP", "KR", "LV", "LT", "LU", "MX", "NL", "NZ", "NO",
    "PL", "PT", "SK", "SI", "ES", "SE", "CH", "TR", "GB", "US"
  )
)

lf_mapped <- lf %>%
  left_join(country_map[, c("iso2c", "id")], by = "iso2c")

# =============================================================================
# FUNZIONE DI NORMALIZZAZIONE
# =============================================================================

normalize_lf <- function(data, name) {
  out <- data %>%
    left_join(lf_mapped[, c("id", "time", "labor_force")], by = c("id", "time")) %>%
    mutate(y = (y / labor_force) * 1000)

  out <- out[, c("time", "y", "x1", "x2", "x3", "x4", "id")]

  cat("===", name, "===\n")
  cat("NAs:", sum(is.na(out)), "\n")
  print(summary(out$y))

  dat_sorted <- out[order(out$id, out$time), ]
  y_by_id <- tapply(dat_sorted$y, dat_sorted$id, function(x) x)
  dy_mat <- t(sapply(y_by_id, diff))
  cat("Rank:", qr(dy_mat)$rank, "/ 37\n\n")

  return(out)
}

# =============================================================================
# APPLICA E SALVA
# =============================================================================

data_37001 <- readRDS("data/output_data/data_balanced_37001.rds")
data_45001 <- readRDS("data/output_data/data_balanced_45001.rds")

data_37001_lf <- normalize_lf(data_37001, "ISO 37001")
data_45001_lf <- normalize_lf(data_45001, "ISO 45001")

saveRDS(data_37001_lf, "data/output_data/data_balanced_37001_lf.rds")
saveRDS(data_45001_lf, "data/output_data/data_balanced_45001_lf.rds")
