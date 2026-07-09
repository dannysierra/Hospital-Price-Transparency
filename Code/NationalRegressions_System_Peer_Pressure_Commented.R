# ============================================================
# NationalRegressions_System_Peer_Pressure.R
# Hospital Price Transparency — National IV Analysis
# Updated: Corrected instrument construction (focal-county fix)
#          County primary: 9-month trailing window (F = 27.1)
#          City robustness: 12-month trailing window (F = 26.7)
#          CBSA dropped (instrument fails under two-way clustering)
#          Two-way clustering throughout: geo-unit + post_month
# ============================================================
#
# PAPER STRUCTURE THIS SCRIPT SUPPORTS:

#
# ── MAIN PAPER ────────────────────────────────────────────────────────────────
#
# Section 4 (Data):
#   Table 1  [tab:summstats]             Summary statistics, county level
#
# Section 5 (Empirical Strategy):
#   Table 2  [tab:first_stage_validity]  First-stage estimates, county panel
#   Figure 1 [fig:first_stage_wald]      First-stage Wald F by service group
#
# Section 6 (Main Results):
#   Table 3  [tab:pooled_iv]             Pooled IV, median price (primary)
#   Table 4  [tab:pooled_outcomes]       Pooled OLS vs IV, all price outcomes
#   Table 5  [tab:service_iv]            Service-level OLS vs IV, median
#
# Section 7 (Shoppability Gradient):
#   Table 6  [tab:metareg]              Shoppability meta-regression (Theory V2)
#   Figure 2 [fig:metareg_combined]     Forest plot + meta-regression scatter
#   Figure 3 [fig:ols_iv_forest]        OLS vs IV forest plot by service
#
# Section 8 (Saturation Dynamics):
#   Figure 4 [fig:saturation_shoppability]  Saturation bin effects
#   Figure 5 [fig:saturation_robust]        Saturation spaghetti across schemes
#
# Section 9 (Heterogeneity):
#   Table 7  [tab:het_main]             Market-level heterogeneity splits
#   Table 8  [tab:het_hospital_type]    Hospital type and proximity splits
#   Table 9  [tab:het_race_2x2]         Race x Medicaid two-by-two
#   Figure 6 [fig:het_forest]           Heterogeneity forest plot
#   Figure 7 [fig:race_2x2]             Race x Medicaid gradient heatmap
#
# Section 10 (Robustness):
#   Figure 8 [fig:rand_inference]       Randomization inference plot
#
# ── APPENDIX ──────────────────────────────────────────────────────────────────
#
# Appendix A (CPT/HCPCS Codes):
#   Table 10 [tab:radiology_codes]       Radiology CPT/HCPCS codes
#   Table 11 [tab:med_surg_codes]        Medicine and surgery codes
#
# Appendix B (Extended Main Results):
#   Table 12 [tab:pooled_ols_iv]         Pooled OLS/IV across all outcomes
#   Table 13 [tab:schemes]               Service classifications, nine schemes
#   Table 14 [tab:geographic_robustness] Two-way geographic comparison (county/city)
#   Table 15 [tab:app_pooled]            Pooled IV by market characteristic
#   Table 16 [tab:app_servicegroup]      Service-group IV by heterogeneity
#   Table 17 [tab:app_schemes]           Shoppability gradient, all nine schemes
#   Table 18 [tab:app_hosptype]          Hospital type: gradient, all nine schemes
#   Table 19 [tab:app_prox]             Geographic proximity: gradient, all nine
#
# Appendix C (Robustness):
#   Table 20 [tab:instrument_variants]   Instrument variant robustness, county
#   Table 21 [tab:loo]                   Leave-one-system-out, county
#   Table 22 [tab:clustering_robustness] Two-way vs one-way clustering
#   Table 23 [tab:placebos]              Placebo and falsification tests
#
# ── KEY NUMBERS (county 9m, two-way clustering) ───────────────────────────────
#   First-stage Wald F (9m, two-way) : 27.1
#   City first-stage Wald F (12m)    : 26.7
#   Instrument coefficient           : TBD (update after running)
#   Wu-Hausman p-value               : TBD
#   Shoppability gradient (V2)       : TBD
#   RI actual F / 95th pctile        : TBD
#   Lead joint F / p                 : TBD
#   LOO gradient range               : TBD
#
# IDENTIFICATION STRATEGY:
#   Primary instrument: system_peer_pressure_county_9m
#     For each county-month, counts distinct out-of-county system
#     peers (from systems PRESENT in the focal county) that posted
#     in the trailing 9-month window. Focal-county restriction
#     eliminates cross-system contamination that inflated the
#     original instrument. Two-way clustering: County_State + post_month.
#   County robustness: system_peer_pressure_county    (6m, F=9.6)
#                      system_peer_pressure_county_12m (12m, F=26.1)
#   City primary:      system_peer_pressure_city_12m  (12m, F=26.7)
#   City robustness:   system_peer_pressure_city_9m   (9m, F=20.0)
#                      system_peer_pressure_city       (6m, F=18.9)
#   CBSA dropped: instrument fails (F<10) under two-way clustering
#                 across all window lengths; 5,305 singleton removals.
#
# PRIMARY MARKET DEFINITION: County
#   County is the standard unit used by CMS, state regulators,
#   and antitrust authorities.
#
# FIXED EFFECTS: market_id (geo x service_group) + post_month
# CLUSTERING:    Two-way: geo-unit + post_month throughout
#
# COLOR SCHEME (FSU):
#   Garnet: #782F40  (primary estimates, significant)
#   Gold:   #CEB888  (reference, robustness, insignificant)
# ============================================================

# ----------------------------------------------------------------------------
# NAVIGATION GUIDE — added organization/comments only
# ----------------------------------------------------------------------------
# This copy preserves the executable code from the original script. The added 
# material is comment-only: section guides, function summaries, and output 
# labels.
# Quick search tags: SECTION 6, OUTPUT_TABLE_03, OUTPUT_FIGURE_01, 
# EXPORT_CSV_results_first_stage, ROBUSTNESS, PLACEBO, HETEROGENEITY.
# Primary spec: county panel, system_peer_pressure_county_9m, two-way 
# clustering by County_State + post_month.
# Use the OUTPUT_* labels below to find the exact code block that creates each
#  table, figure, or exported CSV.
#
# OUTPUT INDEX — TABLES / CONSOLE TABLES
#   OUTPUT_TABLE_01_SUMSTATS                Summary statistics (LaTeX; Section 20)
#   OUTPUT_TABLE_02_FIRST_STAGE             County first stage (etable; Section 6)
#   OUTPUT_TABLE_03_POOLED_OLS_IV           County pooled OLS vs IV, all outcomes
#   OUTPUT_TABLE_05B_SERVICE_OLS_IV         Service-level OLS vs IV comparison
#   OUTPUT_TABLE_06_METAREG                 Shoppability meta-regression
#   OUTPUT_TABLE_14_GEO_COMPARISON          County vs city geographic robustness
#   OUTPUT_TABLE_20_IV_VARIANTS             Instrument window/variant robustness
#   OUTPUT_TABLE_21_LOO                     Leave-one-system-out robustness
#   OUTPUT_TABLE_22_CLUSTERING              Two-way vs one-way clustering
#   OUTPUT_TABLE_23_PLACEBO_LEADS           Instrument lead placebo/pre-trends
#   OUTPUT_TABLE_R1_SYSTEM_FE_TRENDS        System-month FE and system trends
#   OUTPUT_TABLE_R2_ENFORCEMENT             Enforcement-control robustness
#   OUTPUT_TABLE_R3_ROBUSTNESS_SUMMARY      Exclusion-restriction robustness summary
#
# OUTPUT INDEX — FIGURES
#   OUTPUT_FIGURE_01_FIRST_STAGE_WALD       Saved as fig_first_stage_wald.pdf
#   OUTPUT_FIGURE_02_METAREG_COMBINED       Forest plot + meta-regression scatter
#   OUTPUT_FIGURE_03_OLS_IV_FOREST          OLS vs IV forest plot
#   OUTPUT_FIGURE_03B_OLS_IV_SCATTER        OLS vs IV scatter plot
#   OUTPUT_FIGURE_04_SATURATION             Saturation/bin effects
#   OUTPUT_FIGURE_05_SATURATION_ROBUST      Saturation robustness across schemes
#   OUTPUT_FIGURE_06_HETEROGENEITY_FOREST   Market heterogeneity forest plot
#   OUTPUT_FIGURE_07_RACE_MEDICAID_HEATMAP  Race x Medicaid gradient heatmap
#   OUTPUT_FIGURE_08_RANDOMIZATION          Randomization inference histogram
#   OUTPUT_FIGURE_PRETREND_LEADS            Instrument-lead pre-trend plot
#   OUTPUT_FIGURE_CHR_SENSITIVITY           Conley-Hansen-Rossi sensitivity plot
#
# OUTPUT INDEX — CSV EXPORTS
#   Search for EXPORT_CSV_ to jump to each saved results file.
# ============================================================
#
# SCRIPT SECTION NUMBERING (added for navigation):
#   Every main section below is marked with a bold, hash-bordered header:
#     ######## Section N: Title ###########
#   Each of these lines ends in 4+ "#" characters, which RStudio's editor
#   recognizes as a foldable section break — use Edit > Folding > Collapse
#   All, or the folding arrows in the gutter, to collapse the whole script
#   down to just this list of section titles.
#   Sections run 0-23. Two sections were added relative to earlier drafts
#   of this script: Section 18 (additional exclusion-restriction
#   diagnostics — system x month FE, Lee et al. tF, enforcement controls,
#   Conley-Hansen-Rossi bounds) and Sections 22-23 (within-system price
#   dispersion check and the payer-conditional robustness pipeline), which
#   pushed the former Sections 18-20 to 19-21.
# ============================================================


######## Section 0: Libraries, Setup, and Constants ###########
# ----------------------------------------------------------------------------
# Loads packages, defines FSU colors, creates helper constants, and defines 
# the service grouping function used by both county and city panels.

library(data.table)
library(fixest)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(forcats)
library(patchwork)
library(purrr)
library(lubridate)
library(geosphere)
library(stargazer)
library(knitr)
library(kableExtra)
library(hdm)
library(grf)
library(here)

# FSU color constants
FSU_GARNET <- "#782F40"
FSU_GOLD   <- "#CEB888"

# Star notation helper
# ----------------------------------------------------------------------------
# FUNCTION add_stars()
# ----------------------------------------------------------------------------
# Returns conventional significance stars from a p-value. Used only for 
# formatting printed/exported tables.
add_stars <- function(p) {
  case_when(
    is.na(p)  ~ "",
    p < 0.01  ~ "***",
    p < 0.05  ~ "**",
    p < 0.10  ~ "*",
    TRUE      ~ ""
  )
}

# Outcomes list
outcomes_list <- list(
  Median = "ln_median_price",
  Mean   = "ln_mean_price",
  P25    = "ln_p25_price",
  P75    = "ln_p75_price",
  Max    = "ln_max_price",
  Min    = "ln_min_price",
  IQR    = "ln_iqr_price"
)

# Saturation bin cutpoints
BIN1_LOW  <- 1;  BIN1_HIGH <- 3
BIN2_LOW  <- 4;  BIN2_HIGH <- 8

range_levels <- c(
  paste0(BIN1_LOW, "\u2013", BIN1_HIGH),
  paste0(BIN2_LOW, "\u2013", BIN2_HIGH),
  paste0(BIN2_HIGH + 1, "+")
)

# Service group classifications — applied to both geographies
# ----------------------------------------------------------------------------
# FUNCTION add_service_groups()
# ----------------------------------------------------------------------------
# Adds two service-category variables used repeatedly in service-group and 
# shoppability analyses.
add_service_groups <- function(df) {
  df %>%
    mutate(
      SERVICE_GROUP_SHOP = case_when(
        grepl("Ultrasound", SERVICE_GROUP)    ~ "Ultrasound",
        SERVICE_GROUP == "CT Lung"            ~ "CT Lung",
        SERVICE_GROUP == "Mammography"        ~ "Mammography",
        grepl("MRI",    SERVICE_GROUP)        ~ "MRI",
        grepl("CT",     SERVICE_GROUP)        ~ "CT Other",
        grepl("X-Ray",  SERVICE_GROUP)        ~ "X-Ray",
        grepl("Biopsy", SERVICE_GROUP)        ~ "Biopsy",
        SERVICE_GROUP %in%
          c("Colonoscopy", "Endoscopy")       ~ "Endoscopy/Colonoscopy",
        TRUE                                  ~ "Other"
      ),
      SERVICE_GROUP_SHOP2 = case_when(
        grepl("Ultrasound", SERVICE_GROUP) |
          SERVICE_GROUP == "Mammography"      ~ "Ultrasound/Mammography",
        SERVICE_GROUP == "CT Lung"            ~ "CT Lung",
        grepl("MRI",    SERVICE_GROUP)        ~ "MRI",
        grepl("CT",     SERVICE_GROUP)        ~ "CT Other",
        grepl("X-Ray",  SERVICE_GROUP)        ~ "X-Ray",
        grepl("Biopsy", SERVICE_GROUP)        ~ "Biopsy",
        SERVICE_GROUP %in%
          c("Colonoscopy", "Endoscopy")       ~ "Endoscopy/Colonoscopy",
        TRUE                                  ~ "Other"
      )
    )
}


######## Section 1: Load and Prepare County Data (Primary) ###########
# ----------------------------------------------------------------------------
# Reads the county-level cleaned national price file and constructs the 
# primary county entrant-month analysis panel with log price outcomes.
setwd("~/Library/CloudStorage/OneDrive-FloridaStateUniversity/Hospital Price Transparency Paper/Data")
df_county_raw <- fread("../Data/data_clean/National_prices_with_controls_county_final.csv")
cat("County rows:", nrow(df_county_raw), "| Columns:", ncol(df_county_raw), "\n")

df_iv_county <- df_county_raw %>%
  filter(IS_ENTRANT_MONTH_COUNTY == 1) %>%
  mutate(
    post_month      = as.Date(post_month),
    cbsacode        = as.integer(CBSA_CODE),
    County_State    = County_State,
    market_id       = paste0(County_State, "::", SERVICE_GROUP),
    n_prior_posters = as.numeric(N_PRIOR_POSTERS_VERIFIED),
    ln_median_price = log(MEDIAN_PRICE + 1),
    ln_mean_price   = log(MEAN_PRICE   + 1),
    ln_iqr_price    = log(IQR_PRICE    + 1),
    ln_min_price    = log(MIN_PRICE    + 1),
    ln_max_price    = log(MAX_PRICE    + 1),
    ln_p25_price    = log(P25          + 1),
    ln_p75_price    = log(P75          + 1),
    ln_total_beds   = log(TOTAL_BEDS   + 1)
  ) %>%
  filter(cbsacode < 50000) %>%
  add_service_groups()

cat("County: rows =",   nrow(df_iv_county),
    "| counties =",     uniqueN(df_iv_county$County_State),
    "| markets =",      uniqueN(df_iv_county$market_id), "\n")
print(summary(df_iv_county$n_prior_posters))


######## Section 2: County Instrument — Corrected + All Variants ###########
# ----------------------------------------------------------------------------
# KEY FIX applied throughout: only systems WITH hospitals in the
# focal county can exert peer pressure on that county.
# Constructs corrected county instruments for 3m, 6m, 9m, and 12m trailing 
# windows using only health systems present in the focal county.
# This eliminates cross-system contamination

# Step 1: Systems spanning >1 county
system_geo_county <- df_iv_county %>%
  distinct(HOSPITAL_ID, County_State, HEALTH_SYSTEM_NAME) %>%
  filter(!is.na(HEALTH_SYSTEM_NAME), HEALTH_SYSTEM_NAME != "") %>%
  group_by(HEALTH_SYSTEM_NAME) %>%
  summarise(n_counties = n_distinct(County_State), .groups = "drop") %>%
  filter(n_counties > 1)
cat("Systems spanning >1 county:", nrow(system_geo_county), "\n")

# Step 2: Hospital-system-county-month panel
hosp_system_county <- df_iv_county %>%
  distinct(HOSPITAL_ID, County_State, post_month, HEALTH_SYSTEM_NAME) %>%
  filter(!is.na(HEALTH_SYSTEM_NAME), HEALTH_SYSTEM_NAME != "") %>%
  filter(HEALTH_SYSTEM_NAME %in% system_geo_county$HEALTH_SYSTEM_NAME)

# Step 3: KEY FIX — which systems are present in each focal county
focal_county_systems <- df_iv_county %>%
  distinct(County_State, HEALTH_SYSTEM_NAME) %>%
  filter(!is.na(HEALTH_SYSTEM_NAME), HEALTH_SYSTEM_NAME != "") %>%
  filter(HEALTH_SYSTEM_NAME %in% system_geo_county$HEALTH_SYSTEM_NAME)
cat("Focal county-system pairs:", nrow(focal_county_systems), "\n")

# Step 4: Panel of all county-months
panel_county_months <- df_iv_county %>% distinct(County_State, post_month)

# Step 5: Build function (supports window variants)
# ----------------------------------------------------------------------------
# FUNCTION build_county_instrument()
# ----------------------------------------------------------------------------
# Builds the corrected county peer-pressure instrument for a user-specified 
# trailing window; code logic is unchanged.
build_county_instrument <- function(hosp_data, panel_months, focal_systems,
                                    window_months = 6,
                                    var_suffix    = "") {
  hosp_data %>%
    rename(peer_hosp       = HOSPITAL_ID,
           peer_county     = County_State,
           peer_post_month = post_month) %>%
    cross_join(panel_months) %>%
    filter(
      peer_county     != County_State,
      peer_post_month <= post_month,
      peer_post_month >= post_month - months(window_months)
    ) %>%
    inner_join(focal_systems, by = c("HEALTH_SYSTEM_NAME", "County_State")) %>%
    group_by(HEALTH_SYSTEM_NAME, County_State, post_month) %>%
    summarise(n_peers = n_distinct(peer_hosp), .groups = "drop") %>%
    group_by(County_State, post_month) %>%
    summarise(
      !!paste0("system_peer_pressure_county", var_suffix) := sum(n_peers),
      .groups = "drop"
    )
}

cat("Building county instruments...\n")

county_iv_3m  <- build_county_instrument(hosp_system_county, panel_county_months,
                                         focal_county_systems, window_months = 3,
                                         var_suffix = "_3m")   # <-- this was missing

county_iv_6m  <- build_county_instrument(hosp_system_county, panel_county_months,
                                         focal_county_systems, window_months = 6)
county_iv_9m  <- build_county_instrument(hosp_system_county, panel_county_months,
                                         focal_county_systems, window_months = 9,
                                         var_suffix = "_9m")
county_iv_12m <- build_county_instrument(hosp_system_county, panel_county_months,
                                         focal_county_systems, window_months = 12,
                                         var_suffix = "_12m")

df_iv_county <- df_iv_county %>%
  select(-any_of(c("system_peer_pressure_county",
                   "system_peer_pressure_county_3m",
                   "system_peer_pressure_county_9m",
                   "system_peer_pressure_county_12m"))) %>%
  left_join(county_iv_6m,  by = c("County_State", "post_month")) %>%
  left_join(county_iv_3m,  by = c("County_State", "post_month")) %>%
  left_join(county_iv_9m,  by = c("County_State", "post_month")) %>%
  left_join(county_iv_12m, by = c("County_State", "post_month")) %>%
  mutate(across(starts_with("system_peer_pressure_county"), ~replace_na(., 0)))

cat("County instrument summaries:\n")
cat("6m:  "); print(summary(df_iv_county$system_peer_pressure_county))
cat("9m:  "); print(summary(df_iv_county$system_peer_pressure_county_9m))
cat("3m:  "); print(summary(df_iv_county$system_peer_pressure_county_3m))
cat("12m: "); print(summary(df_iv_county$system_peer_pressure_county_12m))

# Confirm first-stage strength across variants
cat("\nCounty first-stage F-stats (two-way clustering):\n")
ivs_county <- c("system_peer_pressure_county",
                "system_peer_pressure_county_3m",
                "system_peer_pressure_county_9m",
                "system_peer_pressure_county_12m")
for (iv in ivs_county) {
  fs <- feols(
    as.formula(paste("n_prior_posters ~", iv,
                     "+ ln_total_beds | market_id + post_month")),
    data    = df_iv_county,
    cluster = ~County_State + post_month
  )
  t_stat <- summary(fs)$coeftable[1, "t value"]
  cat(sprintf("  %-45s  t = %5.2f  F = %6.2f\n", iv, t_stat, t_stat^2))
}


######## Section 3: Load and Prepare City Data (Robustness) ###########
# ----------------------------------------------------------------------------
# Reads the city-level cleaned national price file and constructs the city 
# entrant-month robustness panel.
df_city_raw <- fread("../Data/data_clean/National_prices_with_controls_city_final.csv")

cat("City rows:", nrow(df_city_raw), "| Columns:", ncol(df_city_raw), "\n")

df_iv_city <- df_city_raw %>%
  filter(IS_ENTRANT_MONTH_CITY == 1) %>%
  mutate(
    post_month      = as.Date(post_month),
    cbsacode        = as.integer(CBSA_CODE),
    city_state      = city_state,
    market_id       = paste0(city_state, "::", SERVICE_GROUP),
    n_prior_posters = as.numeric(N_HOSPITALS_CUMULATIVE_CITY - 1),
    ln_median_price = log(MEDIAN_PRICE + 1),
    ln_mean_price   = log(MEAN_PRICE   + 1),
    ln_iqr_price    = log(IQR_PRICE    + 1),
    ln_min_price    = log(MIN_PRICE    + 1),
    ln_max_price    = log(MAX_PRICE    + 1),
    ln_p25_price    = log(P25          + 1),
    ln_p75_price    = log(P75          + 1),
    ln_total_beds   = log(TOTAL_BEDS   + 1)
  ) %>%
  filter(cbsacode < 50000) %>%
  add_service_groups()

cat("City: rows =",  nrow(df_iv_city),
    "| cities =",    uniqueN(df_iv_city$city_state),
    "| markets =",   uniqueN(df_iv_city$market_id), "\n")
print(summary(df_iv_city$n_prior_posters))


######## Section 4: City Instrument — Corrected + All Variants ###########
# ----------------------------------------------------------------------------
# Constructs corrected city instruments for 6m, 9m, and 12m trailing windows 
# using only health systems present in the focal city.

# Step 1: Systems spanning >1 city
system_geo_city <- df_iv_city %>%
  distinct(HOSPITAL_ID, city_state, HEALTH_SYSTEM_NAME) %>%
  filter(!is.na(HEALTH_SYSTEM_NAME), HEALTH_SYSTEM_NAME != "") %>%
  group_by(HEALTH_SYSTEM_NAME) %>%
  summarise(n_cities = n_distinct(city_state), .groups = "drop") %>%
  filter(n_cities > 1)
cat("City: Systems spanning >1 city:", nrow(system_geo_city), "\n")

# Step 2: Hospital-system-city-month panel
hosp_system_city <- df_iv_city %>%
  distinct(HOSPITAL_ID, city_state, post_month, HEALTH_SYSTEM_NAME) %>%
  filter(!is.na(HEALTH_SYSTEM_NAME), HEALTH_SYSTEM_NAME != "") %>%
  filter(HEALTH_SYSTEM_NAME %in% system_geo_city$HEALTH_SYSTEM_NAME)

# Step 3: KEY FIX — which systems are present in each focal city
focal_city_systems <- df_iv_city %>%
  distinct(city_state, HEALTH_SYSTEM_NAME) %>%
  filter(!is.na(HEALTH_SYSTEM_NAME), HEALTH_SYSTEM_NAME != "") %>%
  filter(HEALTH_SYSTEM_NAME %in% system_geo_city$HEALTH_SYSTEM_NAME)
cat("Focal city-system pairs:", nrow(focal_city_systems), "\n")

# Step 4: Panel of all city-months
panel_city_months <- df_iv_city %>% distinct(city_state, post_month)

# Step 5: Build function
# ----------------------------------------------------------------------------
# FUNCTION build_city_instrument()
# ----------------------------------------------------------------------------
# Builds the corrected city peer-pressure instrument for a user-specified 
# trailing window; code logic is unchanged.
build_city_instrument <- function(hosp_data, panel_months, focal_systems,
                                  window_months = 6,
                                  var_suffix    = "") {
  hosp_data %>%
    rename(peer_hosp       = HOSPITAL_ID,
           peer_city       = city_state,
           peer_post_month = post_month) %>%
    cross_join(panel_months) %>%
    filter(
      peer_city       != city_state,
      peer_post_month <= post_month,
      peer_post_month >= post_month - months(window_months)
    ) %>%
    inner_join(focal_systems, by = c("HEALTH_SYSTEM_NAME", "city_state")) %>%
    group_by(HEALTH_SYSTEM_NAME, city_state, post_month) %>%
    summarise(n_peers = n_distinct(peer_hosp), .groups = "drop") %>%
    group_by(city_state, post_month) %>%
    summarise(
      !!paste0("system_peer_pressure_city", var_suffix) := sum(n_peers),
      .groups = "drop"
    )
}

cat("Building city instruments...\n")

city_iv_6m  <- build_city_instrument(hosp_system_city, panel_city_months,
                                     focal_city_systems, window_months = 6)
city_iv_9m  <- build_city_instrument(hosp_system_city, panel_city_months,
                                     focal_city_systems, window_months = 9,
                                     var_suffix = "_9m")
city_iv_12m <- build_city_instrument(hosp_system_city, panel_city_months,
                                     focal_city_systems, window_months = 12,
                                     var_suffix = "_12m")

df_iv_city <- df_iv_city %>%
  select(-any_of(c("system_peer_pressure_city",
                   "system_peer_pressure_city_9m",
                   "system_peer_pressure_city_12m"))) %>%
  left_join(city_iv_6m,  by = c("city_state", "post_month")) %>%
  left_join(city_iv_9m,  by = c("city_state", "post_month")) %>%
  left_join(city_iv_12m, by = c("city_state", "post_month")) %>%
  mutate(across(starts_with("system_peer_pressure_city"), ~replace_na(., 0)))

cat("City instrument summaries:\n")
cat("6m:  "); print(summary(df_iv_city$system_peer_pressure_city))
cat("9m:  "); print(summary(df_iv_city$system_peer_pressure_city_9m))
cat("12m: "); print(summary(df_iv_city$system_peer_pressure_city_12m))

# Confirm first-stage strength across variants
cat("\nCity first-stage F-stats (two-way clustering):\n")
ivs_city <- c("system_peer_pressure_city",
              "system_peer_pressure_city_9m",
              "system_peer_pressure_city_12m")
for (iv in ivs_city) {
  fs <- feols(
    as.formula(paste("n_prior_posters ~", iv,
                     "+ ln_total_beds | market_id + post_month")),
    data    = df_iv_city,
    cluster = ~city_state + post_month
  )
  t_stat <- summary(fs)$coeftable[1, "t value"]
  cat(sprintf("  %-45s  t = %5.2f  F = %6.2f\n", iv, t_stat, t_stat^2))
}


######## Section 5: Shared Estimation Functions ###########
# ----------------------------------------------------------------------------
# Reusable estimation and table-building helpers. These functions run pooled 
# OLS/IV, service-level IV, outcome-by-service IV, and meta-regressions.
# cluster_var accepts a character vector for two-way clustering,
# e.g. c("County_State", "post_month") → ~County_State + post_month

# ----------------------------------------------------------------------------
# FUNCTION extract_first_stage()
# ----------------------------------------------------------------------------
# Extracts first-stage coefficient, standard error, t-stat, and IV Wald F from
#  a fixest IV object.
extract_first_stage <- function(fit, instrument) {
  fs_obj <- tryCatch(fit$iv_first_stage[["n_prior_posters"]], error = function(e) NULL)
  if (is.null(fs_obj)) return(list(fs_coef = NA_real_, fs_se = NA_real_,
                                   fs_t = NA_real_,  fs_f  = NA_real_))
  ct   <- tryCatch(fs_obj$coeftable, error = function(e) NULL)
  coef <- tryCatch(ct[instrument, "Estimate"],   error = function(e) NA_real_)
  se_v <- tryCatch(ct[instrument, "Std. Error"], error = function(e) NA_real_)
  t_v  <- tryCatch(ct[instrument, "t value"],    error = function(e) NA_real_)
  f_v  <- tryCatch(fitstat(fit, "ivwald")[["ivwald1::n_prior_posters"]]$stat,
                   error = function(e) NA_real_)
  list(fs_coef = coef, fs_se = se_v, fs_t = t_v, fs_f = f_v)
}

# ----------------------------------------------------------------------------
# FUNCTION print_iv_diagnostics()
# ----------------------------------------------------------------------------
# Prints a compact diagnostic summary for a pooled IV model: first stage, 
# second stage, p-value, and 95% CI.
print_iv_diagnostics <- function(fit, instrument, label = "") {
  fs <- extract_first_stage(fit, instrument)
  est  <- tryCatch(coef(fit)["fit_n_prior_posters"],   error = function(e) NA_real_)
  se_v <- tryCatch(se(fit)["fit_n_prior_posters"],     error = function(e) NA_real_)
  pval <- tryCatch(pvalue(fit)["fit_n_prior_posters"], error = function(e) NA_real_)
  cat("\n", strrep("=", 60), "\n")
  if (nchar(label) > 0) cat(" ", label, "\n", strrep("-", 60), "\n")
  cat(sprintf("  FIRST STAGE\n"))
  cat(sprintf("    Instrument coef : %.4f\n", fs$fs_coef))
  cat(sprintf("    Instrument SE   : %.4f\n", fs$fs_se))
  cat(sprintf("    Instrument t    : %.2f\n",  fs$fs_t))
  cat(sprintf("    Wald F          : %.1f\n",  fs$fs_f))
  cat(sprintf("  SECOND STAGE\n"))
  cat(sprintf("    n_prior_posters : %.4f\n", est))
  cat(sprintf("    SE              : %.4f\n", se_v))
  cat(sprintf("    p-value         : %.4f\n", pval))
  cat(sprintf("    95%% CI          : [%.4f, %.4f]\n",
              est - 1.96*se_v, est + 1.96*se_v))
  cat(strrep("=", 60), "\n")
}

# ----------------------------------------------------------------------------
# FUNCTION run_ols_pooled()
# ----------------------------------------------------------------------------
# Runs the pooled OLS specification for one outcome with market and month 
# fixed effects.
run_ols_pooled <- function(data, outcome_var = "ln_median_price", cluster_var) {
  feols(
    as.formula(paste0(outcome_var,
                      " ~ n_prior_posters + ln_total_beds | market_id + post_month")),
    data    = data,
    cluster = as.formula(paste0("~", paste(cluster_var, collapse = " + ")))
  )
}

# ----------------------------------------------------------------------------
# FUNCTION run_ols_by_service()
# ----------------------------------------------------------------------------
# Runs the OLS specification separately by service group and returns estimates
#  in percentage-point form.
run_ols_by_service <- function(data, outcome_var = "ln_median_price",
                               cluster_var, min_obs = 100) {
  service_groups <- unique(data$SERVICE_GROUP)
  results <- do.call(rbind, lapply(service_groups, function(sg) {
    df_sub <- data[data$SERVICE_GROUP == sg, ]
    if (nrow(df_sub) < min_obs) return(NULL)
    fit <- tryCatch(
      feols(as.formula(paste0(outcome_var,
                              " ~ n_prior_posters + ln_total_beds | market_id + post_month")),
            data    = df_sub,
            cluster = as.formula(paste0("~", paste(cluster_var, collapse = " + ")))),
      error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    est  <- tryCatch(coef(fit)["n_prior_posters"],   error = function(e) NA_real_)
    se_v <- tryCatch(se(fit)["n_prior_posters"],     error = function(e) NA_real_)
    pval <- tryCatch(pvalue(fit)["n_prior_posters"], error = function(e) NA_real_)
    data.frame(SERVICE_GROUP = sg, n_obs = nrow(df_sub),
               estimate = est, se = se_v, pval = pval,
               estimate_pct = est * 100, se_pct = se_v * 100,
               ci_lo_pct = (est - 1.96 * se_v) * 100,
               ci_hi_pct = (est + 1.96 * se_v) * 100,
               significant = !is.na(est) & !is.na(se_v) &
                 sign((est - 1.96*se_v)) == sign((est + 1.96*se_v)),
               row.names = NULL)
  }))
  results[order(results$estimate_pct), ]
}

# ----------------------------------------------------------------------------
# FUNCTION run_iv_pooled()
# ----------------------------------------------------------------------------
# Runs the pooled IV specification for one outcome and instrument.
run_iv_pooled <- function(data, outcome_var = "ln_median_price",
                          instrument, cluster_var, verbose = FALSE) {
  fit <- feols(
    as.formula(paste0(outcome_var,
                      " ~ ln_total_beds | market_id + post_month | n_prior_posters ~ ",
                      instrument)),
    data    = data,
    cluster = as.formula(paste0("~", paste(cluster_var, collapse = " + ")))
  )
  if (verbose) print_iv_diagnostics(fit, instrument,
                                    label = paste("IV Pooled |", outcome_var))
  fit
}

# ----------------------------------------------------------------------------
# FUNCTION run_iv_by_service()
# ----------------------------------------------------------------------------
# Runs the IV specification separately by service group and returns estimates 
# plus first-stage diagnostics.
run_iv_by_service <- function(data, outcome_var = "ln_median_price",
                              instrument, cluster_var, min_obs = 100) {
  service_groups <- unique(data$SERVICE_GROUP)
  results <- do.call(rbind, lapply(service_groups, function(sg) {
    df_sub <- data[data$SERVICE_GROUP == sg, ]
    if (nrow(df_sub) < min_obs) return(NULL)
    fit <- tryCatch(
      feols(as.formula(paste0(outcome_var,
                              " ~ ln_total_beds | market_id + post_month | n_prior_posters ~ ",
                              instrument)),
            data    = df_sub,
            cluster = as.formula(paste0("~", paste(cluster_var, collapse = " + ")))),
      error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    fs   <- extract_first_stage(fit, instrument)
    est  <- tryCatch(coef(fit)["fit_n_prior_posters"],   error = function(e) NA_real_)
    se_v <- tryCatch(se(fit)["fit_n_prior_posters"],     error = function(e) NA_real_)
    pval <- tryCatch(pvalue(fit)["fit_n_prior_posters"], error = function(e) NA_real_)
    data.frame(SERVICE_GROUP = sg, n_obs = nrow(df_sub),
               estimate = est, se = se_v, pval = pval,
               estimate_pct = est * 100, se_pct = se_v * 100,
               ci_lo = est - 1.96 * se_v, ci_hi = est + 1.96 * se_v,
               ci_lo_pct = (est - 1.96 * se_v) * 100,
               ci_hi_pct = (est + 1.96 * se_v) * 100,
               significant = !is.na(est) & !is.na(se_v) &
                 sign(est - 1.96*se_v) == sign(est + 1.96*se_v),
               fs_coef = fs$fs_coef, fs_se = fs$fs_se,
               fs_t = fs$fs_t, fs_f = fs$fs_f,
               row.names = NULL)
  }))
  results <- results[order(results$estimate_pct), ]
  results$SERVICE_GROUP <- factor(results$SERVICE_GROUP, levels = results$SERVICE_GROUP)
  results
}

# ----------------------------------------------------------------------------
# FUNCTION run_iv_by_large_group()
# ----------------------------------------------------------------------------
# Runs IV estimates by a larger grouping variable such as SERVICE_GROUP_SHOP2.
run_iv_by_large_group <- function(data, group_var,
                                  outcome_var = "ln_median_price",
                                  instrument, cluster_var, min_obs = 100) {
  groups <- unique(data[[group_var]])
  do.call(rbind, lapply(groups, function(grp) {
    df_sub <- data[data[[group_var]] == grp, ]
    if (nrow(df_sub) < min_obs) return(NULL)
    fit <- tryCatch(
      feols(as.formula(paste0(outcome_var,
                              " ~ ln_total_beds | market_id + post_month | n_prior_posters ~ ",
                              instrument)),
            data    = df_sub,
            cluster = as.formula(paste0("~", paste(cluster_var, collapse = " + ")))),
      error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    fs   <- extract_first_stage(fit, instrument)
    est  <- tryCatch(coef(fit)["fit_n_prior_posters"],   error = function(e) NA_real_)
    se_v <- tryCatch(se(fit)["fit_n_prior_posters"],     error = function(e) NA_real_)
    pval <- tryCatch(pvalue(fit)["fit_n_prior_posters"], error = function(e) NA_real_)
    data.frame(group = grp, n_obs = nrow(df_sub),
               estimate_pct = est * 100, se_pct = se_v * 100, pval = pval,
               ci_lo_pct = (est - 1.96 * se_v) * 100,
               ci_hi_pct = (est + 1.96 * se_v) * 100,
               significant = !is.na(est) & !is.na(se_v) &
                 sign((est - 1.96*se_v)) == sign((est + 1.96*se_v)),
               fs_coef = fs$fs_coef, fs_se = fs$fs_se,
               fs_t = fs$fs_t, fs_f = fs$fs_f,
               row.names = NULL)
  }))
}

# ----------------------------------------------------------------------------
# FUNCTION run_iv_by_service_outcomes()
# ----------------------------------------------------------------------------
# Runs service-level IV estimates across the full outcome list: mean, median, 
# percentiles, min, max, and IQR.
run_iv_by_service_outcomes <- function(data, outcomes, instrument,
                                       cluster_var, cum_col, min_obs = 100) {
  service_groups <- unique(data$SERVICE_GROUP)
  do.call(rbind, lapply(service_groups, function(sg) {
    df_sub <- data[data$SERVICE_GROUP == sg, ]
    if (nrow(df_sub) < min_obs) return(NULL)
    do.call(rbind, lapply(names(outcomes), function(out_nm) {
      dep_var <- outcomes[[out_nm]]
      df_use  <- df_sub
      if (dep_var == "ln_iqr_price") {
        df_use <- df_use[df_use[[cum_col]] > 1 & is.finite(df_use$ln_iqr_price), ]
      }
      if (nrow(df_use) < min_obs) return(NULL)
      fit <- tryCatch(
        feols(as.formula(paste0(dep_var,
                                " ~ ln_total_beds | market_id + post_month | n_prior_posters ~ ",
                                instrument)),
              data    = df_use,
              cluster = as.formula(paste0("~", paste(cluster_var, collapse = " + ")))),
        error = function(e) NULL)
      if (is.null(fit)) return(NULL)
      fs   <- extract_first_stage(fit, instrument)
      est  <- tryCatch(coef(fit)["fit_n_prior_posters"],   error = function(e) NA_real_)
      se_v <- tryCatch(se(fit)["fit_n_prior_posters"],     error = function(e) NA_real_)
      pval <- tryCatch(pvalue(fit)["fit_n_prior_posters"], error = function(e) NA_real_)
      data.frame(SERVICE_GROUP = sg, outcome = out_nm,
                 estimate = est, se = se_v, pval = pval,
                 estimate_pct = est * 100, se_pct = se_v * 100,
                 ci_lo_pct = (est - 1.96 * se_v) * 100,
                 ci_hi_pct = (est + 1.96 * se_v) * 100,
                 significant = !is.na(est) & !is.na(se_v) &
                   sign(est - 1.96*se_v) == sign(est + 1.96*se_v),
                 fs_coef = fs$fs_coef, fs_se = fs$fs_se,
                 fs_t = fs$fs_t, fs_f = fs$fs_f,
                 row.names = NULL)
    }))
  }))
}

# ----------------------------------------------------------------------------
# FUNCTION build_meta_data()
# ----------------------------------------------------------------------------
# Converts service-level results into a meta-regression dataset with 
# shoppability indicators and inverse-variance weights.
build_meta_data <- function(res_service) {
  res_service %>%
    mutate(
      SERVICE_GROUP  = as.character(SERVICE_GROUP),
      is_ultrasound  = grepl("Ultrasound", SERVICE_GROUP),
      is_mri         = grepl("MRI",        SERVICE_GROUP),
      is_xray        = grepl("X-Ray",      SERVICE_GROUP),
      is_biopsy      = grepl("Biopsy",     SERVICE_GROUP),
      is_ct_lung     = SERVICE_GROUP == "CT Lung",
      is_mammography = SERVICE_GROUP == "Mammography",
      is_shoppable   = is_ultrasound | is_ct_lung | is_mammography,
      fs_f   = if ("fs_f" %in% names(.)) fs_f else NA_real_,
      weight = ifelse(is.na(se_pct) | se_pct == 0, NA_real_, 1 / (se_pct^2))
    ) %>%
    filter(is.finite(weight))
}

# ----------------------------------------------------------------------------
# FUNCTION run_meta_regressions()
# ----------------------------------------------------------------------------
# Runs WLS and unweighted meta-regressions used to estimate the shoppability 
# gradient.
run_meta_regressions <- function(meta_data) {
  m_simple    <- lm(estimate_pct ~ is_shoppable,
                    data = meta_data, weights = weight)
  m_granular  <- lm(estimate_pct ~ is_ultrasound + is_mri + is_xray +
                      is_biopsy + is_ct_lung + is_mammography,
                    data = meta_data, weights = weight)
  m_simple_uw   <- lm(estimate_pct ~ is_shoppable,   data = meta_data)
  m_granular_uw <- lm(estimate_pct ~ is_ultrasound + is_mri + is_xray +
                        is_biopsy + is_ct_lung + is_mammography, data = meta_data)
  m_fstat <- if ("fs_f" %in% names(meta_data) && sum(!is.na(meta_data$fs_f)) > 3) {
    lm(estimate_pct ~ is_shoppable + fs_f, data = meta_data, weights = weight)
  } else {
    message("  Skipping F-stat sensitivity model (fs_f all NA or too few obs)")
    NULL
  }
  list(simple = m_simple, granular = m_granular,
       simple_uw = m_simple_uw, granular_uw = m_granular_uw,
       fstat_check = m_fstat)
}

# ----------------------------------------------------------------------------
# FUNCTION print_meta_results()
# ----------------------------------------------------------------------------
# Prints the meta-regression model summaries in a consistent order.
print_meta_results <- function(meta_list, label = "") {
  cat("\n", strrep("=", 60), "\n")
  if (nchar(label) > 0) cat(" META-REGRESSION:", label, "\n")
  cat("\n--- (1) Simple shoppable flag (WLS) ---\n");    print(summary(meta_list$simple))
  cat("\n--- (2) Granular service-type dummies (WLS) ---\n"); print(summary(meta_list$granular))
  cat("\n--- (3) Simple shoppable flag (Unweighted) ---\n");  print(summary(meta_list$simple_uw))
  cat("\n--- (4) Granular service-type dummies (Unweighted) ---\n"); print(summary(meta_list$granular_uw))
  if (!is.null(meta_list$fstat_check)) {
    cat("\n--- (5) Shoppable + first-stage F ---\n"); print(summary(meta_list$fstat_check))
  }
  cat(strrep("=", 60), "\n")
}


######## Section 6: County — Primary Results ###########
# ----------------------------------------------------------------------------
# Primary instrument: system_peer_pressure_county_9m
# Clustering: County_State + post_month (two-way)
# Main county results: first stage, Wu-Hausman test, pooled OLS/IV outcomes, 
# service-level IV results, OLS-vs-IV comparison, and county meta-regressions.
# ---------------------------------------------------------------------------
# First stage — stand-alone (two-way clustering)
# ---------------------------------------------------------------------------
fs_county_ctrl <- feols(
  n_prior_posters ~ system_peer_pressure_county_9m + ln_total_beds |
    market_id + post_month,
  data    = df_iv_county,
  cluster = ~County_State + post_month
)

cat("\n=== COUNTY FIRST STAGE (PRIMARY: 9m, two-way clustering) ===\n")
summary(fs_county_ctrl)
cat(sprintf("First-stage t-stat : %.2f\n",
            summary(fs_county_ctrl)$coeftable[1, "t value"]))
cat(sprintf("Implied Wald F     : %.1f\n",
            summary(fs_county_ctrl)$coeftable[1, "t value"]^2))

# Wu-Hausman test
iv_county_ctrl <- feols(
  ln_median_price ~ 1 + ln_total_beds |
    market_id + post_month |
    n_prior_posters ~ system_peer_pressure_county_9m,
  data    = df_iv_county,
  cluster = ~County_State + post_month
)
wh_pval <- fitstat(iv_county_ctrl, "wh")$wh$p
cat(sprintf("Wu-Hausman p-value : %.4e\n", wh_pval))

# First-stage etable
# ----------------------------------------------------------------------------
# OUTPUT_TABLE_02_FIRST_STAGE
# ----------------------------------------------------------------------------
# First-stage estimates for the primary county panel; printed through 
# etable().
etable(
  fs_county_ctrl,
  coefstat   = "se",
  keep       = c("%system_peer_pressure_county_9m", "%ln_total_beds"),
  dict       = c(system_peer_pressure_county_9m = "System Peer Pressure (9m)",
                 ln_total_beds = "ln(Total Beds)",
                 n_prior_posters = "Prior Disclosures (N)"),
  fitstat    = ~ n + r2 + wr2 + f.stat,
  se.below   = TRUE, digits = 4, digits.stats = 3, tex = TRUE
)

# ---------------------------------------------------------------------------
# OLS: pooled — all outcomes
# ---------------------------------------------------------------------------
ols_county_mean   <- run_ols_pooled(df_iv_county, "ln_mean_price",
                                    c("County_State", "post_month"))
ols_county_p25    <- run_ols_pooled(df_iv_county, "ln_p25_price",
                                    c("County_State", "post_month"))
ols_county_median <- run_ols_pooled(df_iv_county, "ln_median_price",
                                    c("County_State", "post_month"))
ols_county_p75    <- run_ols_pooled(df_iv_county, "ln_p75_price",
                                    c("County_State", "post_month"))
ols_county_max    <- run_ols_pooled(df_iv_county, "ln_max_price",
                                    c("County_State", "post_month"))
ols_county_min    <- run_ols_pooled(df_iv_county, "ln_min_price",
                                    c("County_State", "post_month"))
ols_county_iqr    <- feols(
  ln_iqr_price ~ n_prior_posters + ln_total_beds | market_id + post_month,
  data    = df_iv_county %>% filter(N_HOSPITALS_CUMULATIVE_COUNTY > 1,
                                    is.finite(ln_iqr_price)),
  cluster = ~County_State + post_month
)

# ---------------------------------------------------------------------------
# IV: pooled — all outcomes (primary: 9m, two-way clustering)
# ---------------------------------------------------------------------------
iv_county_mean   <- run_iv_pooled(df_iv_county, "ln_mean_price",
                                  "system_peer_pressure_county_9m",
                                  c("County_State", "post_month"), verbose = TRUE)
iv_county_p25    <- run_iv_pooled(df_iv_county, "ln_p25_price",
                                  "system_peer_pressure_county_9m",
                                  c("County_State", "post_month"), verbose = TRUE)
iv_county_pooled <- run_iv_pooled(df_iv_county, "ln_median_price",
                                  "system_peer_pressure_county_9m",
                                  c("County_State", "post_month"), verbose = TRUE)
iv_county_p75    <- run_iv_pooled(df_iv_county, "ln_p75_price",
                                  "system_peer_pressure_county_9m",
                                  c("County_State", "post_month"), verbose = TRUE)
iv_county_max    <- run_iv_pooled(df_iv_county, "ln_max_price",
                                  "system_peer_pressure_county_9m",
                                  c("County_State", "post_month"), verbose = TRUE)
iv_county_min    <- run_iv_pooled(df_iv_county, "ln_min_price",
                                  "system_peer_pressure_county_9m",
                                  c("County_State", "post_month"), verbose = TRUE)
iv_county_iqr    <- feols(
  ln_iqr_price ~ ln_total_beds | market_id + post_month |
    n_prior_posters ~ system_peer_pressure_county_9m,
  data    = df_iv_county %>% filter(N_HOSPITALS_CUMULATIVE_COUNTY > 1,
                                    is.finite(ln_iqr_price)),
  cluster = ~County_State + post_month
)

cat("\n=== IQR IV — First-Stage Wald F ===\n")
cat(sprintf("  Wald F : %.1f\n",
            fitstat(iv_county_iqr, "ivwald")[["ivwald1::n_prior_posters"]]$stat))

# Clustering robustness: two-way (default) vs one-way
iv_county_1way <- feols(
  ln_median_price ~ ln_total_beds | market_id + post_month |
    n_prior_posters ~ system_peer_pressure_county_9m,
  data    = df_iv_county,
  cluster = ~County_State
)

# ----------------------------------------------------------------------------
# OUTPUT_TABLE_22_CLUSTERING
# ----------------------------------------------------------------------------
# Clustering robustness: county two-way clustering versus county-only 
# clustering.
cat("\n=== TABLE 22: Clustering Robustness — County ===\n")
etable(iv_county_pooled, iv_county_1way,
       keep    = "%fit_n_prior_posters",
       headers = c("Two-way cluster (default)", "County-only cluster"))
cat(sprintf("  Two-way Wald F : %.1f\n",
            fitstat(iv_county_pooled, "ivwald")[["ivwald1::n_prior_posters"]]$stat))
cat(sprintf("  One-way Wald F : %.1f\n",
            fitstat(iv_county_1way,   "ivwald")[["ivwald1::n_prior_posters"]]$stat))

# Table: Pooled OLS vs IV — all outcomes
# ----------------------------------------------------------------------------
# OUTPUT_TABLE_03_POOLED_OLS_IV
# ----------------------------------------------------------------------------
# Pooled OLS vs IV estimates across all price outcomes for the county panel.
cat("\n=== TABLE 3: Pooled OLS vs IV — County — All Outcomes ===\n")
etable(
  ols_county_mean,   iv_county_mean,
  ols_county_p25,    iv_county_p25,
  ols_county_median, iv_county_pooled,
  ols_county_p75,    iv_county_p75,
  ols_county_max,    iv_county_max,
  ols_county_min,    iv_county_min,
  ols_county_iqr,    iv_county_iqr,
  keep    = "%n_prior_posters|%fit_n_prior_posters",
  headers = c("OLS","IV","OLS","IV","OLS","IV",
              "OLS","IV","OLS","IV","OLS","IV","OLS","IV")
)

# Wald F across all pooled IV specs
cat("\n=== POOLED IV — Wald F Summary (all outcomes) ===\n")
pooled_iv_list <- list(mean = iv_county_mean, p25 = iv_county_p25,
                       median = iv_county_pooled, p75 = iv_county_p75,
                       max = iv_county_max, min = iv_county_min, iqr = iv_county_iqr)
for (nm in names(pooled_iv_list)) {
  f_val <- tryCatch(fitstat(pooled_iv_list[[nm]], "ivwald")[["ivwald1::n_prior_posters"]]$stat,
                    error = function(e) NA_real_)
  cat(sprintf("  %-8s Wald F = %.1f\n", nm, f_val))
}

# ---------------------------------------------------------------------------
# IV / OLS: by service group (primary spec)
# ---------------------------------------------------------------------------
res_county       <- run_iv_by_service(df_iv_county, "ln_median_price",
                                      "system_peer_pressure_county_9m",
                                      c("County_State", "post_month"))
res_county_shop2 <- run_iv_by_large_group(df_iv_county, "SERVICE_GROUP_SHOP2",
                                          "ln_median_price",
                                          "system_peer_pressure_county_9m",
                                          c("County_State", "post_month"))
res_ols_county   <- run_ols_by_service(df_iv_county, "ln_median_price",
                                       c("County_State", "post_month"))

# Service-level first-stage summary
cat("\n=== SERVICE-LEVEL FIRST-STAGE DIAGNOSTICS (County 9m) ===\n")
cat(sprintf("%-30s  %7s  %7s  %7s  %8s\n",
            "Service Group", "fs_coef", "fs_t", "Wald F", "n_obs"))
cat(strrep("-", 65), "\n")
for (i in seq_len(nrow(res_county))) {
  r    <- res_county[i, ]
  flag <- if (!is.na(r$fs_f) && r$fs_f < 10) " << WEAK" else ""
  cat(sprintf("%-30s  %7.3f  %7.2f  %8.1f  %8d%s\n",
              as.character(r$SERVICE_GROUP),
              r$fs_coef, r$fs_t, r$fs_f, r$n_obs, flag))
}

# Wald F plot by service group
plot_data <- res_county %>%
  arrange(fs_f) %>%
  mutate(SERVICE_GROUP = factor(SERVICE_GROUP, levels = SERVICE_GROUP))

ggplot(plot_data, aes(x = fs_f, y = SERVICE_GROUP)) +
  geom_vline(xintercept = 10, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_segment(aes(x = 0, xend = fs_f, yend = SERVICE_GROUP),
               color = FSU_GARNET, linewidth = 0.7) +
  geom_point(color = FSU_GARNET, size = 2.5) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05)),
                     breaks = c(0, 10, 25, 50, 75, 100)) +
  labs(x = "First-Stage Wald F-Statistic", y = NULL,
       title = "First-Stage Instrument Strength by Service Group",
       subtitle = "County-level panel; instrument: system peer pressure (9-month trailing window)") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(color = "grey93", linewidth = 0.4),
        axis.text.y = element_text(size = 8), plot.title = element_text(size = 11, face = "bold"))

# ----------------------------------------------------------------------------
# OUTPUT_FIGURE_01_FIRST_STAGE_WALD
# ----------------------------------------------------------------------------
# Saved PDF: fig_first_stage_wald.pdf. Shows service-level first-stage Wald 
# F-statistics.
ggsave(filename = "fig_first_stage_wald.pdf",
       path = "~/Library/CloudStorage/OneDrive-FloridaStateUniversity/Hospital Price Transparency Paper/Plots/Paper Plots",
       width = 7, height = 8, device = cairo_pdf)

# OLS vs IV comparison table
ols_iv_compare_county <- res_ols_county %>%
  mutate(SERVICE_GROUP = as.character(SERVICE_GROUP)) %>%
  select(SERVICE_GROUP, ols_est = estimate_pct, ols_se = se_pct,
         ols_pval = pval, ols_sig = significant) %>%
  left_join(res_county %>%
              mutate(SERVICE_GROUP = as.character(SERVICE_GROUP)) %>%
              select(SERVICE_GROUP, iv_est = estimate_pct, iv_se = se_pct,
                     iv_pval = pval, iv_sig = significant,
                     fs_coef, fs_se, fs_t, fs_f),
            by = "SERVICE_GROUP") %>%
  mutate(bias = ols_est - iv_est,
         ols_stars = add_stars(ols_pval), iv_stars = add_stars(iv_pval),
         ols_fmt = sprintf("%.3f%s (%.3f)", ols_est, ols_stars, ols_se),
         iv_fmt  = sprintf("%.3f%s (%.3f)", iv_est,  iv_stars,  iv_se),
         fs_fmt  = sprintf("%.3f (t=%.2f, F=%.1f)", fs_coef, fs_t, fs_f)) %>%
  arrange(iv_est)

# ----------------------------------------------------------------------------
# OUTPUT_TABLE_05B_SERVICE_OLS_IV
# ----------------------------------------------------------------------------
# Console table comparing service-level OLS and IV estimates, including first-
# stage diagnostics.
cat("\n=== TABLE 5B: Service-Level OLS vs IV — County ===\n")
cat(sprintf("%-25s  %-22s  %-22s  %8s  %-28s\n",
            "Service Group", "OLS (SE)", "IV (SE)", "Bias", "First Stage (t, F)"))
cat(strrep("-", 115), "\n")
for (i in seq_len(nrow(ols_iv_compare_county))) {
  r <- ols_iv_compare_county[i, ]
  cat(sprintf("%-25s  %-22s  %-22s  %8.3f  %-28s\n",
              r$SERVICE_GROUP, r$ols_fmt, r$iv_fmt, r$bias, r$fs_fmt))
}
cat("\nBias summary:\n")
cat("  OLS < IV (downward bias):", sum(ols_iv_compare_county$bias < 0, na.rm=TRUE), "\n")
cat("  OLS > IV (upward bias):  ", sum(ols_iv_compare_county$bias > 0, na.rm=TRUE), "\n")
cat("  Mean bias (OLS - IV):    ", round(mean(ols_iv_compare_county$bias, na.rm=TRUE), 3), "%\n")
cat("  Median F across services:", round(median(ols_iv_compare_county$fs_f, na.rm=TRUE), 1), "\n")
cat("  Services with F < 10:   ", sum(ols_iv_compare_county$fs_f < 10, na.rm=TRUE), "\n")

print(res_county_shop2)

# ---------------------------------------------------------------------------
# Meta-regression
# ---------------------------------------------------------------------------
meta_data_county     <- build_meta_data(res_county)
meta_county          <- run_meta_regressions(meta_data_county)
print_meta_results(meta_county, label = "IV — County (9m primary)")

meta_data_ols_county <- build_meta_data(res_ols_county)
meta_county_ols      <- run_meta_regressions(meta_data_ols_county)
print_meta_results(meta_county_ols, label = "OLS — County (comparison)")

if (nrow(res_county_shop2) > 0) {
  cat("\n=== SHOP2 Group Estimates ===\n")
  print(res_county_shop2 %>%
          select(group, n_obs, estimate_pct, se_pct, pval, fs_coef, fs_t, fs_f, significant) %>%
          arrange(estimate_pct))
}

# ---------------------------------------------------------------------------
# All outcomes by service
# ---------------------------------------------------------------------------
res_county_outcomes <- run_iv_by_service_outcomes(
  data        = df_iv_county,
  outcomes    = outcomes_list,
  instrument  = "system_peer_pressure_county_9m",
  cluster_var = c("County_State", "post_month"),
  cum_col     = "N_HOSPITALS_CUMULATIVE_COUNTY"
)

cat("\n=== Significant Results by Outcome — County ===\n")
res_county_outcomes %>%
  group_by(outcome) %>%
  summarise(n_sig = sum(significant, na.rm=TRUE),
            n_sig_neg = sum(significant & estimate < 0, na.rm=TRUE),
            n_sig_pos = sum(significant & estimate > 0, na.rm=TRUE),
            mean_est  = round(mean(estimate_pct, na.rm=TRUE), 3),
            median_f  = round(median(fs_f, na.rm=TRUE), 1),
            pct_weak_f = round(mean(fs_f < 10, na.rm=TRUE) * 100, 1),
            .groups = "drop") %>%
  arrange(desc(n_sig)) %>% as.data.frame() %>% print()

tmp <- res_county_outcomes %>%
  select(SERVICE_GROUP, outcome, estimate_pct, se_pct, pval,
         fs_coef, fs_t, fs_f, significant) %>%
  mutate(across(where(is.numeric), ~round(.x, 3))) %>%
  arrange(outcome, estimate_pct) %>% as.data.frame()
# ----------------------------------------------------------------------------
# EXPORT_CSV_county_outcomes_detail
# ----------------------------------------------------------------------------
# Exports service-by-outcome county IV detail to county_outcomes_detail.csv.
write.csv(tmp, "county_outcomes_detail.csv", row.names = FALSE)
cat("Saved to county_outcomes_detail.csv\n")

######## Section 7: City — Robustness Check ###########
# ----------------------------------------------------------------------------
# Primary city instrument: system_peer_pressure_city_12m
# Clustering: city_state + post_month (two-way)
# City robustness analysis using system_peer_pressure_city_12m with city_state
#  + post_month clustering.
# First stage
fs_city <- feols(
  n_prior_posters ~ system_peer_pressure_city_12m + ln_total_beds |
    market_id + post_month,
  data    = df_iv_city,
  cluster = ~city_state + post_month
)
cat("\n=== CITY FIRST STAGE (PRIMARY: 12m, two-way clustering) ===\n")
summary(fs_city)
cat(sprintf("First-stage t-stat : %.2f\n",
            summary(fs_city)$coeftable[1, "t value"]))
cat(sprintf("Implied Wald F     : %.1f\n",
            summary(fs_city)$coeftable[1, "t value"]^2))

# IV: pooled
iv_city_pooled <- run_iv_pooled(df_iv_city, "ln_median_price",
                                "system_peer_pressure_city_12m",
                                c("city_state", "post_month"), verbose = TRUE)
cat("\n=== CITY IV POOLED ===\n")
summary(iv_city_pooled)
cat("Wald F:", fitstat(iv_city_pooled, "ivwald")[["ivwald1::n_prior_posters"]]$stat, "\n")

# IV: all outcomes
iv_city_mean   <- run_iv_pooled(df_iv_city, "ln_mean_price",
                                "system_peer_pressure_city_12m",
                                c("city_state", "post_month"))
iv_city_p25    <- run_iv_pooled(df_iv_city, "ln_p25_price",
                                "system_peer_pressure_city_12m",
                                c("city_state", "post_month"))
iv_city_p75    <- run_iv_pooled(df_iv_city, "ln_p75_price",
                                "system_peer_pressure_city_12m",
                                c("city_state", "post_month"))
iv_city_max    <- run_iv_pooled(df_iv_city, "ln_max_price",
                                "system_peer_pressure_city_12m",
                                c("city_state", "post_month"))
iv_city_min    <- run_iv_pooled(df_iv_city, "ln_min_price",
                                "system_peer_pressure_city_12m",
                                c("city_state", "post_month"))
iv_city_iqr    <- feols(
  ln_iqr_price ~ ln_total_beds | market_id + post_month |
    n_prior_posters ~ system_peer_pressure_city_12m,
  data    = df_iv_city %>% filter(N_HOSPITALS_CUMULATIVE_CITY > 1,
                                  is.finite(ln_iqr_price)),
  cluster = ~city_state + post_month
)

# IV: by service group
res_city     <- run_iv_by_service(df_iv_city, "ln_median_price",
                                  "system_peer_pressure_city_12m",
                                  c("city_state", "post_month"))
res_city_shop2 <- run_iv_by_large_group(df_iv_city, "SERVICE_GROUP_SHOP2",
                                        "ln_median_price",
                                        "system_peer_pressure_city_12m",
                                        c("city_state", "post_month"))

# City meta-regression
meta_data_city <- build_meta_data(res_city)
meta_city      <- run_meta_regressions(meta_data_city)
cat("\n=== CITY Meta-Regression ===\n")
print(summary(meta_city$simple))
print(summary(meta_city$granular))

# City all outcomes by service
res_city_outcomes <- run_iv_by_service_outcomes(
  data        = df_iv_city,
  outcomes    = outcomes_list,
  instrument  = "system_peer_pressure_city_12m",
  cluster_var = c("city_state", "post_month"),
  cum_col     = "N_HOSPITALS_CUMULATIVE_CITY"
)

cat("\n=== CITY: Significant Results by Outcome ===\n")
res_city_outcomes %>%
  group_by(outcome) %>%
  summarise(n_sig = sum(significant, na.rm=TRUE),
            n_sig_neg = sum(significant & estimate < 0, na.rm=TRUE),
            n_sig_pos = sum(significant & estimate > 0, na.rm=TRUE),
            mean_est  = mean(estimate_pct, na.rm=TRUE),
            .groups = "drop") %>%
  arrange(desc(n_sig)) %>% as.data.frame() %>% print()


######## Section 8: Two-Way Geographic Comparison (Table 14) ###########
# ----------------------------------------------------------------------------
# County (primary) vs City (robustness). Builds county-vs-city comparison
# tables and the shoppability-gradient comparison across geographic
# definitions.

two_way_compare <- res_county %>%
  mutate(SERVICE_GROUP = as.character(SERVICE_GROUP)) %>%
  select(SERVICE_GROUP, county_est = estimate_pct, county_se = se_pct,
         county_pval = pval, county_wald = fs_f) %>%
  left_join(res_city %>%
              mutate(SERVICE_GROUP = as.character(SERVICE_GROUP)) %>%
              select(SERVICE_GROUP, city_est = estimate_pct, city_se = se_pct,
                     city_pval = pval, city_wald = fs_f),
            by = "SERVICE_GROUP") %>%
  mutate(county_fmt = sprintf("%.3f%s (%.3f)", county_est,
                              add_stars(county_pval), county_se),
         city_fmt   = sprintf("%.3f%s (%.3f)", city_est,
                              add_stars(city_pval), city_se)) %>%
  arrange(county_est)

# ----------------------------------------------------------------------------
# OUTPUT_TABLE_14_GEO_COMPARISON
# ----------------------------------------------------------------------------
# Console table comparing county primary estimates to city robustness 
# estimates.
cat("\n=== TABLE 14: Two-Way Geographic Comparison ===\n")
cat(sprintf("%-25s  %-22s  %-8s  %-22s  %-8s\n",
            "Service Group", "County (Primary)", "Wald F",
            "City (Robustness)", "Wald F"))
cat(strrep("-", 95), "\n")
for (i in seq_len(nrow(two_way_compare))) {
  r <- two_way_compare[i,]
  cat(sprintf("%-25s  %-22s  %-8.1f  %-22s  %-8.1f\n",
              r$SERVICE_GROUP, r$county_fmt,
              replace(r$county_wald, is.na(r$county_wald), 0),
              r$city_fmt,
              replace(r$city_wald, is.na(r$city_wald), 0)))
}
cat("*** p<0.01  ** p<0.05  * p<0.10\n")

# Meta-regression shoppability coefficient comparison
cat("\n=== Shoppability gradient: County vs City ===\n")
cat(sprintf("%-10s  %10s  %10s  %8s\n", "Level", "Coef", "SE", "p-value"))
for (nm in c("County", "City")) {
  m  <- switch(nm, County = meta_county$simple, City = meta_city$simple)
  cf <- coef(summary(m))
  cat(sprintf("%-10s  %10.4f  %10.4f  %8.4f\n",
              nm, cf["is_shoppableTRUE","Estimate"],
              cf["is_shoppableTRUE","Std. Error"],
              cf["is_shoppableTRUE","Pr(>|t|)"]))
}


######## Section 9: Figures ###########
# ----------------------------------------------------------------------------
# Creates the first set of paper figures: county IV forest, meta-regression 
# scatter, OLS-vs-IV forest/scatter, heatmap, and county/city figure.
# Figure 1A: Forest plot — county IV by service group
# ----------------------------------------------------------------------------
# OUTPUT_FIGURE_02A_COUNTY_IV_FOREST
# ----------------------------------------------------------------------------
# Figure object fig_forest_county: county IV forest plot by service group.
fig_forest_county <- ggplot(res_county,
                            aes(x = estimate_pct, y = SERVICE_GROUP, color = significant)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5, color = "grey50") +
  geom_errorbarh(aes(xmin = ci_lo_pct, xmax = ci_hi_pct), height = 0, linewidth = 0.45) +
  geom_point(size = 2) +
  scale_color_manual(values = c("FALSE" = FSU_GOLD, "TRUE" = FSU_GARNET),
                     labels = c("Not significant", "Significant (95%)")) +
  scale_x_continuous(labels = function(x) paste0(x, "%")) +
  labs(x = "% change in median negotiated price per additional prior poster",
       y = NULL, color = NULL,
       title = "A. IV Estimates by Service Group (County Markets, 9m instrument)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", axis.text.y = element_text(size = 7.5),
        panel.grid.major.y = element_blank(), panel.grid.minor = element_blank())

# Figure 1B: Meta-regression scatter
# ----------------------------------------------------------------------------
# OUTPUT_FIGURE_02B_METAREG_SCATTER
# ----------------------------------------------------------------------------
# Figure object fig_meta_scatter: shoppability meta-regression scatter.
fig_meta_scatter <- ggplot(meta_data_county,
                           aes(x = as.numeric(is_shoppable), y = estimate_pct,
                               label = SERVICE_GROUP)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(aes(color = is_shoppable), size = 2) +
  geom_smooth(method = "lm", se = TRUE, weights = meta_data_county$weight,
              color = FSU_GARNET, fill = FSU_GOLD, alpha = 0.2) +
  geom_text_repel(size = 2.3, max.overlaps = 20, box.padding = 0.3) +
  scale_color_manual(values = c("FALSE" = FSU_GOLD, "TRUE" = FSU_GARNET),
                     labels = c("Non-shoppable", "Shoppable")) +
  scale_x_continuous(breaks = c(0,1), labels = c("Non-shoppable","Shoppable"),
                     expand = expansion(mult = 0.3)) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(x = NULL, y = "IV Estimate (% change in median price)", color = NULL,
       title = "B. Shoppability Meta-Regression") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")

fig_main <- fig_forest_county + fig_meta_scatter + plot_layout(widths = c(1.8, 1))
# ----------------------------------------------------------------------------
# OUTPUT_FIGURE_02_METAREG_COMBINED
# ----------------------------------------------------------------------------
# Combined paper figure fig_main: county forest plot plus meta-regression 
# scatter.
print(fig_main)

# Figure 2: OLS vs IV forest plot
ols_iv_plot_data <- bind_rows(
  res_ols_county %>%
    mutate(SERVICE_GROUP = factor(as.character(SERVICE_GROUP),
                                  levels = levels(res_county$SERVICE_GROUP)),
           model = "OLS"),
  res_county %>% mutate(se = se_pct, model = "IV (9m)")
) %>%
  mutate(model = factor(model, levels = c("OLS", "IV (9m)")))

# ----------------------------------------------------------------------------
# OUTPUT_FIGURE_03_OLS_IV_FOREST
# ----------------------------------------------------------------------------
# Figure object fig_ols_iv: OLS vs IV forest plot by service group.
fig_ols_iv <- ggplot(ols_iv_plot_data,
                     aes(x = estimate_pct, y = SERVICE_GROUP, color = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5, color = "grey50") +
  geom_errorbarh(aes(xmin = ci_lo_pct, xmax = ci_hi_pct), height = 0,
                 linewidth = 0.45, position = position_dodge(0.6)) +
  geom_point(size = 2, position = position_dodge(0.6)) +
  scale_color_manual(values = c("OLS" = FSU_GOLD, "IV (9m)" = FSU_GARNET)) +
  scale_x_continuous(labels = function(x) paste0(x, "%")) +
  labs(x = "% change in median negotiated price per additional prior poster",
       y = NULL, color = NULL,
       title = "OLS vs IV Estimates by Service Group (County Markets)",
       subtitle = paste0("OLS is more negative than IV for ",
                         sum(ols_iv_compare_county$bias < 0, na.rm=TRUE),
                         " of 46 services")) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", axis.text.y = element_text(size = 7.5),
        panel.grid.major.y = element_blank(), panel.grid.minor = element_blank())
print(fig_ols_iv)


library(stringr)
library(ggrepel)

# ------------------------------------------------------------
# OLS vs IV scatter plot
# IV on x-axis, OLS on y-axis
# ------------------------------------------------------------

ols_iv_scatter_df <- res_ols_county %>%
  select(
    SERVICE_GROUP,
    ols_estimate_pct = estimate_pct,
    ols_se_pct       = se_pct,
    ols_ci_lo_pct    = ci_lo_pct,
    ols_ci_hi_pct    = ci_hi_pct
  ) %>%
  left_join(
    res_county %>%
      select(
        SERVICE_GROUP,
        iv_estimate_pct = estimate_pct,
        iv_se_pct       = se_pct,
        iv_ci_lo_pct    = ci_lo_pct,
        iv_ci_hi_pct    = ci_hi_pct
      ),
    by = "SERVICE_GROUP"
  ) %>%
  mutate(
    SERVICE_GROUP = as.character(SERVICE_GROUP),
    
    # Negative means OLS estimate is more negative than IV estimate
    ols_minus_iv = ols_estimate_pct - iv_estimate_pct,
    ols_more_negative = ols_minus_iv < 0,
    
    service_label = str_wrap(SERVICE_GROUP, width = 22)
  ) %>%
  filter(
    !is.na(ols_estimate_pct),
    !is.na(iv_estimate_pct)
  )

# Count services where OLS is more negative than IV
n_ols_more_negative <- sum(ols_iv_scatter_df$ols_more_negative, na.rm = TRUE)
n_services <- nrow(ols_iv_scatter_df)

# Shared axis limits so the 45-degree line is meaningful
axis_limits <- range(
  c(
    ols_iv_scatter_df$ols_estimate_pct,
    ols_iv_scatter_df$iv_estimate_pct,
    ols_iv_scatter_df$ols_ci_lo_pct,
    ols_iv_scatter_df$ols_ci_hi_pct,
    ols_iv_scatter_df$iv_ci_lo_pct,
    ols_iv_scatter_df$iv_ci_hi_pct
  ),
  na.rm = TRUE
)

axis_pad <- 0.08 * diff(axis_limits)
axis_limits <- c(axis_limits[1] - axis_pad, axis_limits[2] + axis_pad)

# ----------------------------------------------------------------------------
# OUTPUT_FIGURE_03B_OLS_IV_SCATTER
# ----------------------------------------------------------------------------
# Figure object fig_ols_iv_scatter: IV estimate on x-axis, OLS estimate on 
# y-axis.
fig_ols_iv_scatter <- ggplot(
  ols_iv_scatter_df,
  aes(x = iv_estimate_pct, y = ols_estimate_pct)
) +
  
  # Zero reference lines
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.4,
    color = "grey65"
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.4,
    color = "grey65"
  ) +
  
  # 45-degree line: OLS = IV
  geom_abline(
    intercept = 0,
    slope = 1,
    linewidth = 0.75,
    color = "grey35"
  ) +
  
  # Points
  geom_point(
    size = 2.8,
    alpha = 0.85,
    color = FSU_GARNET
  ) +
  
  # Optional: label the largest OLS-IV gaps
  ggrepel::geom_text_repel(
    data = ols_iv_scatter_df %>%
      mutate(abs_gap = abs(ols_minus_iv)) %>%
      slice_max(abs_gap, n = 6),
    aes(label = service_label),
    size = 2.6,
    max.overlaps = 20,
    show.legend = FALSE,
    box.padding = 0.25,
    point.padding = 0.2,
    min.segment.length = 0,
    color = "grey25"
  ) +
  
  annotate(
    "text",
    x = axis_limits[1] + 0.05 * diff(axis_limits),
    y = axis_limits[2] - 0.08 * diff(axis_limits),
    label = "Below line:\nOLS more negative than IV",
    hjust = 0,
    vjust = 1,
    size = 3.1,
    color = "grey35",
    fontface = "italic"
  ) +
  
  coord_equal(
    xlim = axis_limits,
    ylim = axis_limits
  ) +
  
  scale_x_continuous(
    labels = function(x) paste0(x, "%")
  ) +
  
  scale_y_continuous(
    labels = function(x) paste0(x, "%")
  ) +
  
  labs(
    x = "IV estimate: % change in median negotiated price\nper additional prior poster",
    y = "OLS estimate: % change in median negotiated price\nper additional prior poster",
    title = "OLS Estimates Are Systematically More Negative Than IV Estimates",
    subtitle = paste0(
      "OLS is more negative than IV for ",
      n_ols_more_negative, " of ", n_services,
      " service groups"
    ),
    caption = paste0(
      "Notes: Each point is a service group. The 45-degree line represents equality between OLS and IV estimates. ",
      "Points below the line indicate that the OLS estimate is more negative than the IV estimate. ",
      "Estimates are from county-market specifications."
    )
  ) +
  
  theme_bw(base_size = 11) +
  
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9.5, color = "grey30"),
    plot.caption = element_text(
      size = 7.5,
      color = "grey40",
      hjust = 0,
      lineheight = 1.25
    ),
    axis.title.x = element_text(size = 9.5),
    axis.title.y = element_text(size = 9.5),
    panel.grid.minor = element_blank()
  )

print(fig_ols_iv_scatter)






# Figure A1: Heatmap by service and outcome
heat_data_county <- res_county_outcomes %>%
  filter(significant) %>%
  mutate(
    sig_label = case_when(pval < 0.01 ~ "***", pval < 0.05 ~ "**",
                          pval < 0.10 ~ "*",   TRUE ~ ""),
    outcome = factor(outcome, levels = c("Min","P25","Median","Mean","P75","Max","IQR")),
    SERVICE_GROUP = factor(SERVICE_GROUP,
                           levels = res_county_outcomes %>%
                             filter(outcome == "Median") %>%
                             arrange(estimate_pct) %>% pull(SERVICE_GROUP))
  )

# ----------------------------------------------------------------------------
# OUTPUT_FIGURE_A1_OUTCOME_HEATMAP
# ----------------------------------------------------------------------------
# Figure object fig_heatmap_county: significant county effects by service and 
# outcome.
fig_heatmap_county <- ggplot(heat_data_county,
                             aes(x = outcome, y = SERVICE_GROUP, fill = estimate_pct)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sig_label), size = 2.8) +
  scale_fill_gradient2(low = FSU_GARNET, mid = "white", high = FSU_GOLD,
                       midpoint = 0, name = "Est. (%)") +
  labs(x = "Price Outcome", y = NULL,
       title = "IV Estimates by Service Group and Price Outcome (County, 9m)") +
  theme_bw(base_size = 10) +
  theme(axis.text.y = element_text(size = 7.5), panel.grid = element_blank())
print(fig_heatmap_county)

# Figure A2: Two-way geographic comparison forest plot
two_way_plot <- bind_rows(
  res_county %>% mutate(SERVICE_GROUP = as.character(SERVICE_GROUP), level = "County (9m)"),
  res_city   %>% mutate(SERVICE_GROUP = as.character(SERVICE_GROUP), level = "City (12m)")
) %>%
  mutate(
    SERVICE_GROUP = factor(SERVICE_GROUP, levels = levels(res_county$SERVICE_GROUP)),
    level   = factor(level, levels = c("County (9m)", "City (12m)")),
    y_dodge = as.numeric(SERVICE_GROUP) +
      case_when(level == "County (9m)" ~ -0.2, level == "City (12m)" ~ 0.2)
  )

# ----------------------------------------------------------------------------
# OUTPUT_FIGURE_GEO_COMPARISON
# ----------------------------------------------------------------------------
# Figure object fig_two_way: county-vs-city estimate comparison.
fig_two_way <- ggplot(two_way_plot,
                      aes(x = estimate_pct, y = y_dodge, color = level, shape = significant)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, color = "grey50") +
  geom_errorbarh(aes(xmin = ci_lo_pct, xmax = ci_hi_pct), height = 0,
                 linewidth = 0.4, alpha = 0.6) +
  geom_point(size = 1.8, stroke = 0.7) +
  scale_color_manual(values = c("County (9m)" = FSU_GARNET, "City (12m)" = FSU_GOLD),
                     name = "Market Definition") +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                     labels = c("TRUE" = "Significant (95%)", "FALSE" = "Not significant"),
                     name = NULL) +
  scale_y_continuous(breaks = seq_along(levels(res_county$SERVICE_GROUP)),
                     labels = levels(res_county$SERVICE_GROUP),
                     expand = expansion(add = 0.6)) +
  scale_x_continuous(labels = function(x) paste0(x, "%")) +
  labs(x = "% change in median negotiated price per additional prior poster", y = NULL,
       title = "IV Estimates: County vs City Market Definition",
       caption = "County: 9m instrument, two-way clustering. City: 12m instrument, two-way clustering.") +
  theme_bw(base_size = 10) +
  theme(legend.position = "top", axis.text.y = element_text(size = 7),
        panel.grid.major.y = element_blank(), panel.grid.minor = element_blank())
print(fig_two_way)




######## Section 10: Saturation Analysis — Diminishing Returns ###########
# ----------------------------------------------------------------------------
# Saturation analysis: bins prior posters into exposure ranges and checks 
# whether disclosure effects attenuate as markets become more saturated.

cat("=== SATURATION DISTRIBUTION ===\n")
print(summary(df_iv_county$n_prior_posters))
print(quantile(df_iv_county$n_prior_posters,
               probs = c(0.05,0.10,0.25,0.50,0.75,0.90,0.95,0.99), na.rm=TRUE))

df_iv_county <- df_iv_county %>%
  mutate(
    prior_range = factor(case_when(
      n_prior_posters >= BIN1_LOW & n_prior_posters <= BIN1_HIGH ~
        paste0(BIN1_LOW, "\u2013", BIN1_HIGH),
      n_prior_posters >= BIN2_LOW & n_prior_posters <= BIN2_HIGH ~
        paste0(BIN2_LOW, "\u2013", BIN2_HIGH),
      n_prior_posters >= BIN2_HIGH + 1 ~
        paste0(BIN2_HIGH + 1, "+"),
      TRUE ~ NA_character_
    ), levels = range_levels)
  )

cat("\nObs per bin:\n"); print(table(df_iv_county$prior_range, useNA = "always"))

low_label  <- paste0(BIN1_LOW, "-", BIN1_HIGH, " prior posters")
mid_label  <- paste0(BIN2_LOW, "-", BIN2_HIGH, " prior posters")
high_label <- paste0(BIN2_HIGH + 1, "+ prior posters")

# ----------------------------------------------------------------------------
# FUNCTION run_iv_bin()
# ----------------------------------------------------------------------------
# Runs pooled IV within a saturation bin of prior posters.
run_iv_bin <- function(data, label, instrument, cluster_var, min_n = 50) {
  if (nrow(data) < min_n) return(NULL)
  if (length(unique(na.omit(data$n_prior_posters))) < 2) return(NULL)
  instr <- data[[instrument]]
  if (is.null(instr) || length(unique(na.omit(instr))) < 2) return(NULL)
  fit <- tryCatch(
    feols(as.formula(paste0(
      "ln_median_price ~ 1 | market_id + post_month | n_prior_posters ~ ", instrument)),
      data    = data,
      cluster = as.formula(paste0("~", paste(cluster_var, collapse = " + ")))),
    error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  ct <- coeftable(fit)
  if (!"fit_n_prior_posters" %in% rownames(ct)) return(NULL)
  wf <- tryCatch(fitstat(fit, "ivwald")[["ivwald1::n_prior_posters"]]$stat,
                 error = function(e) NA_real_)
  data.frame(bin = label,
             estimate = ct["fit_n_prior_posters", "Estimate"],
             se       = ct["fit_n_prior_posters", "Std. Error"],
             pval     = ct["fit_n_prior_posters", "Pr(>|t|)"],
             wald_f   = wf, n = nrow(data))
}

# Pooled saturation bins (primary: county 9m)
results_3bin_county <- bind_rows(
  run_iv_bin(df_iv_county %>% filter(prior_range == range_levels[1]),
             low_label,  "system_peer_pressure_county_9m", c("County_State","post_month")),
  run_iv_bin(df_iv_county %>% filter(prior_range == range_levels[2]),
             mid_label,  "system_peer_pressure_county_9m", c("County_State","post_month")),
  run_iv_bin(df_iv_county %>% filter(prior_range == range_levels[3]),
             high_label, "system_peer_pressure_county_9m", c("County_State","post_month"))
) %>%
  mutate(pct = 100*(exp(estimate)-1), pct_low = 100*(exp(estimate-1.96*se)-1),
         pct_high = 100*(exp(estimate+1.96*se)-1),
         bin = factor(bin, levels = c(low_label, mid_label, high_label)))

cat("\n=== THREE-BIN SATURATION — County (Pooled) ===\n")
print(results_3bin_county %>% select(bin, n, pct, pct_low, pct_high, wald_f) %>%
        mutate(across(where(is.numeric), ~round(.x, 3))))

# By service group
# ----------------------------------------------------------------------------
# FUNCTION run_iv_bin_by_service()
# ----------------------------------------------------------------------------
# Runs service-level IV within a saturation bin of prior posters.
run_iv_bin_by_service <- function(data, bin_label, instrument,
                                  cluster_var, min_obs = 100) {
  service_groups <- unique(data$SERVICE_GROUP)
  do.call(rbind, lapply(service_groups, function(sg) {
    df_sub <- data[data$SERVICE_GROUP == sg, ]
    if (nrow(df_sub) < min_obs) return(NULL)
    res <- run_iv_bin(df_sub, bin_label, instrument, cluster_var)
    if (is.null(res)) return(NULL)
    res$SERVICE_GROUP <- sg
    res$estimate_pct  <- res$estimate * 100
    res$se_pct        <- res$se * 100
    res$ci_lo_pct     <- (res$estimate - 1.96 * res$se) * 100
    res$ci_hi_pct     <- (res$estimate + 1.96 * res$se) * 100
    res$significant   <- !is.na(res$estimate) & !is.na(res$se) &
      sign(res$estimate - 1.96*res$se) == sign(res$estimate + 1.96*res$se)
    res
  }))
}

sat_service_low  <- run_iv_bin_by_service(
  df_iv_county %>% filter(prior_range == range_levels[1]),
  low_label,  "system_peer_pressure_county_9m", c("County_State","post_month"))
sat_service_mid  <- run_iv_bin_by_service(
  df_iv_county %>% filter(prior_range == range_levels[2]),
  mid_label,  "system_peer_pressure_county_9m", c("County_State","post_month"))
sat_service_high <- run_iv_bin_by_service(
  df_iv_county %>% filter(prior_range == range_levels[3]),
  high_label, "system_peer_pressure_county_9m", c("County_State","post_month"))

sat_service_all <- bind_rows(sat_service_low, sat_service_mid, sat_service_high) %>%
  mutate(estimate_pct = estimate*100, se_pct = se*100,
         ci_lo_pct = (estimate-1.96*se)*100, ci_hi_pct = (estimate+1.96*se)*100,
         significant = !is.na(estimate) & !is.na(se) &
           sign(estimate-1.96*se) == sign(estimate+1.96*se),
         bin = factor(bin, levels = c(low_label, mid_label, high_label)))

# By large service group (SHOP2)
sat_shop2_low  <- run_iv_by_large_group(
  df_iv_county %>% filter(prior_range == range_levels[1]),
  "SERVICE_GROUP_SHOP2", "ln_median_price",
  "system_peer_pressure_county_9m", c("County_State","post_month")) %>%
  mutate(bin = low_label)
sat_shop2_mid  <- run_iv_by_large_group(
  df_iv_county %>% filter(prior_range == range_levels[2]),
  "SERVICE_GROUP_SHOP2", "ln_median_price",
  "system_peer_pressure_county_9m", c("County_State","post_month")) %>%
  mutate(bin = mid_label)
sat_shop2_high <- run_iv_by_large_group(
  df_iv_county %>% filter(prior_range == range_levels[3]),
  "SERVICE_GROUP_SHOP2", "ln_median_price",
  "system_peer_pressure_county_9m", c("County_State","post_month")) %>%
  mutate(bin = high_label)
sat_shop2_all <- bind_rows(sat_shop2_low, sat_shop2_mid, sat_shop2_high) %>%
  mutate(bin = factor(bin, levels = c(low_label, mid_label, high_label)))

cat("\n=== THREE-BIN SATURATION — By Large Service Group (SHOP2) ===\n")
print(sat_shop2_all %>%
        select(bin, group, estimate_pct, se_pct, pval, fs_f, significant) %>%
        mutate(across(where(is.numeric), ~round(.x, 3))) %>% arrange(bin, estimate_pct))

# Saturation meta-regression across classification schemes (18, see n_shop_schemes below)
shoppability_schemes <- list(
  scheme_cms = list(label = "CMS Rule Definition",
                    shoppable = c("Ultrasound","CT Lung","Mammography"), nonshoppable = NULL),
  scheme_theory = list(label = "Theory-Based",
                       shoppable = c("Ultrasound","CT Lung","Mammography"),
                       nonshoppable = c("Biopsy","Colonoscopy","Endoscopy")),
  scheme_broad = list(label = "Broad Shoppable",
                      shoppable = c("Ultrasound","CT Lung","Mammography","X-Ray"),
                      nonshoppable = c("Biopsy","Colonoscopy","Endoscopy")),
  scheme_imaging = list(label = "Imaging vs Procedural",
                        shoppable = c("Ultrasound","CT Lung","Mammography","X-Ray","MRI"),
                        nonshoppable = c("Biopsy","Colonoscopy","Endoscopy")),
  scheme_theory_v2 = list(label = "Theory-Based V2 (MRI as nonshoppable)",
                          shoppable = c("Ultrasound","CT Lung","Mammography"),
                          nonshoppable = c("Biopsy","Colonoscopy","Endoscopy","MRI")),
  scheme_split_nonshop = list(label = "Split Non-Shoppable: MRI vs Procedural",
                              shoppable = c("Ultrasound","CT Lung","Mammography"),
                              nonshoppable_mri  = c("MRI"),
                              nonshoppable_proc = c("Biopsy","Colonoscopy","Endoscopy"),
                              nonshoppable = NULL),
  scheme_cms_statutory = list(label = "CMS Statutory Shoppable List",
                              shoppable = c("Ultrasound","CT Lung","Mammography","Colonoscopy","X-Ray"),
                              nonshoppable = c("Biopsy","MRI")),
  scheme_anatomical = list(label = "High vs Low Shoppability Within Modality",
                           shoppable = c("CT Lung","Mammography","Ultrasound OB","Ultrasound Breast",
                                         "Ultrasound Abdomen","X-Ray Chest","X-Ray Extremity"),
                           nonshoppable = c("Biopsy Pancreas","Biopsy Liver","Biopsy Lung",
                                            "Biopsy Kidney","Biopsy Bone","Colonoscopy","Endoscopy","MRI")),
  scheme_ct_broad = list(label = "CT-Inclusive (all CT modalities shoppable)",
                         shoppable = c("Ultrasound","CT","Mammography"),
                         nonshoppable = c("Biopsy","Colonoscopy","Endoscopy","MRI")),
  scheme_ct_broad_ex_angio = list(label = "CT-Inclusive Except Angio (acute vascular CT excluded)",
                                  shoppable = c("Ultrasound","CT Lung","CT Brain/Head","CT Spine",
                                                "CT Neck","CT Extremity","CT Abdomen","CT Chest",
                                                "CT Pelvis","CT Other","Mammography"),
                                  nonshoppable = c("Biopsy","Colonoscopy","Endoscopy","MRI","CT Angio")),
  
  # --------------------------------------------------------------------------
  # ADDITIONAL SCHEMES 11-18: Eight alternative shoppability frameworks
  # --------------------------------------------------------------------------
  # Added at Danny's request to extend the specification sensitivity analysis
  # beyond the original ten schemes. Each framework classifies the same 46
  # service groups on a different conceptual dimension (clinical invasiveness,
  # CMS legal status, upfront cash-market availability, geography, urgency,
  # staffing, liability exposure, No Surprises Act exposure) rather than
  # re-deriving Theory V2. Services not named in either the shoppable or
  # nonshoppable vector for a given scheme fall into the omitted/baseline
  # category, which the existing modeling functions (build_meta_data_scheme,
  # run_meta_scheme) already treat as "Intermediate" -- consistent with how
  # scheme_theory, scheme_broad, etc. leave MRI/X-Ray unclassified above.
  #
  # NOTE ON SOURCE MATERIAL: these eight frameworks were drafted externally
  # (not derived from this project's own CMS filings or market data), and the
  # source document's own tallies contained several off-by-one arithmetic
  # errors and one substantive inconsistency, corrected here:
  #   1) The source document works from a 45-service total; this project's
  #      locked results (and the SERVICE_GROUP levels actually in the data,
  #      including "X-Ray Spine") use 46. All eight schemes below are coded
  #      against the real 46-group universe, not the source's 45-group count.
  #   2) Several of the source's per-division subtotals (e.g., "24", "31",
  #      "26", "10") did not match the number of items actually itemized
  #      under that heading. Where this happened, the itemized list was
  #      used as the source of truth over the stated count.
  #   3) scheme_div1_operations follows the source document's own
  #      itemization, which -- despite its prose claiming "all 7 X-Rays" --
  #      omits "X-Ray Spine" from the shoppable list for that framework
  #      specifically. That's preserved here rather than silently corrected,
  #      since it's unclear whether the omission was an intended edge case
  #      or a copy-paste error; it leaves "X-Ray Spine" in the baseline
  #      category for this scheme only. Worth a second look before use.
  #   4) The CMS-specific service assignments in scheme_div2_cms_legal
  #      (which services fall in the 70 CMS-specified vs. 230 hospital-
  #      selected list) reflect the source document's assertions and have
  #      NOT been checked against CMS's actual published list of the 70
  #      specified services. The overall 70/230/300 structure is confirmed
  #      accurate against CMS documentation, but the service-level mapping
  #      should be verified against the real CMS table before this scheme
  #      is cited as "CMS-derived" in the paper.
  #   5) scheme_div4_geographic and scheme_alt8_no_surprises produce an
  #      IDENTICAL partition of services (Mammography/Ultrasound/X-Ray
  #      shoppable; Biopsy/Colonoscopy/Endoscopy nonshoppable; CT/MRI
  #      intermediate), despite being framed around different mechanisms
  #      (facility access vs. balance-billing exposure). They are kept as
  #      separate entries below for narrative completeness, but they will
  #      produce numerically identical coefficients -- not independent
  #      robustness checks.
  # --------------------------------------------------------------------------
  
  scheme_div1_operations = list(
    label = "Alt: Core Operations Framework (clinical invasiveness)",
    shoppable = c("Mammography","Ultrasound","CT",
                  "X-Ray Other","X-Ray Skull/Head","X-Ray Chest",
                  "X-Ray Abdomen","X-Ray Pelvis","X-Ray Extremity"),
    # X-Ray Spine deliberately excluded -- see note (3) above. Baseline
    # (Intermediate) = all 10 MRI + X-Ray Spine.
    nonshoppable = c("Biopsy","Colonoscopy","Endoscopy")
  ),
  
  scheme_div2_cms_legal = list(
    label = "Alt: CMS Legal/Regulatory Framework (70 vs 230 vs non-listed)",
    shoppable = c("Mammography","Ultrasound",
                  "X-Ray Other","X-Ray Skull/Head","X-Ray Chest","X-Ray Abdomen",
                  "X-Ray Pelvis","X-Ray Extremity",
                  "CT Brain/Head","CT Abdomen","MRI Brain/Head","MRI Spine",
                  "Colonoscopy","Endoscopy"),
    # Baseline (Intermediate) = CT Lung/Spine/Neck/Extremity/Chest/Pelvis/Other,
    # MRI Other/Pelvis/Chest/Extremity/Breast/Abdomen/Neck, Biopsy Thyroid,
    # Biopsy Breast, Biopsy Lymph Node, X-Ray Spine. See note (4): service-
    # level CMS-list membership unverified against the actual CMS table.
    nonshoppable = c("Biopsy Other","Biopsy Pancreas","Biopsy Kidney","Biopsy Bone",
                     "Biopsy Liver","Biopsy Lung","CT Angio","MRI Angio")
  ),
  
  scheme_div3_mdsave = list(
    label = "Alt: Upfront Cash-Market Framework (MDsave-style voucher availability)",
    shoppable = c("CT","MRI","X-Ray","Ultrasound","Mammography"),
    # Baseline (Intermediate) = Biopsy Thyroid, Biopsy Breast,
    # Biopsy Lymph Node, Colonoscopy, Endoscopy.
    nonshoppable = c("Biopsy Other","Biopsy Pancreas","Biopsy Kidney",
                     "Biopsy Bone","Biopsy Liver","Biopsy Lung")
  ),
  
  scheme_div4_geographic = list(
    label = "Alt: Geographic/Facility Access Framework (retail vs freestanding vs hospital-locked)",
    shoppable = c("Mammography","Ultrasound","X-Ray"),
    # Baseline (Intermediate) = all 10 CT + all 10 MRI.
    nonshoppable = c("Biopsy","Colonoscopy","Endoscopy")
  ),
  
  scheme_div5_urgency = list(
    label = "Alt: Diagnostic Urgency/Lead-Time Framework",
    shoppable = c("Mammography","Ultrasound OB","X-Ray Spine","Colonoscopy"),
    # Baseline (Intermediate) = all CT, all MRI except Spine, all remaining
    # Ultrasound/X-Ray, Endoscopy.
    nonshoppable = c("Biopsy","MRI Spine")
  ),
  
  scheme_alt6_staffing = list(
    label = "Alt: Staffing/Specialist Framework (generalist vs sub-specialist)",
    shoppable = c("CT Lung","CT Brain/Head","CT Spine","CT Neck","CT Extremity",
                  "CT Abdomen","CT Chest","CT Pelvis","CT Other",
                  "MRI Other","MRI Pelvis","MRI Chest","MRI Extremity",
                  "MRI Brain/Head","MRI Abdomen","MRI Spine","MRI Neck",
                  "Mammography","Ultrasound","X-Ray"),
    # Baseline (Intermediate) = CT Angio, MRI Angio, MRI Breast,
    # Biopsy Thyroid, Biopsy Breast, Biopsy Lymph Node.
    nonshoppable = c("Biopsy Other","Biopsy Pancreas","Biopsy Kidney","Biopsy Bone",
                     "Biopsy Liver","Biopsy Lung","Colonoscopy","Endoscopy")
  ),
  
  scheme_alt7_liability = list(
    label = "Alt: Incident Reporting/Liability Framework (risk-based)",
    shoppable = c("Ultrasound","Mammography"),
    # Baseline (Intermediate) = all CT, all MRI, all X-Ray.
    nonshoppable = c("Biopsy","Colonoscopy","Endoscopy")
  ),
  
  scheme_alt8_no_surprises = list(
    label = "Alt: No Surprises Act Legal Framework (balance-billing exposure)",
    # NOTE: identical partition to scheme_div4_geographic -- see note (5).
    shoppable = c("Mammography","Ultrasound","X-Ray"),
    nonshoppable = c("Biopsy","Colonoscopy","Endoscopy")
  )
)

# Number of active classification schemes -- referenced below instead of
# hardcoding "10"/"Ten" so titles/captions/logging stay correct if schemes are
# ever added or removed again. Currently 18 (10 original + 8 added from
# Danny's alternative-framework document).
n_shop_schemes <- length(shoppability_schemes)



# ----------------------------------------------------------------------------
# FUNCTION build_meta_data_scheme()
# ----------------------------------------------------------------------------
# Builds meta-regression data for a specific shoppability classification 
# scheme.
build_meta_data_scheme <- function(res_service, scheme) {
  df <- res_service %>%
    mutate(SERVICE_GROUP = as.character(SERVICE_GROUP)) %>%
    filter(!is.na(se_pct), is.finite(se_pct), se_pct > 0,
           !is.na(estimate_pct), is.finite(estimate_pct))
  if (nrow(df) == 0) return(df)
  df$is_shoppable <- Reduce(`|`, lapply(scheme$shoppable, function(term) {
    grepl(term, df$SERVICE_GROUP, fixed = TRUE) | df$SERVICE_GROUP == term }))
  if (!is.null(scheme$nonshoppable)) {
    df$is_nonshoppable <- Reduce(`|`, lapply(scheme$nonshoppable, function(term) {
      grepl(term, df$SERVICE_GROUP, fixed = TRUE) | df$SERVICE_GROUP == term }))
  } else df$is_nonshoppable <- FALSE
  if (!is.null(scheme$nonshoppable_mri)) {
    df$is_nonshoppable_mri <- Reduce(`|`, lapply(scheme$nonshoppable_mri, function(term) {
      grepl(term, df$SERVICE_GROUP, fixed = TRUE) | df$SERVICE_GROUP == term }))
  } else df$is_nonshoppable_mri <- FALSE
  if (!is.null(scheme$nonshoppable_proc)) {
    df$is_nonshoppable_proc <- Reduce(`|`, lapply(scheme$nonshoppable_proc, function(term) {
      grepl(term, df$SERVICE_GROUP, fixed = TRUE) | df$SERVICE_GROUP == term }))
  } else df$is_nonshoppable_proc <- FALSE
  df$is_nonshoppable_any <- df$is_nonshoppable | df$is_nonshoppable_mri | df$is_nonshoppable_proc
  df$weight <- 1 / (df$se_pct^2)
  df$fs_f   <- if ("fs_f" %in% names(df)) df$fs_f else NA_real_
  df %>% filter(is.finite(weight))
}

# ----------------------------------------------------------------------------
# FUNCTION run_meta_scheme()
# ----------------------------------------------------------------------------
# Runs scheme-specific meta-regressions, including simple and multi-category 
# versions.
run_meta_scheme <- function(res_service, scheme) {
  df       <- build_meta_data_scheme(res_service, scheme)
  m_simple <- lm(estimate_pct ~ is_shoppable, data = df, weights = weight)
  m_three  <- if (any(df$is_nonshoppable)) {
    lm(estimate_pct ~ is_shoppable + is_nonshoppable, data = df, weights = weight)
  } else NULL
  m_split  <- if (any(df$is_nonshoppable_mri) & any(df$is_nonshoppable_proc)) {
    lm(estimate_pct ~ is_shoppable + is_nonshoppable_mri + is_nonshoppable_proc,
       data = df, weights = weight)
  } else NULL
  m_simple_uw <- lm(estimate_pct ~ is_shoppable, data = df)
  m_three_uw  <- if (any(df$is_nonshoppable)) {
    lm(estimate_pct ~ is_shoppable + is_nonshoppable, data = df)
  } else NULL
  list(label = scheme$label, data = df, simple = m_simple, three_way = m_three,
       split = m_split, simple_uw = m_simple_uw, three_uw = m_three_uw)
}

meta_scheme_results <- lapply(shoppability_schemes, function(scheme) {
  run_meta_scheme(res_county, scheme)
})

# Print all scheme results
for (nm in names(meta_scheme_results)) {
  res <- meta_scheme_results[[nm]]
  df  <- res$data
  cat("\n", strrep("=", 65), "\n")
  cat(" SCHEME:", res$label, "\n")
  cat(" Shoppable:", paste(df$SERVICE_GROUP[df$is_shoppable], collapse = ", "), "\n")
  cat("\n--- Model A: Simple shoppable binary (WLS) ---\n"); print(summary(res$simple))
  if (!is.null(res$three_way)) {
    cat("\n--- Model B: Three-way (WLS) ---\n"); print(summary(res$three_way))
  }
  if (!is.null(res$split)) {
    cat("\n--- Model B2: Split Non-Shoppable (WLS) ---\n"); print(summary(res$split))
  }
  cat(strrep("=", 65), "\n")
}

# Summary table: shoppable coefficient across all schemes
cat("\n=== SHOPPABLE COEFFICIENT ACROSS ALL SCHEMES ===\n")
cat(sprintf("%-50s  %8s  %8s  %8s  %6s\n", "Scheme","Coef","SE","p-value","Sig"))
cat(strrep("-", 82), "\n")
for (nm in names(meta_scheme_results)) {
  res <- meta_scheme_results[[nm]]
  ct  <- coef(summary(res$simple))["is_shoppableTRUE", ]
  sig <- ifelse(ct["Pr(>|t|)"] < 0.01, "***", ifelse(ct["Pr(>|t|)"] < 0.05, "**",
                                                     ifelse(ct["Pr(>|t|)"] < 0.10, "*", "")))
  cat(sprintf("%-50s  %8.3f  %8.3f  %8.4f  %6s\n",
              res$label, ct["Estimate"], ct["Std. Error"], ct["Pr(>|t|)"], sig))
}

# Saturation meta-regression across bins and schemes
cat(sprintf("\n=== SATURATION META-REGRESSION ACROSS %d SCHEMES ===\n", n_shop_schemes))
sat_meta_results <- list()
for (bin_nm in c("low","mid","high")) {
  bin_data     <- switch(bin_nm, low = sat_service_low, mid = sat_service_mid, high = sat_service_high)
  bin_label_nm <- switch(bin_nm, low = low_label, mid = mid_label, high = high_label)
  if (is.null(bin_data) || nrow(bin_data) < 5) next
  for (scheme_nm in names(shoppability_schemes)) {
    scheme   <- shoppability_schemes[[scheme_nm]]
    df_meta  <- build_meta_data_scheme(bin_data, scheme)
    if (nrow(df_meta) < 5) next
    m_simple <- tryCatch(lm(estimate_pct ~ is_shoppable, data = df_meta, weights = weight),
                         error = function(e) NULL)
    m_three  <- if (any(df_meta$is_nonshoppable)) {
      tryCatch(lm(estimate_pct ~ is_shoppable + is_nonshoppable, data = df_meta, weights = weight),
               error = function(e) NULL)
    } else NULL
    extract_sat_coefs <- function(model, model_label) {
      if (is.null(model)) return(NULL)
      ct <- as.data.frame(coef(summary(model)))
      ct$term <- rownames(ct); ct$model <- model_label
      ct$scheme <- scheme$label; ct$bin <- bin_label_nm
      ct$r_squared <- summary(model)$r.squared; rownames(ct) <- NULL
      ct %>% select(bin, scheme, model, term,
                    estimate = Estimate, se = `Std. Error`,
                    t_stat = `t value`, p_value = `Pr(>|t|)`, r_squared)
    }
    key <- paste0(bin_nm, "_", scheme_nm)
    sat_meta_results[[paste0(key,"_simple")]] <- extract_sat_coefs(m_simple, "Simple shoppable (WLS)")
    sat_meta_results[[paste0(key,"_three")]]  <- extract_sat_coefs(m_three,  "Three-way (WLS)")
  }
}
sat_meta_df <- do.call(rbind, sat_meta_results) %>%
  mutate(sig = ifelse(p_value < 0.01, "***", ifelse(p_value < 0.05, "**",
                                                    ifelse(p_value < 0.10, "*", ""))),
         across(where(is.numeric), ~round(.x, 4)))

cat("\n--- Summary: Shoppable coefficient by bin x scheme ---\n")
sat_meta_df %>%
  filter(model == "Simple shoppable (WLS)", term == "is_shoppableTRUE") %>%
  arrange(bin, p_value) %>%
  rowwise() %>%
  do({ cat(sprintf("%-20s  %-50s  %8.3f  %8.3f  %8.4f  %4s\n",
                   substr(.$bin,1,20), .$scheme, .$estimate, .$se, .$p_value, .$sig)); data.frame() })

# Saturation figure (Theory V2, primary scheme)
sat_plot_data <- sat_meta_df %>%
  filter(scheme == "Theory-Based V2 (MRI as nonshoppable)",
         term %in% c("is_shoppableTRUE","is_nonshoppableTRUE"),
         model %in% c("Simple shoppable (WLS)","Three-way (WLS)")) %>%
  filter((term == "is_shoppableTRUE"    & model == "Simple shoppable (WLS)") |
           (term == "is_nonshoppableTRUE" & model == "Three-way (WLS)")) %>%
  mutate(
    service_type = factor(ifelse(term=="is_shoppableTRUE","Shoppable","Non-Shoppable"),
                          levels = c("Shoppable","Non-Shoppable")),
    bin_label = factor(bin,
                       levels = c(low_label, mid_label, high_label),
                       labels = c("Low\n(1\u20133 posters)","Moderate\n(4\u20138 posters)",
                                  "High\n(9+ posters)")),
    ci_lo = estimate - 1.96*se, ci_hi = estimate + 1.96*se,
    sig_label = case_when(p_value < 0.01 ~ "***", p_value < 0.05 ~ "**",
                          p_value < 0.10 ~ "*", TRUE ~ ""),
    label_y = ifelse(estimate >= 0, ci_hi + 1.8, ci_lo - 1.8)
  )

# ----------------------------------------------------------------------------
# OUTPUT_FIGURE_04_SATURATION_SHOPPABILITY
# ----------------------------------------------------------------------------
# Figure object fig_saturation_shoppability: saturation-bin effects.
fig_saturation_shoppability <- ggplot(sat_plot_data,
                                      aes(x = bin_label, y = estimate,
                                          color = service_type, fill = service_type, group = service_type)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4, color = "grey50") +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.10, color = NA) +
  geom_line(linewidth = 0.9) + geom_point(size = 3.5) +
  geom_text(aes(y = label_y, label = sig_label), size = 3.5, show.legend = FALSE) +
  scale_color_manual(values = c("Shoppable" = FSU_GARNET, "Non-Shoppable" = FSU_GOLD), name = NULL) +
  scale_fill_manual(values  = c("Shoppable" = FSU_GARNET, "Non-Shoppable" = FSU_GOLD), name = NULL) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(x = "Market Saturation Bin", y = "Meta-regression coefficient (pp)",
       title = "Price Effects by Market Saturation and Service Shoppability",
       subtitle = "Theory V2 scheme; county markets; 9m instrument; two-way clustering.",
       caption = "*** p<0.01  ** p<0.05  * p<0.10. Shaded bands: 95% CIs.") +
  theme_bw(base_size = 11) +
  theme(legend.position = c(0.87, 0.88), panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))
print(fig_saturation_shoppability)

# Saturation robustness across all n_shop_schemes schemes (currently 18)
sat_robust_data <- sat_meta_df %>%
  filter(term %in% c("is_shoppableTRUE","is_nonshoppableTRUE"),
         model %in% c("Simple shoppable (WLS)","Three-way (WLS)")) %>%
  filter((term == "is_shoppableTRUE"    & model == "Simple shoppable (WLS)") |
           (term == "is_nonshoppableTRUE" & model == "Three-way (WLS)")) %>%
  mutate(service_type = factor(ifelse(term=="is_shoppableTRUE","Shoppable","Non-Shoppable"),
                               levels = c("Shoppable","Non-Shoppable"),
                               labels = c("Panel A: Shoppable","Panel B: Non-Shoppable")),
         bin_label = factor(bin,
                            levels = c(low_label, mid_label, high_label),
                            labels = c("Low\n(1\u20133)","Moderate\n(4\u20138)","High\n(9+)")),
         is_theory_v2 = scheme == "Theory-Based V2 (MRI as nonshoppable)")

# ----------------------------------------------------------------------------
# OUTPUT_FIGURE_05_SATURATION_ROBUST
# ----------------------------------------------------------------------------
# Figure object fig_saturation_robust: saturation robustness across 
# classification schemes.
fig_saturation_robust <- ggplot(sat_robust_data,
                                aes(x = bin_label, y = estimate, group = scheme,
                                    color = is_theory_v2, alpha = is_theory_v2, linewidth = is_theory_v2)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4, color = "grey50") +
  geom_line() + geom_point(size = 2) +
  facet_wrap(~service_type, scales = "free_y") +
  scale_color_manual(values = c("TRUE" = FSU_GARNET, "FALSE" = "grey60"), guide = "none") +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.35), guide = "none") +
  scale_linewidth_manual(values = c("TRUE" = 1.2, "FALSE" = 0.5), guide = "none") +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(x = "Market Saturation Bin", y = "Meta-regression coefficient (pp)",
       title = "Saturation Gradient Across All Classification Schemes",
       subtitle = "Garnet = Theory V2 (primary). Grey = alternative schemes. County 9m.",
       caption = "Simple WLS (shoppable) and three-way WLS (non-shoppable).") +
  theme_bw(base_size = 11) + theme(strip.text = element_text(face = "bold", size = 11))
print(fig_saturation_robust)

# Schemes figure: meta-regression forest
schemes_summary <- do.call(rbind, lapply(names(meta_scheme_results), function(nm) {
  res <- meta_scheme_results[[nm]]
  ct  <- coef(summary(res$simple))["is_shoppableTRUE", ]
  data.frame(label = res$label, coef = ct["Estimate"], se = ct["Std. Error"],
             pval = ct["Pr(>|t|)"], row.names = NULL)
})) %>%
  mutate(ci_lo = coef - 1.96*se, ci_hi = coef + 1.96*se, sig = pval < 0.05,
         label_short = case_when(
           grepl("CMS Rule",       label) ~ "CMS Rule Definition",
           grepl("Theory-Based V2",label) ~ "Theory V2 (Primary)",
           grepl("Theory-Based",   label) ~ "Theory-Based",
           grepl("Broad",          label) ~ "Broad Shoppable",
           grepl("Imaging",        label) ~ "Imaging vs Procedural",
           grepl("Split",          label) ~ "Split Non-Shoppable",
           grepl("Statutory",      label) ~ "CMS Statutory List",
           grepl("High vs Low",    label) ~ "High vs Low (Within Modality)",
           grepl("Except Angio",   label) ~ "CT-Inclusive Ex. Angio",
           grepl("CT-Inclusive",   label) ~ "CT-Inclusive (All CT)",
           # -- Short labels for the 8 alternative frameworks added from
           # Danny's classification-scheme document --
           grepl("Core Operations",     label) ~ "Alt: Core Operations",
           grepl("CMS Legal",           label) ~ "Alt: CMS Legal",
           grepl("Upfront Cash-Market", label) ~ "Alt: Cash-Market (MDsave)",
           grepl("Geographic",          label) ~ "Alt: Geographic Access",
           grepl("Urgency",             label) ~ "Alt: Diagnostic Urgency",
           grepl("Staffing",            label) ~ "Alt: Staffing/Specialist",
           grepl("Liability",           label) ~ "Alt: Liability/Risk",
           grepl("No Surprises",        label) ~ "Alt: No Surprises Act",
           TRUE ~ label),
         label_short = factor(label_short, levels = label_short[order(coef)]))

# Refresh results_meta_regression_all_schemes.csv from the live schemes_summary
# object, so panel_b (which reads this file back in further down) reflects the
# current scheme list instead of a stale snapshot from before schemes were
# added or removed.
schemes_summary_csv <- schemes_summary %>%
  transmute(
    scheme   = label,
    model    = "Simple (shoppable)",
    estimate = coef,
    se       = se,
    p_value  = pval,
    stars    = case_when(pval < 0.01 ~ "***", pval < 0.05 ~ "**",
                         pval < 0.10 ~ "*", TRUE ~ "")
  )
write.csv(schemes_summary_csv, "results_meta_regression_all_schemes.csv", row.names = FALSE)

ggplot(schemes_summary, aes(x = coef, y = label_short, color = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.25, linewidth = 0.6) +
  geom_point(size = 3) +
  geom_point(data = filter(schemes_summary, grepl("Primary", label_short)), shape = 18, size = 5) +
  scale_color_manual(values = c("TRUE" = FSU_GARNET, "FALSE" = FSU_GOLD),
                     labels = c("TRUE" = "p < 0.05", "FALSE" = "p \u2265 0.05"), name = NULL) +
  scale_x_continuous(labels = function(x) sprintf("%.1f%%", x)) +
  labs(x = "Shoppable Coefficient (pp)", y = NULL,
       title = paste0("Robustness Across ", n_shop_schemes, " Classification Schemes"),
       subtitle = "Simple WLS meta-regression; whiskers = 95% CIs; diamond = primary spec") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", panel.grid.major.y = element_blank())
# ----------------------------------------------------------------------------
# OUTPUT_FIGURE_METAREG_SCHEMES_FILE
# ----------------------------------------------------------------------------
# Saved PDF: fig_metareg_schemes.pdf. Shows meta-regression results across 
# shoppability schemes.
# NOTE: height increased from 5 to 8 -- with 18 rows (up from the original 10)
# a height of 5 crowds the row labels/whiskers together. Re-check spacing
# after running; adjust further if labels still overlap.
ggsave("fig_metareg_schemes.pdf", width = 8, height = 8, device = cairo_pdf)

# modelsummary for Table 6
library(modelsummary)
m1 <- meta_scheme_results[["scheme_theory_v2"]]$simple
m2 <- meta_scheme_results[["scheme_theory_v2"]]$three_way
m3 <- meta_scheme_results[["scheme_split_nonshop"]]$split
modelsummary(
  list("Simple" = m1, "Three-Way" = m2, "Split Non-Shoppable" = m3),
  estimate = "{estimate}{stars}", statistic = "({std.error})",
  coef_map = c("is_shoppableTRUE" = "Shoppable",
               "is_nonshoppableTRUE" = "Non-Shoppable",
               "is_nonshoppable_mriTRUE" = "Non-Shoppable (MRI)",
               "is_nonshoppable_procTRUE" = "Non-Shoppable (Procedural)",
               "(Intercept)" = "Intercept (Baseline)"),
  gof_map = c("nobs","r.squared","adj.r.squared"),
  stars = c("*"=0.1,"**"=0.05,"***"=0.01), output = "latex"
)





######## Section 11: Heterogeneity — Market Characteristics ###########
# ----------------------------------------------------------------------------
# Market-characteristic heterogeneity: splits markets by ownership, payer mix,
#  poverty, education, insurer concentration, market size, race/ethnicity, and
#  related dimensions.

merged_final <- fread("../Data/data_clean/HHI_CBSA_health_insurance_market_competition.csv")

df_het <- df_iv_county %>%
  left_join(merged_final, by = "cbsacode") %>%
  filter(!is.na(`TOTAL HHI`)) %>%
  mutate(
    hhi_1000       = `TOTAL HHI` / 1000,
    top2_share     = (`Share (%)` + `Share 2 (%)`) / 100,
    high_forprofit = as.integer(Share_ForProfit >= median(Share_ForProfit, na.rm=TRUE)),
    high_college   = as.integer(college_share   >= median(college_share,   na.rm=TRUE)),
    high_elderly   = as.integer(age65plus_share  >= median(age65plus_share, na.rm=TRUE)),
    large_market   = as.integer(population       >= median(population,      na.rm=TRUE)),
    high_poverty   = as.integer(poverty_rate     >= median(poverty_rate,    na.rm=TRUE))
  )

# ----------------------------------------------------------------------------
# FUNCTION run_iv_sub()
# ----------------------------------------------------------------------------
# Runs a pooled IV model on a supplied heterogeneity split/subsample.
run_iv_sub <- function(df, label, instrument, cluster_var, min_n = 50) {
  if (nrow(df) < min_n) return(NULL)
  instr <- df[[instrument]]
  if (is.null(instr) || length(unique(na.omit(instr))) < 2) return(NULL)
  fit <- tryCatch(
    feols(as.formula(paste0("ln_median_price ~ 1 | market_id + post_month | n_prior_posters ~ ",
                            instrument)),
          data    = df,
          cluster = as.formula(paste0("~", paste(cluster_var, collapse = " + ")))),
    error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  ct  <- coeftable(fit)
  if (!"fit_n_prior_posters" %in% rownames(ct)) return(NULL)
  z   <- qnorm(0.975)
  est <- ct["fit_n_prior_posters","Estimate"]
  se  <- ct["fit_n_prior_posters","Std. Error"]
  wf  <- tryCatch(fitstat(fit, "ivwald")[["ivwald1::n_prior_posters"]]$stat,
                  error = function(e) NA_real_)
  data.frame(group = label, estimate = est, se = se,
             pvalue = ct["fit_n_prior_posters","Pr(>|t|)"],
             ci_low = est - z*se, ci_high = est + z*se,
             pct = 100*(exp(est)-1), pct_low = 100*(exp(est-z*se)-1),
             pct_high = 100*(exp(est+z*se)-1),
             n_obs = nrow(df), wald_f = wf, stringsAsFactors = FALSE)
}

# ----------------------------------------------------------------------------
# FUNCTION run_all_results_on_split()
# ----------------------------------------------------------------------------
# For one heterogeneity split, runs pooled IV, service-level IV, SHOP2 
# results, and scheme-specific meta-regressions.
run_all_results_on_split <- function(df_split, label,
                                     instrument  = "system_peer_pressure_county_9m",
                                     cluster_var = c("County_State","post_month"),
                                     min_obs     = 100) {
  pooled    <- run_iv_sub(df_split, label, instrument, cluster_var)
  res_svc   <- run_iv_by_service(df_split, "ln_median_price", instrument, cluster_var, min_obs)
  res_shop2 <- run_iv_by_large_group(df_split, "SERVICE_GROUP_SHOP2", "ln_median_price",
                                     instrument, cluster_var, min_obs)
  meta_list <- list()
  if (!is.null(res_svc) && nrow(res_svc) >= 5) {
    for (scheme_nm in names(shoppability_schemes)) {
      scheme  <- shoppability_schemes[[scheme_nm]]
      df_meta <- build_meta_data_scheme(res_svc, scheme)
      if (nrow(df_meta) < 5) next
      m_simple <- tryCatch(lm(estimate_pct ~ is_shoppable, data=df_meta, weights=weight),
                           error = function(e) NULL)
      m_three  <- if (any(df_meta$is_nonshoppable)) {
        tryCatch(lm(estimate_pct ~ is_shoppable + is_nonshoppable, data=df_meta, weights=weight),
                 error = function(e) NULL)
      } else NULL
      extract_meta <- function(model, model_label) {
        if (is.null(model)) return(NULL)
        ct <- as.data.frame(coef(summary(model)))
        ct$term <- rownames(ct); ct$model <- model_label
        ct$scheme <- scheme$label; ct$split <- label
        ct$r_squared <- summary(model)$r.squared; rownames(ct) <- NULL
        ct %>% select(split, scheme, model, term,
                      estimate = Estimate, se = `Std. Error`,
                      t_stat = `t value`, p_value = `Pr(>|t|)`, r_squared)
      }
      meta_list[[paste0(scheme_nm,"_simple")]] <- extract_meta(m_simple, "Simple shoppable (WLS)")
      meta_list[[paste0(scheme_nm,"_three")]]  <- extract_meta(m_three,  "Three-way (WLS)")
    }
  }
  meta_df <- if (length(meta_list) > 0) {
    do.call(rbind, meta_list) %>%
      mutate(sig = ifelse(p_value<0.01,"***",ifelse(p_value<0.05,"**",ifelse(p_value<0.10,"*",""))),
             across(where(is.numeric), ~round(.x, 4)))
  } else NULL
  list(label = label, pooled = pooled, service = res_svc, shop2 = res_shop2, meta = meta_df)
}

split_defs <- list(
  list(label = "Low For-Profit",      data = df_het %>% filter(high_forprofit == 0)),
  list(label = "High For-Profit",     data = df_het %>% filter(high_forprofit == 1)),
  list(label = "Low College",         data = df_het %>% filter(high_college   == 0)),
  list(label = "High College",        data = df_het %>% filter(high_college   == 1)),
  list(label = "Low HHI",             data = df_het %>% filter(hhi_1000 <  median(hhi_1000,na.rm=TRUE))),
  list(label = "High HHI",            data = df_het %>% filter(hhi_1000 >= median(hhi_1000,na.rm=TRUE))),
  list(label = "Low 65+ Share",       data = df_het %>% filter(high_elderly   == 0)),
  list(label = "High 65+ Share",      data = df_het %>% filter(high_elderly   == 1)),
  list(label = "Small Market",        data = df_het %>% filter(large_market   == 0)),
  list(label = "Large Market",        data = df_het %>% filter(large_market   == 1)),
  list(label = "Low Poverty",         data = df_het %>% filter(high_poverty   == 0)),
  list(label = "High Poverty",        data = df_het %>% filter(high_poverty   == 1)),
  list(label = "Low Uninsured",
       data = df_het %>% filter(uninsured_rate <  median(uninsured_rate,na.rm=TRUE))),
  list(label = "High Uninsured",
       data = df_het %>% filter(uninsured_rate >= median(uninsured_rate,na.rm=TRUE))),
  list(label = "Low Medicaid",
       data = df_het %>% filter(medicaid_share <  median(medicaid_share,na.rm=TRUE))),
  list(label = "High Medicaid",
       data = df_het %>% filter(medicaid_share >= median(medicaid_share,na.rm=TRUE))),
  list(label = "Low Black Share",
       data = df_het %>% filter(black_share <  median(black_share,na.rm=TRUE))),
  list(label = "High Black Share",
       data = df_het %>% filter(black_share >= median(black_share,na.rm=TRUE))),
  list(label = "Low Hispanic Share",
       data = df_het %>% filter(hispanic_share <  median(hispanic_share,na.rm=TRUE))),
  list(label = "High Hispanic Share",
       data = df_het %>% filter(hispanic_share >= median(hispanic_share,na.rm=TRUE)))
)

cat("\nRunning all heterogeneity splits (county 9m, two-way clustering)...\n")
all_het_results <- lapply(split_defs, function(s) {
  cat(sprintf("  Processing: %s (N=%d)\n", s$label, nrow(s$data)))
  run_all_results_on_split(s$data, s$label)
})
names(all_het_results) <- sapply(split_defs, `[[`, "label")




# ============================================================
# MARKET SIZE QUANTILE ANALYSIS
# ============================================================

cat("\n=== MARKET SIZE QUANTILE ANALYSIS: corrected bins ===\n")

N_BINS <- 5

# Step 1: define ONE population value per CBSA
# Use median population within CBSA in case population is repeated or slightly inconsistent.
market_size_bins <- df_het %>%
  filter(!is.na(cbsacode), !is.na(population)) %>%
  group_by(cbsacode) %>%
  summarise(
    market_population = median(population, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(market_population) %>%
  mutate(
    size_quantile = ntile(market_population, N_BINS),
    size_bin = factor(
      size_quantile,
      levels = 1:N_BINS,
      labels = paste0("Q", 1:N_BINS)
    )
  )

# Step 2: check that bins are ordered correctly
market_size_bin_summary <- market_size_bins %>%
  group_by(size_bin) %>%
  summarise(
    n_markets = n(),
    min_pop = min(market_population, na.rm = TRUE),
    median_pop = median(market_population, na.rm = TRUE),
    max_pop = max(market_population, na.rm = TRUE),
    .groups = "drop"
  )

print(market_size_bin_summary)

# Step 3: join bins back to full data
df_het_sizeq <- df_het %>%
  left_join(
    market_size_bins %>%
      select(cbsacode, market_population, size_quantile, size_bin),
    by = "cbsacode"
  ) %>%
  filter(!is.na(size_bin))

# Optional safety check: each CBSA should have exactly one bin
check_bins <- df_het_sizeq %>%
  distinct(cbsacode, size_bin) %>%
  count(cbsacode) %>%
  filter(n > 1)

print(check_bins)

# If check_bins is empty, good.

sizeq_splits <- lapply(levels(df_het_sizeq$size_bin), function(q) {
  list(
    label = paste0("Market Size ", q),
    data = df_het_sizeq %>% filter(size_bin == q)
  )
})

sizeq_results <- lapply(sizeq_splits, function(s) {
  cat(sprintf("  Processing: %s (N=%d)\n", s$label, nrow(s$data)))
  run_all_results_on_split(s$data, s$label)
})

names(sizeq_results) <- sapply(sizeq_splits, `[[`, "label")

sizeq_meta_summary <- do.call(rbind, lapply(sizeq_results, function(r) {
  if (is.null(r$meta)) return(NULL)
  
  r$meta %>%
    filter(
      scheme == "Theory-Based V2 (MRI as nonshoppable)",
      model  == "Simple shoppable (WLS)",
      term   == "is_shoppableTRUE"
    ) %>%
    select(
      split,
      scheme,
      estimate,
      se,
      p_value,
      sig,
      r_squared
    )
}))

sizeq_pooled_summary <- do.call(rbind, lapply(sizeq_results, function(r) {
  if (is.null(r$pooled)) return(NULL)
  
  r$pooled %>%
    transmute(
      split = group,
      wald_f,
      n_obs
    )
}))

sizeq_plot_data <- sizeq_meta_summary %>%
  left_join(sizeq_pooled_summary, by = "split") %>%
  mutate(
    size_quantile = as.integer(gsub("Market Size Q", "", split)),
    size_bin = factor(
      paste0("Q", size_quantile),
      levels = paste0("Q", 1:N_BINS)
    ),
    ci_lo = estimate - 1.96 * se,
    ci_hi = estimate + 1.96 * se,
    weak_fs = !is.na(wald_f) & wald_f < 10,
    sig_label = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      p_value < 0.10  ~ "\u2020",
      TRUE ~ ""
    )
  ) %>%
  left_join(
    market_size_bin_summary %>%
      select(size_bin, n_markets, median_pop),
    by = "size_bin"
  )

print(sizeq_plot_data)



# ------------------------------------------------------------
# Figure: shoppability gradient across market-size quantiles
# ------------------------------------------------------------

# ----------------------------------------------------------------------------
# OUTPUT_FIGURE_MARKET_SIZE_QUANTILES
# ----------------------------------------------------------------------------
# Figure object fig_market_size_quantiles: market-size quantile heterogeneity.
fig_market_size_quantiles <- ggplot(
  sizeq_plot_data,
  aes(x = size_bin, y = estimate, group = 1)
) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.45,
    color = "grey55"
  ) +
  
  geom_line(
    linewidth = 0.85,
    color = FSU_GARNET,
    alpha = 0.85
  ) +
  
  geom_errorbar(
    aes(ymin = ci_lo, ymax = ci_hi),
    width = 0.10,
    linewidth = 0.75,
    color = FSU_GARNET
  ) +
  
  geom_point(
    aes(shape = weak_fs),
    size = 3.6,
    stroke = 1.0,
    color = FSU_GARNET
  ) +
  
  geom_text(
    aes(
      label = sig_label,
      y = ifelse(estimate < 0, ci_lo - 0.08 * diff(range(c(ci_lo, ci_hi), na.rm = TRUE)),
                 ci_hi + 0.08 * diff(range(c(ci_lo, ci_hi), na.rm = TRUE)))
    ),
    size = 3.5,
    fontface = "bold",
    color = FSU_GARNET
  ) +
  
  scale_shape_manual(
    values = c("FALSE" = 16, "TRUE" = 1),
    guide = "none"
  ) +
  
  scale_y_continuous(
    labels = function(x) paste0(x, " pp")
  ) +
  
  labs(
    x = "Market-size quantile",
    y = "Shoppability coefficient\n(pp per additional prior poster)",
    title = "Shoppability Gradient by Market Size Quantile",
    subtitle = "Theory V2 scheme, simple WLS meta-regression; smaller markets show the steepest gradients",
    caption = paste0(
      "Notes: Quantiles are defined at the CBSA level using population. Points show the shoppability coefficient ",
      "estimated from service-level IV estimates within each market-size quantile. Bars show 95% confidence intervals. ",
      "Open circles indicate pooled first-stage Wald F < 10. *** p<0.01, ** p<0.05, * p<0.10, \u2020 p<0.10."
    )
  ) +
  
  theme_bw(base_size = 11) +
  
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9.5, color = "grey30"),
    plot.caption = element_text(
      size = 7.5,
      color = "grey40",
      hjust = 0,
      lineheight = 1.25
    ),
    axis.title.x = element_text(size = 9.5),
    axis.title.y = element_text(size = 9.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

# ----------------------------------------------------------------------------
# OUTPUT_FIGURE_MARKET_SIZE_QUANTILES_PRINT
# ----------------------------------------------------------------------------
# Prints fig_market_size_quantiles.
print(fig_market_size_quantiles)



# Black x Medicaid/Poverty 2x2
cat("\n=== BLACK x MEDICAID / POVERTY INTERACTION ANALYSIS ===\n")
black_med    <- median(df_het$black_share,    na.rm=TRUE)
medicaid_med <- median(df_het$medicaid_share, na.rm=TRUE)
poverty_med  <- median(df_het$poverty_rate,   na.rm=TRUE)

bm_splits <- list(
  list(label="Low Black, Low Medicaid",
       data=df_het%>%filter(black_share<black_med,medicaid_share<medicaid_med)),
  list(label="High Black, Low Medicaid",
       data=df_het%>%filter(black_share>=black_med,medicaid_share<medicaid_med)),
  list(label="Low Black, High Medicaid",
       data=df_het%>%filter(black_share<black_med,medicaid_share>=medicaid_med)),
  list(label="High Black, High Medicaid",
       data=df_het%>%filter(black_share>=black_med,medicaid_share>=medicaid_med)),
  list(label="Low Black, Low Poverty",
       data=df_het%>%filter(black_share<black_med,poverty_rate<poverty_med)),
  list(label="High Black, Low Poverty",
       data=df_het%>%filter(black_share>=black_med,poverty_rate<poverty_med)),
  list(label="Low Black, High Poverty",
       data=df_het%>%filter(black_share<black_med,poverty_rate>=poverty_med)),
  list(label="High Black, High Poverty",
       data=df_het%>%filter(black_share>=black_med,poverty_rate>=poverty_med))
)

bm_results <- lapply(bm_splits, function(s) {
  cat(sprintf("  Processing: %s (N=%d)\n", s$label, nrow(s$data)))
  run_all_results_on_split(s$data, s$label)
})
names(bm_results) <- sapply(bm_splits, `[[`, "label")

print(bm_results)


# ============================================================
# EDUCATION x POVERTY 2x2 INTERACTION ANALYSIS
# ============================================================

cat("\n=== EDUCATION x POVERTY INTERACTION ANALYSIS ===\n")

college_med <- median(df_het$college_share, na.rm = TRUE)
poverty_med <- median(df_het$poverty_rate, na.rm = TRUE)

ep_splits <- list(
  list(
    label = "Low College, High Poverty",
    data  = df_het %>%
      filter(college_share < college_med,
             poverty_rate >= poverty_med)
  ),
  list(
    label = "Low College, Low Poverty",
    data  = df_het %>%
      filter(college_share < college_med,
             poverty_rate < poverty_med)
  ),
  list(
    label = "High College, High Poverty",
    data  = df_het %>%
      filter(college_share >= college_med,
             poverty_rate >= poverty_med)
  ),
  list(
    label = "High College, Low Poverty",
    data  = df_het %>%
      filter(college_share >= college_med,
             poverty_rate < poverty_med)
  )
)

ep_results <- lapply(ep_splits, function(s) {
  cat(sprintf("  Processing: %s (N=%d)\n", s$label, nrow(s$data)))
  run_all_results_on_split(s$data, s$label)
})

names(ep_results) <- sapply(ep_splits, `[[`, "label")

# Extract Theory V2 simple WLS meta-regression estimates
ep_meta_summary <- do.call(rbind, lapply(ep_results, function(r) {
  if (is.null(r$meta)) return(NULL)
  
  r$meta %>%
    filter(
      scheme == "Theory-Based V2 (MRI as nonshoppable)",
      model  == "Simple shoppable (WLS)",
      term   == "is_shoppableTRUE"
    ) %>%
    select(
      split,
      scheme,
      estimate,
      se,
      p_value,
      sig,
      r_squared
    )
}))

# Extract pooled first-stage Wald F and N
ep_pooled_summary <- do.call(rbind, lapply(ep_results, function(r) {
  if (is.null(r$pooled)) return(NULL)
  
  r$pooled %>%
    transmute(
      split = group,
      wald_f,
      n_obs
    )
}))

# Final Education x Poverty table data
ep_2x2_table <- ep_meta_summary %>%
  left_join(ep_pooled_summary, by = "split") %>%
  mutate(
    ci_lo = estimate - 1.96 * se,
    ci_hi = estimate + 1.96 * se,
    
    college_group = ifelse(grepl("^High College", split),
                           "High College Share",
                           "Low College Share"),
    
    poverty_group = ifelse(grepl("High Poverty", split),
                           "High Poverty",
                           "Low Poverty"),
    
    coef_se = sprintf("%.3f%s (%.3f)", estimate, sig, se),
    r2_fmt  = sprintf("%.3f", r_squared),
    f_fmt   = ifelse(
      is.na(wald_f),
      "",
      ifelse(wald_f < 10,
             sprintf("%.1f†", wald_f),
             sprintf("%.1f", wald_f))
    ),
    n_fmt = format(n_obs, big.mark = ",")
  ) %>%
  arrange(
    factor(college_group, levels = c("Low College Share", "High College Share")),
    factor(poverty_group, levels = c("High Poverty", "Low Poverty"))
  )

print(ep_2x2_table)


# Summary tables
dimension_map <- c(
  "Low For-Profit"="Hospital Ownership","High For-Profit"="Hospital Ownership",
  "Low College"="Education","High College"="Education",
  "Low HHI"="Insurer Concentration","High HHI"="Insurer Concentration",
  "Low 65+ Share"="Medicare Exposure","High 65+ Share"="Medicare Exposure",
  "Small Market"="Market Size","Large Market"="Market Size",
  "Low Poverty"="Poverty Rate","High Poverty"="Poverty Rate",
  "Low Uninsured"="Uninsured Rate","High Uninsured"="Uninsured Rate",
  "Low Medicaid"="Medicaid Share","High Medicaid"="Medicaid Share",
  "Low Black Share"="Black Share","High Black Share"="Black Share",
  "Low Hispanic Share"="Hispanic Share","High Hispanic Share"="Hispanic Share"
)

pooled_summary <- do.call(rbind, lapply(all_het_results, function(r) {
  if (is.null(r$pooled)) return(NULL)
  r$pooled %>%
    mutate(dimension = dimension_map[group], sig = add_stars(pvalue),
           fmt = sprintf("%.3f%s (%.3f)", pct, sig, se*100)) %>%
    select(dimension, group, fmt, wald_f, n_obs)
})) %>% mutate(across(where(is.numeric), ~round(.x, 1)))

print(pooled_summary)
# ----------------------------------------------------------------------------
# EXPORT_CSV_het_pooled_summary
# ----------------------------------------------------------------------------
# Exports pooled heterogeneity summary to het_pooled_summary.csv.
write.csv(pooled_summary, "het_pooled_summary.csv", row.names=FALSE)

shop2_summary <- do.call(rbind, lapply(names(all_het_results), function(nm) {
  r <- all_het_results[[nm]]
  if (is.null(r$shop2)) return(NULL)
  r$shop2 %>%
    mutate(split = nm, dimension = dimension_map[nm], sig = add_stars(pval),
           fmt = sprintf("%.3f%s (%.3f)", estimate_pct, sig, se_pct)) %>%
    select(dimension, split, group, fmt, fs_f, n_obs)
}))
# ----------------------------------------------------------------------------
# EXPORT_CSV_het_shop2_summary
# ----------------------------------------------------------------------------
# Exports SHOP2 heterogeneity summary to het_shop2_summary.csv.
write.csv(shop2_summary %>% mutate(across(where(is.numeric), ~round(.x,3))),
          "het_shop2_summary.csv", row.names=FALSE)

meta_summary <- do.call(rbind, lapply(all_het_results, function(r) {
  if (is.null(r$meta)) return(NULL)
  r$meta %>% filter(model=="Simple shoppable (WLS)", term=="is_shoppableTRUE") %>%
    select(split, scheme, estimate, se, p_value, sig, r_squared)
}))
bm_meta_summary <- do.call(rbind, lapply(bm_results, function(r) {
  if (is.null(r$meta)) return(NULL)
  r$meta %>% filter(model=="Simple shoppable (WLS)", term=="is_shoppableTRUE") %>%
    select(split, scheme, estimate, se, p_value, sig, r_squared)
}))
full_meta_summary <- bind_rows(meta_summary, bm_meta_summary) %>%
  mutate(dimension = coalesce(dimension_map[split], "Black \u00d7 Medicaid/Poverty"))
# ----------------------------------------------------------------------------
# EXPORT_CSV_het_meta_summary
# ----------------------------------------------------------------------------
# Exports meta-regression heterogeneity summary to het_meta_summary.csv.
write.csv(full_meta_summary %>% mutate(across(where(is.numeric), ~round(.x,4))),
          "het_meta_summary.csv", row.names=FALSE)

cat("\n--- Meta-regression shoppable coefficient (Theory V2 scheme) ---\n")
full_meta_summary %>%
  filter(scheme == "Theory-Based V2 (MRI as nonshoppable)") %>%
  arrange(dimension, split) %>%
  rowwise() %>%
  do({ cat(sprintf("%-25s  %-22s  %8.3f  %8.3f  %8.4f  %4s\n",
                   substr(.$dimension,1,25), substr(.$split,1,22),
                   .$estimate, .$se, .$p_value, .$sig)); data.frame() })
# ----------------------------------------------------------------------------
# EXPORT_CSV_het_meta_full_matrix
# ----------------------------------------------------------------------------
# Exports full heterogeneity/meta matrix to het_meta_full_matrix.csv.
write.csv(full_meta_summary %>% mutate(across(where(is.numeric),~round(.x,4))) %>%
            arrange(dimension, split, scheme), "het_meta_full_matrix.csv", row.names=FALSE)

# Heterogeneity forest plot
forest_data <- full_meta_summary %>%
  filter(scheme == "Theory-Based V2 (MRI as nonshoppable)",
         !dimension %in% c("Black \u00d7 Medicaid/Poverty")) %>%
  mutate(ci_lo = estimate - 1.96*se, ci_hi = estimate + 1.96*se,
         sig_group = case_when(p_value<0.01~"p < 0.01",p_value<0.05~"p < 0.05",
                               p_value<0.10~"p < 0.10",TRUE~"n.s."),
         mechanism = case_when(
           dimension %in% c("Hospital Ownership","Medicaid Share","Medicare Exposure") ~
             "Panel A: Competitive Incentives",
           dimension %in% c("Poverty Rate","Education","Uninsured Rate") ~
             "Panel B: Patient Capacity",
           dimension %in% c("Insurer Concentration","Market Size","Black Share","Hispanic Share") ~
             "Panel C: Market Structure"),
         weak_fs = split %in% c("Low For-Profit","Low Uninsured"),
         split_label = case_when(
           split=="High For-Profit"~"For-Profit: High",split=="Low For-Profit"~"For-Profit: Low \u2020",
           split=="High Medicaid"~"Medicaid Share: High",split=="Low Medicaid"~"Medicaid Share: Low",
           split=="High 65+ Share"~"Medicare (65+): High",split=="Low 65+ Share"~"Medicare (65+): Low",
           split=="High Poverty"~"Poverty Rate: High",split=="Low Poverty"~"Poverty Rate: Low",
           split=="High College"~"Education: High",split=="Low College"~"Education: Low",
           split=="High Uninsured"~"Uninsured Rate: High",split=="Low Uninsured"~"Uninsured Rate: Low \u2020",
           split=="High HHI"~"Insurer HHI: High",split=="Low HHI"~"Insurer HHI: Low",
           split=="Large Market"~"Market Size: Large",split=="Small Market"~"Market Size: Small",
           split=="High Black Share"~"Black Share: High",split=="Low Black Share"~"Black Share: Low",
           split=="High Hispanic Share"~"Hispanic Share: High",split=="Low Hispanic Share"~"Hispanic Share: Low",
           TRUE~split),
         dim_order = case_when(dimension=="Hospital Ownership"~1,dimension=="Medicaid Share"~2,
                               dimension=="Medicare Exposure"~3,dimension=="Poverty Rate"~4,
                               dimension=="Education"~5,dimension=="Uninsured Rate"~6,
                               dimension=="Insurer Concentration"~7,dimension=="Market Size"~8,
                               dimension=="Black Share"~9,dimension=="Hispanic Share"~10,TRUE~99),
         high_first = grepl("High|Large", split)) %>%
  filter(!is.na(mechanism)) %>%
  arrange(mechanism, dim_order, desc(high_first)) %>%
  mutate(split_label = factor(split_label, levels = rev(unique(split_label))))

# ----------------------------------------------------------------------------
# OUTPUT_FIGURE_06_HETEROGENEITY_FOREST
# ----------------------------------------------------------------------------
# Figure object fig_het_forest: main heterogeneity forest plot.
fig_het_forest <- ggplot(forest_data,
                         aes(x=estimate, y=split_label, color=sig_group, shape=weak_fs)) +
  geom_vline(xintercept=0, linetype="dashed", linewidth=0.4, color="grey50") +
  geom_errorbarh(aes(xmin=ci_lo, xmax=ci_hi), height=0.25, linewidth=0.65) +
  geom_point(size=3) +
  facet_wrap(~mechanism, ncol=1, scales="free_y", strip.position="top") +
  scale_color_manual(values=c("p < 0.01"=FSU_GARNET,"p < 0.05"="#B84040","p < 0.10"=FSU_GOLD,"n.s."="grey60"),
                     name="Significance") +
  scale_shape_manual(values=c("TRUE"=1,"FALSE"=16), guide="none") +
  scale_x_continuous(labels=function(x) paste0(x," pp")) +
  labs(x="Shoppability coefficient (pp per unit, 95% CI)", y=NULL,
       title="Shoppability Gradient by Market Characteristic",
       subtitle="Theory V2 scheme, simple WLS meta-regression. County 9m instrument.",
       caption="\u2020 F-stat subgroup may be weak; interpret with caution.") +
  theme_bw(base_size=11) +
  theme(strip.text=element_text(face="bold",size=10), panel.grid.major.y=element_blank(),
        plot.title=element_text(face="bold"), legend.position="bottom")
# ----------------------------------------------------------------------------
# OUTPUT_FIGURE_06_HETEROGENEITY_FOREST_PRINT
# ----------------------------------------------------------------------------
# Prints fig_het_forest.
print(fig_het_forest)

# Race x Medicaid equity figure
equity_data <- bm_meta_summary %>%
  filter(scheme=="Theory-Based V2 (MRI as nonshoppable)") %>%
  mutate(ci_lo=estimate-1.96*se, ci_hi=estimate+1.96*se,
         black_group=ifelse(grepl("^High Black",split),"High Black Share","Low Black Share"),
         panel=case_when(grepl("Medicaid",split)~"Panel A: \u00d7 Medicaid Share",
                         grepl("Poverty",split) ~"Panel B: \u00d7 Poverty Rate"),
         econ_level=case_when(grepl("High Medic|High Pover",split)~"High\n(more disadvantaged)",
                              grepl("Low Medic|Low Pover",split)  ~"Low\n(less disadvantaged)"),
         weak_fs=split=="High Black, High Medicaid",
         sig_label=case_when(p_value<0.001~"***",p_value<0.01~"**",p_value<0.05~"*",
                             p_value<0.10~"\u2020",TRUE~""),
         label_y=ifelse(estimate>=0, ci_hi+0.45, ci_lo-0.45),
         econ_level=factor(econ_level,levels=c("High\n(more disadvantaged)","Low\n(less disadvantaged)")))

# ----------------------------------------------------------------------------
# OUTPUT_FIGURE_07_RACE_MEDICAID_HEATMAP
# ----------------------------------------------------------------------------
# Figure object fig_race_2x2: race x Medicaid gradient heatmap.
fig_race_2x2 <- ggplot(equity_data,
                       aes(x=econ_level, y=estimate, color=black_group,
                           shape=weak_fs, group=black_group)) +
  geom_hline(yintercept=0, linetype="dashed", linewidth=0.4, color="grey50") +
  geom_line(linewidth=0.85, alpha=0.85) +
  geom_errorbar(aes(ymin=ci_lo, ymax=ci_hi), width=0.10, linewidth=0.7) +
  geom_point(size=3.5) +
  geom_text(aes(y=label_y, label=sig_label), size=3.5, show.legend=FALSE) +
  facet_wrap(~panel, ncol=2) +
  scale_color_manual(values=c("Low Black Share"=FSU_GARNET,"High Black Share"=FSU_GOLD),
                     name="Black population\nshare") +
  scale_shape_manual(values=c("TRUE"=1,"FALSE"=16), guide="none") +
  scale_y_continuous(labels=function(x) paste0(x," pp")) +
  labs(x=NULL, y="Shoppability coefficient (pp per unit)",
       title="Shoppability Gradient: Black Share \u00d7 Economic Disadvantage",
       subtitle="Theory V2 scheme, simple WLS meta-regression. County 9m. 95% CIs shown.",
       caption="Open circle = weak first stage. *** p<0.01, ** p<0.05, * p<0.05, \u2020 p<0.10.") +
  theme_bw(base_size=11) +
  theme(strip.text=element_text(face="bold",size=10), legend.position=c(0.13,0.20),
        legend.background=element_rect(fill="white",color="grey80",linewidth=0.3),
        plot.title=element_text(face="bold"))
# ----------------------------------------------------------------------------
# OUTPUT_FIGURE_07_RACE_MEDICAID_HEATMAP_PRINT
# ----------------------------------------------------------------------------
# Prints fig_race_2x2.
print(fig_race_2x2)


######## Section 11b: Equity Finding — Robustness Across Classification Schemes #####
# ----------------------------------------------------------------------------
# Checks whether the core equity finding -- the shoppability gradient is
# present in low-Medicaid/low-poverty markets and attenuates in high-
# Medicaid/high-poverty markets -- holds across all n_shop_schemes
# classification schemes, not just Theory V2. All of these per-scheme,
# per-split meta-regressions were already computed inside
# run_all_results_on_split() (called once per split in split_defs above) and
# are sitting in full_meta_summary; this section only re-filters and plots
# them, so it adds no new IV estimation.
#
# Local short-label lookup, self-contained on purpose: scheme_short (defined
# later, in the standalone figure-rebuild block) doesn't exist yet at this
# point in the script, so this duplicates the same abbreviations for
# consistency rather than depending on a not-yet-defined object.
scheme_short_local <- c(
  "CMS Rule Definition"                                    = "CMS Rule",
  "Theory-Based"                                           = "Theory-Based",
  "Broad Shoppable"                                        = "Broad",
  "Imaging vs Procedural"                                  = "Imaging vs Proc.",
  "Theory-Based V2 (MRI as nonshoppable)"                  = "Theory V2 \u25C6",
  "Split Non-Shoppable: MRI vs Procedural"                 = "Split Non-Shop.",
  "CMS Statutory Shoppable List"                           = "CMS Statutory",
  "High vs Low Shoppability Within Modality"               = "High vs Low",
  "CT-Inclusive (all CT modalities shoppable)"             = "CT-Inclusive (All CT)",
  "CT-Inclusive Except Angio (acute vascular CT excluded)" = "CT-Inclusive Ex. Angio",
  "Alt: Core Operations Framework (clinical invasiveness)"                    = "Alt: Core Ops.",
  "Alt: CMS Legal/Regulatory Framework (70 vs 230 vs non-listed)"             = "Alt: CMS Legal",
  "Alt: Upfront Cash-Market Framework (MDsave-style voucher availability)"    = "Alt: Cash-Market",
  "Alt: Geographic/Facility Access Framework (retail vs freestanding vs hospital-locked)" = "Alt: Geographic",
  "Alt: Diagnostic Urgency/Lead-Time Framework"                               = "Alt: Urgency",
  "Alt: Staffing/Specialist Framework (generalist vs sub-specialist)"        = "Alt: Staffing",
  "Alt: Incident Reporting/Liability Framework (risk-based)"                 = "Alt: Liability",
  "Alt: No Surprises Act Legal Framework (balance-billing exposure)"         = "Alt: No Surprises"
)

equity_scheme_data <- full_meta_summary %>%
  filter(dimension %in% c("Medicaid Share", "Poverty Rate")) %>%
  mutate(
    disadvantage_level = case_when(
      grepl("^Low",  split) ~ "Low",
      grepl("^High", split) ~ "High",
      TRUE ~ NA_character_
    ),
    scheme_short_lbl = ifelse(is.na(scheme_short_local[scheme]), scheme, scheme_short_local[scheme]),
    ci_lo      = estimate - 1.96 * se,
    ci_hi      = estimate + 1.96 * se,
    is_primary = scheme == "Theory-Based V2 (MRI as nonshoppable)"
  ) %>%
  filter(!is.na(disadvantage_level))

cat("\n=== EQUITY FINDING ROBUSTNESS: Shoppable Coefficient by Medicaid/Poverty Level, All Schemes ===\n")
equity_scheme_data %>%
  arrange(dimension, disadvantage_level, scheme) %>%
  select(dimension, disadvantage_level, scheme_short_lbl, estimate, se, p_value, sig) %>%
  as.data.frame() %>% print()

# ----------------------------------------------------------------------------
# EXPORT_CSV_equity_scheme_robustness
# ----------------------------------------------------------------------------
# Exports the equity-finding-by-scheme comparison to
# equity_scheme_robustness.csv.
write.csv(equity_scheme_data %>% mutate(across(where(is.numeric), ~round(.x, 4))),
          "equity_scheme_robustness.csv", row.names = FALSE)

# ----------------------------------------------------------------------------
# OUTPUT_FIGURE_EQUITY_SCHEME_ROBUSTNESS
# ----------------------------------------------------------------------------
# Figure object fig_equity_scheme_robust: checks whether the shoppability
# gradient's attenuation in high-Medicaid/high-poverty markets (the equity
# finding) holds across all n_shop_schemes classification schemes. If the
# finding is robust, within each scheme the Low-share point should sit
# further from zero (more negative) than the High-share point.
fig_equity_scheme_robust <- ggplot(
  equity_scheme_data,
  aes(x = estimate, y = reorder(scheme_short_lbl, estimate),
      color = disadvantage_level, shape = is_primary)
) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.2, linewidth = 0.6,
                 position = position_dodge(width = 0.6)) +
  geom_point(size = 2.8, position = position_dodge(width = 0.6)) +
  facet_wrap(~dimension, ncol = 1, scales = "free_y") +
  scale_color_manual(values = c("Low" = FSU_GARNET, "High" = FSU_GOLD), name = NULL,
                     labels = c("Low"  = "Low share (less disadvantaged)",
                                "High" = "High share (more disadvantaged)")) +
  scale_shape_manual(values = c("TRUE" = 18, "FALSE" = 16), guide = "none") +
  scale_x_continuous(labels = function(x) paste0(x, "%")) +
  labs(x = "Shoppable coefficient (pp)", y = NULL,
       title = "Equity Finding Robustness Across Classification Schemes",
       subtitle = paste0("Shoppability gradient by Medicaid/poverty level, all ", n_shop_schemes,
                         " schemes. Diamond = Theory V2 (primary)."),
       caption = "Whiskers = 95% CI. Low-share points sitting further from zero than High-share points, within a scheme, supports the equity finding.") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold"),
        panel.grid.major.y = element_blank())
print(fig_equity_scheme_robust)
# ----------------------------------------------------------------------------
# OUTPUT_FIGURE_EQUITY_SCHEME_ROBUSTNESS_SAVE
# ----------------------------------------------------------------------------
# Saved PDF: fig_equity_scheme_robustness.pdf.
ggsave("fig_equity_scheme_robustness.pdf", width = 8, height = 10, device = cairo_pdf)


######## Section 12: Hospital Type Heterogeneity ###########
# ----------------------------------------------------------------------------
# Hospital-type heterogeneity: compares short-term acute care and critical 
# access hospitals, then runs service-level and meta-regression versions.
df_iv_county <- df_iv_county %>%
  mutate(hospital_type_clean = case_when(
    grepl("Short Term Acute Care", HOSPITAL_TYPE, ignore.case=TRUE) ~ "Short Term Acute Care",
    grepl("Critical Access",       HOSPITAL_TYPE, ignore.case=TRUE) ~ "Critical Access",
    grepl("Psychiatric",           HOSPITAL_TYPE, ignore.case=TRUE) ~ "Psychiatric",
    grepl("Rehabilitation",        HOSPITAL_TYPE, ignore.case=TRUE) ~ "Rehabilitation",
    grepl("Children",              HOSPITAL_TYPE, ignore.case=TRUE) ~ "Childrens",
    grepl("Long.Term|Long Term",   HOSPITAL_TYPE, ignore.case=TRUE) ~ "Long Term Acute Care",
    is.na(HOSPITAL_TYPE) | HOSPITAL_TYPE == ""                      ~ NA_character_,
    TRUE ~ "Other"
  ))

cat("\nHospital type distribution:\n")
print(table(df_iv_county$hospital_type_clean, useNA="always"))

hosp_type_pooled <- bind_rows(
  run_iv_sub(df_iv_county %>% filter(hospital_type_clean=="Short Term Acute Care"),
             "Short Term Acute Care", "system_peer_pressure_county_9m",
             c("County_State","post_month")),
  run_iv_sub(df_iv_county %>% filter(hospital_type_clean=="Critical Access"),
             "Critical Access", "system_peer_pressure_county_9m",
             c("County_State","post_month"))
)
cat("\n=== HOSPITAL TYPE: Pooled IV ===\n")
print(hosp_type_pooled %>% select(group,pct,se,pvalue,wald_f,n_obs) %>%
        mutate(across(where(is.numeric),~round(.x,3))))

res_shortterm_svc <- run_iv_by_service(
  df_iv_county %>% filter(hospital_type_clean=="Short Term Acute Care"),
  "ln_median_price","system_peer_pressure_county_9m",c("County_State","post_month"))
res_cah_svc <- run_iv_by_service(
  df_iv_county %>% filter(hospital_type_clean=="Critical Access"),
  "ln_median_price","system_peer_pressure_county_9m",c("County_State","post_month"))
res_shortterm_shop2 <- run_iv_by_large_group(
  df_iv_county %>% filter(hospital_type_clean=="Short Term Acute Care"),
  "SERVICE_GROUP_SHOP2","ln_median_price","system_peer_pressure_county_9m",c("County_State","post_month"))
res_cah_shop2 <- run_iv_by_large_group(
  df_iv_county %>% filter(hospital_type_clean=="Critical Access"),
  "SERVICE_GROUP_SHOP2","ln_median_price","system_peer_pressure_county_9m",c("County_State","post_month"))

cat("\n=== HOSPITAL TYPE: SHOP2 Results ===\n")
bind_rows(res_shortterm_shop2 %>% mutate(type="Short Term Acute Care"),
          res_cah_shop2       %>% mutate(type="Critical Access")) %>%
  select(type,group,estimate_pct,se_pct,pval,fs_f,significant) %>%
  mutate(across(where(is.numeric),~round(.x,3))) %>% print()

hosp_type_meta <- lapply(
  list(list(label="Short Term Acute Care",res=res_shortterm_svc),
       list(label="Critical Access",       res=res_cah_svc)),
  function(ht) {
    if (is.null(ht$res) || nrow(ht$res) < 5) return(NULL)
    meta_rows <- list()
    for (scheme_nm in names(shoppability_schemes)) {
      scheme  <- shoppability_schemes[[scheme_nm]]
      df_meta <- build_meta_data_scheme(ht$res, scheme)
      if (nrow(df_meta) < 5) next
      m <- tryCatch(lm(estimate_pct~is_shoppable,data=df_meta,weights=weight),
                    error=function(e)NULL)
      if (is.null(m)) next
      ct <- coef(summary(m))
      if (!"is_shoppableTRUE" %in% rownames(ct)) next
      meta_rows[[scheme_nm]] <- data.frame(
        hospital_type=ht$label, scheme=scheme$label,
        coef=ct["is_shoppableTRUE","Estimate"], se=ct["is_shoppableTRUE","Std. Error"],
        p_value=ct["is_shoppableTRUE","Pr(>|t|)"], r_squared=summary(m)$r.squared)
    }
    do.call(rbind, meta_rows)
  })
hosp_type_meta_df <- do.call(rbind, hosp_type_meta) %>%
  mutate(sig=add_stars(p_value), across(where(is.numeric),~round(.x,4)))
# ----------------------------------------------------------------------------
# EXPORT_CSV_hosp_type_meta
# ----------------------------------------------------------------------------
# Exports hospital-type meta-regression results to hosp_type_meta.csv.
write.csv(hosp_type_meta_df, "hosp_type_meta.csv", row.names=FALSE)

mod_iv_shortterm <- feols(
  ln_median_price ~ 1 | market_id + post_month |
    n_prior_posters ~ system_peer_pressure_county_9m,
  data    = df_iv_county[df_iv_county$hospital_type_clean=="Short Term Acute Care",],
  cluster = ~County_State + post_month)
mod_iv_cah <- feols(
  ln_median_price ~ 1 | market_id + post_month |
    n_prior_posters ~ system_peer_pressure_county_9m,
  data    = df_iv_county[df_iv_county$hospital_type_clean=="Critical Access",],
  cluster = ~County_State + post_month)
# ----------------------------------------------------------------------------
# OUTPUT_TABLE_HOSPITAL_TYPE_IV
# ----------------------------------------------------------------------------
# Supplemental etable comparing hospital-type IV estimates.
cat("\n=== TABLE: Hospital Type IV (etable) ===\n")
etable(mod_iv_shortterm, mod_iv_cah, fitstat=c("ivwald","n"),
       headers=c("Short Term Acute Care","Critical Access"))


######## Section 13: Geographic Proximity Heterogeneity ###########
# ----------------------------------------------------------------------------
# Geographic-proximity heterogeneity: computes nearest-hospital distances and 
# compares close vs far hospitals within CBSAs.
hosp_coords <- df_iv_county[!is.na(df_iv_county$LATITUDE) & !is.na(df_iv_county$LONGITUDE),] %>%
  distinct(HOSPITAL_ID, cbsacode, LATITUDE, LONGITUDE)

# ----------------------------------------------------------------------------
# FUNCTION compute_min_distance()
# ----------------------------------------------------------------------------
# Computes nearest-neighbor hospital distance within a CBSA using 
# latitude/longitude.
compute_min_distance <- function(df_grp) {
  if (nrow(df_grp) < 2) return(data.frame(HOSPITAL_ID=df_grp$HOSPITAL_ID, dist_nearest_km=NA_real_))
  coords_mat <- as.matrix(df_grp[,c("LONGITUDE","LATITUDE")])
  dist_mat   <- distm(coords_mat, fun=distHaversine) / 1000
  diag(dist_mat) <- NA
  data.frame(HOSPITAL_ID=df_grp$HOSPITAL_ID, dist_nearest_km=apply(dist_mat,1,min,na.rm=TRUE))
}

hosp_dist <- hosp_coords %>%
  group_by(cbsacode) %>%
  group_modify(~compute_min_distance(.x)) %>%
  ungroup()

df_iv_county <- df_iv_county %>%
  left_join(hosp_dist %>% select(HOSPITAL_ID, dist_nearest_km), by="HOSPITAL_ID")

dist_cut    <- quantile(df_iv_county$dist_nearest_km, 0.50, na.rm=TRUE)
label_close <- paste0("Close (\u2264 ",round(dist_cut,1)," km)")
label_far   <- paste0("Far (> ",       round(dist_cut,1)," km)")

df_iv_county <- df_iv_county %>%
  mutate(dist_group = case_when(dist_nearest_km <= dist_cut ~ "close",
                                dist_nearest_km >  dist_cut ~ "far", TRUE ~ NA_character_))

prox_pooled <- bind_rows(
  run_iv_sub(df_iv_county %>% filter(dist_group=="close"), label_close,
             "system_peer_pressure_county_9m", c("County_State","post_month")),
  run_iv_sub(df_iv_county %>% filter(dist_group=="far"),   label_far,
             "system_peer_pressure_county_9m", c("County_State","post_month"))
)
cat("\n=== GEOGRAPHIC PROXIMITY: Pooled IV ===\n")
print(prox_pooled %>% select(group,pct,se,pvalue,wald_f,n_obs) %>%
        mutate(across(where(is.numeric),~round(.x,3))))

res_close_svc <- run_iv_by_service(df_iv_county %>% filter(dist_group=="close"),
                                   "ln_median_price","system_peer_pressure_county_9m",
                                   c("County_State","post_month"))
res_far_svc   <- run_iv_by_service(df_iv_county %>% filter(dist_group=="far"),
                                   "ln_median_price","system_peer_pressure_county_9m",
                                   c("County_State","post_month"))
res_close_shop2 <- run_iv_by_large_group(df_iv_county %>% filter(dist_group=="close"),
                                         "SERVICE_GROUP_SHOP2","ln_median_price",
                                         "system_peer_pressure_county_9m",c("County_State","post_month"))
res_far_shop2   <- run_iv_by_large_group(df_iv_county %>% filter(dist_group=="far"),
                                         "SERVICE_GROUP_SHOP2","ln_median_price",
                                         "system_peer_pressure_county_9m",c("County_State","post_month"))

prox_meta <- lapply(
  list(list(label=label_close,res=res_close_svc), list(label=label_far,res=res_far_svc)),
  function(px) {
    if (is.null(px$res) || nrow(px$res) < 5) return(NULL)
    meta_rows <- list()
    for (scheme_nm in names(shoppability_schemes)) {
      scheme  <- shoppability_schemes[[scheme_nm]]
      df_meta <- build_meta_data_scheme(px$res, scheme)
      if (nrow(df_meta) < 5) next
      m <- tryCatch(lm(estimate_pct~is_shoppable,data=df_meta,weights=weight),
                    error=function(e)NULL)
      if (is.null(m)) next
      ct <- coef(summary(m)); if (!"is_shoppableTRUE" %in% rownames(ct)) next
      meta_rows[[scheme_nm]] <- data.frame(
        proximity=px$label, scheme=scheme$label,
        coef=ct["is_shoppableTRUE","Estimate"], se=ct["is_shoppableTRUE","Std. Error"],
        p_value=ct["is_shoppableTRUE","Pr(>|t|)"], r_squared=summary(m)$r.squared)
    }
    do.call(rbind, meta_rows)
  })
prox_meta_df <- do.call(rbind, prox_meta) %>%
  mutate(sig=add_stars(p_value), across(where(is.numeric),~round(.x,4)))
# ----------------------------------------------------------------------------
# EXPORT_CSV_proximity_meta
# ----------------------------------------------------------------------------
# Exports geographic-proximity meta-regression results to proximity_meta.csv.
write.csv(prox_meta_df, "proximity_meta.csv", row.names=FALSE)

mod_iv_close <- feols(ln_median_price ~ 1 | market_id + post_month |
                        n_prior_posters ~ system_peer_pressure_county_9m,
                      data=df_iv_county %>% filter(dist_group=="close"),
                      cluster=~County_State + post_month)
mod_iv_far   <- feols(ln_median_price ~ 1 | market_id + post_month |
                        n_prior_posters ~ system_peer_pressure_county_9m,
                      data=df_iv_county %>% filter(dist_group=="far"),
                      cluster=~County_State + post_month)
# ----------------------------------------------------------------------------
# OUTPUT_TABLE_PROXIMITY_IV
# ----------------------------------------------------------------------------
# Supplemental etable comparing close vs far hospital proximity estimates.
cat("\n=== TABLE: Geographic Proximity IV (etable) ===\n")
etable(mod_iv_close, mod_iv_far, fitstat=c("ivwald","n"), headers=c(label_close, label_far))


######## Section 14: Service-Level Heterogeneity (Part B & C) ###########
# ----------------------------------------------------------------------------
# Service-level heterogeneity extensions: meta-regression interactions and 
# city-vs-county biopsy amplification check.
cat("\n=== PART B: Meta-Regression Interaction ===\n")
service_mkt_chars <- df_het %>%
  group_by(SERVICE_GROUP) %>%
  summarise(mean_hhi=mean(hhi_1000,na.rm=TRUE), mean_forprofit=mean(Share_ForProfit,na.rm=TRUE),
            mean_pop_log=mean(log(population+1),na.rm=TRUE), mean_college=mean(college_share,na.rm=TRUE),
            .groups="drop")

meta_data_interact <- build_meta_data(res_county) %>%
  left_join(service_mkt_chars, by="SERVICE_GROUP")

mod_interact_hhi  <- lm(estimate_pct~is_shoppable*mean_hhi,       data=meta_data_interact,weights=weight)
mod_interact_fp   <- lm(estimate_pct~is_shoppable*mean_forprofit,  data=meta_data_interact,weights=weight)
mod_interact_pop  <- lm(estimate_pct~is_shoppable*mean_pop_log,    data=meta_data_interact,weights=weight)
mod_interact_full <- lm(estimate_pct~is_shoppable*(mean_hhi+mean_forprofit+mean_pop_log),
                        data=meta_data_interact,weights=weight)
cat("\nB1: Shoppability x HHI:\n");        print(summary(mod_interact_hhi))
cat("\nB2: Shoppability x For-Profit:\n"); print(summary(mod_interact_fp))
cat("\nB3: Shoppability x Market Size:\n");print(summary(mod_interact_pop))
cat("\nB4: Full interaction model:\n");    print(summary(mod_interact_full))

cat("\n=== PART C: City Amplification — Biopsy Services ===\n")
biopsy_compare <- bind_rows(
  res_county %>% mutate(SERVICE_GROUP=as.character(SERVICE_GROUP)) %>%
    filter(grepl("Biopsy",SERVICE_GROUP)) %>% mutate(level="County (9m)"),
  res_city %>% mutate(SERVICE_GROUP=as.character(SERVICE_GROUP)) %>%
    filter(grepl("Biopsy",SERVICE_GROUP)) %>% mutate(level="City (12m)")
) %>% mutate(level=factor(level,levels=c("County (9m)","City (12m)")))

cat("\nBiopsy IV estimates — County vs City:\n")
biopsy_compare %>%
  select(level,SERVICE_GROUP,estimate_pct,se_pct,pval,fs_f) %>%
  mutate(sig=add_stars(pval),fmt=sprintf("%.3f%s (%.3f)",estimate_pct,sig,se_pct)) %>%
  select(level,SERVICE_GROUP,fmt,fs_f) %>% arrange(level,SERVICE_GROUP) %>%
  as.data.frame() %>% print()

biopsy_wide <- biopsy_compare %>%
  select(SERVICE_GROUP,level,estimate_pct) %>%
  mutate(level = case_when(level == "County (9m)" ~ "County", level == "City (12m)" ~ "City")) %>%
  pivot_wider(names_from=level,values_from=estimate_pct) %>%
  mutate(city_minus_county=City-County)
cat(sprintf("\nMean city-county difference: %.3f%%\n",
            mean(biopsy_wide$city_minus_county,na.rm=TRUE)))


######## Section 15: Robustness — Instrument Variants ###########
# ----------------------------------------------------------------------------
# County: 9m is primary; 6m and 12m are window robustness
# Additional: lag (6-12m), large-system exclusion
# City: 12m is primary; 9m and 6m are window robustness
# Instrument-variant robustness: alternative county windows, lagged peer 
# pressure, excluding large systems, and city variants.
# --- County: build lag and large-system-exclusion variants ---
drop_systems_robustness <- c(
  "AdventHealth","HCA Healthcare","Baylor Scott & White Health",
  "HCA Medical City Healthcare",
  "HCA Gulf Coast Division - HCA Houston Healthcare",
  "HCA West Florida Division"
)

# Lag variant: peers posting 6-12 months ago (with focal county fix)
events_lag_county <- hosp_system_county %>%
  rename(peer_hosp=HOSPITAL_ID, peer_county=County_State, peer_post_month=post_month) %>%
  cross_join(panel_county_months) %>%
  filter(peer_county != County_State,
         peer_post_month <= post_month - months(6),
         peer_post_month >= post_month - months(12)) %>%
  inner_join(focal_county_systems, by=c("HEALTH_SYSTEM_NAME","County_State")) %>%
  group_by(HEALTH_SYSTEM_NAME, County_State, post_month) %>%
  summarise(n=n_distinct(peer_hosp), .groups="drop") %>%
  group_by(County_State, post_month) %>%
  summarise(system_peer_pressure_lag=sum(n), .groups="drop")

# Large-system exclusion variant (with focal county fix)
county_month_robust <- hosp_system_county %>%
  filter(!HEALTH_SYSTEM_NAME %in% drop_systems_robustness) %>%
  rename(peer_hosp=HOSPITAL_ID, peer_county=County_State, peer_post_month=post_month) %>%
  cross_join(panel_county_months) %>%
  filter(peer_county != County_State,
         peer_post_month <= post_month,
         peer_post_month >= post_month - months(9)) %>%
  inner_join(focal_county_systems, by=c("HEALTH_SYSTEM_NAME","County_State")) %>%
  group_by(HEALTH_SYSTEM_NAME, County_State, post_month) %>%
  summarise(n=n_distinct(peer_hosp), .groups="drop") %>%
  group_by(County_State, post_month) %>%
  summarise(system_peer_pressure_robust=sum(n), .groups="drop")

df_iv_county <- df_iv_county %>%
  select(-any_of(c("system_peer_pressure_lag","system_peer_pressure_robust"))) %>%
  left_join(events_lag_county,   by=c("County_State","post_month")) %>%
  left_join(county_month_robust, by=c("County_State","post_month")) %>%
  mutate(system_peer_pressure_lag    = replace_na(system_peer_pressure_lag,    0),
         system_peer_pressure_robust = replace_na(system_peer_pressure_robust, 0))

# County F-stat comparison across all variants
variants_county <- list(
  "9m (primary)"         = "system_peer_pressure_county_9m",
  "6m (robustness)"      = "system_peer_pressure_county",
  "12m (robustness)"     = "system_peer_pressure_county_12m",
  "Lag 6-12m"            = "system_peer_pressure_lag",
  "No large systems (9m)"= "system_peer_pressure_robust"
)

# ----------------------------------------------------------------------------
# OUTPUT_TABLE_20_IV_VARIANTS
# ----------------------------------------------------------------------------
# Console table of county instrument variants and first-stage diagnostics.
cat("\n=== TABLE 20: County Instrument Variants (two-way clustering) ===\n")
cat(sprintf("%-30s  %8s  %8s  %8s  %8s\n","Variant","Coef","SE","t-stat","Wald F"))
cat(strrep("-",65),"\n")
results_variants_county <- map_dfr(names(variants_county), function(nm) {
  instr <- variants_county[[nm]]
  if (!instr %in% names(df_iv_county)) return(NULL)
  fit <- tryCatch(feols(as.formula(paste0(
    "ln_median_price ~ ln_total_beds | market_id + post_month | n_prior_posters ~ ",instr)),
    data=df_iv_county, cluster=~County_State+post_month), error=function(e)NULL)
  if (is.null(fit)) return(NULL)
  wald <- tryCatch(fitstat(fit,"ivwald")[[1]]$stat, error=function(e)NA_real_)
  est  <- coef(fit)["fit_n_prior_posters"]
  se_v <- se(fit)["fit_n_prior_posters"]
  cat(sprintf("%-30s  %8.4f  %8.4f  %8.2f  %8.1f\n",
              nm, est*100, se_v*100, est/se_v, wald))
  tibble(variant=nm, estimate=est*100, se=se_v*100,
         pval=pvalue(fit)["fit_n_prior_posters"], wald_f=wald, n_obs=nobs(fit))
})

# City F-stat comparison across window variants
cat("\n=== City Instrument Variants (two-way clustering) ===\n")
cat(sprintf("%-30s  %8s  %8s  %8s\n","Variant","t-stat","Wald F","N"))
cat(strrep("-",55),"\n")
variants_city <- list(
  "12m (primary)"   = "system_peer_pressure_city_12m",
  "9m (robustness)" = "system_peer_pressure_city_9m",
  "6m (robustness)" = "system_peer_pressure_city"
)
for (nm in names(variants_city)) {
  instr <- variants_city[[nm]]
  if (!instr %in% names(df_iv_city)) { cat(sprintf("%-30s  MISSING\n",nm)); next }
  fit <- tryCatch(feols(as.formula(paste0(
    "ln_median_price ~ ln_total_beds | market_id + post_month | n_prior_posters ~ ",instr)),
    data=df_iv_city, cluster=~city_state+post_month), error=function(e)NULL)
  if (is.null(fit)) { cat(sprintf("%-30s  FAILED\n",nm)); next }
  wald <- tryCatch(fitstat(fit,"ivwald")[[1]]$stat, error=function(e)NA_real_)
  t_v  <- se(fit)["fit_n_prior_posters"]
  cat(sprintf("%-30s  %8.2f  %8.1f  %8d\n",nm,coef(fit)["fit_n_prior_posters"]/t_v,wald,nobs(fit)))
}


######## Section 16: Robustness — Leave-One-System-Out ###########
# ----------------------------------------------------------------------------
# Primary instrument: 9m; focal county fix applied to each LOO
# Leave-one-system-out robustness: iteratively drops top systems and re-
# estimates the shoppability gradient.
top_systems_county <- df_iv_county %>%
  distinct(HOSPITAL_ID, HEALTH_SYSTEM_NAME) %>%
  filter(!is.na(HEALTH_SYSTEM_NAME)) %>%
  count(HEALTH_SYSTEM_NAME, sort=TRUE) %>%
  head(15) %>%
  pull(HEALTH_SYSTEM_NAME)

loo_results_county <- map_dfr(top_systems_county, function(sys) {
  events_loo <- hosp_system_county %>%
    filter(HEALTH_SYSTEM_NAME != sys) %>%
    rename(peer_hosp=HOSPITAL_ID, peer_county=County_State, peer_post_month=post_month) %>%
    cross_join(panel_county_months) %>%
    filter(peer_county != County_State,
           peer_post_month <= post_month,
           peer_post_month >= post_month - months(9)) %>%
    inner_join(focal_county_systems, by=c("HEALTH_SYSTEM_NAME","County_State")) %>%
    group_by(HEALTH_SYSTEM_NAME, County_State, post_month) %>%
    summarise(n=n_distinct(peer_hosp), .groups="drop") %>%
    group_by(County_State, post_month) %>%
    summarise(system_peer_pressure_loo=sum(n), .groups="drop")
  
  df_loo <- df_iv_county %>%
    select(-any_of("system_peer_pressure_loo")) %>%
    left_join(events_loo, by=c("County_State","post_month")) %>%
    mutate(system_peer_pressure_loo=replace_na(system_peer_pressure_loo,0))
  
  res_loo <- do.call(rbind, lapply(unique(df_loo$SERVICE_GROUP), function(sg) {
    df_sub <- df_loo[df_loo$SERVICE_GROUP==sg,]
    if (nrow(df_sub) < 100) return(NULL)
    fit <- tryCatch(feols(ln_median_price~1|market_id+post_month|
                            n_prior_posters~system_peer_pressure_loo,
                          data=df_sub, cluster=~County_State+post_month),
                    error=function(e)NULL)
    if (is.null(fit)) return(NULL)
    data.frame(SERVICE_GROUP=sg,
               estimate_pct=tryCatch(coef(fit)["fit_n_prior_posters"]*100,error=function(e)NA_real_),
               se_pct=tryCatch(se(fit)["fit_n_prior_posters"]*100,error=function(e)NA_real_))
  }))
  if (is.null(res_loo) || nrow(res_loo)==0) return(NULL)
  
  meta_loo <- res_loo %>%
    mutate(is_shoppable=grepl("Ultrasound",SERVICE_GROUP)|SERVICE_GROUP%in%c("CT Lung","Mammography"),
           weight=1/se_pct^2) %>% filter(is.finite(weight))
  
  coef_shop <- tryCatch({
    m <- lm(estimate_pct~is_shoppable, data=meta_loo, weights=weight)
    coef(m)["is_shoppableTRUE"]
  }, error=function(e)NA_real_)
  
  tibble(dropped_system=sys, shop_coef=coef_shop)
})

# ----------------------------------------------------------------------------
# OUTPUT_TABLE_21_LEAVE_ONE_SYSTEM_OUT
# ----------------------------------------------------------------------------
# Console table of leave-one-system-out shoppability gradients.
cat("\n=== TABLE 21: Leave-One-System-Out (County 9m) ===\n")
loo_results_county %>%
  arrange(shop_coef) %>%
  mutate(shop_coef=round(shop_coef,3)) %>%
  as.data.frame() %>% print()
loo_results_county %>%
  summarise(min_coef=min(shop_coef,na.rm=TRUE), max_coef=max(shop_coef,na.rm=TRUE),
            mean_coef=mean(shop_coef,na.rm=TRUE), sd_coef=sd(shop_coef,na.rm=TRUE)) %>%
  as.data.frame() %>% print()


######## Section 17: Placebo Testing ###########
# ----------------------------------------------------------------------------
# Placebo testing: instrument leads/pre-trends and randomization inference for
#  first-stage strength.
# Pre-trend test: instrument leads
instrument_leads <- df_iv_county %>%
  select(County_State, post_month, system_peer_pressure_county_9m) %>%
  distinct() %>%
  arrange(County_State, post_month) %>%
  group_by(County_State) %>%
  mutate(spp_lead1 = lead(system_peer_pressure_county_9m, 1),
         spp_lead2 = lead(system_peer_pressure_county_9m, 2),
         spp_lead3 = lead(system_peer_pressure_county_9m, 3)) %>%
  ungroup() %>%
  select(County_State, post_month, spp_lead1, spp_lead2, spp_lead3)

df_leads <- df_iv_county %>%
  left_join(instrument_leads, by=c("County_State","post_month"))

pt_lead1 <- feols(ln_median_price~spp_lead1+ln_total_beds|market_id+post_month,
                  data=df_leads, cluster=~County_State+post_month)
pt_lead2 <- feols(ln_median_price~spp_lead2+ln_total_beds|market_id+post_month,
                  data=df_leads, cluster=~County_State+post_month)
pt_lead3 <- feols(ln_median_price~spp_lead3+ln_total_beds|market_id+post_month,
                  data=df_leads, cluster=~County_State+post_month)
pt_joint  <- feols(ln_median_price~spp_lead1+spp_lead2+spp_lead3+ln_total_beds|market_id+post_month,
                   data=df_leads, cluster=~County_State+post_month)

# ----------------------------------------------------------------------------
# OUTPUT_TABLE_23_PLACEBO_LEADS
# ----------------------------------------------------------------------------
# Pre-trend/placebo lead estimates for the instrument.
cat("\n=== TABLE 23: Pre-Trend Test (Instrument Leads, County 9m) ===\n")
etable(pt_lead1, pt_lead2, pt_lead3, pt_joint,
       headers=c("Lead +1","Lead +2","Lead +3","Joint"),
       keep=c("spp_lead1","spp_lead2","spp_lead3"), fitstat=~f+n)
wald_joint_leads <- wald(pt_joint, c("spp_lead1","spp_lead2","spp_lead3"))
cat(sprintf("\nJoint F-test on all three leads: F = %.3f, p = %.4f\n",
            wald_joint_leads$stat, wald_joint_leads$p))


library(broom)


# ------------------------------------------------------------
# Pre-trend coefficient figure: instrument leads
# Uses the joint specification pt_joint
# ------------------------------------------------------------

# Extract coefficient estimates from the joint model
pretrend_plot_df <- broom::tidy(pt_joint, conf.int = TRUE, conf.level = 0.95) %>%
  filter(term %in% c("spp_lead1", "spp_lead2", "spp_lead3")) %>%
  mutate(
    lead = case_when(
      term == "spp_lead1" ~ 1,
      term == "spp_lead2" ~ 2,
      term == "spp_lead3" ~ 3
    ),
    lead_label = factor(
      lead,
      levels = c(1, 2, 3),
      labels = c("Lead +1", "Lead +2", "Lead +3")
    ),
    sig_5pct = conf.low > 0 | conf.high < 0
  )

# Pull joint Wald test for subtitle/caption
wald_joint_leads <- wald(pt_joint, c("spp_lead1", "spp_lead2", "spp_lead3"))

joint_f <- round(wald_joint_leads$stat, 2)
joint_p <- signif(wald_joint_leads$p, 3)

# Optional: scale to percentage points if outcome is log price.
# Since ln_median_price is in logs, coefficients are approximate log points.
# Multiplying by 100 makes the plot interpretable as approximate percent effects.
pretrend_plot_df <- pretrend_plot_df %>%
  mutate(
    estimate_pp = 100 * estimate,
    ci_low_pp   = 100 * conf.low,
    ci_high_pp  = 100 * conf.high
  )

# ----------------------------------------------------------------------------
# OUTPUT_FIGURE_PRETREND_LEADS
# ----------------------------------------------------------------------------
# Figure object fig_pretrend_leads: coefficients on instrument leads.
fig_pretrend_leads <- ggplot(pretrend_plot_df, aes(x = lead, y = estimate_pp)) +
  
  # Zero reference line
  geom_hline(
    yintercept = 0,
    linewidth = 0.45,
    linetype = "dashed",
    color = "grey55"
  ) +
  
  # Confidence intervals
  geom_errorbar(
    aes(ymin = ci_low_pp, ymax = ci_high_pp),
    width = 0.08,
    linewidth = 0.7,
    color = "#2C3E50"
  ) +
  
  # Coefficient points
  geom_point(
    aes(shape = sig_5pct),
    size = 3.4,
    stroke = 1.0,
    color = "#8B0000"
  ) +
  
  # Connect estimates lightly to emphasize pattern
  geom_line(
    linewidth = 0.55,
    color = "#8B0000",
    alpha = 0.65
  ) +
  
  # Label the one coefficient that barely clears zero, if applicable
  geom_text(
    data = pretrend_plot_df %>% filter(sig_5pct),
    aes(
      label = "*",
      y = ci_high_pp + 0.02 * diff(range(c(ci_low_pp, ci_high_pp), na.rm = TRUE))
    ),
    size = 5,
    fontface = "bold",
    color = "#8B0000"
  ) +
  
  scale_x_continuous(
    breaks = c(1, 2, 3),
    labels = c("Lead +1", "Lead +2", "Lead +3")
  ) +
  
  scale_shape_manual(
    values = c("FALSE" = 16, "TRUE" = 17),
    guide = "none"
  ) +
  
  labs(
    x = NULL,
    y = "Coefficient on Future Peer Pressure\nApprox. percent effect on median price",
    title = "Pre-Trend Test: Future Instrument Values Predict Current Prices",
    subtitle = paste0(
      "Joint test rejects statistically: F = ", joint_f,
      ", p = ", joint_p,
      "; plotted effects are small and wrong-signed"
    ),
    caption = paste0(
      "Notes: Points show coefficients from a single joint regression of log median price on three future values ",
      "of county 9-month system peer pressure, controlling for log beds with market and month fixed effects. ",
      "Bars show 95% confidence intervals clustered by county-state and month."
    )
  ) +
  
  theme_minimal(base_size = 11) +
  
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9.5, color = "grey30"),
    plot.caption = element_text(
      size = 7.5,
      color = "grey40",
      hjust = 0,
      lineheight = 1.25
    ),
    axis.text.x = element_text(size = 10),
    axis.title.y = element_text(size = 9.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

# ----------------------------------------------------------------------------
# OUTPUT_FIGURE_PRETREND_LEADS_PRINT
# ----------------------------------------------------------------------------
# Prints fig_pretrend_leads.
print(fig_pretrend_leads)

# Randomization inference
set.seed(42)
N_PERMS <- 500

instr_spine <- df_iv_county %>%
  select(County_State, post_month, system_peer_pressure_county_9m,
         n_prior_posters, market_id, ln_total_beds) %>%
  distinct()

fs_baseline <- feols(
  n_prior_posters ~ system_peer_pressure_county_9m + ln_total_beds |
    market_id + post_month,
  data=df_iv_county, cluster=~County_State+post_month)
actual_f <- fitstat(fs_baseline, "f")[[1]]$stat

perm_f <- numeric(N_PERMS)
cat(sprintf("\nRunning %d permutations (county 9m, two-way clustering)...\n", N_PERMS))

for (i in seq_len(N_PERMS)) {
  if (i %% 100 == 0) cat(sprintf("  %d / %d\n", i, N_PERMS))
  instr_perm <- instr_spine %>%
    select(County_State, post_month, system_peer_pressure_county_9m) %>%
    distinct() %>%
    group_by(post_month) %>%
    mutate(spp_perm = sample(system_peer_pressure_county_9m, size=n(), replace=FALSE)) %>%
    ungroup() %>%
    select(County_State, post_month, spp_perm)
  df_perm <- instr_spine %>% left_join(instr_perm, by=c("County_State","post_month"))
  fit_perm <- tryCatch(
    feols(n_prior_posters~spp_perm+ln_total_beds|market_id+post_month,
          data=df_perm, cluster=~County_State+post_month),
    error=function(e)NULL)
  perm_f[i] <- tryCatch(fitstat(fit_perm,"f")[[1]]$stat, error=function(e)NA_real_)
}

perm_f_clean <- perm_f[!is.na(perm_f)]
p95_perm  <- quantile(perm_f_clean, 0.95)
ri_pvalue <- mean(perm_f_clean >= actual_f)

cat(sprintf("\n=== Randomization Inference (County 9m) ===\n"))
cat(sprintf("Actual first-stage F       : %7.1f\n", actual_f))
cat(sprintf("Permutation 95th pctile    : %7.1f\n", p95_perm))
cat(sprintf("RI p-value                 : %7.4f\n", ri_pvalue))
cat(sprintf("Permutations completed     : %7d / %d\n", length(perm_f_clean), N_PERMS))

x_cap <- quantile(perm_f_clean, 0.99, na.rm=TRUE)
# ----------------------------------------------------------------------------
# OUTPUT_FIGURE_08_RANDOMIZATION
# ----------------------------------------------------------------------------
# Figure object fig_rand_inference: randomization inference distribution.
fig_rand_inference <- ggplot(tibble(f_stat=perm_f_clean) %>% filter(f_stat<=x_cap), aes(x=f_stat)) +
  geom_histogram(bins=40, fill="grey80", color="white", linewidth=0.25) +
  geom_rug(alpha=0.15, linewidth=0.25) +
  geom_vline(xintercept=p95_perm, color=FSU_GOLD, linewidth=0.9, linetype="dashed") +
  {if(actual_f<=x_cap) geom_vline(xintercept=actual_f,color=FSU_GARNET,linewidth=1.1)
    else geom_vline(xintercept=x_cap,color=FSU_GARNET,linewidth=1.1)} +
  annotate("label", x=Inf, y=Inf, hjust=1.05, vjust=1.2,
           label=sprintf("Actual F = %.1f\nPermutation 95th pctile = %.1f\nRI p-value = %.4f",
                         actual_f, p95_perm, ri_pvalue),
           size=3.2, label.size=0.25, fill="white") +
  scale_x_continuous(labels=scales::comma) +
  labs(x="First-stage F-statistic from permuted instruments",
       y="Number of permutations",
       title="Randomization Inference for First-Stage Strength",
       subtitle=sprintf("%d within-month permutations; county 9m instrument; two-way clustering",
                        length(perm_f_clean)),
       caption=paste0("Gold dashed = permutation 95th percentile. Garnet = actual F",
                      ifelse(actual_f>x_cap," (shown at axis cap).","."),
                      " SE clustered: County_State + post_month.")) +
  theme_bw(base_size=11) +
  theme(panel.grid.minor=element_blank(), plot.title=element_text(face="bold"))
# ----------------------------------------------------------------------------
# OUTPUT_FIGURE_08_RANDOMIZATION_PRINT
# ----------------------------------------------------------------------------
# Prints fig_rand_inference.
print(fig_rand_inference)




######## Section 18: Exclusion-Restriction Diagnostics ###########
# ----------------------------------------------------------------------------
# Additional exclusion-restriction stress tests and robustness diagnostics,
# in four parts: (A) system x month fixed effects and system-specific linear
# trends, (B) Lee et al. (2022) tF diagnostics, (C) CMS enforcement-exposure
# controls, and (D) Conley-Hansen-Rossi sensitivity bounds on the
# meta-regression gradient, ending in a robustness summary table.
# ---- Part A helper: print a compact IV diagnostic row ----
# ----------------------------------------------------------------------------
# FUNCTION iv_row()
# ----------------------------------------------------------------------------
# Extracts one row of IV diagnostics for the exclusion-restriction robustness 
# summary.
iv_row <- function(fit, label, instrument = "system_peer_pressure_county_9m") {
  est  <- tryCatch(coef(fit)["fit_n_prior_posters"],   error = function(e) NA_real_)
  se_v <- tryCatch(se(fit)["fit_n_prior_posters"],     error = function(e) NA_real_)
  pval <- tryCatch(pvalue(fit)["fit_n_prior_posters"], error = function(e) NA_real_)
  f_v  <- tryCatch(
    fitstat(fit, "ivwald")[["ivwald1::n_prior_posters"]]$stat,
    error = function(e) NA_real_)
  fs_ct <- tryCatch(fit$iv_first_stage[["n_prior_posters"]]$coeftable, error = function(e) NULL)
  fs_coef <- tryCatch(fs_ct[instrument, "Estimate"],   error = function(e) NA_real_)
  fs_se   <- tryCatch(fs_ct[instrument, "Std. Error"], error = function(e) NA_real_)
  data.frame(
    spec       = label,
    iv_est_pct = est * 100,
    iv_se_pct  = se_v * 100,
    iv_pval    = pval,
    wald_f     = f_v,
    fs_coef    = fs_coef,
    fs_se      = fs_se,
    row.names  = NULL
  )
}

INSTR <- "system_peer_pressure_county_9m"
CL    <- ~County_State + post_month

# ----------------------------------------------------------------------------
# Part A: System × Month Fixed Effects
# 
# Purpose: Absorbs system-wide pricing shocks, corporate policy
# changes, and any unobservable correlated with system-level 
# timing (threats 2, 3, 4). Identification then comes purely
# from within-system, cross-county timing variation.
#
# Feasibility note: This requires within-system, cross-county
# variation in the instrument CONDITIONAL on system×month.
# If systems are small (few counties each), the first 
# stage may weaken substantially — check Wald F carefully.
# ============================================================

# Check how many systems have multi-county presence
cat("\n=== System Size Distribution ===\n")
df_iv_county %>%
  distinct(HEALTH_SYSTEM_NAME, COUNTY_STATE) %>%
  count(HEALTH_SYSTEM_NAME, name = "n_counties") %>%
  filter(!is.na(HEALTH_SYSTEM_NAME)) %>%
  summarise(
    systems_total      = n(),
    systems_1county    = sum(n_counties == 1),
    systems_multi      = sum(n_counties > 1),
    median_counties    = median(n_counties),
    p90_counties       = quantile(n_counties, 0.9)
  ) %>%
  print()

# Baseline IV (for comparison)
iv_baseline <- feols(
  ln_median_price ~ ln_total_beds | market_id + post_month |
    n_prior_posters ~ system_peer_pressure_county_9m,
  data    = df_iv_county,
  cluster = CL
)

# System × month FE specification
# fixest syntax: HEALTH_SYSTEM_NAME^post_month creates a 
# system-month interaction FE. Note: single-county systems
# will be absorbed into this FE and drop out of identification.
iv_sys_month <- feols(
  ln_median_price ~ ln_total_beds | market_id + HEALTH_SYSTEM_NAME^post_month |
    n_prior_posters ~ system_peer_pressure_county_9m,
  data    = df_iv_county %>% filter(!is.na(HEALTH_SYSTEM_NAME)),
  cluster = CL
)

# System-specific linear trends
# Creates system × linear-time-trend controls.
# Absorbs differential compliance trajectories by system sophistication.
df_iv_county <- df_iv_county %>%
  mutate(month_index = as.integer(factor(post_month, levels = sort(unique(post_month)))))

iv_sys_trends <- feols(
  ln_median_price ~ ln_total_beds + i(HEALTH_SYSTEM_NAME, month_index) |
    market_id + post_month |
    n_prior_posters ~ system_peer_pressure_county_9m,
  data    = df_iv_county %>% filter(!is.na(HEALTH_SYSTEM_NAME)),
  cluster = CL
)

# ----------------------------------------------------------------------------
# OUTPUT_TABLE_R1_SYSTEM_FE_TRENDS
# ----------------------------------------------------------------------------
# Exclusion-restriction robustness: baseline, system-month FE, and system 
# trend specifications.
cat("\n=== TABLE R1: System × Month FE and System Trends ===\n")
etable(
  iv_baseline, iv_sys_month, iv_sys_trends,
  keep    = "%fit_n_prior_posters",
  headers = c("Baseline", "System×Month FE", "System Linear Trends"),
  fitstat = ~ n + ivwald1
)

# Wald F comparison
cat("\n  Baseline      Wald F:", 
    round(fitstat(iv_baseline,  "ivwald")[["ivwald1::n_prior_posters"]]$stat, 1))
cat("\n  Sys×Month FE  Wald F:", 
    round(fitstat(iv_sys_month, "ivwald")[["ivwald1::n_prior_posters"]]$stat, 1))
cat("\n  Sys Trends    Wald F:", 
    round(fitstat(iv_sys_trends,"ivwald")[["ivwald1::n_prior_posters"]]$stat, 1), "\n")



# ----------------------------------------------------------------------------
# Standalone robustness variant: service-level IV with System x Month FE
# Mirrors run_iv_by_service() exactly, but replaces post_month with
# HEALTH_SYSTEM_NAME^post_month in the FE slot. Kept separate so the
# original run_iv_by_service() definition used elsewhere is untouched.
# ----------------------------------------------------------------------------
run_iv_by_service_sysmth <- function(data, outcome_var = "ln_median_price",
                                     instrument, cluster_var, min_obs = 100) {
  service_groups <- unique(data$SERVICE_GROUP)
  results <- do.call(rbind, lapply(service_groups, function(sg) {
    df_sub <- data[data$SERVICE_GROUP == sg, ]
    if (nrow(df_sub) < min_obs) return(NULL)
    fit <- tryCatch(
      feols(as.formula(paste0(outcome_var,
                              " ~ ln_total_beds | market_id + HEALTH_SYSTEM_NAME^post_month | n_prior_posters ~ ",
                              instrument)),
            data    = df_sub,
            cluster = as.formula(paste0("~", paste(cluster_var, collapse = " + ")))),
      error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    fs   <- extract_first_stage(fit, instrument)
    est  <- tryCatch(coef(fit)["fit_n_prior_posters"],   error = function(e) NA_real_)
    se_v <- tryCatch(se(fit)["fit_n_prior_posters"],     error = function(e) NA_real_)
    pval <- tryCatch(pvalue(fit)["fit_n_prior_posters"], error = function(e) NA_real_)
    data.frame(SERVICE_GROUP = sg,
               n_obs        = nrow(df_sub),
               n_used       = nobs(fit),
               estimate = est, se = se_v, pval = pval,
               estimate_pct = est * 100, se_pct = se_v * 100,
               ci_lo = est - 1.96 * se_v, ci_hi = est + 1.96 * se_v,
               ci_lo_pct = (est - 1.96 * se_v) * 100,
               ci_hi_pct = (est + 1.96 * se_v) * 100,
               significant = !is.na(est) & !is.na(se_v) &
                 sign(est - 1.96*se_v) == sign(est + 1.96*se_v),
               fs_coef = fs$fs_coef, fs_se = fs$fs_se,
               fs_t = fs$fs_t, fs_f = fs$fs_f,
               row.names = NULL)
  }))
  results <- results[order(results$estimate_pct), ]
  results$SERVICE_GROUP <- factor(results$SERVICE_GROUP, levels = results$SERVICE_GROUP)
  results
}

# Service-level IV with system x month FEs
res_county_sysmth <- run_iv_by_service_sysmth(
  df_iv_county %>% filter(!is.na(HEALTH_SYSTEM_NAME)),
  outcome_var  = "ln_median_price",
  instrument   = "system_peer_pressure_county_9m",
  cluster_var  = c("County_State", "post_month"),
  min_obs      = 100
)
print(res_county_sysmth)

# Build meta data and run meta-regression
meta_sysmth <- build_meta_data(res_county_sysmth)
meta_results_sysmth <- run_meta_regressions(meta_sysmth)
print_meta_results(meta_results_sysmth, "System x Month FE")

# Compare gradient directly (same weighting scheme in both)
baseline_grad <- coef(
  lm(estimate_pct ~ is_shoppable, data = build_meta_data(res_county),
     weights = weight)
)["is_shoppableTRUE"]

sysmth_grad <- coef(
  lm(estimate_pct ~ is_shoppable, data = meta_sysmth,
     weights = weight)
)["is_shoppableTRUE"]

cat(sprintf("\nBaseline gradient    : %.3f pp\n", baseline_grad))
cat(sprintf("Sys x Month gradient : %.3f pp\n", sysmth_grad))


# ----------------------------------------------------------------------------
# Part B: Lee et al. (2022) tF Diagnostics
# Instrument: system_peer_pressure_county_9m
# Clustering: County_State + post_month (two-way)
# Reference: Lee, McCrary, Moreira & Porter (2022), Journal of
# Econometrics. "Valid t-ratio Inference for IV"
#
# Key idea: With a weak first stage, the conventional 2SLS 
# t-statistic is oversized. The tF procedure adjusts the 
# critical value for the second-stage t-test based on the 
# observed first-stage F. As F → ∞, c(F) → 1.96.
# The correct F to use is the CLUSTERED Wald F (ivwald1),
# not the unclustered ivf reported by fixest.
# ------------------------------------------------------------
# Part B.1: Lee et al. critical value interpolation function
# ------------------------------------------------------------
# Source: Lee et al. (2022), Table 1 (5% significance level,
# one endogenous variable, one instrument)
# As F -> Inf, critical value -> 1.96 (standard normal)

# ----------------------------------------------------------------------------
# FUNCTION lee_critical_value()
# ----------------------------------------------------------------------------
# Returns approximate Lee et al. critical values for weak-instrument 
# sensitivity reporting.
lee_critical_value <- function(F_stat) {
  # Table 1 values from Lee et al. (2022)
  f_grid <- c(0,     5,    10,    15,    20,    25,
              30,    40,    50,   100,   200,   Inf)
  c_grid <- c(Inf, Inf,  3.43,  2.72,  2.40,  2.20,
              2.09,  1.98,  1.93,  1.89,  1.88,  1.96)
  # Linear interpolation; extrapolate using rule=2 at boundaries
  approx(f_grid, c_grid, xout = F_stat, rule = 2)$y
}

# Sanity check the function
cat("=== Lee et al. Critical Value Function Check ===\n")
test_f <- c(10, 15, 20, 25, 27.1, 30, 37.0, 40, 50, 100)
data.frame(
  F_stat = test_f,
  lee_cv = round(sapply(test_f, lee_critical_value), 3)
) %>% print(row.names = FALSE)

# ------------------------------------------------------------
# Part B.2: Extract clustered Wald F from primary specs
# ------------------------------------------------------------
# These are the correct F statistics for Lee et al. purposes.
# fixest's ivwald1 uses clustered SEs; ivf does not.

# ----------------------------------------------------------------------------
# FUNCTION extract_wald_f()
# ----------------------------------------------------------------------------
# Extracts the IV Wald F-statistic from a fixest IV object.
extract_wald_f <- function(fit) {
  tryCatch(
    fitstat(fit, "ivwald")[["ivwald1::n_prior_posters"]]$stat,
    error = function(e) NA_real_
  )
}

# ----------------------------------------------------------------------------
# FUNCTION extract_t_stat()
# ----------------------------------------------------------------------------
# Extracts the first-stage t-statistic from a fixest IV object.
extract_t_stat <- function(fit) {
  tryCatch(
    coef(fit)["fit_n_prior_posters"] / se(fit)["fit_n_prior_posters"],
    error = function(e) NA_real_
  )
}

f_baseline <- extract_wald_f(iv_county_pooled)
f_sysmth   <- extract_wald_f(iv_sys_month)
f_trends   <- extract_wald_f(iv_sys_trends)

t_baseline <- extract_t_stat(iv_county_pooled)
t_sysmth   <- extract_t_stat(iv_sys_month)
t_trends   <- extract_t_stat(iv_sys_trends)

cv_baseline <- lee_critical_value(f_baseline)
cv_sysmth   <- lee_critical_value(f_sysmth)
cv_trends   <- lee_critical_value(f_trends)

cat("\n=== SECTION 2: Pooled Spec — Lee et al. tF Assessment ===\n")
cat(sprintf("%-25s  %8s  %8s  %8s  %8s  %8s\n",
            "Specification", "Wald F", "Lee CV", "2nd-st t", "|t|>CV?", "SY pass?"))
cat(strrep("-", 75), "\n")

specs <- list(
  list(label = "Baseline",        f = f_baseline, cv = cv_baseline, t = t_baseline),
  list(label = "Sys x Month FE",  f = f_sysmth,   cv = cv_sysmth,   t = t_sysmth),
  list(label = "Sys Lin. Trends", f = f_trends,   cv = cv_trends,   t = t_trends)
)

for (s in specs) {
  sy_pass  <- s$f > 16.38
  lee_pass <- abs(s$t) > s$cv
  cat(sprintf("%-25s  %8.2f  %8.3f  %8.3f  %8s  %8s\n",
              s$label, s$f, s$cv, s$t,
              ifelse(lee_pass, "YES", "NO"),
              ifelse(sy_pass,  "YES", "NO")))
}

cat("\nNote: Stock-Yogo (2005) 10% maximal IV size threshold = 16.38\n")
cat("      Lee et al. (2022) CV approaches 1.96 as F → ∞\n")
cat("      Clustered Wald F (ivwald1) is the correct statistic;\n")
cat("      fixest's ivf (unclustered) is NOT appropriate here.\n")

# ------------------------------------------------------------
# Part B.3: Service-level Lee et al. assessment
# ------------------------------------------------------------
# For each service group, check whether:
# (a) First-stage F clears Stock-Yogo (16.38)
# (b) Second-stage |t| clears Lee et al. adjusted CV
# (c) Result remains significant after tF correction

cat("\n=== SECTION 3: Service-Level Lee et al. Assessment ===\n")

service_lee <- res_county %>%
  mutate(
    SERVICE_GROUP = as.character(SERVICE_GROUP),
    t_stat        = estimate / se,
    lee_cv        = sapply(fs_f, lee_critical_value),
    passes_SY     = fs_f > 16.38,
    passes_lee_tF = abs(t_stat) > lee_cv,
    is_shoppable  = SERVICE_GROUP %in% c("CT Lung", "Mammography") |
      grepl("Ultrasound", SERVICE_GROUP),
    sig_conventional = significant
  ) %>%
  select(SERVICE_GROUP, is_shoppable, n_obs,
         estimate_pct, se_pct, t_stat,
         fs_f, lee_cv,
         passes_SY, sig_conventional, passes_lee_tF) %>%
  arrange(fs_f) %>%
  as.data.frame()   # strip tibble/factor attributes causing the print error

# Print with explicit na handling
print(service_lee, na.print = "NA")

# Or use this if the above still errors:
cat(format_df <- capture.output(
  print(service_lee, na.print = "-")
), sep = "\n")

# Most robust option — avoid print entirely:
service_lee %>%
  mutate(across(where(is.numeric), ~round(., 3)),
         across(where(is.logical), ~ifelse(is.na(.), "-", as.character(.)))) %>%
  knitr::kable()

# Summary counts — these don't depend on print working
cat(sprintf("\nTotal service groups            : %d\n", nrow(service_lee)))
cat(sprintf("Pass Stock-Yogo (F > 16.38)     : %d\n", sum(service_lee$passes_SY,        na.rm = TRUE)))
cat(sprintf("Fail Stock-Yogo                 : %d\n", sum(!service_lee$passes_SY,        na.rm = TRUE)))
cat(sprintf("Significant (conventional 1.96) : %d\n", sum(service_lee$sig_conventional, na.rm = TRUE)))
cat(sprintf("Significant (Lee et al. tF)     : %d\n", sum(service_lee$passes_lee_tF,    na.rm = TRUE)))

# Significant services — the key output
cat("\n--- Services significant at conventional threshold ---\n")
service_lee %>%
  filter(sig_conventional == TRUE) %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

# ------------------------------------------------------------
# Part B.4: CT Lung deep dive
# ------------------------------------------------------------
# CT Lung is your only individually significant IV result.
# Check it carefully against Lee et al.

cat("\n=== SECTION 4: CT Lung — Lee et al. Deep Dive ===\n")

ct_row <- service_lee %>% filter(SERVICE_GROUP == "CT Lung")

cat(sprintf("  First-stage Wald F    : %.2f\n", ct_row$fs_f))
cat(sprintf("  Stock-Yogo threshold  : 16.38  → %s\n",
            ifelse(ct_row$passes_SY, "PASS", "FAIL")))
cat(sprintf("  Lee et al. CV         : %.3f\n", ct_row$lee_cv))
cat(sprintf("  Second-stage t-stat   : %.3f\n", ct_row$t_stat))
cat(sprintf("  |t| > Lee CV?         : %s\n",
            ifelse(ct_row$passes_lee_tF, "YES — significant under tF correction",
                   "NO  — loses significance under tF correction")))
cat(sprintf("  Estimate              : %.3f pp\n", ct_row$estimate_pct))

# Also run CT Lung standalone IV to get clean diagnostics
df_ct <- df_iv_county %>% filter(SERVICE_GROUP == "CT Lung")

iv_ct_lung <- feols(
  ln_median_price ~ ln_total_beds | market_id + post_month |
    n_prior_posters ~ system_peer_pressure_county_9m,
  data    = df_ct,
  cluster = ~County_State + post_month
)

cat("\n--- CT Lung standalone IV summary ---\n")
summary(iv_ct_lung)

ct_wald_f  <- extract_wald_f(iv_ct_lung)
ct_t_stat  <- extract_t_stat(iv_ct_lung)
ct_lee_cv  <- lee_critical_value(ct_wald_f)

cat(sprintf("\n  Standalone Wald F  : %.2f\n", ct_wald_f))
cat(sprintf("  Lee et al. CV      : %.3f\n",  ct_lee_cv))
cat(sprintf("  Second-stage t     : %.3f\n",  ct_t_stat))
cat(sprintf("  Passes tF?         : %s\n",
            ifelse(abs(ct_t_stat) > ct_lee_cv, "YES", "NO")))

# ------------------------------------------------------------
# Part B.5: Full spectrum — Wald F vs Lee CV plot
# ------------------------------------------------------------

library(ggplot2)

# Smooth Lee et al. CV curve
f_seq <- seq(10, 120, by = 0.5)
lee_curve <- data.frame(
  F_stat = f_seq,
  lee_cv = sapply(f_seq, lee_critical_value)
)

# Service-level points
plot_data <- service_lee %>%
  filter(!is.na(fs_f), !is.na(t_stat)) %>%
  mutate(
    abs_t       = abs(t_stat),
    above_curve = abs_t > lee_cv,
    label_pt    = ifelse(sig_conventional | is_shoppable, SERVICE_GROUP, NA)
  )

ggplot() +
  # Lee et al. CV curve
  geom_line(data = lee_curve, aes(x = F_stat, y = lee_cv),
            color = "red", linewidth = 1, linetype = "solid") +
  # Stock-Yogo vertical line
  geom_vline(xintercept = 16.38, linetype = "dashed",
             color = "grey50", linewidth = 0.8) +
  # Conventional 1.96 horizontal line
  geom_hline(yintercept = 1.96, linetype = "dotted",
             color = "grey50", linewidth = 0.8) +
  # Service-level points
  geom_point(data = plot_data,
             aes(x = fs_f, y = abs_t,
                 color = is_shoppable,
                 shape = above_curve),
             size = 3, alpha = 0.8) +
  # Labels for notable services
  ggrepel::geom_label_repel(
    data = plot_data %>% filter(!is.na(label_pt)),
    aes(x = fs_f, y = abs_t, label = label_pt, color = is_shoppable),
    size = 2.8, max.overlaps = 15, show.legend = FALSE
  ) +
  scale_color_manual(
    values = c("FALSE" = "steelblue", "TRUE" = "darkred"),
    labels = c("Non-shoppable", "Shoppable"),
    name   = "Service type"
  ) +
  scale_shape_manual(
    values = c("FALSE" = 1, "TRUE" = 16),
    labels = c("Below Lee et al. CV", "Above Lee et al. CV"),
    name   = "tF assessment"
  ) +
  annotate("text", x = 17.5, y = max(plot_data$abs_t, na.rm=TRUE) * 0.95,
           label = "Stock-Yogo\n(F = 16.38)", size = 3, color = "grey40", hjust = 0) +
  annotate("text", x = max(f_seq) * 0.85, y = 2.05,
           label = "Conventional\n1.96 threshold", size = 3, color = "grey40") +
  annotate("text", x = max(f_seq) * 0.6, y = 2.35,
           label = "Lee et al.\ncritical value", size = 3, color = "red") +
  labs(
    x       = "First-stage Wald F (clustered)",
    y       = "Second-stage |t-statistic|",
    title   = "Lee et al. (2022) tF Assessment by Service Group",
    caption = paste0("Red curve: Lee et al. (2022) Table 1 critical values (5% level).\n",
                     "Points above red curve are significant under tF correction.\n",
                     "Instrument: system_peer_pressure_county_9m, two-way clustering.")
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")


# ------------------------------------------------------------
# Part B.6: Summary table for paper
# ------------------------------------------------------------

cat("\n=== SECTION 6: Paper-Ready Summary Table ===\n")
cat("Pooled specifications\n\n")

summary_pooled <- data.frame(
  Specification      = c("Baseline", "System × Month FE", "System Linear Trends"),
  Wald_F             = c(f_baseline, f_sysmth, f_trends),
  SY_pass            = c(f_baseline, f_sysmth, f_trends) > 16.38,
  Lee_CV             = c(cv_baseline, cv_sysmth, cv_trends),
  Second_stage_t     = c(t_baseline, t_sysmth, t_trends),
  Lee_tF_pass        = abs(c(t_baseline, t_sysmth, t_trends)) > 
    c(cv_baseline, cv_sysmth, cv_trends)
) %>%
  mutate(across(where(is.numeric), ~round(., 3)))

print(summary_pooled, row.names = FALSE)

cat("\n\nService-level significant results\n\n")
service_lee %>%
  filter(sig_conventional) %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>%
  select(SERVICE_GROUP, estimate_pct, t_stat, fs_f, lee_cv, 
         passes_SY, sig_conventional, passes_lee_tF) %>%
  print(row.names = FALSE)

cat("\n=== DONE ===\n")




cat("=== CURRENT PRIMARY SPEC — VERIFIED NUMBERS ===\n")
cat(sprintf("Coefficient : %.4f\n", coef(iv_county_pooled)["fit_n_prior_posters"]))
cat(sprintf("SE          : %.4f\n", se(iv_county_pooled)["fit_n_prior_posters"]))
cat(sprintf("t-stat      : %.3f\n", coef(iv_county_pooled)["fit_n_prior_posters"]/se(iv_county_pooled)["fit_n_prior_posters"]))
cat(sprintf("p-value     : %.4f\n", pvalue(iv_county_pooled)["fit_n_prior_posters"]))
cat(sprintf("Wald F      : %.2f\n", fitstat(iv_county_pooled,"ivwald")[["ivwald1::n_prior_posters"]]$stat))
cat(sprintf("N obs       : %d\n", nobs(iv_county_pooled)))


# ----------------------------------------------------------------------------
# Part C: CMS Enforcement Exposure as Control
#
# Purpose: If CMS enforcement drives compliance *and* directly
# affects prices (risk aversion channel), controlling for
# enforcement exposure isolates the peer-cascade mechanism.
# Uses your existing enforcement columns — fine_roll_9m_lag
# is the closest analogue to the 9m trailing instrument window.
# ----------------------------------------------------------------------------

# Control for same-hospital enforcement exposure
iv_enf_control <- feols(
  ln_median_price ~ ln_total_beds + warning_roll_6m |
    market_id + post_month |
    n_prior_posters ~ system_peer_pressure_county_9m,
  data    = df_iv_county,
  cluster = CL
)

# Also control for any system-peer enforcement (county-level cumulative)
# Build system peer enforcement: mean enforcement exposure of same-system
# out-of-county peers, mirroring the instrument construction logic
df_iv_county <- df_iv_county %>%
  group_by(post_month, HEALTH_SYSTEM_NAME) %>%
  mutate(
    system_peer_enforcement = mean(fine_cum[COUNTY_STATE != COUNTY_STATE], na.rm = TRUE),
    system_peer_any_enf     = mean(any_enforcement[COUNTY_STATE != COUNTY_STATE], na.rm = TRUE)
  ) %>%
  ungroup()

# Safer version using lag variables already in data
iv_enf_control2 <- feols(
  ln_median_price ~ ln_total_beds + warning_roll_6m_lag + 
    any_fine_or_warning_roll_9m |
    market_id + post_month |
    n_prior_posters ~ system_peer_pressure_county_9m,
  data    = df_iv_county,
  cluster = CL
)

# ----------------------------------------------------------------------------
# OUTPUT_TABLE_R2_ENFORCEMENT
# ----------------------------------------------------------------------------
# Robustness table adding enforcement controls.
cat("\n=== TABLE R2: Enforcement Controls ===\n")
etable(
  iv_baseline, iv_enf_control, iv_enf_control2,
  keep    = "%fit_n_prior_posters",
  headers = c("Baseline", "+ Enforcement Controls", "+ Any Fine/Warning"),
  fitstat = ~ n + ivwald1
)


# ----------------------------------------------------------------------------
# Part D: Conley-Hansen-Rossi Bounds on the Meta-Regression Gradient
# ----------------------------------------------------------------------------
# 
# The meta-regression estimates: iv_g = alpha + beta * shoppable_g + e_g
# 
# If the instrument has a direct effect delta on ln_price (violating
# exclusion), then each service-group IV estimate is biased:
#   iv_adj_g = (rf_g - delta) / fs_g = iv_g - delta / fs_g
#
# Since fs_g varies by service, the correction delta/fs_g also varies.
# We sweep delta in absolute terms (log units), recompute adjusted IV
# estimates for each service, re-run the meta-regression, and track
# the shoppability gradient.
#
# Key insight: a UNIFORM direct effect (delta constant across services)
# shifts all IV estimates but does NOT necessarily destroy the gradient,
# because the gradient reflects the *differential* between shoppable
# and non-shoppable services. The gradient is only threatened if the
# direct effect is correlated with shoppability — which I test.
# ============================================================

# Step 1: Reconstruct RF and FS for each service group
# I need rf_g and fs_g to compute the CHR correction per service

service_groups_all <- unique(df_iv_county$SERVICE_GROUP)

rf_fs_by_service <- do.call(rbind, lapply(service_groups_all, function(sg) {
  df_sub <- df_iv_county %>% filter(SERVICE_GROUP == sg)
  if (nrow(df_sub) < 100) return(NULL)
  
  # Reduced form
  rf_fit <- tryCatch(
    feols(ln_median_price ~ system_peer_pressure_county_9m + ln_total_beds |
            market_id + post_month,
          data = df_sub, cluster = CL),
    error = function(e) NULL)
  
  # First stage
  fs_fit <- tryCatch(
    feols(n_prior_posters ~ system_peer_pressure_county_9m + ln_total_beds |
            market_id + post_month,
          data = df_sub, cluster = CL),
    error = function(e) NULL)
  
  if (is.null(rf_fit) | is.null(fs_fit)) return(NULL)
  
  rf_c <- tryCatch(coef(rf_fit)["system_peer_pressure_county_9m"], error=function(e) NA_real_)
  rf_s <- tryCatch(se(rf_fit)["system_peer_pressure_county_9m"],   error=function(e) NA_real_)
  fs_c <- tryCatch(coef(fs_fit)["system_peer_pressure_county_9m"], error=function(e) NA_real_)
  fs_s <- tryCatch(se(fs_fit)["system_peer_pressure_county_9m"],   error=function(e) NA_real_)
  
  data.frame(
    SERVICE_GROUP = sg,
    rf_coef = rf_c, rf_se = rf_s,
    fs_coef = fs_c, fs_se = fs_s,
    iv_point = rf_c / fs_c,
    iv_se    = sqrt((rf_s/fs_c)^2 + (rf_c*fs_s/fs_c^2)^2),
    row.names = NULL
  )
}))

cat(sprintf("\nRecovered RF/FS for %d service groups\n", nrow(rf_fs_by_service)))

# Step 2: Join with shoppability flags from meta_data
# (res_county is your existing service-level IV results from run_iv_by_service)
meta_chr <- rf_fs_by_service %>%
  left_join(
    res_county %>% 
      mutate(SERVICE_GROUP = as.character(SERVICE_GROUP)) %>%
      select(SERVICE_GROUP, se_pct, fs_f),
    by = "SERVICE_GROUP"
  ) %>%
  mutate(
    is_ultrasound  = grepl("Ultrasound", SERVICE_GROUP),
    is_mri         = grepl("MRI",        SERVICE_GROUP),
    is_xray        = grepl("X-Ray",      SERVICE_GROUP),
    is_biopsy      = grepl("Biopsy",     SERVICE_GROUP),
    is_ct_lung     = SERVICE_GROUP == "CT Lung",
    is_mammography = SERVICE_GROUP == "Mammography",
    is_shoppable   = is_ultrasound | is_ct_lung | is_mammography,
    weight         = ifelse(is.na(se_pct) | se_pct == 0, NA_real_, 1/(se_pct^2))
  ) %>%
  filter(is.finite(weight), !is.na(fs_coef), fs_coef > 0)

cat(sprintf("Service groups in CHR meta analysis: %d\n", nrow(meta_chr)))

# Baseline meta-regression gradient (should match your existing result ~-2.60)
meta_baseline <- lm(I(iv_point*100) ~ is_shoppable, 
                    data = meta_chr, weights = weight)
baseline_gradient <- coef(meta_baseline)["is_shoppableTRUE"]
cat(sprintf("\nBaseline shoppability gradient: %.3f pp\n", baseline_gradient))

# Step 3: Test whether direct effect is correlated with shoppability
# This is the key diagnostic. If shoppable services have systematically
# larger RF coefficients (in absolute terms), a uniform delta would
# differentially shrink them — threatening the gradient.
cat("\n=== Is the RF correlated with shoppability? ===\n")
rf_shoppable_test <- lm(rf_coef ~ is_shoppable, data = meta_chr, weights = weight)
summary(rf_shoppable_test)

cat("\n=== Is the FS correlated with shoppability? ===\n")
fs_shoppable_test <- lm(fs_coef ~ is_shoppable, data = meta_chr, weights = weight)
summary(fs_shoppable_test)

# Step 4: Sweep delta — uniform direct effect across all services
# Express delta as fraction of the median |RF| to make it interpretable
median_abs_rf <- median(abs(meta_chr$rf_coef), na.rm = TRUE)
cat(sprintf("\nMedian |RF coef| across services: %.5f\n", median_abs_rf))

# Grid: 0 to 100% of median RF (upward-biasing direction)
delta_grid_meta <- seq(0, median_abs_rf, length.out = 100)

gradient_sweep <- lapply(delta_grid_meta, function(delta) {
  meta_adj <- meta_chr %>%
    mutate(
      # Corrected RF: add delta (upward bias from direct effect)
      rf_adj   = rf_coef + delta,
      iv_adj   = rf_adj / fs_coef,
      iv_adj_pct = iv_adj * 100,
      # Recompute SE (same as before — delta is assumed known, not estimated)
      iv_adj_se = sqrt((rf_se/fs_coef)^2 + (rf_coef*fs_se/fs_coef^2)^2),
      weight_adj = 1 / (iv_adj_se * 100)^2
    ) %>%
    filter(is.finite(weight_adj), is.finite(iv_adj_pct))
  
  if (nrow(meta_adj) < 5) return(NULL)
  
  m_adj <- lm(iv_adj_pct ~ is_shoppable, data = meta_adj, weights = weight_adj)
  
  grad     <- coef(m_adj)["is_shoppableTRUE"]
  grad_se  <- summary(m_adj)$coefficients["is_shoppableTRUE", "Std. Error"]
  
  data.frame(
    delta       = delta,
    delta_pct   = delta / median_abs_rf * 100,
    gradient    = grad,
    grad_se     = grad_se,
    ci_lo       = grad - 1.96 * grad_se,
    ci_hi       = grad + 1.96 * grad_se,
    still_neg   = (grad + 1.96 * grad_se) < 0,
    row.names   = NULL
  )
}) %>% bind_rows()

# Thresholds
sign_flip_meta <- gradient_sweep %>% filter(gradient >= 0) %>% slice(1)
ci_flip_meta   <- gradient_sweep %>% filter(!still_neg) %>% slice(1)

cat(sprintf("\n=== Meta-Regression CHR Bounds ===\n"))
cat(sprintf("  Baseline shoppability gradient : %.3f pp\n", baseline_gradient))
cat(sprintf("  Gradient flips sign at delta   = %.1f%% of median |RF|\n",
            sign_flip_meta$delta_pct))
cat(sprintf("  95%% CI crosses zero at delta   = %.1f%% of median |RF|\n",
            ci_flip_meta$delta_pct))

# Selected rows
cat("\n  Selected sensitivity rows:\n")
gradient_sweep %>%
  filter(round(delta_pct) %in% c(0, 10, 20, 30, 40, 50, 75, 100)) %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>%
  select(delta_pct, gradient, grad_se, ci_lo, ci_hi, still_neg) %>%
  print(row.names = FALSE)

# Step 5: Non-uniform delta — proportional to shoppability (worst case)
# This is the most adversarial assumption: direct effect is proportionally
# LARGER for shoppable services (the exact threat to the gradient).
# If the gradient survives even this, the exclusion restriction argument is robust.
cat("\n=== Worst-Case: Proportional delta (larger for shoppable) ===\n")

scale_grid <- seq(0, 1, length.out = 50)  # scale of extra bias on shoppable

worst_case_sweep <- lapply(scale_grid, function(extra_scale) {
  meta_adj <- meta_chr %>%
    mutate(
      # Shoppable services get an extra (extra_scale * median_abs_rf) upward bias
      extra_bias = is_shoppable * extra_scale * median_abs_rf,
      rf_adj     = rf_coef + extra_bias,
      iv_adj     = rf_adj / fs_coef,
      iv_adj_pct = iv_adj * 100,
      iv_adj_se  = sqrt((rf_se/fs_coef)^2 + (rf_coef*fs_se/fs_coef^2)^2),
      weight_adj = 1 / (iv_adj_se * 100)^2
    ) %>%
    filter(is.finite(weight_adj), is.finite(iv_adj_pct))
  
  m_adj    <- lm(iv_adj_pct ~ is_shoppable, data = meta_adj, weights = weight_adj)
  grad     <- coef(m_adj)["is_shoppableTRUE"]
  grad_se  <- summary(m_adj)$coefficients["is_shoppableTRUE", "Std. Error"]
  
  data.frame(
    extra_scale  = extra_scale,
    extra_pct    = extra_scale * 100,
    gradient     = grad,
    grad_se      = grad_se,
    ci_lo        = grad - 1.96 * grad_se,
    ci_hi        = grad + 1.96 * grad_se,
    still_neg    = (grad + 1.96 * grad_se) < 0
  )
}) %>% bind_rows()

wc_sign_flip <- worst_case_sweep %>% filter(gradient >= 0) %>% slice(1)
wc_ci_flip   <- worst_case_sweep %>% filter(!still_neg) %>% slice(1)

wc_ci_flip$extra_pct
cat(sprintf("  Gradient sign flips when shoppable-specific bias = %.1f%% of median |RF|\n",
            wc_sign_flip$extra_pct))
cat(sprintf("  95%% CI crosses zero at shoppable-specific bias  = %.1f%% of median |RF|\n",
            wc_ci_flip$extra_pct))

# Step 6: Combined plot
library(patchwork)

# ----------------------------------------------------------------------------
# OUTPUT_FIGURE_CHR_UNIFORM_PANEL
# ----------------------------------------------------------------------------
# Sensitivity panel p1: uniform direct effect scenario.
p1 <- ggplot(gradient_sweep, aes(x = delta_pct)) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.2, fill = "steelblue") +
  geom_line(aes(y = gradient), color = "steelblue", linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_vline(xintercept = ci_flip_meta$delta_pct, 
             linetype = "dotted", color = "grey40") +
  labs(
    x = "Uniform direct effect (% of median |RF|)",
    y = "Shoppability gradient (pp)",
    title = "Uniform Direct Effect",
    subtitle = sprintf("CI crosses zero at %.0f%% contamination",
                       ci_flip_meta$delta_pct)
  ) +
  theme_minimal()

# ----------------------------------------------------------------------------
# OUTPUT_FIGURE_CHR_ADVERSARIAL_PANEL
# ----------------------------------------------------------------------------
# Sensitivity panel p2: shoppable-specific direct effect scenario.
p2 <- ggplot(worst_case_sweep, aes(x = extra_pct)) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.2, fill = "darkred") +
  geom_line(aes(y = gradient), color = "darkred", linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_vline(xintercept = wc_ci_flip$extra_pct, 
             linetype = "dotted", color = "grey40") +
  labs(
    x = "Shoppable-specific extra bias (% of median |RF|)",
    y = "Shoppability gradient (pp)",
    title = "Worst Case: Shoppable-Specific Direct Effect",
    subtitle = sprintf("CI crosses zero at %.0f%% contamination",
                       wc_ci_flip$extra_pct)
  ) +
  theme_minimal()

p1 + p2 +
  plot_annotation(
    title   = "Conley-Hansen-Rossi Sensitivity: Shoppability Gradient",
    caption = "Left: uniform direct effect on all services. Right: direct effect concentrated on shoppable services (adversarial case)."
  )




# ============================================================
# CONLEY BOUNDS — CT LUNG
# ============================================================

df_ct <- df_iv_county %>% filter(SERVICE_GROUP == "CT Lung")

# Reduced form for CT Lung
rf_ct <- feols(
  ln_median_price ~ system_peer_pressure_county_9m + ln_total_beds |
    market_id + post_month,
  data    = df_ct,
  cluster = CL
)

# First stage for CT Lung
fs_ct <- feols(
  n_prior_posters ~ system_peer_pressure_county_9m + ln_total_beds |
    market_id + post_month,
  data    = df_ct,
  cluster = CL
)

rf_coef_ct  <- coef(rf_ct)["system_peer_pressure_county_9m"]
rf_se_ct    <- se(rf_ct)["system_peer_pressure_county_9m"]
fs_coef_ct  <- coef(fs_ct)["system_peer_pressure_county_9m"]
fs_se_ct    <- se(fs_ct)["system_peer_pressure_county_9m"]
iv_point_ct <- rf_coef_ct / fs_coef_ct

cat(sprintf("\n=== CT Lung: Reduced Form and First Stage ===\n"))
cat(sprintf("  RF coef  : %.5f (SE: %.5f, t: %.2f)\n",
            rf_coef_ct, rf_se_ct, rf_coef_ct/rf_se_ct))
cat(sprintf("  FS coef  : %.5f (SE: %.5f)\n", fs_coef_ct, fs_se_ct))
cat(sprintf("  IV point : %.4f (%.2f%%)\n", iv_point_ct, iv_point_ct*100))

# Conley bounds for CT Lung
# For CT Lung the RF is likely negative and significant, so 
# delta shifts it toward zero. We ask: how much direct positive
# contamination would flip the sign?
# Extend grid to 150% for CT Lung if needed
delta_grid_ct2 <- seq(0, abs(rf_coef_ct) * 1.5, length.out = 200)
# Then re-run chr_ct with this grid and recheck sign_flip_ct
chr_ct2 <- lapply(delta_grid_ct2, function(delta) {
  rf_adj    <- rf_coef_ct + delta
  iv_adj    <- rf_adj / fs_coef_ct
  iv_adj_se <- sqrt((rf_se_ct/fs_coef_ct)^2 + 
                      (rf_coef_ct * fs_se_ct / fs_coef_ct^2)^2)
  data.frame(
    delta_pct       = delta / abs(rf_coef_ct) * 100,
    iv_adj_pct      = iv_adj * 100,
    ci_lo_pct       = (iv_adj - 1.96 * iv_adj_se) * 100,
    ci_hi_pct       = (iv_adj + 1.96 * iv_adj_se) * 100,
    still_neg_ci_hi = (iv_adj + 1.96 * iv_adj_se) < 0
  )
}) %>% bind_rows()

sign_flip_ct2 <- chr_ct2 %>%
  filter(iv_adj_pct >= 0) %>%
  slice(1)

ci_flip_ct2 <- chr_ct2 %>%
  filter(!still_neg_ci_hi) %>%
  slice(1)

cat(sprintf("CT Lung point estimate flips at: %.1f%% of RF\n",
            sign_flip_ct2$delta_pct))

cat(sprintf("  95%% CI crosses zero when delta = %.1f%% of RF\n",
            ci_flip_ct2$delta_pct))

cat("\n  Selected Conley sensitivity rows (CT Lung):\n")
chr_ct2 %>%
  filter(round(delta_pct) %in% c(0, 10, 20, 30, 40, 50, 60, 75)) %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>%
  select(delta_pct, iv_adj_pct, ci_lo_pct, ci_hi_pct, still_neg_ci_hi) %>%
  print(row.names = FALSE)

# Rename/store the CT Lung CHR object consistently
chr_ct <- chr_ct2

# Thresholds
sign_flip_ct <- chr_ct %>%
  filter(iv_adj_pct >= 0) %>%
  slice(1)

ci_flip_ct <- chr_ct %>%
  filter(!still_neg_ci_hi) %>%
  slice(1)

# Selected sensitivity rows
chr_ct %>%
  filter(round(delta_pct) %in% c(0, 10, 20, 30, 40, 50, 60, 75)) %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>%
  select(delta_pct, iv_adj_pct, ci_lo_pct, ci_hi_pct, still_neg_ci_hi) %>%
  print(row.names = FALSE)

# Plot
library(ggplot2)

ggplot(chr_ct, aes(x = delta_pct)) +
  geom_ribbon(
    aes(ymin = ci_lo_pct, ymax = ci_hi_pct),
    alpha = 0.2,
    fill = "steelblue"
  ) +
  geom_line(
    aes(y = iv_adj_pct),
    color = "steelblue",
    linewidth = 1
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "red"
  ) +
  geom_vline(
    xintercept = sign_flip_ct$delta_pct,
    linetype = "dotted",
    color = "grey40"
  ) +
  geom_vline(
    xintercept = ci_flip_ct$delta_pct,
    linetype = "dashed",
    color = "grey40"
  ) +
  labs(
    x = "Direct upward bias in instrument (% of reduced-form magnitude)",
    y = "IV estimate (% price effect, CT Lung)",
    title = "Conley-Hansen-Rossi Bounds: CT Lung",
    subtitle = sprintf(
      "Sign flips at %.0f%% contamination; upper CI crosses zero at %.0f%%",
      sign_flip_ct$delta_pct,
      ci_flip_ct$delta_pct
    ),
    caption = "Direct effect assumed to be upward-biasing (instrument → higher prices)."
  ) +
  theme_minimal()


# ============================================================
# EXTENDED WORST-CASE META-REGRESSION GRID
# ============================================================

# Extend worst-case grid to 300% to find the actual CI flip
scale_grid_ext <- seq(0, 3, length.out = 200)

worst_case_ext <- lapply(scale_grid_ext, function(extra_scale) {
  meta_adj <- meta_chr %>%
    mutate(
      extra_bias = is_shoppable * extra_scale * median_abs_rf,
      rf_adj     = rf_coef + extra_bias,
      iv_adj     = rf_adj / fs_coef,
      iv_adj_pct = iv_adj * 100,
      iv_adj_se  = sqrt((rf_se / fs_coef)^2 +
                          (rf_coef * fs_se / fs_coef^2)^2),
      weight_adj = 1 / (iv_adj_se * 100)^2
    ) %>%
    filter(is.finite(weight_adj), is.finite(iv_adj_pct))
  
  m_adj   <- lm(iv_adj_pct ~ is_shoppable, data = meta_adj, weights = weight_adj)
  grad    <- coef(m_adj)["is_shoppableTRUE"]
  grad_se <- summary(m_adj)$coefficients["is_shoppableTRUE", "Std. Error"]
  
  data.frame(
    extra_scale = extra_scale,
    extra_pct   = extra_scale * 100,
    gradient    = grad,
    grad_se     = grad_se,
    ci_lo       = grad - 1.96 * grad_se,
    ci_hi       = grad + 1.96 * grad_se,
    still_neg   = (grad + 1.96 * grad_se) < 0
  )
}) %>%
  bind_rows()

wc_sign_flip_ext <- worst_case_ext %>%
  filter(gradient >= 0) %>%
  slice(1)

wc_ci_flip_ext <- worst_case_ext %>%
  filter(!still_neg) %>%
  slice(1)

if (nrow(wc_sign_flip_ext) > 0) {
  cat(sprintf(
    "Worst-case gradient sign flip : %.1f%% of median |RF| on shoppable\n",
    wc_sign_flip_ext$extra_pct
  ))
} else {
  cat("Worst-case gradient sign flip : Not reached within 300% of median |RF|\n")
}

if (nrow(wc_ci_flip_ext) > 0) {
  cat(sprintf(
    "Worst-case CI crosses zero    : %.1f%% of median |RF| on shoppable\n",
    wc_ci_flip_ext$extra_pct
  ))
} else {
  cat("Worst-case CI crosses zero    : Not reached within 300% of median |RF|\n")
}

# Print selected rows to see trajectory
worst_case_ext %>%
  filter(round(extra_pct) %in% c(0, 50, 100, 150, 200, 250, 300)) %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>%
  select(extra_pct, gradient, grad_se, ci_lo, ci_hi, still_neg) %>%
  print(row.names = FALSE)


# ============================================================
# FINAL CT LUNG CHR SUMMARY
# ============================================================

# define ct_row if it does not already exist
ct_row <- res_county %>%
  mutate(SERVICE_GROUP = as.character(SERVICE_GROUP)) %>%
  filter(SERVICE_GROUP == "CT Lung") %>%
  slice(1)

cat("\n=== CT Lung CHR Summary — Final Numbers ===\n")
cat(sprintf("  IV point estimate          : %.3f pp\n", ct_row$estimate_pct))
cat(sprintf("  RF t-stat                  : %.2f\n",   rf_coef_ct / rf_se_ct))
cat(sprintf("  Point estimate flips at    : %.1f%% of RF\n", sign_flip_ct$delta_pct))
cat(sprintf("  95%% CI crosses zero at     : %.1f%% of RF\n", ci_flip_ct$delta_pct))



# ============================================================
# CHR SENSITIVITY FIGURE
# fig:chr_sensitivity — Three panels:
#   Panel A: CT Lung bounds
#   Panel B: Meta-regression gradient, uniform contamination
#   Panel C: Meta-regression gradient, adversarial contamination
#
# Required objects (already in session):
#   chr_ct2         — CT Lung CHR bounds, extended grid (0–150% of RF)
#   gradient_sweep  — uniform contamination sweep of meta-regression
#   worst_case_ext  — adversarial contamination sweep (0–300%)
#
# Required packages:
#   ggplot2, patchwork
# ============================================================

# ---- Shared theme ----
chr_theme <- theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border       = element_rect(color = "grey80", fill = NA, linewidth = 0.4),
    axis.title         = element_text(size = 9),
    axis.text          = element_text(size = 8),
    plot.title         = element_text(size = 9, face = "bold", hjust = 0),
    plot.subtitle      = element_text(size = 8, color = "grey40", hjust = 0),
    legend.position    = "none"
  )

line_col   <- "#8B1A1A"   # dark garnet for estimate line
ribbon_col <- "#C97070"   # lighter for CI ribbon

# ============================================================
# PANEL A: CT Lung
# x: delta_pct (contamination as % of RF), range 0–110
# y: iv_adj_pct (IV estimate in %), with CI ribbon
# Reference lines: zero (horizontal), 14.4% CI flip, 100.3% sign flip
# ============================================================

# Subset to 0–110% for clean display
ct_plot_data <- chr_ct2 %>%
  filter(delta_pct <= 110) %>%
  mutate(
    ci_lo_pct = pmax(ci_lo_pct, -16),   # clip ribbon for legibility
    ci_hi_pct = pmin(ci_hi_pct,  10)
  )

# ----------------------------------------------------------------------------
# OUTPUT_FIGURE_CHR_PANEL_A_CT_LUNG
# ----------------------------------------------------------------------------
# CHR sensitivity panel A: CT Lung bounds.
p_ct <- ggplot(ct_plot_data, aes(x = delta_pct)) +
  # CI ribbon
  geom_ribbon(aes(ymin = ci_lo_pct, ymax = ci_hi_pct),
              fill = ribbon_col, alpha = 0.25) +
  # Estimate line
  geom_line(aes(y = iv_adj_pct),
            color = line_col, linewidth = 0.9) +
  # Zero line
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey30", linewidth = 0.5) +
  # CI flip threshold (14.4%)
  geom_vline(xintercept = 14.4, linetype = "dotted",
             color = "grey40", linewidth = 0.5) +
  # Sign flip threshold (100.3%)
  geom_vline(xintercept = 100.3, linetype = "dotted",
             color = "grey40", linewidth = 0.5) +
  # Annotations
  annotate("text", x = 14.4 + 1.5, y = -13.5,
           label = "CI flip\n(14.4%)", size = 2.8,
           color = "grey35", hjust = 0) +
  annotate("text", x = 100.3 + 1.5, y = -13.5,
           label = "Sign flip\n(100.3%)", size = 2.8,
           color = "grey35", hjust = 0) +
  scale_x_continuous(
    limits = c(0, 110),
    breaks = c(0, 25, 50, 75, 100),
    labels = c("0%", "25%", "50%", "75%", "100%")
  ) +
  scale_y_continuous(
    breaks = seq(-12, 4, by = 4),
    labels = function(x) paste0(x, "%")
  ) +
  labs(
    title    = "Panel A: CT Lung",
    subtitle = "IV estimate = \u22125.74%; RF t = \u22122.52",
    x        = "Direct effect (% of reduced-form coefficient)",
    y        = "IV estimate (%)"
  ) +
  chr_theme

# ============================================================
# PANEL B: Meta-Regression Gradient — Uniform Contamination
# x: delta_pct (0–100% of median |RF|)
# y: gradient (pp), with CI ribbon
# Reference line: zero only (CI never crosses)
# Add annotation noting CI never crosses zero
# ============================================================

# ----------------------------------------------------------------------------
# OUTPUT_FIGURE_CHR_PANEL_B_UNIFORM
# ----------------------------------------------------------------------------
# CHR sensitivity panel B: uniform contamination gradient.
p_uniform <- ggplot(gradient_sweep, aes(x = delta_pct)) +
  # CI ribbon
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi),
              fill = ribbon_col, alpha = 0.25) +
  # Estimate line
  geom_line(aes(y = gradient),
            color = line_col, linewidth = 0.9) +
  # Zero line
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey30", linewidth = 0.5) +
  # Annotation: CI never crosses zero
  annotate("text", x = 50, y = -0.6,
           label = "95% CI remains\nbelow zero throughout",
           size = 2.8, color = "grey35", hjust = 0.5) +
  annotate("segment", x = 50, xend = 50,
           y = -0.85, yend = -1.60,
           arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
           color = "grey45", linewidth = 0.4) +
  scale_x_continuous(
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100),
    labels = c("0%", "25%", "50%", "75%", "100%")
  ) +
  scale_y_continuous(
    limits = c(-4.2, 0.5),
    breaks = seq(-4, 0, by = 1),
    labels = function(x) paste0(x, " pp")
  ) +
  labs(
    title    = "Panel B: Gradient \u2014 Uniform Contamination",
    subtitle = "Direct effect uniform across all 46 service groups",
    x        = "Uniform direct effect (% of median |RF|)",
    y        = "Shoppability gradient (pp)"
  ) +
  chr_theme

# ============================================================
# PANEL C: Meta-Regression Gradient — Adversarial Contamination
# x: extra_pct (0–300% of median |RF|, shoppable-targeted)
# y: gradient (pp), with CI ribbon
# Reference lines: zero, 149.2% CI flip, 223.1% sign flip
# ============================================================

# ----------------------------------------------------------------------------
# OUTPUT_FIGURE_CHR_PANEL_C_ADVERSE
# ----------------------------------------------------------------------------
# CHR sensitivity panel C: adversarial shoppable-targeted contamination.
p_adverse <- ggplot(worst_case_ext, aes(x = extra_pct)) +
  # CI ribbon
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi),
              fill = ribbon_col, alpha = 0.25) +
  # Estimate line
  geom_line(aes(y = gradient),
            color = line_col, linewidth = 0.9) +
  # Zero line
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey30", linewidth = 0.5) +
  # CI flip threshold (149.2%)
  geom_vline(xintercept = 149.2, linetype = "dotted",
             color = "grey40", linewidth = 0.5) +
  # Sign flip threshold (223.1%)
  geom_vline(xintercept = 223.1, linetype = "dotted",
             color = "grey40", linewidth = 0.5) +
  # Annotations
  annotate("text", x = 149.2 + 4, y = -3.6,
           label = "CI flip\n(149.2%)", size = 2.8,
           color = "grey35", hjust = 0) +
  annotate("text", x = 223.1 + 4, y = -3.6,
           label = "Sign flip\n(223.1%)", size = 2.8,
           color = "grey35", hjust = 0) +
  scale_x_continuous(
    limits = c(0, 310),
    breaks = c(0, 50, 100, 150, 200, 250, 300),
    labels = c("0%", "50%", "100%", "150%", "200%", "250%", "300%")
  ) +
  scale_y_continuous(
    breaks = seq(-4, 2, by = 1),
    labels = function(x) paste0(x, " pp")
  ) +
  labs(
    title    = "Panel C: Gradient \u2014 Adversarial Contamination",
    subtitle = "Direct effect concentrated exclusively on shoppable services",
    x        = "Shoppable-targeted direct effect (% of median |RF|)",
    y        = "Shoppability gradient (pp)"
  ) +
  chr_theme

# ============================================================
# COMBINE AND EXPORT
# ============================================================

fig_chr <- p_ct | p_uniform | p_adverse

fig_chr_final <- fig_chr +
  plot_annotation(
    caption = paste0(
      "Shaded bands are 95% confidence intervals. ",
      "Red dashed lines mark zero. ",
      "Dotted vertical lines mark CI flip and sign-flip thresholds."
    ),
    theme = theme(
      plot.caption = element_text(size = 7.5, color = "grey40",
                                  hjust = 0, margin = margin(t = 6))
    )
  )

# ----------------------------------------------------------------------------
# OUTPUT_FIGURE_CHR_SENSITIVITY
# ----------------------------------------------------------------------------
# Prints fig_chr_final: three-panel CHR sensitivity figure.
print(fig_chr_final)

# ----------------------------------------------------------------------------
# Part E: Summary Comparison Table — All Robustness Specs
# ----------------------------------------------------------------------------

robustness_summary <- bind_rows(
  iv_row(iv_baseline,   "1. Baseline (market_id + month FE)"),
  iv_row(iv_sys_month,  "2. + System × Month FE"),
  iv_row(iv_sys_trends, "3. + System Linear Trends"),
  iv_row(iv_enf_control,"4. + Enforcement Controls"),
  iv_row(iv_enf_control2,"5. + Any Fine/Warning Control")
)

# ----------------------------------------------------------------------------
# OUTPUT_TABLE_R3_ROBUSTNESS_SUMMARY
# ----------------------------------------------------------------------------
# Robustness summary table across exclusion-restriction stress tests.
cat("\n=== TABLE R3: Robustness Summary ===\n")
robustness_summary %>%
  mutate(
    `IV Est. (%)` = sprintf("%.3f", iv_est_pct),
    `SE (%)`      = sprintf("%.3f", iv_se_pct),
    `p-value`     = sprintf("%.3f", iv_pval),
    `Wald F`      = sprintf("%.1f",  wald_f),
    `FS Coef`     = sprintf("%.4f",  fs_coef)
  ) %>%
  select(spec, `IV Est. (%)`, `SE (%)`, `p-value`, `Wald F`, `FS Coef`) %>%
  print(row.names = FALSE)



######## Section 19: ML Extension — Post-Double-Selection Lasso ###########
# ----------------------------------------------------------------------------
# Now uses county data (df_iv_county) with 9m instrument
# Machine-learning extension: post-double-selection LASSO IV using county data
#  and the 9m instrument.
library(hdm)
controls <- c("ln_total_beds","Share_ForProfit","Share_NonProfit",
              "uninsured_rate","medicaid_share","poverty_rate","log_median_income",
              "age65plus_share","college_share","black_share","hispanic_share",
              "employment_rate","homeowner_rate","population")

X_controls <- df_iv_county %>%
  select(any_of(controls)) %>%
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()

keep_pds <- complete.cases(X_controls) &
  !is.na(df_iv_county$ln_median_price) &
  !is.na(df_iv_county$n_prior_posters) &
  !is.na(df_iv_county$system_peer_pressure_county_9m)

df_pds <- df_iv_county[keep_pds,]
X_pds  <- X_controls[keep_pds,]

pds_result <- tryCatch(
  rlassoIV(x=X_pds, d=as.numeric(df_pds$n_prior_posters),
           y=as.numeric(df_pds$ln_median_price),
           z=as.numeric(df_pds$system_peer_pressure_county_9m),
           select.X=TRUE, select.Z=TRUE),
  error=function(e){ message("PDS LASSO failed: ", conditionMessage(e)); NULL })

if (!is.null(pds_result)) {
  cat("\n=== Post-Double-Selection LASSO IV (County, 9m instrument) ===\n")
  summary(pds_result)
}


######## Section 20: LaTeX Tables for Paper ###########
# ----------------------------------------------------------------------------
# LaTeX and CSV export section. This is the main place to find paper-ready 
# table output and machine-readable result files.
# Table 1: Summary Statistics (county level, primary instrument)
summ_df <- df_iv_county %>%
  select(MEDIAN_PRICE, MEAN_PRICE, P25, P75,
         n_prior_posters, system_peer_pressure_county_9m,
         Share_ForProfit, Share_Government, Share_NonProfit,
         TOTAL_BEDS, population, black_share,
         age65plus_share, uninsured_rate, poverty_rate,
         college_share, employment_rate) %>%
  mutate(across(everything(), as.numeric)) %>%
  as.data.frame()

# ----------------------------------------------------------------------------
# OUTPUT_TABLE_01_SUMSTATS
# ----------------------------------------------------------------------------
# LaTeX Table 1: summary statistics.
stargazer(
  summ_df, type="latex", title="Summary Statistics",
  label="tab:sumstats", digits=2,
  omit.summary.stat=c("min","max"),
  covariate.labels=c(
    "Median Price","Mean Price","P25","P75",
    "Prior Posters","System Peer Pressure (9m)",
    "For-Profit Share","Government Share","Non-Profit Share",
    "Total Beds","Population","Black Share",
    "Age 65+ Share","Uninsured Rate","Poverty Rate",
    "College Share","Employment Rate")
)

# Table 3: Pooled OLS vs IV — all outcomes
# ----------------------------------------------------------------------------
# OUTPUT_TABLE_03_LATEX_POOLED_OLS_IV
# ----------------------------------------------------------------------------
# LaTeX pooled OLS vs IV all-outcomes table.
etable(
  ols_county_mean,   iv_county_mean,
  ols_county_p25,    iv_county_p25,
  ols_county_median, iv_county_pooled,
  ols_county_p75,    iv_county_p75,
  ols_county_max,    iv_county_max,
  ols_county_min,    iv_county_min,
  ols_county_iqr,    iv_county_iqr,
  keep="%n_prior_posters|%fit_n_prior_posters",
  headers=c("OLS","IV","OLS","IV","OLS","IV","OLS","IV","OLS","IV","OLS","IV","OLS","IV"),
  fitstat=c("n","r2","ivwald"), tex=TRUE,
  title="Effect of Price Transparency on Negotiated Prices: OLS and IV Estimates (County, 9m instrument)",
  label="tab:pooled_ols_iv"
)

# Table 6: Meta-regression
meta_simple_coefs   <- coef(summary(meta_county$simple))
meta_granular_coefs <- coef(summary(meta_county$granular))
# ----------------------------------------------------------------------------
# OUTPUT_TABLE_06_METAREG_LATEX
# ----------------------------------------------------------------------------
# LaTeX shoppability meta-regression table.
cat("\n=== LaTeX: Table 6 Meta-Regression (county, 9m) ===\n")
cat("\\begin{table}[h]\n\\centering\n")
cat("\\caption{Shoppability Meta-Regression (County, 9-month instrument)}\n")
cat("\\label{tab:meta_regression}\n")
cat("\\begin{tabular}{lcc}\n\\toprule\n")
cat(" & Simple & Granular \\\\\n\\midrule\n")
all_rows <- union(rownames(meta_simple_coefs), rownames(meta_granular_coefs))
for (r in all_rows) {
  if (r=="(Intercept)") next
  s_est <- if(r%in%rownames(meta_simple_coefs))
    sprintf("%.3f%s (%.3f)",meta_simple_coefs[r,"Estimate"],
            add_stars(meta_simple_coefs[r,"Pr(>|t|)"]),meta_simple_coefs[r,"Std. Error"]) else ""
  g_est <- if(r%in%rownames(meta_granular_coefs))
    sprintf("%.3f%s (%.3f)",meta_granular_coefs[r,"Estimate"],
            add_stars(meta_granular_coefs[r,"Pr(>|t|)"]),meta_granular_coefs[r,"Std. Error"]) else ""
  lab <- gsub("TRUE$","",gsub("is_","",r))
  cat(sprintf("%s & %s & %s \\\\\n",lab,s_est,g_est))
}
cat(sprintf("\\midrule\n$R^2$ & %.3f & %.3f \\\\\n",
            summary(meta_county$simple)$r.squared, summary(meta_county$granular)$r.squared))
cat("$N$ & 46 & 46 \\\\\n")
cat("\\bottomrule\n")
cat("\\multicolumn{3}{l}{\\footnotesize{Weighted by $1/\\text{SE}^2$. *** $p<0.01$, ** $p<0.05$, * $p<0.10$.}}\n")
cat("\\end{tabular}\n\\end{table}\n")



# ── Extract all values for Table: Meta-Regression ──────────────────────────

m1 <- meta_scheme_results[["scheme_theory_v2"]]$simple
m2 <- meta_scheme_results[["scheme_theory_v2"]]$three_way
m3 <- meta_scheme_results[["scheme_split_nonshop"]]$split

# ----------------------------------------------------------------------------
# FUNCTION extract_table_row()
# ----------------------------------------------------------------------------
# Extracts coefficient, standard error, p-value, confidence interval, and 
# stars for custom table construction.
extract_table_row <- function(model, term_name) {
  ct <- coef(summary(model))
  if (!term_name %in% rownames(ct)) return(list(est=NA, se=NA, p=NA))
  list(
    est   = round(ct[term_name, "Estimate"],    4),
    se    = round(ct[term_name, "Std. Error"],  4),
    p     = round(ct[term_name, "Pr(>|t|)"],    4),
    stars = add_stars(ct[term_name, "Pr(>|t|)"])
  )
}

# ----------------------------------------------------------------------------
# OUTPUT_TABLE_VALUES_METAREG
# ----------------------------------------------------------------------------
# Console extraction of meta-regression values for paper table construction.
cat("\n=== TABLE VALUES: Meta-Regression (Theory V2) ===\n\n")
cat("--- Column 1: Simple ---\n")
r <- extract_table_row(m1, "is_shoppableTRUE")
cat(sprintf("  Shoppable:    %8.4f%s (%8.4f)   p = %s\n", r$est, r$stars, r$se, r$p))
r <- extract_table_row(m1, "(Intercept)")
cat(sprintf("  Intercept:    %8.4f%s (%8.4f)   p = %s\n", r$est, r$stars, r$se, r$p))
cat(sprintf("  R2:           %.4f\n", summary(m1)$r.squared))
cat(sprintf("  Adj R2:       %.4f\n", summary(m1)$adj.r.squared))

cat("\n--- Column 2: Three-Way ---\n")
r <- extract_table_row(m2, "is_shoppableTRUE")
cat(sprintf("  Shoppable:    %8.4f%s (%8.4f)   p = %s\n", r$est, r$stars, r$se, r$p))
r <- extract_table_row(m2, "is_nonshoppableTRUE")
cat(sprintf("  Non-Shoppable:%8.4f%s (%8.4f)   p = %s\n", r$est, r$stars, r$se, r$p))
r <- extract_table_row(m2, "(Intercept)")
cat(sprintf("  Intercept:    %8.4f%s (%8.4f)   p = %s\n", r$est, r$stars, r$se, r$p))
cat(sprintf("  R2:           %.4f\n", summary(m2)$r.squared))
cat(sprintf("  Adj R2:       %.4f\n", summary(m2)$adj.r.squared))

cat("\n--- Column 3: Split Non-Shoppable ---\n")
r <- extract_table_row(m3, "is_shoppableTRUE")
cat(sprintf("  Shoppable:        %8.4f%s (%8.4f)   p = %s\n", r$est, r$stars, r$se, r$p))
r <- extract_table_row(m3, "is_nonshoppable_mriTRUE")
cat(sprintf("  Non-Shop (MRI):   %8.4f%s (%8.4f)   p = %s\n", r$est, r$stars, r$se, r$p))
r <- extract_table_row(m3, "is_nonshoppable_procTRUE")
cat(sprintf("  Non-Shop (Proc):  %8.4f%s (%8.4f)   p = %s\n", r$est, r$stars, r$se, r$p))
r <- extract_table_row(m3, "(Intercept)")
cat(sprintf("  Intercept:        %8.4f%s (%8.4f)   p = %s\n", r$est, r$stars, r$se, r$p))
cat(sprintf("  R2:               %.4f\n", summary(m3)$r.squared))
cat(sprintf("  Adj R2:           %.4f\n", summary(m3)$adj.r.squared))

# Table 14: Two-Way Geographic Comparison
latex_2way <- two_way_compare %>%
  mutate(latex_row=sprintf("%-30s & %s & %s \\\\",SERVICE_GROUP,county_fmt,city_fmt)) %>%
  pull(latex_row)
cat(paste(c("\\begin{longtable}{lcc}",
            "\\caption{IV Estimates: County vs City Market Definition} \\\\",
            "\\label{tab:geo_compare} \\\\",
            "\\toprule",
            "Service Group & County 9m (Primary) & City 12m (Robustness) \\\\",
            "\\midrule","\\endfirsthead",
            "\\toprule",
            "Service Group & County 9m (Primary) & City 12m (Robustness) \\\\",
            "\\midrule","\\endhead",
            "\\midrule \\multicolumn{3}{r}{\\textit{Continued...}} \\\\ \\endfoot",
            "\\bottomrule",
            "\\multicolumn{3}{l}{\\footnotesize{Two-way clustering: geo-unit + post\\_month.}} \\\\",
            "\\multicolumn{3}{l}{\\footnotesize{*** $p<0.01$, ** $p<0.05$, * $p<0.10$.}} \\\\",
            "\\endlastfoot",
            latex_2way,
            "\\end{longtable}"), collapse="\n"))

# Table 21: LOO shoppability robustness
loo_results_county %>%
  arrange(shop_coef) %>%
  mutate(shop_coef=round(shop_coef,3)) %>%
  kbl(format="latex", booktabs=TRUE,
      caption="Leave-One-System-Out Shoppability Coefficients (County, 9m instrument)",
      label="tab:loo",
      col.names=c("Dropped System","Shoppability Coef.")) %>%
  kable_styling(latex_options="hold_position",font_size=9) %>%
  cat()

# IQR meta-regression (if res_county_outcomes available)
if (exists("res_county_outcomes")) {
  res_county_iqr_only <- res_county_outcomes %>%
    filter(outcome=="IQR") %>%
    mutate(se_pct=se*100, estimate_pct=estimate*100)
  
  if (nrow(res_county_iqr_only) > 5) {
    cat("\n=== IQR OUTCOME META-REGRESSION ===\n")
    extract_coefs <- function(model, model_label, scheme_label) {
      if (is.null(model)) return(NULL)
      ct <- as.data.frame(coef(summary(model)))
      ct$term <- rownames(ct); ct$model <- model_label; ct$scheme <- scheme_label
      ct$r_squared <- summary(model)$r.squared; ct$adj_r_sq <- summary(model)$adj.r.squared
      rownames(ct) <- NULL
      ct %>% select(scheme,model,term,
                    estimate=Estimate,se=`Std. Error`,t_stat=`t value`,p_value=`Pr(>|t|)`,
                    r_squared,adj_r_sq)
    }
    iqr_results <- list()
    for (scheme_nm in names(shoppability_schemes)) {
      scheme  <- shoppability_schemes[[scheme_nm]]
      df_iqr  <- build_meta_data_scheme(res_county_iqr_only, scheme)
      m_iqr_simple <- tryCatch(lm(estimate_pct~is_shoppable,data=df_iqr,weights=weight),
                               error=function(e)NULL)
      m_iqr_three  <- if(any(df_iqr$is_nonshoppable)) {
        tryCatch(lm(estimate_pct~is_shoppable+is_nonshoppable,data=df_iqr,weights=weight),
                 error=function(e)NULL)
      } else NULL
      iqr_results[[paste0(scheme_nm,"_simple")]] <-
        extract_coefs(m_iqr_simple,"Simple shoppable (WLS)",scheme$label)
      iqr_results[[paste0(scheme_nm,"_three")]]  <-
        extract_coefs(m_iqr_three,"Three-way (WLS)",scheme$label)
    }
    iqr_meta_df <- do.call(rbind,iqr_results) %>%
      mutate(sig=ifelse(p_value<0.01,"***",ifelse(p_value<0.05,"**",ifelse(p_value<0.10,"*",""))),
             across(where(is.numeric),~round(.x,4)))
    # ----------------------------------------------------------------------------
    # EXPORT_CSV_iqr_meta_regression
    # ----------------------------------------------------------------------------
    # Exports IQR meta-regression results to iqr_meta_regression.csv.
    write.csv(iqr_meta_df,"iqr_meta_regression.csv",row.names=FALSE)
    cat("Saved to iqr_meta_regression.csv\n")
  }
}

# Final summary tables
# ----------------------------------------------------------------------------
# OUTPUT_SUMMARY_TABLE_01_HET_POOLED
# ----------------------------------------------------------------------------
# Console summary table: pooled IV by heterogeneity dimension.
cat("\n",strrep("=",70),"\n SUMMARY TABLE 1: POOLED IV BY HETEROGENEITY DIMENSION\n",strrep("=",70),"\n")
cat(sprintf("%-25s %-20s %10s %10s %8s\n", "Dimension", "Group", "Pct Effect", "Wald F", "N"))
pooled_summary %>%
  arrange(dimension,group) %>%
  rowwise() %>%
  do({cat(sprintf("%-25s  %-20s  %10s  %10.1f  %8d\n",
                  substr(.$dimension,1,25),substr(.$group,1,20),
                  .$fmt,.$wald_f,.$n_obs)); data.frame()})

cat("\n",strrep("=",70),"\n META-REGRESSION SHOPPABILITY — Theory V2 — All Splits\n",strrep("=",70),"\n")
full_meta_summary %>%
  filter(scheme=="Theory-Based V2 (MRI as nonshoppable)") %>%
  arrange(dimension,split) %>%
  rowwise() %>%
  do({cat(sprintf("%-25s  %-22s  %8.3f  %8.3f  %8.4f  %4s\n",
                  substr(.$dimension,1,25),substr(.$split,1,22),
                  .$estimate,.$se,.$p_value,.$sig)); data.frame()})

cat("\n",strrep("=",70),"\n HOSPITAL TYPE META-REGRESSION\n",strrep("=",70),"\n")
hosp_type_meta_df %>%
  filter(scheme=="Theory-Based V2 (MRI as nonshoppable)") %>%
  rowwise() %>%
  do({cat(sprintf("%-25s  %8.3f  %8.3f  %8.4f  %4s\n",
                  .$hospital_type,.$coef,.$se,.$p_value,.$sig)); data.frame()})

cat("\n",strrep("=",70),"\n GEOGRAPHIC PROXIMITY META-REGRESSION\n",strrep("=",70),"\n")
prox_meta_df %>%
  filter(scheme=="Theory-Based V2 (MRI as nonshoppable)") %>%
  rowwise() %>%
  do({cat(sprintf("%-30s  %8.3f  %8.3f  %8.4f  %4s\n",
                  .$proximity,.$coef,.$se,.$p_value,.$sig)); data.frame()})


# -----------------------------------------------------------
# 1. FIRST STAGE + WU-HAUSMAN
# -----------------------------------------------------------
first_stage_summary <- data.frame(
  stat  = c("instrument", "coefficient", "std_error", "t_stat",
            "wald_f", "wu_hausman_p", "n_obs",
            "n_market_fes", "n_month_fes",
            "adj_r2", "within_r2"),
  value = c(
    "system_peer_pressure_county_9m",
    summary(fs_county_ctrl)$coeftable[1, "Estimate"],
    summary(fs_county_ctrl)$coeftable[1, "Std. Error"],
    summary(fs_county_ctrl)$coeftable[1, "t value"],
    summary(fs_county_ctrl)$coeftable[1, "t value"]^2,
    wh_pval,
    nobs(fs_county_ctrl),
    uniqueN(df_iv_county$market_id),
    uniqueN(df_iv_county$post_month),
    r2(fs_county_ctrl)["ar2"],   # adjusted R²  (fixest function)
    r2(fs_county_ctrl)["wr2"]    # within R²    (fixest function)
  )
)
# ----------------------------------------------------------------------------
# EXPORT_CSV_results_first_stage
# ----------------------------------------------------------------------------
# Exports first-stage summary to results_first_stage.csv.
write.csv(first_stage_summary, "results_first_stage.csv", row.names = FALSE)
cat("Saved: results_first_stage.csv\n")


# -----------------------------------------------------------
# 2. POOLED IV — ALL OUTCOMES
# -----------------------------------------------------------
# ----------------------------------------------------------------------------
# FUNCTION extract_pooled()
# ----------------------------------------------------------------------------
# Converts a pooled model object into a compact export row for the CSV results
#  tables.
extract_pooled <- function(fit, label, is_iv = TRUE) {
  coef_nm <- if (is_iv) "fit_n_prior_posters" else "n_prior_posters"
  est  <- tryCatch(coef(fit)[coef_nm],   error = function(e) NA_real_)
  se_v <- tryCatch(se(fit)[coef_nm],     error = function(e) NA_real_)
  pval <- tryCatch(pvalue(fit)[coef_nm], error = function(e) NA_real_)
  wald <- if (is_iv) tryCatch(fitstat(fit,"ivwald")[[1]]$stat, error=function(e)NA_real_) else NA_real_
  data.frame(outcome = label, estimator = ifelse(is_iv,"IV","OLS"),
             estimate_pct = est * 100, se_pct = se_v * 100,
             pval = pval, stars = add_stars(pval),
             ci_lo_pct = (est - 1.96*se_v) * 100,
             ci_hi_pct = (est + 1.96*se_v) * 100,
             wald_f = wald, n_obs = nobs(fit))
}

pooled_all_outcomes <- bind_rows(
  extract_pooled(ols_county_mean,   "Mean",   is_iv = FALSE),
  extract_pooled(iv_county_mean,    "Mean",   is_iv = TRUE),
  extract_pooled(ols_county_p25,    "P25",    is_iv = FALSE),
  extract_pooled(iv_county_p25,     "P25",    is_iv = TRUE),
  extract_pooled(ols_county_median, "Median", is_iv = FALSE),
  extract_pooled(iv_county_pooled,  "Median", is_iv = TRUE),
  extract_pooled(ols_county_p75,    "P75",    is_iv = FALSE),
  extract_pooled(iv_county_p75,     "P75",    is_iv = TRUE),
  extract_pooled(ols_county_max,    "Max",    is_iv = FALSE),
  extract_pooled(iv_county_max,     "Max",    is_iv = TRUE),
  extract_pooled(ols_county_min,    "Min",    is_iv = FALSE),
  extract_pooled(iv_county_min,     "Min",    is_iv = TRUE),
  extract_pooled(ols_county_iqr,    "IQR",    is_iv = FALSE),
  extract_pooled(iv_county_iqr,     "IQR",    is_iv = TRUE)
) %>% mutate(across(where(is.numeric), ~round(.x, 4)))

# ----------------------------------------------------------------------------
# EXPORT_CSV_results_pooled_iv_all_outcomes
# ----------------------------------------------------------------------------
# Exports pooled IV estimates across all outcomes.
write.csv(pooled_all_outcomes, "results_pooled_iv_all_outcomes.csv", row.names = FALSE)
cat("Saved: results_pooled_iv_all_outcomes.csv\n")


# -----------------------------------------------------------
# 3. SERVICE-LEVEL OLS vs IV COMPARISON (Table 5B)
# -----------------------------------------------------------
# ----------------------------------------------------------------------------
# EXPORT_CSV_results_service_ols_iv_compare
# ----------------------------------------------------------------------------
# Exports service-level OLS-vs-IV comparison.
write.csv(ols_iv_compare_county %>%
            mutate(across(where(is.numeric), ~round(.x, 4))),
          "results_service_ols_iv_compare.csv", row.names = FALSE)
cat("Saved: results_service_ols_iv_compare.csv\n")


# -----------------------------------------------------------
# 4. SERVICE-LEVEL IV ESTIMATES — FULL DETAIL
# (estimate, SE, CI, significance, first-stage F for each service)
# -----------------------------------------------------------
# ----------------------------------------------------------------------------
# EXPORT_CSV_results_service_iv_detail
# ----------------------------------------------------------------------------
# Exports service-level IV estimate details.
write.csv(res_county %>%
            mutate(SERVICE_GROUP = as.character(SERVICE_GROUP),
                   across(where(is.numeric), ~round(.x, 4))),
          "results_service_iv_detail.csv", row.names = FALSE)
cat("Saved: results_service_iv_detail.csv\n")

# ── Palette ──────────────────────────────────────────────────────
COL_GARNET <- "#8B0000"
COL_GOLD   <- "#8B6914"
COL_SLATE  <- "#2C3E50"
COL_GREY   <- "grey65"

cat_cols <- c(
  "Shoppable"     = COL_GARNET,
  "Intermediate"  = COL_GOLD,
  "Non-Shoppable" = COL_SLATE
)

# ── Theory V2 category assignments ───────────────────────────────
shoppable_v2 <- c(
  "CT Lung", "Mammography",
  "Ultrasound Abdomen", "Ultrasound OB", "Ultrasound Pelvis",
  "Ultrasound Vascular", "Ultrasound Extremity",
  "Ultrasound Other", "Ultrasound Breast"
)
nonshoppable_v2 <- c(
  "MRI Brain/Head", "MRI Spine", "MRI Abdomen", "MRI Pelvis",
  "MRI Neck", "MRI Chest", "MRI Extremity", "MRI Angio",
  "MRI Breast", "MRI Other",
  "Biopsy Lymph Node", "Biopsy Liver", "Biopsy Bone",
  "Biopsy Kidney", "Biopsy Lung", "Biopsy Thyroid",
  "Biopsy Pancreas", "Biopsy Other", "Biopsy Breast",
  "Colonoscopy", "Endoscopy"
)

# ── Panel A data ──────────────────────────────────────────────────
svc <- read.csv("results_service_iv_detail.csv") %>%
  mutate(
    iv_pct   = as.numeric(estimate_pct),
    se_pct   = as.numeric(se_pct),
    sig      = as.numeric(pval) < 0.05,
    weight   = 1 / se_pct^2,
    category = case_when(
      SERVICE_GROUP %in% shoppable_v2    ~ "Shoppable",
      SERVICE_GROUP %in% nonshoppable_v2 ~ "Non-Shoppable",
      TRUE                               ~ "Intermediate"
    ),
    category = factor(category,
                      levels = c("Shoppable", "Intermediate", "Non-Shoppable"))
  )

# Precision-weighted category means
cat_means <- svc %>%
  group_by(category) %>%
  summarise(wmean = weighted.mean(iv_pct, weight), .groups = "drop")

# ── Panel A: scatter plot ─────────────────────────────────────────
set.seed(42)

# ----------------------------------------------------------------------------
# OUTPUT_FIGURE_METAREG_PANEL_A_DRAFT
# ----------------------------------------------------------------------------
# Intermediate panel_a draft; replaced by a cleaner panel_a definition below.
panel_a <- ggplot(svc) +
  aes(x = category, y = iv_pct, color = category) +
  # Zero line
  geom_hline(yintercept = 0, color = COL_GREY, linewidth = 0.4,
             linetype = "dashed") +
  # Jittered points — filled = significant, open = not
  geom_jitter(aes(shape = sig, size = sig),
              width = 0.18, height = 0, alpha = 0.80) +
  # Weighted-mean crossbar for each category
  stat_summary(
    fun = function(x) weighted.mean(x, svc$weight[svc$iv_pct %in% x]),
    geom = "crossbar",
    data = left_join(svc, cat_means, by = "category"),
    aes(y = wmean),
    width = 0.45, linewidth = 0.9, fatten = 0, alpha = 0.75,
    show.legend = FALSE
  ) +
  scale_color_manual(values = cat_cols, name = NULL) +
  scale_shape_manual(values = c("FALSE" = 1, "TRUE" = 16), guide = "none") +
  scale_size_manual(values  = c("FALSE" = 2.0, "TRUE" = 2.5), guide = "none") +
  labs(x = NULL,
       y = "IV Estimate (% per additional prior poster)",
       title = "A  Service-Level Estimates by Category",
       subtitle = "Horizontal bar = precision-weighted mean | filled = p < 0.05") +
  theme_minimal(base_size = 10) +
  theme(
    legend.position  = "none",
    axis.text.x      = element_text(size = 10, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.title         = element_text(face = "bold", size = 11),
    plot.subtitle      = element_text(size = 8, color = "grey45")
  )

# ----------------------------------------------------------------------------
# OUTPUT_FIGURE_METAREG_PANEL_A_FINAL
# ----------------------------------------------------------------------------
# Final panel_a: service-level estimates by shoppability category.
panel_a <- ggplot(svc) +
  aes(x = category, y = iv_pct, color = category) +
  geom_hline(yintercept = 0, color = COL_GREY, linewidth = 0.4,
             linetype = "dashed") +
  # Jittered service dots
  geom_jitter(aes(shape = sig, size = sig),
              width = 0.18, height = 0, alpha = 0.80) +
  # Precision-weighted mean bar (drawn from cat_means)
  geom_crossbar(
    data    = cat_means,
    aes(x   = category, y = wmean,
        ymin = wmean, ymax = wmean),
    width   = 0.45, linewidth = 1.0, fatten = 0,
    show.legend = FALSE
  ) +
  scale_color_manual(values = cat_cols, name = NULL) +
  scale_shape_manual(values = c("FALSE" = 1, "TRUE" = 16), guide = "none") +
  scale_size_manual(values  = c("FALSE" = 2.0, "TRUE" = 2.5), guide = "none") +
  labs(x     = NULL,
       y     = "IV Estimate (% per additional prior poster)",
       title = "A  Service-Level Estimates by Category",
       subtitle = "Horizontal bar = precision-weighted mean") +
  theme_minimal(base_size = 10) +
  theme(
    legend.position    = "none",
    axis.text.x        = element_text(size = 10, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.title         = element_text(face = "bold", size = 11),
    plot.subtitle      = element_text(size = 8, color = "grey45")
  )

# ── Panel B: gradient across schemes (unchanged) ──────────────────
scheme_short <- c(
  "CMS Rule Definition"                                    = "CMS Rule",
  "Theory-Based"                                           = "Theory-Based",
  "Broad Shoppable"                                        = "Broad",
  "Imaging vs Procedural"                                  = "Imaging vs Proc.",
  "Theory-Based V2 (MRI as nonshoppable)"                  = "Theory V2 \u25C6",
  "Split Non-Shoppable: MRI vs Procedural"                 = "Split Non-Shop.",
  "CMS Statutory Shoppable List"                           = "CMS Statutory",
  "High vs Low Shoppability Within Modality"               = "High vs Low",
  "CT-Inclusive (all CT modalities shoppable)"             = "CT-Inclusive (All CT)",
  "CT-Inclusive Except Angio (acute vascular CT excluded)" = "CT-Inclusive Ex. Angio",
  # -- 8 alternative frameworks added from Danny's classification-scheme
  # document; keys must match the `label` field of each scheme exactly --
  "Alt: Core Operations Framework (clinical invasiveness)"                    = "Alt: Core Ops.",
  "Alt: CMS Legal/Regulatory Framework (70 vs 230 vs non-listed)"             = "Alt: CMS Legal",
  "Alt: Upfront Cash-Market Framework (MDsave-style voucher availability)"    = "Alt: Cash-Market",
  "Alt: Geographic/Facility Access Framework (retail vs freestanding vs hospital-locked)" = "Alt: Geographic",
  "Alt: Diagnostic Urgency/Lead-Time Framework"                               = "Alt: Urgency",
  "Alt: Staffing/Specialist Framework (generalist vs sub-specialist)"        = "Alt: Staffing",
  "Alt: Incident Reporting/Liability Framework (risk-based)"                 = "Alt: Liability",
  "Alt: No Surprises Act Legal Framework (balance-billing exposure)"         = "Alt: No Surprises"
)

schemes <- read.csv("results_meta_regression_all_schemes.csv") %>%
  filter(model == "Simple (shoppable)") %>%
  mutate(
    estimate   = as.numeric(estimate),
    se         = as.numeric(se),
    p_value    = as.numeric(p_value),
    label      = scheme_short[scheme],
    label      = ifelse(is.na(label), scheme, label),
    is_primary = (scheme == "Theory-Based V2 (MRI as nonshoppable)")
  ) %>%
  arrange(estimate) %>%
  mutate(label = factor(label, levels = label))

# ----------------------------------------------------------------------------
# OUTPUT_FIGURE_METAREG_PANEL_B
# ----------------------------------------------------------------------------
# panel_b: gradient robustness across the ten shoppability classification 
# schemes.
panel_b <- ggplot(schemes) +
  aes(x = estimate, y = label, shape = is_primary) +
  geom_vline(xintercept = 0, color = COL_GREY, linewidth = 0.4,
             linetype = "dashed") +
  geom_errorbarh(
    aes(xmin = estimate - 1.96 * se,
        xmax = estimate + 1.96 * se),
    color = COL_GARNET, height = 0, linewidth = 0.65
  ) +
  geom_point(color = COL_GARNET, size = 3.2) +
  geom_text(
    aes(x = estimate + 1.96 * se + 0.06, label = stars),
    color = COL_GARNET, size = 3.2, hjust = 0, show.legend = FALSE
  ) +
  scale_shape_manual(values = c("TRUE" = 18, "FALSE" = 16), guide = "none") +
  labs(x     = "Shoppability Gradient (pp per additional prior poster)",
       y     = NULL,
       title = paste0("B  Robustness Across ", n_shop_schemes, " Classification Schemes"),
       subtitle = "Diamond = primary spec (Theory V2)  |  whiskers = 95% CI") +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.y        = element_text(size = 9),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.title         = element_text(face = "bold", size = 11),
    plot.subtitle      = element_text(size = 8, color = "grey45")
  )

# ── Combine and save ──────────────────────────────────────────────
# NOTE: panel_b now has n_shop_schemes (18) rows instead of 10. Whatever
# device/dimensions render this combined figure to PDF/PNG elsewhere in your
# pipeline will likely need more vertical room than before, or panel_b's row
# labels will crowd together -- check the rendered output before finalizing.
fig_metareg_combined <- panel_a + panel_b +
  plot_layout(widths = c(1, 1.4)) &
  theme(plot.background = element_rect(fill = "white", color = NA))

# ----------------------------------------------------------------------------
# OUTPUT_FIGURE_02_METAREG_COMBINED_ALT
# ----------------------------------------------------------------------------
# Prints fig_metareg_combined: combined panel A + panel B meta-regression 
# figure.
print(fig_metareg_combined)


# -----------------------------------------------------------
# 5. SHOP2 GROUP ESTIMATES
# -----------------------------------------------------------
# ----------------------------------------------------------------------------
# EXPORT_CSV_results_shop2_estimates
# ----------------------------------------------------------------------------
# Exports SHOP2 group estimates.
write.csv(res_county_shop2 %>%
            mutate(across(where(is.numeric), ~round(.x, 4))),
          "results_shop2_estimates.csv", row.names = FALSE)
cat("Saved: results_shop2_estimates.csv\n")


# -----------------------------------------------------------
# 6. META-REGRESSION: SHOPPABLE COEFFICIENT ACROSS ALL SCHEMES (n_shop_schemes)
# -----------------------------------------------------------
meta_schemes_export <- do.call(rbind, lapply(names(meta_scheme_results), function(nm) {
  res <- meta_scheme_results[[nm]]
  
  # Simple model (shoppable binary)
  ct_simple <- tryCatch(coef(summary(res$simple))["is_shoppableTRUE", ], error = function(e) NULL)
  
  # Three-way model (nonshoppable coefficient)
  ct_three <- if (!is.null(res$three_way)) {
    tryCatch(coef(summary(res$three_way))["is_nonshoppableTRUE", ], error = function(e) NULL)
  } else NULL
  
  row_simple <- if (!is.null(ct_simple)) {
    data.frame(scheme = res$label, model = "Simple (shoppable)", term = "Shoppable",
               estimate = ct_simple["Estimate"], se = ct_simple["Std. Error"],
               p_value = ct_simple["Pr(>|t|)"], stars = add_stars(ct_simple["Pr(>|t|)"]),
               r_squared = summary(res$simple)$r.squared, n_services = nrow(res$data))
  } else NULL
  
  row_three <- if (!is.null(ct_three)) {
    data.frame(scheme = res$label, model = "Three-way (nonshoppable)", term = "Non-Shoppable",
               estimate = ct_three["Estimate"], se = ct_three["Std. Error"],
               p_value = ct_three["Pr(>|t|)"], stars = add_stars(ct_three["Pr(>|t|)"]),
               r_squared = summary(res$three_way)$r.squared, n_services = nrow(res$data))
  } else NULL
  
  bind_rows(row_simple, row_three)
})) %>% mutate(across(where(is.numeric), ~round(.x, 4)))

# ----------------------------------------------------------------------------
# EXPORT_CSV_results_meta_regression_all_schemes
# ----------------------------------------------------------------------------
# Exports all-scheme shoppability meta-regression results.
write.csv(meta_schemes_export, "results_meta_regression_all_schemes.csv", row.names = FALSE)
cat("Saved: results_meta_regression_all_schemes.csv\n")


# -----------------------------------------------------------
# 7. SATURATION — 3-BIN POOLED RESULTS
# -----------------------------------------------------------
# ----------------------------------------------------------------------------
# EXPORT_CSV_results_saturation_3bin_pooled
# ----------------------------------------------------------------------------
# Exports pooled saturation-bin results.
write.csv(results_3bin_county %>%
            mutate(across(where(is.numeric), ~round(.x, 4))),
          "results_saturation_3bin_pooled.csv", row.names = FALSE)
cat("Saved: results_saturation_3bin_pooled.csv\n")


# -----------------------------------------------------------
# 8. SATURATION — BY SERVICE GROUP
# -----------------------------------------------------------
# ----------------------------------------------------------------------------
# EXPORT_CSV_results_saturation_by_service
# ----------------------------------------------------------------------------
# Exports service-level saturation results.
write.csv(sat_service_all %>%
            mutate(SERVICE_GROUP = as.character(SERVICE_GROUP),
                   across(where(is.numeric), ~round(.x, 4))),
          "results_saturation_by_service.csv", row.names = FALSE)
cat("Saved: results_saturation_by_service.csv\n")


# -----------------------------------------------------------
# 9. SATURATION — SHOP2 BY BIN
# -----------------------------------------------------------
# ----------------------------------------------------------------------------
# EXPORT_CSV_results_saturation_shop2
# ----------------------------------------------------------------------------
# Exports SHOP2 saturation results.
write.csv(sat_shop2_all %>%
            mutate(across(where(is.numeric), ~round(.x, 4))),
          "results_saturation_shop2.csv", row.names = FALSE)
cat("Saved: results_saturation_shop2.csv\n")


# -----------------------------------------------------------
# 10. SATURATION — META-REGRESSION ACROSS ALL SCHEMES BY BIN (n_shop_schemes)
# -----------------------------------------------------------
# ----------------------------------------------------------------------------
# EXPORT_CSV_results_saturation_meta_schemes
# ----------------------------------------------------------------------------
# Exports saturation meta-regressions by classification scheme.
write.csv(sat_meta_df %>%
            filter(term != "(Intercept)") %>%
            mutate(across(where(is.numeric), ~round(.x, 4))),
          "results_saturation_meta_schemes.csv", row.names = FALSE)
cat("Saved: results_saturation_meta_schemes.csv\n")

# ================================================================
# fig:saturation — Shoppable and Non-Shoppable Effects by Bin
# Theory V2 primary scheme
# Bins 2 & 3 visually distinguished as non-causally-identified
# ================================================================

library(ggplot2)
library(dplyr)

TARGET_SCHEME <- "Theory-Based V2 (MRI as nonshoppable)"

sat <- read.csv("results_saturation_meta_schemes.csv") %>%
  filter(scheme == TARGET_SCHEME) %>%
  mutate(
    estimate = as.numeric(estimate),
    se       = as.numeric(se),
    p_value  = as.numeric(p_value),
    ci_lo    = estimate - 1.96 * se,
    ci_hi    = estimate + 1.96 * se,
    
    # Numeric x-position so annotate("rect") works cleanly
    x_pos = case_when(
      bin == "1-3 prior posters" ~ 1,
      bin == "4-8 prior posters" ~ 2,
      bin == "9+ prior posters"  ~ 3,
      TRUE ~ NA_real_
    ),
    
    identified = bin == "1-3 prior posters",
    
    series = case_when(
      model == "Simple shoppable (WLS)" & term == "is_shoppableTRUE"    ~ "Shoppable",
      model == "Three-way (WLS)"        & term == "is_nonshoppableTRUE" ~ "Non-Shoppable",
      TRUE ~ NA_character_
    ),
    
    stars = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ ""
    )
  ) %>%
  filter(!is.na(series), !is.na(x_pos))

# Custom x-axis labels
x_labels <- c(
  "1" = "Bin 1\n1–3 posters\nWald F = 13.0\n[identified]",
  "2" = "Bin 2\n4–8 posters\nWald F = 7.2\n[weak IV]",
  "3" = "Bin 3\n9+ posters\nWald F = 0.006\n[no variation]"
)

COL_GARNET <- "#8B0000"
COL_SLATE  <- "#2C3E50"

# ----------------------------------------------------------------------------
# OUTPUT_FIGURE_04_SATURATION_MAIN
# ----------------------------------------------------------------------------
# Figure object fig_saturation: shoppable and non-shoppable effects by 
# saturation bin.
fig_saturation <- ggplot(sat) +
  
  aes(
    x     = x_pos,
    y     = estimate,
    color = series,
    group = series
  ) +
  
  # Grey shading for non-identified / descriptive region
  annotate(
    "rect",
    xmin = 1.5, xmax = 3.5,
    ymin = -Inf, ymax = Inf,
    fill = "grey92", alpha = 0.65
  ) +
  
  annotate(
    "text",
    x = 2.5, y = Inf,
    label = " ",
    vjust = 1.6,
    size = 3.0,
    color = "grey45",
    fontface = "italic"
  ) +
  
  # Zero reference line
  geom_hline(
    yintercept = 0,
    color = "grey65",
    linewidth = 0.4,
    linetype = "dashed"
  ) +
  
  # Full-sample gradient reference
  geom_hline(
    yintercept = -2.603,
    color = COL_GARNET,
    linewidth = 0.35,
    linetype = "dotted",
    alpha = 0.5
  ) +
  
  annotate(
    "text",
    x = 0.62,
    y = -2.603,
    label = "Full-sample gradient (−2.60 pp)",
    hjust = 0,
    vjust = -0.45,
    size = 2.8,
    color = COL_GARNET,
    fontface = "italic"
  ) +
  
  # Confidence interval ribbon
  geom_ribbon(
    aes(
      ymin = ci_lo,
      ymax = ci_hi,
      fill = series
    ),
    alpha = 0.12,
    color = NA,
    show.legend = FALSE
  ) +
  
  # Descriptive connecting line.
  # Keep dashed throughout to avoid implying Bin 2 is causally identified.
  geom_line(
    linewidth = 0.9,
    linetype = "dashed"
  ) +
  
  # Points: filled = identified Bin 1; open = descriptive Bins 2–3
  geom_point(
    aes(shape = identified),
    size = 3.8,
    stroke = 1.0
  ) +
  
  # Significance stars
  geom_text(
    data = filter(sat, stars != ""),
    aes(
      label = stars,
      y = ci_hi + 1.2
    ),
    size = 4.5,
    show.legend = FALSE,
    fontface = "bold"
  ) +
  
  scale_color_manual(
    values = c(
      "Shoppable"     = COL_GARNET,
      "Non-Shoppable" = COL_SLATE
    ),
    name = NULL
  ) +
  
  scale_fill_manual(
    values = c(
      "Shoppable"     = COL_GARNET,
      "Non-Shoppable" = COL_SLATE
    )
  ) +
  
  scale_shape_manual(
    values = c(
      "TRUE"  = 16,
      "FALSE" = 1
    ),
    guide = "none"
  ) +
  
  scale_x_continuous(
    breaks = 1:3,
    labels = x_labels,
    limits = c(0.6, 3.4)
  ) +
  
  labs(
    x = NULL,
    y = "Meta-Regression Coefficient\n(pp per additional prior poster)",
    title = "Shoppable and Non-Shoppable Price Effects by Market Saturation",
    caption = paste0(
      "Theory V2 scheme. Filled point = Bin 1 (Wald F = 13.0, marginal identification).\n",
      "Open points and shaded region = Bins 2–3, descriptive only. Shaded bands = 95% CIs."
    )
  ) +
  
  theme_minimal(base_size = 11) +
  
  theme(
    legend.position = c(0.15, 0.15),
    legend.background = element_rect(
      fill = "white",
      color = "grey85",
      linewidth = 0.4
    ),
    legend.text = element_text(size = 10),
    axis.text.x = element_text(
      size = 8.5,
      lineheight = 1.25
    ),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.caption = element_text(
      size = 7.5,
      color = "grey40",
      hjust = 0,
      lineheight = 1.3
    ),
    plot.title = element_text(
      face = "bold",
      size = 12
    )
  )

# ----------------------------------------------------------------------------
# OUTPUT_FIGURE_04_SATURATION_MAIN_PRINT
# ----------------------------------------------------------------------------
# Prints fig_saturation.
print(fig_saturation)


# -----------------------------------------------------------
# 11. INSTRUMENT VARIANTS — COUNTY
# -----------------------------------------------------------
# ----------------------------------------------------------------------------
# EXPORT_CSV_results_instrument_variants_county
# ----------------------------------------------------------------------------
# Exports county instrument-variant robustness results.
write.csv(results_variants_county %>%
            mutate(across(where(is.numeric), ~round(.x, 4))),
          "results_instrument_variants_county.csv", row.names = FALSE)
cat("Saved: results_instrument_variants_county.csv\n")


# -----------------------------------------------------------
# 12. INSTRUMENT VARIANTS — CITY
# -----------------------------------------------------------
city_variants_export <- map_dfr(names(variants_city), function(nm) {
  instr <- variants_city[[nm]]
  if (!instr %in% names(df_iv_city)) return(NULL)
  fit <- tryCatch(
    feols(as.formula(paste0(
      "ln_median_price ~ ln_total_beds | market_id + post_month | n_prior_posters ~ ", instr)),
      data = df_iv_city, cluster = ~city_state + post_month),
    error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  wald <- tryCatch(fitstat(fit, "ivwald")[[1]]$stat, error = function(e) NA_real_)
  tibble(variant = nm, instrument = instr,
         estimate_pct = coef(fit)["fit_n_prior_posters"] * 100,
         se_pct       = se(fit)["fit_n_prior_posters"]   * 100,
         pval         = pvalue(fit)["fit_n_prior_posters"],
         stars        = add_stars(pvalue(fit)["fit_n_prior_posters"]),
         wald_f       = wald, n_obs = nobs(fit))
}) %>% mutate(across(where(is.numeric), ~round(.x, 4)))

# ----------------------------------------------------------------------------
# EXPORT_CSV_results_instrument_variants_city
# ----------------------------------------------------------------------------
# Exports city instrument-variant robustness results.
write.csv(city_variants_export, "results_instrument_variants_city.csv", row.names = FALSE)
cat("Saved: results_instrument_variants_city.csv\n")


# -----------------------------------------------------------
# 13. LEAVE-ONE-SYSTEM-OUT
# -----------------------------------------------------------
loo_export <- loo_results_county %>%
  arrange(shop_coef) %>%
  mutate(
    shop_coef  = round(shop_coef, 4),
    vs_full    = round(shop_coef - mean(shop_coef, na.rm = TRUE), 4),
    full_sample_gradient = round(coef(summary(meta_county$simple))["is_shoppableTRUE","Estimate"], 4)
  )
# ----------------------------------------------------------------------------
# EXPORT_CSV_results_loo
# ----------------------------------------------------------------------------
# Exports leave-one-system-out robustness results.
write.csv(loo_export, "results_loo.csv", row.names = FALSE)
cat("Saved: results_loo.csv\n")


# -----------------------------------------------------------
# 14. PRE-TREND LEAD TESTS
# -----------------------------------------------------------
# ----------------------------------------------------------------------------
# FUNCTION extract_lead()
# ----------------------------------------------------------------------------
# Converts one placebo lead model into a compact export row.
extract_lead <- function(fit, lead_name) {
  ct <- tryCatch(coeftable(fit)[lead_name, ], error = function(e) NULL)
  if (is.null(ct)) return(NULL)
  data.frame(lead = lead_name,
             estimate = ct["Estimate"], se = ct["Std. Error"],
             t_stat   = ct["t value"],  p_value = ct["Pr(>|t|)"],
             stars    = add_stars(ct["Pr(>|t|)"]), n_obs = nobs(fit))
}

leads_export <- bind_rows(
  extract_lead(pt_lead1, "spp_lead1"),
  extract_lead(pt_lead2, "spp_lead2"),
  extract_lead(pt_lead3, "spp_lead3")
) %>%
  bind_rows(
    data.frame(lead = "joint_F",
               estimate = wald_joint_leads$stat,
               se = NA, t_stat = NA,
               p_value = wald_joint_leads$p,
               stars = add_stars(wald_joint_leads$p),
               n_obs = nobs(pt_joint))
  ) %>%
  mutate(across(where(is.numeric), ~round(.x, 4)))

# ----------------------------------------------------------------------------
# EXPORT_CSV_results_pre_trend_leads
# ----------------------------------------------------------------------------
# Exports pre-trend/instrument-lead results.
write.csv(leads_export, "results_pre_trend_leads.csv", row.names = FALSE)
cat("Saved: results_pre_trend_leads.csv\n")


# -----------------------------------------------------------
# 15. RANDOMIZATION INFERENCE
# -----------------------------------------------------------
ri_export <- data.frame(
  stat  = c("actual_first_stage_F", "permutation_95th_pctile",
            "ri_p_value", "n_permutations", "instrument", "clustering"),
  value = c(round(actual_f, 2), round(p95_perm, 2),
            round(ri_pvalue, 4), length(perm_f_clean),
            "system_peer_pressure_county_9m",
            "County_State + post_month")
)
# ----------------------------------------------------------------------------
# EXPORT_CSV_results_randomization_inference
# ----------------------------------------------------------------------------
# Exports randomization-inference summary.
write.csv(ri_export, "results_randomization_inference.csv", row.names = FALSE)
cat("Saved: results_randomization_inference.csv\n")

# Also save the full permutation distribution
# ----------------------------------------------------------------------------
# EXPORT_CSV_results_ri_permutation_distribution
# ----------------------------------------------------------------------------
# Exports permutation F-stat distribution.
write.csv(data.frame(perm_f = round(perm_f_clean, 4)),
          "results_ri_permutation_distribution.csv", row.names = FALSE)
cat("Saved: results_ri_permutation_distribution.csv\n")


# -----------------------------------------------------------
# 16. BLACK x MEDICAID / POVERTY 2x2 RESULTS
# -----------------------------------------------------------
bm_pooled_export <- do.call(rbind, lapply(bm_results, function(r) r$pooled)) %>%
  mutate(stars = add_stars(pvalue),
         across(where(is.numeric), ~round(.x, 4)))
# ----------------------------------------------------------------------------
# EXPORT_CSV_results_black_medicaid_2x2_pooled
# ----------------------------------------------------------------------------
# Exports Black share x Medicaid 2x2 pooled estimates.
write.csv(bm_pooled_export, "results_black_medicaid_2x2_pooled.csv", row.names = FALSE)
cat("Saved: results_black_medicaid_2x2_pooled.csv\n")

bm_meta_export <- bm_meta_summary %>%
  mutate(across(where(is.numeric), ~round(.x, 4)))
# ----------------------------------------------------------------------------
# EXPORT_CSV_results_black_medicaid_2x2_meta
# ----------------------------------------------------------------------------
# Exports Black share x Medicaid 2x2 meta-regression estimates.
write.csv(bm_meta_export, "results_black_medicaid_2x2_meta.csv", row.names = FALSE)
cat("Saved: results_black_medicaid_2x2_meta.csv\n")



# -----------------------------------------------------------
# 18. CLUSTERING ROBUSTNESS (two-way vs one-way)
# -----------------------------------------------------------
clustering_export <- bind_rows(
  extract_pooled(iv_county_pooled, "Median (two-way)", is_iv = TRUE) %>%
    mutate(clustering = "County_State + post_month"),
  extract_pooled(iv_county_1way,   "Median (one-way)", is_iv = TRUE) %>%
    mutate(clustering = "County_State only")
) %>% mutate(across(where(is.numeric), ~round(.x, 4)))
# ----------------------------------------------------------------------------
# EXPORT_CSV_results_clustering_robustness
# ----------------------------------------------------------------------------
# Exports clustering robustness comparison.
write.csv(clustering_export, "results_clustering_robustness.csv", row.names = FALSE)
cat("Saved: results_clustering_robustness.csv\n")


# -----------------------------------------------------------
# 19. GEOGRAPHIC COMPARISON: COUNTY vs CITY (all 46 services)
# -----------------------------------------------------------
# ----------------------------------------------------------------------------
# EXPORT_CSV_results_geographic_comparison_county_city
# ----------------------------------------------------------------------------
# Exports county-vs-city geographic comparison.
write.csv(two_way_compare %>%
            mutate(across(where(is.numeric), ~round(.x, 4))),
          "results_geographic_comparison_county_city.csv", row.names = FALSE)
cat("Saved: results_geographic_comparison_county_city.csv\n")


# -----------------------------------------------------------
# 20. PRINT SUMMARY: ALL FILES SAVED
# -----------------------------------------------------------
saved_files <- c(
  "results_first_stage.csv",
  "results_pooled_iv_all_outcomes.csv",
  "results_service_ols_iv_compare.csv",
  "results_service_iv_detail.csv",
  "results_shop2_estimates.csv",
  "results_meta_regression_all_schemes.csv",
  "results_saturation_3bin_pooled.csv",
  "results_saturation_by_service.csv",
  "results_saturation_shop2.csv",
  "results_saturation_meta_schemes.csv",
  "results_instrument_variants_county.csv",
  "results_instrument_variants_city.csv",
  "results_loo.csv",
  "results_pre_trend_leads.csv",
  "results_randomization_inference.csv",
  "results_ri_permutation_distribution.csv",
  "results_black_medicaid_2x2_pooled.csv",
  "results_black_medicaid_2x2_meta.csv",
  "results_causal_forest_ate.csv",
  "results_causal_forest_variable_importance.csv",
  "results_clustering_robustness.csv",
  "results_geographic_comparison_county_city.csv",
  # Previously saved by main script:
  "het_pooled_summary.csv",
  "het_shop2_summary.csv",
  "het_meta_summary.csv",
  "het_meta_full_matrix.csv",
  "hosp_type_meta.csv",
  "proximity_meta.csv",
  "county_outcomes_detail.csv",
  "iqr_meta_regression.csv"
)

cat("\n", strrep("=", 60), "\n")
cat(" ALL RESULTS FILES\n")
cat(strrep("=", 60), "\n")
for (f in saved_files) cat(sprintf("  %s\n", f))
cat(strrep("=", 60), "\n")
cat(sprintf("Total files: %d\n", length(saved_files)))


######## Section 21: ML Extension — Residualized Instrumental Forest ###########
# ----------------------------------------------------------------------------
# Fixed version: keeps X, Y, W, and Z aligned after FE residualization

library(dplyr)
library(fixest)
library(grf)
library(ggplot2)

# ------------------------------------------------------------
# 1. Define heterogeneity variables
# ------------------------------------------------------------

forest_vars <- c(
  "ln_total_beds",
  "Share_ForProfit",
  "Share_NonProfit",
  "uninsured_rate",
  "medicaid_share",
  "poverty_rate",
  "log_median_income",
  "age65plus_share",
  "college_share",
  "black_share",
  "hispanic_share",
  "employment_rate",
  "homeowner_rate",
  "population"
)

forest_var_cols <- intersect(forest_vars, names(df_iv_county))

# ------------------------------------------------------------
# 2. Build clean sample
# ------------------------------------------------------------

keep_forest <- complete.cases(df_iv_county[, forest_var_cols]) &
  !is.na(df_iv_county$ln_median_price) &
  !is.na(df_iv_county$n_prior_posters) &
  !is.na(df_iv_county$system_peer_pressure_county_9m) &
  !is.na(df_iv_county$market_id) &
  !is.na(df_iv_county$post_month)

df_forest <- df_iv_county[keep_forest, ] %>%
  mutate(row_id_forest = row_number())

cat("\nInitial causal forest sample size:", nrow(df_forest), "\n")
cat("Number of forest covariates:", length(forest_var_cols), "\n")

# ------------------------------------------------------------
# 3. Residualize Y, W, and Z using the preferred fixed effects
# ------------------------------------------------------------
# Important:
#   resid(..., na.rm = FALSE) keeps residuals aligned with the original
#   df_forest rows by returning NA for observations dropped by feols.

resid_y_fit <- feols(
  ln_median_price ~ 1 | market_id + post_month,
  data = df_forest
)

resid_w_fit <- feols(
  n_prior_posters ~ 1 | market_id + post_month,
  data = df_forest
)

resid_z_fit <- feols(
  system_peer_pressure_county_9m ~ 1 | market_id + post_month,
  data = df_forest
)

df_forest <- df_forest %>%
  mutate(
    Y_resid = as.numeric(resid(resid_y_fit, na.rm = FALSE)),
    W_resid = as.numeric(resid(resid_w_fit, na.rm = FALSE)),
    Z_resid = as.numeric(resid(resid_z_fit, na.rm = FALSE))
  )

# ------------------------------------------------------------
# 4. Keep only rows with valid residuals
# ------------------------------------------------------------

df_forest_resid <- df_forest %>%
  filter(
    is.finite(Y_resid),
    is.finite(W_resid),
    is.finite(Z_resid)
  )

cat("Residualized causal forest sample size:", nrow(df_forest_resid), "\n")

# ------------------------------------------------------------
# 5. Prepare aligned X, Y, W, and Z
# ------------------------------------------------------------

X_forest <- df_forest_resid %>%
  select(all_of(forest_var_cols)) %>%
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()

Y_resid <- as.numeric(df_forest_resid$Y_resid)
W_resid <- as.numeric(df_forest_resid$W_resid)
Z_resid <- as.numeric(df_forest_resid$Z_resid)

# Remove zero-variance covariates
x_sd <- apply(X_forest, 2, sd, na.rm = TRUE)
X_forest <- X_forest[, x_sd > 0, drop = FALSE]

# Check alignment before running forest
cat("\nAlignment check:\n")
cat("  nrow(X_forest):", nrow(X_forest), "\n")
cat("  length(Y_resid):", length(Y_resid), "\n")
cat("  length(W_resid):", length(W_resid), "\n")
cat("  length(Z_resid):", length(Z_resid), "\n")

stopifnot(
  nrow(X_forest) == length(Y_resid),
  nrow(X_forest) == length(W_resid),
  nrow(X_forest) == length(Z_resid)
)

# ------------------------------------------------------------
# 6. Estimate residualized instrumental forest
# ------------------------------------------------------------

set.seed(42)

iv_forest_resid <- instrumental_forest(
  X = X_forest,
  Y = Y_resid,
  W = W_resid,
  Z = Z_resid,
  num.trees = 4000,
  min.node.size = 10,
  honesty = TRUE,
  seed = 42
)

# ------------------------------------------------------------
# 7. Predicted heterogeneous IV effects
# ------------------------------------------------------------
# Because Z is a count-valued/continuous instrument, grf cannot compute
# average_treatment_effect() for this instrumental forest.
#
# Instead, we use the forest's predicted conditional IV effects tau_hat(x)
# as an exploratory heterogeneity measure.

tau_hat <- predict(iv_forest_resid)$predictions

df_forest_results <- df_forest_resid %>%
  mutate(
    tau_hat     = tau_hat,
    tau_hat_pct = tau_hat * 100
  )

cat("\n=== Residualized Instrumental Forest: Predicted Effects ===\n")
cat("Note: average_treatment_effect() is not used because the instrument is count-valued, not binary.\n")
cat("These are predicted conditional IV effects from the forest.\n\n")

summary(df_forest_results$tau_hat_pct)

cat(sprintf("\nMean predicted effect   : %.3f percent\n",
            mean(df_forest_results$tau_hat_pct, na.rm = TRUE)))
cat(sprintf("Median predicted effect : %.3f percent\n",
            median(df_forest_results$tau_hat_pct, na.rm = TRUE)))
cat(sprintf("10th percentile         : %.3f percent\n",
            quantile(df_forest_results$tau_hat_pct, 0.10, na.rm = TRUE)))
cat(sprintf("90th percentile         : %.3f percent\n",
            quantile(df_forest_results$tau_hat_pct, 0.90, na.rm = TRUE)))

# ------------------------------------------------------------
# 8. Variable importance
# ------------------------------------------------------------

var_imp_df <- data.frame(
  variable   = colnames(X_forest),
  importance = as.numeric(variable_importance(iv_forest_resid))
) %>%
  arrange(desc(importance))

cat("\n=== Variable Importance: Top 10 ===\n")
print(head(var_imp_df, 10), row.names = FALSE)

# ------------------------------------------------------------
# 9. Variable importance
# ------------------------------------------------------------
# Variable importance tells us which covariates the forest used most
# often to detect treatment-effect heterogeneity.

var_imp_df <- data.frame(
  variable   = colnames(X_forest),
  importance = as.numeric(variable_importance(iv_forest_resid))
) %>%
  arrange(desc(importance))

cat("\n=== Variable Importance: Top 10 ===\n")
print(head(var_imp_df, 10), row.names = FALSE)

var_imp_top10 <- var_imp_df %>%
  slice_max(importance, n = 10)

# ------------------------------------------------------------
# 10. Plot variable importance
# ------------------------------------------------------------

fig_cf_varimp <- ggplot(
  var_imp_top10,
  aes(x = importance, y = reorder(variable, importance))
) +
  geom_col(fill = FSU_GARNET, alpha = 0.85) +
  labs(
    title = "Instrumental Forest: Variable Importance",
    subtitle = "Exploratory heterogeneity after residualizing market and month fixed effects",
    x = "Variable importance",
    y = NULL,
    caption = paste0(
      "Notes: Outcome, treatment, and instrument are residualized by market_id and post_month ",
      "before estimating the instrumental forest. Variable importance indicates which covariates ",
      "the forest uses most often to detect heterogeneity; it is not a causal effect of the covariate."
    )
  ) +
  theme_bw(base_size = 11)

print(fig_cf_varimp)

# ------------------------------------------------------------
# 11. Plot distribution of predicted treatment effects
# ------------------------------------------------------------
# tau_hat_pct is the forest-predicted conditional IV effect in percent.

fig_cf_tau_dist <- ggplot(df_forest_results, aes(x = tau_hat_pct)) +
  geom_histogram(bins = 40, fill = FSU_GARNET, alpha = 0.85, color = "white") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  labs(
    title = "Distribution of Instrumental-Forest Predicted Effects",
    subtitle = "Estimated heterogeneous effects of prior posters on median prices",
    x = "Predicted IV effect, percent",
    y = "Number of observations",
    caption = paste0(
      "Notes: Effects are predicted conditional IV effects from a residualized instrumental forest. ",
      "They are exploratory and are not used as baseline causal estimates."
    )
  ) +
  theme_bw(base_size = 11)

print(fig_cf_tau_dist)

# ------------------------------------------------------------
# 12. Forest-implied heterogeneity summaries
# ------------------------------------------------------------
# These summaries compare predicted effects across high/low market
# characteristics. They are descriptive diagnostics, not formal subgroup
# IV estimates.

summarize_tau_split <- function(data, split_var, label) {
  
  cutoff <- median(data[[split_var]], na.rm = TRUE)
  
  data %>%
    mutate(
      group = ifelse(.data[[split_var]] >= cutoff, "High", "Low")
    ) %>%
    group_by(group) %>%
    summarise(
      moderator = label,
      cutoff = cutoff,
      n = n(),
      mean_tau_pct = mean(tau_hat_pct, na.rm = TRUE),
      median_tau_pct = median(tau_hat_pct, na.rm = TRUE),
      p10_tau_pct = quantile(tau_hat_pct, 0.10, na.rm = TRUE),
      p90_tau_pct = quantile(tau_hat_pct, 0.90, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    select(moderator, group, cutoff, n, mean_tau_pct, median_tau_pct,
           p10_tau_pct, p90_tau_pct)
}

heterogeneity_summary_long <- bind_rows(
  summarize_tau_split(df_forest_results, "college_share",  "College share"),
  summarize_tau_split(df_forest_results, "poverty_rate",   "Poverty rate"),
  summarize_tau_split(df_forest_results, "medicaid_share", "Medicaid share"),
  summarize_tau_split(df_forest_results, "black_share",    "Black share"),
  summarize_tau_split(df_forest_results, "population",     "Population")
) %>%
  mutate(across(where(is.numeric), ~round(.x, 3)))

cat("\n=== Forest-Implied Heterogeneity Summary ===\n")
print(heterogeneity_summary_long, row.names = FALSE)

# ------------------------------------------------------------
# 13. High-minus-low differences
# ------------------------------------------------------------

heterogeneity_diff <- heterogeneity_summary_long %>%
  select(moderator, group, mean_tau_pct) %>%
  tidyr::pivot_wider(
    names_from = group,
    values_from = mean_tau_pct
  ) %>%
  mutate(
    high_minus_low = High - Low
  ) %>%
  arrange(high_minus_low)

cat("\n=== Forest-Implied High-minus-Low Differences ===\n")
print(heterogeneity_diff, row.names = FALSE)


# ============================================================
# SCRIPT SUMMARY — PRIMARY AND ROBUSTNESS SPECIFICATIONS
# (Two more sections follow this summary: within-system price
# dispersion, and the payer-conditional robustness pipeline.)
# ============================================================
# Primary specification:   County + system_peer_pressure_county_9m
#   Wald F = 27.1 (two-way: County_State + post_month)
# Robustness check 1:      City  + system_peer_pressure_city_12m
#   Wald F = 26.7 (two-way: city_state + post_month)
# Robustness check 2:      County 6m  (F = 9.6,  window sensitivity)
# Robustness check 3:      County 12m (F = 26.1, window sensitivity)
# Robustness check 4:      City 9m    (F = 20.0, window sensitivity)
# Robustness check 5:      City 6m    (F = 18.9, window sensitivity)
# Robustness check 6:      County lag instrument (6-12m)
# Robustness check 7:      County no large systems
# Robustness check 8:      Leave-one-system-out (Section 16)
# Robustness check 9:      One-way vs two-way clustering (Section 6)
# Placebo tests:           Pre-trend leads (Section 17)
#                          Randomization inference (Section 17)
# ML extensions:           PDS-LASSO (Section 19)
#                          Causal Forest (Section 21)
# CBSA dropped:            All variants fail F < 10 under two-way clustering
# ============================================================

######## Section 22: Within-System Price Dispersion Check ###########
# ----------------------------------------------------------------------------
# Purpose: Test whether hospitals in the same system price
# uniformly across counties. Low dispersion would support the
# centralized contracting concern raised by the referee;
# high dispersion supports local contracting and the exclusion
# restriction of the system peer pressure instrument.

# Step 1: Identify multi-county systems
# Only systems with hospitals in more than one county can
# generate meaningful within-system cross-county dispersion.
# Drop unaffiliated hospitals (missing or blank system name).
multi_county_systems <- df_iv_county %>%
  filter(!is.na(HEALTH_SYSTEM_NAME), HEALTH_SYSTEM_NAME != "") %>%
  group_by(HEALTH_SYSTEM_NAME) %>%
  summarise(n_counties = n_distinct(County_State)) %>%
  filter(n_counties > 1) %>%
  pull(HEALTH_SYSTEM_NAME)

# Step 2: Compute within-system price dispersion by service group
# Each hospital contributes one observation per service group:
# its ln_median_price (log of the median negotiated price across
# payers and CPT codes within that service group). The SD of
# log prices across hospitals within a system-service group cell
# approximates the coefficient of variation and is interpretable
# as approximate percent price variation. Filter to cells with
# at least 2 hospitals since SD is undefined for n = 1.
dispersion <- df_iv_county %>%
  filter(HEALTH_SYSTEM_NAME %in% multi_county_systems) %>%
  group_by(HEALTH_SYSTEM_NAME, SERVICE_GROUP) %>%
  summarise(
    n_hospitals      = n(),
    mean_lnprice     = mean(ln_median_price, na.rm = TRUE),
    sd_across_hosps  = sd(ln_median_price, na.rm = TRUE),
    .groups          = "drop"
  ) %>%
  filter(n_hospitals > 1)

# Step 3: Summarize dispersion across systems by service group
# Reports the median and mean within-system SD of log prices
# for each service group. The median is preferred as the
# summary statistic since the SD distribution is right-skewed.
# An SD of ~0.20 implies roughly 20% price variation within
# the same system across counties — inconsistent with
# centralized uniform pricing.
dispersion %>%
  group_by(SERVICE_GROUP) %>%
  summarise(
    median_sd = median(sd_across_hosps, na.rm = TRUE),
    mean_sd   = mean(sd_across_hosps, na.rm = TRUE)
  ) %>%
  as.data.frame() %>% print()










######## Section 23: Payer-Conditional Robustness Check (Full Pipeline) ###########
# ----------------------------------------------------------------------------
# Excludes "CT Other" -- a residual/catch-all CT category that was
# found to be the single largest driver of the Medicaid and
# MedicareAdv shoppability-gradient sign patterns. Excluded here
# because it aggregates heterogeneous, ambiguously-classified CT
# codes rather than representing a coherent clinical service --
# document this exclusion explicitly in the paper/appendix.
#
# Reuses run_iv_pooled(), run_iv_by_service(), build_meta_data(),
# run_meta_regressions(), print_meta_results() exactly as defined
# earlier in this script -- run this section after sourcing the rest
# of the script (Sections 0-22).

# ---------------------------------------------------------------------
# 0. Load payer-split price panel (from payer_split_prices.sql)
# ---------------------------------------------------------------------
payer_split_raw <- fread("../Data/data_clean/payer_split_prices.csv")

payer_split_raw <- payer_split_raw %>%
  rename(County_State = COUNTY_STATE) %>%             # Snowflake export uppercased this column
  mutate(
    post_month      = as.Date(INGESTED_ON),
    post_month      = lubridate::floor_date(post_month, "month"),
    ln_median_price = log(MEDIAN_PRICE + 1),
    ln_mean_price   = log(MEAN_PRICE   + 1),
    ln_p25_price    = log(P25          + 1),
    ln_p75_price    = log(P75          + 1),
    ln_max_price    = log(MAX_PRICE    + 1),
    ln_min_price    = log(MIN_PRICE    + 1),
    ln_iqr_price    = log(IQR_PRICE    + 1)
  )

payer_buckets <- c("Commercial", "Medicaid", "MedicareAdv")
price_cols <- c("ln_median_price", "ln_mean_price", "ln_p25_price",
                "ln_p75_price", "ln_max_price", "ln_min_price", "ln_iqr_price")

# ---------------------------------------------------------------------
# 1. Build payer-specific panels, EXCLUDING CT Other
# ---------------------------------------------------------------------
build_payer_panel <- function(pb) {
  price_sub <- payer_split_raw %>%
    filter(PAYER_BUCKET == pb, SERVICE_GROUP != "CT Other") %>%
    select(HOSPITAL_ID, County_State, post_month, SERVICE_GROUP, all_of(price_cols))
  
  df_iv_county %>%
    filter(SERVICE_GROUP != "CT Other") %>%
    select(-all_of(price_cols)) %>%
    inner_join(
      price_sub,
      by = c("HOSPITAL_ID", "County_State", "post_month", "SERVICE_GROUP")
    )
}

df_payer_panels <- setNames(lapply(payer_buckets, build_payer_panel), payer_buckets)

for (pb in payer_buckets) {
  cat(sprintf("\n%s panel (CT Other excluded): %d rows | %d hospitals | %d counties | %d services\n",
              pb, nrow(df_payer_panels[[pb]]),
              uniqueN(df_payer_panels[[pb]]$HOSPITAL_ID),
              uniqueN(df_payer_panels[[pb]]$County_State),
              uniqueN(df_payer_panels[[pb]]$SERVICE_GROUP)))
}

# ---------------------------------------------------------------------
# 2. Main pooled IV coefficient, per payer bucket
# ---------------------------------------------------------------------
iv_pooled_by_payer <- lapply(payer_buckets, function(pb) {
  run_iv_pooled(df_payer_panels[[pb]], "ln_median_price",
                "system_peer_pressure_county_9m",
                c("County_State", "post_month"), verbose = TRUE)
})
names(iv_pooled_by_payer) <- payer_buckets

cat("\n=== Main IV coefficient by payer class (CT Other excluded) ===\n")
etable(iv_pooled_by_payer,
       keep    = "%fit_n_prior_posters",
       headers = payer_buckets)

for (pb in payer_buckets) {
  wald_f <- fitstat(iv_pooled_by_payer[[pb]], "ivwald")[["ivwald1::n_prior_posters"]]$stat
  cat(sprintf("  %-12s Wald F: %.1f\n", pb, wald_f))
}

# ---------------------------------------------------------------------
# 3. Shoppability gradient, per payer bucket (CT Other excluded)
#    Stage 1: run_iv_by_service()
#    Stage 2: build_meta_data() + run_meta_regressions()
# ---------------------------------------------------------------------
shop_gradient_by_payer <- lapply(payer_buckets, function(pb) {
  res_service <- run_iv_by_service(
    df_payer_panels[[pb]], "ln_median_price",
    "system_peer_pressure_county_9m",
    c("County_State", "post_month"), min_obs = 100
  )
  meta_data <- build_meta_data(res_service)
  models    <- run_meta_regressions(meta_data)
  list(res_service = res_service, meta_data = meta_data, models = models)
})
names(shop_gradient_by_payer) <- payer_buckets

for (pb in payer_buckets) {
  print_meta_results(shop_gradient_by_payer[[pb]]$models, label = paste(pb, "(CT Other excluded)"))
}

# ---------------------------------------------------------------------
# 4. Summary table: Simple WLS shoppability gradient by payer bucket
# ---------------------------------------------------------------------
shop_gradient_summary <- do.call(rbind, lapply(payer_buckets, function(pb) {
  m <- shop_gradient_by_payer[[pb]]$models$simple
  co <- coef(m)["is_shoppableTRUE"]
  se_v <- summary(m)$coefficients["is_shoppableTRUE", "Std. Error"]
  n_services <- nrow(shop_gradient_by_payer[[pb]]$meta_data)
  data.frame(payer_bucket = pb, gradient_pct = co, se_pct = se_v,
             t_stat = co / se_v, n_services = n_services,
             row.names = NULL)
}))

cat("\n=== Shoppability gradient (Simple WLS), CT Other excluded ===\n")
shop_gradient_summary %>% as.data.frame() %>% print()

write.csv(shop_gradient_summary,
          "EXPORT_CSV_shop_gradient_by_payer_no_ctother.csv", row.names = FALSE)

# ---------------------------------------------------------------------
# 5. Cell-size / thinness diagnostics (mirrors payer_cell_diagnostics.R)
# ---------------------------------------------------------------------
cell_diagnostics_by_payer <- lapply(payer_buckets, function(pb) {
  panel <- df_payer_panels[[pb]]
  res   <- shop_gradient_by_payer[[pb]]$res_service
  
  cell_counts <- panel %>%
    group_by(SERVICE_GROUP) %>%
    summarise(
      n_obs_raw   = n(),
      n_hospitals = n_distinct(HOSPITAL_ID),
      n_counties  = n_distinct(County_State),
      .groups = "drop"
    )
  
  res %>%
    left_join(cell_counts, by = "SERVICE_GROUP") %>%
    mutate(payer_bucket = pb) %>%
    select(payer_bucket, SERVICE_GROUP, n_obs, n_obs_raw, n_hospitals,
           n_counties, fs_f, estimate_pct, se_pct, pval)
})
names(cell_diagnostics_by_payer) <- payer_buckets
cell_diagnostics_all <- do.call(rbind, cell_diagnostics_by_payer)

cat("\n=== Cell size summary by payer bucket (CT Other excluded) ===\n")
cell_diagnostics_all %>%
  group_by(payer_bucket) %>%
  summarise(
    min_n_obs        = min(n_obs, na.rm = TRUE),
    median_n_obs     = median(n_obs, na.rm = TRUE),
    min_fs_f         = min(fs_f, na.rm = TRUE),
    pct_fs_f_below10 = mean(fs_f < 10, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  as.data.frame() %>% print()

write.csv(cell_diagnostics_all,
          "EXPORT_CSV_cell_diagnostics_by_payer_no_ctother.csv", row.names = FALSE)

# ---------------------------------------------------------------------
# 6. LOO stability check (CT Other already excluded, so this shows
#    stability among the remaining 45 services)
# ---------------------------------------------------------------------
loo_gradient_by_payer <- lapply(payer_buckets, function(pb) {
  md <- shop_gradient_by_payer[[pb]]$meta_data
  services <- unique(md$SERVICE_GROUP)
  
  loo_results <- lapply(services, function(sg) {
    md_sub <- md %>% filter(SERVICE_GROUP != sg)
    m <- tryCatch(
      lm(estimate_pct ~ is_shoppable, data = md_sub, weights = weight),
      error = function(e) NULL
    )
    if (is.null(m)) return(NULL)
    data.frame(dropped_service = sg,
               gradient_pct = coef(m)["is_shoppableTRUE"], row.names = NULL)
  })
  do.call(rbind, loo_results) %>% mutate(payer_bucket = pb)
})
names(loo_gradient_by_payer) <- payer_buckets

cat("\n=== LOO gradient stability by payer bucket (CT Other excluded) ===\n")
for (pb in payer_buckets) {
  rng <- range(loo_gradient_by_payer[[pb]]$gradient_pct, na.rm = TRUE)
  full_est <- coef(shop_gradient_by_payer[[pb]]$models$simple)["is_shoppableTRUE"]
  cat(sprintf("  %-12s full estimate: %.3f | LOO range: [%.3f, %.3f]\n",
              pb, full_est, rng[1], rng[2]))
}

# ---------------------------------------------------------------------
# 7. Full service x payer breakdown table + forest plot (CT Other
#    excluded) -- mirrors service_by_payer_breakdown.R
# ---------------------------------------------------------------------
full_breakdown <- lapply(payer_buckets, function(pb) {
  shop_gradient_by_payer[[pb]]$res_service %>% mutate(payer_bucket = pb)
}) %>% bind_rows()

full_breakdown <- full_breakdown %>%
  mutate(
    category = case_when(
      SERVICE_GROUP %in% shoppable_v2    ~ "Shoppable",
      SERVICE_GROUP %in% nonshoppable_v2 ~ "Non-Shoppable",
      TRUE                               ~ "Intermediate"
    ),
    category = factor(category,
                      levels = c("Shoppable", "Intermediate", "Non-Shoppable")),
    payer_bucket = factor(payer_bucket,
                          levels = c("Commercial", "Medicaid", "MedicareAdv")),
    sig = pval < 0.05
  )

cat(sprintf("\n=== Full service x payer breakdown (CT Other excluded): %d rows (expect 45 x 3 = 135) ===\n",
            nrow(full_breakdown)))

full_breakdown_display <- full_breakdown %>%
  select(payer_bucket, SERVICE_GROUP, category, n_obs, fs_f,
         estimate_pct, se_pct, pval, sig) %>%
  arrange(payer_bucket, estimate_pct)

full_breakdown_display %>% as.data.frame() %>% print()

write.csv(full_breakdown_display,
          "EXPORT_CSV_full_service_by_payer_breakdown_no_ctother.csv", row.names = FALSE)

library(kableExtra)
latex_table <- full_breakdown_display %>%
  mutate(
    estimate_pct = round(estimate_pct, 2),
    se_pct       = round(se_pct, 2),
    pval         = signif(pval, 3),
    fs_f         = round(fs_f, 1)
  ) %>%
  select(-sig) %>%
  kbl(
    format = "latex", booktabs = TRUE, longtable = TRUE,
    col.names = c("Payer", "Service Group", "Category", "N", "First-Stage F",
                  "Estimate (pp)", "SE (pp)", "p-value"),
    caption = "Service-level IV estimates by payer class, CT Other excluded",
    label = "tab:service_by_payer_no_ctother"
  ) %>%
  kable_styling(latex_options = c("repeat_header"))

writeLines(as.character(latex_table), "tab_service_by_payer_no_ctother.tex")

fig_service_by_payer_forest <- ggplot(full_breakdown,
                                      aes(x = estimate_pct,
                                          y = reorder(SERVICE_GROUP, estimate_pct),
                                          color = category)) +
  geom_vline(xintercept = 0, color = "grey50", linewidth = 0.4, linetype = "dashed") +
  geom_errorbarh(
    aes(xmin = estimate_pct - 1.96 * se_pct,
        xmax = estimate_pct + 1.96 * se_pct),
    height = 0, linewidth = 0.5
  ) +
  geom_point(aes(shape = sig), size = 2) +
  scale_color_manual(values = c(
    "Shoppable"      = FSU_GARNET,
    "Intermediate"   = "grey40",
    "Non-Shoppable"  = FSU_GOLD
  )) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1), guide = "none") +
  facet_wrap(~ payer_bucket, ncol = 3) +
  labs(
    x = "IV Estimate (% per additional prior poster)",
    y = NULL,
    color = NULL,
    title = "Service-Level IV Estimates by Payer Class",
    subtitle = "CT Other excluded | Filled point = p < 0.05 | Theory V2 shoppability classification"
  ) +
  theme_bw(base_size = 9) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 6.5)
  )

print(fig_service_by_payer_forest)

ggsave("fig_service_by_payer_forest_no_ctother.pdf", fig_service_by_payer_forest,
       width = 12, height = 9)

# ---------------------------------------------------------------------
# 8. Sanity check: mean of service-level estimates vs pooled IV
# ---------------------------------------------------------------------
cat("\n=== Sanity check: mean(service-level) vs pooled IV, CT Other excluded ===\n")
for (pb in payer_buckets) {
  mean_svc_est <- full_breakdown %>%
    filter(payer_bucket == pb) %>%
    summarise(mean_est = mean(estimate_pct, na.rm = TRUE)) %>%
    pull(mean_est)
  pooled_est <- coef(iv_pooled_by_payer[[pb]])["fit_n_prior_posters"] * 100
  cat(sprintf("  %-12s mean(service-level) = %.3f pp | pooled IV = %.3f pp\n",
              pb, mean_svc_est, pooled_est))
}

# ---------------------------------------------------------------------
# 9. Side-by-side comparison: gradient WITH vs WITHOUT CT Other
#    (requires having already run the original all-services version;
#    paste in the original shop_gradient_summary values here, or
#    re-load from EXPORT_CSV_shop_gradient_by_payer.csv if saved)
# ---------------------------------------------------------------------
tryCatch({
  original <- fread("EXPORT_CSV_shop_gradient_by_payer.csv") %>%
    rename(gradient_pct_with_ctother = gradient_pct,
           se_pct_with_ctother = se_pct)
  comparison <- shop_gradient_summary %>%
    select(payer_bucket, gradient_pct_no_ctother = gradient_pct,
           se_pct_no_ctother = se_pct) %>%
    left_join(original %>% select(payer_bucket, gradient_pct_with_ctother,
                                  se_pct_with_ctother),
              by = "payer_bucket")
  cat("\n=== Gradient comparison: with vs without CT Other ===\n")
  comparison %>% as.data.frame() %>% print()
}, error = function(e) {
  cat("\n(Skipping with/without comparison -- original EXPORT_CSV_shop_gradient_by_payer.csv not found in working directory)\n")
})




# ----------------------------------------------------------------------------
# Section 22b: Within-System Price Dispersion, split by payer bucket
# Reuses the same multi_county_systems object and identical logic as Section 22,
# but runs it separately on each payer-split panel from Section 23.
# ----------------------------------------------------------------------------

compute_dispersion <- function(df, multi_county_systems) {
  df %>%
    filter(HEALTH_SYSTEM_NAME %in% multi_county_systems) %>%
    group_by(HEALTH_SYSTEM_NAME, SERVICE_GROUP) %>%
    summarise(
      n_hospitals     = n(),
      mean_lnprice    = mean(ln_median_price, na.rm = TRUE),
      sd_across_hosps = sd(ln_median_price, na.rm = TRUE),
      .groups         = "drop"
    ) %>%
    filter(n_hospitals > 1)
}

dispersion_commercial <- compute_dispersion(df_payer_panels[["Commercial"]], multi_county_systems)
dispersion_medicaid    <- compute_dispersion(df_payer_panels[["Medicaid"]],   multi_county_systems)

# Service-group-level medians, same as your existing Step 3 output
cat("\n=== Commercial: Within-System Dispersion by Service Group ===\n")
dispersion_commercial %>%
  group_by(SERVICE_GROUP) %>%
  summarise(median_sd = median(sd_across_hosps, na.rm = TRUE),
            mean_sd   = mean(sd_across_hosps, na.rm = TRUE)) %>%
  as.data.frame() %>% print()

cat("\n=== Managed Medicaid: Within-System Dispersion by Service Group ===\n")
dispersion_medicaid %>%
  group_by(SERVICE_GROUP) %>%
  summarise(median_sd = median(sd_across_hosps, na.rm = TRUE),
            mean_sd   = mean(sd_across_hosps, na.rm = TRUE)) %>%
  as.data.frame() %>% print()

# ----------------------------------------------------------------------------
# Aggregate to shoppability tier (Theory V2), matching the classification
# already used in build_meta_data() / shop_gradient_by_payer
# ----------------------------------------------------------------------------

assign_tier <- function(df) {
  df %>%
    mutate(tier = case_when(
      SERVICE_GROUP %in% c("CT Lung", "Mammography") |
        grepl("Ultrasound", SERVICE_GROUP) ~ "Shoppable",
      SERVICE_GROUP == "CT Other"          ~ "Intermediate",   # matches Section 23 exclusion note
      grepl("^CT|^X-Ray", SERVICE_GROUP)   ~ "Intermediate",
      TRUE                                  ~ "Non-Shoppable"
    ))
}

tier_summary_commercial <- assign_tier(dispersion_commercial) %>%
  group_by(tier) %>%
  summarise(median_sd = median(sd_across_hosps, na.rm = TRUE), n_cells = n(), .groups = "drop")

tier_summary_medicaid <- assign_tier(dispersion_medicaid) %>%
  group_by(tier) %>%
  summarise(median_sd = median(sd_across_hosps, na.rm = TRUE), n_cells = n(), .groups = "drop")

cat("\n=== Tier-Level Median Dispersion: Commercial ===\n")
print(as.data.frame(tier_summary_commercial))

cat("\n=== Tier-Level Median Dispersion: Managed Medicaid ===\n")
print(as.data.frame(tier_summary_medicaid))

# Pooled "All Payers" figure for the table's bottom row, using the
# already-existing pooled dispersion object from Section 22
tier_summary_pooled <- assign_tier(dispersion) %>%
  group_by(tier) %>%
  summarise(median_sd = median(sd_across_hosps, na.rm = TRUE), n_cells = n(), .groups = "drop")

cat("\n=== Tier-Level Median Dispersion: All Payers (pooled, from Section 22) ===\n")
print(as.data.frame(tier_summary_pooled))


cat(sprintf("Commercial all-service median: %.4f\n",
            median(dispersion_commercial$sd_across_hosps, na.rm = TRUE)))
cat(sprintf("Medicaid all-service median: %.4f\n",
            median(dispersion_medicaid$sd_across_hosps, na.rm = TRUE)))