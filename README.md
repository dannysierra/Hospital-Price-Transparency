# When Transparency Works: Service Shoppability, Market Structure, and the Price Effects of Hospital Disclosure

**Danny Sierra | PhD Candidate, Health Economics & Public Policy | Florida State University**

---

## Overview

This repository contains the replication code for my job market paper, which estimates the causal effects of the **2021 CMS Hospital Price Transparency Rule** on negotiated prices between hospitals and commercial insurers. The rule required hospitals to publicly post machine-readable files of their negotiated rates beginning January 1, 2022.

Rather than testing for a uniform price effect, the paper asks: *when* does transparency work — and for whom? The empirical strategy exploits heterogeneity in service shoppability and local market structure to identify the conditions under which price disclosure translates into competitive pressure on negotiated rates.

**Core findings:**
1. **Pooled null effect** — on average, price transparency had no statistically significant effect on negotiated prices.
2. **Shoppability gradient** — the rule reduced prices for shoppable services (e.g., elective imaging), with no effect for non-shoppable services, across all ten classification schemes tested.
3. **Market structure heterogeneity** — effects are concentrated in commercially active markets and are absent in high-Medicaid-share and high-minority-share markets, raising equity concerns about who benefits from price transparency.

---

## Paper

> **"When Transparency Works: Service Shoppability, Market Structure, and the Price Effects of Hospital Disclosure"**
> Danny Sierra, Florida State University, 2026

*Job market paper. Draft available upon request.*

---

## Data

This project uses **proprietary commercial data** and is not fully replicable without data access. The code is provided for transparency and methodological documentation.

| Source | Description | Access |
|---|---|---|
| **Turquoise Health** | Negotiated rates from hospital machine-readable files | Licensed / commercial |
| **AHA Annual Survey** | Hospital characteristics (`ln_total_beds`, system membership, ownership) | Licensed |
| **County Health Rankings (CHR)** | County-level demographics (Medicaid share, minority share, poverty, college share, HHI) | Public |
| **CMS / IPPS Cost Reports** | Supplemental hospital financial data | Public |

The primary analysis is structured at the **county × service group × month** level, covering 46 radiological and procedural service groups.

---

## Empirical Strategy

### Identification

The paper uses an **instrumental variables (IV)** strategy to address the endogeneity of compliance timing. Hospitals that posted prices early may differ systematically from laggards in ways that also affect pricing.

**Instrument: `system_peer_pressure_county_9m`**

For each county-month, the instrument counts the number of distinct out-of-county hospitals — belonging to health systems *present in the focal county* — that posted prices in a trailing 9-month window. This captures compliance diffusion through administrative channels within multi-hospital systems (corporate mandates, legal guidance, shared compliance infrastructure) rather than local competitive responses to price disclosure.

**Focal-county restriction:** A key correction applied throughout is that only systems with at least one hospital in the focal county can exert peer pressure on that county. This eliminates cross-system contamination that inflated an earlier version of the instrument.

**Two-way clustering:** Standard errors are clustered on `County_State + post_month` throughout.

**CBSA dropped:** All CBSA-level instrument variants fail the relevance threshold (Wald F < 10) under two-way clustering and are excluded from the analysis.

### Primary Specification

Estimated via `fixest::feols`:

```
ln_price_{ict} = α + β · n_prior_posters_{it} + γ · ln_total_beds_{it} + δ_{c×s} + τ_{t} + ε_{ict}
```

where fixed effects `δ_{c×s}` are at the `county × service_group` (market) level and `τ_{t}` are post-month fixed effects. `n_prior_posters` is instrumented by `system_peer_pressure_county_9m`.

---

## Repository Structure

```
.
├── .gitignore
└── NationalRegressions_System_Peer_Pressure_Commented.R
```

---

## Script Structure

The R file is a single self-contained script (~6,300 lines) covering the full estimation pipeline. It is organized into numbered sections with tagged output labels (`OUTPUT_TABLE_*`, `OUTPUT_FIGURE_*`, `EXPORT_CSV_*`) to locate each table, figure, and CSV export in the paper.

### Section 0 — Libraries, Setup, and Constants
Loads packages, defines FSU color constants (Garnet `#782F40`, Gold `#CEB888`), sets saturation bin cutpoints, and defines the `add_service_groups()` function used by both county and city panels.

### Section 1 — Load and Prepare County Data (Primary)
Reads the county-level national price file, filters to entrant-months, constructs log price outcomes (median, mean, P25, P75, min, max, IQR), and builds `market_id` as `County_State :: SERVICE_GROUP`.

### Section 2 — County Instrument Construction (Corrected, All Variants)
Constructs the corrected `system_peer_pressure` instrument for trailing windows of 3, 6, 9, and 12 months using the focal-county restriction. Also builds the lagged variant (6–12 month window) and a large-system-exclusion variant (dropping AdventHealth, HCA, Baylor Scott & White, and HCA regional divisions).

### Section 3 — Estimation Helper Functions
Defines the core reusable functions used throughout the script:
- `run_ols_pooled()` / `run_ols_by_service()` — pooled and service-level OLS
- `run_iv_pooled()` / `run_iv_by_service()` / `run_iv_by_large_group()` / `run_iv_by_service_outcomes()` — IV estimation across outcomes, service groups, and grouping variables
- `extract_first_stage()` / `print_iv_diagnostics()` — first-stage extraction and diagnostics
- `build_meta_data()` / `run_meta_regressions()` / `print_meta_results()` — shoppability meta-regression pipeline

### Section 4 — Load and Prepare City Data (Robustness)
Reads the city-level panel and prepares it in parallel to the county data for geographic robustness checks.

### Section 5 — City Instrument Construction
Builds city-level `system_peer_pressure` variants for 6, 9, and 12-month trailing windows.

### Section 6 — County Primary Results *(Tables 2–5, Figures 2–3)*
- First-stage estimation with two-way clustering; Wald F and Wu-Hausman test
- Pooled OLS and IV across all seven price outcomes (median, mean, P25, P75, min, max, IQR)
- Service-level IV by all 46 service groups
- OLS vs. IV comparison by service (forest plot and scatter)
- County meta-regressions for the shoppability gradient
- LaTeX table output via `etable()` and `modelsummary()`

### Section 7 — City Primary Results (Robustness)
Parallel analysis to Section 6 using the city panel and `system_peer_pressure_city_12m`. Includes city meta-regression and service-level outcome sweep.

### Section 8 — Two-Way Geographic Comparison *(Table 14)*
Side-by-side county vs. city IV estimates by service group, with shoppability gradient comparison across geographic definitions.

### Section 9 — Figures *(Figures 2–3, Appendix)*
Produces all forest plots, the meta-regression scatter, the OLS vs. IV forest and scatter, and a service × outcome heatmap. Figures use `ggplot2` + `patchwork`.

### Section 10 — Shoppability Schemes and Saturation *(Table 6, Figures 4–5)*
- Implements all **ten shoppability classification schemes**: CMS Rule Definition, Theory-Based, Theory-Based V2 (MRI as non-shoppable; primary), Broad Shoppable, Imaging vs. Procedural, Split Non-Shoppable, CMS Statutory List, High vs. Low Within Modality, CT-Inclusive (All CT), CT-Inclusive Ex. Angio
- Runs meta-regressions for each scheme
- Saturation bin analysis (Low: 1–3, Moderate: 4–8, High: 9+ prior posters) by shoppability
- Saturation robustness spaghetti plot across all ten schemes

### Section 11 — Heterogeneity Analysis *(Tables 7–9, Figures 6–7)*
Estimates shoppability gradients on median splits of market characteristics:
- Hospital ownership (for-profit share)
- Education (college share)
- Insurer concentration (HHI)
- Medicare exposure (65+ share)
- Market size (population)
- Poverty rate
- Uninsured rate
- Medicaid share
- Black share
- Hispanic share
- 2×2: Race × Medicaid share *(primary equity finding)*
- 2×2: Education × Poverty rate

Each split runs the full IV and meta-regression pipeline via `run_all_results_on_split()`. Outputs include pooled IV, service-level IV, SHOP2 group IV, and shoppability meta-regressions for each subgroup.

### Section 12 — Hospital Type Heterogeneity *(Table 8)*
Splits by hospital type (teaching, non-teaching, rural, critical access) and runs IV and meta-regressions within each group.

### Section 13 — Geographic Proximity *(Table 19, Appendix)*
Defines close/far hospital pairs by distance and runs IV and shoppability meta-regressions within proximity bins. Includes city-amplification check for biopsy services.

### Section 14 — Service-Level Heterogeneity Extensions
- Meta-regression interactions: shoppability × HHI, × for-profit share, × market size, × full interaction
- City vs. county biopsy amplification comparison

### Section 15 — Instrument Variant Robustness *(Table 20)*
County robustness across window lengths (6m, 9m, 12m), the lagged instrument (6–12m), and the large-system-exclusion variant. City robustness across 6m, 9m, 12m windows.

### Section 16 — Leave-One-System-Out *(Table 21)*
Sequentially drops each major health system and re-estimates the primary county IV. Reports LOO gradient range across systems.

### Section 17 — Placebo and Falsification Tests *(Table 23, Figure 8)*
- **Pre-trend lead test:** Adds forward leads of the instrument to test for pre-trends; reports joint F-test on leads
- **Randomization inference:** Permutes `system_peer_pressure_county_9m` and re-estimates 1,000 times; compares actual Wald F to the permutation distribution
- **Enforcement controls:** Adds CMS enforcement intensity as a control variable

### Section 18 — Clustering Robustness *(Table 22)*
Compares two-way clustering (`County_State + post_month`) to one-way clustering (`County_State` only) and to heteroskedasticity-robust SEs.

### Section 19 — Summary Statistics *(Table 1)*
Constructs the county-level summary statistics table (means, SDs, N) and exports LaTeX output.

### Section 20 — PDS-LASSO (Machine Learning Extension)
Post-double-selection LASSO via the `hdm` package. Selects high-dimensional controls from the CHR covariate set and re-estimates the primary IV with LASSO-selected controls. Reports coefficient stability relative to the primary specification.

### Section 21 — Causal Forest (Machine Learning Extension)
Instrumental forest via the `grf` package to estimate heterogeneous treatment effects. Residualizes outcome, treatment, and instrument on market and month fixed effects before fitting the forest. Outputs: variable importance plot, distribution of predicted effects, and forest-implied heterogeneity summaries (high/low splits on college share, poverty, Medicaid share, Black share, and population).

### Within-System Price Dispersion Check (End of Script)
Computes within-system, cross-county dispersion in log prices by service group to defend the exclusion restriction against centralized-contracting concerns. Finds substantial within-system price variation, inconsistent with uniform pricing across counties.

---

## Output Index

All outputs are tagged in the script for quick navigation.

**Search tags:** `OUTPUT_TABLE_*`, `OUTPUT_FIGURE_*`, `EXPORT_CSV_*`

| Tag | Paper Location | Description |
|---|---|---|
| `OUTPUT_TABLE_01_SUMSTATS` | Table 1 | Summary statistics, county level |
| `OUTPUT_TABLE_02_FIRST_STAGE` | Table 2 | First-stage estimates, county panel |
| `OUTPUT_TABLE_03_POOLED_OLS_IV` | Tables 3–4 | Pooled OLS vs IV, all outcomes |
| `OUTPUT_TABLE_05B_SERVICE_OLS_IV` | Table 5 | Service-level OLS vs IV |
| `OUTPUT_TABLE_06_METAREG` | Table 6 | Shoppability meta-regression |
| `OUTPUT_TABLE_14_GEO_COMPARISON` | Table 14 | County vs city geographic robustness |
| `OUTPUT_TABLE_20_IV_VARIANTS` | Table 20 | Instrument window/variant robustness |
| `OUTPUT_TABLE_21_LOO` | Table 21 | Leave-one-system-out |
| `OUTPUT_TABLE_22_CLUSTERING` | Table 22 | Two-way vs one-way clustering |
| `OUTPUT_TABLE_23_PLACEBO_LEADS` | Table 23 | Pre-trend leads and falsification |
| `OUTPUT_FIGURE_01_FIRST_STAGE_WALD` | Figure 1 | First-stage Wald F by service group |
| `OUTPUT_FIGURE_02_METAREG_COMBINED` | Figure 2 | Forest plot + meta-regression scatter |
| `OUTPUT_FIGURE_03_OLS_IV_FOREST` | Figure 3 | OLS vs IV forest plot |
| `OUTPUT_FIGURE_03B_OLS_IV_SCATTER` | Figure 3 | OLS vs IV scatter |
| `OUTPUT_FIGURE_04_SATURATION_SHOPPABILITY` | Figure 4 | Saturation bin effects |
| `OUTPUT_FIGURE_05_SATURATION_ROBUST` | Figure 5 | Saturation robustness across schemes |
| `OUTPUT_FIGURE_06_HETEROGENEITY_FOREST` | Figure 6 | Heterogeneity forest plot |
| `OUTPUT_FIGURE_07_RACE_MEDICAID_HEATMAP` | Figure 7 | Race × Medicaid gradient heatmap |
| `OUTPUT_FIGURE_08_RANDOMIZATION` | Figure 8 | Randomization inference histogram |

---

## R Packages

```r
# Estimation
library(fixest)       # IV, OLS, two-way FE, two-way clustering
library(hdm)          # Post-double-selection LASSO
library(grf)          # Causal / instrumental forests

# Data manipulation
library(data.table)
library(dplyr)
library(tidyr)
library(purrr)
library(lubridate)
library(here)

# Visualization
library(ggplot2)
library(ggrepel)
library(patchwork)
library(forcats)
library(stringr)

# Table output
library(stargazer)
library(knitr)
library(kableExtra)
library(modelsummary)

# Spatial
library(geosphere)
```

R version: 2026.06.0+242. Package versions available upon request.

---

## Contact

**Danny Sierra**
PhD Candidate, Department of Economics
Florida State University
Ds22c@fsu.edu

---

## Citation

> Sierra, Danny. "When Transparency Works: Service Shoppability, Market Structure, and the Price Effects of Hospital Disclosure." Working Paper, 2026.
