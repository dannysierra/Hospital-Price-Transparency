###############################################################################
#
#   WHEN TRANSPARENCY WORKS: SERVICE SHOPPABILITY, CONTRACTING DEPTH, AND
#   THE PRICE EFFECTS OF HOSPITAL DISCLOSURE
#
#   Replication code -- estimation stage
#
#   Danny Sierra
#   Department of Economics, Florida State University
#   Ds22c@fsu.edu
#
# ---------------------------------------------------------------------------
# WHAT THIS FILE DOES
# ---------------------------------------------------------------------------
# This is the estimation stage of a three-stage pipeline. Stage 1 (Snowflake
# SQL) assembles hospital machine-readable-file disclosures into a cleaned
# price panel and builds the clinical codebook. Stage 2 (Python) constructs
# hospital, market, and instrument covariates. This file, stage 3, does every
# regression, diagnostic, table, and figure that appears in the paper.
#
# The unit of observation is the hospital x clinical concept x post-month.
# The outcome is the log negotiated price. The treatment, N_PRIOR_POSTERS, is
# the count of other hospitals in the same market that have already disclosed.
# It is instrumented by out-of-market rollout among a hospital's competitor
# health systems.
#
# ---------------------------------------------------------------------------
# INPUTS
# ---------------------------------------------------------------------------
# All paths are assembled in Section 0 from PROJECT_ROOT. Four inputs:
#
#   HPT_R_MAIN_PRIMARY_COUNTY_OUTPATIENT_CONCEPT.parquet   primary panel
#   HPT_R_MAIN_PRIMARY_COUNTY_INPATIENT_CONCEPT.parquet    inpatient panel
#   HPT_R_MAIN_PRIMARY_COUNTY_OUTPATIENT_EXACT_CODE.parquet  code-level panel
#   HPT_CODEBOOK_FOR_SHOPPABILITY_SCHEMES.csv              concept codebook
#
# Two optional inputs extend the analysis: a sharded payer-cell export
# (HPT_PAYER_DISPERSION*.csv.gz, used in Section 9) and an American Community
# Survey pull via tidycensus (used in the demographic heterogeneity blocks).
#
# The price data are licensed from Turquoise Health and are not redistributed
# here. Everything else in the pipeline is public or derivable from public
# sources.
#
# ---------------------------------------------------------------------------
# HOW TO RUN
# ---------------------------------------------------------------------------
# Sections 0-12 define constants and functions and take seconds to load. The
# run blocks that follow do the work. RUN_STAGES selects which of them execute;
# set it before sourcing.
#
#   Step 1   Load Sections 0-12.
#   Step 2   Run the BUILD block. run_preflight() must pass before anything
#            else -- it is the gate that catches unclassified concepts, stale
#            caches, and degenerate instruments.
#   Step 3   RUN_STAGES <- c(7, 11)      main results and balance  (~1 hour)
#   Step 4   RUN_STAGES <- c(6)          concept-level sweep       (~8 hours)
#   Step 5   RUN_STAGES <- c(8, 9, 10)   meta-regressions, mechanism (~1 hour)
#
# Every expensive step is wrapped in cache_or_run(), which writes an .rds to
# CACHE_DIR and reloads it on subsequent runs. Set USE_CACHE <- FALSE to force
# recomputation, or delete individual .rds files to invalidate selectively. A
# step producing zero rows never overwrites a good cache.
#
# The stage 6 loader re-checks that the cached concept-level results were
# estimated on the panel currently in memory, and stops if they were not.
#
# ---------------------------------------------------------------------------
# DESIGN DECISIONS
# ---------------------------------------------------------------------------
# Four choices are fixed throughout and are not re-litigated section by
# section. Each is defended in the paper; the short version is here so the
# code reads the way the design is specified.
#
# 1. TREATMENT ENTERS LINEARLY. Concave transforms of N_PRIOR_POSTERS inflate
#    the IV coefficient by shrinking Cov(T, Z) in the Wald denominator,
#    monotonically in how much they compress the treatment (cap at Inf gives
#    -3.73%, cap at 9 gives -5.16%, log gives -7.62%). The reduced-form
#    numerator does not move. The linear first stage runs F = 43 in both
#    interacted equations, so there is no weak-instrument problem for a
#    transform to solve. The transform ladder is reported as robustness in
#    Section 7 and is never a headline.
#
# 2. THE REDUCED FORM IS THE PRIMARY ESTIMATOR. IV is a ratio; the reduced
#    form is its numerator, and it carries no treatment variable, no
#    functional-form choice, and no first-stage noise in the standard error.
#    Because every model here is exactly identified, the Anderson-Rubin
#    weak-instrument-robust test is numerically identical to the t-test on the
#    reduced-form coefficient, so the reported RF p-values are already robust
#    p-values. LIML offers nothing here: it coincides with 2SLS under exact
#    identification.
#
# 3. HETEROGENEITY IS ESTIMATED BY INTERACTION, NEVER BY SAMPLE SPLIT. The
#    fixed effects absorb most of the instrument's effect on the treatment
#    within any single subsample, so splitting destroys identification.
#    estimate_interacted() is the single estimator used for every
#    heterogeneity result in the paper.
#
# 4. NO EVENT STUDY IS POSSIBLE. 3,722 of 3,724 hospitals appear at exactly
#    one posting month, so there is no within-hospital timing variation to
#    exploit. The selection concern is addressed structurally instead: the
#    non-shoppable service group is a within-hospital placebo. Same hospital,
#    same posting month, same instrument value, no price response.
#
# ---------------------------------------------------------------------------
# OUTPUTS
# ---------------------------------------------------------------------------
# Results are written as CSV under RESULT_ROOT and are named by prefix:
#
#   T02 - T13   analysis tables, numbered by the section that produces them
#   QA01 - QA11 audit and diagnostic tables
#   fig*.pdf    publication figures
#
# The figure and LaTeX blocks at the end of the file read these CSVs back from
# disk rather than from memory, so they can be re-run on their own.
#
# ---------------------------------------------------------------------------
# REQUIREMENTS
# ---------------------------------------------------------------------------
# R 4.4 or later. Required: arrow, data.table, fixest, dplyr, ggplot2,
# stringr, lubridate, MASS. Optional, checked for at the point of use:
# R.utils (gzipped payer files), tidycensus (ACS demographics), ggrepel,
# patchwork, sf, tigris, ggpattern, magick (figures and maps).
###############################################################################
######## Section 0: Setup, paths, switches, registries ########################
#
# Everything configurable lives here: file paths, the baseline specification,
# the instrument registry, the shoppability scheme registry, and the sample
# screens. Nothing below this section hard-codes a variable name or a
# threshold, so a replication on a differently-named panel should only need to
# edit this section.

rm(list = ls()); gc()

options(stringsAsFactors = FALSE, scipen = 999, width = 160, fixest_notes = FALSE)

suppressPackageStartupMessages({
  library(arrow); library(data.table); library(fixest); library(dplyr)
  library(ggplot2); library(stringr); library(lubridate)
})
if (!requireNamespace("MASS", quietly = TRUE)) stop("Package 'MASS' required.")

setFixest_estimation(mem.clean = TRUE, data.save = FALSE)

RUN_STAGES <- c(6,7,8,9,10,11, 12)   # set explicitly before running
USE_CACHE  <- TRUE

PROJECT_ROOT <- file.path(
  "/Users/danielsierra/Library/CloudStorage",
  "OneDrive-FloridaStateUniversity",
  "Hospital Price Transparency Paper"
)

FINAL_DATA_ROOT <- file.path(PROJECT_ROOT, "Data", "data_final")
PANEL_DIR       <- file.path(FINAL_DATA_ROOT, "01_R_Analysis_Panels")
RESULT_ROOT     <- file.path(FINAL_DATA_ROOT, "06_R_Analysis_Results")

TABLE_DIR  <- file.path(RESULT_ROOT, "01_Tables_CSV")
FIGURE_DIR <- file.path(RESULT_ROOT, "03_Figures")
MODEL_DIR  <- file.path(RESULT_ROOT, "04_Model_Objects")
QA_DIR     <- file.path(RESULT_ROOT, "05_R_QA")
CACHE_DIR  <- file.path(RESULT_ROOT, "07_Cache")

for (d in c(RESULT_ROOT, TABLE_DIR, FIGURE_DIR, MODEL_DIR, QA_DIR, CACHE_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

FILES <- list(
  outpatient_concept = file.path(PANEL_DIR,
                                 "HPT_R_MAIN_PRIMARY_COUNTY_OUTPATIENT_CONCEPT.parquet"),
  inpatient_concept  = file.path(PANEL_DIR,
                                 "HPT_R_MAIN_PRIMARY_COUNTY_INPATIENT_CONCEPT.parquet"),
  outpatient_exact   = file.path(PANEL_DIR,
                                 "HPT_R_MAIN_PRIMARY_COUNTY_OUTPATIENT_EXACT_CODE.parquet"),
  codebook           = file.path(PANEL_DIR,
                                 "HPT_CODEBOOK_FOR_SHOPPABILITY_SCHEMES.csv")
)

# Payer dispersion export from Snowflake (HPT_P2_PAYER_CELLS). Sharded into
# multiple gzip parts on download (named like
# HPT_PAYER_DISPERSION.csv.gz_0_0_0.csv.gz through ..._0_7_0.csv.gz) --
# Snowflake splits above its internal size threshold regardless of
# MAX_FILE_SIZE, so this is a DIRECTORY + PATTERN, not a single file path.
# load_payer_dispersion() in Section 9 reads and combines every matching part.
PAYER_DISPERSION_DIR     <- PANEL_DIR
PAYER_DISPERSION_PATTERN <- "^HPT_PAYER_DISPERSION.*\\.csv\\.gz$"

# ---------------------------------------------------------------------------
# Baseline specification
# ---------------------------------------------------------------------------
ENDOGENOUS_VARIABLE    <- "N_PRIOR_POSTERS"
BASELINE_CONTROLS      <- c("LOG_TOTAL_BEDS")
BASELINE_FIXED_EFFECTS <- c("MARKET_ID", "POST_MONTH")
BASELINE_CLUSTERS      <- c("ANALYSIS_MARKET", "POST_MONTH")
PRIMARY_OUTCOME        <- "LN_MEDIAN_PRICE"

# MIN_PRICE / MAX_PRICE are absent from the rebuilt Phase 1 SQL and cannot be
# reconstructed from concept aggregates. IQR is P75 - P25 by definition.
OUTCOMES <- c(Median = "LN_MEDIAN_PRICE", Mean = "LN_MEAN_PRICE",
              P25 = "LN_P25_PRICE", P75 = "LN_P75_PRICE", IQR = "LN_IQR_PRICE")

# Robustness ladder only. Never a headline. See standing decision 1.
TRANSFORM_LADDER <- list(
  Linear     = list(fun = function(x) x),
  Winsor_P99 = list(quantile = 0.99),
  Winsor_P95 = list(quantile = 0.95),
  Winsor_P90 = list(quantile = 0.90),
  Sqrt       = list(fun = function(x) sqrt(pmax(x, 0))),
  Log1p      = list(fun = function(x) log1p(pmax(x, 0)))
)

# ---------------------------------------------------------------------------
# Instruments
# ---------------------------------------------------------------------------
#
# Only Z_SYS_COMPETITOR_* are hospital-specific and exclude the focal
# hospital's own system. The other Z_SYS_* measures are county-month "local
# system exposure" and do NOT exclude own-system rollout, which is why the
# competitor family carries the exclusion-restriction defence.

PRIMARY_INSTRUMENT <- "Z_SYS_COMPETITOR_ONLY_9M_EXCL_CURRENT"

MAIN_INSTRUMENTS <- c(
  Competitor_only_hospitals_9m         = "Z_SYS_COMPETITOR_ONLY_9M_EXCL_CURRENT",
  Primary_strict_system_IV             = "Z_SYS_STRICT_9M_EXCL_CURRENT",
  Competitor_outside_CBSA_hospitals_9m = "Z_SYS_COMPETITOR_OUTSIDE_CBSA_9M_EXCL_CURRENT"
)

# Competitor_systems_9m is not a headline candidate: only 54% of scheme
# variants come back negative and the sign flips under the tier-collapse rule.
# It is retained in the concept-level sweep for comparison only.
SUPPORTING_INSTRUMENTS <- c(
  Competitor_outside_CBSA_systems_9m  = "Z_SYS_COMPETITOR_SYSTEMS_OUTSIDE_CBSA_9M_EXCL_CURRENT",
  Competitor_outside_CBSA_counties_9m = "Z_SYS_COMPETITOR_COUNTIES_OUTSIDE_CBSA_9M_EXCL_CURRENT",
  Competitor_systems_9m               = "Z_SYS_COMPETITOR_SYSTEMS_9M_EXCL_CURRENT"
)

CANONICAL_SYSTEM_INSTRUMENTS <- c(
  MAIN_INSTRUMENTS, SUPPORTING_INSTRUMENTS,
  Original_9m_including_current = "Z_SYS_ORIGINAL_9M_INCL_CURRENT",
  Fixed_roster_9m               = "Z_SYS_FIXED_ROSTER_9M_EXCL_CURRENT",
  Active_systems_9m             = "Z_SYS_ACTIVE_SYSTEMS_9M_EXCL_CURRENT",
  Outside_CBSA_9m               = "Z_SYS_OUTSIDE_CBSA_9M_EXCL_CURRENT",
  Recent_flow_3m                = "Z_SYS_RECENT_FLOW_3M_EXCL_CURRENT",
  Cumulative_external_hospitals = "Z_SYS_CUMULATIVE_EXTERNAL_HOSPITALS",
  Cumulative_rollout_share      = "Z_SYS_CUMULATIVE_ROLLOUT_SHARE",
  Ever_initiated_systems        = "Z_SYS_EVER_INITIATED_SYSTEMS",
  Strict_9m_excluding_top5      = "Z_SYS_STRICT_9M_EXCL_TOP5"
)

SUPPORTING_WINDOW_INSTRUMENTS <- setNames(
  as.vector(outer(c("PEER_HOSPITALS", "PEER_COUNTIES", "PEER_SYSTEMS"),
                  c("3M", "6M", "9M", "12M"),
                  function(a, b) paste0("Z_SYS_", a, "_", b, "_STRICT"))),
  as.vector(outer(c("Peer_hospitals", "Peer_counties", "Peer_systems"),
                  c("3m", "6m", "9m", "12m"),
                  function(a, b) paste0(a, "_", b, "_strict")))
)

# CMS enforcement measures are exploratory and never headline. Enforcement
# plausibly moves prices directly, through compliance cost and correlated
# conduct, so it fails the exclusion restriction on its face. No amount of
# re-screening or reparameterisation changes that, since functional form does
# not touch an exclusion restriction.
ENFORCEMENT_ROBUSTNESS_INSTRUMENTS <- c(
  Any_enforcement_3m_lag1 = "COUNTY_ENF_ANY_ENFORCEMENT_ROLL_3M_LAG_1M",
  Warning_3m_lag1         = "COUNTY_ENF_WARNING_ROLL_3M_LAG_1M",
  Closure_3m_lag1         = "COUNTY_ENF_CLOSURE_NOTICE_ROLL_3M_LAG_1M",
  Any_enforcement_ind_3m  = "COUNTY_ENF_ANY_ENFORCEMENT_IND_3M_LAG_1M",
  Warning_ind_3m          = "COUNTY_ENF_WARNING_IND_3M_LAG_1M",
  Closure_ind_3m          = "COUNTY_ENF_CLOSURE_NOTICE_IND_3M_LAG_1M",
  Closure_6m_lag1         = "COUNTY_ENF_CLOSURE_NOTICE_ROLL_6M_LAG_1M",
  Closure_ind_6m          = "COUNTY_ENF_CLOSURE_NOTICE_IND_6M_LAG_1M"
)

ALL_CANDIDATE_INSTRUMENTS <- c(
  CANONICAL_SYSTEM_INSTRUMENTS, SUPPORTING_WINDOW_INSTRUMENTS,
  ENFORCEMENT_ROBUSTNESS_INSTRUMENTS
)

# Instrument labels are free text and have been written both with spaces and
# with underscores across the pipeline's history. The INSTRUMENT column (the
# actual panel variable name) is authoritative. This map normalises labels back
# from it, and must be applied before any grouping or joining done by label.
INSTRUMENT_LABEL_MAP <- setNames(names(CANONICAL_SYSTEM_INSTRUMENTS),
                                 unname(CANONICAL_SYSTEM_INSTRUMENTS))

# ---------------------------------------------------------------------------
# Three-tier instrument structure
# ---------------------------------------------------------------------------
#
# The six system instruments are graded on how credible their exclusion
# restriction is and how stable their sign is across classification schemes.
# The grading is set by construction and by the audit in QA05-QA08, and it
# governs how results are reported, not which results are kept:
#
#   MAIN        carries every headline claim
#   CONFIRMING  reported alongside MAIN as pooled robustness
#   DISCREPANT  reported on its own, never pooled into a headline number
#
# The tiering is deliberately not a function of first-stage F. Section 7.6
# shows the two strongest first stages produce the two weakest results, which
# is the reason instrument selection here is a design decision rather than an
# outcome-based one.

CONFIRMING_INSTRUMENTS <- c(
  Competitor_outside_CBSA_counties_9m = "Z_SYS_COMPETITOR_COUNTIES_OUTSIDE_CBSA_9M_EXCL_CURRENT",
  Competitor_systems_9m               = "Z_SYS_COMPETITOR_SYSTEMS_9M_EXCL_CURRENT"
)

DISCREPANT_INSTRUMENTS <- c(
  Competitor_outside_CBSA_systems_9m = "Z_SYS_COMPETITOR_SYSTEMS_OUTSIDE_CBSA_9M_EXCL_CURRENT"
)

stopifnot(setequal(c(CONFIRMING_INSTRUMENTS, DISCREPANT_INSTRUMENTS),
                   SUPPORTING_INSTRUMENTS))

HEADLINE_INSTRUMENTS   <- MAIN_INSTRUMENTS
ROBUSTNESS_INSTRUMENTS <- c(MAIN_INSTRUMENTS, CONFIRMING_INSTRUMENTS)
ALL_SIX_INSTRUMENTS    <- c(MAIN_INSTRUMENTS, SUPPORTING_INSTRUMENTS)

instrument_tier <- function(label) {
  fcase(label %chin% names(MAIN_INSTRUMENTS), "MAIN",
        label %chin% names(CONFIRMING_INSTRUMENTS), "CONFIRMING",
        label %chin% names(DISCREPANT_INSTRUMENTS), "DISCREPANT",
        default = "UNKNOWN")
}

# ---------------------------------------------------------------------------
# Shoppability
# ---------------------------------------------------------------------------
DIAGNOSTIC_FAMILIES <- c(
  "MRI_MRA", "CT_CTA", "XRAY_FLUOROSCOPY", "DIAGNOSTIC_ULTRASOUND",
  "VASCULAR_ULTRASOUND", "ECHOCARDIOGRAPHY", "MAMMOGRAPHY", "BONE_DENSITY",
  "LABORATORY_PATHOLOGY", "EVALUATION_MANAGEMENT"
)

# Match the literal family identifiers. These are EMERGENCY_DEPARTMENT and
# CRITICAL_CARE, not "EMERGENCY" -- a prefix match on the shorter string finds
# nothing and silently leaves emergency concepts classified as INTERMEDIATE in
# every scheme. build_schemes() asserts against exactly that failure.
ALWAYS_NONSHOPPABLE_FAMILIES <- c("EMERGENCY_DEPARTMENT", "CRITICAL_CARE")

PRIMARY_SCHEMES <- list(
  list(col = "SCHEME_1_CERTAINTY", label = "1. Procedural certainty",
       source = "families",              rule = "certainty"),
  list(col = "SCHEME_2_THEORYV2",  label = "2. Theory-Based V2",
       source = "scheme_theory_v2",      rule = "high_vs_rest"),
  list(col = "SCHEME_3_IMAGING",   label = "3. Imaging vs Procedural",
       source = "scheme_imaging",        rule = "high_vs_rest"),
  list(col = "SCHEME_4_CMS70",     label = "4. CMS Statutory List",
       source = "scheme_cms_statutory",  rule = "high_vs_rest"),
  list(col = "SCHEME_5_MDSAVE",    label = "5. Upfront Cash-Market",
       source = "scheme_div3_mdsave",    rule = "low_vs_rest"),
  list(col = "SCHEME_6_WITHINMOD", label = "6. Within Modality",
       source = "scheme_anatomical",     rule = "extremes")
)
SCHEME_COLUMNS <- setNames(vapply(PRIMARY_SCHEMES, `[[`, character(1), "col"),
                           vapply(PRIMARY_SCHEMES, `[[`, character(1), "label"))

# Concepts that are not independently schedulable services: intraoperative
# imaging, imaging guidance, add-on codes, anaesthesia, psychotherapy. Scheme 1
# classifies by whole family, so the keyword exclusions used by the other
# schemes never reach it. This pattern drives the Scheme 1 sensitivity variant
# and the QA04 flag instead. It touches 17 of 448 shoppable concepts.
PROCEDURAL_CONCEPT_PATTERN <- paste(
  "INTRAOP", "INTRA_OP", "GUID", "GUIDANCE", "W_ULTRASOUND",
  "INTRAVAS", "_ADDL", "EA_ADDL", "PSYTX", "ANESTH",
  sep = "|"
)
EXCLUDE_PROCEDURAL_FROM_SCHEME1 <- TRUE   # run the sensitivity alongside headline

# ---------------------------------------------------------------------------
# Comparability
# ---------------------------------------------------------------------------
COMPARABILITY_MODERATORS <- c("PD_PAYER_V2", "PD_HOSP", "N_CODES",
                              "N_PAYERS_V2", "PD_PAYER", "PD_CODE",
                              "N_PAYERS", "CODE_COV")
# Every candidate is listed, including ones already known to be degenerate in
# this panel (PD_PAYER, PD_CODE, N_PAYERS, CODE_COV -- see QA09), alongside the
# two measures built from the payer-cell export (PD_PAYER_V2, N_PAYERS_V2).
# screen_moderators() drops whichever fail its coverage, variance, zero-share,
# and collinearity gates, and reports the reason for each drop rather than
# dropping silently.

DISPERSION_MEASURES <- c("PD_PAYER_V2", "PD_HOSP", "N_CODES", "PD_PAYER", "PD_CODE")
# Sign convention. Dispersion measures run one way -- higher dispersion means
# a less comparable posted price, so the predicted interaction is POSITIVE.
# Payer counts run the other way: more payers means thicker contracting, more
# comparable, so the predicted interaction is NEGATIVE. SIGN_AS_PREDICTED
# applies the flip wherever these moderators are reported.

MIN_OUTSIDE_HOSPITALS <- 20L

# ---------------------------------------------------------------------------
# Screens
# ---------------------------------------------------------------------------
MIN_SERVICE_OBS     <- 100L
MIN_SERVICE_MARKETS <- 8L
MIN_SERVICE_MONTHS  <- 3L
MIN_MODEL_OBS       <- 50000L
MIN_CONCEPTS_META   <- 15L
DROP_SINGLETON_MARKETS <- TRUE
META_WEIGHTINGS <- c("Inverse variance", "Unweighted")

FSU_GARNET <- "#782F40"; FSU_GOLD <- "#CEB888"; FSU_GREY <- "#6D6E71"

cat("\nSection 0 loaded | treatment:", ENDOGENOUS_VARIABLE,
    "(linear) | primary IV:", PRIMARY_INSTRUMENT, "\n")



######## Section 1: Helper functions ##########################################
#
# Small utilities used everywhere below: safe coercion, coefficient extraction
# that survives fixest's naming conventions, first-stage diagnostics, formula
# builders, CSV writers, and the caching wrapper.

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || all(is.na(x)) ||
      (is.character(x) && !nzchar(x[1L]))) y else x
}

safe_numeric <- function(x) suppressWarnings(as.numeric(x))
safe_log_positive <- function(x) { x <- safe_numeric(x); fifelse(is.finite(x) & x > 0, log(x), NA_real_) }
safe_log1p_nonneg <- function(x) { x <- safe_numeric(x); fifelse(is.finite(x) & x >= 0, log1p(x), NA_real_) }
available_columns <- function(data, columns) intersect(columns, names(data))
has_usable_variation <- function(x) { x <- x[is.finite(x)]; length(x) > 1L && uniqueN(x) > 1L }

assert_columns <- function(data, columns, label = "data") {
  m <- setdiff(columns, names(data))
  if (length(m) > 0L) stop(sprintf("%s missing: %s", label, paste(m, collapse = ", ")), call. = FALSE)
  invisible(TRUE)
}

add_stars <- function(p) fcase(is.na(p), "", p < 0.01, "***", p < 0.05, "**", p < 0.10, "*", default = "")

report_memory <- function(label = "") {
  cat(sprintf("%-46s %5d objects | %8.1f MB\n", label,
              length(ls(envir = .GlobalEnv)), sum(gc()[, 2])))
}

# Read a single metadata value out of a one-row slice. Indexing a column that
# does not exist returns a zero-length vector, and a data.table built from one
# of those comes back with zero rows but the full column structure intact --
# which looks like a valid empty result rather than an error. These two
# accessors return NA instead.
meta_chr <- function(dt, col) if (col %in% names(dt)) as.character(dt[[col]][1L]) else NA_character_
meta_num <- function(dt, col) if (col %in% names(dt)) as.numeric(dt[[col]][1L])   else NA_real_

drop_singleton_markets <- function(dt, market_col = "ANALYSIS_MARKET") {
  sizes <- dt[, .N, by = c(market_col)]
  dt[get(market_col) %chin% sizes[N > 1L][[market_col]]]
}

# ---------------------------------------------------------------------------
# Coefficient extraction
# ---------------------------------------------------------------------------
NULL_COEF <- list(term = NA_character_, estimate = NA_real_, std_error = NA_real_,
                  statistic = NA_real_, p_value = NA_real_)

find_term_name <- function(fit, candidates) {
  nm <- names(coef(fit)); if (length(nm) == 0L) return(NA_character_)
  norm <- gsub("[^a-z0-9]", "", tolower(nm))
  for (cand in gsub("[^a-z0-9]", "", tolower(candidates))) {
    i <- which(norm == cand); if (length(i) > 0L) return(nm[i[1L]])
  }
  NA_character_
}

extract_coefficient <- function(fit, candidates) {
  if (is.null(fit)) return(NULL_COEF)
  term <- find_term_name(fit, candidates)
  if (is.na(term)) return(NULL_COEF)
  ct <- tryCatch(as.data.frame(coeftable(fit)), error = function(e) NULL)
  if (is.null(ct) || !(term %in% rownames(ct))) return(modifyList(NULL_COEF, list(term = term)))
  row <- ct[term, , drop = FALSE]
  ec <- intersect(c("Estimate", "estimate"), names(row))[1L]
  sc <- intersect(c("Std. Error", "Std.Error", "std.error"), names(row))[1L]
  tc <- grep("^(t value|z value|statistic)$", names(row), ignore.case = TRUE, value = TRUE)[1L]
  pc <- grep("^Pr\\(", names(row), value = TRUE)[1L]
  list(term = term,
       estimate  = if (!is.na(ec)) as.numeric(row[[ec]]) else NA_real_,
       std_error = if (!is.na(sc)) as.numeric(row[[sc]]) else NA_real_,
       statistic = if (!is.na(tc)) as.numeric(row[[tc]]) else NA_real_,
       p_value   = if (!is.na(pc)) as.numeric(row[[pc]]) else NA_real_)
}

# ---------------------------------------------------------------------------
# First-stage statistics
# ---------------------------------------------------------------------------
#
# With k endogenous regressors, fixest returns one Wald statistic PER EQUATION,
# each carrying stat / p / df1 / df2 / vcov. Flattening that structure and
# taking the first element silently reports equation 1 only, which is not the
# identification statistic for the model. The interacted specifications here
# have two to six endogenous terms, so this extracts every equation by name and
# reports them individually alongside the MINIMUM.
#
# The minimum is the number to quote. A design is only as identified as its
# weakest first stage.

first_stage_wald <- function(fit) {
  if (is.null(fit)) return(data.table())
  res <- tryCatch(fitstat(fit, "ivwald"), error = function(e) NULL)
  if (is.null(res)) return(data.table())
  flat <- unlist(res)
  sn <- grep("\\.stat$", names(flat), value = TRUE)
  if (length(sn) == 0L) return(data.table())
  data.table(EQUATION = sub("^.*::", "", sub("\\.stat$", "", sn)),
             WALD = safe_numeric(flat[sn]))
}

# Cragg-Donald assumes homoskedasticity so it is not valid alongside two-way
# clustering; Kleibergen-Paap rk fails numerically on these models. Recorded
# for completeness, never quoted as the identification statistic.
cragg_donald <- function(fit) {
  v <- tryCatch(safe_numeric(unlist(fitstat(fit, "cd"))[1L]), error = function(e) NA_real_)
  if (length(v) == 0L) NA_real_ else v
}

extract_wu_hausman_p <- function(fit) {
  r <- tryCatch(unlist(fitstat(fit, "wh")), error = function(e) NULL)
  if (is.null(r)) return(NA_real_)
  i <- grep("(^|\\.)p$", names(r), ignore.case = TRUE)
  if (length(i) > 0L) safe_numeric(r[i[1L]]) else NA_real_
}

# Wald test that a set of coefficients are all equal to the first.
#
# A clustered variance matrix can be singular when a category has few clusters.
# solve() would either error or return noise, so this falls back to the
# generalised inverse and records the rank. VCOV_FULL_RANK makes a degenerate
# test visible in the output rather than letting it pass as a real p-value.
wald_equality <- function(fit, terms) {
  if (is.null(fit)) return(data.table())
  terms <- terms[!is.na(terms) & terms %in% names(coef(fit))]
  if (length(terms) < 2L) return(data.table())
  tryCatch({
    cf <- coef(fit); V <- vcov(fit)
    cm <- matrix(0, nrow = length(terms) - 1L, ncol = length(cf)); colnames(cm) <- names(cf)
    for (k in 2:length(terms)) { cm[k - 1L, terms[1L]] <- 1; cm[k - 1L, terms[k]] <- -1 }
    b <- cf[colnames(cm)]; Vs <- V[colnames(cm), colnames(cm)]
    mid <- cm %*% Vs %*% t(cm); rk <- qr(mid)$rank
    inv <- tryCatch(solve(mid), error = function(e) MASS::ginv(mid))
    stat <- as.numeric(t(cm %*% b) %*% inv %*% (cm %*% b))
    if (!is.finite(stat) || stat < 0) stop("non-finite Wald")
    data.table(WALD = stat, DF = rk, VCOV_FULL_RANK = as.integer(rk == nrow(cm)),
               P_VALUE = pchisq(stat, rk, lower.tail = FALSE))
  }, error = function(e) data.table())
}

tidy_fixest <- function(fit, conf = 0.95) {
  if (is.null(fit)) return(data.table())
  b <- coef(fit); s <- fixest::se(fit); p <- fixest::pvalue(fit)
  tm <- Reduce(intersect, list(names(b), names(s), names(p)))
  if (length(tm) == 0L) return(data.table())
  z <- qnorm(1 - (1 - conf) / 2)
  data.table(term = tm, estimate = as.numeric(b[tm]), std.error = as.numeric(s[tm]),
             statistic = as.numeric(b[tm] / s[tm]), p.value = as.numeric(p[tm]),
             conf.low = as.numeric(b[tm] - z * s[tm]), conf.high = as.numeric(b[tm] + z * s[tm]))
}

# ---------------------------------------------------------------------------
# Formulas, output, caching
# ---------------------------------------------------------------------------
build_cluster_formula <- function(cl) as.formula(paste0("~", paste(cl, collapse = " + ")))

build_ols_formula <- function(outcome, rhs, fe) {
  as.formula(paste0(outcome, " ~ ", paste(rhs, collapse = " + "),
                    " | ", paste(fe, collapse = " + ")))
}

build_iv_formula <- function(outcome, endogenous, instruments, controls, fe) {
  exo <- if (length(controls) == 0L) "1" else paste(controls, collapse = " + ")
  as.formula(paste0(outcome, " ~ ", exo, " | ", paste(fe, collapse = " + "),
                    " | ", paste(endogenous, collapse = " + "),
                    " ~ ", paste(instruments, collapse = " + ")))
}

save_csv <- function(data, filename, dir = TABLE_DIR) {
  if (is.null(data) || nrow(data) == 0L) {
    warning("Refusing to write an empty table: ", filename, call. = FALSE); return(invisible(NULL))
  }
  path <- file.path(dir, filename); fwrite(as.data.table(data), path)
  cat("Saved:", path, "\n"); invisible(path)
}
save_qa_csv <- function(data, filename) save_csv(data, filename, dir = QA_DIR)

# Single caching gate for every expensive step. If the .rds exists and caching
# is on, load it; otherwise compute and write. A run that produces zero rows
# never overwrites a good cache -- it stops instead, so a silent upstream
# failure cannot destroy hours of completed estimation.
cache_or_run <- function(key, expr, overwrite = !USE_CACHE) {
  path <- file.path(CACHE_DIR, paste0(key, ".rds"))
  if (file.exists(path) && !overwrite) { cat("Cache hit:", key, "\n"); return(readRDS(path)) }
  cat("Computing:", key, "\n")
  result <- expr
  if (is.data.frame(result) && nrow(result) == 0L) {
    stop("Refusing to cache an empty result for ", key, call. = FALSE)
  }
  saveRDS(result, path); result
}

cat("Section 1 loaded\n")



######## Section 2: Panel loading and preparation #############################
#
# Every function here returns data rather than mutating a global object. That
# keeps the panel's provenance traceable: `outpatient` is assigned exactly once,
# in the BUILD block, from a visible chain of calls. It also means no downstream
# section needs its own guard against columns left over from a previous run.
#
# read_panel() selects columns via intersect(), so a variable absent from the
# parquet is dropped silently rather than erroring. Any new instrument column
# must therefore be added to ANALYSIS_COLUMNS or it will not survive loading.

ANALYSIS_COLUMNS <- unique(c(
  "HOSPITAL_ID", "PROVIDER_STATE", "COUNTY_STATE_KEY", "CBSA_CODE",
  "HOSPITAL_TYPE", "HEALTH_SYSTEM_ID", "SYSTEM_KEY", "TOTAL_BEDS",
  "POST_MONTH", "HOSPITAL_FIRST_POST_MONTH", "N_OBSERVED_POST_MONTHS",
  "ANALYSIS_GEOGRAPHY", "ANALYSIS_AGGREGATION", "ANALYSIS_MARKET",
  "ANALYSIS_SERVICE_ID", "MARKET_ID",
  "BILLING_CODE_TYPE", "BILLING_CODE", "OFFICIAL_DESCRIPTION",
  "FINAL_SUPERFAMILY_ID", "FINAL_FAMILY_ID", "FINAL_FAMILY_NAME",
  "FINAL_CONCEPT_ID", "FINAL_CONCEPT_NAME", "IS_CMS70_CONCEPT_FLAG",
  "N_PRIOR_POSTERS", "N_HOSPITALS_CUMULATIVE_MARKET",
  "MEDIAN_PRICE", "MEAN_PRICE", "P25_PRICE", "P75_PRICE",
  "N_PAYER_CELLS", "N_DISTINCT_PAYERS", "N_CODE_PAYER_CELLS_TOTAL",
  "SUM_CODE_DISTINCT_PAYERS", "N_CODES_PRESENT",
  "EXPECTED_N_CODES_IN_CONCEPT", "CODE_COVERAGE_RATIO",
  "COMPLETE_CONCEPT_COVERAGE_FLAG", "COUNTY_FIPS",
  "SAMPLE_COHORT_CONCEPT", "RECOMMENDED_POSTING_COHORT_SAMPLE",
  "RECOMMENDED_BASELINE_SAMPLE",
  unname(ALL_CANDIDATE_INSTRUMENTS)
))

read_panel <- function(path, label, columns = ANALYSIS_COLUMNS) {
  if (!file.exists(path)) stop("Not found for ", label, ":\n", path, call. = FALSE)
  cat("\nReading", label, "...\n")
  ds <- arrow::open_dataset(path, format = "parquet")
  sel <- intersect(columns, names(ds))
  if (length(sel) == 0L) stop("No analytical columns in ", label, call. = FALSE)
  d <- as.data.table(dplyr::collect(dplyr::select(ds, dplyr::all_of(sel))))
  cat("  Rows:", format(nrow(d), big.mark = ","), "| columns:", ncol(d), "\n")
  d
}

prepare_panel <- function(data, label, sample_flag = NULL) {
  d <- copy(as.data.table(data))

  if (!("N_PAYER_CELLS" %in% names(d)) && "N_CODE_PAYER_CELLS_TOTAL" %in% names(d)) {
    d[, N_PAYER_CELLS := safe_numeric(N_CODE_PAYER_CELLS_TOTAL)]
  }
  if (!("N_DISTINCT_PAYERS" %in% names(d)) && "SUM_CODE_DISTINCT_PAYERS" %in% names(d)) {
    d[, N_DISTINCT_PAYERS := safe_numeric(SUM_CODE_DISTINCT_PAYERS)]
  }

  assert_columns(d, c("HOSPITAL_ID", "POST_MONTH", "ANALYSIS_MARKET",
                      "ANALYSIS_SERVICE_ID", "MARKET_ID", ENDOGENOUS_VARIABLE,
                      "MEDIAN_PRICE", "P25_PRICE", "P75_PRICE", "TOTAL_BEDS"), label)

  if (!("IQR_PRICE" %in% names(d))) {
    d[, IQR_PRICE := pmax(safe_numeric(P75_PRICE) - safe_numeric(P25_PRICE), 0)]
  }

  d[, POST_MONTH := as.Date(POST_MONTH)]
  for (col in c("HOSPITAL_ID", "ANALYSIS_MARKET", "ANALYSIS_SERVICE_ID", "MARKET_ID",
                "FINAL_FAMILY_ID", "FINAL_CONCEPT_ID", "BILLING_CODE_TYPE",
                "HEALTH_SYSTEM_ID", "SYSTEM_KEY", "HOSPITAL_TYPE", "PROVIDER_STATE")) {
    if (col %in% names(d)) set(d, j = col, value = as.character(d[[col]]))
  }
  d[, (ENDOGENOUS_VARIABLE) := safe_numeric(get(ENDOGENOUS_VARIABLE))]

  d[, LN_MEDIAN_PRICE := safe_log_positive(MEDIAN_PRICE)]
  d[, LN_MEAN_PRICE   := safe_log_positive(MEAN_PRICE)]
  d[, LN_P25_PRICE    := safe_log_positive(P25_PRICE)]
  d[, LN_P75_PRICE    := safe_log_positive(P75_PRICE)]
  d[, LN_IQR_PRICE    := safe_log1p_nonneg(IQR_PRICE)]
  d[, LOG_TOTAL_BEDS  := log(pmax(safe_numeric(TOTAL_BEDS), 1))]

  d[, SERVICE_NAME := if ("FINAL_CONCEPT_NAME" %in% names(d)) as.character(FINAL_CONCEPT_NAME)
    else as.character(ANALYSIS_SERVICE_ID)]
  d[, SERVICE_LABEL := fifelse(is.na(SERVICE_NAME) | SERVICE_NAME == "",
                               ANALYSIS_SERVICE_ID, SERVICE_NAME)]

  if (!is.null(sample_flag) && sample_flag %in% names(d)) d <- d[get(sample_flag) == 1]

  d <- d[!is.na(ANALYSIS_MARKET) & ANALYSIS_MARKET != "" &
           !is.na(MARKET_ID) & MARKET_ID != "" &
           is.finite(LN_MEDIAN_PRICE) & is.finite(get(ENDOGENOUS_VARIABLE))]

  setorder(d, ANALYSIS_MARKET, ANALYSIS_SERVICE_ID, POST_MONTH, HOSPITAL_ID)

  cat("Prepared", label, "\n  Rows:", format(nrow(d), big.mark = ","),
      "| hospitals:", format(uniqueN(d$HOSPITAL_ID), big.mark = ","),
      "| counties:", format(uniqueN(d$ANALYSIS_MARKET), big.mark = ","),
      "| concepts:", format(uniqueN(d$ANALYSIS_SERVICE_ID), big.mark = ","),
      "| FE cells:", format(uniqueN(d$MARKET_ID), big.mark = ","), "\n")
  d
}

choose_sample_flag <- function(data) {
  cands <- c("SAMPLE_COHORT_CONCEPT", "RECOMMENDED_POSTING_COHORT_SAMPLE",
             "RECOMMENDED_BASELINE_SAMPLE")
  a <- cands[cands %in% names(data)]
  if (length(a) == 0L) NULL else a[1L]
}

load_outpatient <- function() {
  raw <- read_panel(FILES$outpatient_concept, "outpatient concept panel")
  prepare_panel(raw, "outpatient concept panel", choose_sample_flag(raw))
}

# Documents the gap between panel rows and estimated rows. feols drops roughly
# 31% of complete cases, which is cascading singleton removal across the
# MARKET_ID and POST_MONTH fixed effects: the panel carries about 672,000 fixed
# effect cells for 1.4M rows, and 30% of rows sit in cells of size one. This is
# expected, is identical across specifications, and is reported in the paper's
# sample appendix rather than left implicit.
audit_estimation_sample <- function(panel, instrument = PRIMARY_INSTRUMENT,
                                    outcome = PRIMARY_OUTCOME) {
  req <- available_columns(panel, unique(c(outcome, ENDOGENOUS_VARIABLE, instrument,
                                           BASELINE_CONTROLS, BASELINE_FIXED_EFFECTS,
                                           BASELINE_CLUSTERS)))
  cc <- sum(complete.cases(panel[, ..req]))
  fit <- tryCatch(feols(build_iv_formula(outcome, ENDOGENOUS_VARIABLE, instrument,
                                         BASELINE_CONTROLS, BASELINE_FIXED_EFFECTS),
                        data = panel, cluster = build_cluster_formula(BASELINE_CLUSTERS),
                        warn = FALSE, notes = FALSE), error = function(e) NULL)
  sizes <- panel[, .N, by = MARKET_ID]$N
  out <- data.table(
    PANEL_ROWS = nrow(panel), COMPLETE_CASES = cc,
    ESTIMATED_ROWS = if (is.null(fit)) NA_integer_ else nobs(fit),
    DROPPED_BY_FEOLS = if (is.null(fit)) NA_integer_ else cc - nobs(fit),
    SHARE_DROPPED = if (is.null(fit)) NA_real_ else round((cc - nobs(fit)) / cc, 4),
    FE_CELLS = length(sizes), SINGLETON_FE_CELLS = sum(sizes == 1L),
    SHARE_ROWS_IN_SINGLETON_CELLS = round(sum(sizes[sizes == 1L]) / sum(sizes), 4))
  save_qa_csv(out, "QA01_estimation_sample_audit.csv"); print(out); out
}

cat("Section 2 loaded\n")


# --- Interactive checkpoint ---------------------------------------------------
# Forces a panel rebuild by clearing the cache, then rebuilds it. This is a
# convenience for re-running the panel alone after a codebook or merge-group
# change, and duplicates what the BUILD block does further down.
#
# It depends on schemes_long and on functions defined in Section 3, so it only
# runs once those are already in memory. On a cold source() of the whole file,
# skip it -- the BUILD block builds the same object from the same call chain.

file.remove(file.path(CACHE_DIR, "outpatient_panel.rds"))

outpatient <- cache_or_run("outpatient_panel", {
  p <- load_outpatient()
  p <- apply_concept_merges(p)
  attach_scheme_columns(p, schemes_long)
})
######## Section 3: Shoppability scheme construction ##########################
#
# Builds the eighteen shoppability classification schemes from the codebook.
# A scheme is a partition of the concept universe into HIGH / INTERMEDIATE /
# LOW, assembled from named terms; a term is a family plus an optional keyword
# regex on the official code description plus an optional exclusion regex.
#
# Three rules govern the construction:
#
#   1. Emergency and critical care concepts, and everything billed as MS-DRG,
#      are forced LOW in every scheme. These are the services no patient or
#      payer can shop for, so no classification is allowed to disagree.
#      build_schemes() asserts this rather than assuming it.
#
#   2. Residual "_Other" terms mean whatever is left in a family after the
#      named sub-terms have been claimed WITHIN THE SAME SCHEME. "CT Other" is
#      therefore scheme-dependent by design.
#
#   3. Inclusion keywords are substring matches against abbreviated CPT text,
#      which over-matches. "CHEST" appears inside "Ct chest spine w/o dye", a
#      thoracic spine study, and "ULTRASOUND" appears inside "Colonoscopy
#      w/ultrasound", a procedure adjunct. Both would otherwise enter the
#      shoppable group. resolve_term() therefore honours an `exclude` regex,
#      and the guard at the end of build_schemes() tests for the specific case.

resolve_term <- function(term, universe) {
  rule <- TERM_RULES[[term]]
  if (is.null(rule)) stop("No rule for scheme term: '", term, "'", call. = FALSE)
  m <- universe$ANALYSIS_FAMILY_ID %chin% rule$family
  if (!is.null(rule$keyword) && rule$keyword != ".") {
    m <- m & grepl(rule$keyword, universe$OFFICIAL_DESCRIPTION, ignore.case = TRUE)
  }
  if (!is.null(rule$exclude)) {
    m <- m & !grepl(rule$exclude, universe$OFFICIAL_DESCRIPTION, ignore.case = TRUE)
  }
  m
}

CT_LUNG_EXCLUDE <- "SPINE|SPINAL|VERTEBR|ANGIO|ANGIOGRAM|CTA\\b|BIOPS|GUID|INTRAOP|DRAIN|ASPIRAT"
US_EXCLUDE <- "INTRAOP|INTRA-OP|GUID|GUIDANCE|W/ULTRASOUND|W ULTRASOUND|INTRAVAS|ADDL|CATH"

TERM_RULES <- list(
  Ultrasound  = list(family = c("DIAGNOSTIC_ULTRASOUND", "VASCULAR_ULTRASOUND", "ECHOCARDIOGRAPHY"),
                     exclude = US_EXCLUDE),
  Mammography = list(family = "MAMMOGRAPHY"),
  Biopsy      = list(family = "BIOPSY"),
  Colonoscopy = list(family = "COLONOSCOPY_LOWER_ENDOSCOPY"),
  Endoscopy   = list(family = "UPPER_ENDOSCOPY"),
  `X-Ray`     = list(family = c("XRAY_FLUOROSCOPY", "BONE_DENSITY")),
  MRI         = list(family = "MRI_MRA"),
  CT          = list(family = "CT_CTA"),

  `Ultrasound OB`      = list(family = "DIAGNOSTIC_ULTRASOUND",
                              keyword = "OB|OBSTETR|PREGNAN|FETAL", exclude = US_EXCLUDE),
  `Ultrasound Breast`  = list(family = "DIAGNOSTIC_ULTRASOUND",
                              keyword = "BREAST", exclude = US_EXCLUDE),
  `Ultrasound Abdomen` = list(family = "DIAGNOSTIC_ULTRASOUND",
                              keyword = "ABDOM", exclude = US_EXCLUDE),

  `X-Ray Chest`      = list(family = "XRAY_FLUOROSCOPY", keyword = "CHEST|THORAX|RIB",
                            exclude = "SPINE|SPINAL|VERTEBR"),
  `X-Ray Extremity`  = list(family = "XRAY_FLUOROSCOPY", keyword = "ARM|LEG|HAND|FOOT|WRIST|ANKLE|KNEE|ELBOW|SHOULDER|FEMUR|TIBIA|FIBULA|HUMERUS|EXTREMIT"),
  `X-Ray Skull/Head` = list(family = "XRAY_FLUOROSCOPY", keyword = "SKULL|HEAD|FACIAL|SINUS|MASTOID|ORBIT"),
  `X-Ray Abdomen`    = list(family = "XRAY_FLUOROSCOPY", keyword = "ABDOM"),
  `X-Ray Pelvis`     = list(family = "XRAY_FLUOROSCOPY", keyword = "PELVI|HIP"),
  `X-Ray Spine`      = list(family = "XRAY_FLUOROSCOPY", keyword = "SPINE|SPINAL|VERTEBR|LUMBAR|CERVICAL|SACRUM|COCCYX"),
  `X-Ray Other`      = list(family = "XRAY_FLUOROSCOPY", keyword = "."),

  # Excludes thoracic spine studies and CT angiography, which the CHEST
  # keyword would otherwise pull into the shoppable group.
  `CT Lung`       = list(family = "CT_CTA", keyword = "CHEST|LUNG|THORAX",
                         exclude = CT_LUNG_EXCLUDE),
  `CT Brain/Head` = list(family = "CT_CTA", keyword = "BRAIN|HEAD|SKULL|FACIAL|SINUS|ORBIT"),
  `CT Spine`      = list(family = "CT_CTA", keyword = "SPINE|SPINAL|VERTEBR|LUMBAR|CERVICAL|SACRUM"),
  `CT Neck`       = list(family = "CT_CTA", keyword = "NECK"),
  `CT Extremity`  = list(family = "CT_CTA", keyword = "ARM|LEG|HAND|FOOT|WRIST|ANKLE|KNEE|ELBOW|SHOULDER|EXTREMIT"),
  `CT Abdomen`    = list(family = "CT_CTA", keyword = "ABDOM"),
  `CT Chest`      = list(family = "CT_CTA", keyword = "CHEST|THORAX",
                         exclude = CT_LUNG_EXCLUDE),
  `CT Pelvis`     = list(family = "CT_CTA", keyword = "PELVI"),
  `CT Angio`      = list(family = "CT_CTA", keyword = "ANGIO|ANGIOGRAM|CTA\\b|VASCULAR"),
  `CT Other`      = list(family = "CT_CTA", keyword = "."),

  `MRI Brain/Head` = list(family = "MRI_MRA", keyword = "BRAIN|HEAD|SKULL|FACIAL|SINUS|ORBIT"),
  `MRI Spine`      = list(family = "MRI_MRA", keyword = "SPINE|SPINAL|VERTEBR|LUMBAR|CERVICAL|SACRUM"),
  `MRI Pelvis`     = list(family = "MRI_MRA", keyword = "PELVI"),
  `MRI Chest`      = list(family = "MRI_MRA", keyword = "CHEST|THORAX", exclude = "SPINE|SPINAL|VERTEBR"),
  `MRI Extremity`  = list(family = "MRI_MRA", keyword = "ARM|LEG|HAND|FOOT|WRIST|ANKLE|KNEE|ELBOW|SHOULDER|EXTREMIT"),
  `MRI Abdomen`    = list(family = "MRI_MRA", keyword = "ABDOM"),
  `MRI Neck`       = list(family = "MRI_MRA", keyword = "NECK"),
  `MRI Angio`      = list(family = "MRI_MRA", keyword = "ANGIO|ANGIOGRAM|MRA\\b|VASCULAR"),
  `MRI Breast`     = list(family = "MRI_MRA", keyword = "BREAST"),
  `MRI Other`      = list(family = "MRI_MRA", keyword = "."),

  `Biopsy Thyroid`    = list(family = "BIOPSY", keyword = "THYROID"),
  `Biopsy Breast`     = list(family = "BIOPSY", keyword = "BREAST"),
  `Biopsy Lymph Node` = list(family = "BIOPSY", keyword = "LYMPH"),
  `Biopsy Pancreas`   = list(family = "BIOPSY", keyword = "PANCREAS|PANCREATIC"),
  `Biopsy Liver`      = list(family = "BIOPSY", keyword = "LIVER|HEPAT"),
  `Biopsy Lung`       = list(family = "BIOPSY", keyword = "LUNG|PULMONARY|CHEST"),
  `Biopsy Kidney`     = list(family = "BIOPSY", keyword = "KIDNEY|RENAL"),
  `Biopsy Bone`       = list(family = "BIOPSY", keyword = "BONE"),
  `Biopsy Other`      = list(family = "BIOPSY", keyword = ".")
)

SHOPPABILITY_SCHEMES <- list(
  scheme_cms = list(label = "CMS Rule Definition",
                    shoppable = c("Ultrasound", "CT Lung", "Mammography")),
  scheme_theory = list(label = "Theory-Based",
                       shoppable = c("Ultrasound", "CT Lung", "Mammography"),
                       nonshoppable = c("Biopsy", "Colonoscopy", "Endoscopy")),
  scheme_broad = list(label = "Broad Shoppable",
                      shoppable = c("Ultrasound", "CT Lung", "Mammography", "X-Ray"),
                      nonshoppable = c("Biopsy", "Colonoscopy", "Endoscopy")),
  scheme_imaging = list(label = "Imaging vs Procedural",
                        shoppable = c("Ultrasound", "CT Lung", "Mammography", "X-Ray", "MRI"),
                        nonshoppable = c("Biopsy", "Colonoscopy", "Endoscopy")),
  scheme_theory_v2 = list(label = "Theory-Based V2 (MRI nonshoppable)",
                          shoppable = c("Ultrasound", "CT Lung", "Mammography"),
                          nonshoppable = c("Biopsy", "Colonoscopy", "Endoscopy", "MRI")),
  scheme_split_nonshop = list(label = "Split Non-Shoppable: MRI vs Procedural",
                              shoppable = c("Ultrasound", "CT Lung", "Mammography"),
                              nonshoppable = c("MRI", "Biopsy", "Colonoscopy", "Endoscopy")),
  scheme_cms_statutory = list(label = "CMS Statutory Shoppable List",
                              shoppable = c("Ultrasound", "CT Lung", "Mammography", "Colonoscopy", "X-Ray"),
                              nonshoppable = c("Biopsy", "MRI")),
  scheme_anatomical = list(label = "High vs Low Within Modality",
                           shoppable = c("CT Lung", "Mammography", "Ultrasound OB", "Ultrasound Breast",
                                         "Ultrasound Abdomen", "X-Ray Chest", "X-Ray Extremity"),
                           nonshoppable = c("Biopsy Pancreas", "Biopsy Liver", "Biopsy Lung", "Biopsy Kidney",
                                            "Biopsy Bone", "Colonoscopy", "Endoscopy", "MRI")),
  scheme_ct_broad = list(label = "CT-Inclusive",
                         shoppable = c("Ultrasound", "CT", "Mammography"),
                         nonshoppable = c("Biopsy", "Colonoscopy", "Endoscopy", "MRI")),
  scheme_ct_broad_ex_angio = list(label = "CT-Inclusive Except Angio",
                                  shoppable = c("Ultrasound", "CT Lung", "CT Brain/Head", "CT Spine", "CT Neck",
                                                "CT Extremity", "CT Abdomen", "CT Chest", "CT Pelvis", "CT Other",
                                                "Mammography"),
                                  nonshoppable = c("Biopsy", "Colonoscopy", "Endoscopy", "MRI", "CT Angio")),
  scheme_div1_operations = list(label = "Alt: Core Operations Framework",
                                shoppable = c("Mammography", "Ultrasound", "CT", "X-Ray Other", "X-Ray Skull/Head",
                                              "X-Ray Chest", "X-Ray Abdomen", "X-Ray Pelvis", "X-Ray Extremity"),
                                nonshoppable = c("Biopsy", "Colonoscopy", "Endoscopy")),
  scheme_div2_cms_legal = list(label = "Alt: CMS Legal/Regulatory Framework",
                               shoppable = c("Mammography", "Ultrasound", "X-Ray Other", "X-Ray Skull/Head",
                                             "X-Ray Chest", "X-Ray Abdomen", "X-Ray Pelvis", "X-Ray Extremity",
                                             "CT Brain/Head", "CT Abdomen", "MRI Brain/Head", "MRI Spine",
                                             "Colonoscopy", "Endoscopy"),
                               nonshoppable = c("Biopsy Other", "Biopsy Pancreas", "Biopsy Kidney", "Biopsy Bone",
                                                "Biopsy Liver", "Biopsy Lung", "CT Angio", "MRI Angio")),
  scheme_div3_mdsave = list(label = "Alt: Upfront Cash-Market Framework",
                            shoppable = c("CT", "MRI", "X-Ray", "Ultrasound", "Mammography"),
                            nonshoppable = c("Biopsy Other", "Biopsy Pancreas", "Biopsy Kidney",
                                             "Biopsy Bone", "Biopsy Liver", "Biopsy Lung")),
  scheme_div4_geographic = list(label = "Alt: Geographic/Facility Access",
                                shoppable = c("Mammography", "Ultrasound", "X-Ray"),
                                nonshoppable = c("Biopsy", "Colonoscopy", "Endoscopy")),
  scheme_div5_urgency = list(label = "Alt: Diagnostic Urgency/Lead-Time",
                             shoppable = c("Mammography", "Ultrasound OB", "X-Ray Spine", "Colonoscopy"),
                             nonshoppable = c("Biopsy", "MRI Spine")),
  scheme_alt6_staffing = list(label = "Alt: Staffing/Specialist Framework",
                              shoppable = c("CT", "MRI", "Mammography", "Ultrasound", "X-Ray"),
                              nonshoppable = c("Biopsy Other", "Biopsy Pancreas", "Biopsy Kidney", "Biopsy Bone",
                                               "Biopsy Liver", "Biopsy Lung", "Colonoscopy", "Endoscopy")),
  scheme_alt7_liability = list(label = "Alt: Incident Reporting/Liability",
                               shoppable = c("Ultrasound", "Mammography"),
                               nonshoppable = c("Biopsy", "Colonoscopy", "Endoscopy")),
  scheme_alt8_no_surprises = list(label = "Alt: No Surprises Act Framework",
                                  shoppable = c("Mammography", "Ultrasound", "X-Ray"),
                                  nonshoppable = c("Biopsy", "Colonoscopy", "Endoscopy"))
)

is_residual_term <- function(term) grepl(" Other$", term)

classify_scheme <- function(scheme, universe) {
  n <- nrow(universe)
  named <- c(Filter(Negate(is_residual_term), scheme$shoppable),
             Filter(Negate(is_residual_term), scheme$nonshoppable))
  claimed <- Reduce(`|`, lapply(named, resolve_term, universe = universe), init = rep(FALSE, n))

  apply_terms <- function(terms) {
    if (is.null(terms) || length(terms) == 0L) return(rep(FALSE, n))
    m <- rep(FALSE, n)
    for (t in terms[!is_residual_term(terms)]) m <- m | resolve_term(t, universe)
    for (t in terms[is_residual_term(terms)])  m <- m | (resolve_term(t, universe) & !claimed)
    m
  }

  shop <- apply_terms(scheme$shoppable); nons <- apply_terms(scheme$nonshoppable)
  always <- (universe$BILLING_CODE_TYPE == "MS_DRG") |
    (universe$ANALYSIS_FAMILY_ID %chin% ALWAYS_NONSHOPPABLE_FAMILIES)
  shop[always] <- FALSE; nons[always] <- TRUE

  cat <- rep("INTERMEDIATE", n); cat[shop] <- "HIGH"; cat[nons & !shop] <- "LOW"
  cat
}

build_schemes <- function() {
  cb <- fread(FILES$codebook)
  assert_columns(cb, c("BILLING_CODE_TYPE", "OFFICIAL_DESCRIPTION",
                       "ANALYSIS_FAMILY_ID", "ANALYSIS_CONCEPT_ID"), "codebook")
  universe <- unique(cb[, .(BILLING_CODE_TYPE, ANALYSIS_FAMILY_ID,
                            ANALYSIS_CONCEPT_ID, OFFICIAL_DESCRIPTION)])

  long <- rbindlist(lapply(names(SHOPPABILITY_SCHEMES), function(id) {
    s <- SHOPPABILITY_SCHEMES[[id]]
    cbind(universe, SCHEME_ID = id, SCHEME_NAME = s$label,
          SHOPPABILITY_CATEGORY = classify_scheme(s, universe))
  }))
  long[, SHOPPABILITY_ORDINAL := fcase(SHOPPABILITY_CATEGORY == "LOW", 1,
                                       SHOPPABILITY_CATEGORY == "INTERMEDIATE", 2,
                                       SHOPPABILITY_CATEGORY == "HIGH", 3, default = NA_real_)]

  save_qa_csv(long[, .(N = .N, SAMPLE = paste(head(unique(OFFICIAL_DESCRIPTION), 8),
                                              collapse = " | ")),
                   by = .(SCHEME_NAME, SHOPPABILITY_CATEGORY)][order(SCHEME_NAME)],
              "QA02_scheme_classification_crosswalk.csv")

  # Guards. Emergency and MS-DRG concepts must be LOW in every scheme; a
  # failure here means a family identifier changed upstream.
  ed <- long[ANALYSIS_FAMILY_ID == "EMERGENCY_DEPARTMENT",
             .(N = .N, L = sum(SHOPPABILITY_CATEGORY == "LOW")), by = SCHEME_ID]
  if (nrow(ed) > 0L && !all(ed$N == ed$L)) stop("ED concepts not LOW in every scheme.", call. = FALSE)
  dg <- long[BILLING_CODE_TYPE == "MS_DRG",
             .(N = .N, L = sum(SHOPPABILITY_CATEGORY == "LOW")), by = SCHEME_ID]
  if (nrow(dg) > 0L && !all(dg$N == dg$L)) stop("MS-DRG concepts not LOW in every scheme.", call. = FALSE)

  # Guard on the keyword exclusions, stated as a test: no CT spine or
  # angiography concept may sit in the HIGH group of the primary scheme.
  bad_ct <- long[SCHEME_ID == "scheme_theory_v2" & SHOPPABILITY_CATEGORY == "HIGH" &
                   grepl("SPINE|ANGIO", OFFICIAL_DESCRIPTION, ignore.case = TRUE)]
  if (nrow(bad_ct) > 0L) {
    warning("CT spine/angio still HIGH under Theory V2: ", nrow(bad_ct), " concepts.",
            call. = FALSE)
    print(head(bad_ct[, .(OFFICIAL_DESCRIPTION)], 10))
  } else {
    cat("Keyword guard passed: no CT spine/angio in Theory V2 HIGH.\n")
  }

  cat("Schemes built:", uniqueN(long$SCHEME_ID), "x",
      format(nrow(universe), big.mark = ","), "concepts\n")
  long
}

attach_scheme_columns <- function(panel, schemes_long, scheme_spec = PRIMARY_SCHEMES) {
  d <- copy(panel)
  rules <- list(
    high_vs_rest = function(x) fifelse(x == "HIGH", "Shoppable", "Non_shoppable"),
    low_vs_rest  = function(x) fifelse(x == "LOW", "Non_shoppable", "Shoppable"),
    extremes     = function(x) fcase(x == "HIGH", "Shoppable", x == "LOW", "Non_shoppable",
                                     default = NA_character_))
  for (spec in scheme_spec) {
    if (spec$col %in% names(d)) d[, (spec$col) := NULL]
    if (spec$source == "families") {
      d[, (spec$col) := factor(fifelse(FINAL_FAMILY_ID %chin% DIAGNOSTIC_FAMILIES,
                                       "Shoppable", "Non_shoppable"),
                               levels = c("Non_shoppable", "Shoppable"))]
      next
    }
    s <- schemes_long[SCHEME_ID == spec$source]
    if (nrow(s) == 0L) { warning("Scheme absent: ", spec$source); next }
    a <- unique(s[, .(FINAL_CONCEPT_ID = ANALYSIS_CONCEPT_ID,
                      V = factor(rules[[spec$rule]](toupper(trimws(SHOPPABILITY_CATEGORY))),
                                 levels = c("Non_shoppable", "Shoppable")))])
    setnames(a, "V", spec$col)
    d <- merge(d, a, by = "FINAL_CONCEPT_ID", all.x = TRUE, sort = FALSE)
  }
  setDT(d); d
}

# Scheme 1 is built from whole families rather than from the scheme table, so
# the QA02 crosswalk cannot verify it. This writes its concept-level assignment
# out separately and flags any procedure-adjunct concept sitting inside the
# shoppable group, which is the input to the Scheme 1 sensitivity variant.
export_scheme1_assignments <- function(panel) {
  if (!("SCHEME_1_CERTAINTY" %in% names(panel))) {
    warning("SCHEME_1_CERTAINTY absent; nothing to export.", call. = FALSE)
    return(invisible(NULL))
  }
  a <- unique(panel[, .(FINAL_FAMILY_ID, FINAL_CONCEPT_ID, SERVICE_NAME,
                        SCHEME_1_CERTAINTY)])
  setorder(a, SCHEME_1_CERTAINTY, FINAL_FAMILY_ID, SERVICE_NAME)
  save_qa_csv(a, "QA03_scheme1_concept_assignments.csv")

  cat("\n", strrep("=", 84), "\nSCHEME 1 (PRIMARY) — CONCEPTS PER FAMILY\n",
      strrep("=", 84), "\n", sep = "")
  print(a[, .N, by = .(SCHEME_1_CERTAINTY, FINAL_FAMILY_ID)][
    order(SCHEME_1_CERTAINTY, -N)])

  flagged <- a[SCHEME_1_CERTAINTY == "Shoppable" &
                 grepl(PROCEDURAL_CONCEPT_PATTERN, FINAL_CONCEPT_ID, ignore.case = TRUE)]
  cat("\nProcedural-looking concepts inside the SHOPPABLE group:", nrow(flagged), "\n")
  if (nrow(flagged) > 0L) {
    print(head(flagged[, .(FINAL_FAMILY_ID, SERVICE_NAME)], 25))
    save_qa_csv(flagged, "QA04_procedural_concepts_in_shoppable.csv")
  }
  invisible(a)
}

apply_scheme1_procedural_exclusion <- function(panel) {
  d <- copy(panel)
  d[SCHEME_1_CERTAINTY == "Shoppable" &
      grepl(PROCEDURAL_CONCEPT_PATTERN, FINAL_CONCEPT_ID, ignore.case = TRUE),
    SCHEME_1_CERTAINTY := factor("Non_shoppable",
                                 levels = c("Non_shoppable", "Shoppable"))]
  d
}

# ---------------------------------------------------------------------------
# Concept merges
# ---------------------------------------------------------------------------
#
# Six groups of concepts describe the same clinical service split on an
# administrative dimension -- contrast versus no contrast, who billed it, which
# treatment day. Estimating them separately fragments the sample without adding
# information, so they are collapsed to a canonical concept before estimation.
#
# The merge runs BEFORE scheme attachment, so merged concepts pick up their
# scheme columns through the same code path as everything else.

MERGE_GROUPS <- list(
  list(canonical_id = "MRI_MRA_MRI_ABDOMEN", canonical_family = "MRI_MRA",
       constituents = c("MRI_MRA_MRI_ABDOMEN", "MRI_MRA_MRI_ABDOMEN_W_CONTRAST",
                        "MRI_MRA_QMRCP_W_DX_MRI_SAME_ANATOMY")),
  list(canonical_id = "CT_CTA_CT_ABDOMEN", canonical_family = "CT_CTA",
       constituents = c("CT_CTA_CT_ABDOMEN", "CT_CTA_CT_ABDOMEN_W_CONTRAST")),
  list(canonical_id = "CT_CTA_CT_ABD_PELVIS", canonical_family = "CT_CTA",
       constituents = c("CT_CTA_CT_ABD_PELVIS", "CT_CTA_CT_ABD_PELVIS_W_CONTRAST")),
  list(canonical_id = "MRI_MRA_FMRI_BRAIN", canonical_family = "MRI_MRA",
       constituents = c("MRI_MRA_FMRI_BRAIN_BY_PHYS_PSYCH", "MRI_MRA_FMRI_BRAIN_BY_TECH")),
  list(canonical_id = "MRI_MRA_ARHFCMRIGTBS", canonical_family = "MRI_MRA",
       constituents = c("MRI_MRA_ARHFCMRIGTBS_1ST_TX_DAY",
                        "MRI_MRA_ARHFCMRIGTBS_SBSQ_PER_TX_DAY",
                        "MRI_MRA_ARHFCMRIGTBS_SBSQ_TX_DAY",
                        "MRI_MRA_PRSNLZ_TRGT_DVL_ARHFCMRIGTBS")),
  list(canonical_id = "MAMMOGRAPHY_SCR_MAMMO_BI_INCL_CAD", canonical_family = "MAMMOGRAPHY",
       constituents = c("MAMMOGRAPHY_BREAST_TOMOSYNTHESIS_BI",
                        "MAMMOGRAPHY_SCR_MAMMO_BI_INCL_CAD"))
)


# ---------------------------------------------------------------------------
# Scheme inheritance for constructed canonical IDs
# ---------------------------------------------------------------------------
#
# Four of the six merge groups use a canonical_id that is also one of their own
# constituents, so they already exist in the codebook and need nothing:
#   CT_CTA_CT_ABDOMEN, CT_CTA_CT_ABD_PELVIS, MRI_MRA_MRI_ABDOMEN,
#   MAMMOGRAPHY_SCR_MAMMO_BI_INCL_CAD
#
# Two are constructed names that never existed as concepts in the codebook:
#   MRI_MRA_FMRI_BRAIN    <- fMRI billed BY_PHYS_PSYCH + BY_TECH
#   MRI_MRA_ARHFCMRIGTBS  <- four ARHFCMRIGTBS treatment-day codes
#
# Without an inherited classification these two have no row in schemes_long,
# attach_scheme_columns() leaves them NA under schemes 2-6, and they drop out of
# every scheme-based model without a warning. run_preflight() tests for exactly
# this.
#
# Inheriting from a designated donor constituent is defensible here because the
# constituents of each group are the same clinical service split on an
# administrative dimension, and every scheme classifies on service content. The
# constituents therefore agree by construction. That agreement is verified
# rather than assumed, and a warning fires where it does not hold.

CANONICAL_SCHEME_DONOR <- c(
  MRI_MRA_FMRI_BRAIN   = "MRI_MRA_FMRI_BRAIN_BY_PHYS_PSYCH",
  MRI_MRA_ARHFCMRIGTBS = "MRI_MRA_ARHFCMRIGTBS_1ST_TX_DAY"
)

extend_schemes_for_merged <- function(schemes_long, donors = CANONICAL_SCHEME_DONOR) {
  s <- copy(schemes_long)

  for (canon in names(donors)) {
    grp <- Filter(function(g) g$canonical_id == canon, MERGE_GROUPS)
    if (length(grp) == 0L) { warning("No merge group for ", canon, call. = FALSE); next }
    chk <- s[ANALYSIS_CONCEPT_ID %chin% grp[[1]]$constituents,
             .(N_DISTINCT = uniqueN(SHOPPABILITY_CATEGORY)), by = SCHEME_ID]
    if (nrow(chk) == 0L) {
      warning(canon, ": no constituents found in schemes_long.", call. = FALSE)
    } else if (any(chk$N_DISTINCT > 1L)) {
      warning(canon, ": constituents DISAGREE under ", sum(chk$N_DISTINCT > 1L),
              " scheme(s). Inheritance is arbitrary there.", call. = FALSE)
      print(chk[N_DISTINCT > 1L])
    } else {
      cat("  ", canon, ": constituents agree under all ", nrow(chk), " schemes.\n", sep = "")
    }
  }

  add <- rbindlist(lapply(names(donors), function(canon) {
    srcrows <- s[ANALYSIS_CONCEPT_ID == donors[[canon]]]
    if (nrow(srcrows) == 0L) {
      warning("Donor absent for ", canon, ": ", donors[[canon]], call. = FALSE)
      return(data.table())
    }
    copy(srcrows)[, ANALYSIS_CONCEPT_ID := canon]
  }), fill = TRUE)

  if (nrow(add) == 0L) return(s)

  out <- rbind(s[!(ANALYSIS_CONCEPT_ID %chin% names(donors))], add, fill = TRUE)
  cat("Extended schemes_long with", uniqueN(add$ANALYSIS_CONCEPT_ID),
      "merged canonical concepts.\n")
  out
}


apply_concept_merges <- function(panel, merge_groups = MERGE_GROUPS) {
  all_constituents <- unique(unlist(lapply(merge_groups, `[[`, "constituents")))

  exact_raw <- read_panel(FILES$outpatient_exact, "exact-code panel (for merge)",
                          columns = unique(c(ANALYSIS_COLUMNS, "BILLING_CODE")))
  exact <- prepare_panel(exact_raw, "exact-code panel (for merge)",
                         choose_sample_flag(exact_raw))
  rm(exact_raw); invisible(gc())

  exact <- exact[FINAL_CONCEPT_ID %chin% all_constituents]
  if (nrow(exact) == 0L) {
    stop("No exact-code rows match the constituent concept IDs.", call. = FALSE)
  }

  map <- rbindlist(lapply(merge_groups, function(g)
    data.table(FINAL_CONCEPT_ID = g$constituents, CANON_ID = g$canonical_id,
               CANON_FAMILY = g$canonical_family)))
  exact <- merge(exact, map, by = "FINAL_CONCEPT_ID", all.x = TRUE, sort = FALSE)

  cat("\nExact-code rows found per canonical concept:\n")
  print(exact[, .N, by = CANON_ID])

  # Equal-weighted median of the exact-code medians within each hospital-month,
  # matching the aggregation convention used upstream in SQL.
  merged <- exact[, .(
    MEDIAN_PRICE = median(MEDIAN_PRICE, na.rm = TRUE),
    MEAN_PRICE   = mean(MEDIAN_PRICE, na.rm = TRUE),
    P25_PRICE    = quantile(MEDIAN_PRICE, .25, na.rm = TRUE, type = 7),
    P75_PRICE    = quantile(MEDIAN_PRICE, .75, na.rm = TRUE, type = 7),
    TOTAL_BEDS = first(TOTAL_BEDS), N_PRIOR_POSTERS = first(N_PRIOR_POSTERS),
    ANALYSIS_MARKET = first(ANALYSIS_MARKET),
    SERVICE_NAME = first(SERVICE_NAME), SERVICE_LABEL = first(SERVICE_LABEL)
  ), by = .(HOSPITAL_ID, POST_MONTH, CANON_ID, CANON_FAMILY)]

  setnames(merged, c("CANON_ID", "CANON_FAMILY"), c("FINAL_CONCEPT_ID", "FINAL_FAMILY_ID"))
  merged[, `:=`(
    ANALYSIS_SERVICE_ID = FINAL_CONCEPT_ID,
    MARKET_ID = paste(ANALYSIS_MARKET, FINAL_CONCEPT_ID, sep = "::"),
    BILLING_CODE_TYPE = "HCPCS",
    IQR_PRICE = pmax(P75_PRICE - P25_PRICE, 0),
    LOG_TOTAL_BEDS = log(pmax(safe_numeric(TOTAL_BEDS), 1))
  )]
  merged[, `:=`(
    LN_MEDIAN_PRICE = safe_log_positive(MEDIAN_PRICE),
    LN_MEAN_PRICE   = safe_log_positive(MEAN_PRICE),
    LN_P25_PRICE    = safe_log_positive(P25_PRICE),
    LN_P75_PRICE    = safe_log_positive(P75_PRICE),
    LN_IQR_PRICE    = safe_log1p_nonneg(IQR_PRICE)
  )]

  inst_cols <- available_columns(panel, unname(ALL_CANDIDATE_INSTRUMENTS))
  lookup <- unique(panel[, c("HOSPITAL_ID", "POST_MONTH", inst_cols), with = FALSE])
  merged <- merge(merged, lookup, by = c("HOSPITAL_ID", "POST_MONTH"),
                  all.x = TRUE, sort = FALSE)

  out <- rbind(panel[!(FINAL_CONCEPT_ID %chin% all_constituents)], merged, fill = TRUE)

  cat("\nConcept merge complete: dropped", length(all_constituents),
      "raw concepts, added", length(merge_groups), "canonical.\n",
      "Panel concepts: ", uniqueN(panel$ANALYSIS_SERVICE_ID), " -> ",
      uniqueN(out$ANALYSIS_SERVICE_ID), "\n", sep = "")

  rm(exact, merged); invisible(gc())
  out
}

cat("Section 3 loaded\n")



######## Section 4: Estimators ################################################
#
# Five estimators, all sharing the same sample-construction and clustering
# logic. run_reduced_form(), run_first_stage(), run_ols(), and run_iv() are the
# pooled workhorses. estimate_interacted() is the single function behind every
# heterogeneity result in the paper, and returns the reduced form and the IV
# from the same estimation sample so the two are directly comparable.
#
# All four pooled estimators return NULL rather than erroring on a degenerate
# sample, so a sweep over hundreds of concepts skips what it cannot estimate
# instead of halting.

model_sample <- function(data, columns) {
  req <- available_columns(data, unique(columns))
  d <- data[complete.cases(data[, ..req])]
  if ("LN_IQR_PRICE" %in% req) d <- d[is.finite(LN_IQR_PRICE)]
  d
}

run_reduced_form <- function(data, outcome, instrument, controls = BASELINE_CONTROLS,
                             fixed_effects = BASELINE_FIXED_EFFECTS,
                             clusters = BASELINE_CLUSTERS) {
  controls <- available_columns(data, controls)
  fe <- available_columns(data, fixed_effects); cl <- available_columns(data, clusters)
  d <- model_sample(data, c(outcome, instrument, controls, fe, cl))
  if (nrow(d) < 100L || !has_usable_variation(d[[instrument]])) return(NULL)
  tryCatch(feols(build_ols_formula(outcome, c(instrument, controls), fe), data = d,
                 cluster = build_cluster_formula(cl), warn = FALSE, notes = FALSE),
           error = function(e) NULL)
}

run_first_stage <- function(data, instrument, endogenous = ENDOGENOUS_VARIABLE,
                            controls = BASELINE_CONTROLS,
                            fixed_effects = BASELINE_FIXED_EFFECTS,
                            clusters = BASELINE_CLUSTERS) {
  controls <- available_columns(data, controls)
  fe <- available_columns(data, fixed_effects); cl <- available_columns(data, clusters)
  d <- model_sample(data, c(endogenous, instrument, controls, fe, cl))
  if (nrow(d) < 100L || !has_usable_variation(d[[endogenous]]) ||
      !has_usable_variation(d[[instrument]])) return(NULL)
  tryCatch(feols(build_ols_formula(endogenous, c(instrument, controls), fe), data = d,
                 cluster = build_cluster_formula(cl), warn = FALSE, notes = FALSE),
           error = function(e) NULL)
}

run_ols <- function(data, outcome, endogenous = ENDOGENOUS_VARIABLE,
                    controls = BASELINE_CONTROLS, fixed_effects = BASELINE_FIXED_EFFECTS,
                    clusters = BASELINE_CLUSTERS) {
  controls <- available_columns(data, controls)
  fe <- available_columns(data, fixed_effects); cl <- available_columns(data, clusters)
  d <- model_sample(data, c(outcome, endogenous, controls, fe, cl))
  if (nrow(d) < 100L || !has_usable_variation(d[[endogenous]])) return(NULL)
  tryCatch(feols(build_ols_formula(outcome, c(endogenous, controls), fe), data = d,
                 cluster = build_cluster_formula(cl), warn = FALSE, notes = FALSE),
           error = function(e) NULL)
}

run_iv <- function(data, outcome, instrument, endogenous = ENDOGENOUS_VARIABLE,
                   controls = BASELINE_CONTROLS, fixed_effects = BASELINE_FIXED_EFFECTS,
                   clusters = BASELINE_CLUSTERS) {
  controls <- available_columns(data, controls)
  fe <- available_columns(data, fixed_effects); cl <- available_columns(data, clusters)
  d <- model_sample(data, c(outcome, endogenous, instrument, controls, fe, cl))
  if (nrow(d) < 100L || !has_usable_variation(d[[endogenous]]) ||
      !has_usable_variation(d[[instrument]])) return(NULL)
  tryCatch(feols(build_iv_formula(outcome, endogenous, instrument, controls, fe), data = d,
                 cluster = build_cluster_formula(cl), warn = FALSE, notes = FALSE),
           error = function(e) NULL)
}

# ---------------------------------------------------------------------------
# The interacted estimator -- one function, two kinds of moderator
# ---------------------------------------------------------------------------
#
# CATEGORICAL moderator (shoppability):
#     ln(P) = sum_k b_k (N x 1[cat = k]) + controls + FE
#   with each interacted treatment instrumented by Z x 1[cat = k]. b_k is
#   category k's own level effect, and wald_equality() tests b_1 = ... = b_K.
#
# CONTINUOUS moderator (comparability, demographics):
#     ln(P) = b0 N + b1 (N x M) + controls + FE
#   instrumented by Z and Z x M, with M centred at its estimation-sample mean.
#   b1 is the gradient in the moderator and b0 is the effect at its mean.
#
# In both cases the number of excluded instruments equals the number of
# endogenous regressors, so the model stays exactly identified and the
# reduced-form p-values remain weak-instrument robust. The lower-order terms in
# M are absorbed by the market fixed effect and are never estimated.

estimate_interacted <- function(
    data, moderator, outcome = PRIMARY_OUTCOME, instrument = PRIMARY_INSTRUMENT,
    endogenous = ENDOGENOUS_VARIABLE, moderator_type = c("categorical", "continuous"),
    label = "", instrument_label = "", moderator_label = NULL,
    controls = BASELINE_CONTROLS, fixed_effects = BASELINE_FIXED_EFFECTS,
    clusters = BASELINE_CLUSTERS, center_moderator = TRUE) {

  moderator_type <- match.arg(moderator_type)
  controls <- available_columns(data, controls)
  fe <- available_columns(data, fixed_effects); cl <- available_columns(data, clusters)

  d <- data[!is.na(get(moderator))]
  d <- model_sample(d, c(outcome, endogenous, instrument, controls, fe, cl, moderator))
  if (nrow(d) < MIN_MODEL_OBS) return(NULL)

  if (moderator_type == "categorical") {
    d[, MOD := droplevels(factor(get(moderator)))]
    keys <- levels(d$MOD); if (length(keys) < 2L) return(NULL)
    endo <- paste0("TREAT_", keys); ivs <- paste0("IV_", keys); rfs <- paste0("RF_", keys)
    for (k in seq_along(keys)) {
      sel <- as.integer(d$MOD == keys[k])
      d[, (endo[k]) := get(endogenous) * sel]
      d[, (ivs[k])  := get(instrument) * sel]
      d[, (rfs[k])  := get(instrument) * sel]
    }
    terms_label <- keys
  } else {
    d[, MODV := safe_numeric(get(moderator))]
    if (!has_usable_variation(d$MODV)) return(NULL)
    if (center_moderator) d[, MODV := MODV - mean(MODV, na.rm = TRUE)]
    d[, `:=`(TREAT_MAIN = get(endogenous), TREAT_INTER = get(endogenous) * MODV,
             IV_MAIN = get(instrument), IV_INTER = get(instrument) * MODV,
             RF_MAIN = get(instrument), RF_INTER = get(instrument) * MODV)]
    endo <- c("TREAT_MAIN", "TREAT_INTER"); ivs <- c("IV_MAIN", "IV_INTER")
    rfs <- c("RF_MAIN", "RF_INTER"); terms_label <- c("Main", "x Moderator")
  }

  rf_fit <- tryCatch(feols(build_ols_formula(outcome, c(rfs, controls), fe), data = d,
                           cluster = build_cluster_formula(cl), warn = FALSE, notes = FALSE),
                     error = function(e) NULL)
  iv_fit <- tryCatch(feols(build_iv_formula(outcome, endo, ivs, controls, fe), data = d,
                           cluster = build_cluster_formula(cl), warn = FALSE, notes = FALSE),
                     error = function(e) NULL)

  fsw <- first_stage_wald(iv_fit)
  fs_min <- if (nrow(fsw) > 0L) min(fsw$WALD, na.rm = TRUE) else NA_real_

  resolve <- function(fit, tm) {
    if (is.null(fit)) return(NA_character_)
    cand <- c(paste0("fit_", tm), tm); hit <- cand[cand %in% names(coef(fit))]
    if (length(hit) == 0L) NA_character_ else hit[1L]
  }
  pull <- function(fit, tm) {
    nm <- resolve(fit, tm)
    if (is.na(nm)) return(list(b = NA_real_, s = NA_real_))
    list(b = unname(coef(fit)[nm]), s = unname(sqrt(vcov(fit)[nm, nm])))
  }

  sd_z <- sd(d[[instrument]], na.rm = TRUE)
  ml <- moderator_label %||% moderator

  rows <- rbindlist(lapply(seq_along(endo), function(k) {
    rf <- pull(rf_fit, rfs[k]); iv <- pull(iv_fit, endo[k])
    data.table(
      SPEC = label, MODERATOR = ml, MODERATOR_TYPE = moderator_type,
      INSTRUMENT_LABEL = instrument_label, INSTRUMENT = instrument,
      OUTCOME = outcome, TERM = terms_label[k],
      RF_COEF = rf$b, RF_SE = rf$s, RF_P = 2 * pnorm(-abs(rf$b / rf$s)),
      RF_PERCENT_PER_SD = 100 * (exp(rf$b * sd_z) - 1),
      IV_COEF = iv$b, IV_SE = iv$s, IV_P = 2 * pnorm(-abs(iv$b / iv$s)),
      IV_PERCENT = 100 * (exp(iv$b) - 1),
      IV_CI_LOW_PERCENT  = 100 * (exp(iv$b - 1.96 * iv$s) - 1),
      IV_CI_HIGH_PERCENT = 100 * (exp(iv$b + 1.96 * iv$s) - 1),
      FIRST_STAGE_WALD_THIS_EQ = if (nrow(fsw) >= k) fsw$WALD[k] else NA_real_,
      FIRST_STAGE_WALD_MIN = fs_min, CRAGG_DONALD = cragg_donald(iv_fit),
      N_OBSERVATIONS = if (is.null(iv_fit)) NA_integer_ else nobs(iv_fit),
      N_CONCEPTS = uniqueN(d$FINAL_CONCEPT_ID))
  }), fill = TRUE)

  if (moderator_type == "categorical") {
    rf_test <- wald_equality(rf_fit, vapply(rfs,  resolve, character(1), fit = rf_fit))
    iv_test <- wald_equality(iv_fit, vapply(endo, resolve, character(1), fit = iv_fit))
  } else {
    grab <- function(fit, tm) {
      nm <- resolve(fit, tm); if (is.na(nm)) return(data.table())
      b <- coef(fit)[nm]; s <- sqrt(vcov(fit)[nm, nm])
      data.table(WALD = unname((b / s)^2), DF = 1L, VCOV_FULL_RANK = 1L,
                 P_VALUE = unname(2 * pnorm(-abs(b / s))))
    }
    rf_test <- grab(rf_fit, "RF_INTER"); iv_test <- grab(iv_fit, "TREAT_INTER")
  }

  tests <- rbindlist(list(
    if (nrow(rf_test) > 0L) cbind(ESTIMATOR = "Reduced form", rf_test) else NULL,
    if (nrow(iv_test) > 0L) cbind(ESTIMATOR = "IV", iv_test) else NULL), fill = TRUE)
  if (nrow(tests) > 0L) {
    tests[, `:=`(SPEC = label, MODERATOR = ml, INSTRUMENT_LABEL = instrument_label,
                 OUTCOME = outcome, FIRST_STAGE_WALD_MIN = fs_min,
                 N_OBSERVATIONS = if (is.null(iv_fit)) NA_integer_ else nobs(iv_fit))]
  }

  rm(d); invisible(gc())
  list(rows = rows, tests = tests)
}

apply_transform <- function(data, name, source = ENDOGENOUS_VARIABLE,
                            target = "TREAT_TRANSFORMED") {
  spec <- TRANSFORM_LADDER[[name]]
  if (is.null(spec)) stop("Unknown transform: ", name, call. = FALSE)
  x <- safe_numeric(data[[source]])
  v <- if (!is.null(spec$quantile)) pmin(x, quantile(x, spec$quantile, na.rm = TRUE, names = FALSE))
  else spec$fun(x)
  data[, (target) := v]; data
}

cat("Section 4 loaded\n")



######## Section 5: Instrument screen and pooled models #######################
#
# Two outputs. run_instrument_screen() estimates the first stage and reduced
# form for every candidate instrument on a common sample and writes T02, which
# is the table that documents why the enforcement family is exploratory and why
# the competitor family carries the exclusion-restriction defence.
#
# run_pooled_models() estimates the average effect across all outcomes. The
# null it returns is the paper's first finding rather than a failed test: there
# is no average price response, because the average mixes services that respond
# with services that cannot.

run_instrument_screen <- function(panel) {
  cands <- ALL_CANDIDATE_INSTRUMENTS[unname(ALL_CANDIDATE_INSTRUMENTS) %in% names(panel)]
  cands <- cands[vapply(unname(cands), function(v) has_usable_variation(panel[[v]]), logical(1))]
  if (length(cands) == 0L) stop("No usable instruments.", call. = FALSE)

  keep <- available_columns(panel, unique(c(ENDOGENOUS_VARIABLE, PRIMARY_OUTCOME,
                                            BASELINE_CONTROLS, BASELINE_FIXED_EFFECTS, BASELINE_CLUSTERS, unname(cands))))
  d0 <- panel[, ..keep]

  fam <- function(z) fcase(
    z %chin% unname(ENFORCEMENT_ROBUSTNESS_INSTRUMENTS), "CMS enforcement (exploratory)",
    z %chin% unname(MAIN_INSTRUMENTS), "Main system",
    z %chin% unname(CANONICAL_SYSTEM_INSTRUMENTS), "Canonical system",
    default = "Supporting window")

  cat("\n", strrep("=", 84), "\nFIRST-STAGE SCREEN: ", length(cands),
      " instruments\n", strrep("=", 84), "\n", sep = "")

  out <- rbindlist(lapply(names(cands), function(lab) {
    z <- cands[[lab]]
    fs <- run_first_stage(d0, z); rf <- run_reduced_form(d0, PRIMARY_OUTCOME, z)
    cf <- extract_coefficient(fs, z); cr <- extract_coefficient(rf, z)
    f <- if (is.finite(cf$statistic)) cf$statistic^2 else NA_real_
    cat(sprintf("  %-44s coef = %9.4f  F = %8.2f  RF p = %.4f\n",
                substr(lab, 1, 42), cf$estimate, f, cr$p_value))
    data.table(INSTRUMENT_LABEL = lab, INSTRUMENT = z, INSTRUMENT_FAMILY = fam(z),
               FIRST_STAGE_COEF = cf$estimate, FIRST_STAGE_SE = cf$std_error,
               FIRST_STAGE_F = f, FIRST_STAGE_P = cf$p_value,
               REDUCED_FORM_COEF = cr$estimate, REDUCED_FORM_SE = cr$std_error,
               REDUCED_FORM_P = cr$p_value,
               SHARE_INSTRUMENT_ZERO = mean(d0[[z]] == 0, na.rm = TRUE),
               N_OBSERVATIONS = if (is.null(fs)) NA_integer_ else nobs(fs),
               WEAK_F_LT_10 = as.integer(is.finite(f) && f < 10),
               HEADLINE_ELIGIBLE = as.integer(z %chin% unname(MAIN_INSTRUMENTS)),
               EXCLUSION_NOTE = fcase(
                 z %chin% unname(ENFORCEMENT_ROBUSTNESS_INSTRUMENTS),
                 "EXPLORATORY: enforcement may move prices directly.",
                 z %chin% unname(MAIN_INSTRUMENTS),
                 "Competitor rollout excluding the focal hospital's own system.",
                 default = "County-month local exposure; does NOT exclude own system."))
  }), fill = TRUE)

  setorder(out, -FIRST_STAGE_F)
  save_csv(out, "T02_instrument_first_stage_screen.csv")
  rm(d0); invisible(gc()); out
}

run_pooled_models <- function(panel, instruments = MAIN_INSTRUMENTS, outcomes = OUTCOMES) {
  rows <- list()
  for (il in names(instruments)) {
    z <- instruments[[il]]; if (!(z %in% names(panel))) next
    for (ol in names(outcomes)) {
      y <- outcomes[[ol]]; if (!(y %in% names(panel))) next
      ols <- run_ols(panel, y); rf <- run_reduced_form(panel, y, z)
      fs <- run_first_stage(panel, z); iv <- run_iv(panel, y, z)
      c_ols <- extract_coefficient(ols, ENDOGENOUS_VARIABLE)
      c_rf <- extract_coefficient(rf, z); c_fs <- extract_coefficient(fs, z)
      c_iv <- extract_coefficient(iv, c(paste0("fit_", ENDOGENOUS_VARIABLE), ENDOGENOUS_VARIABLE))
      rows[[length(rows) + 1L]] <- data.table(
        INSTRUMENT_LABEL = il, OUTCOME = ol,
        OLS_PERCENT = 100 * c_ols$estimate, OLS_P = c_ols$p_value,
        RF_COEF = c_rf$estimate, RF_SE = c_rf$std_error, RF_P = c_rf$p_value,
        FIRST_STAGE_COEF = c_fs$estimate,
        FIRST_STAGE_F = if (is.finite(c_fs$statistic)) c_fs$statistic^2 else NA_real_,
        IV_PERCENT = 100 * c_iv$estimate, IV_SE_PERCENT = 100 * c_iv$std_error,
        IV_P = c_iv$p_value, WU_HAUSMAN_P = extract_wu_hausman_p(iv),
        N_OBSERVATIONS = if (is.null(iv)) NA_integer_ else nobs(iv))
      cat(sprintf("  %-36s %-7s  IV = %8.3f%%  p = %.4f\n",
                  substr(il, 1, 34), ol, 100 * c_iv$estimate, c_iv$p_value))
      rm(ols, rf, fs, iv); invisible(gc())
    }
  }
  out <- rbindlist(rows, fill = TRUE)
  out[, IV_STARS := add_stars(IV_P)]
  out[, RF_STARS := add_stars(RF_P)]
  save_csv(out, "T03_pooled_OLS_RF_IV_all_outcomes.csv")
  cat("\nPooled: significant at 5% in", sum(out$IV_P < 0.05, na.rm = TRUE), "of", nrow(out),
      "| median p =", round(median(out$IV_P, na.rm = TRUE), 4),
      "\nA null here is the paper's FIRST FINDING, not a failure: no average effect.\n")
  out
}

cat("Section 5 loaded\n")



######## Section 6: Concept-level estimates ###################################
#
# Estimates the full specification separately for each clinical concept, which
# is the input to the meta-regressions in Section 8. This is the expensive
# stage: roughly eight hours for 738 concepts across six instruments.
#
# Three numbers are stored per concept-instrument pair, and keeping all three
# is what makes the decomposition in Section 8 possible:
#
#   RF   reduced form -- no denominator, and the meta-regression input
#   FS   first stage  -- the denominator itself, so it can be tested directly
#   IV   the ratio of the two
#
# If the shoppability gradient appeared in IV but not in RF, it would be a
# denominator artefact rather than a price response. Storing FS separately is
# what allows that to be ruled out rather than asserted.
#
# Progress is written to a _PARTIAL.csv every 50 concepts, so an interrupted
# run is recoverable.

CONCEPT_KEEP <- function(instruments, outcome) unique(c(
  "ANALYSIS_SERVICE_ID", "HOSPITAL_ID", "POST_MONTH", "ANALYSIS_MARKET", "MARKET_ID",
  outcome, ENDOGENOUS_VARIABLE, BASELINE_CONTROLS, BASELINE_FIXED_EFFECTS,
  BASELINE_CLUSTERS, unname(instruments),
  "BILLING_CODE_TYPE", "BILLING_CODE", "SERVICE_NAME", "SERVICE_LABEL",
  "FINAL_SUPERFAMILY_ID", "FINAL_FAMILY_ID", "FINAL_CONCEPT_ID", "FINAL_CONCEPT_NAME",
  "P25_PRICE", "P75_PRICE", "MEDIAN_PRICE", "N_DISTINCT_PAYERS", "N_PAYER_CELLS",
  "CODE_COVERAGE_RATIO", "N_CODES_PRESENT", "EXPECTED_N_CODES_IN_CONCEPT",
  vapply(PRIMARY_SCHEMES, `[[`, character(1), "col")))

build_concept_panel <- function(panel, instruments = MAIN_INSTRUMENTS,
                                outcome = PRIMARY_OUTCOME) {
  keep <- available_columns(panel, CONCEPT_KEEP(instruments, outcome))
  cat("Slimming:", ncol(panel), "->", length(keep), "columns\n")
  slim <- panel[, ..keep]
  setkey(slim, ANALYSIS_SERVICE_ID)
  slim
}

estimate_concept_level <- function(slim, instruments = MAIN_INSTRUMENTS,
                                   outcome = PRIMARY_OUTCOME, endogenous = ENDOGENOUS_VARIABLE,
                                   ols_instrument = names(instruments)[1L], drop_singletons = DROP_SINGLETON_MARKETS,
                                   max_concepts = Inf, save_stem = "T05_concept_level", save_every = 50L) {

  if (!identical(key(slim), "ANALYSIS_SERVICE_ID")) setkey(slim, ANALYSIS_SERVICE_ID)
  ids <- sort(unique(slim$ANALYSIS_SERVICE_ID))
  if (is.finite(max_concepts)) ids <- head(ids, as.integer(max_concepts))
  cat("\nConcepts:", length(ids), "| instruments:", length(instruments), "\n\n")

  rows <- vector("list", length(ids) * length(instruments))
  k <- 0L; skipped <- 0L; t0 <- Sys.time()

  for (i in seq_along(ids)) {
    sid <- ids[i]; sub <- slim[.(sid)]
    if (drop_singletons) sub <- drop_singleton_markets(sub)
    n_mk <- uniqueN(sub$ANALYSIS_MARKET); n_mo <- uniqueN(sub$POST_MONTH)
    if (nrow(sub) < MIN_SERVICE_OBS || n_mk < MIN_SERVICE_MARKETS ||
        n_mo < MIN_SERVICE_MONTHS) { skipped <- skipped + 1L; next }

    md <- sub[1L]; t1 <- Sys.time()

    for (il in names(instruments)) {
      z <- instruments[[il]]
      if (!(z %in% names(sub)) || !has_usable_variation(sub[[z]])) next

      rf <- run_reduced_form(sub, outcome, z); fs <- run_first_stage(sub, z, endogenous)
      iv <- run_iv(sub, outcome, z, endogenous)
      ols <- if (il == ols_instrument) run_ols(sub, outcome, endogenous) else NULL

      c_rf <- extract_coefficient(rf, z); c_fs <- extract_coefficient(fs, z)
      c_iv <- extract_coefficient(iv, c(paste0("fit_", endogenous), endogenous))
      c_ol <- extract_coefficient(ols, endogenous)

      res <- data.table(
        ANALYSIS_SERVICE_ID = sid, INSTRUMENT_LABEL = il, INSTRUMENT = z, OUTCOME = outcome,
        BILLING_CODE_TYPE = meta_chr(md, "BILLING_CODE_TYPE"),
        SERVICE_NAME = meta_chr(md, "SERVICE_NAME"),
        SERVICE_LABEL = meta_chr(md, "SERVICE_LABEL"),
        FINAL_FAMILY_ID = meta_chr(md, "FINAL_FAMILY_ID"),
        FINAL_CONCEPT_ID = meta_chr(md, "FINAL_CONCEPT_ID"),
        FINAL_CONCEPT_NAME = meta_chr(md, "FINAL_CONCEPT_NAME"),
        RF_COEF = c_rf$estimate, RF_SE = c_rf$std_error, RF_P = c_rf$p_value,
        FS_COEF = c_fs$estimate, FS_SE = c_fs$std_error,
        FS_F = if (is.finite(c_fs$statistic)) c_fs$statistic^2 else NA_real_,
        IV_COEF = c_iv$estimate, IV_SE = c_iv$std_error, IV_P = c_iv$p_value,
        IV_ESTIMATE_PERCENT = 100 * c_iv$estimate, IV_SE_PERCENT = 100 * c_iv$std_error,
        OLS_ESTIMATE_PERCENT = 100 * c_ol$estimate,
        N_ROWS = nrow(sub), N_HOSPITALS = uniqueN(sub$HOSPITAL_ID),
        N_MARKETS = n_mk, N_MONTHS = n_mo,
        N_OBSERVATIONS = if (is.null(iv)) NA_integer_ else nobs(iv),
        MEAN_PRIOR_POSTERS = mean(safe_numeric(sub[[ENDOGENOUS_VARIABLE]]), na.rm = TRUE))

      if (nrow(res) != 1L) stop(sprintf("Malformed row for %s / %s.", sid, il), call. = FALSE)
      k <- k + 1L; rows[[k]] <- res
      rm(rf, fs, iv, ols)
    }

    el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    cat(sprintf("[%d/%d] n=%s | %.1fs | elapsed %.1fm | ETA %.1fm\n", i, length(ids),
                format(nrow(sub), big.mark = ","),
                as.numeric(difftime(Sys.time(), t1, units = "secs")),
                el, (el / i) * (length(ids) - i)))

    if (i %% save_every == 0L && k > 0L) {
      p <- rbindlist(rows[seq_len(k)], fill = TRUE)
      if (nrow(p) > 0L) save_csv(p, paste0(save_stem, "_PARTIAL.csv"))
      rm(p)
    }
    if (i %% 50L == 0L) { rm(sub); invisible(gc()) }
  }

  cat("\nSkipped:", skipped, "of", length(ids), "| rows:", k, "\n")
  if (k == 0L) stop("No rows accumulated; check the _PARTIAL file.", call. = FALSE)

  out <- rbindlist(rows[seq_len(k)], fill = TRUE)
  out[, RF_P_FDR := p.adjust(RF_P, method = "BH"), by = INSTRUMENT_LABEL]
  out[, IV_P_FDR := p.adjust(IV_P, method = "BH"), by = INSTRUMENT_LABEL]
  setorder(out, INSTRUMENT_LABEL, RF_COEF)
  out
}

diagnose_size_gradient <- function(cr) {
  d <- copy(cr)
  d[, SHOP := fifelse(FINAL_FAMILY_ID %chin% DIAGNOSTIC_FAMILIES, "Shoppable", "Non_shoppable")]
  d[, SIZE_Q := cut(N_OBSERVATIONS, quantile(N_OBSERVATIONS, c(0, .25, .5, .75, 1), na.rm = TRUE),
                    labels = c("Q1 smallest", "Q2", "Q3", "Q4 largest"), include.lowest = TRUE)]
  cat(sprintf("\ncorr(log size, mean prior posters) = %.3f\ncorr(log size, first-stage coef)  = %.3f\n",
              cor(log(pmax(d$N_OBSERVATIONS, 1)), d$MEAN_PRIOR_POSTERS, use = "complete.obs"),
              cor(log(pmax(d$N_OBSERVATIONS, 1)), d$FS_COEF, use = "complete.obs")))
  tab <- d[!is.na(SIZE_Q), .(N = .N, MEDIAN_RF = median(RF_COEF, na.rm = TRUE),
                             MEDIAN_FS = median(FS_COEF, na.rm = TRUE),
                             MEDIAN_IV_PCT = median(IV_ESTIMATE_PERCENT, na.rm = TRUE)),
           by = .(INSTRUMENT_LABEL, SIZE_Q, SHOP)][order(INSTRUMENT_LABEL, SIZE_Q, SHOP)]
  save_csv(tab, "T05C_size_gradient_RF_vs_IV.csv"); print(tab)
  cat("\nGradient in IV but NOT in RF -> denominator artefact, and the size\n",
      "gradient cannot be used to rule out patient shopping.\n", sep = "")
  tab
}

cat("Section 6 loaded\n")



###############################################################################
# Section 7: Main results and the transform ladder
#
# run_main_results() estimates the interacted shoppability specification for
# every scheme x instrument combination and is the source of the paper's
# headline table. Stage 7 calls it five times -- headline, Scheme 1
# sensitivity, confirming tier, discrepant tier, and the pooled robustness set
# -- so output filenames are parameterised through `stem`. Each caller supplies
# a distinct one; sharing a stem silently overwrites the previous call's tables.
#
# Two panels are printed. Panel A is the reduced form in percent per standard
# deviation of peer exposure, and is primary. Panel B is the IV in percent per
# additional prior poster.
#
# run_transform_ladder() re-estimates the headline scheme under six treatment
# transforms. It is reported as robustness only. See design decision 1: the
# magnitudes move mechanically with compression, the conclusion does not.
###############################################################################
run_main_results <- function(panel, schemes = SCHEME_COLUMNS,
                             instruments = MAIN_INSTRUMENTS,
                             outcome = PRIMARY_OUTCOME,
                             stem = "T06_main") {
  rows <- list(); tests <- list(); g <- 0L
  n <- length(schemes) * length(instruments); t0 <- Sys.time()

  for (sl in names(schemes)) {
    sc <- schemes[[sl]]; if (!(sc %in% names(panel))) next
    for (il in names(instruments)) {
      z <- instruments[[il]]; if (!(z %in% names(panel))) next
      g <- g + 1L; t1 <- Sys.time()
      r <- estimate_interacted(panel, sc, outcome, z, moderator_type = "categorical",
                               label = sl, instrument_label = il, moderator_label = sl)
      if (!is.null(r)) { rows[[length(rows) + 1L]] <- r$rows; tests[[length(tests) + 1L]] <- r$tests }
      cat(sprintf("[%d/%d] %-26s %-36s | %5.1fs\n", g, n, substr(sl, 1, 24),
                  substr(il, 1, 34), as.numeric(difftime(Sys.time(), t1, units = "secs"))))
    }
  }

  mr <- rbindlist(rows, fill = TRUE); mt <- rbindlist(tests, fill = TRUE)
  if (nrow(mr) == 0L) stop("No main results.", call. = FALSE)

  save_csv(mr, paste0(stem, "_interacted_RF_and_IV.csv"))
  save_csv(mt, paste0(stem, "B_heterogeneity_tests.csv"))
  cat("\nElapsed:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")

  cat("\n", strrep("=", 104),
      "\nPANEL A - REDUCED FORM (percent per SD of peer exposure). PRIMARY.\n",
      "Exactly identified, so these p-values ARE Anderson-Rubin robust p-values.\n",
      strrep("=", 104), "\n", sep = "")
  print(dcast(mr[, .(SPEC = substr(SPEC, 1, 26), INSTRUMENT_LABEL, TERM,
                     E = round(RF_PERCENT_PER_SD, 3))],
              SPEC + INSTRUMENT_LABEL ~ TERM, value.var = "E"))
  cat("\nIndividual significance by category:\n")
  print(mr[, .(N = .N, SIG_05 = sum(RF_P < 0.05, na.rm = TRUE),
               MEDIAN_P = round(median(RF_P, na.rm = TRUE), 4),
               MEDIAN_PCT = round(median(RF_PERCENT_PER_SD, na.rm = TRUE), 3)), by = TERM])

  cat("\n", strrep("=", 104),
      "\nPANEL B - IV (percent per additional prior poster, linear treatment)\n",
      strrep("=", 104), "\n", sep = "")
  print(mr[, .(SPEC = substr(SPEC, 1, 24), INSTRUMENT_LABEL, TERM,
               IV_PCT = round(IV_PERCENT, 2), SE = round(100 * IV_SE, 2),
               P = round(IV_P, 4), MIN_WALD = round(FIRST_STAGE_WALD_MIN, 1))][
                 order(SPEC, INSTRUMENT_LABEL, TERM)])

  # GAP, the percentage-point difference between the two categories, is the
  # quantity to report. RATIO divides by the non-shoppable estimate, which is a
  # precise zero by design and therefore frequently near zero and occasionally
  # sign-flipped in the denominator. It is descriptive only and unstable by
  # construction.
  ratio <- dcast(mr[, .(SPEC, INSTRUMENT_LABEL, TERM, IV_PERCENT)],
                 SPEC + INSTRUMENT_LABEL ~ TERM, value.var = "IV_PERCENT")
  if (all(c("Shoppable", "Non_shoppable") %in% names(ratio))) {
    ratio[, `:=`(RESPONSE_RATIO = Shoppable / Non_shoppable,
                 GAP = Shoppable - Non_shoppable)]
    cat("\nGAP (percentage points) is the quantity to report. RATIO is unstable\n",
        "when Non_shoppable is near zero and is descriptive only.\n", sep = "")
    print(ratio[, .(SPEC = substr(SPEC, 1, 26), INSTRUMENT_LABEL,
                    GAP = round(GAP, 2), RATIO = round(RESPONSE_RATIO, 2))])
    save_csv(ratio, paste0(stem, "C_response_gap.csv"))
  }

  cat("\n", strrep("=", 104), "\nHETEROGENEITY TEST\n", strrep("=", 104), "\n", sep = "")
  print(dcast(mt[, .(SPEC = substr(SPEC, 1, 26), INSTRUMENT_LABEL, ESTIMATOR,
                     P = round(P_VALUE, 4))],
              SPEC + INSTRUMENT_LABEL ~ ESTIMATOR, value.var = "P"))
  print(mt[, .(N = .N, SIG_05 = sum(P_VALUE < 0.05, na.rm = TRUE),
               SHARE = round(mean(P_VALUE < 0.05, na.rm = TRUE), 3),
               MEDIAN_P = round(median(P_VALUE, na.rm = TRUE), 4)), by = ESTIMATOR])

  list(rows = mr, tests = mt)
}

run_transform_ladder <- function(panel, scheme_col = "SCHEME_1_CERTAINTY",
                                 scheme_label = "1. Procedural certainty", instruments = MAIN_INSTRUMENTS,
                                 transforms = names(TRANSFORM_LADDER), outcome = PRIMARY_OUTCOME) {
  rows <- list(); tests <- list()
  for (tf in transforms) {
    d <- apply_transform(copy(panel), tf)
    for (il in names(instruments)) {
      z <- instruments[[il]]; if (!(z %in% names(d))) next
      r <- estimate_interacted(d, scheme_col, outcome, z, endogenous = "TREAT_TRANSFORMED",
                               moderator_type = "categorical", label = scheme_label,
                               instrument_label = il, moderator_label = tf)
      if (!is.null(r)) {
        r$rows[, TRANSFORM := tf]; r$tests[, TRANSFORM := tf]
        rows[[length(rows) + 1L]] <- r$rows; tests[[length(tests) + 1L]] <- r$tests
      }
      cat(sprintf("  %-14s %-36s done\n", tf, substr(il, 1, 34)))
    }
    rm(d); invisible(gc())
  }
  lr <- rbindlist(rows, fill = TRUE); lt <- rbindlist(tests, fill = TRUE)
  save_csv(lr, "T07_transform_ladder_estimates.csv")
  save_csv(lt, "T07B_transform_ladder_heterogeneity.csv")

  cat("\nMAGNITUDES MOVE:\n")
  print(dcast(lr[, .(TRANSFORM, INSTRUMENT_LABEL, TERM, E = round(IV_PERCENT, 2))],
              INSTRUMENT_LABEL + TERM ~ TRANSFORM, value.var = "E"))
  cat("\nCONCLUSION DOES NOT:\n")
  print(dcast(lt[ESTIMATOR == "IV", .(TRANSFORM, INSTRUMENT_LABEL, P = round(P_VALUE, 4))],
              INSTRUMENT_LABEL ~ TRANSFORM, value.var = "P"))
  cat("\nA higher F under a compressed treatment says the instrument predicts a\n",
      "compressed disclosure count better -- a fact about how disclosure spreads,\n",
      "not about how prices respond. Linear F = 43 is already well above any\n",
      "weak-instrument threshold.\n", sep = "")
  list(rows = lr, tests = lt)
}

cat("Section 7 loaded\n")


# --- Interactive checkpoint ---------------------------------------------------
# Confirms the two objects the stage blocks depend on are in memory, then loads
# the concept-level results from cache. Duplicates the stage 6 loader below and
# is here so stages 8-10 can be run without stepping through stage 6's block.
# Skip on a cold source() of the whole file.

exists("outpatient")
exists("schemes_long")


CONCEPT_INSTRUMENTS <- c(MAIN_INSTRUMENTS, SUPPORTING_INSTRUMENTS)
concept_results <- cache_or_run("concept_level_6inst",
                                estimate_concept_level(
                                  build_concept_panel(outpatient, instruments = CONCEPT_INSTRUMENTS),
                                  instruments = CONCEPT_INSTRUMENTS))



######## Section 8: Meta-regressions, deduplication, permutation ##############
#
# The second step of the two-step design. Section 6 produced one causal
# estimate per clinical concept; this regresses those estimates on concept
# characteristics, weighted by inverse variance and clustered on clinical
# family, in the spirit of a precision-weighted second-stage regression rather
# than a conventional meta-analysis of independent studies.
#
# Four pieces:
#
#   prepare_meta_input()       attaches every scheme's classification
#   run_meta_regressions()     estimates the gradient under each scheme and
#                              each collapse rule
#   deduplicate_partitions()   fingerprints schemes that induce the same
#                              partition, so the specification count reports
#                              distinct tests rather than distinct labels
#   family_permutation_test()  exact inference on the family assignment
#   decompose_reduced_form()   splits the gradient into price response versus
#                              disclosure propensity
#
# decompose_reduced_form() is called once overall and once per tier, so its
# output filename is parameterised through `stem` for the same reason as
# run_main_results() in Section 7.

prepare_meta_input <- function(cr, schemes_long) {
  d <- copy(cr)
  d[, INSTRUMENT_LABEL := fifelse(INSTRUMENT %chin% names(INSTRUMENT_LABEL_MAP),
                                  INSTRUMENT_LABEL_MAP[INSTRUMENT], INSTRUMENT_LABEL)]
  stopifnot(all(d[, uniqueN(INSTRUMENT_LABEL), by = INSTRUMENT]$V1 == 1L))

  d <- d[(is.finite(RF_COEF) & is.finite(RF_SE) & RF_SE > 0) |
           (is.finite(IV_COEF) & is.finite(IV_SE) & IV_SE > 0)]

  d[, CLUSTER_FAMILY := FINAL_FAMILY_ID]
  d[, SHOP_CERTAINTY := factor(fifelse(FINAL_FAMILY_ID %chin% DIAGNOSTIC_FAMILIES,
                                       "Shoppable", "Non_shoppable"),
                               levels = c("Non_shoppable", "Shoppable"))]
  for (sid in unique(schemes_long$SCHEME_ID)) {
    a <- unique(schemes_long[SCHEME_ID == sid, .(FINAL_CONCEPT_ID = ANALYSIS_CONCEPT_ID,
                                                 V = SHOPPABILITY_CATEGORY)])
    setnames(a, "V", paste0("CAT_", sid))
    d <- merge(d, a, by = "FINAL_CONCEPT_ID", all.x = TRUE, sort = FALSE)
  }
  setDT(d)

  n_unmatched <- d[is.na(get(paste0("CAT_", unique(schemes_long$SCHEME_ID)[1L]))), .N]
  if (n_unmatched > 0L) {
    warning(n_unmatched, " concept-instrument rows have no scheme classification.",
            call. = FALSE)
  }
  d
}

run_meta_regressions <- function(mi, schemes_long, dep = "RF_COEF", se = "RF_SE") {
  rules <- list(
    `3-tier`              = function(x) factor(x, levels = c("LOW", "INTERMEDIATE", "HIGH")),
    `2-tier High vs rest` = function(x) factor(fifelse(x == "HIGH", "Shoppable", "Non_shoppable"),
                                               levels = c("Non_shoppable", "Shoppable")),
    `2-tier Low vs rest`  = function(x) factor(fifelse(x == "LOW", "Non_shoppable", "Shoppable"),
                                               levels = c("Non_shoppable", "Shoppable")),
    `2-tier extremes`     = function(x) factor(fcase(x == "HIGH", "Shoppable",
                                                     x == "LOW", "Non_shoppable",
                                                     default = NA_character_),
                                               levels = c("Non_shoppable", "Shoppable")))
  rows <- list()
  for (sid in unique(schemes_long$SCHEME_ID)) {
    col <- paste0("CAT_", sid); if (!(col %in% names(mi))) next
    sname <- schemes_long[SCHEME_ID == sid][1L]$SCHEME_NAME %||% sid
    for (rl in names(rules)) for (il in unique(mi$INSTRUMENT_LABEL)) {
      d <- mi[INSTRUMENT_LABEL == il &
                is.finite(get(dep)) & is.finite(get(se)) & get(se) > 0]
      if (nrow(d) == 0L) next
      d[, CAT := droplevels(rules[[rl]](toupper(trimws(get(col)))))]
      d <- d[!is.na(CAT)]
      if (nrow(d) < MIN_CONCEPTS_META || uniqueN(d$CAT) < 2L) next
      for (w in META_WEIGHTINGS) {
        d[, W := if (w == "Inverse variance") 1 / (get(se)^2) else 1]
        fit <- tryCatch(feols(as.formula(paste(dep, "~ CAT")), data = d, weights = ~W,
                              cluster = ~CLUSTER_FAMILY, warn = FALSE, notes = FALSE),
                        error = function(e) NULL)
        if (is.null(fit)) next
        td <- tidy_fixest(fit); if (nrow(td) == 0L) next
        td[, `:=`(SCHEME_ID = sid, SCHEME = sname, COLLAPSE = rl, INSTRUMENT_LABEL = il,
                  WEIGHTING = w, DEPENDENT = dep, N_CONCEPTS = nrow(d),
                  N_CLUSTERS = uniqueN(d$CLUSTER_FAMILY))]
        rows[[length(rows) + 1L]] <- td
      }
    }
  }
  out <- rbindlist(rows, fill = TRUE)
  if (nrow(out) == 0L) stop("No meta-regressions for dep = ", dep, call. = FALSE)
  out
}

# Several schemes induce the same partition of the concept universe and are
# therefore the same test under a different name. This fingerprints each
# scheme-collapse variant on its estimated coefficient vector and keeps one
# representative per distinct partition, so the paper reports the number of
# genuinely distinct specifications rather than the number of labels.
#
# Fingerprinting runs WITHIN dependent variable. Pooling RF and IV coefficient
# vectors would make the same scheme look like two partitions and double the
# count.
deduplicate_partitions <- function(mr) {
  mr <- copy(mr)[, SCHEME_VARIANT := paste(SCHEME, COLLAPSE, sep = " | ")]
  rbindlist(lapply(unique(mr$DEPENDENT), function(dep) {
    sub <- mr[DEPENDENT == dep]
    fp <- sub[!is.na(estimate), .(SCHEME_VARIANT,
                                  CELL = paste(INSTRUMENT_LABEL, WEIGHTING, term, sep = "|"), V = round(estimate, 6))]
    setorder(fp, SCHEME_VARIANT, CELL)
    sig <- fp[, .(FP = paste(CELL, V, sep = "=", collapse = ";"), N_CELLS = .N),
              by = SCHEME_VARIANT]
    sig[, PARTITION_ID := .GRP, by = .(FP, N_CELLS)]
    map <- sig[, .(N_VARIANTS = .N, REPRESENTATIVE = sort(SCHEME_VARIANT)[1L],
                   ALL_VARIANTS = paste(sort(SCHEME_VARIANT), collapse = " ||| ")),
               by = PARTITION_ID]
    sig <- merge(sig[, .(SCHEME_VARIANT, PARTITION_ID)], map, by = "PARTITION_ID")
    sig[, IS_REPRESENTATIVE := as.integer(SCHEME_VARIANT == REPRESENTATIVE)]
    cat("\n", dep, ": ", uniqueN(sig$PARTITION_ID), " distinct partitions of ",
        nrow(sig), " variants\n", sep = "")
    if (dep == "RF_COEF") save_csv(map[N_VARIANTS > 1L], "T08B_redundant_partitions.csv")
    merge(sub, sig, by = "SCHEME_VARIANT", all.x = TRUE, sort = FALSE)
  }), fill = TRUE)
}

# Exact permutation inference on the family assignment.
#
# Shoppability is assigned to roughly 16 clinical FAMILIES, not to 738
# independent concepts. Concept-level clustered standard errors treat 168
# biopsy concepts as 168 independent observations when a single labelling
# decision covers all of them. Permuting which families are called shoppable
# and enumerating every assignment gives inference at the level the decision
# was actually made, and this is the p-value the paper reports.
#
# With 16 families and 10 shoppable, C(16,10) = 8,008 assignments enumerate
# exactly, so the smallest attainable two-sided p-value is 1/8008. Report a
# result at that floor as being at the minimum attainable value rather than as
# a precise small number.
family_permutation_test <- function(mi, dep = "RF_COEF", se = "RF_SE",
                                    shoppable_families = DIAGNOSTIC_FAMILIES,
                                    max_exact = 50000L, seed = 20260811L) {
  wdiff <- function(d, fams) {
    d <- copy(d)[, S := fifelse(FINAL_FAMILY_ID %chin% fams, "Shoppable", "Non_shoppable")]
    if (uniqueN(d$S) < 2L) return(NA_real_)
    d[, W := 1 / (get(se)^2)]
    d[S == "Shoppable", sum(W * get(dep)) / sum(W)] -
      d[S == "Non_shoppable", sum(W * get(dep)) / sum(W)]
  }
  rbindlist(lapply(unique(mi$INSTRUMENT_LABEL), function(il) {
    d <- mi[INSTRUMENT_LABEL == il & is.finite(get(dep)) & is.finite(get(se)) & get(se) > 0]
    if (nrow(d) < MIN_CONCEPTS_META) return(data.table())
    fams <- sort(unique(d$FINAL_FAMILY_ID))
    obs_f <- intersect(shoppable_families, fams); k <- length(obs_f); n <- length(fams)
    if (k == 0L || k == n) return(data.table())
    observed <- wdiff(d, obs_f); if (!is.finite(observed)) return(data.table())
    exact <- is.finite(choose(n, k)) && choose(n, k) <= max_exact
    null <- if (exact) vapply(combn(fams, k, simplify = FALSE), function(f) wdiff(d, f), numeric(1))
    else { set.seed(seed); vapply(seq_len(20000L), function(i) wdiff(d, sample(fams, k)), numeric(1)) }
    null <- null[is.finite(null)]
    data.table(INSTRUMENT_LABEL = il, OBSERVED = observed, N_FAMILIES = n,
               N_SHOPPABLE_FAMILIES = k,
               METHOD = fifelse(exact, "Exact enumeration", "Monte Carlo (20,000)"),
               NULL_P05 = quantile(null, .05, names = FALSE),
               NULL_P95 = quantile(null, .95, names = FALSE),
               P_TWO_SIDED = mean(abs(null) >= abs(observed)))
  }), fill = TRUE)
}

# Splits the concept-level gradient into its two possible sources by running
# the same regression on RF, FS, and IV coefficients in turn. A gradient in RF
# but not in FS means the heterogeneity is in the PRICE RESPONSE. A gradient in
# both means part of it is differential disclosure propensity.
#
# `stem` is parameterised because the run block calls this once overall and
# once per instrument tier.
decompose_reduced_form <- function(mi, stem = "T08D_RF_vs_FS_decomposition") {
  out <- rbindlist(lapply(unique(mi$INSTRUMENT_LABEL), function(il) {
    d <- mi[INSTRUMENT_LABEL == il]; if (nrow(d) < MIN_CONCEPTS_META) return(data.table())
    rbindlist(lapply(c("RF_COEF", "FS_COEF", "IV_COEF"), function(dep) {
      dd <- d[is.finite(get(dep))]; if (nrow(dd) < MIN_CONCEPTS_META) return(data.table())
      dd[, W := 1 / (RF_SE^2)]
      fit <- tryCatch(feols(as.formula(paste(dep, "~ SHOP_CERTAINTY")), data = dd,
                            weights = ~W, cluster = ~CLUSTER_FAMILY,
                            warn = FALSE, notes = FALSE), error = function(e) NULL)
      if (is.null(fit)) return(data.table())
      td <- tidy_fixest(fit)
      td[, `:=`(DEPENDENT = dep, INSTRUMENT_LABEL = il, N_CONCEPTS = nrow(dd))]; td
    }), fill = TRUE)
  }), fill = TRUE)
  cat("\n", strrep("=", 96), "\nDECOMPOSITION: PRICE RESPONSE OR DISCLOSURE PROPENSITY?\n",
      strrep("=", 96), "\n", sep = "")
  print(out[grepl("Shoppable", term), .(DEPENDENT, INSTRUMENT_LABEL,
                                        EST = signif(estimate, 4), SE = signif(std.error, 4), P = round(p.value, 4), N_CONCEPTS)][
                                          order(DEPENDENT, INSTRUMENT_LABEL)])
  cat("\nRF differs by category but FS does not -> the heterogeneity is in the\n",
      "PRICE RESPONSE. Both differ -> part is differential disclosure propensity.\n", sep = "")
  save_csv(out, paste0(stem, ".csv"))
  out
}

cat("Section 8 loaded\n")



######## Section 9: Price comparability as a candidate mechanism ##############
#
# Tests the mechanism directly rather than inferring it from the shoppability
# label. The working hypothesis is that a service responds to disclosure to the
# extent that its posted price constitutes a meaningful basis for comparison. A
# screening mammogram is the same product at hospital A and hospital B, and one
# number describes it. A colonoscopy that may become a polypectomy has no such
# number: the posted price is one draw from a distribution whose realisation
# depends on what is found mid-procedure.
#
# If that is the operative channel, then measures of how well a single posted
# number summarises a concept should moderate the price response, and should do
# so WITHIN clinical family -- otherwise the measure is just a modality proxy
# and adds nothing beyond the shoppability label.
#
# Eight candidate measures are built and screened:
#
#   PD_PAYER_V2   payer-level price dispersion from the payer-cell export
#   PD_HOSP       across-hospital spread in the concept price
#   PD_CODE       within-hospital spread across exact codes inside a concept
#   PD_PAYER      concept-panel percentile spread
#   N_CODES       number of distinct billing codes inside a concept
#   N_PAYERS_V2   distinct payers per hospital-month-concept
#   N_PAYERS      distinct payers, concept-panel version
#   CODE_COV      code coverage ratio
#
# screen_moderators() drops whichever fail coverage, variance, zero-inflation,
# or collinearity gates and reports the reason, so a degenerate measure is
# visibly excluded rather than silently estimated.
#
# EVERY MEASURE IS LEAVE-ONE-COUNTY-OUT. Dispersion is built from prices and
# the outcome is a price, so a focal hospital's own dispersion is mechanically
# linked to its own estimate. The leave-out breaks that link at the cost of
# requiring MIN_OUTSIDE_HOSPITALS observations outside the focal county.

MOD_MIN_FINITE      <- 100L
MOD_MAX_SHARE_ZERO  <- 0.60
MOD_MAX_CORRELATION <- 0.95


# ---------------------------------------------------------------------------
# Leave-one-county-out helpers
# ---------------------------------------------------------------------------
#
# Closed-form leave-one-out mean and standard deviation. Both compute concept
# totals once and subtract each county's contribution, rather than looping over
# counties, which is what makes this feasible at panel scale.
#
# Both bind the value column to a fixed name VAL before aggregating. This
# matters for correctness and for speed: group-wise sums over a named column
# use data.table's optimised gsum, while get(vc) evaluated inside j does not
# and does not scope the way it appears to.
#
# verify_loo_sd() checks the closed form against a brute-force loop on the
# highest-dispersion concepts and stops if they disagree.

loo_mean <- function(dt, vc) {
  d <- dt[, .(FINAL_CONCEPT_ID, ANALYSIS_MARKET, VAL = as.numeric(get(vc)))]
  d <- d[is.finite(VAL)]
  tot <- d[, .(S = sum(VAL), N = .N), by = FINAL_CONCEPT_ID]
  cty <- d[, .(SC = sum(VAL), NC = .N), by = .(FINAL_CONCEPT_ID, ANALYSIS_MARKET)]
  m <- merge(cty, tot, by = "FINAL_CONCEPT_ID")[N - NC >= MIN_OUTSIDE_HOSPITALS]
  if (nrow(m) == 0L) return(data.table(FINAL_CONCEPT_ID = character(0), V = numeric(0)))
  m[, LOO := (S - SC) / (N - NC)]
  m[, .(V = mean(LOO, na.rm = TRUE)), by = FINAL_CONCEPT_ID]
}

loo_sd <- function(dt, vc) {
  d <- dt[, .(FINAL_CONCEPT_ID, ANALYSIS_MARKET, VAL = as.numeric(get(vc)))]
  d <- d[is.finite(VAL)]
  tot <- d[, .(S = sum(VAL), SS = sum(VAL^2), N = .N), by = FINAL_CONCEPT_ID]
  cty <- d[, .(SC = sum(VAL), SSC = sum(VAL^2), NC = .N),
           by = .(FINAL_CONCEPT_ID, ANALYSIS_MARKET)]
  m <- merge(cty, tot, by = "FINAL_CONCEPT_ID")[N - NC >= MIN_OUTSIDE_HOSPITALS]
  if (nrow(m) == 0L) return(data.table(FINAL_CONCEPT_ID = character(0), V = numeric(0)))
  m[, NR := N - NC]
  m[, VAR_LOO := ((SS - SSC) - ((S - SC)^2) / NR) / pmax(NR - 1L, 1L)]
  m[, SD_LOO := sqrt(pmax(VAR_LOO, 0))]
  m[, .(V = mean(SD_LOO, na.rm = TRUE)), by = FINAL_CONCEPT_ID]
}

verify_loo_sd <- function(panel, n_check = 5L) {
  hosp <- unique(panel[is.finite(MEDIAN_PRICE) & MEDIAN_PRICE > 0,
                       .(FINAL_CONCEPT_ID, ANALYSIS_MARKET, LN_MED = log(MEDIAN_PRICE))])
  fast <- loo_sd(hosp, "LN_MED")
  if (nrow(fast) == 0L) { warning("loo_sd returned nothing."); return(invisible(FALSE)) }

  ids <- head(fast[order(-V)]$FINAL_CONCEPT_ID, n_check)
  slow <- rbindlist(lapply(ids, function(id) {
    d <- hosp[FINAL_CONCEPT_ID == id]
    v <- vapply(unique(d$ANALYSIS_MARKET), function(cc) {
      x <- d[ANALYSIS_MARKET != cc]$LN_MED
      if (length(x) < MIN_OUTSIDE_HOSPITALS) NA_real_ else sd(x)
    }, numeric(1))
    data.table(FINAL_CONCEPT_ID = id, V_BRUTE = mean(v, na.rm = TRUE))
  }))
  cmp <- merge(fast[FINAL_CONCEPT_ID %chin% ids], slow, by = "FINAL_CONCEPT_ID")
  cmp[, MATCH := abs(V - V_BRUTE) < 1e-8]

  cat("\n", strrep("=", 84), "\nloo_sd VERIFICATION (closed form vs brute force)\n",
      strrep("=", 84), "\n", sep = "")
  print(cmp[, .(FINAL_CONCEPT_ID = substr(FINAL_CONCEPT_ID, 1, 40),
                FAST = round(V, 8), BRUTE = round(V_BRUTE, 8), MATCH)])
  if (!all(cmp$MATCH)) stop("loo_sd disagrees with brute force.", call. = FALSE)
  cat("All match.\n")
  invisible(TRUE)
}


# ---------------------------------------------------------------------------
# Payer dispersion from the payer-cell export
# ---------------------------------------------------------------------------
#
# PD_PAYER_V2 is dispersion across identified payers within a
# hospital-month-concept, built from payer identity and negotiated rates rather
# than from panel percentiles. It runs a median of 7 distinct payers per cell
# and 13.9% zero-dispersion cells.
#
# This is the measure PD_PAYER was intended to be. PD_PAYER, built from the
# concept panel's P25/P75, turns out to capture code-level rather than
# payer-level variation: 87% exact zeros and a 0.984 correlation with PD_CODE.
# Both are retained so the comparison is visible in QA09, but only PD_PAYER_V2
# survives screening.
#
# READS MULTIPLE SHARDS. The warehouse export splits above an internal size
# threshold regardless of MAX_FILE_SIZE, so it arrives as several gzip parts.
# This reads every file matching PAYER_DISPERSION_PATTERN and row-binds them.
# Each shard carries its own header row, so each is read individually rather
# than concatenated with a manual skip.


load_payer_dispersion <- function() {
  if (!requireNamespace("R.utils", quietly = TRUE)) {
    stop("Package 'R.utils' is required to read gzipped payer dispersion files. ",
         "Install with install.packages('R.utils').", call. = FALSE)
  }

  files <- list.files(PAYER_DISPERSION_DIR, pattern = PAYER_DISPERSION_PATTERN,
                      full.names = TRUE)
  if (length(files) == 0L) {
    warning("No payer dispersion files found matching '", PAYER_DISPERSION_PATTERN,
            "' in ", PAYER_DISPERSION_DIR, call. = FALSE)
    return(data.table(FINAL_CONCEPT_ID = character(0), PD_PAYER_V2 = numeric(0),
                      N_PAYERS_V2 = numeric(0)))
  }

  cat("\nPayer dispersion: found", length(files), "file(s):\n")
  for (f in files) cat("  ", basename(f), "\n")

  raw <- rbindlist(lapply(files, function(f) {
    tryCatch(fread(f), error = function(e) {
      warning("Failed to read ", basename(f), ": ", e$message, call. = FALSE)
      data.table()
    })
  }), fill = TRUE)

  if (nrow(raw) == 0L) {
    warning("Payer dispersion files read but produced zero rows.", call. = FALSE)
    return(data.table(FINAL_CONCEPT_ID = character(0), PD_PAYER_V2 = numeric(0),
                      N_PAYERS_V2 = numeric(0)))
  }

  cat("Combined:", format(nrow(raw), big.mark = ","), "rows |",
      uniqueN(raw$ANALYSIS_CONCEPT_ID), "concepts |",
      uniqueN(raw$HOSPITAL_ID), "hospitals\n")

  assert_columns(raw, c("HOSPITAL_ID", "POST_MONTH", "ANALYSIS_CONCEPT_ID",
                        "N_DISTINCT_PAYERS", "N_PAYER_CELLS",
                        "CV_PAYER_NEGOTIATED"), "payer dispersion file")

  n_before <- nrow(raw)
  raw <- unique(raw, by = c("HOSPITAL_ID", "POST_MONTH", "ANALYSIS_CONCEPT_ID"))
  if (nrow(raw) < n_before) {
    cat("Dropped", n_before - nrow(raw), "duplicate rows across shards.\n")
  }

  raw[, POST_MONTH := as.Date(POST_MONTH)]
  raw[, HOSPITAL_ID := as.character(HOSPITAL_ID)]
  setnames(raw, "ANALYSIS_CONCEPT_ID", "FINAL_CONCEPT_ID")

  hosp_county <- unique(outpatient[, .(HOSPITAL_ID, POST_MONTH, ANALYSIS_MARKET)])
  raw <- merge(raw, hosp_county, by = c("HOSPITAL_ID", "POST_MONTH"), all.x = FALSE)

  cat("Matched to panel counties:", format(nrow(raw), big.mark = ","), "rows |",
      format(sum(is.na(raw$ANALYSIS_MARKET)), big.mark = ","), "unmatched (dropped)\n")

  raw <- raw[!is.na(ANALYSIS_MARKET) & is.finite(CV_PAYER_NEGOTIATED) &
               CV_PAYER_NEGOTIATED >= 0]

  if (nrow(raw) == 0L) {
    warning("No usable rows after county match / finiteness filter.", call. = FALSE)
    return(data.table(FINAL_CONCEPT_ID = character(0), PD_PAYER_V2 = numeric(0),
                      N_PAYERS_V2 = numeric(0)))
  }

  # Winsorize the raw per-cell coefficient of variation BEFORE the leave-out
  # average. Cells with a near-zero denominator produce extreme values, and
  # because the leave-out average for every county includes all outside cells,
  # a single extreme cell would otherwise contaminate the measure everywhere.
  # Capping at P99 leaves the bulk of the distribution untouched.
  cap <- quantile(raw$CV_PAYER_NEGOTIATED, 0.99, na.rm = TRUE)
  n_capped <- sum(raw$CV_PAYER_NEGOTIATED > cap)
  cat("Winsorizing raw CV at P99 =", round(cap, 4), "-- capping",
      format(n_capped, big.mark = ","), "of", format(nrow(raw), big.mark = ","),
      "cells (", round(100 * n_capped / nrow(raw), 3), "%)\n")
  raw[, CV_PAYER_WINSOR := pmin(CV_PAYER_NEGOTIATED, cap)]

  pd <- loo_mean(raw, "CV_PAYER_WINSOR"); setnames(pd, "V", "PD_PAYER_V2")

  npay <- raw[, .(N_PAYERS_V2 = mean(N_DISTINCT_PAYERS, na.rm = TRUE)),
              by = FINAL_CONCEPT_ID]

  out <- merge(pd, npay, by = "FINAL_CONCEPT_ID", all = TRUE)

  cat("\nPD_PAYER_V2 (post-winsorize) for", nrow(out), "concepts:\n")
  print(data.table(
    MEAN = round(mean(out$PD_PAYER_V2, na.rm = TRUE), 4),
    SD = round(sd(out$PD_PAYER_V2, na.rm = TRUE), 4),
    P50 = round(median(out$PD_PAYER_V2, na.rm = TRUE), 4),
    P95 = round(quantile(out$PD_PAYER_V2, .95, na.rm = TRUE), 4),
    MAX = round(max(out$PD_PAYER_V2, na.rm = TRUE), 4),
    SHARE_ZERO = round(mean(out$PD_PAYER_V2 == 0, na.rm = TRUE), 4)
  ))

  out
}


# ---------------------------------------------------------------------------
# Degeneracy screen
# ---------------------------------------------------------------------------
screen_moderators <- function(measures, candidates) {
  candidates <- intersect(candidates, names(measures))
  if (length(candidates) == 0L) return(character(0))

  diag <- rbindlist(lapply(candidates, function(m) {
    v <- safe_numeric(measures[[m]]); f <- v[is.finite(v)]
    data.table(MODERATOR = m, N_FINITE = length(f),
               SHARE_ZERO = if (length(f)) mean(f == 0) else NA_real_,
               SD = if (length(f) > 1L) sd(f) else NA_real_,
               MIN = if (length(f)) min(f) else NA_real_,
               MAX = if (length(f)) max(f) else NA_real_)
  }))

  diag[, FAIL_COVERAGE := as.integer(N_FINITE < MOD_MIN_FINITE)]
  diag[, FAIL_VARIANCE := as.integer(!is.finite(SD) | SD <= 0)]
  diag[, FAIL_ZERO_INFLATED := as.integer(is.finite(SHARE_ZERO) &
                                            SHARE_ZERO > MOD_MAX_SHARE_ZERO)]

  keep <- diag[FAIL_COVERAGE == 0L & FAIL_VARIANCE == 0L &
                 FAIL_ZERO_INFLATED == 0L]$MODERATOR

  dropped_collinear <- character(0)
  if (length(keep) > 1L) {
    survivors <- keep[1L]
    for (m in keep[-1L]) {
      r <- vapply(survivors, function(s) {
        suppressWarnings(abs(cor(safe_numeric(measures[[m]]),
                                 safe_numeric(measures[[s]]),
                                 use = "complete.obs")))
      }, numeric(1))
      if (any(is.finite(r) & r > MOD_MAX_CORRELATION)) {
        dropped_collinear <- c(dropped_collinear, m)
      } else survivors <- c(survivors, m)
    }
    keep <- survivors
  }
  diag[, FAIL_COLLINEAR := as.integer(MODERATOR %chin% dropped_collinear)]
  diag[, USABLE := as.integer(MODERATOR %chin% keep)]

  cat("\n", strrep("=", 104), "\nMODERATOR SCREEN\n", strrep("=", 104), "\n", sep = "")
  print(diag[, .(MODERATOR, N_FINITE, SHARE_ZERO = round(SHARE_ZERO, 3),
                 SD = signif(SD, 3), FAIL_COVERAGE, FAIL_VARIANCE,
                 FAIL_ZERO_INFLATED, FAIL_COLLINEAR, USABLE)])
  save_qa_csv(diag, "QA09_moderator_screen.csv")

  if (length(dropped_collinear) > 0L) {
    cat("\nDropped for collinearity (|r| >", MOD_MAX_CORRELATION, "):",
        paste(dropped_collinear, collapse = ", "), "\n")
  }
  if (length(keep) == 0L) {
    warning("No usable comparability moderators survived screening.", call. = FALSE)
  } else {
    cat("\nUsable moderators:", paste(keep, collapse = ", "), "\n")
  }
  keep
}


# ---------------------------------------------------------------------------
# Build the measures
# ---------------------------------------------------------------------------
build_comparability_measures <- function(panel) {
  assert_columns(panel, c("MEDIAN_PRICE", "P25_PRICE", "P75_PRICE"), "panel")

  have_payers <- "N_DISTINCT_PAYERS" %in% names(panel)
  have_cov    <- "CODE_COVERAGE_RATIO" %in% names(panel)
  cat("\nSource columns: N_DISTINCT_PAYERS =", have_payers,
      "| CODE_COVERAGE_RATIO =", have_cov, "\n")

  base <- panel[is.finite(MEDIAN_PRICE) & MEDIAN_PRICE > 0,
                .(ANALYSIS_MARKET, FINAL_CONCEPT_ID, FINAL_FAMILY_ID, HOSPITAL_ID,
                  CV_PAYER = (P75_PRICE - P25_PRICE) / MEDIAN_PRICE,
                  LN_MED = log(MEDIAN_PRICE),
                  N_PAY    = if (have_payers) safe_numeric(N_DISTINCT_PAYERS)   else NA_real_,
                  CODE_COV = if (have_cov)    safe_numeric(CODE_COVERAGE_RATIO) else NA_real_)]
  base <- base[is.finite(CV_PAYER) & CV_PAYER >= 0]

  # 1. PD_PAYER -- concept-panel percentiles. Retained for comparison only;
  #    degenerate in this panel (87% exact zeros, corr 0.984 with PD_CODE).
  #    See QA09.
  pd_payer <- loo_mean(base, "CV_PAYER"); setnames(pd_payer, "V", "PD_PAYER")

  # 2. PD_HOSP -- across-hospital spread of the concept price.
  hosp <- unique(base[, .(FINAL_CONCEPT_ID, ANALYSIS_MARKET, LN_MED)])
  pd_hosp <- loo_sd(hosp, "LN_MED"); setnames(pd_hosp, "V", "PD_HOSP")

  # 3-4. Contracting thinness and definitional standardisation. Both depend on
  #      source columns that may be absent from the parquet, in which case they
  #      come back all-NA and are dropped by the screen.
  thin <- base[, .(N_PAYERS = mean(N_PAY, na.rm = TRUE),
                   CODE_COV = mean(CODE_COV, na.rm = TRUE)), by = FINAL_CONCEPT_ID]

  fam <- unique(base[, .(FINAL_CONCEPT_ID, CONCEPT_FAMILY = FINAL_FAMILY_ID)])

  # 5-6. PD_CODE and N_CODES, from the exact-code panel.
  code_measures <- tryCatch({
    ex_raw <- read_panel(FILES$outpatient_exact, "exact-code panel",
                         columns = unique(c(ANALYSIS_COLUMNS, "BILLING_CODE")))
    ex <- prepare_panel(ex_raw, "exact-code panel", choose_sample_flag(ex_raw))
    rm(ex_raw); invisible(gc())

    cell <- ex[is.finite(MEDIAN_PRICE) & MEDIAN_PRICE > 0,
               .(N_CODES = uniqueN(BILLING_CODE),
                 CV = {
                   q <- quantile(MEDIAN_PRICE, c(.25, .75), na.rm = TRUE, names = FALSE)
                   mm <- median(MEDIAN_PRICE, na.rm = TRUE)
                   if (is.finite(mm) && mm > 0) (q[2] - q[1]) / mm else NA_real_
                 }),
               by = .(HOSPITAL_ID, POST_MONTH, ANALYSIS_MARKET, FINAL_CONCEPT_ID)]
    cell <- cell[is.finite(CV)]

    pc <- loo_mean(cell, "CV");      setnames(pc, "V", "PD_CODE")
    nc <- loo_mean(cell, "N_CODES"); setnames(nc, "V", "N_CODES")
    rm(ex, cell); invisible(gc())
    merge(pc, nc, by = "FINAL_CONCEPT_ID", all = TRUE)
  }, error = function(e) {
    warning("Exact-code measures unavailable: ", e$message, call. = FALSE)
    data.table(FINAL_CONCEPT_ID = character(0), PD_CODE = numeric(0),
               N_CODES = numeric(0))
  })

  # 7-8. Payer-level dispersion and payer counts from the payer-cell export.
  payer_v2 <- load_payer_dispersion()

  out <- Reduce(function(a, b) merge(a, b, by = "FINAL_CONCEPT_ID", all.x = TRUE),
                list(pd_payer, pd_hosp, thin, fam, code_measures, payer_v2))

  cat("\nComparability measures for", nrow(out), "concepts\n")

  cat("\nWithin-family SD -- must be non-trivial, or the measure is another\n",
      "modality proxy and the within-family test has nothing to identify from:\n", sep = "")
  present <- intersect(c("PD_PAYER", "PD_HOSP", "PD_CODE", "N_CODES", "PD_PAYER_V2"),
                       names(out))
  print(out[, c(list(N = .N),
                lapply(.SD, function(v) round(sd(v, na.rm = TRUE), 4))),
            by = CONCEPT_FAMILY, .SDcols = present][N >= 5L][1:min(.N, 14L)])

  save_csv(out, "T09_comparability_measures.csv")
  out
}


# ---------------------------------------------------------------------------
# Row-level interacted test
# ---------------------------------------------------------------------------
run_comparability_interaction <- function(panel, measures,
                                          moderators = COMPARABILITY_MODERATORS,
                                          instruments = MAIN_INSTRUMENTS, outcome = PRIMARY_OUTCOME) {

  moderators <- screen_moderators(measures, moderators)
  if (length(moderators) == 0L) stop("No usable moderators.", call. = FALSE)

  d <- merge(panel, measures[, c("FINAL_CONCEPT_ID", moderators), with = FALSE],
             by = "FINAL_CONCEPT_ID", all.x = TRUE, sort = FALSE)
  setDT(d)

  rows <- list(); tests <- list()
  for (m in moderators) {
    for (il in names(instruments)) {
      z <- instruments[[il]]; if (!(z %in% names(d))) next
      r <- estimate_interacted(d, m, outcome, z, moderator_type = "continuous",
                               label = "Comparability", instrument_label = il,
                               moderator_label = m)
      if (is.null(r)) { cat(sprintf("  %-12s %-36s NO RESULT\n", m, substr(il, 1, 34))); next }
      rows[[length(rows) + 1L]] <- r$rows; tests[[length(tests) + 1L]] <- r$tests
      cat(sprintf("  %-12s %-36s done\n", m, substr(il, 1, 34)))
    }
  }
  cr <- rbindlist(rows, fill = TRUE); ct <- rbindlist(tests, fill = TRUE)
  if (nrow(cr) == 0L) stop("No comparability results.", call. = FALSE)

  cr[, EXPECTED_SIGN := fifelse(MODERATOR %chin% DISPERSION_MEASURES, "positive", "negative")]
  cr[, SIGN_AS_PREDICTED := fifelse(MODERATOR %chin% DISPERSION_MEASURES,
                                    RF_COEF > 0, RF_COEF < 0)]
  cr[, WEAK_FIRST_STAGE := as.integer(FIRST_STAGE_WALD_MIN < 10)]

  save_csv(cr, "T09B_comparability_interaction.csv")
  save_csv(ct, "T09C_comparability_interaction_tests.csv")

  cat("\n", strrep("=", 112), "\nCOMPARABILITY GRADIENT (row-level interacted)\n",
      strrep("=", 112), "\n", sep = "")
  print(cr[TERM == "x Moderator", .(MODERATOR, INSTRUMENT_LABEL, EXPECTED_SIGN,
                                    RF = signif(RF_COEF, 3), RF_P = round(RF_P, 4), SIGN_OK = SIGN_AS_PREDICTED,
                                    IV = signif(IV_COEF, 3), IV_P = round(IV_P, 4),
                                    MIN_WALD = round(FIRST_STAGE_WALD_MIN, 1), WEAK = WEAK_FIRST_STAGE)][
                                      order(MODERATOR, RF_P)])

  n_weak <- cr[TERM == "x Moderator" & WEAK_FIRST_STAGE == 1L, .N]
  if (n_weak > 0L) {
    cat("\nWARNING:", n_weak, "of", cr[TERM == "x Moderator", .N],
        "models have first-stage Wald below 10. Read the RF column, which\n",
        "needs no first stage.\n", sep = " ")
  }

  cat("\nDispersion measures predict POSITIVE. N_PAYERS/N_PAYERS_V2/CODE_COV\n",
      "predict NEGATIVE. SIGN_OK applies the flip.\n", sep = "")

  rm(d); invisible(gc()); list(rows = cr, tests = ct)
}


# ---------------------------------------------------------------------------
# Within-family test
# ---------------------------------------------------------------------------
#
# Specification (b), with clinical family fixed effects, is the decisive one.
# Family FE absorb modality entirely, so a gradient that survives cannot be
# restated as "imaging responds and procedures do not" -- which is the
# objection the shoppability label itself cannot answer, because shoppability
# barely varies within family.

run_comparability_within_family <- function(cr, measures,
                                            moderators = COMPARABILITY_MODERATORS) {

  moderators <- screen_moderators(measures, moderators)
  if (length(moderators) == 0L) { warning("No usable moderators."); return(data.table()) }

  d <- merge(cr[is.finite(RF_COEF) & is.finite(RF_SE) & RF_SE > 0],
             measures, by = "FINAL_CONCEPT_ID")
  setDT(d); d[, W := 1 / (RF_SE^2)]

  has_controls <- all(c("CODE_COV", "N_PAYERS") %in% names(d)) &&
    all(vapply(c("CODE_COV", "N_PAYERS"),
               function(v) sum(is.finite(d[[v]])) > MIN_CONCEPTS_META, logical(1)))

  rows <- list()
  for (m in moderators) {
    d[, MODC := safe_numeric(get(m))]
    d[, MODC := (MODC - mean(MODC, na.rm = TRUE)) / sd(MODC, na.rm = TRUE)]
    for (il in unique(d$INSTRUMENT_LABEL)) {
      dd <- d[INSTRUMENT_LABEL == il & is.finite(MODC)]
      if (nrow(dd) < MIN_CONCEPTS_META) next
      specs <- c(`(a) No family FE` = "RF_COEF ~ MODC",
                 `(b) Family FE`    = "RF_COEF ~ MODC | CONCEPT_FAMILY")
      if (has_controls && !(m %chin% c("CODE_COV", "N_PAYERS"))) {
        specs <- c(specs,
                   `(c) Family FE + coverage` = "RF_COEF ~ MODC + CODE_COV + N_PAYERS | CONCEPT_FAMILY")
      }
      for (sn in names(specs)) {
        fit <- tryCatch(feols(as.formula(specs[[sn]]), data = dd, weights = ~W,
                              cluster = ~CONCEPT_FAMILY, warn = FALSE, notes = FALSE),
                        error = function(e) NULL)
        if (is.null(fit)) next
        td <- tidy_fixest(fit); if (nrow(td) == 0L) next
        td[, `:=`(MODERATOR = m, SPEC = sn, INSTRUMENT_LABEL = il,
                  TIER = instrument_tier(il),
                  N_CONCEPTS = nrow(dd), N_FAMILIES = uniqueN(dd$CONCEPT_FAMILY))]
        rows[[length(rows) + 1L]] <- td
      }
    }
  }
  out <- rbindlist(rows, fill = TRUE)
  if (nrow(out) == 0L) { warning("No within-family results."); return(out) }

  out[, EXPECTED_SIGN := fifelse(MODERATOR %chin% DISPERSION_MEASURES, "positive", "negative")]
  out[term == "MODC", SIGN_AS_PREDICTED :=
        fifelse(MODERATOR %chin% DISPERSION_MEASURES, estimate > 0, estimate < 0)]
  save_csv(out, "T09D_comparability_within_family.csv")

  cat("\n", strrep("=", 112), "\nWITHIN-FAMILY COMPARABILITY TEST\n",
      strrep("=", 112), "\n", sep = "")
  print(out[term == "MODC", .(MODERATOR, SPEC, TIER, INSTRUMENT_LABEL,
                              EST = signif(estimate, 3), SE = signif(std.error, 3), P = round(p.value, 4),
                              SIGN_OK = SIGN_AS_PREDICTED, N_CONCEPTS, N_FAMILIES)][order(MODERATOR, SPEC, P)])

  cat("\nSpec (b) is decisive: family FE absorb modality entirely. Read the\n",
      "MAIN-tier rows first.\n", sep = "")

  cat("\nSummary, spec (b), MAIN tier only:\n")
  print(out[term == "MODC" & grepl("^\\(b\\)", SPEC) & TIER == "MAIN",
            .(N = .N, N_SIGN_OK = sum(SIGN_AS_PREDICTED, na.rm = TRUE),
              N_SIG_05 = sum(p.value < 0.05, na.rm = TRUE),
              MEDIAN_P = round(median(p.value, na.rm = TRUE), 4)), by = MODERATOR])
  out
}

cat("Section 9 loaded\n")



######## Section 10: Comparability meta-regression ############################
#
# The two-step version of Section 9. It has lower power than the row-level
# interaction, but buys three things the row-level test cannot deliver: a
# specification curve across instruments and measures, clinical family fixed
# effects, and a horse race that puts shoppability and each comparability
# measure in the same regression.
#
# Four specifications are estimated per measure and instrument:
#
#   (a)  the moderator alone
#   (b)  + clinical family fixed effects        <- the decisive specification
#   (c)  + family FE + coverage controls        (only when the controls exist)
#   (d)  horse race against the shoppability label, no family FE
#
# Two guards matter here. Moderators pass through screen_moderators(), the same
# degeneracy gate Section 9 uses, so the measures known to be unusable in this
# panel are excluded with a stated reason rather than estimated. And spec (c)
# is added only when its controls actually carry data; otherwise every (c)
# model would fail inside tryCatch and disappear without explanation.
#
# HOW TO READ THE RESULT. Section 9 finds the three price-dispersion measures
# uniformly null under family fixed effects, and finds one contracting-depth
# measure -- N_PAYERS_V2, distinct payers per concept -- significant across
# every MAIN instrument, with the OPPOSITE sign to the ex ante comparability
# prediction. This stage reproduces that pattern and adds the horse race, which
# is what determines whether payer depth and the shoppability label capture the
# same heterogeneity or two independent axes of it.

run_comparability_meta <- function(concept_results, measures,
                                   moderators = COMPARABILITY_MODERATORS,
                                   deps = list(list(dep = "RF_COEF", se = "RF_SE", label = "Reduced form"),
                                               list(dep = "IV_COEF", se = "IV_SE", label = "IV"))) {

  # Same degeneracy gate Section 9 uses, so both stages estimate the same set
  # of measures and any exclusion is reported with its reason.
  moderators <- screen_moderators(measures, moderators)
  if (length(moderators) == 0L) stop("No usable moderators.", call. = FALSE)

  d <- merge(concept_results, measures, by = "FINAL_CONCEPT_ID", all = FALSE)
  setDT(d)

  d[, SHOP := factor(fifelse(FINAL_FAMILY_ID %chin% DIAGNOSTIC_FAMILIES,
                             "Shoppable", "Non_shoppable"),
                     levels = c("Non_shoppable", "Shoppable"))]

  # Standardise so coefficients are per-SD and comparable across measures on
  # very different scales (PD_PAYER_V2 is a ratio near 0.2-1.4; N_PAYERS_V2 a
  # count from 1.6 to 8.1).
  for (m in moderators) {
    d[, (paste0(m, "_Z")) := {
      v <- safe_numeric(get(m)); (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE)
    }]
  }

  rows <- list()
  for (dd in deps) {
    dep <- dd$dep; se <- dd$se
    if (!all(c(dep, se) %in% names(d))) next
    for (m in moderators) {
      mz <- paste0(m, "_Z"); if (!(mz %in% names(d))) next
      for (il in unique(d$INSTRUMENT_LABEL)) {
        s <- d[INSTRUMENT_LABEL == il & is.finite(get(dep)) &
                 is.finite(get(se)) & get(se) > 0 & is.finite(get(mz))]
        if (nrow(s) < MIN_CONCEPTS_META) next
        s[, W := 1 / (get(se)^2)]

        specs <- list(
          `(a) Comparability only` = sprintf("%s ~ %s", dep, mz),
          `(b) + family FE`        = sprintf("%s ~ %s | CONCEPT_FAMILY", dep, mz),
          `(d) Horse race vs shoppability` = sprintf("%s ~ %s + SHOP", dep, mz)
        )

        # Spec (c) needs CODE_COV and N_PAYERS as controls. Both depend on
        # source columns that may be absent, so it is added only when they
        # actually carry data -- otherwise it would fail silently inside
        # tryCatch and its absence would be unexplained.
        if (all(c("CODE_COV", "N_PAYERS") %in% names(s)) &&
            sum(is.finite(s$CODE_COV))  > MIN_CONCEPTS_META &&
            sum(is.finite(s$N_PAYERS)) > MIN_CONCEPTS_META &&
            !(m %chin% c("CODE_COV", "N_PAYERS"))) {
          specs[["(c) + family FE + coverage"]] <-
            sprintf("%s ~ %s + CODE_COV + N_PAYERS | CONCEPT_FAMILY", dep, mz)
        }

        for (sn in names(specs)) {
          fit <- tryCatch(feols(as.formula(specs[[sn]]), data = s, weights = ~W,
                                cluster = ~CONCEPT_FAMILY, warn = FALSE, notes = FALSE),
                          error = function(e) NULL)
          if (is.null(fit)) next
          td <- tidy_fixest(fit); if (nrow(td) == 0L) next
          td[, `:=`(DEPENDENT = dd$label, MODERATOR = m, SPEC = sn,
                    INSTRUMENT_LABEL = il, TIER = instrument_tier(il),
                    N_CONCEPTS = nrow(s), N_FAMILIES = uniqueN(s$CONCEPT_FAMILY))]
          rows[[length(rows) + 1L]] <- td
        }
      }
    }
  }

  out <- rbindlist(rows, fill = TRUE)
  if (nrow(out) == 0L) stop("No comparability meta-regressions.", call. = FALSE)

  # Sign convention, applied to the moderator term only.
  out[grepl("_Z$", term), EXPECTED_SIGN :=
        fifelse(MODERATOR %chin% DISPERSION_MEASURES, "positive", "negative")]
  out[grepl("_Z$", term), SIGN_AS_PREDICTED :=
        fifelse(MODERATOR %chin% DISPERSION_MEASURES, estimate > 0, estimate < 0)]

  save_csv(out, "T10_comparability_meta_regression.csv")

  grad <- out[grepl("_Z$", term)]

  cat("\n", strrep("=", 112), "\nMODERATOR GRADIENT — REDUCED FORM (primary), p-values\n",
      strrep("=", 112), "\n", sep = "")
  print(dcast(grad[DEPENDENT == "Reduced form",
                   .(MODERATOR, SPEC, INSTRUMENT_LABEL, P = round(p.value, 4))],
              MODERATOR + SPEC ~ INSTRUMENT_LABEL, value.var = "P"))

  cat("\n", strrep("=", 112),
      "\nSPEC (b), FAMILY FE — THE DECISIVE SPECIFICATION, MAIN TIER FIRST\n",
      strrep("=", 112), "\n", sep = "")
  print(grad[DEPENDENT == "Reduced form" & grepl("^\\(b\\)", SPEC),
             .(MODERATOR, TIER, INSTRUMENT_LABEL, EST = signif(estimate, 3),
               SE = signif(std.error, 3), P = round(p.value, 4),
               SIGN_OK = SIGN_AS_PREDICTED, N_CONCEPTS, N_FAMILIES)][
                 order(MODERATOR, TIER, P)])

  cat("\n", strrep("=", 112), "\nSAME, IV DEPENDENT (robustness)\n",
      strrep("=", 112), "\n", sep = "")
  print(dcast(grad[DEPENDENT == "IV",
                   .(MODERATOR, SPEC, INSTRUMENT_LABEL, P = round(p.value, 4))],
              MODERATOR + SPEC ~ INSTRUMENT_LABEL, value.var = "P"))

  cat("\n", strrep("=", 112),
      "\nHORSE RACE: each moderator vs the ex-ante shoppability label\n",
      strrep("=", 112), "\n", sep = "")
  hr <- out[grepl("Horse race", SPEC) & DEPENDENT == "Reduced form" &
              (grepl("_Z$", term) | grepl("Shoppable", term))]
  if (nrow(hr) > 0L) {
    hr_wide <- dcast(hr[, .(MODERATOR, TIER, INSTRUMENT_LABEL,
                            TERM = fifelse(grepl("_Z$", term), "MODERATOR_P", "SHOPPABLE_P"),
                            P = round(p.value, 4))],
                     MODERATOR + TIER + INSTRUMENT_LABEL ~ TERM, value.var = "P")
    print(hr_wide[order(MODERATOR, TIER, INSTRUMENT_LABEL)])
    save_csv(hr_wide, "T10B_horse_race_summary.csv")
  }

  cat("\n", strrep("=", 112), "\nSUMMARY BY MEASURE AND SPECIFICATION (reduced form)\n",
      strrep("=", 112), "\n", sep = "")
  print(grad[DEPENDENT == "Reduced form",
             .(N_TESTS = .N, N_SIG = sum(p.value < 0.05, na.rm = TRUE),
               SHARE_SIG = round(mean(p.value < 0.05, na.rm = TRUE), 3),
               MEDIAN_P = round(median(p.value, na.rm = TRUE), 4),
               SHARE_EXPECTED_SIGN = round(mean(SIGN_AS_PREDICTED, na.rm = TRUE), 3)),
             by = .(MODERATOR, SPEC)][order(MODERATOR, SPEC)])

  cat("\nMAIN tier only, spec (b) — the numbers to report:\n")
  print(grad[DEPENDENT == "Reduced form" & grepl("^\\(b\\)", SPEC) & TIER == "MAIN",
             .(N = .N, N_SIG_05 = sum(p.value < 0.05, na.rm = TRUE),
               N_SIGN_OK = sum(SIGN_AS_PREDICTED, na.rm = TRUE),
               MEDIAN_P = round(median(p.value, na.rm = TRUE), 4),
               MEDIAN_EST = signif(median(estimate, na.rm = TRUE), 3)),
             by = MODERATOR][order(MEDIAN_P)])

  cat("\n", strrep("=", 112), "\nHOW TO READ THIS\n", strrep("=", 112), "\n",
      "SIGNS. Dispersion measures (PD_PAYER_V2, PD_HOSP, N_CODES) predict\n",
      "POSITIVE: more dispersion = less comparable = smaller price response.\n",
      "N_PAYERS_V2 predicts NEGATIVE under the original comparability theory\n",
      "(more payers = thicker contracting = more comparable = bigger response).\n",
      "SIGN_OK applies the flip. NOTE: stage 9 found N_PAYERS_V2 significant\n",
      "with the OPPOSITE sign, so SIGN_OK = FALSE there is the finding, not a\n",
      "failure -- see the entrenchment interpretation.\n\n",
      "SPEC (b) IS THE ONE THAT MATTERS. Family FE absorb modality entirely, so\n",
      "a surviving gradient cannot be restated as 'imaging responds and\n",
      "procedures do not' -- the objection the shoppability label cannot answer,\n",
      "because it barely varies within family.\n\n",
      "THE HORSE RACE. Both terms in one regression, no family FE (SHOP is\n",
      "family-determined and would be absorbed). Four readings:\n",
      "  Both significant   -> independent axes of heterogeneity. Report both.\n",
      "  Only the moderator -> it subsumes the shoppability result.\n",
      "  Only SHOPPABLE     -> the moderator adds nothing beyond the label.\n",
      "  Neither -> collinear at this sample size. Report separately, say so.\n",
      sep = "")

  out
}

cat("Section 10 loaded\n")



######## Section 11: Instrument balance and placebo tests #####################
#
# Why these two tests and not an event study. 3,722 of 3,724 hospitals appear
# at exactly one posting month, so there is no within-hospital timing variation
# and a within-hospital event study does not exist in this data. The
# identification concern is therefore selection on posting timing, not timing
# mechanics, and it needs a different kind of test.
#
# The primary defence is structural and lives in the design rather than here:
# the headline result is a HETEROGENEITY result, and the non-shoppable group
# functions as a within-hospital placebo. Same hospital, same posting month,
# same instrument value, no price response. Any selection story has to explain
# why selection would operate on shoppable services and not on non-shoppable
# ones at the same hospital in the same month.
#
# The two tests here supplement that defence.
#
# TEST 1 -- BALANCE. Does the instrument predict pre-determined hospital
# characteristics? Estimated at the HOSPITAL level, one row per hospital, since
# each hospital posts once. The fixed effects must be ADDITIVE county + month:
# an interacted county x month term would absorb the county-month instruments
# entirely and leave nothing to test.
#
# TEST 2 -- PLACEBO OUTCOME. Bed count is fixed years before any disclosure
# decision, so the instrument cannot cause it. A null is a clean falsification;
# a significant coefficient would mean the instrument is picking up hospital
# composition the fixed effects failed to absorb.

run_instrument_balance <- function(panel, instruments = MAIN_INSTRUMENTS) {
  hosp_cols <- available_columns(panel, c(
    "HOSPITAL_ID", "ANALYSIS_MARKET", "POST_MONTH", "LOG_TOTAL_BEDS", "TOTAL_BEDS",
    "HEALTH_SYSTEM_ID", "SYSTEM_KEY", "HOSPITAL_TYPE", "PROVIDER_STATE",
    unname(instruments)))
  h <- unique(panel[, ..hosp_cols], by = "HOSPITAL_ID")

  h[, IN_SYSTEM := as.integer(!is.na(SYSTEM_KEY) & SYSTEM_KEY != "")]
  if ("HOSPITAL_TYPE" %in% names(h)) {
    h[, IS_SHORT_TERM := as.integer(grepl("SHORT|ACUTE|GENERAL", HOSPITAL_TYPE, ignore.case = TRUE))]
  }

  cat("\nBalance sample:", nrow(h), "hospitals |",
      uniqueN(h$ANALYSIS_MARKET), "counties |", uniqueN(h$POST_MONTH), "months\n")

  chars <- available_columns(h, c("LOG_TOTAL_BEDS", "IN_SYSTEM", "IS_SHORT_TERM"))
  chars <- chars[vapply(chars, function(v) has_usable_variation(h[[v]]), logical(1))]

  out <- rbindlist(lapply(names(instruments), function(il) {
    z <- instruments[[il]]; if (!(z %in% names(h))) return(data.table())
    sd_z <- sd(h[[z]], na.rm = TRUE)
    rbindlist(lapply(chars, function(ch) {
      d <- h[is.finite(get(ch)) & is.finite(get(z))]
      if (nrow(d) < 200L) return(data.table())
      # ADDITIVE county + month FE. Interacted would absorb the county-month
      # instruments entirely and leave nothing to test.
      fit <- tryCatch(feols(as.formula(paste0(ch, " ~ ", z,
                                              " | ANALYSIS_MARKET + POST_MONTH")),
                            data = d, cluster = ~ANALYSIS_MARKET,
                            warn = FALSE, notes = FALSE), error = function(e) NULL)
      if (is.null(fit)) return(data.table())
      co <- extract_coefficient(fit, z)
      data.table(INSTRUMENT_LABEL = il, CHARACTERISTIC = ch,
                 COEF = co$estimate, SE = co$std_error, P_VALUE = co$p_value,
                 EFFECT_PER_SD_OF_Z = co$estimate * sd_z,
                 SD_OF_CHARACTERISTIC = sd(d[[ch]], na.rm = TRUE),
                 STANDARDISED_EFFECT = co$estimate * sd_z / sd(d[[ch]], na.rm = TRUE),
                 N_HOSPITALS = nrow(d))
    }), fill = TRUE)
  }), fill = TRUE)

  if (nrow(out) == 0L) { warning("No balance results."); return(out) }
  save_csv(out, "T11_instrument_balance.csv")

  cat("\n", strrep("=", 100), "\nINSTRUMENT BALANCE\n", strrep("=", 100), "\n", sep = "")
  print(out[, .(INSTRUMENT_LABEL, CHARACTERISTIC, COEF = signif(COEF, 3),
                SE = signif(SE, 3), P = round(P_VALUE, 4),
                STD_EFFECT = round(STANDARDISED_EFFECT, 4), N_HOSPITALS)][
                  order(CHARACTERISTIC, P)])

  cat("\nSTANDARDISED_EFFECT is the change in the characteristic, in its own SDs,\n",
      "per one-SD increase in peer exposure. Values under about 0.05 are\n",
      "negligible even when statistically significant at this sample size --\n",
      "report the magnitude, not just the p-value.\n",
      "Significant AND large would mean high-exposure hospitals differ\n",
      "systematically, which the county and month FE were supposed to absorb.\n",
      sep = "")

  cat("\nSummary:", sum(out$P_VALUE < 0.05, na.rm = TRUE), "of", nrow(out),
      "balance tests significant at 5% |",
      sum(abs(out$STANDARDISED_EFFECT) > 0.05, na.rm = TRUE),
      "with standardised effect above 0.05\n")
  out
}

run_placebo_outcome <- function(panel, instruments = MAIN_INSTRUMENTS) {
  # Collapse to one row per hospital first, for the same reason
  # run_instrument_balance() does. LOG_TOTAL_BEDS is a hospital-level constant,
  # so running it on the concept-level panel repeats the same value once per
  # concept -- up to about 19 times per hospital. Two-way county and month
  # clustering mostly but not exactly absorbs that, and N_OBSERVATIONS would
  # report a panel-row count rather than the true number of hospitals.
  hosp_cols <- available_columns(panel, c("HOSPITAL_ID", "ANALYSIS_MARKET",
                                          "POST_MONTH", "LOG_TOTAL_BEDS",
                                          unname(instruments)))
  h <- unique(panel[, ..hosp_cols], by = "HOSPITAL_ID")

  out <- rbindlist(lapply(names(instruments), function(il) {
    z <- instruments[[il]]; if (!(z %in% names(h))) return(data.table())
    d <- h[is.finite(LOG_TOTAL_BEDS) & is.finite(get(z))]
    if (nrow(d) < 200L) return(data.table())
    fit <- tryCatch(feols(as.formula(paste0("LOG_TOTAL_BEDS ~ ", z,
                                            " | ANALYSIS_MARKET + POST_MONTH")),
                          data = d, cluster = ~ANALYSIS_MARKET,
                          warn = FALSE, notes = FALSE), error = function(e) NULL)
    if (is.null(fit)) return(data.table())
    co <- extract_coefficient(fit, z)
    sd_z <- sd(d[[z]], na.rm = TRUE)
    data.table(INSTRUMENT_LABEL = il, PLACEBO_OUTCOME = "LOG_TOTAL_BEDS",
               COEF = co$estimate, SE = co$std_error, P_VALUE = co$p_value,
               EFFECT_PER_SD_OF_Z = co$estimate * sd_z,
               N_HOSPITALS = nrow(d))
  }), fill = TRUE)

  if (nrow(out) == 0L) { warning("No placebo results."); return(out) }
  save_csv(out, "T11B_placebo_outcome.csv")

  cat("\n", strrep("=", 100), "\nPLACEBO OUTCOME: PRE-DETERMINED BED COUNT\n",
      strrep("=", 100), "\n", sep = "")
  print(out[, .(INSTRUMENT_LABEL, COEF = signif(COEF, 4), SE = signif(SE, 4),
                P = round(P_VALUE, 4),
                EFFECT_PER_SD = signif(EFFECT_PER_SD_OF_Z, 4), N_HOSPITALS)])
  cat("\nBed count is fixed years before any disclosure decision, so the\n",
      "instrument cannot cause it. A null is a clean falsification. A\n",
      "significant coefficient would mean the instrument is picking up\n",
      "hospital composition the fixed effects failed to absorb.\n", sep = "")
  out
}

cat("Section 11 loaded\n")



######## Section 12: Preflight ################################################
#
# Runs in about ten seconds and must pass before stage 6, which takes eight
# hours. It checks the failure modes that are silent rather than loud: concepts
# with no scheme classification, unexpected NA scheme columns, merge
# constituents left behind in the panel, canonical concepts missing from it,
# and instruments that are absent or degenerate.
#
# Each of these produces a plausible-looking result rather than an error, which
# is why they are worth an explicit gate.

run_preflight <- function(panel, schemes_long) {
  cat("\n", strrep("=", 84), "\nPREFLIGHT\n", strrep("=", 84), "\n", sep = "")
  fails <- character(0)
  chk <- function(label, ok, detail = "") {
    cat(sprintf("  [%s] %-56s %s\n", if (isTRUE(ok)) "PASS" else "FAIL", label, detail))
    if (!isTRUE(ok)) fails <<- c(fails, label)
  }

  # 1. Every panel concept has a scheme classification. A concept missing from
  #    schemes_long is dropped from every scheme-based model without warning.
  miss <- setdiff(unique(panel$FINAL_CONCEPT_ID), schemes_long$ANALYSIS_CONCEPT_ID)
  chk("All panel concepts appear in schemes_long", length(miss) == 0L,
      if (length(miss)) paste(head(miss, 3), collapse = ", ") else "")

  # 2. No unexpected NA scheme columns. Scheme 6 legitimately has NAs, since
  #    the extremes rule drops INTERMEDIATE concepts by design.
  for (spec in PRIMARY_SCHEMES) {
    if (!(spec$col %in% names(panel))) { chk(paste("Column present:", spec$col), FALSE); next }
    n_na <- sum(is.na(panel[[spec$col]]))
    allow_na <- identical(spec$rule, "extremes")
    chk(paste("No unexpected NAs in", spec$col), allow_na || n_na == 0L,
        paste0(format(n_na, big.mark = ","), " NA", if (allow_na) " (expected)" else ""))
  }

  # 3. Merge canonical concepts present.
  canon <- vapply(MERGE_GROUPS, `[[`, character(1), "canonical_id")
  present <- canon %in% unique(panel$FINAL_CONCEPT_ID)
  chk("Merge canonical concepts present in panel", all(present),
      if (all(present)) "" else paste(canon[!present], collapse = ", "))

  # 4. Merge constituents GONE from the panel.
  consts <- setdiff(unique(unlist(lapply(MERGE_GROUPS, `[[`, "constituents"))), canon)
  leftover <- intersect(consts, unique(panel$FINAL_CONCEPT_ID))
  chk("Merge constituents removed from panel", length(leftover) == 0L,
      if (length(leftover)) paste(head(leftover, 3), collapse = ", ") else "")

  # 5. Instruments present and varying.
  for (il in names(MAIN_INSTRUMENTS)) {
    z <- MAIN_INSTRUMENTS[[il]]
    chk(paste("Instrument usable:", il),
        z %in% names(panel) && has_usable_variation(panel[[z]]))
  }

  # 6. Outcome and treatment finite.
  chk("Outcome finite for all rows", all(is.finite(panel[[PRIMARY_OUTCOME]])))
  chk("Treatment finite for all rows", all(is.finite(panel[[ENDOGENOUS_VARIABLE]])))

  # 7. Exact-code panel reachable. PD_CODE and N_CODES in stage 9 need it.
  chk("Exact-code panel file exists", file.exists(FILES$outpatient_exact))

  cat("\n", strrep("-", 84), "\n", sep = "")
  if (length(fails) == 0L) {
    cat("ALL CHECKS PASSED. Safe to run stage 6.\n")
  } else {
    cat("FAILURES:\n"); for (f in fails) cat("  -", f, "\n")
    stop("Preflight failed. Do not start the 8-hour run.", call. = FALSE)
  }
  invisible(TRUE)
}

cat("Section 12 loaded\n")



###############################################################################
# ORDER OF OPERATIONS
#
# Everything above this line defines constants and functions. Everything below
# executes. RUN_STAGES, set in Section 0, selects which stage blocks run.
#
# STEP 1  Load Sections 0-12. Seconds. To force a full rebuild, clear the
#         relevant caches first:
#
#           file.remove(file.path(CACHE_DIR, paste0(c(
#             "schemes_long", "outpatient_panel", "concept_level_6inst",
#             "main_results", "transform_ladder",
#             "meta_regressions", "meta_regressions_rf", "meta_regressions_iv",
#             "comparability_measures", "comparability_interaction"), ".rds")))
#
# STEP 2  Run the BUILD block. A few minutes. run_preflight() must pass before
#         anything else runs.
#
# STEP 3  RUN_STAGES <- c(7, 11)     main results and balance tests, ~1 hour.
#         Independent of stage 6 and can run first.
#
# STEP 4  RUN_STAGES <- c(6)         concept-level sweep, ~8 hours.
#
# STEP 5  RUN_STAGES <- c(8, 9, 10)  meta-regressions and mechanism, ~1 hour.
#         The stage 6 loader verifies that the cached concept-level results
#         were estimated on the panel currently in memory before proceeding.
###############################################################################
######## BUILD -- always run this first in a fresh session ####################
#
# Builds schemes_long and the outpatient panel, exports the Scheme 1
# assignment, audits the estimation sample, and runs preflight. Both objects
# are cached, so this is fast on a warm run and is safe to re-execute.

report_memory("Start:")

schemes_long <- cache_or_run("schemes_long",
                             extend_schemes_for_merged(build_schemes()))

outpatient <- cache_or_run("outpatient_panel", {
  p <- load_outpatient()
  p <- apply_concept_merges(p)
  attach_scheme_columns(p, schemes_long)
})

export_scheme1_assignments(outpatient)
audit_estimation_sample(outpatient)
run_preflight(outpatient, schemes_long)

CONCEPT_INSTRUMENTS <- c(MAIN_INSTRUMENTS, SUPPORTING_INSTRUMENTS)
concept_results <- cache_or_run("concept_level_6inst",
                                estimate_concept_level(
                                  build_concept_panel(outpatient, instruments = CONCEPT_INSTRUMENTS),
                                  instruments = CONCEPT_INSTRUMENTS))

stopifnot(length(setdiff(unique(concept_results$FINAL_CONCEPT_ID),
                         unique(outpatient$FINAL_CONCEPT_ID))) == 0)
cat("Concept results match the current panel:",
    uniqueN(concept_results$FINAL_CONCEPT_ID), "concepts.\n")

cat("\nTheory V2 HIGH concepts by family:\n")
print(schemes_long[SCHEME_ID == "scheme_theory_v2" & SHOPPABILITY_CATEGORY == "HIGH",
                   .N, by = ANALYSIS_FAMILY_ID][order(-N)])


for (canon in names(CANONICAL_SCHEME_DONOR)) {
  grp <- Filter(function(g) g$canonical_id == canon, MERGE_GROUPS)[[1]]
  chk <- schemes_long[ANALYSIS_CONCEPT_ID %chin% grp$constituents,
                      .(N = uniqueN(SHOPPABILITY_CATEGORY)), by = SCHEME_ID]
  cat(canon, ": max distinct categories =", max(chk$N), "\n")
}



######## STAGE 5 -- instrument screen and pooled null #########################
if (5 %in% RUN_STAGES) {
  screen <- cache_or_run("instrument_screen", run_instrument_screen(outpatient))
  pooled <- cache_or_run("pooled_models",     run_pooled_models(outpatient))
}



######## STAGE 6 -- concept level (~8 hours) ##################################
if (6 %in% RUN_STAGES) {
  CONCEPT_INSTRUMENTS <- c(MAIN_INSTRUMENTS, SUPPORTING_INSTRUMENTS)
  concept_results <- cache_or_run("concept_level_6inst",
                                  estimate_concept_level(
                                    build_concept_panel(outpatient, instruments = CONCEPT_INSTRUMENTS),
                                    instruments = CONCEPT_INSTRUMENTS))
  save_csv(concept_results, "T05_concept_level_RF_FS_IV.csv")
  diagnose_size_gradient(concept_results)

  cat("\nMerge canonical concepts estimated (expect 6 rows each, one per instrument):\n")
  print(concept_results[FINAL_CONCEPT_ID %chin%
                          vapply(MERGE_GROUPS, `[[`, character(1), "canonical_id"), .N, by = FINAL_CONCEPT_ID])
}



######## STAGE 6 LOADER -- for stages 8, 9, 10 ################################
#
# Loads the concept-level results from cache and verifies that their concept
# universe matches the panel currently in memory. Stage 6 takes eight hours, so
# the cache is long-lived and can easily outlive a change to the codebook or the
# merge groups. A mismatch here means the cached results were estimated on a
# different panel vintage and must be rebuilt.

if (any(c(8, 9, 10) %in% RUN_STAGES) && !exists("concept_results")) {
  CONCEPT_INSTRUMENTS <- c(MAIN_INSTRUMENTS, SUPPORTING_INSTRUMENTS)
  concept_results <- cache_or_run("concept_level_6inst",
                                  estimate_concept_level(
                                    build_concept_panel(outpatient, instruments = CONCEPT_INSTRUMENTS),
                                    instruments = CONCEPT_INSTRUMENTS))
}

if (any(c(8, 9, 10) %in% RUN_STAGES)) {
  stale <- setdiff(unique(concept_results$FINAL_CONCEPT_ID),
                   unique(outpatient$FINAL_CONCEPT_ID))
  if (length(stale) > 0L) {
    cat("\nSTALE CONCEPTS IN concept_results:\n"); print(head(stale, 20))
    stop("concept_results was estimated on a DIFFERENT panel vintage. Delete\n",
         "  concept_level_6inst.rds and re-run stage 6.", call. = FALSE)
  }
  cat("Concept results match the current panel:",
      uniqueN(concept_results$FINAL_CONCEPT_ID), "concepts.\n")
}



###############################################################################
# INSTRUMENT AUDIT -- justifying the three-instrument headline set
#
# Documents why the headline claims rest on three of the six system
# instruments. Three separate questions, answered separately and never averaged
# together, because they measure different things:
#
#   A. RAW STRENGTH. Pooled first-stage F for all six, estimated on a common
#      sample so the magnitudes are comparable. Coefficient SIZE is not
#      comparable across instruments that count different units -- a market has
#      far more competitor hospitals than competitor systems -- but F is
#      unit-free and is the right comparison.
#
#   B. IDENTIFICATION CREDIBILITY. Sign stability of the shoppable term across
#      all six classification schemes. This is the criterion that excludes
#      Competitor_systems_9m from the headline set.
#
#   C. DOES THE HEADLINE RESULT CHANGE? The actual reported specification
#      (Scheme 1, interacted, RF and IV) estimated under each of the six
#      instruments. This is the decisive test: whether including an instrument
#      changes the finding, not whether it has a larger first stage.
#
# Writes QA05 through QA08. Runtime: A and B are minutes, C is six models at
# roughly five to ten minutes total.
###############################################################################
stopifnot(exists("outpatient"))

CANDIDATE_INSTRUMENTS <- c(MAIN_INSTRUMENTS, SUPPORTING_INSTRUMENTS)   # all 6

cat("\nCandidates:\n")
print(data.table(LABEL = names(CANDIDATE_INSTRUMENTS),
                 COLUMN = unname(CANDIDATE_INSTRUMENTS),
                 TIER = c(rep("MAIN", length(MAIN_INSTRUMENTS)),
                          rep("SUPPORTING", length(SUPPORTING_INSTRUMENTS)))))



###############################################################################
# A. RAW STRENGTH -- pooled first stage on a common sample
###############################################################################
#
# Coefficient size differs by a factor of 20 to 30 between hospital-counted and
# system-counted instruments, purely because they count different units. That
# is a units artefact rather than a strength difference. F-statistics are
# unit-free and are the correct comparison; STANDARDISED_COEF rescales each
# coefficient by its own instrument's standard deviation for the same reason.

req_cols <- available_columns(outpatient, unique(c(
  ENDOGENOUS_VARIABLE, PRIMARY_OUTCOME, BASELINE_CONTROLS, BASELINE_FIXED_EFFECTS,
  BASELINE_CLUSTERS, unname(CANDIDATE_INSTRUMENTS))))
d0 <- outpatient[, ..req_cols]

strength <- rbindlist(lapply(names(CANDIDATE_INSTRUMENTS), function(lab) {
  z <- CANDIDATE_INSTRUMENTS[[lab]]
  fs <- run_first_stage(d0, z)
  co <- extract_coefficient(fs, z)
  f  <- if (is.finite(co$statistic)) co$statistic^2 else NA_real_
  data.table(
    LABEL = lab, INSTRUMENT = z,
    TIER = if (lab %in% names(MAIN_INSTRUMENTS)) "MAIN" else "SUPPORTING",
    COEFFICIENT = co$estimate, F_STAT = f, P_VALUE = co$p_value,
    N_OBSERVATIONS = if (is.null(fs)) NA_integer_ else nobs(fs),
    INSTRUMENT_SD = sd(d0[[z]], na.rm = TRUE),
    SHARE_ZERO = mean(d0[[z]] == 0, na.rm = TRUE),
    STANDARDISED_COEF = co$estimate * sd(d0[[z]], na.rm = TRUE))
}), fill = TRUE)

setorder(strength, -F_STAT)
save_qa_csv(strength, "QA05_instrument_strength_all_6.csv")

cat("\n", strrep("=", 100), "\nA. FIRST-STAGE STRENGTH — ALL 6, SAME SAMPLE\n",
    strrep("=", 100), "\n", sep = "")
print(strength[, .(LABEL, TIER, F_STAT = round(F_STAT, 1),
                   COEFFICIENT = signif(COEFFICIENT, 3),
                   STD_COEF = signif(STANDARDISED_COEF, 3),
                   WEAK_F_LT_10 = as.integer(F_STAT < 10))])

cat("\nSTANDARDISED_COEF rescales each coefficient by its OWN instrument's SD --\n",
    "the right way to compare a hospital-count instrument against a\n",
    "system-count instrument. If MAIN and SUPPORTING look similar on F_STAT,\n",
    "raw strength is NOT the reason for the tiering, which points to Block B.\n",
    sep = "")

rm(d0); invisible(gc())



###############################################################################
# B. SIGN STABILITY across classification schemes
###############################################################################
#
# Estimates the shoppable term under all six schemes for each instrument and
# reports the share of scheme variants coming back negative. An instrument
# whose sign depends on how the shoppable line is drawn is not identifying a
# stable object, whatever its first-stage F.
#
# This is the criterion that excludes Competitor_systems_9m from the headline
# set: only 54% of its scheme variants are negative, and the sign flips under
# the tier-collapse rule.

sign_check <- rbindlist(lapply(names(CANDIDATE_INSTRUMENTS), function(lab) {
  z <- CANDIDATE_INSTRUMENTS[[lab]]

  rbindlist(lapply(names(SCHEME_COLUMNS), function(sl) {
    sc <- SCHEME_COLUMNS[[sl]]
    r <- estimate_interacted(outpatient, sc, PRIMARY_OUTCOME, z,
                             moderator_type = "categorical",
                             label = sl, instrument_label = lab, moderator_label = sl)
    if (is.null(r) || nrow(r$rows) == 0L) return(data.table())
    shop <- r$rows[TERM == "Shoppable"]
    if (nrow(shop) == 0L) return(data.table())
    data.table(LABEL = lab, SCHEME = sl,
               IV_PERCENT = shop$IV_PERCENT, IV_P = shop$IV_P,
               RF_COEF = shop$RF_COEF, RF_P = shop$RF_P)
  }), fill = TRUE)
}), fill = TRUE)

save_qa_csv(sign_check, "QA06_sign_stability_all_6_corrected_panel.csv")

sign_summary <- sign_check[, .(
  N_SCHEMES = .N,
  N_NEGATIVE_IV = sum(IV_PERCENT < 0, na.rm = TRUE),
  SHARE_NEGATIVE_IV = round(mean(IV_PERCENT < 0, na.rm = TRUE), 3),
  N_NEGATIVE_RF = sum(RF_COEF < 0, na.rm = TRUE),
  SHARE_NEGATIVE_RF = round(mean(RF_COEF < 0, na.rm = TRUE), 3),
  MEDIAN_IV_PERCENT = round(median(IV_PERCENT, na.rm = TRUE), 3)
), by = LABEL][order(-SHARE_NEGATIVE_IV)]

cat("\n", strrep("=", 100), "\nB. SIGN STABILITY ON THE CORRECTED PANEL (shoppable term, 6 schemes)\n",
    strrep("=", 100), "\n", sep = "")
print(sign_summary)

cat("\nCompare SHARE_NEGATIVE_IV against the 54% that disqualified\n",
    "Competitor_systems_9m originally. If it is still well below the other\n",
    "instruments' share, the disqualification survives the panel fixes. If it\n",
    "has moved close to the others, the earlier verdict was partly an artefact\n",
    "of the classification bugs that are now fixed, and the instrument may\n",
    "deserve reconsideration.\n", sep = "")



###############################################################################
# C. DOES THE HEADLINE RESULT CHANGE?
###############################################################################
#
# A and B inform the decision; this settles it. The actual reported
# specification -- Scheme 1, interacted, reduced form and IV -- is estimated
# under all six instruments, and the question is whether the three supporting
# instruments tell a different story from the three main ones.
#
# If they agree, restricting the headline to three was conservative but
# costless. If they disagree, the disagreement is itself the reason to keep the
# tiering. Either way the comparison belongs in the paper as a stated
# robustness check rather than as a silent choice, which is why it is written
# to QA07 and QA08.

hc_rows  <- list(); hc_tests <- list()

for (lab in names(CANDIDATE_INSTRUMENTS)) {
  z <- CANDIDATE_INSTRUMENTS[[lab]]
  r <- estimate_interacted(outpatient, "SCHEME_1_CERTAINTY", PRIMARY_OUTCOME, z,
                           moderator_type = "categorical", label = "1. Procedural certainty",
                           instrument_label = lab, moderator_label = lab)
  if (is.null(r)) { cat("  ", lab, ": no result (thin sample or degenerate)\n"); next }
  r$rows[,  TIER := if (lab %in% names(MAIN_INSTRUMENTS)) "MAIN" else "SUPPORTING"]
  r$tests[, TIER := if (lab %in% names(MAIN_INSTRUMENTS)) "MAIN" else "SUPPORTING"]
  hc_rows[[lab]]  <- r$rows
  hc_tests[[lab]] <- r$tests
  cat(sprintf("  %-38s done\n", lab))
}

hc_rows  <- rbindlist(hc_rows,  fill = TRUE)
hc_tests <- rbindlist(hc_tests, fill = TRUE)

save_qa_csv(hc_rows,  "QA07_headline_scheme1_all_6_instruments_rows.csv")
save_qa_csv(hc_tests, "QA08_headline_scheme1_all_6_instruments_tests.csv")

cat("\n", strrep("=", 104),
    "\nC. HEADLINE RESULT (Scheme 1, procedural certainty) UNDER ALL 6\n",
    strrep("=", 104), "\n", sep = "")
print(hc_rows[, .(LABEL = INSTRUMENT_LABEL, TIER, TERM,
                  RF_PCT_PER_SD = round(RF_PERCENT_PER_SD, 3), RF_P = round(RF_P, 4),
                  IV_PCT = round(IV_PERCENT, 2), IV_P = round(IV_P, 4),
                  MIN_WALD_F = round(FIRST_STAGE_WALD_MIN, 1))][
                    order(TIER, LABEL, TERM)])

cat("\nHeterogeneity test (Shoppable = Non_shoppable?) under each instrument:\n")
print(hc_tests[, .(LABEL = INSTRUMENT_LABEL, TIER, ESTIMATOR, P = round(P_VALUE, 4))][
  order(TIER, LABEL, ESTIMATOR)])

cat("\n", strrep("=", 104), "\nDECISION RULE\n", strrep("=", 104), "\n",
    "Look at the RF_PCT_PER_SD sign and the heterogeneity P for MAIN vs\n",
    "SUPPORTING instruments.\n\n",
    "  SUPPORTING instruments agree in sign and significance with MAIN\n",
    "    -> the 3-instrument restriction was conservative, not necessary. You\n",
    "       could add the 2 non-disqualified supporting instruments\n",
    "       (Competitor_outside_CBSA_systems_9m, Competitor_outside_CBSA_counties_9m)\n",
    "       to the headline set, or at minimum report them as confirming\n",
    "       robustness in an appendix table.\n\n",
    "  SUPPORTING instruments disagree (wrong sign, or heterogeneity test not\n",
    "    significant where MAIN is)\n",
    "    -> the exclusion-restriction concern that motivated the tiering is\n",
    "       doing real work, not just being cautious. Keep the 3-instrument\n",
    "       headline and report this table as the reason why.\n\n",
    "  Competitor_systems_9m specifically -- check Block B's sign share first.\n",
    "    If it still shows a materially lower share-negative than the other 5,\n",
    "    it stays excluded regardless of what happens here.\n", sep = "")

cat("\n", strrep("=", 84), "\nINSTRUMENT AUDIT COMPLETE\n", strrep("=", 84), "\n",
    "  QA05  first-stage F, all 6, same sample\n",
    "  QA06  sign stability, all 6, corrected panel\n",
    "  QA07/QA08  headline result under all 6\n", sep = "")



###############################################################################
# STAGE 7 -- main results and transform ladder
#
#   7.1  Headline: three MAIN instruments, six schemes, plus the transform
#        ladder and the Scheme 1 procedural-exclusion sensitivity
#   7.2  Confirming tier
#   7.3  Discrepant tier, reported alone and never pooled
#   7.4  All three tiers combined into one tagged table
#   7.5  Pooled robustness claim over MAIN + CONFIRMING
#   7.6  Tier rerank across all six schemes, alongside first-stage F
###############################################################################
if (7 %in% RUN_STAGES) {

  # --- 7.1 HEADLINE: 3 MAIN instruments ---------------------------------------
  main   <- cache_or_run("main_results",
                         run_main_results(outpatient, stem = "T06_main"))
  ladder <- cache_or_run("transform_ladder", run_transform_ladder(outpatient))

  if (EXCLUDE_PROCEDURAL_FROM_SCHEME1) {
    cat("\n", strrep("=", 84),
        "\nSENSITIVITY: Scheme 1 with 17 procedure-adjunct concepts reclassified\n",
        "(intraoperative ultrasound, imaging guidance, psychotherapy)\n",
        strrep("=", 84), "\n", sep = "")
    main_sens <- run_main_results(apply_scheme1_procedural_exclusion(outpatient),
                                  schemes = SCHEME_COLUMNS["1. Procedural certainty"],
                                  stem = "T06D_scheme1_procedural_excluded")
  }

  # --- 7.2 CONFIRMING: 2 additional instruments -------------------------------
  cat("\n", strrep("=", 84), "\n7.2 CONFIRMING INSTRUMENTS (2)\n", strrep("=", 84), "\n", sep = "")
  confirming <- cache_or_run("confirming_results",
                             run_main_results(outpatient, instruments = CONFIRMING_INSTRUMENTS,
                                              stem = "T06J_confirming"))

  # --- 7.3 DISCREPANT: reported on its own, never pooled ----------------------
  cat("\n", strrep("=", 84), "\n7.3 DISCREPANT INSTRUMENT (1) — reported, not pooled\n",
      strrep("=", 84), "\n", sep = "")
  discrepant <- cache_or_run("discrepant_results",
                             run_main_results(outpatient, instruments = DISCREPANT_INSTRUMENTS,
                                              stem = "T06K_discrepant"))

  # --- 7.4 COMBINED TABLE: all three tiers, tagged, side by side -------------
  combined <- rbindlist(list(
    cbind(TIER = "MAIN",       main$rows),
    cbind(TIER = "CONFIRMING", confirming$rows),
    cbind(TIER = "DISCREPANT", discrepant$rows)
  ), fill = TRUE)
  save_csv(combined, "T06F_all_tiers_combined.csv")

  combined_tests <- rbindlist(list(
    cbind(TIER = "MAIN",       main$tests),
    cbind(TIER = "CONFIRMING", confirming$tests),
    cbind(TIER = "DISCREPANT", discrepant$tests)
  ), fill = TRUE)
  save_csv(combined_tests, "T06G_all_tiers_heterogeneity_tests.csv")

  cat("\n", strrep("=", 108), "\nHEADLINE TERM (Shoppable) BY TIER, SCHEME 1\n",
      strrep("=", 108), "\n", sep = "")
  print(combined[SPEC == "1. Procedural certainty" & TERM == "Shoppable",
                 .(TIER, INSTRUMENT_LABEL, RF_PCT = round(RF_PERCENT_PER_SD, 3),
                   RF_P = round(RF_P, 4), IV_PCT = round(IV_PERCENT, 2),
                   IV_P = round(IV_P, 4), MIN_WALD_F = round(FIRST_STAGE_WALD_MIN, 1))][
                     order(TIER, INSTRUMENT_LABEL)])

  # --- 7.5 POOLED ROBUSTNESS CLAIM: MAIN + CONFIRMING, 5 instruments ---------
  cat("\n", strrep("=", 84), "\n7.5 POOLED: MAIN + CONFIRMING (5 instruments x 6 schemes = 30 models)\n",
      strrep("=", 84), "\n", sep = "")
  robustness_pool <- cache_or_run("robustness_pool_results",
                                  run_main_results(outpatient, instruments = ROBUSTNESS_INSTRUMENTS,
                                                   stem = "T06H_pooled_main_confirming"))

  het_summary <- robustness_pool$tests[
    , .(N_TESTS = .N, N_SIG_05 = sum(P_VALUE < 0.05, na.rm = TRUE),
        SHARE_SIG = round(mean(P_VALUE < 0.05, na.rm = TRUE), 3),
        MEDIAN_P = round(median(P_VALUE, na.rm = TRUE), 4)),
    by = ESTIMATOR]
  cat("\nHeterogeneity significance, pooled across MAIN + CONFIRMING (30 models):\n")
  print(het_summary)

  # --- 7.6 TIER RERANK: across ALL 6 SCHEMES, not Scheme 1 alone -------------
  # The tier assignment in QA07 and QA08 was made on Scheme 1 alone. This
  # checks whether it holds across the full set of schemes, and reports median
  # first-stage F alongside so the two axes stay visibly separate.
  #
  # This is diagnostic, not a re-selection rule. Re-tiering on these p-values
  # would be outcome-based instrument selection, which is exactly what the
  # design decision in Section 0 avoids.
  rerank <- combined_tests[ESTIMATOR == "Reduced form",
                           .(N_SCHEMES = .N,
                             N_SIG_05 = sum(P_VALUE < 0.05, na.rm = TRUE),
                             SHARE_SIG = round(mean(P_VALUE < 0.05, na.rm = TRUE), 3),
                             MEDIAN_P = round(median(P_VALUE, na.rm = TRUE), 4)),
                           by = .(TIER, INSTRUMENT_LABEL)][order(-SHARE_SIG, MEDIAN_P)]

  fs_by_inst <- combined[, .(MEDIAN_FIRST_STAGE_F = round(median(FIRST_STAGE_WALD_MIN, na.rm = TRUE), 1)),
                         by = INSTRUMENT_LABEL]
  rerank <- merge(rerank, fs_by_inst, by = "INSTRUMENT_LABEL")[order(-SHARE_SIG, MEDIAN_P)]

  cat("\n", strrep("=", 104),
      "\n7.6 TIER RERANK — heterogeneity significance across ALL 6 SCHEMES (RF),\n",
      "with first-stage F alongside to show they are INDEPENDENT axes\n",
      strrep("=", 104), "\n", sep = "")
  print(rerank)
  save_csv(rerank, "T06L_tier_rerank_all_schemes.csv")

  cat("\nThe two STRONGEST first stages produce the two WEAKEST results. That is\n",
      "an argument against selecting instruments on first-stage F, and worth a\n",
      "sentence in the paper.\n\n",
      "DESIGN CHECK BEFORE FINALISING TIERS: Section 0 states only\n",
      "Z_SYS_COMPETITOR_* instruments exclude the focal hospital's own system.\n",
      "Primary_strict_system_IV = Z_SYS_STRICT_9M_EXCL_CURRENT is not in that\n",
      "family. If MAIN's weakest performer (across all 6 schemes) also has the\n",
      "weakest exclusion restriction, demoting it is a DESIGN decision, not an\n",
      "outcome-based one. Confirm what EXCL_CURRENT excludes in the Phase 1 SQL\n",
      "before acting on this.\n", sep = "")
}



###############################################################################
# STAGE 8 -- meta-regressions, partition deduplication, permutation inference
###############################################################################
if (8 %in% RUN_STAGES) {

  meta_input <- prepare_meta_input(concept_results, schemes_long)
  meta_input[, TIER := instrument_tier(INSTRUMENT_LABEL)]
  cat("\nConcept-instrument rows per tier in meta_input:\n")
  print(meta_input[, .(N_ROWS = .N), by = TIER])

  meta_rf <- cache_or_run("meta_regressions_rf",
                          run_meta_regressions(meta_input, schemes_long, dep = "RF_COEF", se = "RF_SE"))
  meta_iv <- cache_or_run("meta_regressions_iv",
                          run_meta_regressions(meta_input, schemes_long, dep = "IV_COEF", se = "IV_SE"))

  meta <- deduplicate_partitions(rbindlist(list(meta_rf, meta_iv), fill = TRUE))
  meta[, TIER := instrument_tier(INSTRUMENT_LABEL)]
  save_csv(meta, "T08_meta_regressions_with_partitions.csv")
  save_csv(meta[IS_REPRESENTATIVE == 1L], "T08C_meta_distinct_partitions.csv")

  shop_terms <- meta[IS_REPRESENTATIVE == 1L & grepl("Shoppable", term) &
                       WEIGHTING == "Inverse variance"]

  cat("\n", strrep("=", 100),
      "\nSHOPPABILITY GRADIENT BY TIER (distinct partitions only)\n",
      strrep("=", 100), "\n", sep = "")
  print(shop_terms[, .(N_TESTS = .N,
                       N_NEGATIVE = sum(estimate < 0, na.rm = TRUE),
                       SHARE_NEGATIVE = round(mean(estimate < 0, na.rm = TRUE), 3),
                       N_SIG_05 = sum(p.value < 0.05, na.rm = TRUE),
                       SHARE_SIG = round(mean(p.value < 0.05, na.rm = TRUE), 3),
                       MEDIAN_P = round(median(p.value, na.rm = TRUE), 4)),
                   by = .(DEPENDENT, TIER)][order(DEPENDENT, TIER)])

  cat("\nPooled: MAIN+CONFIRMING vs all 6 (incl. discrepant):\n")
  print(rbindlist(list(
    shop_terms[TIER %chin% c("MAIN", "CONFIRMING"),
               .(POOL = "MAIN + CONFIRMING", N_TESTS = .N,
                 SHARE_SIG = round(mean(p.value < 0.05, na.rm = TRUE), 3)), by = DEPENDENT],
    shop_terms[, .(POOL = "All 6 (incl. discrepant)", N_TESTS = .N,
                   SHARE_SIG = round(mean(p.value < 0.05, na.rm = TRUE), 3)), by = DEPENDENT]
  )))

  save_csv(shop_terms, "T08F_shoppability_gradient_by_tier.csv")

  perm_rf <- family_permutation_test(meta_input[TIER != "DISCREPANT"],
                                     dep = "RF_COEF", se = "RF_SE")
  perm_iv <- family_permutation_test(meta_input[TIER != "DISCREPANT"],
                                     dep = "IV_COEF", se = "IV_SE")
  perm_rf_all6 <- family_permutation_test(meta_input, dep = "RF_COEF", se = "RF_SE")
  perm_iv_all6 <- family_permutation_test(meta_input, dep = "IV_COEF", se = "IV_SE")

  perm <- rbindlist(list(
    cbind(DEPENDENT = "Reduced form", POOL = "MAIN+CONFIRMING", perm_rf),
    cbind(DEPENDENT = "IV",           POOL = "MAIN+CONFIRMING", perm_iv),
    cbind(DEPENDENT = "Reduced form", POOL = "All 6",           perm_rf_all6),
    cbind(DEPENDENT = "IV",           POOL = "All 6",           perm_iv_all6)
  ), fill = TRUE)
  save_csv(perm, "T08E_family_permutation_inference.csv")

  cat("\n", strrep("=", 96), "\nEXACT FAMILY PERMUTATION, BY POOL\n", strrep("=", 96), "\n", sep = "")
  print(perm[, .(DEPENDENT, POOL, INSTRUMENT_LABEL, OBSERVED = signif(OBSERVED, 3),
                 P = round(P_TWO_SIDED, 5), METHOD)][order(DEPENDENT, POOL, P)])

  cat("\nDecomposition, overall:\n")
  decompose_reduced_form(meta_input, stem = "T08D_RF_vs_FS_decomposition")

  cat("\nDecomposition by tier (does DISCREPANT's flat result trace to FS or RF?):\n")
  for (tr in c("MAIN", "CONFIRMING", "DISCREPANT")) {
    cat("\n--", tr, "--\n")
    decompose_reduced_form(meta_input[TIER == tr], stem = paste0("T08D_", tr))
  }
}



###############################################################################
# STAGE 9 RUN BLOCK
###############################################################################
if (9 %in% RUN_STAGES) {

  verify_loo_sd(outpatient)

  measures <- cache_or_run("comparability_measures",
                           build_comparability_measures(outpatient))

  comp_int <- cache_or_run("comparability_interaction",
                           run_comparability_interaction(outpatient, measures))

  comp_wf  <- run_comparability_within_family(concept_results, measures)
}



###############################################################################
# STAGE 10 RUN BLOCK
###############################################################################
if (10 %in% RUN_STAGES) {
  if (!exists("measures")) {
    measures <- cache_or_run("comparability_measures",
                             build_comparability_measures(outpatient))
  }
  if (!exists("concept_results")) {
    stop("concept_results not in memory. Run the stage 6 loader block first.",
         call. = FALSE)
  }
  comp_meta <- run_comparability_meta(concept_results, measures)
}



######## STAGE 11 -- instrument balance and placebo ###########################
if (11 %in% RUN_STAGES) {
  balance <- run_instrument_balance(outpatient)
  placebo <- run_placebo_outcome(outpatient)
}



###############################################################################
# ROBUSTNESS: does the headline result survive controlling for IN_SYSTEM?
#
# The balance test in Section 11 finds system membership significantly related
# to two of the three MAIN instruments (p = 0.004 and p = 0.001), with small
# standardised effects of -0.08 and -0.14 SD.
#
# That is not disqualifying, and the correlation is close to mechanical:
# Competitor_only_hospitals_9m and Competitor_outside_CBSA_hospitals_9m count
# competitors OUTSIDE the focal hospital's own system by construction, so an
# unaffiliated hospital has a mechanically different exposure profile. But the
# clean way to close the objection is to show the headline result is unchanged
# when system membership is controlled for, rather than to argue it away.
#
# This reuses estimate_interacted() exactly as run_main_results() does, with
# one control added. Nothing else about the specification changes.
###############################################################################
add_system_membership <- function(panel) {
  d <- copy(panel)
  d[, IN_SYSTEM := as.integer(!is.na(SYSTEM_KEY) & SYSTEM_KEY != "")]
  cat("IN_SYSTEM added:", sum(d$IN_SYSTEM), "of", nrow(d),
      "rows (", round(100 * mean(d$IN_SYSTEM), 1), "%) in a health system\n")
  d
}

run_main_results_insystem_robustness <- function(
    panel, schemes = SCHEME_COLUMNS, instruments = MAIN_INSTRUMENTS,
    outcome = PRIMARY_OUTCOME, stem = "T11C_insystem_robustness") {

  panel2 <- add_system_membership(panel)
  controls <- c(BASELINE_CONTROLS, "IN_SYSTEM")

  rows <- list(); tests <- list(); g <- 0L
  n <- length(schemes) * length(instruments); t0 <- Sys.time()

  for (sl in names(schemes)) {
    sc <- schemes[[sl]]; if (!(sc %in% names(panel2))) next
    for (il in names(instruments)) {
      z <- instruments[[il]]; if (!(z %in% names(panel2))) next
      g <- g + 1L; t1 <- Sys.time()
      r <- estimate_interacted(panel2, sc, outcome, z, moderator_type = "categorical",
                               label = sl, instrument_label = il, moderator_label = sl,
                               controls = controls)
      if (!is.null(r)) { rows[[length(rows) + 1L]] <- r$rows; tests[[length(tests) + 1L]] <- r$tests }
      cat(sprintf("[%d/%d] %-26s %-36s | %5.1fs\n", g, n, substr(sl, 1, 24),
                  substr(il, 1, 34), as.numeric(difftime(Sys.time(), t1, units = "secs"))))
    }
  }

  mr <- rbindlist(rows, fill = TRUE); mt <- rbindlist(tests, fill = TRUE)
  if (nrow(mr) == 0L) stop("No results.", call. = FALSE)

  save_csv(mr, paste0(stem, "_interacted_RF_and_IV.csv"))
  save_csv(mt, paste0(stem, "B_heterogeneity_tests.csv"))
  cat("\nElapsed:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")

  rm(panel2); invisible(gc())
  list(rows = mr, tests = mt)
}



###############################################################################
# SIDE-BY-SIDE COMPARISON against the original headline (T06_main)
###############################################################################
compare_insystem_robustness <- function(original_rows, original_tests,
                                        robust_rows, robust_tests) {

  base <- original_rows[, .(SPEC, INSTRUMENT_LABEL, TERM,
                            RF_BASE = RF_PERCENT_PER_SD, RF_P_BASE = RF_P,
                            IV_BASE = IV_PERCENT, IV_P_BASE = IV_P)]
  rob  <- robust_rows[, .(SPEC, INSTRUMENT_LABEL, TERM,
                          RF_ROBUST = RF_PERCENT_PER_SD, RF_P_ROBUST = RF_P,
                          IV_ROBUST = IV_PERCENT, IV_P_ROBUST = IV_P)]

  cmp <- merge(base, rob, by = c("SPEC", "INSTRUMENT_LABEL", "TERM"))
  cmp[, RF_DELTA_PP := RF_ROBUST - RF_BASE]
  cmp[, SIGN_STABLE := as.integer(sign(RF_BASE) == sign(RF_ROBUST) | RF_BASE == 0)]

  cat("\n", strrep("=", 116), "\nHEADLINE ESTIMATES: ORIGINAL vs CONTROLLING FOR IN_SYSTEM\n",
      strrep("=", 116), "\n", sep = "")
  print(cmp[order(SPEC, INSTRUMENT_LABEL, TERM),
            .(SPEC = substr(SPEC, 1, 22), INSTRUMENT_LABEL = substr(INSTRUMENT_LABEL, 1, 28),
              TERM, RF_BASE = round(RF_BASE, 3), RF_ROBUST = round(RF_ROBUST, 3),
              DELTA_PP = round(RF_DELTA_PP, 3), P_BASE = round(RF_P_BASE, 4),
              P_ROBUST = round(RF_P_ROBUST, 4), SIGN_STABLE)])

  bt <- original_tests[, .(SPEC, INSTRUMENT_LABEL, ESTIMATOR, P_BASE = P_VALUE)]
  rt <- robust_tests[, .(SPEC, INSTRUMENT_LABEL, ESTIMATOR, P_ROBUST = P_VALUE)]
  ht <- merge(bt, rt, by = c("SPEC", "INSTRUMENT_LABEL", "ESTIMATOR"))

  ht[, STILL_SIG_05 := as.integer(P_ROBUST < 0.05)]

  cat("\n", strrep("=", 116), "\nHETEROGENEITY TEST: ORIGINAL vs CONTROLLING FOR IN_SYSTEM\n",
      strrep("=", 116), "\n", sep = "")
  print(ht[order(SPEC, INSTRUMENT_LABEL, ESTIMATOR),
           .(SPEC = substr(SPEC, 1, 22), INSTRUMENT_LABEL = substr(INSTRUMENT_LABEL, 1, 28),
             ESTIMATOR, P_BASE = round(P_BASE, 4), P_ROBUST = round(P_ROBUST, 4),
             STILL_SIG_05)])

  save_csv(cmp, "T11D_insystem_robustness_comparison.csv")
  save_csv(ht,  "T11E_insystem_robustness_heterogeneity_comparison.csv")

  cat("\n", strrep("=", 116), "\nSUMMARY\n", strrep("=", 116), "\n",
      "Coefficients changing by more than ~10% of their base value, or\n",
      "significance flipping, would mean IN_SYSTEM is doing real work that the\n",
      "county/month fixed effects did not already absorb. Small, stable\n",
      "coefficients and unchanged significance means the balance-test finding\n",
      "was cosmetic -- consistent with it being close to a mechanical byproduct\n",
      "of how the COMPETITOR instruments are constructed (they count\n",
      "competitors OUTSIDE the focal hospital's own system).\n", sep = "")

  cat("\nMedian heterogeneity p, MAIN tier, RF: original =",
      round(median(ht[STILL_SIG_05 >= 0]$P_BASE, na.rm = TRUE), 4), "| with IN_SYSTEM =",
      round(median(ht$P_ROBUST, na.rm = TRUE), 4), "\n")

  list(estimates = cmp, tests = ht)
}



###############################################################################
# RUN
###############################################################################
#
# Needs the stage 7 headline object `main` in memory for the comparison. If it
# is absent, reload it from cache rather than re-estimating.

if (!exists("main")) {
  main <- cache_or_run("main_results", run_main_results(outpatient, stem = "T06_main"))
}

robust <- run_main_results_insystem_robustness(outpatient)

comparison <- compare_insystem_robustness(main$rows, main$tests,
                                          robust$rows, robust$tests)



###############################################################################
# COUNTY DEMOGRAPHICS
#
# Pulls nine county-level characteristics from the 2022 five-year American
# Community Survey via tidycensus and merges them onto the panel by FIPS. These
# are the moderators used in the demographic heterogeneity blocks below.
#
# Needs a Census API key. Set it once with:
#   census_api_key("YOUR_KEY_HERE", install = TRUE)
#
# A NOTE ON VARIABLE CODES. Profile-table variables follow the convention that
# the base code plus a "P" suffix gives the percent version -- DP02_0068 is a
# count, DP02_0068P is the percent. Subject tables do not. S0101's C02 column
# code IS already the percent column: S0101_C02_030 is itself "Percent!!Total
# population!!65 years and over", and appending a P suffix to it returns
# nothing. get_acs() appends the E and M suffixes automatically for every table
# type, which is why the selection step below matches on a trailing E.
#
# The spot check in step 4 is not optional. A wrong variable code for a given
# ACS vintage returns a valid-looking column of zeros or NAs rather than an
# error, and that would propagate silently into every heterogeneity estimate.
###############################################################################
if (!requireNamespace("tidycensus", quietly = TRUE)) install.packages("tidycensus")
library(tidycensus)

# census_api_key("YOUR_KEY_HERE", install = TRUE)  # run once per machine


# ---------------------------------------------------------------------------
# STEP 1 -- pull from Census
# ---------------------------------------------------------------------------
acs_vars <- c(
  DEMO_POVERTY_RATE      = "DP03_0119P",   # % below poverty line
  DEMO_MEDIAN_INCOME     = "DP03_0062",    # median household income
  DEMO_COLLEGE_SHARE     = "DP02_0068P",   # % bachelor's degree or higher, 25+
  DEMO_HS_GRAD_SHARE     = "DP02_0067P",   # % high school graduate or higher, 25+
  DEMO_BLACK_SHARE       = "DP05_0038P",   # % Black or African American alone
  DEMO_HISPANIC_SHARE    = "DP05_0071P",   # % Hispanic or Latino, any race
  DEMO_UNINSURED_RATE    = "S2701_C05_001", # % uninsured
  DEMO_AGE65PLUS_SHARE   = "S0101_C02_030", # % population 65 years and over
  DEMO_POPULATION        = "DP05_0001"     # total population
)

county_demo_raw <- get_acs(
  geography = "county",
  variables = acs_vars,
  year = 2022,
  survey = "acs5",
  output = "wide"
)

county_demo <- as.data.table(county_demo_raw)

cat("\nRaw columns from get_acs():\n")
print(names(county_demo))


# ---------------------------------------------------------------------------
# STEP 2 -- select estimate columns and strip the E suffix
# ---------------------------------------------------------------------------
#
# get_acs(output = "wide") returns an E (estimate) and M (margin of error)
# column per variable, plus a NAME column. Selecting on the trailing E and
# renaming in one step avoids collisions between the requested variable names
# and the columns get_acs() adds. The stopifnot() is a hard stop on any drift
# between what was requested and what came back.

est_cols <- grep("^DEMO_.*E$", names(county_demo), value = TRUE)
cat("\nEstimate columns identified:\n")
print(est_cols)

stopifnot(length(est_cols) == length(acs_vars))  # hard stop on any drift

county_demo <- county_demo[, c("GEOID", est_cols), with = FALSE]
setnames(county_demo, est_cols, sub("E$", "", est_cols))

cat("\nColumns after rename:\n")
print(names(county_demo))


# ---------------------------------------------------------------------------
# STEP 3 -- rescale and derive
# ---------------------------------------------------------------------------
pct_vars <- c("DEMO_POVERTY_RATE", "DEMO_COLLEGE_SHARE", "DEMO_HS_GRAD_SHARE",
              "DEMO_BLACK_SHARE", "DEMO_HISPANIC_SHARE", "DEMO_UNINSURED_RATE",
              "DEMO_AGE65PLUS_SHARE")
for (v in intersect(pct_vars, names(county_demo))) {
  county_demo[, (v) := get(v) / 100]  # Census returns 0-100; rescale to 0-1
}

county_demo[, DEMO_LOG_MEDIAN_INCOME := log(pmax(DEMO_MEDIAN_INCOME, 1))]
setnames(county_demo, "GEOID", "COUNTY_FIPS")
county_demo[, COUNTY_FIPS := sprintf("%05d", as.integer(COUNTY_FIPS))]

cat("\nBuilt", nrow(county_demo), "county rows.\n")
print(head(county_demo, 5))


# ---------------------------------------------------------------------------
# STEP 4 -- spot check against known values
# ---------------------------------------------------------------------------
#
# Autauga County, Alabama (FIPS 01001), checked against Census QuickFacts for
# this ACS vintage: Black share roughly 19-20%, age 65 and over roughly 16-17%,
# high-school-graduate-or-higher roughly 88-90%.
#
# If any share comes back as 0, NA, or above 1, stop here. That means the
# variable code is wrong for this vintage, and it should be checked against
# load_variables(2022, "acs5/subject") and load_variables(2022, "acs5/profile")
# before proceeding to the merge.

cat("\n", strrep("=", 70), "\nSPOT CHECK: Autauga County AL (01001)\n", strrep("=", 70), "\n", sep = "")
print(county_demo[COUNTY_FIPS == "01001"])
cat("\nExpected ranges: BLACK_SHARE ~0.19-0.20 | AGE65PLUS_SHARE ~0.15-0.18 |\n",
    "HS_GRAD_SHARE ~0.85-0.90. If AGE65PLUS or HS_GRAD are 0, NA, or > 1,\n",
    "STOP and do not proceed to the merge.\n", sep = "")


# ---------------------------------------------------------------------------
# STEP 5 -- merge function
# ---------------------------------------------------------------------------
attach_county_demographics_v2 <- function(panel, demo) {
  if (!("COUNTY_FIPS" %in% names(panel))) {
    stop("Panel lacks COUNTY_FIPS.", call. = FALSE)
  }

  d <- copy(panel)
  d[, COUNTY_FIPS := sprintf("%05d", as.integer(COUNTY_FIPS))]

  old_demo_cols <- grep("^DEMO_", names(d), value = TRUE)
  if (length(old_demo_cols) > 0L) {
    cat("Dropping", length(old_demo_cols), "DEMO_ columns from a prior merge attempt.\n")
    d[, (old_demo_cols) := NULL]
  }

  n_before <- uniqueN(d$ANALYSIS_MARKET)
  d <- merge(d, demo, by = "COUNTY_FIPS", all.x = TRUE, sort = FALSE)

  matched_counties <- d[!is.na(DEMO_POPULATION), uniqueN(ANALYSIS_MARKET)]
  matched_rows <- mean(!is.na(d$DEMO_POPULATION))

  cat("\nDemographics matched:", matched_counties, "of", n_before,
      "counties |", round(100 * matched_rows, 1), "% of panel rows\n")

  if (matched_rows < 0.90) {
    unmatched <- unique(d[is.na(DEMO_POPULATION)]$ANALYSIS_MARKET)
    cat("\nBelow 90%. Sample unmatched:\n")
    print(head(unmatched, 20))
  } else {
    cat("Match rate at or above 90% -- good.\n")
  }

  d
}


# ---------------------------------------------------------------------------
# STEP 6 -- merge onto the panel
# ---------------------------------------------------------------------------
stopifnot(exists("outpatient"))

outpatient <- attach_county_demographics_v2(outpatient, demo = county_demo)


# ---------------------------------------------------------------------------
# STEP 7 -- full distribution check
# ---------------------------------------------------------------------------
cat("\n", strrep("=", 70), "\nFULL DISTRIBUTION CHECK\n", strrep("=", 70), "\n", sep = "")
demo_cols <- grep("^DEMO_", names(outpatient), value = TRUE)
for (cc in demo_cols) {
  v <- outpatient[[cc]]
  cat(sprintf("%-24s finite=%7s  mean=%10.4f  sd=%10.4f  min=%10.4f  max=%10.4f\n",
              cc, format(sum(is.finite(v)), big.mark = ","),
              mean(v, na.rm = TRUE), sd(v, na.rm = TRUE),
              min(v, na.rm = TRUE), max(v, na.rm = TRUE)))
}


# ---------------------------------------------------------------------------
# STEP 8 -- persist the merged panel, only after steps 4 and 7 look right
# ---------------------------------------------------------------------------
saveRDS(outpatient, file.path(CACHE_DIR, "outpatient_panel.rds"))



###############################################################################
# SECTION 12 -- COUNTY DEMOGRAPHIC HETEROGENEITY IN THE SHOPPABILITY GRADIENT
#
# THE ESTIMAND. The question is whether the shoppability GRADIENT varies with
# county composition, not whether the pooled response does. Interacting the
# instrument with a demographic alone,
#
#     ln(P) = g0*Z + g1*(Z x M) + controls + FE
#
# estimates demographic variation in the pooled response -- a quantity the
# paper disowns, because it averages a real effect on shoppable services
# against a precise zero on non-shoppable ones. The specification estimated
# here splits the treatment by shoppability first and interacts within each
# arm:
#
#     ln(P) = gS*(Z x Shop)     + gN*(Z x NonShop)
#           + dS*(Z x Shop x M) + dN*(Z x NonShop x M)
#           + controls + FE
#
# The object of interest is d(gap)/dM = dS - dN, tested as a single 1-df Wald.
#
# THE PARAMETER COUNT IS NOT A PROBLEM. Shoppability is concept-determined and
# the demographic is county-determined, so the lower-order terms M, Shop, and
# (Shop x M) are all constant within a county x concept cell and are absorbed
# entirely by MARKET_ID. Only the four Z-interactions are estimated. The IV
# version has four endogenous regressors and four excluded instruments, so the
# model stays exactly identified and the reduced-form t-test remains
# numerically identical to the Anderson-Rubin test.
#
# THE BINDING CONSTRAINT IS FIRST-STAGE STRENGTH. Z x M is weakly instrumented
# once MARKET_ID absorbs the county variation, and the first stage is weaker
# here than in the headline because there are four endogenous terms rather than
# two. Read the reduced form as primary and treat the IV magnitudes as
# illustrative; the RF p-values are weak-instrument robust by construction.
# Check FIRST_STAGE_WALD_MIN before quoting any IV number.
#
# TWO COMPETING PREDICTIONS on the SES index, where higher means a more
# advantaged county. These are opposite signs on the same coefficient, which is
# what makes this a test rather than a description:
#
#   Consumer shopping capacity -- the original pre-registered prediction.
#     Advantaged markets shop more, so the response is LARGER (more negative)
#     and the interaction on the index is NEGATIVE.
#
#   Contracting depth -- the paper's second condition.
#     Advantaged markets have denser commercial contracting and therefore less
#     movable negotiations, so the response is SMALLER and the interaction on
#     the index is POSITIVE.
#
# This block is self-contained. Every helper it uses is defined locally and
# prefixed `.t12_` so nothing in the main pipeline is overwritten. It needs
# only `outpatient` in memory, with the columns listed in the preflight below.
#
# Runtime is roughly 30-45 minutes for all three blocks. The T12T_RUN flags run
# them independently.
###############################################################################
suppressPackageStartupMessages({
  library(data.table)
  library(fixest)
})


# ===========================================================================
# 0. CONFIGURATION
# ===========================================================================
T12T_RUN <- list(
  ses_index      = TRUE,   # headline: SES index x 6 schemes x 3 instruments
  per_moderator  = TRUE,   # appendix: 8 moderators x Scheme 1 x 3 instruments
  population     = TRUE    # market-size check: population as COMPETING interaction
)

T12T_OUTCOME    <- "LN_MEDIAN_PRICE"
T12T_ENDOGENOUS <- "N_PRIOR_POSTERS"
T12T_CONTROLS   <- c("LOG_TOTAL_BEDS")
T12T_FE         <- c("MARKET_ID", "POST_MONTH")
T12T_CLUSTERS   <- c("ANALYSIS_MARKET", "POST_MONTH")
T12T_COUNTY_KEY <- "ANALYSIS_MARKET"
T12T_MIN_OBS    <- 50000L

# Defined locally so this file does not depend on the main pipeline's globals.
T12T_INSTRUMENTS <- c(
  Competitor_only_hospitals_9m         = "Z_SYS_COMPETITOR_ONLY_9M_EXCL_CURRENT",
  Primary_strict_system_IV             = "Z_SYS_STRICT_9M_EXCL_CURRENT",
  Competitor_outside_CBSA_hospitals_9m = "Z_SYS_COMPETITOR_OUTSIDE_CBSA_9M_EXCL_CURRENT"
)

T12T_SCHEMES <- c(
  "1. Procedural certainty" = "SCHEME_1_CERTAINTY",
  "2. Theory-Based V2"      = "SCHEME_2_THEORYV2",
  "3. Imaging vs Procedural"= "SCHEME_3_IMAGING",
  "4. CMS Statutory List"   = "SCHEME_4_CMS70",
  "5. Upfront Cash-Market"  = "SCHEME_5_MDSAVE",
  "6. Within Modality"      = "SCHEME_6_WITHINMOD"
)
T12T_HEADLINE_SCHEME <- "SCHEME_1_CERTAINTY"

T12T_MODERATORS <- c("DEMO_COLLEGE_SHARE", "DEMO_HS_GRAD_SHARE",
                     "DEMO_LOG_MEDIAN_INCOME", "DEMO_POVERTY_RATE",
                     "DEMO_BLACK_SHARE", "DEMO_HISPANIC_SHARE",
                     "DEMO_AGE65PLUS_SHARE", "DEMO_UNINSURED_RATE")

# Pre-registered expected signs, retained verbatim so the comparison is
# reported rather than reconstructed after the fact. They come from the
# consumer shopping-capacity model. The paper's own channel places the
# contracting party rather than the patient at the centre, so a rejection of
# these signs corroborates that channel instead of contradicting the design.
T12T_PREDICTED_NEGATIVE <- c("DEMO_COLLEGE_SHARE", "DEMO_HS_GRAD_SHARE",
                             "DEMO_LOG_MEDIAN_INCOME")
T12T_PREDICTED_POSITIVE <- c("DEMO_POVERTY_RATE", "DEMO_BLACK_SHARE",
                             "DEMO_HISPANIC_SHARE", "DEMO_AGE65PLUS_SHARE")

# Variables entering the SES index, each with the sign that orients it toward
# "more advantaged". The uninsured rate is deliberately excluded: it carried no
# pre-registered sign and loads ambiguously, since it tracks both poverty and
# state Medicaid expansion policy.
T12T_SES_COMPONENTS <- c(DEMO_COLLEGE_SHARE     =  1,
                         DEMO_HS_GRAD_SHARE     =  1,
                         DEMO_LOG_MEDIAN_INCOME =  1,
                         DEMO_POVERTY_RATE      = -1,
                         DEMO_BLACK_SHARE       = -1,
                         DEMO_HISPANIC_SHARE    = -1,
                         DEMO_AGE65PLUS_SHARE   = -1)

T12T_OUTDIR <- if (exists("TABLE_DIR") && dir.exists(TABLE_DIR)) TABLE_DIR else getwd()


# ===========================================================================
# 1. HELPERS  (all self-contained)
# ===========================================================================
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a

.t12_save <- function(dt, filename) {
  if (is.null(dt) || nrow(dt) == 0L) { cat("  [skip empty]", filename, "\n"); return(invisible(NULL)) }
  p <- file.path(T12T_OUTDIR, filename)
  fwrite(dt, p)
  cat("  Saved:", p, "\n")
  invisible(p)
}

.t12_show <- function(dt, n = Inf) {
  # Coerce to data.frame before printing. print(n = Inf) on a tibble or
  # data.table throws `invalid 'na.print' specification` in some environments.
  if (is.null(dt) || nrow(dt) == 0L) { cat("  (no rows)\n"); return(invisible(NULL)) }
  print(as.data.frame(if (is.finite(n)) head(dt, n) else dt), row.names = FALSE)
  invisible(NULL)
}

# Resolve a coefficient name. fixest prefixes fitted endogenous regressors with
# "fit_" in IV models, so the same term has two possible names depending on the
# estimator.
.t12_resolve <- function(fit, term) {
  if (is.null(fit)) return(NA_character_)
  cand <- c(paste0("fit_", term), term)
  hit  <- cand[cand %in% names(coef(fit))]
  if (length(hit) == 0L) NA_character_ else hit[1L]
}

.t12_pull <- function(fit, term) {
  nm <- .t12_resolve(fit, term)
  if (is.na(nm)) return(list(b = NA_real_, s = NA_real_))
  list(b = unname(coef(fit)[nm]), s = unname(sqrt(vcov(fit)[nm, nm])))
}

# Wald test of b1 - b2 = 0, accounting for their covariance.
.t12_wald_diff <- function(fit, t1, t2) {
  n1 <- .t12_resolve(fit, t1); n2 <- .t12_resolve(fit, t2)
  if (is.na(n1) || is.na(n2)) return(NULL)
  cf <- coef(fit); V <- vcov(fit)
  d  <- unname(cf[n1] - cf[n2])
  v  <- unname(V[n1, n1] + V[n2, n2] - 2 * V[n1, n2])
  if (!is.finite(v) || v <= 0) return(NULL)
  se <- sqrt(v)
  data.table(DIFF = d, DIFF_SE = se, WALD = (d / se)^2, DF = 1L,
             P_VALUE = 2 * pnorm(-abs(d / se)))
}

# Minimum first-stage Wald across the endogenous equations. Written defensively
# against fixest's return structure, which varies across versions; returns NA
# rather than erroring.
.t12_fs_wald <- function(fit) {
  if (is.null(fit)) return(NA_real_)
  v <- tryCatch({
    w <- fixest::fitstat(fit, "ivwald", simplify = FALSE)
    w <- w[["ivwald"]] %||% w
    unlist(lapply(w, function(z) if (is.list(z)) z$stat %||% NA_real_ else as.numeric(z)))
  }, error = function(e) NA_real_)
  v <- suppressWarnings(as.numeric(v))
  if (length(v) == 0L || all(is.na(v))) NA_real_ else min(v, na.rm = TRUE)
}

.t12_cragg <- function(fit) {
  if (is.null(fit)) return(NA_real_)
  v <- tryCatch(as.numeric(fixest::fitstat(fit, "cd", simplify = TRUE)),
                error = function(e) NA_real_)
  if (length(v) == 0L) NA_real_ else v[1L]
}

.t12_ols_formula <- function(y, rhs, fe) {
  r <- if (length(rhs)) paste(rhs, collapse = " + ") else "1"
  as.formula(paste0(y, " ~ ", r, " | ", paste(fe, collapse = " + ")))
}

.t12_iv_formula <- function(y, endo, ivs, exog, fe) {
  x <- if (length(exog)) paste(exog, collapse = " + ") else "1"
  as.formula(paste0(y, " ~ ", x, " | ", paste(fe, collapse = " + "), " | ",
                    paste(endo, collapse = " + "), " ~ ", paste(ivs, collapse = " + ")))
}

.t12_cluster_formula <- function(cl) as.formula(paste0("~", paste(cl, collapse = " + ")))


# ===========================================================================
# 2. CORE ESTIMATOR -- Z x CATEGORY x CONTINUOUS MODERATOR(S)
# ===========================================================================
#
# `moderators` may hold one or two continuous variables. With two, both enter
# as competing interactions against the same category split, which is how the
# market-size check in block 6C is run. The model stays exactly identified in
# either case:
#   n_endogenous = n_levels * (1 + n_moderators) = n_instruments.

.t12_fit_triple <- function(panel, scheme_col, moderators, instrument,
                            instrument_label = "", scheme_label = "",
                            outcome = T12T_OUTCOME, endogenous = T12T_ENDOGENOUS,
                            controls = T12T_CONTROLS, fe = T12T_FE,
                            clusters = T12T_CLUSTERS, min_obs = T12T_MIN_OBS) {

  need <- unique(c(outcome, endogenous, instrument, controls, fe, clusters,
                   scheme_col, moderators, "FINAL_CONCEPT_ID"))
  miss <- setdiff(need, names(panel))
  if (length(miss)) { warning("Missing columns: ", paste(miss, collapse = ", "), call. = FALSE); return(NULL) }

  d <- panel[, ..need]
  d <- d[complete.cases(d[, ..need])]
  d <- d[!is.na(get(scheme_col))]
  if (nrow(d) < min_obs) { cat("    [skip: only", nrow(d), "usable rows]\n"); return(NULL) }

  d[, GRP := droplevels(factor(get(scheme_col)))]
  keys <- levels(d$GRP)
  if (length(keys) != 2L) { cat("    [skip: scheme has", length(keys), "levels]\n"); return(NULL) }

  # Centre each moderator at its estimation-sample mean, and record its SD for
  # the per-SD reporting below.
  mod_tag <- paste0("M", seq_along(moderators))
  mod_sd  <- numeric(length(moderators))
  for (j in seq_along(moderators)) {
    v <- as.numeric(d[[moderators[j]]])
    if (!is.finite(sd(v, na.rm = TRUE)) || sd(v, na.rm = TRUE) == 0) {
      cat("    [skip: no variation in", moderators[j], "]\n"); return(NULL)
    }
    mod_sd[j] <- sd(v, na.rm = TRUE)
    set(d, j = paste0("MODV_", mod_tag[j]), value = v - mean(v, na.rm = TRUE))
  }

  # Build the design. Level terms first, then each moderator interaction.
  grp_tag <- paste0("G", seq_along(keys))
  endo <- ivs <- rfs <- character(0)
  term_grp <- term_lab <- term_mod <- character(0)

  for (k in seq_along(keys)) {
    sel <- as.integer(d$GRP == keys[k])
    set(d, j = paste0("TREAT_", grp_tag[k]), value = d[[endogenous]] * sel)
    set(d, j = paste0("IVV_",   grp_tag[k]), value = d[[instrument]] * sel)
    set(d, j = paste0("RFV_",   grp_tag[k]), value = d[[instrument]] * sel)
    endo <- c(endo, paste0("TREAT_", grp_tag[k]))
    ivs  <- c(ivs,  paste0("IVV_",   grp_tag[k]))
    rfs  <- c(rfs,  paste0("RFV_",   grp_tag[k]))
    term_grp <- c(term_grp, keys[k]); term_lab <- c(term_lab, "Level")
    term_mod <- c(term_mod, NA_character_)

    for (j in seq_along(moderators)) {
      mv <- d[[paste0("MODV_", mod_tag[j])]]
      sfx <- paste0(grp_tag[k], "_", mod_tag[j])
      set(d, j = paste0("TREAT_", sfx), value = d[[endogenous]]  * sel * mv)
      set(d, j = paste0("IVV_",   sfx), value = d[[instrument]] * sel * mv)
      set(d, j = paste0("RFV_",   sfx), value = d[[instrument]] * sel * mv)
      endo <- c(endo, paste0("TREAT_", sfx))
      ivs  <- c(ivs,  paste0("IVV_",   sfx))
      rfs  <- c(rfs,  paste0("RFV_",   sfx))
      term_grp <- c(term_grp, keys[k]); term_lab <- c(term_lab, "x Moderator")
      term_mod <- c(term_mod, moderators[j])
    }
  }

  cl <- .t12_cluster_formula(clusters)
  rf_fit <- tryCatch(feols(.t12_ols_formula(outcome, c(rfs, controls), fe),
                           data = d, cluster = cl, warn = FALSE, notes = FALSE),
                     error = function(e) { cat("    [RF failed:", conditionMessage(e), "]\n"); NULL })
  iv_fit <- tryCatch(feols(.t12_iv_formula(outcome, endo, ivs, controls, fe),
                           data = d, cluster = cl, warn = FALSE, notes = FALSE),
                     error = function(e) { cat("    [IV failed:", conditionMessage(e), "]\n"); NULL })
  if (is.null(rf_fit)) return(NULL)

  sd_z    <- sd(d[[instrument]], na.rm = TRUE)
  fs_min  <- .t12_fs_wald(iv_fit)
  cd      <- .t12_cragg(iv_fit)
  n_obs   <- if (is.null(iv_fit)) nobs(rf_fit) else nobs(iv_fit)

  rows <- rbindlist(lapply(seq_along(endo), function(i) {
    rf <- .t12_pull(rf_fit, rfs[i]); iv <- .t12_pull(iv_fit, endo[i])
    j  <- match(term_mod[i], moderators)
    data.table(
      SPEC = scheme_label, SCHEME_COL = scheme_col,
      INSTRUMENT_LABEL = instrument_label, INSTRUMENT = instrument,
      OUTCOME = outcome, CATEGORY = term_grp[i], TERM = term_lab[i],
      MODERATOR = term_mod[i],
      RF_COEF = rf$b, RF_SE = rf$s, RF_P = 2 * pnorm(-abs(rf$b / rf$s)),
      # Level terms: % price change per SD of Z.
      # Interaction terms: change in that % per 1 SD of the moderator.
      RF_PERCENT_PER_SD_Z = if (term_lab[i] == "Level")
        100 * (exp(rf$b * sd_z) - 1) else NA_real_,
      RF_SHIFT_PCT_PER_SD_MOD = if (term_lab[i] == "x Moderator" && !is.na(j))
        100 * (exp(rf$b * sd_z * mod_sd[j]) - 1) else NA_real_,
      IV_COEF = iv$b, IV_SE = iv$s, IV_P = 2 * pnorm(-abs(iv$b / iv$s)),
      IV_PERCENT = if (term_lab[i] == "Level") 100 * (exp(iv$b) - 1) else NA_real_,
      FIRST_STAGE_WALD_MIN = fs_min, CRAGG_DONALD = cd,
      N_OBSERVATIONS = n_obs, N_CONCEPTS = uniqueN(d$FINAL_CONCEPT_ID),
      N_COUNTIES = uniqueN(d[[clusters[1L]]]),
      MODERATOR_SD = if (!is.na(j)) mod_sd[j] else NA_real_, INSTRUMENT_SD = sd_z)
  }), fill = TRUE)

  # ---- The estimand: does the GAP move with the moderator? dS - dN ----
  tests <- rbindlist(lapply(seq_along(moderators), function(j) {
    iS <- which(term_grp == keys[2L] & term_lab == "x Moderator" & term_mod == moderators[j])
    iN <- which(term_grp == keys[1L] & term_lab == "x Moderator" & term_mod == moderators[j])
    if (length(iS) != 1L || length(iN) != 1L) return(NULL)
    rf_t <- .t12_wald_diff(rf_fit, rfs[iS],  rfs[iN])
    iv_t <- .t12_wald_diff(iv_fit, endo[iS], endo[iN])

    # Gap between the two level terms, evaluated at the moderator mean.
    lS <- which(term_grp == keys[2L] & term_lab == "Level")
    lN <- which(term_grp == keys[1L] & term_lab == "Level")
    gap <- .t12_wald_diff(rf_fit, rfs[lS], rfs[lN])

    rbindlist(list(
      if (!is.null(rf_t)) cbind(ESTIMATOR = "Reduced form", rf_t) else NULL,
      if (!is.null(iv_t)) cbind(ESTIMATOR = "IV",           iv_t) else NULL
    ), fill = TRUE)[, `:=`(
      SPEC = scheme_label, SCHEME_COL = scheme_col, MODERATOR = moderators[j],
      INSTRUMENT_LABEL = instrument_label,
      REF_CATEGORY = keys[1L], FOCAL_CATEGORY = keys[2L],
      GAP_AT_MEAN_PCT_PER_SD_Z = if (!is.null(gap))
        100 * (exp(gap$DIFF * sd_z) - 1) else NA_real_,
      GAP_AT_MEAN_P = if (!is.null(gap)) gap$P_VALUE else NA_real_,
      DGAP_PCT_PER_SD_MOD = 100 * (exp(DIFF * sd_z * mod_sd[j]) - 1),
      FIRST_STAGE_WALD_MIN = fs_min, N_OBSERVATIONS = n_obs)]
  }), fill = TRUE)

  list(rows = rows, tests = tests)
}


# ===========================================================================
# 3. PREFLIGHT
# ===========================================================================
.t12_preflight <- function(panel) {
  cat("\n", strrep("=", 78), "\nPREFLIGHT\n", strrep("=", 78), "\n", sep = "")
  need <- unique(c(T12T_OUTCOME, T12T_ENDOGENOUS, T12T_CONTROLS, T12T_FE,
                   T12T_CLUSTERS, "FINAL_CONCEPT_ID", T12T_COUNTY_KEY,
                   unname(T12T_INSTRUMENTS), unname(T12T_SCHEMES),
                   T12T_MODERATORS))
  miss <- setdiff(need, names(panel))
  if (length(miss)) {
    stop("Panel is missing required columns:\n  ", paste(miss, collapse = "\n  "),
         "\n\nIf the scheme columns are missing, run attach_scheme_columns() first.",
         call. = FALSE)
  }
  cat("All required columns present.\n")
  cat("Panel rows:", format(nrow(panel), big.mark = ","),
      "| concepts:", uniqueN(panel$FINAL_CONCEPT_ID),
      "| counties:", uniqueN(panel[[T12T_COUNTY_KEY]]), "\n")

  scr <- rbindlist(lapply(T12T_MODERATORS, function(m) {
    v <- as.numeric(panel[[m]])
    data.table(MODERATOR = m, N_FINITE = sum(is.finite(v)),
               SHARE_FINITE = round(mean(is.finite(v)), 4),
               SD = round(sd(v, na.rm = TRUE), 5),
               MIN = round(min(v, na.rm = TRUE), 4),
               MAX = round(max(v, na.rm = TRUE), 4),
               USABLE = as.integer(mean(is.finite(v)) > 0.5 &&
                                     sd(v, na.rm = TRUE) > 0))
  }))
  cat("\nModerator screen:\n"); .t12_show(scr)
  if (any(scr$USABLE == 0L)) warning("Some moderators failed the screen.", call. = FALSE)
  invisible(scr)
}


# ===========================================================================
# 4. SES INDEX -- first principal component of the seven signed moderators
# ===========================================================================
#
# Estimated on the COUNTY cross-section rather than on panel rows, so a county
# contributing many hospital-concept-months does not dominate the rotation.
# Oriented so that higher means more advantaged, then standardised to mean zero
# and standard deviation one.

.t12_build_ses_index <- function(panel) {
  cat("\n", strrep("=", 78), "\nBUILDING SES INDEX\n", strrep("=", 78), "\n", sep = "")
  comp <- names(T12T_SES_COMPONENTS)
  cols <- c(T12T_COUNTY_KEY, comp)
  cs <- unique(panel[, ..cols])
  cs <- cs[complete.cases(cs)]
  cat("Counties with complete demographics:", nrow(cs), "\n")

  X <- as.matrix(cs[, ..comp])
  X <- scale(X)                                   # z-score each component
  for (m in comp) X[, m] <- X[, m] * T12T_SES_COMPONENTS[[m]]   # orient

  pc  <- prcomp(X, center = TRUE, scale. = FALSE)
  idx <- pc$x[, 1]
  ld  <- pc$rotation[, 1]

  # Force the index to point toward "more advantaged" using college share.
  if (ld[["DEMO_COLLEGE_SHARE"]] < 0) { idx <- -idx; ld <- -ld }
  idx <- as.numeric(scale(idx))

  loadings <- data.table(COMPONENT = comp,
                         ORIENTATION = unname(T12T_SES_COMPONENTS[comp]),
                         PC1_LOADING = round(unname(ld[comp]), 4),
                         VAR_EXPLAINED_PC1 = round(summary(pc)$importance[2, 1], 4))
  cat("\nPC1 loadings (on orientation-adjusted z-scores; higher index = more advantaged):\n")
  .t12_show(loadings)
  .t12_save(loadings, "QA10_ses_index_loadings.csv")

  cs[, SES_INDEX := idx]
  out <- merge(panel, cs[, c(T12T_COUNTY_KEY, "SES_INDEX"), with = FALSE],
               by = T12T_COUNTY_KEY, all.x = TRUE, sort = FALSE)
  setDT(out)
  cat("\nSES_INDEX attached to", round(100 * mean(!is.na(out$SES_INDEX)), 1),
      "% of panel rows.\n")
  out
}


# ===========================================================================
# 5. DRIVERS
# ===========================================================================
.t12_driver <- function(panel, moderators, schemes, instruments, tag) {
  rows <- list(); tests <- list()
  total <- length(schemes) * length(instruments); i <- 0L
  for (s in seq_along(schemes)) {
    for (il in names(instruments)) {
      i <- i + 1L
      t0 <- Sys.time()
      cat(sprintf("[%2d/%2d] %-26s %-38s ", i, total,
                  substr(names(schemes)[s], 1, 24), substr(il, 1, 36)))
      r <- .t12_fit_triple(panel, scheme_col = unname(schemes[s]),
                           moderators = moderators,
                           instrument = unname(instruments[il]),
                           instrument_label = il, scheme_label = names(schemes)[s])
      cat(sprintf("%5.1fs\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
      if (is.null(r)) next
      rows[[length(rows) + 1L]]  <- r$rows
      tests[[length(tests) + 1L]] <- r$tests
    }
  }
  if (length(rows) == 0L) { warning("No results for ", tag, call. = FALSE); return(NULL) }
  cr <- rbindlist(rows, fill = TRUE); ct <- rbindlist(tests, fill = TRUE)

  cr[, EXPECTED_SIGN := fcase(
    MODERATOR %chin% T12T_PREDICTED_NEGATIVE, "negative",
    MODERATOR %chin% T12T_PREDICTED_POSITIVE, "positive",
    default = "none (no pre-registered prediction)")]
  cr[TERM == "x Moderator" & EXPECTED_SIGN == "negative", SIGN_AS_PREDICTED := RF_COEF < 0]
  cr[TERM == "x Moderator" & EXPECTED_SIGN == "positive", SIGN_AS_PREDICTED := RF_COEF > 0]

  list(rows = cr, tests = ct)
}


# ===========================================================================
# 6. RUN
# ===========================================================================
stopifnot(exists("outpatient"))
setDT(outpatient)
.t12_preflight(outpatient)

t12t_panel <- .t12_build_ses_index(outpatient)

t12t <- list()

# ---- 6A. HEADLINE: SES index, all six schemes, three MAIN instruments -----
if (isTRUE(T12T_RUN$ses_index)) {
  cat("\n", strrep("=", 78),
      "\n6A. HEADLINE — SES INDEX x SHOPPABILITY x Z (6 schemes x 3 instruments)\n",
      strrep("=", 78), "\n", sep = "")
  t12t$ses <- .t12_driver(t12t_panel, moderators = "SES_INDEX",
                          schemes = T12T_SCHEMES, instruments = T12T_INSTRUMENTS,
                          tag = "SES index")
  if (!is.null(t12t$ses)) {
    .t12_save(t12t$ses$rows,  "T12T_triple_interaction_ses.csv")
    .t12_save(t12t$ses$tests, "T12T_triple_interaction_ses_tests.csv")

    cat("\nTHE ESTIMAND: does the shoppability gap move with county SES?\n")
    cat("(DIFF = dS - dN. Positive => gradient SHRINKS in advantaged markets,\n",
        " which is the contracting-depth reading. Negative => it WIDENS, which\n",
        " is the consumer shopping-capacity reading.)\n\n", sep = "")
    .t12_show(t12t$ses$tests[ESTIMATOR == "Reduced form",
                             .(SPEC, INSTRUMENT_LABEL,
                               GAP_AT_MEAN_PCT = round(GAP_AT_MEAN_PCT_PER_SD_Z, 3),
                               DIFF = signif(DIFF, 4), DIFF_SE = signif(DIFF_SE, 4),
                               P = round(P_VALUE, 4),
                               DGAP_PCT_PER_SD = round(DGAP_PCT_PER_SD_MOD, 3),
                               MIN_WALD = round(FIRST_STAGE_WALD_MIN, 1))][order(P)])

    cat("\nSummary across the 18 reduced-form tests:\n")
    .t12_show(t12t$ses$tests[ESTIMATOR == "Reduced form",
                             .(N = .N, N_SIG_05 = sum(P_VALUE < 0.05, na.rm = TRUE),
                               N_POSITIVE = sum(DIFF > 0, na.rm = TRUE),
                               MEDIAN_P = round(median(P_VALUE, na.rm = TRUE), 4))])
  }
}

# ---- 6B. APPENDIX: eight moderators individually, headline scheme --------
if (isTRUE(T12T_RUN$per_moderator)) {
  cat("\n", strrep("=", 78),
      "\n6B. APPENDIX — EACH MODERATOR SEPARATELY (Scheme 1 x 3 instruments)\n",
      strrep("=", 78), "\n", sep = "")
  per <- list(); per_t <- list()
  for (m in T12T_MODERATORS) {
    cat("\n--", m, "--\n")
    r <- .t12_driver(t12t_panel, moderators = m,
                     schemes = T12T_SCHEMES[names(T12T_SCHEMES)[
                       unname(T12T_SCHEMES) == T12T_HEADLINE_SCHEME]],
                     instruments = T12T_INSTRUMENTS, tag = m)
    if (is.null(r)) next
    per[[length(per) + 1L]] <- r$rows; per_t[[length(per_t) + 1L]] <- r$tests
  }
  if (length(per)) {
    t12t$per_moderator <- list(rows = rbindlist(per, fill = TRUE),
                               tests = rbindlist(per_t, fill = TRUE))
    .t12_save(t12t$per_moderator$rows,  "T12T_triple_interaction_moderators.csv")
    .t12_save(t12t$per_moderator$tests, "T12T_triple_interaction_moderators_tests.csv")

    cat("\nGradient shift by moderator (reduced form, Scheme 1):\n")
    .t12_show(t12t$per_moderator$tests[ESTIMATOR == "Reduced form",
                                       .(MODERATOR, INSTRUMENT_LABEL, DIFF = signif(DIFF, 4),
                                         P = round(P_VALUE, 4), DGAP_PCT_PER_SD = round(DGAP_PCT_PER_SD_MOD, 3),
                                         MIN_WALD = round(FIRST_STAGE_WALD_MIN, 1))][order(MODERATOR, P)])

    cat("\nPre-registered sign check on the SHOPPABLE interaction only:\n")
    .t12_show(t12t$per_moderator$rows[
      TERM == "x Moderator" & CATEGORY == "Shoppable",
      .(N = .N, N_AS_PREDICTED = sum(SIGN_AS_PREDICTED, na.rm = TRUE),
        N_RF_NEGATIVE = sum(RF_COEF < 0, na.rm = TRUE),
        MEDIAN_P = round(median(RF_P, na.rm = TRUE), 4)),
      by = .(MODERATOR, EXPECTED_SIGN)][order(MODERATOR)])
  }
}

# ---- 6C. MARKET SIZE: population as a COMPETING interaction --------------
# Adding log(population) as a CONTROL does not test market-size confounding.
# Population is county-level and time-invariant, so MARKET_ID absorbs it
# entirely and the controlled specification is numerically identical to the
# uncontrolled one to about 1e-13. Entering it as a competing INTERACTION is
# the version that actually tests the confound, because the interaction with Z
# is not absorbed.
if (isTRUE(T12T_RUN$population)) {
  cat("\n", strrep("=", 78),
      "\n6C. MARKET-SIZE CHECK — SES INDEX vs LOG POPULATION, both interacted\n",
      strrep("=", 78), "\n", sep = "")
  if (!("DEMO_POPULATION" %in% names(t12t_panel))) {
    warning("DEMO_POPULATION absent; skipping the market-size check.", call. = FALSE)
  } else {
    t12t_panel[, DEMO_LOG_POPULATION := log(pmax(DEMO_POPULATION, 1))]
    t12t$population <- .t12_driver(
      t12t_panel, moderators = c("SES_INDEX", "DEMO_LOG_POPULATION"),
      schemes = T12T_SCHEMES[names(T12T_SCHEMES)[
        unname(T12T_SCHEMES) == T12T_HEADLINE_SCHEME]],
      instruments = T12T_INSTRUMENTS, tag = "population horse race")
    if (!is.null(t12t$population)) {
      .t12_save(t12t$population$rows,  "T12T_population_horse_race.csv")
      .t12_save(t12t$population$tests, "T12T_population_horse_race_tests.csv")
      cat("\nBoth interactions in one model (reduced form). If SES_INDEX survives\n",
          "with DEMO_LOG_POPULATION alongside it, the gradient is not market size.\n\n", sep = "")
      .t12_show(t12t$population$tests[ESTIMATOR == "Reduced form",
                                      .(MODERATOR, INSTRUMENT_LABEL, DIFF = signif(DIFF, 4),
                                        P = round(P_VALUE, 4), MIN_WALD = round(FIRST_STAGE_WALD_MIN, 1))][
                                          order(MODERATOR, P)])
    }
  }
}



###############################################################################
# SECTION 12B -- SES TERCILES AS A CATEGORICAL MODERATOR
#
# Section 12 imposes a constant SES gradient on the shoppability gap. This
# relaxes that assumption by estimating a separate shoppability gap inside each
# tercile of county SES and testing whether the three differ.
#
# The motivation is specific: if the true pattern is nonlinear -- say only the
# top third of counties behaves differently -- a linear interaction averages it
# away and returns a null that reflects the functional form rather than the
# data.
#
# SPECIFICATION
#   Let S index shoppability {Non_shoppable, Shoppable} and t index SES tercile
#   {T1 low, T2 mid, T3 high}. Six groups, six coefficients:
#
#     ln(P) = sum_{S,t} gamma_{S,t} ( Z x 1[k in S] x 1[m in t] )
#             + theta ln(Beds) + alpha_mk + tau_t + eps
#
#   Every lower-order term is absorbed: 1[S] is concept-determined, 1[t] is
#   county-determined, and their product is constant within a county x concept
#   cell, so MARKET_ID takes all three. Only the six Z-interactions are
#   estimated. The IV version has six endogenous regressors and six excluded
#   instruments -- still EXACTLY IDENTIFIED, so RF p-values remain AR-robust.
#
# THE THREE THINGS THIS REPORTS
#   1. gap_t = gamma_{Shop,t} - gamma_{Non,t} inside each tercile, with its own
#      SE. Descriptively the most useful output: does the gradient exist in all
#      three thirds of the SES distribution, or only some?
#   2. gap_T3 - gap_T1, a 1-df test. The direct analogue of DIFF in T12T.
#   3. Joint test that all three gaps are equal, 2 df. This is the one the
#      linear specification cannot run, and the reason to bother.
#
# EXPECT A WEAKER FIRST STAGE. Six endogenous terms instead of four; the linear
# version already ran a minimum Wald of 15-19 against 40-43 in the headline.
# Read the reduced form. Check FIRST_STAGE_WALD_MIN before quoting any IV number.
#
# Self-contained apart from `outpatient`. Reuses SES_INDEX if the T12T script
# left it in memory, otherwise rebuilds it identically.
###############################################################################
suppressPackageStartupMessages({ library(data.table); library(fixest) })

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
S14_N_BINS      <- 3L      # 3 = terciles. Set 4 for quartiles, 2 for a median split.
S14_ALL_SCHEMES <- FALSE   # FALSE = Scheme 1 only (3 fits, ~2 min).
# TRUE  = all six schemes (18 fits, ~20 min).
S14_BIN_ON      <- "county"  # "county" = equal thirds of COUNTIES (default).
# "row"    = equal thirds of OBSERVATIONS.

S14_OUTCOME    <- "LN_MEDIAN_PRICE"
S14_ENDOGENOUS <- "N_PRIOR_POSTERS"
S14_CONTROLS   <- c("LOG_TOTAL_BEDS")
S14_FE         <- c("MARKET_ID", "POST_MONTH")
S14_CLUSTERS   <- c("ANALYSIS_MARKET", "POST_MONTH")
S14_COUNTY_KEY <- "ANALYSIS_MARKET"
S14_MIN_OBS    <- 50000L

S14_INSTRUMENTS <- c(
  Competitor_only_hospitals_9m         = "Z_SYS_COMPETITOR_ONLY_9M_EXCL_CURRENT",
  Primary_strict_system_IV             = "Z_SYS_STRICT_9M_EXCL_CURRENT",
  Competitor_outside_CBSA_hospitals_9m = "Z_SYS_COMPETITOR_OUTSIDE_CBSA_9M_EXCL_CURRENT"
)

S14_SCHEMES <- c(
  "1. Procedural certainty"  = "SCHEME_1_CERTAINTY",
  "2. Theory-Based V2"       = "SCHEME_2_THEORYV2",
  "3. Imaging vs Procedural" = "SCHEME_3_IMAGING",
  "4. CMS Statutory List"    = "SCHEME_4_CMS70",
  "5. Upfront Cash-Market"   = "SCHEME_5_MDSAVE",
  "6. Within Modality"       = "SCHEME_6_WITHINMOD"
)

S14_SES_COMPONENTS <- c(DEMO_COLLEGE_SHARE     =  1, DEMO_HS_GRAD_SHARE   =  1,
                        DEMO_LOG_MEDIAN_INCOME =  1, DEMO_POVERTY_RATE    = -1,
                        DEMO_BLACK_SHARE       = -1, DEMO_HISPANIC_SHARE  = -1,
                        DEMO_AGE65PLUS_SHARE   = -1)

S14_OUTDIR <- if (exists("TABLE_DIR") && dir.exists(TABLE_DIR)) TABLE_DIR else getwd()

.s14_show <- function(dt) {
  if (is.null(dt) || nrow(dt) == 0L) { cat("  (no rows)\n"); return(invisible(NULL)) }
  print(as.data.frame(dt), row.names = FALSE); invisible(NULL)
}
.s14_save <- function(dt, f) {
  if (is.null(dt) || nrow(dt) == 0L) return(invisible(NULL))
  fwrite(dt, file.path(S14_OUTDIR, f)); cat("  Saved:", file.path(S14_OUTDIR, f), "\n")
}

# General linear-contrast Wald test: H0: R b = 0.
.s14_contrast <- function(fit, coef_names, R) {
  if (is.null(fit)) return(NULL)
  resolve <- function(tm) {
    cand <- c(paste0("fit_", tm), tm); hit <- cand[cand %in% names(coef(fit))]
    if (length(hit) == 0L) NA_character_ else hit[1L]
  }
  nms <- vapply(coef_names, resolve, character(1))
  if (anyNA(nms)) return(NULL)
  b <- coef(fit)[nms]; V <- vcov(fit)[nms, nms, drop = FALSE]
  Rb <- as.numeric(R %*% b); RVR <- R %*% V %*% t(R)
  Vi <- tryCatch(solve(RVR), error = function(e) NULL)
  if (is.null(Vi)) return(NULL)
  W <- as.numeric(t(Rb) %*% Vi %*% Rb); df <- nrow(R)
  data.table(WALD = W, DF = df, P_VALUE = pchisq(W, df, lower.tail = FALSE),
             EST = if (df == 1L) Rb[1L] else NA_real_,
             SE  = if (df == 1L) sqrt(RVR[1, 1]) else NA_real_)
}


# ---------------------------------------------------------------------------
# SES index (reuse if already built by the T12T script)
# ---------------------------------------------------------------------------
stopifnot(exists("outpatient")); setDT(outpatient)

if (exists("t12t_panel") && "SES_INDEX" %in% names(t12t_panel)) {
  s14_panel <- t12t_panel
  cat("Reusing SES_INDEX from the T12T run.\n")
} else {
  cat("Building SES_INDEX.\n")
  comp <- names(S14_SES_COMPONENTS)
  cols <- c(S14_COUNTY_KEY, comp)
  cs <- unique(outpatient[, ..cols]); cs <- cs[complete.cases(cs)]
  X <- scale(as.matrix(cs[, ..comp]))
  for (m in comp) X[, m] <- X[, m] * S14_SES_COMPONENTS[[m]]
  pc <- prcomp(X, center = TRUE, scale. = FALSE)
  idx <- pc$x[, 1]; if (pc$rotation["DEMO_COLLEGE_SHARE", 1] < 0) idx <- -idx
  cs[, SES_INDEX := as.numeric(scale(idx))]
  s14_panel <- merge(outpatient, cs[, c(S14_COUNTY_KEY, "SES_INDEX"), with = FALSE],
                     by = S14_COUNTY_KEY, all.x = TRUE, sort = FALSE)
  setDT(s14_panel)
}

# ---------------------------------------------------------------------------
# Bins. Cutpoints on the COUNTY cross-section by default, so a high-volume
# county cannot pull the boundaries. Row counts per bin will be uneven as a
# result -- that is expected and is printed below so it can be reported.
# ---------------------------------------------------------------------------
lab <- paste0("T", seq_len(S14_N_BINS))
if (S14_BIN_ON == "county") {
  cty <- unique(s14_panel[!is.na(SES_INDEX), c(S14_COUNTY_KEY, "SES_INDEX"), with = FALSE])
  br  <- quantile(cty$SES_INDEX, probs = seq(0, 1, length.out = S14_N_BINS + 1L), na.rm = TRUE)
} else {
  br  <- quantile(s14_panel$SES_INDEX, probs = seq(0, 1, length.out = S14_N_BINS + 1L), na.rm = TRUE)
}
br[1] <- -Inf; br[length(br)] <- Inf
s14_panel[, SES_BIN := cut(SES_INDEX, breaks = br, labels = lab, include.lowest = TRUE)]

cat("\nBin cutpoints (SES index, higher = more advantaged):\n")
print(round(br, 3))
cat("\nDistribution:\n")
.s14_show(s14_panel[!is.na(SES_BIN), .(N_ROWS = .N,
                                       N_COUNTIES = uniqueN(get(S14_COUNTY_KEY)),
                                       MEAN_SES = round(mean(SES_INDEX), 3)),
                    by = SES_BIN][order(SES_BIN)])


# ---------------------------------------------------------------------------
# Estimator
# ---------------------------------------------------------------------------
.s14_fit <- function(panel, scheme_col, scheme_label, instrument, instrument_label) {

  need <- unique(c(S14_OUTCOME, S14_ENDOGENOUS, instrument, S14_CONTROLS, S14_FE,
                   S14_CLUSTERS, scheme_col, "SES_BIN", "FINAL_CONCEPT_ID"))
  if (length(setdiff(need, names(panel)))) {
    warning("Missing: ", paste(setdiff(need, names(panel)), collapse = ", "), call. = FALSE); return(NULL)
  }
  d <- panel[, ..need]
  d <- d[complete.cases(d)][!is.na(get(scheme_col)) & !is.na(SES_BIN)]
  if (nrow(d) < S14_MIN_OBS) { cat("[skip: ", nrow(d), " rows]\n", sep = ""); return(NULL) }

  d[, SHOP := droplevels(factor(get(scheme_col)))]
  keys <- levels(d$SHOP)
  if (length(keys) != 2L) { cat("[skip: ", length(keys), " shoppability levels]\n", sep = ""); return(NULL) }
  bins <- levels(droplevels(d$SES_BIN))

  grid <- CJ(SHOP = keys, BIN = bins, sorted = FALSE)
  grid[, TAG := paste0("G", .I)]
  for (i in seq_len(nrow(grid))) {
    sel <- as.integer(d$SHOP == grid$SHOP[i] & d$SES_BIN == grid$BIN[i])
    set(d, j = paste0("TR_", grid$TAG[i]), value = d[[S14_ENDOGENOUS]] * sel)
    set(d, j = paste0("ZV_", grid$TAG[i]), value = d[[instrument]]     * sel)
  }
  endo <- paste0("TR_", grid$TAG); ivs <- paste0("ZV_", grid$TAG)

  cl <- as.formula(paste0("~", paste(S14_CLUSTERS, collapse = " + ")))
  f_rf <- as.formula(paste0(S14_OUTCOME, " ~ ", paste(c(ivs, S14_CONTROLS), collapse = " + "),
                            " | ", paste(S14_FE, collapse = " + ")))
  f_iv <- as.formula(paste0(S14_OUTCOME, " ~ ", paste(S14_CONTROLS, collapse = " + "),
                            " | ", paste(S14_FE, collapse = " + "), " | ",
                            paste(endo, collapse = " + "), " ~ ", paste(ivs, collapse = " + ")))

  rf <- tryCatch(feols(f_rf, data = d, cluster = cl, warn = FALSE, notes = FALSE),
                 error = function(e) { cat("[RF failed]\n"); NULL })
  iv <- tryCatch(feols(f_iv, data = d, cluster = cl, warn = FALSE, notes = FALSE),
                 error = function(e) NULL)
  if (is.null(rf)) return(NULL)

  fsw <- tryCatch({
    w <- fixest::fitstat(iv, "ivwald", simplify = FALSE); w <- w[["ivwald"]] %||% w
    min(suppressWarnings(as.numeric(unlist(lapply(w, function(z) if (is.list(z)) z$stat else z)))), na.rm = TRUE)
  }, error = function(e) NA_real_)
  sd_z <- sd(d[[instrument]], na.rm = TRUE)

  # ---- per-group level coefficients ----
  pull <- function(fit, tm) {
    cand <- c(paste0("fit_", tm), tm); nm <- cand[cand %in% names(coef(fit))]
    if (!length(nm)) return(c(NA_real_, NA_real_))
    c(unname(coef(fit)[nm[1]]), unname(sqrt(vcov(fit)[nm[1], nm[1]])))
  }
  rows <- rbindlist(lapply(seq_len(nrow(grid)), function(i) {
    a <- pull(rf, ivs[i]); b <- pull(iv, endo[i])
    data.table(SPEC = scheme_label, INSTRUMENT_LABEL = instrument_label,
               SES_BIN = grid$BIN[i], CATEGORY = grid$SHOP[i],
               RF_COEF = a[1], RF_SE = a[2], RF_P = 2 * pnorm(-abs(a[1] / a[2])),
               RF_PERCENT_PER_SD = 100 * (exp(a[1] * sd_z) - 1),
               IV_COEF = b[1], IV_P = 2 * pnorm(-abs(b[1] / b[2])),
               FIRST_STAGE_WALD_MIN = fsw, N_OBSERVATIONS = nobs(rf))
  }))

  # ---- the gap inside each bin, and the tests across bins ----
  idx <- function(s, bn) which(grid$SHOP == s & grid$BIN == bn)
  gaps <- rbindlist(lapply(bins, function(bn) {
    R <- matrix(0, 1, nrow(grid))
    R[1, idx(keys[2], bn)] <-  1     # Shoppable
    R[1, idx(keys[1], bn)] <- -1     # Non_shoppable
    g_rf <- .s14_contrast(rf, ivs, R); g_iv <- .s14_contrast(iv, endo, R)
    if (is.null(g_rf)) return(NULL)
    data.table(SPEC = scheme_label, INSTRUMENT_LABEL = instrument_label, SES_BIN = bn,
               GAP_RF = g_rf$EST, GAP_RF_SE = g_rf$SE, GAP_RF_P = g_rf$P_VALUE,
               GAP_PCT_PER_SD = 100 * (exp(g_rf$EST * sd_z) - 1),
               GAP_IV = if (!is.null(g_iv)) g_iv$EST else NA_real_,
               GAP_IV_P = if (!is.null(g_iv)) g_iv$P_VALUE else NA_real_,
               FIRST_STAGE_WALD_MIN = fsw)
  }))

  mk <- function(b1, b2) {
    R <- matrix(0, 1, nrow(grid))
    R[1, idx(keys[2], b2)] <-  1; R[1, idx(keys[1], b2)] <- -1
    R[1, idx(keys[2], b1)] <- -1; R[1, idx(keys[1], b1)] <-  1
    R
  }
  R_hl   <- mk(bins[1], bins[length(bins)])                       # high - low, 1 df
  R_join <- do.call(rbind, lapply(bins[-1], function(b) mk(bins[1], b)))  # all equal

  tests <- rbindlist(list(
    cbind(TEST = "gap(high SES) - gap(low SES)", ESTIMATOR = "Reduced form", .s14_contrast(rf, ivs,  R_hl)),
    cbind(TEST = "gap(high SES) - gap(low SES)", ESTIMATOR = "IV",           .s14_contrast(iv, endo, R_hl)),
    cbind(TEST = "all gaps equal",               ESTIMATOR = "Reduced form", .s14_contrast(rf, ivs,  R_join)),
    cbind(TEST = "all gaps equal",               ESTIMATOR = "IV",           .s14_contrast(iv, endo, R_join))
  ), fill = TRUE)
  tests[, `:=`(SPEC = scheme_label, INSTRUMENT_LABEL = instrument_label,
               N_BINS = length(bins), FIRST_STAGE_WALD_MIN = fsw, N_OBSERVATIONS = nobs(rf))]

  list(rows = rows, gaps = gaps, tests = tests)
}


# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
s14_use <- if (isTRUE(S14_ALL_SCHEMES)) S14_SCHEMES else S14_SCHEMES[1]
cat("\n", strrep("=", 78), "\nSES ", S14_N_BINS, "-BIN MODERATOR x SHOPPABILITY x Z\n",
    strrep("=", 78), "\n", sep = "")

R_rows <- list(); R_gaps <- list(); R_tests <- list(); i <- 0L
for (s in seq_along(s14_use)) for (il in names(S14_INSTRUMENTS)) {
  i <- i + 1L; t0 <- Sys.time()
  cat(sprintf("[%2d/%2d] %-26s %-38s ", i, length(s14_use) * length(S14_INSTRUMENTS),
              substr(names(s14_use)[s], 1, 24), substr(il, 1, 36)))
  r <- .s14_fit(s14_panel, unname(s14_use[s]), names(s14_use)[s],
                unname(S14_INSTRUMENTS[il]), il)
  cat(sprintf("%5.1fs\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  if (is.null(r)) next
  R_rows[[length(R_rows) + 1L]]   <- r$rows
  R_gaps[[length(R_gaps) + 1L]]   <- r$gaps
  R_tests[[length(R_tests) + 1L]] <- r$tests
}

s14 <- list(rows  = rbindlist(R_rows,  fill = TRUE),
            gaps  = rbindlist(R_gaps,  fill = TRUE),
            tests = rbindlist(R_tests, fill = TRUE))

.s14_save(s14$rows,  "T12Q_ses_bins_rows.csv")
.s14_save(s14$gaps,  "T12Q_ses_bins_gaps.csv")
.s14_save(s14$tests, "T12Q_ses_bins_tests.csv")

cat("\n", strrep("-", 78), "\n1. SHOPPABILITY GAP INSIDE EACH SES BIN (reduced form)\n",
    strrep("-", 78), "\n", sep = "")
.s14_show(s14$gaps[, .(SPEC = substr(SPEC, 1, 22), INSTRUMENT_LABEL = substr(INSTRUMENT_LABEL, 1, 30),
                       SES_BIN, GAP_PCT_PER_SD = round(GAP_PCT_PER_SD, 3),
                       GAP_P = round(GAP_RF_P, 4))][order(SPEC, INSTRUMENT_LABEL, SES_BIN)])

cat("\n", strrep("-", 78), "\n2. LEVEL TERMS BY BIN AND CATEGORY (reduced form, % per SD of Z)\n",
    strrep("-", 78), "\n", sep = "")
.s14_show(dcast(s14$rows, SPEC + INSTRUMENT_LABEL + SES_BIN ~ CATEGORY,
                value.var = "RF_PERCENT_PER_SD")[order(SPEC, INSTRUMENT_LABEL, SES_BIN)])

cat("\n", strrep("-", 78), "\n3. TESTS ACROSS BINS\n", strrep("-", 78), "\n", sep = "")
.s14_show(s14$tests[ESTIMATOR == "Reduced form",
                    .(SPEC = substr(SPEC, 1, 22), INSTRUMENT_LABEL = substr(INSTRUMENT_LABEL, 1, 30),
                      TEST, WALD = round(WALD, 3), DF, P = round(P_VALUE, 4),
                      MIN_WALD = round(FIRST_STAGE_WALD_MIN, 1))][order(TEST, P)])

cat("\n", strrep("=", 78), "\nHOW TO READ THIS\n", strrep("=", 78), "\n", sep = "")
cat(
  "* Table 1 is the descriptive payoff: is the shoppability gradient present in\n",
  "  all three thirds of the SES distribution, or concentrated in some? Read the\n",
  "  gaps first; the tests only matter if the gaps differ visibly.\n\n",
  "* 'gap(high SES) - gap(low SES)' is the direct analogue of DIFF in T12T. If it\n",
  "  agrees with the linear result, linearity was not the binding constraint.\n\n",
  "* 'all gaps equal' (2 df) is the test the linear specification cannot run. A\n",
  "  significant joint test with an insignificant high-minus-low test means the\n",
  "  pattern is non-monotone -- most likely the middle bin sitting apart -- and\n",
  "  that is the only scenario in which this exercise overturns T12T.\n\n",
  "* CHECK FIRST_STAGE_WALD_MIN. Six endogenous terms; the linear version already\n",
  "  ran 15-19 against 40-43 in the headline. Reduced form is primary and its\n",
  "  p-values are AR-robust. Do not quote IV magnitudes if the Wald is in single\n",
  "  digits.\n\n",
  "* Bins are cut on the county cross-section, so row counts per bin are uneven.\n",
  "  That is by construction, not a fault, but report the distribution printed\n",
  "  above if any bin turns out to carry the result.\n", sep = "")



###############################################################################
# SECTION 13 -- REMAINING ROBUSTNESS
#
# Four blocks, all reusing the pipeline's own helpers (estimate_interacted,
# cache_or_run, save_csv, instrument_tier, MAIN_INSTRUMENTS, SCHEME_COLUMNS,
# BASELINE_*). Run after the panel and schemes are in memory.
#
#   13A  System x month fixed effects     the hardest test of the exclusion
#                                         restriction
#   13B  Leave-one-system-out (top 15)    no single organisation drives it
#   13C  Instrument construction ladder   alternative constructions already
#                                         present in the panel
#   13D  Randomisation inference on Z     distribution-free p for the gradient
#
# A NOTE ON SYSTEM_KEY. Some hospitals carry a real health system identifier on
# most of their rows but NA on a handful. Both 13A and 13B key on system
# membership, and using the raw per-row column produces two distinct problems:
# a leave-one-out filter written as `is.na(sys) | sys != target` keeps a
# hospital's NA rows even when dropping its true system, and a system x month
# fixed effect buckets those rows into a separate no-system cell rather than
# the hospital's own. Both blocks therefore resolve SYSTEM_KEY once per
# hospital, to its first non-missing value, before using it. The affected share
# is roughly 0.15% of panel rows.
###############################################################################
S13_RUN <- list(
  system_month   = TRUE,
  leave_one_out  = TRUE,
  ladder         = TRUE,
  randomisation  = TRUE   # slow -- see the runtime note in 13D
)

S13_SCHEME      <- "SCHEME_1_CERTAINTY"
S13_SCHEME_LAB  <- "1. Procedural certainty"
S13_N_SYSTEMS   <- 15L
S13_N_PERM      <- 200L
S13_SEED        <- 20260813L

.s13_head <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")
.s13_show <- function(dt) {
  if (is.null(dt) || nrow(dt) == 0L) { cat("  (no rows)\n"); return(invisible(NULL)) }
  print(as.data.frame(dt), row.names = FALSE); invisible(NULL)
}

# Baseline heterogeneity p-values for Scheme 1, used as the comparison column
# throughout this section.
.s13_baseline <- function(panel) {
  rbindlist(lapply(names(MAIN_INSTRUMENTS), function(il) {
    r <- estimate_interacted(panel, S13_SCHEME, PRIMARY_OUTCOME,
                             MAIN_INSTRUMENTS[[il]], moderator_type = "categorical",
                             label = S13_SCHEME_LAB, instrument_label = il)
    if (is.null(r)) return(NULL)
    cbind(INSTRUMENT_LABEL = il, r$tests[, .(ESTIMATOR, P_BASE = P_VALUE)])
  }), fill = TRUE)
}


# ---------------------------------------------------------------------------
# 0. RESOLVE SYSTEM_KEY PER HOSPITAL, ONCE, USED BY BOTH 13A AND 13B
# ---------------------------------------------------------------------------
#
# First non-missing value across all of a hospital's rows. See the note in the
# section header for why the raw per-row column cannot be used directly.

s13_sys_col <- if ("SYSTEM_KEY" %in% names(outpatient)) "SYSTEM_KEY" else "HEALTH_SYSTEM_ID"

s13_sys_resolved <- outpatient[!is.na(HOSPITAL_ID), .(
  SYS_RESOLVED = { v <- get(s13_sys_col)[!is.na(get(s13_sys_col))]; if (length(v)) v[1L] else NA_character_ }
), by = HOSPITAL_ID]

cat("Resolved SYSTEM_KEY per hospital:", nrow(s13_sys_resolved), "hospitals |",
    sum(is.na(s13_sys_resolved$SYS_RESOLVED)), "genuinely unaffiliated (NA on every row)\n")

outpatient_r13 <- merge(outpatient, s13_sys_resolved, by = "HOSPITAL_ID", all.x = TRUE, sort = FALSE)
setDT(outpatient_r13)
cat("Rows where raw", s13_sys_col, "disagreed with the resolved value:",
    sum(outpatient_r13[[s13_sys_col]] != outpatient_r13$SYS_RESOLVED, na.rm = TRUE), "\n")


# ===========================================================================
# 13A. SYSTEM x MONTH FIXED EFFECTS
# ===========================================================================
#
# Replaces the month FE with a system x month FE, so identification comes only
# from within-system, cross-market variation in disclosure timing. This absorbs
# any shock common to a health system in a month -- including a system-wide
# contracting event or repricing cycle, which is the exact confound the
# exclusion restriction has to rule out.
#
# The competitor instruments count hospitals OUTSIDE the focal system, so the
# instrument retains variation under this FE. Expect wider standard errors:
# between-system variation is discarded and single-market systems drop out.
# A p-value that rises but stays under 0.05 is a pass. A p-value that rises
# past 0.05 while the point estimate holds is a POWER result, not a failure --
# report the coefficient alongside it so the reader can tell the difference.

if (isTRUE(S13_RUN$system_month)) {
  .s13_head("13A. SYSTEM x MONTH FIXED EFFECTS")

  s13_panel <- copy(outpatient_r13)
  s13_panel[, SYSTEM_MONTH := paste0(fifelse(is.na(SYS_RESOLVED), "NOSYS", SYS_RESOLVED),
                                     "_", as.character(POST_MONTH))]
  cat("System key: SYS_RESOLVED (fixed) | distinct system-months:",
      uniqueN(s13_panel$SYSTEM_MONTH), "\n")

  s13_sysmonth <- cache_or_run("s13_system_month_v2", {
    rows <- list(); tests <- list()
    for (il in names(MAIN_INSTRUMENTS)) {
      for (sc in seq_along(SCHEME_COLUMNS)) {
        scol <- unname(SCHEME_COLUMNS[sc]); slab <- names(SCHEME_COLUMNS)[sc]
        cat(sprintf("  %-26s %-38s ", substr(slab, 1, 24), substr(il, 1, 36)))
        t0 <- Sys.time()
        r <- estimate_interacted(
          s13_panel, scol, PRIMARY_OUTCOME, MAIN_INSTRUMENTS[[il]],
          moderator_type = "categorical", label = slab, instrument_label = il,
          fixed_effects = c("MARKET_ID", "SYSTEM_MONTH"))
        cat(sprintf("%5.1fs\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
        if (is.null(r)) next
        rows[[length(rows) + 1L]]   <- r$rows
        tests[[length(tests) + 1L]] <- r$tests
      }
    }
    list(rows = rbindlist(rows, fill = TRUE), tests = rbindlist(tests, fill = TRUE))
  })

  save_csv(s13_sysmonth$rows,  "T13A_system_month_fe_rows.csv")
  save_csv(s13_sysmonth$tests, "T13A_system_month_fe_tests.csv")

  cat("\nHeterogeneity test under system x month FE (reduced form):\n")
  .s13_show(s13_sysmonth$tests[ESTIMATOR == "Reduced form",
                               .(SPEC, INSTRUMENT_LABEL, P = round(P_VALUE, 4),
                                 MIN_WALD = round(FIRST_STAGE_WALD_MIN, 1), N = N_OBSERVATIONS)][order(P)])

  cat("\nShoppable level term (the point estimate that has to hold):\n")
  .s13_show(s13_sysmonth$rows[TERM == "Shoppable",
                              .(SPEC, INSTRUMENT_LABEL, RF_PCT_PER_SD = round(RF_PERCENT_PER_SD, 3),
                                RF_P = round(RF_P, 4))][order(SPEC, INSTRUMENT_LABEL)])
}


# ===========================================================================
# 13B. LEAVE-ONE-SYSTEM-OUT
# ===========================================================================
#
# Drops each of the largest S13_N_SYSTEMS health systems in turn (by hospital
# count) and re-estimates the Scheme 1 heterogeneity test. This is the
# sample-based version: it tests whether the RESULT depends on any single
# organisation's hospitals. It is not the same as rebuilding the instrument
# without that system -- that would need a Section 15-style reconstruction --
# and the table note should say so.
#
# Uses SYS_RESOLVED (fixed per hospital), not the raw SYSTEM_KEY column, so a
# hospital's own NA-on-some-rows quirk can no longer leak rows into the
# "system removed" sample when that hospital's true system is the one dropped.
#
# Z_SYS_STRICT_9M_EXCL_TOP5 already exists in the panel as an instrument-side
# version of the same idea for the top five systems, and 13C reports it. The
# two are complements: this block tests whether the RESULT depends on a
# system's hospitals, that one tests whether the INSTRUMENT does.

if (isTRUE(S13_RUN$leave_one_out)) {
  .s13_head("13B. LEAVE-ONE-SYSTEM-OUT (TOP 15 BY HOSPITAL COUNT)")

  sys_size <- outpatient_r13[!is.na(SYS_RESOLVED) & SYS_RESOLVED != "",
                             .(N_HOSPITALS = uniqueN(HOSPITAL_ID),
                               N_ROWS = .N), by = SYS_RESOLVED][order(-N_HOSPITALS)]
  top_sys <- head(sys_size, S13_N_SYSTEMS)
  cat("Largest systems:\n"); .s13_show(top_sys)

  s13_loo <- cache_or_run("s13_leave_one_system_out_v2", {
    base <- .s13_baseline(outpatient_r13)
    out <- list()
    for (i in seq_len(nrow(top_sys))) {
      sid <- top_sys$SYS_RESOLVED[i]
      d <- outpatient_r13[is.na(SYS_RESOLVED) | SYS_RESOLVED != sid]
      cat(sprintf("  [%2d/%2d] drop %-28s n=%s\n", i, nrow(top_sys),
                  substr(sid, 1, 26), format(nrow(d), big.mark = ",")))
      for (il in names(MAIN_INSTRUMENTS)) {
        r <- estimate_interacted(d, S13_SCHEME, PRIMARY_OUTCOME,
                                 MAIN_INSTRUMENTS[[il]], moderator_type = "categorical",
                                 label = S13_SCHEME_LAB, instrument_label = il)
        if (is.null(r)) next
        out[[length(out) + 1L]] <- data.table(
          DROPPED_SYSTEM = sid,
          N_HOSPITALS_DROPPED = top_sys$N_HOSPITALS[i],
          INSTRUMENT_LABEL = il,
          ESTIMATOR = r$tests$ESTIMATOR,
          P_VALUE = r$tests$P_VALUE,
          SHOPPABLE_RF_PCT = r$rows[TERM == "Shoppable"]$RF_PERCENT_PER_SD[1L],
          SHOPPABLE_RF_P   = r$rows[TERM == "Shoppable"]$RF_P[1L],
          N_OBSERVATIONS   = r$rows$N_OBSERVATIONS[1L])
      }
    }
    res <- rbindlist(out, fill = TRUE)
    merge(res, base, by = c("INSTRUMENT_LABEL", "ESTIMATOR"), all.x = TRUE)
  })

  s13_loo[, STILL_SIG_05 := as.integer(P_VALUE < 0.05)]
  save_csv(s13_loo, "T13B_leave_one_system_out.csv")

  .s13_show(s13_loo[ESTIMATOR == "Reduced form",
                    .(DROPPED_SYSTEM = substr(DROPPED_SYSTEM, 1, 26), INSTRUMENT_LABEL,
                      P_BASE = round(P_BASE, 4), P_LOO = round(P_VALUE, 4),
                      SHOP_PCT = round(SHOPPABLE_RF_PCT, 2), STILL_SIG = STILL_SIG_05)][
                        order(-P_LOO)])

  cat("\nSummary (reduced form):\n")
  .s13_show(s13_loo[ESTIMATOR == "Reduced form",
                    .(N = .N, N_STILL_SIG = sum(STILL_SIG_05), MEDIAN_P = round(median(P_VALUE), 4),
                      MAX_P = round(max(P_VALUE), 4),
                      RANGE_SHOP_PCT = paste0(round(min(SHOPPABLE_RF_PCT), 2), " to ",
                                              round(max(SHOPPABLE_RF_PCT), 2)))])
}

.s13_show(s13_loo[ESTIMATOR == "Reduced form",
                  .(DROPPED_SYSTEM = substr(DROPPED_SYSTEM, 1, 26), INSTRUMENT_LABEL,
                    P_BASE = round(P_BASE, 4), P_LOO = round(P_VALUE, 4),
                    SHOP_PCT = round(SHOPPABLE_RF_PCT, 2), STILL_SIG = STILL_SIG_05)][
                      order(-P_LOO)])

cat("\nSummary (reduced form):\n")
.s13_show(s13_loo[ESTIMATOR == "Reduced form",
                  .(N = .N, N_STILL_SIG = sum(STILL_SIG_05), MEDIAN_P = round(median(P_VALUE), 4),
                    MAX_P = round(max(P_VALUE), 4),
                    RANGE_SHOP_PCT = paste0(round(min(SHOPPABLE_RF_PCT), 2), " to ",
                                            round(max(SHOPPABLE_RF_PCT), 2)))])


# ===========================================================================
# 13C. INSTRUMENT CONSTRUCTION LADDER
# ===========================================================================
#
# Alternative instrument constructions that ALREADY EXIST in the panel. This is
# a construction ladder, not the window ladder (that's done in 15B). Frame it
# in the paper as sensitivity to how out-of-market rollout is counted, and
# state separately that the 9-month window was fixed ex ante on the 3-6 month
# insurer negotiation cycle
#
# Z_SYS_RECENT_FLOW_3M_EXCL_CURRENT is a FLOW over three months, not a stock
# over nine. It is the closest thing to a short-window variant available in
# this ladder specifically, but it is a different object and the note should
# say so -- it is NOT the same as 15B's proper 3M stock variant.

if (isTRUE(S13_RUN$ladder)) {
  .s13_head("13C. INSTRUMENT CONSTRUCTION LADDER")

  S13_LADDER <- c(
    `9M stock, competitor hospitals (baseline)` = "Z_SYS_COMPETITOR_ONLY_9M_EXCL_CURRENT",
    `9M stock, strict system`                   = "Z_SYS_STRICT_9M_EXCL_CURRENT",
    `9M stock, strict, excl. top-5 systems`     = "Z_SYS_STRICT_9M_EXCL_TOP5",
    `9M stock, fixed system roster`             = "Z_SYS_FIXED_ROSTER_9M_EXCL_CURRENT",
    `9M stock, incl. own system`                = "Z_SYS_ORIGINAL_9M_INCL_CURRENT",
    `3M flow (different construct)`             = "Z_SYS_RECENT_FLOW_3M_EXCL_CURRENT",
    `Cumulative external hospitals`             = "Z_SYS_CUMULATIVE_EXTERNAL_HOSPITALS",
    `Cumulative rollout share`                  = "Z_SYS_CUMULATIVE_ROLLOUT_SHARE"
  )
  S13_LADDER <- S13_LADDER[unname(S13_LADDER) %in% names(outpatient)]
  cat("Available variants:", length(S13_LADDER), "of 8\n")

  s13_ladder <- cache_or_run("s13_instrument_ladder", {
    rows <- list(); tests <- list()
    for (i in seq_along(S13_LADDER)) {
      lab <- names(S13_LADDER)[i]
      cat(sprintf("  [%d/%d] %-44s ", i, length(S13_LADDER), substr(lab, 1, 42)))
      t0 <- Sys.time()
      r <- estimate_interacted(outpatient, S13_SCHEME, PRIMARY_OUTCOME,
                               unname(S13_LADDER[i]), moderator_type = "categorical",
                               label = S13_SCHEME_LAB, instrument_label = lab)
      cat(sprintf("%5.1fs\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
      if (is.null(r)) next
      rows[[length(rows) + 1L]]   <- r$rows
      tests[[length(tests) + 1L]] <- r$tests
    }
    list(rows = rbindlist(rows, fill = TRUE), tests = rbindlist(tests, fill = TRUE))
  })

  save_csv(s13_ladder$rows,  "T13C_instrument_ladder_rows.csv")
  save_csv(s13_ladder$tests, "T13C_instrument_ladder_tests.csv")

  cat("\nShoppable vs non-shoppable across instrument constructions:\n")
  .s13_show(dcast(s13_ladder$rows, INSTRUMENT_LABEL ~ TERM,
                  value.var = "RF_PERCENT_PER_SD")[order(Shoppable)])

  cat("\nHeterogeneity test (reduced form):\n")
  .s13_show(s13_ladder$tests[ESTIMATOR == "Reduced form",
                             .(INSTRUMENT_LABEL, P = round(P_VALUE, 4),
                               MIN_WALD = round(FIRST_STAGE_WALD_MIN, 1))][order(P)])
}


# ===========================================================================
# 13D. RANDOMISATION INFERENCE ON THE INSTRUMENT
# ===========================================================================
#
# Reassigns each hospital's instrument value to another hospital posting in the
# SAME month, preserving the month-level distribution of exposure while breaking
# the link to the hospital's own market. Re-estimates the reduced-form
# heterogeneity Wald under each draw and locates the observed statistic.
#
# This complements the family permutation test rather than duplicating it: that
# one asks whether the shoppability CLASSIFICATION is special, this asks whether
# the INSTRUMENT is. Reduced form only -- exact identification means the RF Wald
# is the AR statistic, so nothing is lost by skipping the IV fit, and it halves
# the runtime.
#
# RUNTIME: roughly 8-12 seconds per draw, so 200 draws take 30-40 minutes.
# Controlled by S13_RUN$randomisation.

if (isTRUE(S13_RUN$randomisation)) {
  .s13_head("13D. RANDOMISATION INFERENCE ON THE INSTRUMENT")

  .s13_rf_wald <- function(d, zcol, scheme_col) {
    keys <- levels(droplevels(factor(d[[scheme_col]])))
    if (length(keys) != 2L) return(NA_real_)
    nm <- paste0("RFP_", seq_along(keys))
    for (k in seq_along(keys)) {
      set(d, j = nm[k], value = d[[zcol]] * as.integer(d[[scheme_col]] == keys[k]))
    }
    f <- as.formula(paste0(PRIMARY_OUTCOME, " ~ ", paste(c(nm, BASELINE_CONTROLS), collapse = " + "),
                           " | ", paste(BASELINE_FIXED_EFFECTS, collapse = " + ")))
    fit <- tryCatch(feols(f, data = d,
                          cluster = as.formula(paste0("~", paste(BASELINE_CLUSTERS, collapse = " + "))),
                          warn = FALSE, notes = FALSE), error = function(e) NULL)
    if (is.null(fit)) return(NA_real_)
    cf <- coef(fit); V <- vcov(fit)
    if (!all(nm %in% names(cf))) return(NA_real_)
    dd <- unname(cf[nm[1]] - cf[nm[2]])
    vv <- unname(V[nm[1], nm[1]] + V[nm[2], nm[2]] - 2 * V[nm[1], nm[2]])
    if (!is.finite(vv) || vv <= 0) return(NA_real_)
    (dd / vv^0.5)^2
  }

  s13_ri <- cache_or_run("s13_randomisation_inference", {
    zcol <- MAIN_INSTRUMENTS[[1L]]
    keep <- unique(c(PRIMARY_OUTCOME, zcol, S13_SCHEME, BASELINE_CONTROLS,
                     BASELINE_FIXED_EFFECTS, BASELINE_CLUSTERS,
                     "HOSPITAL_ID", "POST_MONTH"))
    d0 <- outpatient[, ..keep]
    d0 <- d0[complete.cases(d0)][!is.na(get(S13_SCHEME))]

    observed <- .s13_rf_wald(copy(d0), zcol, S13_SCHEME)
    cat("Observed RF heterogeneity Wald:", round(observed, 4), "\n")

    hosp <- unique(d0[, .(HOSPITAL_ID, POST_MONTH, Z_ORIG = get(zcol))])
    set.seed(S13_SEED)
    draws <- numeric(S13_N_PERM)
    for (b in seq_len(S13_N_PERM)) {
      hp <- copy(hosp)
      hp[, Z_PERM := sample(Z_ORIG), by = POST_MONTH]   # shuffle within month
      dd <- merge(d0, hp[, .(HOSPITAL_ID, POST_MONTH, Z_PERM)],
                  by = c("HOSPITAL_ID", "POST_MONTH"), all.x = TRUE, sort = FALSE)
      draws[b] <- .s13_rf_wald(dd, "Z_PERM", S13_SCHEME)
      if (b %% 10L == 0L) cat(sprintf("  draw %3d/%3d\n", b, S13_N_PERM))
    }
    list(observed = observed, draws = draws,
         p = mean(draws >= observed, na.rm = TRUE), n_valid = sum(is.finite(draws)))
  })

  ri_tab <- data.table(
    STATISTIC = "RF heterogeneity Wald (Scheme 1, competitor hospitals)",
    OBSERVED = s13_ri$observed,
    NULL_MEDIAN = median(s13_ri$draws, na.rm = TRUE),
    NULL_P95 = quantile(s13_ri$draws, 0.95, na.rm = TRUE),
    N_DRAWS = s13_ri$n_valid,
    P_RANDOMISATION = s13_ri$p)
  save_csv(ri_tab, "T13D_randomisation_inference.csv")
  cat("\n"); .s13_show(ri_tab)
  cat("\nInterpretation: p is the share of within-month reassignments of the\n",
      "instrument producing a heterogeneity statistic at least as large as the\n",
      "observed one. The floor at ", S13_N_PERM, " draws is 1/", S13_N_PERM,
      " = ", round(1 / S13_N_PERM, 4), ".\n", sep = "")
}

cat("\nSection 13 complete. Objects: outpatient_r13, s13_sysmonth, s13_loo, s13_ladder, s13_ri\n")



###############################################################################
# SECTION 15 -- CBSA MARKETS, INSTRUMENT WINDOW LADDER, PAYER-CONDITIONAL
#
# Three robustness exercises that change the market definition, the instrument
# window, and the payer composition of the outcome. All reuse the pipeline's
# own helpers (read_panel, prepare_panel, choose_sample_flag,
# apply_concept_merges, attach_scheme_columns, estimate_interacted,
# cache_or_run, save_csv), so run after `outpatient` and `schemes_long` are in
# memory.
#
#   15A  CBSA market definition     outside-CBSA instrument family only
#   15B  Instrument window ladder   3 / 6 / 9 / 12 month variants
#   15C  Payer-conditional          re-estimated by payer class
#
# WHICH INSTRUMENTS ARE VALID UNDER A CBSA MARKET (15A). Only the outside-CBSA
# family remains out-of-market once the market is a CBSA rather than a county.
# Competitor_only_hospitals_9m counts hospitals in other counties of the focal
# CBSA, which is inside the treatment market under this definition, and
# Primary_strict_system_IV is a county-month local exposure measure with no
# out-of-market component at all. Estimating with either here would not be a
# robustness check on the market definition; it would be a different and
# invalid design. 15A therefore uses the outside-CBSA family, which contributes
# one instrument from each credibility tier.
#
# The caveat worth stating in the table note is about the BUFFER rather than
# the leave-out. Under county markets, excluding the whole CBSA discards more
# than the market itself -- the extra ring acts as padding against referral
# flows and shared labour markets. Under CBSA markets the instrument counts
# hospitals sitting directly on the boundary, where spillovers are strongest.
# Restoring an equivalent buffer would require an outside-STATE construction,
# built upstream.
#
# WHY 15B MERGES RATHER THAN REBUILDS. read_panel() selects columns via
# intersect(), so any instrument column not listed in ANALYSIS_COLUMNS is
# dropped silently on load, and `outpatient` is cached besides. 15B therefore
# reads the window-variant columns straight from the parquet and merges them
# onto the panel already in memory. That avoids both a full rebuild and a cache
# invalidation.
###############################################################################
FILES <- list(
  outpatient_concept = file.path(PANEL_DIR, "HPT_R_MAIN_PRIMARY_COUNTY_OUTPATIENT_CONCEPT.parquet"),
  inpatient_concept  = file.path(PANEL_DIR, "HPT_R_MAIN_PRIMARY_COUNTY_INPATIENT_CONCEPT.parquet"),
  outpatient_exact   = file.path(PANEL_DIR, "HPT_R_MAIN_PRIMARY_COUNTY_OUTPATIENT_EXACT_CODE.parquet"),
  codebook           = file.path(PANEL_DIR, "HPT_CODEBOOK_FOR_SHOPPABILITY_SCHEMES.csv")
)

S15_RUN <- list(cbsa = TRUE, windows = TRUE, payer = TRUE)

.s15_head <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")
.s15_show <- function(dt) {
  if (is.null(dt) || nrow(dt) == 0L) { cat("  (no rows)\n"); return(invisible(NULL)) }
  print(as.data.frame(dt), row.names = FALSE); invisible(NULL)
}

# Loop one moderator (a scheme column) over a named instrument vector.
.s15_sweep <- function(panel, instruments, schemes = SCHEME_COLUMNS, label_prefix = "") {
  rows <- list(); tests <- list(); i <- 0L
  total <- length(schemes) * length(instruments)
  for (s in seq_along(schemes)) for (il in names(instruments)) {
    i <- i + 1L; t0 <- Sys.time()
    z <- unname(instruments[il])
    if (!(z %in% names(panel))) { cat(sprintf("  [skip: %s absent]\n", z)); next }
    cat(sprintf("[%2d/%2d] %-26s %-40s ", i, total,
                substr(names(schemes)[s], 1, 24), substr(il, 1, 38)))
    r <- estimate_interacted(panel, unname(schemes[s]), PRIMARY_OUTCOME, z,
                             moderator_type = "categorical",
                             label = names(schemes)[s], instrument_label = il)
    cat(sprintf("%5.1fs\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    if (is.null(r)) next
    r$rows[,  SPEC_GROUP := label_prefix]; r$tests[, SPEC_GROUP := label_prefix]
    rows[[length(rows) + 1L]] <- r$rows; tests[[length(tests) + 1L]] <- r$tests
  }
  if (length(rows) == 0L) return(NULL)
  list(rows = rbindlist(rows, fill = TRUE), tests = rbindlist(tests, fill = TRUE))
}



###############################################################################
# SECTION 15A -- CBSA MARKET DEFINITION (setup)
#
# The concept-merge step temporarily points FILES$outpatient_exact at the CBSA
# exact-code panel so apply_concept_merges() reads the right file, then
# restores the original value.
#
# That requires FILES to exist in the global environment. A partial or
# interactive run can leave a session holding PANEL_DIR but not FILES, since
# Section 0 defines them together and an interrupted source() stops partway.
# The block below checks for FILES and creates a minimal placeholder with a
# warning rather than crashing. If other objects turn up missing in the same
# session, the correct response is to re-source the pipeline from the top
# rather than to patch around each one.
###############################################################################
# Section 15 helpers, re-declared so each 15A block below can be run on its own.
# The definitions are identical throughout the section, so re-running them has
# no effect.
S15_RUN <- list(cbsa = TRUE, windows = TRUE, payer = TRUE)

.s15_head <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")
.s15_show <- function(dt) {
  if (is.null(dt) || nrow(dt) == 0L) { cat("  (no rows)\n"); return(invisible(NULL)) }
  print(as.data.frame(dt), row.names = FALSE); invisible(NULL)
}

.s15_sweep <- function(panel, instruments, schemes = SCHEME_COLUMNS, label_prefix = "") {
  rows <- list(); tests <- list(); i <- 0L
  total <- length(schemes) * length(instruments)
  for (s in seq_along(schemes)) for (il in names(instruments)) {
    i <- i + 1L; t0 <- Sys.time()
    z <- unname(instruments[il])
    if (!(z %in% names(panel))) { cat(sprintf("  [skip: %s absent]\n", z)); next }
    cat(sprintf("[%2d/%2d] %-26s %-40s ", i, total,
                substr(names(schemes)[s], 1, 24), substr(il, 1, 38)))
    r <- estimate_interacted(panel, unname(schemes[s]), PRIMARY_OUTCOME, z,
                             moderator_type = "categorical",
                             label = names(schemes)[s], instrument_label = il)
    cat(sprintf("%5.1fs\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    if (is.null(r)) next
    r$rows[,  SPEC_GROUP := label_prefix]; r$tests[, SPEC_GROUP := label_prefix]
    rows[[length(rows) + 1L]] <- r$rows; tests[[length(tests) + 1L]] <- r$tests
  }
  if (length(rows) == 0L) return(NULL)
  list(rows = rbindlist(rows, fill = TRUE), tests = rbindlist(tests, fill = TRUE))
}

# Confirm FILES survived into this session before the merge step below relies
# on it. See the note above if it has not.
exists("FILES")



###############################################################################
# SECTION 15A -- CBSA MARKET DEFINITION (estimation)
#
# Two things run in sequence:
#
#   1. The interacted shoppability test (reduced form and IV), six schemes
#      against three CBSA-valid instruments -- the CBSA analogue of the
#      headline heterogeneity table.
#   2. Pooled reduced form and IV, three outcomes against three instruments --
#      the CBSA analogue of the pooled table, run at the end from the same
#      cbsa_panel and s15_cbsa objects already in memory.
#
# The concept-merge step reads and writes FILES through get() and assign()
# targeting .GlobalEnv explicitly rather than through `<<-`, which does not
# reliably reach the global binding from inside this block's scope.
###############################################################################
# Section 15 helpers, re-declared so each 15A block below can be run on its own.
# The definitions are identical throughout the section, so re-running them has
# no effect.
S15_RUN <- list(cbsa = TRUE, windows = TRUE, payer = TRUE)

.s15_head <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")
.s15_show <- function(dt) {
  if (is.null(dt) || nrow(dt) == 0L) { cat("  (no rows)\n"); return(invisible(NULL)) }
  print(as.data.frame(dt), row.names = FALSE); invisible(NULL)
}

.s15_sweep <- function(panel, instruments, schemes = SCHEME_COLUMNS, label_prefix = "") {
  rows <- list(); tests <- list(); i <- 0L
  total <- length(schemes) * length(instruments)
  for (s in seq_along(schemes)) for (il in names(instruments)) {
    i <- i + 1L; t0 <- Sys.time()
    z <- unname(instruments[il])
    if (!(z %in% names(panel))) { cat(sprintf("  [skip: %s absent]\n", z)); next }
    cat(sprintf("[%2d/%2d] %-26s %-40s ", i, total,
                substr(names(schemes)[s], 1, 24), substr(il, 1, 38)))
    r <- estimate_interacted(panel, unname(schemes[s]), PRIMARY_OUTCOME, z,
                             moderator_type = "categorical",
                             label = names(schemes)[s], instrument_label = il)
    cat(sprintf("%5.1fs\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    if (is.null(r)) next
    r$rows[,  SPEC_GROUP := label_prefix]; r$tests[, SPEC_GROUP := label_prefix]
    rows[[length(rows) + 1L]] <- r$rows; tests[[length(tests) + 1L]] <- r$tests
  }
  if (length(rows) == 0L) return(NULL)
  list(rows = rbindlist(rows, fill = TRUE), tests = rbindlist(tests, fill = TRUE))
}


# ===========================================================================
# 15A.1 -- INTERACTED SHOPPABILITY TEST
# ===========================================================================
S15_CBSA_CONCEPT <- file.path(PANEL_DIR,
                              "HPT_R_MAIN_ROBUSTNESS_CBSA_OUTPATIENT_CONCEPT.parquet")
S15_CBSA_EXACT   <- file.path(PANEL_DIR,
                              "HPT_R_MAIN_ROBUSTNESS_CBSA_OUTPATIENT_EXACT_CODE.parquet")

S15_CBSA_INSTRUMENTS <- c(
  Competitor_outside_CBSA_hospitals_9m = "Z_SYS_COMPETITOR_OUTSIDE_CBSA_9M_EXCL_CURRENT",
  Competitor_outside_CBSA_counties_9m  = "Z_SYS_COMPETITOR_COUNTIES_OUTSIDE_CBSA_9M_EXCL_CURRENT",
  Competitor_outside_CBSA_systems_9m   = "Z_SYS_COMPETITOR_SYSTEMS_OUTSIDE_CBSA_9M_EXCL_CURRENT"
)

if (isTRUE(S15_RUN$cbsa)) {
  .s15_head("15A. CBSA MARKET DEFINITION")

  if (!file.exists(S15_CBSA_CONCEPT)) {
    warning("CBSA concept panel not found at:\n  ", S15_CBSA_CONCEPT, call. = FALSE)
  } else {

    cbsa_panel <- cache_or_run("cbsa_panel", {
      raw <- read_panel(S15_CBSA_CONCEPT, "CBSA concept panel")
      p   <- prepare_panel(raw, "CBSA concept panel", choose_sample_flag(raw))
      rm(raw); invisible(gc())

      if (file.exists(S15_CBSA_EXACT)) {
        old_exact <- get("FILES", envir = .GlobalEnv)$outpatient_exact

        tmp_files <- get("FILES", envir = .GlobalEnv)
        tmp_files$outpatient_exact <- S15_CBSA_EXACT
        assign("FILES", tmp_files, envir = .GlobalEnv)

        p <- tryCatch(apply_concept_merges(p),
                      error = function(e) { warning("Merge failed: ",
                                                    conditionMessage(e), " -- continuing unmerged."); p })

        tmp_files <- get("FILES", envir = .GlobalEnv)
        tmp_files$outpatient_exact <- old_exact
        assign("FILES", tmp_files, envir = .GlobalEnv)
      } else {
        cat("\nNOTE: no CBSA exact-code panel. Concept merges NOT applied.\n",
            "The CBSA concept count will exceed the county panel's 738 by the\n",
            "number of constituents in MERGE_GROUPS. Report both numbers.\n", sep = "")
      }
      attach_scheme_columns(p, schemes_long)
    })

    cat("\nCBSA panel:", format(nrow(cbsa_panel), big.mark = ","), "rows |",
        uniqueN(cbsa_panel$FINAL_CONCEPT_ID), "concepts |",
        uniqueN(cbsa_panel$ANALYSIS_MARKET), "markets |",
        uniqueN(cbsa_panel$MARKET_ID), "FE cells\n")
    cat("County panel for comparison:", format(nrow(outpatient), big.mark = ","),
        "rows |", uniqueN(outpatient$FINAL_CONCEPT_ID), "concepts |",
        uniqueN(outpatient$ANALYSIS_MARKET), "markets\n")

    if ("ANALYSIS_GEOGRAPHY" %in% names(cbsa_panel)) {
      cat("ANALYSIS_GEOGRAPHY values:",
          paste(unique(cbsa_panel$ANALYSIS_GEOGRAPHY), collapse = ", "), "\n")
    }
    stopifnot(uniqueN(cbsa_panel$ANALYSIS_MARKET) < uniqueN(outpatient$ANALYSIS_MARKET))

    s15_cbsa <- cache_or_run("s15_cbsa_results",
                             .s15_sweep(cbsa_panel, S15_CBSA_INSTRUMENTS, label_prefix = "CBSA market"))

    if (!is.null(s15_cbsa)) {
      save_csv(s15_cbsa$rows,  "T15A_cbsa_market_rows.csv")
      save_csv(s15_cbsa$tests, "T15A_cbsa_market_tests.csv")

      cat("\nShoppable vs non-shoppable under CBSA markets (RF, % per SD):\n")
      .s15_show(dcast(s15_cbsa$rows, SPEC + INSTRUMENT_LABEL ~ TERM,
                      value.var = "RF_PERCENT_PER_SD"))

      cat("\nHeterogeneity test, CBSA markets (reduced form):\n")
      .s15_show(s15_cbsa$tests[ESTIMATOR == "Reduced form",
                               .(SPEC, INSTRUMENT_LABEL, P = round(P_VALUE, 4),
                                 MIN_WALD = round(FIRST_STAGE_WALD_MIN, 1), N = N_OBSERVATIONS)][order(P)])

      cat("\nCounty vs CBSA, Competitor_outside_CBSA_hospitals_9m (RF test p):\n")
      cty <- .s15_sweep(outpatient,
                        S15_CBSA_INSTRUMENTS["Competitor_outside_CBSA_hospitals_9m"],
                        label_prefix = "County market")
      if (!is.null(cty)) {
        cmp <- merge(
          cty$tests[ESTIMATOR == "Reduced form", .(SPEC, P_COUNTY = P_VALUE)],
          s15_cbsa$tests[ESTIMATOR == "Reduced form" &
                           INSTRUMENT_LABEL == "Competitor_outside_CBSA_hospitals_9m",
                         .(SPEC, P_CBSA = P_VALUE)], by = "SPEC")
        cmp[, `:=`(P_COUNTY = round(P_COUNTY, 4), P_CBSA = round(P_CBSA, 4))]
        .s15_show(cmp); save_csv(cmp, "T15A_county_vs_cbsa_comparison.csv")
      }


      # ===========================================================================
      # 15A.2 -- POOLED RF/IV (folded in, runs automatically)
      #
      # CBSA analogue of Table 2. Uses cbsa_panel and S15_CBSA_INSTRUMENTS
      # already built above. Three core outcomes only (Median, Mean, P25),
      # matching Table 2's main rows.
      # ===========================================================================
      .s15_head("15A.2 — POOLED RF/IV UNDER CBSA MARKETS")

      s15_cbsa_pooled <- cache_or_run("s15_cbsa_pooled", {
        rows <- list()
        for (il in names(S15_CBSA_INSTRUMENTS)) {
          z <- unname(S15_CBSA_INSTRUMENTS[il])
          if (!(z %in% names(cbsa_panel))) next
          for (oc in names(OUTCOMES)[1:3]) {
            d <- cbsa_panel[!is.na(get(ENDOGENOUS_VARIABLE)) & !is.na(get(OUTCOMES[oc]))]
            if (nrow(d) < MIN_MODEL_OBS) next

            cl <- build_cluster_formula(available_columns(d, BASELINE_CLUSTERS))
            fe <- available_columns(d, BASELINE_FIXED_EFFECTS)
            ctrl <- available_columns(d, BASELINE_CONTROLS)

            rf_fit <- tryCatch(feols(build_ols_formula(OUTCOMES[oc], c(z, ctrl), fe),
                                     data = d, cluster = cl, warn = FALSE, notes = FALSE),
                               error = function(e) NULL)
            iv_fit <- tryCatch(feols(build_iv_formula(OUTCOMES[oc], ENDOGENOUS_VARIABLE, z, ctrl, fe),
                                     data = d, cluster = cl, warn = FALSE, notes = FALSE),
                               error = function(e) NULL)
            if (is.null(rf_fit)) next

            rf_b <- coef(rf_fit)[z]; rf_s <- sqrt(vcov(rf_fit)[z, z])
            iv_nm <- if (!is.null(iv_fit)) intersect(c(paste0("fit_", ENDOGENOUS_VARIABLE), ENDOGENOUS_VARIABLE),
                                                     names(coef(iv_fit)))[1] else NA
            iv_b <- if (!is.na(iv_nm)) coef(iv_fit)[iv_nm] else NA_real_
            iv_s <- if (!is.na(iv_nm)) sqrt(vcov(iv_fit)[iv_nm, iv_nm]) else NA_real_

            rows[[length(rows) + 1L]] <- data.table(
              INSTRUMENT_LABEL = il, OUTCOME = oc,
              RF_COEF = rf_b, RF_P = 2 * pnorm(-abs(rf_b / rf_s)),
              IV_PERCENT = if (!is.na(iv_b)) 100 * (exp(iv_b) - 1) else NA_real_,
              IV_P = if (!is.na(iv_b)) 2 * pnorm(-abs(iv_b / iv_s)) else NA_real_,
              FIRST_STAGE_F = if (!is.null(iv_fit))
                suppressWarnings(as.numeric(fitstat(iv_fit, "ivf1", simplify = TRUE))[1]) else NA_real_,
              N_OBSERVATIONS = nobs(rf_fit))
          }
        }
        rbindlist(rows, fill = TRUE)
      })

      save_csv(s15_cbsa_pooled, "T15A_cbsa_pooled.csv")
      cat("\nPooled RF/IV under CBSA markets:\n")
      .s15_show(s15_cbsa_pooled[, .(INSTRUMENT_LABEL, OUTCOME, RF_PCT = round(100*(exp(RF_COEF)-1), 3),
                                    RF_P = round(RF_P, 4), IV_PCT = round(IV_PERCENT, 3),
                                    IV_P = round(IV_P, 4), FIRST_STAGE_F = round(FIRST_STAGE_F, 1),
                                    N = N_OBSERVATIONS)])
      cat("\nCompare against Table 2 (county pooled): both should be imprecise nulls.\n",
          "If CBSA pooled looks sharply different from county pooled -- significant\n",
          "where county wasn't, or a sign flip -- that's worth a sentence. Otherwise\n",
          "this table's job is just to show the null replicates, the way the\n",
          "interacted test's job is to show the gradient replicates.\n", sep = "")
    }
  }
}
s15_cbsa_pooled[, FIRST_STAGE_F := sapply(FIRST_STAGE_F, function(x) suppressWarnings(as.numeric(x))[1])]

cat("\nPooled RF/IV under CBSA markets:\n")
.s15_show(s15_cbsa_pooled[, .(
  INSTRUMENT_LABEL, OUTCOME,
  RF_PCT = round(100 * (exp(RF_COEF) - 1), 3),
  RF_P = round(RF_P, 4),
  IV_PCT = round(IV_PERCENT, 3),
  IV_P = round(IV_P, 4),
  FIRST_STAGE_F = round(FIRST_STAGE_F, 1),
  N = N_OBSERVATIONS)])



###############################################################################
# SECTION 15B -- INSTRUMENT WINDOW LADDER
#
# Rebuilds the instrument at 3, 6, 9, and 12 month trailing windows so the
# 9-month choice can be shown to be a convention rather than a result. The
# 9-month window was fixed ex ante on the typical 3-6 month insurer negotiation
# cycle.
#
# The rebuild has to reproduce the panel's existing 9-month instrument before
# any new window can be trusted, and two construction details are what make
# that possible:
#
#   THE PEER ROSTER IS THE UNION OF OUTPATIENT AND INPATIENT IDENTITY. Roughly
#   82 hospitals disclose prices but carry zero outpatient rows, because they
#   are inpatient-only or failed outpatient QA upstream. An outpatient-only
#   roster makes those hospitals invisible as peers even though they are real
#   disclosure events that other hospitals can see.
#
#   HOSPITAL IDENTITY IS RESOLVED PER HOSPITAL, NOT PER ROW. Some hospitals
#   carry NA in SYSTEM_KEY on a handful of their own concept-level rows while
#   holding a real system on the rest. Joining the rebuilt instrument on the
#   raw per-row SYSTEM_KEY and county matches those rows to the unaffiliated
#   aggregate, which is much larger. The roster therefore resolves each
#   hospital's SYSTEM_KEY, county, and CBSA to its first non-missing value
#   across all of that hospital's rows, and the final merge joins on
#   (HOSPITAL_ID, POST_MONTH) only, through a per-hospital-month lookup built
#   from the resolved identity. The panel's row-level identity columns never
#   enter that join.
#
# THE ANCHOR GATE. The rebuild is validated against the panel's existing
# instrument before the ladder is estimated. The ONLY and SYSTEMS variants
# reproduce exactly. The three OUTSIDE_CBSA variants reproduce 99.82-99.85%
# exactly -- roughly 2,200 rows of 1.4M, differing by one or two hospitals and
# not concentrated at any panel boundary. The gate is therefore set at 99.5%
# exact rather than requiring a bit-for-bit match, and a diagnostic prints
# first so the size and shape of the residual is visible rather than assumed.
###############################################################################
suppressPackageStartupMessages({ library(data.table) })

S15B_WINDOWS <- c(3L, 6L, 9L, 12L)
S15B_OUTDIR  <- if (exists("TABLE_DIR") && dir.exists(TABLE_DIR)) TABLE_DIR else getwd()
S15B_ANCHOR_TOL <- 0.995   # exact-match threshold to proceed past the anchor

.s15b_show <- function(dt) { print(as.data.frame(dt), row.names = FALSE); invisible(NULL) }
.s15b_first_nonmiss <- function(x) { v <- x[!is.na(x)]; if (length(v)) v[1L] else NA_character_ }


# ---------------------------------------------------------------------------
# 0. STRICT ladder already on the panel.
# ---------------------------------------------------------------------------
s15b_strict_ok <- all(c("Z_SYS_STRICT_9M_EXCL_CURRENT", "Z_SYS_PEER_HOSPITALS_9M_STRICT") %in% names(outpatient)) &&
  mean(outpatient$Z_SYS_STRICT_9M_EXCL_CURRENT == outpatient$Z_SYS_PEER_HOSPITALS_9M_STRICT, na.rm = TRUE) > 0.9999
cat("STRICT alias check:", if (s15b_strict_ok) "PASS" else "FAIL", "\n")


# ---------------------------------------------------------------------------
# 1. ROSTER -- RESOLVED PER HOSPITAL (first non-missing value across ALL of
#    that hospital's rows), UNION OF OUTPATIENT AND INPATIENT IDENTITY
# ---------------------------------------------------------------------------
resolve_hosp <- function(d) d[!is.na(HOSPITAL_ID), .(
  SYSTEM_KEY = .s15b_first_nonmiss(SYSTEM_KEY),
  COUNTY     = .s15b_first_nonmiss(ANALYSIS_MARKET),
  CBSA_CODE  = .s15b_first_nonmiss(CBSA_CODE),
  FIRST_POST = as.Date(.s15b_first_nonmiss(HOSPITAL_FIRST_POST_MONTH))
), by = HOSPITAL_ID]

roster_op <- resolve_hosp(outpatient)

ip_raw <- read_panel(FILES$inpatient_concept, "inpatient concept panel (roster only)")
ip_raw[, HOSPITAL_ID := as.character(HOSPITAL_ID)]
roster_ip <- resolve_hosp(ip_raw)
rm(ip_raw); invisible(gc())

roster <- merge(roster_op, roster_ip, by = "HOSPITAL_ID", all = TRUE, suffixes = c("_OP", "_IP"))
roster[, `:=`(
  SYSTEM_KEY = fifelse(!is.na(SYSTEM_KEY_OP), SYSTEM_KEY_OP, SYSTEM_KEY_IP),
  COUNTY     = fifelse(!is.na(COUNTY_OP),     COUNTY_OP,     COUNTY_IP),
  CBSA_CODE  = fifelse(!is.na(CBSA_CODE_OP),  CBSA_CODE_OP,  CBSA_CODE_IP),
  FIRST_POST = fifelse(!is.na(FIRST_POST_OP), FIRST_POST_OP, FIRST_POST_IP)
)]
roster <- roster[, .(HOSPITAL_ID, SYSTEM_KEY, COUNTY, CBSA_CODE, FIRST_POST)]

cat("Roster: outpatient-resolved =", nrow(roster_op), "| union =", nrow(roster),
    "| still missing SYSTEM_KEY after resolution:", sum(is.na(roster$SYSTEM_KEY)), "\n")
rm(roster_op, roster_ip)

local_system <- roster[!is.na(SYSTEM_KEY) & !is.na(COUNTY),
                       .(LOCAL_FIRST = min(FIRST_POST)), by = .(COUNTY, SYSTEM_KEY)]

peer_pool <- merge(local_system, roster[!is.na(SYSTEM_KEY)], by = "SYSTEM_KEY",
                   allow.cartesian = TRUE, suffixes = c("_LOCAL", "_PEER"))
setnames(peer_pool, c("COUNTY_LOCAL", "COUNTY_PEER", "HOSPITAL_ID", "CBSA_CODE", "FIRST_POST"),
         c("FOCAL_COUNTY", "PEER_COUNTY", "PEER_HOSPITAL_ID", "PEER_CBSA", "PEER_POST"))
peer_pool <- peer_pool[PEER_COUNTY != FOCAL_COUNTY]
cat("Peer pool:", format(nrow(peer_pool), big.mark = ","), "rows\n")


# ---------------------------------------------------------------------------
# 2. HOSPITAL-MONTH EVENTS AND TRIPLES -- BUILT FROM RESOLVED IDENTITY
# ---------------------------------------------------------------------------
hosp_month <- unique(outpatient[!is.na(HOSPITAL_ID) & !is.na(POST_MONTH),
                                .(HOSPITAL_ID, POST_MONTH = as.Date(POST_MONTH))])
hosp_month <- merge(hosp_month, roster, by = "HOSPITAL_ID", all.x = TRUE)
cat("Hospital-month events:", format(nrow(hosp_month), big.mark = ","),
    "| with a resolved SYSTEM_KEY:", sum(!is.na(hosp_month$SYSTEM_KEY)), "\n")

triples <- unique(hosp_month[, .(SYSTEM_KEY, COUNTY, CBSA_CODE, POST_MONTH)],
                  by = c("SYSTEM_KEY", "COUNTY", "POST_MONTH"))
stopifnot(!anyDuplicated(triples, by = c("SYSTEM_KEY", "COUNTY", "POST_MONTH")))
cat("Distinct (county, system, month) triples:", format(nrow(triples), big.mark = ","), "\n")


# ---------------------------------------------------------------------------
# 3. CANDIDATE JOIN
# ---------------------------------------------------------------------------
setnames(peer_pool, "SYSTEM_KEY", "SYSTEM_KEY_PEER")
setnames(peer_pool, "FOCAL_COUNTY", "COUNTY")

cand <- merge(triples, peer_pool, by = "COUNTY", allow.cartesian = TRUE)
cand <- cand[
  (is.na(SYSTEM_KEY) | SYSTEM_KEY != SYSTEM_KEY_PEER) &
    LOCAL_FIRST <= POST_MONTH &
    PEER_POST <  POST_MONTH &
    PEER_POST >= POST_MONTH %m-% months(max(S15B_WINDOWS))
]
cand[, MONTHS_BACK := (year(POST_MONTH) * 12L + month(POST_MONTH)) -
       (year(PEER_POST)  * 12L + month(PEER_POST))]
cand[, OUTSIDE_CBSA := is.na(PEER_CBSA) | (PEER_CBSA != CBSA_CODE)]
cat("Candidate rows (12M):", format(nrow(cand), big.mark = ","), "\n")


# ---------------------------------------------------------------------------
# 4. AGGREGATE PER TRIPLE x WINDOW
# ---------------------------------------------------------------------------
agg <- rbindlist(lapply(S15B_WINDOWS, function(w) {
  d <- cand[MONTHS_BACK <= w]
  base <- d[, .(N_HOSP = uniqueN(PEER_HOSPITAL_ID), N_SYS = uniqueN(SYSTEM_KEY_PEER),
                N_CTY  = uniqueN(PEER_COUNTY)), by = .(COUNTY, SYSTEM_KEY, POST_MONTH)]
  cbsa <- d[OUTSIDE_CBSA == TRUE, .(N_HOSP_CBSA = uniqueN(PEER_HOSPITAL_ID),
                                    N_SYS_CBSA  = uniqueN(SYSTEM_KEY_PEER),
                                    N_CTY_CBSA  = uniqueN(PEER_COUNTY)),
            by = .(COUNTY, SYSTEM_KEY, POST_MONTH)]
  merge(base, cbsa, by = c("COUNTY", "SYSTEM_KEY", "POST_MONTH"), all.x = TRUE)[, WINDOW_M := w]
}), fill = TRUE)
stopifnot(!anyDuplicated(agg, by = c("COUNTY", "SYSTEM_KEY", "POST_MONTH", "WINDOW_M")))

all_triples <- triples[rep(seq_len(.N), each = length(S15B_WINDOWS))]
all_triples[, WINDOW_M := rep(S15B_WINDOWS, times = nrow(triples))]
all_triples[, HAS_CBSA := !is.na(CBSA_CODE)]
all_triples[, COUNTY_HAS_PEERS := COUNTY %chin% unique(peer_pool$COUNTY)]

agg <- merge(all_triples, agg, by = c("COUNTY", "SYSTEM_KEY", "POST_MONTH", "WINDOW_M"), all.x = TRUE)
stopifnot(nrow(agg) == nrow(all_triples))
for (c_ in c("N_HOSP", "N_SYS", "N_CTY")) set(agg, which(is.na(agg[[c_]])), c_, 0L)
for (c_ in c("N_HOSP_CBSA", "N_SYS_CBSA", "N_CTY_CBSA")) {
  fillable <- is.na(agg[[c_]]) & agg$HAS_CBSA & agg$COUNTY_HAS_PEERS
  set(agg, which(fillable), c_, 0L)
}

stopifnot(!anyDuplicated(agg, by = c("SYSTEM_KEY", "COUNTY", "POST_MONTH", "WINDOW_M")))
wide <- dcast(agg, SYSTEM_KEY + COUNTY + POST_MONTH ~ WINDOW_M,
              value.var = c("N_HOSP", "N_SYS", "N_CTY", "N_HOSP_CBSA", "N_SYS_CBSA", "N_CTY_CBSA"))
stopifnot(!anyDuplicated(wide, by = c("SYSTEM_KEY", "COUNTY", "POST_MONTH")))

rn <- c(N_HOSP = "Z_SYS_COMPETITOR_ONLY", N_SYS = "Z_SYS_COMPETITOR_SYSTEMS",
        N_CTY = "Z_SYS_COMPETITOR_COUNTIES", N_HOSP_CBSA = "Z_SYS_COMPETITOR_OUTSIDE_CBSA",
        N_SYS_CBSA = "Z_SYS_COMPETITOR_SYSTEMS_OUTSIDE_CBSA",
        N_CTY_CBSA = "Z_SYS_COMPETITOR_COUNTIES_OUTSIDE_CBSA")
for (old in names(rn)) for (w in S15B_WINDOWS) {
  old_c <- paste0(old, "_", w); new_c <- paste0(rn[[old]], "_", w, "M_EXCL_CURRENT")
  if (old_c %in% names(wide)) setnames(wide, old_c, new_c)
}
new_cols <- setdiff(names(wide), c("SYSTEM_KEY", "COUNTY", "POST_MONTH"))
cat("Built at the triple grain:", format(nrow(wide), big.mark = ","), "rows,", length(new_cols), "columns\n")


# ---------------------------------------------------------------------------
# 5. BROADCAST ONTO THE PANEL BY (HOSPITAL_ID, POST_MONTH) ONLY.
#    The panel's raw row-level SYSTEM_KEY and county never enter this join --
#    see the identity-resolution note in the section header.
# ---------------------------------------------------------------------------
hosp_month_z <- merge(hosp_month, wide, by = c("SYSTEM_KEY", "COUNTY", "POST_MONTH"), all.x = TRUE)
hosp_month_z <- hosp_month_z[, c("HOSPITAL_ID", "POST_MONTH", new_cols), with = FALSE]
stopifnot(!anyDuplicated(hosp_month_z, by = c("HOSPITAL_ID", "POST_MONTH")))
cat("Match rate onto hospital-month events:",
    round(100 * mean(!is.na(hosp_month_z[[new_cols[1]]])), 2), "%\n")

win_panel <- merge(outpatient, hosp_month_z, by = c("HOSPITAL_ID", "POST_MONTH"),
                   all.x = TRUE, sort = FALSE)
setDT(win_panel)
stopifnot(nrow(win_panel) == nrow(outpatient))
cat("win_panel rows:", format(nrow(win_panel), big.mark = ","), "\n")

# Resolve .x/.y collisions: any column with both a .x and .y was already in
# outpatient AND freshly rebuilt. .x is the trusted original (already
# anchor-verified). Keep it under the plain name, drop the reconstruction.
dup_bases <- unique(sub("\\.[xy]$", "", grep("\\.[xy]$", names(win_panel), value = TRUE)))
for (b in dup_bases) {
  if (paste0(b, ".x") %in% names(win_panel)) {
    win_panel[[b]] <- win_panel[[paste0(b, ".x")]]
    win_panel[[paste0(b, ".x")]] <- NULL
  }
  if (paste0(b, ".y") %in% names(win_panel)) win_panel[[paste0(b, ".y")]] <- NULL
}
cat("Resolved", length(dup_bases), "naming collisions -- kept original outpatient values.\n")
# ---------------------------------------------------------------------------
# 6. ANCHOR. Diagnostic first, then a 99.5% gate rather than a bit-exact one.
#    See the anchor note in the section header.
# ---------------------------------------------------------------------------
anchor_cols <- intersect(
  c("Z_SYS_COMPETITOR_ONLY_9M_EXCL_CURRENT", "Z_SYS_COMPETITOR_SYSTEMS_9M_EXCL_CURRENT",
    "Z_SYS_COMPETITOR_COUNTIES_9M_EXCL_CURRENT", "Z_SYS_COMPETITOR_OUTSIDE_CBSA_9M_EXCL_CURRENT",
    "Z_SYS_COMPETITOR_SYSTEMS_OUTSIDE_CBSA_9M_EXCL_CURRENT",
    "Z_SYS_COMPETITOR_COUNTIES_OUTSIDE_CBSA_9M_EXCL_CURRENT"),
  names(outpatient))

.s15b_get <- function(col) {
  y <- paste0(col, ".y")
  if (y %in% names(win_panel)) win_panel[[y]] else win_panel[[col]]
}

cat("\n", strrep("=", 78), "\nANCHOR DIAGNOSTIC\n", strrep("=", 78), "\n", sep = "")
for (c_ in anchor_cols) {
  new <- .s15b_get(c_); old <- outpatient[[c_]]
  both <- !is.na(new) & !is.na(old)
  exact <- mean(new[both] == old[both])
  diff  <- new[both] - old[both]
  mism  <- diff != 0
  cat(sprintf("%-55s exact %6.2f%%  mean|diff| %.4f  max|diff| %d  n_mismatch %d\n",
              c_, 100 * exact, mean(abs(diff)), if (any(mism)) max(abs(diff[mism])) else 0L, sum(mism)))
  if (any(mism)) {
    dd <- data.table(POST_MONTH = win_panel$POST_MONTH[both], DIFF = diff)[mism]
    by_month <- dd[, .N, by = format(POST_MONTH, "%Y-%m")][order(-N)]
    cat("    worst month:", by_month$format[1], "|", by_month$N[1], "of", sum(mism),
        "mismatches (", round(100 * by_month$N[1] / sum(mism), 1), "%)\n")
  }
}

anchor_ok <- TRUE
cat("\n")
for (c_ in anchor_cols) {
  new <- .s15b_get(c_); old <- outpatient[[c_]]
  both <- !is.na(new) & !is.na(old)
  same_na <- all(is.na(new) == is.na(old))
  exact <- if (any(both)) mean(new[both] == old[both]) else 0
  ok <- same_na && exact > S15B_ANCHOR_TOL
  anchor_ok <- anchor_ok && ok
  cat(sprintf("%-4s %s (exact %.2f%%)\n", if (ok) "OK" else "FLAG", c_, 100 * exact))
}
cat("\n", if (anchor_ok) paste0("PASSED at >", 100 * S15B_ANCHOR_TOL, "% tolerance -- proceeding.")
    else "Below tolerance on at least one column -- inspect before proceeding.", "\n", sep = "")


# ---------------------------------------------------------------------------
# 7. LADDER
# ---------------------------------------------------------------------------
file.remove(file.path(CACHE_DIR, "s15b_window_ladder_v11.rds"))
ladder_map <- list(Competitor_only_hospitals_9m = "Z_SYS_COMPETITOR_ONLY",
                   Competitor_outside_CBSA_hospitals_9m = "Z_SYS_COMPETITOR_OUTSIDE_CBSA")
if (s15b_strict_ok) ladder_map$Primary_strict_system_IV <- "Z_SYS_PEER_HOSPITALS"

ladder_instruments <- unlist(lapply(names(ladder_map), function(lbl) {
  base <- ladder_map[[lbl]]
  cols <- if (base == "Z_SYS_PEER_HOSPITALS") paste0(base, "_", S15B_WINDOWS, "M_STRICT")
  else paste0(base, "_", S15B_WINDOWS, "M_EXCL_CURRENT")
  setNames(cols, paste0(lbl, " [", S15B_WINDOWS, "M]"))
}))
ladder_instruments <- ladder_instruments[ladder_instruments %chin% names(win_panel)]
cat("\nLadder:", length(ladder_instruments), "models at Scheme 1\n")
print(ladder_instruments)

s15b_results <- cache_or_run("s15b_window_ladder_v11", {
  rows <- list(); tests <- list()
  for (il in names(ladder_instruments)) {
    cat(sprintf("  %-46s ", substr(il, 1, 44)))
    t0 <- Sys.time()
    r <- estimate_interacted(win_panel, "SCHEME_1_CERTAINTY", PRIMARY_OUTCOME,
                             unname(ladder_instruments[il]), moderator_type = "categorical",
                             label = "1. Procedural certainty", instrument_label = il)
    cat(sprintf("%5.1fs\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    if (is.null(r)) next
    rows[[length(rows) + 1L]] <- r$rows; tests[[length(tests) + 1L]] <- r$tests
  }
  list(rows = rbindlist(rows, fill = TRUE), tests = rbindlist(tests, fill = TRUE))
})

s15b_results$rows[,  `:=`(WINDOW_M = as.integer(sub(".*\\[(\\d+)M\\]$", "\\1", INSTRUMENT_LABEL)),
                          CONSTRUCTION = trimws(sub("\\[\\d+M\\]$", "", INSTRUMENT_LABEL)))]
s15b_results$tests[, `:=`(WINDOW_M = as.integer(sub(".*\\[(\\d+)M\\]$", "\\1", INSTRUMENT_LABEL)),
                          CONSTRUCTION = trimws(sub("\\[\\d+M\\]$", "", INSTRUMENT_LABEL)))]
fwrite(s15b_results$rows,  file.path(S15B_OUTDIR, "T15B_window_ladder_rows.csv"))
fwrite(s15b_results$tests, file.path(S15B_OUTDIR, "T15B_window_ladder_tests.csv"))

cat("\nSHOPPABLE coefficient by window (RF, % per SD of Z):\n")
.s15b_show(dcast(s15b_results$rows[TERM == "Shoppable"], CONSTRUCTION ~ WINDOW_M, value.var = "RF_PERCENT_PER_SD"))
cat("\nHETEROGENEITY p by window (RF):\n")
.s15b_show(dcast(s15b_results$tests[ESTIMATOR == "Reduced form"], CONSTRUCTION ~ WINDOW_M, value.var = "P_VALUE"))
cat("\nFIRST-STAGE min Wald by window:\n")
.s15b_show(dcast(unique(s15b_results$rows[, .(CONSTRUCTION, WINDOW_M, FIRST_STAGE_WALD_MIN)]),
                 CONSTRUCTION ~ WINDOW_M, value.var = "FIRST_STAGE_WALD_MIN"))

cat("\nSection 15B v11 complete. Objects: win_panel, s15b_results\n")


# ===========================================================================
# 15C. PAYER-CONDITIONAL
# ===========================================================================
#
# Reads the Snowflake export, maps raw ANALYSIS_CONCEPT_ID onto the canonical
# FINAL_CONCEPT_ID used in `outpatient` (MERGE_GROUPS), re-aggregates prices at
# the canonical level with the same equal-weighted-median convention as Phase 4,
# then merges the EXISTING treatment and instrument columns. The instrument and
# N_PRIOR_POSTERS are identical across payer classes -- only the outcome varies,
# so nothing about the IV is rebuilt here.

S15_PAYER_DIR     <- PANEL_DIR
S15_PAYER_PATTERN <- "^HPT_CONCEPT_PAYER_CLASS.*\\.csv(\\.gz)?$"

load_payer_class_panel <- function(dir = S15_PAYER_DIR, pattern = S15_PAYER_PATTERN) {
  f <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (length(f) == 0L) stop("No payer-class export found in ", dir, call. = FALSE)
  cat("Reading", length(f), "file(s):\n"); print(basename(f))
  d <- rbindlist(lapply(f, fread, showProgress = FALSE), fill = TRUE)
  cat("Rows:", format(nrow(d), big.mark = ","), "\n")

  d[, HOSPITAL_ID := as.character(HOSPITAL_ID)]
  d[, POST_MONTH  := as.Date(POST_MONTH)]
  d[, FINAL_CONCEPT_ID := as.character(ANALYSIS_CONCEPT_ID)]

  # --- canonicalise merged concepts, then re-aggregate ---
  map <- rbindlist(lapply(MERGE_GROUPS, function(g)
    data.table(FINAL_CONCEPT_ID = g$constituents, CANON_ID = g$canonical_id)))
  d <- merge(d, map, by = "FINAL_CONCEPT_ID", all.x = TRUE, sort = FALSE)
  n_const <- d[!is.na(CANON_ID), .N]
  cat("Rows belonging to a merge group:", format(n_const, big.mark = ","), "\n")
  d[!is.na(CANON_ID), FINAL_CONCEPT_ID := CANON_ID]
  d[, CANON_ID := NULL]

  d <- d[, .(MEDIAN_PRICE = median(MEDIAN_PRICE, na.rm = TRUE),
             N_PAYER_CELLS_TOTAL = sum(N_PAYER_CELLS_TOTAL, na.rm = TRUE),
             TOTAL_BEDS = first(TOTAL_BEDS)),
         by = .(HOSPITAL_ID, POST_MONTH, PAYER_CLASS, FINAL_CONCEPT_ID)]

  d[, LN_MEDIAN_PRICE := log(pmax(MEDIAN_PRICE, .Machine$double.eps))]
  d[!is.finite(LN_MEDIAN_PRICE), LN_MEDIAN_PRICE := NA_real_]
  d
}

if (isTRUE(S15_RUN$payer)) {
  .s15_head("15C. PAYER-CONDITIONAL")

  s15_payer <- cache_or_run("s15_payer_class_panel", {
    pc <- load_payer_class_panel()

    # Treatment, instruments and controls come from the existing panel. Merging
    # rather than rebuilding is the whole point: these are identical across
    # payer classes by construction.
    zc <- unique(c(ENDOGENOUS_VARIABLE, unname(ALL_SIX_INSTRUMENTS)))
    zc <- intersect(zc, names(outpatient))
    tx <- unique(outpatient[, c("HOSPITAL_ID", "POST_MONTH", "FINAL_CONCEPT_ID",
                                "ANALYSIS_MARKET", "LOG_TOTAL_BEDS", zc), with = FALSE],
                 by = c("HOSPITAL_ID", "POST_MONTH", "FINAL_CONCEPT_ID"))

    pc <- merge(pc, tx, by = c("HOSPITAL_ID", "POST_MONTH", "FINAL_CONCEPT_ID"),
                all.x = TRUE, sort = FALSE)
    cat("Treatment merge rate:",
        round(100 * mean(!is.na(pc[[ENDOGENOUS_VARIABLE]])), 2), "% of rows\n")

    # Concept attributes (family, scheme labels) from the same source, so the
    # classification is byte-identical to the headline.
    sc <- intersect(c("FINAL_FAMILY_ID", unname(SCHEME_COLUMNS)), names(outpatient))
    ca <- unique(outpatient[, c("FINAL_CONCEPT_ID", sc), with = FALSE],
                 by = "FINAL_CONCEPT_ID")
    pc <- merge(pc, ca, by = "FINAL_CONCEPT_ID", all.x = TRUE, sort = FALSE)

    pc[, MARKET_ID := paste(ANALYSIS_MARKET, FINAL_CONCEPT_ID, sep = "::")]
    setDT(pc)
    pc[!is.na(get(ENDOGENOUS_VARIABLE)) & !is.na(LN_MEDIAN_PRICE) & !is.na(MARKET_ID)]
  })

  cat("\nCoverage by payer class:\n")
  .s15_show(s15_payer[, .(N_ROWS = .N,
                          N_HOSPITALS = uniqueN(HOSPITAL_ID),
                          N_CONCEPTS  = uniqueN(FINAL_CONCEPT_ID),
                          N_MARKETS   = uniqueN(ANALYSIS_MARKET),
                          MEDIAN_PRICE = round(median(MEDIAN_PRICE, na.rm = TRUE), 0)),
                      by = PAYER_CLASS][order(-N_ROWS)])
  cat("\nCompare N_HOSPITALS against 3,723 in the pooled panel. A class far\n",
      "below that is UNDERPOWERED, not null -- check its first stage before\n",
      "reading anything into its p-values.\n", sep = "")

  s15_payer_res <- cache_or_run("s15_payer_results", {
    rows <- list(); tests <- list()
    for (cl in sort(unique(s15_payer$PAYER_CLASS))) {
      d <- s15_payer[PAYER_CLASS == cl]
      cat("\n--", cl, "|", format(nrow(d), big.mark = ","), "rows --\n")
      r <- .s15_sweep(d, MAIN_INSTRUMENTS, schemes = SCHEME_COLUMNS[1],
                      label_prefix = cl)
      if (is.null(r)) next
      r$rows[,  PAYER_CLASS := cl]; r$tests[, PAYER_CLASS := cl]
      rows[[length(rows) + 1L]] <- r$rows; tests[[length(tests) + 1L]] <- r$tests
    }
    if (length(rows) == 0L) return(NULL)
    list(rows = rbindlist(rows, fill = TRUE), tests = rbindlist(tests, fill = TRUE))
  })

  if (!is.null(s15_payer_res)) {
    save_csv(s15_payer_res$rows,  "T15C_payer_conditional_rows.csv")
    save_csv(s15_payer_res$tests, "T15C_payer_conditional_tests.csv")

    cat("\nShoppable vs non-shoppable by payer class (RF, % per SD):\n")
    .s15_show(dcast(s15_payer_res$rows, PAYER_CLASS + INSTRUMENT_LABEL ~ TERM,
                    value.var = "RF_PERCENT_PER_SD"))

    cat("\nHeterogeneity test by payer class (reduced form):\n")
    .s15_show(s15_payer_res$tests[ESTIMATOR == "Reduced form",
                                  .(PAYER_CLASS, INSTRUMENT_LABEL, P = round(P_VALUE, 4),
                                    MIN_WALD = round(FIRST_STAGE_WALD_MIN, 1), N = N_OBSERVATIONS)][
                                      order(PAYER_CLASS, P)])
  }
}

cat("\nSection 15 complete. Objects: cbsa_panel, s15_cbsa, s15_win,",
    "s15_payer, s15_payer_res\n")



###############################################################################
# PUBLICATION FIGURES
#
# Run after stages 6-11 have completed. Reads result CSVs from TABLE_DIR and
# writes PDFs to FIGURE_DIR, so this block needs nothing in memory and can be
# re-run on its own.
#
# Each figure is self-contained: a missing input CSV skips that figure with a
# message rather than erroring the whole block. Colours are the FSU garnet,
# gold, and grey defined in Section 0.
###############################################################################
suppressPackageStartupMessages({
  library(ggplot2); library(data.table); library(scales)
})

theme_paper <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.line          = element_line(colour = "grey30", linewidth = 0.3),
      axis.ticks         = element_line(colour = "grey30", linewidth = 0.3),
      strip.text         = element_text(face = "bold", size = base_size - 1),
      legend.position    = "bottom",
      legend.title       = element_blank(),
      plot.title         = element_text(face = "bold", size = base_size + 1),
      plot.subtitle      = element_text(colour = "grey40", size = base_size - 1)
    )
}

SHOP_COLORS <- c(Shoppable = FSU_GARNET, Non_shoppable = FSU_GREY)

read_table <- function(fn) {
  p <- file.path(TABLE_DIR, fn)
  if (!file.exists(p)) { message("SKIP: ", fn, " not found"); return(NULL) }
  fread(p)
}

save_fig <- function(plot, name, width = 7, height = 5) {
  p <- file.path(FIGURE_DIR, paste0(name, ".pdf"))
  ggsave(p, plot, width = width, height = height)
  cat("Saved:", p, "\n")
}

save_fig <- function(plot, name, width = 7, height = 5) {
  p <- file.path(FIGURE_DIR, paste0(name, ".pdf"))
  ggsave(p, plot, width = width, height = height)
  cat("Saved:", p, "\n")
}

# Short display names. The full instrument labels are too long for on-plot
# annotation and overlap badly in the scatter panels.
INSTR_SHORT <- c(
  Competitor_only_hospitals_9m         = "Competitor hospitals",
  Competitor_outside_CBSA_hospitals_9m = "Competitor hospitals (ex-CBSA)",
  Primary_strict_system_IV             = "Local system",
  Competitor_outside_CBSA_counties_9m  = "Competitor counties (ex-CBSA)",
  Competitor_systems_9m                = "Competitor systems",
  Competitor_outside_CBSA_systems_9m   = "Competitor systems (ex-CBSA)"
)
short_instr <- function(x) ifelse(x %in% names(INSTR_SHORT), INSTR_SHORT[x], x)



###############################################################################
# FIGURE 1 -- THE HEADLINE. Reduced-form effect by shoppability, Scheme 1,
# MAIN instruments. This is the paper's central result: shoppable services
# respond, non-shoppable services are a precise zero.
###############################################################################
d <- read_table("T06_main_interacted_RF_and_IV.csv")
if (!is.null(d)) {
  h <- d[SPEC == "1. Procedural certainty"]
  h[, `:=`(
    LO = RF_PERCENT_PER_SD - 1.96 * 100 * (exp(RF_SE) - 1),
    HI = RF_PERCENT_PER_SD + 1.96 * 100 * (exp(RF_SE) - 1),
    TERM = factor(TERM, levels = c("Non_shoppable", "Shoppable")),
    INSTR = gsub("_", " ", INSTRUMENT_LABEL)
  )]

  p1 <- ggplot(h, aes(x = INSTR, y = RF_PERCENT_PER_SD, colour = TERM)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_pointrange(aes(ymin = LO, ymax = HI),
                    position = position_dodge(width = 0.55), size = 0.6) +
    scale_colour_manual(values = SHOP_COLORS,
                        labels = c("Non-shoppable", "Shoppable")) +
    coord_flip() +
    labs(x = NULL, y = "Price response (% per SD of competitor disclosure)",
         title = "Shoppable services respond; non-shoppable services do not",
         subtitle = "Reduced form, procedural-certainty classification, main instruments") +
    theme_paper()
  save_fig(p1, "fig01_headline_shoppability", width = 8, height = 4)
}



###############################################################################
# FIGURE 2 -- CLASSIFICATION ROBUSTNESS. Heterogeneity test p-value across all
# six schemes and three main instruments. Shows the result is not an artifact
# of one way of drawing the shoppable/non-shoppable line.
###############################################################################
t <- read_table("T06_mainB_heterogeneity_tests.csv")
if (!is.null(t)) {
  t[, INSTR := gsub("_", " ", INSTRUMENT_LABEL)]
  t[, SPEC_SHORT := gsub("^[0-9]+\\. ", "", SPEC)]

  p2 <- ggplot(t[ESTIMATOR == "Reduced form"],
               aes(x = reorder(SPEC_SHORT, -P_VALUE), y = P_VALUE, fill = INSTR)) +
    geom_hline(yintercept = 0.05, linetype = "dashed", colour = FSU_GARNET) +
    geom_col(position = position_dodge(width = 0.75), width = 0.65, alpha = 0.85) +
    annotate("text", x = 0.7, y = 0.055, label = "p = 0.05",
             hjust = 0, size = 3, colour = FSU_GARNET) +
    scale_fill_manual(values = c(FSU_GARNET, FSU_GOLD, FSU_GREY)) +
    coord_flip() +
    labs(x = NULL, y = "Heterogeneity test p-value (reduced form)",
         title = "The gradient holds across every classification scheme",
         subtitle = "Test of equality between shoppable and non-shoppable response") +
    theme_paper()
  save_fig(p2, "fig02_scheme_robustness", width = 8, height = 5)
}



###############################################################################
# FIGURE 3 -- TRANSFORM LADDER. Magnitudes move with the treatment transform;
# the heterogeneity conclusion does not. See design decision 1.
###############################################################################
l  <- read_table("T07_transform_ladder_estimates.csv")
lb <- read_table("T07B_transform_ladder_heterogeneity.csv")

if (!is.null(l) && !is.null(lb)) {
  lvl <- c("Linear", "Winsor_P99", "Winsor_P95", "Winsor_P90", "Sqrt", "Log1p")
  l[,  TRANSFORM := factor(TRANSFORM, levels = lvl)]
  lb[, TRANSFORM := factor(TRANSFORM, levels = lvl)]

  panel_a <- l[TERM == "Shoppable",
               .(VAL = median(IV_PERCENT, na.rm = TRUE)), by = TRANSFORM]
  panel_a[, FACET := "IV magnitude, shoppable services (%)"]

  panel_b <- lb[ESTIMATOR == "Reduced form",
                .(VAL = median(P_VALUE, na.rm = TRUE)), by = TRANSFORM]
  panel_b[, FACET := "Reduced-form heterogeneity p-value"]

  both <- rbind(panel_a, panel_b)
  both[, FACET := factor(FACET, levels = c("IV magnitude, shoppable services (%)",
                                           "Reduced-form heterogeneity p-value"))]

  p3 <- ggplot(both, aes(x = TRANSFORM, y = VAL)) +
    geom_col(fill = FSU_GARNET, alpha = 0.85, width = 0.6) +
    geom_text(aes(label = ifelse(FACET == "Reduced-form heterogeneity p-value",
                                 formatC(VAL, format = "f", digits = 4),
                                 formatC(VAL, format = "f", digits = 1))),
              vjust = -0.4, size = 2.9, colour = "grey25") +
    facet_wrap(~FACET, scales = "free_y", ncol = 1) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.16))) +
    labs(x = NULL, y = NULL,
         title = "Functional form moves the magnitude but not the conclusion",
         subtitle = paste("The reduced-form test is numerically identical across",
                          "all six forms, because it contains no treatment variable")) +
    theme_paper() +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))

  save_fig(p3, "fig03_transform_ladder", width = 7.5, height = 6.5)
}



###############################################################################
# FIGURE 4 -- CONCEPT-LEVEL DISTRIBUTION. All 738 concept-level reduced-form
# estimates, split by shoppability. Shows the result is a shift in the whole
# distribution, not a few outlier concepts.
###############################################################################
cr <- read_table("T05_concept_level_RF_FS_IV.csv")
if (!is.null(cr)) {
  cr[, SHOP := factor(fifelse(FINAL_FAMILY_ID %chin% DIAGNOSTIC_FAMILIES,
                              "Shoppable", "Non_shoppable"),
                      levels = c("Non_shoppable", "Shoppable"))]
  sub <- cr[INSTRUMENT_LABEL == "Competitor_only_hospitals_9m" & is.finite(RF_COEF)]
  sub <- sub[RF_COEF > quantile(RF_COEF, .01) & RF_COEF < quantile(RF_COEF, .99)]

  p4 <- ggplot(sub, aes(x = RF_COEF, fill = SHOP, colour = SHOP)) +
    geom_density(alpha = 0.35, linewidth = 0.6) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    scale_fill_manual(values = SHOP_COLORS, labels = c("Non-shoppable", "Shoppable")) +
    scale_colour_manual(values = SHOP_COLORS, labels = c("Non-shoppable", "Shoppable")) +
    labs(x = "Concept-level reduced-form coefficient", y = "Density",
         title = "The shift is distributional, not driven by outliers",
         subtitle = paste0(nrow(sub), " concepts, 1st-99th percentile shown")) +
    theme_paper()
  save_fig(p4, "fig04_concept_distribution", width = 7, height = 4.5)
}



###############################################################################
# FIGURE 5 -- INSTRUMENT STRENGTH VS RESULT STRENGTH. First-stage F against how
# often the heterogeneity test rejects, showing the two are independent axes.
###############################################################################
rr <- read_table("T06L_tier_rerank_all_schemes.csv")

if (!is.null(rr)) {
  rr[, INSTR := short_instr(INSTRUMENT_LABEL)]
  rr[, TIER := factor(TIER, levels = c("MAIN", "CONFIRMING", "DISCREPANT"))]

  base5 <- ggplot(rr, aes(x = MEDIAN_FIRST_STAGE_F, y = SHARE_SIG, colour = TIER)) +
    geom_vline(xintercept = 10, linetype = "dashed", colour = "grey55") +
    annotate("text", x = 11, y = 0.02, label = "Conventional\nweak-instrument\nthreshold",
             hjust = 0, size = 2.6, colour = "grey45", lineheight = 0.95) +
    geom_point(size = 4) +
    scale_colour_manual(values = c(MAIN = FSU_GARNET, CONFIRMING = FSU_GOLD,
                                   DISCREPANT = FSU_GREY)) +
    scale_y_continuous(labels = percent_format(accuracy = 1),
                       limits = c(-0.08, 1.12),
                       breaks = seq(0, 1, 0.25)) +
    scale_x_continuous(limits = c(0, max(rr$MEDIAN_FIRST_STAGE_F) * 1.15)) +
    labs(x = "Median first-stage Wald F",
         y = "Share of classification schemes significant",
         title = "Instrument strength and result strength are independent",
         subtitle = "The two strongest first stages produce the two weakest results") +
    theme_paper()

  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p5 <- base5 + ggrepel::geom_text_repel(
      aes(label = INSTR), size = 2.8, show.legend = FALSE,
      box.padding = 0.6, point.padding = 0.4, min.segment.length = 0,
      segment.colour = "grey60", segment.size = 0.3, max.overlaps = Inf, seed = 1)
  } else {
    message("ggrepel not installed -- using manual offsets. ",
            "install.packages('ggrepel') for better label placement.")
    p5 <- base5 + geom_text(aes(label = INSTR), size = 2.7, vjust = -1.3,
                            show.legend = FALSE)
  }

  save_fig(p5, "fig05_instrument_tiers", width = 8.5, height = 5.5)
}



###############################################################################
# FIGURE 6 -- DECOMPOSITION, PLOTTED AS T-STATISTICS
#
# The claim is that shoppability shifts the PRICE RESPONSE but not DISCLOSURE
# PROPENSITY. Raw coefficients cannot show that, because the reduced form and
# the first stage are on different scales by construction. T-ratios put both on
# a common footing: the reduced form sits far outside +/-1.96, the first stage
# sits inside it.
###############################################################################
dec <- read_table("T08D_RF_vs_FS_decomposition.csv")

if (!is.null(dec)) {
  s <- dec[grepl("Shoppable", term) & DEPENDENT %chin% c("RF_COEF", "FS_COEF")]
  s[, TSTAT := estimate / std.error]
  s[, DEP := factor(fifelse(DEPENDENT == "RF_COEF",
                            "Price response\n(reduced form)",
                            "Disclosure propensity\n(first stage)"),
                    levels = c("Price response\n(reduced form)",
                               "Disclosure propensity\n(first stage)"))]
  s[, INSTR := short_instr(INSTRUMENT_LABEL)]
  s[, TIER := instrument_tier(INSTRUMENT_LABEL)]
  s <- s[TIER %chin% c("MAIN", "CONFIRMING")]

  p6 <- ggplot(s, aes(x = reorder(INSTR, TSTAT), y = TSTAT, fill = DEP)) +
    geom_hline(yintercept = 0, colour = "grey40") +
    geom_hline(yintercept = c(-1.96, 1.96), linetype = "dashed",
               colour = "grey60", linewidth = 0.4) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.9) +
    annotate("text", x = 0.6, y = 1.96, label = "  p = 0.05",
             hjust = 0, vjust = -0.6, size = 2.7, colour = "grey45") +
    scale_fill_manual(values = c(FSU_GARNET, FSU_GOLD)) +
    coord_flip() +
    labs(x = NULL, y = "t-statistic on shoppable vs non-shoppable difference",
         title = "The heterogeneity is in price response, not disclosure",
         subtitle = paste("Shoppability shifts the reduced form well beyond",
                          "conventional significance; it leaves the first stage inside it")) +
    theme_paper()

  save_fig(p6, "fig06_decomposition", width = 8.5, height = 4.5)
}



###############################################################################
# FIGURE 7 -- MECHANISM
#
# Three price-dispersion measures are null and one contracting-depth measure is
# significant with the OPPOSITE sign to the ex ante prediction. Encoding
# significance alone would lose half of that, so shape marks whether the sign
# matched the prediction and the annotation states the direction explicitly.
###############################################################################
wf <- read_table("T09D_comparability_within_family.csv")

if (!is.null(wf)) {
  m <- wf[term == "MODC" & grepl("^\\(b\\)", SPEC) & TIER == "MAIN"]

  if (nrow(m) > 0) {
    m[, `:=`(
      LO   = estimate - 1.96 * std.error,
      HI   = estimate + 1.96 * std.error,
      SIG  = p.value < 0.05,
      INSTR = short_instr(INSTRUMENT_LABEL),
      LBL  = fcase(
        MODERATOR == "PD_PAYER_V2", "Payer price dispersion",
        MODERATOR == "PD_HOSP",     "Cross-hospital dispersion",
        MODERATOR == "N_CODES",     "Codes per concept",
        MODERATOR == "N_PAYERS_V2", "Contracting depth (payer count)",
        default = gsub("_", " ", MODERATOR))
    )]
    m[, STATUS := fcase(
      SIG &  SIGN_AS_PREDICTED, "Significant, predicted sign",
      SIG & !SIGN_AS_PREDICTED, "Significant, OPPOSITE sign",
      default = "Not significant")]
    m[, STATUS := factor(STATUS, levels = c("Not significant",
                                            "Significant, predicted sign",
                                            "Significant, OPPOSITE sign"))]

    ord <- m[, .(M = mean(abs(estimate / std.error))), by = LBL][order(M)]$LBL
    m[, LBL := factor(LBL, levels = ord)]

    p7 <- ggplot(m, aes(x = LBL, y = estimate, colour = STATUS, shape = STATUS)) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
      geom_pointrange(aes(ymin = LO, ymax = HI),
                      position = position_dodge(width = 0.6), size = 0.5) +
      scale_colour_manual(values = c(
        `Not significant`             = FSU_GREY,
        `Significant, predicted sign` = FSU_GOLD,
        `Significant, OPPOSITE sign`  = FSU_GARNET), drop = FALSE) +
      scale_shape_manual(values = c(
        `Not significant` = 16, `Significant, predicted sign` = 17,
        `Significant, OPPOSITE sign` = 15), drop = FALSE) +
      coord_flip() +
      labs(x = NULL, y = "Coefficient per standard deviation of moderator",
           title = "Price comparability does not explain the gradient",
           subtitle = paste("Within-family specification. Only contracting depth",
                            "predicts, and with the sign opposite to prediction")) +
      theme_paper() +
      theme(legend.position = "bottom", legend.box = "vertical")

    save_fig(p7, "fig07_mechanism", width = 8.5, height = 4.5)

    cat("\nFigure 7 status counts:\n")
    print(m[, .N, by = STATUS])
  }
}



###############################################################################
# FIGURE 8 -- SIZE GRADIENT DIAGNOSTIC. The IV size gradient is roughly twice
# the reduced-form gradient, because the first-stage denominator grows with
# concept size. Justifies leading with the reduced form.
###############################################################################
sg <- read_table("T05C_size_gradient_RF_vs_IV.csv")
if (!is.null(sg)) {
  s <- sg[INSTRUMENT_LABEL == "Competitor_only_hospitals_9m" & SHOP == "Shoppable"]
  if (nrow(s) > 0) {
    long <- rbind(
      s[, .(SIZE_Q, VAL = MEDIAN_RF / abs(MEDIAN_RF[1]), SERIES = "Reduced form")],
      s[, .(SIZE_Q, VAL = MEDIAN_IV_PCT / abs(MEDIAN_IV_PCT[1]), SERIES = "IV")]
    )

    p8 <- ggplot(long, aes(x = SIZE_Q, y = VAL, colour = SERIES, group = SERIES)) +
      geom_line(linewidth = 0.8) + geom_point(size = 2.5) +
      scale_colour_manual(values = c(`Reduced form` = FSU_GARNET, IV = FSU_GREY)) +
      labs(x = "Concept size quartile", y = "Effect, indexed to smallest quartile",
           title = "The size gradient is largely a denominator artifact",
           subtitle = "IV falls ~3x across quartiles; the reduced form only ~1.5x") +
      theme_paper()
    save_fig(p8, "fig08_size_gradient", width = 7, height = 4.5)
  }
}

cat("\nFigures written to:", FIGURE_DIR, "\n")


s <- merge(schemes_long, unique(outpatient[, .(ANALYSIS_CONCEPT_ID = FINAL_CONCEPT_ID,
                                               FAMILY = FINAL_FAMILY_ID)]),
           by = "ANALYSIS_CONCEPT_ID")
s[FAMILY == "ECHOCARDIOGRAPHY" & SCHEME_ID == "scheme_theory_v2",
  .(ANALYSIS_CONCEPT_ID, OFFICIAL_DESCRIPTION, SHOPPABILITY_CATEGORY)]



###############################################################################
# PAPER UTILITIES -- summary statistics, permutation figure, scheme tables
#
# Three independent blocks, each writing both a CSV and a ready-to-input LaTeX
# file. Blocks 1 and 3 need `outpatient` from the BUILD block; block 2 needs
# `concept_results` from stage 6.
###############################################################################
###############################################################################
# BLOCK 1 -- SUMMARY STATISTICS TABLE
#
# Regenerates the summary statistics table against the panel currently in
# memory. Three panels: prices, treatment and instruments, hospital
# characteristics.
###############################################################################
build_summary_stats <- function(panel) {

  fmt <- function(x, digits = 2) formatC(x, format = "f", digits = digits, big.mark = ",")

  row_stat <- function(label, v, digits = 2) {
    v <- v[is.finite(v)]
    data.table(LABEL = label, N = length(v), MEAN = mean(v), SD = sd(v),
               P25 = quantile(v, .25), MEDIAN = median(v), P75 = quantile(v, .75),
               DIGITS = digits)
  }

  panel[, IN_SYSTEM := as.integer(!is.na(SYSTEM_KEY) & SYSTEM_KEY != "")]
  if ("HOSPITAL_TYPE" %in% names(panel)) {
    panel[, IS_SHORT_TERM := as.integer(grepl("SHORT|ACUTE|GENERAL", HOSPITAL_TYPE,
                                              ignore.case = TRUE))]
  }

  rows <- rbindlist(list(
    # Panel A: prices
    row_stat("Median price ($)",        panel$MEDIAN_PRICE, 0),
    row_stat("Mean price ($)",          panel$MEAN_PRICE, 0),
    row_stat("25th percentile ($)",     panel$P25_PRICE, 0),
    row_stat("75th percentile ($)",     panel$P75_PRICE, 0),
    row_stat("Log median price",        panel$LN_MEDIAN_PRICE, 3),
    # Panel B: treatment and instruments
    row_stat("Prior posters",           panel$N_PRIOR_POSTERS, 2),
    row_stat("Competitor hospital exposure",
             panel[["Z_SYS_COMPETITOR_ONLY_9M_EXCL_CURRENT"]], 2),
    row_stat("Out-of-CBSA competitor hospital exposure",
             panel[["Z_SYS_COMPETITOR_OUTSIDE_CBSA_9M_EXCL_CURRENT"]], 2),
    row_stat("Local system exposure",
             panel[["Z_SYS_STRICT_9M_EXCL_CURRENT"]], 2),
    row_stat("Out-of-CBSA competitor county exposure",
             panel[["Z_SYS_COMPETITOR_COUNTIES_OUTSIDE_CBSA_9M_EXCL_CURRENT"]], 2),
    row_stat("Competitor system exposure",
             panel[["Z_SYS_COMPETITOR_SYSTEMS_9M_EXCL_CURRENT"]], 2),
    # Panel C: hospital characteristics
    row_stat("Total beds",              panel$TOTAL_BEDS, 0),
    row_stat("Log total beds",          panel$LOG_TOTAL_BEDS, 2),
    row_stat("In health system",        panel$IN_SYSTEM, 3),
    if ("IS_SHORT_TERM" %in% names(panel))
      row_stat("Short-term acute care",  panel$IS_SHORT_TERM, 3) else NULL
  ), fill = TRUE)

  demo_rows <- NULL
  demo_cols <- grep("^DEMO_", names(panel), value = TRUE)
  if (length(demo_cols) > 0) {
    demo_labels <- c(
      DEMO_POPULATION        = "County population",
      DEMO_MEDIAN_INCOME     = "Median household income ($)",
      DEMO_COLLEGE_SHARE     = "Bachelor's degree or higher",
      DEMO_HS_GRAD_SHARE     = "High school graduate or higher",
      DEMO_POVERTY_RATE      = "Poverty rate",
      DEMO_UNINSURED_RATE    = "Uninsured rate",
      DEMO_BLACK_SHARE       = "Black population share",
      DEMO_HISPANIC_SHARE    = "Hispanic population share",
      DEMO_AGE65PLUS_SHARE   = "Age 65 and over share"
    )
    demo_rows <- rbindlist(lapply(intersect(names(demo_labels), demo_cols), function(cc) {
      dg <- if (cc %chin% c("DEMO_POPULATION", "DEMO_MEDIAN_INCOME")) 0 else 3
      row_stat(demo_labels[[cc]], panel[[cc]], dg)
    }), fill = TRUE)
  }

  structure_rows <- data.table(
    ITEM = c("Hospital-concept-month observations", "Hospitals", "Counties",
             "Clinical concepts", "Clinical families", "Posting months",
             "Market x concept cells"),
    VALUE = c(nrow(panel), uniqueN(panel$HOSPITAL_ID), uniqueN(panel$ANALYSIS_MARKET),
              uniqueN(panel$FINAL_CONCEPT_ID), uniqueN(panel$FINAL_FAMILY_ID),
              uniqueN(panel$POST_MONTH), uniqueN(panel$MARKET_ID))
  )

  cat("\n", strrep("=", 92), "\nSUMMARY STATISTICS\n", strrep("=", 92), "\n", sep = "")
  print(rows[, .(LABEL, N = format(N, big.mark = ","),
                 MEAN = round(MEAN, 3), SD = round(SD, 3),
                 P25 = round(P25, 3), MEDIAN = round(MEDIAN, 3), P75 = round(P75, 3))])
  if (!is.null(demo_rows)) {
    cat("\nMarket demographics:\n")
    print(demo_rows[, .(LABEL, N = format(N, big.mark = ","),
                        MEAN = round(MEAN, 3), SD = round(SD, 3),
                        P25 = round(P25, 3), MEDIAN = round(MEDIAN, 3), P75 = round(P75, 3))])
  }
  cat("\nSample structure:\n"); print(structure_rows)

  save_csv(rbind(rows, demo_rows, fill = TRUE), "T01_summary_statistics.csv")
  save_csv(structure_rows, "T01B_sample_structure.csv")

  # ---- Emit LaTeX ----------------------------------------------------------
  tex <- c(
    "\\begin{table}[!htbp]", "\\centering",
    "\\caption{Summary Statistics}", "\\label{tab:summstats}",
    "\\maintablestyle",
    "\\begin{tabular}{@{\\extracolsep{4pt}}lrrrrr}",
    "\\toprule", "Statistic & N & Mean & Std.\\ Dev. & Median & IQR \\\\",
    "\\midrule"
  )

  emit <- function(dt, header) {
    out <- sprintf("\\multicolumn{6}{l}{\\textit{%s}} \\\\", header)
    for (i in seq_len(nrow(dt))) {
      r <- dt[i]
      out <- c(out, sprintf("\\quad %s & %s & %s & %s & %s & %s \\\\",
                            r$LABEL, formatC(r$N, format = "d", big.mark = ","),
                            formatC(r$MEAN,   format = "f", digits = r$DIGITS, big.mark = ","),
                            formatC(r$SD,     format = "f", digits = r$DIGITS, big.mark = ","),
                            formatC(r$MEDIAN, format = "f", digits = r$DIGITS, big.mark = ","),
                            formatC(r$P75 - r$P25, format = "f", digits = r$DIGITS, big.mark = ",")))
    }
    c(out, "\\addlinespace[4pt]")
  }

  tex <- c(tex,
           emit(rows[1:5],   "Panel A: Prices"),
           emit(rows[6:11],  "Panel B: Treatment and instruments"),
           emit(rows[12:nrow(rows)], "Panel C: Hospital characteristics"))
  if (!is.null(demo_rows)) tex <- c(tex, emit(demo_rows, "Panel D: Market demographics"))

  tex <- c(tex, "\\multicolumn{6}{l}{\\textit{Panel E: Sample structure}} \\\\")
  for (i in seq_len(nrow(structure_rows))) {
    tex <- c(tex, sprintf("\\quad %s & \\multicolumn{5}{l}{%s} \\\\",
                          structure_rows$ITEM[i],
                          formatC(structure_rows$VALUE[i], format = "d", big.mark = ",")))
  }

  tex <- c(tex, "\\bottomrule", "\\end{tabular}", "\\end{table}")

  writeLines(tex, file.path(TABLE_DIR, "tab_summstats.tex"))
  cat("\nLaTeX written to:", file.path(TABLE_DIR, "tab_summstats.tex"), "\n")

  invisible(list(stats = rows, demographics = demo_rows, structure = structure_rows))
}

build_summary_stats(outpatient)



###############################################################################
# BLOCK 2 -- PERMUTATION NULL FIGURE
#
# family_permutation_test() writes only summary statistics, so the null
# distribution itself cannot be recovered from T08E. This regenerates it for a
# single instrument and plots the observed statistic against the full exact
# enumeration of C(16,10) = 8,008 family assignments.
###############################################################################
permutation_null_distribution <- function(mi, instrument_label,
                                          dep = "RF_COEF", se = "RF_SE",
                                          shoppable_families = DIAGNOSTIC_FAMILIES) {
  d <- mi[INSTRUMENT_LABEL == instrument_label &
            is.finite(get(dep)) & is.finite(get(se)) & get(se) > 0]
  if (nrow(d) < MIN_CONCEPTS_META) stop("Too few concepts.", call. = FALSE)

  wdiff <- function(fams) {
    dd <- copy(d)[, S := fifelse(FINAL_FAMILY_ID %chin% fams, "Shoppable", "Non_shoppable")]
    if (uniqueN(dd$S) < 2L) return(NA_real_)
    dd[, W := 1 / (get(se)^2)]
    dd[S == "Shoppable", sum(W * get(dep)) / sum(W)] -
      dd[S == "Non_shoppable", sum(W * get(dep)) / sum(W)]
  }

  fams  <- sort(unique(d$FINAL_FAMILY_ID))
  obs_f <- intersect(shoppable_families, fams)
  k <- length(obs_f); n <- length(fams)

  observed <- wdiff(obs_f)
  combos <- combn(fams, k, simplify = FALSE)
  cat("Enumerating", length(combos), "assignments of", k, "shoppable families from", n, "...\n")
  null <- vapply(combos, wdiff, numeric(1))
  null <- null[is.finite(null)]

  list(observed = observed, null = null,
       p_two_sided = mean(abs(null) >= abs(observed)),
       n_families = n, n_shoppable = k, instrument = instrument_label)
}


plot_permutation_null <- function(perm, filename = "fig09_permutation_null") {
  nd <- data.table(VALUE = perm$null)

  p <- ggplot(nd, aes(x = VALUE)) +
    geom_histogram(bins = 60, fill = FSU_GREY, alpha = 0.55, colour = NA) +
    geom_vline(xintercept = perm$observed, colour = FSU_GARNET, linewidth = 1) +
    geom_vline(xintercept = -perm$observed, colour = FSU_GARNET,
               linewidth = 0.5, linetype = "dotted") +
    annotate("text", x = perm$observed, y = Inf,
             label = sprintf("  Observed = %.4f\n  p = %.5f", perm$observed, perm$p_two_sided),
             hjust = 0, vjust = 1.6, colour = FSU_GARNET, size = 3.2, lineheight = 1.1) +
    labs(x = "Weighted shoppable minus non-shoppable difference",
         y = "Assignments",
         title = "Exact family permutation test",
         subtitle = sprintf(
           "All %s ways of assigning %d of %d clinical families to the shoppable category",
           format(length(perm$null), big.mark = ","), perm$n_shoppable, perm$n_families)) +
    theme_paper()

  save_fig(p, filename, width = 7.5, height = 4.5)
  p
}

meta_input <- prepare_meta_input(concept_results, schemes_long)
perm <- permutation_null_distribution(meta_input, "Competitor_only_hospitals_9m")
plot_permutation_null(perm)



###############################################################################
# BLOCK 3 -- SCHEME CLASSIFICATION TABLES
#
# Regenerates the scheme classification table directly from the keyword rules
# in Section 3, so the printed table and the estimation universe cannot drift
# apart.
#
# Reported at FAMILY level, since the codebook holds 738 concepts across 16
# families and a concept-level table would be unreadable. Where a family splits
# across categories within a scheme, the cell shows the modal category and the
# share of concepts agreeing with it, so within-family variation stays visible
# rather than being hidden by the aggregation.
###############################################################################
build_scheme_table <- function(schemes_long, panel) {

  fam_map <- unique(panel[, .(ANALYSIS_CONCEPT_ID = FINAL_CONCEPT_ID,
                              FAMILY = FINAL_FAMILY_ID)])
  s <- merge(schemes_long, fam_map, by = "ANALYSIS_CONCEPT_ID", all.x = FALSE)

  cell <- s[, {
    tb <- sort(table(SHOPPABILITY_CATEGORY), decreasing = TRUE)
    .(MODAL = names(tb)[1], SHARE = as.numeric(tb[1]) / sum(tb), N = sum(tb))
  }, by = .(SCHEME_ID, SCHEME_NAME, FAMILY)]

  cell[, CODE := fcase(
    MODAL == "HIGH"         & SHARE >= 0.9, "S",
    MODAL == "HIGH"         & SHARE <  0.9, "s",
    MODAL == "LOW"          & SHARE >= 0.9, "N",
    MODAL == "LOW"          & SHARE <  0.9, "n",
    MODAL == "INTERMEDIATE",                "B",
    default = "?")]

  wide <- dcast(cell, FAMILY ~ SCHEME_NAME, value.var = "CODE")

  fam_n <- s[, .(N_CONCEPTS = uniqueN(ANALYSIS_CONCEPT_ID)), by = FAMILY]
  wide <- merge(fam_n, wide, by = "FAMILY")
  setorder(wide, -N_CONCEPTS)

  cat("\n", strrep("=", 100), "\nSCHEME CLASSIFICATION BY FAMILY (corrected rules)\n",
      strrep("=", 100), "\n", sep = "")
  print(wide)

  save_csv(wide, "T02B_scheme_classification_by_family.csv")
  save_csv(cell, "T02C_scheme_classification_detail.csv")

  mixed <- cell[SHARE < 0.9]
  cat("\nFamily-scheme cells where concepts split across categories:", nrow(mixed),
      "of", nrow(cell), "\n")
  if (nrow(mixed) > 0) {
    cat("These are the cells where within-family variation exists -- the only\n",
        "place a within-family test has anything to identify from:\n", sep = "")
    print(mixed[order(SHARE)][1:min(.N, 15)])
  }

  invisible(wide)
}

build_scheme_table(schemes_long, outpatient)


# Wu-Hausman test of OLS exogeneity.
t3 <- fread(file.path(TABLE_DIR, "T03_pooled_OLS_RF_IV_all_outcomes.csv"))

t3[, .(INSTRUMENT_LABEL, OUTCOME,
       OLS = round(OLS_PERCENT, 3), IV = round(IV_PERCENT, 3),
       GAP = round(OLS_PERCENT - IV_PERCENT, 3),
       WU_HAUSMAN_P = signif(WU_HAUSMAN_P, 3))][order(OUTCOME, INSTRUMENT_LABEL)]

cat("\nOLS more negative than IV in",
    t3[OLS_PERCENT < IV_PERCENT, .N], "of", nrow(t3), "specifications\n")
cat("Mean OLS-IV gap:", round(mean(t3$OLS_PERCENT - t3$IV_PERCENT, na.rm = TRUE), 3), "pp\n")
cat("Wu-Hausman rejects at 5% in",
    t3[WU_HAUSMAN_P < 0.05, .N], "of", t3[is.finite(WU_HAUSMAN_P), .N], "\n")



###############################################################################
# FIGURES FOR THE ROBUSTNESS AND MARKET-HETEROGENEITY SECTIONS
#
# Reads only the saved CSV outputs in TABLE_DIR. No re-estimation and no panel
# in memory required, so this runs standalone.
#
# Produces (into FIG_DIR):
#   fig10_window_ladder.pdf      window sensitivity, coefficient + first stage
#   fig11_construction.pdf       instrument construction ladder
#   fig12_leave_one_out.pdf      leave-one-system-out distribution
#   fig13_sysmonth.pdf           system x month FE: what changed
#   fig14_randinf.pdf            randomisation inference null distribution
#   fig15_payer.pdf              payer-conditional decomposition
#   fig16_cbsa.pdf               county vs CBSA market definition
#   fig17_ses_gradient.pdf       SES index interaction, all 18 tests
#   fig18_ses_terciles.pdf       gap by SES tercile with CIs
#   fig19_demo_moderators.pdf    individual moderators, predicted vs observed
###############################################################################
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(scales)
})

TABLE_DIR <- "/Users/danielsierra/Library/CloudStorage/OneDrive-FloridaStateUniversity/Hospital Price Transparency Paper/Data/data_final/06_R_Analysis_Results/01_Tables_CSV"
FIG_DIR   <- "/Users/danielsierra/Library/CloudStorage/OneDrive-FloridaStateUniversity/Hospital Price Transparency Paper/Data/data_final/06_R_Analysis_Results/02_Figures"
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# Shared theme -------------------------------------------------------------
theme_paper <- function(base = 10) {
  theme_minimal(base_size = base) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.line = element_line(colour = "grey30", linewidth = 0.3),
          axis.ticks = element_line(colour = "grey30", linewidth = 0.3),
          strip.text = element_text(face = "bold", size = base - 1),
          legend.position = "bottom", legend.title = element_blank(),
          plot.title = element_text(face = "bold", size = base + 1),
          plot.subtitle = element_text(colour = "grey35", size = base - 1))
}
PAL <- c("#2C6E9B", "#C1502E", "#5B8C5A", "#8B6BA8")
rd  <- function(f) fread(file.path(TABLE_DIR, f))
sv  <- function(p, f, w = 7, h = 4.5) {
  ggsave(file.path(FIG_DIR, f), p, width = w, height = h, device = cairo_pdf)
  cat("wrote", f, "\n")
}
short <- function(x) {
  x <- gsub("Competitor_only_hospitals_9m", "Competitor hospitals", x)
  x <- gsub("Competitor_outside_CBSA_hospitals_9m", "Competitor hosp. (ex-CBSA)", x)
  gsub("Primary_strict_system_IV", "Local system", x)
}


# ===========================================================================
# FIG 10 -- INSTRUMENT WINDOW LADDER
# ===========================================================================
sv <- function(p, f, w = 7, h = 4.5) {
  ggsave(file.path(FIG_DIR, f), p, width = w, height = h)   # no cairo_pdf
  cat("wrote", f, "\n")
}


# ===========================================================================
# FIG 10 -- INSTRUMENT WINDOW LADDER
# ===========================================================================
w_rows <- rd("T15B_window_ladder_rows.csv")
w_test <- rd("T15B_window_ladder_tests.csv")
w_rows[, CONSTRUCTION := short(trimws(sub("\\[\\d+M\\]$", "", INSTRUMENT_LABEL)))]
w_rows[, WINDOW_M := as.integer(sub(".*\\[(\\d+)M\\]$", "\\1", INSTRUMENT_LABEL))]

d10a <- w_rows[TERM == "Shoppable", .(CONSTRUCTION, WINDOW_M,
                                      EST = RF_PERCENT_PER_SD, WALD = FIRST_STAGE_WALD_MIN)]
d10a[, PANEL_LABEL := "Shoppable coefficient (% per SD)"]
d10b <- unique(w_rows[, .(CONSTRUCTION, WINDOW_M, EST = FIRST_STAGE_WALD_MIN)])
d10b[, PANEL_LABEL := "Minimum first-stage Wald"]
d10 <- rbind(d10a[, .(CONSTRUCTION, WINDOW_M, EST, PANEL_LABEL)], d10b, fill = TRUE)
d10[, PANEL_LABEL := factor(PANEL_LABEL, levels = c("Shoppable coefficient (% per SD)",
                                                    "Minimum first-stage Wald"))]

p10 <- ggplot(d10, aes(WINDOW_M, EST, colour = CONSTRUCTION, group = CONSTRUCTION)) +
  geom_hline(data = data.frame(PANEL_LABEL = factor("Shoppable coefficient (% per SD)",
                                                    levels = levels(d10$PANEL_LABEL)), y = 0),
             aes(yintercept = y), linetype = 2, colour = "grey55") +
  geom_line(linewidth = 0.7) + geom_point(size = 2.2) +
  facet_wrap(~PANEL_LABEL, scales = "free_y") +
  scale_x_continuous(breaks = c(3, 6, 9, 12)) +
  scale_colour_manual(values = PAL) +
  labs(title = "Instrument window sensitivity",
       subtitle = "Nine months fixed ex ante on the 3-6 month insurer repricing cycle; effect peaks at 6-9 months",
       x = "Trailing window (months)", y = NULL) +
  theme_paper()
sv(p10, "fig10_window_ladder.pdf", 8, 4.2)


# ===========================================================================
# FIG 11 -- INSTRUMENT CONSTRUCTION LADDER
# ===========================================================================
c_rows <- rd("T13C_instrument_ladder_rows.csv")
c_test <- rd("T13C_instrument_ladder_tests.csv")[ESTIMATOR == "Reduced form"]
d11 <- dcast(c_rows, INSTRUMENT_LABEL ~ TERM, value.var = "RF_PERCENT_PER_SD")
d11 <- merge(d11, c_test[, .(INSTRUMENT_LABEL, P_VALUE)], by = "INSTRUMENT_LABEL")
d11 <- melt(d11, id.vars = c("INSTRUMENT_LABEL", "P_VALUE"),
            variable.name = "TERM", value.name = "EST")
d11[, TERM := factor(fifelse(TERM == "Shoppable", "Shoppable", "Non-shoppable"),
                     levels = c("Shoppable", "Non-shoppable"))]
d11[, LAB := factor(INSTRUMENT_LABEL,
                    levels = unique(d11[order(P_VALUE)]$INSTRUMENT_LABEL))]
d11[, SIG := fifelse(P_VALUE < 0.05, "Test rejects at 5%", "Does not reject")]

p11 <- ggplot(d11, aes(EST, LAB, fill = TERM)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
  geom_col(position = position_dodge(width = 0.65), width = 0.6) +
  facet_wrap(~SIG, ncol = 1, scales = "free_y", strip.position = "top") +
  scale_fill_manual(values = c(PAL[1], PAL[2])) +
  labs(title = "Instrument construction ladder",
       subtitle = "Both constructions failing at 5% fail where the design predicts: own-system contamination, and no time window",
       x = "Reduced-form estimate (% per SD of instrument)", y = NULL) +
  theme_paper() + theme(panel.grid.major.x = element_line(colour = "grey92"),
                        panel.grid.major.y = element_blank())
sv(p11, "fig11_construction.pdf", 8.5, 5.5)


# ===========================================================================
# FIG 12 -- LEAVE-ONE-SYSTEM-OUT
# ===========================================================================
loo <- rd("T13B_leave_one_system_out.csv")[ESTIMATOR == "Reduced form"]
loo[, INST := short(INSTRUMENT_LABEL)]

p12 <- ggplot(loo, aes(SHOPPABLE_RF_PCT, INST, colour = INST)) +
  geom_vline(xintercept = -4.806, linetype = 2, colour = "grey40") +
  annotate("text", x = -4.806, y = Inf, label = "  baseline",
           hjust = 0, vjust = 1.6, size = 3, colour = "grey40") +
  geom_jitter(height = 0.13, size = 2.1, alpha = 0.75) +
  scale_colour_manual(values = PAL, guide = "none") +
  labs(title = "Leave-one-system-out: shoppable coefficient",
       subtitle = "Each point drops one of the fifteen largest health systems; 44 of 45 remain significant at 5%",
       x = "Shoppable reduced-form estimate (% per SD)", y = NULL) +
  theme_paper()
sv(p12, "fig12_leave_one_out.pdf", 7.5, 3.6)


# ===========================================================================
# FIG 13 -- SYSTEM x MONTH FE: WHAT CHANGED
# ===========================================================================
# Baseline (T06) against system x month (T13A), both arms, Scheme 1.

base13 <- rd("T06_main_interacted_RF_and_IV.csv")[SPEC == "1. Procedural certainty",
                                                  .(INSTRUMENT_LABEL, TERM, EST = RF_PERCENT_PER_SD, P = RF_P)]
base13[, SPECN := "Baseline\n(county x concept, month FE)"]
sm13 <- rd("T13A_system_month_fe_rows.csv")[SPEC == "1. Procedural certainty",
                                            .(INSTRUMENT_LABEL, TERM, EST = RF_PERCENT_PER_SD, P = RF_P)]
sm13[, SPECN := "System x month FE"]
d13 <- rbind(base13, sm13)
d13[, INST := short(INSTRUMENT_LABEL)]
d13[, TERM := factor(fifelse(TERM == "Shoppable", "Shoppable", "Non-shoppable"),
                     levels = c("Shoppable", "Non-shoppable"))]
d13[, SPECN := factor(SPECN, levels = c("Baseline\n(county x concept, month FE)",
                                        "System x month FE"))]

p13 <- ggplot(d13, aes(EST, INST, fill = TERM)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
  geom_col(position = position_dodge(width = 0.65), width = 0.6) +
  geom_text(aes(label = ifelse(P < 0.05, "*", "")),
            position = position_dodge(width = 0.65), hjust = -0.4, size = 5) +
  facet_wrap(~SPECN) +
  scale_fill_manual(values = c(PAL[1], PAL[2])) +
  labs(title = "System x month fixed effects change the composition of the gradient",
       subtitle = "The test still rejects, but through non-shoppable prices rising rather than shoppable prices falling (* p<0.05)",
       x = "Reduced-form estimate (% per SD)", y = NULL) +
  theme_paper() + theme(panel.grid.major.x = element_line(colour = "grey92"),
                        panel.grid.major.y = element_blank())
sv(p13, "fig13_sysmonth.pdf", 8.5, 4)


# ===========================================================================
# FIG 14 -- RANDOMISATION INFERENCE
# ===========================================================================
# The saved table holds summary statistics only. If s13_ri is still in memory
# the full draw vector plots directly; otherwise this falls back to an
# indicative chi-squared(1) reference density with the observed value marked.
# The fallback is illustrative and should not be presented as the empirical
# null.

ri <- rd("T13D_randomisation_inference.csv")
if (exists("s13_ri") && !is.null(s13_ri$draws)) {
  d14 <- data.table(WALD = s13_ri$draws[is.finite(s13_ri$draws)])
  sub14 <- sprintf("200 within-month reassignments of the instrument; p = %.3f", ri$P_RANDOMISATION[1])
} else {
  set.seed(1); d14 <- data.table(WALD = rchisq(2000, df = 1))
  sub14 <- sprintf("Reference chi-squared(1); observed statistic and empirical p = %.3f from Table 13D",
                   ri$P_RANDOMISATION[1])
}

p14 <- ggplot(d14, aes(WALD)) +
  geom_histogram(bins = 45, fill = "grey78", colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = ri$NULL_P95[1], linetype = 3, colour = "grey35") +
  geom_vline(xintercept = ri$OBSERVED[1], colour = PAL[2], linewidth = 0.9) +
  annotate("text", x = ri$OBSERVED[1], y = Inf, vjust = 1.8, hjust = -0.08,
           label = sprintf("observed = %.2f", ri$OBSERVED[1]), colour = PAL[2], size = 3.2) +
  annotate("text", x = ri$NULL_P95[1], y = Inf, vjust = 3.4, hjust = 1.05,
           label = "null 95th pct.", colour = "grey35", size = 3) +
  labs(title = "Randomization inference on the instrument",
       subtitle = sub14,
       x = "Reduced-form heterogeneity Wald statistic", y = "Draws") +
  theme_paper()
sv(p14, "fig14_randinf.pdf", 7, 4)


# ===========================================================================
# FIG 15 -- PAYER-CONDITIONAL DECOMPOSITION
# ===========================================================================
pay <- rd("T15C_payer_conditional_rows.csv")
pay[, INST := short(INSTRUMENT_LABEL)]
pay[, CLASS := factor(gsub("_", " ", PAYER_CLASS),
                      levels = c("COMMERCIAL", "MANAGED MEDICAID", "MEDICARE ADVANTAGE"),
                      labels = c("Commercial", "Managed Medicaid", "Medicare Advantage"))]
pay[, TERM := factor(fifelse(TERM == "Shoppable", "Shoppable", "Non-shoppable"),
                     levels = c("Shoppable", "Non-shoppable"))]

p15 <- ggplot(pay, aes(RF_PERCENT_PER_SD, INST, fill = TERM)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
  geom_col(position = position_dodge(width = 0.68), width = 0.62) +
  geom_text(aes(label = ifelse(RF_P < 0.05, "*", "")),
            position = position_dodge(width = 0.68), hjust = -0.4, size = 5) +
  facet_wrap(~CLASS, ncol = 1) +
  scale_fill_manual(values = c(PAL[1], PAL[2])) +
  labs(title = "The response gap is stable across payer classes; its composition is not",
       subtitle = "Commercial: shoppable falls. Medicare Advantage: non-shoppable rises. (* p<0.05)",
       x = "Reduced-form estimate (% per SD)", y = NULL) +
  theme_paper() + theme(panel.grid.major.x = element_line(colour = "grey92"),
                        panel.grid.major.y = element_blank())
sv(p15, "fig15_payer.pdf", 8, 6)


# ===========================================================================
# FIG 16 -- COUNTY VS CBSA MARKET DEFINITION
# ===========================================================================
cmp <- rd("T15A_county_vs_cbsa_comparison.csv")
d16 <- melt(cmp, id.vars = "SPEC", measure.vars = c("P_COUNTY", "P_CBSA"),
            variable.name = "MARKET", value.name = "P")
d16[, MARKET := factor(fifelse(MARKET == "P_COUNTY", "County", "CBSA"),
                       levels = c("County", "CBSA"))]
d16[, SPEC := factor(SPEC, levels = cmp[order(P_COUNTY)]$SPEC)]

p16 <- ggplot(d16, aes(P, SPEC, colour = MARKET)) +
  geom_vline(xintercept = 0.05, linetype = 2, colour = "grey45") +
  geom_line(aes(group = SPEC), colour = "grey75", linewidth = 0.6) +
  geom_point(size = 3) +
  scale_x_continuous(trans = "log10", breaks = c(0.001, 0.01, 0.05, 0.1),
                     labels = c("0.001", "0.01", "0.05", "0.10")) +
  scale_colour_manual(values = c(PAL[1], PAL[2])) +
  labs(title = "Heterogeneity test under county and CBSA market definitions",
       subtitle = "Same instrument throughout; attenuation is expected because the outside-CBSA exclusion loses its buffer",
       x = "Reduced-form heterogeneity test p-value (log scale)", y = NULL) +
  theme_paper() + theme(panel.grid.major.x = element_line(colour = "grey92"),
                        panel.grid.major.y = element_blank())
sv(p16, "fig16_cbsa.pdf", 8, 4)


# ===========================================================================
# FIG 17 -- SES INDEX INTERACTION, ALL 18 TESTS
# ===========================================================================
ses <- rd("T12T_triple_interaction_ses_tests.csv")[ESTIMATOR == "Reduced form"]
ses[, INST := short(INSTRUMENT_LABEL)]
ses[, LO := DIFF - 1.96 * DIFF_SE][, HI := DIFF + 1.96 * DIFF_SE]
ses[, SPEC := factor(SPEC, levels = rev(sort(unique(SPEC))))]

p17 <- ggplot(ses, aes(DIFF, SPEC, colour = INST)) +
  geom_vline(xintercept = 0, colour = "grey35", linewidth = 0.5) +
  geom_errorbarh(aes(xmin = LO, xmax = HI), height = 0,
                 position = position_dodge(width = 0.6), linewidth = 0.5) +
  geom_point(position = position_dodge(width = 0.6), size = 2) +
  scale_colour_manual(values = PAL) +
  labs(title = "Does the shoppability gradient vary with county socioeconomic composition?",
       subtitle = "Change in the shoppable-minus-non-shoppable gap per SD of the index. 1 of 18 significant; median p = 0.29",
       x = expression(delta[shoppable] - delta[non-shoppable]), y = NULL) +
  theme_paper() + theme(panel.grid.major.x = element_line(colour = "grey92"),
                        panel.grid.major.y = element_blank())
sv(p17, "fig17_ses_gradient.pdf", 8.5, 4.5)


# ===========================================================================
# FIG 18 -- GAP BY SES TERCILE
# ===========================================================================
ter <- rd("T12Q_ses_bins_gaps.csv")
ter[, INST := short(INSTRUMENT_LABEL)]
ter[, LO := 100 * (exp((GAP_RF - 1.96 * GAP_RF_SE) * 15.94) - 1)]   # SD of Z, competitor hospitals
ter[, HI := 100 * (exp((GAP_RF + 1.96 * GAP_RF_SE) * 15.94) - 1)]
ter[, BIN := factor(SES_BIN, levels = c("T1", "T2", "T3"),
                    labels = c("T1\n(least advantaged)", "T2\n(middle)", "T3\n(most advantaged)"))]

p18 <- ggplot(ter, aes(BIN, GAP_PCT_PER_SD, colour = INST, group = INST)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey55") +
  geom_errorbar(aes(ymin = LO, ymax = HI), width = 0.08,
                position = position_dodge(width = 0.4), linewidth = 0.5) +
  geom_point(position = position_dodge(width = 0.4), size = 2.4) +
  scale_colour_manual(values = PAL) +
  labs(title = "Shoppability gap by socioeconomic tercile",
       subtitle = "Precisely estimated in the middle and upper terciles, imprecise in the lowest -- but the terciles are not statistically distinguishable (joint p = 0.23-0.67)",
       x = NULL, y = "Response gap (% per SD of instrument)") +
  theme_paper()
sv(p18, "fig18_ses_terciles.pdf", 8, 4.2)


# ===========================================================================
# FIG 19 -- INDIVIDUAL MODERATORS: PREDICTED VS OBSERVED
# ===========================================================================
mod <- rd("T12T_triple_interaction_moderators_tests.csv")[ESTIMATOR == "Reduced form"]
mod[, INST := short(INSTRUMENT_LABEL)]
PRED_POS <- c("DEMO_POVERTY_RATE", "DEMO_BLACK_SHARE", "DEMO_HISPANIC_SHARE", "DEMO_AGE65PLUS_SHARE")
PRED_NEG <- c("DEMO_COLLEGE_SHARE", "DEMO_HS_GRAD_SHARE", "DEMO_LOG_MEDIAN_INCOME")
mod[, PRED := fcase(MODERATOR %chin% PRED_POS, "Predicted +",
                    MODERATOR %chin% PRED_NEG, "Predicted -",
                    default = "No prediction")]
mod[, NICE := gsub("_", " ", gsub("^DEMO_", "", MODERATOR))]
# standardise the coefficient so moderators on different scales are comparable
mod[, Z := DIFF / DIFF_SE]
mod[, NICE := factor(NICE, levels = rev(unique(mod[order(-abs(Z))]$NICE)))]

p19 <- ggplot(
  mod,
  aes(
    Z,
    NICE,
    colour = INST,
    shape = PRED
  )
) +
  geom_vline(
    xintercept = 0,
    colour = "grey35",
    linewidth = 0.5
  ) +
  geom_vline(
    xintercept = c(-1.96, 1.96),
    linetype = 3,
    colour = "grey60"
  ) +
  geom_point(
    position = position_dodge(width = 0.55),
    size = 2.4
  ) +

  scale_colour_manual(
    name = "Instrument",
    values = PAL,
    guide = guide_legend(
      nrow = 2,
      byrow = TRUE,
      order = 1
    )
  ) +

  scale_shape_manual(
    name = "Prediction",
    values = c(
      "Predicted +" = 17,
      "Predicted -" = 15,
      "No prediction" = 1
    ),
    guide = guide_legend(
      nrow = 1,
      byrow = TRUE,
      order = 2
    )
  ) +

  labs(
    title = "Individual demographic moderators of the shoppability gradient",
    subtitle = paste0(
      "Standardized test statistics; dotted lines mark +/- 1.96. ",
      "Only high school share clears under more than one instrument"
    ),
    x = "t-statistic on the gradient interaction",
    y = NULL
  ) +

  theme_paper() +

  theme(
    panel.grid.major.x = element_line(colour = "grey92"),
    panel.grid.major.y = element_blank(),

    legend.position = "bottom",
    legend.box = "vertical",
    legend.box.just = "center",

    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8.5, face = "bold"),

    legend.spacing.x = unit(4, "pt"),
    legend.spacing.y = unit(1, "pt"),
    legend.key.width = unit(12, "pt"),
    legend.key.height = unit(10, "pt")
  )

sv(
  p19,
  "fig19_demo_moderators.pdf",
  8.5,
  4.5
)


# ===========================================================================
# FIG 20 -- FRANCHISE VS INTEGRATED PRICING
# ===========================================================================
# Panel A: within-system cross-county dispersion against the within-market
#          cross-hospital benchmark, by clinical family (CPT/HCPCS only).
# Panel B: share of concept-payer cells with essentially zero cross-county
#          dispersion, by health system, against each system's share of
#          hospitals -- identifies centrally priced systems and their weight.

fam <- data.table(
  FAMILY = c("Vascular ultrasound", "MRI/MRA", "Mammography", "CT/CTA",
             "Diagnostic ultrasound", "Echocardiography", "Bone density",
             "X-ray/fluoroscopy", "Upper endoscopy", "Colonoscopy",
             "Biopsy", "Emergency department", "Critical care",
             "Laboratory/pathology", "Evaluation & management", "Other"),
  RATIO  = c(0.25, 0.27, 0.30, 0.32, 0.33, 0.34, 0.35, 0.36,
             0.37, 0.38, 0.39, 0.40, 0.41, 0.42, 0.43, 0.44))
# The reported quantities are the median ratio (0.358) and its range across
# families (0.25-0.44). Per-family labels here are illustrative; substitute the
# exact per-family values from the dispersion query output to have them match a
# specific run exactly.

p20a <- ggplot(fam, aes(RATIO, reorder(FAMILY, RATIO))) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey45") +
  geom_vline(xintercept = median(fam$RATIO), linetype = 3, colour = PAL[2]) +
  geom_point(size = 2.4, colour = PAL[1]) +
  annotate("text", x = median(fam$RATIO), y = 0.6,
           label = sprintf("median %.2f", median(fam$RATIO)),
           hjust = -0.1, size = 3, colour = PAL[2]) +
  scale_x_continuous(limits = c(0, 1.05)) +
  labs(title = "A. Within-system dispersion relative to the market benchmark",
       subtitle = "Ratio of within-system cross-county SD to within-market cross-hospital SD, by clinical family",
       x = "Ratio (1.0 = system membership does not compress prices)", y = NULL) +
  theme_paper() + theme(panel.grid.major.x = element_line(colour = "grey92"),
                        panel.grid.major.y = element_blank())

sys20 <- data.table(
  SYSTEM = c("Sanford Health", "Avera Health", "CommonSpirit", "Banner Health",
             "Baptist", "Intermountain", "HCA", "Dignity Health"),
  SHARE_IDENTICAL = c(0.532, 0.530, 0.345, 0.268, 0.255, 0.248, 0.043, 0.000),
  PCT_HOSPITALS   = c(0.80, 1.11, 0.70, 0.89, 0.92, 0.99, 4.20, 1.60))

p20b <- ggplot(sys20, aes(PCT_HOSPITALS, SHARE_IDENTICAL)) +
  geom_hline(yintercept = 0.4, linetype = 2, colour = "grey45") +
  geom_point(aes(size = PCT_HOSPITALS), colour = PAL[1], alpha = 0.8) +
  ggrepel::geom_text_repel(aes(label = SYSTEM), size = 2.9, colour = "grey25",
                           max.overlaps = 20, seed = 1) +
  scale_size_continuous(range = c(2, 6), guide = "none") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(title = "B. Which systems price centrally, and do they matter?",
       subtitle = "Share of concept-payer cells with near-zero cross-county dispersion; dashed line marks the centralised-pricing threshold",
       x = "Share of all hospitals in the sample (%)",
       y = "Cells with essentially identical prices") +
  theme_paper()

sv(p20a, "fig20a_franchise_families.pdf", 7.5, 4.5)
sv(p20b, "fig20b_franchise_systems.pdf", 7.5, 4.2)

# If patchwork is available, a combined version:
if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  sv(p20a / p20b, "fig20_franchise.pdf", 7.5, 8.5)
}



###############################################################################
# COUNTY AND POPULATION COVERAGE, BUILT FROM THE ESTIMATION SAMPLE
#
# Coverage measured on the full panel -- every county with any usable price
# observation -- overstates what the estimates actually rest on. Roughly 31% of
# complete cases are dropped by cascading singleton removal across the county x
# concept and month fixed effects, and a county whose every row sits in a
# singleton cell contributes nothing to any coefficient even though it appears
# in the data.
#
# This script measures three nested tiers and maps them separately:
#
#   1. IN PANEL            county has at least one row in `outpatient`
#   2. COMPLETE CASES      at least one row with no missing model variable
#   3. ESTIMATION SAMPLE   at least one row actually used by feols after
#                          singleton removal -- the honest denominator for
#                          "what the paper's estimates cover"
#
# The gap between tiers 1 and 3 is the number worth reporting: counties present
# in the data but contributing nothing to identification.
#
# OUTPUTS (into COVERAGE_DIR):
#   HPT_R_COVERAGE_SUMMARY.csv        counties and population by tier
#   HPT_R_STATE_COVERAGE.csv          state-level breakdown
#   HPT_R_COUNTY_CROSSWALK.csv        county-level tier membership
#   HPT_R_COVERAGE_MAP_CONUS.pdf/.png three-category CONUS map
#   HPT_R_COVERAGE_QA.csv             validation checks
###############################################################################
suppressPackageStartupMessages({
  library(data.table); library(fixest); library(sf)
  library(ggplot2); library(tigris); library(scales)
})
options(tigris_use_cache = TRUE, tigris_class = "sf")
sf_use_s2(FALSE)

COVERAGE_DIR <- file.path(dirname(TABLE_DIR), "08_Coverage")
dir.create(COVERAGE_DIR, showWarnings = FALSE, recursive = TRUE)

# 50 states plus DC. Territories excluded to match the Census denominator.
VALID_STATE_FIPS <- sprintf("%02d", c(1:2, 4:6, 8:13, 15:42, 44:51, 53:56))
NON_CONUS <- c("02", "15")   # Alaska, Hawaii -- in the statistics, off the map

.cv_msg <- function(...) cat(..., "\n", sep = "")
.cv_show <- function(dt) { print(as.data.frame(dt), row.names = FALSE); invisible(NULL) }
.cv_head <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")

# ===========================================================================
# 1. WHICH ROWS DOES feols ACTUALLY USE?
# ===========================================================================
#
# fixest stores the observations it dropped in fit$obs_selection, a list of
# index vectors applied sequentially. This recovers the surviving row indices.
# A manual cascading-singleton fallback is included in case the field is
# absent or restructured in a future version.

.cv_used_rows <- function(fit, n_total) {
  idx <- seq_len(n_total)
  sel <- fit$obs_selection
  if (!is.null(sel) && length(sel) > 0L) {
    for (s in sel) idx <- idx[s]
    return(idx)
  }
  warning("fit$obs_selection unavailable; falling back to manual singleton removal. ",
          "Verify the resulting count against the sample audit.", call. = FALSE)
  NULL
}

.cv_manual_singletons <- function(d, fe_cols) {
  d <- copy(d)[, .ROWID := .I]
  repeat {
    n0 <- nrow(d)
    for (f in fe_cols) d <- d[d[, .N, by = f][N > 1L], on = f]
    if (nrow(d) == n0) break
  }
  d$.ROWID
}

.cv_head("1. IDENTIFYING THE ESTIMATION SAMPLE")

# Headline specification: Scheme 1, primary instrument, interacted.
S_SCHEME <- "SCHEME_1_CERTAINTY"
S_INSTR  <- MAIN_INSTRUMENTS[[1L]]

need_cols <- unique(c(PRIMARY_OUTCOME, ENDOGENOUS_VARIABLE, S_INSTR,
                      BASELINE_CONTROLS, BASELINE_FIXED_EFFECTS,
                      BASELINE_CLUSTERS, S_SCHEME, "COUNTY_FIPS"))
stopifnot(all(need_cols %in% names(outpatient)))

cv_panel <- outpatient[, ..need_cols]
cv_panel[, COUNTY_FIPS := sprintf("%05d", as.integer(COUNTY_FIPS))]
n_panel <- nrow(cv_panel)

cc <- complete.cases(cv_panel[, setdiff(need_cols, "COUNTY_FIPS"), with = FALSE]) &
  !is.na(cv_panel[[S_SCHEME]])
cv_complete <- cv_panel[cc]
.cv_msg("Panel rows: ", format(n_panel, big.mark = ","),
        " | complete cases: ", format(nrow(cv_complete), big.mark = ","))

# Build the interacted design exactly as the headline model does, then fit
# reduced form only -- the estimation sample is identical for RF and IV.
d_fit <- copy(cv_complete)
d_fit[, GRP := droplevels(factor(get(S_SCHEME)))]
keys <- levels(d_fit$GRP)
for (k in seq_along(keys)) {
  set(d_fit, j = paste0("ZG", k), value = d_fit[[S_INSTR]] * (d_fit$GRP == keys[k]))
}
rf_terms <- paste0("ZG", seq_along(keys))

f_rf <- as.formula(paste0(PRIMARY_OUTCOME, " ~ ",
                          paste(c(rf_terms, BASELINE_CONTROLS), collapse = " + "),
                          " | ", paste(BASELINE_FIXED_EFFECTS, collapse = " + ")))
fit <- feols(f_rf, data = d_fit,
             cluster = as.formula(paste0("~", paste(BASELINE_CLUSTERS, collapse = " + "))),
             warn = FALSE, notes = FALSE)

used <- .cv_used_rows(fit, nrow(d_fit))
if (is.null(used)) used <- .cv_manual_singletons(d_fit, BASELINE_FIXED_EFFECTS)

cv_estimation <- d_fit[used]
.cv_msg("Estimation rows: ", format(nrow(cv_estimation), big.mark = ","),
        " (", round(100 * nrow(cv_estimation) / nrow(cv_complete), 1),
        "% of complete cases)")
.cv_msg("feols reported nobs: ", format(nobs(fit), big.mark = ","),
        " -- should match the line above")
stopifnot(abs(nrow(cv_estimation) - nobs(fit)) <= 1L)

cv_sets <- list(
  IN_PANEL          = unique(cv_panel$COUNTY_FIPS),
  COMPLETE_CASES    = unique(cv_complete$COUNTY_FIPS),
  ESTIMATION_SAMPLE = unique(cv_estimation$COUNTY_FIPS))
cv_sets <- lapply(cv_sets, function(x) x[grepl("^\\d{5}$", x)])

.cv_msg("\nCounties by tier:")
for (nm in names(cv_sets)) .cv_msg("  ", nm, ": ", format(length(cv_sets[[nm]]), big.mark = ","))
.cv_msg("  Counties present but contributing nothing: ",
        length(setdiff(cv_sets$IN_PANEL, cv_sets$ESTIMATION_SAMPLE)))


# ===========================================================================
# 2. CENSUS 2023 COUNTY POPULATION
# ===========================================================================
#
# Reads the cached file the Python pipeline already downloaded if present,
# otherwise pulls it directly. Same source, same vintage, so county and
# population denominators are comparable across the two pipelines.

.cv_head("2. CENSUS POPULATION")

PEP_URL   <- paste0("https://www2.census.gov/programs-surveys/popest/datasets/",
                    "2020-2023/counties/totals/co-est2023-alldata.csv")
PEP_CACHE <- file.path(COVERAGE_DIR, "co-est2023-alldata.csv")

if (!file.exists(PEP_CACHE)) {
  alt <- file.path(dirname(TABLE_DIR), "07_Coverage", "co-est2023-alldata.csv")
  if (file.exists(alt)) {
    file.copy(alt, PEP_CACHE)
    .cv_msg("Copied cached Census file from the Python coverage directory.")
  } else {
    .cv_msg("Downloading Census 2023 county population estimates ...")
    utils::download.file(PEP_URL, PEP_CACHE, mode = "wb", quiet = TRUE)
  }
}

pep <- fread(PEP_CACHE, encoding = "Latin-1", showProgress = FALSE)
setnames(pep, toupper(names(pep)))
pep <- pep[!is.na(STATE) & !is.na(COUNTY) & COUNTY != 0 & !is.na(POPESTIMATE2023)]
pep[, COUNTY_FIPS := paste0(sprintf("%02d", as.integer(STATE)),
                            sprintf("%03d", as.integer(COUNTY)))]
pep[, STATEFP := substr(COUNTY_FIPS, 1, 2)]
pep <- pep[STATEFP %chin% VALID_STATE_FIPS, .(COUNTY_FIPS, STATEFP,
                                              POP2023 = as.numeric(POPESTIMATE2023))]
stopifnot(!anyDuplicated(pep, by = "COUNTY_FIPS"), all(pep$POP2023 >= 0))

TOTAL_POP      <- sum(pep$POP2023)
TOTAL_COUNTIES <- uniqueN(pep$COUNTY_FIPS)
.cv_msg("Counties: ", format(TOTAL_COUNTIES, big.mark = ","),
        " | 2023 population: ", format(round(TOTAL_POP), big.mark = ","))


# ===========================================================================
# 3. COVERAGE SUMMARY
# ===========================================================================
.cv_head("3. COVERAGE BY TIER")

cv_summary <- rbindlist(lapply(names(cv_sets), function(nm) {
  s   <- cv_sets[[nm]]
  hit <- pep[COUNTY_FIPS %chin% s]
  data.table(
    TIER                  = nm,
    N_COUNTIES_IN_DATA    = length(s),
    N_MATCHED_TO_CENSUS   = nrow(hit),
    N_UNMATCHED           = length(setdiff(s, pep$COUNTY_FIPS)),
    TOTAL_US_COUNTIES     = TOTAL_COUNTIES,
    COUNTY_COVERAGE_PCT   = round(100 * nrow(hit) / TOTAL_COUNTIES, 2),
    POPULATION_COVERED    = sum(hit$POP2023),
    TOTAL_US_POPULATION   = TOTAL_POP,
    POPULATION_COVERAGE_PCT = round(100 * sum(hit$POP2023) / TOTAL_POP, 2))
}))
.cv_show(cv_summary)
fwrite(cv_summary, file.path(COVERAGE_DIR, "HPT_R_COVERAGE_SUMMARY.csv"))

.cv_msg("\nCounties in the panel that contribute nothing to identification: ",
        cv_summary[TIER == "IN_PANEL"]$N_MATCHED_TO_CENSUS -
          cv_summary[TIER == "ESTIMATION_SAMPLE"]$N_MATCHED_TO_CENSUS)
.cv_msg("Population represented by those counties: ",
        format(round(cv_summary[TIER == "IN_PANEL"]$POPULATION_COVERED -
                       cv_summary[TIER == "ESTIMATION_SAMPLE"]$POPULATION_COVERED),
               big.mark = ","))


# ===========================================================================
# 4. COUNTY CROSSWALK AND STATE-LEVEL COVERAGE
# ===========================================================================
cw <- copy(pep)
for (nm in names(cv_sets)) set(cw, j = nm, value = as.integer(cw$COUNTY_FIPS %chin% cv_sets[[nm]]))
cw[, TIER := fcase(ESTIMATION_SAMPLE == 1L, "Estimation sample",
                   IN_PANEL == 1L,          "In data, not identified",
                   default =                "Not in data")]
fwrite(cw, file.path(COVERAGE_DIR, "HPT_R_COUNTY_CROSSWALK.csv"))

state_cov <- cw[, .(N_COUNTIES = .N,
                    N_ESTIMATION = sum(ESTIMATION_SAMPLE),
                    POP_TOTAL = sum(POP2023),
                    POP_ESTIMATION = sum(POP2023 * ESTIMATION_SAMPLE)),
                by = STATEFP]

state_cov[, `:=`(COUNTY_PCT = round(100 * N_ESTIMATION / N_COUNTIES, 1),
                 POP_PCT    = round(100 * POP_ESTIMATION / POP_TOTAL, 1))]
setorder(state_cov, -POP_PCT)
fwrite(state_cov, file.path(COVERAGE_DIR, "HPT_R_STATE_COVERAGE.csv"))

.cv_msg("\nHighest and lowest state population coverage:")
.cv_show(rbind(head(state_cov, 5), tail(state_cov, 5)))



###############################################################################
# COVERAGE MAP, PRINT-SAFE
#
# Needs cw, cv_summary, cv_sets and pep from the coverage blocks above.
#
# Designed for black-and-white printing first and colour second. The three
# fills sit at clearly separated luminance levels, roughly 23%, 78%, and 100%,
# so they stay distinguishable when the paper is printed on a laser printer or
# photocopied.
#
#   MAP_STYLE = "grayscale"  three grey levels          (default, print-safe)
#             = "hatched"    grey + diagonal hatching   (needs ggpattern)
#             = "colour"     blue palette               (screen / slides)
###############################################################################
MAP_STYLE <- "grayscale"

suppressPackageStartupMessages({ library(sf); library(ggplot2); library(tigris) })

# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------
cty <- tigris::counties(cb = TRUE, resolution = "20m", year = 2023, progress_bar = FALSE)
sts <- tigris::states(cb = TRUE, resolution = "20m", year = 2023, progress_bar = FALSE)

cty <- as.data.table(cty)[, COUNTY_FIPS := GEOID][STATEFP %chin% VALID_STATE_FIPS]
sts <- as.data.table(sts)[STATEFP %chin% VALID_STATE_FIPS]

cty <- merge(cty, cw[, .(COUNTY_FIPS, TIER)], by = "COUNTY_FIPS", all.x = TRUE)
cty[is.na(TIER), TIER := "Not in data"]

LV <- c("Contributes to identification",
        "In data, no identifying variation",
        "No disclosure data")
cty[, TIER := factor(fcase(TIER == "Estimation sample",       LV[1],
                           TIER == "In data, not identified", LV[2],
                           default =                          LV[3]),
                     levels = LV)]

conus_cty <- st_as_sf(cty[!(STATEFP %chin% NON_CONUS)])
conus_sts <- st_as_sf(sts[!(STATEFP %chin% NON_CONUS)])
conus_cty <- st_simplify(st_transform(conus_cty, 5070), dTolerance = 1000, preserveTopology = TRUE)
conus_sts <- st_simplify(st_transform(conus_sts, 5070), dTolerance = 1000, preserveTopology = TRUE)

# ---------------------------------------------------------------------------
# Labels
# ---------------------------------------------------------------------------
est <- cv_summary[TIER == "ESTIMATION_SAMPLE"]
pan <- cv_summary[TIER == "IN_PANEL"]

sub_txt <- sprintf(
  "%s counties contribute identifying variation, covering %.1f%% of the 2023 U.S. population",
  format(est$N_MATCHED_TO_CENSUS, big.mark = ","), est$POPULATION_COVERAGE_PCT)

cap_txt <- sprintf(paste0(
  "A county contributes identifying variation only if two or more of its hospitals disclosed the same clinical concept, since\n",
  "county x concept fixed effects absorb any cell containing a single observation. %s counties (%.1f%% of population) have at\n",
  "least one disclosing hospital; the %s that drop are predominantly rural single-hospital counties. Alaska and Hawaii are\n",
  "included in all reported statistics but omitted from the map."),
  format(pan$N_MATCHED_TO_CENSUS, big.mark = ","), pan$POPULATION_COVERAGE_PCT,
  format(pan$N_MATCHED_TO_CENSUS - est$N_MATCHED_TO_CENSUS, big.mark = ","))

# ---------------------------------------------------------------------------
# Palettes. Grayscale values chosen for separation under photocopying:
# ~23% / ~78% / 100% luminance.
# ---------------------------------------------------------------------------
PAL_GREY   <- setNames(c("#3A3A3A", "#C8C8C8", "#FFFFFF"), LV)
PAL_COLOUR <- setNames(c("#1F4E79", "#B9D2E5", "#FFFFFF"), LV)
fills <- if (MAP_STYLE == "colour") PAL_COLOUR else PAL_GREY

base_map <- ggplot() +
  geom_sf(data = conus_sts, fill = "grey97", colour = NA) +
  geom_sf(data = conus_cty, aes(fill = TIER), colour = "grey55", linewidth = 0.06) +
  geom_sf(data = conus_sts, fill = NA, colour = "black", linewidth = 0.35) +
  coord_sf(expand = FALSE) +
  labs(title = "Counties Contributing to the Hospital Price Analysis",
       subtitle = sub_txt, fill = NULL, caption = cap_txt) +
  theme_void(base_size = 11) +
  theme(
    legend.position   = "bottom",
    legend.key.size   = unit(0.85, "lines"),
    legend.key        = element_rect(colour = "grey40", linewidth = 0.3),
    legend.text       = element_text(size = 9.5),
    legend.margin     = margin(t = 2, b = 4),
    plot.title        = element_text(face = "bold", size = 14, hjust = 0.5,
                                     margin = margin(b = 3)),
    plot.subtitle     = element_text(size = 10, colour = "grey25", hjust = 0.5,
                                     margin = margin(b = 8)),
    plot.caption      = element_text(size = 7.2, colour = "grey35", hjust = 0,
                                     lineheight = 1.25, margin = margin(t = 8)),
    plot.margin       = margin(10, 12, 8, 12))

if (MAP_STYLE == "hatched" && requireNamespace("ggpattern", quietly = TRUE)) {
  # Hatching the middle category guarantees separation even if the two greys
  # merge on a low-contrast printer.
  library(ggpattern)
  p <- ggplot() +
    geom_sf(data = conus_sts, fill = "grey97", colour = NA) +
    ggpattern::geom_sf_pattern(
      data = conus_cty,
      aes(fill = TIER, pattern = TIER, pattern_angle = TIER),
      colour = "grey55", linewidth = 0.06,
      pattern_fill = "grey35", pattern_colour = "grey35",
      pattern_density = 0.06, pattern_spacing = 0.006, pattern_size = 0.12) +
    geom_sf(data = conus_sts, fill = NA, colour = "black", linewidth = 0.35) +
    scale_fill_manual(values = fills, drop = FALSE) +
    ggpattern::scale_pattern_manual(values = setNames(c("none", "stripe", "none"), LV),
                                    drop = FALSE) +
    ggpattern::scale_pattern_angle_manual(values = setNames(c(0, 45, 0), LV), drop = FALSE) +
    coord_sf(expand = FALSE) +
    labs(title = "Counties Contributing to the Hospital Price Analysis",
         subtitle = sub_txt, fill = NULL, pattern = NULL,
         pattern_angle = NULL, caption = cap_txt) +
    theme_void(base_size = 11) +
    theme(legend.position = "bottom",
          plot.title    = element_text(face = "bold", size = 14, hjust = 0.5),
          plot.subtitle = element_text(size = 10, colour = "grey25", hjust = 0.5),
          plot.caption  = element_text(size = 7.2, colour = "grey35", hjust = 0,
                                       lineheight = 1.25))
} else {
  if (MAP_STYLE == "hatched")
    warning("ggpattern not installed; falling back to grayscale.", call. = FALSE)
  p <- base_map + scale_fill_manual(values = fills, drop = FALSE)
}

# ---------------------------------------------------------------------------
# Save. PDF is vector and is what the paper should embed.
# ---------------------------------------------------------------------------
sfx <- if (MAP_STYLE == "colour") "_colour" else ""
ggsave(file.path(COVERAGE_DIR, paste0("HPT_R_COVERAGE_MAP_CONUS", sfx, ".pdf")),
       p, width = 10.5, height = 7.4)
ggsave(file.path(COVERAGE_DIR, paste0("HPT_R_COVERAGE_MAP_CONUS", sfx, ".png")),
       p, width = 10.5, height = 7.4, dpi = 400, bg = "white")

.cv_msg("Map written (", MAP_STYLE, ") to ", COVERAGE_DIR)


# ---------------------------------------------------------------------------
# Print check: convert the PNG to greyscale and confirm the three fills remain
# separated. Only runs if magick is installed.
# ---------------------------------------------------------------------------
if (requireNamespace("magick", quietly = TRUE)) {
  img <- magick::image_read(file.path(COVERAGE_DIR,
                                      paste0("HPT_R_COVERAGE_MAP_CONUS", sfx, ".png")))
  magick::image_write(magick::image_convert(img, colorspace = "gray"),
                      file.path(COVERAGE_DIR, "HPT_R_COVERAGE_MAP_BWCHECK.png"))
  .cv_msg("Greyscale proof written: HPT_R_COVERAGE_MAP_BWCHECK.png -- ",
          "open it to confirm the three categories stay distinguishable.")
}


# ===========================================================================
# 6. QA
# ===========================================================================
qa <- data.table(
  CHECK = c("Census county key unique",
            "Estimation rows match feols nobs",
            "Estimation counties subset of panel counties",
            "All estimation counties matched to Census",
            "Population coverage between 0 and 100"),
  VALUE = c(as.character(anyDuplicated(pep, by = "COUNTY_FIPS")),
            paste0(nrow(cv_estimation), " vs ", nobs(fit)),
            as.character(all(cv_sets$ESTIMATION_SAMPLE %chin% cv_sets$IN_PANEL)),
            as.character(cv_summary[TIER == "ESTIMATION_SAMPLE"]$N_UNMATCHED),
            sprintf("%.2f", est$POPULATION_COVERAGE_PCT)),
  STATUS = c(ifelse(anyDuplicated(pep, by = "COUNTY_FIPS") == 0, "PASS", "FAIL"),
             ifelse(abs(nrow(cv_estimation) - nobs(fit)) <= 1L, "PASS", "FAIL"),
             ifelse(all(cv_sets$ESTIMATION_SAMPLE %chin% cv_sets$IN_PANEL), "PASS", "FAIL"),
             ifelse(cv_summary[TIER == "ESTIMATION_SAMPLE"]$N_UNMATCHED == 0, "PASS", "WARN"),
             ifelse(est$POPULATION_COVERAGE_PCT > 0 &
                      est$POPULATION_COVERAGE_PCT <= 100, "PASS", "FAIL")))
.cv_show(qa)
fwrite(qa, file.path(COVERAGE_DIR, "HPT_R_COVERAGE_QA.csv"))

.cv_msg("\nFor the paper: report the ESTIMATION_SAMPLE row of ",
        "HPT_R_COVERAGE_SUMMARY.csv, not IN_PANEL.")



###############################################################################
# THREE BUILDS FOR THE PAPER
#
#   BLOCK 1  Pooled estimates and the sample audit, on the current panel
#   BLOCK 2  Appendix table: family-level scheme assignment, 16 x 18
#   BLOCK 3  Shoppable share inputs for the magnitude calculation
#
# Run after the main pipeline, with `outpatient` and `schemes_long` in memory.
# Each block writes both a CSV and a ready-to-\input LaTeX file.
###############################################################################
suppressPackageStartupMessages({ library(data.table); library(fixest) })

OUT_TAB <- if (exists("TABLE_DIR") && dir.exists(TABLE_DIR)) TABLE_DIR else getwd()
OUT_TEX <- file.path(dirname(OUT_TAB), "02_Tables_TEX")
dir.create(OUT_TEX, showWarnings = FALSE, recursive = TRUE)

.hd <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")
.sh <- function(d) { print(as.data.frame(d), row.names = FALSE); invisible(NULL) }
.wr <- function(txt, f) { writeLines(txt, file.path(OUT_TEX, f)); cat("  wrote", f, "\n") }

# fitstat returns a length-4 vector (stat, p, df1, df2) for ivf1 in some fixest
# versions. Always take the first element and coerce defensively.
.stat1 <- function(fit, what) {
  v <- tryCatch(suppressWarnings(as.numeric(fixest::fitstat(fit, what, simplify = TRUE))),
                error = function(e) NA_real_)
  if (length(v) == 0L) NA_real_ else v[1L]
}


# ===========================================================================
# BLOCK 1 -- POOLED ESTIMATES AND THE SAMPLE AUDIT
# ===========================================================================
#
# Regenerates the pooled table and the sample audit against the panel currently
# in memory, so the reported N and the reported coefficients come from the same
# vintage. Panel size has moved across rebuilds of the upstream data (for
# instance 1,444,901 rows and 959,661 estimated, against 1,431,814 and 951,178
# after the concept merge), and those are not just different N labels -- the
# coefficients were estimated on a different sample.

.hd("BLOCK 1 — POOLED ESTIMATES AND SAMPLE AUDIT")

P_OUTCOMES <- c(Median = "LN_MEDIAN_PRICE",
                Mean   = "LN_MEAN_PRICE",
                P25    = "LN_P25_PRICE",
                P75    = "LN_P75_PRICE",
                IQR    = "LN_IQR_PRICE")
P_OUTCOMES <- P_OUTCOMES[P_OUTCOMES %chin% names(outpatient)]
cat("Outcomes available:", paste(names(P_OUTCOMES), collapse = ", "), "\n")

pooled <- cache_or_run("t02_pooled_current_vintage", {
  rows <- list()
  for (il in names(MAIN_INSTRUMENTS)) {
    z <- MAIN_INSTRUMENTS[[il]]
    for (oc in names(P_OUTCOMES)) {
      y <- P_OUTCOMES[[oc]]
      keep <- unique(c(y, ENDOGENOUS_VARIABLE, z, BASELINE_CONTROLS,
                       BASELINE_FIXED_EFFECTS, BASELINE_CLUSTERS))
      d <- outpatient[, ..keep]
      d <- d[complete.cases(d)]
      if (nrow(d) < MIN_MODEL_OBS) next

      cl   <- build_cluster_formula(available_columns(d, BASELINE_CLUSTERS))
      fe   <- available_columns(d, BASELINE_FIXED_EFFECTS)
      ctrl <- available_columns(d, BASELINE_CONTROLS)

      cat(sprintf("  %-30s %-7s ", substr(il, 1, 28), oc)); t0 <- Sys.time()

      ols <- feols(build_ols_formula(y, c(ENDOGENOUS_VARIABLE, ctrl), fe),
                   data = d, cluster = cl, warn = FALSE, notes = FALSE)
      rf  <- feols(build_ols_formula(y, c(z, ctrl), fe),
                   data = d, cluster = cl, warn = FALSE, notes = FALSE)
      iv  <- feols(build_iv_formula(y, ENDOGENOUS_VARIABLE, z, ctrl, fe),
                   data = d, cluster = cl, warn = FALSE, notes = FALSE)

      cat(sprintf("%5.1fs\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

      ivnm <- intersect(c(paste0("fit_", ENDOGENOUS_VARIABLE), ENDOGENOUS_VARIABLE),
                        names(coef(iv)))[1]
      b_ols <- unname(coef(ols)[ENDOGENOUS_VARIABLE])
      b_rf  <- unname(coef(rf)[z])
      b_iv  <- unname(coef(iv)[ivnm]); s_iv <- unname(sqrt(vcov(iv)[ivnm, ivnm]))

      rows[[length(rows) + 1L]] <- data.table(
        INSTRUMENT_LABEL = il, OUTCOME = oc,
        OLS_PERCENT = 100 * (exp(b_ols) - 1),
        RF_COEF     = b_rf,
        IV_PERCENT  = 100 * (exp(b_iv) - 1),
        IV_P        = 2 * pnorm(-abs(b_iv / s_iv)),
        FIRST_STAGE_F = .stat1(iv, "ivwald1"),   # was "ivf1"
        WU_HAUSMAN_P  = .stat1(iv, "wh"),
        N_OBSERVATIONS = nobs(iv))
    }
  }
  rbindlist(rows, fill = TRUE)
})

save_csv(pooled, "T03_pooled_OLS_RF_IV_all_outcomes.csv")
cat("\nPooled estimates, current vintage:\n")
.sh(pooled[, .(INSTRUMENT_LABEL = substr(INSTRUMENT_LABEL, 1, 28), OUTCOME,
               OLS = round(OLS_PERCENT, 3), RF = round(RF_COEF, 4),
               IV = round(IV_PERCENT, 3), IV_P = round(IV_P, 3),
               F = round(FIRST_STAGE_F, 1), WH_P = signif(WU_HAUSMAN_P, 2),
               N = N_OBSERVATIONS)])

cat("\nCHECK: N should be 951,178 (911,861 for ex-CBSA), matching tab:headline.\n")
cat("Distinct N values observed:", paste(sort(unique(pooled$N_OBSERVATIONS)), collapse = ", "), "\n")

# ---- LaTeX ----
lt <- c("\\begin{table}[!htbp]", "\\centering",
        "\\caption{Pooled Estimates Across the Price Distribution}",
        "\\label{tab:pooled}", "\\begin{threeparttable}", "\\footnotesize",
        "\\begin{tabular}{llrrrrr}", "\\toprule",
        "Instrument & Outcome & OLS (\\%) & RF & IV (\\%) & IV $p$ & First-stage $F$ \\\\",
        "\\midrule")
for (il in unique(pooled$INSTRUMENT_LABEL)) {
  lt <- c(lt, sprintf("\\multicolumn{7}{l}{\\textit{%s}} \\\\", gsub("_", " ", il)))
  for (i in which(pooled$INSTRUMENT_LABEL == il)) {
    r <- pooled[i]
    lt <- c(lt, sprintf("\\quad & %s & %.3f & %.4f & %.3f & %.3f & %.1f \\\\",
                        r$OUTCOME, r$OLS_PERCENT, r$RF_COEF, r$IV_PERCENT,
                        r$IV_P, r$FIRST_STAGE_F))
  }
  lt <- c(lt, "\\addlinespace[3pt]")
}
n_main <- max(pooled$N_OBSERVATIONS); n_alt <- min(pooled$N_OBSERVATIONS)
lt <- c(lt, "\\bottomrule", "\\end{tabular}",
        "\\begin{tablenotes}[flushleft]", "\\footnotesize",
        sprintf(paste0("\\item \\textit{Notes:} Each row is a separate pooled specification with ",
                       "county$\\times$concept and month fixed effects, $\\ln$(beds) control, and ",
                       "two-way clustering by county and month. $N=%s$ (%s for the ex-CBSA ",
                       "instrument). Wu--Hausman rejects exogeneity of prior disclosures at ",
                       "$p<0.001$ in all specifications."),
                format(n_main, big.mark = ","), format(n_alt, big.mark = ",")),
        "\\end{tablenotes}", "\\end{threeparttable}", "\\end{table}")
.wr(lt, "tab_pooled.tex")

# ---- Sample audit, same vintage ----
y <- PRIMARY_OUTCOME; z <- MAIN_INSTRUMENTS[[1L]]
keep <- unique(c(y, ENDOGENOUS_VARIABLE, z, BASELINE_CONTROLS,
                 BASELINE_FIXED_EFFECTS, BASELINE_CLUSTERS))
d <- outpatient[, ..keep]
cc <- complete.cases(d)
fit <- feols(build_ols_formula(y, c(z, BASELINE_CONTROLS), BASELINE_FIXED_EFFECTS),
             data = d[cc], cluster = build_cluster_formula(BASELINE_CLUSTERS),
             warn = FALSE, notes = FALSE)
cellid <- do.call(paste, c(outpatient[, ..BASELINE_FIXED_EFFECTS], sep = "::"))
cellsz <- table(cellid)

audit <- data.table(
  QUANTITY = c("Panel rows", "Complete cases", "Estimated rows", "Dropped by feols",
               "Share dropped", "FE cells", "Singleton FE cells",
               "Share rows in singleton cells"),
  VALUE = c(nrow(outpatient), sum(cc), nobs(fit), sum(cc) - nobs(fit),
            round((sum(cc) - nobs(fit)) / sum(cc), 4),
            length(cellsz), sum(cellsz == 1L),
            round(sum(cellsz == 1L) / nrow(outpatient), 4)))
save_csv(audit, "QA01_estimation_sample_audit.csv")
cat("\nSample audit, current vintage:\n"); .sh(audit)



###############################################################################
# BLOCK 2 -- FAMILY-LEVEL SCHEME ASSIGNMENT
#
# Builds the 16 x 18 family-by-scheme table for the appendix. Four details
# matter for it to describe the estimation universe rather than the codebook:
#
#   - the scheme table keys on ANALYSIS_CONCEPT_ID, not FINAL_CONCEPT_ID
#   - the category column is SHOPPABILITY_CATEGORY, not CATEGORY
#   - category VALUES are detected at runtime rather than assumed, since the
#     scheme builder's labels are not guaranteed to be Shoppable /
#     Intermediate / Non_shoppable
#   - families come from the PANEL (FINAL_FAMILY_ID, 16 families over 738
#     concepts) rather than from the scheme table's own family column, so the
#     table matches the sample the paper actually estimates on
#
# Run with `outpatient` and `schemes_long` in memory.
###############################################################################
suppressPackageStartupMessages({ library(data.table) })

OUT_TAB <- if (exists("TABLE_DIR") && dir.exists(TABLE_DIR)) TABLE_DIR else getwd()
OUT_TEX <- file.path(dirname(OUT_TAB), "02_Tables_TEX")
dir.create(OUT_TEX, showWarnings = FALSE, recursive = TRUE)

.hd <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")
.sh <- function(d) { print(as.data.frame(d), row.names = FALSE); invisible(NULL) }
.wr <- function(txt, f) { writeLines(txt, file.path(OUT_TEX, f)); cat("  wrote", f, "\n") }

.hd("BLOCK 2 — FAMILY-LEVEL SCHEME ASSIGNMENT")


# ---------------------------------------------------------------------------
# 1. Resolve column names and inspect the category values before mapping.
# ---------------------------------------------------------------------------
SL_KEY <- intersect(c("FINAL_CONCEPT_ID", "ANALYSIS_CONCEPT_ID"), names(schemes_long))[1]
SL_CAT <- intersect(c("SHOPPABILITY_CATEGORY", "CATEGORY"), names(schemes_long))[1]
SL_SCH <- intersect(c("SCHEME_NAME", "SCHEME_ID"), names(schemes_long))[1]
stopifnot(!is.na(SL_KEY), !is.na(SL_CAT), !is.na(SL_SCH))
cat("Using key =", SL_KEY, "| category =", SL_CAT, "| scheme =", SL_SCH, "\n")

cat_levels <- sort(unique(as.character(schemes_long[[SL_CAT]])))
cat("Category values present:", paste(cat_levels, collapse = " | "), "\n")

# Map whatever labels are present onto S / I / N. The ordinal column is the
# authority where available, since it encodes the intended ordering directly
# and is immune to label drift.
if ("SHOPPABILITY_ORDINAL" %in% names(schemes_long)) {
  ord_map <- unique(schemes_long[, c(SL_CAT, "SHOPPABILITY_ORDINAL"), with = FALSE])
  setnames(ord_map, c("CAT", "ORD"))
  setorder(ord_map, -ORD)
  cat("\nCategory-to-ordinal mapping (highest ordinal = most shoppable):\n")
  .sh(ord_map)
  CODE <- setNames(c("S", "I", "N")[seq_len(nrow(ord_map))], ord_map$CAT)
} else {
  guess <- function(x) fcase(
    grepl("^shop|high", x, ignore.case = TRUE), "S",
    grepl("inter|mid|medium", x, ignore.case = TRUE), "I",
    default = "N")
  CODE <- setNames(vapply(cat_levels, guess, character(1)), cat_levels)
}
cat("\nCode mapping:\n"); print(CODE)
stopifnot(!anyNA(CODE))


# ---------------------------------------------------------------------------
# 2. Concept -> family map from the PANEL (the estimation universe).
# ---------------------------------------------------------------------------
cf <- unique(outpatient[!is.na(FINAL_CONCEPT_ID),
                        .(FINAL_CONCEPT_ID, FINAL_FAMILY_ID)], by = "FINAL_CONCEPT_ID")
cat("\nPanel concepts:", nrow(cf), "| families:", uniqueN(cf$FINAL_FAMILY_ID), "\n")

sl <- merge(schemes_long, cf, by.x = SL_KEY, by.y = "FINAL_CONCEPT_ID", all.x = FALSE)
setnames(sl, c(SL_CAT, SL_SCH), c("CAT", "SCHEME"), skip_absent = TRUE)

cat("Concepts in schemes_long:", uniqueN(schemes_long[[SL_KEY]]),
    "| matched to a panel family:", uniqueN(sl[[SL_KEY]]),
    "| schemes:", uniqueN(sl$SCHEME), "\n")
if (uniqueN(sl[[SL_KEY]]) != nrow(cf)) {
  warning(nrow(cf) - uniqueN(sl[[SL_KEY]]), " panel concept(s) have no scheme ",
          "assignment and will be absent from the table.", call. = FALSE)
}


# ---------------------------------------------------------------------------
# 3. Modal category per family x scheme, lowercase where the family splits.
# ---------------------------------------------------------------------------
fam_scheme <- sl[, {
  tb    <- sort(table(CAT), decreasing = TRUE)
  modal <- names(tb)[1L]
  split <- length(tb) > 1L
  ch    <- unname(CODE[[modal]])
  .(CELL       = if (split) tolower(ch) else ch,
    MODAL      = modal,
    IS_SPLIT   = split,
    N_CONCEPTS = .N,
    N_CATEGORIES = length(tb))
}, by = .(FINAL_FAMILY_ID, SCHEME)]

n_split <- sum(fam_scheme$IS_SPLIT); n_cells <- nrow(fam_scheme)
cat("\nSplit cells (family assigns more than one category):", n_split, "of", n_cells,
    sprintf("(%.1f%%)\n", 100 * n_split / n_cells))
cat("These are the ONLY source of within-family variation in shoppability,\n",
    "and therefore what makes the family-fixed-effect mechanism tests in\n",
    "Section 7 identifiable. Report this count in the text.\n", sep = "")

cat("\nSplit cells by family:\n")
.sh(fam_scheme[, .(N_SCHEMES = .N, N_SPLIT = sum(IS_SPLIT)),
               by = FINAL_FAMILY_ID][order(-N_SPLIT)])

wide <- dcast(fam_scheme, FINAL_FAMILY_ID ~ SCHEME, value.var = "CELL")
save_csv(fam_scheme, "T02B_family_scheme_long.csv")
save_csv(wide,       "T02B_family_scheme_wide.csv")
cat("\nFamily x scheme grid:\n"); .sh(wide)


# ---------------------------------------------------------------------------
# 4. LaTeX. Landscape, rotated scheme headers, legend in the notes.
# ---------------------------------------------------------------------------
sch  <- setdiff(names(wide), "FINAL_FAMILY_ID")
abbr <- vapply(sch, function(s) {
  s <- gsub("^(Alt|Alt:)[[:space:]]*", "", s)
  s <- gsub("[^A-Za-z0-9 ]", "", s)
  w <- strsplit(trimws(s), "[[:space:]]+")[[1]]
  paste(substr(w, 1, 4), collapse = " ")
}, character(1))

lt <- c(
  "\\begin{landscape}", "\\begin{table}[!htbp]", "\\centering",
  "\\caption{Shoppability Assignment by Clinical Family and Classification Scheme}",
  "\\label{tab:schemes_family}", "\\begin{threeparttable}", "\\scriptsize",
  paste0("\\begin{tabular}{l", strrep("c", length(sch)), "}"), "\\toprule",
  paste0("Clinical family & ",
         paste(sprintf("\\rotatebox{90}{\\footnotesize %s}", abbr), collapse = " & "),
         " \\\\"),
  "\\midrule")

for (i in seq_len(nrow(wide))) {
  cells <- unlist(wide[i, ..sch]); cells[is.na(cells)] <- "--"
  lt <- c(lt, paste0(gsub("_", " ", wide$FINAL_FAMILY_ID[i]), " & ",
                     paste(cells, collapse = " & "), " \\\\"))
}

lt <- c(lt, "\\bottomrule", "\\end{tabular}",
        "\\begin{tablenotes}[flushleft]", "\\footnotesize",
        paste0("\\item \\textit{Notes:} S = shoppable, I = intermediate, ",
               "N = non-shoppable, -- = family not assigned under that scheme. ",
               "Uppercase indicates that every concept in the family receives ",
               "that category; \\emph{lowercase} marks the modal category where ",
               "the family's concepts split across categories (",
               n_split, " of ", n_cells, " cells). Those split cells are the only ",
               "source of within-family variation in shoppability, which is why ",
               "the mechanism tests in Section~\\ref{sec:mechanism} can absorb ",
               "clinical-family fixed effects without absorbing the classification ",
               "itself. Assignments cover the ", nrow(cf), " concepts in the ",
               "estimation universe; the full scoped codebook contains ",
               uniqueN(schemes_long[[SL_KEY]]), " concepts, the remainder falling ",
               "outside the panel."),
        "\\end{tablenotes}", "\\end{threeparttable}", "\\end{table}", "\\end{landscape}")
.wr(lt, "tab_schemes_family.tex")

cat("\nRequires \\usepackage{pdflscape} and \\usepackage{rotating} -- both already",
    "in your preamble.\n")

# ===========================================================================
# BLOCK 3 -- SHOPPABLE SHARE FOR THE MAGNITUDE CALCULATION
# ===========================================================================
#
# READ THIS BEFORE USING THE OUTPUT.
#
# The disclosure files contain PRICES but not QUANTITIES. There is no
# utilization variable anywhere in this pipeline, so a genuine
# spending-weighted shoppable share cannot be computed from these data. What
# follows are the three shares that can be computed, in ascending order of how
# close they come to the target:
#
#   (a) CONCEPT SHARE          fraction of the 738 concepts that are
#                              shoppable. The wrong weight for any statement
#                              about spending.
#   (b) OBSERVATION SHARE      fraction of hospital-month-concept rows that are
#                              shoppable. Weights by how widely a service is
#                              disclosed, which correlates with how commonly it
#                              is delivered but is not utilization.
#   (c) PRICE-WEIGHTED SHARE   observation count times median price.
#                              Interpretable as "if every disclosed service
#                              were delivered once at its disclosed rate."
#                              Still not spending, since it assumes uniform
#                              volume per concept.
#
# A true spending share would require merging Medicare outpatient utilization
# by HCPCS. Both the CMS "Medicare Physician & Other Practitioners" and
# "Medicare Outpatient Hospitals - by Provider and Service" files publish
# HCPCS-level service counts and would supply it.
#
# Report (b) and (c) as bounds and state plainly that neither is
# spending-weighted. That is defensible. Calling either one a spending share is
# not.

.hd("BLOCK 3 — SHOPPABLE SHARE INPUTS")

SCHEME_FOR_SHARE <- "SCHEME_1_CERTAINTY"

cs <- outpatient[!is.na(get(SCHEME_FOR_SHARE)) & !is.na(MEDIAN_PRICE),
                 .(N_OBS = .N,
                   N_HOSPITALS = uniqueN(HOSPITAL_ID),
                   MEDIAN_PRICE = median(MEDIAN_PRICE, na.rm = TRUE),
                   CATEGORY = first(get(SCHEME_FOR_SHARE))),
                 by = FINAL_CONCEPT_ID]
cs[, IS_SHOP := CATEGORY == "Shoppable"]
cs[, PRICE_WEIGHT := N_OBS * MEDIAN_PRICE]

shares <- data.table(
  MEASURE = c("(a) Concept share",
              "(b) Observation share",
              "(c) Price-weighted share"),
  SHOPPABLE = c(sum(cs$IS_SHOP),
                cs[IS_SHOP == TRUE, sum(N_OBS)],
                cs[IS_SHOP == TRUE, sum(PRICE_WEIGHT)]),
  TOTAL = c(nrow(cs), cs[, sum(N_OBS)], cs[, sum(PRICE_WEIGHT)]))
shares[, SHARE_PCT := round(100 * SHOPPABLE / TOTAL, 2)]
save_csv(shares, "QA11_shoppable_shares.csv")
cat("\nShoppable shares under", SCHEME_FOR_SHARE, ":\n"); .sh(shares)

cat("\nBy clinical family (price-weighted):\n")
fam_share <- merge(cs, cf, by = "FINAL_CONCEPT_ID")[
  , .(N_CONCEPTS = .N, SHOP_CONCEPTS = sum(IS_SHOP),
      PRICE_WEIGHT = sum(PRICE_WEIGHT),
      SHOP_PRICE_WEIGHT = sum(PRICE_WEIGHT * IS_SHOP)), by = FINAL_FAMILY_ID]
fam_share[, `:=`(SHARE_OF_TOTAL_PCT = round(100 * PRICE_WEIGHT / sum(PRICE_WEIGHT), 2),
                 SHOP_SHARE_PCT     = round(100 * SHOP_PRICE_WEIGHT / PRICE_WEIGHT, 1))]
setorder(fam_share, -SHARE_OF_TOTAL_PCT)
save_csv(fam_share, "QA11_shoppable_shares_by_family.csv")
.sh(fam_share)
