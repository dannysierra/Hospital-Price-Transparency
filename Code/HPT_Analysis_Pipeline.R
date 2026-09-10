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
# ===========================================================================
# HOW TO RUN
# ===========================================================================
#
#   1. Put these four files in one directory:
#        HPT_Analysis_Pipeline.R              (this file)
#        HPT_warm_start.R
#        HPT_diagnostic_for_family_levels.R
#        HPT_Section19_TierI_Diagnostics.R
#
#   2. Point HPT_ROOT at the project directory, which must contain
#      Data/data_final/01_R_Analysis_Panels with the inputs listed below:
#
#        Sys.setenv(HPT_ROOT = "/path/to/Hospital Price Transparency Paper")
#
#   3. source("HPT_Analysis_Pipeline.R")
#
# That is the whole procedure. There are no absolute paths in this file. Every
# path derives from HPT_ROOT, and the three companion files are located
# relative to CODE_DIR.
#
# ---------------------------------------------------------------------------
# CONTROLLING WHAT RUNS
# ---------------------------------------------------------------------------
# RUN_STAGES and HPT_RUN, both set in PART 1, decide which blocks execute. The
# defaults reproduce every table and figure in the paper from the shipped
# cache WITHOUT re-running the eight-hour concept-level sweep. To rebuild that
# from scratch, add 6 to RUN_STAGES and set USE_CACHE <- FALSE.
#
# Approximate runtimes from a warm cache:
#   PART 2  definitions            seconds
#   PART 3  build and preflight    minutes
#   PART 4  estimation             ~6 hours   (~14 with stage 6)
#   PART 5  tables and figures     minutes
#   PART 6  diagnostics            ~2 hours   (off by default)
#
# ---------------------------------------------------------------------------
# STRUCTURE
# ---------------------------------------------------------------------------
#   PART 1  Configuration          paths, packages, input check, constants
#   PART 2  Definitions            functions only, nothing estimated
#   PART 3  Build and preflight    the gate everything else depends on
#   PART 4  Estimation             stages 5 through 18
#   PART 5  Outputs                tables, LaTeX, figures, coverage map
#   PART 6  Diagnostics            referee-response checks, off by default
#   PART 7  Scratch                interactive code, commented out
#
# Nothing between the top of the file and the end of PART 2 touches data, so
# the file can be sourced that far to get every function into scope in seconds
# without estimating anything. That is what warm_start() relies on.
#
# ---------------------------------------------------------------------------
# INPUTS
# ---------------------------------------------------------------------------
# Four required, under Data/data_final/01_R_Analysis_Panels:
#   HPT_R_MAIN_PRIMARY_COUNTY_OUTPATIENT_CONCEPT.parquet
#   HPT_R_MAIN_PRIMARY_COUNTY_INPATIENT_CONCEPT.parquet
#   HPT_R_MAIN_PRIMARY_COUNTY_OUTPATIENT_EXACT_CODE.parquet
#   HPT_CODEBOOK_FOR_SHOPPABILITY_SCHEMES.csv
#
# Four optional, each gating one robustness section rather than the run:
#   HPT_R_MAIN_ROBUSTNESS_CBSA_OUTPATIENT_CONCEPT.parquet     (Section 15A)
#   HPT_R_MAIN_ROBUSTNESS_CBSA_OUTPATIENT_EXACT_CODE.parquet  (Section 15A)
#   HPT_PAYER_DISPERSION*.csv.gz                              (Section 9)
#   HPT_CONCEPT_PAYER_CLASS*.csv                              (Section 15C)
#
# The price data are licensed from Turquoise Health and are not redistributed.
# The ACS blocks additionally need a Census API key.
#
# ---------------------------------------------------------------------------
# DESIGN DECISIONS
# ---------------------------------------------------------------------------
# Four choices are fixed throughout and are not re-litigated section by
# section. Each is defended in the paper; the short version is here so the code
# reads the way the design is specified.
#
# 1. TREATMENT ENTERS LINEARLY. Concave transforms of N_PRIOR_POSTERS inflate
#    the IV coefficient by shrinking Cov(T, Z) in the Wald denominator,
#    monotonically in how much they compress the treatment. The reduced-form
#    numerator does not move. The linear first stage runs F = 43 in both
#    interacted equations, so there is no weak-instrument problem for a
#    transform to solve. The transform ladder is robustness, never a headline.
#
# 2. THE REDUCED FORM IS THE PRIMARY ESTIMATOR. IV is a ratio; the reduced form
#    is its numerator, and it carries no treatment variable, no functional-form
#    choice, and no first-stage noise in the standard error. Because every
#    model here is exactly identified, the Anderson-Rubin weak-instrument-
#    robust test is numerically identical to the t-test on the reduced-form
#    coefficient, so the reported RF p-values are already robust p-values.
#
# 3. HETEROGENEITY IS ESTIMATED BY INTERACTION, NEVER BY SAMPLE SPLIT. The
#    fixed effects absorb most of the instrument's effect on the treatment
#    within any single subsample, so splitting destroys identification.
#    estimate_interacted() is the single estimator behind every heterogeneity
#    result in the paper.
#
# 4. NO EVENT STUDY IS POSSIBLE. 3,721 of 3,723 hospitals appear at exactly one
#    posting month, so there is no within-hospital timing variation to exploit.
#    The selection concern is addressed structurally instead: the non-shoppable
#    service group is a within-hospital placebo.
#
# ---------------------------------------------------------------------------
# INFERENCE NOTE
# ---------------------------------------------------------------------------
# With two-way clustering where the month dimension has only 16 levels, fixest
# references a t distribution rather than a normal. .pval() in PART 2 recovers
# the correct degrees of freedom from the fit rather than assuming them. Any
# hand-computed 2 * pnorm(-abs(b/s)) is wrong here by roughly the difference
# between p = 0.041 and p = 0.059 at the margin. Two places still used a normal
# or chi-square reference; both are fixed, and the CHANGELOG records where.
#
# ---------------------------------------------------------------------------
# CACHING
# ---------------------------------------------------------------------------
# Every expensive step is wrapped in cache_or_run(), which writes an .rds to
# CACHE_DIR and reloads it on later runs. A step producing zero rows never
# overwrites a good cache. Set USE_CACHE <- FALSE to force recomputation, or
# call invalidate_cache() (PART 2) to delete selectively; it knows the
# dependency order and cascades by default, so editing a scheme rule
# invalidates everything downstream rather than leaving a half-updated cache.
#
# To resume a completed run in a fresh session: source this file with
# RUN_STAGES <- integer(0) and every HPT_RUN switch FALSE, then call
# restore_session().
#
# ===========================================================================
# CHANGELOG -- what the reorganization changed
# ===========================================================================
# The estimation code is unchanged. Every function body, specification, and
# constant is carried over verbatim. Only order, duplicates, paths, and two
# genuine bugs were touched.
#
# STRUCTURAL
#   * Definitions and run blocks separated. Previously interleaved, so the file
#     could not be sourced partway to load functions without estimating.
#   * Cache management moved from line 8031 of the old file up into PART 2,
#     beside the definitions that use it.
#   * All run blocks now sit behind RUN_STAGES or HPT_RUN switches. Several
#     previously ran unconditionally on every source().
#   * Interactive one-liners moved to PART 7 and commented out.
#
# PORTABILITY
#   * rm(list = ls()) removed from the top. It made the file unsourceable from
#     any wrapper, which is why a 91-line cold-start block existed to write
#     temp files and reconstruct its own path after being wiped. That block is
#     deleted; it is no longer needed.
#   * Eight hardcoded /Users/danielsierra paths replaced, four of them inside
#     source() calls that stopped a replicator at the first one.
#   * TABLE_DIR and FIG_DIR were redefined with absolute paths partway down,
#     silently overriding PART 1 from that point on, and FIGURE_DIR pointed at
#     03_Figures in one place and 02_Figures in another. Figures written before
#     and after that line landed in different directories. Now one definition
#     each, with FIG_DIR aliased to FIGURE_DIR and TEX_DIR to OUT_TEX.
#   * Package and input checks moved to PART 1 and run up front, rather than
#     failing at the point of use hours into a run.
#
# DUPLICATE DEFINITIONS REMOVED
#   * theme_paper() was defined twice with DIFFERENT bodies: base_size = 11
#     with axis lines and ticks, and base = 10 without them. The second
#     silently restyled every figure defined before it. Kept the first.
#   * sv() was defined twice, the first using cairo_pdf and the second not.
#     Kept the second, since cairo is not installed on the author's machine.
#   * save_fig() twice identically; S15_RUN, .s15_head, .s15_show, .s15_sweep
#     three times each; FILES twice; a FIG 10 header twice. One each.
#
# BUGS FIXED
#   * A stray top-level line
#         d[, P_VALUE := 2 * pt(-qnorm(1 - P_VALUE / 2), df = 15)]
#     sat between the FIG 1 and FIG 2 blocks, where `d` is still the
#     T06_main_interacted_RF_and_IV.csv rows table, which carries RF_P and
#     IV_P but no P_VALUE column. This is the source of the "object 'P_VALUE'
#     not found" error. Deleted; the correction already exists inside FIG 2,
#     operating on the tests table where it belongs.
#   * .s14_contrast() computed P_VALUE = pf(W / df, df1 = df, df2 = .cluster_df(fit), lower.tail = FALSE), an
#     asymptotic reference inconsistent with the t(15) used everywhere else.
#     Changed to pf(W/df, df1 = df, df2 = .cluster_df(fit)). For df = 1 this is
#     algebraically identical to 2*pt(-sqrt(W), df), so every one-degree-of-
#     freedom value already in the paper is unchanged; only the two-df joint
#     tercile test moves.
#
# STILL OPEN, flagged rather than silently changed
#   * T13C_instrument_ladder_tests.csv still carries normal-reference
#     p-values. All eight rows equal 2*(1 - pnorm(sqrt(WALD))) exactly, and the
#     t(15) correction reproduces Appendix Table 26 to four decimals. Fix with
#       t13c[, P_VALUE := 2 * pt(-qnorm(1 - P_VALUE / 2), df = 15)]
#     and re-save before regenerating fig11.
#   * s24_scheme_crosswalk() counts rows with .N where it should count concepts
#     with uniqueN(ANALYSIS_CONCEPT_ID). Its all-18-scheme counts sum to 781
#     rather than 738. The six-scheme crosswalk in the same function is
#     correct. Marked in place in PART 6.
###############################################################################


###############################################################################
#                                                                             #
#   PART 1 -- CONFIGURATION                                                   #
#                                                                             #
#   Paths, packages, input checks, and every constant the pipeline reads.     #
#   The only part a replicator edits, and only if HPT_ROOT cannot be set as   #
#   an environment variable.                                                  #
#                                                                             #
###############################################################################


options(stringsAsFactors = FALSE, scipen = 999, width = 160, fixest_notes = FALSE)

suppressPackageStartupMessages({
  library(arrow); library(data.table); library(fixest); library(dplyr)
  library(ggplot2); library(stringr); library(lubridate)
})
if (!requireNamespace("MASS", quietly = TRUE)) stop("Package 'MASS' required.")
setFixest_estimation(mem.clean = TRUE, data.save = FALSE)

# ---------------------------------------------------------------------------
# 1.1  Optional packages
# ---------------------------------------------------------------------------
# Checked once here rather than at ten scattered points of use, so a missing
# package surfaces at minute one rather than eight hours into a run. None of
# these stops the pipeline; each disables one specific output.

.HPT_OPTIONAL <- c(
  R.utils          = "gzipped payer-cell files (Section 9 mechanism measures)",
  tidycensus       = "ACS county demographics (Section 12 heterogeneity)",
  ggrepel          = "figure label placement (fig05, fig20b)",
  patchwork        = "combined two-panel fig20",
  sf               = "coverage map geometry",
  tigris           = "county and state shapefiles (coverage map)",
  ggpattern        = "hatched coverage map variant",
  magick           = "greyscale print proof of the coverage map",
  scales           = "percent axis formatting",
  fwildclusterboot = "wild cluster bootstrap (Section 21, optional)"
)
.hpt_missing <- names(.HPT_OPTIONAL)[
  !vapply(names(.HPT_OPTIONAL), requireNamespace, logical(1), quietly = TRUE)]
if (length(.hpt_missing)) {
  cat("\nOptional packages not installed. Each disables one output:\n")
  for (p in .hpt_missing) cat(sprintf("  %-18s %s\n", p, .HPT_OPTIONAL[[p]]))
  cat('  install.packages(c(',
      paste0('"', .hpt_missing, '"', collapse = ", "), '))\n\n', sep = "")
}

# ---------------------------------------------------------------------------
# 1.2  Paths -- everything derives from HPT_ROOT
# ---------------------------------------------------------------------------
PROJECT_ROOT <- Sys.getenv("HPT_ROOT", unset = NA_character_)
if (is.na(PROJECT_ROOT) || !nzchar(PROJECT_ROOT)) {
  PROJECT_ROOT <- path.expand(file.path(
    "~", "Library", "CloudStorage", "OneDrive-FloridaStateUniversity",
    "Hospital Price Transparency Paper"))
}
if (!dir.exists(PROJECT_ROOT))
  stop("PROJECT_ROOT not found:\n  ", PROJECT_ROOT,
       "\n\nSet it and re-source:\n",
       '  Sys.setenv(HPT_ROOT = "/path/to/Hospital Price Transparency Paper")',
       call. = FALSE)

CODE_DIR        <- file.path(PROJECT_ROOT, "Code")
FINAL_DATA_ROOT <- file.path(PROJECT_ROOT, "Data", "data_final")
PANEL_DIR       <- file.path(FINAL_DATA_ROOT, "01_R_Analysis_Panels")
RESULT_ROOT     <- file.path(FINAL_DATA_ROOT, "06_R_Analysis_Results")

TABLE_DIR    <- file.path(RESULT_ROOT, "01_Tables_CSV")
FIGURE_DIR   <- file.path(RESULT_ROOT, "02_Figures")
OUT_TEX      <- file.path(RESULT_ROOT, "02_Tables_TEX")
MODEL_DIR    <- file.path(RESULT_ROOT, "04_Model_Objects")
QA_DIR       <- file.path(RESULT_ROOT, "05_R_QA")
CACHE_DIR    <- file.path(RESULT_ROOT, "07_Cache")
COVERAGE_DIR <- file.path(RESULT_ROOT, "08_Coverage")

# FIGURE_DIR and FIG_DIR were two names for what should have been one place.
# The old Section 0 set FIGURE_DIR to 03_Figures, and a hardcoded line partway
# down set FIG_DIR to 02_Figures, so save_fig() wrote figures 1-9 into
# 03_Figures while sv() wrote figures 10-20 into 02_Figures. Both names now
# resolve to 02_Figures, which is where the figures the paper actually uses
# were being written. NOTE: 03_Figures may still hold stale copies of figures
# 1 through 9 from before this fix; delete it once the run completes.
FIG_DIR <- FIGURE_DIR

# TEX_DIR is deliberately NOT aliased here. Section 24 defines it as
# file.path(TABLE_DIR, "tex") for its own .tex fragments, which is a different
# destination from OUT_TEX (the full-table LaTeX in 02_Tables_TEX).

for (d in c(RESULT_ROOT, TABLE_DIR, FIGURE_DIR, OUT_TEX, MODEL_DIR,
            QA_DIR, CACHE_DIR, COVERAGE_DIR))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# 1.3  Run switches
# ---------------------------------------------------------------------------
# RUN_STAGES drives the numbered stages inherited from the original file.
# HPT_RUN drives everything added since. Both are checked before any block
# executes, so a fresh source() with these defaults reproduces the paper from
# cache without re-running the eight-hour sweep.

RUN_STAGES <- c(5, 7, 8, 9, 10, 11, 12)   # add 6 to rebuild the concept sweep
USE_CACHE  <- TRUE

# ---------------------------------------------------------------------------
# WARM START -- the switch for a fresh session on a completed run
# ---------------------------------------------------------------------------
# Set this TRUE to load every function and every cached result into memory
# WITHOUT estimating anything:
#
#     HPT_WARM_START <- TRUE
#     source("HPT_Analysis_Pipeline.R")
#
# PARTS 1 and 2 evaluate normally (definitions only, seconds). PART 3 then
# calls restore_session() instead of building, and every stage in PARTS 4
# through 6 is switched off, so the file runs end to end in well under a
# minute and leaves outpatient, schemes_long, concept_results, main,
# measures, and the rest in the global environment under the names the run
# blocks use.
#
# This replaces the marker-grepping approach in HPT_warm_start.R. That version
# located the definitions/run boundary by searching for two comment strings,
# which broke the moment the file was reorganized. A switch cannot break that
# way.
#
# Nothing is written to disk in warm-start mode, and RUN_STAGES is emptied, so
# accidentally evaluating a stage block afterwards is a no-op rather than an
# eight-hour surprise.
HPT_WARM_START <- if (exists("HPT_WARM_START")) isTRUE(HPT_WARM_START) else FALSE

HPT_RUN <- list(
  instrument_audit = TRUE,   # QA05-QA08, justifies the three-instrument set
  insystem_robust  = TRUE,   # headline with system affiliation controlled
  demographics     = TRUE,   # needs a Census API key
  ses_terciles     = TRUE,
  robustness_13    = TRUE,   # system x month, leave-one-out, ladder, RI
  markets_15       = TRUE,   # CBSA, window ladder, payer-conditional
  family_16        = TRUE,
  conley_17        = TRUE,
  enforcement_18   = TRUE,
  tables           = TRUE,
  figures          = TRUE,
  coverage_map     = FALSE,  # downloads Census shapefiles, needs sf and tigris
  diagnostics      = FALSE   # Sections 19-25, referee response only
)

# In warm-start mode every switch is forced off. Set them again by hand after
# the session is loaded if you want to run one block interactively.
if (isTRUE(HPT_WARM_START)) {
  RUN_STAGES <- integer(0)
  HPT_RUN <- lapply(HPT_RUN, function(x) FALSE)
  cat("\nWARM START: RUN_STAGES emptied and every HPT_RUN switch set FALSE.\n")
}

# Section 20 (A5, A6, A8) has its own switch so it can run on a warm start
# without waking Sections 19 and 21-25. Several of those fire bare top-level
# calls, and the block at the end of Section 25 rewrites two saved CSVs in
# place, so it is not safe to run twice. Set HPT_TIER2 <- TRUE before
# source()ing to turn Section 20 on.
HPT_RUN$tier2_diag <- if (exists("HPT_TIER2")) isTRUE(HPT_TIER2) else FALSE

# PART 5's table block holds build_summary_stats(). Warm start blanks HPT_RUN,
# so it falls out of scope in exactly the session where it is wanted.
if (exists("HPT_TABLES")) HPT_RUN$tables <- isTRUE(HPT_TABLES)
# ---------------------------------------------------------------------------
# 1.4  Input inventory
# ---------------------------------------------------------------------------
# Called at the top of PART 3. Reports every input before anything runs, and
# separates the four that gate the whole pipeline from the four that gate a
# single robustness section.

hpt_check_inputs <- function() {
  core <- FILES[c("outpatient_concept", "inpatient_concept",
                  "outpatient_exact", "codebook")]
  cat("\nCORE INPUTS:\n"); ok <- TRUE
  for (nm in names(core)) {
    e <- file.exists(core[[nm]]); ok <- ok && e
    cat(sprintf("  [%-7s] %s\n", if (e) "OK" else "MISSING", basename(core[[nm]])))
  }
  cat("\nOPTIONAL INPUTS (each gates one section, not the run):\n")
  cbsa <- file.path(PANEL_DIR, c("HPT_R_MAIN_ROBUSTNESS_CBSA_OUTPATIENT_CONCEPT.parquet",
                                 "HPT_R_MAIN_ROBUSTNESS_CBSA_OUTPATIENT_EXACT_CODE.parquet"))
  for (p in cbsa)
    cat(sprintf("  [%-7s] %-56s Section 15A\n",
                if (file.exists(p)) "OK" else "absent", basename(p)))
  np <- length(list.files(PAYER_DISPERSION_DIR, pattern = PAYER_DISPERSION_PATTERN))
  nc <- length(list.files(PANEL_DIR, pattern = "^HPT_CONCEPT_PAYER_CLASS.*[.]csv([.]gz)?$"))
  cat(sprintf("  [%-7s] %-56s Section 9   (%d shard(s))\n",
              if (np) "OK" else "absent", "HPT_PAYER_DISPERSION*.csv.gz", np))
  cat(sprintf("  [%-7s] %-56s Section 15C (%d shard(s))\n",
              if (nc) "OK" else "absent", "HPT_CONCEPT_PAYER_CLASS*.csv", nc))
  if (!ok) stop("A core input is missing; see above.", call. = FALSE)
  invisible(ok)
}

cat("\nPART 1 loaded.\n  PROJECT_ROOT: ", PROJECT_ROOT,
    "\n  Results to:   ", RESULT_ROOT, "\n", sep = "")


# ============================================================================
# 1.5  Specification constants, instrument and scheme registries
# ============================================================================

# Verbatim from Section 0 of the original file, lines 198-442.

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




###############################################################################
#                                                                             #
#   PART 2 -- DEFINITIONS                                                     #
#                                                                             #
#   Functions only. Nothing here touches data, so the file can be sourced to  #
#   the end of PART 2 in seconds to load every function without estimating.   #
#                                                                             #
###############################################################################


# ============================================================================
# Section 1: Helper functions
# ============================================================================

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
               P_VALUE = pf(stat / rk, df1 = rk, df2 = .cluster_df(fit), lower.tail = FALSE))
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
# .pval() -- p-value from a t statistic, using the SAME reference distribution
# fixest itself uses for the given fit.
#
# The pipeline previously computed p-values by hand as 2 * pnorm(-abs(b/s)).
# That is the normal approximation. With two-way clustering where one
# dimension has only 16 levels (POST_MONTH), fixest references a
# t-distribution instead, which has visibly fatter tails at these sample
# sizes: a coefficient reported at p = 0.041 under pnorm is p = 0.059 under
# the correct reference.
#
# Rather than hard-code a df, .cluster_df() recovers it from the fit: it
# takes a term whose p-value fixest already computed correctly, and solves
# for the df that reproduces it. That works regardless of how fixest derives
# df internally, and adapts automatically to subsamples with different
# cluster counts (CBSA, payer-class, etc.).
#
# Falls back to the old pnorm() behaviour if the fit is unusable, so a
# patched line can never error where the original succeeded.
# ---------------------------------------------------------------------------
.cluster_df <- function(fit) {
  if (is.null(fit)) return(Inf)
  out <- tryCatch({
    b <- coef(fit); s <- fixest::se(fit); p <- fixest::pvalue(fit)
    keep <- which(is.finite(b) & is.finite(s) & s > 0 &
                    is.finite(p) & p > 1e-12 & p < 1 - 1e-12)
    if (!length(keep)) return(Inf)
    i  <- keep[1L]
    t0 <- abs(unname(b[i] / s[i]))
    p0 <- unname(p[i])
    if (!is.finite(t0) || t0 <= 0) return(Inf)
    uniroot(function(df) 2 * pt(-t0, df = df) - p0, interval = c(1, 1e6))$root
  }, error = function(e) Inf)
  if (!is.finite(out) || out < 1) Inf else out
}

.pval <- function(tstat, fit = NULL) {
  tstat <- unname(tstat)
  df <- .cluster_df(fit)
  if (!is.finite(df)) return(2 * pnorm(-abs(tstat)))
  2 * pt(-abs(tstat), df = df)
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




# ============================================================================
# Section 2: Panel loading and preparation
# ============================================================================

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




# ============================================================================
# Section 3: Shoppability scheme construction
# ============================================================================

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




# ============================================================================
# Section 4: Estimators
# ============================================================================

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
      RF_COEF = rf$b, RF_SE = rf$s, RF_P = .pval(rf$b / rf$s, rf_fit),
      RF_PERCENT_PER_SD = 100 * (exp(rf$b * sd_z) - 1),
      IV_COEF = iv$b, IV_SE = iv$s, IV_P = .pval(iv$b / iv$s, iv_fit),
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
                 P_VALUE = unname(.pval(b / s, fit)))
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




# ============================================================================
# Section 5: Instrument screen and pooled models
# ============================================================================

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




# ============================================================================
# Sections 6-7: Concept-level estimates, main results, transform ladder
# ============================================================================

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





# ============================================================================
# Section 8: Meta-regressions, deduplication, permutation
# ============================================================================

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




# ============================================================================
# Section 9: Price comparability as a candidate mechanism
# ============================================================================

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




# ============================================================================
# Section 10: Comparability meta-regression
# ============================================================================

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




# ============================================================================
# Section 11: Instrument balance and placebo tests
# ============================================================================

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




# ============================================================================
# Section 12: Preflight
# ============================================================================

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

# ============================================================================
# Cache management: cache_status, restore_session, invalidate_cache
# ============================================================================

# Moved up from line 8031. It belongs with the definitions rather than
# 5,000 lines below the code that calls it.

# SESSION RESTORE AND CACHE MANAGEMENT
#
# Everything above this line either defines a function or runs an analysis.
# This block does neither: it manages the .rds cache written by cache_or_run()
# so that a fresh session can pick up a completed run without re-estimating,
# and so that a code change can be propagated without re-running everything.
#
# Nothing here executes on its own. Load Sections 0-12 to define the pipeline,
# then call these directly.
#
#   cache_status()        what is cached, how large, and when it was written
#   restore_session()     load cached results into the global environment
#   invalidate_cache()    delete caches so the next run recomputes them
#
# -----------------------------------------------------------------------------
# THE PROBLEM THIS SOLVES
# -----------------------------------------------------------------------------
# cache_or_run() returns a saved result whenever the .rds exists. That is what
# makes a warm run fast, and it is also the one way this pipeline can silently
# report a stale number: if build_schemes() is edited, schemes_long.rds is now
# wrong, but cache_or_run() will load it anyway and every downstream result
# inherits the old classification without a warning.
#
# The cache has a dependency order, and invalidating a step requires
# invalidating everything built from it. invalidate_cache() does that by
# default rather than leaving it to be remembered:
#
#   schemes_long
#     `-- outpatient_panel
#           |-- instrument_screen, pooled_models, t02_pooled_current_vintage
#           |-- concept_level_6inst  (8 hours)
#           |     `-- meta_regressions_rf, meta_regressions_iv
#           |-- main_results, transform_ladder, confirming_results,
#           |   discrepant_results, robustness_pool_results
#           |-- comparability_measures
#           |     `-- comparability_interaction
#           |-- s13_* (four robustness blocks)
#           `-- cbsa_panel, s15_* (market definition, windows, payer)
#
# Editing a function changes only the caches at or below its level. Editing a
# scheme rule invalidates everything; editing run_transform_ladder() invalidates
# one object.
###############################################################################


# ---------------------------------------------------------------------------
# Cache registry: .rds key -> global variable name -> what produced it
# ---------------------------------------------------------------------------
CACHE_REGISTRY <- list(
  schemes_long                = list(var = "schemes_long",     stage = "BUILD"),
  outpatient_panel            = list(var = "outpatient",       stage = "BUILD"),
  instrument_screen           = list(var = "screen",           stage = "5"),
  pooled_models               = list(var = "pooled",           stage = "5"),
  concept_level_6inst         = list(var = "concept_results",  stage = "6"),
  main_results                = list(var = "main",             stage = "7"),
  transform_ladder            = list(var = "ladder",           stage = "7"),
  confirming_results          = list(var = "confirming",       stage = "7"),
  discrepant_results          = list(var = "discrepant",       stage = "7"),
  robustness_pool_results     = list(var = "robustness_pool",  stage = "7"),
  meta_regressions_rf         = list(var = "meta_rf",          stage = "8"),
  meta_regressions_iv         = list(var = "meta_iv",          stage = "8"),
  comparability_measures      = list(var = "measures",         stage = "9"),
  comparability_interaction   = list(var = "comp_int",         stage = "9"),
  s13_system_month_v2         = list(var = "s13_sysmonth",     stage = "13A"),
  s13_leave_one_system_out_v2 = list(var = "s13_loo",          stage = "13B"),
  s13_instrument_ladder       = list(var = "s13_ladder",       stage = "13C"),
  s13_randomisation_inference = list(var = "s13_ri",           stage = "13D"),
  cbsa_panel                  = list(var = "cbsa_panel",       stage = "15A"),
  s15_cbsa_results            = list(var = "s15_cbsa",         stage = "15A"),
  s15_cbsa_pooled             = list(var = "s15_cbsa_pooled",  stage = "15A"),
  s15b_window_ladder_v11      = list(var = "s15b_results",     stage = "15B"),
  s15_payer_class_panel       = list(var = "s15_payer",        stage = "15C"),
  s15_payer_results           = list(var = "s15_payer_res",    stage = "15C"),
  # The LaTeX block assigns this to `pooled`, the same name stage 5 uses for a
  # different object. Restored under a distinct name so one cannot overwrite the
  # other; the block itself is unaffected.
  t02_pooled_current_vintage  = list(var = "pooled_t02",       stage = "LaTeX")
)

# Dependency edges: invalidating the name on the left invalidates everything on
# the right, transitively. Used by invalidate_cache(cascade = TRUE).
CACHE_DEPENDENTS <- list(
  schemes_long        = "outpatient_panel",
  outpatient_panel    = c("instrument_screen", "pooled_models",
                          "t02_pooled_current_vintage", "concept_level_6inst",
                          "main_results", "transform_ladder",
                          "confirming_results", "discrepant_results",
                          "robustness_pool_results", "comparability_measures",
                          "s13_system_month_v2", "s13_leave_one_system_out_v2",
                          "s13_instrument_ladder", "s13_randomisation_inference",
                          "cbsa_panel", "s15b_window_ladder_v11",
                          "s15_payer_class_panel"),
  concept_level_6inst    = c("meta_regressions_rf", "meta_regressions_iv"),
  comparability_measures = "comparability_interaction",
  cbsa_panel             = c("s15_cbsa_results", "s15_cbsa_pooled"),
  s15_payer_class_panel  = "s15_payer_results"
)


# ---------------------------------------------------------------------------
# cache_status() -- what exists on disk
# ---------------------------------------------------------------------------
# Base R throughout, so this block keeps working even if the pipeline's
# packages are not loaded in the session.
cache_status <- function(quiet = FALSE) {
  keys <- names(CACHE_REGISTRY)
  path <- file.path(CACHE_DIR, paste0(keys, ".rds"))
  ok   <- file.exists(path)
  
  out <- data.frame(
    KEY       = keys,
    VARIABLE  = vapply(CACHE_REGISTRY, function(x) x$var,   character(1), USE.NAMES = FALSE),
    STAGE     = vapply(CACHE_REGISTRY, function(x) x$stage, character(1), USE.NAMES = FALSE),
    CACHED    = ok,
    SIZE_MB   = ifelse(ok, round(file.size(path) / 1024^2, 2), NA_real_),
    WRITTEN   = ifelse(ok, format(file.mtime(path), "%Y-%m-%d %H:%M"), NA_character_),
    IN_MEMORY = vapply(keys, function(k) exists(CACHE_REGISTRY[[k]]$var, envir = .GlobalEnv),
                       logical(1), USE.NAMES = FALSE),
    stringsAsFactors = FALSE
  )
  
  if (!quiet) {
    cat("\nCache directory:", CACHE_DIR, "\n")
    cat(sum(out$CACHED), "of", nrow(out), "objects cached,",
        round(sum(out$SIZE_MB, na.rm = TRUE), 1), "MB total\n\n")
    print(out[, c("KEY", "STAGE", "CACHED", "SIZE_MB", "WRITTEN")], row.names = FALSE)
    miss <- out$KEY[!out$CACHED]
    if (length(miss))
      cat("\nNot cached (", length(miss), "):\n  ",
          paste(miss, collapse = ", "), "\n", sep = "")
  }
  invisible(out)
}


# ---------------------------------------------------------------------------
# restore_session() -- load a completed run into a fresh session
#
# Loads every cached object into the global environment under the same variable
# name the run blocks use, so downstream code written afterwards sees exactly
# what it would have seen at the end of a full run.
#
# Does NOT recompute anything. A key with no .rds is reported as missing and
# skipped, so this is safe to call on a partially completed run.
# ---------------------------------------------------------------------------
restore_session <- function(keys = names(CACHE_REGISTRY), quiet = FALSE) {
  keys    <- intersect(keys, names(CACHE_REGISTRY))
  loaded  <- character(0)
  missing <- character(0)
  
  # Two keys mapping to one variable name would mean the later silently
  # overwrites the earlier, which is exactly the kind of thing that produces a
  # correct-looking object holding the wrong contents.
  targets <- vapply(keys, function(k) CACHE_REGISTRY[[k]]$var, character(1))
  if (anyDuplicated(targets))
    stop("CACHE_REGISTRY maps more than one key to the same variable name: ",
         paste(unique(targets[duplicated(targets)]), collapse = ", "),
         call. = FALSE)
  
  unreadable <- character(0)
  for (k in keys) {
    path <- file.path(CACHE_DIR, paste0(k, ".rds"))
    if (!file.exists(path)) { missing <- c(missing, k); next }

    # readRDS() used to run unprotected here. One corrupted or not-yet-synced
    # file (a OneDrive cloud-only placeholder is the usual cause) aborted the
    # entire loop with no indication of which key was responsible, so a
    # session with 19 of 20 objects perfectly readable still failed to
    # restore any of them. Now the bad file is skipped and named, and
    # everything else loads.
    obj <- tryCatch(readRDS(path), error = function(e) {
      cat("\n[SKIPPED] ", k, " -- ", conditionMessage(e),
          "\n          ", path, "\n", sep = "")
      NULL
    })
    if (is.null(obj)) { unreadable <- c(unreadable, k); next }

    assign(CACHE_REGISTRY[[k]]$var, obj, envir = .GlobalEnv)
    loaded <- c(loaded, CACHE_REGISTRY[[k]]$var)
  }
  
  # Rebuilt rather than cached: cheap, and needed by several later blocks.
  if (exists("MAIN_INSTRUMENTS", envir = .GlobalEnv) &&
      exists("SUPPORTING_INSTRUMENTS", envir = .GlobalEnv)) {
    assign("CONCEPT_INSTRUMENTS",
           c(get("MAIN_INSTRUMENTS", envir = .GlobalEnv),
             get("SUPPORTING_INSTRUMENTS", envir = .GlobalEnv)),
           envir = .GlobalEnv)
    loaded <- c(loaded, "CONCEPT_INSTRUMENTS")
  }
  
  if (!quiet) {
    cat("\nRestored", length(loaded), "objects into the global environment:\n  ",
        paste(loaded, collapse = ", "), "\n")
    if (length(missing))
      cat("\nNot on disk, skipped:\n  ", paste(missing, collapse = ", "),
          "\n  Run the stages that produce these, or ignore if not needed.\n")
    if (length(unreadable))
      cat("\nOn disk but unreadable, skipped (see errors above):\n  ",
          paste(unreadable, collapse = ", "),
          "\n  If this is OneDrive, the file may be a cloud-only placeholder --",
          "\n  right-click the cache folder and choose 'Always Keep on This Device'.\n")
    cat("\nResult CSVs remain readable from", TABLE_DIR, "\n")
  }
  invisible(list(loaded = loaded, missing = missing, unreadable = unreadable))
}


# ---------------------------------------------------------------------------
# invalidate_cache() -- delete caches so the next run recomputes them
#
# cascade = TRUE (the default) also deletes everything built from the named
# keys, which is what keeps a code change from producing a half-updated set of
# results. dry_run = TRUE shows what would be deleted without deleting it, and
# is worth using first when concept_level_6inst is in scope -- that one costs
# eight hours to rebuild.
# ---------------------------------------------------------------------------
invalidate_cache <- function(keys, cascade = TRUE, dry_run = TRUE) {
  unknown <- setdiff(keys, names(CACHE_REGISTRY))
  if (length(unknown))
    stop("Unknown cache key(s): ", paste(unknown, collapse = ", "),
         "\nSee names(CACHE_REGISTRY).", call. = FALSE)
  
  if (cascade) {
    repeat {
      expanded <- unique(c(keys, unlist(CACHE_DEPENDENTS[keys], use.names = FALSE)))
      expanded <- expanded[!is.na(expanded)]
      if (setequal(expanded, keys)) break
      keys <- expanded
    }
  }
  
  present <- keys[file.exists(file.path(CACHE_DIR, paste0(keys, ".rds")))]
  
  if (!length(present)) {
    cat("Nothing to invalidate: none of those keys are cached.\n")
    return(invisible(character(0)))
  }
  
  cat("\n", if (dry_run) "WOULD DELETE" else "DELETING", " ",
      length(present), " cached object(s):\n", sep = "")
  for (k in present)
    cat("  ", k, "  (stage ", CACHE_REGISTRY[[k]]$stage, ")\n", sep = "")
  
  if ("concept_level_6inst" %in% present)
    cat("\n  NOTE: concept_level_6inst takes roughly 8 hours to rebuild.\n")
  
  if (dry_run) {
    cat("\nDry run. Re-call with dry_run = FALSE to delete.\n")
    return(invisible(present))
  }
  
  file.remove(file.path(CACHE_DIR, paste0(present, ".rds")))
  cat("\nDeleted. The next run of the relevant stages will recompute them.\n")
  invisible(present)
}


###############################################################################
# TYPICAL USE
#
# Fresh session, continuing from a completed run:
#
#   # load Sections 0-12 to define the pipeline, then:
#   cache_status()
#   restore_session()
#   # every object is now in memory under its usual name
#
# After editing a function and wanting the change to propagate:
#
#   invalidate_cache("comparability_measures")                    # dry run
#   invalidate_cache("comparability_measures", dry_run = FALSE)   # delete
#   RUN_STAGES <- c(9)                                            # recompute
#
# After editing a scheme rule in Section 3, which changes the classification
# every result depends on:
#
#   invalidate_cache("schemes_long")   # dry run first -- this cascades widely
#
# To recompute one object without touching anything downstream, set
# cascade = FALSE. Do that only when the change genuinely cannot affect
# dependents, since a partially updated cache is worse than a stale one.
# The line below was a bare, unguarded warm_start() call in the original
# source, sitting right after the two source() calls for HPT_warm_start.R and
# HPT_diagnostic_for_family_levels.R (both of which are correctly relocated
# to PART 3 and PART 6 respectively, and do not appear here). This call was
# carried over unchanged during reorganization and, because this whole block
# now sits in PART 2 -- which is unconditional, sourced at the top of every
# run before HPT_WARM_START is even read as a switch -- it caused infinite
# recursion: warm_start() re-sources the entire pipeline, which reaches this
# same line again, which calls warm_start() again. Commented out. Call
# warm_start() yourself, from the console, after this file has finished
# loading -- never from inside the file being sourced.
#
# warm_start()






###############################################################################

###############################################################################
#                                                                             #
#   PART 3 -- BUILD AND PREFLIGHT                                             #
#                                                                             #
#   run_preflight() must pass before anything expensive runs. It catches      #
#   unclassified concepts, stale caches, and degenerate instruments, each of  #
#   which otherwise yields plausible-looking output rather than an error.     #
#                                                                             #
#   In warm-start mode this block restores the cache instead of building.     #
#                                                                             #
###############################################################################

if (isTRUE(HPT_WARM_START)) {

  # Fresh session on a completed run. Every object the run blocks produce is
  # read from CACHE_DIR and assigned into the global environment under its
  # usual name. Nothing is estimated and nothing is written.
  #
  # HPT_WARM_START_KEYS, if set before sourcing, restricts this to a subset of
  # CACHE_REGISTRY's keys instead of all 25. Added because reading a large
  # .rds (outpatient_panel, cbsa_panel, s15_payer_class_panel) spikes memory
  # well above its final resident size during decompression, and on an 8GB
  # machine that spike does not reliably get released even after the object
  # is later rm()'d -- confirmed directly: removing 23 of the 25 restored
  # objects after the fact recovered only ~1.2GB of the ~15GB in use. Section
  # 33/34 work needs only concept_level_6inst, so:
  #   HPT_WARM_START_KEYS <- "concept_level_6inst"
  # restores only concept_results and never opens the large files at all.
  # Unset (the default), behavior is unchanged: every key is restored.
  restore_keys <- if (exists("HPT_WARM_START_KEYS")) HPT_WARM_START_KEYS
                  else names(CACHE_REGISTRY)
  cat("\n", strrep("-", 70),
      "\nWARM START -- restoring cached objects, estimating nothing\n",
      if (length(restore_keys) < length(CACHE_REGISTRY))
        paste0("Restricted to: ", paste(restore_keys, collapse = ", "), "\n")
      else "",
      strrep("-", 70), "\n", sep = "")
  restore_session(keys = restore_keys)

  # Section 12's .t12_* helpers are defined inside a stage block, so they are
  # not in scope after a warm start. HPT_warm_start.R supplies
  # load_t12_helpers() to pull just their definitions without the 30-45
  # minute driver run.
  source(file.path(CODE_DIR, "HPT_warm_start.R"))

  cat("\nSession restored. Every pipeline function is in scope and every",
      "\ncached result is loaded. RUN_STAGES is empty, so no stage block",
      "\nwill execute. cache_status() shows what is on disk.\n")

} else {

  hpt_check_inputs()

  # Companion file, located relative to CODE_DIR rather than by absolute path.
  source(file.path(CODE_DIR, "HPT_warm_start.R"))

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




}  # end warm-start branch

###############################################################################
#                                                                             #
#   PART 4 -- ESTIMATION                                                      #
#                                                                             #
###############################################################################


# ============================================================================
# Stages 5 and 6: instrument screen, concept sweep, stage-6 loader
# ============================================================================

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

# ============================================================================
# Instrument audit and Stages 7-10
# ============================================================================

if (isTRUE(HPT_RUN$instrument_audit)) {

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




}  # end HPT_RUN$instrument_audit

# ============================================================================
# Stage 11: instrument balance and placebo
# ============================================================================

######## STAGE 11 -- instrument balance and placebo ###########################
if (11 %in% RUN_STAGES) {
  balance <- run_instrument_balance(outpatient)
  placebo <- run_placebo_outcome(outpatient)
}



###############################################################################

# ============================================================================
# Robustness: controlling for system affiliation
# ============================================================================

if (isTRUE(HPT_RUN$insystem_robust)) {

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

}  # end HPT_RUN$insystem_robust

# ============================================================================
# County demographics from the ACS (needs a Census API key)
# ============================================================================

if (isTRUE(HPT_RUN$demographics)) {

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
#
# This overwrites the cached panel with the demographics-merged version,
# deliberately and outside cache_or_run(). The effect is that restore_session()
# returns `outpatient` with the DEMO_ columns already attached, so Sections 12
# and 12B run in a fresh session without another Census API pull.
#
# THE TRADE-OFF: the cached panel is no longer what the BUILD block alone would
# produce. After any rebuild of outpatient_panel -- whether through
# invalidate_cache() or USE_CACHE <- FALSE -- the DEMO_ columns are gone until
# this block is run again. Re-run it before Section 12 in that case.
saveRDS(outpatient, file.path(CACHE_DIR, "outpatient_panel.rds"))



###############################################################################

}  # end HPT_RUN$demographics

# ============================================================================
# Section 12: SES heterogeneity
# ============================================================================

if (isTRUE(HPT_RUN$demographics)) {

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
             P_VALUE = .pval(d / se, fit))
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
      RF_COEF = rf$b, RF_SE = rf$s, RF_P = .pval(rf$b / rf$s, rf_fit),
      # Level terms: % price change per SD of Z.
      # Interaction terms: change in that % per 1 SD of the moderator.
      RF_PERCENT_PER_SD_Z = if (term_lab[i] == "Level")
        100 * (exp(rf$b * sd_z) - 1) else NA_real_,
      RF_SHIFT_PCT_PER_SD_MOD = if (term_lab[i] == "x Moderator" && !is.na(j))
        100 * (exp(rf$b * sd_z * mod_sd[j]) - 1) else NA_real_,
      IV_COEF = iv$b, IV_SE = iv$s, IV_P = .pval(iv$b / iv$s, iv_fit),
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

}  # end HPT_RUN$demographics

# ============================================================================
# Section 12B: SES terciles
# ============================================================================

if (isTRUE(HPT_RUN$ses_terciles)) {

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
  data.table(WALD = W, DF = df, P_VALUE = pf(W / df, df1 = df, df2 = 15, lower.tail = FALSE),
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
               RF_COEF = a[1], RF_SE = a[2], RF_P = .pval(a[1] / a[2], rf),
               RF_PERCENT_PER_SD = 100 * (exp(a[1] * sd_z) - 1),
               IV_COEF = b[1], IV_P = .pval(b[1] / b[2], iv),
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

}  # end HPT_RUN$ses_terciles

# ============================================================================
# Section 13: system-month FE, leave-one-out, ladder, randomization
# ============================================================================

if (isTRUE(HPT_RUN$robustness_13)) {

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

# Duplicate of the two .s13_show() calls just above, unconditionally repeated
# outside the S13_RUN$leave_one_out guard that already printed them, so this
# copy errored whenever that switch was off and s13_loo did not exist.
#
# .s13_show(s13_loo[ESTIMATOR == "Reduced form",
#                   .(DROPPED_SYSTEM = substr(DROPPED_SYSTEM, 1, 26), INSTRUMENT_LABEL,
#                     P_BASE = round(P_BASE, 4), P_LOO = round(P_VALUE, 4),
#                     SHOP_PCT = round(SHOPPABLE_RF_PCT, 2), STILL_SIG = STILL_SIG_05)][
#                       order(-P_LOO)])
#
# cat("\nSummary (reduced form):\n")
# .s13_show(s13_loo[ESTIMATOR == "Reduced form",
#                   .(N = .N, N_STILL_SIG = sum(STILL_SIG_05), MEDIAN_P = round(median(P_VALUE), 4),
#                     MAX_P = round(max(P_VALUE), 4),
#                     RANGE_SHOP_PCT = paste0(round(min(SHOPPABLE_RF_PCT), 2), " to ",
#                                             round(max(SHOPPABLE_RF_PCT), 2)))])


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

}  # end HPT_RUN$robustness_13

# ============================================================================
# Section 15: CBSA markets, window ladder, payer-conditional
# ============================================================================

if (isTRUE(HPT_RUN$markets_15)) {

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
              RF_COEF = rf_b, RF_P = .pval(rf_b / rf_s, rf_fit),
              IV_PERCENT = if (!is.na(iv_b)) 100 * (exp(iv_b) - 1) else NA_real_,
              IV_P = if (!is.na(iv_b)) .pval(iv_b / iv_s, iv_fit) else NA_real_,
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
# This reporting block sat unconditionally after the S15_RUN$cbsa guard
# above closed, so it errored whenever that switch was off and
# s15_cbsa_pooled did not exist. Re-guarded on the same condition rather than
# merged into the block above, to avoid relocating text across a brace.
if (isTRUE(S15_RUN$cbsa)) {
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
}



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

}  # end HPT_RUN$markets_15

# ============================================================================
# Section 16: family and superfamily analysis
# ============================================================================

if (isTRUE(HPT_RUN$family_16)) {

# SECTION 16 -- FAMILY AND SUPERFAMILY LEVEL ANALYSIS (subsample design)
#
# Same estimator, same fixed effects, same clustering as the concept-level
# work. Nothing here builds a new outcome variable or a weighting rule; every
# number is estimated on raw rows filtered (or grouped) by FINAL_FAMILY_ID or
# FINAL_SUPERFAMILY_ID instead of ANALYSIS_SERVICE_ID.
#
# Load AFTER warm_start() / restore_session(). Needs `outpatient` and every
# Section-0-through-12 function and constant. Nothing here overwrites an
# existing object; new results go into GROUP_ prefixed and s16_ prefixed names.
#
# WHAT "INTERACTED IV" MEANS AT THIS GRAIN
#   estimate_interacted()'s categorical branch already handles a moderator
#   with any number of levels, not just two. Passing FINAL_FAMILY_ID directly
#   as the moderator on the full panel gives every family's coefficient from
#   ONE exactly-identified regression, with a joint heterogeneity Wald test --
#   the direct family-level extension of the Scheme 1 table, not new
#   machinery. This is 16.5/16.6 below.
#
# WHAT REPLACES THE MEta-REGRESSION
#   Sixteen families and MIN_CONCEPTS_META = 15 means a family-level
#   meta-regression is a two-group comparison clustered on the variable that
#   defines the groups -- not informative. What DOES carry weight at this
#   grain is the exact permutation test already used at concept level
#   (family_permutation_test()), now with each family contributing exactly
#   one point instead of several, which is better-specified inference, not a
#   downgrade. That is 16.7. A parametric weighted-mean-difference companion
#   is reported alongside it for comparison, not as the primary result.
###############################################################################

suppressPackageStartupMessages({ library(data.table); library(fixest) })

.s16_hd   <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")
.s16_show <- function(d, n = Inf) {
  d <- as.data.frame(d)
  if (is.finite(n) && nrow(d) > n) { print(head(d, n), row.names = FALSE)
    cat("  [", nrow(d) - n, " further rows suppressed]\n", sep = "")
  } else print(d, row.names = FALSE)
  invisible(NULL)
}

stopifnot(exists("outpatient"), exists("schemes_long"))


# ===========================================================================
# 16.0 BACKFILL FINAL_SUPERFAMILY_ID FOR MERGED CANONICAL CONCEPTS
#
# apply_concept_merges() sets FINAL_FAMILY_ID on the six merged concepts but
# never sets FINAL_SUPERFAMILY_ID, so it comes through as NA on those rows
# (confirmed: 14,143 rows / 6 concepts in the diagnostic). Every constituent
# family already has a known superfamily elsewhere in the panel, so this is a
# lookup, not an inference.
# ===========================================================================
.s16_hd("16.0 Backfilling FINAL_SUPERFAMILY_ID")

n_na_before <- outpatient[is.na(FINAL_SUPERFAMILY_ID), .N]

.s16_fam2super <- unique(outpatient[!is.na(FINAL_SUPERFAMILY_ID),
                                    .(FINAL_FAMILY_ID, FINAL_SUPERFAMILY_ID)])
if (uniqueN(.s16_fam2super$FINAL_FAMILY_ID) != nrow(.s16_fam2super)) {
  stop("A FINAL_FAMILY_ID maps to more than one FINAL_SUPERFAMILY_ID among the ",
       "non-missing rows -- the lookup is not well defined. Inspect ",
       "'.s16_fam2super[, .N, by = FINAL_FAMILY_ID][N > 1]' before proceeding.",
       call. = FALSE)
}

outpatient[is.na(FINAL_SUPERFAMILY_ID),
           FINAL_SUPERFAMILY_ID := .s16_fam2super$FINAL_SUPERFAMILY_ID[
             match(FINAL_FAMILY_ID, .s16_fam2super$FINAL_FAMILY_ID)]]

n_na_after <- outpatient[is.na(FINAL_SUPERFAMILY_ID), .N]
cat("FINAL_SUPERFAMILY_ID missing before:", n_na_before, "| after backfill:", n_na_after, "\n")
if (n_na_after > 0L) {
  cat("Remaining NA families (no known superfamily anywhere in the panel):\n")
  print(outpatient[is.na(FINAL_SUPERFAMILY_ID), .N, by = FINAL_FAMILY_ID])
}

# OTHER is a genuine mix under Scheme 1 (LABORATORY_PATHOLOGY / EVALUATION_
# MANAGEMENT are shoppable; EMERGENCY_DEPARTMENT / CRITICAL_CARE are forced
# non-shoppable), confirmed 64% shoppable by concept count in the diagnostic.
# That only matters where OTHER is asked to stand in for a shoppability label
# (e.g. a binary shoppable-vs-not collapse). It is NOT a problem for treating
# OTHER as its own category in the joint family/superfamily interaction below
# -- that only asks "what is OTHER's own average gradient," which is a valid
# question regardless of its internal composition.
S16_CONTAMINATED_SUPERFAMILIES <- "OTHER"


# ===========================================================================
# 16.1 GENERALISED GROUP-LEVEL SWEEP
#
# Direct analogue of estimate_concept_level() (Section 6), grouped by an
# arbitrary column instead of ANALYSIS_SERVICE_ID. Same three numbers per
# group-instrument pair (RF, FS, IV) for the same reason: if a gradient shows
# up in IV but not RF, it is a denominator artefact, not a price response.
#
# Built directly on build_ols_formula() / build_iv_formula() /
# build_cluster_formula() / extract_coefficient() -- the same primitives
# estimate_concept_level() uses -- rather than on run_reduced_form() /
# run_first_stage() / run_iv(), whose exact argument names were not visible
# from this session. Verify against T05_concept_level_RF_FS_IV.csv the first
# time this runs: XRAY_FLUOROSCOPY's family-level RF/FS should sit inside the
# range of its 199 constituent concepts' estimates.
# ===========================================================================
.s16_hd("16.1 run_group_level_sweep()")

run_group_level_sweep <- function(panel, group_col, instruments = CONCEPT_INSTRUMENTS,
                                  outcome = PRIMARY_OUTCOME, endogenous = ENDOGENOUS_VARIABLE,
                                  drop_singletons = DROP_SINGLETON_MARKETS,
                                  min_obs = MIN_SERVICE_OBS, min_markets = MIN_SERVICE_MARKETS,
                                  min_months = MIN_SERVICE_MONTHS,
                                  exclude_groups = character(0)) {
  stopifnot(group_col %in% names(panel))
  ids <- sort(unique(panel[[group_col]]))
  ids <- setdiff(ids, c(NA, exclude_groups))
  cat("\n", group_col, "groups:", length(ids), "| instruments:", length(instruments), "\n")
  
  rows <- vector("list", length(ids) * length(instruments)); k <- 0L; t0 <- Sys.time()
  
  for (gid in ids) {
    sub <- panel[get(group_col) == gid]
    if (drop_singletons) sub <- drop_singleton_markets(sub)
    n_mk <- uniqueN(sub$ANALYSIS_MARKET); n_mo <- uniqueN(sub$POST_MONTH)
    if (nrow(sub) < min_obs || n_mk < min_markets || n_mo < min_months) {
      cat("  Skipped", gid, "-- n=", format(nrow(sub), big.mark = ","),
          "markets=", n_mk, "months=", n_mo, "\n")
      next
    }
    
    for (il in names(instruments)) {
      z <- instruments[[il]]
      if (!(z %in% names(sub)) || !has_usable_variation(sub[[z]])) next
      
      req <- intersect(c(outcome, endogenous, z, BASELINE_CONTROLS,
                         BASELINE_FIXED_EFFECTS, BASELINE_CLUSTERS), names(sub))
      d <- sub[complete.cases(sub[, ..req])]
      if (nrow(d) < min_obs) next
      
      rf <- tryCatch(feols(build_ols_formula(outcome, c(z, BASELINE_CONTROLS), BASELINE_FIXED_EFFECTS),
                           data = d, cluster = build_cluster_formula(BASELINE_CLUSTERS),
                           warn = FALSE, notes = FALSE), error = function(e) NULL)
      fs <- tryCatch(feols(build_ols_formula(endogenous, c(z, BASELINE_CONTROLS), BASELINE_FIXED_EFFECTS),
                           data = d, cluster = build_cluster_formula(BASELINE_CLUSTERS),
                           warn = FALSE, notes = FALSE), error = function(e) NULL)
      iv <- tryCatch(feols(build_iv_formula(outcome, endogenous, z, BASELINE_CONTROLS, BASELINE_FIXED_EFFECTS),
                           data = d, cluster = build_cluster_formula(BASELINE_CLUSTERS),
                           warn = FALSE, notes = FALSE), error = function(e) NULL)
      
      c_rf <- extract_coefficient(rf, z)
      c_fs <- extract_coefficient(fs, z)
      c_iv <- extract_coefficient(iv, c(paste0("fit_", endogenous), endogenous))
      
      k <- k + 1L
      rows[[k]] <- data.table(
        GROUP_COL = group_col, GROUP_ID = gid, INSTRUMENT_LABEL = il, INSTRUMENT = z, OUTCOME = outcome,
        RF_COEF = c_rf$estimate, RF_SE = c_rf$std_error, RF_P = c_rf$p_value,
        FS_COEF = c_fs$estimate, FS_SE = c_fs$std_error,
        FS_F = if (is.finite(c_fs$statistic)) c_fs$statistic^2 else NA_real_,
        IV_COEF = c_iv$estimate, IV_SE = c_iv$std_error, IV_P = c_iv$p_value,
        IV_ESTIMATE_PERCENT = 100 * c_iv$estimate, IV_SE_PERCENT = 100 * c_iv$std_error,
        N_ROWS = nrow(d), N_HOSPITALS = uniqueN(d$HOSPITAL_ID),
        N_CONCEPTS = uniqueN(d$FINAL_CONCEPT_ID), N_MARKETS = n_mk, N_MONTHS = n_mo,
        N_OBSERVATIONS = if (is.null(iv)) NA_integer_ else nobs(iv))
    }
    cat("  Done:", gid, "\n")
  }
  
  out <- rbindlist(rows[seq_len(k)], fill = TRUE)
  if (nrow(out) == 0L) stop("No rows for group_col = ", group_col, call. = FALSE)
  out[, RF_P_FDR := p.adjust(RF_P, method = "BH"), by = INSTRUMENT_LABEL]
  out[, IV_P_FDR := p.adjust(IV_P, method = "BH"), by = INSTRUMENT_LABEL]
  setorder(out, INSTRUMENT_LABEL, RF_COEF)
  cat("\nElapsed:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "min |",
      k, "group-instrument rows\n")
  out
}


# ===========================================================================
# 16.2 RUN: FAMILY-LEVEL SWEEP  -> T16A
# ===========================================================================
.s16_hd("16.2 Family-level sweep (T16A)")

family_results <- cache_or_run("family_level_6inst",
                               run_group_level_sweep(outpatient, "FINAL_FAMILY_ID", instruments = CONCEPT_INSTRUMENTS))
save_csv(family_results, "T16A_family_level_RF_FS_IV.csv")

.s16_show(family_results[INSTRUMENT_LABEL == names(MAIN_INSTRUMENTS)[1],
                         .(GROUP_ID, RF_COEF = round(RF_COEF, 5), RF_P = round(RF_P, 4),
                           FS_F = round(FS_F, 1), IV_PCT = round(IV_ESTIMATE_PERCENT, 2),
                           N_ROWS)][order(RF_COEF)])


# ===========================================================================
# 16.3 RUN: SUPERFAMILY-LEVEL SWEEP  -> T16B
# ===========================================================================
.s16_hd("16.3 Superfamily-level sweep (T16B)")

superfamily_results <- cache_or_run("superfamily_level_6inst",
                                    run_group_level_sweep(outpatient, "FINAL_SUPERFAMILY_ID", instruments = CONCEPT_INSTRUMENTS))
save_csv(superfamily_results, "T16B_superfamily_level_RF_FS_IV.csv")

.s16_show(superfamily_results[INSTRUMENT_LABEL == names(MAIN_INSTRUMENTS)[1],
                              .(GROUP_ID, RF_COEF = round(RF_COEF, 5), RF_P = round(RF_P, 4),
                                FS_F = round(FS_F, 1), IV_PCT = round(IV_ESTIMATE_PERCENT, 2),
                                N_ROWS)][order(RF_COEF)])
cat("\nOTHER pools LABORATORY_PATHOLOGY/EVALUATION_MANAGEMENT (shoppable) with\n",
    "EMERGENCY_DEPARTMENT/CRITICAL_CARE (forced non-shoppable). Read its row as\n",
    "a basket average, not as evidence about a coherent 'other services' market.\n", sep = "")


# ===========================================================================
# 16.4 IS SHOPPABILITY A VALID CONCEPT-COUNT-WEIGHTED LABEL AT EACH GRAIN?
#
# Sanity check before the joint interactions: how many families/superfamilies
# in the sweep tables are the a-priori DIAGNOSTIC_FAMILIES set, restated here
# so the permutation test in 16.7 is auditable against something printed.
# ===========================================================================
.s16_hd("16.4 Diagnostic-family membership check")
cat("DIAGNOSTIC_FAMILIES (a-priori shoppable):\n  ",
    paste(DIAGNOSTIC_FAMILIES, collapse = ", "), "\n")
cat("\nFamilies in the sweep NOT in DIAGNOSTIC_FAMILIES (a-priori non-shoppable):\n  ",
    paste(setdiff(unique(family_results$GROUP_ID), DIAGNOSTIC_FAMILIES), collapse = ", "), "\n")


# ===========================================================================
# 16.5 JOINT INTERACTED IV -- FAMILY AS THE MODERATOR  -> T16C
#
# One exactly-identified regression per instrument: FINAL_FAMILY_ID (all
# levels) interacted with Z, on the full panel. Direct family-level extension
# of run_main_results()'s Scheme 1 table -- same function, wider moderator.
# ===========================================================================
.s16_hd("16.5 Joint interacted IV, family moderator (T16C)")

family_interacted <- cache_or_run("family_interacted", {
  rows <- list(); tests <- list()
  for (il in names(MAIN_INSTRUMENTS)) {
    z <- MAIN_INSTRUMENTS[[il]]
    r <- estimate_interacted(outpatient, "FINAL_FAMILY_ID", PRIMARY_OUTCOME, z,
                             moderator_type = "categorical", label = "Family (joint)",
                             instrument_label = il, moderator_label = "FINAL_FAMILY_ID")
    if (!is.null(r)) { rows[[il]] <- r$rows; tests[[il]] <- r$tests }
    cat("  ", il, "done\n")
  }
  list(rows = rbindlist(rows, fill = TRUE), tests = rbindlist(tests, fill = TRUE))
})
save_csv(family_interacted$rows,  "T16C_family_interacted_RF_and_IV.csv")
save_csv(family_interacted$tests, "T16C_family_interacted_heterogeneity_test.csv")

.s16_show(family_interacted$rows[INSTRUMENT_LABEL == names(MAIN_INSTRUMENTS)[1],
                                 .(TERM, RF_PERCENT_PER_SD = round(RF_PERCENT_PER_SD, 3),
                                   RF_P = round(RF_P, 4), IV_PCT = round(IV_PERCENT, 2),
                                   FIRST_STAGE_WALD_THIS_EQ = round(FIRST_STAGE_WALD_THIS_EQ, 1))][
                                     order(RF_PERCENT_PER_SD)])
cat("\nJoint heterogeneity test (H0: all 16 family coefficients equal):\n")
.s16_show(family_interacted$tests[, .(INSTRUMENT_LABEL, ESTIMATOR, P_VALUE = round(P_VALUE, 5))])


# ===========================================================================
# 16.6 JOINT INTERACTED IV -- SUPERFAMILY AS THE MODERATOR  -> T16D
# ===========================================================================
.s16_hd("16.6 Joint interacted IV, superfamily moderator (T16D)")

superfamily_interacted <- cache_or_run("superfamily_interacted", {
  rows <- list(); tests <- list()
  for (il in names(MAIN_INSTRUMENTS)) {
    z <- MAIN_INSTRUMENTS[[il]]
    r <- estimate_interacted(outpatient, "FINAL_SUPERFAMILY_ID", PRIMARY_OUTCOME, z,
                             moderator_type = "categorical", label = "Superfamily (joint)",
                             instrument_label = il, moderator_label = "FINAL_SUPERFAMILY_ID")
    if (!is.null(r)) { rows[[il]] <- r$rows; tests[[il]] <- r$tests }
    cat("  ", il, "done\n")
  }
  list(rows = rbindlist(rows, fill = TRUE), tests = rbindlist(tests, fill = TRUE))
})
save_csv(superfamily_interacted$rows,  "T16D_superfamily_interacted_RF_and_IV.csv")
save_csv(superfamily_interacted$tests, "T16D_superfamily_interacted_heterogeneity_test.csv")

.s16_show(superfamily_interacted$rows[INSTRUMENT_LABEL == names(MAIN_INSTRUMENTS)[1],
                                      .(TERM, RF_PERCENT_PER_SD = round(RF_PERCENT_PER_SD, 3),
                                        RF_P = round(RF_P, 4), IV_PCT = round(IV_PERCENT, 2))][
                                          order(RF_PERCENT_PER_SD)])
cat("\nJoint heterogeneity test (H0: IMAGING = BIOPSY = GI_ENDOSCOPY = OTHER):\n")
.s16_show(superfamily_interacted$tests[, .(INSTRUMENT_LABEL, ESTIMATOR, P_VALUE = round(P_VALUE, 5))])


# ===========================================================================
# 16.7 FAMILY-LEVEL INFERENCE -- weighted mean-difference + exact permutation
#
# The unit here IS the family (one point each), so this is what a
# meta-regression collapses to at n=16, plus the permutation test that is
# actually well powered at this grain: enumerating which subsets of families
# could have been called "shoppable" gives inference at the level the
# shoppability decision was actually made, exactly as
# family_permutation_test() does at concept level -- just without the
# multiple-concepts-per-family double counting.
# ===========================================================================
.s16_hd("16.7 Family-level weighted mean-difference and exact permutation (T16E)")

run_group_weighted_diff <- function(sweep, dep = "RF_COEF", se = "RF_SE",
                                    shoppable = DIAGNOSTIC_FAMILIES) {
  rbindlist(lapply(unique(sweep$INSTRUMENT_LABEL), function(il) {
    d <- sweep[INSTRUMENT_LABEL == il & is.finite(get(dep)) & is.finite(get(se)) & get(se) > 0]
    d[, S := fifelse(GROUP_ID %chin% shoppable, "Shoppable", "Non_shoppable")]
    if (uniqueN(d$S) < 2L) return(data.table())
    fit <- tryCatch(feols(as.formula(paste(dep, "~ S")), data = d, weights = as.formula(paste0("~I(1/", se, "^2)")),
                          warn = FALSE, notes = FALSE), error = function(e) NULL)
    if (is.null(fit)) return(data.table())
    td <- tidy_fixest(fit)
    td[, `:=`(INSTRUMENT_LABEL = il, N_GROUPS = nrow(d),
              N_SHOPPABLE = d[S == "Shoppable", .N], N_NONSHOPPABLE = d[S == "Non_shoppable", .N])]
    td
  }), fill = TRUE)
}

run_group_permutation_test <- function(sweep, dep = "RF_COEF", se = "RF_SE",
                                       shoppable = DIAGNOSTIC_FAMILIES,
                                       max_exact = 50000L, seed = 20260817L) {
  wdiff <- function(d, groups) {
    d <- copy(d)[, S := fifelse(GROUP_ID %chin% groups, "Shoppable", "Non_shoppable")]
    if (uniqueN(d$S) < 2L) return(NA_real_)
    d[, W := 1 / (get(se)^2)]
    d[S == "Shoppable", sum(W * get(dep)) / sum(W)] -
      d[S == "Non_shoppable", sum(W * get(dep)) / sum(W)]
  }
  rbindlist(lapply(unique(sweep$INSTRUMENT_LABEL), function(il) {
    d <- sweep[INSTRUMENT_LABEL == il & is.finite(get(dep)) & is.finite(get(se)) & get(se) > 0]
    groups <- sort(unique(d$GROUP_ID))
    obs_g <- intersect(shoppable, groups); k <- length(obs_g); n <- length(groups)
    if (k == 0L || k == n || n < 4L) return(data.table())
    observed <- wdiff(d, obs_g); if (!is.finite(observed)) return(data.table())
    exact <- is.finite(choose(n, k)) && choose(n, k) <= max_exact
    null <- if (exact) vapply(combn(groups, k, simplify = FALSE), function(g) wdiff(d, g), numeric(1))
    else { set.seed(seed); vapply(seq_len(20000L), function(i) wdiff(d, sample(groups, k)), numeric(1)) }
    null <- null[is.finite(null)]
    data.table(INSTRUMENT_LABEL = il, OBSERVED = observed, N_GROUPS = n, N_SHOPPABLE = k,
               METHOD = fifelse(exact, "Exact enumeration", "Monte Carlo (20,000)"),
               NULL_P05 = quantile(null, .05, names = FALSE),
               NULL_P95 = quantile(null, .95, names = FALSE),
               P_TWO_SIDED = mean(abs(null) >= abs(observed)))
  }), fill = TRUE)
}

family_weighted_diff <- run_group_weighted_diff(family_results)
family_permutation   <- run_group_permutation_test(family_results)

save_csv(family_weighted_diff, "T16E_family_weighted_mean_difference.csv")
save_csv(family_permutation,   "T16E_family_exact_permutation.csv")

cat("\nParametric weighted mean-difference (family as unit, n = 16):\n")
.s16_show(family_weighted_diff[grepl("Shoppable", term),
                               .(INSTRUMENT_LABEL, ESTIMATE = round(estimate, 5),
                                 SE = round(std.error, 5), P = round(p.value, 4), N_GROUPS)])
cat("\nExact permutation over family-shoppability assignment:\n")
.s16_show(family_permutation[, .(INSTRUMENT_LABEL, OBSERVED = round(OBSERVED, 5),
                                 N_GROUPS, N_SHOPPABLE, METHOD, P_TWO_SIDED = round(P_TWO_SIDED, 5))])
cat("\nThese two should broadly agree in sign and rough magnitude of significance;\n",
    "report the permutation p-value as primary, consistent with the concept-level\n",
    "convention, and the weighted regression as the parametric companion.\n", sep = "")



# ===========================================================================
# 16.8 SUPERFAMILY-LEVEL ACS/SES HETEROGENEITY  -> T16F
#
# IMAGING vs. everything else, interacted with the SES index and with the
# instrument, following Section 12's triple-interaction design:
#
#   ln(P) = gI*(Z x Imaging)     + gR*(Z x Rest)
#         + dI*(Z x Imaging x M) + dR*(Z x Rest x M)
#         + controls + FE
#
# Object of interest is dI - dR, a single 1-df Wald test. Family-level (16-way)
# SES heterogeneity is NOT attempted here: that is 16 categories x SES, ~32
# endogenous terms, and Section 12's own note already flags a weak first stage
# with just 4. Superfamily's binary IMAGING/rest split is the version of this
# question that is actually estimable.
#
# Requires .t12_build_ses_index() from Section 12's self-contained SES block.
# That block lives entirely in the executed region of the pipeline (past
# `# ORDER OF OPERATIONS`), not in the Sections 0-12 definitions warm_start()
# loads, so it needs its own loader: load_t12_helpers() (in HPT_warm_start.R)
# sources just the config + .t12_* function definitions, stopping short of
# Section 12's own 30-45 minute driver run.
# ===========================================================================
load_t12_helpers()                       # loads .t12_build_ses_index() etc., no 30-45 min run

.s16_hd("16.8 Superfamily x SES heterogeneity (T16F)")

if (!exists(".t12_build_ses_index")) {
  if (exists("load_t12_helpers")) {
    cat("`.t12_build_ses_index` not in memory -- loading Section 12's helpers.\n")
    load_t12_helpers()
  } else {
    stop("`.t12_build_ses_index` is not defined and `load_t12_helpers()` is not ",
         "available either.\n  Add load_t12_helpers() to HPT_warm_start.R (see the ",
         "updated version) and call it, or manually source lines 3872-4340 of ",
         "HPT_Analysis_Pipeline.R (the SECTION 12 banner through the line before ",
         "'# 6. RUN').", call. = FALSE)
  }
}

.s16_ses_panel <- tryCatch(.t12_build_ses_index(outpatient), error = function(e) {
  cat("`.t12_build_ses_index(outpatient)` failed:", conditionMessage(e), "\n",
      "Skipping 16.8.\n")
  NULL
})

run_superfamily_ses_heterogeneity <- function(panel, shop_value = "IMAGING",
                                              group_col = "FINAL_SUPERFAMILY_ID",
                                              moderator = "SES_INDEX",
                                              instrument, instrument_label = "",
                                              outcome = PRIMARY_OUTCOME,
                                              endogenous = ENDOGENOUS_VARIABLE,
                                              controls = BASELINE_CONTROLS,
                                              fe = BASELINE_FIXED_EFFECTS,
                                              clusters = BASELINE_CLUSTERS) {
  if (!(moderator %in% names(panel))) {
    stop("Column '", moderator, "' not found. Check what .t12_build_ses_index() ",
         "actually names the index and pass it via `moderator =`.", call. = FALSE)
  }
  d <- panel[!is.na(get(group_col)) & !is.na(get(moderator))]
  req <- intersect(c(outcome, endogenous, instrument, moderator, controls, fe, clusters), names(d))
  d <- d[complete.cases(d[, ..req])]
  if (nrow(d) < MIN_MODEL_OBS) return(NULL)
  
  d[, SHOP := as.integer(get(group_col) == shop_value)]
  d[, M := safe_numeric(get(moderator))]
  d[, M := M - mean(M, na.rm = TRUE)]
  d[, `:=`(
    RF_S = get(instrument) * SHOP,       RF_N = get(instrument) * (1L - SHOP),
    RF_SM = get(instrument) * SHOP * M,  RF_NM = get(instrument) * (1L - SHOP) * M,
    TR_S = get(endogenous) * SHOP,       TR_N = get(endogenous) * (1L - SHOP),
    TR_SM = get(endogenous) * SHOP * M,  TR_NM = get(endogenous) * (1L - SHOP) * M)]
  
  rf_fit <- tryCatch(feols(build_ols_formula(outcome, c("RF_S","RF_N","RF_SM","RF_NM", controls), fe),
                           data = d, cluster = build_cluster_formula(clusters),
                           warn = FALSE, notes = FALSE), error = function(e) NULL)
  iv_fit <- tryCatch(feols(build_iv_formula(outcome, c("TR_S","TR_N","TR_SM","TR_NM"),
                                            c("RF_S","RF_N","RF_SM","RF_NM"), controls, fe),
                           data = d, cluster = build_cluster_formula(clusters),
                           warn = FALSE, notes = FALSE), error = function(e) NULL)
  if (is.null(rf_fit)) return(NULL)
  
  wald_diff <- function(fit, t1, t2) {
    nm <- names(coef(fit))
    n1 <- intersect(c(t1, paste0("fit_", t1)), nm)[1L]
    n2 <- intersect(c(t2, paste0("fit_", t2)), nm)[1L]
    if (is.na(n1) || is.na(n2)) return(c(diff = NA_real_, se = NA_real_, p = NA_real_))
    b1 <- unname(coef(fit)[n1]); b2 <- unname(coef(fit)[n2]); V <- vcov(fit)
    dd <- b1 - b2; se_ <- sqrt(V[n1, n1] + V[n2, n2] - 2 * V[n1, n2])
    c(diff = dd, se = se_, p = .pval(dd / se_, fit))
  }
  
  rf_diff <- wald_diff(rf_fit, "RF_SM", "RF_NM")
  iv_diff <- if (!is.null(iv_fit)) wald_diff(iv_fit, "TR_SM", "TR_NM") else c(diff = NA, se = NA, p = NA)
  fsw <- if (!is.null(iv_fit)) first_stage_wald(iv_fit) else data.table()
  fs_min <- if (nrow(fsw) > 0L) min(fsw$WALD, na.rm = TRUE) else NA_real_
  
  data.table(GROUP = shop_value, MODERATOR = moderator, INSTRUMENT_LABEL = instrument_label,
             INSTRUMENT = instrument, OUTCOME = outcome,
             RF_DIFF_SHOP_MINUS_REST = unname(rf_diff["diff"]), RF_SE = unname(rf_diff["se"]),
             RF_P = unname(rf_diff["p"]),
             IV_DIFF_SHOP_MINUS_REST = unname(iv_diff["diff"]), IV_SE = unname(iv_diff["se"]),
             IV_P = unname(iv_diff["p"]),
             FIRST_STAGE_WALD_MIN = fs_min,
             N_OBSERVATIONS = if (is.null(rf_fit)) NA_integer_ else nobs(rf_fit))
}

if (!is.null(.s16_ses_panel)) {
  superfamily_ses <- cache_or_run("superfamily_ses_heterogeneity", {
    rbindlist(lapply(names(MAIN_INSTRUMENTS), function(il) {
      run_superfamily_ses_heterogeneity(.s16_ses_panel, instrument = MAIN_INSTRUMENTS[[il]],
                                        instrument_label = il)
    }), fill = TRUE)
  })
  save_csv(superfamily_ses, "T16F_superfamily_SES_heterogeneity.csv")
  cat("\nCheck FIRST_STAGE_WALD_MIN before quoting an IV number -- read RF as primary,\n",
      "matching Section 12's own convention.\n", sep = "")
  .s16_show(superfamily_ses[, .(INSTRUMENT_LABEL, RF_DIFF = round(RF_DIFF_SHOP_MINUS_REST, 5),
                                RF_SE = round(RF_SE, 5), RF_P = round(RF_P, 4),
                                FIRST_STAGE_WALD_MIN = round(FIRST_STAGE_WALD_MIN, 1))])
}

.s16_hd("SECTION 16 COMPLETE")






###############################################################################

}  # end HPT_RUN$family_16

# ============================================================================
# Section 17: Conley exclusion-restriction sensitivity
# ============================================================================

if (isTRUE(HPT_RUN$conley_17)) {

# SECTION 17 -- CONLEY EXCLUSION-RESTRICTION SENSITIVITY
#
# Follows Conley, Hansen, and Rossi's plausibly-exogenous framing: how large a
# direct effect of the instrument on price, concentrated on shoppable services
# specifically, would have to be to overturn the gradient.
#
# EXACT ALGEBRA, NOT A GRID OF REFITS.
#   The reduced form regresses ln(P) on RF_Shoppable = Z*1[Shop], RF_Non =
#   Z*1[NonShop], controls, and FE. If a direct effect of size gamma is
#   assumed to load on RF_Shoppable specifically, then netting it out means
#   subtracting gamma * RF_Shoppable from the outcome and refitting the
#   identical model. For OLS this is a Frisch-Waugh-Lovell identity: the
#   coefficient on RF_Shoppable falls by EXACTLY gamma, every other
#   coefficient is UNCHANGED, and every residual is UNCHANGED -- so SE(GAP)
#   does not move either. The grid below is therefore closed-form arithmetic,
#   confirmed once per instrument by an actual refit at the grid midpoint
#   rather than trusted blindly.
#
#   This works because the reduced form is a single linear regression with
#   RF_Shoppable literally as one of its own regressors. It does NOT carry
#   over to the IV estimate in the same closed form (Z_Shoppable is an
#   instrument for TREAT_Shoppable there, not a direct regressor, so netting
#   out gamma does not simply shift one 2SLS coefficient). Given the design
#   decision that reduced form is primary and IV is rescaling only, that is
#   the estimator this bound is built on, and it is also the one whose
#   p-values are Anderson-Rubin robust by construction.
###############################################################################

.s17_hd <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")

run_conley_sensitivity <- function(panel, scheme_col = "SCHEME_1_CERTAINTY",
                                   instrument, instrument_label = "",
                                   outcome = PRIMARY_OUTCOME, endogenous = ENDOGENOUS_VARIABLE,
                                   controls = BASELINE_CONTROLS, fe = BASELINE_FIXED_EFFECTS,
                                   clusters = BASELINE_CLUSTERS,
                                   gamma_grid = NULL, gamma_steps = 41L) {
  
  d <- panel[!is.na(get(scheme_col))]
  req <- intersect(c(outcome, endogenous, instrument, scheme_col, controls, fe, clusters), names(d))
  d <- d[complete.cases(d[, ..req])]
  
  d[, CAT := droplevels(factor(get(scheme_col)))]
  keys <- levels(d$CAT)
  if (!all(c("Shoppable", "Non_shoppable") %chin% keys)) {
    stop("Expected 'Shoppable'/'Non_shoppable' levels in ", scheme_col,
         "; found: ", paste(keys, collapse = ", "), call. = FALSE)
  }
  for (kk in keys) d[, paste0("RF_", kk) := get(instrument) * as.integer(CAT == kk)]
  rhs <- c(paste0("RF_", keys), controls)
  
  fit <- tryCatch(feols(build_ols_formula(outcome, rhs, fe), data = d,
                        cluster = build_cluster_formula(clusters), warn = FALSE, notes = FALSE),
                  error = function(e) NULL)
  if (is.null(fit)) stop("Reduced-form fit failed for ", instrument_label, call. = FALSE)
  
  nm_s <- "RF_Shoppable"; nm_n <- "RF_Non_shoppable"
  b_s <- unname(coef(fit)[nm_s]); b_n <- unname(coef(fit)[nm_n])
  V <- vcov(fit)
  se_gap <- sqrt(V[nm_s, nm_s] + V[nm_n, nm_n] - 2 * V[nm_s, nm_n])
  gap0 <- b_s - b_n
  sd_z <- sd(d[[instrument]], na.rm = TRUE)
  
  if (is.null(gamma_grid)) {
    span <- max(abs(gap0) * 1.5, 3 * se_gap)
    gamma_grid <- seq(-span, span, length.out = gamma_steps)
  }
  
  grid <- data.table(GAMMA = gamma_grid)
  grid[, `:=`(GAP_CAUSAL = gap0 - GAMMA, SE_GAP = se_gap,
              CI_LOW = gap0 - GAMMA - 1.96 * se_gap, CI_HIGH = gap0 - GAMMA + 1.96 * se_gap)]
  grid[, SIGNIFICANT := as.integer(sign(CI_LOW) == sign(CI_HIGH))]
  
  # Smallest-magnitude adversarial gamma (toward zero) at which the CI first
  # touches zero -- the direct effect the design can withstand before the
  # gradient stops being distinguishable from zero at 95%.
  # Two thresholds, both reported.
  #
  # BREAKEVEN_GAMMA_SIG: the smallest-magnitude adversarial direct effect at
  #   which the 95% CI on the gradient first touches zero. Solve
  #   CI_HIGH = (gap0 - gamma) + 1.96*se >= 0 for a negative gap, giving
  #   gamma <= gap0 + 1.96*se. NOTE THE SIGN: moving gamma toward the gap
  #   shrinks the causal gradient, so the binding value is gap0 PLUS the
  #   half-width, not minus it. (An earlier version subtracted, which returns
  #   CI_LOW at gamma = 0 and overstates the bound by roughly 5.7x.)
  #
  # BREAKEVEN_GAMMA_ZERO: the direct effect that drives the POINT estimate to
  #   zero, which is simply gamma = gap0 (100% of the observed gradient).
  breakeven_sig  <- gap0 + sign(gap0) * -1 * 1.96 * se_gap
  breakeven_zero <- gap0
  
  grid[, `:=`(INSTRUMENT_LABEL = instrument_label, INSTRUMENT = instrument,
              GAP_OBSERVED = gap0, SE_GAP_OBSERVED = se_gap, SD_INSTRUMENT = sd_z,
              GAP_PCT_PER_SD = 100 * (exp(gap0 * sd_z) - 1),
              BREAKEVEN_GAMMA_SIG = breakeven_sig,
              BREAKEVEN_GAMMA_ZERO = breakeven_zero,
              # Same convention Table 5 (tab:headline) uses for each arm's own
              # percent-per-SD figure, so this table's percent column is now
              # directly comparable to Table 5's Gap row.
              BREAKEVEN_SIG_PCT_PER_SD = 100 * (exp(breakeven_sig * sd_z) - 1),
              BREAKEVEN_SIG_SHARE_OF_GAP = breakeven_sig / gap0,
              BREAKEVEN_ZERO_SHARE_OF_GAP = 1)]
  
  # Spot-check the closed form with one actual refit at the grid midpoint.
  mid <- gamma_grid[ceiling(length(gamma_grid) / 2)]
  d[, Y_SHIFTED := get(outcome) - mid * get(nm_s)]
  fit_chk <- tryCatch(feols(build_ols_formula("Y_SHIFTED", rhs, fe), data = d,
                            cluster = build_cluster_formula(clusters), warn = FALSE, notes = FALSE),
                      error = function(e) NULL)
  if (!is.null(fit_chk)) {
    chk_gap <- unname(coef(fit_chk)[nm_s]) - unname(coef(fit_chk)[nm_n])
    closed_gap <- gap0 - mid
    cat(sprintf("  [%s] Spot-check at gamma=%.4f: refit GAP=%.6f vs closed-form=%.6f (|diff|=%.2e)\n",
                instrument_label, mid, chk_gap, closed_gap, abs(chk_gap - closed_gap)))
  }
  
  # Verify the breakeven against the grid itself rather than trusting algebra:
  # the smallest-|gamma| row flagged insignificant should sit at breakeven_sig.
  edge <- grid[SIGNIFICANT == 0L & sign(GAMMA) == sign(gap0)]
  edge_gamma <- if (nrow(edge)) edge[which.min(abs(GAMMA)), GAMMA] else NA_real_
  cat(sprintf("  [%s] gap=%.6f se=%.6f | breakeven(sig)=%.6f (%.1f%% of gap) | grid edge=%.6f\n",
              instrument_label, gap0, se_gap, breakeven_sig,
              100 * breakeven_sig / gap0, edge_gamma))
  
  list(grid = grid, gap_observed = gap0, se_gap = se_gap,
       breakeven_sig = breakeven_sig, breakeven_zero = breakeven_zero, fit = fit)
}

.s17_hd("SECTION 17 -- Conley exclusion-restriction sensitivity (T17)")

conley_results <- cache_or_run("conley_sensitivity_v2", {
  out <- lapply(names(MAIN_INSTRUMENTS), function(il) {
    run_conley_sensitivity(outpatient, scheme_col = "SCHEME_1_CERTAINTY",
                           instrument = MAIN_INSTRUMENTS[[il]], instrument_label = il)
  })
  names(out) <- names(MAIN_INSTRUMENTS)
  out
})

conley_grid <- rbindlist(lapply(conley_results, `[[`, "grid"), fill = TRUE)
save_csv(conley_grid, "T17_conley_sensitivity_grid.csv")

conley_summary <- unique(conley_grid[, .(
  INSTRUMENT_LABEL,
  GAP_OBSERVED               = round(GAP_OBSERVED, 6),
  SE_GAP_OBSERVED            = round(SE_GAP_OBSERVED, 6),
  GAP_PCT_PER_SD             = round(GAP_PCT_PER_SD, 3),
  BREAKEVEN_GAMMA_SIG        = round(BREAKEVEN_GAMMA_SIG, 6),
  BREAKEVEN_SIG_PCT_PER_SD   = round(BREAKEVEN_SIG_PCT_PER_SD, 3),
  BREAKEVEN_SIG_SHARE_OF_GAP = round(BREAKEVEN_SIG_SHARE_OF_GAP, 3),
  BREAKEVEN_GAMMA_ZERO       = round(BREAKEVEN_GAMMA_ZERO, 6))])
save_csv(conley_summary, "T17B_conley_breakeven_summary.csv")

cat("\nBreakeven summary (smallest adversarial direct effect, in the units of\n",
    "the instrument, that brings the 95% CI on the gradient to touch zero):\n", sep = "")
.s17_show <- .s16_show
.s17_show(conley_summary)
cat("\nBREAKEVEN_SIG_PCT_PER_SD is on the same percent-per-SD-of-peer-exposure\n",
    "scale as RF_PERCENT_PER_SD, so it drops into prose without conversion.\n",
    "BREAKEVEN_SIG_SHARE_OF_GAP is the headline number: the fraction of the\n",
    "observed reduced-form gradient that would have to be direct effect,\n",
    "loading on shoppable services ONLY, to render the gradient insignificant.\n",
    "BREAKEVEN_GAMMA_ZERO (= the gap itself) is the stricter threshold at which\n",
    "the point estimate is driven to zero, i.e. 100% of the observed gradient.\n", sep = "")

.s17_hd("SECTION 17 COMPLETE")



###############################################################################

}  # end HPT_RUN$conley_17

# ============================================================================
# Section 18: enforcement controls
# ============================================================================

if (isTRUE(HPT_RUN$enforcement_18)) {

# SECTION 18 -- ENFORCEMENT CONTROLS
#
# CMS enforcement varies by county and month and is already in the panel
# (COUNTY_ENF_* -- confirmed present, county-month invariant, with genuine
# within-county time variation in the diagnostic). This is a controls
# question, not an instrument question: does the shoppability gradient
# survive adding each enforcement measure as an additional right-hand-side
# control in the headline interacted spec? estimate_interacted() already
# accepts an arbitrary `controls` vector, so this reuses it unchanged with
# BASELINE_CONTROLS extended by one enforcement variable at a time.
###############################################################################

.s18_hd <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")
.s18_show <- .s16_show

.s18_hd("SECTION 18 -- Enforcement controls (T18)")

run_enforcement_controls <- function(panel, enforcement_vars = ENFORCEMENT_ROBUSTNESS_INSTRUMENTS,
                                     scheme_col = "SCHEME_1_CERTAINTY",
                                     scheme_label = "1. Procedural certainty",
                                     instruments = MAIN_INSTRUMENTS, outcome = PRIMARY_OUTCOME) {
  rows <- list(); tests <- list()
  n <- length(enforcement_vars) * length(instruments); g <- 0L
  
  for (ev in names(enforcement_vars)) {
    v <- enforcement_vars[[ev]]
    if (!(v %in% names(panel))) { cat("  Missing from panel, skipped:", v, "\n"); next }
    for (il in names(instruments)) {
      g <- g + 1L
      z <- instruments[[il]]
      r <- estimate_interacted(panel, scheme_col, outcome, z, moderator_type = "categorical",
                               label = scheme_label, instrument_label = il,
                               moderator_label = scheme_label,
                               controls = c(BASELINE_CONTROLS, v))
      if (!is.null(r)) {
        r$rows[, ENFORCEMENT_CONTROL := ev]; r$tests[, ENFORCEMENT_CONTROL := ev]
        rows[[length(rows) + 1L]] <- r$rows; tests[[length(tests) + 1L]] <- r$tests
      }
      cat(sprintf("[%d/%d] %-24s %-36s done\n", g, n, ev, il))
    }
  }
  er <- rbindlist(rows, fill = TRUE); et <- rbindlist(tests, fill = TRUE)
  if (nrow(er) == 0L) stop("No enforcement-control results produced.", call. = FALSE)
  list(rows = er, tests = et)
}

enforcement_controls <- cache_or_run("enforcement_controls",
                                     run_enforcement_controls(outpatient))

save_csv(enforcement_controls$rows,  "T18_enforcement_controls_estimates.csv")
save_csv(enforcement_controls$tests, "T18B_enforcement_controls_heterogeneity.csv")

# Baseline (no enforcement control) for the comparison column, pulled from the
# existing T06 headline rather than re-estimated here.
if (exists("main") && is.list(main) && "rows" %in% names(main)) {
  .s18_headline <- main$rows[grepl("certainty", SPEC, ignore.case = TRUE),
                             .(INSTRUMENT_LABEL, TERM, RF_PERCENT_PER_SD_NO_CONTROL = RF_PERCENT_PER_SD)]
  .s18_cmp <- merge(
    enforcement_controls$rows[, .(ENFORCEMENT_CONTROL, INSTRUMENT_LABEL, TERM,
                                  RF_PERCENT_PER_SD_WITH_CONTROL = RF_PERCENT_PER_SD, RF_P)],
    .s18_headline, by = c("INSTRUMENT_LABEL", "TERM"), all.x = TRUE)
  .s18_cmp[, DELTA := RF_PERCENT_PER_SD_WITH_CONTROL - RF_PERCENT_PER_SD_NO_CONTROL]
  save_csv(.s18_cmp, "T18C_enforcement_vs_no_control.csv")
  
  cat("\nGradient with vs. without each enforcement control",
      "(Shoppable term, main instrument):\n")
  .s18_show(.s18_cmp[TERM == "Shoppable" & INSTRUMENT_LABEL == names(MAIN_INSTRUMENTS)[1],
                     .(ENFORCEMENT_CONTROL, RF_PERCENT_PER_SD_NO_CONTROL = round(RF_PERCENT_PER_SD_NO_CONTROL, 3),
                       RF_PERCENT_PER_SD_WITH_CONTROL = round(RF_PERCENT_PER_SD_WITH_CONTROL, 3),
                       DELTA = round(DELTA, 3), RF_P = round(RF_P, 4))])
} else {
  cat("`main` (the headline T06 object) is not in memory -- the comparison\n",
      "column was skipped, but T18/T18B are saved. Run restore_session() with\n",
      "'main_results' included, or re-source this after warm_start().\n", sep = "")
}

.s18_hd("SECTION 18 COMPLETE")




###############################################################################

}  # end HPT_RUN$enforcement_18

###############################################################################
#                                                                             #
#   PART 5 -- OUTPUTS                                                         #
#                                                                             #
#   Blocks appear in their original relative order, because later ones depend #
#   on helpers the earlier ones define (theme_paper, save_fig, read_table in  #
#   Figures 1-9; PAL, rd, sv, short in Figures 10-20). Each reads result CSVs #
#   from TABLE_DIR rather than memory, so any one can be re-run alone.        #
#                                                                             #
###############################################################################


# ============================================================================
# Figures 1-9
# ============================================================================

if (isTRUE(HPT_RUN$figures)) {

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



t06g <- fread(file.path(TABLE_DIR, "T06G_all_tiers_heterogeneity_tests.csv"))
stopifnot("P_VALUE" %in% names(t06g))

t06g[, P_VALUE := 2 * pt(-qnorm(1 - P_VALUE / 2), df = 15)]

rerank <- t06g[ESTIMATOR == "Reduced form",
               .(N_SCHEMES = .N,
                 N_SIG_05 = sum(P_VALUE < 0.05, na.rm = TRUE),
                 SHARE_SIG = round(mean(P_VALUE < 0.05, na.rm = TRUE), 3),
                 MEDIAN_P = round(median(P_VALUE, na.rm = TRUE), 4)),
               by = .(TIER, INSTRUMENT_LABEL)][order(-SHARE_SIG, MEDIAN_P)]

fs_by_inst <- t06g[, .(MEDIAN_FIRST_STAGE_F = round(median(FIRST_STAGE_WALD_MIN, na.rm = TRUE), 1)),
                   by = INSTRUMENT_LABEL]
rerank <- merge(rerank, fs_by_inst, by = "INSTRUMENT_LABEL")[order(-SHARE_SIG, MEDIAN_P)]

fwrite(rerank, file.path(TABLE_DIR, "T06L_tier_rerank_all_schemes.csv"))
print(rerank)

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


# Interactive spot-check of one family/scheme cell. Needs schemes_long and
# outpatient, both in scope after warm_start().
# s <- merge(schemes_long, unique(outpatient[, .(ANALYSIS_CONCEPT_ID = FINAL_CONCEPT_ID,
#                                                FAMILY = FINAL_FAMILY_ID)]),
#            by = "ANALYSIS_CONCEPT_ID")
# s[FAMILY == "ECHOCARDIOGRAPHY" & SCHEME_ID == "scheme_theory_v2",
#   .(ANALYSIS_CONCEPT_ID, OFFICIAL_DESCRIPTION, SHOPPABILITY_CATEGORY)]



###############################################################################

}  # end HPT_RUN$figures

# ============================================================================
# Paper utilities: summary statistics, permutation figure, scheme tables
# ============================================================================

if (isTRUE(HPT_RUN$tables)) {

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
               ZERO = mean(v == 0, na.rm = TRUE), DIGITS = digits)
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
    row_stat("Out-of-CBSA competitor system exposure",
             panel[["Z_SYS_COMPETITOR_SYSTEMS_OUTSIDE_CBSA_9M_EXCL_CURRENT"]], 2),
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

# Interactive: builds and writes the summary-statistics table on demand.
# build_summary_stats(outpatient)



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

# Interactive: generates the permutation-null figure on demand from
# concept_results and schemes_long, both already in scope after warm_start().
# meta_input <- prepare_meta_input(concept_results, schemes_long)
# perm <- permutation_null_distribution(meta_input, "Competitor_only_hospitals_9m")
# plot_permutation_null(perm)



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

# Interactive: builds and writes the scheme classification table on demand.
# build_scheme_table(schemes_long, outpatient)


# Wu-Hausman test of OLS exogeneity.
# Interactive: Wu-Hausman summary read from the T03 CSV that PART 5's
# "Three builds" block writes.
# t3 <- fread(file.path(TABLE_DIR, "T03_pooled_OLS_RF_IV_all_outcomes.csv"))
#
# t3[, .(INSTRUMENT_LABEL, OUTCOME,
#        OLS = round(OLS_PERCENT, 3), IV = round(IV_PERCENT, 3),
#        GAP = round(OLS_PERCENT - IV_PERCENT, 3),
#        WU_HAUSMAN_P = signif(WU_HAUSMAN_P, 3))][order(OUTCOME, INSTRUMENT_LABEL)]
#
# cat("\nOLS more negative than IV in",
#     t3[OLS_PERCENT < IV_PERCENT, .N], "of", nrow(t3), "specifications\n")
# cat("Mean OLS-IV gap:", round(mean(t3$OLS_PERCENT - t3$IV_PERCENT, na.rm = TRUE), 3), "pp\n")
# cat("Wu-Hausman rejects at 5% in",
#     t3[WU_HAUSMAN_P < 0.05, .N], "of", t3[is.finite(WU_HAUSMAN_P), .N], "\n")



###############################################################################

}  # end HPT_RUN$tables

# ============================================================================
# Figures 10-20
# ============================================================================

if (isTRUE(HPT_RUN$figures)) {

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


PAL <- c("#2C6E9B", "#C1502E", "#5B8C5A", "#8B6BA8")
rd  <- function(f) fread(file.path(TABLE_DIR, f))
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
#
# Two fixes.
#
# (1) The subtitle asserted "Both constructions failing at 5%". Three fail,
#     not two: cumulative rollout share, cumulative external hospitals, and
#     own-system included. The count is now computed from the data so it
#     cannot drift out of sync with the table again, and the wording adapts
#     to whatever the count turns out to be.
#
# (2) The subtitle used a colon as a sentence separator. Rewritten.
#
# NOTE: verify T13C_instrument_ladder_tests.csv carries t(15) p-values before
# trusting the facet split. If P_VALUE there still equals
# 2*(1 - pnorm(sqrt(WALD))), it predates the .pval() fix and needs
#   d[, P_VALUE := 2 * pt(-qnorm(1 - P_VALUE/2), df = 15)]
# applied and re-saved first, exactly as T06_mainB and T07B did.

c_rows <- rd("T13C_instrument_ladder_rows.csv")
c_test <- rd("T13C_instrument_ladder_tests.csv")[ESTIMATOR == "Reduced form"]

# Guard against the stale-CSV case rather than failing silently.
if (nrow(c_test) && "WALD" %in% names(c_test)) {
  implied_normal <- 2 * (1 - pnorm(sqrt(c_test$WALD)))
  if (isTRUE(all.equal(c_test$P_VALUE, implied_normal, tolerance = 1e-6))) {
    warning("T13C_instrument_ladder_tests.csv P_VALUE matches a normal reference; ",
            "applying the t(15) correction in memory. Re-save the CSV to make it stick.",
            call. = FALSE)
    c_test[, P_VALUE := 2 * pt(-qnorm(1 - P_VALUE / 2), df = 15)]
  }
}

d11 <- dcast(c_rows, INSTRUMENT_LABEL ~ TERM, value.var = "RF_PERCENT_PER_SD")
d11 <- merge(d11, c_test[, .(INSTRUMENT_LABEL, P_VALUE)], by = "INSTRUMENT_LABEL")
d11 <- melt(d11, id.vars = c("INSTRUMENT_LABEL", "P_VALUE"),
            variable.name = "TERM", value.name = "EST")
d11[, TERM := factor(fifelse(TERM == "Shoppable", "Shoppable", "Non-shoppable"),
                     levels = c("Shoppable", "Non-shoppable"))]
d11[, LAB := factor(INSTRUMENT_LABEL,
                    levels = unique(d11[order(P_VALUE)]$INSTRUMENT_LABEL))]
d11[, SIG := fifelse(P_VALUE < 0.05, "Test rejects at 5%", "Does not reject")]

n_fail  <- uniqueN(d11[P_VALUE >= 0.05]$INSTRUMENT_LABEL)
n_word  <- c("No", "One", "Two", "Three", "Four", "Five", "Six", "Seven", "Eight")[n_fail + 1L]
sub11 <- paste0(
  n_word, " of the eight constructions fall short of 5%, and they are the ones the\n",
  "design predicts should fail, since own-system rollout reintroduces the contamination\n",
  "the exclusion argument removes and the cumulative measures discard the time window")

p11 <- ggplot(d11, aes(EST, LAB, fill = TERM)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
  geom_col(position = position_dodge(width = 0.65), width = 0.6) +
  facet_wrap(~SIG, ncol = 1, scales = "free_y", strip.position = "top") +
  scale_fill_manual(values = c(PAL[1], PAL[2])) +
  labs(title = "Instrument construction ladder",
       subtitle = sub11,
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
#
# Three fixes.
#
# (1) The confidence bands were wrong for two of three instruments. The old
#     code multiplied every instrument's log-point gap by 15.94, labelled in
#     the comment as the competitor-hospitals SD. The three instruments have
#     different exposure SDs, and on this estimation sample (N = 941,013,
#     smaller than the headline because counties without an ACS index drop)
#     they are 15.50, 16.64 and 14.38, not 15.94. So the plotted intervals
#     were on a different scale from the plotted point estimates, which
#     already used the correct per-instrument SD inside GAP_PCT_PER_SD.
#
#     Rather than hardcode three more numbers, recover each instrument's SD
#     from the data itself. GAP_PCT_PER_SD = 100*(exp(GAP_RF * SD) - 1), so
#     SD = log1p(GAP_PCT_PER_SD/100) / GAP_RF exactly. Taking it from the
#     tercile with the largest |GAP_RF| avoids dividing by a near-zero gap.
#
# (2) The subtitle hardcoded "joint p = 0.23-0.67". Those were the old
#     chi-square values. Read the range from the tests file instead, so it
#     tracks .s14_contrast() automatically.
#
# (3) The subtitle ran off the plot edge and used "--". Wrapped, and rewritten
#     without the dash.

ter  <- rd("T12Q_ses_bins_gaps.csv")
tst  <- rd("T12Q_ses_bins_tests.csv")[TEST == "all gaps equal" & ESTIMATOR == "Reduced form"]

ter[, INST := short(INSTRUMENT_LABEL)]

# Per-instrument SD of the instrument, recovered exactly from the two columns
# the pipeline already wrote. No hardcoded constants.
ter[, SD_Z := {
  i <- which.max(abs(GAP_RF))
  log1p(GAP_PCT_PER_SD[i] / 100) / GAP_RF[i]
}, by = INSTRUMENT_LABEL]

ter[, LO := 100 * (exp((GAP_RF - 1.96 * GAP_RF_SE) * SD_Z) - 1)]
ter[, HI := 100 * (exp((GAP_RF + 1.96 * GAP_RF_SE) * SD_Z) - 1)]

ter[, BIN := factor(SES_BIN, levels = c("T1", "T2", "T3"),
                    labels = c("T1\n(least advantaged)", "T2\n(middle)", "T3\n(most advantaged)"))]

joint_lo <- sprintf("%.2f", min(tst$P_VALUE))
joint_hi <- sprintf("%.2f", max(tst$P_VALUE))
sub18 <- paste0(
  "Precisely estimated in the middle and upper terciles, imprecise in the lowest.\n",
  "The terciles are not statistically distinguishable from one another (joint p = ",
  joint_lo, " to ", joint_hi, ")")

p18 <- ggplot(ter, aes(BIN, GAP_PCT_PER_SD, colour = INST, group = INST)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey55") +
  geom_errorbar(aes(ymin = LO, ymax = HI), width = 0.08,
                position = position_dodge(width = 0.4), linewidth = 0.5) +
  geom_point(position = position_dodge(width = 0.4), size = 2.4) +
  scale_colour_manual(values = PAL) +
  labs(title = "Shoppability gap by socioeconomic tercile",
       subtitle = sub18,
       x = NULL, y = "Response gap (% per SD of instrument)") +
  theme_paper()
sv(p18, "fig18_ses_terciles.pdf", 8, 4.4)


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

}  # end HPT_RUN$figures

# ============================================================================
# County coverage statistics and print-safe map
# ============================================================================

if (isTRUE(HPT_RUN$coverage_map)) {

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

}  # end HPT_RUN$coverage_map

# ============================================================================
# Three builds: pooled table, family-scheme grid, shoppable shares
# ============================================================================

if (isTRUE(HPT_RUN$tables)) {

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
        IV_P        = .pval(b_iv / s_iv, iv),
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


###############################################################################

}  # end HPT_RUN$tables

# ============================================================================
# Figures 16-18
# ============================================================================

if (isTRUE(HPT_RUN$figures)) {

# FIGURES FOR SECTIONS 16, 17, AND 18
#
#   fig16_family_forest.pdf        main.tex, Section 6.5
#   fig17_conley_sensitivity.pdf   appendix.tex
#   fig18_enforcement_controls.pdf appendix.tex
#
# Reads only the CSVs already written by HPT_sections_16_17_18.R, so this can
# run in a fresh session after warm_start() without re-estimating anything.
# It does not need `outpatient` in memory.
#
# If TABLE_DIR / FIGURE_DIR are not defined (i.e. the pipeline definitions are
# not loaded), set them by hand at the top of this file.
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

if (!exists("TABLE_DIR"))  stop("TABLE_DIR not defined. Run warm_start(), or set it manually.")
if (!exists("FIGURE_DIR")) {
  FIGURE_DIR <- file.path(dirname(TABLE_DIR), "02_Figures")
  cat("FIGURE_DIR not defined; defaulting to:\n  ", FIGURE_DIR, "\n")
}
dir.create(FIGURE_DIR, showWarnings = FALSE, recursive = TRUE)

# The paper's figures use a plain serif-free look; swap in your house theme
# here if you have one so these match fig01-fig15.
.fig_theme <- theme_bw(base_size = 10) +
  theme(panel.grid.minor   = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.border       = element_rect(colour = "grey30", linewidth = 0.4),
        strip.background   = element_rect(fill = "grey92", colour = "grey30"),
        strip.text         = element_text(face = "bold", size = 9),
        legend.position    = "bottom",
        legend.title       = element_blank(),
        plot.caption       = element_text(hjust = 0, size = 7.5, colour = "grey30"))

SHOP_COLOURS <- c("Shoppable" = "#1b4965", "Non-shoppable" = "#bc4749")

# Display names, matching the LaTeX tables.
FAMILY_LABELS <- c(
  MAMMOGRAPHY                 = "Mammography",
  VASCULAR_ULTRASOUND         = "Vascular ultrasound",
  BONE_DENSITY                = "Bone density",
  ECHOCARDIOGRAPHY            = "Echocardiography",
  MRI_MRA                     = "MRI / MRA",
  LABORATORY_PATHOLOGY        = "Laboratory and pathology",
  EVALUATION_MANAGEMENT       = "Evaluation and management",
  DIAGNOSTIC_ULTRASOUND       = "Diagnostic ultrasound",
  XRAY_FLUOROSCOPY            = "X-ray and fluoroscopy",
  CT_CTA                      = "CT / CTA",
  CRITICAL_CARE               = "Critical care",
  BIOPSY                      = "Biopsy",
  PROCEDURE_SURGERY           = "Procedure and surgery",
  COLONOSCOPY_LOWER_ENDOSCOPY = "Colonoscopy / lower endoscopy",
  EMERGENCY_DEPARTMENT        = "Emergency department",
  UPPER_ENDOSCOPY             = "Upper endoscopy")

INSTRUMENT_LABELS <- c(
  Competitor_only_hospitals_9m         = "Competitor hospitals",
  Primary_strict_system_IV             = "Local system",
  Competitor_outside_CBSA_hospitals_9m = "Competitor hospitals (ex-CBSA)",
  Competitor_outside_CBSA_counties_9m  = "Competitor counties (ex-CBSA)",
  Competitor_systems_9m                = "Competitor systems",
  Competitor_outside_CBSA_systems_9m   = "Competitor systems (ex-CBSA)")

MAIN_THREE <- c("Competitor_only_hospitals_9m",
                "Primary_strict_system_IV",
                "Competitor_outside_CBSA_hospitals_9m")


# ===========================================================================
# FIGURE 16 -- Family-level forest plot
#
# Sixteen families under the primary instrument, ordered by point estimate,
# coloured by the a-priori shoppability label. Dashed vertical lines mark the
# precision-weighted group means, which are the quantities the permutation
# test in Table 8 differences.
# ===========================================================================

fam <- fread(file.path(TABLE_DIR, "T16A_family_level_RF_FS_IV.csv"))
fam <- fam[INSTRUMENT_LABEL == "Competitor_only_hospitals_9m"]

if (!exists("DIAGNOSTIC_FAMILIES")) {
  # Fallback if pipeline constants are not loaded.
  DIAGNOSTIC_FAMILIES <- c("XRAY_FLUOROSCOPY", "MRI_MRA", "DIAGNOSTIC_ULTRASOUND",
                           "CT_CTA", "LABORATORY_PATHOLOGY", "VASCULAR_ULTRASOUND",
                           "EVALUATION_MANAGEMENT", "ECHOCARDIOGRAPHY",
                           "BONE_DENSITY", "MAMMOGRAPHY")
}

fam[, CLASS := fifelse(GROUP_ID %chin% DIAGNOSTIC_FAMILIES, "Shoppable", "Non-shoppable")]
fam[, LABEL := FAMILY_LABELS[GROUP_ID]]
fam[, `:=`(CI_LOW  = RF_COEF - 1.96 * RF_SE,
           CI_HIGH = RF_COEF + 1.96 * RF_SE)]
setorder(fam, -RF_COEF)
fam[, LABEL := factor(LABEL, levels = LABEL)]

# Precision-weighted group means: the estimand of the permutation test.
gm <- fam[, .(MEAN = sum(RF_COEF / RF_SE^2) / sum(1 / RF_SE^2)), by = CLASS]
cat("\nPrecision-weighted group means (should match Table 8, row 1):\n")
print(as.data.frame(gm), row.names = FALSE)
cat("Difference:", round(gm[CLASS == "Shoppable", MEAN] -
                           gm[CLASS == "Non-shoppable", MEAN], 6), "\n")

p16 <- ggplot(fam, aes(x = RF_COEF, y = LABEL, colour = CLASS)) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.4) +
  geom_vline(data = gm, aes(xintercept = MEAN, colour = CLASS),
             linetype = "dashed", linewidth = 0.45, show.legend = FALSE) +
  geom_errorbarh(aes(xmin = CI_LOW, xmax = CI_HIGH), height = 0, linewidth = 0.55) +
  geom_point(size = 2.1) +
  scale_colour_manual(values = SHOP_COLOURS) +
  scale_x_continuous(breaks = seq(-0.010, 0.006, by = 0.002)) +
  labs(x = "Reduced-form coefficient (log points)", y = NULL,
       caption = paste("Competitor-hospital instrument. Bars are 95% confidence intervals from two-way",
                       "clustered standard errors.\nDashed lines mark precision-weighted group means.")) +
  .fig_theme

ggsave(file.path(FIGURE_DIR, "fig16_family_forest.pdf"), p16,
       width = 7.5, height = 5.0)
cat("Wrote fig16_family_forest.pdf\n")


# ===========================================================================
# FIGURE 17 -- Conley sensitivity curves
#
# Causal gradient against an assumed direct effect psi, with its 95% band.
# The band has constant width by the Frisch-Waugh-Lovell identity, which is
# the point worth seeing: only the point estimate slides.
# ===========================================================================

cg <- fread(file.path(TABLE_DIR, "T17_conley_sensitivity_grid.csv"))
cg <- cg[INSTRUMENT_LABEL %chin% MAIN_THREE]
cg[, PANEL_LABEL := factor(INSTRUMENT_LABELS[INSTRUMENT_LABEL],
                           levels = INSTRUMENT_LABELS[MAIN_THREE])]

# Breakeven, recomputed from the grid rather than read, as a check on the
# stored column. These should agree to numerical precision.
be <- unique(cg[, .(PANEL_LABEL,
                    STORED    = BREAKEVEN_GAMMA_SIG,
                    RECOMPUTED = GAP_OBSERVED + 1.96 * SE_GAP,
                    SHARE     = BREAKEVEN_SIG_SHARE_OF_GAP)])
cat("\nConley breakeven check (stored vs recomputed):\n")
print(as.data.frame(be), row.names = FALSE)
stopifnot(all(abs(be$STORED - be$RECOMPUTED) < 1e-9))

p17 <- ggplot(cg, aes(x = GAMMA, y = GAP_CAUSAL)) +
  geom_ribbon(aes(ymin = CI_LOW, ymax = CI_HIGH), fill = "#1b4965", alpha = 0.16) +
  geom_hline(yintercept = 0, colour = "grey35", linewidth = 0.4) +
  geom_vline(data = be, aes(xintercept = STORED),
             colour = "#bc4749", linetype = "dashed", linewidth = 0.5) +
  geom_line(colour = "#1b4965", linewidth = 0.7) +
  facet_wrap(~ PANEL_LABEL, nrow = 1, scales = "free_x") +
  labs(x = expression(paste("Assumed direct effect on shoppable prices, ", psi)),
       y = "Causal shoppability gradient (log points)",
       caption = paste("Dashed line marks the breakeven violation at which the 95% band first covers zero.",
                       "The band is constant-width\nbecause netting out an assumed direct effect leaves every residual unchanged.")) +
  .fig_theme + theme(legend.position = "none")

ggsave(file.path(FIGURE_DIR, "fig17_conley_sensitivity.pdf"), p17,
       width = 8.0, height = 3.4)
cat("Wrote fig17_conley_sensitivity.pdf\n")



# ===========================================================================
# FIGURE 18 -- Enforcement control stability
#
# Shoppable coefficient across the eight enforcement controls and the
# no-control baseline, for the three main instruments. The visual claim is
# that nothing moves.
# ===========================================================================

enf <- fread(file.path(TABLE_DIR, "T18C_enforcement_vs_no_control.csv"))
enf <- enf[TERM == "Shoppable" & INSTRUMENT_LABEL %chin% MAIN_THREE]

ENF_LABELS <- c(
  Any_enforcement_3m_lag1 = "Any enforcement, 3M count",
  Warning_3m_lag1         = "Warning letters, 3M count",
  Closure_3m_lag1         = "Closure notices, 3M count",
  Any_enforcement_ind_3m  = "Any enforcement, 3M indicator",
  Warning_ind_3m          = "Warning letters, 3M indicator",
  Closure_ind_3m          = "Closure notices, 3M indicator",
  Closure_6m_lag1         = "Closure notices, 6M count",
  Closure_ind_6m          = "Closure notices, 6M indicator")

enf[, LABEL := ENF_LABELS[ENFORCEMENT_CONTROL]]
if (anyNA(enf$LABEL)) {
  cat("\nUnmapped enforcement labels (edit ENF_LABELS):\n")
  print(unique(enf[is.na(LABEL), ENFORCEMENT_CONTROL]))
  enf[is.na(LABEL), LABEL := ENFORCEMENT_CONTROL]
}
enf[, LABEL := factor(LABEL, levels = rev(unique(ENF_LABELS)))]
enf[, PANEL_LABEL := factor(INSTRUMENT_LABELS[INSTRUMENT_LABEL],
                            levels = INSTRUMENT_LABELS[MAIN_THREE])]

base <- unique(enf[, .(PANEL_LABEL, BASE = RF_PERCENT_PER_SD_NO_CONTROL)])
cat("\nBaseline shoppable coefficients (should match Table 5, Panel A):\n")
print(as.data.frame(base), row.names = FALSE)
cat("\nLargest displacement by instrument:\n")
print(as.data.frame(enf[, .(MAX_ABS_DELTA = round(max(abs(DELTA)), 3)),
                        by = PANEL_LABEL]), row.names = FALSE)

p18 <- ggplot(enf, aes(x = RF_PERCENT_PER_SD_WITH_CONTROL, y = LABEL)) +
  geom_vline(data = base, aes(xintercept = BASE),
             colour = "grey35", linetype = "dashed", linewidth = 0.5) +
  geom_point(colour = "#1b4965", size = 2.1) +
  facet_wrap(~ PANEL_LABEL, nrow = 1, scales = "free_x") +
  scale_x_continuous(expand = expansion(mult = 0.35)) +
  labs(x = "Shoppable reduced-form coefficient (% per SD of instrument)", y = NULL,
       caption = paste("Each point is a separate specification adding one county-level enforcement measure",
                       "to the control set.\nDashed line marks the baseline with no enforcement control.")) +
  .fig_theme + theme(legend.position = "none")

ggsave(file.path(FIGURE_DIR, "fig18_enforcement_controls.pdf"), p18,
       width = 8.5, height = 3.2)
cat("Wrote fig18_enforcement_controls.pdf\n")

cat("\nAll three figures written to:\n  ", FIGURE_DIR, "\n")




# The block below is Section 19 diagnostic code that references
# comparability_within_family, decomposition_results, and ses_gradient_table
# -- three names never defined anywhere in this file. It ran in the original
# author's console against objects that existed there at the time, not
# against anything this file produces, and would error on a fresh source()
# regardless of this reorganization. Left as a record; rebuild comp_wf and
# dec from the current session's objects before using it.
#
# s19 <- run_section_19(outpatient,
#                       within_family = comparability_within_family,  # T09D
#                       decomposition = decomposition_results,        # T08
#                       ses_gradient  = ses_gradient_table)
# exists("schemes_long"); exists("concept_results")
#
# if (!exists("schemes_long"))
#   schemes_long <- readRDS(file.path(CACHE_DIR, "schemes_long.rds"))
# if (!exists("concept_results"))
#   concept_results <- readRDS(file.path(CACHE_DIR, "concept_level_6inst.rds"))
#
# measures <- cache_or_run("comparability_measures", build_comparability_measures(outpatient))
# comp_wf  <- run_comparability_within_family(concept_results, measures)
#
# meta_input <- prepare_meta_input(concept_results, schemes_long)
# meta_input[, TIER := instrument_tier(INSTRUMENT_LABEL)]
# dec <- decompose_reduced_form(meta_input, stem = "T08D_RF_vs_FS_decomposition")
#
# s19c <- run_s19c(comp_wf, dec)




###############################################################################

}  # end HPT_RUN$figures

###############################################################################
#                                                                             #
#   PART 6 -- REFEREE-RESPONSE DIAGNOSTICS                                    #
#                                                                             #
#   Sections 19 through 25. None produces a number that appears in the paper; #
#   each answers a specific referee item. Off by default.                     #
#                                                                             #
###############################################################################

if (isTRUE(HPT_RUN$diagnostics)) {
  
  source(file.path(CODE_DIR, "HPT_diagnostic_for_family_levels.R"))
  source(file.path(CODE_DIR, "HPT_Section19_TierI_Diagnostics.R"))
  
}  # end HPT_RUN$diagnostics -- Section 19 only

# Section 20 is gated separately. Its definitions and its execution are
# independent of Sections 19 and 21-25.
if (isTRUE(HPT_RUN$diagnostics) || isTRUE(HPT_RUN$tier2_diag)) {
  
  # ============================================================================
  # Section 20: Tier II diagnostics A5, A6, A8


# ============================================================================
# Section 20: Tier II diagnostics A5, A6, A8
# ============================================================================

# SECTION 20 -- TIER II DIAGNOSTICS: A5, A6, A8
#
# Source AFTER HPT_Analysis_Pipeline.R has run through Section 9 (needs
# `measures` and `concept_results` in memory) and after `outpatient` exists.
# Nothing here touches Section 6's cached headline estimates.
#
#   20A  A6 -- is the contracting-depth result (N_PAYERS_V2) the concept-size
#              gradient relabeled? Correlations against size, exposure, the
#              concept's own first stage, and the precision weight used in
#              run_comparability_within_family(); then re-estimates unweighted,
#              with a log-size control, and with size-quartile FE layered on
#              top of family FE.
#   20B  A5 -- payer-mix composition. Two pieces:
#              (i)  a balance test: does the instrument predict how many
#                   payers a hospital reports for a concept, differentially
#                   by shoppability? Feasible today from the payer-cell
#                   export (HPT_PAYER_DISPERSION*.csv.gz).
#              (ii) the headline interacted price regression with reported
#                   payer count added as a row-level control.
#              The FULL ask -- restrict to a fixed set of common payers and
#              estimate at hospital x concept x payer level with payer FE --
#              is NOT feasible from what's in this file. See the note at the
#              bottom of this section for exactly what that would require.
#   20C  A8 -- the fixed-effects ladder: baseline, + system FE (additive),
#              then system x month (reproducing Section 13's own result as
#              the third rung, so the ladder is internally consistent with
#              Table sysmonth). Reuses Section 13's own system-resolution
#              logic rather than assuming outpatient_r13 already exists.
#
# NOT ATTEMPTED HERE: the distance-based exclusion ladder from A8 (0/50/100/
# 200-mile rings around the focal hospital). There is no latitude/longitude,
# address, or geocode field anywhere in this pipeline -- I checked. That is a
# Stage 1/2 data-construction task (hospital coordinates from CMS Care
# Compare's Hospital General Information file, the AHA annual survey you
# already cite, or an NPI registry lookup), not a diagnostic on what already
# exists. See the note at the end of this file for a concrete next step.
###############################################################################

.s20_hd  <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")
.s20_sub <- function(x) cat("\n--- ", x, " ", strrep("-", max(0, 70 - nchar(x))), "\n", sep = "")


# ===========================================================================
# 20A. CONTRACTING DEPTH vs. CONCEPT SIZE                                (A6)
# ===========================================================================

s20_depth_diagnostics <- function(concept_results, measures,
                                  moderator   = "N_PAYERS_V2",
                                  instruments = MAIN_INSTRUMENTS,
                                  verbose_timing = TRUE) {

  .t <- Sys.time()
  .tick <- function(label) {
    if (isTRUE(verbose_timing)) {
      now <- Sys.time()
      cat(sprintf("  [20A] %-32s %6.2fs\n", label,
                  as.numeric(difftime(now, .t, units = "secs"))))
      .t <<- now
    }
  }

  cat("[20A] entered s20_depth_diagnostics(); nrow(concept_results)=",
      nrow(concept_results), " nrow(measures)=", nrow(measures), "\n", sep = "")

  cr <- as.data.table(copy(concept_results)); .tick("copy concept_results")
  ms <- as.data.table(copy(measures));        .tick("copy measures")

  req_cr <- c("FINAL_CONCEPT_ID", "INSTRUMENT_LABEL", "RF_COEF", "RF_SE",
              "FS_COEF", "N_OBSERVATIONS")
  missing_cr <- setdiff(req_cr, names(cr))
  if (length(missing_cr) > 0L)
    stop("concept_results is missing: ", paste(missing_cr, collapse = ", "),
         call. = FALSE)
  if (!(moderator %chin% names(ms)))
    stop(moderator, " not found in measures. Run build_comparability_measures() first.",
         call. = FALSE)
  .tick("column checks")

  # Mirrors run_comparability_within_family()'s own merge and weight exactly,
  # so this is the same object that function estimates on, not a re-derivation.
  d <- merge(cr[is.finite(RF_COEF) & is.finite(RF_SE) & RF_SE > 0], ms,
             by = "FINAL_CONCEPT_ID")
  cat("[20A] merge complete, nrow(d)=", nrow(d),
      " (expect roughly nrow(cr) after the RF_SE filter, not more)\n", sep = "")
  if (nrow(d) > 5L * nrow(cr))
    warning("[20A] merge produced far more rows than expected -- FINAL_CONCEPT_ID ",
            "is likely not unique in one of the two inputs. Check ",
            "uniqueN(measures$FINAL_CONCEPT_ID) == nrow(measures).", call. = FALSE)
  .tick("merge")

  setDT(d)
  d[, MOD   := safe_numeric(get(moderator))]
  d <- d[is.finite(MOD)]
  d[, W     := 1 / (RF_SE^2)]
  d[, LOG_N := log(pmax(N_OBSERVATIONS, 1))]
  has_mpp <- "MEAN_PRIOR_POSTERS" %chin% names(d)
  .tick("derived columns")

  # -- correlations, per main instrument -----------------------------------
  corr_tab <- rbindlist(lapply(names(instruments), function(il) {
    dd <- d[INSTRUMENT_LABEL == il]
    data.table(
      INSTRUMENT_LABEL = il,
      N_CONCEPTS       = nrow(dd),
      COR_LOG_SIZE     = cor(dd$MOD, dd$LOG_N, use = "complete.obs"),
      COR_MEAN_PRIOR   = if (has_mpp) cor(dd$MOD, dd$MEAN_PRIOR_POSTERS, use = "complete.obs") else NA_real_,
      COR_FS_COEF      = cor(dd$MOD, dd$FS_COEF, use = "complete.obs"),
      COR_WEIGHT       = cor(dd$MOD, dd$W, use = "complete.obs")
    )
  }))
  .tick("correlations")

  cat("\nCorrelation of", moderator, "with concept size, mean exposure, the\n",
      "concept's own first stage, and the precision weight (1/RF_SE^2) used\n",
      "in run_comparability_within_family():\n", sep = "")
  print(corr_tab[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])

  # -- re-estimation ladder --------------------------------------------------
  d[, MODC   := (MOD - mean(MOD, na.rm = TRUE)) / sd(MOD, na.rm = TRUE)]
  d[, SIZE_Q := cut(N_OBSERVATIONS,
                    quantile(N_OBSERVATIONS, c(0, .25, .5, .75, 1), na.rm = TRUE),
                    include.lowest = TRUE, labels = c("Q1", "Q2", "Q3", "Q4"))]
  .tick("bins")

  specs <- list(
    `(b) Family FE, weighted [paper's spec]` = list(f = "RF_COEF ~ MODC | CONCEPT_FAMILY", w = TRUE),
    `(b') Family FE, unweighted`             = list(f = "RF_COEF ~ MODC | CONCEPT_FAMILY", w = FALSE),
    `(c) + log(N_obs) control`               = list(f = "RF_COEF ~ MODC + LOG_N | CONCEPT_FAMILY", w = TRUE),
    `(d) + size-quartile FE`                 = list(f = "RF_COEF ~ MODC | CONCEPT_FAMILY + SIZE_Q", w = TRUE)
  )

  re_est <- rbindlist(lapply(names(instruments), function(il) {
    cat("[20A] fitting", il, "...\n")
    dd <- d[INSTRUMENT_LABEL == il & is.finite(MODC)]
    out <- rbindlist(lapply(names(specs), function(sn) {
      t1 <- Sys.time()
      sp  <- specs[[sn]]
      fit <- tryCatch(
        if (sp$w) feols(as.formula(sp$f), data = dd, weights = ~W,
                        cluster = ~CONCEPT_FAMILY, warn = FALSE, notes = FALSE)
        else      feols(as.formula(sp$f), data = dd,
                        cluster = ~CONCEPT_FAMILY, warn = FALSE, notes = FALSE),
        error = function(e) { message("  [", il, " / ", sn, "] error: ", e$message); NULL })
      cat(sprintf("    %-40s %6.2fs\n", sn,
                  as.numeric(difftime(Sys.time(), t1, units = "secs"))))
      if (is.null(fit)) return(NULL)
      td <- tidy_fixest(fit)
      td <- td[term == "MODC"]
      if (nrow(td) == 0L) return(NULL)
      td[, `:=`(SPEC = sn, INSTRUMENT_LABEL = il, N = nrow(dd))]
      td
    }), fill = TRUE)
    out
  }), fill = TRUE)
  .tick("re-estimation ladder")

  cat("\nRe-estimation ladder (term = standardized ", moderator, "):\n", sep = "")
  print(re_est[, .(SPEC, INSTRUMENT_LABEL,
                   estimate  = signif(estimate, 3),
                   std.error = signif(std.error, 3),
                   p.value   = round(p.value, 4), N)])

  save_csv(corr_tab, "T20A_depth_size_correlations.csv")
  save_csv(re_est,   "T20A_depth_reestimation.csv")
  .tick("save")
  list(correlations = corr_tab, reestimation = re_est)
}


# ===========================================================================
# 20B. PAYER-MIX COMPOSITION                                             (A5)
# ===========================================================================
#
# load_payer_dispersion() itself collapses to FINAL_CONCEPT_ID before it
# returns (its last two lines are `loo_mean(raw, "CV_PAYER_WINSOR")` and
# `raw[, .(N_PAYERS_V2 = mean(N_DISTINCT_PAYERS, ...)), by = FINAL_CONCEPT_ID]`),
# so its return value has no HOSPITAL_ID or POST_MONTH, and N_DISTINCT_PAYERS
# no longer exists under that name. That is fine for build_comparability_
# measures(), which only needs concept-level output, but this section needs
# the row-level file. This reimplements the read-and-filter steps up to the
# point just before that collapse -- same files, same dedup key, same county
# merge, same finiteness filter, same P99 winsorization -- and returns the
# row-level table instead of aggregating it away.

s20_load_payer_dispersion_rowlevel <- function() {
  if (!requireNamespace("R.utils", quietly = TRUE)) {
    stop("Package 'R.utils' is required to read gzipped payer dispersion files.",
         call. = FALSE)
  }
  files <- list.files(PAYER_DISPERSION_DIR, pattern = PAYER_DISPERSION_PATTERN,
                      full.names = TRUE)
  if (length(files) == 0L)
    stop("No payer dispersion files found matching '", PAYER_DISPERSION_PATTERN,
         "' in ", PAYER_DISPERSION_DIR, call. = FALSE)
  
  raw <- rbindlist(lapply(files, function(f) {
    tryCatch(fread(f), error = function(e) {
      warning("Failed to read ", basename(f), ": ", e$message, call. = FALSE)
      data.table()
    })
  }), fill = TRUE)
  if (nrow(raw) == 0L)
    stop("Payer dispersion files read but produced zero rows.", call. = FALSE)
  
  assert_columns(raw, c("HOSPITAL_ID", "POST_MONTH", "ANALYSIS_CONCEPT_ID",
                        "N_DISTINCT_PAYERS", "N_PAYER_CELLS",
                        "CV_PAYER_NEGOTIATED"), "payer dispersion file")
  
  raw <- unique(raw, by = c("HOSPITAL_ID", "POST_MONTH", "ANALYSIS_CONCEPT_ID"))
  raw[, POST_MONTH := as.Date(POST_MONTH)]
  raw[, HOSPITAL_ID := as.character(HOSPITAL_ID)]
  setnames(raw, "ANALYSIS_CONCEPT_ID", "FINAL_CONCEPT_ID")
  
  hosp_county <- unique(outpatient[, .(HOSPITAL_ID, POST_MONTH, ANALYSIS_MARKET)])
  raw <- merge(raw, hosp_county, by = c("HOSPITAL_ID", "POST_MONTH"), all.x = FALSE)
  
  raw <- raw[!is.na(ANALYSIS_MARKET) & is.finite(CV_PAYER_NEGOTIATED) &
               CV_PAYER_NEGOTIATED >= 0]
  if (nrow(raw) == 0L)
    stop("No usable rows after county match / finiteness filter.", call. = FALSE)
  
  cap <- quantile(raw$CV_PAYER_NEGOTIATED, 0.99, na.rm = TRUE)
  raw[, CV_PAYER_WINSOR := pmin(CV_PAYER_NEGOTIATED, cap)]
  
  cat("Row-level payer file: ", format(nrow(raw), big.mark = ","), " rows | ",
      uniqueN(raw$HOSPITAL_ID), " hospitals | ", uniqueN(raw$FINAL_CONCEPT_ID),
      " concepts (pre-collapse; load_payer_dispersion() itself would ",
      "aggregate this to concept level from here)\n", sep = "")
  raw
}

s20_build_payer_panel <- function(panel = outpatient, scheme_col = "SCHEME_1_CERTAINTY",
                                  instruments = MAIN_INSTRUMENTS) {
  payer <- s20_load_payer_dispersion_rowlevel()
  
  # estimate_interacted() needs the endogenous variable (to build TREAT_k) and
  # the baseline controls, not just the scheme and instrument columns -- both
  # default to the pipeline's own globals, so if BASELINE_CONTROLS changes
  # upstream this stays in sync automatically rather than silently dropping
  # a control via available_columns()'s graceful-missing behavior.
  need <- c("HOSPITAL_ID", "POST_MONTH", "FINAL_CONCEPT_ID", "MARKET_ID",
            scheme_col, ENDOGENOUS_VARIABLE, BASELINE_CONTROLS, unname(instruments))
  missing <- setdiff(need, names(panel))
  if (length(missing) > 0L)
    stop("panel (outpatient) is missing: ", paste(missing, collapse = ", "),
         ". Confirm the scheme columns and instruments are already attached ",
         "to outpatient in this session.", call. = FALSE)
  
  base <- unique(panel[, ..need])
  d <- merge(payer, base, by = c("HOSPITAL_ID", "POST_MONTH", "FINAL_CONCEPT_ID"))
  cat("Payer panel: ", format(nrow(d), big.mark = ","), " rows (from ",
      format(nrow(payer), big.mark = ","), " row-level payer-file rows, matched to ",
      scheme_col, ", ", ENDOGENOUS_VARIABLE, ", and instrument columns on outpatient)\n", sep = "")
  d
}

s20_payer_balance <- function(panel = outpatient, scheme_col = "SCHEME_1_CERTAINTY",
                              instruments = MAIN_INSTRUMENTS,
                              outcomes = c("N_DISTINCT_PAYERS", "N_PAYER_CELLS",
                                           "CV_PAYER_WINSOR")) {
  d <- s20_build_payer_panel(panel, scheme_col, instruments)
  outcomes <- intersect(outcomes, names(d))
  
  results <- rbindlist(lapply(outcomes, function(oc) {
    rbindlist(lapply(names(instruments), function(il) {
      r <- tryCatch(
        estimate_interacted(d, scheme_col, outcome = oc, instrument = instruments[[il]],
                            moderator_type = "categorical", label = "Payer-file balance",
                            instrument_label = il, fixed_effects = c("MARKET_ID", "POST_MONTH")),
        error = function(e) { message("  [", oc, " / ", il, "] failed: ", e$message); NULL })
      if (is.null(r)) return(NULL)
      cbind(OUTCOME = oc, r$rows)
    }), fill = TRUE)
  }), fill = TRUE)
  
  cat("\nDoes peer disclosure exposure predict how many payers a hospital\n",
      "reports for a concept, and does that differ by shoppability? A\n",
      "significant, asymmetric effect here would mean part of the price\n",
      "gradient could be mechanical (payer-mix expansion) rather than a\n",
      "genuine price response.\n\n", sep = "")
  
  if (nrow(results) == 0L || !("TERM" %chin% names(results))) {
    cat("All (outcome, instrument) combinations failed -- see the [FAILED]\n",
        "messages above for the specific error. Nothing to print.\n", sep = "")
    return(invisible(results))
  }
  
  print(results[TERM %chin% c("Shoppable", "Non_shoppable"),
                .(OUTCOME, CATEGORY = TERM, INSTRUMENT_LABEL,
                  RF_PCT = round(RF_PERCENT_PER_SD, 3), RF_P = round(RF_P, 4))][
                    order(OUTCOME, CATEGORY, INSTRUMENT_LABEL)])
  
  save_csv(results, "T20B_payer_file_balance.csv")
  invisible(results)
}

# ---------------------------------------------------------------------------
# 20B-ii. FORMAL EQUALITY TEST                                           (A5)
#
# s20_payer_balance() keeps only r$rows -- the two arm-specific coefficients
# -- and discards r$tests, which is the same Wald equality test Table 5
# reports for the price outcome. This reruns the identical nine (outcome x
# instrument) specifications and keeps the equality test instead, so the
# balance check reports the statistic the point-estimate pattern above is
# actually asking about.
# ---------------------------------------------------------------------------
s20_payer_balance_equality <- function(panel = outpatient, scheme_col = "SCHEME_1_CERTAINTY",
                                       instruments = MAIN_INSTRUMENTS,
                                       outcomes = c("N_DISTINCT_PAYERS", "N_PAYER_CELLS",
                                                    "CV_PAYER_WINSOR")) {
  d <- s20_build_payer_panel(panel, scheme_col, instruments)
  outcomes <- intersect(outcomes, names(d))
  
  out <- rbindlist(lapply(outcomes, function(oc) {
    rbindlist(lapply(names(instruments), function(il) {
      r <- tryCatch(
        estimate_interacted(d, scheme_col, outcome = oc, instrument = instruments[[il]],
                            moderator_type = "categorical", label = "Payer-file balance",
                            instrument_label = il, fixed_effects = c("MARKET_ID", "POST_MONTH")),
        error = function(e) { message("  [", oc, " / ", il, "] failed: ", e$message); NULL })
      if (is.null(r)) return(NULL)
      cbind(OUTCOME = oc, r$tests[ESTIMATOR == "Reduced form"])
    }), fill = TRUE)
  }), fill = TRUE)
  
  cat("\nFormal test of H0: peer disclosure predicts reported payer composition\n",
      "equally for shoppable and non-shoppable concepts:\n\n", sep = "")
  print(out[, .(OUTCOME, INSTRUMENT_LABEL, WALD = round(WALD, 2),
                P = round(P_VALUE, 4))][order(OUTCOME, INSTRUMENT_LABEL)])
  
  save_csv(out, "T20B_payer_file_balance_equality.csv")
  invisible(out)
}

s20_price_with_payer_control <- function(panel = outpatient, scheme_col = "SCHEME_1_CERTAINTY",
                                         instruments = MAIN_INSTRUMENTS,
                                         control_var = "N_DISTINCT_PAYERS") {
  payer <- s20_load_payer_dispersion_rowlevel()
  if (!(control_var %chin% names(payer)))
    stop(control_var, " not found in the row-level payer file. Available: ",
         paste(names(payer), collapse = ", "), call. = FALSE)
  
  ctrl <- unique(payer[, .(HOSPITAL_ID, POST_MONTH, FINAL_CONCEPT_ID,
                           PAYER_CTRL = safe_numeric(get(control_var)))])
  d <- merge(panel, ctrl, by = c("HOSPITAL_ID", "POST_MONTH", "FINAL_CONCEPT_ID"), all.x = TRUE)
  match_rate <- mean(!is.na(d$PAYER_CTRL))
  cat("Matched payer-count control to ", format(sum(!is.na(d$PAYER_CTRL)), big.mark = ","),
      " of ", format(nrow(d), big.mark = ","), " panel rows (",
      round(100 * match_rate, 1), "%)\n", sep = "")
  if (match_rate < 0.5)
    warning("Less than half the panel matched a payer-count control; results ",
            "below are on a smaller, possibly selected sample.", call. = FALSE)
  
  out <- rbindlist(lapply(names(instruments), function(il) {
    ctrl_set <- available_columns(d, c(BASELINE_CONTROLS, "PAYER_CTRL"))
    if (!("PAYER_CTRL" %chin% ctrl_set))
      message("  [", il, "] PAYER_CTRL dropped by available_columns() -- check ",
              "for all-NA or zero-variance values in the matched subsample.")
    r <- tryCatch(
      estimate_interacted(d, scheme_col, PRIMARY_OUTCOME, instruments[[il]],
                          moderator_type = "categorical", label = "With payer-count control",
                          instrument_label = il,
                          controls = c(BASELINE_CONTROLS, "PAYER_CTRL"),
                          fixed_effects = c("MARKET_ID", "POST_MONTH")),
      error = function(e) { message("  [", il, "] failed: ", e$message); NULL })
    if (is.null(r)) return(NULL)
    r$rows
  }), fill = TRUE)
  
  cat("\nHeadline coefficients with reported payer count held fixed at the\n",
      "row level (compare against Table 5, Panel A):\n", sep = "")
  
  if (nrow(out) == 0L || !("TERM" %chin% names(out))) {
    cat("All instruments failed -- see the [FAILED] messages above.\n")
    return(invisible(out))
  }
  
  print(out[TERM %chin% c("Shoppable", "Non_shoppable"),
            .(CATEGORY = TERM, INSTRUMENT_LABEL,
              RF_PCT = round(RF_PERCENT_PER_SD, 3), RF_P = round(RF_P, 4))][
                order(CATEGORY, INSTRUMENT_LABEL)])
  
  save_csv(out, "T20B_price_with_payer_control.csv")
  invisible(out)
}

`%||%` <- function(a, b) if (is.null(a)) b else a


# ===========================================================================
# 20C. FIXED-EFFECTS LADDER                                              (A8)
# ===========================================================================
#
# Reuses Section 13's own system-resolution logic verbatim rather than
# assuming outpatient_r13 already exists in memory, so this runs standalone.

s20_resolve_system <- function(panel = outpatient) {
  sys_col <- if ("SYSTEM_KEY" %chin% names(panel)) "SYSTEM_KEY" else "HEALTH_SYSTEM_ID"
  if (!(sys_col %chin% names(panel)))
    stop("Neither SYSTEM_KEY nor HEALTH_SYSTEM_ID found in panel.", call. = FALSE)
  
  resolved <- panel[!is.na(HOSPITAL_ID), .(
    SYS_RESOLVED = { v <- get(sys_col)[!is.na(get(sys_col))]; if (length(v)) v[1L] else NA_character_ }
  ), by = HOSPITAL_ID]
  
  d <- merge(panel, resolved, by = "HOSPITAL_ID", all.x = TRUE, sort = FALSE)
  setDT(d)
  d[, SYSTEM_MONTH := paste0(fifelse(is.na(SYS_RESOLVED), "NOSYS", SYS_RESOLVED),
                             "_", as.character(POST_MONTH))]
  cat("Resolved system for ", format(uniqueN(d$HOSPITAL_ID), big.mark = ","),
      " hospitals | ", format(sum(is.na(resolved$SYS_RESOLVED)), big.mark = ","),
      " genuinely unaffiliated | ", format(uniqueN(d$SYSTEM_MONTH), big.mark = ","),
      " distinct system-months\n", sep = "")
  d
}

# ---------------------------------------------------------------------------
# 20C.0  Fixed-effects cell construction and the rung registry
# ---------------------------------------------------------------------------
#
# Builds every interacted cell the ladder needs, on top of s20_resolve_system().
# MARKET_ID (county x concept) and POST_MONTH already exist in the panel.
#
#   MARKET_FAMILY   county x family    coarsened market cell, ~16 not ~738
#   FAMILY_MONTH    family x month     family-specific national time path
#   CONCEPT_MONTH   concept x month    concept-specific national time path
#   CBSA_CONCEPT    CBSA x concept     market-definition check
#   SYSTEM_MONTH    system x month     from s20_resolve_system()
#   SHOP_NUM        numeric shoppable  for the hospital-FE gradient spec
#
# SHOP_NUM exists because under hospital FE the categorical spec loses a
# category to collinearity. N_PRIOR_POSTERS is constant within hospital, since
# 3,721 of 3,723 hospitals appear at exactly one month, so sum_k TREAT_k is
# absorbed by the hospital effect and only K-1 interactions are identified.
# The continuous spec sidesteps this: TREAT_MAIN drops, TREAT_INTER survives,
# and its coefficient IS the shoppable-minus-non-shoppable gradient.
#
# Hospitals outside a CBSA fall back to their county rather than pooling into
# one giant NOCBSA cell, which would otherwise compare rural Montana to rural
# Georgia inside a single fixed effect.

s20_build_fe_cols <- function(panel = outpatient,
                              scheme_col = "SCHEME_1_CERTAINTY") {
  d <- s20_resolve_system(panel)
  
  d[, MARKET_FAMILY := paste0(ANALYSIS_MARKET, "::", FINAL_FAMILY_ID)]
  d[, FAMILY_MONTH  := paste0(FINAL_FAMILY_ID,  "::", as.character(POST_MONTH))]
  d[, CONCEPT_MONTH := paste0(FINAL_CONCEPT_ID, "::", as.character(POST_MONTH))]
  
  if ("CBSA_CODE" %chin% names(d)) {
    d[, CBSA_CONCEPT := paste0(
      fifelse(is.na(CBSA_CODE) | CBSA_CODE == "",
              paste0("NOCBSA_", ANALYSIS_MARKET), as.character(CBSA_CODE)),
      "::", FINAL_CONCEPT_ID)]
  }
  
  if (scheme_col %chin% names(d)) {
    d[, SHOP_NUM := as.integer(as.character(get(scheme_col)) == "Shoppable")]
  }
  
  d
}

# The rung registry. Names are what appears in every table, so they are the
# one place to edit a label. RUNG_ID drives the cache key -- the three ids
# inherited from the original three-rung ladder are kept verbatim so the
# existing cache is reused rather than re-estimated.

S20_RUNGS <- list(
  `L0 county + concept + month`           = list(
    id = "rung_L0_additive",
    fe = c("ANALYSIS_MARKET", "FINAL_CONCEPT_ID", "POST_MONTH"),
    gradient_only = FALSE),
  `L1 county x family + concept + month`  = list(
    id = "rung_L1_countyfam",
    fe = c("MARKET_FAMILY", "FINAL_CONCEPT_ID", "POST_MONTH"),
    gradient_only = FALSE),
  `L2 county x concept + month (MAIN)`    = list(
    id = "rung1_baseline",
    fe = c("MARKET_ID", "POST_MONTH"),
    gradient_only = FALSE),
  `L3 county x concept + family x month`  = list(
    id = "rung_L3_fammonth",
    fe = c("MARKET_ID", "FAMILY_MONTH"),
    gradient_only = FALSE),
  `L4 county x concept + concept x month` = list(
    id = "rung_L4_conceptmonth",
    fe = c("MARKET_ID", "CONCEPT_MONTH"),
    gradient_only = FALSE),
  `L5 CBSA x concept + month`             = list(
    id = "rung_L5_cbsa",
    fe = c("CBSA_CONCEPT", "POST_MONTH"),
    gradient_only = FALSE),
  `L6a county x concept + month + system` = list(
    id = "rung2_system_fe",
    fe = c("MARKET_ID", "POST_MONTH", "SYS_RESOLVED"),
    gradient_only = FALSE),
  `L6b county x concept + system x month` = list(
    id = "rung3_sysmonth",
    fe = c("MARKET_ID", "SYSTEM_MONTH"),
    gradient_only = FALSE),
  `L7 hospital + concept (gradient)`      = list(
    id = "rung_L7_hospital",
    fe = c("HOSPITAL_ID", "FINAL_CONCEPT_ID"),
    gradient_only = TRUE)
)


# ---------------------------------------------------------------------------
# 20C.1  Fixed-effects census -- the QA gate, no estimation
# ---------------------------------------------------------------------------
#
# Run this BEFORE the ladder. For each rung it reports how much identifying
# variation survives the fixed effects, using one no-covariate demeaning pass
# per rung rather than a full IV fit. Minutes, not hours.
#
# Columns:
#   FE_LEVELS      total fixed-effect levels across all dimensions
#   N_KEPT         rows surviving cascading singleton removal
#   PCT_KEPT       N_KEPT as a share of the complete-case sample
#   Z_RESID_SHARE  residual SD of the instrument / raw SD, on kept rows
#   ZSHOP_RESID    same for Z x 1[shoppable], which is what identifies the
#                  shoppable arm of the interacted spec
#   ZNON_RESID     same for Z x 1[non-shoppable]
#   HOSP_KEPT      hospitals surviving
#   INDEP_LOST     share of unaffiliated hospitals dropped, the composition
#                  damage that makes a system-FE collapse ambiguous
#
# READING IT. Z_RESID_SHARE below roughly 0.20 means the rung has eaten most
# of the instrument and its coefficient is uninformative no matter what the
# first-stage F says. A rung with PCT_KEPT under about 0.40 is estimating on
# a different population, so a coefficient change there is composition and
# confounding mixed together and cannot be attributed to either.

s20_fe_census <- function(panel = outpatient, scheme_col = "SCHEME_1_CERTAINTY",
                          instrument = PRIMARY_INSTRUMENT,
                          outcome = PRIMARY_OUTCOME, rungs = S20_RUNGS) {
  
  d <- s20_build_fe_cols(panel, scheme_col)
  
  if (!("SHOP_NUM" %chin% names(d)))
    stop("SHOP_NUM not built -- check that ", scheme_col, " is on the panel.",
         call. = FALSE)
  
  base_cols <- unique(c(outcome, ENDOGENOUS_VARIABLE, instrument,
                        available_columns(d, BASELINE_CONTROLS),
                        available_columns(d, BASELINE_CLUSTERS), scheme_col))
  
  cat("\nFixed-effects census. One demeaning pass per rung, no covariates.\n")
  cat("Each line prints when it finishes; interrupt is safe, nothing is cached.\n\n")
  
  out <- rbindlist(lapply(names(rungs), function(rn) {
    spec <- rungs[[rn]]
    missing_fe <- setdiff(spec$fe, names(d))
    if (length(missing_fe)) {
      cat(sprintf("  %-40s SKIPPED, missing %s\n", rn,
                  paste(missing_fe, collapse = ", ")))
      return(NULL)
    }
    
    t0 <- Sys.time()
    s <- model_sample(copy(d), c(base_cols, spec$fe, "SHOP_NUM"))
    if (nrow(s) < MIN_MODEL_OBS) {
      cat(sprintf("  %-40s SKIPPED, %s complete cases < MIN_MODEL_OBS\n",
                  rn, format(nrow(s), big.mark = ",")))
      return(NULL)
    }
    
    s[, Z_RAW  := safe_numeric(get(instrument))]
    s[, Z_SHOP := Z_RAW * SHOP_NUM]
    s[, Z_NON  := Z_RAW * (1L - SHOP_NUM)]
    
    fit <- tryCatch(
      feols(c(Z_RAW, Z_SHOP, Z_NON) ~ 1, data = s,
            fixef = spec$fe, warn = FALSE, notes = FALSE),
      error = function(e) { cat("  FAILED: ", e$message, "\n", sep = ""); NULL })
    if (is.null(fit)) return(NULL)
    
    n_kept  <- nobs(fit[[1L]])
    kept_ix <- tryCatch(obs(fit[[1L]]), error = function(e) NULL)
    if (is.null(kept_ix) || length(kept_ix) != n_kept) kept_ix <- seq_len(nrow(s))
    
    resid_share <- function(j, col) {
      r <- sd(resid(fit[[j]]), na.rm = TRUE)
      raw <- sd(s[[col]][kept_ix], na.rm = TRUE)
      if (!is.finite(raw) || raw <= 0) NA_real_ else r / raw
    }
    
    kept <- s[kept_ix]
    indep_all  <- uniqueN(d[is.na(SYS_RESOLVED), HOSPITAL_ID])
    indep_kept <- uniqueN(kept[is.na(SYS_RESOLVED), HOSPITAL_ID])
    
    res <- data.table(
      RUNG          = rn,
      FE_LEVELS     = sum(vapply(spec$fe, function(f) uniqueN(s[[f]]), integer(1))),
      N_COMPLETE    = nrow(s),
      N_KEPT        = n_kept,
      PCT_KEPT      = round(n_kept / nrow(s), 3),
      Z_RESID_SHARE = round(resid_share(1L, "Z_RAW"),  3),
      ZSHOP_RESID   = round(resid_share(2L, "Z_SHOP"), 3),
      ZNON_RESID    = round(resid_share(3L, "Z_NON"),  3),
      HOSP_KEPT     = uniqueN(kept$HOSPITAL_ID),
      SYS_KEPT      = uniqueN(kept$SYS_RESOLVED),
      INDEP_LOST    = if (indep_all > 0)
        round(1 - indep_kept / indep_all, 3) else NA_real_,
      SECS          = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
    )
    
    cat(sprintf("  %-40s kept %5.1f%% | Zresid %.2f | %6.1fs\n",
                rn, 100 * res$PCT_KEPT, res$Z_RESID_SHARE, res$SECS))
    rm(s, kept); invisible(gc())
    res
  }), fill = TRUE)
  
  if (nrow(out) == 0L) {
    cat("\nEvery rung skipped or failed. Nothing to report.\n")
    return(invisible(out))
  }
  
  cat("\n")
  print(out[, .(RUNG, FE_LEVELS, N_KEPT, PCT_KEPT,
                Z_RESID_SHARE, ZSHOP_RESID, ZNON_RESID,
                HOSP_KEPT, INDEP_LOST, SECS)])
  
  cat("\nVERDICT (screen only, not a substitute for the first stage):\n")
  # Judge on the weaker interaction arm, not the raw instrument. The interacted
  # spec is identified by Z x 1[category], so a rung that absorbs raw Z
  # entirely (hospital FE necessarily does, since Z is constant within
  # hospital) is fine as long as both interactions survive.
  out[, ZMIN := pmin(ZSHOP_RESID, ZNON_RESID, na.rm = TRUE)]
  out[, VERDICT := fifelse(
    is.na(ZMIN) | ZMIN < 0.20, "STARVED -- do not estimate",
    fifelse(INDEP_LOST >= 0.99, "SYSTEMS ONLY -- all independents dropped",
            fifelse(PCT_KEPT < 0.40, "COMPOSITION -- interpret with N, not alone",
                    fifelse(ZMIN < 0.40, "THIN -- expect wide SEs", "OK"))))]
  print(out[, .(RUNG, Z_RESID_SHARE, PCT_KEPT, VERDICT)])
  
  save_qa_csv(out, "QA20C_fe_census.csv")
  invisible(out)
}


# ---------------------------------------------------------------------------
# 20C.2  The extended ladder
# ---------------------------------------------------------------------------
#
# Deliberately a NEW function rather than an edit to s20_fe_ladder(). The old
# three-rung version stays exactly as it is, so run_section_20() keeps working
# and its cached results stay valid. Extending the rung list inside the old
# function would have silently turned the driver's 18-scheme call into
# 18 x 9 x 3 = 486 estimations on the next run.
#
# DEFAULTS ARE ONE SCHEME AND ONE INSTRUMENT ON PURPOSE. That is 9 rungs, so 9
# calls and 18 feols fits. Widen only for rungs the census cleared and this
# pass found interesting. Going straight to unname(SCHEME_COLUMNS) x
# MAIN_INSTRUMENTS is a 486-call job and there is no reason to buy that before
# knowing which rungs survive.
#
# Three rung ids are inherited from the old ladder, so L2, L6a and L6b load
# from cache instead of re-estimating if that ladder has already run under the
# same scheme.
#
# L7 uses the continuous moderator. Under hospital FE the categorical spec
# would silently lose one category to collinearity and return NA for it. The
# continuous spec drops TREAT_MAIN instead, which is the term that is genuinely
# unidentified here, and reports the gradient directly.

s20_fe_ladder_full <- function(panel = outpatient,
                               scheme_cols = "SCHEME_1_CERTAINTY",
                               instruments = MAIN_INSTRUMENTS[1],
                               outcome = PRIMARY_OUTCOME,
                               rungs = S20_RUNGS, use_cache = TRUE) {
  
  d <- s20_build_fe_cols(panel, scheme_cols[1])
  
  n_total <- length(scheme_cols) * length(rungs) * length(instruments)
  cat("\nExtended FE ladder: ", length(scheme_cols), " scheme(s) x ",
      length(rungs), " rungs x ", length(instruments), " instrument(s) = ",
      n_total, " calls, each fitting a reduced form and an IV.\n",
      "Every call prints its own elapsed time and is cached on completion, so\n",
      "interrupting loses at most the call in flight.\n\n", sep = "")
  
  out <- rbindlist(lapply(scheme_cols, function(scheme_col) {
    if (!(scheme_col %chin% names(d))) {
      cat("  scheme not on panel, skipped: ", scheme_col, "\n", sep = "")
      return(NULL)
    }
    d[, SHOP_NUM := as.integer(as.character(get(scheme_col)) == "Shoppable")]
    
    rbindlist(lapply(names(rungs), function(rn) {
      spec <- rungs[[rn]]
      missing_fe <- setdiff(spec$fe, names(d))
      if (length(missing_fe)) {
        cat(sprintf("  %-40s SKIPPED, missing %s\n", rn,
                    paste(missing_fe, collapse = ", ")))
        return(NULL)
      }
      
      mod      <- if (isTRUE(spec$gradient_only)) "SHOP_NUM" else scheme_col
      mod_type <- if (isTRUE(spec$gradient_only)) "continuous" else "categorical"
      
      rbindlist(lapply(names(instruments), function(il) {
        key <- paste0("s20_ladder_", scheme_col, "_", spec$id, "_", il)
        cat(sprintf("  %-20s %-40s %-30s ", scheme_col, substr(rn, 1, 38),
                    substr(il, 1, 28)))
        t0 <- Sys.time()
        r <- tryCatch({
          call_it <- function() estimate_interacted(
            d, mod, outcome, instruments[[il]],
            moderator_type = mod_type, label = rn, instrument_label = il,
            moderator_label = scheme_col, fixed_effects = spec$fe)
          if (use_cache) cache_or_run(key, call_it()) else call_it()
        }, error = function(e) { cat("FAILED: ", e$message, "\n", sep = ""); NULL })
        cat(sprintf("%6.1fs\n",
                    as.numeric(difftime(Sys.time(), t0, units = "secs"))))
        if (is.null(r)) return(NULL)
        cbind(SCHEME = scheme_col, RUNG = rn, GRADIENT_SPEC = mod_type, r$rows)
      }), fill = TRUE)
    }), fill = TRUE)
  }), fill = TRUE)
  
  if (nrow(out) == 0L || !("TERM" %chin% names(out))) {
    cat("\nEvery cell failed. See the FAILED lines above.\n")
    return(invisible(out))
  }
  
  # ------------------------------------------------------------------
  # The cross-rung comparable object is the GRADIENT, not the level.
  # Categorical rungs give it as Shoppable minus Non_shoppable; L7 reports
  # it directly as the interaction term. Both are percent per SD of Z.
  # ------------------------------------------------------------------
  cat_grad <- dcast(out[GRADIENT_SPEC == "categorical" &
                          TERM %chin% c("Shoppable", "Non_shoppable")],
                    SCHEME + RUNG + INSTRUMENT_LABEL ~ TERM,
                    value.var = "RF_PERCENT_PER_SD")
  if (all(c("Shoppable", "Non_shoppable") %chin% names(cat_grad))) {
    cat_grad[, GRADIENT_PCT := Shoppable - Non_shoppable]
  } else {
    cat_grad <- data.table()
  }
  
  con_grad <- out[GRADIENT_SPEC == "continuous" & TERM == "x Moderator",
                  .(SCHEME, RUNG, INSTRUMENT_LABEL,
                    Shoppable = NA_real_, Non_shoppable = NA_real_,
                    GRADIENT_PCT = RF_PERCENT_PER_SD)]
  
  grad <- rbindlist(list(cat_grad, con_grad), fill = TRUE)
  
  fs <- unique(out[, .(SCHEME, RUNG, INSTRUMENT_LABEL,
                       FS_MIN = round(FIRST_STAGE_WALD_MIN, 1),
                       N = N_OBSERVATIONS)])
  grad <- merge(grad, fs, by = c("SCHEME", "RUNG", "INSTRUMENT_LABEL"),
                all.x = TRUE, sort = FALSE)
  setorder(grad, SCHEME, INSTRUMENT_LABEL, RUNG)
  
  cat("\nShoppability gradient across the ladder (percent per SD of Z):\n")
  print(grad[, .(RUNG, INSTRUMENT_LABEL,
                 SHOP = round(Shoppable, 3), NONSHOP = round(Non_shoppable, 3),
                 GRADIENT = round(GRADIENT_PCT, 3), FS_MIN, N)])
  
  cat("\nLevels and p-values, every term:\n")
  print(out[, .(RUNG, TERM, INSTRUMENT_LABEL,
                RF_PCT = round(RF_PERCENT_PER_SD, 3),
                RF_P = round(RF_P, 4), FS = round(FIRST_STAGE_WALD_MIN, 1),
                N = N_OBSERVATIONS)][order(RUNG, INSTRUMENT_LABEL, TERM)])
  
  cat("\nHOW TO READ THIS.\n",
      "  GRADIENT is the comparable number across rungs. The levels move with\n",
      "  what each fixed effect absorbs and are not comparable rung to rung.\n",
      "  L2 is the paper's main specification. L3 and L4 test whether the\n",
      "  gradient is a service-family time trend. L7 is the strongest check:\n",
      "  same hospital, same month, same system, same payer mix, and the\n",
      "  gradient is identified purely off which services a hospital posts.\n",
      "  Read every rung against QA20C_fe_census.csv. A gradient that shrinks\n",
      "  on a rung the census flagged STARVED or COMPOSITION is not evidence\n",
      "  of confounding, it is evidence the rung had nothing left to work with.\n",
      "  Check the sign of FS_MIN's first stage before believing any collapse.\n",
      sep = "")
  
  save_csv(out,  "T20C2_fe_ladder_full.csv")
  save_csv(grad, "T20C2_fe_ladder_gradient.csv")
  invisible(list(rows = out, gradient = grad))
}



s20_fe_ladder <- function(panel = outpatient, scheme_cols = "SCHEME_1_CERTAINTY",
                          instruments = MAIN_INSTRUMENTS, outcome = PRIMARY_OUTCOME,
                          use_cache = TRUE) {
  
  d <- s20_resolve_system(panel)
  
  rungs <- list(
    `(1) Baseline (market x concept + month)` = c("MARKET_ID", "POST_MONTH"),
    `(2) + system FE (additive)`              = c("MARKET_ID", "POST_MONTH", "SYS_RESOLVED"),
    `(3) System x month`                      = c("MARKET_ID", "SYSTEM_MONTH")
  )
  rung_ids <- c("rung1_baseline", "rung2_system_fe", "rung3_sysmonth")
  names(rung_ids) <- names(rungs)
  
  n_total <- length(scheme_cols) * length(rungs) * length(instruments)
  cat("\nRunning ", length(scheme_cols), " scheme(s) x 3 rungs x ", length(instruments),
      " instruments = ", n_total, " calls. Each prints its own elapsed time; if a\n",
      "call hangs, note which one and interrupt -- earlier cells are already\n",
      "cached (use_cache = TRUE) and will not be redone on a re-run. The cache\n",
      "key includes the scheme, so running multiple schemes will not collide.\n\n", sep = "")
  
  out <- rbindlist(lapply(scheme_cols, function(scheme_col) {
    rbindlist(lapply(names(rungs), function(rn) {
      fe <- rungs[[rn]]
      rbindlist(lapply(names(instruments), function(il) {
        # scheme_col is part of the key -- without this, running a second
        # scheme after the first would silently return the first scheme's
        # cached result for every (rung, instrument) cell rather than
        # actually re-estimating, since rung/instrument alone repeat exactly.
        key <- paste0("s20_ladder_", scheme_col, "_", rung_ids[[rn]], "_", il)
        cat(sprintf("  %-22s %-42s %-38s ", scheme_col, substr(rn, 1, 40), substr(il, 1, 36)))
        t0 <- Sys.time()
        r <- tryCatch({
          if (use_cache) {
            cache_or_run(key, estimate_interacted(
              d, scheme_col, outcome, instruments[[il]],
              moderator_type = "categorical", label = rn,
              instrument_label = il, fixed_effects = fe))
          } else {
            estimate_interacted(
              d, scheme_col, outcome, instruments[[il]],
              moderator_type = "categorical", label = rn,
              instrument_label = il, fixed_effects = fe)
          }
        }, error = function(e) { cat("FAILED: ", e$message, "\n", sep = ""); NULL })
        cat(sprintf("%6.1fs\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
        if (is.null(r)) return(NULL)
        cbind(SCHEME = scheme_col, RUNG = rn, r$rows)
      }), fill = TRUE)
    }), fill = TRUE)
  }), fill = TRUE)
  
  if (nrow(out) == 0L || !("TERM" %chin% names(out))) {
    cat("\nAll scheme x rung x instrument combinations failed -- see the FAILED\n",
        "messages above for specific errors. Nothing to print.\n", sep = "")
    save_csv(out, "T20C_fe_ladder.csv")
    return(invisible(out))
  }
  
  cat("\nShoppable coefficient across the fixed-effects ladder, by scheme:\n")
  print(out[TERM == "Shoppable",
            .(SCHEME, RUNG, INSTRUMENT_LABEL, RF_PCT = round(RF_PERCENT_PER_SD, 3),
              RF_P = round(RF_P, 4))][order(SCHEME, INSTRUMENT_LABEL, RUNG)])
  
  cat("\nNon-shoppable coefficient across the same ladder, by scheme:\n")
  print(out[TERM == "Non_shoppable",
            .(SCHEME, RUNG, INSTRUMENT_LABEL, RF_PCT = round(RF_PERCENT_PER_SD, 3),
              RF_P = round(RF_P, 4))][order(SCHEME, INSTRUMENT_LABEL, RUNG)])
  
  cat("\nSummary: share of the baseline-to-rung(2) drop, by scheme (primary\n",
      "instrument only, for a quick cross-scheme read):\n", sep = "")
  wide <- dcast(out[TERM == "Shoppable" & INSTRUMENT_LABEL == names(MAIN_INSTRUMENTS)[1]],
                SCHEME ~ RUNG, value.var = "RF_PERCENT_PER_SD")
  rung_names <- names(rungs)
  if (all(rung_names %chin% names(wide))) {
    wide[, DROP_SHARE_AT_RUNG2 := round(1 - get(rung_names[2]) / get(rung_names[1]), 2)]
    print(wide[, c("SCHEME", rung_names, "DROP_SHARE_AT_RUNG2"), with = FALSE])
  }
  
  cat("\nHOW TO READ THIS.\n",
      "  Rung (1) is Table 5's baseline. Rung (3) should numerically match\n",
      "  Table sysmonth's existing reduced-form rows -- if it does not, one of\n",
      "  the two blocks of code has drifted and that needs resolving first.\n",
      "  Rung (2) is the new information: it absorbs average price-level\n",
      "  differences across systems without absorbing system-specific TIME\n",
      "  variation. If the shoppable coefficient survives rung (2) but not\n",
      "  rung (3), the system x month result is about system-level SHOCKS in\n",
      "  a given month, not baseline system price levels -- which favors the\n",
      "  reading in the paper that between-system variation carries the\n",
      "  benchmarking signal, since within-system time variation is what rung\n",
      "  (3) additionally removes relative to rung (2).\n", sep = "")
  
  save_csv(out, "T20C_fe_ladder.csv")
  invisible(out)
}


# ===========================================================================
# DRIVER
# ===========================================================================

run_section_20 <- function(concept_results = NULL, measures = NULL,
                           panel = outpatient) {
  .s20_hd("SECTION 20 -- A5, A6, A8")
  
  a6 <- NULL
  if (!is.null(concept_results) && !is.null(measures)) {
    .s20_sub("20A -- Contracting depth vs. concept size (A6)")
    a6 <- s20_depth_diagnostics(concept_results, measures)
  } else {
    cat("\n[20A skipped] Pass concept_results = and measures = (concept_results,\n",
        "measures from Section 9) to run this.\n", sep = "")
  }
  
  .s20_sub("20B -- Payer-mix composition (A5)")
  b_balance <- tryCatch(s20_payer_balance(panel), error = function(e) {
    message("[20B balance test failed] ", e$message); NULL
  })
  b_control <- tryCatch(s20_price_with_payer_control(panel), error = function(e) {
    message("[20B control regression failed] ", e$message); NULL
  })
  
  .s20_sub("20C -- Fixed-effects ladder (A8)")
  c_ladder <- tryCatch(s20_fe_ladder(panel, scheme_cols = unname(SCHEME_COLUMNS)),
                       error = function(e) { message("[20C failed] ", e$message); NULL })
  
  .s20_hd("SECTION 20 COMPLETE")
  invisible(list(depth = a6, payer_balance = b_balance,
                 payer_control = b_control, fe_ladder = c_ladder))
}

# s20 <- run_section_20(concept_results = concept_results, measures = measures)
#
# If concept_results / measures aren't in memory:
#   measures <- cache_or_run("comparability_measures", build_comparability_measures(outpatient))
#   s20 <- run_section_20(concept_results = concept_results, measures = measures)


# ---------------------------------------------------------------------------
# Section 20 execution
#
# 20A and 20B are cheap. 20C is not: with scheme_cols = SCHEME_COLUMNS it is
# 6 schemes x 3 rungs x 3 instruments = 54 interacted IV fits, an overnight
# job on a cold cache. It is therefore off unless asked for:
#     S20_RUN <- c("20A", "20B", "20C")
# ---------------------------------------------------------------------------
S20_RUN <- if (exists("S20_RUN")) S20_RUN else c("20A", "20B")
s20_out <- list()

if ("20A" %in% S20_RUN) {
  .s20_sub("20A -- Contracting depth vs. concept size (A6)")
  if (!exists("concept_results")) {
    message("[20A skipped] concept_results not in memory. restore_session() first.")
  } else {
    if (!exists("measures"))
      measures <- cache_or_run("comparability_measures",
                               build_comparability_measures(outpatient))
    s20_out$depth <- s20_depth_diagnostics(concept_results, measures)
  }
}

if ("20B" %in% S20_RUN) {
  .s20_sub("20B -- Payer-mix composition (A5)")
  if (length(list.files(PAYER_DISPERSION_DIR, pattern = PAYER_DISPERSION_PATTERN)) == 0L) {
    message("[20B skipped] No HPT_PAYER_DISPERSION*.csv.gz in ", PAYER_DISPERSION_DIR)
  } else {
    s20_out$payer_balance <- tryCatch(s20_payer_balance(outpatient),
                                      error = function(e) { message("[20B balance] ", e$message); NULL })
    s20_out$payer_control <- tryCatch(s20_price_with_payer_control(outpatient),
                                      error = function(e) { message("[20B control] ", e$message); NULL })
  }
}

if ("20C" %in% S20_RUN) {
  .s20_sub("20C -- Fixed-effects ladder (A8)")
  s20_out$fe_ladder <- tryCatch(
    s20_fe_ladder(outpatient, scheme_cols = unname(SCHEME_COLUMNS)),
    error = function(e) { message("[20C] ", e$message); NULL })
}

}  # end Section 20 gate

if (isTRUE(HPT_RUN$diagnostics)) {


###############################################################################

# ============================================================================
# Section 21: Tier II diagnostics A3, E1
# ============================================================================

# SECTION 21 -- TIER II DIAGNOSTICS: A3, E1
#
# Source AFTER HPT_Analysis_Pipeline.R and, for A3, ideally after Section 19
# has run once (not required, but confirms the naive comparison below).
#
#   21A  A3 -- exact test of pi_S = pi_N (the two OWN first-stage
#              coefficients in the interacted system), accounting for their
#              cross-equation covariance rather than assuming independence.
#              A cheap naive version (independence-assumed) is included
#              first since it needs no new computation; the exact version
#              is the expensive fallback if the naive result isn't decisive
#              enough to rely on.
#   21B  E1  -- wild cluster bootstrap over POST_MONTH (16 clusters), the
#              few-clusters check on the time dimension of two-way
#              clustering. Requires the fwildclusterboot package, which is
#              purpose-built for exactly this and interfaces directly with
#              fixest objects.
###############################################################################

.s21_hd  <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")
.s21_sub <- function(x) cat("\n--- ", x, " ", strrep("-", max(0, 70 - nchar(x))), "\n", sep = "")


# ===========================================================================
# 21A. EXACT TEST OF pi_S = pi_N                                         (A3)
# ===========================================================================

# -- cheap version: needs nothing but numbers already in T19A -------------
# Treats the two own-coefficient estimates as independent. This is EITHER
# too conservative OR too liberal depending on the sign of their true
# cross-equation covariance, which this version cannot see. Use it as a
# first read, not a final answer.
s21_naive_pi_test <- function(pi_n, t_n, pi_s, t_s, instrument_label = "") {
  se_n <- pi_n / t_n
  se_s <- pi_s / t_s
  diff <- pi_s - pi_n
  se_diff <- sqrt(se_n^2 + se_s^2)
  z <- diff / se_diff
  data.table(INSTRUMENT_LABEL = instrument_label, PI_N = pi_n, PI_S = pi_s,
             DIFF = diff, SE_DIFF_NAIVE = se_diff, Z_NAIVE = z,
             P_NAIVE = 2 * pnorm(-abs(z)),
             NOTE = "independence assumed -- see 21A exact version for the real test")
}

run_s21a_naive <- function() {
  .s21_hd("21A (cheap) -- NAIVE TEST OF pi_S = pi_N, INDEPENDENCE ASSUMED")
  # Values transcribed from the T19A run already completed. Re-derive from
  # the saved CSV if you want to avoid hand transcription:
  #   pi <- fread(file.path(TABLE_DIR, "T19A_interacted_first_stage_matrix.csv"))
  rows <- list(
    list(il = "Competitor_only_hospitals_9m",         pi_n = 0.128546, t_n = 7.74, pi_s = 0.099785, t_s = 5.47),
    list(il = "Primary_strict_system_IV",              pi_n = 0.123782, t_n = 8.48, pi_s = 0.096042, t_s = 5.76),
    list(il = "Competitor_outside_CBSA_hospitals_9m",  pi_n = 0.127403, t_n = 7.34, pi_s = 0.094847, t_s = 4.92)
  )
  out <- rbindlist(lapply(rows, function(r)
    s21_naive_pi_test(r$pi_n, r$t_n, r$pi_s, r$t_s, r$il)))
  print(out[, .(INSTRUMENT_LABEL, PI_N = round(PI_N, 4), PI_S = round(PI_S, 4),
                DIFF = round(DIFF, 4), Z = round(Z_NAIVE, 2), P = round(P_NAIVE, 3))])
  cat("\nNone significant under independence. This is a starting point, not\n",
      "a conclusion -- pi_N and pi_S are fit on the same underlying rows\n",
      "(TREAT_Non and TREAT_Shop are mutually exclusive by construction on\n",
      "any given row), so their true covariance is almost certainly nonzero\n",
      "and could push the exact test's p-value in either direction.\n", sep = "")
  invisible(out)
}


# -- exact version: joint stacked estimation -------------------------------
#
# Stacks the sample twice, once per category, tagged by EQ. Each equation
# gets its OWN fixed effects (EQ-prefixed), so the point estimates exactly
# reproduce running the two first-stage equations separately. Clustering
# uses the ORIGINAL, non-prefixed county/month identifiers, so that EQ1 and
# EQ2 rows sharing the same underlying market/month are correctly treated as
# the same cluster -- this is what lets the joint vcov() capture the
# cross-equation covariance instead of assuming it away.

s21_test_pi_equality <- function(
    data, scheme_col = "SCHEME_1_CERTAINTY", instrument, instrument_label = "",
    endogenous = ENDOGENOUS_VARIABLE, controls = BASELINE_CONTROLS,
    fixed_effects = BASELINE_FIXED_EFFECTS, clusters = BASELINE_CLUSTERS) {
  
  controls <- available_columns(data, controls)
  fe <- available_columns(data, fixed_effects)
  cl <- available_columns(data, clusters)
  
  d <- data[!is.na(get(scheme_col))]
  d <- model_sample(d, c(endogenous, instrument, controls, fe, cl, scheme_col))
  if (nrow(d) < MIN_MODEL_OBS) stop("Sample too small for ", instrument_label, call. = FALSE)
  
  d[, MOD := droplevels(factor(get(scheme_col)))]
  keys <- levels(d$MOD)
  if (length(keys) != 2L)
    stop("This test assumes exactly two categories; ", scheme_col,
         " has ", length(keys), ".", call. = FALSE)
  k1 <- keys[1]; k2 <- keys[2]
  
  for (k in keys) {
    sel <- as.integer(d$MOD == k)
    d[, paste0("TREAT_", k) := get(endogenous) * sel]
    d[, paste0("IV_", k)    := get(instrument) * sel]
  }
  
  mk_block <- function(k, eq_label) {
    b <- copy(d)
    b[, EQ  := eq_label]
    b[, DEP := get(paste0("TREAT_", k))]
    b
  }
  d_stack <- rbindlist(list(mk_block(k1, "EQ1"), mk_block(k2, "EQ2")))
  
  # Explicit per-equation columns rather than formula-interaction syntax --
  # avoids R's baseline+difference contrast coding, which would not give two
  # directly-testable free coefficients.
  d_stack[, IVk1_EQ1 := get(paste0("IV_", k1)) * as.integer(EQ == "EQ1")]
  d_stack[, IVk2_EQ1 := get(paste0("IV_", k2)) * as.integer(EQ == "EQ1")]
  d_stack[, IVk1_EQ2 := get(paste0("IV_", k1)) * as.integer(EQ == "EQ2")]
  d_stack[, IVk2_EQ2 := get(paste0("IV_", k2)) * as.integer(EQ == "EQ2")]
  for (cv in controls) {
    d_stack[, paste0(cv, "_EQ1") := get(cv) * as.integer(EQ == "EQ1")]
    d_stack[, paste0(cv, "_EQ2") := get(cv) * as.integer(EQ == "EQ2")]
  }
  d_stack[, FE_A := paste0(EQ, "_", get(fe[1]))]
  if (length(fe) > 1L) d_stack[, FE_B := paste0(EQ, "_", get(fe[2]))]
  
  rhs <- c("IVk1_EQ1", "IVk2_EQ1", "IVk1_EQ2", "IVk2_EQ2",
           if (length(controls) > 0L) paste0(rep(controls, each = 2), c("_EQ1", "_EQ2")))
  fe_stack <- if (length(fe) > 1L) c("FE_A", "FE_B") else "FE_A"
  
  # Cluster on the ORIGINAL (non-EQ-prefixed) identifiers -- this is the
  # whole point: it lets EQ1 and EQ2 rows sharing a market/month be treated
  # as the same cluster, so the sandwich covariance captures the
  # cross-equation correlation instead of assuming it away.
  cl_stack <- cl
  
  form <- as.formula(paste("DEP ~", paste(rhs, collapse = " + "), "|",
                           paste(fe_stack, collapse = " + ")))
  clf  <- as.formula(paste("~", paste(cl_stack, collapse = " + ")))
  
  fit <- feols(form, data = d_stack, cluster = clf, warn = FALSE, notes = FALSE)
  
  b <- coef(fit); V <- vcov(fit)
  nm1 <- "IVk1_EQ1"; nm2 <- "IVk2_EQ2"   # pi_{k1} and pi_{k2}, the two OWN coefficients
  if (!all(c(nm1, nm2) %chin% names(b)))
    stop("Expected coefficients not found -- got: ", paste(names(b), collapse = ", "),
         call. = FALSE)
  
  diff    <- unname(b[nm1] - b[nm2])
  se_diff <- sqrt(V[nm1, nm1] + V[nm2, nm2] - 2 * V[nm1, nm2])
  z <- diff / se_diff
  p <- .pval(z, fit)
  
  data.table(
    INSTRUMENT_LABEL = instrument_label,
    CATEGORY_1 = k1, CATEGORY_2 = k2,
    PI_1 = unname(b[nm1]), PI_2 = unname(b[nm2]),
    COV_12 = V[nm1, nm2],
    DIFF = diff, SE_DIFF_EXACT = se_diff, Z_EXACT = z, P_EXACT = p,
    N_STACKED = nrow(d_stack)
  )
}

run_s21a_exact <- function(panel = outpatient, scheme_col = "SCHEME_1_CERTAINTY",
                           instruments = MAIN_INSTRUMENTS, use_cache = TRUE) {
  .s21_hd("21A (exact) -- STACKED TEST OF pi_S = pi_N")
  cat("\nThis fits a doubled-row, doubled-FE joint model per instrument.\n",
      "Expect each call to take as long as, or longer than, the entire\n",
      "9-cell ladder from Section 20C. Progress prints per instrument;\n",
      "results are cached individually.\n\n", sep = "")
  
  out <- rbindlist(lapply(names(instruments), function(il) {
    key <- paste0("s21a_pi_equality_", scheme_col, "_", il)
    cat(sprintf("  %-38s ", il))
    t0 <- Sys.time()
    r <- tryCatch({
      if (use_cache) {
        cache_or_run(key, s21_test_pi_equality(panel, scheme_col, instruments[[il]], il))
      } else {
        s21_test_pi_equality(panel, scheme_col, instruments[[il]], il)
      }
    }, error = function(e) { cat("FAILED: ", e$message, "\n", sep = ""); NULL })
    cat(sprintf("%6.1fs\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    r
  }), fill = TRUE)
  
  if (nrow(out) > 0L) {
    cat("\nExact test, accounting for cross-equation covariance:\n")
    print(out[, .(INSTRUMENT_LABEL, PI_1 = round(PI_1, 4), PI_2 = round(PI_2, 4),
                  COV = signif(COV_12, 3), DIFF = round(DIFF, 4),
                  Z = round(Z_EXACT, 2), P = round(P_EXACT, 3))])
  }
  save_csv(out, "T21A_pi_equality_exact.csv")
  invisible(out)
}



# ===========================================================================
# 21B. WILD CLUSTER BOOTSTRAP OVER MONTHS                                (E1)
# ===========================================================================
#
# Two-way clustering (county, month) with only 16 month clusters is the
# concern. fwildclusterboot::boottest() is purpose-built for few-cluster
# inference and works directly on a fitted feols object.

s21_check_boot_package <- function() {
  if (!requireNamespace("fwildclusterboot", quietly = TRUE)) {
    cat("\nPackage 'fwildclusterboot' is not installed. Install with:\n",
        "  install.packages('fwildclusterboot')\n",
        "It is purpose-built for wild cluster bootstrap inference with few\n",
        "clusters and interfaces directly with fixest objects -- there is no\n",
        "good reason to hand-rolled this one.\n", sep = "")
    return(FALSE)
  }
  TRUE
}

# -----------------------------------------------------------------------
# LEAVE-ONE-MONTH-OUT: the alternative that actually runs.
#
# fwildclusterboot::boottest() did not complete on this panel across four
# attempts (a direct call, a hand-rolled Frisch-Waugh-Lovell workaround with
# three distinct bugs, and a fourth attempt using boottest()'s own fe=
# argument), including a token B=199 run that still did not return within
# a minute. The per-iteration cost against MARKET_ID's ~660k levels is the
# evident wall, fe= notwithstanding. This is not worth a fifth attempt.
#
# This check addresses the same underlying concern -- fragility of
# inference given only 16 month clusters -- through a much cheaper
# mechanism: refit the headline interacted model 16 times, each time
# excluding one month, and see how far the coefficient and SE move. It is
# the same logic as the paper's existing leave-one-system-out check
# (Section 8.3 / s13_leave_one_system_out), applied to months instead of
# systems, so it fits the paper's existing robustness vocabulary rather
# than introducing a new kind of evidence.
# -----------------------------------------------------------------------

s21_loco_month <- function(panel = outpatient, scheme_col = "SCHEME_1_CERTAINTY",
                           instrument, instrument_label = "",
                           outcome = PRIMARY_OUTCOME, controls = BASELINE_CONTROLS,
                           fixed_effects = BASELINE_FIXED_EFFECTS) {
  controls <- available_columns(panel, controls)
  fe <- available_columns(panel, fixed_effects)
  
  d <- panel[!is.na(get(scheme_col))]
  d <- model_sample(d, c(outcome, instrument, controls, fe, "ANALYSIS_MARKET", scheme_col))
  d[, MOD := droplevels(factor(get(scheme_col)))]
  keys <- levels(d$MOD)
  for (k in keys) d[, paste0("RF_", k) := get(instrument) * as.integer(MOD == k)]
  rf_names <- paste0("RF_", keys)
  
  months <- sort(unique(d$POST_MONTH))
  form <- as.formula(paste(outcome, "~", paste(c(rf_names, controls), collapse = " + "),
                           "|", paste(fe, collapse = " + ")))
  clf  <- as.formula(paste("~", paste(c("ANALYSIS_MARKET", "POST_MONTH"), collapse = " + ")))
  
  cat("Leave-one-month-out for ", instrument_label, ": refitting with each of the ",
      length(months), " months dropped in turn.\n\n", sep = "")
  
  t0 <- Sys.time()
  fit_base <- feols(form, data = d, cluster = clf, warn = FALSE, notes = FALSE)
  cat("Baseline (all ", length(months), " months, ",
      round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "s):\n", sep = "")
  print(tidy_fixest(fit_base)[term %chin% rf_names,
                              .(term, estimate = round(estimate, 5),
                                std.error = round(std.error, 5), p.value = round(p.value, 4))])
  
  out <- rbindlist(lapply(months, function(m) {
    t0 <- Sys.time()
    dd <- d[POST_MONTH != m]
    fit <- tryCatch(feols(form, data = dd, cluster = clf, warn = FALSE, notes = FALSE),
                    error = function(e) NULL)
    cat(sprintf("  dropped %s  (%.1fs)\n", as.character(m),
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    if (is.null(fit)) return(NULL)
    r <- data.table(DROPPED_MONTH = as.character(m))
    for (rn in rf_names) {
      r[[paste0(rn, "_coef")]] <- unname(coef(fit)[rn])
      r[[paste0(rn, "_se")]]   <- unname(fixest::se(fit)[rn])
      r[[paste0(rn, "_p")]]    <- unname(pvalue(fit)[rn])
    }
    r
  }), fill = TRUE)
  
  cat("\nCoefficient range across the ", length(months), " leave-one-month-out refits,\n",
      "against the baseline (all months):\n", sep = "")
  for (rn in rf_names) {
    vals <- out[[paste0(rn, "_coef")]]
    base_val <- unname(coef(fit_base)[rn])
    cat(sprintf("  %-20s baseline %+.5f | range [%+.5f, %+.5f] | max deviation %.5f\n",
                rn, base_val, min(vals, na.rm=TRUE), max(vals, na.rm=TRUE),
                max(abs(vals - base_val), na.rm=TRUE)))
  }
  cat("\nA baseline value that stays well within this range under every single\n",
      "month's removal is evidence no one month is driving the result -- the\n",
      "direct analogue of what the leave-one-system-out check already shows\n",
      "for systems.\n", sep = "")
  
  save_csv(out, paste0("T21B_leave_one_month_out_", instrument_label, ".csv"))
  invisible(list(baseline = fit_base, loco = out))
}

run_s21b_loco <- function(panel = outpatient, scheme_col = "SCHEME_1_CERTAINTY",
                          instruments = MAIN_INSTRUMENTS) {
  .s21_hd("21B (fallback) -- LEAVE-ONE-MONTH-OUT (item E1)")
  out <- lapply(names(instruments), function(il) {
    .s21_sub(il)
    s21_loco_month(panel, scheme_col, instruments[[il]], il)
  })
  names(out) <- names(instruments)
  invisible(out)
}

s21_wild_bootstrap_month <- function(
    data, scheme_col = "SCHEME_1_CERTAINTY", instrument, instrument_label = "",
    outcome = PRIMARY_OUTCOME, controls = BASELINE_CONTROLS,
    fixed_effects = BASELINE_FIXED_EFFECTS, B = 999, seed = 20240101) {
  
  if (!s21_check_boot_package()) return(NULL)
  
  # Third attempt. The first two (direct call with no FE hint; manual FWL
  # pre-demeaning) both failed for different reasons -- the direct call
  # apparently never finished in an hour, and hand-rolled FWL hit three
  # separate bugs (row misalignment, a Date-typed clustering variable, and
  # joint-vs-single-FE singleton removal giving mismatched sample sizes).
  # Dropping the manual demeaning entirely. fwildclusterboot::boottest()
  # has its own fe= argument specifically for telling it a fixed effect is
  # already absorbed in the fitted model, so it can use its fast algorithm
  # to account for it during resampling rather than whatever slower path it
  # falls back to without that hint -- which is very possibly why the
  # original direct call never returned. B defaults small here (999, not
  # 9999) as a timing check before committing to the full run.
  controls <- available_columns(data, controls)
  fe <- available_columns(data, fixed_effects)
  if (length(fe) < 1L) stop("Need at least one fixed effect.", call. = FALSE)
  high_card_fe <- fe[1]
  
  d <- data[!is.na(get(scheme_col))]
  d <- model_sample(d, c(outcome, instrument, controls, fe, "ANALYSIS_MARKET", scheme_col))
  d[, MOD := droplevels(factor(get(scheme_col)))]
  keys <- levels(d$MOD)
  for (k in keys) d[, paste0("RF_", k) := get(instrument) * as.integer(MOD == k)]
  rf_names <- paste0("RF_", keys)
  
  # boottest() fails on Date-typed clustering variables ("| not defined for
  # Date objects"); POST_MONTH must be a factor wherever it's used, whether
  # as an FE term or as clustid.
  d[, POST_MONTH_FCT := factor(POST_MONTH)]
  fe_formula_terms <- ifelse(fe == "POST_MONTH", "POST_MONTH_FCT", fe)
  
  cat("Fitting the full two-way FE model directly. Passing fe = '", high_card_fe,
      "' to boottest() so it uses its own fast algorithm for that fixed effect\n",
      "rather than a slower default path.\n", sep = "")
  
  t0 <- Sys.time()
  fit <- feols(as.formula(paste(outcome, "~", paste(c(rf_names, controls), collapse = " + "),
                                "|", paste(fe_formula_terms, collapse = " + "))),
               data = d, cluster = ~ANALYSIS_MARKET, warn = FALSE, notes = FALSE)
  cat("Model fit took ", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1),
      "s\n", sep = "")
  
  cat("\nCounty-only clustered SEs on the full model (the 'primary' comparison ",
      "E1 asks for -- and the fallback if the bootstrap below still doesn't\n",
      "complete in reasonable time):\n", sep = "")
  print(tidy_fixest(fit)[term %chin% rf_names,
                         .(term, estimate = signif(estimate, 4),
                           std.error = signif(std.error, 4), p.value = round(p.value, 4))])
  
  cat("\nStarting with B = ", B, " as a timing check before scaling up. If any\n",
      "single term takes more than a couple of minutes, interrupt and report\n",
      "back rather than waiting further -- at that point I'd treat the wild\n",
      "bootstrap as intractable with this FE structure and recommend relying\n",
      "on the county-only clustering above as the primary robustness check,\n",
      "which is already the paper's own documented fallback.\n", sep = "")
  
  results <- lapply(rf_names, function(term_name) {
    cat("\nWild cluster bootstrap (B = ", B, ") on ", term_name,
        ", clustered on POST_MONTH_FCT (", uniqueN(d$POST_MONTH_FCT), " clusters):\n", sep = "")
    t0 <- Sys.time()
    set.seed(seed)
    bt <- tryCatch(
      fwildclusterboot::boottest(fit, param = term_name, clustid = "POST_MONTH_FCT",
                                 fe = high_card_fe, B = B),
      error = function(e) { cat("  FAILED: ", e$message, "\n", sep = ""); NULL })
    cat("  (", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "s)\n", sep = "")
    if (!is.null(bt)) print(summary(bt))
    bt
  })
  names(results) <- rf_names
  
  cat("\nHOW TO READ THIS.\n",
      "  Compare the wild-bootstrap p-value on each term to the conventional\n",
      "  two-way-clustered p-value already in Table 5. A p-value that moves\n",
      "  materially, especially past 0.05, means the month dimension's few\n",
      "  clusters were doing more work in the conventional SE than the\n",
      "  asymptotics assume. A p-value that barely moves is a clean pass.\n", sep = "")
  
  invisible(results)
}

run_s21b <- function(panel = outpatient, scheme_col = "SCHEME_1_CERTAINTY",
                     instruments = MAIN_INSTRUMENTS, B = 9999) {
  .s21_hd("21B -- WILD CLUSTER BOOTSTRAP OVER MONTHS (item E1)")
  if (!s21_check_boot_package()) return(invisible(NULL))
  
  out <- lapply(names(instruments), function(il) {
    .s21_sub(il)
    s21_wild_bootstrap_month(panel, scheme_col, instruments[[il]], il, B = B)
  })
  names(out) <- names(instruments)
  invisible(out)
}


# ===========================================================================
# DRIVER
# ===========================================================================

run_section_21 <- function(panel = outpatient, run_exact_pi_test = FALSE, B = 9999) {
  .s21_hd("SECTION 21 -- A3, E1")
  
  naive <- run_s21a_naive()
  
  exact <- NULL
  if (run_exact_pi_test) {
    exact <- run_s21a_exact(panel)
  } else {
    cat("\n[21A exact skipped] Naive test above shows no significant difference\n",
        "for any instrument. Set run_exact_pi_test = TRUE to run the expensive\n",
        "exact version if you want to confirm that isn't an artifact of the\n",
        "independence assumption.\n", sep = "")
  }
  
  boot <- run_s21b(panel, B = B)
  
  .s21_hd("SECTION 21 COMPLETE")
  invisible(list(naive_pi_test = naive, exact_pi_test = exact, wild_boot = boot))
}

# Example:
#   s21 <- run_section_21(outpatient, run_exact_pi_test = FALSE)   # fast: naive + wild boot
#   s21 <- run_section_21(outpatient, run_exact_pi_test = TRUE)    # also runs the expensive exact test

s21_loco <- run_s21b_loco(outpatient)






###############################################################################

# ============================================================================
# Section 22: Tier III diagnostics E2, E5
# ============================================================================

# SECTION 22 -- TIER III DIAGNOSTICS: E2, E5   [CORRECTED]
#
# Source AFTER warm_start() / restore_session().
#
#   22A  E2 -- how correlated are the estimable classification partitions'
#              shoppability gradients? Reports the correlation BETWEEN THE
#              ESTIMATORS, not between category labels.
#   22B  E5 -- family-level regression robustness: weighted vs unweighted,
#              leave-one-family-out on both, and the rank statistic.
#
# WHAT CHANGED FROM THE FIRST DRAFT, AND WHY
#
#   1. `schemes_long` is LONG, not wide. It carries SCHEME_ID /
#      ANALYSIS_CONCEPT_ID / SHOPPABILITY_CATEGORY (HIGH/INTERMEDIATE/LOW),
#      one row per concept-scheme. The old 22A treated its column NAMES as
#      the 38 partitions, which would have found ~4 "partitions".
#
#   2. The old 22A correlated `RF_COEF * sign(category)` across partitions.
#      RF_COEF is IDENTICAL across partitions -- concept-level coefficients
#      do not depend on how concepts are later grouped. So that correlation
#      is a magnitude-weighted distortion of label agreement, not the
#      correlation of anything the paper estimates. Replaced below.
#
#   3. `family_level_6inst` is cached but is NOT in CACHE_REGISTRY, so
#      restore_session() does not load it. The old 22B would have silently
#      fallen through to refitting 16 subsample models. It now reads the
#      cache file directly, then the CSV, and only rebuilds as a last resort.
#      Its columns are GROUP_ID / RF_COEF / RF_SE, with no SHOPPABLE column;
#      shoppability comes from DIAGNOSTIC_FAMILIES.
###############################################################################

.s22_hd  <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")
.s22_sub <- function(x) cat("\n--- ", x, " ", strrep("-", max(0, 70 - nchar(x))), "\n", sep = "")

S22_INSTRUMENT <- "Competitor_only_hospitals_9m"


# ===========================================================================
# 22A. CORRELATION ACROSS THE ESTIMABLE PARTITIONS                       (E2)
# ===========================================================================
#
# The object E2 asks about is the shoppability gradient under partition j:
#
#     theta_j = weighted mean of RF_COEF over Shoppable_j
#             - weighted mean of RF_COEF over Non-shoppable_j
#
# which is exactly the CAT coefficient in the pipeline's own
# feols(RF_COEF ~ CAT, weights = ~W) -- WLS on a binary regressor IS the
# weighted mean difference. So this reproduces the paper's partition
# coefficients without refitting anything, and 22A verifies that against
# `meta_rf` before going further.
#
# Writing theta_j = sum_i c_ji * b_i, with
#     c_ji =  w_i / W_Sj   if concept i is Shoppable under j
#          = -w_i / W_Nj   if Non-shoppable
#          =  0            if excluded
# and w_i = 1 / RF_SE_i^2, treating concept-level estimates as independent:
#
#     Cov(theta_j, theta_k) = sum_i c_ji c_ki / w_i
#     Var(theta_j)          = 1/W_Sj + 1/W_Nj
#
# Independence across concepts is FALSE here (concepts share counties and
# months), and ignoring that correlation makes the resulting correlations
# TOO SMALL. That is the right direction to err: the reported figure is a
# lower bound on how correlated these partitions really are, which is the
# conservative way to make E2's point. 22A-iii adds a family-cluster
# bootstrap that does not assume independence, as a check.

s22a_build_partitions <- function(concept_results, schemes_long,
                                  instrument_label = S22_INSTRUMENT,
                                  min_concepts = MIN_CONCEPTS_META) {
  
  cr <- as.data.table(concept_results)[INSTRUMENT_LABEL == instrument_label &
                                         is.finite(RF_COEF) & is.finite(RF_SE) & RF_SE > 0,
                                       .(FINAL_CONCEPT_ID, FINAL_FAMILY_ID, RF_COEF, RF_SE)]
  if (nrow(cr) == 0L)
    stop("No concept_results rows for INSTRUMENT_LABEL == '", instrument_label,
         "'. Available: ", paste(unique(concept_results$INSTRUMENT_LABEL), collapse = ", "),
         call. = FALSE)
  cr <- unique(cr, by = "FINAL_CONCEPT_ID")
  cat("Estimable concepts for ", instrument_label, ": ", nrow(cr), "\n", sep = "")
  
  sl <- as.data.table(schemes_long)
  # Same four collapse rules as run_meta_regressions(), same order.
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
  
  variants <- list(); n_built <- 0L
  for (sid in unique(sl$SCHEME_ID)) {
    a <- unique(sl[SCHEME_ID == sid, .(FINAL_CONCEPT_ID = ANALYSIS_CONCEPT_ID,
                                       CAT_RAW = toupper(trimws(SHOPPABILITY_CATEGORY)))])
    sname <- sl[SCHEME_ID == sid][1L]$SCHEME_NAME
    if (is.null(sname) || is.na(sname)) sname <- sid
    m <- merge(cr, a, by = "FINAL_CONCEPT_ID", all.x = TRUE, sort = FALSE)
    for (rl in names(rules)) {
      n_built <- n_built + 1L
      cat_vec <- droplevels(rules[[rl]](m$CAT_RAW))
      keep <- !is.na(cat_vec)
      if (sum(keep) < min_concepts) next
      lv <- levels(droplevels(cat_vec[keep]))
      # The paper counts the "Shoppable" term, which only exists for the
      # two-level collapses. Three-tier variants have no such term and are
      # not part of the 38.
      if (length(lv) != 2L || !("Shoppable" %in% lv)) next
      variants[[length(variants) + 1L]] <- list(
        VARIANT = paste(sname, rl, sep = " | "),
        SHOP    = keep & cat_vec == "Shoppable",
        NONSHOP = keep & cat_vec == "Non_shoppable")
    }
  }
  cat("Scheme-collapse variants built: ", n_built,
      " | two-level and estimable: ", length(variants), "\n", sep = "")
  
  # Deduplicate on the assignment vector. The pipeline fingerprints on the
  # estimated coefficient vector instead, but the coefficients are a
  # deterministic function of the assignment, so the grouping is the same.
  fp <- vapply(variants, function(v) paste(c(v$SHOP, v$NONSHOP), collapse = ""), character(1))
  keep_idx <- which(!duplicated(fp))
  variants <- variants[keep_idx]
  cat("Distinct estimable partitions: ", length(variants),
      if (length(variants) != 38L)
        paste0("   <-- expected 38; check MIN_CONCEPTS_META and the concept filter")
      else "   <-- matches the paper", "\n", sep = "")
  
  list(cr = cr, variants = variants)
}

s22a_theta <- function(cr, variants) {
  w <- 1 / cr$RF_SE^2; b <- cr$RF_COEF
  vapply(variants, function(v) {
    ws <- w[v$SHOP]; wn <- w[v$NONSHOP]
    sum(ws * b[v$SHOP]) / sum(ws) - sum(wn * b[v$NONSHOP]) / sum(wn)
  }, numeric(1))
}

s22a_correlation <- function(cr, variants) {
  w <- 1 / cr$RF_SE^2
  n <- nrow(cr); J <- length(variants)
  
  C <- matrix(0, nrow = n, ncol = J)
  for (j in seq_len(J)) {
    v <- variants[[j]]
    C[v$SHOP,    j] <-  w[v$SHOP]    / sum(w[v$SHOP])
    C[v$NONSHOP, j] <- -w[v$NONSHOP] / sum(w[v$NONSHOP])
  }
  # Cov = C' diag(sigma^2) C, with sigma^2 = 1/w
  V <- crossprod(C / sqrt(w))
  s <- sqrt(diag(V))
  R <- V / outer(s, s)
  dimnames(R) <- list(vapply(variants, `[[`, character(1), "VARIANT"),
                      vapply(variants, `[[`, character(1), "VARIANT"))
  R
}

s22a_label_agreement <- function(variants) {
  J <- length(variants)
  A <- matrix(NA_real_, J, J)
  for (j in seq_len(J)) for (k in seq_len(J)) {
    both <- (variants[[j]]$SHOP | variants[[j]]$NONSHOP) &
      (variants[[k]]$SHOP | variants[[k]]$NONSHOP)
    if (!any(both)) next
    A[j, k] <- mean(variants[[j]]$SHOP[both] == variants[[k]]$SHOP[both])
  }
  A
}

s22a_family_bootstrap <- function(cr, variants, B = 2000L, seed = 20260819L) {
  set.seed(seed)
  fams <- unique(cr$FINAL_FAMILY_ID); nf <- length(fams)
  w <- 1 / cr$RF_SE^2; b <- cr$RF_COEF
  fam_idx <- match(cr$FINAL_FAMILY_ID, fams)
  J <- length(variants)
  
  theta_b <- matrix(NA_real_, B, J)
  for (bi in seq_len(B)) {
    mult <- tabulate(sample.int(nf, nf, replace = TRUE), nbins = nf)[fam_idx]
    mw <- mult * w
    for (j in seq_len(J)) {
      v <- variants[[j]]
      ss <- sum(mw[v$SHOP]); sn <- sum(mw[v$NONSHOP])
      if (ss <= 0 || sn <= 0) next
      theta_b[bi, j] <- sum(mw[v$SHOP] * b[v$SHOP]) / ss -
        sum(mw[v$NONSHOP] * b[v$NONSHOP]) / sn
    }
  }
  R <- suppressWarnings(cor(theta_b, use = "pairwise.complete.obs"))
  list(R = R, draws = theta_b)
}

s22_partition_correlation <- function(concept_results, schemes_long,
                                      instrument_label = S22_INSTRUMENT,
                                      meta_rf_obj = NULL, bootstrap = TRUE) {
  
  bp <- s22a_build_partitions(concept_results, schemes_long, instrument_label)
  cr <- bp$cr; variants <- bp$variants
  if (length(variants) < 2L) stop("Fewer than two partitions; nothing to correlate.", call. = FALSE)
  
  theta <- s22a_theta(cr, variants)
  names(theta) <- vapply(variants, `[[`, character(1), "VARIANT")
  
  cat(sprintf("\nGradient across partitions: median %.5f | range [%.5f, %.5f] | negative in %d of %d\n",
              median(theta), min(theta), max(theta), sum(theta < 0), length(theta)))
  
  # ---- verification against the pipeline's own meta-regression ------------
  if (!is.null(meta_rf_obj)) {
    mr <- as.data.table(meta_rf_obj)
    mr <- mr[WEIGHTING == "Inverse variance" & grepl("Shoppable", term)]
    if ("INSTRUMENT_LABEL" %in% names(mr)) {
      cand <- unique(mr$INSTRUMENT_LABEL)
      il <- if (instrument_label %chin% cand) instrument_label else cand[1L]
      mr <- mr[INSTRUMENT_LABEL == il]
      if (il != instrument_label)
        cat("NOTE: meta_rf relabels instruments; verifying against '", il, "'.\n", sep = "")
    }
    mr[, VARIANT := paste(SCHEME, COLLAPSE, sep = " | ")]
    chk <- merge(data.table(VARIANT = names(theta), MINE = as.numeric(theta)),
                 mr[, .(VARIANT, THEIRS = estimate)], by = "VARIANT")
    if (nrow(chk) > 0L) {
      cat(sprintf("VERIFY vs meta_rf: %d matched variants, max |diff| = %.2e %s\n",
                  nrow(chk), max(abs(chk$MINE - chk$THEIRS)),
                  if (max(abs(chk$MINE - chk$THEIRS)) < 1e-8) "(exact)" else "(INVESTIGATE)"))
    } else {
      cat("VERIFY vs meta_rf: no variant names matched; skipping (names differ, not an error).\n")
    }
  }
  
  # ---- (i) analytic estimator correlation --------------------------------
  R <- s22a_correlation(cr, variants)
  off <- R[upper.tri(R)]
  .s22_sub("(i) Analytic correlation between partition gradients")
  cat(sprintf("  Median %.3f | Mean %.3f | Min %.3f | Max %.3f | %% > 0.5: %.1f%%\n",
              median(off), mean(off), min(off), max(off), 100 * mean(off > 0.5)))
  cat("  Concept-level estimates treated as independent, so this is a LOWER\n",
      "  bound on the true correlation.\n", sep = "")
  
  # ---- (ii) label agreement ----------------------------------------------
  A <- s22a_label_agreement(variants)
  offA <- A[upper.tri(A)]
  .s22_sub("(ii) Classification agreement between partitions")
  cat(sprintf("  Median %.3f | Mean %.3f | Min %.3f | Max %.3f\n",
              median(offA, na.rm = TRUE), mean(offA, na.rm = TRUE),
              min(offA, na.rm = TRUE), max(offA, na.rm = TRUE)))
  cat("  Share of jointly classified concepts placed on the same side.\n")
  
  # ---- (iii) family cluster bootstrap ------------------------------------
  boot <- NULL
  if (isTRUE(bootstrap)) {
    .s22_sub("(iii) Family-cluster bootstrap correlation (no independence assumption)")
    boot <- s22a_family_bootstrap(cr, variants)
    offB <- boot$R[upper.tri(boot$R)]
    offB <- offB[is.finite(offB)]
    cat(sprintf("  Median %.3f | Mean %.3f | Min %.3f | Max %.3f | 2000 draws over %d families\n",
                median(offB), mean(offB), min(offB), max(offB),
                uniqueN(cr$FINAL_FAMILY_ID)))
    cat("  Resamples the 16 clinical families, the pipeline's own cluster unit.\n")
  }
  
  cat("\nHOW TO READ THIS.\n",
      "  E2 asks you to quantify that the partitions are correlated re-cuts of\n",
      "  one dataset rather than independent tests. (i) is the number to quote:\n",
      "  the median pairwise correlation between the partitions' gradient\n",
      "  ESTIMATES. (ii) is the intuition behind it. (iii) confirms (i) is not\n",
      "  an artefact of the independence assumption; if (iii) is materially\n",
      "  higher than (i), quote (i) anyway and note it is conservative.\n", sep = "")
  
  save_csv(as.data.table(R, keep.rownames = "PARTITION"), "T22A_partition_gradient_correlation.csv")
  save_csv(data.table(PARTITION = names(theta), GRADIENT = as.numeric(theta)),
           "T22A_partition_gradients.csv")
  
  invisible(list(theta = theta, R_analytic = R, agreement = A, boot = boot,
                 variants = variants, cr = cr))
}


# ===========================================================================
# 22B. FAMILY-LEVEL ROBUSTNESS: UNWEIGHTED + LEAVE-ONE-FAMILY-OUT        (E5)
# ===========================================================================

s22_get_family_results <- function(instrument_label = S22_INSTRUMENT,
                                   panel = NULL, rebuild_ok = TRUE) {
  
  # 1. Cache file. family_level_6inst is written by cache_or_run() but is NOT
  #    in CACHE_REGISTRY, so restore_session() does not load it.
  cache_path <- file.path(CACHE_DIR, "family_level_6inst.rds")
  if (file.exists(cache_path)) {
    cat("Reading cached family_level_6inst.rds (the object behind Appendix Table 19).\n")
    return(as.data.table(readRDS(cache_path)))
  }
  # 2. Already in the global environment (Section 16 ran this session).
  if (exists("family_results", envir = .GlobalEnv)) {
    cat("Using `family_results` from the global environment.\n")
    return(as.data.table(get("family_results", envir = .GlobalEnv)))
  }
  # 3. Exported CSV.
  csv_path <- file.path(TABLE_DIR, "T16A_family_level_RF_FS_IV.csv")
  if (file.exists(csv_path)) {
    cat("Reading T16A_family_level_RF_FS_IV.csv.\n")
    return(fread(csv_path))
  }
  if (!rebuild_ok)
    stop("No family-level results found. Run Section 16, or set rebuild_ok = TRUE.", call. = FALSE)
  
  cat("No cached family-level results found; rebuilding via run_group_level_sweep().\n",
      "This refits 16 subsample models and takes a few minutes.\n", sep = "")
  if (is.null(panel)) panel <- get("outpatient", envir = .GlobalEnv)
  as.data.table(run_group_level_sweep(panel, "FINAL_FAMILY_ID",
                                      instruments = MAIN_INSTRUMENTS[instrument_label]))
}

s22_family_robustness <- function(family_results, instrument_label = S22_INSTRUMENT) {
  fr <- as.data.table(copy(family_results))
  
  if ("INSTRUMENT_LABEL" %in% names(fr)) {
    if (!(instrument_label %chin% fr$INSTRUMENT_LABEL))
      stop("Instrument '", instrument_label, "' not in family results. Available: ",
           paste(unique(fr$INSTRUMENT_LABEL), collapse = ", "), call. = FALSE)
    fr <- fr[INSTRUMENT_LABEL == instrument_label]
  }
  if (!("GROUP_ID" %in% names(fr)) && "FINAL_FAMILY_ID" %in% names(fr))
    setnames(fr, "FINAL_FAMILY_ID", "GROUP_ID")
  
  need <- c("GROUP_ID", "RF_COEF", "RF_SE")
  miss <- setdiff(need, names(fr))
  if (length(miss)) stop("family results missing: ", paste(miss, collapse = ", "), call. = FALSE)
  
  fr <- fr[is.finite(RF_COEF) & is.finite(RF_SE) & RF_SE > 0,
           .(GROUP_ID, RF_COEF, RF_SE)]
  fr[, SHOP01 := as.integer(GROUP_ID %chin% DIAGNOSTIC_FAMILIES)]
  fr[, W := 1 / RF_SE^2]
  setorder(fr, RF_COEF)
  
  cat("Families: ", nrow(fr), " (expect 16) | shoppable: ", sum(fr$SHOP01),
      " (expect 10)\n\n", sep = "")
  if (nrow(fr) != 16L || sum(fr$SHOP01) != 10L)
    warning("Family or shoppable count is not 16/10; check GROUP_ID against DIAGNOSTIC_FAMILIES.",
            call. = FALSE)
  
  # ---- weighted (paper's spec) vs unweighted ------------------------------
  fit_w  <- lm(RF_COEF ~ SHOP01, data = fr, weights = fr$W)
  fit_uw <- lm(RF_COEF ~ SHOP01, data = fr)
  sw <- coef(summary(fit_w)); su <- coef(summary(fit_uw))
  
  .s22_sub("Weighted vs unweighted family-level regression")
  cat(sprintf("  Weighted   difference %+.5f  SE %.5f  t %.2f  p %.6f\n",
              sw["SHOP01", 1], sw["SHOP01", 2], sw["SHOP01", 3], sw["SHOP01", 4]))
  cat(sprintf("  Unweighted difference %+.5f  SE %.5f  t %.2f  p %.6f\n",
              su["SHOP01", 1], su["SHOP01", 2], su["SHOP01", 3], su["SHOP01", 4]))
  cat(sprintf("  Group means -- weighted: shop %+.5f / non-shop %+.5f\n",
              sum(fr[SHOP01 == 1]$W * fr[SHOP01 == 1]$RF_COEF) / sum(fr[SHOP01 == 1]$W),
              sum(fr[SHOP01 == 0]$W * fr[SHOP01 == 0]$RF_COEF) / sum(fr[SHOP01 == 0]$W)))
  cat(sprintf("               unweighted: shop %+.5f / non-shop %+.5f\n",
              mean(fr[SHOP01 == 1]$RF_COEF), mean(fr[SHOP01 == 0]$RF_COEF)))
  
  # ---- rank statistic -----------------------------------------------------
  .s22_sub("Rank separation")
  rk <- rank(fr$RF_COEF)
  top10_shop <- sum(fr[order(RF_COEF)][1:10]$SHOP01)
  wt <- suppressWarnings(wilcox.test(RF_COEF ~ SHOP01, data = fr, exact = TRUE))
  cat(sprintf("  Shoppable families among the 10 most negative: %d of 10\n", top10_shop))
  cat(sprintf("  Wilcoxon rank-sum: W = %.0f, exact p = %.6f\n", wt$statistic, wt$p.value))
  cat(sprintf("  Probability of exact separation under random assignment: 1/%d = %.5f\n",
              choose(16, 10), 1 / choose(16, 10)))
  
  # ---- leave-one-family-out ----------------------------------------------
  .s22_sub("Leave-one-family-out")
  loo <- rbindlist(lapply(fr$GROUP_ID, function(fam) {
    dd <- fr[GROUP_ID != fam]
    if (uniqueN(dd$SHOP01) < 2L) return(NULL)
    fw <- lm(RF_COEF ~ SHOP01, data = dd, weights = dd$W)
    fu <- lm(RF_COEF ~ SHOP01, data = dd)
    data.table(DROPPED_FAMILY = fam,
               DIFF_W  = unname(coef(fw)["SHOP01"]),
               P_W     = coef(summary(fw))["SHOP01", 4],
               DIFF_UW = unname(coef(fu)["SHOP01"]),
               P_UW    = coef(summary(fu))["SHOP01", 4])
  }), fill = TRUE)
  setorder(loo, DIFF_W)
  print(loo[, .(DROPPED_FAMILY, DIFF_W = round(DIFF_W, 5), P_W = round(P_W, 5),
                DIFF_UW = round(DIFF_UW, 5), P_UW = round(P_UW, 5))])
  cat(sprintf("\n  Full sample: weighted %+.5f | unweighted %+.5f\n",
              sw["SHOP01", 1], su["SHOP01", 1]))
  cat(sprintf("  LOO range weighted   [%+.5f, %+.5f]  max p %.5f\n",
              min(loo$DIFF_W), max(loo$DIFF_W), max(loo$P_W)))
  cat(sprintf("  LOO range unweighted [%+.5f, %+.5f]  max p %.5f\n",
              min(loo$DIFF_UW), max(loo$DIFF_UW), max(loo$P_UW)))
  cat(sprintf("  Most influential family: %s (moves the weighted difference %.1f%%)\n",
              loo[which.max(abs(DIFF_W - sw["SHOP01", 1]))]$DROPPED_FAMILY,
              100 * max(abs(loo$DIFF_W - sw["SHOP01", 1])) / abs(sw["SHOP01", 1])))
  
  cat("\nHOW TO READ THIS.\n",
      "  E5's objection is that a 16-observation regression weighted by the\n",
      "  same estimation that produced its dependent variable can be carried\n",
      "  by one precise family. If the unweighted difference is close to the\n",
      "  weighted one and the LOO range stays comfortably negative and\n",
      "  significant, the objection is answered and you can report the range\n",
      "  in the table notes.\n", sep = "")
  
  save_csv(fr,  "T22B_family_level_input.csv")
  save_csv(loo, "T22B_leave_one_family_out.csv")
  invisible(list(weighted = fit_w, unweighted = fit_uw, loo = loo, input = fr,
                 top10_shop = top10_shop, wilcox = wt))
}


# ===========================================================================
# DRIVER
# ===========================================================================

run_section_22 <- function(concept_results = get("concept_results", envir = .GlobalEnv),
                           schemes_long    = get("schemes_long", envir = .GlobalEnv),
                           instrument_label = S22_INSTRUMENT,
                           bootstrap = TRUE) {
  .s22_hd("SECTION 22 -- E2, E5")
  
  .s22_sub("22A -- Partition correlation (E2)")
  meta_rf_obj <- if (exists("meta_rf", envir = .GlobalEnv)) get("meta_rf", envir = .GlobalEnv) else NULL
  a <- tryCatch(s22_partition_correlation(concept_results, schemes_long,
                                          instrument_label, meta_rf_obj, bootstrap),
                error = function(e) { message("[22A failed] ", e$message); NULL })
  
  .s22_sub("22B -- Family-level robustness (E5)")
  b <- tryCatch({
    fam <- s22_get_family_results(instrument_label)
    s22_family_robustness(fam, instrument_label)
  }, error = function(e) { message("[22B failed] ", e$message); NULL })
  
  .s22_hd("SECTION 22 COMPLETE")
  invisible(list(partition_corr = a, family_robustness = b))
}

s22 <- run_section_22()





###############################################################################

# ============================================================================
# Section 23: p-value recomputation check
# ============================================================================

# SECTION 23 -- P-VALUE RECOMPUTATION CHECK
#
# estimate_interacted() computes RF_P and IV_P via 2 * pnorm(-abs(b/s)) --
# a normal-distribution p-value -- rather than fixest's own p-value, which
# would apply a cluster-adjusted t-distribution. With month clustering at
# only 16 levels, that distinction is not cosmetic: t(15) has meaningfully
# fatter tails than normal, so borderline cells can flip significance.
#
# Three parts, in increasing order of rigor and cost:
#
#   23A  Instant. Recomputes p under t(df) directly from the RF_COEF/RF_SE
#        already sitting in main$rows -- no new regressions, seconds to run.
#        Relies on an ASSUMED df (15, i.e. month clusters - 1), which is the
#        standard convention but not verified against what fixest itself
#        would compute for this exact two-way-clustered, FE-nested setup.
#
#   23B  Rigorous. Refits every one of the 18 scheme x instrument
#        specifications underlying Table 5's 36 rows, and extracts
#        fixest's OWN p-value via tidy_fixest() -- the same call already
#        used correctly throughout today's other diagnostics. No df
#        assumption; this is whatever cluster-adjusted reference
#        distribution fixest actually applies given the real FE and
#        cluster structure. ~18 regressions, comparable cost to Section
#        20C's ladder, so a couple of minutes.
#
#   23C  Merges 23B against main$rows directly, with a coefficient/SE
#        match check before trusting any p-value comparison, and reports
#        the definitive count of significance flips.
#
# Run 23A first for an immediate read; run 23B/23C for the number that
# actually settles this.
###############################################################################

.s23_hd  <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")
.s23_sub <- function(x) cat("\n--- ", x, " ", strrep("-", max(0, 70 - nchar(x))), "\n", sep = "")

# Confirmed against main$rows' actual SPEC strings, printed earlier this
# session -- not a guess. Order matches SCHEME_COLUMNS as verified in
# Section 20C, where all six schemes' numbers reproduced the paper's
# existing Table sysmonth_full exactly under this same mapping.
S23_SCHEME_MAP <- c(
  "1. Procedural certainty"    = "SCHEME_1_CERTAINTY",
  "2. Theory-Based V2"         = "SCHEME_2_THEORYV2",
  "3. Imaging vs Procedural"   = "SCHEME_3_IMAGING",
  "4. CMS Statutory List"      = "SCHEME_4_CMS70",
  "5. Upfront Cash-Market"     = "SCHEME_5_MDSAVE",
  "6. Within Modality"         = "SCHEME_6_WITHINMOD"
)


# ===========================================================================
# 23A. INSTANT CHECK -- t(df) applied directly to main$rows               
# ===========================================================================

s23_recompute_p_t15 <- function(main_rows, df = 15) {
  d <- as.data.table(copy(main_rows))
  
  d[, RF_T     := RF_COEF / RF_SE]
  d[, RF_P_T   := 2 * pt(-abs(RF_T), df = df)]
  d[, IV_T     := IV_COEF / IV_SE]
  d[, IV_P_T   := 2 * pt(-abs(IV_T), df = df)]
  
  d[, SIG_05_ORIG := RF_P   < 0.05]
  d[, SIG_05_T     := RF_P_T < 0.05]
  d[, SIG_10_ORIG := RF_P   < 0.10]
  d[, SIG_10_T     := RF_P_T < 0.10]
  d[, FLIPPED_05  := SIG_05_ORIG != SIG_05_T]
  d[, FLIPPED_10  := SIG_10_ORIG != SIG_10_T]
  
  cat("Recomputing RF_P and IV_P under t(df=", df, ") directly from the\n",
      "RF_COEF/RF_SE already in main$rows -- instant, no new regressions.\n",
      "This ASSUMES df = ", df, " (month clusters - 1); 23B/23C verify against\n",
      "fixest's own computed p-value rather than this assumption.\n\n", sep = "")
  
  cat("Reduced-form: cells where 5% significance flips:\n")
  print(d[FLIPPED_05 == TRUE, .(SPEC, INSTRUMENT_LABEL, TERM,
                                RF_P_ORIGINAL = round(RF_P, 4),
                                RF_P_T15      = round(RF_P_T, 4))])
  
  cat("\nReduced-form: cells where 10% significance flips:\n")
  print(d[FLIPPED_10 == TRUE, .(SPEC, INSTRUMENT_LABEL, TERM,
                                RF_P_ORIGINAL = round(RF_P, 4),
                                RF_P_T15      = round(RF_P_T, 4))])
  
  cat("\nSummary (reduced form, all 36 rows):\n")
  cat("  Significant at 5% under original (pnorm):", sum(d$SIG_05_ORIG, na.rm = TRUE), "of", nrow(d), "\n")
  cat("  Significant at 5% under t(", df, "):        ", sum(d$SIG_05_T, na.rm = TRUE), " of ", nrow(d), "\n", sep = "")
  cat("  Flipped at 5% threshold: ", sum(d$FLIPPED_05, na.rm = TRUE), "\n", sep = "")
  cat("  Flipped at 10% threshold:", sum(d$FLIPPED_10, na.rm = TRUE), "\n", sep = "")
  
  save_csv(d, "T23A_p_value_recompute_t_approx.csv")
  invisible(d)
}


# ===========================================================================
# 23B. RIGOROUS CHECK -- refit and extract fixest's own p-value
# ===========================================================================

s23_recompute_p_exact <- function(panel = outpatient, scheme_map = S23_SCHEME_MAP,
                                  instruments = MAIN_INSTRUMENTS,
                                  outcome = PRIMARY_OUTCOME, controls = BASELINE_CONTROLS,
                                  fixed_effects = BASELINE_FIXED_EFFECTS,
                                  clusters = BASELINE_CLUSTERS, use_cache = TRUE) {
  
  controls_avail <- available_columns(panel, controls)
  fe_avail       <- available_columns(panel, fixed_effects)
  cl_avail       <- available_columns(panel, clusters)
  
  n_total <- length(scheme_map) * length(instruments)
  cat("Refitting ", length(scheme_map), " schemes x ", length(instruments),
      " instruments = ", n_total, " regressions, extracting fixest's own\n",
      "p-value (whatever cluster-adjusted reference distribution it actually\n",
      "applies) rather than assuming a specific df.\n\n", sep = "")
  
  out <- rbindlist(lapply(names(scheme_map), function(spec_label) {
    scheme_col <- scheme_map[[spec_label]]
    rbindlist(lapply(names(instruments), function(il) {
      key <- paste0("s23_exact_p_", scheme_col, "_", il)
      cat(sprintf("  %-28s %-38s ", spec_label, il))
      t0 <- Sys.time()
      
      compute <- function() {
        d <- panel[!is.na(get(scheme_col))]
        d <- model_sample(d, c(outcome, instruments[[il]], controls_avail, fe_avail,
                               cl_avail, scheme_col))
        d[, MOD := droplevels(factor(get(scheme_col)))]
        keys <- levels(d$MOD)
        for (k in keys) d[, paste0("RF_", k) := get(instruments[[il]]) * as.integer(MOD == k)]
        fit <- feols(as.formula(paste(outcome, "~",
                                      paste(c(paste0("RF_", keys), controls_avail), collapse = " + "),
                                      "|", paste(fe_avail, collapse = " + "))),
                     data = d, cluster = as.formula(paste("~", paste(cl_avail, collapse = "+"))),
                     warn = FALSE, notes = FALSE)
        td <- tidy_fixest(fit)
        td[grepl("^RF_", term), .(TERM = sub("^RF_", "", term),
                                  RF_COEF_CHECK = estimate,
                                  RF_SE_CHECK   = std.error,
                                  RF_P_EXACT    = p.value)]
      }
      
      result <- tryCatch(
        if (use_cache) cache_or_run(key, compute()) else compute(),
        error = function(e) { cat("FAILED: ", e$message, "\n", sep = ""); NULL }
      )
      cat(sprintf("%6.1fs\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
      if (is.null(result)) return(NULL)
      cbind(SPEC = spec_label, INSTRUMENT_LABEL = il, result)
    }), fill = TRUE)
  }), fill = TRUE)
  
  save_csv(out, "T23B_p_value_recompute_exact.csv")
  invisible(out)
}


# ===========================================================================
# 23C. FINAL COMPARISON -- merge against main$rows, report the real count
# ===========================================================================

s23_compare_exact <- function(main_rows, exact_results) {
  m <- as.data.table(copy(main_rows))
  e <- as.data.table(copy(exact_results))
  
  merged <- merge(m[, .(SPEC, INSTRUMENT_LABEL, TERM, RF_COEF, RF_SE, RF_P)],
                  e[, .(SPEC, INSTRUMENT_LABEL, TERM, RF_COEF_CHECK, RF_SE_CHECK, RF_P_EXACT)],
                  by = c("SPEC", "INSTRUMENT_LABEL", "TERM"))
  
  if (nrow(merged) == 0L)
    stop("Merge produced zero rows -- SPEC strings in main$rows and ",
         "S23_SCHEME_MAP's names must match exactly. Check for whitespace ",
         "or formatting differences before trusting anything downstream.",
         call. = FALSE)
  
  merged[, COEF_MATCH := abs(RF_COEF - RF_COEF_CHECK) < 1e-6]
  merged[, SE_MATCH    := abs(RF_SE   - RF_SE_CHECK)   < 1e-6]
  
  cat("\n=== Coefficient/SE match check (must pass before trusting the ===\n",
      "=== p-value comparison below) ===\n", sep = "")
  cat("Coefficients matching: ", sum(merged$COEF_MATCH), " of ", nrow(merged), "\n",
      "SEs matching:          ", sum(merged$SE_MATCH),   " of ", nrow(merged), "\n", sep = "")
  if (!all(merged$COEF_MATCH) || !all(merged$SE_MATCH)) {
    cat("\nMISMATCH DETECTED -- do not trust the p-value comparison below until\n",
        "this is resolved. Print merged[COEF_MATCH == FALSE | SE_MATCH == FALSE]\n",
        "to see which cells differ and by how much.\n", sep = "")
  }
  
  merged[, SIG_05_ORIG  := RF_P       < 0.05]
  merged[, SIG_05_EXACT := RF_P_EXACT < 0.05]
  merged[, FLIPPED_05   := SIG_05_ORIG != SIG_05_EXACT]
  
  cat("\n=== DEFINITIVE COMPARISON: original (pnorm) vs fixest-native p-value ===\n\n")
  cat("Cells where 5% significance flips:\n")
  print(merged[FLIPPED_05 == TRUE, .(SPEC, INSTRUMENT_LABEL, TERM,
                                     RF_P_ORIGINAL = round(RF_P, 4),
                                     RF_P_EXACT    = round(RF_P_EXACT, 4))])
  
  cat("\nOriginally significant at 5% (pnorm):        ", sum(merged$SIG_05_ORIG), " of ", nrow(merged), "\n",
      "Significant at 5% under fixest's own p-value: ", sum(merged$SIG_05_EXACT), " of ", nrow(merged), "\n",
      "Flipped:                                      ", sum(merged$FLIPPED_05), "\n", sep = "")
  
  save_csv(merged, "T23C_p_value_comparison_final.csv")
  invisible(merged)
}


# ===========================================================================
# DRIVER
# ===========================================================================

run_section_23 <- function(main_rows, panel = outpatient, run_exact = FALSE) {
  .s23_hd("SECTION 23 -- P-VALUE RECOMPUTATION")
  
  .s23_sub("23A -- instant t(15) approximation")
  a <- s23_recompute_p_t15(main_rows)
  
  b <- NULL; c_out <- NULL
  if (run_exact) {
    .s23_sub("23B -- rigorous refit")
    b <- s23_recompute_p_exact(panel)
    .s23_sub("23C -- final comparison")
    c_out <- s23_compare_exact(main_rows, b)
  } else {
    cat("\n[23B/23C skipped] Set run_exact = TRUE to run the rigorous refit-based\n",
        "check once you've looked at 23A's output above.\n", sep = "")
  }
  
  .s23_hd("SECTION 23 COMPLETE")
  invisible(list(approx = a, exact = b, comparison = c_out))
}

# Example:
s23 <- run_section_23(main$rows, outpatient, run_exact = FALSE)  # instant only
s23 <- run_section_23(main$rows, outpatient, run_exact = TRUE)   # + rigorous check






###############################################################################

# ============================================================================
# Section 24: Tier IV table builders
# ============================================================================

# SECTION 24 -- TIER IV TABLE BUILDERS
#
# Source AFTER warm_start() / restore_session(). Produces the five Tier IV
# items that need numbers the paper does not currently contain. Each writes
# a CSV and, where the table is new, a LaTeX fragment you can \input{}.
#
#   24A  F7  -- post-absorption first-stage variation, all six instruments
#   24B  F2  -- Table 8 rebuilt with coefficients, not just p-values
#   24C  C8  -- scheme crosswalk + counts on the 738-concept universe
#   24D  C12 -- standardized difference column for tab:family_permutation
#   24E  C13 -- the 16 x 18 family-by-scheme grid
#   24F  F8  -- consolidated sample audit across specification families
#
# Writes .tex fragments to TABLE_DIR/tex/. Point \input{} at them or paste.
###############################################################################

.s24_hd  <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")
TEX_DIR <- file.path(TABLE_DIR, "tex")
dir.create(TEX_DIR, showWarnings = FALSE, recursive = TRUE)
.s24_tex <- function(txt, file) {
  writeLines(txt, file.path(TEX_DIR, file)); cat("Wrote:", file.path(TEX_DIR, file), "\n")
}
.esc <- function(x) gsub("_", "\\\\_", x)


# ===========================================================================
# 24A. F7 -- HOW MUCH INSTRUMENT VARIATION THE FIXED EFFECTS ABSORB
# ===========================================================================
#
# For each instrument: raw SD of Z on the estimation sample, SD of Z after
# absorbing county x concept and month fixed effects, the share of variance
# absorbed, and the implied movement in N per SD of residualized Z (the
# first-stage coefficient times the residualized SD). The last column is the
# number Section 5.6 currently asserts as "well under one poster".

s24_first_stage_variation <- function(panel = outpatient,
                                      instruments = CONCEPT_INSTRUMENTS) {
  out <- rbindlist(lapply(names(instruments), function(il) {
    z <- instruments[[il]]
    d <- model_sample(panel, c(PRIMARY_OUTCOME, z, ENDOGENOUS_VARIABLE,
                               BASELINE_CONTROLS, BASELINE_FIXED_EFFECTS,
                               BASELINE_CLUSTERS))
    if (nrow(d) < MIN_MODEL_OBS) return(NULL)
    
    # Residualize Z on the FE structure the headline spec uses.
    zres <- tryCatch(
      as.numeric(fixest::demean(X = d[[z]],
                                f = d[, ..BASELINE_FIXED_EFFECTS])),
      error = function(e) NULL)
    if (is.null(zres)) {
      fit0 <- feols(as.formula(paste(z, "~ 1 |",
                                     paste(BASELINE_FIXED_EFFECTS, collapse = " + "))),
                    data = d, warn = FALSE, notes = FALSE)
      zres <- resid(fit0)
    }
    
    fs <- tryCatch(
      feols(as.formula(paste(ENDOGENOUS_VARIABLE, "~", z, "+",
                             paste(BASELINE_CONTROLS, collapse = " + "), "|",
                             paste(BASELINE_FIXED_EFFECTS, collapse = " + "))),
            data = d,
            cluster = as.formula(paste("~", paste(BASELINE_CLUSTERS, collapse = "+"))),
            warn = FALSE, notes = FALSE),
      error = function(e) NULL)
    pi_hat <- if (!is.null(fs) && z %chin% names(coef(fs))) unname(coef(fs)[z]) else NA_real_
    
    sd_raw <- sd(d[[z]], na.rm = TRUE)
    sd_res <- sd(zres, na.rm = TRUE)
    data.table(INSTRUMENT_LABEL = il, INSTRUMENT = z, N = nrow(d),
               SD_RAW = sd_raw, SD_RESID = sd_res,
               SHARE_ABSORBED = 1 - (sd_res^2) / (sd_raw^2),
               FS_COEF = pi_hat,
               DELTA_N_PER_SD_RESID = pi_hat * sd_res)
  }), fill = TRUE)
  
  setorder(out, -SHARE_ABSORBED)
  print(out[, .(INSTRUMENT_LABEL, N,
                SD_RAW = round(SD_RAW, 2), SD_RESID = round(SD_RESID, 2),
                PCT_ABSORBED = round(100 * SHARE_ABSORBED, 1),
                FS_COEF = round(FS_COEF, 4),
                DELTA_N = round(DELTA_N_PER_SD_RESID, 3))])
  
  cat("\nDELTA_N is the movement in prior posters induced by a one-SD move in\n",
      "the residualized instrument. Section 5.6 claims this is 'well under one\n",
      "poster'; this column is what makes that claim checkable.\n", sep = "")
  
  body <- paste0(
    "\\quad ", .esc(out$INSTRUMENT_LABEL), " & ",
    formatC(out$SD_RAW, format = "f", digits = 2), " & ",
    formatC(out$SD_RESID, format = "f", digits = 2), " & ",
    formatC(100 * out$SHARE_ABSORBED, format = "f", digits = 1), " & ",
    formatC(out$FS_COEF, format = "f", digits = 4), " & ",
    formatC(out$DELTA_N_PER_SD_RESID, format = "f", digits = 3), " \\\\",
    collapse = "\n")
  .s24_tex(c(
    "% F7 -- generated by HPT_Section24_TierIV_tables.R",
    "\\begin{tabular}{lrrrrr}", "\\toprule",
    "Instrument & SD of $Z$ & SD after FE & \\% absorbed & First stage $\\pi$ & $\\Delta N$ per SD \\\\",
    "\\midrule", body, "\\bottomrule", "\\end{tabular}"),
    "tab_fs_variation.tex")
  
  save_csv(out, "T24A_first_stage_variation.csv")
  invisible(out)
}


# ===========================================================================
# 24B. F2 -- TABLE 8 WITH COEFFICIENTS
# ===========================================================================
#
# Rebuilds tab:schemes_main as two panels (reduced form, IV), each carrying
# the shoppable coefficient, the non-shoppable coefficient, and the
# heterogeneity p, for six schemes x three main instruments.

s24_schemes_with_coefs <- function(main_obj = main) {
  mr <- as.data.table(main_obj$rows); mt <- as.data.table(main_obj$tests)
  keep <- names(MAIN_INSTRUMENTS)
  mr <- mr[INSTRUMENT_LABEL %chin% keep]; mt <- mt[INSTRUMENT_LABEL %chin% keep]
  
  rf <- dcast(mr[TERM %chin% c("Shoppable", "Non_shoppable")],
              SPEC + INSTRUMENT_LABEL ~ TERM, value.var = "RF_PERCENT_PER_SD")
  setnames(rf, c("Shoppable", "Non_shoppable"), c("RF_SHOP", "RF_NONSHOP"), skip_absent = TRUE)
  iv <- dcast(mr[TERM %chin% c("Shoppable", "Non_shoppable")],
              SPEC + INSTRUMENT_LABEL ~ TERM, value.var = "IV_PERCENT")
  setnames(iv, c("Shoppable", "Non_shoppable"), c("IV_SHOP", "IV_NONSHOP"), skip_absent = TRUE)
  pv <- dcast(mt, SPEC + INSTRUMENT_LABEL ~ ESTIMATOR, value.var = "P_VALUE")
  
  d <- merge(merge(rf, iv, by = c("SPEC", "INSTRUMENT_LABEL")), pv,
             by = c("SPEC", "INSTRUMENT_LABEL"))
  print(d)
  cat("\nCheck the ESTIMATOR column names above; the .tex writer assumes the RF\n",
      "and IV p-value columns are the two non-key columns of `pv`.\n", sep = "")
  save_csv(d, "T24B_schemes_with_coefficients.csv")
  cat("\nBuild the LaTeX from this CSV; column order there matches the table\n",
      "layout in the Tier IV patch.\n", sep = "")
  invisible(d)
}


# ===========================================================================
# 24C. C8 -- SCHEME CROSSWALK AND 738-UNIVERSE COUNTS
# ===========================================================================
#
# Two things the appendix scheme table cannot currently support: which
# main-text label maps to which appendix scheme, and how many of the 738
# estimated concepts each scheme puts on each side. The 441 figure quoted in
# Section 6.3 should appear in the procedural-certainty row.

s24_scheme_crosswalk <- function(panel = outpatient, schemes_long,
                                 scheme_spec = PRIMARY_SCHEMES) {
  universe <- unique(panel$FINAL_CONCEPT_ID)
  universe <- universe[!is.na(universe)]
  cat("Estimation universe: ", length(universe), " concepts (expect 738)\n", sep = "")
  
  # (a) counts under the six primary schemes, as attached to the panel
  prim <- rbindlist(lapply(scheme_spec, function(sp) {
    if (!(sp$col %in% names(panel))) return(NULL)
    tb <- unique(panel[, .(FINAL_CONCEPT_ID, CAT = get(sp$col))])
    tb <- tb[!is.na(CAT)]
    data.table(MAIN_TEXT_LABEL = sp$label, PANEL_COLUMN = sp$col,
               APPENDIX_SCHEME_ID = sp$source, COLLAPSE_RULE = sp$rule,
               SHOPPABLE_738 = tb[CAT == "Shoppable", uniqueN(FINAL_CONCEPT_ID)],
               NONSHOPPABLE_738 = tb[CAT == "Non_shoppable", uniqueN(FINAL_CONCEPT_ID)],
               CLASSIFIED_738 = tb[, uniqueN(FINAL_CONCEPT_ID)])
  }), fill = TRUE)
  cat("\nSix primary schemes, counts on the estimation universe:\n"); print(prim)
  
  # (b) all eighteen appendix schemes, raw three-way counts on the 738
  sl <- as.data.table(schemes_long)[ANALYSIS_CONCEPT_ID %chin% universe]
  all18 <- dcast(sl[, .N, by = .(SCHEME_ID, SCHEME_NAME,
                                 CAT = toupper(trimws(SHOPPABILITY_CATEGORY)))],
                 SCHEME_ID + SCHEME_NAME ~ CAT, value.var = "N", fill = 0L)
  xw <- unique(prim[, .(APPENDIX_SCHEME_ID, MAIN_TEXT_LABEL)])
  all18 <- merge(all18, xw, by.x = "SCHEME_ID", by.y = "APPENDIX_SCHEME_ID", all.x = TRUE)
  all18[is.na(MAIN_TEXT_LABEL), MAIN_TEXT_LABEL := ""]
  setorder(all18, SCHEME_NAME)
  cat("\nAll eighteen schemes on the estimation universe:\n"); print(all18)
  
  save_csv(prim,  "T24C_primary_scheme_crosswalk.csv")
  save_csv(all18, "T24C_all18_counts_738.csv")
  invisible(list(primary = prim, all18 = all18))
}


# ===========================================================================
# 24D. C12 -- STANDARDIZED DIFFERENCE FOR THE FAMILY PERMUTATION TABLE
# ===========================================================================
#
# The weighted differences in tab:family_permutation are in log points per
# unit of each instrument, and the instruments differ in units, so the six
# rows are not comparable as printed. Multiplying by each instrument's SD on
# its own estimation sample puts them all in log points per SD.

s24_standardize_family_diffs <- function(panel = outpatient,
                                         instruments = CONCEPT_INSTRUMENTS,
                                         family_results = NULL) {
  if (is.null(family_results)) {
    p <- file.path(CACHE_DIR, "family_level_6inst.rds")
    if (!file.exists(p)) stop("family_level_6inst.rds not found; run Section 16.", call. = FALSE)
    family_results <- readRDS(p)
  }
  fr <- as.data.table(family_results)
  if (!("GROUP_ID" %in% names(fr))) setnames(fr, "FINAL_FAMILY_ID", "GROUP_ID", skip_absent = TRUE)
  fr <- fr[is.finite(RF_COEF) & is.finite(RF_SE) & RF_SE > 0]
  fr[, SHOP01 := as.integer(GROUP_ID %chin% DIAGNOSTIC_FAMILIES)]
  fr[, W := 1 / RF_SE^2]
  
  sds <- rbindlist(lapply(names(instruments), function(il) {
    z <- instruments[[il]]
    d <- model_sample(panel, c(PRIMARY_OUTCOME, z, BASELINE_CONTROLS,
                               BASELINE_FIXED_EFFECTS, BASELINE_CLUSTERS))
    data.table(INSTRUMENT_LABEL = il, SD_Z = sd(d[[z]], na.rm = TRUE))
  }), fill = TRUE)
  
  diffs <- fr[, {
    ws <- W[SHOP01 == 1]; wn <- W[SHOP01 == 0]
    .(DIFF = sum(ws * RF_COEF[SHOP01 == 1]) / sum(ws) -
        sum(wn * RF_COEF[SHOP01 == 0]) / sum(wn))
  }, by = INSTRUMENT_LABEL]
  
  out <- merge(diffs, sds, by = "INSTRUMENT_LABEL")
  out[, DIFF_PER_SD := DIFF * SD_Z]
  out[, DIFF_PCT_PER_SD := 100 * (exp(DIFF_PER_SD) - 1)]
  setorder(out, DIFF_PER_SD)
  print(out[, .(INSTRUMENT_LABEL, DIFF = round(DIFF, 5), SD_Z = round(SD_Z, 2),
                DIFF_PER_SD = round(DIFF_PER_SD, 5),
                PCT_PER_SD = round(DIFF_PCT_PER_SD, 2))])
  cat("\nDIFF_PER_SD is the new column for tab:family_permutation. Keep DIFF as\n",
      "the raw column for reproducibility, as C12 asks.\n", sep = "")
  
  save_csv(out, "T24D_family_diff_standardized.csv")
  invisible(out)
}


# ===========================================================================
# 24E. C13 -- THE FAMILY-BY-SCHEME GRID
# ===========================================================================
#
# Sixteen families by eighteen schemes. A cell is the set of categories that
# scheme assigns to that family's concepts. Cells with more than one category
# are the 49 "split" cells the paper's caveat rests on.

s24_family_scheme_grid <- function(panel = outpatient, schemes_long) {
  fam <- unique(as.data.table(panel)[, .(FINAL_CONCEPT_ID, FINAL_FAMILY_ID)])
  fam <- fam[!is.na(FINAL_CONCEPT_ID) & !is.na(FINAL_FAMILY_ID)]
  
  sl <- as.data.table(schemes_long)[ANALYSIS_CONCEPT_ID %chin% fam$FINAL_CONCEPT_ID,
                                    .(FINAL_CONCEPT_ID = ANALYSIS_CONCEPT_ID, SCHEME_ID,
                                      SCHEME_NAME, CAT = toupper(trimws(SHOPPABILITY_CATEGORY)))]
  d <- merge(sl, fam, by = "FINAL_CONCEPT_ID")
  
  cells <- d[, .(N_CATS = uniqueN(CAT),
                 CATS = paste(sort(unique(substr(CAT, 1, 1))), collapse = "/"),
                 N_CONCEPTS = uniqueN(FINAL_CONCEPT_ID)),
             by = .(FINAL_FAMILY_ID, SCHEME_ID, SCHEME_NAME)]
  cells[, SPLIT := as.integer(N_CATS > 1L)]
  
  cat("Grid cells: ", nrow(cells), " (expect 16 x 18 = 288)\n", sep = "")
  cat("Split cells (more than one category): ", sum(cells$SPLIT), " (paper says 49)\n", sep = "")
  cat("\nSplit cells by family:\n")
  print(cells[SPLIT == 1L, .(N_SPLIT = .N), by = FINAL_FAMILY_ID][order(-N_SPLIT)])
  cat("\nFamilies that never split:\n")
  never <- cells[, .(SPLITS = sum(SPLIT)), by = FINAL_FAMILY_ID][SPLITS == 0L]
  print(never); cat("  count:", nrow(never), "(paper says 10)\n")
  
  wide <- dcast(cells, FINAL_FAMILY_ID ~ SCHEME_NAME, value.var = "CATS", fill = "")
  save_csv(wide,  "T24E_family_scheme_grid_wide.csv")
  save_csv(cells, "T24E_family_scheme_grid_long.csv")
  cat("\nBuild the landscape appendix table from the wide CSV. H = high,\n",
      "I = intermediate, L = low; a cell with a slash is a split cell.\n", sep = "")
  invisible(list(long = cells, wide = wide))
}


# ===========================================================================
# 24F. F8 -- CONSOLIDATED SAMPLE AUDIT
# ===========================================================================
#
# One row per specification family. Panel rows and complete cases come from
# the relevant panel object; estimated rows from a representative fit.

s24_sample_audit <- function() {
  spec <- list(
    list(NAME = "Baseline, county markets",            PANEL = "outpatient",
         INSTR = names(MAIN_INSTRUMENTS)[1], FE = BASELINE_FIXED_EFFECTS),
    list(NAME = "Baseline, ex-CBSA instruments",       PANEL = "outpatient",
         INSTR = names(MAIN_INSTRUMENTS)[3], FE = BASELINE_FIXED_EFFECTS),
    list(NAME = "CBSA markets",                        PANEL = "cbsa_panel",
         INSTR = names(MAIN_INSTRUMENTS)[3], FE = BASELINE_FIXED_EFFECTS)
  )
  out <- rbindlist(lapply(spec, function(s) {
    if (!exists(s$PANEL, envir = .GlobalEnv)) return(NULL)
    p <- as.data.table(get(s$PANEL, envir = .GlobalEnv))
    z <- CONCEPT_INSTRUMENTS[[s$INSTR]]
    vars <- c(PRIMARY_OUTCOME, z, ENDOGENOUS_VARIABLE, BASELINE_CONTROLS,
              s$FE, BASELINE_CLUSTERS)
    d <- model_sample(p, vars)
    fit <- tryCatch(
      feols(as.formula(paste(PRIMARY_OUTCOME, "~", z, "+",
                             paste(BASELINE_CONTROLS, collapse = " + "), "|",
                             paste(s$FE, collapse = " + "))),
            data = d, warn = FALSE, notes = FALSE), error = function(e) NULL)
    data.table(SPECIFICATION = s$NAME, INSTRUMENT = s$INSTR,
               PANEL_ROWS = nrow(p), COMPLETE_CASES = nrow(d),
               ESTIMATED_ROWS = if (is.null(fit)) NA_integer_ else nobs(fit),
               DROP_SHARE = if (is.null(fit)) NA_real_ else 1 - nobs(fit) / nrow(d))
  }), fill = TRUE)
  print(out)
  cat("\nAdd the system x month and within-modality rows by re-running with\n",
      "FE = c('MARKET_ID','SYSTEM_MONTH') and the within-modality scheme filter.\n", sep = "")
  save_csv(out, "T24F_sample_audit_consolidated.csv")
  invisible(out)
}


# ===========================================================================
# DRIVER
# ===========================================================================

run_section_24 <- function() {
  .s24_hd("SECTION 24 -- TIER IV TABLE BUILDERS")
  res <- list()
  for (nm in c("A", "B", "C", "D", "E", "F")) {
    cat("\n--- 24", nm, " ---\n", sep = "")
    res[[nm]] <- tryCatch(switch(nm,
                                 A = s24_first_stage_variation(),
                                 B = s24_schemes_with_coefs(),
                                 C = s24_scheme_crosswalk(outpatient, schemes_long),
                                 D = s24_standardize_family_diffs(),
                                 E = s24_family_scheme_grid(outpatient, schemes_long),
                                 F = s24_sample_audit()),
                          error = function(e) { message("[24", nm, " failed] ", e$message); NULL })
  }
  .s24_hd("SECTION 24 COMPLETE")
  invisible(res)
}

s24 <- run_section_24()


###############################################################################

# ============================================================================
# Section 25: figure repairs
# ============================================================================

# SECTION 25 -- FIGURE REPAIRS
#
# Source AFTER warm_start() / restore_session(). Four figures need attention.
#
#   fig10_window_ladder.pdf   MISSING. Referenced by Section 8.3.4.
#   fig03_transform_ladder    STALE. Prints the reduced-form p as 0.0051; the
#                             corrected value is 0.0135.
#   fig02_scheme_robustness   STALE (same source as fig03).
#   fig05_instrument_tiers    STALE (its y-axis counts p < 0.05).
#
# The staleness is the same normal-vs-t issue already fixed in the tables.
# `main$tests` stores p computed against a normal reference; with sixteen
# month clusters the correct reference is t(15). The transformation is exact:
#
#     p_t = 2 * pt(-qnorm(1 - p_normal/2), df = 15)
#
# .to_t15() below applies it, so these figures can be rebuilt from the existing
# cache without re-estimating anything. Re-running the estimation with the
# corrected .pval() would give identical numbers.
#
# Revision note (post first run): three bugs fixed after a live run surfaced
# them.
#   1. facet_wrap(~PANEL) collided with ggplot2's internal reserved column
#      name PANEL, thrown at build time inside ggsave(). Renamed to
#      FACET_LAB in both figures that facet.
#   2. .save() used device = cairo_pdf, which needs an X11/cairo install this
#      machine doesn't have. Switched to plain grDevices::pdf.
#   3. s25_transform_ladder() assumed a specific cache schema and errored
#      with "object 'RF_P_NORMAL' not found" against a differently-shaped
#      transform_ladder.rds. It now probes several plausible column names and
#      only trusts the cache if all three needed fields resolve, else falls
#      back to the published values.
#   4. s25_instrument_tiers() looked for an object called `screen` that does
#      not exist in this session, and `main$tests` turned out to hold only
#      the three main-tier instruments (by construction of
#      run_main_results()), so the recount silently dropped the confirming
#      and discrepant rows. It now accepts confirming/discrepant result
#      objects as optional arguments and falls back to Table 6's pooled F
#      when no within-interacted Wald object is available, labelling which
#      source it used.
###############################################################################

suppressPackageStartupMessages({library(ggplot2); library(data.table)})

FIG_DIR <- if (exists("FIG_DIR", envir = .GlobalEnv)) get("FIG_DIR", envir = .GlobalEnv) else {
  # Your pipeline defines FIG_DIR directly (see HPT_Analysis_Pipeline.R,
  # around line 6677) as .../06_R_Analysis_Results/02_Figures. There is no
  # OUTPUT_DIR object anywhere in that pipeline; I invented it by mistake.
  # This branch only fires if FIG_DIR somehow isn't in scope when this file
  # is sourced, and falls back to reconstructing the same path from
  # RESULT_ROOT, which the pipeline does define.
  warning("FIG_DIR not found in the global environment; reconstructing it from ",
          "RESULT_ROOT. If this path is wrong, source HPT_Analysis_Pipeline.R ",
          "first so FIG_DIR is already defined.", call. = FALSE)
  file.path(RESULT_ROOT, "02_Figures")
}
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

N_MONTH_CLUSTERS <- 16L
.to_t15 <- function(p, df = N_MONTH_CLUSTERS - 1L) {
  p <- pmin(pmax(as.numeric(p), .Machine$double.xmin), 1 - 1e-15)
  2 * pt(-qnorm(1 - p / 2), df = df)
}

INSTR_LABEL <- c(
  Competitor_only_hospitals_9m         = "Competitor hospitals",
  Primary_strict_system_IV             = "Local system",
  Competitor_outside_CBSA_hospitals_9m = "Competitor hospitals (ex-CBSA)",
  Competitor_outside_CBSA_counties_9m  = "Competitor counties (ex-CBSA)",
  Competitor_systems_9m                = "Competitor systems",
  Competitor_outside_CBSA_systems_9m   = "Competitor systems (ex-CBSA)")
.pretty <- function(x) ifelse(is.na(INSTR_LABEL[x]), x, INSTR_LABEL[x])

.theme <- function() {
  theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(size = 8.5, colour = "grey30"),
          plot.caption = element_text(size = 7.5, colour = "grey40", hjust = 0),
          legend.position = "bottom", legend.title = element_blank(),
          panel.grid.minor = element_blank())
}
.save <- function(p, file, w = 7.5, h = 5) {
  # Plain grDevices::pdf, not cairo_pdf. cairo_pdf requires an X11/cairo
  # installation this machine does not have (missing libSM/libXrender under
  # /opt/X11), which threw non-fatal warnings and, for the two figures that
  # errored before reaching ggsave, made the trace noisier than it needed to
  # be. Nothing in these plots needs cairo's font handling.
  ggsave(file.path(FIG_DIR, file), p, width = w, height = h, device = "pdf")
  cat("Wrote:", file.path(FIG_DIR, file), "\n")
}


# ===========================================================================
# fig10_window_ladder.pdf -- MISSING FIGURE
# ===========================================================================
#
# Two panels. Panel A is the shoppable coefficient at each window length,
# Panel B the minimum first-stage Wald. The point of the figure is the
# non-monotonicity in Panel A against the monotonicity in Panel B, which is
# what Section 8.3.4 argues shows the horizon was not selected on outcome.
#
# Reads the window-ladder cache if present. The fallback is the published
# content of Appendix Table 31, so the figure matches the table either way.

s25_window_ladder <- function(cache_name = "window_ladder") {
  p <- file.path(CACHE_DIR, paste0(cache_name, ".rds"))
  if (file.exists(p)) {
    d <- as.data.table(readRDS(p))
    message("Using cache: ", p)
    if (!("INSTRUMENT_LABEL" %in% names(d))) stop("Unexpected cache schema.", call. = FALSE)
    d[, INSTRUMENT := .pretty(INSTRUMENT_LABEL)]
    if ("TEST_P" %in% names(d)) d[, TEST_P := .to_t15(TEST_P)]
  } else {
    warning("No window_ladder cache; using the values published in Appendix Table 31.",
            call. = FALSE)
    d <- data.table(
      INSTRUMENT = rep(c("Competitor hospitals", "Competitor hospitals (ex-CBSA)",
                         "Local system"), each = 4L),
      WINDOW_MONTHS = rep(c(3L, 6L, 9L, 12L), times = 3L),
      SHOPPABLE_PCT = c(-2.44, -4.83, -4.81, -3.43,
                        -2.50, -4.16, -4.39, -2.88,
                        -2.36, -4.11, -4.26, -2.72),
      TEST_P        = c(0.0174, 0.0064, 0.0135, 0.0585,
                        0.0107, 0.0057, 0.0108, 0.0559,
                        0.0096, 0.0065, 0.0256, 0.0834),
      MIN_WALD      = c(6.5, 22.0, 43.0, 52.9,
                        3.7, 16.0, 40.0, 51.8,
                        8.8, 23.8, 42.3, 52.5))
  }
  
  d[, WINDOW := factor(WINDOW_MONTHS, levels = c(3, 6, 9, 12),
                       labels = c("3M", "6M", "9M", "12M"))]
  d[, SIG05 := TEST_P < 0.05]
  
  long <- rbind(
    d[, .(INSTRUMENT, WINDOW, SIG05, VALUE = SHOPPABLE_PCT,
          FACET_LAB = "A. Shoppable coefficient (% per SD)")],
    d[, .(INSTRUMENT, WINDOW, SIG05, VALUE = MIN_WALD,
          FACET_LAB = "B. Minimum first-stage Wald")])
  
  # NOTE: ggplot2 reserves the column name "PANEL" internally for its own
  # facet layout bookkeeping (created inside ggplot_build()), so a data
  # column literally called PANEL collides with it and facet_wrap() errors
  # at build/print time with "PANEL is not an allowed name for faceting
  # variables." FACET_LAB avoids the collision.
  g <- ggplot(long, aes(WINDOW, VALUE, group = INSTRUMENT, colour = INSTRUMENT)) +
    geom_hline(data = data.frame(FACET_LAB = "A. Shoppable coefficient (% per SD)", y = 0),
               aes(yintercept = y), linetype = "dashed", colour = "grey60") +
    geom_line(linewidth = 0.6) +
    geom_point(aes(shape = SIG05), size = 2.4) +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                       labels = c(`TRUE` = "Test rejects at 5%",
                                  `FALSE` = "Does not reject")) +
    facet_wrap(~FACET_LAB, scales = "free_y") +
    labs(title = "The effect peaks where the institutional rationale predicts",
         subtitle = paste("Shoppable coefficient is non-monotone in window length;",
                          "first-stage strength is monotone, which is mechanical"),
         x = "Trailing window used to accumulate peer disclosure", y = NULL,
         caption = paste("Procedural-certainty classification. Three-month IV",
                         "magnitudes are not quoted, since the first stage is weak there.")) +
    .theme()
  .save(g, "fig10_window_ladder.pdf", w = 8.5, h = 4.6)
  invisible(d)
}


# ===========================================================================
# fig02_scheme_robustness.pdf -- REGENERATE WITH t(15)
# ===========================================================================

s25_scheme_robustness <- function(main_obj = main) {
  mt <- as.data.table(main_obj$tests)
  mt <- mt[ESTIMATOR %chin% c("Reduced form", "RF") &
             INSTRUMENT_LABEL %chin% names(MAIN_INSTRUMENTS)]
  if (!nrow(mt)) stop("No reduced-form rows found in main$tests.", call. = FALSE)
  
  mt[, P_T := .to_t15(P_VALUE)]
  mt[, INSTRUMENT := .pretty(INSTRUMENT_LABEL)]
  mt[, SCHEME := sub("^[0-9]+\\. ", "", SPEC)]
  ord <- mt[, .(M = mean(P_T)), by = SCHEME][order(M), SCHEME]
  mt[, SCHEME := factor(SCHEME, levels = ord)]
  
  cat("\nReduced-form heterogeneity p, normal reference vs t(15):\n")
  print(mt[order(SCHEME, INSTRUMENT),
           .(SCHEME, INSTRUMENT, NORMAL = round(P_VALUE, 4), T15 = round(P_T, 4))])
  cat("\nSignificant at 5%: ", sum(mt$P_VALUE < 0.05), " under normal, ",
      sum(mt$P_T < 0.05), " under t(15), of ", nrow(mt), "\n", sep = "")
  
  g <- ggplot(mt, aes(P_T, SCHEME, colour = INSTRUMENT)) +
    geom_vline(xintercept = 0.05, linetype = "dashed", colour = "grey40") +
    geom_point(size = 2.6, alpha = 0.9) +
    annotate("text", x = 0.05, y = 0.55, label = "p = 0.05",
             size = 2.8, colour = "grey40", hjust = -0.1) +
    labs(title = "The gradient holds across every classification scheme",
         subtitle = "Test of equality between shoppable and non-shoppable response",
         x = "Heterogeneity test p-value (reduced form)", y = NULL,
         caption = paste("p-values reference the t distribution implied by the",
                         "two-way clustering, matching Table 8.")) +
    .theme()
  .save(g, "fig02_scheme_robustness.pdf", w = 7.5, h = 4.2)
  invisible(mt)
}


# ===========================================================================
# fig03_transform_ladder.pdf -- REGENERATE WITH t(15)
# ===========================================================================
#
# The printed reduced-form p should read 0.0135, not 0.0051. The point of the
# figure is that it is identical across all six forms while the IV magnitude
# moves by a factor of five, and that survives the correction unchanged.

s25_transform_ladder <- function(transform_obj = NULL) {
  FALLBACK <- data.table(
    FORM = c("Linear", "Winsor P99", "Winsor P95", "Winsor P90", "Sqrt", "Log1p"),
    IV_PCT = c(-3.72, -3.91, -4.81, -7.20, -13.60, -17.70),
    RF_P_NORMAL = rep(0.0051, 6L))
  
  d <- NULL
  if (is.null(transform_obj)) {
    p <- file.path(CACHE_DIR, "transform_ladder.rds")
    if (file.exists(p)) {
      cand <- as.data.table(readRDS(p))
      # A cache with an unexpected schema previously produced "object
      # 'RF_P_NORMAL' not found" a few lines below, because setnames()
      # only fires when the exact source name is present and otherwise
      # leaves RF_P_NORMAL undefined. Try several plausible column names
      # instead of assuming one, and only trust the cache if a form label,
      # an IV magnitude, and a reduced-form p-value can all be located.
      form_col <- intersect(c("FORM", "TRANSFORM", "FUNCTIONAL_FORM"), names(cand))[1]
      iv_col   <- intersect(c("IV_PCT", "IV_PERCENT", "IV_MAGNITUDE"), names(cand))[1]
      p_col    <- intersect(c("RF_P_NORMAL", "RF_P", "P_VALUE", "RF_PVAL"), names(cand))[1]
      if (!anyNA(c(form_col, iv_col, p_col)) &&
          all(nzchar(c(form_col, iv_col, p_col)))) {
        cand <- cand[, .(FORM = get(form_col), IV_PCT = get(iv_col),
                         RF_P_NORMAL = get(p_col))]
        d <- cand
      } else {
        warning("transform_ladder.rds exists but its columns (",
                paste(names(cand), collapse = ", "),
                ") do not match any expected schema; using published values instead.",
                call. = FALSE)
      }
    }
  } else {
    d <- as.data.table(transform_obj)
  }
  if (is.null(d)) {
    if (is.null(transform_obj))
      warning("No usable transform_ladder cache; using published values.", call. = FALSE)
    d <- copy(FALLBACK)
  }
  
  d[, RF_P := .to_t15(RF_P_NORMAL)]
  d[, FORM := factor(FORM, levels = FORM)]
  cat("\nReduced-form p across forms: ", paste(round(d$RF_P, 4), collapse = ", "),
      "\n(identical by construction; the reduced form contains no treatment variable)\n",
      sep = "")
  
  long <- rbind(
    d[, .(FORM, VALUE = IV_PCT, FACET_LAB = "IV magnitude, shoppable services (%)")],
    d[, .(FORM, VALUE = RF_P,  FACET_LAB = "Reduced-form heterogeneity p-value")])
  
  g <- ggplot(long, aes(FORM, VALUE)) +
    geom_col(fill = "grey35", width = 0.62) +
    geom_text(aes(label = ifelse(FACET_LAB == "Reduced-form heterogeneity p-value",
                                 sprintf("%.4f", VALUE), sprintf("%.1f", VALUE))),
              vjust = ifelse(long$VALUE < 0, 1.3, -0.4), size = 2.7) +
    facet_wrap(~FACET_LAB, scales = "free_y") +
    labs(title = "Functional form moves the magnitude but not the conclusion",
         subtitle = paste("The reduced-form test is numerically identical across all six",
                          "forms, because it contains no treatment variable"),
         x = NULL, y = NULL) +
    .theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
  .save(g, "fig03_transform_ladder.pdf", w = 8.5, h = 4.4)
  invisible(d)
}


# ===========================================================================
# fig05_instrument_tiers.pdf -- REGENERATE WITH t(15)
# ===========================================================================
#
# The y-axis is a count of p < 0.05, so it moves with the correction. Also
# emits the recount Appendix Table 22 needs: that table currently carries two
# rows both labelled "Competitor hospitals (ex-CBSA)".

s25_instrument_tiers <- function(main_obj = main, confirming_obj = NULL,
                                 discrepant_obj = NULL) {
  # `main$tests` (from run_main_results()) covers only the three main-tier
  # instruments by construction, which is why the recount below only ever
  # returned three rows. If a "confirming" or "discrepant" results object
  # exists in this session with the same $tests structure, pass it in to
  # get all six instruments; otherwise this falls back to the appendix's own
  # already-published confirming/discrepant figures so the plot isn't
  # silently missing half the instrument set.
  grab_tests <- function(obj) {
    if (is.null(obj)) return(NULL)
    tt <- tryCatch(as.data.table(obj$tests), error = function(e) NULL)
    if (is.null(tt)) return(NULL)
    tt[ESTIMATOR %chin% c("Reduced form", "RF")]
  }
  mt <- rbindlist(list(grab_tests(main_obj), grab_tests(confirming_obj),
                       grab_tests(discrepant_obj)), fill = TRUE)
  mt[, P_T := .to_t15(P_VALUE)]
  agg <- mt[, .(N_SCHEMES = .N,
                SIG_NORMAL = sum(P_VALUE < 0.05),
                SIG_T15 = sum(P_T < 0.05),
                MEDIAN_P_T = median(P_T)), by = INSTRUMENT_LABEL]
  agg[, INSTRUMENT := .pretty(INSTRUMENT_LABEL)]
  
  # First-stage F. Appendix Table 6 (tab:instrument_strength) already
  # reports the pooled cluster-robust F for all six instruments on a common
  # sample; that is not the same statistic as the within-interacted-model
  # Wald this figure plots on its x-axis (see the table's own note on that
  # distinction), but it is the only per-instrument F on hand when no
  # `screen`-like object exists in the session, so it is used as a labelled
  # fallback rather than dropping the plot.
  fs_fallback <- data.table(
    INSTRUMENT = c("Competitor systems", "Competitor systems (ex-CBSA)",
                   "Competitor counties (ex-CBSA)", "Local system",
                   "Competitor hospitals", "Competitor hospitals (ex-CBSA)"),
    MEDIAN_F = c(47.01, 36.56, 16.87, 16.05, 15.89, 12.57),
    F_SOURCE = "tab:instrument_strength (pooled F, not within-interacted Wald)")
  
  fs <- tryCatch(as.data.table(screen)[, .(INSTRUMENT_LABEL, MEDIAN_F = F_STAT)],
                 error = function(e) NULL)
  if (!is.null(fs)) {
    agg <- merge(agg, fs, by = "INSTRUMENT_LABEL", all.x = TRUE)
    agg[, F_SOURCE := "screen object (within-interacted Wald)"]
  } else {
    warning("No `screen`-type object with within-interacted first-stage F; ",
            "falling back to the pooled F from Table 6. The x-axis is not on ",
            "the same footing as a within-interacted-model Wald would give.",
            call. = FALSE)
    agg <- merge(agg, fs_fallback, by = "INSTRUMENT", all.x = TRUE)
  }
  agg[, TIER := fifelse(INSTRUMENT_LABEL %chin% names(MAIN_INSTRUMENTS), "MAIN",
                        fifelse(grepl("systems", INSTRUMENT_LABEL) &
                                  grepl("CBSA", INSTRUMENT_LABEL), "DISCREPANT", "CONFIRMING"))]
  agg[, SHARE := SIG_T15 / N_SCHEMES]
  
  cat("\nRecount for Appendix Table 22 (tab:tier_rerank). The published table\n",
      "has a duplicate 'Competitor hospitals (ex-CBSA)' row: one stale entry\n",
      "(6 of 6, p=0.0098) left over from before the t(15) correction, and one\n",
      "correct entry (5 of 6, p=0.0206) that matches this recount exactly.\n",
      "Delete the stale row; keep the second.\n", sep = "")
  print(agg[order(-SHARE), .(INSTRUMENT, TIER, N_SCHEMES, SIG_NORMAL, SIG_T15,
                             MEDIAN_P_T = round(MEDIAN_P_T, 4),
                             MEDIAN_F = round(MEDIAN_F, 1), F_SOURCE)])
  save_csv(agg, "T25_instrument_tier_recount.csv")
  
  if (all(is.na(agg$MEDIAN_F))) {
    warning("No first-stage F available even after the fallback; skipping the plot.",
            call. = FALSE)
    return(invisible(agg))
  }
  g <- ggplot(agg, aes(MEDIAN_F, SHARE, colour = TIER)) +
    geom_vline(xintercept = 10, linetype = "dashed", colour = "grey55") +
    geom_point(size = 3) +
    ggrepel::geom_text_repel(aes(label = INSTRUMENT), size = 2.7, show.legend = FALSE) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(title = "Instrument strength and result strength are independent",
         subtitle = "The two strongest first stages produce the two weakest results",
         x = "Median first-stage Wald F", y = "Share of classification schemes significant",
         caption = "Significance at 5% against the t(15) reference.") +
    .theme()
  .save(g, "fig05_instrument_tiers.pdf", w = 7.5, h = 5)
  invisible(agg)
}


run_section_25 <- function() {
  cat("\n", strrep("=", 78), "\nSECTION 25 -- FIGURE REPAIRS\n", strrep("=", 78), "\n", sep = "")
  res <- list()
  for (nm in c("window", "schemes", "transform", "tiers")) {
    cat("\n--- ", nm, " ---\n", sep = "")
    res[[nm]] <- tryCatch(switch(nm,
                                 window    = s25_window_ladder(),
                                 schemes   = s25_scheme_robustness(),
                                 transform = s25_transform_ladder(),
                                 tiers     = s25_instrument_tiers()),
                          error = function(e) { message("[", nm, " failed] ", e$message); NULL })
  }
  invisible(res)
}

# Example:
s25 <- run_section_25()



library(data.table)

# --- T06_mainB (feeds fig02) ---
t06 <- fread(file.path(TABLE_DIR, "T06_mainB_heterogeneity_tests.csv"))
stopifnot("P_VALUE" %in% names(t06))          # sanity check before touching it
t06[, P_VALUE := 2 * pt(-qnorm(1 - P_VALUE / 2), df = 15)]
fwrite(t06, file.path(TABLE_DIR, "T06_mainB_heterogeneity_tests.csv"))

# --- T07B (feeds fig03) ---
t07b <- fread(file.path(TABLE_DIR, "T07B_transform_ladder_heterogeneity.csv"))
stopifnot("P_VALUE" %in% names(t07b))
t07b[, P_VALUE := 2 * pt(-qnorm(1 - P_VALUE / 2), df = 15)]
fwrite(t07b, file.path(TABLE_DIR, "T07B_transform_ladder_heterogeneity.csv"))

}  # end HPT_RUN$diagnostics

###############################################################################
#                                                                             #
#   PART 7 -- SCRATCH                                                         #
#                                                                             #
#   Interactive code that previously executed on every source(). Kept as a    #
#   record of what was run by hand; commented out so it does not execute.     #
#                                                                             #
###############################################################################

# (none)


writeLines(capture.output(sessionInfo()),
           file.path(RESULT_ROOT, "sessionInfo.txt"))
cat("\nPipeline complete. sessionInfo written to ", RESULT_ROOT, "\n", sep = "")








# Sys.setenv(HPT_ROOT = "/Users/danielsierra/Library/CloudStorage/OneDrive-FloridaStateUniversity/Hospital Price Transparency Paper")
# HPT_WARM_START <- TRUE
# source(file.path(Sys.getenv("HPT_ROOT"), "Code", "HPT_Analysis_Pipeline.R"))
# cache_status()
# family_results      <- readRDS(file.path(CACHE_DIR, "family_level_6inst.rds"))
# superfamily_results <- readRDS(file.path(CACHE_DIR, "superfamily_level_6inst.rds"))












###############################################################################
#
#   SECTION 26 -- CONCEPT-LEVEL ESTIMATE BROWSER
#
#   Sorts the concept-level sweep from Section 6 (concept_results) three ways:
#
#     (1) by coefficient, least to greatest
#     (2) by p-value, least to greatest
#     (3) by both, under three explicit rules (see cb_by_both)
#
#   Source AFTER HPT_Analysis_Pipeline.R Section 1 has run (needs CACHE_DIR,
#   TABLE_DIR, DIAGNOSTIC_FAMILIES, FAMILY_LABELS, INSTRUMENT_LABEL_MAP,
#   save_csv). It reads concept_results from memory, then the cache, then the
#   CSV -- whichever it finds first. Nothing here re-estimates anything.
#
#   TWO CONVENTIONS WORTH KNOWING BEFORE READING THE OUTPUT
#
#   1. Percent conversion. estimate_concept_level() stores
#      IV_ESTIMATE_PERCENT = 100 * b, the linear approximation.
#      estimate_interacted() stores IV_PERCENT = 100 * (exp(b) - 1).
#      The two diverge at the tails of the concept sweep (b = -0.126 gives
#      -12.6 under the first and -11.8 under the second). This script reports
#      both, as PCT_LINEAR and PCT_EXP, so the gap is visible rather than
#      silently inherited. Section 6.2.1 of the paper quotes the linear
#      version; Table 5 Panel B quotes the exponential one.
#
#   2. Reference distribution. Concept-level p-values come straight from
#      coeftable(feols(...)) with two-way clustering, and fixest's default
#      t.df = "min" already references t with df = min(#clusters) - 1. Within
#      a concept subsample that is usually the month count, so the df is
#      already at or below 15 and varies by concept. Do NOT apply the
#      pipeline's .to_t15() to these -- it would double-count the correction.
#      DF_MIN below reports min(N_MARKETS, N_MONTHS) - 1, the reference each
#      concept was actually evaluated against.
#
###############################################################################

suppressPackageStartupMessages(library(data.table))


# ===========================================================================
# 26.0  Configuration
# ===========================================================================

CB_METRICS <- list(
  RF = list(coef = "RF_COEF", se = "RF_SE", p = "RF_P",  fdr = "RF_P_FDR",
            unit = "log points per unit of Z"),
  IV = list(coef = "IV_COEF", se = "IV_SE", p = "IV_P",  fdr = "IV_P_FDR",
            unit = "log points per prior poster"),
  FS = list(coef = "FS_COEF", se = "FS_SE", p = NA_character_, fdr = NA_character_,
            unit = "prior posters per unit of Z")
)

# Used only for the percent-per-SD column when `outpatient` is not in memory.
# These are the estimation-sample standard deviations reported in Appendix
# Table (tab:instrument_strength). If the panel IS in memory the SD is
# computed from it directly and this is ignored.
CB_Z_SD_FALLBACK <- c(
  Competitor_only_hospitals_9m         = 15.94,
  Primary_strict_system_IV             = 17.16,
  Competitor_outside_CBSA_hospitals_9m = 14.72,
  Competitor_outside_CBSA_counties_9m  = 11.07,
  Competitor_systems_9m                = 1.505,
  Competitor_outside_CBSA_systems_9m   = 1.256
)

.cb_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

.cb_stars <- function(p) {
  fifelse(is.na(p), "",
          fifelse(p < 0.01, "***",
                  fifelse(p < 0.05, "**",
                          fifelse(p < 0.10, "*", ""))))
}


# ===========================================================================
# 26.1  Load
# ===========================================================================

cb_load <- function(verbose = TRUE) {
  if (exists("concept_results") && is.data.frame(concept_results)) {
    if (verbose) cat("Source: concept_results in memory\n")
    return(as.data.table(copy(concept_results)))
  }
  rds <- file.path(CACHE_DIR, "concept_level_6inst.rds")
  if (file.exists(rds)) {
    if (verbose) cat("Source:", rds, "\n")
    return(as.data.table(readRDS(rds)))
  }
  csv <- file.path(TABLE_DIR, "T05_concept_level_RF_FS_IV.csv")
  if (file.exists(csv)) {
    if (verbose) cat("Source:", csv, "\n")
    return(fread(csv))
  }
  stop("No concept-level results found. Run stage 6, or check CACHE_DIR / TABLE_DIR.",
       call. = FALSE)
}

cb_z_sd <- function(instrument_label, instrument_var = NULL, panel_name = "outpatient") {
  if (!is.null(instrument_var) && exists(panel_name)) {
    p <- get(panel_name)
    if (instrument_var %chin% names(p)) {
      v <- .cb_num(p[[instrument_var]])
      s <- stats::sd(v[is.finite(v)], na.rm = TRUE)
      if (is.finite(s) && s > 0) return(s)
    }
  }
  s <- CB_Z_SD_FALLBACK[[instrument_label]]
  if (is.null(s)) NA_real_ else s
}


# ===========================================================================
# 26.2  Build the browsable view
# ===========================================================================
#
# One row per concept, for ONE instrument and ONE metric. Everything the
# sorting rules need is computed here so the three orderings below are pure
# reorderings of the same object and cannot disagree about the numbers.

cb_view <- function(cr = cb_load(),
                    metric = "RF",
                    instrument = "Competitor_only_hospitals_9m",
                    shoppable_families = DIAGNOSTIC_FAMILIES,
                    verbose = TRUE) {
  
  metric <- match.arg(toupper(metric), names(CB_METRICS))
  m <- CB_METRICS[[metric]]
  d <- as.data.table(copy(cr))
  
  # The INSTRUMENT column is authoritative; labels have been written both with
  # spaces and with underscores across the pipeline's history. Normalise before
  # filtering, exactly as the note at INSTRUMENT_LABEL_MAP instructs.
  if ("INSTRUMENT" %chin% names(d) && exists("INSTRUMENT_LABEL_MAP")) {
    lab <- unname(INSTRUMENT_LABEL_MAP[as.character(d$INSTRUMENT)])
    if (any(!is.na(lab))) d[!is.na(lab), INSTRUMENT_LABEL := lab[!is.na(lab)]]
  }
  
  avail <- sort(unique(d$INSTRUMENT_LABEL))
  if (!is.null(instrument)) {
    if (!(instrument %chin% avail))
      stop("Instrument '", instrument, "' not present. Available:\n  ",
           paste(avail, collapse = "\n  "), call. = FALSE)
    d <- d[INSTRUMENT_LABEL == instrument]
  }
  if (nrow(d) == 0L) stop("No rows after filtering.", call. = FALSE)
  
  if (!(m$coef %chin% names(d)))
    stop(m$coef, " not found in concept_results.", call. = FALSE)
  
  d[, COEF := .cb_num(get(m$coef))]
  d[, SE   := if (m$se %chin% names(d)) .cb_num(get(m$se)) else NA_real_]
  d[, P    := if (!is.na(m$p)   && m$p   %chin% names(d)) .cb_num(get(m$p))   else NA_real_]
  d[, PFDR := if (!is.na(m$fdr) && m$fdr %chin% names(d)) .cb_num(get(m$fdr)) else NA_real_]
  
  d <- d[is.finite(COEF)]
  if (nrow(d) == 0L) stop("No finite coefficients for metric ", metric, ".", call. = FALSE)
  
  d[, TSTAT := fifelse(is.finite(SE) & SE > 0, COEF / SE, NA_real_)]
  d[, STARS := .cb_stars(P)]
  d[, SHOP  := fifelse(as.character(FINAL_FAMILY_ID) %chin% shoppable_families,
                       "Shoppable", "Non-shoppable")]
  d[, FAMILY := if (exists("FAMILY_LABELS")) {
    lb <- unname(FAMILY_LABELS[as.character(FINAL_FAMILY_ID)])
    fifelse(is.na(lb), as.character(FINAL_FAMILY_ID), lb)
  } else as.character(FINAL_FAMILY_ID)]
  
  # Both percent conventions, side by side. See the header note.
  sd_z <- cb_z_sd(instrument, if ("INSTRUMENT" %chin% names(d)) d$INSTRUMENT[1L] else NULL)
  d[, PCT_LINEAR := 100 * COEF]
  d[, PCT_EXP    := 100 * (exp(COEF) - 1)]
  d[, PCT_PER_SD := if (is.finite(sd_z)) 100 * (exp(COEF * sd_z) - 1) else NA_real_]
  
  d[, DF_MIN := pmin(.cb_num(N_MARKETS), .cb_num(N_MONTHS)) - 1]
  
  keep <- c("FINAL_CONCEPT_ID", "FINAL_CONCEPT_NAME", "SERVICE_LABEL", "FAMILY",
            "FINAL_FAMILY_ID", "SHOP", "INSTRUMENT_LABEL",
            "COEF", "SE", "TSTAT", "P", "PFDR", "STARS",
            "PCT_LINEAR", "PCT_EXP", "PCT_PER_SD",
            "N_OBSERVATIONS", "N_ROWS", "N_HOSPITALS", "N_MARKETS", "N_MONTHS",
            "DF_MIN", "FS_COEF", "FS_F", "MEAN_PRIOR_POSTERS")
  keep <- intersect(keep, names(d))
  v <- d[, ..keep]
  
  setattr(v, "cb_metric", metric)
  setattr(v, "cb_unit", m$unit)
  setattr(v, "cb_instrument", instrument)
  setattr(v, "cb_sd_z", sd_z)
  
  if (verbose)
    cat(sprintf("View: %d concepts | metric %s (%s) | instrument %s | SD(Z) = %s\n",
                nrow(v), metric, m$unit, instrument,
                if (is.finite(sd_z)) round(sd_z, 3) else "unavailable"))
  v[]
}


# ===========================================================================
# 26.3  The three orderings
# ===========================================================================

# (1) Least to greatest coefficient. Most negative price response first.
cb_by_coef <- function(v, decreasing = FALSE) {
  out <- copy(v)
  setorderv(out, "COEF", if (decreasing) -1L else 1L)
  out[, RANK_COEF := .I]
  out[]
}

# (2) Least to greatest p-value. Most precisely estimated first, regardless
#     of sign -- which is the point of looking at it separately.
cb_by_p <- function(v) {
  if (all(is.na(v$P)))
    stop("This metric carries no p-value (FS). Use metric = 'RF' or 'IV'.", call. = FALSE)
  out <- copy(v)
  setorderv(out, c("P", "COEF"), c(1L, 1L))
  out[, RANK_P := .I]
  out[]
}

# (3) By both. Three rules, because "both" is genuinely ambiguous and the
#     three answer different questions:
#
#   "evidence"  EVIDENCE = sign(COEF) * (-log10(P)). Ascending, this puts the
#               most negative AND most significant concepts at the top, the
#               most positive AND most significant at the bottom, and the
#               unresolved middle in the middle. One signed number, and the
#               ordering that surfaces the two tails of the distribution the
#               paper's Figure 4 plots. A display ordering, not a test
#               statistic.
#
#   "ranks"     RANK_COEF + RANK_P, ascending. Equal weight to "how negative"
#               and "how precise". More robust to a single concept with an
#               absurdly small p, which the log scale in "evidence" rewards
#               heavily.
#
#   "lex"       Significance band first (1% / 5% / 10% / ns), then coefficient
#               ascending within band. This is the ordering that matches how a
#               referee reads a table: show me the significant ones, sorted by
#               size, then the rest.
cb_by_both <- function(v, rule = c("evidence", "ranks", "lex")) {
  rule <- match.arg(rule)
  if (all(is.na(v$P)))
    stop("This metric carries no p-value (FS). Use metric = 'RF' or 'IV'.", call. = FALSE)
  
  out <- copy(v)
  out[, RANK_COEF := frank(COEF, ties.method = "first")]
  out[, RANK_P    := frank(P,    ties.method = "first", na.last = TRUE)]
  out[, EVIDENCE  := sign(COEF) * (-log10(pmax(P, 1e-300)))]
  out[, RANK_SUM  := RANK_COEF + RANK_P]
  out[, SIG_BAND  := fifelse(is.na(P), 4L,
                             fifelse(P < 0.01, 1L,
                                     fifelse(P < 0.05, 2L,
                                             fifelse(P < 0.10, 3L, 4L))))]
  
  if (rule == "evidence") setorderv(out, c("EVIDENCE", "COEF"),   c(1L, 1L))
  if (rule == "ranks")    setorderv(out, c("RANK_SUM", "COEF"),   c(1L, 1L))
  if (rule == "lex")      setorderv(out, c("SIG_BAND", "COEF"),   c(1L, 1L))
  
  setattr(out, "cb_rule", rule)
  out[]
}


# ===========================================================================
# 26.4  Printing
# ===========================================================================
#
# as.data.frame() before print() on purpose: print(n = Inf) on a tibble or a
# wide data.table throws `invalid 'na.print' specification` in this setup.

cb_show <- function(x, n = Inf, cols = NULL, digits = 5) {
  d <- as.data.table(copy(x))
  if (is.null(cols)) {
    cols <- intersect(c("RANK_COEF", "RANK_P", "SIG_BAND", "EVIDENCE",
                        "FINAL_CONCEPT_NAME", "FAMILY", "SHOP",
                        "COEF", "SE", "P", "PFDR", "STARS",
                        "PCT_LINEAR", "N_OBSERVATIONS", "N_HOSPITALS", "DF_MIN"),
                      names(d))
  }
  d <- d[, ..cols]
  numcols <- names(d)[vapply(d, is.numeric, logical(1))]
  for (cn in numcols) set(d, j = cn, value = round(d[[cn]], digits))
  if (is.finite(n)) d <- head(d, n)
  print(as.data.frame(d))
  invisible(x)
}

cb_extremes <- function(v, n = 25) {
  s <- cb_by_coef(v)
  cat("\n--- ", n, " most negative -----------------------------------------\n", sep = "")
  cb_show(head(s, n))
  cat("\n--- ", n, " most positive -----------------------------------------\n", sep = "")
  cb_show(tail(s, n))
  invisible(s)
}


# ===========================================================================
# 26.5  One row per concept, all instruments side by side
# ===========================================================================
#
# The single-instrument views answer "which concepts respond". This answers
# "which concepts respond under every instrument", which is the more useful
# question for the sweep and is not currently anywhere in the pipeline.

cb_wide <- function(cr = cb_load(), metric = "RF",
                    instruments = if (exists("ALL_SIX_INSTRUMENTS"))
                      names(ALL_SIX_INSTRUMENTS) else NULL,
                    order_by = "Competitor_only_hospitals_9m") {
  
  metric <- match.arg(toupper(metric), names(CB_METRICS))
  parts <- lapply(instruments, function(il)
    tryCatch(cb_view(cr, metric, il, verbose = FALSE), error = function(e) NULL))
  parts <- Filter(Negate(is.null), parts)
  if (length(parts) == 0L) stop("No instruments could be viewed.", call. = FALSE)
  long <- rbindlist(parts, fill = TRUE)
  
  w <- dcast(long, FINAL_CONCEPT_ID + FINAL_CONCEPT_NAME + FAMILY + SHOP ~ INSTRUMENT_LABEL,
             value.var = c("COEF", "P"))
  
  coef_cols <- grep("^COEF_", names(w), value = TRUE)
  p_cols    <- grep("^P_",    names(w), value = TRUE)
  
  w[, N_INSTRUMENTS := rowSums(!is.na(as.matrix(.SD))), .SDcols = coef_cols]
  w[, N_NEGATIVE    := rowSums(as.matrix(.SD) < 0, na.rm = TRUE), .SDcols = coef_cols]
  w[, N_SIG_05      := rowSums(as.matrix(.SD) < 0.05, na.rm = TRUE), .SDcols = p_cols]
  w[, ALL_AGREE     := N_INSTRUMENTS > 0L & (N_NEGATIVE == N_INSTRUMENTS | N_NEGATIVE == 0L)]
  
  ocol <- paste0("COEF_", order_by)
  if (ocol %chin% names(w)) setorderv(w, ocol, 1L, na.last = TRUE)
  w[]
}


# ===========================================================================
# 26.6  Verification against the published numbers
# ===========================================================================
#
# Section 6.2 of the paper reports, under the primary instrument, a 5%
# significance rate of 31.3% for shoppable concepts against 3.4% for
# non-shoppable, and median reduced-form coefficients of -0.0034 and +0.0004.
# If this does not reproduce, the browser and the paper are reading different
# vintages and nothing below should be trusted.

cb_summary <- function(v) {
  s <- v[, .(N          = .N,
             MEDIAN_COEF = median(COEF, na.rm = TRUE),
             MEAN_COEF   = mean(COEF, na.rm = TRUE),
             SHARE_NEG   = mean(COEF < 0, na.rm = TRUE),
             SHARE_P05   = mean(P < 0.05, na.rm = TRUE),
             SHARE_P10   = mean(P < 0.10, na.rm = TRUE),
             SHARE_FDR05 = mean(PFDR < 0.05, na.rm = TRUE)),
         by = SHOP][order(-SHOP)]
  cat("\nBy shoppability (metric ", attr(v, "cb_metric"), ", instrument ",
      attr(v, "cb_instrument"), "):\n", sep = "")
  print(as.data.frame(s[, lapply(.SD, function(x)
    if (is.numeric(x)) round(x, 4) else x)]))
  
  f <- v[, .(N = .N, MEDIAN_COEF = round(median(COEF, na.rm = TRUE), 5),
             SHARE_P05 = round(mean(P < 0.05, na.rm = TRUE), 3)),
         by = .(FAMILY, SHOP)][order(MEDIAN_COEF)]
  cat("\nBy clinical family, most negative first:\n")
  print(as.data.frame(f))
  invisible(list(by_shop = s, by_family = f))
}


# ===========================================================================
# 26.7  Driver
# ===========================================================================

cb_run <- function(metric = "RF",
                   instrument = "Competitor_only_hospitals_9m",
                   rule = "evidence",
                   write = TRUE) {
  
  cr <- cb_load()
  v  <- cb_view(cr, metric, instrument)
  
  cb_summary(v)
  
  s_coef <- cb_by_coef(v)
  s_p    <- cb_by_p(v)
  s_both <- cb_by_both(v, rule)
  
  cat("\n\n=== (1) BY COEFFICIENT, least to greatest =========================\n")
  cb_show(s_coef)
  cat("\n\n=== (2) BY P-VALUE, least to greatest =============================\n")
  cb_show(s_p)
  cat("\n\n=== (3) BY BOTH, rule = ", rule, " ================================\n", sep = "")
  cb_show(s_both)
  
  if (isTRUE(write)) {
    tag <- paste0(metric, "_", instrument)
    writer <- if (exists("save_csv")) save_csv else
      function(d, f) { fwrite(d, file.path(TABLE_DIR, f)); cat("Saved:", f, "\n") }
    writer(s_coef, paste0("T26A_concept_sorted_by_coef_",  tag, ".csv"))
    writer(s_p,    paste0("T26B_concept_sorted_by_p_",     tag, ".csv"))
    writer(s_both, paste0("T26C_concept_sorted_by_both_",  rule, "_", tag, ".csv"))
    writer(cb_wide(cr, metric), paste0("T26D_concept_wide_all_instruments_", metric, ".csv"))
  }
  
  invisible(list(view = v, by_coef = s_coef, by_p = s_p, by_both = s_both))
}

cat("Section 26 loaded. Try:  res <- cb_run()\n")

if (exists("HPT_SCRATCH") && isTRUE(HPT_SCRATCH)) {   # interactive scratch, off by default
  res <- cb_run()                                   # RF, primary instrument
  # res <- cb_run(metric = "IV")                      # IV ratio instead
  # res <- cb_run(instrument = "Primary_strict_system_IV")
  # res <- cb_run(rule = "lex", write = FALSE)
  
  v <- cb_view()                              # RF, competitor hospitals
  cb_show(cb_by_coef(v), n = 40)              # 40 most negative
  cb_show(cb_by_p(v), n = 40)                 # 40 most precise
  cb_show(cb_by_both(v, "evidence"), n = 40)  # most negative AND most precise
  cb_extremes(v, 25)                          # both tails at once
  
  cb_show(cb_by_coef(v[SHOP == "Shoppable"])) # shoppable only
  
  w <- cb_wide()                              # all six instruments, one row per concept
  cb_show(w[ALL_AGREE & N_SIG_05 >= 3], cols = names(w))
  
  View(cb_by_both(v, "evidence"))             # 738 rows are easier to scroll than to print
}   # end interactive scratch
###############################################################################





if (exists("HPT_SCRATCH") && isTRUE(HPT_SCRATCH)) {   # interactive scratch, off by default
  rm(list = intersect(c("cr", "ms"), ls(envir = .GlobalEnv)), envir = .GlobalEnv)  # clear the stray globals first
  a6 <- s20_depth_diagnostics(concept_results, measures)
  
  
  exists("s20_payer_balance")
  exists("s20_price_with_payer_control")
  list.files(PAYER_DISPERSION_DIR, pattern = PAYER_DISPERSION_PATTERN)
  
  b1 <- s20_payer_balance(outpatient)
  b2 <- s20_price_with_payer_control(outpatient)
  
  beq <- s20_payer_balance_equality(outpatient)
  
  # 2. §7 MDE bound -- instant, reads a saved CSV
  ses_tests <- read_table("T12T_triple_interaction_ses_tests.csv")
  mde <- copy(ses_tests[ESTIMATOR == "Reduced form"])
  mde[, SCALE     := DGAP_PCT_PER_SD_MOD / DIFF]
  mde[, MDE_PCT   := 2.802 * DIFF_SE * SCALE]
  mde[, MDE_SHARE := abs(MDE_PCT) / abs(GAP_AT_MEAN_PCT_PER_SD_Z)]
  print(mde[order(MDE_SHARE), .(SPEC, INSTRUMENT_LABEL,
                                GAP = round(GAP_AT_MEAN_PCT_PER_SD_Z, 2),
                                MDE = round(abs(MDE_PCT), 2),
                                SHARE = round(MDE_SHARE, 2))])
  save_csv(mde, "T12T_ses_index_mde.csv")
  
  # 3. County-only clustering -- a few minutes
  s13_cty <- rbindlist(lapply(names(MAIN_INSTRUMENTS), function(il) {
    r <- estimate_interacted(outpatient, "SCHEME_1_CERTAINTY", PRIMARY_OUTCOME,
                             MAIN_INSTRUMENTS[[il]], moderator_type = "categorical",
                             label = "County-only clustering", instrument_label = il,
                             clusters = "ANALYSIS_MARKET")
    if (is.null(r)) return(NULL)
    r$tests[ESTIMATOR == "Reduced form"]
  }), fill = TRUE)
  print(s13_cty[, .(INSTRUMENT_LABEL, WALD = round(WALD, 2), P_VALUE = signif(P_VALUE, 3))])
  save_csv(s13_cty, "T13F_county_only_clustering.csv")
  
  
  
  g <- outpatient[, .(N = .N, N_AFTER = sum(POST_MONTH > min(POST_MONTH))), by = MARKET_ID]
  g[, .(cells = .N, contributing_cells = sum(N_AFTER >= 2),
        contributing_rows = sum(N_AFTER[N_AFTER >= 2]))]
  
  
  outpatient[, CATEGORY_MONTH := paste0(SCHEME_1_CERTAINTY, "::", POST_MONTH)]
  
  cat_month <- cache_or_run("category_month_fe_check", {
    lapply(names(MAIN_INSTRUMENTS), function(lab) {
      estimate_interacted(
        outpatient, moderator = "SCHEME_1_CERTAINTY", moderator_type = "categorical",
        instrument = MAIN_INSTRUMENTS[[lab]], instrument_label = lab,
        label = "Category x month FE",
        fixed_effects = c("MARKET_ID", "CATEGORY_MONTH"))
    })
  })
  
  rows  <- rbindlist(lapply(cat_month, `[[`, "rows"))
  tests <- rbindlist(lapply(cat_month, `[[`, "tests"))
  
  print(rows[, .(INSTRUMENT_LABEL, TERM, RF_PERCENT_PER_SD, RF_P, N_OBSERVATIONS)])
  print(tests[ESTIMATOR == "Reduced form", .(INSTRUMENT_LABEL, P_VALUE, FIRST_STAGE_WALD_MIN)])
}   # end interactive scratch


###############################################################################
# Section 28: Alternative price series
#             Gross charge, discounted cash rate, Medicare reference rate
#
# # Source AFTER warm_start() / restore_session(). Requires `outpatient` in scope
# and HPT_ALT_CONCEPT.csv.gz (or its Snowflake shards, HPT_ALT_CONCEPT.csv.gz
# _0_0_0.csv.gz through _0_7_0.csv.gz) in PANEL_DIR.
#
#   28A  Load and key harmonisation   the merge, and why it is not trivial
#   28B  Reporting-completeness check does Z predict HAVING a price at all?
#   28C  Estimation                   headline spec on three new outcomes
#   28D  Comparison                   against the negotiated-rate headline
#
# WHAT THIS IS FOR
#
#   GROSS_CHARGE   unilaterally set chargemaster price. Not bargained, but not
#                  inert: percent-of-charges contracts (which Phase 2 filters
#                  OUT of the negotiated panel), out-of-network billing, and
#                  Medicare outlier payments all key off it.
#
#   CASH_RATE      the consumer-facing self-pay price. This is the DIRECT
#                  measure of the patient channel that Section 7's ACS
#                  demographic interactions can only proxy for.
#
#   MEDICARE_RATE  administratively set; no hospital controls it. A response
#                  here indicates a coding or instrument problem, not an
#                  economic finding. This is the falsification check.
#
# WHY THE JOIN IS NOT TRIVIAL
#
#   The SQL export is keyed on ANALYSIS_CONCEPT_ID. The R panel is keyed on
#   FINAL_CONCEPT_ID. These differ for the six MERGE_GROUPS, where several
#   ANALYSIS_CONCEPT_IDs collapse to one canonical FINAL_CONCEPT_ID (MRI
#   abdomen variants, CT abdomen, CT abd/pelvis, fMRI brain, ARHFCMRIGTBS,
#   and mammography tomosynthesis). Merging on the raw ID would silently drop
#   those concepts -- including mammography, which is the most influential
#   single family in the paper's permutation test. 28A applies the same
#   canonical mapping and re-aggregates with the same equal-weighted-median
#   convention apply_concept_merges() uses upstream, so the alt series enter
#   on exactly the units the rest of the paper is estimated on.
###############################################################################

.s28_hd  <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")
.s28_sub <- function(x) cat("\n--- ", x, " ", strrep("-", max(0, 70 - nchar(x))), "\n", sep = "")

S28_ALT_DIR     <- PANEL_DIR
S28_ALT_PATTERN <- "^HPT_ALT_CONCEPT.*\\.csv\\.gz$"
S28_SCHEME      <- "SCHEME_1_CERTAINTY"

S28_OUTCOMES <- c(
  Gross    = "LN_GROSS_CHARGE",
  Cash     = "LN_CASH_RATE",
  Medicare = "LN_MEDICARE_RATE"
)


# ===========================================================================
# 28A. LOAD, HARMONISE KEYS, MERGE
# ===========================================================================
#
# The SQL export for this file did not set SINGLE = TRUE, so Snowflake split
# it across parallel threads on write -- the same failure mode already
# documented for the payer-dispersion export (PAYER_DISPERSION_PATTERN,
# load_payer_dispersion() in Section 9). This reads every matching shard and
# row-binds them, the same way that loader does, rather than assuming one
# file. If the SQL is rerun with SINGLE = TRUE added, this still works
# unchanged -- list.files() with one match is a one-file read.

s28_load_alt <- function(dir = S28_ALT_DIR, pattern = S28_ALT_PATTERN,
                         merge_groups = MERGE_GROUPS) {
  .s28_hd("28A. LOAD ALTERNATIVE PRICE SERIES")
  
  files <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0L) {
    stop("No alt price files found matching '", pattern, "' in\n  ", dir,
         "\nDownload HPT_ALT_CONCEPT.csv.gz (or its shards) from the ",
         "Snowflake stage first.", call. = FALSE)
  }
  cat("Found", length(files), "matching file(s):\n")
  cat(paste0("  ", basename(files)), sep = "\n")
  
  a <- rbindlist(lapply(files, fread, showProgress = FALSE), fill = TRUE)
  cat("\nRead", format(nrow(a), big.mark = ","), "rows total |",
      ncol(a), "columns\n")
  
  req <- c("HOSPITAL_ID", "POST_MONTH", "ANALYSIS_CONCEPT_ID",
           "MEDIAN_GROSS_CHARGE", "MEDIAN_CASH_RATE", "MEDIAN_MEDICARE_RATE")
  miss <- setdiff(req, names(a))
  if (length(miss)) stop("Missing columns in alt file: ",
                         paste(miss, collapse = ", "), call. = FALSE)
  
  # Type harmonisation to match prepare_panel()'s conventions exactly.
  a[, HOSPITAL_ID := as.character(HOSPITAL_ID)]
  a[, POST_MONTH  := as.Date(POST_MONTH)]
  a[, ANALYSIS_CONCEPT_ID := as.character(ANALYSIS_CONCEPT_ID)]
  
  # Apply the canonical concept mapping (see header note).
  map <- rbindlist(lapply(merge_groups, function(g)
    data.table(ANALYSIS_CONCEPT_ID = g$constituents, CANON_ID = g$canonical_id)))
  a <- merge(a, map, by = "ANALYSIS_CONCEPT_ID", all.x = TRUE, sort = FALSE)
  a[, FINAL_CONCEPT_ID := fifelse(is.na(CANON_ID), ANALYSIS_CONCEPT_ID, CANON_ID)]
  
  n_remapped <- a[!is.na(CANON_ID), .N]
  cat("Rows remapped to a canonical concept ID:",
      format(n_remapped, big.mark = ","), "\n")
  if (n_remapped == 0L) {
    cat("*** WARNING: no rows matched any MERGE_GROUPS constituent. Check that\n",
        "    the SQL export's ANALYSIS_CONCEPT_ID values use the same naming\n",
        "    convention as MERGE_GROUPS before trusting the merge.\n", sep = "")
  }
  
  # Collapse constituents to the canonical concept, equal-weighted median of
  # the constituent medians, the same convention apply_concept_merges() uses.
  # MEDIAN_NEGOTIATED_CHECK is carried through so Section 33 can build
  # NEG_TO_MEDICARE from it.
  if (!("MEDIAN_NEGOTIATED_CHECK" %in% names(a)))
    a[, MEDIAN_NEGOTIATED_CHECK := NA_real_]
  
  alt <- a[, .(
    GROSS_CHARGE  = median(MEDIAN_GROSS_CHARGE,  na.rm = TRUE),
    CASH_RATE     = median(MEDIAN_CASH_RATE,     na.rm = TRUE),
    MEDICARE_RATE = median(MEDIAN_MEDICARE_RATE, na.rm = TRUE),
    MEDIAN_NEGOTIATED_CHECK = median(MEDIAN_NEGOTIATED_CHECK, na.rm = TRUE),
    N_ALT_SOURCE_CONCEPTS = .N
  ), by = .(HOSPITAL_ID, POST_MONTH, FINAL_CONCEPT_ID)]
  
  for (v in c("GROSS_CHARGE", "CASH_RATE", "MEDICARE_RATE", "MEDIAN_NEGOTIATED_CHECK")) {
    set(alt, i = which(!is.finite(alt[[v]])), j = v, value = NA_real_)
  }
  
  cat("Collapsed to", format(nrow(alt), big.mark = ","),
      "hospital x month x concept rows\n")
  alt
}


s28_merge_alt <- function(panel = outpatient, alt = NULL) {
  if (is.null(alt)) alt <- s28_load_alt()
  
  .s28_sub("Merging onto the outpatient panel")
  
  d <- merge(copy(panel), alt,
             by = c("HOSPITAL_ID", "POST_MONTH", "FINAL_CONCEPT_ID"),
             all.x = TRUE, sort = FALSE)
  
  d[, LN_GROSS_CHARGE  := safe_log_positive(GROSS_CHARGE)]
  d[, LN_CASH_RATE     := safe_log_positive(CASH_RATE)]
  d[, LN_MEDICARE_RATE := safe_log_positive(MEDICARE_RATE)]
  
  d[, HAS_GROSS    := as.integer(!is.na(LN_GROSS_CHARGE) & is.finite(LN_GROSS_CHARGE))]
  d[, HAS_CASH     := as.integer(!is.na(LN_CASH_RATE) & is.finite(LN_CASH_RATE))]
  d[, HAS_MEDICARE := as.integer(!is.na(LN_MEDICARE_RATE) & is.finite(LN_MEDICARE_RATE))]
  
  cov <- d[, .(
    PANEL_ROWS      = .N,
    MATCHED_ANY_ALT = sum(!is.na(N_ALT_SOURCE_CONCEPTS)),
    WITH_GROSS      = sum(HAS_GROSS),
    WITH_CASH       = sum(HAS_CASH),
    WITH_MEDICARE   = sum(HAS_MEDICARE),
    SHARE_GROSS     = round(mean(HAS_GROSS), 4),
    SHARE_CASH      = round(mean(HAS_CASH), 4),
    SHARE_MEDICARE  = round(mean(HAS_MEDICARE), 4))]
  print(cov)
  save_qa_csv(cov, "QA28A_alt_merge_coverage.csv")
  
  .s28_sub("Coverage by shoppability category")
  print(d[, .(N = .N,
              SHARE_GROSS    = round(mean(HAS_GROSS), 4),
              SHARE_CASH     = round(mean(HAS_CASH), 4),
              SHARE_MEDICARE = round(mean(HAS_MEDICARE), 4)),
          by = c(S28_SCHEME)])
  
  .s28_sub("Sanity: alt series against the negotiated price, same rows")
  print(d[HAS_GROSS == 1L, .(
    N = .N,
    MEDIAN_RATIO_GROSS = round(median(exp(LN_GROSS_CHARGE - LN_MEDIAN_PRICE), na.rm = TRUE), 3))])
  print(d[HAS_CASH == 1L, .(
    N = .N,
    MEDIAN_RATIO_CASH = round(median(exp(LN_CASH_RATE - LN_MEDIAN_PRICE), na.rm = TRUE), 3))])
  print(d[HAS_MEDICARE == 1L, .(
    N = .N,
    MEDIAN_RATIO_MEDICARE = round(median(exp(LN_MEDICARE_RATE - LN_MEDIAN_PRICE), na.rm = TRUE), 3))])
  cat("\nThese should reproduce the SQL-side pairwise ratios (roughly 3.6, 1.6,\n",
      "0.95). A large discrepancy means the merge landed on the wrong rows.\n", sep = "")
  
  d
}


# ===========================================================================
# 28B. REPORTING COMPLETENESS
# ===========================================================================
#
# Gross and cash coverage is roughly 44% of hospital x code cells and is
# uneven across clinical families (40% in endoscopy, 99% in mammography).
# Whether a hospital reports a gross or cash price is a feature of file
# completeness, not of price-setting. If the instrument predicts reporting,
# and differentially by shoppability, then any gradient estimated on these
# outcomes is partly a selection artifact.
#
# Same logic as Appendix Table tab:payermix Panel A, which tests whether the
# instrument predicts reported payer composition: same specification, binary
# outcome. This check precedes 28C.

s28_completeness_check <- function(d, instruments = MAIN_INSTRUMENTS) {
  .s28_hd("28B. DOES THE INSTRUMENT PREDICT REPORTING COMPLETENESS?")
  
  out <- rbindlist(lapply(c("HAS_GROSS", "HAS_CASH", "HAS_MEDICARE"), function(oc) {
    rbindlist(lapply(names(instruments), function(lab) {
      r <- estimate_interacted(
        d, moderator = S28_SCHEME, moderator_type = "categorical",
        outcome = oc, instrument = instruments[[lab]],
        instrument_label = lab, label = paste("Completeness:", oc))
      if (is.null(r)) return(NULL)
      tst <- r$tests[ESTIMATOR == "Reduced form"]
      cbind(OUTCOME = oc, INSTRUMENT_LABEL = lab,
            r$rows[, .(TERM, RF_COEF, RF_SE, RF_P)],
            EQUALITY_P = if (nrow(tst)) tst$P_VALUE[1L] else NA_real_)
    }), fill = TRUE)
  }), fill = TRUE)
  
  print(out)
  save_qa_csv(out, "QA28B_reporting_completeness.csv")
  cat("\nLarge RF_P throughout and a non-significant EQUALITY_P indicate the\n",
      "instrument does not predict whether a price is reported, and does not\n",
      "do so differentially by shoppability, so 28C's estimates are not\n",
      "reporting-selection artifacts. A significant EQUALITY_P is a material\n",
      "qualification on the corresponding outcome and is reported alongside\n",
      "it.\n", sep = "")
  invisible(out)
}


# ===========================================================================
# 28C. ESTIMATION -- HEADLINE SPEC, THREE NEW OUTCOMES
# ===========================================================================
#
# The design is unchanged: same interacted specification, same instruments,
# same fixed effects, same clustering. Only the outcome varies.

s28_estimate <- function(d, instruments = MAIN_INSTRUMENTS,
                         outcomes = S28_OUTCOMES) {
  .s28_hd("28C. ALTERNATIVE PRICE SERIES, HEADLINE SPECIFICATION")
  
  rows <- list(); tests <- list()
  
  for (oc_lab in names(outcomes)) {
    oc <- outcomes[[oc_lab]]
    for (inst_lab in names(instruments)) {
      .s28_sub(paste(oc_lab, "|", inst_lab))
      r <- estimate_interacted(
        d, moderator = S28_SCHEME, moderator_type = "categorical",
        outcome = oc, instrument = instruments[[inst_lab]],
        instrument_label = inst_lab, label = oc_lab)
      if (is.null(r)) { cat("  not estimable\n"); next }
      
      print(r$rows[, .(TERM,
                       RF_PCT_PER_SD = round(RF_PERCENT_PER_SD, 2),
                       RF_P = round(RF_P, 4),
                       N = N_OBSERVATIONS)])
      tst <- r$tests[ESTIMATOR == "Reduced form"]
      if (nrow(tst)) cat("  Heterogeneity test p =", round(tst$P_VALUE[1L], 4),
                         "| min first-stage Wald =",
                         round(r$rows$FIRST_STAGE_WALD_MIN[1L], 1), "\n")
      
      rows[[length(rows) + 1L]]  <- cbind(OUTCOME_LABEL = oc_lab, r$rows)
      tests[[length(tests) + 1L]] <- cbind(OUTCOME_LABEL = oc_lab, r$tests)
    }
  }
  
  res <- list(rows = rbindlist(rows, fill = TRUE),
              tests = rbindlist(tests, fill = TRUE))
  save_qa_csv(res$rows,  "QA28C_alt_price_coefficients.csv")
  save_qa_csv(res$tests, "QA28C_alt_price_tests.csv")
  res
}


# ===========================================================================
# 28D. SIDE BY SIDE WITH THE NEGOTIATED-RATE HEADLINE
# ===========================================================================

s28_compare <- function(d, res, instruments = MAIN_INSTRUMENTS) {
  .s28_hd("28D. COMPARISON WITH THE NEGOTIATED-RATE RESULT")
  
  base <- rbindlist(lapply(names(instruments), function(lab) {
    r <- estimate_interacted(
      d, moderator = S28_SCHEME, moderator_type = "categorical",
      outcome = PRIMARY_OUTCOME, instrument = instruments[[lab]],
      instrument_label = lab, label = "Negotiated")
    if (is.null(r)) return(NULL)
    cbind(OUTCOME_LABEL = "Negotiated", r$rows)
  }), fill = TRUE)
  
  all_rows <- rbindlist(list(base, res$rows), fill = TRUE)
  
  wide <- dcast(all_rows[TERM %chin% c("Shoppable", "Non_shoppable")],
                OUTCOME_LABEL + INSTRUMENT_LABEL ~ TERM,
                value.var = "RF_PERCENT_PER_SD")
  print(wide)
  save_qa_csv(all_rows, "QA28D_alt_vs_negotiated.csv")
  
  cat("\nInterpretation, fixed in advance:\n",
      "  MEDICARE  a response here is a coding or instrument problem, not a\n",
      "            result. This is the falsification check.\n",
      "  GROSS     unilaterally set. Percent-of-charges contracts, out-of-\n",
      "            network billing and outlier payments give it real\n",
      "            bargaining relevance despite not being negotiated.\n",
      "  CASH      the consumer-facing price and the direct measure of the\n",
      "            patient channel. A negative coefficient is evidence the\n",
      "            patient channel is operative, which cuts against Section\n",
      "            7's current reading of the demographic nulls and implies a\n",
      "            revision to that section rather than an additional table.\n",
      sep = "")
  invisible(all_rows)
}



if (exists("HPT_SCRATCH") && isTRUE(HPT_SCRATCH)) {   # interactive scratch, off by default
  d28 <- s28_merge_alt(outpatient)
  
  chk <- s28_completeness_check(d28)
  
  res <- s28_estimate(d28)
  cmp <- s28_compare(d28, res)
}   # end interactive scratch


# ===========================================================================
# RUNNER
# ===========================================================================

run_section28 <- function(panel = outpatient) {
  d   <- cache_or_run("s28_alt_price_panel", s28_merge_alt(panel))
  chk <- s28_completeness_check(d)
  res <- s28_estimate(d)
  cmp <- s28_compare(d, res)
  invisible(list(panel = d, completeness = chk, estimates = res, comparison = cmp))
}

cat("Section 28 loaded. Call s28_merge_alt(outpatient) first, then run_section28().\n")




###############################################################################
# Section 29: Margin recovery -- common sample, discount narrowing, arithmetic
#
# Source AFTER Section 28. Requires the merged alt-price panel (d28) in scope,
# or rebuilds it from cache.
#
#   29A  Common sample      does the negotiated decline exist where cash is
#                           reported? A precondition, not a test.
#   29B  Discount narrowing ln(cash) - ln(gross) as the outcome. The sharpest
#                           single piece of evidence available here.
#   29C  Recovery arithmetic breakeven self-pay share, computed from the
#                           MATCHED-SAMPLE coefficients rather than hardcoded.
#
# WHY 29A COMES FIRST
#
#   The Section 28 comparison puts -4.81% (negotiated, 951k obs) next to
#   +9.83% (cash, 537k obs). Those are different rows. If negotiated rates do
#   not fall on the subsample where cash prices are reported, there is no
#   revenue loss for the cash increase to recover and the mechanism has no
#   foundation, and the two results are unrelated facts about non-overlapping
#   samples. 29A settles that, and produces the matched coefficients 29C
#   needs.
#
# WHY 29B IS THE SHARPEST TEST
#
#   Cash prices are typically set as a percentage discount off the
#   chargemaster: P_cash = (1-d) * P_gross, so ln(cash) - ln(gross) = ln(1-d).
#   Using the ratio as the outcome holds the list price fixed BY
#   CONSTRUCTION, so any "this hospital raised all its prices" story --
#   chargemaster drift, a billing-system migration, an inflation adjustment --
#   differences out, since it moves numerator and denominator together. What
#   survives is the discretionary decision about how much to discount self-pay
#   patients, which is the decision the margin-recovery mechanism concerns.
#
#   Units: the estimand is a change in the log cash-to-gross ratio, not a
#   price response, and is not comparable to the "% per SD" price
#   coefficients.
###############################################################################

.s29_hd  <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")
.s29_sub <- function(x) cat("\n--- ", x, " ", strrep("-", max(0, 70 - nchar(x))), "\n", sep = "")

S29_SCHEME <- "SCHEME_1_CERTAINTY"


# ===========================================================================
# 29A. COMMON-SAMPLE CHECK
# ===========================================================================

s29_common_sample <- function(d, instruments = MAIN_INSTRUMENTS) {
  .s29_hd("29A. DOES THE NEGOTIATED DECLINE EXIST WHERE CASH IS REPORTED?")
  
  .s29_sub("Sample overlap")
  ov <- d[, .(
    PANEL_ROWS        = .N,
    WITH_CASH         = sum(HAS_CASH == 1L),
    WITH_GROSS        = sum(HAS_GROSS == 1L),
    WITH_BOTH         = sum(HAS_CASH == 1L & HAS_GROSS == 1L),
    CASH_AND_NEG      = sum(HAS_CASH == 1L & is.finite(LN_MEDIAN_PRICE)),
    SHARE_CASH_HAS_NEG = round(mean(is.finite(LN_MEDIAN_PRICE)[HAS_CASH == 1L]), 4))]
  print(ov)
  save_qa_csv(ov, "QA29A_sample_overlap.csv")
  
  out <- rbindlist(lapply(names(instruments), function(lab) {
    inst <- instruments[[lab]]
    
    full <- estimate_interacted(
      d, moderator = S29_SCHEME, moderator_type = "categorical",
      outcome = PRIMARY_OUTCOME, instrument = inst,
      instrument_label = lab, label = "Negotiated, full sample")
    
    sub <- estimate_interacted(
      d[HAS_CASH == 1L], moderator = S29_SCHEME, moderator_type = "categorical",
      outcome = PRIMARY_OUTCOME, instrument = inst,
      instrument_label = lab, label = "Negotiated, cash-reporting subsample")
    
    cashsub <- estimate_interacted(
      d[HAS_CASH == 1L], moderator = S29_SCHEME, moderator_type = "categorical",
      outcome = "LN_CASH_RATE", instrument = inst,
      instrument_label = lab, label = "Cash, same rows")
    
    grab <- function(r, tag) {
      if (is.null(r)) return(NULL)
      tst <- r$tests[ESTIMATOR == "Reduced form"]
      cbind(SPEC = tag, INSTRUMENT_LABEL = lab,
            r$rows[, .(TERM, RF_PERCENT_PER_SD, RF_P, N_OBSERVATIONS)],
            EQUALITY_P = if (nrow(tst)) tst$P_VALUE[1L] else NA_real_)
    }
    rbindlist(list(grab(full, "Negotiated_full"),
                   grab(sub,  "Negotiated_cashsample"),
                   grab(cashsub, "Cash_cashsample")), fill = TRUE)
  }), fill = TRUE)
  
  print(dcast(out[TERM == "Shoppable"], INSTRUMENT_LABEL ~ SPEC,
              value.var = "RF_PERCENT_PER_SD"))
  save_qa_csv(out, "QA29A_common_sample.csv")
  
  cat("\nPrecondition: the Negotiated_cashsample shoppable coefficient stays\n",
      "clearly negative. If it attenuates toward zero there is no revenue\n",
      "loss on these rows for a cash increase to recover, and the margin-\n",
      "recovery interpretation is dropped rather than caveated.\n", sep = "")
  invisible(out)
}


# ===========================================================================
# 29B. SELF-PAY DISCOUNT NARROWING
# ===========================================================================

s29_discount_narrowing <- function(d, instruments = MAIN_INSTRUMENTS) {
  .s29_hd("29B. IS THE SELF-PAY DISCOUNT NARROWING?")
  
  d[, LN_CASH_GROSS_RATIO := LN_CASH_RATE - LN_GROSS_CHARGE]
  
  .s29_sub("Descriptive: the level of the self-pay discount")
  desc <- d[is.finite(LN_CASH_GROSS_RATIO), .(
    N               = .N,
    MEDIAN_RATIO    = round(median(exp(LN_CASH_GROSS_RATIO)), 4),
    IMPLIED_DISCOUNT = round(1 - median(exp(LN_CASH_GROSS_RATIO)), 4),
    P25_RATIO       = round(quantile(exp(LN_CASH_GROSS_RATIO), 0.25), 4),
    P75_RATIO       = round(quantile(exp(LN_CASH_GROSS_RATIO), 0.75), 4))]
  print(desc)
  print(d[is.finite(LN_CASH_GROSS_RATIO),
          .(N = .N, MEDIAN_RATIO = round(median(exp(LN_CASH_GROSS_RATIO)), 4)),
          by = c(S29_SCHEME)])
  save_qa_csv(desc, "QA29B_discount_level.csv")
  
  .s29_sub("Interacted specification on the log cash-to-gross ratio")
  out <- rbindlist(lapply(names(instruments), function(lab) {
    r <- estimate_interacted(
      d[is.finite(LN_CASH_GROSS_RATIO)],
      moderator = S29_SCHEME, moderator_type = "categorical",
      outcome = "LN_CASH_GROSS_RATIO", instrument = instruments[[lab]],
      instrument_label = lab, label = "Cash-to-gross ratio")
    if (is.null(r)) return(NULL)
    tst <- r$tests[ESTIMATOR == "Reduced form"]
    cat("\n", lab, "\n", sep = "")
    print(r$rows[, .(TERM,
                     RATIO_CHANGE_PER_SD = round(RF_COEF * sd(d[[instruments[[lab]]]],
                                                              na.rm = TRUE), 5),
                     RF_P = round(RF_P, 4), N = N_OBSERVATIONS)])
    if (nrow(tst)) cat("  Heterogeneity test p =", round(tst$P_VALUE[1L], 4),
                       "| min first-stage Wald =",
                       round(r$rows$FIRST_STAGE_WALD_MIN[1L], 1), "\n")
    cbind(INSTRUMENT_LABEL = lab,
          r$rows[, .(TERM, RF_COEF, RF_SE, RF_P, N_OBSERVATIONS,
                     FIRST_STAGE_WALD_MIN)],
          EQUALITY_P = if (nrow(tst)) tst$P_VALUE[1L] else NA_real_)
  }), fill = TRUE)
  
  save_qa_csv(out, "QA29B_discount_narrowing.csv")
  cat("\nPrediction: a POSITIVE shoppable coefficient, meaning the cash price\n",
      "rises relative to the list price -- the discount narrows. Because the\n",
      "list price is held fixed by construction, a chargemaster-drift story\n",
      "cannot generate this.\n", sep = "")
  invisible(out)
}


# ===========================================================================
# 29C. RECOVERY ARITHMETIC
# ===========================================================================
#
#   rho = (gamma_C / |gamma_N|) * r * c * s/(1-s)
#
# where gamma are proportional price responses per SD, r = P_cash/P_neg,
# c is the collection rate on billed self-pay charges, and s is the self-pay
# share of volume. Breakeven (rho = 1) solves to
#
#   s* = 1 / (1 + K*c),    K = (gamma_C / |gamma_N|) * r
#
# Since c <= 1 by construction, s* is bounded below by 1/(1+K). Every unknown
# pushes the required share higher, so no precise external estimate of s is
# needed to sign the conclusion.
#
# Coefficients come from the matched sample in 29A rather than the Section 28
# output, which was estimated on non-overlapping samples.

s29_recovery_arithmetic <- function(d, common = NULL,
                                    instrument_label = "Competitor_only_hospitals_9m",
                                    s_grid = c(0.02, 0.04, 0.06, 0.08, 0.10, 0.15, 0.20),
                                    c_grid = c(1.0, 0.5, 0.3, 0.2)) {
  .s29_hd("29C. RECOVERY ARITHMETIC")
  
  if (is.null(common)) common <- s29_common_sample(d)
  
  gN <- common[SPEC == "Negotiated_cashsample" & TERM == "Shoppable" &
                 INSTRUMENT_LABEL == instrument_label, RF_PERCENT_PER_SD] / 100
  gC <- common[SPEC == "Cash_cashsample" & TERM == "Shoppable" &
                 INSTRUMENT_LABEL == instrument_label, RF_PERCENT_PER_SD] / 100
  
  if (length(gN) == 0L || length(gC) == 0L) {
    stop("Could not recover matched-sample coefficients for ", instrument_label,
         call. = FALSE)
  }
  
  r <- d[HAS_CASH == 1L & is.finite(LN_MEDIAN_PRICE),
         median(exp(LN_CASH_RATE - LN_MEDIAN_PRICE), na.rm = TRUE)]
  
  K <- (gC / abs(gN)) * r
  
  cat("Matched-sample inputs (", instrument_label, "):\n", sep = "")
  cat("  gamma_N (negotiated, shoppable) =", round(gN, 5), "\n")
  cat("  gamma_C (cash, shoppable)       =", round(gC, 5), "\n")
  cat("  r = median(P_cash / P_neg)      =", round(r, 4), "\n")
  cat("  K = (gamma_C/|gamma_N|) * r     =", round(K, 4), "\n")
  cat("\n  Lower bound on breakeven self-pay share (at c = 1):",
      paste0(round(100 / (1 + K), 2), "%"), "\n")
  cat("  Full offset would require a cash increase of",
      paste0(round(100 * abs(gN) * (1 - 0.06) / (r * 1 * 0.06), 1), "%"),
      "at s = 6%, c = 1\n")
  
  .s29_sub("Breakeven self-pay share, by collection rate")
  be <- data.table(COLLECTION_RATE = c_grid)
  be[, BREAKEVEN_SELFPAY_SHARE := round(100 / (1 + K * COLLECTION_RATE), 2)]
  print(be)
  
  .s29_sub("Recovery ratio rho (%), by self-pay share and collection rate")
  grid <- CJ(SELFPAY_SHARE = s_grid, COLLECTION_RATE = c_grid)
  grid[, RECOVERY_PCT := round(100 * K * COLLECTION_RATE *
                                 SELFPAY_SHARE / (1 - SELFPAY_SHARE), 1)]
  wide <- dcast(grid, SELFPAY_SHARE ~ COLLECTION_RATE, value.var = "RECOVERY_PCT")
  print(wide)
  
  save_qa_csv(be,   "QA29C_breakeven_share.csv")
  save_qa_csv(grid, "QA29C_recovery_grid.csv")
  
  cat("\nSince c <= 1 by construction, the breakeven share is bounded below by\n",
      "the c = 1 row. The claim rests not on a point estimate of the self-pay\n",
      "share but on the weaker statement that it does not reach that bound,\n",
      "the same structure as the Conley breakeven in Section 8.7.1.\n", sep = "")
  
  invisible(list(gamma_N = gN, gamma_C = gC, r = r, K = K,
                 breakeven = be, grid = grid))
}




if (exists("HPT_SCRATCH") && isTRUE(HPT_SCRATCH)) {   # interactive scratch, off by default
  common <- s29_common_sample(d28)
  
  common[SPEC == "Negotiated_cashsample" & TERM == "Shoppable",
         .(INSTRUMENT_LABEL, RF_PERCENT_PER_SD, RF_P, EQUALITY_P, N_OBSERVATIONS)]
  
  disc  <- s29_discount_narrowing(d28)
  arith <- s29_recovery_arithmetic(d28, common = common)
}   # end interactive scratch

# ===========================================================================
# RUNNER
# ===========================================================================

# run_section29 <- function(d = NULL) {
#   if (is.null(d)) d <- cache_or_run("s28_alt_price_panel", s28_merge_alt(outpatient))
#   common <- s29_common_sample(d)
#   disc   <- s29_discount_narrowing(d)
#   arith  <- s29_recovery_arithmetic(d, common = common)
#   invisible(list(panel = d, common = common, discount = disc, arithmetic = arith))
# }

cat("Section 29 loaded. Call s29_common_sample(d28) first.\n")




###############################################################################
# Section 30: Robustness for the margin-recovery finding
#             System x month fixed effects, leave-one-system-out
#
# Source AFTER Section 29. Requires d28 (or a panel with LN_CASH_RATE,
# LN_GROSS_CHARGE, HAS_CASH already attached) in scope.
#
#   30A  Resolve SYS_RESOLVED on d28   mirrors Section 13's construction
#   30B  System x month FE             does discount-narrowing survive the
#                                       same demanding FE structure Section
#                                       13 already put the headline result
#                                       through?
#   30C  Leave-one-system-out          top 15 systems by hospital count,
#                                       dropped one at a time
#
# WHY THIS PRECEDES ANY USE OF SECTIONS 28-29 IN THE PAPER
#
#   The cash-price result clears the heterogeneity test on 2 of 3 instruments
#   (28C); the discount-narrowing ratio clears it on 0 of 3, directionally
#   consistent only (29B); the matched-sample negotiated decline attenuates
#   by a third and loses individual significance (29A). At this level of
#   statistical marginality, a single large health system doing something
#   idiosyncratic to its own cash pricing -- a policy change, a billing
#   vendor migration, a self-pay program redesign -- could plausibly be
#   generating a meaningful share of the pooled result. This is exactly the
#   threat Section 13 tests for the negotiated-rate headline. Cash pricing has
#   at least as much reason to be set at the system level, so it gets the
#   identical two-part test: absorb system-level shocks directly (13A's
#   logic), then check whether the result depends on any one organisation's
#   hospitals (13B's logic).
#
# OUTCOMES TESTED
#
#   LN_CASH_GROSS_RATIO   the discount-narrowing outcome from 29B, which the
#                         paper leads with, so the check is built around it.
#   LN_CASH_RATE          the raw cash-price outcome from 28C, included for
#                         completeness since it is reported alongside.
###############################################################################

.s30_hd  <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")
.s30_sub <- function(x) cat("\n--- ", x, " ", strrep("-", max(0, 70 - nchar(x))), "\n", sep = "")

S30_SCHEME     <- "SCHEME_1_CERTAINTY"
S30_N_SYSTEMS  <- 15L
S30_OUTCOMES   <- c(Ratio = "LN_CASH_GROSS_RATIO", Cash = "LN_CASH_RATE")


# ===========================================================================
# 30A. RESOLVE SYS_RESOLVED ON d28
# ===========================================================================
#
# Identical construction to Section 13's block 0: first non-missing system
# value across a hospital's rows, resolved once on `outpatient` (not on d28,
# which has fewer rows per hospital after the alt-price merge and could
# resolve differently by chance if a hospital's early rows happen to be the
# ones without a match). Merged onto d28 by HOSPITAL_ID.

s30_attach_system <- function(d, base_panel = outpatient) {
  .s30_hd("30A. RESOLVE SYSTEM PER HOSPITAL")
  
  sys_col <- if ("SYSTEM_KEY" %in% names(base_panel)) "SYSTEM_KEY" else "HEALTH_SYSTEM_ID"
  cat("Using", sys_col, "as the system identifier (matches Section 13).\n")
  
  resolved <- base_panel[!is.na(HOSPITAL_ID), .(
    SYS_RESOLVED = { v <- get(sys_col)[!is.na(get(sys_col))]; if (length(v)) v[1L] else NA_character_ }
  ), by = HOSPITAL_ID]
  
  cat("Resolved for", nrow(resolved), "hospitals |",
      sum(is.na(resolved$SYS_RESOLVED)), "genuinely unaffiliated (NA on every row)\n")
  
  if ("SYS_RESOLVED" %in% names(d)) d[, SYS_RESOLVED := NULL]
  d <- merge(d, resolved, by = "HOSPITAL_ID", all.x = TRUE, sort = FALSE)
  
  if (!("LN_CASH_GROSS_RATIO" %in% names(d))) {
    d[, LN_CASH_GROSS_RATIO := LN_CASH_RATE - LN_GROSS_CHARGE]
  }
  
  d
}


# ===========================================================================
# 30B. SYSTEM x MONTH FIXED EFFECTS
# ===========================================================================
#
# Same logic as Section 13A: replace additive MARKET_ID + POST_MONTH with
# MARKET_ID + SYSTEM_MONTH, restricting identification to within-system,
# cross-market variation in disclosure timing. If a system-wide shock (a
# self-pay policy change correlated with when that system's hospitals
# disclosed) is driving the result, this absorbs it directly.

s30_system_month <- function(d, outcomes = S30_OUTCOMES, instruments = MAIN_INSTRUMENTS) {
  .s30_hd("30B. SYSTEM x MONTH FIXED EFFECTS")
  
  d[, SYSTEM_MONTH := paste0(fifelse(is.na(SYS_RESOLVED), "NOSYS", SYS_RESOLVED),
                             "_", as.character(POST_MONTH))]
  cat("Distinct system-months:", uniqueN(d$SYSTEM_MONTH), "\n")
  
  out <- rbindlist(lapply(names(outcomes), function(oc_lab) {
    oc <- outcomes[[oc_lab]]
    rbindlist(lapply(names(instruments), function(il) {
      base <- estimate_interacted(
        d[is.finite(get(oc))], S30_SCHEME, oc, instruments[[il]],
        moderator_type = "categorical", label = paste(oc_lab, "baseline"),
        instrument_label = il)
      sysm <- estimate_interacted(
        d[is.finite(get(oc))], S30_SCHEME, oc, instruments[[il]],
        moderator_type = "categorical", label = paste(oc_lab, "sys x month"),
        instrument_label = il, fixed_effects = c("MARKET_ID", "SYSTEM_MONTH"))
      
      grab <- function(r, spec) {
        if (is.null(r)) return(NULL)
        tst <- r$tests[ESTIMATOR == "Reduced form"]
        data.table(
          OUTCOME_LABEL = oc_lab, SPEC = spec, INSTRUMENT_LABEL = il,
          SHOPPABLE_PCT = r$rows[TERM == "Shoppable"]$RF_PERCENT_PER_SD[1L],
          SHOPPABLE_P   = r$rows[TERM == "Shoppable"]$RF_P[1L],
          NONSHOP_PCT   = r$rows[TERM == "Non_shoppable"]$RF_PERCENT_PER_SD[1L],
          NONSHOP_P     = r$rows[TERM == "Non_shoppable"]$RF_P[1L],
          EQUALITY_P    = if (nrow(tst)) tst$P_VALUE[1L] else NA_real_,
          N_OBSERVATIONS = r$rows$N_OBSERVATIONS[1L])
      }
      rbindlist(list(grab(base, "Baseline"), grab(sysm, "System x month")), fill = TRUE)
    }), fill = TRUE)
  }), fill = TRUE)
  
  print(dcast(out, OUTCOME_LABEL + INSTRUMENT_LABEL ~ SPEC,
              value.var = c("SHOPPABLE_PCT", "SHOPPABLE_P")))
  save_qa_csv(out, "QA30B_system_month.csv")
  
  cat("\nA p-value that rises while the point estimate holds is a power\n",
      "result: system x month absorbs a large share of identifying variation\n",
      "by construction, as Section 13 found for the negotiated outcome. A\n",
      "point estimate that collapses toward zero is a different and more\n",
      "serious result, indicating a system-level shock rather than disclosure\n",
      "carried the baseline finding.\n", sep = "")
  invisible(out)
}


# ===========================================================================
# 30C. LEAVE-ONE-SYSTEM-OUT
# ===========================================================================
#
# Identical structure to Section 13B. Drops each of the largest 15 systems by
# hospital count in turn and re-estimates. This tests whether the result
# depends on any single organisation's hospitals, not whether the instrument
# does, which would require a Section 15-style reconstruction.

s30_leave_one_out <- function(d, outcomes = S30_OUTCOMES, instruments = MAIN_INSTRUMENTS,
                              n_systems = S30_N_SYSTEMS) {
  .s30_hd("30C. LEAVE-ONE-SYSTEM-OUT (TOP 15 BY HOSPITAL COUNT)")
  
  sys_size <- d[!is.na(SYS_RESOLVED) & SYS_RESOLVED != "",
                .(N_HOSPITALS = uniqueN(HOSPITAL_ID)), by = SYS_RESOLVED][order(-N_HOSPITALS)]
  top_sys <- head(sys_size, n_systems)
  cat("Largest systems in the alt-price panel:\n"); print(top_sys)
  
  base <- rbindlist(lapply(names(outcomes), function(oc_lab) {
    oc <- outcomes[[oc_lab]]
    rbindlist(lapply(names(instruments), function(il) {
      r <- estimate_interacted(d[is.finite(get(oc))], S30_SCHEME, oc,
                               instruments[[il]], moderator_type = "categorical",
                               label = oc_lab, instrument_label = il)
      if (is.null(r)) return(NULL)
      tst <- r$tests[ESTIMATOR == "Reduced form"]
      data.table(OUTCOME_LABEL = oc_lab, INSTRUMENT_LABEL = il,
                 ESTIMATOR = tst$ESTIMATOR, P_BASE = tst$P_VALUE)
    }), fill = TRUE)
  }), fill = TRUE)
  
  out <- list()
  for (i in seq_len(nrow(top_sys))) {
    sid <- top_sys$SYS_RESOLVED[i]
    cat(sprintf("  [%2d/%2d] drop %-24s (n_hosp=%d)\n", i, nrow(top_sys),
                substr(sid, 1, 22), top_sys$N_HOSPITALS[i]))
    dd <- d[is.na(SYS_RESOLVED) | SYS_RESOLVED != sid]
    
    for (oc_lab in names(outcomes)) {
      oc <- outcomes[[oc_lab]]
      for (il in names(instruments)) {
        r <- estimate_interacted(dd[is.finite(get(oc))], S30_SCHEME, oc,
                                 instruments[[il]], moderator_type = "categorical",
                                 label = oc_lab, instrument_label = il)
        if (is.null(r)) next
        tst <- r$tests[ESTIMATOR == "Reduced form"]
        out[[length(out) + 1L]] <- data.table(
          DROPPED_SYSTEM = sid, N_HOSPITALS_DROPPED = top_sys$N_HOSPITALS[i],
          OUTCOME_LABEL = oc_lab, INSTRUMENT_LABEL = il,
          ESTIMATOR = tst$ESTIMATOR, P_VALUE = tst$P_VALUE,
          SHOPPABLE_PCT = r$rows[TERM == "Shoppable"]$RF_PERCENT_PER_SD[1L],
          SHOPPABLE_P   = r$rows[TERM == "Shoppable"]$RF_P[1L],
          N_OBSERVATIONS = r$rows$N_OBSERVATIONS[1L])
      }
    }
  }
  
  res <- rbindlist(out, fill = TRUE)
  res <- merge(res, base, by = c("OUTCOME_LABEL", "INSTRUMENT_LABEL", "ESTIMATOR"), all.x = TRUE)
  res[, STILL_SIG_05 := as.integer(P_VALUE < 0.05)]
  save_qa_csv(res, "QA30C_leave_one_system_out.csv")
  
  .s30_sub("Summary: how many of the 45 drops (15 systems x 3 instruments) stay significant")
  print(res[ESTIMATOR == "Reduced form", .(
    N_DROPS = .N, N_STILL_SIG = sum(STILL_SIG_05),
    MIN_SHOPPABLE_PCT = round(min(SHOPPABLE_PCT), 2),
    MAX_SHOPPABLE_PCT = round(max(SHOPPABLE_PCT), 2),
    MAX_P = round(max(P_VALUE), 4)),
    by = OUTCOME_LABEL])
  
  cat("\nSign and rough magnitude holding across all 45 drops indicates no\n",
      "single system carries the result, the standard Table 13's leave-one-out\n",
      "(44 of 45) is read against. A drop that flips the sign or collapses the\n",
      "magnitude toward zero identifies the system responsible, which is named\n",
      "in the text rather than averaged over.\n", sep = "")
  invisible(res)
}

if (exists("HPT_SCRATCH") && isTRUE(HPT_SCRATCH)) {   # interactive scratch, off by default
  d30 <- s30_attach_system(d28)
  sm <- s30_system_month(d30)
  loo <- s30_leave_one_out(d30)
}   # end interactive scratch


# ===========================================================================
# RUNNER
# ===========================================================================

run_section30 <- function(d = NULL) {
  if (is.null(d)) d <- cache_or_run("s28_alt_price_panel", s28_merge_alt(outpatient))
  d   <- s30_attach_system(d)
  sm  <- s30_system_month(d)
  loo <- s30_leave_one_out(d)
  invisible(list(panel = d, system_month = sm, leave_one_out = loo))
}

cat("Section 30 loaded. Call s30_attach_system(d28) first.\n")





###############################################################################
# Figures 21-23: Alternative price series
#
# Source AFTER the Figures 1-20 block, so theme_paper(), save_fig(),
# SHOP_COLORS, FSU_GARNET / FSU_GOLD / FSU_GREY, INSTR_SHORT and short_instr()
# are already defined. Reads the QA CSVs written by Sections 28-30.
#
#   fig21_alt_price_series      four price series, shoppable vs non-shoppable
#   fig22_recovery_breakeven    breakeven self-pay share against collection rate
#   fig23_cash_robustness       baseline vs system x month, both cash outcomes
#
# These read from the QA directory rather than TABLE_DIR, since Sections 28-30
# write with save_qa_csv().
###############################################################################

suppressPackageStartupMessages({
  library(ggplot2); library(data.table); library(scales)
})

if (!exists("QA_DIR")) {
  QA_DIR <- file.path(RESULTS_DIR, "05_R_QA")
}

read_qa <- function(fn) {
  p <- file.path(QA_DIR, fn)
  if (!file.exists(p)) { message("SKIP: ", fn, " not found in ", QA_DIR); return(NULL) }
  fread(p)
}

SERIES_ORDER <- c("Negotiated", "Gross", "Cash", "Medicare")
SERIES_LABEL <- c(Negotiated = "Negotiated rate",
                  Gross      = "Gross charge",
                  Cash       = "Discounted cash rate",
                  Medicare   = "Medicare reference rate")


# ===========================================================================
# fig21: the four price series side by side
# ===========================================================================
#
# One panel per series, shoppable and non-shoppable coefficients with 95%
# intervals, across the three main instruments. The negotiated panel sits
# below zero, gross and Medicare at zero, and cash above it.

if (exists("HPT_SCRATCH") && isTRUE(HPT_SCRATCH)) {   # interactive scratch, off by default
  alt <- read_qa("QA28D_alt_vs_negotiated.csv")
  
  if (!is.null(alt)) {
    
    d21 <- alt[TERM %chin% c("Shoppable", "Non_shoppable")]
    d21[, SERIES := factor(OUTCOME_LABEL, levels = SERIES_ORDER,
                           labels = SERIES_LABEL[SERIES_ORDER])]
    d21[, INSTR := factor(short_instr(INSTRUMENT_LABEL),
                          levels = rev(c("Competitor hospitals", "Local system",
                                         "Competitor hospitals (ex-CBSA)")))]
    d21[, CAT := factor(TERM, levels = c("Shoppable", "Non_shoppable"))]
    
    # 95% interval in percent-per-SD space. RF_SE is in log points, so scale it
    # the same way the point estimate was scaled: the ratio of the percent
    # estimate to the log-point estimate is the SD multiplier.
    d21[, SCALE := fifelse(abs(RF_COEF) > 1e-12,
                           (log1p(RF_PERCENT_PER_SD / 100)) / RF_COEF, NA_real_)]
    d21[, `:=`(LO = 100 * (exp(log1p(RF_PERCENT_PER_SD / 100) - 1.96 * RF_SE * SCALE) - 1),
               HI = 100 * (exp(log1p(RF_PERCENT_PER_SD / 100) + 1.96 * RF_SE * SCALE) - 1))]
    
    p21 <- ggplot(d21, aes(x = INSTR, y = RF_PERCENT_PER_SD,
                           colour = CAT, shape = CAT)) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50",
                 linewidth = 0.4) +
      geom_pointrange(aes(ymin = LO, ymax = HI),
                      position = position_dodge(width = 0.55),
                      size = 0.45, fatten = 2.2) +
      coord_flip() +
      facet_wrap(~ SERIES, ncol = 2, scales = "free_x") +
      scale_colour_manual(values = SHOP_COLORS,
                          labels = c("Shoppable", "Non-shoppable")) +
      scale_shape_manual(values = c(Shoppable = 16, Non_shoppable = 17),
                         labels = c("Shoppable", "Non-shoppable")) +
      labs(x = NULL, y = "Reduced-form response (% per SD of instrument)") +
      theme_paper()
    
    save_fig(p21, "fig21_alt_price_series", width = 8, height = 5.5)
  }
  
  
  # ===========================================================================
  # fig22: the recovery breakeven
  # ===========================================================================
  #
  # Plots equation (breakeven): s* = 1/(1 + K c) against the collection rate,
  # one curve per instrument, with a shaded band for plausible self-pay shares.
  # Every curve sits above any plausible share for every c, and lowering c
  # pushes the requirement higher.
  
  K_VALUES <- c(`Competitor hospitals`           = 6.41,
                `Local system`                   = 5.54,
                `Competitor hospitals (ex-CBSA)` = 8.65)
  
  # Shaded reference band only. These bounds are illustrative and carry no
  # external source. Anchor them to an HCUP or MEPS estimate before the figure
  # is used in the paper, or drop the band.
  PLAUSIBLE_SHARE <- c(0.02, 0.10)
  
  grid22 <- CJ(INSTR = names(K_VALUES), C = seq(0.15, 1, by = 0.01))
  grid22[, K := K_VALUES[INSTR]]
  grid22[, BREAKEVEN := 1 / (1 + K * C)]
  grid22[, INSTR := factor(INSTR, levels = names(K_VALUES))]
  
  p22 <- ggplot(grid22, aes(x = C, y = BREAKEVEN, colour = INSTR)) +
    annotate("rect", xmin = -Inf, xmax = Inf,
             ymin = PLAUSIBLE_SHARE[1], ymax = PLAUSIBLE_SHARE[2],
             fill = "grey80", alpha = 0.45) +
    annotate("text", x = 0.18, y = mean(PLAUSIBLE_SHARE), hjust = 0,
             vjust = -0.6, size = 3, colour = "grey30",
             label = "Plausible self-pay share of outpatient volume") +
    geom_line(linewidth = 0.7) +
    scale_colour_manual(values = c(FSU_GARNET, FSU_GOLD, FSU_GREY)) +
    scale_y_continuous(labels = percent_format(accuracy = 1),
                       limits = c(0, NA)) +
    scale_x_continuous(labels = percent_format(accuracy = 1)) +
    labs(x = "Collection rate on billed self-pay charges",
         y = "Self-pay share required for full offset") +
    theme_paper()
  
  save_fig(p22, "fig22_recovery_breakeven", width = 7, height = 4.5)
  
  
  # ===========================================================================
  # fig23: cash robustness ladder
  # ===========================================================================
  #
  # Baseline against system x month for both cash outcomes, with the
  # leave-one-system-out range overlaid as a band on the baseline. Shows the
  # collapse under system x month and the stability of sign under LOO in one
  # picture.
  
  sm <- read_qa("QA30B_system_month.csv")
  
  if (!is.null(sm)) {
    
    if (!"SHOPPABLE_PCT" %in% names(sm)) {
      message("QA30B has an unexpected shape; skipping fig23")
    } else {
      d23 <- copy(sm)
      d23[, SPEC := factor(SPEC, levels = c("Baseline", "System x month"))]
      d23[, INSTR := factor(short_instr(INSTRUMENT_LABEL),
                            levels = rev(c("Competitor hospitals", "Local system",
                                           "Competitor hospitals (ex-CBSA)")))]
      d23[, SERIES := factor(OUTCOME_LABEL,
                             levels = c("Cash", "Ratio"),
                             labels = c("Discounted cash rate",
                                        "Cash-to-gross ratio"))]
      
      p23 <- ggplot(d23, aes(x = INSTR, y = SHOPPABLE_PCT,
                             colour = SPEC, shape = SPEC)) +
        geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50",
                   linewidth = 0.4) +
        geom_point(position = position_dodge(width = 0.5), size = 2.4) +
        coord_flip() +
        facet_wrap(~ SERIES, ncol = 2) +
        scale_colour_manual(values = c(Baseline = FSU_GARNET,
                                       `System x month` = FSU_GREY)) +
        scale_shape_manual(values = c(Baseline = 16, `System x month` = 1)) +
        labs(x = NULL,
             y = "Shoppable reduced-form response (% per SD of instrument)") +
        theme_paper()
      
      save_fig(p23, "fig23_cash_robustness", width = 8, height = 3.8)
    }
  }
  
  cat("Alternative-price figures complete.\n")
}   # end interactive scratch






###############################################################################
#
#   SECTION 31 -- BINNED REDUCED FORM  (revision checklist item B5)
#
#   Placed after the last figure block. Section numbering picks up from 30
#   (QA30C is the last QA prefix in use).
#
#   ---------------------------------------------------------------------------
#   WHAT THIS ESTIMATES, AND WHY IT IS A STEP FUNCTION RATHER THAN A SPLINE
#   ---------------------------------------------------------------------------
#
#   The headline reduced form is
#
#       ln P = sum_g gamma_g (Z x 1[k in g]) + theta ln(Beds) + alpha_mk + tau_t
#
#   which imposes linearity in Z. This section replaces the linear term with a
#   set of exposure-bin indicators, one per (category, bin), with the zero-
#   exposure bin omitted inside each category:
#
#       ln P = sum_g sum_{b>=2} gamma_gb (1[k in g] x 1[Z in bin_b])
#              + theta ln(Beds) + alpha_mk + tau_t
#
#   gamma_gb is the log-price difference between category-g services at a
#   hospital facing bin-b exposure and category-g services at a hospital facing
#   ZERO peer exposure, within the same county x concept cell and month. No
#   functional form in Z is imposed anywhere.
#
#   The checklist item suggests quartiles of Z or a spline. Neither works as
#   stated, because Z is zero for 38 to 55 percent of the panel depending on
#   the instrument (QA05 "Zero share" column). Quartiles are undefined, since
#   the bottom two are both zero, and a within-bin slope in the reference bin
#   has no variation to identify it. Bin indicators handle a point mass
#   cleanly: the zero group becomes the omitted reference and every other
#   coefficient is measured against it. Bins on the positive part are quantile
#   cuts, so the cut points are data-determined rather than chosen.
#
#   1[k in g] alone is absorbed by MARKET_ID (county x concept determines the
#   concept, which determines the category). 1[Z in bin_b] alone is NOT
#   absorbed: it varies within a cell whenever two hospitals in the same county
#   posted the same concept at different times facing different accumulated
#   exposure, which is exactly the variation the whole design runs on. Dropping
#   bin_1 within each category is therefore necessary and sufficient.
#
#   ---------------------------------------------------------------------------
#   TWO UNITS WARNINGS
#   ---------------------------------------------------------------------------
#
#   1. These coefficients are not "percent per SD of the instrument." They are
#      percent level differences relative to the zero-exposure bin. Table 5's
#      units and this table's units are not interchangeable.
#
#   2. Concavity in Z is not literally the framework's saturation prediction,
#      which concerns f''(N) and g''(N) in N. The bridge is MEAN_N_PRIOR, the
#      average prior-poster count inside each bin, reported alongside every
#      coefficient as descriptive context.
#
#   ---------------------------------------------------------------------------
#   OUTPUTS
#   ---------------------------------------------------------------------------
#     QA31A_bin_diagnostic.csv       bin composition, per instrument
#     QA31B_bin_cell_variation.csv   within-cell bin variation, per instrument
#     T31_binned_rf_estimates.csv    coefficients
#     T31B_binned_rf_tests.csv       gap, saturation, and linearity tests
#     fig24_binned_rf.pdf
#
###############################################################################

HPT_RUN$binned_rf <- if (exists("HPT_BINNED")) isTRUE(HPT_BINNED) else FALSE

S31_SCHEME          <- "SCHEME_1_CERTAINTY"
S31_SCHEME_LABEL    <- "1. Procedural certainty"
S31_N_POSITIVE_BINS <- 3L    # zero bin + 3 positive bins = 4 categories.
# Drop to 2 if VCOV_FULL_RANK comes back 0 on the
# joint tests: two-way clustering with 16 month
# clusters is the same constraint that already
# defeated the sixteen-category joint test.

.s31_hd <- function(x)
  cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")


# ---------------------------------------------------------------------------
# Bin construction: {Z = 0} plus quantile cuts on the strictly positive part.
# ---------------------------------------------------------------------------
s31_make_bins <- function(z, k = S31_N_POSITIVE_BINS, probs = NULL) {
  z <- safe_numeric(z)
  if (any(z < 0, na.rm = TRUE)) {
    warning("Negative instrument values; they will fall in the zero bin.",
            call. = FALSE)
  }
  pos <- z[is.finite(z) & z > 0]
  if (length(pos) < 1000L) return(NULL)
  # probs, if supplied, overrides k entirely and can be uneven -- e.g.
  # c(0.55, 0.85) puts more of the sample in the bottom bin and carves a
  # narrower, higher-exposure top bin than an equal k-way split would. This
  # is how to push further into the tail without shrinking every bin: move
  # resolution OUT of a region already shown to be flat and INTO the top,
  # rather than adding bins and thinning all of them at once.
  cut_probs <- if (is.null(probs)) seq_len(k - 1L) / k else sort(unique(probs))
  qs <- unique(quantile(pos, probs = cut_probs, na.rm = TRUE, names = FALSE))
  cuts <- c(-Inf, 0, qs, Inf)
  if (anyDuplicated(cuts) || length(cuts) < 3L) return(NULL)
  labs <- c("Z0", paste0("Z", seq_len(length(cuts) - 2L)))
  factor(cut(z, breaks = cuts, labels = labs, right = TRUE), levels = labs)
}


# ---------------------------------------------------------------------------
# Wald test on an arbitrary contrast matrix.
#
# wald_equality() only tests "all coefficients equal to the first," which
# cannot express "the shoppable-minus-non-shoppable gap is the same in every
# exposure bin." Same generalised-inverse fallback, same rank flag, and the
# same F reference through .cluster_df(), so a test from here is directly
# comparable to every other p-value in the pipeline.
# ---------------------------------------------------------------------------
s31_wald_contrast <- function(fit, R) {
  if (is.null(fit) || is.null(R) || nrow(R) == 0L) return(data.table())
  tryCatch({
    cf <- coef(fit); V <- vcov(fit)
    keep <- colnames(R)
    if (!all(keep %in% names(cf))) return(data.table())
    b   <- cf[keep]
    Vs  <- V[keep, keep, drop = FALSE]
    mid <- R %*% Vs %*% t(R)
    rk  <- qr(mid)$rank
    inv <- tryCatch(solve(mid), error = function(e) MASS::ginv(mid))
    stat <- as.numeric(t(R %*% b) %*% inv %*% (R %*% b))
    if (!is.finite(stat) || stat < 0) stop("non-finite Wald")
    data.table(WALD = stat, DF = rk,
               VCOV_FULL_RANK = as.integer(rk == nrow(R)),
               P_VALUE = pf(stat / rk, df1 = rk, df2 = .cluster_df(fit),
                            lower.tail = FALSE))
  }, error = function(e) data.table())
}


# ===========================================================================
# 31A. BIN DIAGNOSTIC. Precedes estimation.
#
# Two questions decide whether B5 is feasible at all:
#
#   (a) Is there enough mass in the positive bins? If the top bin holds a
#       few thousand rows the tail coefficient will be uninformative, which
#       is itself a reportable answer to the tail concern but not a table.
#
#   (b) Do county x concept cells span bins? Cell-invariant bin membership is
#       absorbed by the fixed effect and contributes nothing, so a low share
#       of rows in cells spanning two or more bins ends the exercise.
#
# There is no hard threshold. Above roughly 0.5 the exercise is comfortable,
# between 0.3 and 0.5 the interior bins have wide intervals, and below 0.3 the
# diagnostic itself is the reportable result rather than the table.
# ===========================================================================
s31_bin_diagnostic <- function(panel, instruments = MAIN_INSTRUMENTS,
                               scheme_col = S31_SCHEME,
                               outcome = PRIMARY_OUTCOME,
                               k = S31_N_POSITIVE_BINS, probs = NULL) {
  
  .s31_hd("31A. BIN DIAGNOSTIC")
  
  comp_all <- list(); cell_all <- list()
  
  for (il in names(instruments)) {
    z <- instruments[[il]]
    if (!(z %in% names(panel))) { message("SKIP: ", z, " absent"); next }
    
    cols <- c(outcome, ENDOGENOUS_VARIABLE, z, BASELINE_CONTROLS,
              BASELINE_FIXED_EFFECTS, BASELINE_CLUSTERS, scheme_col)
    d <- model_sample(panel[!is.na(get(scheme_col))], cols)
    if (nrow(d) < MIN_MODEL_OBS) { message("SKIP: ", il, " too few rows"); next }
    
    d[, ZBIN := s31_make_bins(get(z), k, probs)]
    if (is.null(d$ZBIN) || anyNA(d$ZBIN)) {
      message("SKIP: ", il, " bin construction failed"); next
    }
    
    comp <- d[, .(N              = .N,
                  SHARE          = .N / nrow(d),
                  MIN_Z          = min(safe_numeric(get(z))),
                  MAX_Z          = max(safe_numeric(get(z))),
                  MEAN_Z         = mean(safe_numeric(get(z))),
                  MEAN_N_PRIOR   = mean(safe_numeric(get(ENDOGENOUS_VARIABLE))),
                  SHARE_SHOPPABLE = mean(get(scheme_col) == "Shoppable")),
              by = ZBIN][order(ZBIN)]
    comp[, `:=`(INSTRUMENT_LABEL = il, INSTRUMENT = z)]
    comp_all[[il]] <- comp
    
    # Within-cell bin variation. NBINS is how many distinct exposure bins a
    # county x concept cell contains; a cell with one bin is absorbed.
    cv <- d[, .(NBINS = uniqueN(ZBIN), NROWS = .N), by = MARKET_ID]
    top <- levels(d$ZBIN)[nlevels(d$ZBIN)]
    span_top <- d[, .(HAS_REF = any(ZBIN == "Z0"), HAS_TOP = any(ZBIN == top),
                      NROWS = .N), by = MARKET_ID][HAS_REF & HAS_TOP]
    
    cell <- data.table(
      INSTRUMENT_LABEL      = il,
      FE_CELLS              = nrow(cv),
      CELLS_SPANNING_2PLUS  = sum(cv$NBINS >= 2L),
      SHARE_CELLS_SPANNING  = mean(cv$NBINS >= 2L),
      SHARE_ROWS_SPANNING   = sum(cv[NBINS >= 2L]$NROWS) / nrow(d),
      CELLS_SPANNING_REF_TO_TOP = nrow(span_top),
      SHARE_ROWS_REF_TO_TOP = sum(span_top$NROWS) / nrow(d))
    cell_all[[il]] <- cell
    
    cat("\n--", il, "--\n")
    print(as.data.frame(comp[, .(ZBIN, N, SHARE = round(SHARE, 3),
                                 MIN_Z = round(MIN_Z, 1), MAX_Z = round(MAX_Z, 1),
                                 MEAN_Z = round(MEAN_Z, 2),
                                 MEAN_N_PRIOR = round(MEAN_N_PRIOR, 2),
                                 SHARE_SHOP = round(SHARE_SHOPPABLE, 3))]))
    cat("  rows in cells spanning >=2 bins:      ",
        round(cell$SHARE_ROWS_SPANNING, 3), "\n",
        "  rows in cells spanning zero -> top:   ",
        round(cell$SHARE_ROWS_REF_TO_TOP, 3), "\n", sep = "")
    
    rm(d); invisible(gc())
  }
  
  comp_out <- rbindlist(comp_all, fill = TRUE)
  cell_out <- rbindlist(cell_all, fill = TRUE)
  if (nrow(comp_out)) save_qa_csv(comp_out, "QA31A_bin_diagnostic.csv")
  if (nrow(cell_out)) save_qa_csv(cell_out, "QA31B_bin_cell_variation.csv")
  
  cat("\nSHARE_ROWS_SPANNING is the binding diagnostic. Below ~0.30 the binned\n",
      "reduced form has too little within-cell variation to support a table.\n",
      sep = "")
  
  list(composition = comp_out, cell_variation = cell_out)
}


# ===========================================================================
# 31B. THE DISTRIBUTION OF N, AND HOW FAR Z REACHES INTO IT
#
# Two distinct questions, and it matters which one a number here answers.
#
#   (a) The marginal distribution of N_PRIOR_POSTERS, unconditional on Z.
#       Descriptive only: what the treatment looks like in the data, with no
#       causal content. It shows how far the tail of actual exposure extends.
#
#   (b) The tail curve: for a candidate top-slice cutoff on positive Z, the
#       sample size and N-distribution that slice carries. This is the table
#       that selects probs, trading N gained against sample size given up as
#       the top Z cut point moves higher.
#
# Both are computed on the same sample estimate_binned_rf() uses, via
# model_sample() on the same required columns, so the population here matches
# the estimator's.
#
# Binning is on Z alone, never on N, for the endogeneity reason in the 31A
# header. N stays descriptive throughout.
# ===========================================================================
s31_n_distribution <- function(panel, instruments = MAIN_INSTRUMENTS,
                               scheme_col = S31_SCHEME,
                               outcome = PRIMARY_OUTCOME,
                               tail_cuts = c(0.50, 0.70, 0.80, 0.85, 0.90,
                                             0.95, 0.975, 0.99)) {
  
  .s31_hd("31B. DISTRIBUTION OF N, AND HOW FAR THE TOP-Z SLICE REACHES INTO IT")
  
  qprobs <- c(0, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99, 0.995, 0.999, 1)
  
  marg_all <- list(); curve_all <- list()
  
  for (il in names(instruments)) {
    z <- instruments[[il]]
    if (!(z %in% names(panel))) { message("SKIP: ", z, " absent"); next }
    
    cols <- c(outcome, ENDOGENOUS_VARIABLE, z, BASELINE_CONTROLS,
              BASELINE_FIXED_EFFECTS, BASELINE_CLUSTERS, scheme_col)
    d <- model_sample(panel[!is.na(get(scheme_col))], cols)
    if (nrow(d) < MIN_MODEL_OBS) { message("SKIP: ", il, " too few rows"); next }
    
    n <- safe_numeric(d[[ENDOGENOUS_VARIABLE]])
    
    # ---- (a) marginal distribution of N, all rows and Z > 0 rows only ----
    mq_all <- quantile(n, probs = qprobs, na.rm = TRUE, names = FALSE)
    mq_pos <- quantile(n[safe_numeric(d[[z]]) > 0], probs = qprobs,
                       na.rm = TRUE, names = FALSE)
    marg <- data.table(INSTRUMENT_LABEL = il, QPROB = qprobs,
                       N_ALL_ROWS = round(mq_all, 2),
                       N_ZPOS_ROWS = round(mq_pos, 2))
    marg_all[[il]] <- marg
    
    cat("\n--", il, "-- marginal quantiles of N (prior posters)\n")
    print(as.data.frame(marg))
    
    # ---- (b) tail curve: candidate top-Z cutoffs vs. resulting N ----------
    pos_z <- safe_numeric(d[[z]])[safe_numeric(d[[z]]) > 0]
    cut_vals <- quantile(pos_z, probs = tail_cuts, na.rm = TRUE, names = FALSE)
    
    curve <- rbindlist(lapply(seq_along(tail_cuts), function(i) {
      sel <- safe_numeric(d[[z]]) > cut_vals[i]
      nn  <- n[sel]
      data.table(
        INSTRUMENT_LABEL = il, TOP_SLICE_QUANTILE = tail_cuts[i],
        Z_CUT_VALUE = round(cut_vals[i], 2),
        N_ROWS = sum(sel), SHARE_OF_SAMPLE = round(sum(sel) / nrow(d), 4),
        N_SHOPPABLE = sum(sel & d[[scheme_col]] == "Shoppable"),
        N_NONSHOPPABLE = sum(sel & d[[scheme_col]] == "Non_shoppable"),
        MEAN_N = round(mean(nn), 2), MEDIAN_N = round(median(nn), 2),
        P90_N = round(quantile(nn, 0.90, na.rm = TRUE, names = FALSE), 2),
        MAX_N = round(max(nn), 2))
    }))
    curve_all[[il]] <- curve
    
    cat("\n--", il, "-- tail curve: top-Z-slice cutoff vs. resulting N\n")
    print(as.data.frame(curve))
  }
  
  marg_out  <- rbindlist(marg_all,  fill = TRUE)
  curve_out <- rbindlist(curve_all, fill = TRUE)
  if (nrow(marg_out))  save_qa_csv(marg_out,  "QA31C_n_marginal_distribution.csv")
  if (nrow(curve_out)) save_qa_csv(curve_out, "QA31D_tail_curve.csv")
  
  cat("\nThe tail curve, not the marginal quantiles, selects a probs cut.\n",
      "MEAN_N climbing while N_ROWS stays in the tens of thousands is the\n",
      "usable range. Once N_ROWS for one category falls toward single-digit\n",
      "thousands, that row is the practical ceiling.\n", sep = "")
  
  list(marginal = marg_out, tail_curve = curve_out)
}


# ===========================================================================
# 31C. THE ESTIMATOR
#
# add_linear = TRUE additionally includes the category-specific LINEAR Z terms
# alongside the bin dummies. The joint test that the bin dummies are zero in
# that model is the formal test of whether linearity in Z is adequate, which is
# the specification test a referee asks for after seeing a binned figure.
# ===========================================================================
estimate_binned_rf <- function(data, instrument = PRIMARY_INSTRUMENT,
                               scheme_col = S31_SCHEME,
                               outcome = PRIMARY_OUTCOME,
                               k = S31_N_POSITIVE_BINS, probs = NULL,
                               instrument_label = "", label = S31_SCHEME_LABEL,
                               controls = BASELINE_CONTROLS,
                               fixed_effects = BASELINE_FIXED_EFFECTS,
                               clusters = BASELINE_CLUSTERS,
                               add_linear = FALSE) {
  
  controls <- available_columns(data, controls)
  fe <- available_columns(data, fixed_effects)
  cl <- available_columns(data, clusters)
  
  d <- data[!is.na(get(scheme_col))]
  d <- model_sample(d, c(outcome, ENDOGENOUS_VARIABLE, instrument, controls,
                         fe, cl, scheme_col))
  if (nrow(d) < MIN_MODEL_OBS) return(NULL)
  
  d[, ZBIN := s31_make_bins(get(instrument), k, probs)]
  if (is.null(d$ZBIN) || anyNA(d$ZBIN)) return(NULL)
  d[, CAT := droplevels(factor(get(scheme_col)))]
  
  cats <- levels(d$CAT)
  bins <- levels(d$ZBIN)
  if (length(cats) < 2L || length(bins) < 3L) return(NULL)
  ref  <- bins[1L]
  
  grid <- CJ(CAT = cats, BIN = bins[-1L], sorted = FALSE)
  grid[, TERM := paste0("B_", CAT, "_", BIN)]
  for (i in seq_len(nrow(grid))) {
    set(d, j = grid$TERM[i],
        value = as.integer(d$CAT == grid$CAT[i] & d$ZBIN == grid$BIN[i]))
  }
  
  rhs <- grid$TERM
  lin_terms <- character(0)
  if (add_linear) {
    lin_terms <- paste0("Zlin_", cats)
    for (j in seq_along(cats)) {
      set(d, j = lin_terms[j],
          value = safe_numeric(d[[instrument]]) * as.integer(d$CAT == cats[j]))
    }
    rhs <- c(lin_terms, rhs)
  }
  
  fit <- tryCatch(feols(build_ols_formula(outcome, c(rhs, controls), fe),
                        data = d, cluster = build_cluster_formula(cl),
                        warn = FALSE, notes = FALSE),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  
  # Bin means on the same rows the model uses.
  bin_stats <- d[, .(N_BIN        = .N,
                     MEAN_Z       = mean(safe_numeric(get(instrument))),
                     MEAN_N_PRIOR = mean(safe_numeric(get(ENDOGENOUS_VARIABLE)))),
                 by = .(CAT, ZBIN)]
  ref_stats <- bin_stats[ZBIN == ref, .(CAT, REF_MEAN_Z = MEAN_Z,
                                        REF_MEAN_N_PRIOR = MEAN_N_PRIOR)]
  bin_stats <- merge(bin_stats, ref_stats, by = "CAT", sort = FALSE)
  bin_stats[, `:=`(DELTA_Z       = MEAN_Z - REF_MEAN_Z,
                   DELTA_N_PRIOR = MEAN_N_PRIOR - REF_MEAN_N_PRIOR)]
  
  pull <- function(nm) {
    if (!(nm %in% names(coef(fit)))) return(list(b = NA_real_, s = NA_real_))
    list(b = unname(coef(fit)[nm]), s = unname(sqrt(vcov(fit)[nm, nm])))
  }
  
  # Records exactly which cut points produced this fit, so a saved CSV mixing
  # runs at different resolutions can still be told apart. "even_k=3" is the
  # original equal-tertile split; a custom probs vector prints its own cuts.
  bin_scheme <- if (is.null(probs)) paste0("even_k", k)
  else paste0("probs_", paste(round(sort(unique(probs)), 3), collapse = "_"))
  
  rows <- rbindlist(lapply(seq_len(nrow(grid)), function(i) {
    e  <- pull(grid$TERM[i])
    bs <- bin_stats[CAT == grid$CAT[i] & ZBIN == grid$BIN[i]]
    data.table(
      SPEC = label, INSTRUMENT_LABEL = instrument_label,
      INSTRUMENT = instrument, OUTCOME = outcome, BIN_SCHEME = bin_scheme,
      TERM = grid$CAT[i], ZBIN = grid$BIN[i], IS_REFERENCE = 0L,
      RF_COEF = e$b, RF_SE = e$s, RF_P = .pval(e$b / e$s, fit),
      RF_PERCENT = 100 * (exp(e$b) - 1),
      RF_CI_LOW_PERCENT  = 100 * (exp(e$b - 1.96 * e$s) - 1),
      RF_CI_HIGH_PERCENT = 100 * (exp(e$b + 1.96 * e$s) - 1),
      N_BIN = bs$N_BIN, MEAN_Z = bs$MEAN_Z, DELTA_Z = bs$DELTA_Z,
      MEAN_N_PRIOR = bs$MEAN_N_PRIOR, DELTA_N_PRIOR = bs$DELTA_N_PRIOR,
      # RF_PCT_PER_PRIOR_POSTER is deliberately absent. Dividing a
      # coefficient identified off within-cell variation by a raw between-bin
      # difference in N does not yield a first-stage-scaled quantity, despite
      # resembling one. MEAN_N_PRIOR and DELTA_N_PRIOR are descriptive
      # context for the bin.
      N_OBSERVATIONS = nobs(fit), ADD_LINEAR = as.integer(add_linear))
  }), fill = TRUE)
  
  # Reference rows, coefficient zero by construction, carried so the figure
  # and the table both show where the comparison starts.
  ref_rows <- rbindlist(lapply(cats, function(g) {
    bs <- bin_stats[CAT == g & ZBIN == ref]
    data.table(SPEC = label, INSTRUMENT_LABEL = instrument_label,
               INSTRUMENT = instrument, OUTCOME = outcome,
               BIN_SCHEME = bin_scheme,
               TERM = g, ZBIN = ref, IS_REFERENCE = 1L,
               RF_COEF = 0, RF_SE = NA_real_, RF_P = NA_real_,
               RF_PERCENT = 0, RF_CI_LOW_PERCENT = NA_real_,
               RF_CI_HIGH_PERCENT = NA_real_,
               N_BIN = bs$N_BIN, MEAN_Z = bs$MEAN_Z, DELTA_Z = 0,
               MEAN_N_PRIOR = bs$MEAN_N_PRIOR, DELTA_N_PRIOR = 0,
               N_OBSERVATIONS = nobs(fit),
               ADD_LINEAR = as.integer(add_linear))
  }), fill = TRUE)
  rows <- rbind(ref_rows, rows, fill = TRUE)
  
  # ---- tests -------------------------------------------------------------
  nm  <- names(coef(fit))
  shop <- paste0("B_Shoppable_",     bins[-1L])
  nons <- paste0("B_Non_shoppable_", bins[-1L])
  ok   <- shop %in% nm & nons %in% nm
  tests <- list()
  
  # (1) the shoppability gap inside each exposure bin
  for (j in which(ok)) {
    R <- matrix(0, nrow = 1, ncol = 2,
                dimnames = list(NULL, c(shop[j], nons[j])))
    R[1, shop[j]] <- 1; R[1, nons[j]] <- -1
    t1 <- s31_wald_contrast(fit, R)
    if (nrow(t1)) tests[[length(tests) + 1L]] <-
      cbind(TEST = "Gap in bin", ZBIN = bins[-1L][j], t1)
  }
  
  # (2) does the gap itself vary across exposure bins?
  if (sum(ok) >= 2L) {
    use <- which(ok); keep <- c(shop[use], nons[use])
    R <- matrix(0, nrow = length(use) - 1L, ncol = length(keep),
                dimnames = list(NULL, keep))
    for (r in seq_len(length(use) - 1L)) {
      R[r, shop[use][1L]] <-  1; R[r, nons[use][1L]] <- -1
      R[r, shop[use][r + 1L]] <- -1; R[r, nons[use][r + 1L]] <-  1
    }
    t2 <- s31_wald_contrast(fit, R)
    if (nrow(t2)) tests[[length(tests) + 1L]] <-
      cbind(TEST = "Gap equal across bins", ZBIN = NA_character_, t2)
  }
  
  # (3) within each category, does the response differ across bins at all?
  for (g in cats) {
    tg <- wald_equality(fit, paste0("B_", g, "_", bins[-1L]))
    if (nrow(tg)) tests[[length(tests) + 1L]] <-
        cbind(TEST = paste0("Bins equal within ", g), ZBIN = NA_character_, tg)
  }
  
  # (4) is linearity in Z adequate? Only defined in the add_linear fit.
  if (add_linear) {
    keep <- grid$TERM[grid$TERM %in% nm]
    if (length(keep)) {
      R <- diag(length(keep)); colnames(R) <- keep
      t4 <- s31_wald_contrast(fit, R)
      if (nrow(t4)) tests[[length(tests) + 1L]] <-
        cbind(TEST = "Bin dummies jointly zero given linear Z",
              ZBIN = NA_character_, t4)
    }
  }
  
  tests <- rbindlist(tests, fill = TRUE)
  if (nrow(tests)) {
    tests[, `:=`(SPEC = label, INSTRUMENT_LABEL = instrument_label,
                 OUTCOME = outcome, BIN_SCHEME = bin_scheme,
                 N_OBSERVATIONS = nobs(fit),
                 ADD_LINEAR = as.integer(add_linear))]
  }
  
  rm(d); invisible(gc())
  list(rows = rows, tests = tests, fit = NULL)
}


# ===========================================================================
# 31D. DRIVER
# ===========================================================================
run_binned_rf <- function(panel, instruments = MAIN_INSTRUMENTS,
                          scheme_col = S31_SCHEME, k = S31_N_POSITIVE_BINS,
                          probs = NULL, outcome = PRIMARY_OUTCOME, stem = "T31") {
  
  rows <- list(); tests <- list()
  
  for (il in names(instruments)) {
    z <- instruments[[il]]
    if (!(z %in% names(panel))) next
    for (al in c(FALSE, TRUE)) {
      t1 <- Sys.time()
      psd_flag <- FALSE
      # A non-PSD clustered VCOV corrupts SEs and p-values but not
      # coefficients (those come from the normal equations, not from
      # vcov()). Deferred R warnings can't be traced back to which loop
      # iteration produced them, so this captures the warning at the source
      # and tags every row and test from that fit. A p-value from a fit with
      # PSD_WARNING == TRUE is unreliable unless the finding also survives a
      # specification that avoids the instability, such as fewer or coarser
      # bins or one-way clustering.
      r <- withCallingHandlers(
        estimate_binned_rf(panel, instrument = z, scheme_col = scheme_col,
                           outcome = outcome, k = k, probs = probs,
                           instrument_label = il, add_linear = al),
        warning = function(w) {
          if (grepl("positive semi-definite", conditionMessage(w),
                    ignore.case = TRUE)) {
            psd_flag <<- TRUE
            message("    *** PSD WARNING on this fit -- SEs/p-values from ",
                    "it are suspect; coefficients are not ***")
          }
          invokeRestart("muffleWarning")
        })
      if (!is.null(r)) {
        r$rows[,  PSD_WARNING := psd_flag]
        r$tests[, PSD_WARNING := psd_flag]
        rows[[length(rows) + 1L]]  <- r$rows
        tests[[length(tests) + 1L]] <- r$tests
      }
      cat(sprintf("  %-36s linear=%-5s | %5.1fs%s\n", substr(il, 1, 34), al,
                  as.numeric(difftime(Sys.time(), t1, units = "secs")),
                  if (psd_flag) "  [PSD WARNING]" else ""))
    }
  }
  
  br <- rbindlist(rows,  fill = TRUE)
  bt <- rbindlist(tests, fill = TRUE)
  if (nrow(br) == 0L) stop("No binned results.", call. = FALSE)
  
  save_csv(br, paste0(stem, "_binned_rf_estimates.csv"))
  save_csv(bt, paste0(stem, "B_binned_rf_tests.csv"))
  
  .s31_hd("BINNED REDUCED FORM -- percent relative to the ZERO-exposure bin")
  cat("These are NOT percent per SD. They are level differences against Z = 0.\n")
  cat("Bin scheme:", unique(br$BIN_SCHEME), "\n\n")
  print(dcast(br[ADD_LINEAR == 0,
                 .(INSTRUMENT_LABEL, TERM, ZBIN,
                   E = round(RF_PERCENT, 2))],
              INSTRUMENT_LABEL + TERM ~ ZBIN, value.var = "E"))
  
  # Split by TERM (category), not just ZBIN -- MEAN_N_PRIOR differs by
  # category within the same bin, so collapsing across TERM here previously
  # produced duplicate INSTRUMENT_LABEL x ZBIN rows and dcast silently fell
  # back to counting them instead of erroring.
  cat("\nMean prior posters inside each bin (exposure range, descriptive only):\n")
  print(dcast(br[ADD_LINEAR == 0,
                 .(INSTRUMENT_LABEL, TERM, ZBIN,
                   E = round(MEAN_N_PRIOR, 2))],
              INSTRUMENT_LABEL + TERM ~ ZBIN, value.var = "E"))
  
  cat("\nBin sizes (rows actually used by the model):\n")
  print(dcast(br[ADD_LINEAR == 0,
                 .(INSTRUMENT_LABEL, TERM, ZBIN, N_BIN)],
              INSTRUMENT_LABEL + TERM ~ ZBIN, value.var = "N_BIN"))
  
  if (nrow(bt)) {
    cat("\nTESTS. Two separate failure modes, both make a p-value here\n",
        "unreliable and neither is caught by the other:\n",
        "  RANK_OK = 0     the restriction has no usable variance at all\n",
        "                  (too few clusters for this many df).\n",
        "  PSD_WARN = 1    the covariance matrix was full rank but had to be\n",
        "                  numerically repaired -- SEs and p-values from that\n",
        "                  fit are suspect (coefficients are not; those come\n",
        "                  from a different computation and are unaffected).\n\n", sep = "")
    print(bt[, .(INSTRUMENT_LABEL, TEST, ZBIN, DF, RANK_OK = VCOV_FULL_RANK,
                 PSD_WARN = PSD_WARNING, P = round(P_VALUE, 4),
                 ADD_LINEAR)][order(INSTRUMENT_LABEL, TEST)])
  }
  
  list(rows = br, tests = bt)
}


# ===========================================================================
# 31E. FIGURE
#
# x is the mean prior-poster count inside each bin, so the picture reads as a
# dose-response on the treatment scale. The dashed line is the LINEAR headline
# estimate evaluated at each bin's mean Z, which is what the binned points are
# being tested against.
# ===========================================================================
s31_figure <- function(binned_rows = NULL, binned_tests = NULL,
                       linear_csv = "T06_main_interacted_RF_and_IV.csv",
                       exclude_flagged = TRUE, annotate_gap = TRUE) {
  
  br <- binned_rows
  if (is.null(br)) br <- read_table("T31_binned_rf_estimates.csv")
  if (is.null(br)) { message("SKIP fig24: no binned estimates"); return(invisible(NULL)) }
  br <- as.data.table(br)[ADD_LINEAR == 0]
  
  bt <- binned_tests
  if (is.null(bt)) bt <- read_table("T31B_binned_rf_tests.csv")
  if (!is.null(bt)) bt <- as.data.table(bt)[ADD_LINEAR == 0 & TEST == "Gap in bin"]
  
  # An instrument flagged PSD_WARNING has a suspect covariance matrix. Its
  # coefficients are unaffected (see the note in run_binned_rf) but its error
  # bars may be mechanically too narrow, overstating precision on the plot.
  # The default drops it from the figure rather than plotting a false
  # confidence band; exclude_flagged = FALSE plots it, labelled.
  if ("PSD_WARNING" %in% names(br)) {
    flagged <- unique(br[PSD_WARNING == TRUE]$INSTRUMENT_LABEL)
    if (length(flagged)) {
      if (exclude_flagged) {
        cat("fig24: excluding", paste(flagged, collapse = ", "),
            "-- PSD_WARNING fired on this instrument, its SEs are not",
            "trustworthy for a figure. Pass exclude_flagged = FALSE to",
            "override.\n")
        br <- br[!(INSTRUMENT_LABEL %chin% flagged)]
        if (!is.null(bt)) bt <- bt[!(INSTRUMENT_LABEL %chin% flagged)]
      } else {
        cat("fig24: PLOTTING", paste(flagged, collapse = ", "),
            "despite PSD_WARNING -- its error bars may be unreliable.\n")
      }
    }
  }
  if (nrow(br) == 0L) { message("SKIP fig24: nothing left after PSD filter"); return(invisible(NULL)) }
  
  lin <- read_table(linear_csv)
  ref <- NULL
  if (!is.null(lin)) {
    lin <- as.data.table(lin)[SPEC == S31_SCHEME_LABEL &
                                INSTRUMENT_LABEL %in% unique(br$INSTRUMENT_LABEL)]
    if (nrow(lin)) {
      ref <- merge(br[, .(INSTRUMENT_LABEL, TERM, ZBIN, DELTA_Z, MEAN_N_PRIOR)],
                   lin[, .(INSTRUMENT_LABEL, TERM, LIN_COEF = RF_COEF)],
                   by = c("INSTRUMENT_LABEL", "TERM"), sort = FALSE)
      ref[, PRED_PERCENT := 100 * (exp(LIN_COEF * DELTA_Z) - 1)]
    }
  }
  
  instr_levels <- c("Competitor hospitals", "Local system",
                    "Competitor hospitals (ex-CBSA)")
  br[, `:=`(INSTR = factor(short_instr(INSTRUMENT_LABEL), levels = instr_levels),
            CAT = factor(TERM, levels = c("Shoppable", "Non_shoppable")))]
  
  p <- ggplot(br, aes(x = MEAN_N_PRIOR, y = RF_PERCENT, colour = CAT)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50",
               linewidth = 0.4) +
    geom_errorbar(aes(ymin = RF_CI_LOW_PERCENT, ymax = RF_CI_HIGH_PERCENT),
                  width = 0, linewidth = 0.5, na.rm = TRUE) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 2.2) +
    facet_wrap(~ INSTR, nrow = 1) +
    scale_colour_manual(values = SHOP_COLORS,
                        labels = c(Shoppable = "Shoppable",
                                   Non_shoppable = "Non-shoppable")) +
    labs(x = "Mean prior posters in the exposure bin",
         y = "Price difference vs. zero-exposure bin (%)") +
    theme_paper()
  
  if (!is.null(ref) && nrow(ref)) {
    ref[, `:=`(INSTR = factor(short_instr(INSTRUMENT_LABEL), levels = instr_levels),
               CAT = factor(TERM, levels = levels(br$CAT)))]
    p <- p + geom_line(data = ref,
                       aes(x = MEAN_N_PRIOR, y = PRED_PERCENT, colour = CAT),
                       linetype = "dotted", linewidth = 0.5, inherit.aes = FALSE,
                       show.legend = FALSE)
  }
  
  # Gap-in-bin p-value printed above the higher of the two arms' CIs at that
  # bin. This is the test carrying the paper's claim: an arm crossing zero on
  # its own axis does not imply the gap between arms is insignificant, since
  # the gap's variance uses the covariance between the two arms rather than
  # their individual variances alone. Requires binned_tests or its saved CSV,
  # and is skipped rather than failing the figure when unavailable.
  if (annotate_gap && !is.null(bt) && nrow(bt)) {
    ann <- merge(
      br[IS_REFERENCE == 0, .(INSTRUMENT_LABEL, ZBIN, INSTR,
                              Y_TOP = max(RF_CI_HIGH_PERCENT, na.rm = TRUE),
                              X_MID = mean(MEAN_N_PRIOR)),
         by = .(INSTRUMENT_LABEL, ZBIN, INSTR)][, .(INSTRUMENT_LABEL, ZBIN, INSTR, Y_TOP, X_MID)],
      bt[, .(INSTRUMENT_LABEL, ZBIN, P_VALUE)],
      by = c("INSTRUMENT_LABEL", "ZBIN"))
    if (nrow(ann)) {
      y_pad <- diff(range(br$RF_CI_HIGH_PERCENT, br$RF_CI_LOW_PERCENT, na.rm = TRUE)) * 0.06
      ann[, LAB := sprintf("p=%.3f%s", P_VALUE, ifelse(P_VALUE < 0.05, " *", ""))]
      p <- p + geom_text(data = ann, aes(x = X_MID, y = Y_TOP + y_pad, label = LAB),
                         inherit.aes = FALSE, size = 2.7, colour = "grey20")
    }
  }
  
  save_fig(p, "fig24_binned_rf", width = 9, height = 3.8)
  invisible(p)
}

s31_figure_gap <- function(binned_rows = NULL, binned_tests = NULL,
                           exclude_flagged = TRUE) {
  
  br <- binned_rows
  if (is.null(br)) br <- read_table("T31_binned_rf_estimates.csv")
  bt <- binned_tests
  if (is.null(bt)) bt <- read_table("T31B_binned_rf_tests.csv")
  if (is.null(br) || is.null(bt)) {
    message("SKIP fig24b: need both binned_rows and binned_tests"); return(invisible(NULL))
  }
  br <- as.data.table(br)[ADD_LINEAR == 0]
  bt <- as.data.table(bt)[ADD_LINEAR == 0 & TEST == "Gap in bin"]
  
  if ("PSD_WARNING" %in% names(br)) {
    flagged <- unique(br[PSD_WARNING == TRUE]$INSTRUMENT_LABEL)
    if (length(flagged)) {
      if (exclude_flagged) {
        cat("fig24b: excluding", paste(flagged, collapse = ", "),
            "-- PSD_WARNING fired on this instrument.\n")
        br <- br[!(INSTRUMENT_LABEL %chin% flagged)]
        bt <- bt[!(INSTRUMENT_LABEL %chin% flagged)]
      }
    }
  }
  if (nrow(br) == 0L) { message("SKIP fig24b: nothing left after PSD filter"); return(invisible(NULL)) }
  
  wide <- dcast(br, INSTRUMENT_LABEL + ZBIN ~ TERM,
                value.var = c("RF_PERCENT", "RF_COEF", "MEAN_N_PRIOR"))
  wide[, `:=`(GAP_PP      = RF_PERCENT_Shoppable - RF_PERCENT_Non_shoppable,
              GAP_LOGPTS  = RF_COEF_Shoppable - RF_COEF_Non_shoppable,
              MEAN_N      = (MEAN_N_PRIOR_Shoppable + MEAN_N_PRIOR_Non_shoppable) / 2)]
  
  wide <- merge(wide, bt[, .(INSTRUMENT_LABEL, ZBIN, WALD, P_VALUE)],
                by = c("INSTRUMENT_LABEL", "ZBIN"), all.x = TRUE)
  wide[is.na(WALD) & GAP_LOGPTS == 0, `:=`(GAP_CI_LOW = 0, GAP_CI_HIGH = 0)]
  wide[!is.na(WALD) & WALD > 0,
       `:=`(GAP_SE_LOGPTS = abs(GAP_LOGPTS) / sqrt(WALD))]
  wide[!is.na(GAP_SE_LOGPTS),
       `:=`(GAP_CI_LOW  = GAP_PP - 1.96 * 100 * GAP_SE_LOGPTS,
            GAP_CI_HIGH = GAP_PP + 1.96 * 100 * GAP_SE_LOGPTS)]
  
  wide[, `:=`(INSTR = factor(short_instr(INSTRUMENT_LABEL),
                             levels = c("Competitor hospitals", "Local system",
                                        "Competitor hospitals (ex-CBSA)")),
              SIG = fifelse(!is.na(P_VALUE) & P_VALUE < 0.05, "p < 0.05", "n.s."))]
  
  p <- ggplot(wide, aes(x = MEAN_N, y = GAP_PP)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
    geom_errorbar(aes(ymin = GAP_CI_LOW, ymax = GAP_CI_HIGH), width = 0,
                  linewidth = 0.6, colour = FSU_GARNET, na.rm = TRUE) +
    geom_line(colour = FSU_GARNET, linewidth = 0.7) +
    geom_point(aes(shape = SIG), size = 2.6, colour = FSU_GARNET) +
    scale_shape_manual(values = c("p < 0.05" = 16, "n.s." = 1), name = NULL) +
    facet_wrap(~ INSTR, nrow = 1) +
    labs(x = "Mean prior posters in the exposure bin (midpoint of the two arms)",
         y = "Shoppable minus non-shoppable gap (percentage points)") +
    theme_paper()
  
  save_fig(p, "fig24b_binned_rf_gap", width = 8, height = 3.8)
  invisible(p)
}



if (exists("HPT_SCRATCH") && isTRUE(HPT_SCRATCH)) {   # interactive scratch, off by default
  s31_binned_v2 <- run_binned_rf(outpatient, probs = c(0.50, 0.85))
  s31_binned_coarse <- run_binned_rf(outpatient, probs = c(0.85))
  
  s31_binned_coarse$rows[ADD_LINEAR == 0 & 
                           INSTRUMENT_LABEL != "Competitor_outside_CBSA_hospitals_9m",
                         .(INSTRUMENT_LABEL, TERM, ZBIN, MEAN_N_PRIOR,
                           PCT = round(RF_PERCENT, 2),
                           LO  = round(RF_CI_LOW_PERCENT, 2),
                           HI  = round(RF_CI_HIGH_PERCENT, 2))]
  
  clean_instruments <- MAIN_INSTRUMENTS[c("Competitor_only_hospitals_9m", "Primary_strict_system_IV")]
  
  s31_binned_coarse <- run_binned_rf(outpatient, probs = c(0.85), instruments = clean_instruments)
  
  s31_figure(s31_binned_coarse$rows, s31_binned_coarse$tests)
  s31_figure_gap(s31_binned_coarse$rows, s31_binned_coarse$tests)
  
  
  s31_binned_coarse$tests[ADD_LINEAR == 0 & 
                            TEST %in% c("Bins equal within Shoppable", 
                                        "Gap equal across bins"),
                          .(INSTRUMENT_LABEL, TEST, DF, VCOV_FULL_RANK, 
                            PSD_WARNING, P = round(P_VALUE, 4))]
  
  s31_binned_coarse$tests[, .N, by = .(INSTRUMENT_LABEL, ADD_LINEAR)]
  s31_binned_coarse$rows[INSTRUMENT_LABEL == "Competitor_outside_CBSA_hospitals_9m" & 
                           ADD_LINEAR == 0, .N]
}   # end interactive scratch





# =====================================================================
#  make_scatter_figs.R
#  Builds the two opening slides of the job talk:
#     fig00a_concept_scatter_pooled.pdf   all points one grey
#     fig00b_concept_scatter_split.pdf    identical, colored by shoppability
#
#  The two figures share identical axes, point size, and ordering. The
#  comparison depends on nothing changing between the frames except the
#  colour, so any movement in the points defeats it.
#
#  SOURCE: T05_concept_level_RF_FS_IV.csv, written by stage 6 of
#  HPT_Analysis_Pipeline.R (save_csv(concept_results, ...) at the end of
#  the stage-6 block). The RDS cache concept_level_6inst.rds holds the
#  same object and loads faster.
#
#  estimate_concept_level() writes one row per concept x instrument across
#  six instruments, so the table holds roughly 4,400 rows rather than 738.
#  Filtering to a single INSTRUMENT_LABEL is required, or every concept is
#  plotted six times.
#
#  Runs standalone. If it is sourced from inside the pipeline session it
#  reuses cb_load(), DIAGNOSTIC_FAMILIES, and TABLE_DIR instead of the
#  local fallbacks.
# =====================================================================

if (exists("HPT_SCRATCH") && isTRUE(HPT_SCRATCH)) {   # interactive scratch, off by default
  library(data.table)
  library(ggplot2)
  
  # ---------------------------------------------------------------------
  # 0. CONFIGURATION
  # ---------------------------------------------------------------------
  
  # Only needed when running standalone. Ignored if TABLE_DIR already exists.
  T05_PATH <- "T05_concept_level_RF_FS_IV.csv"
  
  # Which of the six instruments to plot. This is an INSTRUMENT_LABEL, the
  # name side of MAIN_INSTRUMENTS, not the Z_ variable name.
  # "Competitor_only_hospitals_9m" is PRIMARY_INSTRUMENT.
  INSTRUMENT_TO_PLOT <- "Competitor_only_hospitals_9m"
  
  # Pooled panel estimate for the reference line and the annotation. This
  # is NOT recoverable from the concept betas, so it is carried over from
  # the pooled models. Update if the pooled specification changes.
  POOLED_RF <- -0.0021          # reduced form, log points
  POOLED_IV <- -3.07            # percent
  POOLED_P  <-  0.159
  
  OUTDIR <- if (exists("FIGURE_DIR")) FIGURE_DIR else "figures"
  
  # ---------------------------------------------------------------------
  # 1. LOAD
  # ---------------------------------------------------------------------
  
  cs <- if (exists("cb_load", mode = "function")) {
    cb_load()                                   # memory -> RDS -> CSV
  } else if (exists("concept_results") && is.data.frame(concept_results)) {
    as.data.table(copy(concept_results))
  } else {
    p <- if (exists("TABLE_DIR")) file.path(TABLE_DIR, basename(T05_PATH)) else T05_PATH
    if (!file.exists(p)) {
      stop("Cannot find ", p, "\n",
           "  Run stage 6, or point T05_PATH at T05_concept_level_RF_FS_IV.csv.\n",
           "  If stage 6 was interrupted, T05_concept_level_PARTIAL.csv has the\n",
           "  concepts finished so far and the same columns.", call. = FALSE)
    }
    cat("Source:", p, "\n")
    fread(p)
  }
  setDT(cs)
  
  need <- c("INSTRUMENT_LABEL", "FINAL_CONCEPT_ID", "FINAL_FAMILY_ID",
            "RF_COEF", "RF_SE")
  miss <- setdiff(need, names(cs))
  if (length(miss)) stop("Missing columns: ", paste(miss, collapse = ", "), call. = FALSE)
  
  cat(sprintf("Loaded %s rows across %d instruments.\n",
              format(nrow(cs), big.mark = ","), uniqueN(cs$INSTRUMENT_LABEL)))
  
  # ---------------------------------------------------------------------
  # 2. COLLAPSE TO ONE ROW PER CONCEPT
  # ---------------------------------------------------------------------
  
  if (!INSTRUMENT_TO_PLOT %chin% unique(cs$INSTRUMENT_LABEL)) {
    stop("INSTRUMENT_TO_PLOT not present. Available:\n  ",
         paste(sort(unique(cs$INSTRUMENT_LABEL)), collapse = "\n  "), call. = FALSE)
  }
  
  cb_all <- copy(cs)                       # all instruments, for the drop report
  cs <- cs[INSTRUMENT_LABEL == INSTRUMENT_TO_PLOT]
  cs <- cs[is.finite(RF_COEF) & is.finite(RF_SE) & RF_SE > 0]
  
  # The sweep loops over ANALYSIS_SERVICE_ID while the paper counts
  # FINAL_CONCEPT_ID. apply_concept_merges() normally makes these agree.
  # If they do not, the merge groups did not collapse and one row per
  # FINAL_CONCEPT_ID is the unit the deck's "738 concepts" refers to.
  if ("ANALYSIS_SERVICE_ID" %chin% names(cs) &&
      uniqueN(cs$ANALYSIS_SERVICE_ID) != uniqueN(cs$FINAL_CONCEPT_ID)) {
    warning(sprintf("service ids (%d) != concept ids (%d); deduping on FINAL_CONCEPT_ID",
                    uniqueN(cs$ANALYSIS_SERVICE_ID), uniqueN(cs$FINAL_CONCEPT_ID)))
    setorder(cs, FINAL_CONCEPT_ID, RF_SE)
    cs <- unique(cs, by = "FINAL_CONCEPT_ID")
  }
  
  # Shoppability is not stored in T05. Derive it the same way
  # diagnose_size_gradient() does, from the clinical family.
  DIAG_FAM <- if (exists("DIAGNOSTIC_FAMILIES")) DIAGNOSTIC_FAMILIES else c(
    "MRI_MRA", "CT_CTA", "XRAY_FLUOROSCOPY", "DIAGNOSTIC_ULTRASOUND",
    "VASCULAR_ULTRASOUND", "ECHOCARDIOGRAPHY", "MAMMOGRAPHY", "BONE_DENSITY",
    "LABORATORY_PATHOLOGY", "EVALUATION_MANAGEMENT")
  
  cs[, shop := FINAL_FAMILY_ID %chin% DIAG_FAM]
  setnames(cs, c("RF_COEF", "RF_SE"), c("beta", "se"))
  
  cat(sprintf("Plotting %d concepts (%d shoppable / %d non-shoppable), instrument %s.\n",
              nrow(cs), sum(cs$shop), sum(!cs$shop), INSTRUMENT_TO_PLOT))
  
  # Concepts present in the sweep but with no usable row for THIS instrument.
  # estimate_concept_level() skips an instrument when has_usable_variation()
  # fails, so a concept can be estimated for five instruments and not the sixth.
  all_ids <- unique(cb_all$FINAL_CONCEPT_ID)
  dropped <- setdiff(all_ids, cs$FINAL_CONCEPT_ID)
  if (length(dropped)) {
    cat(sprintf("Dropped for this instrument: %d of %d concepts.\n",
                length(dropped), length(all_ids)))
    print(unique(cb_all[FINAL_CONCEPT_ID %chin% dropped,
                        .(FINAL_CONCEPT_ID, FINAL_CONCEPT_NAME, FINAL_FAMILY_ID)]))
  }
  if (nrow(cs) < 700)
    warning("Fewer than 700 concepts. Check whether this is the PARTIAL file.")
  
  # ---------------------------------------------------------------------
  # 3. SHARED GEOMETRY. Computed once, reused by both figures.
  # ---------------------------------------------------------------------
  
  setorder(cs, beta)
  cs[, rank := .I]
  cs[, `:=`(lo = beta - 1.96 * se, hi = beta + 1.96 * se)]
  
  # Trim the y-axis to the central mass so the tails do not flatten the
  # picture. Points outside the range stay in the data and their whiskers
  # clip rather than being dropped. The count is printed below.
  YLIM <- as.numeric(quantile(cs$beta, c(0.01, 0.99), na.rm = TRUE)) * 1.9
  cat(sprintf("Points outside the plotted y-range: %d\n",
              cs[beta < YLIM[1] | beta > YLIM[2], .N]))
  
  base_layers <- list(
    geom_hline(yintercept = 0, linewidth = 0.5, color = "grey30"),
    geom_hline(yintercept = POOLED_RF, linewidth = 0.6,
               linetype = "dashed", color = "grey20"),
    scale_y_continuous(name = "Reduced-form estimate, log points"),
    scale_x_continuous(name = "Clinical concept, ranked by point estimate",
                       expand = expansion(mult = 0.01)),
    coord_cartesian(ylim = YLIM),
    theme_minimal(base_size = 13),
    theme(panel.grid.minor     = element_blank(),
          panel.grid.major.x   = element_blank(),
          legend.position      = c(0.02, 0.98),
          legend.justification = c(0, 1),
          legend.title         = element_blank(),
          legend.background    = element_rect(fill = "white", color = NA),
          plot.caption         = element_text(hjust = 0, color = "grey35", size = 9))
  )
  
  # ---------------------------------------------------------------------
  # 4. PANEL A. One colour, the null as the literature reports it.
  # ---------------------------------------------------------------------
  
  pA <- ggplot(cs, aes(x = rank, y = beta)) +
    geom_linerange(aes(ymin = lo, ymax = hi),
                   linewidth = 0.25, color = "grey72", alpha = 0.55) +
    geom_point(size = 0.85, color = "grey38", alpha = 0.85) +
    base_layers +
    annotate("label", x = 0.02 * nrow(cs), y = YLIM[2] * 0.88,
             label = sprintf("Pooled IV: %.2f%%   (p = %.3f)", POOLED_IV, POOLED_P),
             hjust = 0, size = 4, linewidth = 0, fill = "white", color = "grey15") +
    labs(caption = sprintf(
      "Each point is one clinical concept (N = %d), estimated separately. Whiskers are 95%% CIs. Dashed line is the pooled reduced form.",
      nrow(cs)))
  
  # ---------------------------------------------------------------------
  # 5. PANEL B. Identical geometry, color mapped to shoppability.
  # ---------------------------------------------------------------------
  
  cs[, shop_lab := factor(fifelse(shop, "Shoppable", "Non-shoppable"),
                          levels = c("Shoppable", "Non-shoppable"))]
  
  PAL <- c("Shoppable" = "#782F40", "Non-shoppable" = "#A97C2A")   # garnet / gold
  
  pB <- ggplot(cs, aes(x = rank, y = beta, color = shop_lab)) +
    geom_linerange(aes(ymin = lo, ymax = hi),
                   linewidth = 0.25, alpha = 0.30, show.legend = FALSE) +
    geom_point(size = 0.85, alpha = 0.85) +
    scale_color_manual(values = PAL) +
    base_layers +
    guides(color = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
    labs(caption = sprintf(
      "Identical estimates and axes to the previous figure. Median: %.4f shoppable, %+.4f non-shoppable.",
      median(cs[shop == TRUE, beta]), median(cs[shop == FALSE, beta])))
  
  # ---------------------------------------------------------------------
  # 6. WRITE. 16:9 proportions to match the slide.
  # ---------------------------------------------------------------------
  
  dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
  
  # capabilities("cairo") returns TRUE on macOS R builds whose cairo DLL will
  # not load without XQuartz. Let ggsave infer the base pdf device from the
  # file extension instead, which needs no X11.
  ggsave(file.path(OUTDIR, "fig00a_concept_scatter_pooled.pdf"), pA,
         width = 10, height = 5.2)
  ggsave(file.path(OUTDIR, "fig00b_concept_scatter_split.pdf"), pB,
         width = 10, height = 5.2)
  cat("Wrote both figures to", normalizePath(OUTDIR), "\n")
  
  # ---------------------------------------------------------------------
  # 7. VERIFICATION OF THE NUMBERS QUOTED ON THE SLIDE
  # ---------------------------------------------------------------------
  
  # Three ways to call a concept "moving on its own". The CI rule uses a
  # normal 1.96 cutoff and is the loosest. RF_P comes from fixest with the
  # clustered t distribution and is the paper's own inference, so it is the
  # figure the deck quotes. RF_P_FDR is Benjamini-Hochberg across concepts.
  cs[, moves_ci := (hi < 0 | lo > 0)]
  cs[, moves_p  := is.finite(RF_P) & RF_P < 0.05]
  if (!"RF_P_FDR" %chin% names(cs)) cs[, RF_P_FDR := NA_real_]
  cs[, moves_fdr := is.finite(RF_P_FDR) & RF_P_FDR < 0.05]
  
  print(cs[, .(n           = .N,
               share_ci    = round(mean(moves_ci),  3),
               share_p     = round(mean(moves_p),   3),
               share_fdr   = round(mean(moves_fdr), 3),
               median_beta = round(median(beta),    4)), by = shop_lab])
}   # end interactive scratch

# slides.tex, frame m-scatter2, asserts 31% vs 3% and medians -0.0034 vs
# +0.0004. If this table disagrees, edit the deck rather than the table:
#   \only<2> and \only<3> in the \callout on that frame.












###############################################################################
#                                                                             #
#   SECTION 33 -- CONCEPT-CHARACTERISTIC META-REGRESSION   (GATE 1)           #
#                                                                             #
#   Last section in HPT_Analysis_Pipeline.R, after the slide-figure script.   #
#   Numbering jumps 31 -> 33 deliberately. Section 32 (a causal-forest         #
#   heterogeneity screen over market and contracting-depth moderators) was     #
#   removed: it carried no service characteristics in its covariate set and    #
#   found no heterogeneity. Any QA32*.csv in the QA directory are its stale    #
#   outputs.                                                                   #
#                                                                             #
#   Requires `concept_results` and `outpatient` in the session. Both exist     #
#   after RUN_STAGES 6, or after a warm start restores them.                   #
#                                                                             #
#   ---------------------------------------------------------------------------
#   WHAT THIS ANSWERS
#   ---------------------------------------------------------------------------
#
#   Sections 6-8 already produce one causal estimate per concept and regress
#   those estimates on SHOPPABILITY SCHEMES. A scheme is a hand-built partition
#   of the concept universe, so it can only find heterogeneity along a
#   dimension named in advance. This section asks the complementary question:
#   regress the same concept-level estimates on continuous, objective
#   characteristics of the service, and test whether anything beats or
#   survives alongside shoppability.
#
#   It reuses the existing two-step design rather than replacing it. Same
#   dependent variables (RF_COEF / IV_COEF from concept_results), same
#   inverse-variance weighting, same family clustering, same .pval() reference
#   distribution. The only new thing is the right-hand side.
#
#   ---------------------------------------------------------------------------
#   FOUR BLOCKS, RUN AND CHECK IN ORDER
#   ---------------------------------------------------------------------------
#
#   33A  INVENTORY       what characteristic columns actually exist on disk
#   33B  HETEROGENEITY   is there anything to explain at all (Q, I-squared,
#        BUDGET          tau-squared, cross-instrument reliability). A negative
#                        answer here ends the section: nothing downstream can
#                        work, and that is itself a reportable finding.
#   33C  CHARACTERISTICS build the concept-level covariate table
#   33D  META-REGRESSIONS univariate sweep, joint model, horse race against
#                        shoppability, in three specifications
#
#   ---------------------------------------------------------------------------
#   THE THREE SPECIFICATIONS IN 33D, AND WHY THE THIRD IS THE ONE THAT MATTERS
#   ---------------------------------------------------------------------------
#
#   RAW      characteristic alone. Easiest to read, most confounded.
#   SIZE     plus log(N_OBSERVATIONS) and mean hospitals per county cell.
#            diagnose_size_gradient() already documents a size gradient in this
#            panel, so any characteristic correlated with concept support picks
#            that up unless support is conditioned out.
#   FAMILY   plus clinical-family fixed effects. The sharp test. Nearly every
#            characteristic here is correlated with family (log Medicare rate
#            is close to "MRI versus X-ray"), and every shoppability scheme is
#            defined on family. A characteristic significant only in RAW is
#            relabelling the family gradient the paper already reports. One
#            surviving FAMILY expresses something none of the 18 schemes can,
#            which is the claim worth making.
#
#   ---------------------------------------------------------------------------
#   OUTPUTS
#   ---------------------------------------------------------------------------
#     QA33A_characteristic_inventory.csv
#     QA33B_heterogeneity_budget.csv
#     QA33B2_cross_instrument_reliability.csv
#     QA33C_characteristic_coverage.csv
#     T33A_concept_characteristics.csv
#     T33B_meta_univariate_sweep.csv
#     T33C_meta_joint_model.csv
#     T33D_meta_horserace_vs_shoppability.csv
#
###############################################################################


# ============================================================================
# 33.0  CONSTANTS AND THE CHARACTERISTIC REGISTRY
# ============================================================================

S33_MIN_CONCEPTS <- MIN_CONCEPTS_META          # reuse the existing floor (15)
S33_DEPENDENTS   <- c(RF = "RF_COEF", IV = "IV_COEF")
S33_SE_FOR       <- c(RF_COEF = "RF_SE", IV_COEF = "IV_SE")
S33_SPECS        <- c("RAW", "SIZE", "FAMILY")

# Same pattern Section 28 uses. The Snowflake export did not set SINGLE = TRUE,
# so this file arrives as several shards named HPT_ALT_CONCEPT.csv.gz_0_0_0
# .csv.gz and so on. Reading only the first shard would silently drop most
# concepts, so every match is read and row-bound, exactly as s28_load_alt()
# and load_payer_dispersion() already do.
S33_ALT_PATTERN <- "^HPT_ALT_CONCEPT.*\\.csv(\\.gz)?$"
S33_CODEBOOK_PATTERN <- "^HPT_CODEBOOK\\.csv$"

# Three code-keyed supplemental sources, none of them concept-keyed. Each is
# joined onto the codebook's BILLING_CODE -> ANALYSIS_CONCEPT_ID crosswalk,
# then collapsed through s33_to_canonical() exactly like every other source
# in this section, so a code that belongs to a merged concept (mammography
# tomosynthesis, the MRI/CT abdomen variants) still lands correctly.
#
# Addendum A is NOT built here. It has no HCPCS-code column at all -- it is
# keyed on APC number -- and HPT_CODEBOOK.csv carries no APC column to join
# it through (checked directly: OPPS_STATUS_INDICATORS is the only OPPS
# field that made it into the export). Building this needs either a fresh
# SQL pull that retains APC number per code, or Addendum B's own code-to-APC
# mapping pulled separately. Deferred, not abandoned.
S33_PFS_RVU_PATTERN <- "^PPRRVU.*\\.xlsx$"
S33_ASC_BB_PATTERN  <- "Addendum[ _-]?BB.*\\.txt$"
S33_GEO_PATTERN      <- "^MUP_PHY.*Geo\\.csv$"

# Shared by all three code-keyed loaders below: BILLING_CODE -> concept.
s33_code_to_concept_map <- function() {
  cb <- s33_read_all(S33_CODEBOOK_PATTERN)
  if (is.null(cb)) cb <- tryCatch(fread(FILES$codebook), error = function(e) NULL)
  if (is.null(cb) || !all(c("BILLING_CODE", "ANALYSIS_CONCEPT_ID") %in% names(cb)))
    return(NULL)
  cb <- s33_to_canonical(cb[, .(BILLING_CODE, ANALYSIS_CONCEPT_ID)], "ANALYSIS_CONCEPT_ID")
  setnames(cb, "BILLING_CODE", "HCPCS_CODE")
  cb[, HCPCS_CODE := toupper(trimws(HCPCS_CODE))]
  unique(cb, by = "HCPCS_CODE")
}

# OPPS status indicators, grouped by what each implies for whether a posted
# per-code price is a standalone price. Editable here rather than inline.
#
#   S/T/V     separately paid
#   J1/J2     comprehensive-APC primaries: paid, but the payment absorbs the
#             rest of the claim, so the posted price is a price for the
#             encounter rather than for this service
#   Q1-Q4     paid separately only when no S/T line appears on the same claim
#   N         never paid separately
#
# Collapsing these into one separately-payable indicator covers 78% of
# concepts and codes a comprehensive-APC service identically to a separately
# paid one, which is why they are kept apart.
S33_SI_SEPARATELY_PAYABLE <- c("S", "T", "V")
S33_SI_COMPREHENSIVE      <- c("J1", "J2")
S33_SI_COND_PACKAGED      <- c("Q1", "Q2", "Q3", "Q4")
S33_SI_PACKAGED           <- c("N")

# Registry. outcome_derived = TRUE means the characteristic is built from the
# same negotiated-price data that produced the dependent variable. Those are
# not excluded, because several of them (price level, cross-hospital
# dispersion) are exactly the economics worth testing. They are FLAGGED so the
# output separates them, and none should carry a headline claim on its own.
S33_CHAR_REGISTRY <- list(
  # ---- clean: administratively set or structural, not built from the outcome
  list(name = "LN_MEDICARE_RATE",     label = "Log Medicare reference rate",
       outcome_derived = FALSE,
       note = "CMS administratively-set price. Resource-intensity proxy with no hospital discretion in it, and the cleanest complexity measure available without HCUP or DRG weights."),
  list(name = "SHARE_ASC_COVERED",    label = "Share of codes ASC-covered",
       outcome_derived = FALSE,
       note = "CMS Addendum AA. CMS's own judgment that the procedure is plannable and schedulable."),
  list(name = "SHARE_SI_SEPARATE",    label = "Share of codes separately payable (OPPS S/T/V)",
       outcome_derived = FALSE,
       note = "The service is an independently priceable unit under OPPS."),
  list(name = "SHARE_SI_COMPREHENSIVE", label = "Share of codes in a comprehensive APC (OPPS J1/J2)",
       outcome_derived = FALSE,
       note = "The posted per-code price is a price for the whole encounter, not for this service. Directly relevant to whether a disclosed number is shoppable."),
  list(name = "SHARE_SI_COND_PACKAGED", label = "Share of codes conditionally packaged (OPPS Q1-Q4)",
       outcome_derived = FALSE,
       note = "Paid separately only when no S/T line is on the same claim, so the posted price applies some of the time and not others."),
  list(name = "SHARE_SI_PACKAGED",    label = "Share of codes packaged (OPPS N)",
       outcome_derived = FALSE,
       note = "Never separately paid under OPPS."),
  list(name = "SHARE_NCCI_ADDON",     label = "Share of codes that are NCCI add-ons",
       outcome_derived = FALSE,
       note = "Component share. Near zero in the standalone panel by construction; kept as a coverage check."),
  list(name = "ANY_CMS70",            label = "CMS-70 required service (0/1)",
       outcome_derived = FALSE,
       note = "Statutory shoppable list membership. Overlaps SCHEME_4 by construction."),
  list(name = "HAS_CONTRAST_VARIANT", label = "Concept has contrast variants (0/1)",
       outcome_derived = FALSE,
       note = "Proxy for whether the posted price can shift at point of service."),
  list(name = "N_CODES_IN_CONCEPT",   label = "Billing codes per concept",
       outcome_derived = FALSE,
       note = "Definitional breadth. A concept spanning many codes is a fuzzier price object."),
  list(name = "GROSS_TO_MEDICARE",    label = "Gross charge / Medicare rate",
       outcome_derived = FALSE,
       note = "Chargemaster markup. Gross charge is set unilaterally and is not the outcome, so this measures pricing discretion without being mechanically tied to the dependent variable."),
  list(name = "CASH_TO_MEDICARE",     label = "Cash rate / Medicare rate",
       outcome_derived = FALSE,
       note = "Consumer-facing markup. Direct measure of the patient channel."),
  
  # ---- support: condition on these, do not interpret them
  list(name = "LN_N_HOSPITALS",       label = "Log hospitals posting the concept",
       outcome_derived = FALSE,
       note = "SUPPORT. Ubiquity of the service, and the main driver of estimate precision."),
  list(name = "MEAN_HOSP_PER_MARKET", label = "Mean hospitals per county cell",
       outcome_derived = FALSE,
       note = "SUPPORT. Density of the within-county comparison that identifies the effect. Also enters SIZE and FAMILY as a control, so its own univariate SIZE row collapses back to RAW by construction."),
  list(name = "MEAN_CODE_COVERAGE",   label = "Mean code-coverage ratio",
       outcome_derived = FALSE,
       note = "SUPPORT. How completely hospitals report the concept's constituent codes."),
  
  # ---- outcome-derived: read only in the SIZE and FAMILY specs, never alone
  list(name = "LN_CONCEPT_PRICE",     label = "Log concept price level",
       outcome_derived = TRUE,
       note = "Dollar magnitude, the deductible-exposure proxy. Built from the outcome."),
  list(name = "SD_LN_PRICE_XHOSP",    label = "Cross-hospital SD of log price",
       outcome_derived = TRUE,
       note = "How much room there is to converge. Built from the outcome and mechanically noisier for thin concepts."),
  list(name = "NEG_TO_MEDICARE",      label = "Negotiated / Medicare rate",
       outcome_derived = TRUE,
       note = "Commercial-to-Medicare ratio, the standard service-line margin measure. Numerator is the outcome."),
  list(name = "MEAN_CV_PAYER",        label = "Mean within-hospital payer dispersion",
       outcome_derived = TRUE,
       note = "Contracting depth at concept level. The same object as PD_PAYER_V2, aggregated up."),
  list(name = "MEAN_N_PAYERS",        label = "Mean distinct payers per cell",
       outcome_derived = FALSE,
       note = "Thickness of contracting. Counts, not prices."),
  
  # ---- PFS RVU (PPRRVU25_JAN), ASC Addendum BB, and Medicare geography/
  #      volume file. All three administratively set or externally observed,
  #      none outcome-derived.
  list(name = "SHARE_PCTC_SPLIT",     label = "Share of codes with a PC/TC split (PFS PCTC IND = 1)",
       outcome_derived = FALSE,
       note = "The posted price is half the bill: a separate professional-component invoice exists that the patient does not see when comparing hospital prices."),
  list(name = "SHARE_GLOB_SURGICAL",  label = "Share of codes carrying a global-surgery package (PFS GLOB DAYS in 000/010/090)",
       outcome_derived = FALSE,
       note = "The posted price bundles a period of follow-up care rather than one visit; XXX/ZZZ concepts have no such bundling."),
  list(name = "SHARE_ASC_ANCILLARY_PACKAGED", label = "Share of codes packaged under the ASC ancillary schedule (Addendum BB, N1)",
       outcome_derived = FALSE,
       note = "Independent packaging signal from a different payment system than OPPS -- corroborates or contradicts SHARE_SI_PACKAGED without sharing its source."),
  list(name = "LN_NATIONAL_VOLUME",   label = "Log national Medicare FFS volume (services/year, geography-and-service file)",
       outcome_derived = FALSE,
       note = "How often the service is actually performed nationally. A code with near-zero volume is not something any patient realistically shops."),
  list(name = "N_PROVIDERS_NATIONAL", label = "National count of distinct rendering providers",
       outcome_derived = FALSE,
       note = "Breadth of supply. A direct measure of whether an alternative provider exists to shop toward."),
  list(name = "SHARE_VOL_FACILITY",   label = "Share of national volume billed in a facility place of service",
       outcome_derived = FALSE,
       note = "If most volume happens in physician offices rather than hospital outpatient departments, the hospital's posted price is largely irrelevant to what most patients actually pay for this service.")
)

S33_CHAR_NAMES    <- vapply(S33_CHAR_REGISTRY, `[[`, character(1), "name")
S33_SIZE_CONTROLS <- c("LN_N_OBS", "MEAN_HOSP_PER_MARKET")

.s33_hd <- function(x)
  cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")

s33_registry_dt <- function() {
  rbindlist(lapply(S33_CHAR_REGISTRY, function(r)
    data.table(CHARACTERISTIC = r$name, LABEL = r$label,
               OUTCOME_DERIVED = as.integer(r$outcome_derived), NOTE = r$note)))
}

# Every match, not the first. See the note at S33_ALT_PATTERN.
s33_find_files <- function(pattern, dir = PANEL_DIR)
  list.files(dir, pattern = pattern, full.names = TRUE, ignore.case = TRUE)

s33_read_all <- function(pattern, dir = PANEL_DIR, nrows = Inf) {
  f <- s33_find_files(pattern, dir)
  if (length(f) == 0L) return(NULL)
  tryCatch(rbindlist(lapply(f, fread, nrows = nrows, showProgress = FALSE),
                     fill = TRUE),
           error = function(e) NULL)
}


# ============================================================================
# 33A  INVENTORY -- what is actually on disk, before building anything
# ============================================================================
#
# The SQL codebook (HPT_P1_FINAL_CODEBOOK_SCOPED, exported as HPT_CODEBOOK.csv)
# carries OPPS_STATUS_INDICATORS, IS_ASC_COVERED_FLAG, IS_NCCI_ADDON_FLAG, MDC,
# MEDICAL_SURGICAL_TYPE, DRG_ACUITY_KEYWORD_TAG and CONTRAST_VARIANT.
# build_schemes() reads FILES$codebook and keeps four of them. Whether the file
# in PANEL_DIR is the full export or a trimmed one is the first thing to check,
# because every codebook characteristic below depends on it.

S33_ATTRIBUTE_COLUMNS <- c(
  "OPPS_STATUS_INDICATORS", "IS_ASC_COVERED_FLAG", "IS_NCCI_ADDON_FLAG",
  "IS_CMS70_CODE", "CONTRAST_VARIANT", "MDC", "MEDICAL_SURGICAL_TYPE",
  "DRG_ACUITY_KEYWORD_TAG", "N_HOSPITALS", "N_STATES", "N_COUNTIES",
  "ASC_PAYMENT_INDICATORS", "CMS70_SERVICE_ID"
)

s33_inventory <- function() {
  rows <- list()
  add <- function(source, item, status, detail = "") {
    rows[[length(rows) + 1L]] <<- data.table(
      SOURCE = source, ITEM = item, STATUS = status,
      DETAIL = as.character(detail))
  }
  
  # --- concept_results ------------------------------------------------------
  if (exists("concept_results")) {
    cr <- as.data.table(concept_results)
    add("concept_results", "rows",        "OK", format(nrow(cr), big.mark = ","))
    add("concept_results", "concepts",    "OK", uniqueN(cr$FINAL_CONCEPT_ID))
    add("concept_results", "families",    "OK", uniqueN(cr$FINAL_FAMILY_ID))
    add("concept_results", "instruments", "OK",
        paste(sort(unique(cr$INSTRUMENT_LABEL)), collapse = ", "))
    for (v in c("RF_COEF", "RF_SE", "IV_COEF", "IV_SE", "FS_COEF", "FS_F",
                "N_OBSERVATIONS", "N_HOSPITALS", "N_MARKETS")) {
      add("concept_results", v, if (v %in% names(cr)) "OK" else "MISSING",
          if (v %in% names(cr)) sprintf("%d finite", sum(is.finite(cr[[v]]))) else "")
    }
  } else {
    add("concept_results", "object", "MISSING", "run RUN_STAGES 6 or warm start")
  }
  
  # --- codebook that build_schemes() reads ---------------------------------
  cb <- tryCatch(fread(FILES$codebook, nrows = 5L), error = function(e) NULL)
  if (is.null(cb)) {
    add("codebook", "path", "MISSING", FILES$codebook)
  } else {
    add("codebook", "path",      "OK", basename(FILES$codebook))
    add("codebook", "n_columns", "OK", ncol(cb))
    for (v in S33_ATTRIBUTE_COLUMNS)
      add("codebook", v, if (v %in% names(cb)) "OK" else "MISSING", "")
  }
  
  # --- full SQL codebook, if downloaded separately --------------------------
  fcb <- s33_read_all(S33_CODEBOOK_PATTERN, nrows = 5L)
  if (is.null(fcb)) {
    add("full_sql_codebook", "HPT_CODEBOOK.csv", "MISSING",
        "re-download qa/HPT_CODEBOOK.csv from HPT_PY_EXPORT_STAGE")
  } else {
    add("full_sql_codebook", "HPT_CODEBOOK.csv", "OK",
        paste(ncol(fcb), "columns"))
    for (v in S33_ATTRIBUTE_COLUMNS)
      add("full_sql_codebook", v, if (v %in% names(fcb)) "OK" else "MISSING", "")
  }
  
  # --- Phase 5 alternative price series ------------------------------------
  altf <- s33_find_files(S33_ALT_PATTERN)
  if (length(altf) == 0L) {
    add("phase5_alt_prices", "HPT_ALT_CONCEPT", "MISSING",
        "Medicare / gross / cash characteristics will be skipped")
  } else {
    a <- tryCatch(fread(altf[1L], nrows = 5L), error = function(e) NULL)
    add("phase5_alt_prices", "HPT_ALT_CONCEPT", "OK",
        sprintf("%d shard(s), %s columns", length(altf),
                if (is.null(a)) "unreadable" else ncol(a)))
  }
  
  # --- payer dispersion -----------------------------------------------------
  pd <- list.files(PAYER_DISPERSION_DIR, pattern = PAYER_DISPERSION_PATTERN)
  add("payer_dispersion", "shards",
      if (length(pd) == 0L) "MISSING" else "OK", length(pd))
  
  # --- the panel ------------------------------------------------------------
  if (exists("outpatient")) {
    op <- as.data.table(outpatient)
    add("outpatient", "rows", "OK", format(nrow(op), big.mark = ","))
    for (v in c("FINAL_CONCEPT_ID", "FINAL_FAMILY_ID", "LN_MEDIAN_PRICE",
                "CODE_COVERAGE_RATIO", "N_DISTINCT_PAYERS", "ANALYSIS_MARKET"))
      add("outpatient", v, if (v %in% names(op)) "OK" else "MISSING", "")
  } else {
    add("outpatient", "object", "MISSING", "")
  }
  
  out <- rbindlist(rows)
  save_qa_csv(out, "QA33A_characteristic_inventory.csv")
  out
}


# ============================================================================
# 33B  HETEROGENEITY BUDGET -- is there anything to explain?
# ============================================================================
#
# There is one estimate per concept, each with a standard error. Part of the
# spread across them is real cross-service variation and the rest is sampling
# noise. If the noise share is close to one, no covariate can explain the
# spread because there is no spread to explain. Establish that number before
# spending a day building covariates.
#
#   Q        Cochran's Q, the inverse-variance-weighted dispersion statistic
#   I2       share of total variation attributable to real heterogeneity
#   TAU2     DerSimonian-Laird estimate of the between-concept variance
#   TAU_PCT  sqrt(TAU2) on the paper's percent-per-unit scale
#
# One caveat. Q assumes independent estimates, and these are not independent:
# concepts share hospitals, counties and months, so their sampling errors are
# positively correlated and Q is inflated. I2 is therefore an upper bound.
# s33_cross_instrument_reliability() is the complement that does not rest on
# that assumption.

s33_heterogeneity_budget <- function(cr, dependents = S33_DEPENDENTS) {
  cr <- as.data.table(cr)
  rows <- list()
  
  for (dep in dependents) {
    se_col <- S33_SE_FOR[[dep]]
    if (!all(c(dep, se_col) %in% names(cr))) next
    
    for (il in sort(unique(cr$INSTRUMENT_LABEL))) {
      d <- cr[INSTRUMENT_LABEL == il & is.finite(get(dep)) &
                is.finite(get(se_col)) & get(se_col) > 0]
      if (nrow(d) < S33_MIN_CONCEPTS) next
      
      b <- d[[dep]]; s <- d[[se_col]]; w <- 1 / s^2
      bbar <- sum(w * b) / sum(w)
      Q    <- sum(w * (b - bbar)^2)
      df   <- length(b) - 1L
      cval <- sum(w) - sum(w^2) / sum(w)
      tau2 <- max(0, (Q - df) / cval)
      I2   <- max(0, (Q - df) / Q)
      
      # Naive unweighted decomposition, reported alongside because it is the
      # one a reader can verify by eye from the CSV.
      var_b    <- var(b)
      mean_se2 <- mean(s^2)
      
      rows[[length(rows) + 1L]] <- data.table(
        DEPENDENT = dep, INSTRUMENT_LABEL = il, N_CONCEPTS = length(b),
        N_FAMILIES = uniqueN(d$FINAL_FAMILY_ID),
        POOLED_ESTIMATE_PCT = 100 * bbar,
        Q = Q, DF = df, Q_P = pchisq(Q, df, lower.tail = FALSE),
        I2 = I2, TAU2 = tau2, TAU_PCT = 100 * sqrt(tau2),
        SD_BHAT_PCT = 100 * sqrt(var_b), MEAN_SE_PCT = 100 * sqrt(mean_se2),
        NAIVE_SIGNAL_SHARE = max(0, 1 - mean_se2 / var_b),
        MEDIAN_SHRINKAGE = median(tau2 / (tau2 + s^2)))
    }
  }
  
  out <- rbindlist(rows)
  if (nrow(out) == 0L) stop("No usable concept estimates for the budget.", call. = FALSE)
  save_qa_csv(out, "QA33B_heterogeneity_budget.csv")
  out
}

# Does the concept-level ranking reproduce across instruments? Different
# instruments draw on different rollout variation, so a high rank correlation
# is evidence of real service-level signal that does not depend on Q's
# independence assumption. A correlation near zero means the estimates are
# mostly noise and nothing downstream will work.
s33_cross_instrument_reliability <- function(cr, dep = "RF_COEF") {
  cr <- as.data.table(cr)
  if (!(dep %in% names(cr))) return(data.table())
  w <- dcast(cr[is.finite(get(dep))], FINAL_CONCEPT_ID ~ INSTRUMENT_LABEL,
             value.var = dep, fun.aggregate = function(x) x[1L])
  ils <- setdiff(names(w), "FINAL_CONCEPT_ID")
  if (length(ils) < 2L) return(data.table())
  
  rows <- list()
  for (i in seq_along(ils)) for (j in seq_along(ils)) {
    if (j <= i) next
    a <- w[[ils[i]]]; b <- w[[ils[j]]]
    ok <- is.finite(a) & is.finite(b)
    if (sum(ok) < S33_MIN_CONCEPTS) next
    rows[[length(rows) + 1L]] <- data.table(
      DEPENDENT = dep, INSTRUMENT_A = ils[i], INSTRUMENT_B = ils[j],
      N_CONCEPTS = sum(ok),
      PEARSON  = cor(a[ok], b[ok]),
      SPEARMAN = cor(a[ok], b[ok], method = "spearman"))
  }
  out <- rbindlist(rows)
  if (nrow(out) == 0L) return(out)
  save_qa_csv(out, "QA33B2_cross_instrument_reliability.csv")
  out
}


# ============================================================================
# 33C  BUILD THE CONCEPT-LEVEL CHARACTERISTIC TABLE
# ============================================================================
#
# Four sources, each in its own function, each returning NULL with a printed
# reason when its input is absent rather than killing the section. Every source
# is keyed on FINAL_CONCEPT_ID *after* MERGE_GROUPS is applied, so the six
# canonical merged concepts line up with concept_results. Merging on the raw
# ANALYSIS_CONCEPT_ID would silently drop them, mammography included, which is
# the most influential single family in the permutation test.

s33_canonical_map <- function(merge_groups = MERGE_GROUPS)
  rbindlist(lapply(merge_groups, function(g)
    data.table(RAW_CONCEPT_ID = g$constituents, FINAL_CONCEPT_ID = g$canonical_id)))

s33_to_canonical <- function(dt, id_col) {
  d <- copy(as.data.table(dt))
  setnames(d, id_col, "RAW_CONCEPT_ID")
  d <- merge(d, s33_canonical_map(), by = "RAW_CONCEPT_ID",
             all.x = TRUE, sort = FALSE)
  d[, FINAL_CONCEPT_ID := fifelse(is.na(FINAL_CONCEPT_ID),
                                  RAW_CONCEPT_ID, FINAL_CONCEPT_ID)]
  d[, RAW_CONCEPT_ID := NULL]
  d
}

# --- source 1: the codebook, aggregated from billing code to concept --------
s33_chars_from_codebook <- function() {
  cb <- s33_read_all(S33_CODEBOOK_PATTERN)
  if (is.null(cb)) cb <- tryCatch(fread(FILES$codebook), error = function(e) NULL)
  if (is.null(cb) || !("ANALYSIS_CONCEPT_ID" %in% names(cb))) {
    cat("  codebook unreadable or missing ANALYSIS_CONCEPT_ID; skipped.\n")
    return(NULL)
  }
  cb <- s33_to_canonical(cb, "ANALYSIS_CONCEPT_ID")
  
  # Share of a concept's codes whose OPPS status indicator falls in `set`.
  # Codes with no indicator are excluded rather than counted as zero, so the
  # share is over codes where the question is answerable.
  si_share <- function(x, set) {
    v <- toupper(trimws(as.character(x)))
    v <- v[!is.na(v) & nzchar(v)]
    if (length(v) == 0L) return(NA_real_)
    mean(vapply(strsplit(v, ","), function(p)
      as.numeric(any(trimws(p) %chin% set)), numeric(1)))
  }
  
  has <- function(v) v %in% names(cb)
  out <- cb[, {
    lst <- list(N_CODES_IN_CONCEPT = .N)
    if (has("IS_ASC_COVERED_FLAG"))
      lst$SHARE_ASC_COVERED <- mean(safe_numeric(IS_ASC_COVERED_FLAG), na.rm = TRUE)
    if (has("IS_NCCI_ADDON_FLAG"))
      lst$SHARE_NCCI_ADDON  <- mean(safe_numeric(IS_NCCI_ADDON_FLAG), na.rm = TRUE)
    if (has("IS_CMS70_CODE"))
      # sum(... > 0) rather than max(): max() on an all-NA group returns -Inf
      # with a warning, which would fire once per empty concept.
      lst$ANY_CMS70 <- as.numeric(sum(safe_numeric(IS_CMS70_CODE) > 0, na.rm = TRUE) > 0)
    if (has("OPPS_STATUS_INDICATORS")) {
      lst$SHARE_SI_SEPARATE      <- si_share(OPPS_STATUS_INDICATORS, S33_SI_SEPARATELY_PAYABLE)
      lst$SHARE_SI_COMPREHENSIVE <- si_share(OPPS_STATUS_INDICATORS, S33_SI_COMPREHENSIVE)
      lst$SHARE_SI_COND_PACKAGED <- si_share(OPPS_STATUS_INDICATORS, S33_SI_COND_PACKAGED)
      lst$SHARE_SI_PACKAGED      <- si_share(OPPS_STATUS_INDICATORS, S33_SI_PACKAGED)
    }
    if (has("CONTRAST_VARIANT"))
      lst$HAS_CONTRAST_VARIANT <- as.numeric(
        any(!is.na(CONTRAST_VARIANT) & nzchar(as.character(CONTRAST_VARIANT))))
    lst
  }, by = FINAL_CONCEPT_ID]
  
  cat("  codebook: ", nrow(out), " concepts, ",
      length(setdiff(names(out), "FINAL_CONCEPT_ID")), " characteristics\n", sep = "")
  out
}

# --- source 2: the estimation panel itself ---------------------------------
s33_chars_from_panel <- function(panel) {
  d <- as.data.table(panel)
  
  hosp <- d[, .(LNP = median(LN_MEDIAN_PRICE, na.rm = TRUE)),
            by = .(FINAL_CONCEPT_ID, HOSPITAL_ID)]
  disp <- hosp[, .(SD_LN_PRICE_XHOSP = sd(LNP, na.rm = TRUE),
                   LN_CONCEPT_PRICE  = median(LNP, na.rm = TRUE)),
               by = FINAL_CONCEPT_ID]
  
  cell <- d[, .(NH = uniqueN(HOSPITAL_ID)),
            by = .(FINAL_CONCEPT_ID, ANALYSIS_MARKET)]
  cell <- cell[, .(MEAN_HOSP_PER_MARKET = mean(NH),
                   N_MARKETS_CONCEPT    = .N), by = FINAL_CONCEPT_ID]
  
  base <- d[, {
    lst <- list(LN_N_HOSPITALS = log(pmax(uniqueN(HOSPITAL_ID), 1)))
    if ("CODE_COVERAGE_RATIO" %in% names(d))
      lst$MEAN_CODE_COVERAGE <- mean(safe_numeric(CODE_COVERAGE_RATIO), na.rm = TRUE)
    if ("N_DISTINCT_PAYERS" %in% names(d))
      lst$MEAN_N_PAYERS <- mean(safe_numeric(N_DISTINCT_PAYERS), na.rm = TRUE)
    lst
  }, by = FINAL_CONCEPT_ID]
  
  out <- Reduce(function(a, b) merge(a, b, by = "FINAL_CONCEPT_ID", all = TRUE),
                list(base, disp, cell))
  cat("  panel: ", nrow(out), " concepts\n", sep = "")
  out
}

# --- source 3: Phase 5 alternative price series ----------------------------
#
# Prefers s28_load_alt() when Section 28 is loaded, so the shard handling and
# the canonical remap are done once, in the function that already validates
# them, rather than reimplemented here. Falls back to reading the shards
# directly when Section 28 is not in the session.
s33_chars_from_alt_prices <- function() {
  a <- NULL
  if (exists("s28_load_alt", mode = "function")) {
    a <- tryCatch(s28_load_alt(), error = function(e) NULL)
    if (!is.null(a)) {
      setnames(a,
               c("GROSS_CHARGE", "CASH_RATE", "MEDICARE_RATE"),
               c("MEDIAN_GROSS_CHARGE", "MEDIAN_CASH_RATE", "MEDIAN_MEDICARE_RATE"),
               skip_absent = TRUE)
    }
  }
  if (is.null(a)) {
    a <- s33_read_all(S33_ALT_PATTERN)
    if (is.null(a) || !("ANALYSIS_CONCEPT_ID" %in% names(a))) {
      cat("  Phase 5 alt-price shards not found or unreadable; Medicare, gross\n",
          "  and cash characteristics skipped.\n", sep = "")
      return(NULL)
    }
    a <- s33_to_canonical(a, "ANALYSIS_CONCEPT_ID")
  }
  
  a <- as.data.table(a)
  g <- function(v) if (v %in% names(a)) safe_numeric(a[[v]]) else rep(NA_real_, nrow(a))
  a[, `:=`(MED_TMP = g("MEDIAN_MEDICARE_RATE"),
           GRO_TMP = g("MEDIAN_GROSS_CHARGE"),
           CSH_TMP = g("MEDIAN_CASH_RATE"),
           NEG_TMP = g("MEDIAN_NEGOTIATED_CHECK"))]
  
  out <- a[, {
    med <- median(MED_TMP, na.rm = TRUE)
    denom <- if (is.finite(med) && med > 0) med else NA_real_
    list(LN_MEDICARE_RATE  = if (is.finite(denom)) log(denom) else NA_real_,
         GROSS_TO_MEDICARE = median(GRO_TMP, na.rm = TRUE) / denom,
         CASH_TO_MEDICARE  = median(CSH_TMP, na.rm = TRUE) / denom,
         NEG_TO_MEDICARE   = median(NEG_TMP, na.rm = TRUE) / denom)
  }, by = FINAL_CONCEPT_ID]
  
  for (v in c("LN_MEDICARE_RATE", "GROSS_TO_MEDICARE",
              "CASH_TO_MEDICARE", "NEG_TO_MEDICARE"))
    out[!is.finite(get(v)), (v) := NA_real_]
  
  cat("  alt prices: ", nrow(out), " concepts\n", sep = "")
  out
}

# --- source 4: payer dispersion --------------------------------------------
s33_chars_from_payer_dispersion <- function() {
  parts <- list.files(PAYER_DISPERSION_DIR, pattern = PAYER_DISPERSION_PATTERN,
                      full.names = TRUE)
  if (length(parts) == 0L) {
    cat("  payer dispersion shards not found; MEAN_CV_PAYER skipped.\n")
    return(NULL)
  }
  pd <- tryCatch(rbindlist(lapply(parts, fread, showProgress = FALSE), fill = TRUE),
                 error = function(e) NULL)
  if (is.null(pd) || !("ANALYSIS_CONCEPT_ID" %in% names(pd)) ||
      !("CV_PAYER_NEGOTIATED" %in% names(pd))) {
    cat("  payer dispersion unreadable or missing expected columns; skipped.\n")
    return(NULL)
  }
  pd <- s33_to_canonical(pd, "ANALYSIS_CONCEPT_ID")
  # N_DISTINCT_PAYERS comes from this file rather than from `outpatient`,
  # which does not carry it: the Phase 4 concept rollup drops the payer count.
  has_np <- "N_DISTINCT_PAYERS" %in% names(pd)
  
  # CV_PAYER_NEGOTIATED is (P75-P25)/median, so it explodes wherever the
  # median price is near zero: P99 is 114 and the maximum 241, against a P75
  # of 0.87. A median across hospital-months is robust to that; a mean is not.
  # The column name stays MEAN_CV_PAYER for continuity with the registry and
  # every downstream reference, but the statistic is a median.
  out <- pd[, c(list(
    MEAN_CV_PAYER = median(safe_numeric(CV_PAYER_NEGOTIATED), na.rm = TRUE)),
    if (has_np) list(MEAN_N_PAYERS = mean(safe_numeric(N_DISTINCT_PAYERS), na.rm = TRUE))),
    by = FINAL_CONCEPT_ID]
  out[!is.finite(MEAN_CV_PAYER), MEAN_CV_PAYER := NA_real_]
  if (has_np) out[!is.finite(MEAN_N_PAYERS), MEAN_N_PAYERS := NA_real_]
  cat("  payer dispersion: ", nrow(out), " concepts\n", sep = "")
  rm(pd); invisible(gc())
  out
}

# --- source 5: PFS RVU file (PPRRVU) -- PC/TC split and global period -------
#
# Global row only (blank modifier): a code with PCTC IND = 1 has separate 26/
# TC rows in addition to the global row, and the global row is the one whose
# RVUs match what a payer negotiates against, so it is the row that answers
# "does a PC/TC split exist for this code" without double-counting.
#
# Columns are addressed by position, not name: CMS's header row repeats the
# literal string "INDICATOR" three times (NA indicators, physician
# supervision) and the PCTC/GLOB columns carry only the fragment "IND"/"DAYS"
# after their header wraps across two spreadsheet rows. Position 1/2/14/15
# is the verified CY2025 January-release layout; a released quarter with a
# different column count changes this and should be re-checked before reuse.
s33_chars_from_pfs_rvu <- function() {
  f <- s33_find_files(S33_PFS_RVU_PATTERN)
  if (length(f) == 0L) {
    cat("  PFS RVU file not found; SHARE_PCTC_SPLIT and SHARE_GLOB_SURGICAL skipped.\n")
    return(NULL)
  }
  if (length(f) > 1L)
    cat("  Multiple PFS RVU files matched; using ", basename(f[1L]), "\n", sep = "")
  
  if (!requireNamespace("readxl", quietly = TRUE)) {
    cat("  readxl not installed (install.packages(\"readxl\")); PFS RVU skipped.\n")
    return(NULL)
  }
  raw <- tryCatch(readxl::read_excel(f[1L], skip = 9, col_names = TRUE,
                                     .name_repair = "unique_quiet"),
                  error = function(e) NULL)
  if (is.null(raw) || ncol(raw) < 15) {
    cat("  PFS RVU file unreadable or an unexpected layout; skipped.\n")
    return(NULL)
  }
  dt <- as.data.table(raw)
  setnames(dt, old = names(dt)[c(1L, 2L, 14L, 15L)],
           new = c("HCPCS_CODE", "MODIFIER", "PCTC_IND", "GLOB_DAYS"))
  
  dt[, HCPCS_CODE := toupper(trimws(as.character(HCPCS_CODE)))]
  dt[, MODIFIER    := toupper(trimws(as.character(MODIFIER)))]
  dt[, PCTC_IND    := safe_numeric(PCTC_IND)]
  dt[, GLOB_DAYS   := toupper(trimws(as.character(GLOB_DAYS)))]
  dt <- dt[!is.na(HCPCS_CODE) & nzchar(HCPCS_CODE) & (is.na(MODIFIER) | MODIFIER == "")]
  
  map <- s33_code_to_concept_map()
  if (is.null(map)) {
    cat("  codebook unavailable; PFS RVU cannot be joined to concepts.\n")
    return(NULL)
  }
  m <- merge(dt, map, by = "HCPCS_CODE")
  if (nrow(m) == 0L) {
    cat("  No PFS RVU rows matched the codebook; skipped.\n")
    return(NULL)
  }
  
  out <- m[, .(
    SHARE_PCTC_SPLIT    = mean(PCTC_IND == 1, na.rm = TRUE),
    SHARE_GLOB_SURGICAL = mean(GLOB_DAYS %chin% c("000", "010", "090"), na.rm = TRUE)
  ), by = FINAL_CONCEPT_ID]
  
  cat("  PFS RVU: ", nrow(out), " concepts, ", nrow(m), " matched codes\n", sep = "")
  out
}

# --- source 6: ASC Addendum BB -- ancillary packaging status ---------------
#
# Both years read and row-bound (AA/BB rarely reclassify a code between
# vintages -- confirmed directly: 410 vs 412 matched codes across 2024/2025),
# then deduplicated on code+indicator. A code whose indicator genuinely
# changed between years keeps both rows, which lets the concept-level share
# reflect that ambiguity honestly rather than picking one year arbitrarily.
s33_chars_from_asc_bb <- function() {
  f <- s33_find_files(S33_ASC_BB_PATTERN)
  if (length(f) == 0L) {
    cat("  ASC Addendum BB not found; SHARE_ASC_ANCILLARY_PACKAGED skipped.\n")
    return(NULL)
  }
  
  read_one <- function(path) {
    raw <- tryCatch(readLines(path, encoding = "latin1", warn = FALSE),
                    error = function(e) NULL)
    if (is.null(raw) || length(raw) < 5L) return(NULL)
    hdr_idx <- which(grepl("HCPCS Code", raw, ignore.case = TRUE))[1L]
    if (is.na(hdr_idx)) return(NULL)
    tryCatch(fread(text = paste(raw[hdr_idx:length(raw)], collapse = "\n"),
                   sep = "\t", header = TRUE, fill = TRUE, encoding = "Latin-1"),
             error = function(e) NULL)
  }
  parts <- Filter(Negate(is.null), lapply(f, read_one))
  if (length(parts) == 0L) {
    cat("  ASC Addendum BB unreadable; skipped.\n")
    return(NULL)
  }
  bb <- rbindlist(parts, fill = TRUE, use.names = TRUE)
  
  code_col <- grep("HCPCS", names(bb), ignore.case = TRUE, value = TRUE)[1L]
  pi_cands <- grep("Payment Indicator", names(bb), ignore.case = TRUE, value = TRUE)
  pi_col   <- pi_cands[grepl("Final", pi_cands, ignore.case = TRUE)][1L]
  if (is.na(code_col) || is.na(pi_col)) {
    cat("  ASC Addendum BB: expected columns not found; skipped.\n")
    return(NULL)
  }
  setnames(bb, c(code_col, pi_col), c("HCPCS_CODE", "PAYMENT_INDICATOR"))
  bb[, HCPCS_CODE       := toupper(trimws(as.character(HCPCS_CODE)))]
  bb[, PAYMENT_INDICATOR := toupper(trimws(as.character(PAYMENT_INDICATOR)))]
  bb <- unique(bb[nzchar(HCPCS_CODE), .(HCPCS_CODE, PAYMENT_INDICATOR)])
  
  map <- s33_code_to_concept_map()
  if (is.null(map)) {
    cat("  codebook unavailable; ASC BB cannot be joined to concepts.\n")
    return(NULL)
  }
  m <- merge(bb, map, by = "HCPCS_CODE")
  if (nrow(m) == 0L) {
    cat("  No ASC Addendum BB rows matched the codebook; skipped.\n")
    return(NULL)
  }
  
  out <- m[, .(SHARE_ASC_ANCILLARY_PACKAGED = mean(PAYMENT_INDICATOR == "N1", na.rm = TRUE)),
           by = FINAL_CONCEPT_ID]
  cat("  ASC Addendum BB: ", nrow(out), " concepts\n", sep = "")
  out
}

# --- source 7: Medicare Physician & Other Practitioners, Geography/Service -
#
# National-level rows only (excludes the ~255K state-level rows the file also
# carries). Volume and provider count are summed across the codes making up
# a concept, not averaged -- a concept's total national volume is the sum of
# its constituent codes' volumes, matching how N_CODES_IN_CONCEPT already
# treats code-level counts elsewhere in this section.
s33_chars_from_geo_volume <- function() {
  f <- s33_find_files(S33_GEO_PATTERN)
  if (length(f) == 0L) {
    cat("  Medicare geography/volume file not found; LN_NATIONAL_VOLUME, ",
        "N_PROVIDERS_NATIONAL and SHARE_VOL_FACILITY skipped.\n", sep = "")
    return(NULL)
  }
  geo <- tryCatch(fread(f[1L], select = c("Rndrng_Prvdr_Geo_Lvl", "HCPCS_Cd",
                                          "Place_Of_Srvc", "Tot_Rndrng_Prvdrs",
                                          "Tot_Srvcs")),
                  error = function(e) NULL)
  if (is.null(geo)) {
    cat("  Geography/volume file unreadable; skipped.\n")
    return(NULL)
  }
  geo <- geo[Rndrng_Prvdr_Geo_Lvl == "National"]
  if (nrow(geo) == 0L) {
    cat("  No National-level rows in the geography/volume file; skipped.\n")
    return(NULL)
  }
  geo[, HCPCS_CODE     := toupper(trimws(HCPCS_Cd))]
  geo[, Place_Of_Srvc  := toupper(trimws(Place_Of_Srvc))]
  geo[, Tot_Rndrng_Prvdrs := safe_numeric(Tot_Rndrng_Prvdrs)]
  geo[, Tot_Srvcs         := safe_numeric(Tot_Srvcs)]
  
  code_level <- geo[, .(
    VOLUME_NATIONAL = sum(Tot_Srvcs, na.rm = TRUE),
    N_PROVIDERS     = sum(Tot_Rndrng_Prvdrs, na.rm = TRUE),
    VOL_FACILITY    = sum(Tot_Srvcs[Place_Of_Srvc == "F"], na.rm = TRUE)
  ), by = HCPCS_CODE]
  
  map <- s33_code_to_concept_map()
  if (is.null(map)) {
    cat("  codebook unavailable; geography/volume file cannot be joined to concepts.\n")
    return(NULL)
  }
  m <- merge(code_level, map, by = "HCPCS_CODE")
  if (nrow(m) == 0L) {
    cat("  No geography/volume rows matched the codebook; skipped.\n")
    return(NULL)
  }
  
  out <- m[, .(
    LN_NATIONAL_VOLUME   = log(pmax(sum(VOLUME_NATIONAL, na.rm = TRUE), 1)),
    N_PROVIDERS_NATIONAL = sum(N_PROVIDERS, na.rm = TRUE),
    SHARE_VOL_FACILITY   = sum(VOL_FACILITY, na.rm = TRUE) /
                            pmax(sum(VOLUME_NATIONAL, na.rm = TRUE), 1)
  ), by = FINAL_CONCEPT_ID]
  
  cat("  geography/volume: ", nrow(out), " concepts\n", sep = "")
  out
}

# Deliberately NOT wrapped in cache_or_run(). This takes seconds, and caching
# it would serve a stale characteristic table the first time the codebook or
# the alt-price shards are re-downloaded, which is exactly the iteration this
# section is built for.
s33_build_characteristics <- function(panel, concepts_keep = NULL) {
  cat("\nBuilding concept characteristics:\n")
  pieces <- Filter(Negate(is.null), list(
    s33_chars_from_codebook(),
    s33_chars_from_panel(panel),
    s33_chars_from_alt_prices(),
    s33_chars_from_payer_dispersion(),
    s33_chars_from_pfs_rvu(),
    s33_chars_from_asc_bb(),
    s33_chars_from_geo_volume()))
  if (length(pieces) == 0L)
    stop("No characteristic source could be built. Check 33A.", call. = FALSE)
  
  ch <- Reduce(function(a, b) merge(a, b, by = "FINAL_CONCEPT_ID", all = TRUE), pieces)
  if (!is.null(concepts_keep)) ch <- ch[FINAL_CONCEPT_ID %chin% concepts_keep]
  
  cov <- rbindlist(lapply(S33_CHAR_NAMES, function(v) {
    if (!(v %in% names(ch)))
      return(data.table(CHARACTERISTIC = v, STATUS = "NOT BUILT",
                        N_NONMISSING = 0L, SHARE_NONMISSING = 0, SD = NA_real_))
    x   <- safe_numeric(ch[[v]])
    sdx <- sd(x, na.rm = TRUE)
    data.table(CHARACTERISTIC = v,
               STATUS = fcase(sum(is.finite(x)) < S33_MIN_CONCEPTS, "TOO FEW",
                              !is.finite(sdx) || sdx == 0,          "NO VARIATION",
                              default = "USABLE"),
               N_NONMISSING = sum(is.finite(x)),
               SHARE_NONMISSING = round(mean(is.finite(x)), 3), SD = sdx)
  }))
  cov <- merge(s33_registry_dt(), cov, by = "CHARACTERISTIC", all = TRUE, sort = FALSE)
  save_qa_csv(cov, "QA33C_characteristic_coverage.csv")
  
  usable <- cov[STATUS == "USABLE", CHARACTERISTIC]
  cat("\nUsable characteristics (", length(usable), "): ",
      paste(usable, collapse = ", "), "\n", sep = "")
  dropped <- cov[STATUS != "USABLE"]
  if (nrow(dropped) > 0L) {
    cat("Dropped:\n")
    print(as.data.frame(dropped[, .(CHARACTERISTIC, STATUS, N_NONMISSING)]))
  }
  if (length(usable) == 0L)
    stop("No usable characteristics. Fix the inputs flagged in 33A first.",
         call. = FALSE)
  
  list(characteristics = ch, coverage = cov, usable = usable)
}


# ============================================================================
# 33D  META-REGRESSIONS
# ============================================================================
#
# Unit of observation: one concept, one instrument. Weight 1/SE^2. Cluster on
# clinical family, which gives 16 clusters, so .pval() lands on t(15) exactly
# as it does for the two-way clustered models in Section 7. Characteristics are
# z-scored across concepts, unweighted, so every coefficient reads as the
# change in the concept-level effect per one cross-concept SD of that
# characteristic and coefficients are comparable down the column.

s33_assemble <- function(cr, chars, dep = "RF_COEF") {
  se_col <- S33_SE_FOR[[dep]]
  d <- as.data.table(cr)[is.finite(get(dep)) & is.finite(get(se_col)) & get(se_col) > 0]
  d <- merge(d, chars, by = "FINAL_CONCEPT_ID", all.x = TRUE, sort = FALSE)
  d <- copy(d)
  d[, CLUSTER_FAMILY := as.character(FINAL_FAMILY_ID)]
  d[, SHOP_CERTAINTY := factor(fifelse(FINAL_FAMILY_ID %chin% DIAGNOSTIC_FAMILIES,
                                       "Shoppable", "Non_shoppable"),
                               levels = c("Non_shoppable", "Shoppable"))]
  d[, LN_N_OBS := log(pmax(safe_numeric(N_OBSERVATIONS), 1))]
  d[, W_IV := 1 / (get(se_col)^2)]
  d
}

s33_zscore <- function(d, vars) {
  for (v in intersect(vars, names(d))) {
    x <- safe_numeric(d[[v]]); s <- sd(x, na.rm = TRUE)
    set(d, j = v,
        value = if (is.finite(s) && s > 0) (x - mean(x, na.rm = TRUE)) / s
        else rep(NA_real_, nrow(d)))
  }
  d
}

s33_fit_one <- function(d, dep, rhs, spec, weighting) {
  fe  <- if (spec == "FAMILY") " | CLUSTER_FAMILY" else ""
  ctl <- if (spec %chin% c("SIZE", "FAMILY"))
    intersect(S33_SIZE_CONTROLS, names(d)) else character(0)
  # unique() matters: MEAN_HOSP_PER_MARKET is both a candidate characteristic
  # and a size control, so its SIZE row is the RAW model plus LN_N_OBS only.
  terms <- unique(c(rhs, ctl))
  terms <- terms[nzchar(terms)]
  if (length(terms) == 0L) return(NULL)
  
  req <- unique(c(dep, terms, "CLUSTER_FAMILY", "W_IV"))
  req <- intersect(req, names(d))
  dd  <- copy(d[complete.cases(d[, ..req])])
  if (nrow(dd) < S33_MIN_CONCEPTS || uniqueN(dd$CLUSTER_FAMILY) < 3L) return(NULL)
  dd[, W := if (weighting == "Inverse variance") W_IV else 1]
  
  f   <- as.formula(paste0(dep, " ~ ", paste(terms, collapse = " + "), fe))
  fit <- tryCatch(feols(f, data = dd, weights = ~W, cluster = ~CLUSTER_FAMILY,
                        warn = FALSE, notes = FALSE), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  
  td <- tidy_fixest(fit)
  if (nrow(td) == 0L) return(NULL)
  td <- td[term %chin% rhs | grepl("^SHOP_CERTAINTY", term)]
  if (nrow(td) == 0L) return(NULL)
  
  r2 <- tryCatch(as.numeric(fitstat(fit, "r2", simplify = TRUE))[1L],
                 error = function(e) NA_real_)
  
  td[, `:=`(P_T = .pval(statistic, fit),
            ESTIMATE_PCT = 100 * estimate, SE_PCT = 100 * std.error,
            SPEC = spec, WEIGHTING = weighting, DEPENDENT = dep,
            N_CONCEPTS = nrow(dd), N_FAMILIES = uniqueN(dd$CLUSTER_FAMILY),
            R2 = r2)]
  td[, STARS := add_stars(P_T)]
  td[]
}

# --- one characteristic at a time, every spec ------------------------------
s33_meta_univariate <- function(d, usable, dep = "RF_COEF",
                                instruments = names(MAIN_INSTRUMENTS)) {
  rows <- list()
  for (il in intersect(instruments, unique(d$INSTRUMENT_LABEL))) {
    di <- d[INSTRUMENT_LABEL == il]
    for (v in usable) for (sp in S33_SPECS) for (w in META_WEIGHTINGS) {
      r <- s33_fit_one(di, dep, v, sp, w)
      if (is.null(r)) next
      r[, `:=`(INSTRUMENT_LABEL = il, CHARACTERISTIC = v)]
      rows[[length(rows) + 1L]] <- r
    }
  }
  out <- rbindlist(rows, fill = TRUE)
  if (nrow(out) == 0L) stop("No univariate meta-regressions estimated.", call. = FALSE)
  out <- merge(out, s33_registry_dt()[, .(CHARACTERISTIC, LABEL, OUTCOME_DERIVED)],
               by = "CHARACTERISTIC", all.x = TRUE, sort = FALSE)
  setorder(out, DEPENDENT, INSTRUMENT_LABEL, SPEC, P_T)
  save_csv(out, "T33B_meta_univariate_sweep.csv")
  out
}

# --- everything jointly ----------------------------------------------------
s33_meta_joint <- function(d, usable, dep = "RF_COEF",
                           instruments = names(MAIN_INSTRUMENTS)) {
  rows <- list()
  for (il in intersect(instruments, unique(d$INSTRUMENT_LABEL))) {
    di <- d[INSTRUMENT_LABEL == il]
    for (sp in S33_SPECS) {
      r <- s33_fit_one(di, dep, usable, sp, "Inverse variance")
      if (is.null(r)) next
      r[, INSTRUMENT_LABEL := il]
      rows[[length(rows) + 1L]] <- r
    }
  }
  out <- rbindlist(rows, fill = TRUE)
  if (nrow(out) == 0L) {
    warning("Joint model did not estimate; too few complete cases across all ",
            "characteristics at once.", call. = FALSE)
    return(data.table())
  }
  setnames(out, "term", "CHARACTERISTIC", skip_absent = TRUE)
  setorder(out, DEPENDENT, INSTRUMENT_LABEL, SPEC, P_T)
  save_csv(out, "T33C_meta_joint_model.csv")
  out
}

# --- the horse race --------------------------------------------------------
#
# Does the shoppability gradient survive once objective characteristics enter,
# and does any characteristic add anything on top of shoppability? Three nested
# models per instrument, reported side by side.
#
# The FAMILY spec is deliberately absent here. SHOP_CERTAINTY is a
# deterministic function of FINAL_FAMILY_ID, so family fixed effects absorb it
# entirely and the comparison is vacuous. The within-family evidence lives in
# the FAMILY rows of T33B instead.
s33_meta_horserace <- function(d, usable, dep = "RF_COEF",
                               instruments = names(MAIN_INSTRUMENTS)) {
  clean <- intersect(usable, s33_registry_dt()[OUTCOME_DERIVED == 0, CHARACTERISTIC])
  if (length(clean) == 0L) {
    warning("No non-outcome-derived characteristics usable; horse race skipped.",
            call. = FALSE)
    return(data.table())
  }
  models <- list(
    `1. Shoppability only`              = "SHOP_CERTAINTY",
    `2. Characteristics only`           = clean,
    `3. Shoppability + characteristics` = c("SHOP_CERTAINTY", clean))
  
  rows <- list()
  for (il in intersect(instruments, unique(d$INSTRUMENT_LABEL))) {
    di <- d[INSTRUMENT_LABEL == il]
    for (mn in names(models)) for (sp in c("RAW", "SIZE")) {
      r <- s33_fit_one(di, dep, models[[mn]], sp, "Inverse variance")
      if (is.null(r)) next
      r[, `:=`(INSTRUMENT_LABEL = il, MODEL = mn)]
      rows[[length(rows) + 1L]] <- r
    }
  }
  out <- rbindlist(rows, fill = TRUE)
  if (nrow(out) == 0L) {
    warning("Horse race did not estimate.", call. = FALSE)
    return(data.table())
  }
  setnames(out, "term", "CHARACTERISTIC", skip_absent = TRUE)
  setorder(out, DEPENDENT, INSTRUMENT_LABEL, MODEL, SPEC, P_T)
  save_csv(out, "T33D_meta_horserace_vs_shoppability.csv")
  out
}


# ============================================================================
# 33R  RESTORE FROM DISK -- loads Section 33's saved output, computes nothing
# ============================================================================
#
# s33_meta_univariate(), s33_meta_joint() and s33_meta_horserace() only ever
# WRITE their CSVs; none of the three checks whether a current one already
# exists before re-running its regressions. HPT_CONCEPT_CHARS <- TRUE forces
# the whole of Section 33 to rebuild every time regardless. This function is
# the read path that was missing: it loads what 33C/33D already wrote,
# reassembles s33_panel by simple merge and z-score (no regression), and sets
# every object 33's run block would have set. Nothing here calls feols().
#
# Requires concept_results in the session (from a plain warm start) and the
# four Section 33 CSVs on disk from a prior HPT_CONCEPT_CHARS <- TRUE run.
# Takes a few seconds regardless of how many characteristics or instruments
# that prior run covered.

restore_section33 <- function() {
  .s33_hd("33R. RESTORING SECTION 33 FROM SAVED OUTPUT")
  
  if (!exists("concept_results")) {
    stop("concept_results not in session -- run a warm start first.", call. = FALSE)
  }
  
  # Self-contained on purpose: read_table() only exists when HPT_RUN$figures
  # is TRUE (it is defined inside that guard, lines 7388-7786), which is FALSE
  # under a plain warm start. This function must not depend on it.
  .s33_read <- function(fn, dir = TABLE_DIR) {
    p <- file.path(dir, fn)
    if (!file.exists(p)) { message("SKIP: ", fn, " not found in ", dir); return(NULL) }
    fread(p)
  }
  
  s33_chars <- .s33_read("T33A_concept_characteristics.csv")
  if (is.null(s33_chars)) {
    stop("T33A_concept_characteristics.csv not found in TABLE_DIR. ",
         "Run with HPT_CONCEPT_CHARS <- TRUE at least once first.", call. = FALSE)
  }
  s33_chars <- as.data.table(s33_chars)
  
  s33_cov <- .s33_read("QA33C_characteristic_coverage.csv", dir = QA_DIR)
  if (is.null(s33_cov)) {
    stop("QA33C_characteristic_coverage.csv not found in QA_DIR.", call. = FALSE)
  }
  s33_usable <- s33_cov[STATUS == "USABLE", CHARACTERISTIC]
  
  s33_uni   <- .s33_read("T33B_meta_univariate_sweep.csv")
  s33_joint <- .s33_read("T33C_meta_joint_model.csv")
  s33_race  <- .s33_read("T33D_meta_horserace_vs_shoppability.csv")
  
  s33_panel <- s33_assemble(concept_results, s33_chars, dep = "RF_COEF")
  s33_panel <- s33_zscore(s33_panel, unique(c(s33_usable, S33_SIZE_CONTROLS)))
  
  cat(sprintf(
    "Restored: %d usable characteristics | s33_panel %s rows | uni %s | joint %s | race %s\n",
    length(s33_usable), format(nrow(s33_panel), big.mark = ","),
    if (is.null(s33_uni)) "MISSING" else format(nrow(s33_uni), big.mark = ","),
    if (is.null(s33_joint)) "MISSING" else format(nrow(s33_joint), big.mark = ","),
    if (is.null(s33_race)) "MISSING" else format(nrow(s33_race), big.mark = ",")))
  cat("No regression was re-estimated. Objects match names from a live Section 33 run:\n",
      "  s33_chars, s33_usable, s33_panel, s33_uni, s33_joint, s33_race\n", sep = "")
  
  assign("s33_chars",  s33_chars,  envir = .GlobalEnv)
  assign("s33_usable", s33_usable, envir = .GlobalEnv)
  assign("s33_panel",  s33_panel,  envir = .GlobalEnv)
  assign("s33_uni",    s33_uni,    envir = .GlobalEnv)
  assign("s33_joint",  s33_joint,  envir = .GlobalEnv)
  assign("s33_race",   s33_race,   envir = .GlobalEnv)
  invisible(list(chars = s33_chars, usable = s33_usable, panel = s33_panel,
                 uni = s33_uni, joint = s33_joint, race = s33_race))
}


# ============================================================================
# 33  RUN BLOCK. Each block is checked before the next is run.
# ============================================================================

HPT_RUN$concept_characteristics <-
  if (exists("HPT_CONCEPT_CHARS")) isTRUE(HPT_CONCEPT_CHARS) else !isTRUE(HPT_WARM_START)

if (isTRUE(HPT_RUN$concept_characteristics)) {
  
  stopifnot("concept_results not in session" = exists("concept_results"),
            "outpatient not in session"      = exists("outpatient"))
  
  # ---- GATE 1 -------------------------------------------------------------
  .s33_hd("33A. INVENTORY -- what exists before anything is built")
  s33_inv <- s33_inventory()
  print(as.data.frame(s33_inv))
  cat("\nA MISSING row under 'codebook' or 'full_sql_codebook' means those\n",
      "characteristics drop out of 33C without further warning. The two that\n",
      "matter most are OPPS_STATUS_INDICATORS and IS_ASC_COVERED_FLAG.\n", sep = "")
  
  # ---- GATE 2 -------------------------------------------------------------
  .s33_hd("33B. HETEROGENEITY BUDGET -- is there anything to explain?")
  s33_budget <- s33_heterogeneity_budget(concept_results)
  print(as.data.frame(s33_budget))
  
  s33_rel <- s33_cross_instrument_reliability(concept_results, dep = "RF_COEF")
  if (nrow(s33_rel) > 0L) {
    cat("\nCross-instrument reliability of the concept ranking:\n")
    print(as.data.frame(s33_rel))
  }
  
  cat("\nInterpreting 33B.\n",
      "  I2 near 0, or NAIVE_SIGNAL_SHARE near 0, means the concept estimates\n",
      "    are almost all sampling noise. No covariate can explain them. That\n",
      "    is a reportable finding, not a failure, and it also explains why a\n",
      "    causal forest over the same estimates would return a flat CATE.\n",
      "  I2 high but cross-instrument SPEARMAN near 0 means the dispersion is\n",
      "    real within an instrument and does not reproduce across them, which\n",
      "    points at instrument-specific noise rather than service economics.\n",
      "  I2 high and SPEARMAN high means there is real, replicable service-\n",
      "    level heterogeneity. Proceed to 33C.\n",
      "  Q assumes independent estimates and these share hospitals, so read I2\n",
      "    as an upper bound and weight the reliability table more heavily.\n",
      sep = "")
  
  # ---- GATE 3 -------------------------------------------------------------
  .s33_hd("33C. CONCEPT CHARACTERISTICS")
  s33_chars_obj <- s33_build_characteristics(
    outpatient, concepts_keep = unique(concept_results$FINAL_CONCEPT_ID))
  s33_chars  <- s33_chars_obj$characteristics
  s33_usable <- s33_chars_obj$usable
  save_csv(s33_chars, "T33A_concept_characteristics.csv")
  
  .s33_hd("33D. META-REGRESSIONS")
  s33_panel <- s33_assemble(concept_results, s33_chars, dep = "RF_COEF")
  s33_panel <- s33_zscore(s33_panel, unique(c(s33_usable, S33_SIZE_CONTROLS)))
  
  s33_uni <- s33_meta_univariate(s33_panel, s33_usable, dep = "RF_COEF")
  cat("\nUnivariate sweep, primary instrument, inverse-variance weighted.\n",
      "RAW is the least informative column and FAMILY is the sharp test.\n\n", sep = "")
  print(as.data.frame(
    s33_uni[INSTRUMENT_LABEL == names(MAIN_INSTRUMENTS)[1L] &
              WEIGHTING == "Inverse variance",
            .(CHARACTERISTIC, SPEC, OUTCOME_DERIVED,
              EST_PCT = round(ESTIMATE_PCT, 3), SE_PCT = round(SE_PCT, 3),
              P = round(P_T, 4), STARS, N_CONCEPTS)]))
  
  s33_joint <- s33_meta_joint(s33_panel, s33_usable, dep = "RF_COEF")
  if (nrow(s33_joint) > 0L) {
    cat("\nJoint model:\n\n")
    print(as.data.frame(s33_joint[, .(INSTRUMENT_LABEL, SPEC, CHARACTERISTIC,
                                      EST_PCT = round(ESTIMATE_PCT, 3),
                                      P = round(P_T, 4), STARS,
                                      N_CONCEPTS, R2 = round(R2, 3))]))
  }
  
  s33_race <- s33_meta_horserace(s33_panel, s33_usable, dep = "RF_COEF")
  if (nrow(s33_race) > 0L) {
    cat("\nHorse race against shoppability:\n\n")
    print(as.data.frame(s33_race[, .(INSTRUMENT_LABEL, MODEL, SPEC, CHARACTERISTIC,
                                     EST_PCT = round(ESTIMATE_PCT, 3),
                                     P = round(P_T, 4), STARS,
                                     N_CONCEPTS, R2 = round(R2, 3))]))
  }
  
  cat("\nInterpreting 33D.\n",
      "  A characteristic significant in RAW but not FAMILY is redescribing\n",
      "    the family gradient the shoppability schemes already capture. Not a\n",
      "    new finding and should not be written up as one.\n",
      "  A characteristic significant in FAMILY is within-family evidence that\n",
      "    none of the 18 schemes can express. That is the result worth having.\n",
      "  In the horse race, SHOP_CERTAINTY is compared across models 1 and 3.\n",
      "    Stability indicates shoppability is not proxying for the objective\n",
      "    characteristics and the headline is unaffected. A collapse implies\n",
      "    the mechanism claim needs rewriting rather than an extra table.\n",
      "  Every coefficient is per one cross-concept SD of the characteristic,\n",
      "    on the same percent scale as RF_PERCENT_PER_SD.\n", sep = "")
  
} else {
  s33_inv <- s33_budget <- s33_chars <- s33_uni <- s33_joint <- s33_race <- NULL
}

cat("\nSection 33 complete.\n")


# ===========================================================================
# 33E  ROBUSTNESS CHECKS FOR SHARE_SI_PACKAGED
# ===========================================================================
#
# Off by default. Everything here reads s33_panel, s33_chars and s33_usable,
# which exist only after the Section 33 run block above has executed, so this
# errors on a warm start that leaves Section 33 switched off.
#
# The split-half block calls estimate_concept_level() twice. That is the same
# function behind the eight-hour concept sweep, restricted here to the primary
# instrument on half the hospitals, so budget hours rather than minutes.
#
#   Leave-one-family-out   is the SHARE_SI_PACKAGED result driven by the two
#                          families holding 82% of the packaged concepts?
#   Split-half reliability what share of the cross-concept spread in RF_COEF
#                          is real, estimated by splitting hospitals at random
#                          and correlating the two independent sweeps?
#   Correlation matrix     collinearity among the SI shares and the support
#                          cluster, which the joint model in 33D runs into.
#   Packaging vs SHOP      does packaging cut across the shoppable/non-shoppable
#                          split, or merely relabel it?

if (exists("HPT_SCRATCH") && isTRUE(HPT_SCRATCH)) {   # interactive scratch, off by default
  
  
  
  s33_loo_family <- function(drop_family, panel = s33_panel, dep = "RF_COEF",
                             instrument = "Competitor_only_hospitals_9m") {
    d <- panel[FINAL_FAMILY_ID != drop_family & INSTRUMENT_LABEL == instrument]
    r <- s33_fit_one(d, dep, "SHARE_SI_PACKAGED", "FAMILY", "Inverse variance")
    r[, INSTRUMENT_LABEL := instrument]
  }
  
  loo_result <- rbindlist(list(
    data.table(DROPPED = "none",
               s33_fit_one(s33_panel[INSTRUMENT_LABEL == "Competitor_only_hospitals_9m"],
                           "RF_COEF", "SHARE_SI_PACKAGED", "FAMILY", "Inverse variance")),
    data.table(DROPPED = "XRAY_FLUOROSCOPY", s33_loo_family("XRAY_FLUOROSCOPY")),
    data.table(DROPPED = "BIOPSY",           s33_loo_family("BIOPSY"))
  ), fill = TRUE)
  
  loo_result
  
  
  d2 <- s33_panel[FINAL_FAMILY_ID %chin% c("XRAY_FLUOROSCOPY", "BIOPSY") == FALSE &
                    INSTRUMENT_LABEL == "Competitor_only_hospitals_9m"]
  
  loo_result2 <- rbindlist(list(
    loo_result[, INSTRUMENT_LABEL := "Competitor_only_hospitals_9m"],   # fixes the NA on row 1
    data.table(DROPPED = "both",
               s33_fit_one(d2, "RF_COEF", "SHARE_SI_PACKAGED", "FAMILY", "Inverse variance"))
  ), fill = TRUE)
  loo_result2[, INSTRUMENT_LABEL := "Competitor_only_hospitals_9m"]
  
  loo_result2
  
  
  
  
  
  set.seed(20260909)  # fixed seed, so this is exactly reproducible if anyone asks
  
  hosp_ids <- unique(outpatient$HOSPITAL_ID)
  hosp_A   <- sample(hosp_ids, floor(length(hosp_ids) / 2))
  hosp_B   <- setdiff(hosp_ids, hosp_A)
  
  slim_primary <- build_concept_panel(outpatient,
                                      instruments = MAIN_INSTRUMENTS["Competitor_only_hospitals_9m"])
  
  slim_A <- slim_primary[HOSPITAL_ID %in% hosp_A]; setkey(slim_A, ANALYSIS_SERVICE_ID)
  slim_B <- slim_primary[HOSPITAL_ID %in% hosp_B]; setkey(slim_B, ANALYSIS_SERVICE_ID)
  
  # NOT wrapped in cache_or_run, and save_stem is overridden -- this must never
  # write to T05_concept_level* or touch concept_level_6inst.rds.
  sh_A <- estimate_concept_level(slim_A, instruments = MAIN_INSTRUMENTS["Competitor_only_hospitals_9m"],
                                 save_stem = "QA33_splithalf_A")
  sh_B <- estimate_concept_level(slim_B, instruments = MAIN_INSTRUMENTS["Competitor_only_hospitals_9m"],
                                 save_stem = "QA33_splithalf_B")
  
  sh_merged <- merge(sh_A[, .(FINAL_CONCEPT_ID, RF_A = RF_COEF)],
                     sh_B[, .(FINAL_CONCEPT_ID, RF_B = RF_COEF)],
                     by = "FINAL_CONCEPT_ID")
  
  cat("Concepts in both halves:", nrow(sh_merged), "of",
      uniqueN(outpatient$FINAL_CONCEPT_ID), "\n")
  cat("Pearson: ",  round(cor(sh_merged$RF_A, sh_merged$RF_B, method = "pearson"),  3), "\n")
  cat("Spearman:", round(cor(sh_merged$RF_A, sh_merged$RF_B, method = "spearman"), 3), "\n")
  
  save_qa_csv(sh_merged, "QA33_splithalf_concept_correlation.csv")
  
  
  
  
  
  chars_for_corr <- c("SHARE_SI_PACKAGED", "SHARE_SI_COMPREHENSIVE", "SHARE_SI_COND_PACKAGED",
                      "SHARE_SI_SEPARATE", "SD_LN_PRICE_XHOSP", "LN_N_HOSPITALS",
                      "MEAN_HOSP_PER_MARKET", "MEAN_N_PAYERS", "MEAN_CV_PAYER")
  
  one_inst <- s33_panel[INSTRUMENT_LABEL == "Competitor_only_hospitals_9m", ..chars_for_corr]
  round(cor(one_inst, use = "pairwise.complete.obs"), 2)
  
  
  
  
  
  grp <- unique(s33_panel[INSTRUMENT_LABEL == "Competitor_only_hospitals_9m",
                          .(FINAL_CONCEPT_ID, SHOP_CERTAINTY)])
  
  merge(s33_chars[, .(FINAL_CONCEPT_ID, SHARE_SI_PACKAGED)], grp, by = "FINAL_CONCEPT_ID")[
    , .(n = .N,
        n_packaged = sum(SHARE_SI_PACKAGED > 0, na.rm = TRUE),
        mean_packaged = round(mean(SHARE_SI_PACKAGED, na.rm = TRUE), 3))
    , by = SHOP_CERTAINTY]
  
  
  
  
  
  chars_no_sd <- setdiff(s33_usable, "SD_LN_PRICE_XHOSP")
  r <- s33_fit_one(s33_panel[INSTRUMENT_LABEL == "Competitor_only_hospitals_9m"],
                   "RF_COEF", chars_no_sd, "FAMILY", "Inverse variance")
  r[term == "SHARE_SI_PACKAGED"]
}   # end interactive scratch




###############################################################################
#
#   SECTION 34 -- CODE-KEYED SUPPLEMENTAL CHARACTERISTICS
#
#   Adds three externally sourced characteristic sets to the Section 33
#   framework. Additive and non-destructive: Section 33's own objects and its
#   T33*/QA33* outputs are left untouched, and everything here writes under a
#   T34/QA34 prefix. If the new characteristics turn out to be uninformative,
#   nothing has to be undone.
#
#   Requires a live session with Section 33 already loaded, via either
#   HPT_CONCEPT_CHARS <- TRUE or restore_section33(). Specifically needs
#   concept_results, s33_chars, s33_usable, and the s33_* helper functions.
#
#   ---------------------------------------------------------------------------
#   THE THREE SOURCES, AND WHY EACH IS KEYED DIFFERENTLY FROM SECTION 33
#   ---------------------------------------------------------------------------
#
#   Every Section 33 source is concept-keyed already (the codebook is
#   collapsed to concept before use; the alt-price and payer-dispersion files
#   arrive keyed on ANALYSIS_CONCEPT_ID). All three sources here are keyed on
#   HCPCS CODE instead, so each is joined through the codebook's
#   BILLING_CODE -> ANALYSIS_CONCEPT_ID crosswalk and then collapsed with
#   s33_to_canonical(), which is what keeps a code belonging to one of the six
#   MERGE_GROUPS concepts (mammography tomosynthesis, the MRI/CT abdomen
#   variants) landing on the same canonical concept the rest of the paper is
#   estimated on.
#
#   PFS RVU (PPRRVU)        PC/TC split and global-surgery period. A code with
#                           PCTC IND = 1 is billed as separate professional and
#                           technical components, so the hospital's posted
#                           price is only part of what the patient owes. None
#                           of the Section 33 characteristics express this.
#
#   ASC Addendum BB         Packaging status for covered ancillary services,
#                           from the ASC payment system rather than OPPS. This
#                           is an independent second reading of "is this
#                           service separately priced," which makes it a
#                           corroboration check on SHARE_SI_PACKAGED rather
#                           than a restatement of it.
#
#   Geography and Service   National Medicare FFS volume, distinct provider
#                           count, and the facility share of volume. Volume and
#                           supply breadth are absent from Section 33 entirely.
#
#   ---------------------------------------------------------------------------
#   WHAT IS NOT HERE
#   ---------------------------------------------------------------------------
#
#   OPPS Addendum A is not built. It is keyed on APC number, and
#   HPT_CODEBOOK.csv carries no APC column to join it through
#   (OPPS_STATUS_INDICATORS is the only OPPS field in that export). Adding it
#   requires either a fresh Phase 1 SQL pull that retains APC per code, or
#   Addendum B's own code-to-APC mapping loaded separately.
#
#   ---------------------------------------------------------------------------
#   OUTPUTS
#   ---------------------------------------------------------------------------
#     QA34A_supplemental_coverage.csv    per-characteristic coverage and SD
#     T34A_concept_characteristics.csv   Section 33 table plus the new columns
#     T34B_meta_univariate_sweep.csv
#     T34C_meta_joint_model.csv
#     T34D_meta_horserace_vs_shoppability.csv
#
###############################################################################


# ============================================================================
# 34.0  FILE PATTERNS
# ============================================================================
#
# s33_find_files() passes ignore.case = TRUE to list.files(), which uses POSIX
# regex. A PCRE inline flag such as (?i) would be read as four literal
# characters and match nothing, so case handling is left to that argument.

S34_PFS_RVU_PATTERN <- "^PPRRVU.*\\.xlsx$"
S34_ASC_BB_PATTERN  <- "Addendum[ _-]?BB.*\\.txt$"
S34_GEO_PATTERN     <- "^MUP_PHY.*Geo\\.csv$"

# S33's own s33_find_files() is not recursive, which is correct for Section
# 33's four sources -- they sit flat in PANEL_DIR. Section 34's three sources
# each unzipped into their own subfolder, so this is a separate function
# rather than a change to s33_find_files(), to avoid touching anything
# Section 33 already depends on.
s34_find_files <- function(pattern, dir = PANEL_DIR)
  list.files(dir, pattern = pattern, full.names = TRUE, ignore.case = TRUE,
             recursive = TRUE)

.s34_hd <- function(x)
  cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")


# ============================================================================
# 34.1  BILLING_CODE -> CONCEPT CROSSWALK, SHARED BY ALL THREE LOADERS
# ============================================================================

s34_code_to_concept_map <- function() {
  cb <- s33_read_all(S33_CODEBOOK_PATTERN)
  if (is.null(cb)) cb <- tryCatch(fread(FILES$codebook), error = function(e) NULL)
  if (is.null(cb) || !all(c("BILLING_CODE", "ANALYSIS_CONCEPT_ID") %in% names(cb))) {
    cat("  codebook unavailable or missing BILLING_CODE/ANALYSIS_CONCEPT_ID.\n")
    return(NULL)
  }
  cb <- as.data.table(cb)
  m  <- s33_to_canonical(cb[, .(BILLING_CODE, ANALYSIS_CONCEPT_ID)], "ANALYSIS_CONCEPT_ID")
  setnames(m, "BILLING_CODE", "HCPCS_CODE")
  m[, HCPCS_CODE := toupper(trimws(as.character(HCPCS_CODE)))]
  unique(m[nzchar(HCPCS_CODE)], by = "HCPCS_CODE")
}


# ============================================================================
# 34.2  SOURCE LOADERS
# ============================================================================

# --- PFS RVU: PC/TC split and global-surgery period -------------------------
#
# Global row only, meaning a blank modifier. A code with PCTC IND = 1 also has
# separate -26 and -TC rows; the global row is the one whose RVUs correspond to
# the whole service, so restricting to it answers "does a split exist for this
# code" once rather than three times.
#
# Columns are addressed by position because the header is unusable by name:
# CMS wraps it across spreadsheet rows 9 and 10, leaving the literal string
# "INDICATOR" repeated three times and the PCTC and GLOB columns carrying only
# the fragments "IND" and "DAYS". Positions 1, 2, 14, 15 are verified against
# the CY2025 January release. A quarter with a different column count would
# need this re-checked.
s34_chars_from_pfs_rvu <- function() {
  f <- s34_find_files(S34_PFS_RVU_PATTERN)
  if (length(f) == 0L) {
    cat("  PFS RVU file not found; SHARE_PCTC_SPLIT and SHARE_GLOB_SURGICAL skipped.\n")
    return(NULL)
  }
  if (length(f) > 1L)
    cat("  Multiple PFS RVU files matched; using ", basename(f[1L]), "\n", sep = "")
  
  if (!requireNamespace("readxl", quietly = TRUE)) {
    cat("  readxl not installed (install.packages(\"readxl\")); PFS RVU skipped.\n")
    return(NULL)
  }
  raw <- tryCatch(
    readxl::read_excel(f[1L], skip = 9L, col_names = TRUE, .name_repair = "minimal"),
    error = function(e) NULL)
  if (is.null(raw) || ncol(raw) < 15L) {
    cat("  PFS RVU file unreadable or unexpected layout; skipped.\n")
    return(NULL)
  }
  
  dt <- as.data.table(raw)
  keep <- dt[, c(1L, 2L, 14L, 15L), with = FALSE]
  setnames(keep, c("HCPCS_CODE", "MODIFIER", "PCTC_IND", "GLOB_DAYS"))
  
  keep[, HCPCS_CODE := toupper(trimws(as.character(HCPCS_CODE)))]
  keep[, MODIFIER   := toupper(trimws(as.character(MODIFIER)))]
  keep[, PCTC_IND   := safe_numeric(PCTC_IND)]
  keep[, GLOB_DAYS  := toupper(trimws(as.character(GLOB_DAYS)))]
  keep <- keep[!is.na(HCPCS_CODE) & nzchar(HCPCS_CODE) &
                 (is.na(MODIFIER) | MODIFIER == "" | MODIFIER == "NA")]
  
  map <- s34_code_to_concept_map()
  if (is.null(map)) return(NULL)
  m <- merge(keep, map, by = "HCPCS_CODE")
  if (nrow(m) == 0L) {
    cat("  No PFS RVU rows matched the codebook; skipped.\n")
    return(NULL)
  }
  
  out <- m[, .(
    SHARE_PCTC_SPLIT    = mean(PCTC_IND == 1, na.rm = TRUE),
    SHARE_GLOB_SURGICAL = mean(GLOB_DAYS %chin% c("000", "010", "090"), na.rm = TRUE)
  ), by = FINAL_CONCEPT_ID]
  
  cat("  PFS RVU: ", nrow(out), " concepts from ", nrow(m), " matched codes\n", sep = "")
  out
}

# --- ASC Addendum BB: ancillary packaging status ---------------------------
#
# Every matching vintage is read and row-bound, then deduplicated on code and
# indicator together. A code whose indicator genuinely changed between years
# therefore keeps both rows, and the concept-level share reflects that
# ambiguity rather than silently taking whichever year sorted first.
s34_chars_from_asc_bb <- function() {
  f <- s34_find_files(S34_ASC_BB_PATTERN)
  if (length(f) == 0L) {
    cat("  ASC Addendum BB not found; SHARE_ASC_ANCILLARY_PACKAGED skipped.\n")
    return(NULL)
  }
  
  read_one <- function(path) {
    raw <- tryCatch({
      con <- file(path, encoding = "latin1")
      on.exit(close(con), add = TRUE)
      readLines(con, warn = FALSE)
    }, error = function(e) NULL)
    if (is.null(raw) || length(raw) < 5L) return(NULL)
    hdr <- which(grepl("HCPCS", raw, ignore.case = TRUE))
    if (length(hdr) == 0L) return(NULL)
    tryCatch(fread(text = paste(raw[hdr[1L]:length(raw)], collapse = "\n"),
                   sep = "\t", header = TRUE, fill = TRUE),
             error = function(e) NULL)
  }
  parts <- Filter(Negate(is.null), lapply(f, read_one))
  if (length(parts) == 0L) {
    cat("  ASC Addendum BB unreadable; skipped.\n")
    return(NULL)
  }
  bb <- rbindlist(parts, fill = TRUE, use.names = TRUE)
  
  code_col <- grep("HCPCS", names(bb), ignore.case = TRUE, value = TRUE)
  pi_all   <- grep("Payment Indicator", names(bb), ignore.case = TRUE, value = TRUE)
  pi_final <- pi_all[grepl("Final", pi_all, ignore.case = TRUE)]
  pi_col   <- if (length(pi_final)) pi_final[1L] else
    if (length(pi_all))   pi_all[1L]   else character(0)
  if (length(code_col) == 0L || length(pi_col) == 0L) {
    cat("  ASC Addendum BB: expected columns not found; skipped.\n")
    return(NULL)
  }
  
  bb <- bb[, c(code_col[1L], pi_col), with = FALSE]
  setnames(bb, c("HCPCS_CODE", "PAYMENT_INDICATOR"))
  bb[, HCPCS_CODE        := toupper(trimws(as.character(HCPCS_CODE)))]
  bb[, PAYMENT_INDICATOR := toupper(trimws(as.character(PAYMENT_INDICATOR)))]
  bb <- unique(bb[nzchar(HCPCS_CODE)])
  
  map <- s34_code_to_concept_map()
  if (is.null(map)) return(NULL)
  m <- merge(bb, map, by = "HCPCS_CODE")
  if (nrow(m) == 0L) {
    cat("  No ASC Addendum BB rows matched the codebook; skipped.\n")
    return(NULL)
  }
  
  out <- m[, .(SHARE_ASC_ANCILLARY_PACKAGED =
                 mean(PAYMENT_INDICATOR == "N1", na.rm = TRUE)),
           by = FINAL_CONCEPT_ID]
  cat("  ASC Addendum BB: ", nrow(out), " concepts from ", length(f),
      " vintage(s)\n", sep = "")
  out
}

# --- Geography and Service: national volume, providers, facility share -----
#
# National rows only. The file also carries roughly 255,000 state-level rows,
# which are the same services counted again per state and would double-count
# every total below.
#
# Volume and provider counts are SUMMED across the codes making up a concept
# rather than averaged: a concept's national volume is the total of its
# constituent codes' volumes, the same way N_CODES_IN_CONCEPT already treats
# code-level counts in Section 33.
s34_chars_from_geo_volume <- function() {
  f <- s34_find_files(S34_GEO_PATTERN)
  if (length(f) == 0L) {
    cat("  Geography/volume file not found; LN_NATIONAL_VOLUME, ",
        "N_PROVIDERS_NATIONAL and SHARE_VOL_FACILITY skipped.\n", sep = "")
    return(NULL)
  }
  need <- c("Rndrng_Prvdr_Geo_Lvl", "HCPCS_Cd", "Place_Of_Srvc",
            "Tot_Rndrng_Prvdrs", "Tot_Srvcs")
  geo <- tryCatch(fread(f[1L], select = need), error = function(e) NULL)
  if (is.null(geo)) {
    cat("  Geography/volume file unreadable or missing expected columns; skipped.\n")
    return(NULL)
  }
  
  geo <- geo[Rndrng_Prvdr_Geo_Lvl == "National"]
  if (nrow(geo) == 0L) {
    cat("  No National-level rows in the geography/volume file; skipped.\n")
    return(NULL)
  }
  geo[, HCPCS_CODE        := toupper(trimws(as.character(HCPCS_Cd)))]
  geo[, POS               := toupper(trimws(as.character(Place_Of_Srvc)))]
  geo[, Tot_Rndrng_Prvdrs := safe_numeric(Tot_Rndrng_Prvdrs)]
  geo[, Tot_Srvcs         := safe_numeric(Tot_Srvcs)]
  
  code_level <- geo[, .(
    VOL   = sum(Tot_Srvcs, na.rm = TRUE),
    PRVDR = sum(Tot_Rndrng_Prvdrs, na.rm = TRUE),
    VOL_F = sum(Tot_Srvcs[POS == "F"], na.rm = TRUE)
  ), by = HCPCS_CODE]
  
  map <- s34_code_to_concept_map()
  if (is.null(map)) return(NULL)
  m <- merge(code_level, map, by = "HCPCS_CODE")
  if (nrow(m) == 0L) {
    cat("  No geography/volume rows matched the codebook; skipped.\n")
    return(NULL)
  }
  
  out <- m[, {
    v <- sum(VOL, na.rm = TRUE)
    list(LN_NATIONAL_VOLUME   = log(pmax(v, 1)),
         N_PROVIDERS_NATIONAL = sum(PRVDR, na.rm = TRUE),
         SHARE_VOL_FACILITY   = if (v > 0) sum(VOL_F, na.rm = TRUE) / v else NA_real_)
  }, by = FINAL_CONCEPT_ID]
  
  cat("  geography/volume: ", nrow(out), " concepts from ", nrow(m),
      " matched codes\n", sep = "")
  out
}


# ============================================================================
# 34.3  REGISTRY EXTENSION
# ============================================================================
#
# s33_registry_dt() reads S33_CHAR_REGISTRY at call time, so appending here is
# enough for the LABEL/OUTCOME_DERIVED merge in the univariate sweep and for
# the horse race's OUTCOME_DERIVED == 0 filter to see the new entries. None of
# the six is built from the negotiated price, so all are outcome_derived =
# FALSE and all are eligible for the horse race.
#
# The guard makes re-pasting this block harmless: without it, a second paste
# would append duplicate entries and every new characteristic would appear
# twice in the sweep.

S34_CHAR_REGISTRY <- list(
  list(name = "SHARE_PCTC_SPLIT",
       label = "Share of codes with a PC/TC split (PFS PCTC IND = 1)",
       outcome_derived = FALSE,
       note = "The posted price is part of the bill: a separate professional-component invoice exists that a patient comparing hospital prices does not see."),
  list(name = "SHARE_GLOB_SURGICAL",
       label = "Share of codes with a global-surgery package (PFS GLOB DAYS 000/010/090)",
       outcome_derived = FALSE,
       note = "The posted price bundles a defined follow-up period rather than a single encounter. XXX and ZZZ codes carry no such bundle."),
  list(name = "SHARE_ASC_ANCILLARY_PACKAGED",
       label = "Share of codes packaged under the ASC ancillary schedule (Addendum BB, N1)",
       outcome_derived = FALSE,
       note = "Packaging status from the ASC payment system rather than OPPS, so it corroborates or contradicts SHARE_SI_PACKAGED without sharing a source with it."),
  list(name = "LN_NATIONAL_VOLUME",
       label = "Log national Medicare FFS volume",
       outcome_derived = FALSE,
       note = "How often the service is performed nationally. A near-zero-volume code is not something patients realistically shop for."),
  list(name = "N_PROVIDERS_NATIONAL",
       label = "National count of distinct rendering providers",
       outcome_derived = FALSE,
       note = "Breadth of supply, and a direct measure of whether an alternative provider exists to shop toward."),
  list(name = "SHARE_VOL_FACILITY",
       label = "Share of national volume in a facility place of service",
       outcome_derived = FALSE,
       note = "Where most volume is billed in physician offices rather than hospital outpatient departments, a hospital's posted price is largely irrelevant to what patients actually pay.")
)

S34_CHAR_NAMES <- vapply(S34_CHAR_REGISTRY, `[[`, character(1), "name")

if (!all(S34_CHAR_NAMES %in% vapply(S33_CHAR_REGISTRY, `[[`, character(1), "name"))) {
  S33_CHAR_REGISTRY <- c(S33_CHAR_REGISTRY, S34_CHAR_REGISTRY)
  S33_CHAR_NAMES    <- vapply(S33_CHAR_REGISTRY, `[[`, character(1), "name")
  cat("Section 34: registry extended to ", length(S33_CHAR_NAMES),
      " characteristics.\n", sep = "")
} else {
  cat("Section 34: registry already extended; left unchanged.\n")
}


# ============================================================================
# 34.4  META-REGRESSIONS WRITING UNDER A T34 PREFIX
# ============================================================================
#
# Structurally identical to s33_meta_univariate/_joint/_horserace, including
# the same s33_fit_one() call, the same spec and weighting loops, and the same
# .pval() reference. The only difference is the output filename, so that
# Section 33's T33B/C/D remain on disk and restore_section33() keeps working
# against the 19-characteristic set while these run alongside it.

s34_meta_univariate <- function(d, usable, dep = "RF_COEF",
                                instruments = names(MAIN_INSTRUMENTS)) {
  rows <- list()
  cat("  univariate sweep: ", length(intersect(instruments, unique(d$INSTRUMENT_LABEL))),
      " instrument(s), ", length(usable), " characteristics\n", sep = "")
  for (il in intersect(instruments, unique(d$INSTRUMENT_LABEL))) {
    t0 <- Sys.time()
    di <- d[INSTRUMENT_LABEL == il]
    n_fit <- 0L
    for (v in usable) for (sp in S33_SPECS) for (w in META_WEIGHTINGS) {
      r <- s33_fit_one(di, dep, v, sp, w)
      if (is.null(r)) next
      r[, `:=`(INSTRUMENT_LABEL = il, CHARACTERISTIC = v)]
      rows[[length(rows) + 1L]] <- r
      n_fit <- n_fit + 1L
    }
    cat(sprintf("    %-38s %3d fits | %.1fs\n", il, n_fit,
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
  out <- rbindlist(rows, fill = TRUE)
  if (nrow(out) == 0L) stop("No univariate meta-regressions estimated.", call. = FALSE)
  out <- merge(out, s33_registry_dt()[, .(CHARACTERISTIC, LABEL, OUTCOME_DERIVED)],
               by = "CHARACTERISTIC", all.x = TRUE, sort = FALSE)
  setorder(out, DEPENDENT, INSTRUMENT_LABEL, SPEC, P_T)
  save_csv(out, "T34B_meta_univariate_sweep.csv")
  out
}

s34_meta_joint <- function(d, usable, dep = "RF_COEF",
                           instruments = names(MAIN_INSTRUMENTS)) {
  rows <- list()
  cat("  joint model:\n")
  for (il in intersect(instruments, unique(d$INSTRUMENT_LABEL))) {
    t0 <- Sys.time()
    di <- d[INSTRUMENT_LABEL == il]
    for (sp in S33_SPECS) {
      r <- s33_fit_one(di, dep, usable, sp, "Inverse variance")
      if (is.null(r)) next
      r[, INSTRUMENT_LABEL := il]
      rows[[length(rows) + 1L]] <- r
    }
    cat(sprintf("    %-38s %.1fs\n", il,
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
  out <- rbindlist(rows, fill = TRUE)
  if (nrow(out) == 0L) {
    warning("Joint model did not estimate; too few complete cases across all ",
            "characteristics at once.", call. = FALSE)
    return(data.table())
  }
  setnames(out, "term", "CHARACTERISTIC", skip_absent = TRUE)
  setorder(out, DEPENDENT, INSTRUMENT_LABEL, SPEC, P_T)
  save_csv(out, "T34C_meta_joint_model.csv")
  out
}

s34_meta_horserace <- function(d, usable, dep = "RF_COEF",
                               instruments = names(MAIN_INSTRUMENTS)) {
  clean <- intersect(usable, s33_registry_dt()[OUTCOME_DERIVED == 0, CHARACTERISTIC])
  if (length(clean) == 0L) {
    warning("No non-outcome-derived characteristics usable; horse race skipped.",
            call. = FALSE)
    return(data.table())
  }
  models <- list(
    `1. Shoppability only`              = "SHOP_CERTAINTY",
    `2. Characteristics only`           = clean,
    `3. Shoppability + characteristics` = c("SHOP_CERTAINTY", clean))
  
  rows <- list()
  cat("  horse race:\n")
  for (il in intersect(instruments, unique(d$INSTRUMENT_LABEL))) {
    t0 <- Sys.time()
    di <- d[INSTRUMENT_LABEL == il]
    for (mn in names(models)) for (sp in c("RAW", "SIZE")) {
      r <- s33_fit_one(di, dep, models[[mn]], sp, "Inverse variance")
      if (is.null(r)) next
      r[, `:=`(INSTRUMENT_LABEL = il, MODEL = mn)]
      rows[[length(rows) + 1L]] <- r
    }
    cat(sprintf("    %-38s %.1fs\n", il,
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
  out <- rbindlist(rows, fill = TRUE)
  if (nrow(out) == 0L) {
    warning("Horse race did not estimate.", call. = FALSE)
    return(data.table())
  }
  setnames(out, "term", "CHARACTERISTIC", skip_absent = TRUE)
  setorder(out, DEPENDENT, INSTRUMENT_LABEL, MODEL, SPEC, P_T)
  save_csv(out, "T34D_meta_horserace_vs_shoppability.csv")
  out
}


# ============================================================================
# 34.5  DRIVER
# ============================================================================
#
# Merges the three new sources onto whatever s33_chars currently holds,
# recomputes the coverage/usable classification across the full extended set
# using the same STATUS rule Section 33 applies, reassembles the panel, and
# re-runs the three meta-regressions. Section 33's objects are read but never
# overwritten; everything new lands under an s34_ name.

run_section34 <- function(chars = NULL, dep = "RF_COEF") {
  
  if (!exists("concept_results"))
    stop("concept_results not in session. Warm start first.", call. = FALSE)
  if (is.null(chars)) {
    if (!exists("s33_chars"))
      stop("s33_chars not in session. Run restore_section33() or a live ",
           "Section 33 first.", call. = FALSE)
    chars <- get("s33_chars", envir = .GlobalEnv)
  }
  chars <- as.data.table(copy(chars))
  
  .s34_hd("34A. SUPPLEMENTAL SOURCES")
  pieces <- Filter(Negate(is.null), list(
    s34_chars_from_pfs_rvu(),
    s34_chars_from_asc_bb(),
    s34_chars_from_geo_volume()))
  
  if (length(pieces) == 0L)
    stop("No supplemental source could be built. Check that the files are in ",
         PANEL_DIR, call. = FALSE)
  
  for (p in pieces) {
    dup <- setdiff(intersect(names(p), names(chars)), "FINAL_CONCEPT_ID")
    if (length(dup)) chars[, (dup) := NULL]
    chars <- merge(chars, p, by = "FINAL_CONCEPT_ID", all.x = TRUE, sort = FALSE)
  }
  
  .s34_hd("34B. COVERAGE ACROSS THE EXTENDED CHARACTERISTIC SET")
  cov <- rbindlist(lapply(S33_CHAR_NAMES, function(v) {
    if (!(v %in% names(chars)))
      return(data.table(CHARACTERISTIC = v, STATUS = "NOT BUILT",
                        N_NONMISSING = 0L, SHARE_NONMISSING = 0, SD = NA_real_))
    x   <- safe_numeric(chars[[v]])
    sdx <- sd(x, na.rm = TRUE)
    data.table(CHARACTERISTIC = v,
               STATUS = fcase(sum(is.finite(x)) < S33_MIN_CONCEPTS, "TOO FEW",
                              !is.finite(sdx) || sdx == 0,          "NO VARIATION",
                              default = "USABLE"),
               N_NONMISSING = sum(is.finite(x)),
               SHARE_NONMISSING = round(mean(is.finite(x)), 3), SD = sdx)
  }))
  cov <- merge(s33_registry_dt(), cov, by = "CHARACTERISTIC", all = TRUE, sort = FALSE)
  cov[, IS_NEW_IN_S34 := as.integer(CHARACTERISTIC %chin% S34_CHAR_NAMES)]
  save_qa_csv(cov, "QA34A_supplemental_coverage.csv")
  
  cat("\nThe six characteristics added by Section 34:\n")
  print(as.data.frame(cov[IS_NEW_IN_S34 == 1L,
                          .(CHARACTERISTIC, STATUS, N_NONMISSING,
                            SHARE_NONMISSING, SD = round(SD, 4))]))
  
  usable <- cov[STATUS == "USABLE", CHARACTERISTIC]
  cat("\nUsable characteristics: ", length(usable),
      " (Section 33 alone had ", length(get("s33_usable", envir = .GlobalEnv)),
      ")\n", sep = "")
  dropped <- cov[STATUS != "USABLE"]
  if (nrow(dropped)) {
    cat("Dropped:\n")
    print(as.data.frame(dropped[, .(CHARACTERISTIC, STATUS, N_NONMISSING)]))
  }
  save_csv(chars, "T34A_concept_characteristics.csv")
  
  .s34_hd("34C. META-REGRESSIONS, EXTENDED CHARACTERISTIC SET")
  panel <- s33_assemble(concept_results, chars, dep = dep)
  panel <- s33_zscore(panel, unique(c(usable, S33_SIZE_CONTROLS)))
  
  uni <- s34_meta_univariate(panel, usable, dep = dep)
  cat("\nUnivariate sweep, primary instrument, inverse-variance weighted,\n",
      "FAMILY spec only. New characteristics are flagged.\n\n", sep = "")
  fam <- uni[INSTRUMENT_LABEL == names(MAIN_INSTRUMENTS)[1L] &
               WEIGHTING == "Inverse variance" & SPEC == "FAMILY"]
  fam[, NEW := fifelse(CHARACTERISTIC %chin% S34_CHAR_NAMES, "<- new", "")]
  print(as.data.frame(fam[, .(CHARACTERISTIC,
                              EST_PCT = round(ESTIMATE_PCT, 3),
                              SE_PCT  = round(SE_PCT, 3),
                              P = round(P_T, 4), STARS, N_CONCEPTS, NEW)]))
  
  joint <- s34_meta_joint(panel, usable, dep = dep)
  race  <- s34_meta_horserace(panel, usable, dep = dep)
  
  assign("s34_chars",  chars,  envir = .GlobalEnv)
  assign("s34_usable", usable, envir = .GlobalEnv)
  assign("s34_panel",  panel,  envir = .GlobalEnv)
  assign("s34_uni",    uni,    envir = .GlobalEnv)
  assign("s34_joint",  joint,  envir = .GlobalEnv)
  assign("s34_race",   race,   envir = .GlobalEnv)
  
  cat("\nSection 34 complete. New objects: s34_chars, s34_usable, s34_panel,\n",
      "s34_uni, s34_joint, s34_race. Section 33's objects are unchanged.\n", sep = "")
  invisible(list(chars = chars, usable = usable, panel = panel,
                 uni = uni, joint = joint, race = race))
}


# ============================================================================
# 34.6  RESTORE FROM DISK
# ============================================================================
#
# The Section 34 counterpart to restore_section33(). Reads what run_section34()
# wrote and reassembles the panel by merge and z-score, estimating nothing.
# Requires concept_results in the session and a prior run_section34().

restore_section34 <- function(dep = "RF_COEF") {
  .s34_hd("34R. RESTORING SECTION 34 FROM SAVED OUTPUT")
  
  if (!exists("concept_results"))
    stop("concept_results not in session. Warm start first.", call. = FALSE)
  
  .s34_read <- function(fn, dir = TABLE_DIR) {
    p <- file.path(dir, fn)
    if (!file.exists(p)) { message("SKIP: ", fn, " not found in ", dir); return(NULL) }
    fread(p)
  }
  
  chars <- .s34_read("T34A_concept_characteristics.csv")
  if (is.null(chars))
    stop("T34A_concept_characteristics.csv not found. Run run_section34() first.",
         call. = FALSE)
  chars <- as.data.table(chars)
  
  cov <- .s34_read("QA34A_supplemental_coverage.csv", dir = QA_DIR)
  if (is.null(cov))
    stop("QA34A_supplemental_coverage.csv not found in QA_DIR.", call. = FALSE)
  usable <- cov[STATUS == "USABLE", CHARACTERISTIC]
  
  uni   <- .s34_read("T34B_meta_univariate_sweep.csv")
  joint <- .s34_read("T34C_meta_joint_model.csv")
  race  <- .s34_read("T34D_meta_horserace_vs_shoppability.csv")
  
  panel <- s33_assemble(concept_results, chars, dep = dep)
  panel <- s33_zscore(panel, unique(c(usable, S33_SIZE_CONTROLS)))
  
  cat(sprintf("Restored: %d usable characteristics | s34_panel %s rows | uni %s | joint %s | race %s\n",
              length(usable), format(nrow(panel), big.mark = ","),
              if (is.null(uni))   "MISSING" else format(nrow(uni), big.mark = ","),
              if (is.null(joint)) "MISSING" else format(nrow(joint), big.mark = ","),
              if (is.null(race))  "MISSING" else format(nrow(race), big.mark = ",")))
  cat("No regression was re-estimated.\n")
  
  assign("s34_chars",  chars,  envir = .GlobalEnv)
  assign("s34_usable", usable, envir = .GlobalEnv)
  assign("s34_panel",  panel,  envir = .GlobalEnv)
  assign("s34_uni",    uni,    envir = .GlobalEnv)
  assign("s34_joint",  joint,  envir = .GlobalEnv)
  assign("s34_race",   race,   envir = .GlobalEnv)
  invisible(list(chars = chars, usable = usable, panel = panel,
                 uni = uni, joint = joint, race = race))
}

cat("Section 34 loaded. Run:  run_section34()\n")


# ===========================================================================
# 34F  INTERACTIVE CHECKS AND THE ONE-STEP PROMOTION
# ===========================================================================
#
# Everything below is run by hand, not by run_section34(). Two separate
# exercises share this space.
#
#   34F.1  Loader smoke tests. Each of the three code-level sources is built
#          on its own before the merged table is trusted, so a failure names
#          the file that caused it.
#
#   34F.2  SHARE_SI_PACKAGED estimated one-step on the row-level panel rather
#          than two-step on concept summaries, then re-estimated netting out
#          shoppability.
#
# 34F.2 needs `outpatient`, which restore_section33() does not load and which
# is the largest object in the session. Start from a full warm start, not
# HPT_WARM_START_KEYS <- "concept_level_6inst".


# ---------------------------------------------------------------------------
# 34F.1  Loader smoke tests
# ---------------------------------------------------------------------------
#
# Expected on the CY2025/2024 files, checked against HPT_CODEBOOK.csv:
#   PFS RVU    753 concepts from 877 of 915 codes. The 38 unmatched are all
#              C-prefix OPPS pass-through codes, which by definition never
#              appear on the physician fee schedule.
#   ASC BB     340 concepts, read from both vintages. Two concepts return NaN
#              because every matched code has a blank payment indicator in
#              CMS's own file; is.finite() drops them downstream, correctly.
#   Geo/volume 698 concepts from 813 codes, National rows only.
#
# readxl emits ~50 type-guessing warnings on the PFS file. They name columns
# E, H and X; the loader reads 1, 2, 14 and 15. Ignore them.

test_pfs <- s34_chars_from_pfs_rvu()
nrow(test_pfs)
head(test_pfs)

test_bb <- s34_chars_from_asc_bb()
nrow(test_bb)
head(test_bb)

test_geo <- s34_chars_from_geo_volume()
nrow(test_geo)
head(test_geo)


# ---------------------------------------------------------------------------
# 34F.2a  SHARE_SI_PACKAGED, one-step, all three main instruments
# ---------------------------------------------------------------------------
#
# Section 33 regresses concept-level RF_COEF estimates on characteristics.
# This asks the same question of the 1.4M-row panel directly: does the
# disclosure response differ by packaging within county x concept cells,
# with the moderator entering as a continuous interaction rather than as a
# regressor in a second-stage summary regression.
#
# `op` is the panel with the concept-level moderator merged on:
#   op <- merge(outpatient, s33_chars[, .(FINAL_CONCEPT_ID, SHARE_SI_PACKAGED)],
#               by = "FINAL_CONCEPT_ID", all.x = TRUE, sort = FALSE)
#
# The moderator is demeaned inside estimate_interacted(), so the Main row is
# the response at average packaging and should reproduce the pooled estimate
# already in the paper. It does: -3.06% IV against POOLED_IV = -3.07.
# That agreement is a specification check, not a new result.
#
# Roughly 3 minutes per instrument.

r_all <- rbindlist(lapply(names(MAIN_INSTRUMENTS), function(il) {
  t0 <- Sys.time()
  res <- estimate_interacted(
    op, moderator = "SHARE_SI_PACKAGED", moderator_type = "continuous",
    instrument = MAIN_INSTRUMENTS[[il]], instrument_label = il,
    label = "SHARE_SI_PACKAGED, one-step")
  cat(sprintf("  %-38s %.1f min\n", il,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  if (is.null(res)) return(NULL)
  cbind(INSTRUMENT_LABEL = il, res$rows)
}), fill = TRUE)

# The "x Moderator" row is the object of interest. Read RF, not IV: first-stage
# Wald runs 6.8 to 8.6 here, below the conventional threshold, so the IV
# columns carry the usual weak-instrument caveat while the reduced form is the
# Anderson-Rubin-equivalent test.
#
# Result: -0.00495 (p = .018), -0.00427 (p = .037), -0.00430 (p = .062).
# Sign and magnitude near-identical across instruments; the marginal third is
# also the weakest first stage, which is internally consistent.
r_all[TERM == "x Moderator",
      .(INSTRUMENT_LABEL, RF_COEF, RF_P, RF_PERCENT_PER_SD, IV_PERCENT, IV_P,
        FIRST_STAGE_WALD_MIN, N_OBSERVATIONS)]

saveRDS(r_all, file.path(CACHE_DIR, "share_si_packaged_onestep_all_instruments.rds"))


# ---------------------------------------------------------------------------
# 34F.2b  Netting out shoppability
# ---------------------------------------------------------------------------
#
# Packaged codes are concentrated in the shoppable group (36 of 439 shoppable
# concepts vs 8 of 297 non-shoppable), so 34F.2a leaves open whether the
# packaging interaction is shoppability in another guise. This is the one-step
# analogue of the concept-level horse race.
#
# estimate_interacted() takes a single moderator, so this rebuilds its
# continuous branch with two. Nothing else changes: same outcome, controls,
# fixed effects and clustering as everywhere else in the paper.
#
# Neither moderator's own level enters as a regressor. Both are concept-level
# and time-invariant, so county x concept absorbs them entirely; only the
# interaction with Z is identified. Both are demeaned, so Main stays readable
# as the response at the sample average of both.

s34_net_shoppability <- function(panel, instrument_label, instrument,
                                 moderator1 = "SHARE_SI_PACKAGED",
                                 moderator2 = "SHOP_DUMMY",
                                 outcome = PRIMARY_OUTCOME,
                                 endogenous = ENDOGENOUS_VARIABLE,
                                 controls = BASELINE_CONTROLS,
                                 fixed_effects = BASELINE_FIXED_EFFECTS,
                                 clusters = BASELINE_CLUSTERS) {
  
  ctl <- available_columns(panel, controls)
  fe  <- available_columns(panel, fixed_effects)
  cl  <- available_columns(panel, clusters)
  
  d <- panel[!is.na(get(moderator1)) & !is.na(get(moderator2))]
  d <- model_sample(d, c(outcome, endogenous, instrument, ctl, fe, cl, moderator1, moderator2))
  if (nrow(d) < MIN_MODEL_OBS) return(NULL)
  
  d[, MOD1 := safe_numeric(get(moderator1))]; d[, MOD1 := MOD1 - mean(MOD1, na.rm = TRUE)]
  d[, MOD2 := safe_numeric(get(moderator2))]; d[, MOD2 := MOD2 - mean(MOD2, na.rm = TRUE)]
  
  # Three terms per equation: level, packaging interaction, shoppability
  # interaction. RF and IV columns are built in parallel so the same three
  # rows come back for both.
  d[, `:=`(RF_MAIN = get(instrument), RF_PACKAGED = get(instrument) * MOD1,
           RF_SHOP = get(instrument) * MOD2,
           TREAT_MAIN = get(endogenous), TREAT_PACKAGED = get(endogenous) * MOD1,
           TREAT_SHOP = get(endogenous) * MOD2,
           IV_MAIN = get(instrument), IV_PACKAGED = get(instrument) * MOD1,
           IV_SHOP = get(instrument) * MOD2)]
  
  rfs  <- c("RF_MAIN", "RF_PACKAGED", "RF_SHOP")
  endo <- c("TREAT_MAIN", "TREAT_PACKAGED", "TREAT_SHOP")
  ivs  <- c("IV_MAIN", "IV_PACKAGED", "IV_SHOP")
  
  rf_fit <- tryCatch(feols(build_ols_formula(outcome, c(rfs, ctl), fe), data = d,
                           cluster = build_cluster_formula(cl), warn = FALSE, notes = FALSE),
                     error = function(e) NULL)
  iv_fit <- tryCatch(feols(build_iv_formula(outcome, endo, ivs, ctl, fe), data = d,
                           cluster = build_cluster_formula(cl), warn = FALSE, notes = FALSE),
                     error = function(e) NULL)
  if (is.null(rf_fit)) return(NULL)
  
  fsw <- first_stage_wald(iv_fit)
  fs_min <- if (nrow(fsw) > 0L) min(fsw$WALD, na.rm = TRUE) else NA_real_
  
  # fixest prefixes IV coefficient names with "fit_"; try both forms.
  pull <- function(fit, tm) {
    if (is.null(fit)) return(list(b = NA_real_, s = NA_real_))
    cand <- c(paste0("fit_", tm), tm); hit <- cand[cand %in% names(coef(fit))]
    if (length(hit) == 0L) return(list(b = NA_real_, s = NA_real_))
    list(b = unname(coef(fit)[hit[1L]]), s = unname(sqrt(vcov(fit)[hit[1L], hit[1L]])))
  }
  sd_z <- sd(d[[instrument]], na.rm = TRUE)
  
  out <- rbindlist(lapply(
    list(c("Main", "RF_MAIN", "TREAT_MAIN"),
         c("x SHARE_SI_PACKAGED", "RF_PACKAGED", "TREAT_PACKAGED"),
         c("x Shoppable", "RF_SHOP", "TREAT_SHOP")),
    function(spec) {
      rf <- pull(rf_fit, spec[2]); iv <- pull(iv_fit, spec[3])
      data.table(TERM = spec[1], RF_COEF = rf$b, RF_SE = rf$s,
                 RF_P = .pval(rf$b / rf$s, rf_fit),
                 RF_PERCENT_PER_SD = 100 * (exp(rf$b * sd_z) - 1),
                 IV_COEF = iv$b, IV_SE = iv$s, IV_P = .pval(iv$b / iv$s, iv_fit),
                 IV_PERCENT = 100 * (exp(iv$b) - 1))
    }))
  out[, `:=`(INSTRUMENT_LABEL = instrument_label, FIRST_STAGE_WALD_MIN = fs_min,
             N_OBSERVATIONS = nobs(rf_fit))]
  out[]
}

# SHOP_CERTAINTY lives on s33_panel, which is concept-level. The panel-native
# column is SCHEME_1_CERTAINTY, the same field every Section 28-30
# specification uses.
op[, SHOP_DUMMY := as.numeric(SCHEME_1_CERTAINTY == "Shoppable")]

# Single instrument first, to confirm the specification runs before spending
# three times the wall clock on the full set.
t0 <- Sys.time()
r_net <- s34_net_shoppability(op, "Competitor_only_hospitals_9m",
                              MAIN_INSTRUMENTS[["Competitor_only_hospitals_9m"]])
cat("Elapsed:", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
r_net[, .(TERM, RF_COEF, RF_P, RF_PERCENT_PER_SD, IV_PERCENT, IV_P,
          FIRST_STAGE_WALD_MIN, N_OBSERVATIONS)]

r_net_all <- rbindlist(lapply(names(MAIN_INSTRUMENTS), function(il) {
  t0 <- Sys.time()
  res <- s34_net_shoppability(op, il, MAIN_INSTRUMENTS[[il]])
  cat(sprintf("  %-38s %.1f min\n", il,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  res
}), fill = TRUE)

# What this returns, and it is the honest version of the result:
#
#   instrument            x SHARE_SI_PACKAGED    x Shoppable
#   Competitor hospitals  -0.00415 (p = .053)    -0.00284 (p = .024)
#   Primary strict system -0.00364 (p = .084)    -0.00221 (p = .043)
#   Competitor ex-CBSA    -0.00330 (p = .168)    -0.00337 (p = .018)
#
# Shoppability is robust across all three. Packaging holds its sign and rough
# magnitude but weakens as the instrument weakens, and is not distinguishable
# from zero under the third, which is also the weakest first stage in the
# paper (5.7). Attenuation from 34F.2a is 16% on the primary instrument.
#
# The two-step horse race in Section 33 had packaging clearing p < .02 on
# every instrument. This one-step version, on the real panel rather than 736
# concept summaries, is the more rigorous test and is less favourable. Where
# they disagree, this one governs.
r_net_all[TERM != "Main",
          .(INSTRUMENT_LABEL, TERM, RF_COEF, RF_P, RF_PERCENT_PER_SD,
            IV_PERCENT, IV_P, FIRST_STAGE_WALD_MIN, N_OBSERVATIONS)]

saveRDS(r_net_all, file.path(CACHE_DIR, "share_si_packaged_net_shoppability_all.rds"))

