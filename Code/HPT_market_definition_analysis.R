###############################################################################
#
#   WHEN TRANSPARENCY WORKS: SERVICE SHOPPABILITY, CONTRACTING DEPTH, AND
#   THE PRICE EFFECTS OF HOSPITAL DISCLOSURE
#
#   Replication code: market definition and the exclusion restriction
#
#   Danny Sierra
#   Department of Economics, Florida State University
#   Ds22c@fsu.edu
#
###############################################################################
#
# HOW TO RUN
# -----------------------------------------------------------------------------
# This file is a companion to HPT_Analysis_Pipeline.R and uses the pipeline's
# functions, constants, paths, and cached results. Source the pipeline first,
# in the same R session. A warm start loads its functions and cached results
# without estimating anything:
#   Sys.setenv(HPT_ROOT = "/path/to/Hospital Price Transparency Paper")
#   HPT_WARM_START <- TRUE
#   source(file.path(Sys.getenv("HPT_ROOT"), "Code", "HPT_Analysis_Pipeline.R"))
# Then source this file. A session in which the pipeline has completed a full
# run also works. Either way, source() of the pipeline stops in its Section
# 34F with "object 'op' not found" (see the pipeline's KNOWN ISSUES); every
# object this file uses is defined before that point. This file has no run
# switches: source() runs every part in order, and cache_or_run() follows the
# pipeline's USE_CACHE.
#
# Part 0 checks the pipeline objects listed in MD_REQUIRED_OBJECTS by name and
# stops if any is missing. Sourcing the pipeline only through PART 2 does not
# pass this check: CONCEPT_INSTRUMENTS is defined in the pipeline's PART 3
# (BUILD block) or by restore_session(), which a warm start calls. Two cached
# results of the main county run must also be in memory: concept_results
# (cache key concept_level_6inst) from Part 9 (9C) onward, and schemes_long in
# Part 15E. A warm start loads both; if HPT_WARM_START_KEYS restricts it, the
# keys must include concept_level_6inst and schemes_long.
#
# Inputs. Part 0 stops if any of these is missing from PANEL_DIR
# (Data/data_final/01_R_Analysis_Panels):
#   HPT_R_MAIN_PRIMARY_COUNTY_OUTPATIENT_CONCEPT.parquet
#   HPT_R_MAIN_ROBUSTNESS_<LEVEL>_OUTPATIENT_CONCEPT.parquet, for LEVEL = CITY,
#     CBSA, HSA, and HRR
#   HPT_RING_EXPOSURE.csv, HPT_ADJACENCY_DECOMP.csv, HPT_HOSPITAL_GEO.csv
# Read later in the file:
#   HPT_HOSPITAL_DISCLOSURE_EVENTS.parquet, in
#     Data/data_final/02_Treatments_and_Instruments (10.1)
#   the Dartmouth HSA boundary shapefile, in PANEL_DIR (11.1)
#   HPT_COUNTY_ADJACENCY_PAIRS.csv, in PANEL_DIR (11.2)
# Every call to md_load_panel() also reads the pipeline's county exact-code
# panel (FILES$outpatient_exact) through apply_concept_merges().
#
# WHAT THE FILE ESTIMATES
# -----------------------------------------------------------------------------
# Two questions about the headline shoppability result.
#
# 1. Market definition. The headline gradient, estimated on county markets,
#    is re-estimated with markets defined by city, CBSA, HSA, and HRR, with
#    the number of prior posters within 30 miles as the treatment, in
#    rural/urban and multi-hospital subsamples, and at the concept level at
#    HSA and HRR.
#
# 2. Exclusion restriction. Instruments IV01-IV14 apply a "different county"
#    filter and, in the CBSA variants, an "outside the CBSA" filter. Neither
#    excludes adjacent counties, and metropolitan areas span adjacent
#    counties. IV15 (Z_SYS_COMPETITOR_EXCL_ADJACENT_9M_EXCL_CURRENT) is IV06
#    (PRIMARY_INSTRUMENT, Z_SYS_COMPETITOR_ONLY_9M_EXCL_CURRENT) with peers in
#    adjacent counties excluded. The file compares the two and estimates the
#    own-county and adjacent-county channels separately. Part 17 adds
#    out-of-state and ownership-dated versions of IV06.
#
# Unless a part says otherwise, models use the pipeline's baseline
# specification: LOG_TOTAL_BEDS as control, MARKET_ID and POST_MONTH fixed
# effects, and two-way clustering by ANALYSIS_MARKET and POST_MONTH, where
# ANALYSIS_MARKET is the market of the panel in use (0.2).
#
# Its tables are CSV files whose names begin with MD, written to TABLE_DIR
# (MD06 to QA_DIR). One block also writes under the pipeline's file names:
# 15E.7 calls the pipeline's run_comparability_within_family() and
# run_comparability_meta(), which save T09D_comparability_within_family.csv,
# T10_comparability_meta_regression.csv, and T10B_horse_race_summary.csv to
# TABLE_DIR and QA09_moderator_screen.csv to QA_DIR (see KNOWN ISSUES).
# Parts 15E and 16 also write PDF figures to FIGURE_DIR and
# tab_schemes_geography_body.tex to TABLE_DIR. Its own caches use keys that
# begin with md_, p12_, or p15_. It also reads three of the pipeline's caches:
# schemes_long (Part 2), meta_regressions_iv (15E.3), and
# comparability_measures (15E.7). When a cache file is missing or USE_CACHE
# is FALSE, schemes_long is rebuilt and saved with the pipeline's own call,
# and the other two fail (see KNOWN ISSUES). The file also changes the
# session: 0.1 extends ANALYSIS_COLUMNS, 11.1 sets sf_use_s2(FALSE), and 15.0
# removes `outpatient` and other large pipeline objects from memory.
#
# STRUCTURE
# -----------------------------------------------------------------------------
# The parts in the order they run:
#   Part 0    Preflight: required objects, column list, input files
#   Part 1    Hospital-level support data (rural/urban, rings, adjacency)
#   Part 2    Panel loader, md_load_panel()
#   Part 3    Market definition ladder
#   Part 4    Ring exposure and the spatial decay diagnostic
#   Part 5    Own-county against adjacent-county decomposition
#   Part 6    IV15 against IV06
#   Part 7    Rural/urban and multi-hospital subsamples
#   Part 9    Concept-level analysis at HRR, one instrument
#   Part 8    Assembled summary of Parts 3-7
#   Part 10   Competitor instruments rebuilt in R and validated; six at HRR
#   Part 11   HSA adjacency from boundary polygons; nine HSA instruments
#   Part 12   Concept-level estimation and meta-regressions at HSA and HRR
#   Part 13   Between-concept heterogeneity across county, HSA, and HRR
#   Part 14   Fixed-effect feasibility
#   Part 15   Six-instrument concept runs and the restricted county run
#   Part 15E  Cross-geography comparison tables and figures
#   Part 16   Distance-banded instruments and the ring decay IV
#   Part 17   Out-of-state and ownership-dated constructions
#
# Part 9 is placed before Part 8 in the file and runs first; Part 8 uses only
# the results of Parts 3-7. The part numbers are kept because the paper and
# its appendix refer to them.
#
# RUNTIME
# -----------------------------------------------------------------------------
# The concept-level loops in Parts 9, 12, and 15 take most of the time: about
# 11-12 hours per single-instrument geography and 8-12 hours per
# six-instrument geography. Each loop runs inside cache_or_run(), so with
# USE_CACHE = TRUE a completed run is loaded instead of repeated.
# estimate_concept_level() saves a _PARTIAL.csv every 50 concepts as a record
# of the concepts finished. It does not resume from that file, so an
# interrupted loop starts again from the first concept.
#
# Two costs remain when every result is cached. The Part 3 ladder is cached
# with overwrite = TRUE, so its five models are re-estimated on every run.
# md_load_panel() runs for every panel a block needs (five times in Part 3,
# once each in Parts 4, 9, 10.6, and 11.5, twice per run in Part 12, and again
# in Parts 15 and 17.4), also when the estimates come from the cache, and each
# call reads the county exact-code panel as well.
#
#
###############################################################################

library(arrow)
library(data.table)
library(fixest)
library(ggplot2)
library(lubridate)
library(sf)


# =============================================================================
# Part 0: Preflight
# =============================================================================
#
# Every name in MD_REQUIRED_OBJECTS is defined by HPT_Analysis_Pipeline.R, not
# by this file (see HOW TO RUN). The check stops with the list of missing
# names, so a missing prerequisite produces an immediate error instead of a
# NULL model hours into a concept loop. The list holds the main functions and
# constants; other pipeline objects used below, such as wald_equality(),
# SCHEME_COLUMNS, and FINAL_DATA_ROOT, are defined in PARTS 1 and 2 of the
# pipeline.

MD_REQUIRED_OBJECTS <- c(
  "estimate_interacted", "prepare_panel", "read_panel", "build_schemes",
  "extend_schemes_for_merged", "attach_scheme_columns", "apply_concept_merges",
  "choose_sample_flag", "cache_or_run", "save_csv", "save_qa_csv",
  "available_columns", "model_sample", "build_iv_formula",
  "build_cluster_formula", "first_stage_wald", ".pval",
  "estimate_concept_level", "build_concept_panel", "prepare_meta_input",
  "run_meta_regressions",
  "ENDOGENOUS_VARIABLE", "BASELINE_CONTROLS", "BASELINE_FIXED_EFFECTS",
  "BASELINE_CLUSTERS", "PRIMARY_OUTCOME", "PRIMARY_INSTRUMENT",
  "MAIN_INSTRUMENTS", "SUPPORTING_INSTRUMENTS", "CONCEPT_INSTRUMENTS",
  "MIN_MODEL_OBS", "PANEL_DIR", "ANALYSIS_COLUMNS"
)

md_missing <- MD_REQUIRED_OBJECTS[!vapply(MD_REQUIRED_OBJECTS, exists, logical(1))]
if (length(md_missing) > 0L) {
  stop("Definitions from HPT_Analysis_Pipeline.R are not loaded.\n",
       "Missing: ", paste(md_missing, collapse = ", "), "\n\n",
       "Run Parts 1-2 of the main pipeline first, or restore_session().",
       call. = FALSE)
}

cat("\n=== MARKET DEFINITION ANALYSIS ===\n")
cat("Endogenous:", ENDOGENOUS_VARIABLE, "| Primary IV:", PRIMARY_INSTRUMENT, "\n")


# -----------------------------------------------------------------------------
# 0.1  Extend ANALYSIS_COLUMNS so the IV15 family survives read_panel()
# -----------------------------------------------------------------------------
# read_panel() keeps only intersect(ANALYSIS_COLUMNS, names(ds)). The three
# IV15-family columns in IV_EXCLUDE_ADJACENT are not in
# ALL_CANDIDATE_INSTRUMENTS, so without this extension they are dropped at
# load with no warning and any model that uses them returns NULL. HSA_NUM,
# HRR_NUM, and RUCC_2023 are added so that they are kept where a panel has
# them. The extended ANALYSIS_COLUMNS stays in effect for the rest of the
# session.

IV_EXCLUDE_ADJACENT <- c(
  Competitor_excl_adjacent_hospitals_9m = "Z_SYS_COMPETITOR_EXCL_ADJACENT_9M_EXCL_CURRENT",
  Competitor_excl_adjacent_systems_9m   = "Z_SYS_COMPETITOR_EXCL_ADJACENT_SYSTEMS_9M_EXCL_CURRENT",
  Competitor_excl_adjacent_counties_9m  = "Z_SYS_COMPETITOR_EXCL_ADJACENT_COUNTIES_9M_EXCL_CURRENT"
)

ANALYSIS_COLUMNS <- unique(c(ANALYSIS_COLUMNS,
                             unname(IV_EXCLUDE_ADJACENT),
                             "HSA_NUM", "HRR_NUM", "RUCC_2023"))

cat("ANALYSIS_COLUMNS extended with", length(IV_EXCLUDE_ADJACENT),
    "adjacency instruments.\n")


# -----------------------------------------------------------------------------
# 0.2  Market definitions and their panel files
# -----------------------------------------------------------------------------
# The headline specification uses the concept panel
# (FILES$outpatient_concept), not the exact-code panel, so every rung of the
# ladder reads the concept panel for its geography. md_panel_path() returns
# the file in PANEL_DIR: HPT_R_MAIN_PRIMARY_COUNTY_OUTPATIENT_CONCEPT.parquet
# for COUNTY and HPT_R_MAIN_ROBUSTNESS_<LEVEL>_OUTPATIENT_CONCEPT.parquet for
# the other levels.
#
# Fixed effects and clustering follow the panel. Each exported panel renames
# its geography column to ANALYSIS_MARKET and builds MARKET_ID as
# ANALYSIS_MARKET::FINAL_CONCEPT_ID. BASELINE_FIXED_EFFECTS (MARKET_ID,
# POST_MONTH) and BASELINE_CLUSTERS (ANALYSIS_MARKET, POST_MONTH) use these
# generic names, so changing the panel changes the fixed effects and the
# cluster level together. Setting them by hand for one geography would put
# the fixed effects out of step with the treatment without any error.

MD_GEOGRAPHIES <- c("COUNTY", "CITY", "CBSA", "HSA", "HRR")

md_panel_path <- function(level) {
  role <- if (level == "COUNTY") "PRIMARY" else "ROBUSTNESS"
  file.path(PANEL_DIR,
            sprintf("HPT_R_MAIN_%s_%s_OUTPATIENT_CONCEPT.parquet", role, level))
}

# Hospital-level support files in PANEL_DIR: ring exposure, the own- and
# adjacent-county decomposition, and hospital geography. The two loops below
# stop if any support file or any of the five panels is missing.
MD_SUPPORT_FILES <- list(
  ring      = file.path(PANEL_DIR, "HPT_RING_EXPOSURE.csv"),
  adjacency = file.path(PANEL_DIR, "HPT_ADJACENCY_DECOMP.csv"),
  geo       = file.path(PANEL_DIR, "HPT_HOSPITAL_GEO.csv")
)

for (nm in names(MD_SUPPORT_FILES)) {
  if (!file.exists(MD_SUPPORT_FILES[[nm]]))
    stop("Missing support file (", nm, "): ", MD_SUPPORT_FILES[[nm]], call. = FALSE)
}
for (g in MD_GEOGRAPHIES) {
  if (!file.exists(md_panel_path(g)))
    stop("Missing panel for ", g, ": ", md_panel_path(g), call. = FALSE)
}

MD_SCHEME <- "SCHEME_1_CERTAINTY"   # headline classification


# =============================================================================
# Part 1: Hospital-level support data
# =============================================================================
#
# Builds md_hospital_attrs, one row per hospital, with RUCC_2023 and
# METRO_STATUS (RUCC 1-3 metro, 4-9 nonmetro), prior-poster counts in
# distance rings around the hospital (0-15, 15-30, 30-60, and 0-30 miles)
# and RING_TOTAL_0_30, and the own-county and adjacent-county prior-poster
# counts. md_load_panel() (Part 2) merges it onto every panel by HOSPITAL_ID.
#
# unique() keeps one row per HOSPITAL_ID in each support file. Four hospitals
# (La Salle IL, La Porte IN, and two Valdez-Cordova AK records) appear twice
# in HPT_HOSPITAL_GEO because the join in Section 3 of the Snowflake build
# matches them through both HPT_REF_GEO_OVERRIDES and the name crosswalk. The
# duplicate rows are byte-identical, so keeping the first loses nothing.

md_hosp_geo <- unique(
  fread(MD_SUPPORT_FILES$geo,
        colClasses = list(character = c("HOSPITAL_ID", "COUNTY_FIPS",
                                        "HSA_NUM", "HRR_NUM"))),
  by = "HOSPITAL_ID")
md_hosp_geo[, HOSPITAL_ID := as.character(HOSPITAL_ID)]

md_ring <- unique(fread(MD_SUPPORT_FILES$ring), by = "HOSPITAL_ID")
md_ring[, HOSPITAL_ID := as.character(HOSPITAL_ID)]

md_adj <- unique(fread(MD_SUPPORT_FILES$adjacency), by = "HOSPITAL_ID")
md_adj[, HOSPITAL_ID := as.character(HOSPITAL_ID)]

md_hospital_attrs <- Reduce(
  function(x, y) merge(x, y, by = "HOSPITAL_ID", all.x = TRUE),
  list(
    md_hosp_geo[, .(HOSPITAL_ID, RUCC_2023, COUNTY_FIPS_GEO = COUNTY_FIPS)],
    md_ring[, .(HOSPITAL_ID,
                RING_PRIOR_0_15  = N_PRIOR_POSTERS_0_15MI,
                RING_PRIOR_15_30 = N_PRIOR_POSTERS_15_30MI,
                RING_PRIOR_30_60 = N_PRIOR_POSTERS_30_60MI,
                RING_PRIOR_0_30  = N_PRIOR_POSTERS_0_30MI,
                RING_TOTAL_0_30  = N_TOTAL_0_30MI)],
    md_adj[, .(HOSPITAL_ID,
               ADJ_PRIOR_OWN      = N_PRIOR_POSTERS_OWN_COUNTY,
               ADJ_PRIOR_ADJACENT = N_PRIOR_POSTERS_ADJACENT_COUNTY)]
  ))

md_hospital_attrs[, METRO_STATUS := fifelse(
  RUCC_2023 %in% 1:3, "Metro",
  fifelse(RUCC_2023 %in% 4:9, "Nonmetro", NA_character_))]

cat("Hospital attributes assembled:",
    format(nrow(md_hospital_attrs), big.mark = ","), "hospitals\n")


# =============================================================================
# Part 2: Panel loader
# =============================================================================
#
# md_load_panel(level) builds a geography's panel with the steps the
# pipeline's BUILD block uses for `outpatient`: read_panel() and
# prepare_panel() as in load_outpatient(), then apply_concept_merges() and
# attach_scheme_columns(), applied to the file from md_panel_path(level). It
# then merges the Part 1 hospital attributes (stops if the row count changes)
# and sets ANALYSIS_GEOGRAPHY_LABEL. Following the same sequence keeps every
# rung comparable with the headline. md_schemes_long is read from the
# pipeline's schemes_long cache; if that file is missing or USE_CACHE is
# FALSE, it is rebuilt and saved with the pipeline's own call.
#
# Market keys of the merged concepts. apply_concept_merges() builds the six
# merged concepts from the county exact-code panel at every geography, so the
# merged rows carry a county ANALYSIS_MARKET. At CITY, CBSA, HSA, and HRR
# this replaces the market key with a county string in roughly 12,000 to
# 14,000 rows per geography. The loader therefore saves ANALYSIS_MARKET by
# HOSPITAL_ID, POST_MONTH, and FINAL_CONCEPT_ID before the merge, writes it
# back afterwards, and rebuilds MARKET_ID; apply_concept_merges() is used as
# the pipeline defines it. A stop() checks that the write-back leaves the row
# count unchanged, and a warning() flags keys that still look like county
# keys.
#
# A merged row whose combination of HOSPITAL_ID, POST_MONTH, and
# FINAL_CONCEPT_ID did not occur before the merge gets ANALYSIS_MARKET = NA,
# at every geography including COUNTY. This covers every row of the two
# constructed canonical concepts (MRI_MRA_FMRI_BRAIN, MRI_MRA_ARHFCMRIGTBS)
# and the rows of the other four merged concepts in hospital-months that had
# no row for the canonical concept itself. ANALYSIS_MARKET is a clustering
# variable, so these rows drop out of every model in this file, while the
# pipeline's `outpatient` panel keeps them.

md_schemes_long <- cache_or_run("schemes_long",
                                extend_schemes_for_merged(build_schemes()))

md_load_panel <- function(level) {
  raw <- read_panel(md_panel_path(level), paste(level, "concept panel"))
  cat("  raw panel size:", format(object.size(raw), units = "GB"), "\n")
  
  p <- prepare_panel(raw, paste(level, "concept panel"), choose_sample_flag(raw))
  rm(raw); invisible(gc())
  
  market_key_backup <- p[, .(HOSPITAL_ID, POST_MONTH, FINAL_CONCEPT_ID,
                             ANALYSIS_MARKET_CORRECT = ANALYSIS_MARKET)]
  
  p <- apply_concept_merges(p)
  invisible(gc())
  
  n_before <- nrow(p)
  p <- merge(p, market_key_backup,
             by = c("HOSPITAL_ID", "POST_MONTH", "FINAL_CONCEPT_ID"),
             all.x = TRUE)
  if (nrow(p) != n_before)
    stop("Market key restore changed row count for ", level,
         " (", n_before, " -> ", nrow(p), "). apply_concept_merges() may ",
         "have altered FINAL_CONCEPT_ID for some rows.", call. = FALSE)
  
  n_restored <- sum(p$ANALYSIS_MARKET != p$ANALYSIS_MARKET_CORRECT, na.rm = TRUE)
  cat("  Corrected", format(n_restored, big.mark = ","),
      "rows with a leaked county-level market key.\n")
  
  p[, ANALYSIS_MARKET := ANALYSIS_MARKET_CORRECT]
  p[, ANALYSIS_MARKET_CORRECT := NULL]
  p[, MARKET_ID := paste(ANALYSIS_MARKET, FINAL_CONCEPT_ID, sep = "::")]
  
  # A restored key longer than 10 characters at CITY, CBSA, HSA, or HRR is
  # taken to be a county key that was not restored, and triggers a warning.
  # At COUNTY every key is a county key, so the check is skipped. Keys set to
  # NA above are ignored by max(na.rm = TRUE).
  max_len <- max(nchar(p$ANALYSIS_MARKET), na.rm = TRUE)
  if (level != "COUNTY" && max_len > 10)
    warning(level, ": ANALYSIS_MARKET max length is ", max_len,
            " after restore -- verify the fix held.", call. = FALSE)
  
  p <- attach_scheme_columns(p, md_schemes_long)
  cat("  after scheme attach:", format(object.size(p), units = "GB"), "\n")
  invisible(gc())
  
  n_before <- nrow(p)
  p <- merge(p, md_hospital_attrs, by = "HOSPITAL_ID", all.x = TRUE)
  if (nrow(p) != n_before)
    stop("Hospital attribute merge changed row count for ", level,
         " (", n_before, " -> ", nrow(p), ").", call. = FALSE)
  
  p[, ANALYSIS_GEOGRAPHY_LABEL := level]
  invisible(gc())
  p
}


# =============================================================================
# Part 3: Market definition ladder
# =============================================================================
#
# Re-estimates the headline specification (estimate_interacted() with the
# Scheme 1 categories in MD_SCHEME, PRIMARY_OUTCOME, and PRIMARY_INSTRUMENT)
# once for each geography in MD_GEOGRAPHIES. MD01 has one row per category
# and geography, with the number of markets, fixed-effect cells, hospitals,
# and panel rows behind it; MD02 has the equality tests.
#
# The instrument is county-based at every rung. Every Z_SYS_* instrument is
# built at county level in the Python pipeline, so when the treatment and the
# fixed effects move to HRR the instrument does not. This is a choice rather
# than a consequence of swapping the panel: the instrument's role is to
# predict the posting decision, not to define the market. Parts 10 and 11
# build instruments at HRR and HSA.
#
# cache_or_run() is called with overwrite = TRUE, so the ladder is
# re-estimated, and its five panels reloaded, every time this part runs.

md_run_ladder_row <- function(level) {
  cat("\n=== Ladder:", level, "===\n")
  p <- md_load_panel(level)
  
  res <- estimate_interacted(
    data             = p,
    moderator        = MD_SCHEME,
    moderator_type   = "categorical",
    outcome          = PRIMARY_OUTCOME,
    instrument       = PRIMARY_INSTRUMENT,
    label            = paste0("Market=", level),
    instrument_label = names(MAIN_INSTRUMENTS)[
      match(PRIMARY_INSTRUMENT, MAIN_INSTRUMENTS)],
    moderator_label  = MD_SCHEME
  )
  
  if (is.null(res)) {
    cat("  Model returned NULL (below MIN_MODEL_OBS or no variation).\n")
    rm(p); invisible(gc())
    return(NULL)
  }
  
  rows <- copy(res$rows)
  rows[, `:=`(MARKET_DEFINITION = level,
              N_MARKETS  = uniqueN(p$ANALYSIS_MARKET),
              N_FE_CELLS = uniqueN(p$MARKET_ID),
              N_HOSPITALS = uniqueN(p$HOSPITAL_ID),
              PANEL_ROWS  = nrow(p))]
  
  tests <- if (!is.null(res$tests) && nrow(res$tests) > 0L) {
    t <- copy(res$tests); t[, MARKET_DEFINITION := level]; t
  } else NULL
  
  # Progress report for this geography (rows, hospitals, markets, the IV
  # estimates, and gc() memory use), printed before the next one is loaded.
  cat(sprintf("  DONE %-6s | rows: %s | hospitals: %s | markets: %s\n",
              level, format(nrow(p), big.mark = ","),
              format(uniqueN(p$HOSPITAL_ID), big.mark = ","),
              format(uniqueN(p$ANALYSIS_MARKET), big.mark = ",")))
  print(rows[, .(TERM, IV_PERCENT = round(IV_PERCENT, 3),
                 IV_P = round(IV_P, 4), N_OBSERVATIONS)])
  print(gc())
  
  rm(p); invisible(gc())
  list(rows = rows, tests = tests)
}

md_ladder <- cache_or_run("md_market_definition_ladder", {
  out <- lapply(MD_GEOGRAPHIES, md_run_ladder_row)
  list(
    rows  = rbindlist(lapply(out, `[[`, "rows"),  fill = TRUE),
    tests = rbindlist(lapply(out, `[[`, "tests"), fill = TRUE)
  )
}, overwrite = TRUE)

save_csv(md_ladder$rows,  "MD01_market_definition_ladder_rows.csv")
save_csv(md_ladder$tests, "MD02_market_definition_ladder_tests.csv")

cat("\n=== LADDER SUMMARY (corrected market keys) ===\n")
print(md_ladder$rows[, .(MARKET_DEFINITION, TERM, N_MARKETS,
                         IV_PERCENT = round(IV_PERCENT, 3),
                         IV_P = round(IV_P, 4),
                         FIRST_STAGE_WALD_MIN = round(FIRST_STAGE_WALD_MIN, 1),
                         N_OBSERVATIONS)])


# =============================================================================
# Part 4: Ring exposure and the spatial decay diagnostic
# =============================================================================
#
# md_county_panel is loaded once here and stays in memory for the rest of the
# file. It is used in Parts 4-7, as the focal panel for the instruments built
# in Parts 10, 11, 16, and 17, and in Parts 14 and 15.

md_county_panel <- md_load_panel("COUNTY")


# -----------------------------------------------------------------------------
# 4.1  Ring exposure as a treatment
# -----------------------------------------------------------------------------
# The treatment is replaced by the number of prior posters within 30 miles
# (RING_PRIOR_0_30). The instrument stays PRIMARY_INSTRUMENT, and
# ANALYSIS_MARKET and MARKET_ID stay county-based because a ring has no
# discrete market key. This rung changes the treatment without changing the
# fixed effects, so it is reported as a distance-based treatment check
# (MD03_ring_treatment_rows.csv), not as a sixth market definition.

md_ring_panel <- copy(md_county_panel)
md_ring_panel[, N_PRIOR_POSTERS_COUNTY_ORIG := get(ENDOGENOUS_VARIABLE)]
md_ring_panel[, (ENDOGENOUS_VARIABLE) := safe_numeric(RING_PRIOR_0_30)]
md_ring_panel <- md_ring_panel[is.finite(get(ENDOGENOUS_VARIABLE))]

md_ring_result <- estimate_interacted(
  data             = md_ring_panel,
  moderator        = MD_SCHEME,
  moderator_type   = "categorical",
  outcome          = PRIMARY_OUTCOME,
  instrument       = PRIMARY_INSTRUMENT,
  label            = "Market=RING_0_30MI",
  instrument_label = "Competitor_only_hospitals_9m",
  moderator_label  = MD_SCHEME
)

if (!is.null(md_ring_result)) {
  ring_rows <- copy(md_ring_result$rows)
  ring_rows[, `:=`(MARKET_DEFINITION = "RING_0_30MI",
                   N_MARKETS = NA_integer_,
                   N_FE_CELLS = uniqueN(md_ring_panel$MARKET_ID),
                   N_HOSPITALS = uniqueN(md_ring_panel$HOSPITAL_ID),
                   PANEL_ROWS = nrow(md_ring_panel))]
  save_csv(ring_rows, "MD03_ring_treatment_rows.csv")
  print(ring_rows[, .(MARKET_DEFINITION, TERM,
                      IV_PERCENT = round(IV_PERCENT, 3),
                      IV_P = round(IV_P, 4), N_OBSERVATIONS)])
}


# -----------------------------------------------------------------------------
# 4.2  Spatial decay, OLS
# -----------------------------------------------------------------------------
# One OLS regression of the outcome on the prior-poster counts in all three
# bands (0-15, 15-30, and 30-60 miles), with the baseline controls, fixed
# effects, and clusters and no instrument; the printed heading calls it the
# reduced form. It does not define a market. It estimates the distance at
# which the peer-posting relationship fades. Writes
# MD04_ring_spatial_decay.csv.
#
# The share of ring peers that posted first is nearly flat across the three
# bands (0.398, 0.402, and 0.407 for 0-15, 15-30, and 30-60 miles), so flat
# coefficients describe the setting, in which peer posting is not spatially
# local, rather than a data problem.

md_decay_fit <- tryCatch(
  feols(
    build_ols_formula(
      PRIMARY_OUTCOME,
      c("RING_PRIOR_0_15", "RING_PRIOR_15_30", "RING_PRIOR_30_60",
        available_columns(md_county_panel, BASELINE_CONTROLS)),
      available_columns(md_county_panel, BASELINE_FIXED_EFFECTS)),
    data    = md_county_panel[is.finite(RING_PRIOR_0_15) &
                                is.finite(RING_PRIOR_15_30) &
                                is.finite(RING_PRIOR_30_60)],
    cluster = build_cluster_formula(
      available_columns(md_county_panel, BASELINE_CLUSTERS)),
    warn = FALSE, notes = FALSE),
  error = function(e) { cat("  Decay model failed:", conditionMessage(e), "\n"); NULL })

if (!is.null(md_decay_fit)) {
  md_decay <- data.table(
    BAND  = c("0-15mi", "15-30mi", "30-60mi"),
    TERM  = c("RING_PRIOR_0_15", "RING_PRIOR_15_30", "RING_PRIOR_30_60"))
  md_decay[, COEF := vapply(TERM, function(t) unname(coef(md_decay_fit)[t]), numeric(1))]
  md_decay[, SE   := vapply(TERM, function(t) unname(sqrt(vcov(md_decay_fit)[t, t])), numeric(1))]
  md_decay[, `:=`(T_STAT = COEF / SE,
                  P_VALUE = .pval(COEF / SE, md_decay_fit),
                  N_OBSERVATIONS = nobs(md_decay_fit))]
  save_csv(md_decay, "MD04_ring_spatial_decay.csv")
  cat("\n=== SPATIAL DECAY (reduced form, all bands together) ===\n")
  print(md_decay[, .(BAND, COEF = round(COEF, 6), SE = round(SE, 6),
                     P_VALUE = round(P_VALUE, 4))])
}


# -----------------------------------------------------------------------------
# 4.3  Band contrasts and rings interacted with shoppability
# -----------------------------------------------------------------------------
# A claim of decay is a claim about differences between bands, which 4.2 does
# not test. If the selection bias that attenuates OLS is similar across
# bands, it cancels in those differences, so the contrasts are less exposed
# to it than the band coefficients. The interacted model asks whether the
# shoppability gradient itself is local: each band count is interacted with
# the Scheme 1 shoppable and non-shoppable indicators (six terms). Both models
# are OLS; Part 16 builds the three excluded instruments that a banded IV
# needs and compares that IV with OLS.
#
# md_lincomb() returns a linear combination of coefficients with its
# delta-method standard error and stops if a term is missing from the model;
# md_run_contrasts() applies it to a named list of weight vectors. The block
# stops if the levels of MD_SCHEME are not MD_SHOP_LEVEL and
# MD_NONSHOP_LEVEL. Writes MD04B_ring_band_contrasts.csv,
# MD04C_ring_shop_coefficients.csv, and MD04D_ring_shop_contrasts.csv.

MD_SHOP_LEVEL    <- "Shoppable"
MD_NONSHOP_LEVEL <- "Non_shoppable"

md_lincomb <- function(fit, w) {
  b <- coef(fit); V <- vcov(fit)
  stopifnot("contrast term missing from the model" = all(names(w) %in% names(b)))
  est <- sum(w * b[names(w)])
  se  <- sqrt(drop(t(w) %*% V[names(w), names(w)] %*% w))
  data.table(ESTIMATE = est, SE = se, T_STAT = est / se,
             P_VALUE = .pval(est / se, fit))
}

md_run_contrasts <- function(fit, contrasts) {
  rbindlist(lapply(names(contrasts), function(nm)
    cbind(CONTRAST = nm, md_lincomb(fit, contrasts[[nm]]))))
}

if (!is.null(md_decay_fit)) {
  md_decay_diff <- md_run_contrasts(md_decay_fit, list(
    "0-15 minus 15-30" = c(RING_PRIOR_0_15 = 1, RING_PRIOR_15_30 = -1),
    "0-15 minus 30-60" = c(RING_PRIOR_0_15 = 1, RING_PRIOR_30_60 = -1)))
  save_csv(md_decay_diff, "MD04B_ring_band_contrasts.csv")
  cat("\n=== NEAR MINUS FAR, POOLED (log points per poster) ===\n")
  print(md_decay_diff[, .(CONTRAST, ESTIMATE = round(ESTIMATE, 6),
                          SE = round(SE, 6), P_VALUE = round(P_VALUE, 4))])
}

md_ring_cols <- unique(c(
  PRIMARY_OUTCOME, MD_SCHEME,
  "RING_PRIOR_0_15", "RING_PRIOR_15_30", "RING_PRIOR_30_60",
  available_columns(md_county_panel,
                    c(BASELINE_CONTROLS, BASELINE_FIXED_EFFECTS, BASELINE_CLUSTERS))))

md_ring_i <- md_county_panel[is.finite(RING_PRIOR_0_15) & is.finite(RING_PRIOR_15_30) &
                               is.finite(RING_PRIOR_30_60) & !is.na(get(MD_SCHEME)),
                             ..md_ring_cols]

md_scheme_levels <- sort(unique(md_ring_i[[MD_SCHEME]]))
cat("Levels of", MD_SCHEME, ":", paste(md_scheme_levels, collapse = " | "), "\n")
stopifnot("scheme labels differ -- set MD_SHOP_LEVEL / MD_NONSHOP_LEVEL to the printed levels" =
            setequal(md_scheme_levels, c(MD_SHOP_LEVEL, MD_NONSHOP_LEVEL)))

md_ring_i[, SHOP := as.integer(get(MD_SCHEME) == MD_SHOP_LEVEL)]
for (bd in c("0_15", "15_30", "30_60")) {
  src <- paste0("RING_PRIOR_", bd)
  md_ring_i[, (paste0("RS_", bd)) := get(src) * SHOP]
  md_ring_i[, (paste0("RN_", bd)) := get(src) * (1L - SHOP)]
}
RING_I_TERMS <- c("RS_0_15", "RS_15_30", "RS_30_60", "RN_0_15", "RN_15_30", "RN_30_60")

md_ring_i_fit <- tryCatch(
  feols(
    build_ols_formula(PRIMARY_OUTCOME,
                      c(RING_I_TERMS, available_columns(md_ring_i, BASELINE_CONTROLS)),
                      available_columns(md_ring_i, BASELINE_FIXED_EFFECTS)),
    data    = md_ring_i,
    cluster = build_cluster_formula(available_columns(md_ring_i, BASELINE_CLUSTERS)),
    warn = FALSE, notes = FALSE),
  error = function(e) {
    cat("  Interacted ring model failed:", conditionMessage(e), "\n"); NULL })

if (!is.null(md_ring_i_fit)) {
  md_ring_i_coefs <- data.table(
    TERM = RING_I_TERMS,
    ARM  = rep(c(MD_SHOP_LEVEL, MD_NONSHOP_LEVEL), each = 3L),
    BAND = rep(c("0-15mi", "15-30mi", "30-60mi"), 2L))
  md_ring_i_coefs[, `:=`(COEF = unname(coef(md_ring_i_fit)[TERM]),
                         SE   = unname(sqrt(diag(vcov(md_ring_i_fit))[TERM])))]
  md_ring_i_coefs[, P_VALUE := .pval(COEF / SE, md_ring_i_fit)]
  
  md_ring_i_diff <- md_run_contrasts(md_ring_i_fit, list(
    "Gap at 0-15 (shoppable minus non-shoppable)" = c(RS_0_15 = 1, RN_0_15 = -1),
    "Gap at 15-30"                   = c(RS_15_30 = 1, RN_15_30 = -1),
    "Gap at 30-60"                   = c(RS_30_60 = 1, RN_30_60 = -1),
    "Shoppable, 0-15 minus 15-30"    = c(RS_0_15 = 1, RS_15_30 = -1),
    "Shoppable, 0-15 minus 30-60"    = c(RS_0_15 = 1, RS_30_60 = -1),
    "Gap at 0-15 minus gap at 15-30" = c(RS_0_15 = 1, RN_0_15 = -1,
                                         RS_15_30 = -1, RN_15_30 = 1),
    "Gap at 0-15 minus gap at 30-60" = c(RS_0_15 = 1, RN_0_15 = -1,
                                         RS_30_60 = -1, RN_30_60 = 1)))
  
  save_csv(md_ring_i_coefs, "MD04C_ring_shop_coefficients.csv")
  save_csv(md_ring_i_diff,  "MD04D_ring_shop_contrasts.csv")
  
  cat("\n=== RINGS BY SHOPPABILITY (OLS, log points per poster) ===\n")
  print(md_ring_i_coefs[, .(ARM, BAND, COEF = round(COEF, 6), SE = round(SE, 6),
                            P_VALUE = round(P_VALUE, 4))])
  cat("\n=== CONTRASTS ===\n")
  print(md_ring_i_diff[, .(CONTRAST, ESTIMATE = round(ESTIMATE, 6),
                           SE = round(SE, 6), P_VALUE = round(P_VALUE, 4))])
  cat("N:", format(nobs(md_ring_i_fit), big.mark = ","), "\n")
}
rm(md_ring_i); invisible(gc())


# =============================================================================
# Part 5: Own-county against adjacent-county decomposition
# =============================================================================
#
# OLS of the outcome on the own-county and adjacent-county prior-poster counts
# (ADJ_PRIOR_OWN, ADJ_PRIOR_ADJACENT) as two regressors, with the baseline
# controls, fixed effects, and clusters, followed by a Wald test that the two
# coefficients are equal (wald_equality(), pipeline Section 1); the printed
# heading calls it the reduced form. A near-zero adjacent-county coefficient
# supports the county boundary; a large one measures the cross-border
# channel. Writes MD05_adjacency_decomposition.csv to TABLE_DIR and
# MD06_adjacency_equality_test.csv to QA_DIR.
#
# The two counts are not collinear: the own-county mean is 2.343 and the
# adjacent-county mean 5.711, with a correlation of 0.468, low enough to
# identify both coefficients.

md_adj_panel <- md_county_panel[is.finite(ADJ_PRIOR_OWN) & is.finite(ADJ_PRIOR_ADJACENT)]

md_adj_fit <- tryCatch(
  feols(
    build_ols_formula(
      PRIMARY_OUTCOME,
      c("ADJ_PRIOR_OWN", "ADJ_PRIOR_ADJACENT",
        available_columns(md_adj_panel, BASELINE_CONTROLS)),
      available_columns(md_adj_panel, BASELINE_FIXED_EFFECTS)),
    data    = md_adj_panel,
    cluster = build_cluster_formula(
      available_columns(md_adj_panel, BASELINE_CLUSTERS)),
    warn = FALSE, notes = FALSE),
  error = function(e) { cat("  Adjacency model failed:", conditionMessage(e), "\n"); NULL })

if (!is.null(md_adj_fit)) {
  md_adj_out <- data.table(
    CHANNEL = c("Own county", "Adjacent counties"),
    TERM    = c("ADJ_PRIOR_OWN", "ADJ_PRIOR_ADJACENT"))
  md_adj_out[, COEF := vapply(TERM, function(t) unname(coef(md_adj_fit)[t]), numeric(1))]
  md_adj_out[, SE   := vapply(TERM, function(t) unname(sqrt(vcov(md_adj_fit)[t, t])), numeric(1))]
  md_adj_out[, `:=`(P_VALUE = .pval(COEF / SE, md_adj_fit),
                    N_OBSERVATIONS = nobs(md_adj_fit))]
  save_csv(md_adj_out, "MD05_adjacency_decomposition.csv")
  cat("\n=== ADJACENCY DECOMPOSITION (reduced form) ===\n")
  print(md_adj_out[, .(CHANNEL, COEF = round(COEF, 6), SE = round(SE, 6),
                       P_VALUE = round(P_VALUE, 4))])
  
  eq_test <- wald_equality(md_adj_fit, c("ADJ_PRIOR_OWN", "ADJ_PRIOR_ADJACENT"))
  if (nrow(eq_test) > 0L) {
    save_qa_csv(eq_test, "MD06_adjacency_equality_test.csv")
    cat("\nTest that own == adjacent:\n"); print(eq_test)
  }
}


# =============================================================================
# Part 6: IV15 against IV06
# =============================================================================
#
# The headline specification on the county panel, estimated once with IV06
# (PRIMARY_INSTRUMENT) and once with IV15
# (Competitor_excl_adjacent_hospitals_9m), each on its own full sample. Writes
# the coefficient rows to MD07_iv15_vs_iv06_rows.csv.
#
# The adjacency filter removes part of the peer set but leaves the instrument
# largely intact: IV15 drops 19.6% of IV06's peers, moving the mean from
# 7.168 to 5.763, and the two correlate at 0.9806.
#
# IV15 is NaN for 0.78% of hospital-months, almost all in Connecticut, where
# the 2022 planning-region reorganization leaves no usable county FIPS code.

md_iv_compare <- rbindlist(lapply(
  c(Competitor_only_hospitals_9m = PRIMARY_INSTRUMENT,
    IV_EXCLUDE_ADJACENT["Competitor_excl_adjacent_hospitals_9m"]),
  function(iv_col) {
    iv_name <- names(which(c(
      Competitor_only_hospitals_9m = PRIMARY_INSTRUMENT,
      IV_EXCLUDE_ADJACENT) == iv_col))[1]
    cat("\n--- Instrument:", iv_name, "---\n")
    
    if (!(iv_col %in% names(md_county_panel))) {
      cat("  ABSENT from panel. Check the ANALYSIS_COLUMNS extension in Part 0.\n")
      return(NULL)
    }
    
    res <- estimate_interacted(
      data             = md_county_panel,
      moderator        = MD_SCHEME,
      moderator_type   = "categorical",
      outcome          = PRIMARY_OUTCOME,
      instrument       = iv_col,
      label            = paste0("IVcheck=", iv_name),
      instrument_label = iv_name,
      moderator_label  = MD_SCHEME)
    
    if (is.null(res)) return(NULL)
    r <- copy(res$rows); r[, INSTRUMENT_TESTED := iv_name]; r
  }), fill = TRUE)

if (nrow(md_iv_compare) > 0L) {
  save_csv(md_iv_compare, "MD07_iv15_vs_iv06_rows.csv")
  cat("\n=== IV06 VS IV15 (headline spec, county market) ===\n")
  print(md_iv_compare[, .(INSTRUMENT_TESTED, TERM,
                          IV_PERCENT = round(IV_PERCENT, 3),
                          IV_P = round(IV_P, 4),
                          FIRST_STAGE_WALD_MIN = round(FIRST_STAGE_WALD_MIN, 1),
                          N_OBSERVATIONS)])
}


# -----------------------------------------------------------------------------
# 6.1  The same comparison in reduced form
# -----------------------------------------------------------------------------
# Part 6 keeps only the coefficient rows and prints the IV columns. The other
# robustness checks in the paper are read off the reduced form and its
# equality test, so this block re-estimates the same two models, caches them
# (md_iv15_vs_iv06_rf_2inst), and writes the rows
# (MD07B_iv15_vs_iv06_rf_rows.csv) and the equality tests
# (MD07C_iv15_vs_iv06_tests.csv). The calls match Part 6, so IV_PERCENT
# reproduces MD07 row for row, which checks the two blocks against each
# other.

md_iv_compare_rf <- cache_or_run("md_iv15_vs_iv06_rf_2inst", {
  ivs <- c(Competitor_only_hospitals_9m = PRIMARY_INSTRUMENT,
           Competitor_excl_adjacent_hospitals_9m =
             unname(IV_EXCLUDE_ADJACENT["Competitor_excl_adjacent_hospitals_9m"]))
  out <- lapply(names(ivs), function(nm) {
    res <- estimate_interacted(
      data = md_county_panel, moderator = MD_SCHEME, moderator_type = "categorical",
      outcome = PRIMARY_OUTCOME, instrument = ivs[[nm]],
      label = paste0("IVcheckRF=", nm), instrument_label = nm,
      moderator_label = MD_SCHEME)
    if (is.null(res)) return(NULL)
    list(rows  = cbind(INSTRUMENT_TESTED = nm, res$rows),
         tests = cbind(INSTRUMENT_TESTED = nm, res$tests))
  })
  list(rows  = rbindlist(lapply(out, `[[`, "rows"),  fill = TRUE),
       tests = rbindlist(lapply(out, `[[`, "tests"), fill = TRUE))
})

save_csv(md_iv_compare_rf$rows,  "MD07B_iv15_vs_iv06_rf_rows.csv")
save_csv(md_iv_compare_rf$tests, "MD07C_iv15_vs_iv06_tests.csv")

cat("\n=== PANEL B IN REDUCED FORM ===\n")
print(md_iv_compare_rf$rows[, .(INSTRUMENT_TESTED, TERM,
                                RF_PCT_PER_SD = round(RF_PERCENT_PER_SD, 2),
                                RF_P = round(RF_P, 4),
                                IV_PERCENT = round(IV_PERCENT, 2),
                                FIRST_STAGE_WALD_MIN = round(FIRST_STAGE_WALD_MIN, 1),
                                N_OBSERVATIONS)])
print(md_iv_compare_rf$tests[, .(INSTRUMENT_TESTED, ESTIMATOR,
                                 P_VALUE = round(P_VALUE, 4))])


# =============================================================================
# Part 7: Rural/urban and multi-hospital subsamples
# =============================================================================
#
# The headline specification on the county panel in five samples: full,
# metro, nonmetro (METRO_STATUS), and markets with at least two or at least
# three hospitals. Writes MD08_subsample_rows.csv. Both restrictions leave
# estimable samples: 2,551 metro against 1,138 nonmetro hospitals, and 2,541
# of 3,724 hospitals (68.2%) are in counties with at least two hospitals.
#
# The multi-hospital restriction bears most directly on the market
# definition. In a single-hospital county there is no local competitor for
# the treatment to vary against, so those counties contribute level
# differences but no within-market identifying variation.
#
# The number of hospitals per market is counted in md_county_panel and merged
# onto it as N_HOSP_IN_MARKET. Running this part a second time in the same
# session duplicates that column (N_HOSP_IN_MARKET.x and .y), and the two
# hospital-count filters then fail.

md_market_sizes <- unique(md_county_panel, by = "HOSPITAL_ID")[
  !is.na(ANALYSIS_MARKET), .(N_HOSP_IN_MARKET = .N), by = ANALYSIS_MARKET]
md_county_panel <- merge(md_county_panel, md_market_sizes,
                         by = "ANALYSIS_MARKET", all.x = TRUE)

md_subsamples <- list(
  Full            = function(d) d,
  Metro           = function(d) d[METRO_STATUS == "Metro"],
  Nonmetro        = function(d) d[METRO_STATUS == "Nonmetro"],
  Multi_hospital  = function(d) d[N_HOSP_IN_MARKET >= 2],
  Three_plus_hosp = function(d) d[N_HOSP_IN_MARKET >= 3]
)

md_subsample_results <- rbindlist(lapply(names(md_subsamples), function(nm) {
  cat("\n--- Subsample:", nm, "---\n")
  d <- md_subsamples[[nm]](md_county_panel)
  cat("  Rows:", format(nrow(d), big.mark = ","),
      "| hospitals:", uniqueN(d$HOSPITAL_ID),
      "| markets:", uniqueN(d$ANALYSIS_MARKET), "\n")
  
  if (nrow(d) < MIN_MODEL_OBS) {
    cat("  Below MIN_MODEL_OBS (", format(MIN_MODEL_OBS, big.mark = ","), ").\n")
    return(NULL)
  }
  
  res <- estimate_interacted(
    data             = d,
    moderator        = MD_SCHEME,
    moderator_type   = "categorical",
    outcome          = PRIMARY_OUTCOME,
    instrument       = PRIMARY_INSTRUMENT,
    label            = paste0("Subsample=", nm),
    instrument_label = "Competitor_only_hospitals_9m",
    moderator_label  = MD_SCHEME)
  
  if (is.null(res)) return(NULL)
  r <- copy(res$rows)
  r[, `:=`(SUBSAMPLE = nm,
           N_HOSPITALS = uniqueN(d$HOSPITAL_ID),
           N_MARKETS   = uniqueN(d$ANALYSIS_MARKET),
           PANEL_ROWS  = nrow(d))]
  r
}), fill = TRUE)

if (nrow(md_subsample_results) > 0L) {
  save_csv(md_subsample_results, "MD08_subsample_rows.csv")
  cat("\n=== SUBSAMPLES ===\n")
  print(md_subsample_results[, .(SUBSAMPLE, TERM, N_HOSPITALS, N_MARKETS,
                                 IV_PERCENT = round(IV_PERCENT, 3),
                                 IV_P = round(IV_P, 4),
                                 N_OBSERVATIONS)])
}


# =============================================================================
# Part 9: Concept-level analysis at HRR, one instrument
# =============================================================================
#
# Runs the pipeline's concept-level estimator (estimate_concept_level(),
# Section 6) and meta-regressions (run_meta_regressions(), Section 8) on the
# HRR panel with one instrument, Competitor_only_hospitals_9m
# (CONCEPT_INSTRUMENTS_MD), and compares the concept estimates with the
# county run.
#
# HSA is not run here. In a 30-concept sample a median 73.4% of MARKET_ID
# fixed-effect cells at HSA are singletons, against 20.1% at HRR, so this
# panel does not support per-concept estimation at HSA. Part 12 estimates HSA
# on a restricted sample.
#
# One instrument instead of six. estimate_concept_level() fits a reduced
# form, a first stage, and an IV for each instrument plus one OLS, so six
# instruments mean nineteen fits per concept and one instrument four.
# Benchmarks with system.time() on this fixed-effect structure put three
# instruments at roughly 28-30 hours for 738 concepts and one instrument at
# roughly 11-12 hours. The comparison of IV06 with IV15 in Part 6 does not
# depend on this choice, and Part 15 runs six instruments at HRR and HSA.
#
# md_run_concept_level() loads the panel, then runs the concept loop and the
# reduced-form and IV meta-regressions under the cache keys
# md_concept_level_<geo>_<n>inst, md_meta_rf_<geo>_<n>inst, and
# md_meta_iv_<geo>_<n>inst, and writes MD13_concept_level, MD14_meta_rf, and
# MD15_meta_iv with the suffix _<geo>_<n>inst. The keys record the number of
# instruments, not their names: a run with a different single instrument
# loads the cached Competitor_only_hospitals_9m results unless those cache
# files are deleted.

CONCEPT_INSTRUMENTS_MD <- c(Competitor_only_hospitals_9m = PRIMARY_INSTRUMENT)

md_key <- function(stem, level, instruments = CONCEPT_INSTRUMENTS_MD) {
  sprintf("%s_%s_%dinst", stem, tolower(level), length(instruments))
}

md_run_concept_level <- function(level, instruments = CONCEPT_INSTRUMENTS_MD) {
  cat("\n", strrep("=", 70), "\n", sep = "")
  cat("CONCEPT-LEVEL:", level, "| instruments:", length(instruments), "\n")
  cat("  ", paste(names(instruments), collapse = "\n   "), "\n", sep = "")
  cat("  started:", format(Sys.time(), "%H:%M:%S"), "\n")
  cat(strrep("=", 70), "\n", sep = "")
  
  t_start <- Sys.time()
  p <- md_load_panel(level)
  
  cr <- cache_or_run(
    md_key("md_concept_level", level, instruments),
    estimate_concept_level(
      build_concept_panel(p, instruments = instruments),
      instruments = instruments,
      save_stem   = paste0("MD13_concept_level_", tolower(level),
                           "_", length(instruments), "inst")
    ))
  
  cat("\nConcept-level done for", level, "--",
      format(nrow(cr), big.mark = ","), "rows,",
      format(uniqueN(cr$FINAL_CONCEPT_ID), big.mark = ","), "concepts |",
      sprintf("%.1f min\n", as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
  
  cr[, `:=`(MARKET_DEFINITION = level, N_INSTRUMENTS_USED = length(instruments))]
  
  meta_rf_g <- cache_or_run(md_key("md_meta_rf", level, instruments),
                            run_meta_regressions(prepare_meta_input(cr, md_schemes_long),
                                                 md_schemes_long, dep = "RF_COEF", se = "RF_SE"))
  meta_iv_g <- cache_or_run(md_key("md_meta_iv", level, instruments),
                            run_meta_regressions(prepare_meta_input(cr, md_schemes_long),
                                                 md_schemes_long, dep = "IV_COEF", se = "IV_SE"))
  
  meta_rf_g[, `:=`(MARKET_DEFINITION = level, N_INSTRUMENTS_USED = length(instruments))]
  meta_iv_g[, `:=`(MARKET_DEFINITION = level, N_INSTRUMENTS_USED = length(instruments))]
  
  suffix <- sprintf("_%s_%dinst.csv", tolower(level), length(instruments))
  save_csv(cr,        paste0("MD13_concept_level", suffix))
  save_csv(meta_rf_g, paste0("MD14_meta_rf",       suffix))
  save_csv(meta_iv_g, paste0("MD15_meta_iv",       suffix))
  
  rm(p); invisible(gc())
  list(concept = cr, meta_rf = meta_rf_g, meta_iv = meta_iv_g)
}


# -----------------------------------------------------------------------------
# 9A  HRR concept level
# -----------------------------------------------------------------------------
# Roughly 11-12 hours: 738 concepts at four fits each. Progress is saved to
# MD13_concept_level_hrr_1inst_PARTIAL.csv every 50 concepts, and a completed
# run is loaded from the cache.
md_hrr <- md_run_concept_level("HRR")

cat("\n=== HRR META-REGRESSION, IV (all schemes) ===\n")
print(md_hrr$meta_iv)

cat("\n=== HRR META-REGRESSION, REDUCED FORM (all schemes) ===\n")
print(md_hrr$meta_rf)


# -----------------------------------------------------------------------------
# 9B  Summary of the HRR concept-level estimates
# -----------------------------------------------------------------------------
# Number of concepts, median IV estimate in percent, share of negative IV
# coefficients, share significant at 5%, and median first-stage F
# (MD19_hrr_concept_summary_1inst.csv).
md_hrr_summary <- md_hrr$concept[, .(
  N_CONCEPTS      = .N,
  MEDIAN_IV_PCT   = round(median(IV_ESTIMATE_PERCENT, na.rm = TRUE), 3),
  SHARE_NEGATIVE  = round(mean(IV_COEF < 0, na.rm = TRUE), 3),
  SHARE_SIG_5PCT  = round(mean(IV_P < 0.05, na.rm = TRUE), 3),
  MEDIAN_FS_F     = round(median(FS_F, na.rm = TRUE), 1)
), by = INSTRUMENT_LABEL]

save_csv(md_hrr_summary, "MD19_hrr_concept_summary_1inst.csv")
cat("\n=== HRR CONCEPT-LEVEL SUMMARY ===\n")
print(md_hrr_summary)


# -----------------------------------------------------------------------------
# 9C  County against HRR, concept by concept
# -----------------------------------------------------------------------------
# The county side is concept_results, the pipeline's cached six-instrument
# county run, restricted to Competitor_only_hospitals_9m, the instrument used
# at HRR; the block stops if concept_results is not in memory. The delta is
# the county IV coefficient minus the HRR coefficient
# (IV_DELTA_LOCAL_MINUS_HRR), and MD16 holds it for every matched concept.
# Concepts with a first-stage F of at least 10 on both sides are then
# summarized by category under each of the 18 schemes in md_schemes_long
# (MD17), by the slope of the delta on the ordinal shoppability rank (MD17B),
# and by clinical family (MD18). SHARE_LOCAL_STRONGER is the share of
# negative deltas.
if (!exists("concept_results"))
  stop("concept_results (cached county run) not in scope. ",
       "Re-run warm_start().", call. = FALSE)

md_delta <- merge(
  as.data.table(concept_results)[
    INSTRUMENT_LABEL == "Competitor_only_hospitals_9m",
    .(FINAL_CONCEPT_ID, INSTRUMENT_LABEL, FINAL_FAMILY_ID, SERVICE_LABEL,
      IV_LOCAL = IV_COEF, IV_SE_LOCAL = IV_SE, IV_P_LOCAL = IV_P,
      FS_F_LOCAL = FS_F, N_OBS_LOCAL = N_OBSERVATIONS)],
  md_hrr$concept[
    , .(FINAL_CONCEPT_ID, INSTRUMENT_LABEL,
        IV_HRR = IV_COEF, IV_SE_HRR = IV_SE, IV_P_HRR = IV_P,
        FS_F_HRR = FS_F, N_OBS_HRR = N_OBSERVATIONS)],
  by = c("FINAL_CONCEPT_ID", "INSTRUMENT_LABEL"))

md_delta[, IV_DELTA_LOCAL_MINUS_HRR := IV_LOCAL - IV_HRR]

md_delta_ok <- md_delta[is.finite(IV_DELTA_LOCAL_MINUS_HRR) &
                          FS_F_LOCAL >= 10 & FS_F_HRR >= 10]

cat(sprintf("\nConcepts matched: %d | F>=10 both sides: %d\n",
            nrow(md_delta), nrow(md_delta_ok)))

md_delta_schemes <- copy(md_delta_ok)
for (sid in unique(md_schemes_long$SCHEME_ID)) {
  a <- unique(md_schemes_long[SCHEME_ID == sid,
                              .(FINAL_CONCEPT_ID = ANALYSIS_CONCEPT_ID,
                                CAT = SHOPPABILITY_CATEGORY,
                                ORD = SHOPPABILITY_ORDINAL)])
  setnames(a, c("CAT", "ORD"), paste0(c("CAT_", "ORD_"), sid))
  md_delta_schemes <- merge(md_delta_schemes, a, by = "FINAL_CONCEPT_ID", all.x = TRUE)
}

md_crossover_test <- rbindlist(lapply(unique(md_schemes_long$SCHEME_ID), function(sid) {
  ccol <- paste0("CAT_", sid)
  if (!(ccol %in% names(md_delta_schemes))) return(NULL)
  d <- md_delta_schemes[!is.na(get(ccol))]
  if (uniqueN(d[[ccol]]) < 2L || nrow(d) < MIN_CONCEPTS_META) return(NULL)
  
  d[, .(SCHEME_ID = sid, N_CONCEPTS = .N,
        MEAN_DELTA   = round(mean(IV_DELTA_LOCAL_MINUS_HRR, na.rm = TRUE), 4),
        MEDIAN_DELTA = round(median(IV_DELTA_LOCAL_MINUS_HRR, na.rm = TRUE), 4),
        SHARE_LOCAL_STRONGER = round(mean(IV_DELTA_LOCAL_MINUS_HRR < 0, na.rm = TRUE), 3)),
    by = c(CATEGORY = ccol)]
}), fill = TRUE)

md_crossover_slope <- rbindlist(lapply(unique(md_schemes_long$SCHEME_ID), function(sid) {
  ocol <- paste0("ORD_", sid)
  if (!(ocol %in% names(md_delta_schemes))) return(NULL)
  d <- md_delta_schemes[is.finite(get(ocol)) & is.finite(IV_DELTA_LOCAL_MINUS_HRR)]
  if (nrow(d) < MIN_CONCEPTS_META || !has_usable_variation(d[[ocol]])) return(NULL)
  
  y <- d$IV_DELTA_LOCAL_MINUS_HRR; x <- d[[ocol]]
  fit <- tryCatch(lm(y ~ x), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  s <- summary(fit)$coefficients
  if (nrow(s) < 2L) return(NULL)
  
  data.table(SCHEME_ID = sid, N_CONCEPTS = nrow(d),
             SLOPE = round(s[2, 1], 5), SE = round(s[2, 2], 5),
             P_VALUE = round(s[2, 4], 4))
}), fill = TRUE)

save_csv(md_delta,           "MD16_concept_local_hrr_delta_1inst.csv")
save_csv(md_crossover_test,  "MD17_crossover_by_scheme_1inst.csv")
save_csv(md_crossover_slope, "MD17B_crossover_slope_by_scheme_1inst.csv")

cat("\n=== CROSS-OVER: mean COUNTY-minus-HRR delta by scheme category ===\n")
print(md_crossover_test)
cat("\n=== CROSS-OVER SLOPE ===\n")
print(md_crossover_slope)

md_family_delta <- md_delta_ok[, .(
  N_PAIRS      = .N,
  MEDIAN_DELTA = round(median(IV_DELTA_LOCAL_MINUS_HRR, na.rm = TRUE), 4),
  SHARE_LOCAL_STRONGER = round(mean(IV_DELTA_LOCAL_MINUS_HRR < 0, na.rm = TRUE), 3)
), by = FINAL_FAMILY_ID][order(MEDIAN_DELTA)]

save_csv(md_family_delta, "MD18_family_local_hrr_delta_1inst.csv")
cat("\n=== BY CLINICAL FAMILY ===\n")
print(md_family_delta)


# =============================================================================
# Part 8: Assembled summary of Parts 3 to 7
# =============================================================================
#
# Stacks the IV rows of the ladder (Part 3), the ring treatment (4.1), the
# IV06 against IV15 comparison (Part 6), and the subsamples (Part 7) into
# MD09_committee_response_summary.csv. Part 8 sits after Part 9 in the file
# and runs after it, but uses only results from Parts 3-7. ring_rows enters
# only if it exists, and md_iv_compare and md_subsample_results only if they
# have rows. The "Shoppability x market" block refers to md_crossover, which
# is not defined in this file or in the pipeline, so its exists() guard
# skips it.

md_summary <- rbindlist(list(
  md_ladder$rows[, .(BLOCK = "Market definition", VARIANT = MARKET_DEFINITION,
                     TERM, IV_PERCENT, IV_P, FIRST_STAGE_WALD_MIN, N_OBSERVATIONS)],
  if (exists("ring_rows"))
    ring_rows[, .(BLOCK = "Distance treatment", VARIANT = MARKET_DEFINITION,
                  TERM, IV_PERCENT, IV_P, FIRST_STAGE_WALD_MIN, N_OBSERVATIONS)],
  if (nrow(md_iv_compare) > 0L)
    md_iv_compare[, .(BLOCK = "Instrument", VARIANT = INSTRUMENT_TESTED,
                      TERM, IV_PERCENT, IV_P, FIRST_STAGE_WALD_MIN, N_OBSERVATIONS)],
  if (exists("md_crossover") && nrow(md_crossover) > 0L)
    md_crossover[, .(BLOCK = "Shoppability x market", 
                     VARIANT = paste(MARKET_DEFINITION, SHOPPABILITY_TIER, sep = " / "),
                     TERM = SHOPPABILITY_TIER, IV_PERCENT, IV_P,
                     FIRST_STAGE_WALD_MIN = FIRST_STAGE_WALD, N_OBSERVATIONS)],
  if (nrow(md_subsample_results) > 0L)
    md_subsample_results[, .(BLOCK = "Subsample", VARIANT = SUBSAMPLE,
                             TERM, IV_PERCENT, IV_P, FIRST_STAGE_WALD_MIN, N_OBSERVATIONS)]
), fill = TRUE)

save_csv(md_summary, "MD09_committee_response_summary.csv")

cat("\n\n=== COMMITTEE RESPONSE SUMMARY ===\n")
print(md_summary[, .(BLOCK, VARIANT, TERM,
                     IV_PERCENT = round(IV_PERCENT, 3),
                     IV_P = round(IV_P, 4),
                     N_OBSERVATIONS)])

cat("\nWritten to:", TABLE_DIR, "\n")
cat("  MD01/MD02  market definition ladder\n")
cat("  MD03       ring treatment\n")
cat("  MD04       spatial decay diagnostic\n")
cat("  MD05/MD06  adjacency decomposition\n")
cat("  MD07       IV15 vs IV06\n")
cat("  MD08       subsamples\n")
cat("  MD09       assembled summary\n")
cat("\n=== DONE ===\n")


# =============================================================================
# Part 10: Competitor instruments at HRR
# =============================================================================
#
# Parts 10 and 11 rebuild the Python pipeline's competitor instrument in R,
# check the rebuild against the instrument columns in the county panel, and
# use the same builder to construct instrument sets at HRR and HSA, which the
# Python pipeline does not build. 10.3 and 11.2 stop unless every tested
# rebuild matches its panel column exactly in at least 98% of hospital-months.


# -----------------------------------------------------------------------------
# 10.1  Hospital roster
# -----------------------------------------------------------------------------
# Peers come from HPT_HOSPITAL_DISCLOSURE_EVENTS.parquet in INSTRUMENT_DIR,
# which the Python pipeline assembles from all three exact-code frames
# (outpatient, inpatient, and component). It covers 3,805 hospitals, while
# md_county_panel (outpatient concepts only) covers 3,723. A roster built from
# the panel would omit inpatient-only and component-only hospitals, which are
# valid peers, and would understate every instrument. md_event_roster has one
# row per hospital with the first non-missing SYSTEM_KEY, posting month
# (PEER_POST_MONTH), COUNTY_STATE_KEY, and COUNTY_FIPS.
INSTRUMENT_DIR <- file.path(FINAL_DATA_ROOT, "02_Treatments_and_Instruments")

first_non_na <- function(x) {
  v <- x[!is.na(x)]
  if (length(v)) v[1] else x[NA_integer_][1]
}

md_events <- as.data.table(read_parquet(
  file.path(INSTRUMENT_DIR, "HPT_HOSPITAL_DISCLOSURE_EVENTS.parquet")))
md_events[, HOSPITAL_ID := as.character(HOSPITAL_ID)]

md_event_roster <- md_events[, .(
  SYSTEM_KEY       = first_non_na(SYSTEM_KEY),
  PEER_POST_MONTH  = as.Date(first_non_na(HOSPITAL_FIRST_POST_MONTH)),
  COUNTY_STATE_KEY = first_non_na(COUNTY_STATE_KEY),
  COUNTY_FIPS      = first_non_na(COUNTY_FIPS)
), by = HOSPITAL_ID]

stopifnot("duplicate HOSPITAL_ID in event roster" =
            !any(duplicated(md_event_roster$HOSPITAL_ID)))

cat(sprintf("Event roster: %d hospitals, %d missing SYSTEM_KEY, %d missing post month\n",
            nrow(md_event_roster), sum(is.na(md_event_roster$SYSTEM_KEY)),
            sum(is.na(md_event_roster$PEER_POST_MONTH))))


# -----------------------------------------------------------------------------
# 10.2  Instrument builder
# -----------------------------------------------------------------------------
# build_dynamic_competitor_instrument() is an R port of the Python pipeline's
# build_competitor_instrument() for any geography; geo_col names the roster
# column that defines it. For each hospital-month in focal_panel it counts
# the peers that meet all of these conditions:
#   * the peer's system has a hospital in the focal hospital's geography that
#     posted no later than the focal month (LOCAL_SYSTEM_FIRST_POST_MONTH <=
#     POST_MONTH);
#   * the peer is in a different geography;
#   * the peer posted in the lookback_months (default 9) before the focal
#     month, not counting the focal month itself;
#   * with exclude_own_system = TRUE, the peer's system differs from the focal
#     hospital's.
# It returns one row per hospital-month with the number of qualifying peer
# hospitals, systems, and geographies in <out_prefix>_HOSPITALS, _SYSTEMS, and
# _GEOS. The counts are zero where no peer qualifies, including hospitals with
# no geography code.
#
#   exclude_own_system = TRUE    competitor family (the IV06 analog)
#   exclude_own_system = FALSE   strict family (the IV02 analog)
#   adjacency = <pair table>     drops peers in geographies adjacent to the
#                                focal one (columns GEO_A and GEO_B; the IV15
#                                analog); NULL skips the filter
#   presence = "posted"          default; applies the first condition in full
#   presence = "owned"           drops the date in the first condition: a
#                                rival system counts whenever it has a
#                                hospital (with a recorded posting month) in
#                                the focal geography
#   return_pairs = TRUE          returns the qualifying hospital-peer pairs
#                                instead of the counts; no call in this file
#                                uses it
#
# The first condition, the LOCAL_SYSTEM_FIRST_POST_MONTH gate, is the core of
# the construction: without it the rebuilt counts are off by a factor of 30
# to 1700. 10.3 and 11.2 validate the defaults (presence = "posted",
# return_pairs = FALSE), and Part 17 uses presence = "owned". The
# distance-band instruments in Part 16 come from a separate function,
# build_ring_instruments(), that applies the same conditions.
build_dynamic_competitor_instrument <- function(focal_panel, roster, geo_col,
                                                out_prefix, lookback_months = 9,
                                                exclude_own_system = TRUE,
                                                adjacency = NULL,
                                                presence = c("posted", "owned"),
                                                return_pairs = FALSE) {
  
  presence <- match.arg(presence)
  
  r <- copy(roster); setnames(r, geo_col, "GEO")
  r <- r[, .(HOSPITAL_ID, SYSTEM_KEY, PEER_POST_MONTH, GEO)]
  stopifnot("duplicate HOSPITAL_ID in roster" = !any(duplicated(r$HOSPITAL_ID)))
  
  cat(sprintf("  [%s] roster %d hospitals, %d missing SYSTEM_KEY, %d missing geo\n",
              out_prefix, nrow(r), sum(is.na(r$SYSTEM_KEY)), sum(is.na(r$GEO))))
  
  local_system_roster <- r[!is.na(GEO) & !is.na(SYSTEM_KEY) & !is.na(PEER_POST_MONTH),
                           .(LOCAL_SYSTEM_FIRST_POST_MONTH = min(PEER_POST_MONTH)),
                           by = .(FOCAL_GEO = GEO, SYSTEM_KEY)]
  
  peer_events <- r[!is.na(SYSTEM_KEY) & !is.na(GEO) & !is.na(PEER_POST_MONTH),
                   .(PEER_HOSPITAL_ID = HOSPITAL_ID, SYSTEM_KEY,
                     PEER_GEO = GEO, PEER_POST_MONTH)]
  
  peer_universe <- merge(local_system_roster, peer_events, by = "SYSTEM_KEY",
                         allow.cartesian = TRUE)[FOCAL_GEO != PEER_GEO]
  
  if (!is.null(adjacency)) {
    n_pre <- nrow(peer_universe)
    peer_universe <- peer_universe[!adjacency,
                                   on = .(FOCAL_GEO = GEO_A, PEER_GEO = GEO_B)]
    cat(sprintf("  [%s] adjacency filter %s to %s peer rows (%.1f%% dropped)\n",
                out_prefix, format(n_pre, big.mark = ","),
                format(nrow(peer_universe), big.mark = ","),
                100 * (1 - nrow(peer_universe) / n_pre)))
  }
  
  cat(sprintf("  [%s] peer universe %s rows, %d focal geographies, %d systems\n",
              out_prefix, format(nrow(peer_universe), big.mark = ","),
              uniqueN(peer_universe$FOCAL_GEO), uniqueN(peer_universe$SYSTEM_KEY)))
  
  focal_events <- unique(focal_panel[, .(HOSPITAL_ID, POST_MONTH)])
  n_before <- nrow(focal_events)
  focal_events <- merge(focal_events, r[, .(HOSPITAL_ID, GEO, SYSTEM_KEY)],
                        by = "HOSPITAL_ID", all.x = TRUE)
  stopifnot("roster merge changed focal row count" = nrow(focal_events) == n_before)
  setnames(focal_events, c("GEO", "SYSTEM_KEY"), c("FOCAL_GEO", "FOCAL_SYSTEM"))
  
  joined <- merge(focal_events, peer_universe, by = "FOCAL_GEO", allow.cartesian = TRUE)
  # presence = "posted" keeps a peer only once a hospital of its system in the
  # focal geography has posted; "owned" drops that condition.
  joined[, DYNAMIC_OK := if (presence == "owned") TRUE else
    LOCAL_SYSTEM_FIRST_POST_MONTH <= POST_MONTH]
  joined[, LOWER := POST_MONTH %m-% months(lookback_months)]
  joined[, STRICT := PEER_POST_MONTH >= LOWER & PEER_POST_MONTH < POST_MONTH]
  joined[, COMPETITOR := if (exclude_own_system)
    is.na(FOCAL_SYSTEM) | (SYSTEM_KEY != FOCAL_SYSTEM) else TRUE]
  
  if (isTRUE(return_pairs))
    return(unique(joined[DYNAMIC_OK & STRICT & COMPETITOR,
                         .(HOSPITAL_ID, POST_MONTH, PEER_HOSPITAL_ID,
                           SYSTEM_KEY, PEER_GEO)]))
  
  counts <- joined[DYNAMIC_OK & STRICT & COMPETITOR,
                   .(N_HOSPITALS = uniqueN(PEER_HOSPITAL_ID),
                     N_SYSTEMS   = uniqueN(SYSTEM_KEY),
                     N_GEOS      = uniqueN(PEER_GEO)),
                   by = .(HOSPITAL_ID, POST_MONTH)]
  
  result <- merge(focal_events[, .(HOSPITAL_ID, POST_MONTH)], counts,
                  by = c("HOSPITAL_ID", "POST_MONTH"), all.x = TRUE)
  result[is.na(N_HOSPITALS), `:=`(N_HOSPITALS = 0L, N_SYSTEMS = 0L, N_GEOS = 0L)]
  stopifnot("duplicate (HOSPITAL_ID, POST_MONTH) in result" =
              !any(duplicated(result, by = c("HOSPITAL_ID", "POST_MONTH"))))
  
  setnames(result, c("N_HOSPITALS", "N_SYSTEMS", "N_GEOS"),
           paste0(out_prefix, c("_HOSPITALS", "_SYSTEMS", "_GEOS")))
  result[]
}


# -----------------------------------------------------------------------------
# 10.3  Validation against the Python-generated columns
# -----------------------------------------------------------------------------
# Rebuilds the county instruments on COUNTY_STATE_KEY with both values of
# exclude_own_system and compares them with the panel's columns: competitor
# hospitals, systems, and counties with Z_SYS_COMPETITOR_ONLY_9M_EXCL_CURRENT,
# Z_SYS_COMPETITOR_SYSTEMS_9M_EXCL_CURRENT, and
# Z_SYS_COMPETITOR_COUNTIES_9M_EXCL_CURRENT, and strict hospitals with
# Z_SYS_STRICT_9M_EXCL_CURRENT. The run stops unless both branches are tested
# and every tested variant matches exactly in at least 98% of hospital-months.
# The adjacency branch is validated in 11.2 against IV15, and 17.2 checks
# that presence = "owned" never yields fewer peers than "posted".
#
# Z_SYS_COMPETITOR_COUNTIES_9M_EXCL_CURRENT is not in
# ALL_CANDIDATE_INSTRUMENTS, so it is not in ANALYSIS_COLUMNS and read_panel()
# does not load it; that check reports "column not loaded" instead of a match
# rate. Its N_GEOS count uses the same logic as N_SYSTEMS and N_HOSPITALS,
# both of which reproduce their columns exactly. Adding the column to
# ANALYSIS_COLUMNS in 0.1 would make it testable.
cty_comp_rebuild <- build_dynamic_competitor_instrument(
  md_county_panel, md_event_roster, "COUNTY_STATE_KEY", "CTY_COMP",
  exclude_own_system = TRUE)

cty_strict_rebuild <- build_dynamic_competitor_instrument(
  md_county_panel, md_event_roster, "COUNTY_STATE_KEY", "CTY_STRICT",
  exclude_own_system = FALSE)

truth_cols <- list(
  T_COMP_H   = "Z_SYS_COMPETITOR_ONLY_9M_EXCL_CURRENT",
  T_COMP_S   = "Z_SYS_COMPETITOR_SYSTEMS_9M_EXCL_CURRENT",
  T_COMP_C   = "Z_SYS_COMPETITOR_COUNTIES_9M_EXCL_CURRENT",
  T_STRICT_H = "Z_SYS_STRICT_9M_EXCL_CURRENT")

avail <- truth_cols[vapply(truth_cols,
                           function(x) x %in% names(md_county_panel), logical(1))]
cat(sprintf("\nValidation columns available: %d of %d\n", length(avail), length(truth_cols)))

truth <- unique(md_county_panel[, c("HOSPITAL_ID", "POST_MONTH", unlist(avail)),
                                with = FALSE])
setnames(truth, unlist(avail), names(avail))

v <- Reduce(function(a, b) merge(a, b, by = c("HOSPITAL_ID", "POST_MONTH")),
            list(cty_comp_rebuild, cty_strict_rebuild, truth))

val_pairs <- list(
  list("competitor hospitals", "CTY_COMP_HOSPITALS",   "T_COMP_H",   "exclude_own_system=TRUE"),
  list("competitor systems",   "CTY_COMP_SYSTEMS",     "T_COMP_S",   "exclude_own_system=TRUE"),
  list("competitor counties",  "CTY_COMP_GEOS",        "T_COMP_C",   "exclude_own_system=TRUE"),
  list("strict hospitals",     "CTY_STRICT_HOSPITALS", "T_STRICT_H", "exclude_own_system=FALSE"))

val_report <- rbindlist(lapply(val_pairs, function(p) {
  if (!(p[[3]] %in% names(v)))
    return(data.table(CHECK = p[[1]], BRANCH = p[[4]], MATCH = NA_real_,
                      MEAN_ABS_DIFF = NA_real_, NOTE = "column not loaded"))
  data.table(CHECK = p[[1]], BRANCH = p[[4]],
             MATCH = round(mean(v[[p[[2]]]] == v[[p[[3]]]], na.rm = TRUE), 4),
             MEAN_ABS_DIFF = round(mean(abs(v[[p[[2]]]] - v[[p[[3]]]]), na.rm = TRUE), 4),
             NOTE = "")
}), fill = TRUE)

cat("\n=== COUNTY VALIDATION ===\n"); print(val_report)

tested <- val_report[!is.na(MATCH)]
stopifnot("only one branch of exclude_own_system was tested" =
            uniqueN(tested$BRANCH) == 2)
stopifnot("a tested variant failed to reproduce" = all(tested$MATCH >= 0.98))
cat(sprintf("Validated: %d variants across %d branches.\n",
            nrow(tested), uniqueN(tested$BRANCH)))


# -----------------------------------------------------------------------------
# 10.4  Six HRR instruments
# -----------------------------------------------------------------------------
# Adds HRR_NUM from md_hosp_geo to the roster and crosses the competitor
# filter (on or off) with the unit counted (peer hospitals, peer systems, or
# peer HRRs). HRR_INSTRUMENTS lists the six columns of hrr_six.
#
# The outside-CBSA variants of the county set are not built. CBSA lies between
# county and HRR in size (929 CBSAs against 306 HRRs), so "outside the focal
# CBSA" is largely implied by "outside the focal HRR". The adjacency variant
# is not built either: HRRs are delineated from observed patient travel for
# major cardiovascular surgery and neurosurgery, so excluding adjacent HRRs
# would drop valid variation without strengthening the exclusion restriction.
md_event_roster_hrr <- merge(md_event_roster,
                             unique(md_hosp_geo[, .(HOSPITAL_ID, HRR_NUM)]),
                             by = "HOSPITAL_ID", all.x = TRUE)

hrr_comp <- build_dynamic_competitor_instrument(
  md_county_panel, md_event_roster_hrr, "HRR_NUM", "HRR_COMP",
  exclude_own_system = TRUE)

hrr_strict <- build_dynamic_competitor_instrument(
  md_county_panel, md_event_roster_hrr, "HRR_NUM", "HRR_STRICT",
  exclude_own_system = FALSE)

hrr_six <- merge(hrr_comp, hrr_strict, by = c("HOSPITAL_ID", "POST_MONTH"))

HRR_INSTRUMENTS <- c(
  HRR_competitor_hospitals = "HRR_COMP_HOSPITALS",
  HRR_competitor_systems   = "HRR_COMP_SYSTEMS",
  HRR_competitor_hrrs      = "HRR_COMP_GEOS",
  HRR_strict_hospitals     = "HRR_STRICT_HOSPITALS",
  HRR_strict_systems       = "HRR_STRICT_SYSTEMS",
  HRR_strict_hrrs          = "HRR_STRICT_GEOS")


# -----------------------------------------------------------------------------
# 10.5  Descriptives and the correlation matrix
# -----------------------------------------------------------------------------
# iv_descriptives() returns the mean, median, SD, share of zeros, 90th
# percentile, and maximum of each instrument and prints the roster size for
# scale (MD24_hrr_six_descriptives.csv). The correlation matrix is printed,
# not saved. Instruments that correlate above roughly 0.95 measure the same
# variation and are not independent evidence.
iv_descriptives <- function(dt, instruments, roster_n) {
  out <- rbindlist(lapply(names(instruments), function(nm) {
    x <- dt[[instruments[[nm]]]]
    data.table(INSTRUMENT = nm, MEAN = round(mean(x), 3), MEDIAN = median(x),
               SD = round(sd(x), 3), SHARE_ZERO = round(mean(x == 0), 3),
               P90 = quantile(x, 0.90), MAX = max(x))
  }))
  cat(sprintf("\nRoster size for scale reference: %d hospitals\n", roster_n))
  out
}

hrr_desc <- iv_descriptives(hrr_six, HRR_INSTRUMENTS, nrow(md_event_roster))
save_csv(hrr_desc, "MD24_hrr_six_descriptives.csv")
cat("\n=== SIX HRR INSTRUMENTS ===\n"); print(hrr_desc)

iv_mat <- as.matrix(hrr_six[, unname(HRR_INSTRUMENTS), with = FALSE])
colnames(iv_mat) <- names(HRR_INSTRUMENTS)
cat("\n=== CORRELATION AMONG THE SIX ===\n")
print(round(cor(iv_mat, use = "complete.obs"), 3))


# -----------------------------------------------------------------------------
# 10.6  Pooled second stage with an explicit first stage
# -----------------------------------------------------------------------------
# For each instrument, run_pooled_with_first_stage() fits the pooled first
# stage of the treatment on the instrument (baseline controls, fixed effects,
# and clusters) and reports its coefficient, t-statistic, within-R2, and
# effect per SD of the instrument; it then fits the interacted shoppability
# model with estimate_interacted(). Instruments absent from the panel or
# without usable variation are skipped with a message. shoppability_gap()
# returns the shoppable minus non-shoppable IV_PERCENT for each instrument,
# in percentage points.
#
# The HRR panel is loaded again as md_hrr_pooled, and the six HRR instruments
# are merged on by hospital-month. PRIMARY_INSTRUMENT (IV06, built at county
# level) is added as a benchmark under the name County_gated_IV06_benchmark.
# Writes MD25_hrr_six_instruments_results.csv and
# MD26_hrr_shoppability_gap_by_instrument.csv. md_hrr_pooled stays in memory
# until 15.0 removes it.
run_pooled_with_first_stage <- function(panel, instruments, geo_tag) {
  rbindlist(lapply(names(instruments), function(nm) {
    iv <- instruments[[nm]]
    cat("\n---", nm, "---\n")
    
    if (!(iv %in% names(panel)) || !has_usable_variation(panel[[iv]])) {
      cat("  absent or no usable variation, skipped\n"); return(NULL)
    }
    
    fs <- tryCatch(feols(
      build_ols_formula(ENDOGENOUS_VARIABLE,
                        c(iv, available_columns(panel, BASELINE_CONTROLS)),
                        available_columns(panel, BASELINE_FIXED_EFFECTS)),
      data = panel, cluster = build_cluster_formula(BASELINE_CLUSTERS),
      warn = FALSE, notes = FALSE), error = function(e) NULL)
    if (is.null(fs)) { cat("  first stage failed\n"); return(NULL) }
    
    fs_b  <- unname(coef(fs)[iv])
    fs_se <- unname(sqrt(vcov(fs)[iv, iv]))
    fs_r2 <- tryCatch(fitstat(fs, "wr2")$wr2, error = function(e) NA_real_)
    iv_sd <- sd(panel[[iv]], na.rm = TRUE)
    cat(sprintf("  first stage: coef %.5f, t %.2f, within-R2 %.4f, per-SD %.3f\n",
                fs_b, fs_b / fs_se, fs_r2, fs_b * iv_sd))
    
    res <- estimate_interacted(
      data = panel, moderator = MD_SCHEME, moderator_type = "categorical",
      outcome = PRIMARY_OUTCOME, instrument = iv,
      label = paste0(geo_tag, "/", nm), instrument_label = nm,
      moderator_label = MD_SCHEME)
    if (is.null(res)) { cat("  second stage returned NULL\n"); return(NULL) }
    
    r <- copy(res$rows)
    r[, `:=`(IV_NAME = nm, IV_COLUMN = iv, FS_COEF = fs_b, FS_SE = fs_se,
             FS_T = fs_b / fs_se, FS_WITHIN_R2 = fs_r2,
             FS_EFFECT_PER_SD = fs_b * iv_sd, IV_SD = iv_sd)]
    r
  }), fill = TRUE)
}

shoppability_gap <- function(results) {
  g <- dcast(results, IV_NAME ~ TERM, value.var = "IV_PERCENT")
  if (!all(c("Shoppable", "Non_shoppable") %in% names(g))) return(NULL)
  g[, GAP_PP := round(Shoppable - Non_shoppable, 3)]
  g[order(GAP_PP)]
}

md_hrr_pooled <- md_load_panel("HRR")
md_hrr_pooled <- merge(md_hrr_pooled, hrr_six,
                       by = c("HOSPITAL_ID", "POST_MONTH"), all.x = TRUE)

hrr_results <- run_pooled_with_first_stage(
  md_hrr_pooled,
  c(HRR_INSTRUMENTS, County_gated_IV06_benchmark = PRIMARY_INSTRUMENT),
  "HRR")

save_csv(hrr_results, "MD25_hrr_six_instruments_results.csv")

cat("\n=== FIRST STAGE, HRR ===\n")
print(unique(hrr_results[, .(IV_NAME, FS_COEF = round(FS_COEF, 5),
                             FS_T = round(FS_T, 2),
                             FS_WITHIN_R2 = round(FS_WITHIN_R2, 4),
                             FS_PER_SD = round(FS_EFFECT_PER_SD, 3),
                             IV_SD = round(IV_SD, 2))], by = "IV_NAME"))

cat("\n=== POOLED SECOND STAGE, HRR ===\n")
print(hrr_results[, .(IV_NAME, TERM, IV_PERCENT = round(IV_PERCENT, 3),
                      IV_P = round(IV_P, 4),
                      WALD = round(FIRST_STAGE_WALD_THIS_EQ, 1), N_OBSERVATIONS)])

hrr_gap <- shoppability_gap(hrr_results)
if (!is.null(hrr_gap)) {
  cat("\n=== SHOPPABILITY GAP, HRR (pp) ===\n"); print(hrr_gap)
  save_csv(hrr_gap, "MD26_hrr_shoppability_gap_by_instrument.csv")
}


# =============================================================================
# Part 11: HSA instruments
# =============================================================================
#
# HSA adjacency from the Dartmouth boundary file, a check of the builder's
# adjacency branch against IV15, nine HSA instruments, and the pooled model
# at HSA.


# -----------------------------------------------------------------------------
# 11.1  HSA adjacency from Dartmouth boundary polygons
# -----------------------------------------------------------------------------
# Adjacency is computed from the HSA polygons (st_touches()) rather than
# inferred from county adjacency and the hospital-to-HSA mapping, which would
# miss HSAs that connect only through a county with no sample hospital. The
# block stops if the shapefile (HSA_SHP) is missing, if a geometry is still
# invalid after st_make_valid(), or if the adjacency is not symmetric. Writes
# the pairs to MD30_hsa_adjacency_pairs.csv.
#
# sf_use_s2(FALSE) switches sf to planar geometry for the rest of the session.
# The spherical engine (s2) rejects these 1993 boundaries with an
# edge-crossing error even after st_make_valid(); the shapefile is not
# corrupt. Adjacency is topological at this spatial scale, so planar geometry
# does not affect the result.
sf_use_s2(FALSE)

HSA_SHP <- file.path(PANEL_DIR, "HSA_Bdry__AK_HI_unmodified", "hsa-shapefile",
                     "HsaBdry_AK_HI_unmodified.shp")
stopifnot("HSA shapefile not found" = file.exists(HSA_SHP))

hsa_sf <- st_read(HSA_SHP, quiet = TRUE)
n_invalid <- sum(!st_is_valid(hsa_sf))
if (n_invalid > 0) {
  cat(sprintf("Repairing %d invalid geometries.\n", n_invalid))
  hsa_sf <- st_make_valid(hsa_sf)
  stopifnot("geometry still invalid after repair" = all(st_is_valid(hsa_sf)))
}

nb <- st_touches(hsa_sf, sparse = TRUE)
hsa_adjacency <- rbindlist(lapply(seq_along(nb), function(i) {
  if (!length(nb[[i]])) return(NULL)
  data.table(GEO_A = as.character(hsa_sf$HSA93[i]),
             GEO_B = as.character(hsa_sf$HSA93[nb[[i]]]))
}))
setkey(hsa_adjacency, GEO_A, GEO_B)

adj_deg <- hsa_adjacency[, .N, by = GEO_A]
cat(sprintf("\nHSA adjacency: %s pairs, %d HSAs with at least one neighbour\n",
            format(nrow(hsa_adjacency), big.mark = ","), nrow(adj_deg)))
cat(sprintf("Neighbours per HSA: mean %.1f, median %d, max %d. Islands: %d\n",
            mean(adj_deg$N), median(adj_deg$N), max(adj_deg$N),
            nrow(hsa_sf) - nrow(adj_deg)))

rev_pairs <- hsa_adjacency[, .(GEO_A = GEO_B, GEO_B = GEO_A)]
stopifnot("adjacency is not symmetric" = nrow(fsetdiff(hsa_adjacency, rev_pairs)) == 0)

# HSA-to-HRR nesting, on which the HSA-against-HRR comparison in Part 12
# rests: prints the number of HSAs that map to more than one HRR in the
# boundary file (no check).
nest <- as.data.table(st_drop_geometry(hsa_sf))[, .(N_HRR = uniqueN(HRR93)), by = HSA93]
cat(sprintf("HSAs mapping to more than one HRR: %d of %d. Distinct HRRs: %d\n",
            nest[N_HRR > 1, .N], nrow(nest),
            uniqueN(st_drop_geometry(hsa_sf)$HRR93)))

# Empty string is a missing value that would otherwise become its own market.
md_hosp_geo[HSA_NUM == "", HSA_NUM := NA_character_]

save_csv(hsa_adjacency, "MD30_hsa_adjacency_pairs.csv")


# -----------------------------------------------------------------------------
# 11.2  Validation of the adjacency branch
# -----------------------------------------------------------------------------
# IV15 (Z_SYS_COMPETITOR_EXCL_ADJACENT_9M_EXCL_CURRENT) is the county
# competitor-hospitals count with peers in adjacent counties excluded. The
# builder is run on five-digit county FIPS codes with the county pairs in
# HPT_COUNTY_ADJACENCY_PAIRS.csv as `adjacency`, and the run stops unless it
# matches IV15 exactly in at least 98% of hospital-months.
#
# IV15 is NA wherever the focal hospital has no county FIPS code, almost
# entirely in Connecticut after the 2022 planning-region reorganization;
# those rows are left out of the comparison. formatC() turns a missing
# COUNTY_FIPS into "   NA", so roster hospitals without a code form one
# county in this rebuild (17.1 avoids this with md_fips5()), and it pads a
# character code shorter than five characters with spaces, not zeros.
cty_adj <- fread(file.path(PANEL_DIR, "HPT_COUNTY_ADJACENCY_PAIRS.csv"),
                 colClasses = list(character = c("COUNTY_GEOID", "NEIGHBOR_GEOID")))
cty_adj_pairs <- unique(cty_adj[, .(
  GEO_A = formatC(COUNTY_GEOID,   width = 5, flag = "0"),
  GEO_B = formatC(NEIGHBOR_GEOID, width = 5, flag = "0"))])
setkey(cty_adj_pairs, GEO_A, GEO_B)

md_event_roster_fips <- copy(md_event_roster)
md_event_roster_fips[, COUNTY_FIPS := formatC(COUNTY_FIPS, width = 5, flag = "0")]

cty_excl_adj <- build_dynamic_competitor_instrument(
  md_county_panel, md_event_roster_fips, "COUNTY_FIPS", "CTY_EXCLADJ",
  exclude_own_system = TRUE, adjacency = cty_adj_pairs)

va <- merge(cty_excl_adj,
            unique(md_county_panel[, .(HOSPITAL_ID, POST_MONTH,
                                       T_IV15 = Z_SYS_COMPETITOR_EXCL_ADJACENT_9M_EXCL_CURRENT)]),
            by = c("HOSPITAL_ID", "POST_MONTH"))[!is.na(T_IV15)]

adj_match <- mean(va$CTY_EXCLADJ_HOSPITALS == va$T_IV15)
cat(sprintf("\nAdjacency validation against IV15: %s rows, exact match %.4f, mean abs diff %.4f\n",
            format(nrow(va), big.mark = ","), adj_match,
            mean(abs(va$CTY_EXCLADJ_HOSPITALS - va$T_IV15))))
stopifnot("adjacency branch does not reproduce IV15" = adj_match >= 0.98)


# -----------------------------------------------------------------------------
# 11.3  Nine HSA instruments
# -----------------------------------------------------------------------------
# The six instruments of the HRR design (competitor and strict families,
# counting hospitals, systems, or HSAs) plus the competitor family with peers
# in adjacent HSAs excluded. The adjacency variant is built here, unlike at
# HRR, because HSAs are small enough that spillover from adjacent HSAs is a
# concern. Writes MD27_hsa_nine_descriptives.csv and prints the correlation
# matrix.
md_event_roster_hsa <- merge(md_event_roster,
                             unique(md_hosp_geo[!is.na(HSA_NUM), .(HOSPITAL_ID, HSA_NUM)]),
                             by = "HOSPITAL_ID", all.x = TRUE)

hsa_comp    <- build_dynamic_competitor_instrument(
  md_county_panel, md_event_roster_hsa, "HSA_NUM", "HSA_COMP",
  exclude_own_system = TRUE)

hsa_strict  <- build_dynamic_competitor_instrument(
  md_county_panel, md_event_roster_hsa, "HSA_NUM", "HSA_STRICT",
  exclude_own_system = FALSE)

hsa_excladj <- build_dynamic_competitor_instrument(
  md_county_panel, md_event_roster_hsa, "HSA_NUM", "HSA_EXCLADJ",
  exclude_own_system = TRUE, adjacency = hsa_adjacency)

hsa_nine <- Reduce(function(a, b) merge(a, b, by = c("HOSPITAL_ID", "POST_MONTH")),
                   list(hsa_comp, hsa_strict, hsa_excladj))

HSA_INSTRUMENTS <- c(
  HSA_competitor_hospitals    = "HSA_COMP_HOSPITALS",
  HSA_competitor_systems      = "HSA_COMP_SYSTEMS",
  HSA_competitor_hsas         = "HSA_COMP_GEOS",
  HSA_strict_hospitals        = "HSA_STRICT_HOSPITALS",
  HSA_strict_systems          = "HSA_STRICT_SYSTEMS",
  HSA_strict_hsas             = "HSA_STRICT_GEOS",
  HSA_excl_adjacent_hospitals = "HSA_EXCLADJ_HOSPITALS",
  HSA_excl_adjacent_systems   = "HSA_EXCLADJ_SYSTEMS",
  HSA_excl_adjacent_hsas      = "HSA_EXCLADJ_GEOS")

hsa_desc <- iv_descriptives(hsa_nine, HSA_INSTRUMENTS, nrow(md_event_roster))
save_csv(hsa_desc, "MD27_hsa_nine_descriptives.csv")
cat("\n=== NINE HSA INSTRUMENTS ===\n"); print(hsa_desc)

iv_mat <- as.matrix(hsa_nine[, unname(HSA_INSTRUMENTS), with = FALSE])
colnames(iv_mat) <- names(HSA_INSTRUMENTS)
cat("\n=== CORRELATION AMONG THE NINE ===\n")
print(round(cor(iv_mat, use = "complete.obs"), 3))


# -----------------------------------------------------------------------------
# 11.4  Where the HSA instrument has variation
# -----------------------------------------------------------------------------
# Dartmouth delineated HSAs around individual hospital catchments, so most
# contain a single hospital. In a single-hospital HSA the only system that can
# be present is the focal hospital's own, which the competitor filter
# excludes, so the competitor instruments are zero there. The block prints the
# number of HSAs by hospital count (from md_hosp_geo), the share of
# county-panel rows with a positive HSA_COMP_HOSPITALS by HSA size, and the
# HSA size and metro share of rows with a positive and a zero value. This
# sets the sample restriction in Part 12 and the scope condition for every
# HSA estimate.
hsa_size <- md_hosp_geo[!is.na(HSA_NUM), .(N_HOSP_IN_HSA = .N), by = HSA_NUM]
hsa_size[, BUCKET := fifelse(N_HOSP_IN_HSA == 1, "1",
                             fifelse(N_HOSP_IN_HSA == 2, "2",
                                     fifelse(N_HOSP_IN_HSA == 3, "3", "4+")))]

md_hsa_check <- merge(md_county_panel[, .(HOSPITAL_ID, POST_MONTH)],
                      hsa_nine, by = c("HOSPITAL_ID", "POST_MONTH"))
md_hsa_check <- merge(md_hsa_check,
                      unique(md_hosp_geo[, .(HOSPITAL_ID, HSA_NUM, RUCC_2023)]),
                      by = "HOSPITAL_ID", all.x = TRUE)
md_hsa_check <- merge(md_hsa_check, hsa_size, by = "HSA_NUM", all.x = TRUE)
md_hsa_check[, NONZERO := HSA_COMP_HOSPITALS > 0]

cat("\n=== SAMPLE HSAs BY HOSPITAL COUNT ===\n")
print(hsa_size[, .N, by = BUCKET][order(BUCKET)])

cat("\n=== SHARE WITH NONZERO INSTRUMENT, BY HSA SIZE ===\n")
print(md_hsa_check[, .(N = .N, SHARE_NONZERO = round(mean(NONZERO), 3)),
                   by = BUCKET][order(BUCKET)])

cat("\n=== COMPOSITION OF THE NONZERO SUBSAMPLE ===\n")
print(md_hsa_check[, .(N = .N,
                       MEAN_HOSP_IN_HSA = round(mean(N_HOSP_IN_HSA, na.rm = TRUE), 2),
                       SHARE_METRO = round(mean(RUCC_2023 %in% 1:3, na.rm = TRUE), 3)),
                   by = NONZERO])


# -----------------------------------------------------------------------------
# 11.5  Pooled second stage at HSA
# -----------------------------------------------------------------------------
# The HSA panel is loaded as md_hsa_pooled, the nine HSA instruments are
# merged on, and run_pooled_with_first_stage() (10.6) is run with IV06 as a
# benchmark (MD28_hsa_instruments_results.csv, MD29_hsa_shoppability_gap.csv).
# The HSA_COMP_HOSPITALS model is then re-estimated on the rows where that
# instrument is positive, so that the estimate and the population it
# describes refer to the same units; that result is printed, not saved.
# md_hsa_pooled stays in memory until 15.0 removes it.
md_hsa_pooled <- md_load_panel("HSA")
md_hsa_pooled <- merge(md_hsa_pooled, hsa_nine,
                       by = c("HOSPITAL_ID", "POST_MONTH"), all.x = TRUE)

hsa_results <- run_pooled_with_first_stage(
  md_hsa_pooled,
  c(HSA_INSTRUMENTS, County_gated_IV06_benchmark = PRIMARY_INSTRUMENT),
  "HSA")

save_csv(hsa_results, "MD28_hsa_instruments_results.csv")

cat("\n=== FIRST STAGE, HSA ===\n")
print(unique(hsa_results[, .(IV_NAME, FS_COEF = round(FS_COEF, 5),
                             FS_T = round(FS_T, 2),
                             FS_WITHIN_R2 = round(FS_WITHIN_R2, 4),
                             FS_PER_SD = round(FS_EFFECT_PER_SD, 3),
                             IV_SD = round(IV_SD, 2))], by = "IV_NAME"))

cat("\n=== POOLED SECOND STAGE, HSA ===\n")
print(hsa_results[, .(IV_NAME, TERM, IV_PERCENT = round(IV_PERCENT, 3),
                      IV_P = round(IV_P, 4),
                      WALD = round(FIRST_STAGE_WALD_THIS_EQ, 1), N_OBSERVATIONS)])

hsa_gap <- shoppability_gap(hsa_results)
if (!is.null(hsa_gap)) {
  cat("\n=== SHOPPABILITY GAP, HSA (pp) ===\n"); print(hsa_gap)
  save_csv(hsa_gap, "MD29_hsa_shoppability_gap.csv")
}

# Nonzero subsample: the rows whose own HSA_COMP_HOSPITALS is positive. Part
# 12 restricts by market instead, keeping every row of an HSA in which the
# instrument is positive for at least one row.
md_hsa_nz <- md_hsa_pooled[HSA_COMP_HOSPITALS > 0]
cat(sprintf("\nNonzero subsample: %s rows, %d hospitals, %d HSAs\n",
            format(nrow(md_hsa_nz), big.mark = ","),
            uniqueN(md_hsa_nz$HOSPITAL_ID), uniqueN(md_hsa_nz$ANALYSIS_MARKET)))

res_nz <- estimate_interacted(
  data = md_hsa_nz, moderator = MD_SCHEME, moderator_type = "categorical",
  outcome = PRIMARY_OUTCOME, instrument = "HSA_COMP_HOSPITALS",
  label = "HSA/competitor_hospitals, nonzero only",
  instrument_label = "HSA_competitor_hospitals_nonzero", moderator_label = MD_SCHEME)

if (!is.null(res_nz))
  print(res_nz$rows[, .(TERM, IV_PERCENT = round(IV_PERCENT, 3),
                        IV_P = round(IV_P, 4),
                        WALD = round(FIRST_STAGE_WALD_THIS_EQ, 1), N_OBSERVATIONS)])


# =============================================================================
# Part 12: Concept-level analysis at HSA and HRR
# =============================================================================
#
# The county concept-level machinery applied at HSA and HRR. For each concept,
# estimate_concept_level() gives reduced-form, first-stage, and IV estimates,
# and run_meta_regressions() covers every SCHEME_ID in md_schemes_long (all 18
# classification schemes) in one call. P12_RUNS defines three runs, each with
# one instrument:
#   hsa_comp      HSA, HSA_COMP_HOSPITALS, restricted sample      (12.A)
#   hsa_excladj   HSA, HSA_EXCLADJ_HOSPITALS, restricted sample   (12.B)
#   hrr_comp      HRR, HRR_COMP_HOSPITALS, full sample            (12.C)
#
# Instrument choice. HRR_strict_systems has the strongest HRR first stage,
# with a within-R2 of 0.425, but the strict family keeps own-system peers.
# Only the competitor family excludes own-system rollout, which is the basis
# of the exclusion-restriction argument in the county specification, so the
# HRR and HSA estimates use competitor variants and rest on the same argument
# as the county headline.
#
# HSA sample restriction. The HSA runs keep only HSAs in which the instrument
# is positive for at least one row. Single-hospital HSAs cannot generate
# instrument variation (11.4), so the restriction follows from the geography.
# It also makes concept-level estimation feasible, cutting the singleton
# fixed-effect share from 0.734 to 0.257. HSA estimates therefore describe
# metropolitan multi-hospital service areas.
#
# Fixed effects and clustering are the baseline ones at every geography (see
# 0.2).
#
# Runtime: roughly two hours per geography per instrument. Each run has its
# own cache keys (p12_concept_<run>_1inst, p12_meta_rf_<run>_1inst, and
# p12_meta_iv_<run>_1inst), saves MD31_concept_<run>_PARTIAL.csv every 50
# concepts, and writes MD31_concept_level_<run>.csv, MD32_meta_rf_<run>.csv,
# and MD33_meta_iv_<run>.csv.

P12_RUNS <- list(
  list(key = "hsa_comp",    geo = "HSA", iv_col = "HSA_COMP_HOSPITALS",
       iv_label = "HSA_competitor_hospitals",    restrict_nonzero = TRUE),
  list(key = "hsa_excladj", geo = "HSA", iv_col = "HSA_EXCLADJ_HOSPITALS",
       iv_label = "HSA_excl_adjacent_hospitals", restrict_nonzero = TRUE),
  list(key = "hrr_comp",    geo = "HRR", iv_col = "HRR_COMP_HOSPITALS",
       iv_label = "HRR_competitor_hospitals",    restrict_nonzero = FALSE))

p12_key <- function(stem, run) sprintf("%s_%s_1inst", stem, run$key)

# p12_build_panel() loads the panel for run$geo, merges the run's instrument
# from hsa_nine or hrr_six by hospital-month (stops if the row count
# changes), and reports the rows left without an instrument value. With
# restrict_nonzero = TRUE it keeps only the markets in which the instrument is
# positive for at least one row. It prints the resulting fixed-effect cell
# structure.
#
# The instruments are built with focal_panel = md_county_panel, so their
# hospital-months come from the county panel (3,723 hospitals) rather than
# the HSA or HRR panel (3,722). Rows of the HSA or HRR panel without a match
# have NA for the instrument and drop out of the models; the unmatched count
# shows how many there are.
p12_build_panel <- function(run) {
  cat("\n", strrep("=", 70), "\n", sep = "")
  cat("PANEL:", run$geo, "| instrument:", run$iv_label, "\n")
  cat(strrep("=", 70), "\n", sep = "")
  
  p <- md_load_panel(run$geo)
  iv_source <- if (run$geo == "HSA") hsa_nine else hrr_six
  stopifnot("instrument column missing" = run$iv_col %in% names(iv_source))
  
  keep <- c("HOSPITAL_ID", "POST_MONTH", run$iv_col)
  n_before <- nrow(p)
  p <- merge(p, iv_source[, ..keep], by = c("HOSPITAL_ID", "POST_MONTH"), all.x = TRUE)
  stopifnot("instrument merge changed row count" = nrow(p) == n_before)
  
  n_unmatched <- sum(is.na(p[[run$iv_col]]))
  cat(sprintf("  instrument merged, unmatched rows %s (%.3f%%)\n",
              format(n_unmatched, big.mark = ","), 100 * n_unmatched / nrow(p)))
  
  if (isTRUE(run$restrict_nonzero)) {
    nz_markets <- unique(p[get(run$iv_col) > 0, ANALYSIS_MARKET])
    n_pre <- nrow(p)
    p <- p[ANALYSIS_MARKET %in% nz_markets]
    cat(sprintf("  restricted to markets with instrument variation: %s to %s rows, %d markets\n",
                format(n_pre, big.mark = ","), format(nrow(p), big.mark = ","),
                uniqueN(p$ANALYSIS_MARKET)))
  }
  
  fe_sizes <- p[, .N, by = MARKET_ID]$N
  cat(sprintf("  FE cells %s, singleton share %.3f, mean rows per cell %.2f\n",
              format(length(fe_sizes), big.mark = ","),
              mean(fe_sizes == 1), mean(fe_sizes)))
  cat(sprintf("  rows %s, hospitals %s, markets %s, concepts %s\n",
              format(nrow(p), big.mark = ","),
              format(uniqueN(p$HOSPITAL_ID), big.mark = ","),
              format(uniqueN(p$ANALYSIS_MARKET), big.mark = ","),
              format(uniqueN(p$FINAL_CONCEPT_ID), big.mark = ",")))
  p
}

p12_run_concept <- function(run) {
  t0 <- Sys.time()
  cat("\nStarted:", format(t0, "%Y-%m-%d %H:%M:%S"), "\n")
  
  p <- p12_build_panel(run)
  instruments <- setNames(run$iv_col, run$iv_label)
  
  cr <- cache_or_run(
    p12_key("p12_concept", run),
    estimate_concept_level(
      build_concept_panel(p, instruments = instruments),
      instruments = instruments,
      save_stem = paste0("MD31_concept_", run$key)))
  
  cr[, `:=`(GEOGRAPHY = run$geo, IV_LABEL = run$iv_label, RUN_KEY = run$key)]
  
  cat(sprintf("\nDone %s: %s rows, %s concepts, %.1f min\n",
              run$key, format(nrow(cr), big.mark = ","),
              format(uniqueN(cr$FINAL_CONCEPT_ID), big.mark = ","),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  
  save_csv(cr, paste0("MD31_concept_level_", run$key, ".csv"))
  rm(p); invisible(gc())
  cr
}

p12_run_meta <- function(cr, run) {
  mi <- prepare_meta_input(cr, md_schemes_long)
  
  m_rf <- cache_or_run(p12_key("p12_meta_rf", run),
                       run_meta_regressions(mi, md_schemes_long, dep = "RF_COEF", se = "RF_SE"))
  m_iv <- cache_or_run(p12_key("p12_meta_iv", run),
                       run_meta_regressions(mi, md_schemes_long, dep = "IV_COEF", se = "IV_SE"))
  
  m_rf[, `:=`(GEOGRAPHY = run$geo, IV_LABEL = run$iv_label, RUN_KEY = run$key)]
  m_iv[, `:=`(GEOGRAPHY = run$geo, IV_LABEL = run$iv_label, RUN_KEY = run$key)]
  
  save_csv(m_rf, paste0("MD32_meta_rf_", run$key, ".csv"))
  save_csv(m_iv, paste0("MD33_meta_iv_", run$key, ".csv"))
  list(rf = m_rf, iv = m_iv)
}

# p12_run_pooled_all_schemes() loads the run's panel with p12_build_panel()
# and fits the pooled interacted model (estimate_interacted(), categorical)
# once for each of the six primary scheme columns in SCHEME_COLUMNS. Writes
# MD34_pooled_all_schemes_<run>.csv. Called in 12.F.
p12_run_pooled_all_schemes <- function(run) {
  p <- p12_build_panel(run)
  
  out <- rbindlist(lapply(names(SCHEME_COLUMNS), function(scheme_label) {
    sc <- SCHEME_COLUMNS[[scheme_label]]
    if (!(sc %in% names(p))) return(NULL)
    if (uniqueN(p[[sc]][!is.na(p[[sc]])]) < 2L) return(NULL)
    
    res <- estimate_interacted(
      data = p, moderator = sc, moderator_type = "categorical",
      outcome = PRIMARY_OUTCOME, instrument = run$iv_col,
      label = paste0(run$key, "/", sc),
      instrument_label = run$iv_label, moderator_label = scheme_label)
    if (is.null(res)) return(NULL)
    
    r <- copy(res$rows)
    r[, `:=`(SCHEME_LABEL = scheme_label, SCHEME_COLUMN = sc,
             GEOGRAPHY = run$geo, RUN_KEY = run$key)]
    r
  }), fill = TRUE)
  
  if (nrow(out) > 0L)
    save_csv(out, paste0("MD34_pooled_all_schemes_", run$key, ".csv"))
  
  rm(p); invisible(gc())
  out
}


# -----------------------------------------------------------------------------
# 12.A  HSA, competitor instrument
# -----------------------------------------------------------------------------
run_hsa_comp  <- P12_RUNS[[1]]
cr_hsa_comp   <- p12_run_concept(run_hsa_comp)
meta_hsa_comp <- p12_run_meta(cr_hsa_comp, run_hsa_comp)
cat("\n=== HSA competitor, meta-regression (IV) ===\n"); print(meta_hsa_comp$iv)


# -----------------------------------------------------------------------------
# 12.B  HSA, adjacency-excluded instrument
# -----------------------------------------------------------------------------
run_hsa_adj  <- P12_RUNS[[2]]
cr_hsa_adj   <- p12_run_concept(run_hsa_adj)
meta_hsa_adj <- p12_run_meta(cr_hsa_adj, run_hsa_adj)
cat("\n=== HSA adjacency-excluded, meta-regression (IV) ===\n"); print(meta_hsa_adj$iv)


# -----------------------------------------------------------------------------
# 12.C  HRR, competitor instrument
# -----------------------------------------------------------------------------
run_hrr_comp  <- P12_RUNS[[3]]
cr_hrr_comp   <- p12_run_concept(run_hrr_comp)
meta_hrr_comp <- p12_run_meta(cr_hrr_comp, run_hrr_comp)
cat("\n=== HRR competitor, meta-regression (IV) ===\n"); print(meta_hrr_comp$iv)


# -----------------------------------------------------------------------------
# 12.D  Concept-level summaries across the three runs
# -----------------------------------------------------------------------------
# Stacks the three runs (MD35_concept_level_all_runs.csv) and summarizes each:
# number of concepts, median IV estimate in percent, share negative, share
# significant at 5%, median first-stage F, and share with a first-stage F of
# at least 10 (MD36_concept_summary_all_runs.csv).
p12_all_concepts <- rbindlist(list(cr_hsa_comp, cr_hsa_adj, cr_hrr_comp), fill = TRUE)
save_csv(p12_all_concepts, "MD35_concept_level_all_runs.csv")

p12_summary <- p12_all_concepts[, .(
  N_CONCEPTS     = .N,
  MEDIAN_IV_PCT  = round(median(IV_ESTIMATE_PERCENT, na.rm = TRUE), 3),
  SHARE_NEGATIVE = round(mean(IV_COEF < 0, na.rm = TRUE), 3),
  SHARE_SIG_5PCT = round(mean(IV_P < 0.05, na.rm = TRUE), 3),
  MEDIAN_FS_F    = round(median(FS_F, na.rm = TRUE), 1),
  SHARE_FS_F_10  = round(mean(FS_F >= 10, na.rm = TRUE), 3)
), by = .(RUN_KEY, GEOGRAPHY, IV_LABEL)]

save_csv(p12_summary, "MD36_concept_summary_all_runs.csv")
cat("\n=== CONCEPT-LEVEL SUMMARY ===\n"); print(p12_summary)


# -----------------------------------------------------------------------------
# 12.E  HSA against HRR, concept by concept
# -----------------------------------------------------------------------------
# p12_crossover() takes, for each concept, the HSA IV coefficient minus the
# HRR coefficient from cr_hrr_comp, once for each HSA run. A negative delta
# means the HSA market disciplines price more for that concept, the pattern
# predicted for routine services delivered within a local catchment; a
# positive delta is the referral-catchment pattern. Concepts with a
# first-stage F of at least 10 on both sides are summarized by category under
# each of the 18 schemes, by the slope of the delta on the ordinal
# shoppability rank, and by clinical family. Writes MD37 (every matched
# concept), MD38, MD39, and MD40.
#
# The comparison is descriptive. The two sides come from non-nested models on
# different samples (the HSA runs are restricted, the HRR run is not), so the
# per-concept difference has no p-value. The slope on ordinal shoppability
# has one, but it tests whether the delta trends with shoppability, which is
# a weaker claim than equality across models. Parts 15 and 15E repeat the
# cross-geography comparison with six instruments at each geography, and
# 15E.11 checks the robustness of the slope from this block.
p12_crossover <- function(cr_hsa, hsa_label) {
  d <- merge(
    cr_hsa[, .(FINAL_CONCEPT_ID, FINAL_FAMILY_ID, SERVICE_LABEL,
               IV_HSA = IV_COEF, IV_SE_HSA = IV_SE, IV_P_HSA = IV_P,
               FS_F_HSA = FS_F, N_OBS_HSA = N_OBSERVATIONS)],
    cr_hrr_comp[, .(FINAL_CONCEPT_ID,
                    IV_HRR = IV_COEF, IV_SE_HRR = IV_SE, IV_P_HRR = IV_P,
                    FS_F_HRR = FS_F, N_OBS_HRR = N_OBSERVATIONS)],
    by = "FINAL_CONCEPT_ID")
  
  d[, IV_DELTA_HSA_MINUS_HRR := IV_HSA - IV_HRR]
  d[, HSA_INSTRUMENT := hsa_label]
  
  ok <- d[is.finite(IV_DELTA_HSA_MINUS_HRR) & FS_F_HSA >= 10 & FS_F_HRR >= 10]
  cat(sprintf("\n[%s] concepts matched %d, first-stage F at least 10 on both sides %d\n",
              hsa_label, nrow(d), nrow(ok)))
  
  ds <- copy(ok)
  for (sid in unique(md_schemes_long$SCHEME_ID)) {
    a <- unique(md_schemes_long[SCHEME_ID == sid,
                                .(FINAL_CONCEPT_ID = ANALYSIS_CONCEPT_ID,
                                  CAT = SHOPPABILITY_CATEGORY,
                                  ORD = SHOPPABILITY_ORDINAL)])
    setnames(a, c("CAT", "ORD"), paste0(c("CAT_", "ORD_"), sid))
    ds <- merge(ds, a, by = "FINAL_CONCEPT_ID", all.x = TRUE)
  }
  
  by_cat <- rbindlist(lapply(unique(md_schemes_long$SCHEME_ID), function(sid) {
    ccol <- paste0("CAT_", sid)
    if (!(ccol %in% names(ds))) return(NULL)
    dd <- ds[!is.na(get(ccol))]
    if (uniqueN(dd[[ccol]]) < 2L || nrow(dd) < MIN_CONCEPTS_META) return(NULL)
    dd[, .(SCHEME_ID = sid, HSA_INSTRUMENT = hsa_label, N_CONCEPTS = .N,
           MEAN_DELTA   = round(mean(IV_DELTA_HSA_MINUS_HRR, na.rm = TRUE), 4),
           MEDIAN_DELTA = round(median(IV_DELTA_HSA_MINUS_HRR, na.rm = TRUE), 4),
           SHARE_HSA_STRONGER = round(mean(IV_DELTA_HSA_MINUS_HRR < 0, na.rm = TRUE), 3)),
       by = c(CATEGORY = ccol)]
  }), fill = TRUE)
  
  slopes <- rbindlist(lapply(unique(md_schemes_long$SCHEME_ID), function(sid) {
    ocol <- paste0("ORD_", sid)
    if (!(ocol %in% names(ds))) return(NULL)
    dd <- ds[is.finite(get(ocol)) & is.finite(IV_DELTA_HSA_MINUS_HRR)]
    if (nrow(dd) < MIN_CONCEPTS_META || !has_usable_variation(dd[[ocol]])) return(NULL)
    fit <- tryCatch(lm(dd$IV_DELTA_HSA_MINUS_HRR ~ dd[[ocol]]), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    s <- summary(fit)$coefficients
    if (nrow(s) < 2L) return(NULL)
    data.table(SCHEME_ID = sid, HSA_INSTRUMENT = hsa_label, N_CONCEPTS = nrow(dd),
               SLOPE = round(s[2, 1], 5), SE = round(s[2, 2], 5),
               P_VALUE = round(s[2, 4], 4))
  }), fill = TRUE)
  
  by_family <- ok[, .(N_CONCEPTS = .N,
                      MEDIAN_DELTA = round(median(IV_DELTA_HSA_MINUS_HRR, na.rm = TRUE), 4),
                      SHARE_HSA_STRONGER = round(mean(IV_DELTA_HSA_MINUS_HRR < 0, na.rm = TRUE), 3)),
                  by = .(FINAL_FAMILY_ID, HSA_INSTRUMENT)][order(MEDIAN_DELTA)]
  
  list(delta = d, by_cat = by_cat, slopes = slopes, by_family = by_family)
}

xo_comp <- p12_crossover(cr_hsa_comp, "HSA_competitor")
xo_adj  <- p12_crossover(cr_hsa_adj,  "HSA_excl_adjacent")

save_csv(rbindlist(list(xo_comp$delta,     xo_adj$delta),     fill = TRUE),
         "MD37_concept_hsa_hrr_delta.csv")
save_csv(rbindlist(list(xo_comp$by_cat,    xo_adj$by_cat),    fill = TRUE),
         "MD38_crossover_by_scheme.csv")
save_csv(rbindlist(list(xo_comp$slopes,    xo_adj$slopes),    fill = TRUE),
         "MD39_crossover_slopes.csv")
save_csv(rbindlist(list(xo_comp$by_family, xo_adj$by_family), fill = TRUE),
         "MD40_crossover_by_family.csv")

cat("\n=== CROSS-OVER SLOPE ON ORDINAL SHOPPABILITY ===\n")
cat("A negative slope means more shoppable concepts favour HSA.\n\n")
print(rbindlist(list(xo_comp$slopes, xo_adj$slopes), fill = TRUE))

cat("\n=== CROSS-OVER BY CLINICAL FAMILY ===\n")
print(rbindlist(list(xo_comp$by_family, xo_adj$by_family), fill = TRUE))


# -----------------------------------------------------------------------------
# 12.F  Pooled interacted models, six primary schemes
# -----------------------------------------------------------------------------
# Runs p12_run_pooled_all_schemes() for the three runs
# (MD41_pooled_all_schemes_all_runs.csv) and tabulates IV_PERCENT by run,
# scheme, and category (MD42_pooled_gap_by_scheme.csv). This is fast relative
# to the concept-level loop. Pooling averages over roughly 738 concepts whose
# effects run in opposite directions, so a pooled coefficient near zero is
# expected even where a concept-level gradient exists. These models show the
# spread across schemes and are not read as significance tests.
p12_pooled_all <- rbindlist(list(
  p12_run_pooled_all_schemes(run_hsa_comp),
  p12_run_pooled_all_schemes(run_hsa_adj),
  p12_run_pooled_all_schemes(run_hrr_comp)), fill = TRUE)

save_csv(p12_pooled_all, "MD41_pooled_all_schemes_all_runs.csv")

pooled_gap <- dcast(p12_pooled_all, RUN_KEY + SCHEME_LABEL ~ TERM,
                    value.var = "IV_PERCENT", fun.aggregate = mean)
cat("\n=== POOLED SHOPPABILITY GAP BY SCHEME AND RUN ===\n"); print(pooled_gap)
save_csv(pooled_gap, "MD42_pooled_gap_by_scheme.csv")


# =============================================================================
# Part 13: Between-concept heterogeneity
# =============================================================================
#
# A shoppability gradient in the meta-regression requires variation in the
# concept-level coefficients beyond sampling error. Part 13 measures that
# variation at county, HSA, and HRR with the same function, het_decomp(), and
# the competitor-hospitals instrument at each geography: the
# Competitor_only_hospitals_9m rows of the cached county run
# (concept_results) and the Part 12 runs cr_hsa_comp and cr_hrr_comp. The
# stopifnot() below requires all three objects.
#
# Cochran's Q compares the observed dispersion with what sampling error alone
# would produce; its expected value without heterogeneity is df. I-squared is
# the share of total variance that lies between concepts. tau is the
# between-concept standard deviation in coefficient units (DerSimonian-Laird).
# TAU_OVER_SE, tau divided by the median standard error, is the summary
# measure: below about 0.3, the spread is mostly estimation noise.
#
# Median first-stage F rises from county to HSA to HRR, the same order in
# which I-squared falls. A weaker instrument gives noisier coefficients and
# more apparent heterogeneity, so the pattern partly reflects instrument
# strength rather than geography alone. For this reason MD44 repeats the
# decomposition under first-stage screens of 0, 5, 10, and 20. The
# per-concept cross-over in 12.E compares the same concepts across
# geographies and does not rely on a variance ratio.
#
# Writes MD43 (decomposition), MD44 (sensitivity to the first-stage screen),
# and MD45 (family medians); prints the variance share explained by family.

stopifnot("cr_hrr_comp not in scope"  = exists("cr_hrr_comp"),
          "cr_hsa_comp not in scope"  = exists("cr_hsa_comp"),
          "concept_results not in scope" = exists("concept_results"))

# het_decomp() summarizes one set of concept-level IV estimates: concepts with
# finite IV_COEF and IV_SE, IV_SE > 0, and FS_F >= fs_min (default 10), or
# NULL if fewer than 10 remain. It returns one row: the inverse-variance pooled
# coefficient (also in percent), Q with its df and chi-square p-value,
# I-squared, tau, the median SE, TAU_OVER_SE, the raw SD of the coefficients,
# and the median first-stage F. Parts 13, 15.D, and 15E all use it; 15.D and
# 15E call it through het_decomp_compact() and het_decomp_minimal(), which
# keep fewer columns.
het_decomp <- function(d, label, fs_min = 10) {
  d <- as.data.table(d)
  x <- d[is.finite(IV_COEF) & is.finite(IV_SE) & IV_SE > 0 & FS_F >= fs_min]
  if (nrow(x) < 10) return(NULL)
  
  w     <- 1 / x$IV_SE^2
  theta <- sum(w * x$IV_COEF) / sum(w)
  Q     <- sum(w * (x$IV_COEF - theta)^2)
  df    <- nrow(x) - 1L
  I2    <- max(0, (Q - df) / Q)
  tau2  <- max(0, (Q - df) / (sum(w) - sum(w^2) / sum(w)))
  
  data.table(
    GEOGRAPHY   = label,
    N_CONCEPTS  = nrow(x),
    POOLED_COEF = round(theta, 6),
    POOLED_PCT  = round(100 * (exp(theta) - 1), 4),
    Q           = round(Q, 1),
    DF          = df,
    Q_MINUS_DF  = round(Q - df, 1),
    Q_P_VALUE   = signif(pchisq(Q, df, lower.tail = FALSE), 3),
    I_SQUARED   = round(I2, 4),
    TAU         = round(sqrt(tau2), 6),
    MEDIAN_SE   = round(median(x$IV_SE), 6),
    TAU_OVER_SE = round(sqrt(tau2) / median(x$IV_SE), 3),
    OBSERVED_SD = round(sd(x$IV_COEF), 6),
    MEDIAN_FS_F = round(median(x$FS_F), 1))
}

# The cached county run has one row per concept and instrument. Only the
# Competitor_only_hospitals_9m rows are kept, to match the single-instrument
# HSA and HRR runs of Part 12.
cty_comp <- as.data.table(concept_results)[
  INSTRUMENT_LABEL == "Competitor_only_hospitals_9m"]

decomp <- rbindlist(list(
  het_decomp(cty_comp,    "COUNTY (competitor IV)"),
  het_decomp(cr_hsa_comp, "HSA (competitor IV)"),
  het_decomp(cr_hrr_comp, "HRR (competitor IV)")), fill = TRUE)

cat("\n=== VARIANCE DECOMPOSITION ===\n"); print(decomp)
save_csv(decomp, "MD43_heterogeneity_decomposition.csv")

sens <- rbindlist(lapply(c(0, 5, 10, 20), function(f) {
  rbindlist(list(
    het_decomp(cty_comp,    sprintf("COUNTY (FS_F>=%d)", f), fs_min = f),
    het_decomp(cr_hsa_comp, sprintf("HSA (FS_F>=%d)", f),    fs_min = f),
    het_decomp(cr_hrr_comp, sprintf("HRR (FS_F>=%d)", f),    fs_min = f)), fill = TRUE)
}), fill = TRUE)

cat("\n=== SENSITIVITY TO THE FIRST-STAGE SCREEN ===\n")
print(sens[, .(GEOGRAPHY, N_CONCEPTS, I_SQUARED, TAU, MEDIAN_SE, TAU_OVER_SE)])
save_csv(sens, "MD44_heterogeneity_sensitivity_fs_screen.csv")

# fam_table() summarizes the raw concept-level IV coefficients (FS_F >=
# fs_min) by clinical family, for families with at least min_concepts
# concepts. The medians are not shrunk: where I-squared is zero, empirical
# Bayes would return the pooled mean for every concept. Where I-squared is
# zero, differences between family medians are within sampling noise.
fam_table <- function(d, label, fs_min = 10, min_concepts = 5) {
  x <- as.data.table(d)[is.finite(IV_COEF) & FS_F >= fs_min]
  x[, .(GEOGRAPHY = label, N_CONCEPTS = .N,
        MEDIAN_IV      = round(median(IV_COEF), 5),
        MEDIAN_IV_PCT  = round(100 * (exp(median(IV_COEF)) - 1), 3),
        IQR_IV         = round(IQR(IV_COEF), 5),
        SHARE_NEGATIVE = round(mean(IV_COEF < 0), 3),
        SHARE_SIG_5PCT = round(mean(IV_P < 0.05, na.rm = TRUE), 3),
        MEDIAN_SE      = round(median(IV_SE), 5),
        MEDIAN_FS_F    = round(median(FS_F), 1)),
    by = FINAL_FAMILY_ID][N_CONCEPTS >= min_concepts][order(MEDIAN_IV)]
}

fam_all <- rbindlist(list(
  fam_table(cty_comp,    "COUNTY"),
  fam_table(cr_hsa_comp, "HSA"),
  fam_table(cr_hrr_comp, "HRR")), fill = TRUE)

cat("\n=== FAMILY MEDIANS BY GEOGRAPHY ===\n"); print(fam_all)
save_csv(fam_all, "MD45_family_medians_by_geography.csv")

# Share of the variance of the raw coefficients explained by clinical family
# (R-squared of IV_COEF on family dummies, FS_F >= 10), printed for each
# geography as a second check on the decomposition above.
for (nm in c("COUNTY", "HSA", "HRR")) {
  d <- switch(nm, COUNTY = cty_comp, HSA = as.data.table(cr_hsa_comp),
              HRR = as.data.table(cr_hrr_comp))
  d <- d[is.finite(IV_COEF) & FS_F >= 10]
  fit <- tryCatch(lm(IV_COEF ~ factor(FINAL_FAMILY_ID), data = d), error = function(e) NULL)
  if (is.null(fit)) next
  s <- summary(fit)
  cat(sprintf("%s: family explains R2 = %.4f (adj %.4f), F = %.2f, p = %.3g, n = %d\n",
              nm, s$r.squared, s$adj.r.squared, s$fstatistic[1],
              pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3], lower.tail = FALSE),
              nrow(d)))
}


# =============================================================================
# Part 14: Fixed-effect feasibility
# =============================================================================
#
# The MARKET_ID (market-by-concept) fixed effect absorbs level differences
# but is constant over the study window, so it cannot absorb a shock to the
# local competitive environment that arrives partway through. Part 14
# considers three finer alternatives: market by month, hospital by month, and
# system by month. For each geography, the loop below reports the share of
# market-month cells with a single distinct treatment value and the number
# of hospitals and systems with no within-unit instrument variation. These
# three diagnostics decide which of the finer fixed effects can be estimated.
# Uses md_county_panel and the pooled HSA and HRR panels of 11.5 and 10.6
# (md_hsa_pooled, md_hrr_pooled).
#
# An interacted fixed effect must exist as a column. estimate_interacted()
# passes fixed_effects through available_columns(), which keeps only names
# found in names(data), so a fixest interaction such as
# "ANALYSIS_MARKET^POST_MONTH" would be dropped without a warning and the
# model would run on BASELINE_FIXED_EFFECTS alone.

for (nm in c("COUNTY", "HSA", "HRR")) {
  p  <- switch(nm, COUNTY = md_county_panel, HSA = md_hsa_pooled, HRR = md_hrr_pooled)
  iv <- switch(nm, COUNTY = PRIMARY_INSTRUMENT,
               HSA = "HSA_COMP_HOSPITALS", HRR = "HRR_COMP_HOSPITALS")
  
  mkt <- p[, .(N = uniqueN(get(ENDOGENOUS_VARIABLE))), by = .(ANALYSIS_MARKET, POST_MONTH)]
  hos <- p[, .(SD = sd(get(iv), na.rm = TRUE)), by = HOSPITAL_ID]
  sys <- p[!is.na(SYSTEM_KEY), .(SD = sd(get(iv), na.rm = TRUE)), by = SYSTEM_KEY]
  
  cat(sprintf("\n%s\n  market-month cells with one distinct treatment value: %.3f\n",
              nm, mean(mkt$N == 1)))
  cat(sprintf("  hospitals with no within-hospital instrument variation: %d of %d\n",
              sum(hos$SD == 0, na.rm = TRUE), nrow(hos)))
  cat(sprintf("  systems with no within-system instrument variation: %d of %d\n",
              sum(sys$SD == 0, na.rm = TRUE), nrow(sys)))
}

# Market by month is not estimable: the treatment is defined at the market
# level and is near-constant within a market-month cell, so the fixed effect
# absorbs the regressor. Hospital by month is not estimable either, since the
# instrument changes at most a few times per hospital across the window.
# System by month is estimable at all three geographies; SYSTEM_MONTH_FE is
# added to the three panels in place.
#
# run_sysmonth() fits the headline interacted model twice, with the baseline
# fixed effects and with SYSTEM_MONTH_FE added. The second model also
# clusters by SYSTEM_KEY, so model_sample() drops rows with a missing
# SYSTEM_KEY from it and its sample can be smaller than the baseline's.
# Writes MD51_system_x_month_fe.csv.
md_county_panel[, SYSTEM_MONTH_FE := paste(SYSTEM_KEY, POST_MONTH, sep = "::")]
md_hsa_pooled[,   SYSTEM_MONTH_FE := paste(SYSTEM_KEY, POST_MONTH, sep = "::")]
md_hrr_pooled[,   SYSTEM_MONTH_FE := paste(SYSTEM_KEY, POST_MONTH, sep = "::")]

FE_SYSTEM_MONTH <- c(BASELINE_FIXED_EFFECTS, "SYSTEM_MONTH_FE")

run_sysmonth <- function(panel, iv_col, iv_label, geo_label) {
  t0 <- Sys.time()
  cat(sprintf("\n%s system-by-month starting %s\n", geo_label, format(t0, "%H:%M:%S")))
  
  base <- estimate_interacted(
    data = panel, moderator = MD_SCHEME, moderator_type = "categorical",
    outcome = PRIMARY_OUTCOME, instrument = iv_col,
    label = paste0(geo_label, "/baseline"), instrument_label = iv_label,
    moderator_label = MD_SCHEME)
  
  sysm <- estimate_interacted(
    data = panel, moderator = MD_SCHEME, moderator_type = "categorical",
    outcome = PRIMARY_OUTCOME, instrument = iv_col,
    label = paste0(geo_label, "/system_x_month"), instrument_label = iv_label,
    moderator_label = MD_SCHEME,
    fixed_effects = FE_SYSTEM_MONTH,
    clusters = available_columns(panel, c(BASELINE_CLUSTERS, "SYSTEM_KEY")))
  
  cat(sprintf("%s done in %.1fs\n", geo_label,
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  
  rbindlist(list(
    cbind(base$rows, SPEC = "baseline"),
    cbind(sysm$rows, SPEC = "system x month FE")), fill = TRUE)[, GEOGRAPHY := geo_label]
}

sysmonth_all <- rbindlist(list(
  run_sysmonth(md_county_panel, PRIMARY_INSTRUMENT,
               "Competitor_only_hospitals_9m", "COUNTY"),
  run_sysmonth(md_hsa_pooled, "HSA_COMP_HOSPITALS",
               "HSA_competitor_hospitals", "HSA"),
  run_sysmonth(md_hrr_pooled, "HRR_COMP_HOSPITALS",
               "HRR_competitor_hospitals", "HRR")), fill = TRUE)

save_csv(sysmonth_all, "MD51_system_x_month_fe.csv")
cat("\n=== BASELINE AGAINST SYSTEM BY MONTH ===\n")
print(sysmonth_all[, .(GEOGRAPHY, SPEC, TERM, IV_PERCENT = round(IV_PERCENT, 3),
                       IV_P = round(IV_P, 4),
                       WALD = round(FIRST_STAGE_WALD_THIS_EQ, 1), N_OBSERVATIONS)])


# =============================================================================
# Part 15: Six-instrument concept runs and the restricted county run
# =============================================================================
#
# Four analyses that extend Parts 12 and 13:
#
#   15.A  County concept level, restricted to markets with instrument
#         variation. The HSA runs are restricted this way and the cached
#         county run is not, so the Part 13 comparison mixes market
#         definition with sample composition. 15.A isolates the restriction.
#   15.B  HRR concept level with six instruments. The county gradient is a
#         six-instrument result, while the Part 12 HRR run uses one.
#   15.C  HSA concept level with six instruments, for the same reason.
#   15.D  The Part 13 decomposition on the three six-instrument runs, with
#         the instrument count held constant.
#
# Each run goes through p15_run() (15.1), which caches its results under keys
# that include the instrument count and writes MD52-MD54; 15.D writes MD55 and
# MD56.
#
# Runtime: six instruments cost roughly five times a single-instrument run,
# so about 8 to 12 hours per six-instrument geography and 2 to 4 hours for
# 15.A. A completed run is loaded from its cache and not repeated, and an
# interruption costs only the run in progress.


# -----------------------------------------------------------------------------
# 15.0  Working set
# -----------------------------------------------------------------------------
# Removes large objects these analyses do not need. On an 8 GB machine a few
# leftover panels are enough to push R into swap, which slows the concept
# loop by roughly an order of magnitude. md_hsa_pooled and md_hrr_pooled are
# rebuilt in 15.B and 15.C with the six-instrument sets.
#
# `outpatient` is removed here, but 15E.7 passes it to
# build_comparability_measures() inside a cache_or_run() call with the key
# "comparability_measures". That call works when
# 07_Cache/comparability_measures.rds exists, because the argument is then
# never evaluated. Without the cache file, or with USE_CACHE = FALSE, it
# stops with "object 'outpatient' not found".
rm(list = intersect(ls(), c("Dec16", "Jun21", "outpatient", "cbsa_panel",
                            "s15_payer_class_panel", "s15_payer",
                            "test_panel", "md_hsa_pooled", "md_hrr_pooled")))
invisible(gc())
cat("Largest objects after cleanup:\n")
print(head(sort(sapply(ls(), function(x) object.size(get(x))), decreasing = TRUE), 5))


# The six-instrument HSA set is the competitor and strict families built in
# Part 11, without the three adjacency-excluded variants. Both source objects
# are already in scope, so this is a merge rather than a rebuild.
hsa_six <- merge(hsa_comp, hsa_strict, by = c("HOSPITAL_ID", "POST_MONTH"))

HSA_SIX <- c(
  HSA_competitor_hospitals = "HSA_COMP_HOSPITALS",
  HSA_competitor_systems   = "HSA_COMP_SYSTEMS",
  HSA_competitor_hsas      = "HSA_COMP_GEOS",
  HSA_strict_hospitals     = "HSA_STRICT_HOSPITALS",
  HSA_strict_systems       = "HSA_STRICT_SYSTEMS",
  HSA_strict_hsas          = "HSA_STRICT_GEOS")

# HRR_INSTRUMENTS from Part 10.4 is already the six-instrument HRR set and is
# reused directly in 15.B.

cat(sprintf("\nHRR_COMP_HOSPITALS mean %.3f, share zero %.3f\n",
            mean(hrr_six$HRR_COMP_HOSPITALS), mean(hrr_six$HRR_COMP_HOSPITALS == 0)))
cat(sprintf("HSA_COMP_HOSPITALS mean %.3f, share zero %.3f\n",
            mean(hsa_six$HSA_COMP_HOSPITALS), mean(hsa_six$HSA_COMP_HOSPITALS == 0)))

# Stops before the long runs if any listed object, from this file or from the
# main pipeline, is missing.
md_p15_required <- c("md_load_panel", "md_schemes_long", "md_hospital_attrs",
                     "md_event_roster", "build_dynamic_competitor_instrument",
                     "hrr_six", "hsa_six", "md_county_panel", "MD_SCHEME",
                     "estimate_concept_level", "build_concept_panel",
                     "prepare_meta_input", "run_meta_regressions",
                     "cache_or_run", "CONCEPT_INSTRUMENTS", "MAIN_INSTRUMENTS",
                     "SUPPORTING_INSTRUMENTS")
md_p15_missing <- md_p15_required[!vapply(md_p15_required, exists, logical(1))]
if (length(md_p15_missing))
  stop("Part 15 prerequisites missing: ",
       paste(md_p15_missing, collapse = ", "), call. = FALSE)
cat("\nPart 15 prerequisites present.\n"); gc()


# -----------------------------------------------------------------------------
# 15.1  Configuration
# -----------------------------------------------------------------------------
# Cache and output keys encode the instrument count, following the
# concept_level_6inst convention in the main pipeline. cache_or_run() keys on
# the string alone, so without the count a later run with a different
# instrument set would return the earlier object and report a cache hit.
p15_key <- function(stem, tag, n_inst) sprintf("%s_%s_%dinst", stem, tag, n_inst)

# p15_run() prints the panel's size and fixed-effect cell structure, fits the
# concept-level models with estimate_concept_level(), and runs the
# reduced-form and IV meta-regressions on the result. It returns
# list(concept, meta_rf, meta_iv). Cache keys: p15_concept_<tag>_<n>inst,
# p15_meta_rf_<tag>_<n>inst, and p15_meta_iv_<tag>_<n>inst. Writes
# MD52_concept_level, MD53_meta_rf, and MD54_meta_iv, each with the suffix
# _<tag>_<n>inst.csv; partial files use the stem MD52_concept_<tag>_<n>inst.
p15_run <- function(panel, instruments, tag, geo_label) {
  t0 <- Sys.time()
  cat("\n", strrep("=", 70), "\n", sep = "")
  cat(sprintf("%s | %d instruments | started %s\n",
              geo_label, length(instruments), format(t0, "%H:%M:%S")))
  cat(strrep("=", 70), "\n", sep = "")
  
  fe_sizes <- panel[, .N, by = MARKET_ID]$N
  cat(sprintf("  rows %s, hospitals %s, markets %s, concepts %s\n",
              format(nrow(panel), big.mark = ","),
              format(uniqueN(panel$HOSPITAL_ID), big.mark = ","),
              format(uniqueN(panel$ANALYSIS_MARKET), big.mark = ","),
              format(uniqueN(panel$FINAL_CONCEPT_ID), big.mark = ",")))
  cat(sprintf("  FE cells %s, singleton share %.3f, mean rows per cell %.2f\n",
              format(length(fe_sizes), big.mark = ","),
              mean(fe_sizes == 1), mean(fe_sizes)))
  
  n_inst <- length(instruments)
  cr <- cache_or_run(
    p15_key("p15_concept", tag, n_inst),
    estimate_concept_level(
      build_concept_panel(panel, instruments = instruments),
      instruments = instruments,
      save_stem = sprintf("MD52_concept_%s_%dinst", tag, n_inst)))
  
  cr[, `:=`(GEOGRAPHY = geo_label, RUN_TAG = tag, N_INSTRUMENTS_USED = n_inst)]
  
  cat(sprintf("\nDone: %s rows, %s concepts, %.1f min\n",
              format(nrow(cr), big.mark = ","),
              format(uniqueN(cr$FINAL_CONCEPT_ID), big.mark = ","),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  
  mi <- prepare_meta_input(cr, md_schemes_long)
  m_rf <- cache_or_run(p15_key("p15_meta_rf", tag, n_inst),
                       run_meta_regressions(mi, md_schemes_long, dep = "RF_COEF", se = "RF_SE"))
  m_iv <- cache_or_run(p15_key("p15_meta_iv", tag, n_inst),
                       run_meta_regressions(mi, md_schemes_long, dep = "IV_COEF", se = "IV_SE"))
  
  m_rf[, `:=`(GEOGRAPHY = geo_label, RUN_TAG = tag, N_INSTRUMENTS_USED = n_inst)]
  m_iv[, `:=`(GEOGRAPHY = geo_label, RUN_TAG = tag, N_INSTRUMENTS_USED = n_inst)]
  
  suffix <- sprintf("_%s_%dinst.csv", tag, n_inst)
  save_csv(cr,   paste0("MD52_concept_level", suffix))
  save_csv(m_rf, paste0("MD53_meta_rf",       suffix))
  save_csv(m_iv, paste0("MD54_meta_iv",       suffix))
  
  list(concept = cr, meta_rf = m_rf, meta_iv = m_iv)
}


# -----------------------------------------------------------------------------
# 15.A  County, restricted to markets with instrument variation
# -----------------------------------------------------------------------------
# Uses CONCEPT_INSTRUMENTS, the same six instruments as the cached county run
# (concept_level_6inst), so apart from the merged-concept rows that Part 2
# leaves without a market key, the restriction is the only difference from
# that run. A change in the heterogeneity decomposition then reflects sample
# composition rather than instrument count or market definition.
#
# The restriction rule matches the HSA runs: a market is kept if
# PRIMARY_INSTRUMENT is nonzero for at least one of its hospital-months. Run
# tag county_nz, geography label COUNTY_RESTRICTED.

nz_county_markets <- unique(
  md_county_panel[get(PRIMARY_INSTRUMENT) > 0, ANALYSIS_MARKET])

md_county_nz <- md_county_panel[ANALYSIS_MARKET %in% nz_county_markets]

cat(sprintf("\nCounty restriction: %s to %s rows, %d to %d markets, %d to %d hospitals\n",
            format(nrow(md_county_panel), big.mark = ","),
            format(nrow(md_county_nz), big.mark = ","),
            uniqueN(md_county_panel$ANALYSIS_MARKET),
            uniqueN(md_county_nz$ANALYSIS_MARKET),
            uniqueN(md_county_panel$HOSPITAL_ID),
            uniqueN(md_county_nz$HOSPITAL_ID)))

county_nz <- p15_run(md_county_nz, CONCEPT_INSTRUMENTS, "county_nz", "COUNTY_RESTRICTED")
cat("\n=== COUNTY RESTRICTED, META-REGRESSION (IV) ===\n"); print(county_nz$meta_iv)

rm(md_county_nz); invisible(gc())


# -----------------------------------------------------------------------------
# 15.B  HRR, six instruments
# -----------------------------------------------------------------------------
# HRR_INSTRUMENTS is the six-instrument HRR set defined in Part 10.4. The HRR
# panel is not restricted. The run stops if any panel row has no instrument
# value after the merge.

md_hrr_pooled <- merge(md_load_panel("HRR"), hrr_six,
                       by = c("HOSPITAL_ID", "POST_MONTH"), all.x = TRUE)
stopifnot("unmatched instrument rows at HRR" =
            sum(is.na(md_hrr_pooled$HRR_COMP_HOSPITALS)) == 0)

hrr_6inst <- p15_run(md_hrr_pooled, HRR_INSTRUMENTS, "hrr", "HRR")
cat("\n=== HRR SIX INSTRUMENTS, META-REGRESSION (IV) ===\n"); print(hrr_6inst$meta_iv)

rm(md_hrr_pooled); invisible(gc())


# -----------------------------------------------------------------------------
# 15.C  HSA, six instruments
# -----------------------------------------------------------------------------
# Uses HSA_SIX, defined after 15.0. The panel is restricted to markets where
# HSA_COMP_HOSPITALS is nonzero for at least one hospital-month, matching
# Part 12. Single-hospital HSAs cannot generate instrument variation, so the
# restriction is forced by the geography rather than chosen.

md_hsa_pooled <- merge(md_load_panel("HSA"), hsa_six,
                       by = c("HOSPITAL_ID", "POST_MONTH"), all.x = TRUE)
stopifnot("unmatched instrument rows at HSA" =
            sum(is.na(md_hsa_pooled$HSA_COMP_HOSPITALS)) == 0)

nz_hsa_markets <- unique(md_hsa_pooled[HSA_COMP_HOSPITALS > 0, ANALYSIS_MARKET])
md_hsa_nz <- md_hsa_pooled[ANALYSIS_MARKET %in% nz_hsa_markets]

cat(sprintf("\nHSA restriction: %s to %s rows, %d to %d markets\n",
            format(nrow(md_hsa_pooled), big.mark = ","),
            format(nrow(md_hsa_nz), big.mark = ","),
            uniqueN(md_hsa_pooled$ANALYSIS_MARKET),
            uniqueN(md_hsa_nz$ANALYSIS_MARKET)))

hsa_6inst <- p15_run(md_hsa_nz, HSA_SIX, "hsa", "HSA")
cat("\n=== HSA SIX INSTRUMENTS, META-REGRESSION (IV) ===\n"); print(hsa_6inst$meta_iv)

rm(md_hsa_pooled, md_hsa_nz); invisible(gc())


# -----------------------------------------------------------------------------
# 15.D  Heterogeneity across the six-instrument runs
# -----------------------------------------------------------------------------
# Repeats the Part 13 decomposition on the three six-instrument runs. The
# instrument count is the same at every geography and the county run is
# restricted as the HSA run is, so the comparison isolates market definition
# more closely than Part 13 does. Each run has one row per concept and
# instrument, so the decomposition is computed for each instrument separately
# rather than pooled across them. Writes MD55_heterogeneity_six_instrument.csv
# and, when concept_results is in memory,
# MD56_heterogeneity_full_vs_restricted.csv.

# Selects the column set these tables report from het_decomp() in Part 13.
het_decomp_compact <- function(d, label, fs_min = 10) {
  r <- het_decomp(d, label, fs_min)
  if (is.null(r)) return(NULL)
  r[, .(GEOGRAPHY, N_CONCEPTS, POOLED_PCT, Q, DF, I_SQUARED, TAU,
        MEDIAN_SE, TAU_OVER_SE, MEDIAN_FS_F)]
}

p15_decomp <- rbindlist(lapply(
  list(list(county_nz$concept, "COUNTY_RESTRICTED"),
       list(hrr_6inst$concept, "HRR"),
       list(hsa_6inst$concept, "HSA")),
  function(z) {
    d <- as.data.table(z[[1]])
    rbindlist(lapply(unique(d$INSTRUMENT_LABEL), function(il)
      het_decomp_compact(d[INSTRUMENT_LABEL == il],
                         sprintf("%s / %s", z[[2]], il))), fill = TRUE)
  }), fill = TRUE)

save_csv(p15_decomp, "MD55_heterogeneity_six_instrument.csv")
cat("\n=== HETEROGENEITY, SIX-INSTRUMENT RUNS ===\n")
print(p15_decomp[, .(GEOGRAPHY, N_CONCEPTS, I_SQUARED, TAU, MEDIAN_SE,
                     TAU_OVER_SE, MEDIAN_FS_F)])

# Against the cached county run on the full sample. A difference here is
# attributable to the sample restriction alone, since instrument set,
# estimator, and market definition are all held fixed.
if (exists("concept_results")) {
  cty_full <- as.data.table(concept_results)
  full_decomp <- rbindlist(lapply(unique(cty_full$INSTRUMENT_LABEL), function(il)
    het_decomp_compact(cty_full[INSTRUMENT_LABEL == il],
                       sprintf("COUNTY_FULL / %s", il))), fill = TRUE)
  
  cat("\n=== COUNTY, FULL AGAINST RESTRICTED ===\n")
  print(rbindlist(list(
    full_decomp,
    p15_decomp[GEOGRAPHY %like% "COUNTY_RESTRICTED"]), fill = TRUE)[
      , .(GEOGRAPHY, N_CONCEPTS, I_SQUARED, TAU_OVER_SE, MEDIAN_FS_F)])
  
  save_csv(rbindlist(list(full_decomp, p15_decomp), fill = TRUE),
           "MD56_heterogeneity_full_vs_restricted.csv")
}

cat("\n=== PART 15 COMPLETE ===\n")


# =============================================================================
# Part 15E: Cross-geography comparison tables and figures
# =============================================================================
#
# Reads the Part 15 results back from TABLE_DIR (MD52-MD56) and builds the
# cross-geography comparisons: a forest plot of the shoppability gradient
# under each market definition, the concept-level distribution, the
# shoppable and non-shoppable levels stated separately, and the appendix
# table.
#
# Because the Part 15 results are read from CSV, this part can run in a fresh
# session once Part 15 has completed. It also needs the main pipeline's
# functions with concept_results and schemes_long (a warm start provides
# them), het_decomp() from Part 13, md_schemes_long (Part 2), and xo_comp
# from 12.E (for 15E.11).
#
# County comparator. 15E.1 and 15E.2 use the restricted county run
# (county_nz, 15.A), which is the comparator for the restricted-against-full
# heterogeneity comparison; MD60 is written from it. It is not the paper's
# county estimate. 15E.3 replaces the county side with the full
# six-instrument county run (concept_results and the cached
# meta_regressions_iv) that every other table in the paper uses, and 15E.4
# onward use that version.
#
# Two known problems stop this part when it is sourced (15E.7, wf_hsa; 15E.9,
# abs_anchor_dedup), and two calls work only with existing cache files
# (15E.3, 15E.7); see the comments there.

read_tbl <- function(f) as.data.table(fread(file.path(TABLE_DIR, f)))

county_nz_concept <- read_tbl("MD52_concept_level_county_nz_6inst.csv")
county_nz_meta_rf <- read_tbl("MD53_meta_rf_county_nz_6inst.csv")
county_nz_meta_iv <- read_tbl("MD54_meta_iv_county_nz_6inst.csv")

hrr_concept <- read_tbl("MD52_concept_level_hrr_6inst.csv")
hrr_meta_rf <- read_tbl("MD53_meta_rf_hrr_6inst.csv")
hrr_meta_iv <- read_tbl("MD54_meta_iv_hrr_6inst.csv")

hsa_concept <- read_tbl("MD52_concept_level_hsa_6inst.csv")
hsa_meta_rf <- read_tbl("MD53_meta_rf_hsa_6inst.csv")
hsa_meta_iv <- read_tbl("MD54_meta_iv_hsa_6inst.csv")

het_six        <- read_tbl("MD55_heterogeneity_six_instrument.csv")
het_full_vs_nz <- read_tbl("MD56_heterogeneity_full_vs_restricted.csv")

cat("=== LOAD CHECK ===\n")
for (nm in c("county_nz_concept", "hrr_concept", "hsa_concept")) {
  d <- get(nm)
  cat(sprintf("%-20s %6s rows | %3d concepts | %d instruments\n",
              nm, format(nrow(d), big.mark = ","),
              uniqueN(d$FINAL_CONCEPT_ID), uniqueN(d$INSTRUMENT_LABEL)))
}

cat("County instruments:\n"); print(unique(county_nz_concept$INSTRUMENT_LABEL))
cat("\nHRR instruments:\n");    print(unique(hrr_concept$INSTRUMENT_LABEL))
cat("\nHSA instruments:\n");    print(unique(hsa_concept$INSTRUMENT_LABEL))


cat("=== SIX-INSTRUMENT HETEROGENEITY, ALL THREE GEOGRAPHIES ===\n")
print(het_six[, .(GEOGRAPHY, N_CONCEPTS, POOLED_PCT, I_SQUARED, TAU,
                  TAU_OVER_SE, MEDIAN_FS_F)])

cat("\n=== COUNTY: FULL SAMPLE vs RESTRICTED ===\n")
print(het_full_vs_nz[GEOGRAPHY %like% "COUNTY",
                     .(GEOGRAPHY, N_CONCEPTS, I_SQUARED, TAU_OVER_SE, MEDIAN_FS_F)])


# -----------------------------------------------------------------------------
# 15E.1  Anchor instrument only
# -----------------------------------------------------------------------------
# One instrument per geography, the competitor-hospitals variant at each, so
# the three columns rest on the same exclusion-restriction argument. Prints
# the decomposition for each geography.
ANCHOR <- data.table(
  GEOGRAPHY = c("COUNTY_RESTRICTED", "HRR", "HSA"),
  INSTRUMENT_LABEL = c("Competitor_only_hospitals_9m",
                       "HRR_competitor_hospitals",
                       "HSA_competitor_hospitals"))

anchor_concept <- rbindlist(list(
  county_nz_concept[INSTRUMENT_LABEL == ANCHOR[1, INSTRUMENT_LABEL]],
  hrr_concept[INSTRUMENT_LABEL == ANCHOR[2, INSTRUMENT_LABEL]],
  hsa_concept[INSTRUMENT_LABEL == ANCHOR[3, INSTRUMENT_LABEL]]
), fill = TRUE)

# Selects the column set this table reports from het_decomp() in Part 13.
het_decomp_minimal <- function(d, label, fs_min = 10) {
  r <- het_decomp(d, label, fs_min)
  if (is.null(r)) return(NULL)
  r[, .(GEOGRAPHY, N_CONCEPTS, I_SQUARED, TAU_OVER_SE, MEDIAN_FS_F)]
}

anchor_decomp <- rbindlist(lapply(ANCHOR$GEOGRAPHY, function(g)
  het_decomp_minimal(anchor_concept[GEOGRAPHY == g], g)), fill = TRUE)

cat("\n=== ANCHOR INSTRUMENT ONLY ===\n")
print(anchor_decomp)


# -----------------------------------------------------------------------------
# 15E.2  Meta-regression terms, county_nz comparator
# -----------------------------------------------------------------------------
# extract_top_term() keeps, for each scheme and collapse, the non-intercept
# row with the most negative estimate. The minimum is used because the name
# of the shoppable-category term differs by collapse (CATHIGH for 3-tier,
# CATShoppable for the 2-tier collapses). The meta-regression tables hold
# both weightings (inverse variance and unweighted) and the selection does
# not filter on WEIGHTING, so the row kept is the more negative of the two;
# for the 3-tier collapse it can also be the INTERMEDIATE term. MD60 is
# written from these rows.
#
# extract_intercept() returns the intercept, which is the mean of the
# reference category (Non_shoppable for the 2-tier collapses, LOW for
# 3-tier), with its standard error from the same fit; one row per scheme,
# collapse, and weighting.
extract_top_term <- function(meta_iv_tbl, instrument_label) {
  d <- meta_iv_tbl[INSTRUMENT_LABEL == instrument_label & term != "(Intercept)"]
  d[, .SD[which.min(estimate)], by = .(SCHEME_ID, SCHEME, COLLAPSE)]
}

top_county <- extract_top_term(county_nz_meta_iv, ANCHOR[1, INSTRUMENT_LABEL])[, GEOGRAPHY := "COUNTY_RESTRICTED"]
top_hrr    <- extract_top_term(hrr_meta_iv,       ANCHOR[2, INSTRUMENT_LABEL])[, GEOGRAPHY := "HRR"]
top_hsa    <- extract_top_term(hsa_meta_iv,       ANCHOR[3, INSTRUMENT_LABEL])[, GEOGRAPHY := "HSA"]

meta_compare <- rbindlist(list(top_county, top_hrr, top_hsa), fill = TRUE)

meta_wide <- dcast(meta_compare, SCHEME + COLLAPSE ~ GEOGRAPHY,
                   value.var = c("estimate", "p.value"))

cat("\n=== SHOPPABILITY EFFECT BY SCHEME, ANCHOR INSTRUMENT, ALL THREE GEOGRAPHIES ===\n")
print(meta_wide)
save_csv(meta_wide, "MD60_anchor_meta_comparison.csv")

extract_intercept <- function(meta_iv_tbl, instrument_label) {
  d <- meta_iv_tbl[INSTRUMENT_LABEL == instrument_label & term == "(Intercept)"]
  d[, .(SCHEME_ID, SCHEME, COLLAPSE, estimate, std.error, p.value)]
}

int_county <- extract_intercept(county_nz_meta_iv, ANCHOR[1, INSTRUMENT_LABEL])[, GEOGRAPHY := "COUNTY_RESTRICTED"]
int_hrr    <- extract_intercept(hrr_meta_iv,       ANCHOR[2, INSTRUMENT_LABEL])[, GEOGRAPHY := "HRR"]
int_hsa    <- extract_intercept(hsa_meta_iv,       ANCHOR[3, INSTRUMENT_LABEL])[, GEOGRAPHY := "HSA"]

intercept_compare <- rbindlist(list(int_county, int_hrr, int_hsa), fill = TRUE)

intercept_sign <- intercept_compare[, .(
  N_SCHEMES  = .N,
  N_POSITIVE = sum(estimate > 0),
  N_SIG_POS  = sum(estimate > 0 & p.value < 0.05),
  N_SIG_NEG  = sum(estimate < 0 & p.value < 0.05)
), by = GEOGRAPHY]

cat("=== REFERENCE-CATEGORY (roughly 'non-shoppable') SIGN, PROPER TEST ===\n")
print(intercept_sign)


# -----------------------------------------------------------------------------
# 15E.3  County side rebuilt on the full sample
# -----------------------------------------------------------------------------
# Replaces the county_nz objects above with the full six-instrument county
# run for every table and figure that follows. ANCHOR is redefined with
# GEOGRAPHY = "COUNTY" in place of "COUNTY_RESTRICTED" to match.
#
# meta_regressions_iv is the main pipeline's stage 8 cache, computed from
# prepare_meta_input(concept_results, schemes_long). The call below passes
# concept_results itself, so it works only on a cache hit. On a cache miss,
# or with USE_CACHE = FALSE, run_meta_regressions() finds no CAT_ columns and
# stops with "No meta-regressions for dep = IV_COEF".
meta_iv_full_county <- cache_or_run("meta_regressions_iv",
                                    run_meta_regressions(concept_results, schemes_long,
                                                         dep = "IV_COEF", se = "IV_SE"))

ANCHOR <- data.table(
  GEOGRAPHY = c("COUNTY", "HRR", "HSA"),
  INSTRUMENT_LABEL = c("Competitor_only_hospitals_9m",
                       "HRR_competitor_hospitals",
                       "HSA_competitor_hospitals"))

top_county <- extract_top_term(meta_iv_full_county, ANCHOR[1, INSTRUMENT_LABEL])[, GEOGRAPHY := "COUNTY"]
top_hrr    <- extract_top_term(hrr_meta_iv,          ANCHOR[2, INSTRUMENT_LABEL])[, GEOGRAPHY := "HRR"]
top_hsa    <- extract_top_term(hsa_meta_iv,          ANCHOR[3, INSTRUMENT_LABEL])[, GEOGRAPHY := "HSA"]
meta_compare <- rbindlist(list(top_county, top_hrr, top_hsa), fill = TRUE)

int_county <- extract_intercept(meta_iv_full_county, ANCHOR[1, INSTRUMENT_LABEL])[, GEOGRAPHY := "COUNTY"]
int_hrr    <- extract_intercept(hrr_meta_iv,          ANCHOR[2, INSTRUMENT_LABEL])[, GEOGRAPHY := "HRR"]
int_hsa    <- extract_intercept(hsa_meta_iv,          ANCHOR[3, INSTRUMENT_LABEL])[, GEOGRAPHY := "HSA"]
intercept_compare <- rbindlist(list(int_county, int_hrr, int_hsa), fill = TRUE)

# anchor_concept is rebuilt the same way, with the full county run
# (concept_results) in place of county_nz_concept.
anchor_concept <- rbindlist(list(
  as.data.table(concept_results)[INSTRUMENT_LABEL == ANCHOR[1, INSTRUMENT_LABEL]][, GEOGRAPHY := "COUNTY"],
  hrr_concept[INSTRUMENT_LABEL == ANCHOR[2, INSTRUMENT_LABEL]][, GEOGRAPHY := "HRR"],
  hsa_concept[INSTRUMENT_LABEL == ANCHOR[3, INSTRUMENT_LABEL]][, GEOGRAPHY := "HSA"]
), fill = TRUE)

GEO_LEVELS  <- c("COUNTY", "HSA", "HRR")
GEO_COLOURS <- c(COUNTY = "#782F40", HSA = "#4472A8", HRR = "#8A8A8A")


# -----------------------------------------------------------------------------
# 15E.4  Five-scheme forest
# -----------------------------------------------------------------------------
# The shoppable-category effect from meta_compare (15E.3) for five main
# scheme-collapse pairs (MAIN_SCHEMES) at each geography, with intervals of
# 1.96 standard errors. The plot is printed; no file is written.
MAIN_SCHEMES <- data.table(
  SCHEME   = c("Theory-Based V2 (MRI nonshoppable)", "Imaging vs Procedural",
               "CMS Statutory Shoppable List", "Alt: Upfront Cash-Market Framework",
               "High vs Low Within Modality"),
  COLLAPSE = c("2-tier High vs rest", "2-tier High vs rest",
               "2-tier High vs rest", "2-tier Low vs rest", "2-tier extremes"),
  LABEL    = c("Theory-Based V2", "Imaging vs Procedural",
               "CMS Statutory List", "Upfront Cash-Market", "Within Modality"))

meta_main <- merge(meta_compare, MAIN_SCHEMES, by = c("SCHEME", "COLLAPSE"))
meta_main[, GEOGRAPHY := factor(GEOGRAPHY, levels = GEO_LEVELS)]

p_forest <- ggplot(meta_main, aes(x = LABEL, y = estimate, colour = GEOGRAPHY)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_pointrange(aes(ymin = estimate - 1.96 * std.error,
                      ymax = estimate + 1.96 * std.error),
                  position = position_dodge(width = 0.5), linewidth = 0.7, size = 0.5) +
  coord_flip() +
  scale_colour_manual(values = GEO_COLOURS, name = NULL) +
  labs(x = NULL, y = "Shoppable-category effect (log points)",
       title = "Shoppability gradient by market definition, five main schemes") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")
print(p_forest)


# -----------------------------------------------------------------------------
# 15E.5  Concept-level distribution
# -----------------------------------------------------------------------------
# Violin and box plots of the concept-level IV coefficients for the anchor
# instrument at each geography, concepts with FS_F >= 10. The y-axis is
# limited to the 2nd to 98th percentiles with coord_cartesian(), which drops
# no data. The plot is printed; no file is written.
anchor_fs10 <- anchor_concept[is.finite(IV_COEF) & FS_F >= 10]

p_dist <- ggplot(anchor_fs10, aes(x = factor(GEOGRAPHY, levels = GEO_LEVELS), y = IV_COEF)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_violin(fill = "grey92", colour = "grey45", linewidth = 0.4) +
  geom_boxplot(width = 0.10, outlier.shape = NA, linewidth = 0.5) +
  coord_cartesian(ylim = quantile(anchor_fs10$IV_COEF, c(0.02, 0.98), na.rm = TRUE)) +
  labs(x = NULL, y = "Concept-level IV coefficient (log points)",
       title = "Spread of concept-level effects, anchor instrument") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())
print(p_dist)


# -----------------------------------------------------------------------------
# 15E.6  All-schemes forest and slope chart
# -----------------------------------------------------------------------------
# One 2-tier collapse per scheme: "2-tier High vs rest", or "2-tier Low vs
# rest" for schemes without a High-vs-rest row. Writes
# explore_forest_all_schemes.pdf and explore_slope_schemes.pdf to FIGURE_DIR.
ALL_SCHEMES_2TIER <- meta_compare[COLLAPSE == "2-tier High vs rest" |
                                    (COLLAPSE == "2-tier Low vs rest" &
                                       !SCHEME %in% meta_compare[COLLAPSE == "2-tier High vs rest", unique(SCHEME)])]
ALL_SCHEMES_2TIER[, GEOGRAPHY := factor(GEOGRAPHY, levels = GEO_LEVELS)]
ALL_SCHEMES_2TIER[, SCHEME_SHORT := gsub("^(Alt: |CMS |High vs Low )", "", SCHEME)]

p_forest_all <- ggplot(ALL_SCHEMES_2TIER,
                       aes(x = reorder(SCHEME_SHORT, estimate), y = estimate, colour = GEOGRAPHY)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_pointrange(aes(ymin = estimate - 1.96 * std.error,
                      ymax = estimate + 1.96 * std.error),
                  position = position_dodge(width = 0.6), linewidth = 0.5, size = 0.35) +
  coord_flip() +
  scale_colour_manual(values = GEO_COLOURS, name = NULL) +
  labs(x = NULL, y = "Category effect (log points)") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "top",
        axis.text.y = element_text(size = 8))
print(p_forest_all)
ggsave(file.path(FIGURE_DIR, "explore_forest_all_schemes.pdf"), p_forest_all,
       width = 8, height = 10)

p_slope <- ggplot(ALL_SCHEMES_2TIER,
                  aes(x = GEOGRAPHY, y = estimate, group = SCHEME_SHORT)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_line(alpha = 0.5, colour = "grey40") +
  geom_point(aes(colour = p.value < 0.05), size = 2.2) +
  scale_colour_manual(values = c(`TRUE` = "#782F40", `FALSE` = "grey70"), name = "p < .05") +
  labs(x = NULL, y = "Category effect (log points)",
       title = "Shoppability category effect, by scheme, across market definitions") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")
print(p_slope)
ggsave(file.path(FIGURE_DIR, "explore_slope_schemes.pdf"), p_slope, width = 6.5, height = 5)


# -----------------------------------------------------------------------------
# 15E.7  Contracting-depth mechanism tests at HRR and HSA
# -----------------------------------------------------------------------------
# Runs the main pipeline's comparability tests on the six-instrument HRR and
# HSA concept results: run_comparability_within_family() (Section 9) and
# run_comparability_meta() (Section 10). Both functions save their tables
# under fixed names in TABLE_DIR, so these calls overwrite the main
# pipeline's county files T09D_comparability_within_family.csv,
# T10_comparability_meta_regression.csv, and T10B_horse_race_summary.csv, and
# through screen_moderators() QA09_moderator_screen.csv in QA_DIR; the files
# left on disk hold the HSA results.
#
# build_comparability_measures() needs `outpatient`, which 15.0 removes, so
# the next line works only with a cached comparability_measures.rds (see
# 15.0).

measures <- cache_or_run("comparability_measures", build_comparability_measures(outpatient))

wf_hrr <- run_comparability_within_family(hrr_concept, measures)
wf_hsa <- run_comparability_within_family(hsa_concept, measures)

# TIER is still "UNKNOWN" for every HRR and HSA row here (relabel_tier() is
# applied below), so this loop prints empty tables.
cat("\n=== CONTRACTING DEPTH, WITHIN-FAMILY, SPEC (b), MAIN TIER ===\n")
for (nm in c("wf_hrr", "wf_hsa")) {
  d <- get(nm)[term == "MODC" & grepl("^\\(b\\)", SPEC) & TIER == "MAIN" &
                 MODERATOR == "N_PAYERS_V2"]
  cat(sprintf("\n%s:\n", toupper(nm)))
  print(d[, .(INSTRUMENT_LABEL, estimate, std.error, p.value, N_CONCEPTS, N_FAMILIES)])
}

meta_hrr <- run_comparability_meta(hrr_concept, measures)
meta_hsa <- run_comparability_meta(hsa_concept, measures)


# instrument_tier() (main pipeline, PART 1.5) assigns tiers only to the six
# county instrument labels in MAIN_INSTRUMENTS, CONFIRMING_INSTRUMENTS, and
# DISCREPANT_INSTRUMENTS. Every HRR_* and HSA_* label gets the default
# "UNKNOWN", so a TIER == "MAIN" filter matches no row. relabel_tier() assigns
# tiers to those labels by the unit counted: hospitals MAIN, systems
# CONFIRMING, HRRs or HSAs OTHER_UNIT. The pipeline's lists are not changed.
relabel_tier <- function(dt) {
  dt <- copy(dt)
  dt[grepl("^HRR_|^HSA_", INSTRUMENT_LABEL), TIER :=
       fifelse(grepl("_competitor_hospitals$|_strict_hospitals$", INSTRUMENT_LABEL), "MAIN",
               fifelse(grepl("_competitor_systems$|_strict_systems$",     INSTRUMENT_LABEL), "CONFIRMING",
                       "OTHER_UNIT"))]
  dt
}

wf_hrr <- relabel_tier(wf_hrr)
wf_hsa <- relabel_tier

# The line above assigns the function relabel_tier to wf_hsa instead of the
# result of relabel_tier(wf_hsa). wf_hsa then holds a function, not the HSA
# results, and the wf_hsa pass of the last loop in 15E.7 stops with "object
# 'term' not found".


meta_hrr <- relabel_tier(meta_hrr)
meta_hsa <- relabel_tier(meta_hsa)

# The spec (b), MAIN-tier summary by moderator that run_comparability_meta()
# prints is empty for HRR and HSA labels; it is recomputed here after
# relabel_tier(). Unlike that summary, these counts pool the reduced-form and
# IV rows (no filter on DEPENDENT).
for (nm in c("meta_hrr", "meta_hsa")) {
  d <- get(nm)[grepl("_Z$", term) & SPEC == "(b) + family FE" & TIER == "MAIN"]
  cat(sprintf("\n=== %s, MAIN TIER, SPEC (b) ===\n", toupper(nm)))
  print(d[, .(N = .N, N_SIG_05 = sum(p.value < 0.05, na.rm = TRUE),
              N_SIGN_OK = sum(SIGN_AS_PREDICTED, na.rm = TRUE),
              MEDIAN_P = round(median(p.value, na.rm = TRUE), 4)), by = MODERATOR])
}

for (nm in c("wf_hrr", "wf_hsa")) {
  d <- get(nm)[term == "MODC" & grepl("^\\(b\\)", SPEC) & TIER == "MAIN" & MODERATOR == "N_PAYERS_V2"]
  cat(sprintf("\n%s, contracting depth, MAIN tier:\n", toupper(nm)))
  print(d[, .(INSTRUMENT_LABEL, estimate, std.error, p.value, N_CONCEPTS, N_FAMILIES)])
}


# -----------------------------------------------------------------------------
# 15E.8  Shoppable and non-shoppable levels stated separately
# -----------------------------------------------------------------------------
# run_meta_regressions_absolute() follows run_meta_regressions() (main
# pipeline, Section 8): the same formula (dep ~ CAT), weightings (inverse
# variance and unweighted, META_WEIGHTINGS), clustering by CLUSTER_FAMILY,
# sample screens (MIN_CONCEPTS_META, at least two categories), and 2-tier
# collapse rules. It omits the 3-tier collapse. Each model is fitted twice,
# with Non_shoppable and then Shoppable as the reference level, so each
# level's mean is the intercept of one fit and its standard error comes from
# that fit; no intercept-slope covariance is needed. Returns one row per
# scheme, collapse, instrument, and weighting with the NONSHOP_* and SHOP_*
# estimates, standard errors, and p-values.
run_meta_regressions_absolute <- function(mi, schemes_long, dep = "IV_COEF", se = "IV_SE") {
  rules <- list(
    `2-tier High vs rest` = function(x) factor(fifelse(x == "HIGH", "Shoppable", "Non_shoppable"),
                                               levels = c("Non_shoppable", "Shoppable")),
    `2-tier Low vs rest`  = function(x) factor(fifelse(x == "LOW", "Non_shoppable", "Shoppable"),
                                               levels = c("Non_shoppable", "Shoppable")),
    `2-tier extremes`     = function(x) factor(fcase(x == "HIGH", "Shoppable",
                                                     x == "LOW", "Non_shoppable",
                                                     default = NA_character_),
                                               levels = c("Non_shoppable", "Shoppable")))
  rows <- list()
  t0 <- Sys.time()
  schemes <- unique(schemes_long$SCHEME_ID)
  
  for (sid in schemes) {
    col <- paste0("CAT_", sid); if (!(col %in% names(mi))) next
    sname <- schemes_long[SCHEME_ID == sid][1L]$SCHEME_NAME %||% sid
    
    for (rl in names(rules)) for (il in unique(mi$INSTRUMENT_LABEL)) {
      d <- mi[INSTRUMENT_LABEL == il & is.finite(get(dep)) & is.finite(get(se)) & get(se) > 0]
      if (nrow(d) == 0L) next
      d[, CAT := droplevels(rules[[rl]](toupper(trimws(get(col)))))]
      d <- d[!is.na(CAT)]
      if (nrow(d) < MIN_CONCEPTS_META || uniqueN(d$CAT) < 2L) next
      
      for (w in META_WEIGHTINGS) {
        d[, W := if (w == "Inverse variance") 1 / (get(se)^2) else 1]
        
        fit_ns <- tryCatch(feols(as.formula(paste(dep, "~ CAT")), data = d, weights = ~W,
                                 cluster = ~CLUSTER_FAMILY, warn = FALSE, notes = FALSE),
                           error = function(e) NULL)
        d2 <- copy(d); d2[, CAT := relevel(CAT, ref = "Shoppable")]
        fit_s <- tryCatch(feols(as.formula(paste(dep, "~ CAT")), data = d2, weights = ~W,
                                cluster = ~CLUSTER_FAMILY, warn = FALSE, notes = FALSE),
                          error = function(e) NULL)
        if (is.null(fit_ns) || is.null(fit_s)) next
        
        ns <- tidy_fixest(fit_ns)[term == "(Intercept)"]
        sh <- tidy_fixest(fit_s)[term == "(Intercept)"]
        if (nrow(ns) == 0L || nrow(sh) == 0L) next
        
        rows[[length(rows) + 1L]] <- data.table(
          SCHEME_ID = sid, SCHEME = sname, COLLAPSE = rl, INSTRUMENT_LABEL = il,
          WEIGHTING = w, DEPENDENT = dep,
          NONSHOP_EST = ns$estimate, NONSHOP_SE = ns$std.error, NONSHOP_P = ns$p.value,
          SHOP_EST    = sh$estimate, SHOP_SE    = sh$std.error, SHOP_P    = sh$p.value,
          N_CONCEPTS = nrow(d), N_CLUSTERS = uniqueN(d$CLUSTER_FAMILY))
      }
    }
    cat(sprintf("  scheme %d/%d (%s) | rows so far %d | %.1fs elapsed\n",
                which(schemes == sid), length(schemes), sname, length(rows),
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
  rbindlist(rows, fill = TRUE)
}

# prepare_meta_input() rebuilds the meta-input at each geography (a merge, not
# a refit). The absolute-level regressions run on the IV coefficients for
# every instrument and both weightings; abs_anchor keeps the anchor
# instrument with inverse-variance weights and is written to MD61.
mi_county <- prepare_meta_input(concept_results, schemes_long)
mi_hrr    <- prepare_meta_input(hrr_concept,      schemes_long)
mi_hsa    <- prepare_meta_input(hsa_concept,       schemes_long)

abs_county <- run_meta_regressions_absolute(mi_county, schemes_long, dep = "IV_COEF", se = "IV_SE")[, GEOGRAPHY := "COUNTY"]
abs_hrr    <- run_meta_regressions_absolute(mi_hrr,    schemes_long, dep = "IV_COEF", se = "IV_SE")[, GEOGRAPHY := "HRR"]
abs_hsa    <- run_meta_regressions_absolute(mi_hsa,    schemes_long, dep = "IV_COEF", se = "IV_SE")[, GEOGRAPHY := "HSA"]

abs_all <- rbindlist(list(abs_county, abs_hrr, abs_hsa), fill = TRUE)

ANCHOR_INST <- c(COUNTY = "Competitor_only_hospitals_9m",
                 HRR    = "HRR_competitor_hospitals",
                 HSA    = "HSA_competitor_hospitals")

abs_anchor <- abs_all[WEIGHTING == "Inverse variance" &
                        INSTRUMENT_LABEL == ANCHOR_INST[GEOGRAPHY]]

save_csv(abs_anchor, "MD61_shoppable_nonshoppable_by_geography.csv")
cat("Rows:", nrow(abs_anchor), "| schemes:", uniqueN(abs_anchor$SCHEME),
    "| geographies:", uniqueN(abs_anchor$GEOGRAPHY), "\n")


# -----------------------------------------------------------------------------
# 15E.9  Two-panel forest, shoppable against non-shoppable
# -----------------------------------------------------------------------------
# abs_anchor_dedup, used here and in 15E.10, is not created in this file or
# in the main pipeline, so this block stops with "object 'abs_anchor_dedup'
# not found". The code expects one row per SCHEME and GEOGRAPHY (see the
# check on scheme_order below); abs_anchor, the nearest object, has one row
# per scheme, collapse, and geography. Writes
# fig26_shoppable_nonshoppable_geography.pdf to FIGURE_DIR.
abs_long <- rbindlist(list(
  abs_anchor_dedup[, .(SCHEME, COLLAPSE, GEOGRAPHY, CATEGORY = "Shoppable",
                       estimate = SHOP_EST, std.error = SHOP_SE, p.value = SHOP_P)],
  abs_anchor_dedup[, .(SCHEME, COLLAPSE, GEOGRAPHY, CATEGORY = "Non-shoppable",
                       estimate = NONSHOP_EST, std.error = NONSHOP_SE, p.value = NONSHOP_P)]
), fill = TRUE)

abs_long[, SCHEME_SHORT := gsub("^(Alt: |CMS |High vs Low )", "", SCHEME)]
abs_long[, GEOGRAPHY := factor(GEOGRAPHY, levels = GEO_LEVELS)]
abs_long[, CATEGORY  := factor(CATEGORY, levels = c("Non-shoppable", "Shoppable"))]

# Checks that stripping the prefixes did not give two different schemes the
# same SCHEME_SHORT. If it did, the full SCHEME names are used instead.
dup_check <- abs_long[, uniqueN(SCHEME), by = SCHEME_SHORT][V1 > 1]
if (nrow(dup_check) > 0) {
  cat("WARNING: these short labels still collide across different schemes:\n")
  print(dup_check)
  cat("Using the full SCHEME name instead of SCHEME_SHORT for ordering.\n")
  abs_long[, SCHEME_SHORT := SCHEME]
}

scheme_order <- abs_long[GEOGRAPHY == "COUNTY" & CATEGORY == "Shoppable"][order(estimate), SCHEME_SHORT]
stopifnot("scheme_order still has duplicates" = !anyDuplicated(scheme_order))
abs_long[, SCHEME_SHORT := factor(SCHEME_SHORT, levels = scheme_order)]

p_forest_split <- ggplot(abs_long, aes(x = SCHEME_SHORT, y = estimate, colour = GEOGRAPHY)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_pointrange(aes(ymin = estimate - 1.96 * std.error, ymax = estimate + 1.96 * std.error),
                  position = position_dodge(width = 0.6), linewidth = 0.5, size = 0.35) +
  coord_flip() +
  facet_wrap(~ CATEGORY, ncol = 2) +
  scale_colour_manual(values = GEO_COLOURS, name = NULL) +
  labs(x = NULL, y = "Concept-level IV coefficient (log points)") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "top",
        axis.text.y = element_text(size = 8), strip.text = element_text(face = "bold"))
print(p_forest_split)

ggsave(file.path(FIGURE_DIR, "fig26_shoppable_nonshoppable_geography.pdf"),
       p_forest_split, width = 9, height = 7.5)


# -----------------------------------------------------------------------------
# 15E.10  Appendix table and summary counts
# -----------------------------------------------------------------------------
# Writes the body rows of the appendix table to tab_schemes_geography_body.tex
# in TABLE_DIR (not OUT_TEX): scheme, geography, shoppable estimate and
# p-value, non-shoppable estimate and p-value, the gap, and the number of
# concepts. Stars mark p < 0.10, 0.05, and 0.01. Uses abs_anchor_dedup (see
# 15E.9).
stars <- function(p) fifelse(is.na(p), "",
                             fifelse(p < 0.01, "$^{***}$",
                                     fifelse(p < 0.05, "$^{**}$",
                                             fifelse(p < 0.10, "$^{*}$", ""))))

tex <- copy(abs_anchor_dedup)
tex[, GEOGRAPHY := factor(GEOGRAPHY, levels = c("COUNTY", "HSA", "HRR"))]
tex[, GAP := SHOP_EST - NONSHOP_EST]
setorder(tex, SCHEME, GEOGRAPHY)
tex[, SCHEME_SHOWN := fifelse(duplicated(SCHEME), "", SCHEME)]

body <- tex[, sprintf(
  "%s & %s & %.5f%s & %.4f & %.5f%s & %.4f & %.5f & %d \\\\",
  gsub("&", "\\\\&", SCHEME_SHOWN), as.character(GEOGRAPHY),
  SHOP_EST, stars(SHOP_P), SHOP_P,
  NONSHOP_EST, stars(NONSHOP_P), NONSHOP_P,
  GAP, N_CONCEPTS)]

writeLines(body, file.path(TABLE_DIR, "tab_schemes_geography_body.tex"))
cat("Wrote", length(body), "rows to tab_schemes_geography_body.tex\n")

# Summary counts by geography -------------------------------------------------
cat("\n=== SUMMARY FOR THE WRITE-UP ===\n")
print(tex[, .(N_SCHEMES      = .N,
              SHOP_NEGATIVE  = sum(SHOP_EST < 0),
              SHOP_SIG_NEG   = sum(SHOP_EST < 0 & SHOP_P < 0.05),
              NONSHOP_SIG    = sum(NONSHOP_P < 0.05),
              GAP_NEGATIVE   = sum(GAP < 0),
              MEDIAN_SHOP    = round(median(SHOP_EST), 5),
              MEDIAN_NONSHOP = round(median(NONSHOP_EST), 5),
              MEDIAN_GAP     = round(median(GAP), 5)), by = GEOGRAPHY])

cat("\n=== SERVICE COUNTS BY COLLAPSE (the caveat worth checking) ===\n")
print(tex[GEOGRAPHY == "COUNTY", .(SCHEME, COLLAPSE, N_CONCEPTS)][order(N_CONCEPTS)])


# -----------------------------------------------------------------------------
# 15E.11  Robustness on the cross-over slope
# -----------------------------------------------------------------------------
# Checks the 12.E cross-over slope (HSA minus HRR coefficient on ordinal
# shoppability) for the HSA competitor run. Uses xo_comp$delta (12.E) and
# md_schemes_long (Part 2); no concept-level model is refitted, so this runs
# in seconds. For each of the first six SCHEME_IDs, on concepts with FS_F >=
# 10 at both HSA and HRR (at least 30 concepts), it refits the slope leaving
# out one clinical family at a time and computes a two-sided permutation
# p-value from 2,000 random shuffles of the ordinal rank (set.seed(1)).
# Writes MD57 (summary) and MD58 (leave-one-family-out detail).
crossover_robustness <- function(delta_table, hsa_label, schemes_to_check) {
  
  d <- as.data.table(delta_table)[FS_F_HSA >= 10 & FS_F_HRR >= 10 &
                                    is.finite(IV_DELTA_HSA_MINUS_HRR)]
  
  scheme_rows <- list()
  lofo_all    <- list()
  
  for (sid in schemes_to_check) {
    ord <- unique(md_schemes_long[SCHEME_ID == sid,
                                  .(FINAL_CONCEPT_ID = ANALYSIS_CONCEPT_ID,
                                    ORD = SHOPPABILITY_ORDINAL)])
    dd <- merge(d, ord, by = "FINAL_CONCEPT_ID")[is.finite(ORD)]
    if (nrow(dd) < 30 || !has_usable_variation(dd$ORD)) next
    
    fit0 <- lm(IV_DELTA_HSA_MINUS_HRR ~ ORD, data = dd)
    slope_obs <- unname(coef(fit0)[2])
    
    fam_check <- rbindlist(lapply(unique(dd$FINAL_FAMILY_ID), function(fam) {
      sub <- dd[FINAL_FAMILY_ID != fam]
      if (nrow(sub) < 30 || !has_usable_variation(sub$ORD)) return(NULL)
      fit <- lm(IV_DELTA_HSA_MINUS_HRR ~ ORD, data = sub)
      s <- summary(fit)$coefficients
      data.table(SCHEME_ID = sid, HSA_INSTRUMENT = hsa_label,
                 DROPPED_FAMILY = fam, N_REMAINING = nrow(sub),
                 SLOPE = round(s[2, 1], 5), P_VALUE = round(s[2, 4], 4))
    }), fill = TRUE)
    
    set.seed(1)
    n_perm <- 2000
    perm_slopes <- vapply(seq_len(n_perm), function(i) {
      dd_p <- copy(dd)
      dd_p[, ORD := sample(ORD)]
      unname(coef(lm(IV_DELTA_HSA_MINUS_HRR ~ ORD, data = dd_p))[2])
    }, numeric(1))
    perm_p <- mean(abs(perm_slopes) >= abs(slope_obs))
    
    scheme_rows[[sid]] <- data.table(
      SCHEME_ID = sid, HSA_INSTRUMENT = hsa_label, N_CONCEPTS = nrow(dd),
      OBSERVED_SLOPE = round(slope_obs, 5), PERMUTATION_P = round(perm_p, 4),
      MIN_LOFO_SLOPE = round(min(fam_check$SLOPE), 5),
      MAX_LOFO_SLOPE = round(max(fam_check$SLOPE), 5),
      SHARE_LOFO_STILL_SIG = round(mean(fam_check$P_VALUE < 0.05), 3))
    lofo_all[[sid]] <- fam_check
  }
  
  list(summary = rbindlist(scheme_rows, fill = TRUE),
       lofo    = rbindlist(lofo_all, fill = TRUE))
}

key_schemes <- unique(md_schemes_long$SCHEME_ID)[1:6]

rob_comp <- crossover_robustness(xo_comp$delta, "HSA_competitor", key_schemes)

save_csv(rob_comp$summary, "MD57_crossover_robustness_summary.csv")
save_csv(rob_comp$lofo,    "MD58_crossover_robustness_lofo_detail.csv")

cat("\n=== CROSS-OVER ROBUSTNESS: HSA_competitor ===\n")
print(rob_comp$summary)


# =============================================================================
# Part 16: Distance-banded instruments and the ring decay IV
# =============================================================================
#
# The three-band ring decomposition of Part 4 (4.2 and 4.3) is estimated by
# OLS, because three endogenous band terms need three excluded instruments
# and the design supplies one. Part 16 builds the three band instruments,
# estimates the three-band IV next to OLS (16.4), and plots the OLS estimates
# (16.5).
#
# Construction. The instrument's exclusion argument requires peers outside
# the focal market, while the ring decomposition is about local reach, and
# most hospitals within fifteen miles of a focal hospital are in its own
# county. The band instruments are therefore county-gated and out-of-county.
# A peer counts for a focal hospital-month if its system had already posted
# in the focal hospital's county (the LOCAL_SYSTEM_FIRST_POST_MONTH gate of
# IV06), it is in a different county, it belongs to a different system, it
# posted within the trailing nine months, and it lies in the distance band.
# The near band is thinner than the far bands by construction; 16.3 reports
# how thin before anything is estimated. If the 0-15 band instrument is zero
# for most hospital-months, the three-band IV is not identified and the OLS
# decomposition is the one to report.
#
# Treatment. The endogenous variables are the SQL ring counts
# (RING_PRIOR_0_15, RING_PRIOR_15_30, RING_PRIOR_30_60), which count all
# prior posters in the band, including same-county and same-system ones.
# This parallels the headline design, where the treatment counts all prior
# posters in the county and the instrument counts only out-of-county
# competitors.
#
# Checks that stop the run: coordinate columns found (16.1); no duplicate
# roster or output rows (16.2); the band sum never exceeds IV06 (16.3); the
# instrument merge keeps the row count and the ring counts are present
# (16.4); three OLS rows (16.5).
#
# Needs md_county_panel (with the RING_PRIOR_* columns from Part 1),
# md_event_roster (10.1), and md_hosp_geo (Part 1). Writes
# MD59_ring_decay_iv.csv, MD59_ring_decay_ols.csv, and fig25_ring_decay.pdf.

MILES_TO_M <- 1609.344
RING_CUTS  <- c(0, 15, 30, 60)          # miles
RING_LABS  <- c("0_15", "15_30", "30_60")


# -----------------------------------------------------------------------------
# 16.1  Hospital coordinates and pairwise distances
# -----------------------------------------------------------------------------
# Candidate pairs within 60 miles come from st_is_within_distance() on points
# projected to Albers Equal Area (EPSG 5070), which sf can index. Distances
# for the candidate pairs are then recomputed geodesically with s2 and
# assigned to the bands [0, 15), [15, 30), and [30, 60] miles. EPSG 5070
# distorts Alaska and Hawaii, so the candidate search there is approximate:
# the geodesic pass corrects the distances of the candidates but cannot add
# pairs the projected search missed. s2 is switched on only for the distance
# step (points raise no validity problems) and off again afterwards, the
# setting used for the polygon work in Part 11.
#
# Column names in HPT_HOSPITAL_GEO vary by export vintage, so the latitude
# and longitude columns are detected rather than assumed.
lat_col <- grep("^(LAT|LATITUDE|HQ_LATITUDE)$",  names(md_hosp_geo), value = TRUE)[1]
lon_col <- grep("^(LON|LONG|LONGITUDE|HQ_LONGITUDE)$", names(md_hosp_geo), value = TRUE)[1]
stopifnot("no latitude column found in md_hosp_geo"  = !is.na(lat_col),
          "no longitude column found in md_hosp_geo" = !is.na(lon_col))
cat(sprintf("Coordinates: %s / %s\n", lat_col, lon_col))

hosp_xy <- unique(md_hosp_geo[, .(HOSPITAL_ID,
                                  LAT = as.numeric(get(lat_col)),
                                  LON = as.numeric(get(lon_col)),
                                  COUNTY_STATE_KEY_GEO = COUNTY_FIPS)],
                  by = "HOSPITAL_ID")
hosp_xy <- hosp_xy[is.finite(LAT) & is.finite(LON)]
cat(sprintf("Hospitals with usable coordinates: %d of %d\n",
            nrow(hosp_xy), uniqueN(md_hosp_geo$HOSPITAL_ID)))

# Candidate pairs within 60 miles on the projected CRS.
sf_use_s2(FALSE)
pts_proj <- st_transform(
  st_as_sf(hosp_xy, coords = c("LON", "LAT"), crs = 4326, remove = FALSE),
  5070)

cat("Finding candidate pairs within 60 miles ...\n")
nb <- st_is_within_distance(pts_proj, dist = max(RING_CUTS) * MILES_TO_M)

pairs <- rbindlist(lapply(seq_along(nb), function(i) {
  j <- nb[[i]][nb[[i]] != i]
  if (!length(j)) return(NULL)
  data.table(FOCAL_HOSPITAL_ID = hosp_xy$HOSPITAL_ID[i],
             PEER_HOSPITAL_ID  = hosp_xy$HOSPITAL_ID[j])
}))
cat(sprintf("Candidate pairs: %s\n", format(nrow(pairs), big.mark = ",")))

# Geodesic distances for the candidate pairs (s2 on for this step only), then
# the band assignment.
pairs <- merge(pairs, hosp_xy[, .(FOCAL_HOSPITAL_ID = HOSPITAL_ID,
                                  F_LON = LON, F_LAT = LAT)],
               by = "FOCAL_HOSPITAL_ID")
pairs <- merge(pairs, hosp_xy[, .(PEER_HOSPITAL_ID = HOSPITAL_ID,
                                  P_LON = LON, P_LAT = LAT)],
               by = "PEER_HOSPITAL_ID")

sf_use_s2(TRUE)
pairs[, DIST_MI := as.numeric(st_distance(
  st_as_sf(data.table(x = F_LON, y = F_LAT), coords = c("x", "y"), crs = 4326),
  st_as_sf(data.table(x = P_LON, y = P_LAT), coords = c("x", "y"), crs = 4326),
  by_element = TRUE)) / MILES_TO_M]
sf_use_s2(FALSE)

pairs <- pairs[DIST_MI <= max(RING_CUTS)]
pairs[, BAND := cut(DIST_MI, breaks = RING_CUTS, labels = RING_LABS,
                    include.lowest = TRUE, right = FALSE)]
pairs <- pairs[!is.na(BAND)]
pairs[, c("F_LON", "F_LAT", "P_LON", "P_LAT") := NULL]

cat(sprintf("Pairs within 60 miles: %s\n", format(nrow(pairs), big.mark = ",")))
print(pairs[, .(PAIRS = .N, MEDIAN_MI = round(median(DIST_MI), 1)), by = BAND][order(BAND)])

invisible(gc())


# -----------------------------------------------------------------------------
# 16.2  Band-specific competitor instruments
# -----------------------------------------------------------------------------
# build_ring_instruments() applies the gates of
# build_dynamic_competitor_instrument() (Part 10.2) with presence = "posted"
# and exclude_own_system = TRUE. The peer universe comes from the distance
# pairs instead of a geography, so the out-of-market restriction is applied
# on county (COUNTY_STATE_KEY) and the band on distance. The four stages are
# marked below. Returns one row per focal hospital-month with Z_RING_0_15,
# Z_RING_15_30, and Z_RING_30_60, the number of distinct peer hospitals in
# each band (zero where there are none).
build_ring_instruments <- function(focal_panel, roster, pair_table,
                                   lookback_months = 9) {
  
  r <- copy(roster)[, .(HOSPITAL_ID, SYSTEM_KEY, PEER_POST_MONTH,
                        GEO = COUNTY_STATE_KEY)]
  stopifnot("duplicate HOSPITAL_ID in roster" = !any(duplicated(r$HOSPITAL_ID)))
  
  # Stage 1: earliest month each system posted in each county. Identical to
  # the county-level builder.
  local_system_roster <- r[!is.na(GEO) & !is.na(SYSTEM_KEY) & !is.na(PEER_POST_MONTH),
                           .(LOCAL_SYSTEM_FIRST_POST_MONTH = min(PEER_POST_MONTH)),
                           by = .(FOCAL_GEO = GEO, SYSTEM_KEY)]
  
  # Stage 2: attach roster attributes to both sides of each pair, then drop
  # same-county pairs. This is the out-of-market restriction.
  p <- merge(pair_table,
             r[, .(FOCAL_HOSPITAL_ID = HOSPITAL_ID, FOCAL_GEO = GEO,
                   FOCAL_SYSTEM = SYSTEM_KEY)],
             by = "FOCAL_HOSPITAL_ID")
  p <- merge(p,
             r[, .(PEER_HOSPITAL_ID = HOSPITAL_ID, PEER_GEO = GEO,
                   PEER_SYSTEM = SYSTEM_KEY, PEER_POST_MONTH)],
             by = "PEER_HOSPITAL_ID")
  
  n_pre <- nrow(p)
  p <- p[!is.na(FOCAL_GEO) & !is.na(PEER_GEO) & FOCAL_GEO != PEER_GEO]
  cat(sprintf("  out-of-county restriction: %s to %s pairs (%.1f%% dropped)\n",
              format(n_pre, big.mark = ","), format(nrow(p), big.mark = ","),
              100 * (1 - nrow(p) / n_pre)))
  
  # Stage 3: the local-presence gate. Inner join, so a peer whose system has
  # never posted in the focal county is excluded entirely.
  p <- merge(p, local_system_roster,
             by.x = c("FOCAL_GEO", "PEER_SYSTEM"),
             by.y = c("FOCAL_GEO", "SYSTEM_KEY"))
  cat(sprintf("  after local-system gate: %s pairs, %d focal hospitals\n",
              format(nrow(p), big.mark = ","), uniqueN(p$FOCAL_HOSPITAL_ID)))
  
  # Stage 4: expand to focal hospital-months and apply the time and
  # competitor gates.
  focal_events <- unique(focal_panel[, .(FOCAL_HOSPITAL_ID = HOSPITAL_ID, POST_MONTH)])
  j <- merge(focal_events, p, by = "FOCAL_HOSPITAL_ID", allow.cartesian = TRUE)
  
  j[, DYNAMIC_OK := LOCAL_SYSTEM_FIRST_POST_MONTH <= POST_MONTH]
  j[, LOWER      := POST_MONTH %m-% months(lookback_months)]
  j[, STRICT     := PEER_POST_MONTH >= LOWER & PEER_POST_MONTH < POST_MONTH]
  j[, COMPETITOR := is.na(FOCAL_SYSTEM) | (PEER_SYSTEM != FOCAL_SYSTEM)]
  
  counts <- j[DYNAMIC_OK & STRICT & COMPETITOR,
              .(N = uniqueN(PEER_HOSPITAL_ID)),
              by = .(HOSPITAL_ID = FOCAL_HOSPITAL_ID, POST_MONTH, BAND)]
  
  wide <- dcast(counts, HOSPITAL_ID + POST_MONTH ~ BAND,
                value.var = "N", fill = 0L)
  setnames(wide, RING_LABS, paste0("Z_RING_", RING_LABS))
  
  out <- merge(unique(focal_panel[, .(HOSPITAL_ID, POST_MONTH)]), wide,
               by = c("HOSPITAL_ID", "POST_MONTH"), all.x = TRUE)
  for (b in paste0("Z_RING_", RING_LABS))
    if (!(b %in% names(out))) out[, (b) := 0L] else out[is.na(get(b)), (b) := 0L]
  
  stopifnot("duplicate output rows" =
              !any(duplicated(out, by = c("HOSPITAL_ID", "POST_MONTH"))))
  out[]
}

cat("\n=== Building distance-banded instruments ===\n")
ring_iv <- build_ring_instruments(md_county_panel, md_event_roster, pairs)


# -----------------------------------------------------------------------------
# 16.3  Diagnostics before estimating
# -----------------------------------------------------------------------------
# The near band is thin by construction, since most hospitals within fifteen
# miles share the focal hospital's county and are excluded. Reports the
# distribution of each band instrument and their correlations, checks that
# the bands nest within IV06's peers, and flags any band that is zero in
# more than 95% of rows.
ring_desc <- rbindlist(lapply(RING_LABS, function(b) {
  x <- ring_iv[[paste0("Z_RING_", b)]]
  data.table(BAND = b, MEAN = round(mean(x), 3), MEDIAN = median(x),
             SD = round(sd(x), 3), SHARE_ZERO = round(mean(x == 0), 3),
             P90 = quantile(x, 0.90), MAX = max(x))
}))
cat("\n=== BAND INSTRUMENTS ===\n"); print(ring_desc)

cat("\n=== CORRELATION AMONG BAND INSTRUMENTS ===\n")
print(round(cor(as.matrix(ring_iv[, paste0("Z_RING_", RING_LABS), with = FALSE])), 3))

# Construction check: the bands partition a subset of IV06's peers, so their
# sum cannot exceed IV06. The run stops if it does.
chk <- merge(ring_iv,
             unique(md_county_panel[, .(HOSPITAL_ID, POST_MONTH,
                                        IV06 = Z_SYS_COMPETITOR_ONLY_9M_EXCL_CURRENT)]),
             by = c("HOSPITAL_ID", "POST_MONTH"))
chk[, BAND_SUM := Z_RING_0_15 + Z_RING_15_30 + Z_RING_30_60]

cat(sprintf("\nBand sum never exceeds IV06: %s (violations: %d)\n",
            all(chk$BAND_SUM <= chk$IV06, na.rm = TRUE),
            sum(chk$BAND_SUM > chk$IV06, na.rm = TRUE)))
cat(sprintf("Mean band sum %.3f against IV06 mean %.3f (%.1f%% of IV06 within 60 mi)\n",
            mean(chk$BAND_SUM), mean(chk$IV06, na.rm = TRUE),
            100 * mean(chk$BAND_SUM) / mean(chk$IV06, na.rm = TRUE)))
cat(sprintf("Correlation of band sum with IV06: %.4f\n",
            cor(chk$BAND_SUM, chk$IV06, use = "complete.obs")))

stopifnot("band sum exceeds IV06, peer universes do not nest" =
            all(chk$BAND_SUM <= chk$IV06, na.rm = TRUE))

# Three endogenous terms need three instruments with variation, and a band
# that is zero in more than 95% of rows cannot support its own first stage.
# This check only prints a message; the run continues and 16.4 still
# estimates the IV.
share_zero <- ring_desc$SHARE_ZERO
names(share_zero) <- ring_desc$BAND
if (any(share_zero > 0.95)) {
  cat("\n", strrep("!", 70), "\n", sep = "")
  cat("At least one band instrument is zero for more than 95% of rows.\n")
  cat("The three-band IV is not identified. Report the OLS decomposition.\n")
  cat(strrep("!", 70), "\n", sep = "")
  print(share_zero)
}


# -----------------------------------------------------------------------------
# 16.4  Three-band IV and the OLS comparison
# -----------------------------------------------------------------------------
# PRIMARY_OUTCOME on the three ring counts with the baseline controls and
# fixed effects, by OLS and by IV with Z_RING_0_15, Z_RING_15_30, and
# Z_RING_30_60 as instruments, clustered by BASELINE_CLUSTERS. Sample: county
# panel rows with complete outcome, ring counts, band instruments, and
# controls. Writes MD59_ring_decay_iv.csv (OLS and IV rows) and prints the
# first-stage Wald statistics.
#
# This reassigns md_ring_panel, the name of the Part 4.1 ring-treatment
# panel, which is not used after Part 4.
md_ring_panel <- merge(md_county_panel, ring_iv,
                       by = c("HOSPITAL_ID", "POST_MONTH"), all.x = TRUE)
stopifnot("instrument merge changed row count" =
            nrow(md_ring_panel) == nrow(md_county_panel))

ENDO_RING <- paste0("RING_PRIOR_", RING_LABS)
INST_RING <- paste0("Z_RING_",     RING_LABS)
stopifnot("ring treatment columns missing from the panel" =
            all(ENDO_RING %in% names(md_ring_panel)))

md_ring_est <- md_ring_panel[
  complete.cases(md_ring_panel[, c(PRIMARY_OUTCOME, ENDO_RING, INST_RING,
                                   available_columns(md_ring_panel, BASELINE_CONTROLS)),
                               with = FALSE])]
cat(sprintf("\nEstimation sample: %s rows, %d hospitals, %d markets\n",
            format(nrow(md_ring_est), big.mark = ","),
            uniqueN(md_ring_est$HOSPITAL_ID), uniqueN(md_ring_est$ANALYSIS_MARKET)))

# The IV formula is assembled with sprintf() in the fixest form
# y ~ controls | FE | endo1 + endo2 + endo3 ~ z1 + z2 + z3.
ctrl <- available_columns(md_ring_est, BASELINE_CONTROLS)
fe   <- available_columns(md_ring_est, BASELINE_FIXED_EFFECTS)

f_iv <- as.formula(sprintf("%s ~ %s | %s | %s ~ %s",
                           PRIMARY_OUTCOME, paste(ctrl, collapse = " + "),
                           paste(fe, collapse = " + "),
                           paste(ENDO_RING, collapse = " + "),
                           paste(INST_RING, collapse = " + ")))
f_ols <- as.formula(sprintf("%s ~ %s + %s | %s",
                            PRIMARY_OUTCOME, paste(ENDO_RING, collapse = " + "),
                            paste(ctrl, collapse = " + "),
                            paste(fe, collapse = " + ")))

cat("\nIV formula:\n"); print(f_iv)

cl <- build_cluster_formula(available_columns(md_ring_est, BASELINE_CLUSTERS))

fit_ols <- feols(f_ols, data = md_ring_est, cluster = cl, warn = FALSE, notes = FALSE)
fit_iv  <- tryCatch(feols(f_iv, data = md_ring_est, cluster = cl,
                          warn = FALSE, notes = FALSE),
                    error = function(e) { cat("IV failed:", conditionMessage(e), "\n"); NULL })

extract_bands <- function(fit, label, prefix = "") {
  if (is.null(fit)) return(NULL)
  rbindlist(lapply(RING_LABS, function(b) {
    nm <- paste0(prefix, "RING_PRIOR_", b)
    if (!(nm %in% names(coef(fit)))) return(NULL)
    est <- unname(coef(fit)[nm]); se <- unname(sqrt(vcov(fit)[nm, nm]))
    data.table(ESTIMATOR = label, BAND = b, COEF = est, SE = se,
               T_STAT = est / se, P_VALUE = .pval(est / se, fit),
               CI_LO = est - 1.96 * se, CI_HI = est + 1.96 * se,
               N_OBSERVATIONS = nobs(fit))
  }), fill = TRUE)
}

ring_results <- rbindlist(list(
  extract_bands(fit_ols, "OLS"),
  extract_bands(fit_iv,  "IV", prefix = "fit_")), fill = TRUE)

# Fallback for IV coefficient names without the fit_ prefix. It runs only
# when ring_results is empty, so not when the OLS rows were extracted.
if (nrow(ring_results) == 0L)
  ring_results <- extract_bands(fit_iv, "IV")

save_csv(ring_results, "MD59_ring_decay_iv.csv")
cat("\n=== RING DECAY: OLS AND IV ===\n")
print(ring_results[, .(ESTIMATOR, BAND, COEF = round(COEF, 6),
                       SE = round(SE, 6), P_VALUE = round(P_VALUE, 4),
                       N_OBSERVATIONS)])

if (!is.null(fit_iv)) {
  fs <- tryCatch(fitstat(fit_iv, "ivwald"), error = function(e) NULL)
  cat("\n=== FIRST STAGE, ALL THREE BANDS ===\n")
  if (!is.null(fs)) print(fs) else cat("  ivwald unavailable; see summary(fit_iv, stage = 1)\n")
  cat("\nA design is only as identified as its weakest first stage. If the\n")
  cat("0 to 15 mile equation is weak, that coefficient is not interpretable\n")
  cat("even if the far bands are strong.\n")
}


# -----------------------------------------------------------------------------
# 16.5  Figure, OLS only
# -----------------------------------------------------------------------------
# The figure shows the OLS estimates only. The three-band IV system is weakly
# identified: its standard errors are two to five times the OLS ones and the
# sign pattern reverses across bands, as expected with near-collinear
# instruments in a system of three endogenous terms and three instruments.
# In Table 8, OLS is attenuated toward zero at the pooled level (-1.05%
# against an instrumented -3.07%), so a significant near-band OLS
# coefficient is, if anything, a lower bound on the local effect.
#
# Writes fig25_ring_decay.pdf (FIGURE_DIR) and MD59_ring_decay_ols.csv.
# Intervals are 1.96 standard errors, while the filled points mark P_VALUE <
# 0.05 from .pval() (t reference), so the two can disagree near the margin.

ring_ols <- ring_results[ESTIMATOR == "OLS"]
stopifnot("OLS rows missing from ring_results" = nrow(ring_ols) == 3)

band_levels <- c("0 to 15", "15 to 30", "30 to 60")
ring_plot <- copy(ring_ols)
ring_plot[, BAND_LAB := factor(gsub("_", " to ", BAND), levels = band_levels)]
ring_plot[, SIG := P_VALUE < 0.05]

cat("\n=== RING DECAY, OLS (figure data) ===\n")
print(ring_plot[, .(BAND_LAB, COEF, SE, P_VALUE, CI_LO, CI_HI, SIG)])

p_ring <- ggplot(ring_plot, aes(x = BAND_LAB, y = COEF)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey40", linewidth = 0.4) +
  geom_errorbar(aes(ymin = CI_LO, ymax = CI_HI), width = 0.10,
                linewidth = 0.6, colour = "grey25") +
  geom_point(aes(fill = SIG), shape = 21, size = 3.8,
             stroke = 0.7, colour = "grey15") +
  scale_fill_manual(values = c(`TRUE` = "#782F40", `FALSE` = "white"),
                    guide = "none") +
  scale_y_continuous(labels = function(x) sprintf("%.3f", x)) +
  labs(x = "Distance from focal hospital (miles)",
       y = "Change in log price per prior poster in band") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.line          = element_line(colour = "grey20", linewidth = 0.3),
    axis.ticks         = element_line(colour = "grey20", linewidth = 0.3),
    axis.title.x       = element_text(margin = margin(t = 8)),
    axis.title.y       = element_text(margin = margin(r = 8)),
    plot.margin        = margin(6, 10, 6, 6))

ggsave(file.path(FIGURE_DIR, "fig25_ring_decay.pdf"),
       p_ring, width = 6.0, height = 3.6, device = "pdf")
cat("\nFigure written to fig25_ring_decay.pdf\n")

save_csv(ring_ols[, .(BAND, COEF, SE, T_STAT, P_VALUE, CI_LO, CI_HI, N_OBSERVATIONS)],
         "MD59_ring_decay_ols.csv")

rm(pairs, md_ring_panel, fit_iv, ring_iv); invisible(gc())
cat("\n=== PART 16 COMPLETE (OLS reported; IV attempted and rejected on weak-system grounds) ===\n")


# =============================================================================
# Part 17: Out-of-state and ownership-dated constructions
# =============================================================================
#
# Two alternative constructions of the competitor instrument, each aimed at
# one channel through which the exclusion restriction could fail.
#
# 17A, multi-market contracting. A carrier that contracts with a system in
# several markets, or a system that negotiates one contract covering several
# hospitals, could carry disclosure elsewhere into local prices without
# passing through local peers. Both channels are densest within a state.
# CTY_OOS is IV06 with every peer in the focal hospital's own state removed.
#
# 17B, simultaneity. IV06 counts a rival system only once that system's own
# hospital in the focal county has posted, so the instrument's support
# depends on that hospital's timing. CTY_OWN counts the rival system
# throughout the window instead, through presence = "owned" in the Part 10.2
# builder.
#
# Both use build_dynamic_competitor_instrument() at county geography.
# Passing the same-state county pairs as `adjacency` removes in-state peers
# in other counties (CTY_OOS); passing the cross-state pairs removes
# out-of-state peers (CTY_INS). The builder always drops same-county peers,
# so CTY_OOS + CTY_INS must equal the unfiltered count CTY_ALL (checked in
# 17.2).
#
# Needs Parts 0-2, 4, and 10, with md_county_panel and md_event_roster in
# scope. Writes MD_P17A_construction_descriptives.csv,
# MD_P17B_construction_rows.csv, and MD_P17C_construction_tests.csv.


# -----------------------------------------------------------------------------
# 17.1  Clean FIPS roster and the state pair tables
# -----------------------------------------------------------------------------
# Missing FIPS codes stay NA. formatC() turns NA into "   NA", which the
# builder would treat as one county shared by every hospital with a missing
# code; md_fips5() pads only the codes that are present.
md_fips5 <- function(x) {
  v   <- suppressWarnings(as.integer(trimws(as.character(x))))
  out <- rep(NA_character_, length(v))
  out[!is.na(v)] <- formatC(v[!is.na(v)], width = 5, flag = "0")
  out
}

md_roster_st <- copy(md_event_roster)
md_roster_st[, COUNTY_FIPS := md_fips5(COUNTY_FIPS)]
md_roster_st[, STATE_FIPS  := substr(COUNTY_FIPS, 1L, 2L)]
cat(sprintf("Roster %d hospitals | missing county FIPS %d | states %d\n",
            nrow(md_roster_st), sum(is.na(md_roster_st$COUNTY_FIPS)),
            uniqueN(md_roster_st$STATE_FIPS, na.rm = TRUE)))

md_cty_state <- unique(md_roster_st[!is.na(COUNTY_FIPS), .(GEO = COUNTY_FIPS, STATE_FIPS)])
stopifnot("a county maps to two states" = !anyDuplicated(md_cty_state$GEO))

md_pairs <- CJ(GEO_A = md_cty_state$GEO, GEO_B = md_cty_state$GEO)[GEO_A != GEO_B]
md_pairs[md_cty_state, ST_A := i.STATE_FIPS, on = .(GEO_A = GEO)]
md_pairs[md_cty_state, ST_B := i.STATE_FIPS, on = .(GEO_B = GEO)]
md_same_state_pairs  <- md_pairs[ST_A == ST_B, .(GEO_A, GEO_B)]
md_cross_state_pairs <- md_pairs[ST_A != ST_B, .(GEO_A, GEO_B)]
setkey(md_same_state_pairs,  GEO_A, GEO_B)
setkey(md_cross_state_pairs, GEO_A, GEO_B)
cat(sprintf("County pairs: %s same-state | %s cross-state\n",
            format(nrow(md_same_state_pairs),  big.mark = ","),
            format(nrow(md_cross_state_pairs), big.mark = ",")))
rm(md_pairs); invisible(gc())


# -----------------------------------------------------------------------------
# 17.2  Build the four counts and check them
# -----------------------------------------------------------------------------
# CTY_ALL (IV06 rebuilt on the clean FIPS codes), CTY_OOS (out-of-state
# peers), CTY_INS (in-state peers), and CTY_OWN (ownership-dated). Three
# checks stop the run on failure: in-state and out-of-state peers partition
# the unfiltered set; dating presence by ownership can only add peers, never
# remove them; and CTY_ALL reproduces the panel's IV06 in at least 98% of
# hospital-months.
md_p17 <- Reduce(function(a, b) merge(a, b, by = c("HOSPITAL_ID", "POST_MONTH")), list(
  build_dynamic_competitor_instrument(md_county_panel, md_roster_st, "COUNTY_FIPS",
                                      "CTY_ALL", exclude_own_system = TRUE),
  build_dynamic_competitor_instrument(md_county_panel, md_roster_st, "COUNTY_FIPS",
                                      "CTY_OOS", exclude_own_system = TRUE,
                                      adjacency = md_same_state_pairs),
  build_dynamic_competitor_instrument(md_county_panel, md_roster_st, "COUNTY_FIPS",
                                      "CTY_INS", exclude_own_system = TRUE,
                                      adjacency = md_cross_state_pairs),
  build_dynamic_competitor_instrument(md_county_panel, md_roster_st, "COUNTY_FIPS",
                                      "CTY_OWN", exclude_own_system = TRUE,
                                      presence = "owned")))

# Check 1: in-state and out-of-state peers partition the unfiltered set.
stopifnot("in-state + out-of-state != unfiltered count" =
            all(md_p17$CTY_ALL_HOSPITALS ==
                  md_p17$CTY_OOS_HOSPITALS + md_p17$CTY_INS_HOSPITALS))

# Check 2: dating presence by ownership can only add peers.
stopifnot("ownership-dated count below posting-dated count" =
            all(md_p17$CTY_OWN_HOSPITALS >= md_p17$CTY_ALL_HOSPITALS))

# A focal hospital without a county FIPS cannot be placed in a state; its four
# counts are set to NA.
md_p17_cols <- c("CTY_ALL_HOSPITALS", "CTY_OOS_HOSPITALS",
                 "CTY_INS_HOSPITALS", "CTY_OWN_HOSPITALS")
md_p17[HOSPITAL_ID %in% md_roster_st[is.na(COUNTY_FIPS), HOSPITAL_ID],
       (md_p17_cols) := NA_integer_]

# Check 3: the clean-FIPS rebuild reproduces the panel's IV06 (exact match in
# at least 98% of hospital-months; the first mismatches are printed if not).
md_p17_check <- merge(
  md_p17[!is.na(CTY_ALL_HOSPITALS), .(HOSPITAL_ID, POST_MONTH, CTY_ALL_HOSPITALS)],
  unique(md_county_panel[, .(HOSPITAL_ID, POST_MONTH, T_IV06 = get(PRIMARY_INSTRUMENT))]),
  by = c("HOSPITAL_ID", "POST_MONTH"))[!is.na(T_IV06)]
md_p17_match <- mean(md_p17_check$CTY_ALL_HOSPITALS == md_p17_check$T_IV06)
cat(sprintf("\nClean-FIPS rebuild vs panel IV06: exact match %.4f on %s hospital-months\n",
            md_p17_match, format(nrow(md_p17_check), big.mark = ",")))
if (md_p17_match < 0.98) print(head(md_p17_check[CTY_ALL_HOSPITALS != T_IV06], 10))
stopifnot("rebuild does not reproduce IV06 -- do not estimate" = md_p17_match >= 0.98)

md_p17_desc <- md_p17[!is.na(CTY_ALL_HOSPITALS), .(
  HOSPITAL_MONTHS          = .N,
  MEAN_IV06_REBUILD        = mean(CTY_ALL_HOSPITALS),
  MEAN_OUT_OF_STATE        = mean(CTY_OOS_HOSPITALS),
  MEAN_OWNERSHIP_DATED     = mean(CTY_OWN_HOSPITALS),
  SHARE_PEERS_OUT_OF_STATE = sum(CTY_OOS_HOSPITALS) / sum(CTY_ALL_HOSPITALS),
  ZERO_IV06                = mean(CTY_ALL_HOSPITALS == 0L),
  ZERO_OUT_OF_STATE        = mean(CTY_OOS_HOSPITALS == 0L),
  ZERO_OWNERSHIP_DATED     = mean(CTY_OWN_HOSPITALS == 0L),
  CORR_OOS_WITH_IV06       = cor(CTY_OOS_HOSPITALS, CTY_ALL_HOSPITALS),
  CORR_OWN_WITH_IV06       = cor(CTY_OWN_HOSPITALS, CTY_ALL_HOSPITALS))]
save_csv(md_p17_desc, "MD_P17A_construction_descriptives.csv")
print(t(md_p17_desc))


# -----------------------------------------------------------------------------
# 17.3  Headline specification under each construction, common sample
# -----------------------------------------------------------------------------
# The headline interacted model with each of three instruments: IV06
# (PRIMARY_INSTRUMENT), CTY_OOS_HOSPITALS, and CTY_OWN_HOSPITALS. All three
# are estimated on the county panel rows where the out-of-state count is
# defined, so the comparison is not confounded by sample changes. The update
# join adds CTY_OOS_HOSPITALS and CTY_OWN_HOSPITALS to md_county_panel in
# place. Cached as md_p17_constructions_3inst; writes MD_P17B (rows) and
# MD_P17C (tests).
md_county_panel[md_p17, `:=`(CTY_OOS_HOSPITALS = i.CTY_OOS_HOSPITALS,
                             CTY_OWN_HOSPITALS = i.CTY_OWN_HOSPITALS),
                on = .(HOSPITAL_ID, POST_MONTH)]

P17_IVS <- c(Competitor_only_hospitals_9m         = PRIMARY_INSTRUMENT,
             Competitor_out_of_state_hospitals_9m = "CTY_OOS_HOSPITALS",
             Competitor_ownership_dated_9m        = "CTY_OWN_HOSPITALS")

md_p17_res <- cache_or_run("md_p17_constructions_3inst", {
  d <- md_county_panel[!is.na(CTY_OOS_HOSPITALS)]
  out <- lapply(names(P17_IVS), function(nm) {
    cat("\n---", nm, "---\n")
    res <- estimate_interacted(
      data = d, moderator = MD_SCHEME, moderator_type = "categorical",
      outcome = PRIMARY_OUTCOME, instrument = P17_IVS[[nm]],
      label = paste0("P17=", nm), instrument_label = nm,
      moderator_label = MD_SCHEME)
    if (is.null(res)) return(NULL)
    list(rows  = cbind(INSTRUMENT_TESTED = nm, res$rows),
         tests = cbind(INSTRUMENT_TESTED = nm, res$tests))
  })
  rm(d); invisible(gc())
  list(rows  = rbindlist(lapply(out, `[[`, "rows"),  fill = TRUE),
       tests = rbindlist(lapply(out, `[[`, "tests"), fill = TRUE))
})

save_csv(md_p17_res$rows,  "MD_P17B_construction_rows.csv")
save_csv(md_p17_res$tests, "MD_P17C_construction_tests.csv")

cat("\n=== OUT-OF-STATE AND OWNERSHIP-DATED CONSTRUCTIONS (common sample) ===\n")
print(md_p17_res$rows[, .(INSTRUMENT_TESTED, TERM,
                          RF_PCT_PER_SD = round(RF_PERCENT_PER_SD, 2),
                          RF_P = round(RF_P, 4),
                          IV_PERCENT = round(IV_PERCENT, 2),
                          FIRST_STAGE_WALD_MIN = round(FIRST_STAGE_WALD_MIN, 1),
                          N_OBSERVATIONS)])
print(md_p17_res$tests[, .(INSTRUMENT_TESTED, ESTIMATOR, P_VALUE = round(P_VALUE, 4))])


# -----------------------------------------------------------------------------
# 17.4  HSA size distribution in the estimation sample
# -----------------------------------------------------------------------------
# Reports how many sample HSAs contain one, two, three, or four or more
# hospitals. This is the descriptive behind the HSA scope condition used in
# Parts 11.4 and 12. Reloads the HSA panel with md_load_panel().
p_hsa <- md_load_panel("HSA")
hsa_n <- unique(p_hsa[!is.na(ANALYSIS_MARKET), .(HOSPITAL_ID, ANALYSIS_MARKET)])[
  , .(N_HOSP = .N), by = ANALYSIS_MARKET]
hsa_n[, BUCKET := fcase(N_HOSP == 1L, "1", N_HOSP == 2L, "2",
                        N_HOSP == 3L, "3", default = "4+")]
print(hsa_n[, .(SAMPLE_HSAS = .N, SHARE = round(.N / nrow(hsa_n), 3)),
            by = BUCKET][order(BUCKET)])
cat("Sample HSAs:", nrow(hsa_n), "| hospitals:", sum(hsa_n$N_HOSP), "\n")
rm(p_hsa); invisible(gc())

cat("\n=== MARKET DEFINITION ANALYSIS COMPLETE ===\n")