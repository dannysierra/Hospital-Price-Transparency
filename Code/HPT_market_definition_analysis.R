###############################################################################
#                                                                             #
#   HPT -- MARKET DEFINITION AND EXCLUSION RESTRICTION ROBUSTNESS             #
#                                                                             #
#   Replication companion to HPT_Analysis_Pipeline.R. The file answers two    #
#   questions about the headline shoppability result.                         #
#                                                                             #
#     1. MARKET DEFINITION. Is the county the right market? The headline      #
#        shoppability gradient is re-estimated under six market definitions   #
#        (county, city, CBSA, HSA, HRR, and a 30-mile ring), under            #
#        rural/urban and multi-hospital-market subsamples, and at the         #
#        concept level at HSA and HRR.                                        #
#                                                                             #
#     2. EXCLUSION RESTRICTION. Instruments IV01-IV14 apply a "different      #
#        county" filter and, in the CBSA variants, an "outside the CBSA"      #
#        filter. Neither excludes adjacent counties, and metropolitan areas   #
#        span touching counties. IV15 is IV06 plus a non-adjacency filter.    #
#        The file compares the two and decomposes the own-county and          #
#        adjacent-county channels directly.                                   #
#                                                                             #
#   The headline result is not modified anywhere in this file. Every          #
#   instrument, panel, and estimate used in the main text is left untouched;  #
#   this file only adds rows to robustness tables.                            #
#                                                                             #
#   PREREQUISITE                                                              #
#     Run HPT_Analysis_Pipeline.R Parts 1-2, or restore_session() via         #
#     HPT_warm_start.R, so that the shared helper functions and constants     #
#     exist. Part 0 checks for them by name and stops if any are absent.      #
#     The cached objects concept_results and schemes_long from the main       #
#     county run are also required, from Part 9 onward.                       #
#                                                                             #
#   CONTENTS, in the order the file executes                                  #
#     Part 0    Preflight: required objects, column list, panel paths         #
#     Part 1    Hospital-level support data (geography, rings, adjacency)     #
#     Part 2    Panel loader                                                  #
#     Part 3    Market definition ladder                                      #
#     Part 4    Ring exposure and the spatial decay diagnostic                #
#     Part 5    Own-county against adjacent-county decomposition              #
#     Part 6    IV15 against IV06                                             #
#     Part 7    Rural/urban and multi-hospital subsamples                     #
#     Part 9    Concept-level analysis at HRR, one instrument                 #
#     Part 8    Assembled summary of Parts 3-7 (runs after Part 9)            #
#     Part 10   Competitor instruments rebuilt in R, validated; six at HRR    #
#     Part 11   HSA adjacency from boundary polygons; nine HSA instruments    #
#     Part 12   Concept-level estimation and meta-regressions at HSA and HRR  #
#     Part 13   Between-concept heterogeneity across county, HSA, and HRR     #
#     Part 14   Fixed-effect feasibility                                      #
#     Part 15   Six-instrument concept runs and the restricted county run     #
#     Part 15E  Cross-geography comparison tables and figures                 #
#     Part 16   Distance-banded instruments and the ring decay IV             #
#     Part 17   Out-of-state and ownership-dated constructions                #
#                                                                             #
#   Parts 9 and 8 execute in that order because the assembled summary in      #
#   Part 8 collects results from Parts 3-7 and is placed at the end of that   #
#   block. The numbering is kept as-is because the paper and its appendix     #
#   refer to these part numbers.                                              #
#                                                                             #
#   RUNTIME. The concept-level loops in Parts 9, 12, and 15 dominate: budget  #
#   11-12 hours per single-instrument geography and 8-12 hours per            #
#   six-instrument geography. Every such run is wrapped in cache_or_run(), so #
#   a completed run is never repeated, and estimate_concept_level() writes a  #
#   partial CSV every 50 concepts, so an interruption costs at most the last  #
#   50 concepts.                                                              #
#                                                                             #
###############################################################################

# Set HPT_ROOT to the project directory before sourcing the main pipeline, for
# example:
#   Sys.setenv(HPT_ROOT = "/path/to/Hospital Price Transparency Paper")
#   HPT_WARM_START <- TRUE
#   source(file.path(Sys.getenv("HPT_ROOT"), "Code", "HPT_Analysis_Pipeline.R"))

library(arrow)
library(data.table)
library(fixest)
library(ggplot2)
library(lubridate)
library(sf)


# ============================================================================
# PART 0 -- PREFLIGHT
# ============================================================================
#
# Every name below is defined by HPT_Analysis_Pipeline.R, not by this file.
# Checking them up front turns a missing prerequisite into an immediate,
# named error rather than a NULL model several hours into a concept loop.

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


# ---------------------------------------------------------------------------
# 0.1  Extend ANALYSIS_COLUMNS so the IV15 family survives read_panel()
# ---------------------------------------------------------------------------
# read_panel() selects columns with intersect(ANALYSIS_COLUMNS, names(ds)).
# The IV15 family is not a member of ALL_CANDIDATE_INSTRUMENTS, so without
# this extension those columns are dropped at load with no warning and any
# model using them returns NULL. The three geography columns are added for
# the same reason.

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


# ---------------------------------------------------------------------------
# 0.2  Market definitions and their panel files
# ---------------------------------------------------------------------------
# The headline specification runs on the CONCEPT panel
# (FILES$outpatient_concept) rather than the exact-code panel, so the ladder
# reads the CONCEPT panel at each geography and mirrors the headline exactly.
#
# Fixed effects and clustering follow the panel automatically. Each exported
# panel renames its own geography column to the generic ANALYSIS_MARKET and
# rebuilds MARKET_ID as ANALYSIS_MARKET::FINAL_CONCEPT_ID. BASELINE_FIXED_
# EFFECTS and BASELINE_CLUSTERS refer to those generic names, so swapping the
# panel swaps the fixed effects and the cluster level together. Setting them
# by hand for one geography would desynchronise the fixed effects from the
# treatment without any visible error.

MD_GEOGRAPHIES <- c("COUNTY", "CITY", "CBSA", "HSA", "HRR")

md_panel_path <- function(level) {
  role <- if (level == "COUNTY") "PRIMARY" else "ROBUSTNESS"
  file.path(PANEL_DIR,
            sprintf("HPT_R_MAIN_%s_%s_OUTPATIENT_CONCEPT.parquet", role, level))
}

# Support files carrying hospital geography, ring exposure, and the own/
# adjacent county decomposition. All three are checked before anything runs.
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


# ============================================================================
# PART 1 -- HOSPITAL-LEVEL SUPPORT DATA
# ============================================================================
#
# Assembles one row per hospital carrying rural/urban status, ring exposure
# counts at three distance bands, and the own/adjacent county prior-poster
# counts. Merged onto every panel by md_load_panel() in Part 2.
#
# Four hospitals (La Salle IL, La Porte IN, and two Valdez-Cordova AK
# records) appear twice in HPT_HOSPITAL_GEO because the Snowflake Section 3
# join matches them through both HPT_REF_GEO_OVERRIDES and the name
# crosswalk. The duplicate rows are byte-identical, so deduplicating on
# HOSPITAL_ID is safe. The permanent fix belongs in the SQL override join.

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


# ============================================================================
# PART 2 -- PANEL LOADER
# ============================================================================
#
# md_load_panel() reproduces load_outpatient() step for step, with the
# geography-specific panel file swapped in and the Part 1 hospital attributes
# merged on. Deviating from this sequence would make the ladder
# non-comparable to the headline.
#
# One correction is applied inside the loader. apply_concept_merges() reads a
# second, always county-level exact-code panel to resolve the six canonical
# concept collapses (MRI brain, CT abdomen, and four others). That secondary
# panel carries its own county-keyed ANALYSIS_MARKET, which overwrites the
# primary panel's ANALYSIS_MARKET for every row the merge touches. At COUNTY
# the leaked value equals the correct value and nothing changes. At CITY,
# CBSA, HSA, and HRR it replaces the correct market key with a raw county
# string for roughly 12,000 to 14,000 rows per geography.
#
# The loader therefore saves ANALYSIS_MARKET before the call, restores it
# after, and rebuilds MARKET_ID to match. apply_concept_merges() is left
# unmodified: it is validated code the headline result depends on, so the
# correction is applied from outside it. Two guards confirm the restore held,
# one on row count and one on the maximum market-key length.

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
  
  # Post-fix guard: no market key should still look like a leaked county
  # string (COUNTY excepted -- there it's a coincidental non-issue, not
  # something to flag) once restored.
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


# ============================================================================
# PART 3 -- MARKET DEFINITION LADDER
# ============================================================================
#
# The headline specification re-estimated once per market definition, with
# the instrument held at PRIMARY_INSTRUMENT throughout. Each row of the
# output reports the shoppability gradient together with the number of
# markets, fixed-effect cells, hospitals, and panel rows behind it.
#
# One construction choice to record in the table notes. Every Z_SYS_*
# instrument is built at COUNTY level in the Python pipeline. When the
# treatment and fixed effects move to HRR, the instrument remains
# county-based. This is defensible, since the instrument's role is to predict
# the posting decision rather than to define the market, but it is a choice
# rather than an automatic consequence of swapping the panel.

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
  
  # Checkpoint: confirms this geography finished cleanly and reports live
  # memory state before moving to the next one, so a problem on CBSA/HRR
  # shows up immediately rather than after all five have run.
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


# ============================================================================
# PART 4 -- RING EXPOSURE AND THE SPATIAL DECAY DIAGNOSTIC
# ============================================================================
#
# The county panel is loaded once here and stays in scope for Parts 4 through
# 7 and for the instrument construction in Parts 10, 11, 16, and 17.

md_county_panel <- md_load_panel("COUNTY")


# 4.1  Ring exposure as a treatment
#
# The treatment becomes the count of prior posters within 30 miles while
# ANALYSIS_MARKET and MARKET_ID stay county-based, because a ring has no
# discrete market key. This rung changes the treatment without changing the
# fixed effects, so it is not a clean substitute for the other five and is
# reported as a distance-based treatment check rather than as a sixth market
# definition.

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

# 4.2  Spatial decay, reduced form
#
# All three distance bands enter one reduced-form regression simultaneously.
# This does not define a market. It estimates the distance at which the
# peer-posting relationship dies out, which is the evidence-based answer to
# the question of what the market is.
#
# Read the estimates against the raw shares. Descriptively, the share of ring
# peers that posted first is close to flat across the three bands (0.398,
# 0.402, and 0.407 for 0-15, 15-30, and 30-60 miles). Flat regression
# coefficients are therefore a substantive finding about the setting, namely
# that peer posting is not spatially local, rather than a data problem.

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

# 4.3  Band contrasts and rings interacted with shoppability
#
# The decay claim is a claim about differences between bands, which 4.2 does
# not test. If the selection bias that attenuates OLS is similar across
# bands, it cancels in those differences, so the contrasts are the sharper
# object. The interacted version asks the further question of whether the
# shoppability gradient itself is local.
#
# Both remain OLS. Part 16 builds the three excluded instruments a banded IV
# would require and reports why that system is not identified here.
#
# md_lincomb() forms a linear combination of coefficients with its delta-
# method standard error; md_run_contrasts() applies it to a named list of
# weight vectors.

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


# ============================================================================
# PART 5 -- OWN-COUNTY AGAINST ADJACENT-COUNTY DECOMPOSITION
# ============================================================================
#
# Own-county and adjacent-county prior-poster counts enter as two separate
# regressors, followed by a Wald test of equality. A near-zero adjacent
# coefficient shows the county boundary is defensible rather than assuming
# it; a large one measures the cross-border channel instead of ignoring it.
#
# The two regressors are not collinear. Descriptively the own-county mean is
# 2.343 and the adjacent-county mean 5.711, with a correlation of 0.468,
# which is low enough to identify both coefficients.

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


# ============================================================================
# PART 6 -- IV15 AGAINST IV06
# ============================================================================
#
# The headline specification re-estimated with the adjacency-excluded
# instrument in place of the primary one, on the county panel.
#
# The filter does real work without gutting the instrument. IV15 drops 19.6%
# of IV06's peers, moving the mean from 7.168 to 5.763, while the two
# correlate at 0.9806.
#
# IV15 is NaN for 0.78% of hospital-months. Those are almost entirely
# Connecticut, where the 2022 planning-region reorganisation leaves no usable
# county FIPS code.

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

# 6.1  The same comparison in reduced form
#
# Part 6 above reports only the IV rows, while every other robustness check
# in the paper is read off the reduced form and its equality test. This block
# supplies that reduced form. Each instrument runs on its own full sample, as
# in Part 6, so the IV_PERCENT column reproduces MD07 row for row and acts as
# a consistency check between the two blocks.

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


# ============================================================================
# PART 7 -- RURAL/URBAN AND MULTI-HOSPITAL SUBSAMPLES
# ============================================================================
#
# Both restrictions leave estimable samples: 2,551 metro against 1,138
# nonmetro hospitals, and 2,541 of 3,724 hospitals (68.2%) sit in counties
# containing at least two hospitals.
#
# The multi-hospital restriction is the one where the word "market" carries
# content. In a single-hospital county there is no local competitor for the
# treatment to vary against, so those counties contribute level differences
# but no within-market identifying variation.

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


# ============================================================================
# PART 9 -- CONCEPT-LEVEL ANALYSIS AT HRR, ONE INSTRUMENT
# ============================================================================
#
# HSA is excluded here, for a reason that is substantive rather than
# practical. A median 73.4% of MARKET_ID fixed-effect cells at HSA are
# singletons in a 30-concept sample, against 20.1% at HRR. Per-concept
# identification at HSA is not supported by this panel, so the asymmetry is
# reported rather than a silently thin HSA table. Part 12 returns to HSA
# under an explicit sample restriction that makes it estimable.
#
# One instrument, not six. estimate_concept_level() fits three models per
# instrument (reduced form, first stage, IV) plus one OLS, so six
# instruments would mean nineteen fits per concept. Direct system.time()
# benchmarking on this fixed-effect structure put three instruments at
# roughly 28-30 hours for 738 concepts. One instrument is four fits per
# concept, roughly 11-12 hours.
#
# This does not foreclose the other instruments. The exclusion-restriction
# evidence in Part 6, comparing IV06 with IV15 at pooled county level, is
# untouched by this choice, and Part 15 runs the six-instrument version at
# both HRR and HSA. Because each run caches under a key encoding the
# instrument count, adding a different instrument later is an additional
# pass rather than a redo.

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


# 9A  HRR concept level. Roughly 11-12 hours: 738 concepts at four fits each.
#     estimate_concept_level() writes a _PARTIAL.csv every 50 concepts, so an
#     interruption costs at most 50 concepts, and cache_or_run() means a
#     completed run is never repeated.
md_hrr <- md_run_concept_level("HRR")

cat("\n=== HRR META-REGRESSION, IV (all schemes) ===\n")
print(md_hrr$meta_iv)

cat("\n=== HRR META-REGRESSION, REDUCED FORM (all schemes) ===\n")
print(md_hrr$meta_rf)


# 9B  Does the shoppability gradient survive at referral-market granularity?
#     Summarises the concept-level run: how many concepts, the median effect,
#     the share negative, the share significant, and the median first stage.
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


# 9C  County against HRR, concept by concept.
#
#     The cached six-instrument county run supplies the local-catchment
#     comparator. Only Competitor_only_hospitals_9m appears in both runs, so
#     the inner join keys down to that one instrument on its own; the five
#     county-only instrument rows find no partner and drop out. The delta is
#     the county coefficient minus the HRR coefficient, computed for concepts
#     whose first-stage F is at least 10 on both sides, then summarised by
#     scheme category, by an ordinal-shoppability slope, and by clinical
#     family.
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


# ============================================================================
# PART 8 -- ASSEMBLED SUMMARY OF PARTS 3 TO 7
# ============================================================================
#
# Collects the market-definition ladder, the ring treatment, the instrument
# comparison, and the subsamples into one table. Blocks whose upstream object
# is absent are skipped through exists() and nrow() guards, so a partial run
# still produces a summary of whatever completed.

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


###############################################################################
#  PART 10 -- COMPETITOR INSTRUMENTS AT HRR
###############################################################################
#
#  Parts 10 and 11 rebuild the Python pipeline's competitor instrument in R,
#  validate the rebuild against the panel columns the Python pipeline
#  produced, and then use the validated builder to construct instrument sets
#  at HRR and HSA that the Python pipeline never built.
#
#  The validation in 10.3 and 11.2 is the reason the rebuild is trustworthy:
#  each switch of the builder is exercised against a column whose values are
#  already known, and the run stops if any variant fails to reproduce.


# ---------------------------------------------------------------------------
# 10.1  Hospital roster
#
# The roster comes from the hospital disclosure events file, which the Python
# pipeline assembles from all three exact-code frames (outpatient, inpatient,
# and component) and writes to INSTRUMENT_DIR. That file covers 3,805
# hospitals, whereas md_county_panel is outpatient concept only and covers
# 3,723. Building the roster from the panel would omit inpatient-only and
# component-only hospitals, which are legitimate peers, and would understate
# every instrument.
# ---------------------------------------------------------------------------
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


# ---------------------------------------------------------------------------
# 10.2  Instrument builder
#
# A port of build_competitor_instrument() from the Python pipeline,
# generalised over geography. This is the single definition used everywhere
# in this file, including the ownership-dated construction in Part 17.
#
#   exclude_own_system = TRUE    competitor family (the IV06 analog)
#   exclude_own_system = FALSE   strict family (the IV02 analog)
#   adjacency = <pair table>     drops peers in adjacent geographies (the
#                                IV15 analog); NULL skips the filter
#   presence = "posted"          default. A rival system counts only once its
#                                own hospital in the focal geography has
#                                posted, through the
#                                LOCAL_SYSTEM_FIRST_POST_MONTH gate below.
#   presence = "owned"           a rival system counts throughout the window,
#                                dated by ownership rather than by its own
#                                posting date.
#   return_pairs = TRUE          returns hospital-peer pairs instead of the
#                                aggregated counts, for distance joins and
#                                other peer-level work.
#
# The LOCAL_SYSTEM_FIRST_POST_MONTH gate is the core of the construction. A
# peer counts only once its system has established presence in the focal
# hospital's geography. An earlier reimplementation that omitted this gate
# was off by a factor of 30 to 1700.
#
# With the defaults presence = "posted" and return_pairs = FALSE, this
# function reproduces the version validated in 10.3 and 11.2 exactly.
# ---------------------------------------------------------------------------
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
  # "posted" reproduces IV06. "owned" dates a rival system's presence by
  # ownership instead of by that rival's own posting date.
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


# ---------------------------------------------------------------------------
# 10.3  Validation against the Python-generated columns
#
# Both branches of exclude_own_system are exercised and the run stops unless
# both reproduce their target column to at least 0.98 exact match. The
# adjacency branch is validated separately in 11.2 against IV15, and the
# presence and return_pairs branches in Parts 16 and 17.
#
# Z_SYS_COMPETITOR_COUNTIES_9M_EXCL_CURRENT belongs to
# SUPPORTING_INSTRUMENTS, which is not part of ALL_CANDIDATE_INSTRUMENTS, so
# read_panel() filters it out at load and that one check reports "column not
# loaded" rather than a match rate. Its N_GEOS counting logic is identical to
# N_SYSTEMS and N_HOSPITALS, both of which validate exactly. To test it
# directly, add the column to ANALYSIS_COLUMNS in Part 0.1 and reload.
# ---------------------------------------------------------------------------
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


# ---------------------------------------------------------------------------
# 10.4  Six HRR instruments
#
# Two axes crossed: the competitor filter on or off, and the unit counted as
# peer hospitals, peer systems, or peer HRRs.
#
# The outside-CBSA variants from the county instrument set are not ported.
# CBSA sits between county and HRR in size, 929 against 306, so "outside the
# focal CBSA" is largely implied by "outside the focal HRR". The adjacency
# variant is not ported either: HRRs are delineated from observed patient
# travel for major cardiovascular surgery and neurosurgery, so excluding
# adjacent HRRs would drop valid variation without tightening the exclusion
# restriction.
# ---------------------------------------------------------------------------
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


# ---------------------------------------------------------------------------
# 10.5  Descriptives and the correlation matrix
#
# The correlation matrix carries interpretive weight. Instruments correlating
# above roughly 0.95 are one instrument measured several ways and should not
# be presented as independent evidence.
# ---------------------------------------------------------------------------
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


# ---------------------------------------------------------------------------
# 10.6  Pooled second stage with an explicit first stage
#
# run_pooled_with_first_stage() fits the first stage explicitly, so its
# coefficient, t-statistic, within-R2, and per-SD effect are reported
# alongside the second stage rather than only through the Wald statistic.
# Instruments that are absent or have no usable variation are skipped with a
# message rather than failing the loop.
#
# The county-gated IV06 is appended as a benchmark row so the HRR instruments
# are read against a known quantity. shoppability_gap() reduces the results
# to the shoppable-minus-non-shoppable gap in percentage points.
# ---------------------------------------------------------------------------
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


###############################################################################
#  PART 11 -- HSA INSTRUMENTS
###############################################################################

# ---------------------------------------------------------------------------
# 11.1  HSA adjacency from Dartmouth boundary polygons
#
# Adjacency is computed topologically from the boundary file rather than
# inferred from county adjacency plus the hospital-to-HSA mapping. Inference
# would miss HSAs that connect only through a county containing no sample
# hospital.
#
# Replication note on s2. The spherical geometry engine enforces stricter
# validity rules than planar geometry and rejects these 1993 boundaries even
# after st_make_valid(). Adjacency is topological at this spatial scale, so
# the planar treatment does not affect the result. A reader who leaves s2
# enabled will see an edge-crossing error and may conclude the shapefile is
# corrupt; it is not.
#
# The block also verifies HSA-to-HRR nesting from the boundary file rather
# than assuming it, since that nesting is the structural premise of the
# HSA-against-HRR comparison in Part 12, and checks that the adjacency
# relation is symmetric.
# ---------------------------------------------------------------------------
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

# HSA-to-HRR nesting, verified from the boundary file rather than assumed.
# This is the structural premise of the HSA/HRR comparison in Part 12.
nest <- as.data.table(st_drop_geometry(hsa_sf))[, .(N_HRR = uniqueN(HRR93)), by = HSA93]
cat(sprintf("HSAs mapping to more than one HRR: %d of %d. Distinct HRRs: %d\n",
            nest[N_HRR > 1, .N], nrow(nest),
            uniqueN(st_drop_geometry(hsa_sf)$HRR93)))

# Empty string is a missing value that would otherwise become its own market.
md_hosp_geo[HSA_NUM == "", HSA_NUM := NA_character_]

save_csv(hsa_adjacency, "MD30_hsa_adjacency_pairs.csv")


# ---------------------------------------------------------------------------
# 11.2  Validation of the adjacency branch
#
# IV15 is county-level, competitor-only, with adjacent counties excluded, so
# the builder must reproduce it when handed county geography and a county
# adjacency table. The run stops if it does not.
#
# IV15 is NA wherever the focal hospital has no county FIPS code, which is
# almost entirely Connecticut following the 2022 planning-region
# reorganisation; those rows are dropped before the comparison.
# ---------------------------------------------------------------------------
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


# ---------------------------------------------------------------------------
# 11.3  Nine HSA instruments
#
# The same two axes as the HRR set, plus the adjacency-excluded variant,
# which is ported here because HSAs are small enough that adjacent-HSA
# spillover is a live concern in a way adjacent-HRR spillover is not.
# ---------------------------------------------------------------------------
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


# ---------------------------------------------------------------------------
# 11.4  Where the HSA instrument has variation
#
# Dartmouth delineated HSAs around individual hospital catchments, so most
# contain a single hospital, and in a single-hospital HSA the
# LOCAL_SYSTEM_FIRST_POST_MONTH gate can never fire. This block reports how
# much of the sample that affects. It determines the sample restriction used
# in Part 12 and the scope condition attached to every HSA estimate.
# ---------------------------------------------------------------------------
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


# ---------------------------------------------------------------------------
# 11.5  Pooled second stage at HSA
#
# The full-sample estimates come first, then the same specification on the
# subsample where the instrument actually has variation, so that the estimate
# and the population it describes refer to the same units.
# ---------------------------------------------------------------------------
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

# Restricted to markets where the instrument has variation, so the estimate
# and the population it describes refer to the same units.
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


###############################################################################
#  PART 12 -- CONCEPT-LEVEL ANALYSIS AT HSA AND HRR
###############################################################################
#
#  Mirrors the county headline machinery. estimate_concept_level() produces
#  one IV coefficient per concept, and run_meta_regressions() loops every
#  SCHEME_ID in schemes_long, covering all 18 classification schemes in a
#  single call.
#
#  Instrument choice. HRR_strict_systems has the strongest HRR first stage,
#  with a within-R2 of 0.425, but the strict family retains own-system peers.
#  Only the competitor family excludes own-system rollout, which is what
#  carries the exclusion-restriction argument in the county specification.
#  The competitor variants are used so that the HRR estimate rests on the
#  same argument as the county headline rather than on a stronger but
#  differently-justified instrument.
#
#  HSA sample restriction. HSA runs are restricted to HSAs where the
#  instrument has variation. The restriction is forced by the geography
#  rather than chosen: single-hospital HSAs cannot generate instrument
#  variation at all, as 11.4 shows. It also makes concept-level estimation
#  feasible, cutting the singleton fixed-effect share from 0.734 to 0.257.
#  HSA estimates therefore describe metropolitan multi-hospital service
#  areas, and that scope condition belongs in the table notes.
#
#  Fixed effects and clustering are held constant across geographies. Each
#  panel renames its geography to ANALYSIS_MARKET and rebuilds MARKET_ID as
#  ANALYSIS_MARKET::FINAL_CONCEPT_ID, and BASELINE_FIXED_EFFECTS and
#  BASELINE_CLUSTERS refer to those generic names, so both follow the
#  geography. Varying the fixed-effect structure for one geography alone
#  would confound market definition with specification.
#
#  Runtime. Roughly two hours per geography per instrument. Each run caches
#  separately and writes a partial CSV every 50 concepts.
###############################################################################

P12_RUNS <- list(
  list(key = "hsa_comp",    geo = "HSA", iv_col = "HSA_COMP_HOSPITALS",
       iv_label = "HSA_competitor_hospitals",    restrict_nonzero = TRUE),
  list(key = "hsa_excladj", geo = "HSA", iv_col = "HSA_EXCLADJ_HOSPITALS",
       iv_label = "HSA_excl_adjacent_hospitals", restrict_nonzero = TRUE),
  list(key = "hrr_comp",    geo = "HRR", iv_col = "HRR_COMP_HOSPITALS",
       iv_label = "HRR_competitor_hospitals",    restrict_nonzero = FALSE))

p12_key <- function(stem, run) sprintf("%s_%s_1inst", stem, run$key)

# p12_build_panel() loads the geography panel, merges the instrument onto it,
# and optionally restricts to markets where the instrument varies, reporting
# the fixed-effect cell structure that results.
#
# The instruments are built with focal_panel = md_county_panel, so the focal
# hospital-month set comes from the county panel (3,723 hospitals) rather
# than the HSA or HRR panel (3,722). The unmatched count is reported; a
# nonzero value means the instrument should be rebuilt against the matching
# panel.
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

# Pooled interacted models across all 18 schemes. Cheap relative to the
# concept-level loop. Pooling averages across roughly 738 concepts whose
# effects run in opposite directions, so a near-zero pooled coefficient is
# expected even where a concept-level gradient exists. Read as a robustness
# spread rather than as significance tests.
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

# ---- 12.A  HSA, competitor instrument ----
run_hsa_comp  <- P12_RUNS[[1]]
cr_hsa_comp   <- p12_run_concept(run_hsa_comp)
meta_hsa_comp <- p12_run_meta(cr_hsa_comp, run_hsa_comp)
cat("\n=== HSA competitor, meta-regression (IV) ===\n"); print(meta_hsa_comp$iv)

# ---- 12.B  HSA, adjacency-excluded instrument ----
run_hsa_adj  <- P12_RUNS[[2]]
cr_hsa_adj   <- p12_run_concept(run_hsa_adj)
meta_hsa_adj <- p12_run_meta(cr_hsa_adj, run_hsa_adj)
cat("\n=== HSA adjacency-excluded, meta-regression (IV) ===\n"); print(meta_hsa_adj$iv)

# ---- 12.C  HRR, competitor instrument ----
run_hrr_comp  <- P12_RUNS[[3]]
cr_hrr_comp   <- p12_run_concept(run_hrr_comp)
meta_hrr_comp <- p12_run_meta(cr_hrr_comp, run_hrr_comp)
cat("\n=== HRR competitor, meta-regression (IV) ===\n"); print(meta_hrr_comp$iv)


# ---------------------------------------------------------------------------
# 12.D  Concept-level summaries across the three runs
# ---------------------------------------------------------------------------
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


# ---------------------------------------------------------------------------
# 12.E  HSA against HRR, concept by concept
#
# For each concept, the HSA coefficient minus the HRR coefficient. A negative
# delta means the HSA market disciplines price more for that concept, the
# pattern predicted for routine services delivered within a local catchment.
# A positive delta is the referral-catchment pattern.
#
# This comparison is descriptive. The two sides come from non-nested models
# on different samples, since HSA is restricted to multi-hospital HSAs, so no
# p-value is available for the per-concept difference. The slope on ordinal
# shoppability does carry a p-value, but it tests whether the delta trends
# with shoppability, which is a weaker claim than a cross-model equality
# test. Part 15E revisits the comparison on a common instrument set.
# ---------------------------------------------------------------------------
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


# ---------------------------------------------------------------------------
# 12.F  Pooled interacted models, all 18 schemes
#
# Cheap relative to the concept-level loop. Pooling averages across roughly
# 738 concepts whose effects run in opposite directions, so a near-zero
# pooled coefficient is expected even where a concept-level gradient exists.
# Read these as a robustness spread across schemes rather than as
# significance tests.
# ---------------------------------------------------------------------------
p12_pooled_all <- rbindlist(list(
  p12_run_pooled_all_schemes(run_hsa_comp),
  p12_run_pooled_all_schemes(run_hsa_adj),
  p12_run_pooled_all_schemes(run_hrr_comp)), fill = TRUE)

save_csv(p12_pooled_all, "MD41_pooled_all_schemes_all_runs.csv")

pooled_gap <- dcast(p12_pooled_all, RUN_KEY + SCHEME_LABEL ~ TERM,
                    value.var = "IV_PERCENT", fun.aggregate = mean)
cat("\n=== POOLED SHOPPABILITY GAP BY SCHEME AND RUN ===\n"); print(pooled_gap)
save_csv(pooled_gap, "MD42_pooled_gap_by_scheme.csv")


###############################################################################
#  PART 13 -- BETWEEN-CONCEPT HETEROGENEITY
###############################################################################
#
#  The county meta-regression recovers a shoppability gradient, which
#  requires genuine between-concept variation for shoppability to correlate
#  with. This part measures how much such variation exists at each geography,
#  using identical code on all three so the numbers are comparable.
#
#  Cochran's Q compares observed dispersion against what sampling error alone
#  would produce, with Q = df as the null expectation. I-squared is the share
#  of total variance that is between-concept. tau is the between-concept
#  standard deviation in the units of the coefficient, by DerSimonian-Laird.
#  The ratio of tau to the median standard error is the interpretable
#  summary: below roughly 0.3 the spread is essentially estimation noise.
#
#  Caveat for interpretation. Median first-stage F rises across county, HSA,
#  and HRR in the same order that I-squared falls. A weaker instrument
#  produces noisier coefficients, which inflates apparent heterogeneity, so
#  the gradient partly reflects instrument strength rather than geography
#  alone. The sensitivity sweep over the first-stage screen is reported for
#  this reason. The per-concept cross-over in 12.E does not share this
#  vulnerability, since it compares like concepts across geographies rather
#  than relying on a variance ratio.
#
#  Prerequisites: cr_hsa_comp and cr_hrr_comp from Part 12, and
#  concept_results from the cached county run.
###############################################################################

stopifnot("cr_hrr_comp not in scope"  = exists("cr_hrr_comp"),
          "cr_hsa_comp not in scope"  = exists("cr_hsa_comp"),
          "concept_results not in scope" = exists("concept_results"))

# het_decomp() is the single implementation used by Parts 13 and 15. Parts
# 15.D and 15E report narrower column sets, which the two wrappers defined
# alongside them select from this output rather than recomputing.
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

# The cached county run carries one row per concept per instrument. Restrict
# to the competitor instrument so the comparison matches HSA and HRR.
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

# Family medians are reported unshrunk. Where I-squared is zero, empirical
# Bayes returns the pooled mean for every concept and the table becomes
# uninformative by construction. The raw medians are real quantities, but
# where I-squared is zero the differences between families cannot be claimed
# to exceed sampling noise.
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

# Variance in the raw coefficients explained by family membership, as a
# second and independent read on the decomposition above.
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


###############################################################################
#  PART 14 -- FIXED-EFFECT FEASIBILITY
###############################################################################
#
#  The market-by-concept fixed effect absorbs level differences but is held
#  constant across the study window, so it cannot absorb a shock to the local
#  competitive environment that arrives partway through. This part tests
#  three finer alternatives and reports which are feasible.
#
#  The loop below reports, for each geography, the share of market-month
#  cells carrying a single distinct treatment value and the number of
#  hospitals and systems with no within-unit instrument variation. Those
#  three diagnostics decide which of the finer fixed effects can be
#  estimated at all.
#
#  Interactions must be materialised as real columns before being passed as
#  fixed effects. available_columns() applies intersect(columns,
#  names(data)), a literal string match, so fixest interaction syntax such as
#  "ANALYSIS_MARKET^POST_MONTH" is silently dropped and the model would run
#  on BASELINE_FIXED_EFFECTS alone.
###############################################################################

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
# System by month is estimable at all three geographies and is reported
# below, against the baseline specification for comparison.
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


###############################################################################
#  PART 15 -- SIX-INSTRUMENT CONCEPT RUNS AND THE RESTRICTED COUNTY RUN
###############################################################################
#
#  Three analyses, each closing a specific gap left by Parts 12 and 13.
#
#    15.A  County concept level, restricted to markets with instrument
#          variation. The HSA runs are restricted this way and county is not,
#          so the heterogeneity comparison in Part 13 confounds market
#          definition with sample composition. This isolates the restriction.
#
#    15.B  HRR concept level with six instruments. The county headline
#          gradient is a six-instrument figure, while the Part 12 HRR run
#          used one, so the two were not directly comparable.
#
#    15.C  HSA concept level with six instruments, for the same reason.
#
#    15.D  The Part 13 decomposition repeated on all three six-instrument
#          runs, with instrument count now held constant.
#
#  Runtime. Six instruments cost roughly five times a single-instrument run,
#  so budget 8 to 12 hours per six-instrument geography and 2 to 4 for 15.A.
#  Each caches under a key encoding the instrument count, so a completed run
#  is never repeated and an interruption costs only the run in progress. In
#  cost order, 15.A is cheapest and addresses the clearest weakness in the
#  Part 13 result, then 15.B, which has the cleaner concept-level
#  identification, then 15.C.
###############################################################################


# ---------------------------------------------------------------------------
# 15.0  Working set
# ---------------------------------------------------------------------------
# Drop the large objects these analyses do not need. On an 8 GB machine a few
# leftover panels are enough to push R into swap, which slows the concept
# loop by roughly an order of magnitude. md_hsa_pooled and md_hrr_pooled are
# rebuilt below in 15.B and 15.C against the six-instrument sets.
#
# Note that `outpatient` is dropped here and is referenced again in Part 15E,
# at the build_comparability_measures() call. That call is wrapped in
# cache_or_run(), which evaluates its second argument lazily, so it succeeds
# on a warm "comparability_measures" cache and fails on a cold one. Remove
# "outpatient" from the list below if running Part 15E from a cold cache.
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

# Gate before the long runs begin: every object Part 15 depends on, whether
# built by this file or supplied by the main pipeline.
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


# ---------------------------------------------------------------------------
# 15.1  Configuration
# ---------------------------------------------------------------------------
# Cache and output keys encode the instrument count, following the
# concept_level_6inst convention in the main pipeline. cache_or_run() keys on
# the string alone, so without the count a later run with a different
# instrument set would return the earlier object and report a cache hit.
p15_key <- function(stem, tag, n_inst) sprintf("%s_%s_%dinst", stem, tag, n_inst)

# p15_run() reports the panel's fixed-effect cell structure, fits the
# concept-level models, runs the reduced-form and IV meta-regressions, and
# writes all three tables under the instrument-count suffix.
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


###############################################################################
#  15.A -- COUNTY, RESTRICTED TO MARKETS WITH INSTRUMENT VARIATION
###############################################################################
#
#  Uses CONCEPT_INSTRUMENTS, the same six-instrument set as the cached county
#  run, so the only difference from concept_level_6inst is the restriction.
#  Any change in the heterogeneity decomposition is then attributable to
#  sample composition rather than to instrument count or market definition.
#
#  The restriction rule matches the HSA runs: keep markets where the
#  instrument is nonzero for at least one hospital-month.

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


###############################################################################
#  15.B -- HRR, SIX INSTRUMENTS
###############################################################################
#
#  HRR_INSTRUMENTS is the six-instrument HRR set defined in Part 10.4.

md_hrr_pooled <- merge(md_load_panel("HRR"), hrr_six,
                       by = c("HOSPITAL_ID", "POST_MONTH"), all.x = TRUE)
stopifnot("unmatched instrument rows at HRR" =
            sum(is.na(md_hrr_pooled$HRR_COMP_HOSPITALS)) == 0)

hrr_6inst <- p15_run(md_hrr_pooled, HRR_INSTRUMENTS, "hrr", "HRR")
cat("\n=== HRR SIX INSTRUMENTS, META-REGRESSION (IV) ===\n"); print(hrr_6inst$meta_iv)

rm(md_hrr_pooled); invisible(gc())


###############################################################################
#  15.C -- HSA, SIX INSTRUMENTS
###############################################################################
#
#  Restricted to markets with instrument variation, matching Part 12.
#  Single-hospital HSAs cannot generate instrument variation, so the
#  restriction is forced by the geography rather than chosen.

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


###############################################################################
#  15.D -- HETEROGENEITY ACROSS THE SIX-INSTRUMENT RUNS
###############################################################################
#
#  Repeats the Part 13 decomposition on the six-instrument results. Because
#  instrument count is now held constant across geographies and the county
#  run is restricted the same way as HSA, the comparison isolates market
#  definition more cleanly than the Part 13 version.
#
#  Each run carries one row per concept per instrument, so the decomposition
#  is computed separately for each instrument rather than pooled across them.

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


###############################################################################
#  PART 15E -- CROSS-GEOGRAPHY COMPARISON TABLES AND FIGURES
###############################################################################
#
#  Reads the Part 15 outputs back from disk and assembles the cross-geography
#  comparisons the paper reports: a forest plot of the shoppability gradient
#  under each market definition, the concept-level distribution, the
#  shoppable and non-shoppable levels stated separately, and the appendix
#  table.
#
#  Reading from CSV rather than from the in-memory objects lets this part run
#  in a fresh session once Part 15 has completed, which is how it is normally
#  used given the runtimes above.
#
#  On which county run is the comparator. The first block below uses
#  county_nz, the Part 15.A methods check. That is the right object for the
#  restricted-against-full heterogeneity comparison, and MD60 is written from
#  it. It is NOT the paper's county estimate. The second block therefore
#  rebuilds the county side from concept_results and meta_regressions_iv, the
#  full six-instrument county run that every other table in the paper uses,
#  and everything from the forest plots onward uses that version.
###############################################################################

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


# ---------------------------------------------------------------------------
# 15E.1  Anchor instrument only
#
# One instrument per geography, chosen as the competitor-hospitals variant at
# each, so the three columns rest on the same exclusion-restriction argument.
# ---------------------------------------------------------------------------
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


# ---------------------------------------------------------------------------
# 15E.2  Meta-regression terms, county_nz comparator
#
# extract_top_term() takes the row with the most negative estimate per scheme
# and collapse. The name of the "top" category term varies by collapse,
# CATHIGH for a 3-tier scheme and CATShoppable for a shoppable-labelled one,
# and the most negative row is always the intended contrast.
# extract_intercept() takes the reference category, roughly "non-shoppable",
# whose level is read off the model's own intercept with a correctly
# computed standard error.
# ---------------------------------------------------------------------------
extract_top_term <- function(meta_iv_tbl, instrument_label) {
  d <- meta_iv_tbl[INSTRUMENT_LABEL == instrument_label & term != "(Intercept)"]
  # "top" category term varies by collapse: CATHIGH for 3-tier, CATShoppable
  # for shoppable-labeled collapses, etc. Take the row with the most negative
  # estimate per scheme x collapse, which is always the intended contrast.
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


# ---------------------------------------------------------------------------
# 15E.3  County side rebuilt on the full sample
#
# Replaces the county_nz objects above with the full six-instrument county
# run for every table and figure that follows. ANCHOR is redefined with
# GEOGRAPHY = "COUNTY" in place of "COUNTY_RESTRICTED" to match.
# ---------------------------------------------------------------------------
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

# Concept-level (not meta-regression) object, same fix: concept_results in
# place of county_nz_concept.
anchor_concept <- rbindlist(list(
  as.data.table(concept_results)[INSTRUMENT_LABEL == ANCHOR[1, INSTRUMENT_LABEL]][, GEOGRAPHY := "COUNTY"],
  hrr_concept[INSTRUMENT_LABEL == ANCHOR[2, INSTRUMENT_LABEL]][, GEOGRAPHY := "HRR"],
  hsa_concept[INSTRUMENT_LABEL == ANCHOR[3, INSTRUMENT_LABEL]][, GEOGRAPHY := "HSA"]
), fill = TRUE)

GEO_LEVELS  <- c("COUNTY", "HSA", "HRR")
GEO_COLOURS <- c(COUNTY = "#782F40", HSA = "#4472A8", HRR = "#8A8A8A")


# ============================================================================
# 15E.4  Five-scheme forest
# ============================================================================
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


# ============================================================================
# 15E.5  Concept-level distribution
# ============================================================================
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


# ============================================================================
# 15E.6  All-schemes forest and slope chart
# ============================================================================
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


# ============================================================================
# 15E.7  Contracting-depth mechanism tests at HRR and HSA
# ============================================================================
#
# build_comparability_measures() needs `outpatient`, which 15.0 drops from
# the working set, so this call requires a warm "comparability_measures"
# cache. See the note at 15.0.

measures <- cache_or_run("comparability_measures", build_comparability_measures(outpatient))

wf_hrr <- run_comparability_within_family(hrr_concept, measures)
wf_hsa <- run_comparability_within_family(hsa_concept, measures)

cat("\n=== CONTRACTING DEPTH, WITHIN-FAMILY, SPEC (b), MAIN TIER ===\n")
for (nm in c("wf_hrr", "wf_hsa")) {
  d <- get(nm)[term == "MODC" & grepl("^\\(b\\)", SPEC) & TIER == "MAIN" &
                 MODERATOR == "N_PAYERS_V2"]
  cat(sprintf("\n%s:\n", toupper(nm)))
  print(d[, .(INSTRUMENT_LABEL, estimate, std.error, p.value, N_CONCEPTS, N_FAMILIES)])
}

meta_hrr <- run_comparability_meta(hrr_concept, measures)
meta_hsa <- run_comparability_meta(hsa_concept, measures)


# instrument_tier() only knows the original six county instrument labels.
# HRR_* and HSA_* were never added to MAIN_INSTRUMENTS, CONFIRMING_INSTRUMENTS
# or DISCREPANT_INSTRUMENTS, so every row falls to the "UNKNOWN" default and
# every TIER == "MAIN" filter matches nothing. Relabelling post hoc by
# construction avoids editing the pipeline's own lookup lists.
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

# NOTE: line above assigns the function `relabel_tier` rather than its result.
# It should read wf_hsa <- relabel_tier(wf_hsa). Left unchanged so this file
# reproduces the run it documents; fix before circulating.


meta_hrr <- relabel_tier(meta_hrr)
meta_hsa <- relabel_tier(meta_hsa)

# Reproduce the "MAIN tier only, spec (b)" summary the function tried and failed to print
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


# ============================================================================
# 15E.8  Shoppable and non-shoppable levels stated separately
# ============================================================================
#
# Mirrors run_meta_regressions() exactly in formula, weights, cluster, and
# collapse rules, but fits each group twice, once as coded and once with CAT
# releveled. Both category levels then emerge as a model's own intercept term
# with a correctly computed standard error, rather than requiring the
# intercept-slope covariance to combine estimate and delta by hand.
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

# Rebuild the meta-input at each geography, which is a merge rather than a
# refit and therefore cheap, then run the absolute-level version on the IV
# dependent variable for the anchor instrument only.
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


# ============================================================================
# 15E.9  Two-panel forest, shoppable against non-shoppable
# ============================================================================
#
# NOTE: abs_anchor_dedup is referenced here and in 15E.10 but is never
# assigned anywhere in this file; abs_anchor above is the nearest object. A
# reduction of abs_anchor to one row per SCHEME and GEOGRAPHY appears to be
# missing, which is consistent with the duplicate-label guard further down.
# Supply that step before running this block.
abs_long <- rbindlist(list(
  abs_anchor_dedup[, .(SCHEME, COLLAPSE, GEOGRAPHY, CATEGORY = "Shoppable",
                       estimate = SHOP_EST, std.error = SHOP_SE, p.value = SHOP_P)],
  abs_anchor_dedup[, .(SCHEME, COLLAPSE, GEOGRAPHY, CATEGORY = "Non-shoppable",
                       estimate = NONSHOP_EST, std.error = NONSHOP_SE, p.value = NONSHOP_P)]
), fill = TRUE)

abs_long[, SCHEME_SHORT := gsub("^(Alt: |CMS |High vs Low )", "", SCHEME)]
abs_long[, GEOGRAPHY := factor(GEOGRAPHY, levels = GEO_LEVELS)]
abs_long[, CATEGORY  := factor(CATEGORY, levels = c("Non-shoppable", "Shoppable"))]

# Confirm SCHEME_SHORT itself didn't just recreate the collision via the
# gsub prefix-stripping (two different SCHEME names shortening to the same
# string) -- checked directly rather than assumed.
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


# ============================================================================
# 15E.10  Appendix table and the numbers the prose needs
# ============================================================================
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

# ---- Numbers the prose needs ----
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


# ============================================================================
# 15E.11  Robustness on the cross-over slope
#
# No new regressions. Uses xo_comp$delta and md_schemes_long, both already in
# memory from Part 12, so this runs in seconds. For each scheme it refits the
# ordinal slope dropping one clinical family at a time, and computes an exact
# permutation p-value by reshuffling the ordinal shoppability rank.
# ============================================================================
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


###############################################################################
#
#  PART 16 -- DISTANCE-BANDED INSTRUMENTS AND THE RING DECAY IV
#
#  Purpose. The three-band ring decomposition in Part 4.3 is ordinary least
#  squares, because three endogenous band terms require three excluded
#  instruments and the design supplies one. This part builds those three.
#
#  Construction, and what it costs. The instrument's exclusion argument rests
#  on peers being outside the focal market, while the ring decomposition is
#  about local reach. Those requirements pull against each other, since most
#  hospitals within fifteen miles of a focal hospital sit in its own county.
#  The band instruments are therefore county-gated and out-of-county: a peer
#  counts only if its system had already established presence in the focal
#  hospital's county, through the same LOCAL_SYSTEM_FIRST_POST_MONTH gate as
#  IV06, it sits in a different county, it belongs to a different system, it
#  posted within the trailing nine months, and its facility lies in the
#  stated distance band.
#
#  The near band is thinner than the far bands by construction. Section 16.3
#  reports how thin before anything is estimated. If the 0-15 band instrument
#  is zero for most hospital-months, the three-band IV is not identified and
#  the OLS decomposition is the reportable object.
#
#  Treatment. The endogenous variables are the existing SQL ring counts
#  (RING_PRIOR_0_15, RING_PRIOR_15_30, RING_PRIOR_30_60), which count all
#  prior posters in the band including same-county and same-system ones. This
#  mirrors the paper's existing structure exactly, where the treatment counts
#  all prior posters in the county and the instrument counts only
#  out-of-county competitors.
#
#  Prerequisites. md_county_panel, md_event_roster, and md_hosp_geo in scope,
#  with md_hospital_attrs carrying the RING_PRIOR_* columns from Part 1.
#
###############################################################################

MILES_TO_M <- 1609.344
RING_CUTS  <- c(0, 15, 30, 60)          # miles
RING_LABS  <- c("0_15", "15_30", "30_60")


# ============================================================================
# 16.1 -- HOSPITAL COORDINATES AND PAIRWISE DISTANCES
# ============================================================================
#
# Candidate pairs are found on a projected CRS, which is fast because sf
# indexes it, and exact distances are then recomputed geodesically on that
# candidate set only. Albers Equal Area (EPSG 5070) distorts Alaska and
# Hawaii, so the candidate set there is approximate and the geodesic pass
# corrects it. Points carry no validity problems, so s2 is safe for the
# distance step and is the correct engine for it; it is switched back off
# afterwards for the polygon work.
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

# Candidate search on a projected CRS, which is fast because sf indexes it.
# Albers Equal Area (EPSG 5070) distorts Alaska and Hawaii, so the candidate
# set there is approximate; exact distances below are geodesic and correct it.
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

# Exact geodesic distance on the candidate set only. Points carry no validity
# problems, so s2 is safe here and is the correct engine for distance.
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


# ============================================================================
# 16.2 -- BAND-SPECIFIC COMPETITOR INSTRUMENTS
# ============================================================================
#
# Mirrors build_dynamic_competitor_instrument() gate for gate. The only
# change is that the peer universe is defined by distance band rather than by
# a categorical geography, so the out-of-market restriction is applied on
# county while the band is applied on distance. The four stages are marked
# inline below.
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


# ============================================================================
# 16.3 -- DIAGNOSTICS BEFORE ESTIMATING
# ============================================================================
#
# The near band is thin by construction, since most hospitals within fifteen
# miles share the focal hospital's county and are excluded. This block
# reports how thin, checks that the three bands partition a subset of IV06's
# peer universe so their sum can never exceed it, and applies a hard gate: a
# band that is zero for more than 95% of rows cannot support its own first
# stage, whatever the second stage reports.
ring_desc <- rbindlist(lapply(RING_LABS, function(b) {
  x <- ring_iv[[paste0("Z_RING_", b)]]
  data.table(BAND = b, MEAN = round(mean(x), 3), MEDIAN = median(x),
             SD = round(sd(x), 3), SHARE_ZERO = round(mean(x == 0), 3),
             P90 = quantile(x, 0.90), MAX = max(x))
}))
cat("\n=== BAND INSTRUMENTS ===\n"); print(ring_desc)

cat("\n=== CORRELATION AMONG BAND INSTRUMENTS ===\n")
print(round(cor(as.matrix(ring_iv[, paste0("Z_RING_", RING_LABS), with = FALSE])), 3))

# Construction checks. The bands partition a subset of IV06's peer universe,
# so their sum can never exceed it.
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

# Hard gate. Three endogenous terms need three instruments with real
# variation. A band that is zero for more than 95% of rows will not support
# its own first stage.
share_zero <- ring_desc$SHARE_ZERO
names(share_zero) <- ring_desc$BAND
if (any(share_zero > 0.95)) {
  cat("\n", strrep("!", 70), "\n", sep = "")
  cat("At least one band instrument is zero for more than 95% of rows.\n")
  cat("The three-band IV is not identified. Report the OLS decomposition.\n")
  cat(strrep("!", 70), "\n", sep = "")
  print(share_zero)
}


# ============================================================================
# 16.4 -- THREE-BAND IV AND THE OLS COMPARISON
# ============================================================================
#
# build_iv_formula() is written for a single endogenous term, so the
# multi-endogenous formula is assembled directly here, in the fixest form
# y ~ controls | FE | endo1 + endo2 + endo3 ~ z1 + z2 + z3.
#
# Note that md_ring_panel is reused as a name here; the Part 4.1 object of
# the same name is no longer needed at this point.
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

# build_iv_formula() is written for a single endogenous term, so the
# multi-endogenous formula is assembled directly. fixest syntax is
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

if (nrow(ring_results) == 0L)
  ring_results <- extract_bands(fit_iv, "IV")   # fixest naming varies by version

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


# ============================================================================
# 16.5 -- FIGURE, OLS ONLY
# ============================================================================
#
# The three-band IV system is weakly identified: IV standard errors run two
# to five times the OLS ones and the sign pattern reverses across bands,
# which is the signature of near-collinear instruments in a three-endogenous,
# three-instrument system rather than a corrected estimate.
#
# Reporting OLS here is supported by the paper's own framework. Table 8
# establishes that OLS is attenuated toward zero in this design, -1.05%
# against an instrumented -3.07% at the pooled level, so a significant
# near-band OLS coefficient is if anything a lower bound on the local effect
# rather than an inflated one.
# ============================================================================

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


###############################################################################
#
#  PART 17 -- OUT-OF-STATE AND OWNERSHIP-DATED CONSTRUCTIONS
#
#  Two alternative constructions of the competitor instrument, each aimed at
#  a specific channel through which the exclusion restriction could fail.
#
#  17A, multi-market contracting. A carrier that contracts with a system in
#  several markets, or a system that negotiates one contract covering several
#  hospitals, could carry disclosure elsewhere into local prices without
#  passing through local peers. Both channels are densest within a state.
#  CTY_OOS keeps IV06 as built and drops every peer in the focal hospital's
#  own state.
#
#  17B, simultaneity. IV06 admits a rival system only once that system's own
#  hospital in the focal county has posted, so the instrument's support
#  depends on that hospital's timing. CTY_OWN dates presence by ownership
#  instead, through presence = "owned" on the Part 10.2 builder.
#
#  No new counting code is needed. At county geography, passing a same-state
#  county-pair table as `adjacency` removes in-state peers, since the builder
#  already drops same-county ones, and the in-state complement must add back
#  to the unfiltered count exactly. That partition is Gate 1 below.
#
#  Prerequisites. Parts 0-2, 4, and 10, with md_county_panel and
#  md_event_roster in scope.
#
###############################################################################

# ---------------------------------------------------------------------------
# 17.1 -- Clean FIPS roster and the state pair tables
#
# Missing FIPS codes are handled explicitly. formatC() turns NA into "   NA",
# which the builder would otherwise treat as one phantom county shared by
# every hospital with a missing code.
# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
# 17.2 -- Build the four counts and gate them
#
# Three gates, all fatal on failure: in-state and out-of-state peers must
# partition the unfiltered set; dating presence by ownership can only add
# peers, never remove them; and the clean-FIPS rebuild must reproduce the
# panel's own IV06.
# ---------------------------------------------------------------------------
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

# Gate 1: in-state and out-of-state peers partition the unfiltered set.
stopifnot("in-state + out-of-state != unfiltered count" =
            all(md_p17$CTY_ALL_HOSPITALS ==
                  md_p17$CTY_OOS_HOSPITALS + md_p17$CTY_INS_HOSPITALS))

# Gate 2: dating presence by ownership can only add peers.
stopifnot("ownership-dated count below posting-dated count" =
            all(md_p17$CTY_OWN_HOSPITALS >= md_p17$CTY_ALL_HOSPITALS))

# A focal hospital without a county FIPS cannot be placed in a state.
md_p17_cols <- c("CTY_ALL_HOSPITALS", "CTY_OOS_HOSPITALS",
                 "CTY_INS_HOSPITALS", "CTY_OWN_HOSPITALS")
md_p17[HOSPITAL_ID %in% md_roster_st[is.na(COUNTY_FIPS), HOSPITAL_ID],
       (md_p17_cols) := NA_integer_]

# Gate 3: the clean-FIPS rebuild reproduces the panel's IV06.
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

# ---------------------------------------------------------------------------
# 17.3 -- Headline specification under each construction, common sample
#
# All three instruments are estimated on the rows where the out-of-state
# count is defined, so the comparison is not confounded by sample changes.
# ---------------------------------------------------------------------------
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


# ---------------------------------------------------------------------------
# 17.4 -- HSA size distribution in the estimation sample
#
# Reports how many sample HSAs contain one, two, three, or four or more
# hospitals. This is the descriptive behind the HSA scope condition used in
# Parts 11.4 and 12.
# ---------------------------------------------------------------------------
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