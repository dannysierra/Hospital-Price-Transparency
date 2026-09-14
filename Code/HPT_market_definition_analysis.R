###############################################################################
#                                                                             #
#   HPT -- MARKET DEFINITION AND EXCLUSION RESTRICTION ROBUSTNESS             #
#                                                                             #
#   Companion to HPT_Analysis_Pipeline.R. Addresses two committee concerns:   #
#                                                                             #
#     1. MARKET DEFINITION (Katie). Is county the right market? Estimates     #
#        the headline shoppability gradient under six market definitions      #
#        (county, city, CBSA, HSA, HRR, 30-mile ring) plus rural/urban and    #
#        multi-hospital-market subsamples.                                    #
#                                                                             #
#     2. EXCLUSION RESTRICTION (Nick / adjacency). None of IV01-IV14 ever     #
#        excluded adjacent counties -- the only geographic filter every       #
#        variant inherits is "different county," and the CBSA variants add    #
#        "outside the CBSA," which is not the same thing since metros span    #
#        touching counties. IV15 is IV06 plus a genuine non-adjacency         #
#        filter. This file compares them and runs the own/adjacent            #
#        decomposition that tests whether the channel matters at all.        #
#                                                                             #
#   NOTHING HERE MODIFIES THE HEADLINE. Every existing instrument, panel,     #
#   and result is untouched. This file only adds rows to robustness tables.   #
#                                                                             #
#   PREREQUISITE. Run HPT_Analysis_Pipeline.R Parts 1-2 first (or            #
#   restore_session() via HPT_warm_start.R) so the helper functions and       #
#   constants exist. This file deliberately does not redefine them.          #
#                                                                             #
###############################################################################

# Sys.setenv(HPT_ROOT = "/Users/danielsierra/Library/CloudStorage/OneDrive-FloridaStateUniversity/Hospital Price Transparency Paper")
# HPT_WARM_START <- TRUE
# source(file.path(Sys.getenv("HPT_ROOT"), "Code", "HPT_Analysis_Pipeline.R"))

# ============================================================================
# PART 0 -- PREFLIGHT
# ============================================================================

MD_REQUIRED_OBJECTS <- c(
  "estimate_interacted", "prepare_panel", "read_panel", "build_schemes",
  "extend_schemes_for_merged", "attach_scheme_columns", "apply_concept_merges",
  "choose_sample_flag", "cache_or_run", "save_csv", "save_qa_csv",
  "available_columns", "model_sample", "build_iv_formula",
  "build_cluster_formula", "first_stage_wald", ".pval",
  "ENDOGENOUS_VARIABLE", "BASELINE_CONTROLS", "BASELINE_FIXED_EFFECTS",
  "BASELINE_CLUSTERS", "PRIMARY_OUTCOME", "PRIMARY_INSTRUMENT",
  "MAIN_INSTRUMENTS", "MIN_MODEL_OBS", "PANEL_DIR", "ANALYSIS_COLUMNS"
)

md_missing <- MD_REQUIRED_OBJECTS[!vapply(MD_REQUIRED_OBJECTS, exists, logical(1))]
if (length(md_missing) > 0L) {
  stop("Definitions from HPT_Analysis_Pipeline.R are not loaded.\n",
       "Missing: ", paste(md_missing, collapse = ", "), "\n\n",
       "Run Parts 1-2 of the main pipeline first, or restore_session().",
       call. = FALSE)
}

library(arrow)
library(data.table)
library(fixest)

cat("\n=== MARKET DEFINITION ANALYSIS ===\n")
cat("Endogenous:", ENDOGENOUS_VARIABLE, "| Primary IV:", PRIMARY_INSTRUMENT, "\n")


# ---------------------------------------------------------------------------
# 0.1  CRITICAL: extend ANALYSIS_COLUMNS so IV15 survives read_panel()
# ---------------------------------------------------------------------------
# read_panel() filters columns with intersect(ANALYSIS_COLUMNS, names(ds)).
# The IV15 family is not in ALL_CANDIDATE_INSTRUMENTS, so without this the
# new columns are SILENTLY DROPPED at load -- no error, the column simply is
# not there, and any model using it returns NULL. This is the single most
# likely way for this file to appear to "work" while doing nothing.

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
# NOTE ON THE UNIT OF ANALYSIS. The headline runs on the CONCEPT panel
# (FILES$outpatient_concept), not the exact-code panel, so the ladder uses
# the CONCEPT panels to mirror the headline specification exactly.

MD_GEOGRAPHIES <- c("COUNTY", "CITY", "CBSA", "HSA", "HRR")

md_panel_path <- function(level) {
  role <- if (level == "COUNTY") "PRIMARY" else "ROBUSTNESS"
  file.path(PANEL_DIR,
            sprintf("HPT_R_MAIN_%s_%s_OUTPATIENT_CONCEPT.parquet", role, level))
}

# FIXED EFFECTS AND CLUSTERING MOVE AUTOMATICALLY.
# Each exported panel renames its own geography to the generic ANALYSIS_MARKET
# and rebuilds MARKET_ID as ANALYSIS_MARKET::FINAL_CONCEPT_ID. Because
# BASELINE_FIXED_EFFECTS and BASELINE_CLUSTERS reference those generic names,
# swapping the panel swaps the FE and the cluster level together. No manual
# adjustment is needed and none should be added -- doing so would silently
# desynchronise the FE from the treatment.

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

# Four hospitals (La Salle IL, La Porte IN, 2x Valdez-Cordova AK) appear twice
# in HPT_HOSPITAL_GEO because the Snowflake Section 3 join matched them via
# both HPT_REF_GEO_OVERRIDES and the name crosswalk once the regex fix made
# the normal path work. The duplicate rows are byte-identical, so dedup is
# safe. The permanent fix belongs in the SQL override join.
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
# Mirrors load_outpatient() exactly, with the geography file swapped and the
# hospital attributes merged on. Any deviation from this sequence would make
# the ladder non-comparable to the headline.
#
# BUG FOUND AND FIXED (this session): apply_concept_merges() pulls in a
# second, always-county-level exact-code panel internally to resolve the six
# canonical concept collapses (MRI brain, CT abdomen, etc.). That secondary
# panel's own county-keyed ANALYSIS_MARKET was leaking into the PRIMARY
# panel's ANALYSIS_MARKET for every row touched by that merge -- invisible
# for the COUNTY geography (leaked value == correct value, so nothing looked
# wrong) but real corruption for CITY/CBSA/HSA/HRR, where it silently
# overwrote the correct market key with a raw county string for ~14,000 rows
# per geography. Confirmed directly: corrupted row count matched the merge's
# own "processed N groups" log line exactly, and every corrupted concept was
# one of the six canonical merge targets and no others.
#
# Fix: save ANALYSIS_MARKET before calling apply_concept_merges(), restore it
# after, rebuild MARKET_ID to match. apply_concept_merges() itself is
# untouched -- it is existing, validated code the headline result depends on,
# and editing it blind under deadline is a worse risk than working around it
# from outside.

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
# Headline specification, re-estimated once per market definition. Instrument
# held at PRIMARY_INSTRUMENT throughout.
#
# ONE CHOICE TO DOCUMENT. Every Z_SYS_* instrument is constructed at COUNTY
# level in the Python pipeline. When the treatment moves to HRR, the
# instrument stays county-based. That is defensible -- the instrument's job is
# to predict the posting decision, not to define the market -- but it is a
# choice, not an automatic consequence, and it belongs in the table notes.

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
# PART 4 -- RING EXPOSURE AS A SIXTH RUNG, AND THE DECAY DIAGNOSTIC
# ============================================================================

md_county_panel <- md_load_panel("COUNTY")

# 4.1  Ring as a market definition: swap the treatment, keep county FE.
#
# The ring has no discrete market key, so ANALYSIS_MARKET and MARKET_ID stay
# county-based while the treatment becomes the 30-mile prior-poster count.
# This rung is therefore NOT a clean substitute for the others -- it changes
# the treatment without changing the FE -- and should be reported as a
# distance-based treatment check rather than as a sixth market definition.

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

# 4.2  Spatial decay diagnostic.
#
# Reduced form with all three distance bands entered simultaneously. This does
# not define a market -- it estimates where the peer-posting relationship dies
# out, which is the evidence-based answer to "what is a market."
#
# READ THIS AGAINST THE RAW SHARES. Descriptively, the share of ring peers
# that posted first is essentially flat across bands (0.398 / 0.402 / 0.407
# for 0-15, 15-30, 30-60 miles). If the regression bands are also flat, that
# is a substantive finding about the setting, not a data problem: peer posting
# is not spatially local. Frame it explicitly rather than burying it.

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


# ============================================================================
# PART 5 -- ADJACENCY DECOMPOSITION
# ============================================================================
#
# Own-county and adjacent-county prior-poster counts entered as two separate
# regressors. If the adjacent coefficient is near zero, the county boundary is
# defensible and you have shown it rather than assumed it. If it is large, the
# cross-border channel is real and measured.
#
# Descriptively these are NOT redundant: own-county mean 2.343, adjacent-county
# mean 5.711, correlation 0.468 -- low enough to identify both.

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
# PART 6 -- IV15 VERSUS IV06
# ============================================================================
#
# The headline specification, re-estimated with the adjacency-excluded
# instrument in place of the primary one. Descriptively IV15 drops 19.6% of
# IV06's peers (mean 7.168 -> 5.763) while correlating 0.9806, so the filter
# does real work without gutting the instrument.
#
# IV15 is NaN for 0.78% of hospital-months, essentially all Connecticut, where
# the 2022 planning-region reorganisation leaves no usable county FIPS.

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


# ============================================================================
# PART 7 -- SUBSAMPLES: RURAL/URBAN AND MULTI-HOSPITAL MARKETS
# ============================================================================
#
# Katie's second and third points. 2,551 metro against 1,138 nonmetro
# hospitals, and 2,541 of 3,724 hospitals (68.2%) sit in counties with at
# least two hospitals, so both restrictions leave estimable samples.
#
# The multi-hospital restriction is where "market" has content: in a
# single-hospital county the concept is empty and the treatment has no
# within-market variation to give.

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
# PART 9 -- CONCEPT-LEVEL ANALYSIS AT HSA AND HRR
#           (replaces the earlier two-bucket cross-over, which was not a
#            valid test -- see note below)
# ============================================================================
#
# WHY THE EARLIER VERSION WAS WRONG. Part 9 originally split the panel into
# two coarse buckets (Shoppable / Non_shoppable under SCHEME_1 only) and
# compared pooled coefficients across geographies. That cannot test the
# HSA/HRR hypothesis, for a structural reason: HRRs are built by AGGREGATING
# HSAs -- every HSA nests inside exactly one HRR. Averaging hundreds of
# concepts of wildly different complexity into two buckets washes out exactly
# the concept-specific variation the hypothesis is about. Its apparent
# "HRR favoured for both tiers" result is uninformative, not evidence against.
#
# WHAT THIS DOES INSTEAD. Runs the SAME concept-level estimator the headline
# uses (estimate_concept_level -> run_meta_regressions) once per geography.
# Two things come for free from reusing that machinery:
#   - every concept gets its own IV coefficient per instrument, so roughly
#     738 concepts x 3 instruments ~ 2,200 rows per geography
#   - run_meta_regressions() loops over EVERY SCHEME_ID in schemes_long, so
#     all 18 shoppability schemes are covered, not just SCHEME_1
#
# THE ACTUAL TEST (Part 9C). For each concept, compute its HSA coefficient
# minus its HRR coefficient. Negative delta = HSA disciplines price more =
# behaves like a routine, local-catchment service. Positive = HRR more =
# behaves like a referral-catchment service. Then check whether that delta
# tracks shoppability. That is the cross-over, at the right granularity.

# Instruments: MAIN_INSTRUMENTS (3), not the 6 used in the headline concept
# run (CONCEPT_INSTRUMENTS = MAIN + SUPPORTING). Market definition and
# instrument choice are separate axes, and Part 6 already validated IV06 vs
# IV15 at county level, so re-running all 6 at every geography answers a
# question not being asked here. Three covers the credibility-tier "Main"
# instruments at roughly half the cost. NOTE for table notes: the headline
# -2.603pp gradient is a 6-instrument number, so these are not directly
# comparable to it without that caveat stated.
CONCEPT_INSTRUMENTS_MD <- MAIN_INSTRUMENTS

# Cache and output keys encode the instrument COUNT, mirroring the original
# pipeline's "concept_level_6inst". cache_or_run() keys purely on the string,
# so without this a later run with a different instrument set would return
# tonight's 3-instrument object and announce "Cache hit" as though correct.
# The output CSVs carry the same suffix so two vintages never overwrite each
# other unnoticed.
md_key <- function(stem, level, instruments = CONCEPT_INSTRUMENTS_MD) {
  sprintf("%s_%s_%dinst", stem, tolower(level), length(instruments))
}

md_run_concept_level <- function(level, instruments = CONCEPT_INSTRUMENTS_MD) {
  cat("\n", strrep("=", 70), "\n", sep = "")
  cat("CONCEPT-LEVEL:", level, "| instruments:", length(instruments), "\n")
  cat("  ", paste(names(instruments), collapse = "\n   "), "\n", sep = "")
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
  
  cr[, MARKET_DEFINITION := level]
  cr[, N_INSTRUMENTS_USED := length(instruments)]
  
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


# ---------------------------------------------------------------------------
# 9A -- HSA. Run this alone first and check it before starting HRR.
# ---------------------------------------------------------------------------
md_hsa <- md_run_concept_level("HSA")

cat("\n=== HSA META-REGRESSION (all schemes) ===\n")
print(md_hsa$meta_iv)


# ---------------------------------------------------------------------------
# 9B -- HRR. Only after 9A looks sane.
# ---------------------------------------------------------------------------
md_hrr <- md_run_concept_level("HRR")

cat("\n=== HRR META-REGRESSION (all schemes) ===\n")
print(md_hrr$meta_iv)


# ---------------------------------------------------------------------------
# 9C -- THE CROSS-OVER TEST: per-concept HSA vs HRR delta
# ---------------------------------------------------------------------------
#
# Both sides carry one row per concept x instrument, so the merge is keyed on
# BOTH, not on concept alone. Merging on concept alone would fan out to
# 3 x 3 = 9 rows per concept and silently corrupt every downstream mean.

md_delta <- merge(
  md_hsa$concept[, .(FINAL_CONCEPT_ID, INSTRUMENT_LABEL, FINAL_FAMILY_ID, SERVICE_LABEL,
                     IV_HSA = IV_COEF, IV_SE_HSA = IV_SE, IV_P_HSA = IV_P,
                     FS_F_HSA = FS_F, N_OBS_HSA = N_OBSERVATIONS)],
  md_hrr$concept[, .(FINAL_CONCEPT_ID, INSTRUMENT_LABEL,
                     IV_HRR = IV_COEF, IV_SE_HRR = IV_SE, IV_P_HRR = IV_P,
                     FS_F_HRR = FS_F, N_OBS_HRR = N_OBSERVATIONS)],
  by = c("FINAL_CONCEPT_ID", "INSTRUMENT_LABEL"))

# Negative delta: HSA gives a more negative (stronger) price response than
# HRR for this concept -> routine, local-catchment behaviour.
md_delta[, IV_DELTA_HSA_MINUS_HRR := IV_HSA - IV_HRR]

# Keep only concepts where BOTH geographies produced a usable first stage.
# A delta built on a weak first stage on either side is noise, not signal.
md_delta_ok <- md_delta[is.finite(IV_DELTA_HSA_MINUS_HRR) &
                          FS_F_HSA >= 10 & FS_F_HRR >= 10]

cat(sprintf("\nConcept x instrument pairs in both geographies: %d | F>=10 both sides: %d\n",
            nrow(md_delta), nrow(md_delta_ok)))
cat("Distinct concepts surviving:", uniqueN(md_delta_ok$FINAL_CONCEPT_ID), "\n")
print(md_delta_ok[, .N, by = INSTRUMENT_LABEL])

# Attach every scheme's classification, then test the delta against each.
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
  
  d[, .(SCHEME_ID = sid,
        N_CONCEPTS = .N,
        MEAN_DELTA = round(mean(IV_DELTA_HSA_MINUS_HRR, na.rm = TRUE), 4),
        MEDIAN_DELTA = round(median(IV_DELTA_HSA_MINUS_HRR, na.rm = TRUE), 4),
        SHARE_HSA_STRONGER = round(mean(IV_DELTA_HSA_MINUS_HRR < 0, na.rm = TRUE), 3)),
    by = c(CATEGORY = ccol, INSTRUMENT = "INSTRUMENT_LABEL")]
}), fill = TRUE)

# The sharper version: regress the per-concept delta on the ordinal
# shoppability score, once per scheme x instrument. A negative slope means
# more-shoppable concepts have more negative deltas, i.e. HSA disciplines them
# more than HRR does -- which is the cross-over the hypothesis predicts.
md_crossover_slope <- rbindlist(lapply(unique(md_schemes_long$SCHEME_ID), function(sid) {
  ocol <- paste0("ORD_", sid)
  if (!(ocol %in% names(md_delta_schemes))) return(NULL)
  
  rbindlist(lapply(unique(md_delta_schemes$INSTRUMENT_LABEL), function(il) {
    d <- md_delta_schemes[INSTRUMENT_LABEL == il &
                            is.finite(get(ocol)) & is.finite(IV_DELTA_HSA_MINUS_HRR)]
    if (nrow(d) < MIN_CONCEPTS_META || !has_usable_variation(d[[ocol]])) return(NULL)
    
    y <- d$IV_DELTA_HSA_MINUS_HRR; x <- d[[ocol]]
    fit <- tryCatch(lm(y ~ x), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    s <- summary(fit)$coefficients
    if (nrow(s) < 2L) return(NULL)
    
    data.table(SCHEME_ID = sid, INSTRUMENT = il, N_CONCEPTS = nrow(d),
               SLOPE = round(s[2, 1], 5), SE = round(s[2, 2], 5),
               P_VALUE = round(s[2, 4], 4))
  }), fill = TRUE)
}), fill = TRUE)

isuf <- sprintf("_%dinst.csv", length(CONCEPT_INSTRUMENTS_MD))
save_csv(md_delta,           paste0("MD16_concept_hsa_hrr_delta", isuf))
save_csv(md_crossover_test,  paste0("MD17_crossover_by_scheme", isuf))
save_csv(md_crossover_slope, paste0("MD17B_crossover_slope_by_scheme", isuf))

cat("\n=== CROSS-OVER: mean HSA-minus-HRR delta by scheme category ===\n")
cat("Negative = HSA stronger (routine/local). Positive = HRR stronger (referral).\n\n")
print(md_crossover_test)

cat("\n=== CROSS-OVER SLOPE: delta regressed on ordinal shoppability ===\n")
cat("Negative slope = more shoppable concepts favour HSA. That is the prediction.\n\n")
print(md_crossover_slope)

# Family-level view: which clinical families lean local vs referral?
md_family_delta <- md_delta_ok[, .(
  N_PAIRS      = .N,
  MEAN_DELTA   = round(mean(IV_DELTA_HSA_MINUS_HRR, na.rm = TRUE), 4),
  MEDIAN_DELTA = round(median(IV_DELTA_HSA_MINUS_HRR, na.rm = TRUE), 4),
  SHARE_HSA_STRONGER = round(mean(IV_DELTA_HSA_MINUS_HRR < 0, na.rm = TRUE), 3)
), by = .(FINAL_FAMILY_ID, INSTRUMENT_LABEL)][order(MEDIAN_DELTA)]

save_csv(md_family_delta, paste0("MD18_family_hsa_hrr_delta", isuf))
cat("\n=== BY CLINICAL FAMILY (sorted: most HSA-favouring first) ===\n")
print(md_family_delta)



# ============================================================================
# PART 8 -- ASSEMBLED SUMMARY
# ============================================================================

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