###############################################################################
# DIAGNOSTIC BLOCK
#
# Run AFTER Sections 0-12 are loaded, with `outpatient` and `schemes_long` in
# memory. Read-only: touches nothing, defines nothing outside `.dx_*`.
# Runtime under a minute. Paste the entire console output back.
#
# Purpose: establish (a) whether FINAL_SUPERFAMILY_ID survived loading and what
# its values are, (b) whether family and superfamily subsets clear the existing
# screens, (c) whether shoppability is well defined at those levels or has to be
# collapsed by a rule, (d) whether treatment and instruments are invariant within
# hospital x month, (e) what the enforcement columns look like, (f) the headline
# gradient numbers that set the scale for the Conley grid.
###############################################################################

.dx_hd <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")
.dx_show <- function(d, n = Inf, label = NULL) {
  if (!is.null(label)) cat("\n--- ", label, " ---\n", sep = "")
  d <- as.data.frame(d)
  if (is.finite(n) && nrow(d) > n) {
    print(head(d, n), row.names = FALSE)
    cat("  [", nrow(d) - n, " further rows suppressed]\n", sep = "")
  } else {
    print(d, row.names = FALSE)
  }
  invisible(NULL)
}


# ===========================================================================
# 0. SESSION AND OBJECTS
# ===========================================================================
.dx_hd("0. SESSION AND OBJECTS")

cat("R", paste(R.version$major, R.version$minor, sep = "."),
    "| fixest", as.character(packageVersion("fixest")),
    "| data.table", as.character(packageVersion("data.table")), "\n")

cat("In memory:", paste(intersect(
  c("outpatient", "schemes_long", "concept_results", "main", "mi", "mr", "perm"),
  ls(envir = .GlobalEnv)), collapse = ", "), "\n")

if (!exists("outpatient")) {
  stop("`outpatient` is not in memory.\n",
       "  In a fresh session, load a completed run without re-estimating:\n",
       "    source(\"~/path/to/HPT_warm_start.R\"); warm_start()\n",
       "  Or manually: run lines 1-2956 of HPT_Analysis_Pipeline.R (definitions),\n",
       "  then lines 7884-end (restore utilities), then restore_session().\n",
       "  Do NOT run the BUILD block -- restore_session() supplies `outpatient`.",
       call. = FALSE)
}
if (!exists("schemes_long")) {
  warning("`schemes_long` is not in memory. Checks 1-9 still run; only the\n",
          "  scheme cross-tabs would use it.", call. = FALSE)
}
cat("outpatient:", format(nrow(outpatient), big.mark = ","), "rows,",
    ncol(outpatient), "columns\n")


# ===========================================================================
# 1. COLUMN INVENTORY
# ===========================================================================
.dx_hd("1. COLUMN INVENTORY")
print(sort(names(outpatient)))


# ===========================================================================
# 2. FIRST FIVE ROWS, KEY COLUMNS, TRANSPOSED
# ===========================================================================
.dx_hd("2. FIRST 5 ROWS (transposed so long IDs stay readable)")

.dx_key <- intersect(c(
  "HOSPITAL_ID", "POST_MONTH", "ANALYSIS_MARKET", "COUNTY_FIPS", "CBSA_CODE",
  "MARKET_ID", "ANALYSIS_SERVICE_ID", "ANALYSIS_AGGREGATION",
  "FINAL_SUPERFAMILY_ID", "FINAL_FAMILY_ID", "FINAL_FAMILY_NAME",
  "FINAL_CONCEPT_ID", "FINAL_CONCEPT_NAME", "BILLING_CODE_TYPE",
  "MEDIAN_PRICE", "LN_MEDIAN_PRICE", "LOG_TOTAL_BEDS", "TOTAL_BEDS",
  "N_PRIOR_POSTERS", PRIMARY_INSTRUMENT,
  unname(SCHEME_COLUMNS)), names(outpatient))

print(t(as.data.frame(head(outpatient[, .dx_key, with = FALSE], 5))))


# ===========================================================================
# 3. SUPERFAMILY -- DOES IT EXIST, AND WHAT IS IN IT
# ===========================================================================
.dx_hd("3. SUPERFAMILY (candidate major-family level)")

if (!("FINAL_SUPERFAMILY_ID" %in% names(outpatient))) {
  cat("FINAL_SUPERFAMILY_ID IS NOT IN THE PANEL.\n",
      "read_panel() drops absent columns silently, so it is either missing from\n",
      "the parquet or named differently. Candidates present:\n", sep = "")
  print(grep("SUPER|GROUP|DOMAIN|MODALITY|CATEGORY", names(outpatient),
             value = TRUE, ignore.case = TRUE))
} else {
  .dx_show(outpatient[, .(
    N_ROWS     = .N,
    N_FAMILIES = uniqueN(FINAL_FAMILY_ID),
    N_CONCEPTS = uniqueN(ANALYSIS_SERVICE_ID),
    N_HOSP     = uniqueN(HOSPITAL_ID),
    N_MARKETS  = uniqueN(ANALYSIS_MARKET),
    N_MONTHS   = uniqueN(POST_MONTH),
    N_FE_CELLS = uniqueN(MARKET_ID)
  ), by = FINAL_SUPERFAMILY_ID][order(-N_ROWS)],
  label = "Superfamily totals")

  .dx_show(unique(outpatient[, .(FINAL_SUPERFAMILY_ID, FINAL_FAMILY_ID)])[
    order(FINAL_SUPERFAMILY_ID, FINAL_FAMILY_ID)],
  label = "Family -> superfamily crosswalk")
}


# ===========================================================================
# 4. FAMILY TAXONOMY AND SCREEN FEASIBILITY
# ===========================================================================
.dx_hd("4. FAMILY TAXONOMY AND SCREEN FEASIBILITY")

.dx_by <- intersect(c("FINAL_SUPERFAMILY_ID", "FINAL_FAMILY_ID"), names(outpatient))

.dx_fam <- outpatient[, .(
  N_ROWS     = .N,
  N_CONCEPTS = uniqueN(ANALYSIS_SERVICE_ID),
  N_HOSP     = uniqueN(HOSPITAL_ID),
  N_MARKETS  = uniqueN(ANALYSIS_MARKET),
  N_MONTHS   = uniqueN(POST_MONTH),
  N_FE_CELLS = uniqueN(MARKET_ID),
  MEAN_TREAT = round(mean(safe_numeric(get(ENDOGENOUS_VARIABLE)), na.rm = TRUE), 2),
  SD_Z       = round(sd(safe_numeric(get(PRIMARY_INSTRUMENT)), na.rm = TRUE), 3)
), by = .dx_by][order(-N_ROWS)]

.dx_fam[, IS_DIAGNOSTIC_FAMILY := as.integer(FINAL_FAMILY_ID %chin% DIAGNOSTIC_FAMILIES)]
.dx_fam[, PASSES_SCREENS := as.integer(
  N_ROWS >= MIN_SERVICE_OBS & N_MARKETS >= MIN_SERVICE_MARKETS &
    N_MONTHS >= MIN_SERVICE_MONTHS)]

.dx_show(.dx_fam)

cat("\nScreens in force: MIN_SERVICE_OBS =", MIN_SERVICE_OBS,
    "| MIN_SERVICE_MARKETS =", MIN_SERVICE_MARKETS,
    "| MIN_SERVICE_MONTHS =", MIN_SERVICE_MONTHS,
    "\nMIN_MODEL_OBS (interacted) =", MIN_MODEL_OBS,
    "| MIN_CONCEPTS_META =", MIN_CONCEPTS_META, "\n")


# ===========================================================================
# 5. IS SHOPPABILITY WELL DEFINED AT FAMILY AND SUPERFAMILY LEVEL?
#
# A family-level estimate needs a family-level label. If a scheme splits a
# family internally, the label has to come from a collapse rule (modal
# category, unanimity, or drop-if-mixed) and that rule has to be stated in the
# paper. This reports how often the question arises.
# ===========================================================================
.dx_hd("5. WITHIN-FAMILY SCHEME UNIFORMITY")

.dx_sc <- intersect(unname(SCHEME_COLUMNS), names(outpatient))
.dx_cc <- unique(outpatient[, c(.dx_by, "ANALYSIS_SERVICE_ID", .dx_sc), with = FALSE])
cat("Distinct concepts carrying scheme labels:", nrow(.dx_cc), "\n")

for (s in .dx_sc) {
  a <- copy(.dx_cc)
  a[, CAT := toupper(trimws(as.character(get(s))))]
  a <- a[!is.na(CAT) & CAT != ""]
  if (nrow(a) == 0L) { cat("\n", s, ": no non-missing labels\n", sep = ""); next }

  tab <- a[, .N, by = c("FINAL_FAMILY_ID", "CAT")]
  fam <- tab[order(FINAL_FAMILY_ID, -N),
             .(N_CONCEPTS = sum(N), N_CATS = .N, MODE = CAT[1L],
               MODE_SHARE = round(N[1L] / sum(N), 3)), by = FINAL_FAMILY_ID]

  cat("\n", s, ": ", sum(fam$N_CATS == 1L), " of ", nrow(fam),
      " families uniform | mixed families listed below\n", sep = "")
  .dx_show(fam[N_CATS > 1L][order(-N_CONCEPTS)], n = 12)

  if ("FINAL_SUPERFAMILY_ID" %in% names(a)) {
    tabs <- a[, .N, by = c("FINAL_SUPERFAMILY_ID", "CAT")]
    sup <- tabs[order(FINAL_SUPERFAMILY_ID, -N),
                .(N_CONCEPTS = sum(N), N_CATS = .N, MODE = CAT[1L],
                  MODE_SHARE = round(N[1L] / sum(N), 3)), by = FINAL_SUPERFAMILY_ID]
    .dx_show(sup[order(-N_CONCEPTS)], label = paste0(s, " at superfamily level"))
  }
}


# ===========================================================================
# 6. ARE TREATMENT AND INSTRUMENTS INVARIANT WITHIN HOSPITAL x MONTH?
#
# If they are, a family-level or superfamily-level unit can be formed without
# any decision about how to aggregate the right-hand side. If they are not,
# every aggregation needs a weighting rule and the exercise gets much harder.
# ===========================================================================
.dx_hd("6. VARIATION OF RHS WITHIN HOSPITAL x MONTH")

.dx_iv <- intersect(unname(ALL_SIX_INSTRUMENTS), names(outpatient))
.dx_rhs <- intersect(c(ENDOGENOUS_VARIABLE, "LOG_TOTAL_BEDS", .dx_iv), names(outpatient))

.dx_u <- outpatient[, lapply(.SD, uniqueN), by = .(HOSPITAL_ID, POST_MONTH),
                    .SDcols = .dx_rhs]

.dx_show(data.table(
  VARIABLE = .dx_rhs,
  SHARE_HOSP_MONTH_CELLS_WITH_MULTIPLE_VALUES = vapply(
    .dx_rhs, function(v) round(mean(.dx_u[[v]] > 1L), 5), numeric(1)),
  MAX_DISTINCT_IN_A_CELL = vapply(
    .dx_rhs, function(v) max(.dx_u[[v]]), numeric(1))
))

cat("\nHospital x month cells:", format(nrow(.dx_u), big.mark = ","), "\n")
rm(.dx_u); invisible(gc())


# ===========================================================================
# 7. ENFORCEMENT VARIABLES
# ===========================================================================
.dx_hd("7. ENFORCEMENT VARIABLES")

.dx_enf <- intersect(unname(ENFORCEMENT_ROBUSTNESS_INSTRUMENTS), names(outpatient))
cat("Present in panel:", length(.dx_enf), "of",
    length(ENFORCEMENT_ROBUSTNESS_INSTRUMENTS), "\n")
if (length(.dx_enf) < length(ENFORCEMENT_ROBUSTNESS_INSTRUMENTS)) {
  cat("Absent:", paste(setdiff(unname(ENFORCEMENT_ROBUSTNESS_INSTRUMENTS),
                               .dx_enf), collapse = ", "), "\n")
}

if (length(.dx_enf) > 0L) {
  .dx_show(rbindlist(lapply(.dx_enf, function(v) {
    x <- safe_numeric(outpatient[[v]])
    data.table(VARIABLE = v,
               SHARE_NA    = round(mean(is.na(x)), 4),
               MEAN        = round(mean(x, na.rm = TRUE), 4),
               SD          = round(sd(x, na.rm = TRUE), 4),
               SHARE_ZERO  = round(mean(x == 0, na.rm = TRUE), 4),
               MAX         = suppressWarnings(max(x, na.rm = TRUE)),
               N_DISTINCT  = uniqueN(x))
  })))

  # Is enforcement constant within county x month? It should be, by construction.
  .dx_cm <- unique(outpatient[, c("ANALYSIS_MARKET", "POST_MONTH", .dx_enf),
                              with = FALSE])
  .dx_dup <- .dx_cm[, .N, by = .(ANALYSIS_MARKET, POST_MONTH)][N > 1L, .N]
  cat("\nDistinct county-month cells:", format(nrow(.dx_cm), big.mark = ","),
      "| cells carrying MORE THAN ONE enforcement vector:", .dx_dup,
      "\n(should be 0 if enforcement is truly county-month)\n")

  # Within-county time variation is what makes enforcement identifiable
  # alongside POST_MONTH fixed effects.
  .dx_show(rbindlist(lapply(.dx_enf, function(v) {
    w <- .dx_cm[, .(SDV = sd(safe_numeric(get(v)), na.rm = TRUE),
                    ANY = as.integer(any(safe_numeric(get(v)) > 0, na.rm = TRUE))),
                by = ANALYSIS_MARKET]
    data.table(VARIABLE = v,
               N_COUNTIES = nrow(w),
               SHARE_COUNTIES_EVER_NONZERO = round(mean(w$ANY == 1L), 4),
               SHARE_COUNTIES_WITH_TIME_VARIATION =
                 round(mean(w$SDV > 0, na.rm = TRUE), 4))
  })), label = "Within-county time variation")
  rm(.dx_cm); invisible(gc())
}


# ===========================================================================
# 8. CONCEPT COVERAGE WITHIN HOSPITAL x MONTH x FAMILY
#
# Only matters if the family-level design turns out to need an aggregated price
# index rather than a subsample. If the concept set moves around inside a
# hospital-family over time, a raw median-of-medians index confounds price
# changes with composition changes.
# ===========================================================================
.dx_hd("8. CONCEPT COVERAGE WITHIN HOSPITAL x MONTH x FAMILY")

.dx_cov <- outpatient[, .(N_CONCEPTS = .N),
                      by = .(HOSPITAL_ID, POST_MONTH, FINAL_FAMILY_ID)]

.dx_show(.dx_cov[, .(CELLS = .N,
                     MEAN_CONCEPTS = round(mean(N_CONCEPTS), 2),
                     MEDIAN_CONCEPTS = as.numeric(median(N_CONCEPTS)),
                     P90 = quantile(N_CONCEPTS, .90, names = FALSE),
                     MAX = max(N_CONCEPTS)),
                 by = FINAL_FAMILY_ID][order(-CELLS)])

.dx_bal <- .dx_cov[, .(SDN = sd(N_CONCEPTS), NM = .N),
                   by = .(HOSPITAL_ID, FINAL_FAMILY_ID)]
cat("\nHospital-family cells observed in >1 month:",
    format(.dx_bal[NM > 1L, .N], big.mark = ","),
    "\nOf those, share whose concept COUNT never changes across months:",
    round(.dx_bal[NM > 1L, mean(SDN == 0, na.rm = TRUE)], 4), "\n")
rm(.dx_cov, .dx_bal); invisible(gc())


# ===========================================================================
# 9. HEADLINE GRADIENT -- sets the scale for the Conley grid
#
# The Conley bound is reported as a fraction of the observed reduced-form
# response on shoppable services, so I need those coefficients and the
# interacted first-stage Wald to write sensible defaults.
# ===========================================================================
.dx_hd("9. HEADLINE GRADIENT")

.dx_t06 <- file.path(TABLE_DIR, "T06_main_interacted_RF_and_IV.csv")
if (file.exists(.dx_t06)) {
  .dx_m <- fread(.dx_t06)
  .dx_cols <- intersect(c("SPEC", "INSTRUMENT_LABEL", "TERM", "RF_COEF", "RF_SE",
                          "RF_PERCENT_PER_SD", "IV_COEF", "IV_SE", "IV_PERCENT",
                          "FIRST_STAGE_WALD_MIN", "N_OBSERVATIONS"),
                        names(.dx_m))
  .dx_show(.dx_m[grepl("certainty", SPEC, ignore.case = TRUE),
                 .dx_cols, with = FALSE])
} else {
  cat("T06 not on disk at:", .dx_t06, "\n",
      "If `main` is in memory, run: .dx_show(main$rows)\n")
}

cat("\nSD of the primary instrument in the full panel:",
    round(sd(safe_numeric(outpatient[[PRIMARY_INSTRUMENT]]), na.rm = TRUE), 5), "\n")


# ===========================================================================
# 10. CACHE STATE
# ===========================================================================
.dx_hd("10. CACHE STATE")
if (exists("cache_status")) invisible(cache_status()) else
  print(sort(list.files(CACHE_DIR, pattern = "\\.rds$")))

cat("\n", strrep("=", 78), "\nDIAGNOSTIC COMPLETE\n", strrep("=", 78), "\n", sep = "")
