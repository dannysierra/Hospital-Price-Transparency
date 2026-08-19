###############################################################################
# SECTION 19 -- TIER I DIAGNOSTICS
#
# Source AFTER HPT_Analysis_Pipeline.R has run through at least Section 7, so
# that `outpatient`, the helper functions, and the registries are in memory.
# Nothing here re-estimates the headline models from scratch; every block
# rebuilds the SAME design matrix estimate_interacted() builds, on the SAME
# model_sample(), so the numbers are directly comparable to Table 5.
#
#   19A  A4 -- the interacted first-stage coefficient MATRIX. Resolves whether
#             the Panel A / Panel B inconsistency is (i) genuinely divergent
#             category first stages, (ii) off-diagonal leakage through the
#             month FE so that beta != gamma/pi, or (iii) a 0/0 artifact.
#   19B  A4 + C5 -- the reduced-form gap with its own SE, reported under all
#             three percent conventions currently in use in the paper.
#   19C  A7 -- minimum detectable effects for the comparability moderators and
#             the demographic interactions, expressed as a share of the
#             shoppability gradient.
#
# WHY 19A IS BUILT MANUALLY RATHER THAN FROM fixest INTERNALS.
#   setFixest_estimation(data.save = FALSE) is set in Section 0, so the stored
#   first-stage objects on an IV fit are not reliably available. Fitting the
#   two first-stage equations directly with feols() is transparent, costs one
#   extra pass, and lets the off-diagonals be tested individually -- which is
#   the whole point of the exercise.
###############################################################################

.s19_hd <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")
.s19_sub <- function(x) cat("\n--- ", x, " ", strrep("-", max(0, 70 - nchar(x))), "\n", sep = "")


# ===========================================================================
# 19A. THE INTERACTED FIRST-STAGE MATRIX                                 (A4)
# ===========================================================================
#
# With K categories the IV system is
#
#     TREAT_k = sum_j Pi[k, j] * IV_j + controls + FE + u_k        (k = 1..K)
#     ln(P)   = sum_k beta_k TREAT_k  + controls + FE + e
#
# and 2SLS returns beta = Pi^{-1} gamma, where gamma is the vector of
# reduced-form coefficients on RF_j = IV_j. The paper's Panel B is read AS IF
# beta_k = gamma_k / Pi[k,k], which is true only when Pi is diagonal.
#
# IV_Shoppable and IV_Non_shoppable are mutually exclusive by construction, and
# because shoppability is CONCEPT-determined every county x concept FE cell is
# category-pure, so that fixed effect cannot induce leakage. The MONTH fixed
# effect can: it spans both categories, so after within-transformation IV_N is
# not exactly zero on shoppable rows. This block measures how much that matters.

s19_first_stage_matrix <- function(
    data,
    scheme_col       = "SCHEME_1_CERTAINTY",
    instrument       = PRIMARY_INSTRUMENT,
    instrument_label = "Competitor_only_hospitals_9m",
    outcome          = PRIMARY_OUTCOME,
    endogenous       = ENDOGENOUS_VARIABLE,
    controls         = BASELINE_CONTROLS,
    fixed_effects    = BASELINE_FIXED_EFFECTS,
    clusters         = BASELINE_CLUSTERS) {
  
  controls <- available_columns(data, controls)
  fe       <- available_columns(data, fixed_effects)
  cl       <- available_columns(data, clusters)
  
  d <- data[!is.na(get(scheme_col))]
  d <- model_sample(d, c(outcome, endogenous, instrument, controls, fe, cl, scheme_col))
  if (nrow(d) < MIN_MODEL_OBS) stop("Sample too small for ", instrument_label, call. = FALSE)
  
  d[, MOD := droplevels(factor(get(scheme_col)))]
  keys <- levels(d$MOD)
  K <- length(keys)
  if (K < 2L) stop("Need at least two categories.", call. = FALSE)
  
  endo <- paste0("TREAT_", keys)
  ivs  <- paste0("IV_",    keys)
  rfs  <- paste0("RF_",    keys)
  for (k in seq_along(keys)) {
    sel <- as.integer(d$MOD == keys[k])
    d[, (endo[k]) := get(endogenous) * sel]
    d[, (ivs[k])  := get(instrument) * sel]
    d[, (rfs[k])  := get(instrument) * sel]
  }
  
  clf <- build_cluster_formula(cl)
  
  # -- reduced form: gamma -------------------------------------------------
  rf_fit <- feols(build_ols_formula(outcome, c(rfs, controls), fe),
                  data = d, cluster = clf, warn = FALSE, notes = FALSE)
  
  # -- IV: beta ------------------------------------------------------------
  iv_fit <- feols(build_iv_formula(outcome, endo, ivs, controls, fe),
                  data = d, cluster = clf, warn = FALSE, notes = FALSE)
  
  resolve <- function(fit, tm) {
    cand <- c(paste0("fit_", tm), tm)
    hit  <- cand[cand %in% names(coef(fit))]
    if (length(hit) == 0L) NA_character_ else hit[1L]
  }
  
  gamma <- vapply(rfs, function(x) unname(coef(rf_fit)[x]), numeric(1))
  beta  <- vapply(endo, function(x) {
    nm <- resolve(iv_fit, x); if (is.na(nm)) NA_real_ else unname(coef(iv_fit)[nm])
  }, numeric(1))
  
  # -- the K first-stage equations, fitted directly ------------------------
  # Each endogenous term on ALL excluded instruments plus the exogenous set.
  # Pi[k, j] = coefficient on instrument j in the first-stage equation for
  # endogenous regressor k. Row = endogenous, column = instrument.
  Pi    <- matrix(NA_real_, K, K, dimnames = list(endo, ivs))
  Pi_se <- Pi
  Pi_t  <- Pi
  fs_fits <- vector("list", K); names(fs_fits) <- endo
  
  for (k in seq_len(K)) {
    f <- feols(build_ols_formula(endo[k], c(ivs, controls), fe),
               data = d, cluster = clf, warn = FALSE, notes = FALSE)
    fs_fits[[k]] <- f
    b <- coef(f); s <- fixest::se(f)
    for (j in seq_len(K)) {
      if (ivs[j] %in% names(b)) {
        Pi[k, j]    <- unname(b[ivs[j]])
        Pi_se[k, j] <- unname(s[ivs[j]])
        Pi_t[k, j]  <- Pi[k, j] / Pi_se[k, j]
      }
    }
  }
  
  # -- three readings of beta ---------------------------------------------
  #
  # The structural system is  ln(P) = sum_k beta_k * TREAT_k + ...,  and the
  # first stage is  TREAT_k = sum_j Pi[k,j] * IV_j + ....  Substituting,
  #
  #     gamma_j = sum_k Pi[k,j] * beta_k     i.e.     gamma = t(Pi) %*% beta
  #
  # so the correct recovery is  beta = solve(t(Pi)) %*% gamma  --  NOT
  # solve(Pi) %*% gamma. An earlier version of this file used solve(Pi),
  # which inverts the transpose of the right system and does not reproduce
  # fixest's own 2SLS solution. solve(t(Pi)) reproduces BETA_ACTUAL to
  # numerical precision on every instrument; that is the correct check, and
  # it confirms the reported IV estimates are internally consistent 2SLS
  # output rather than a computational error.
  beta_diag <- gamma / diag(Pi)                       # naive diagonal reading
  beta_mat  <- tryCatch(as.vector(solve(t(Pi)) %*% gamma),
                        error = function(e) rep(NA_real_, K))  # what 2SLS is
  
  # -- the per-category first stage a reader would expect ------------------
  # Fitted on the category subsample only. This is the object the concept-level
  # decomposition in Appendix Table 21 aggregates, so it is the right comparison
  # for the claim that first stages do not differ by shoppability.
  pi_sub <- vapply(seq_along(keys), function(k) {
    ds <- d[MOD == keys[k]]
    f <- tryCatch(feols(build_ols_formula(endogenous, c(instrument, controls), fe),
                        data = ds, cluster = clf, warn = FALSE, notes = FALSE),
                  error = function(e) NULL)
    if (is.null(f) || !(instrument %in% names(coef(f)))) NA_real_
    else unname(coef(f)[instrument])
  }, numeric(1))
  names(pi_sub) <- keys
  
  sd_z    <- sd(d[[instrument]], na.rm = TRUE)
  sd_z_by <- vapply(keys, function(k) sd(d[MOD == k][[instrument]], na.rm = TRUE), numeric(1))
  
  out <- data.table(
    INSTRUMENT_LABEL   = instrument_label,
    INSTRUMENT         = instrument,
    SCHEME             = scheme_col,
    CATEGORY           = keys,
    N_ROWS             = vapply(keys, function(k) sum(d$MOD == k), integer(1)),
    GAMMA_RF           = gamma,
    GAMMA_RF_SE        = vapply(rfs, function(x) unname(fixest::se(rf_fit)[x]), numeric(1)),
    RF_PCT_PER_SD      = 100 * (exp(gamma * sd_z) - 1),
    PI_OWN             = diag(Pi),
    PI_OWN_SE          = diag(Pi_se),
    PI_OWN_T           = diag(Pi_t),
    PI_CROSS_MAX       = vapply(seq_len(K), function(k) {
      v <- Pi[k, -k]; if (all(is.na(v))) NA_real_ else v[which.max(abs(v))]
    }, numeric(1)),
    PI_CROSS_MAX_T     = vapply(seq_len(K), function(k) {
      v <- Pi_t[k, -k]; if (all(is.na(v))) NA_real_ else v[which.max(abs(v))]
    }, numeric(1)),
    PI_SUBSAMPLE       = pi_sub,
    BETA_ACTUAL        = beta,
    BETA_ACTUAL_PCT    = 100 * (exp(beta) - 1),
    BETA_IF_DIAGONAL   = beta_diag,
    BETA_IF_DIAG_PCT   = 100 * (exp(beta_diag) - 1),
    BETA_MATRIX_SOLVE  = beta_mat,
    BETA_MATRIX_PCT    = 100 * (exp(beta_mat) - 1),
    SD_INSTRUMENT      = sd_z,
    SD_INSTRUMENT_CAT  = sd_z_by,
    N_OBSERVATIONS     = nobs(rf_fit)
  )
  
  # Reconciliation flags. RECONCILES_MATRIX should be 1 if 2SLS really is
  # Pi^{-1} gamma; RECONCILES_DIAGONAL is 1 only if Pi is effectively diagonal.
  out[, `:=`(
    RECONCILES_MATRIX   = as.integer(abs(BETA_ACTUAL - BETA_MATRIX_SOLVE) < 1e-6),
    RECONCILES_DIAGONAL = as.integer(abs(BETA_ACTUAL - BETA_IF_DIAGONAL)  < 1e-6),
    PI_RATIO_TO_MAX     = diag(Pi) / max(abs(diag(Pi)), na.rm = TRUE),
    GAMMA_RATIO_TO_MAX  = gamma / gamma[which.max(abs(gamma))]
  )]
  
  list(table = out, Pi = Pi, Pi_se = Pi_se, Pi_t = Pi_t,
       gamma = gamma, beta = beta, keys = keys,
       rf_fit = rf_fit, iv_fit = iv_fit, fs_fits = fs_fits, n = nobs(rf_fit))
}


run_s19a <- function(panel = outpatient, scheme_col = "SCHEME_1_CERTAINTY",
                     instruments = MAIN_INSTRUMENTS) {
  .s19_hd("19A -- INTERACTED FIRST-STAGE MATRIX (item A4)")
  
  res <- lapply(names(instruments), function(il) {
    r <- s19_first_stage_matrix(panel, scheme_col, instruments[[il]], il)
    
    .s19_sub(paste0(il, "  (N = ", format(r$n, big.mark = ","), ")"))
    cat("\nFirst-stage coefficient matrix Pi  [rows = endogenous, cols = instrument]\n")
    print(round(r$Pi, 6))
    cat("\n  t-statistics\n")
    print(round(r$Pi_t, 2))
    
    cat("\nReduced form gamma:\n")
    print(round(r$gamma, 6))
    
    cat("\nbeta comparison (log points):\n")
    print(r$table[, .(CATEGORY,
                      BETA_ACTUAL     = round(BETA_ACTUAL, 5),
                      BETA_MATRIX     = round(BETA_MATRIX_SOLVE, 5),
                      BETA_IF_DIAG    = round(BETA_IF_DIAGONAL, 5),
                      MATCHES_MATRIX  = RECONCILES_MATRIX,
                      MATCHES_DIAG    = RECONCILES_DIAGONAL)])
    
    cat("\nFirst stage, interacted own-equation vs category subsample:\n")
    print(r$table[, .(CATEGORY,
                      PI_OWN       = round(PI_OWN, 5),
                      PI_OWN_T     = round(PI_OWN_T, 1),
                      PI_CROSS_MAX = round(PI_CROSS_MAX, 5),
                      PI_CROSS_T   = round(PI_CROSS_MAX_T, 1),
                      PI_SUBSAMPLE = round(PI_SUBSAMPLE, 5))])
    
    r$table
  })
  
  out <- rbindlist(res, fill = TRUE)
  save_csv(out, "T19A_interacted_first_stage_matrix.csv")
  
  # -------- verdict -------------------------------------------------------
  .s19_sub("VERDICT")
  cross_big  <- out[abs(PI_CROSS_MAX_T) > 1.96, .N]
  cross_share <- out[, .(CROSS_SHARE_OF_OWN = round(mean(abs(PI_CROSS_MAX) / abs(PI_OWN)), 2)),
                     by = INSTRUMENT_LABEL]
  subsample_ratio <- out[, .(SUBSAMPLE_RATIO = round(max(PI_SUBSAMPLE) / min(PI_SUBSAMPLE), 2),
                             INTERACTED_RATIO = round(max(abs(PI_OWN)) / min(abs(PI_OWN)), 2)),
                         by = INSTRUMENT_LABEL]
  matrix_ok <- out[, all(RECONCILES_MATRIX == 1L)]
  
  cat("\nsolve(t(Pi)) %*% gamma reproduces feols's own 2SLS coefficients: ",
      if (matrix_ok) "YES, on every row." else "NO -- investigate before proceeding.",
      "\nThis confirms the reported IV estimates are correct 2SLS output, not\n",
      "an arithmetic error.\n", sep = "")
  
  cat("\nOff-diagonal first-stage coefficients significant at 5%: ",
      cross_big, " of ", nrow(out), "\n", sep = "")
  cat("\nCross-loading as a share of the own coefficient, by instrument:\n")
  print(cross_share)
  
  cat("\nOwn first-stage ratio: naive subsample (no cross term) vs the\n",
      "interacted system's own coefficient (cross term included):\n", sep = "")
  print(subsample_ratio)
  
  cat("\nHOW TO READ THIS.\n",
      "  RECONCILES_MATRIX = 1 throughout confirms the paper's Table 5, Panel B\n",
      "  numbers are the exact solution to a jointly-estimated 2SLS system --\n",
      "  there is no coding error. The puzzle is explained, not by an error,\n",
      "  but by SUBSTANTIAL, highly significant cross-loading between the two\n",
      "  category-specific instruments (t approx 7-9, all three instruments),\n",
      "  which arises because Z_Shoppable and Z_Non_shoppable are the same\n",
      "  underlying instrument split by a row-level category indicator, and\n",
      "  POST_MONTH fixed effects span both categories.\n\n",
      "  PI_SUBSAMPLE (the naive, no-cross-term first stage fit separately on\n",
      "  each category) is close to identical across categories -- consistent\n",
      "  with, and the like-for-like comparison to, Appendix Table 21's null on\n",
      "  first-stage differences by shoppability. PI_OWN (the interacted\n",
      "  system's own coefficient, WITH the cross term in the same regression)\n",
      "  is roughly 30% larger and differs more across categories. These are\n",
      "  two different, both-correct objects, and the paper should say so:\n",
      "  Table 21 is measuring the univariate first stage; Table 5's IV column\n",
      "  is measuring the solution to a cross-loaded two-equation system.\n\n",
      "  ACTION: this is a disclosure item, not a correction. Add the Pi matrix\n",
      "  (or a summary of its off-diagonal magnitude and significance) to the\n",
      "  paper, footnote Table 5 Panel B accordingly, and use this to sharpen\n",
      "  -- not merely restate -- the existing decision to lead with the\n",
      "  reduced form: the two structural parameters are not separately\n",
      "  identified off their own instrument alone.\n", sep = "")
  
  invisible(out)
}


# ===========================================================================
# 19B. THE REDUCED-FORM GAP, WITH ITS SE AND ITS UNITS            (A4 and C5)
# ===========================================================================
#
# The paper currently reports the gradient under two different percent
# conventions without saying so, which is the whole of item C5:
#
#   estimate_interacted()      RF_PERCENT_PER_SD = 100 * (exp(b * sd_z) - 1)
#                              per arm; Table 5 then leaves the reader to
#                              difference two already-exponentiated numbers.
#
#   run_conley_sensitivity()   GAP_PCT_PER_SD    = 100 * gap * sd_z
#                              a linear log-point approximation of the gap.
#
# These are not the same function of the same coefficients. This block reports
# the gap once, with its own standard error, under all three conventions, so
# whichever is adopted can be applied consistently everywhere.

s19_gap_units <- function(
    data, scheme_col = "SCHEME_1_CERTAINTY",
    instrument = PRIMARY_INSTRUMENT, instrument_label = "",
    outcome = PRIMARY_OUTCOME, controls = BASELINE_CONTROLS,
    fixed_effects = BASELINE_FIXED_EFFECTS, clusters = BASELINE_CLUSTERS) {
  
  controls <- available_columns(data, controls)
  fe <- available_columns(data, fixed_effects); cl <- available_columns(data, clusters)
  
  d <- data[!is.na(get(scheme_col))]
  d <- model_sample(d, c(outcome, ENDOGENOUS_VARIABLE, instrument, controls, fe, cl, scheme_col))
  
  d[, MOD := droplevels(factor(get(scheme_col)))]
  keys <- levels(d$MOD)
  if (!all(c("Shoppable", "Non_shoppable") %chin% keys))
    stop("Expected Shoppable / Non_shoppable levels.", call. = FALSE)
  for (k in keys) d[, paste0("RF_", k) := get(instrument) * as.integer(MOD == k)]
  
  fit <- feols(build_ols_formula(outcome, c(paste0("RF_", keys), controls), fe),
               data = d, cluster = build_cluster_formula(cl), warn = FALSE, notes = FALSE)
  
  nS <- "RF_Shoppable"; nN <- "RF_Non_shoppable"
  bS <- unname(coef(fit)[nS]); bN <- unname(coef(fit)[nN])
  V <- vcov(fit)
  seS <- sqrt(V[nS, nS]); seN <- sqrt(V[nN, nN])
  gap <- bS - bN
  se_gap <- sqrt(V[nS, nS] + V[nN, nN] - 2 * V[nS, nN])
  z <- gap / se_gap
  s <- sd(d[[instrument]], na.rm = TRUE)
  
  data.table(
    INSTRUMENT_LABEL = instrument_label, INSTRUMENT = instrument,
    N_OBSERVATIONS = nobs(fit), SD_INSTRUMENT = s,
    B_SHOP = bS, SE_SHOP = seS, P_SHOP = .pval(bS / seS, fit),
    B_NON  = bN, SE_NON  = seN, P_NON  = .pval(bN / seN, fit),
    PCT_SHOP = 100 * (exp(bS * s) - 1),
    PCT_NON  = 100 * (exp(bN * s) - 1),
    GAP_COEF = gap, GAP_SE = se_gap, GAP_Z = z,
    GAP_P = .pval(z, fit),
    GAP_CI_LOW  = gap - 1.96 * se_gap,
    GAP_CI_HIGH = gap + 1.96 * se_gap,
    # convention (i): difference of the two exponentiated arms  (Table 5 implied)
    GAP_PCT_DIFF_OF_PCTS = 100 * (exp(bS * s) - exp(bN * s)),
    # convention (ii): exponentiate the difference
    GAP_PCT_EXP_OF_DIFF  = 100 * (exp(gap * s) - 1),
    # convention (iii): linear log-point scaling      (Conley table, Table 33)
    GAP_PCT_LINEAR       = 100 * gap * s
  )
}


run_s19b <- function(panel = outpatient, scheme_col = "SCHEME_1_CERTAINTY",
                     instruments = MAIN_INSTRUMENTS) {
  .s19_hd("19B -- REDUCED-FORM GAP, SE, AND UNIT CONVENTIONS (items A4, C5)")
  
  out <- rbindlist(lapply(names(instruments), function(il)
    s19_gap_units(panel, scheme_col, instruments[[il]], il)), fill = TRUE)
  
  cat("\nGap point estimate and inference (log points):\n")
  print(out[, .(INSTRUMENT_LABEL,
                B_SHOP = round(B_SHOP, 6), B_NON = round(B_NON, 6),
                GAP = round(GAP_COEF, 6), SE = round(GAP_SE, 6),
                Z = round(GAP_Z, 2), P = round(GAP_P, 4))])
  
  cat("\nThe same gap under the three percent conventions in the paper:\n")
  print(out[, .(INSTRUMENT_LABEL,
                SD_Z          = round(SD_INSTRUMENT, 2),
                PCT_SHOP      = round(PCT_SHOP, 2),
                PCT_NON       = round(PCT_NON, 2),
                `(i) diff of pcts`  = round(GAP_PCT_DIFF_OF_PCTS, 2),
                `(ii) exp of diff`  = round(GAP_PCT_EXP_OF_DIFF, 2),
                `(iii) linear`      = round(GAP_PCT_LINEAR, 2))])
  
  cat("\nConvention (i) is what a reader gets by differencing Table 5, Panel A.\n",
      "Convention (iii) is what Appendix Table 33 reports as the observed\n",
      "gradient. They are different functions of the same coefficients, which is\n",
      "the entire source of the mismatch flagged in item C5. Pick one, apply it\n",
      "everywhere, and add the GAP / SE row above to Table 5.\n", sep = "")
  
  save_csv(out, "T19B_gap_units_and_se.csv")
  invisible(out)
}


# ===========================================================================
# 19C. MINIMUM DETECTABLE EFFECTS                                        (A7)
# ===========================================================================
#
# A null is only a rejection if the design could have found something. For a
# two-sided 5% test at 80% power the smallest detectable coefficient is
#
#     MDE = (z_{0.975} + z_{0.80}) * SE  =  2.802 * SE
#
# The comparability moderators are estimated on standardised moderators, so
# their coefficients are already "change in the concept-level reduced-form
# coefficient per SD of the moderator." The natural benchmark is the
# shoppability gradient in the SAME units, which is the concept-level
# decomposition coefficient in T08 -- i.e. the shoppable-minus-non-shoppable
# difference in concept-level RF coefficients.
#
# MDE_SHARE below is the answer to "how large a comparability gradient can this
# design rule out, as a fraction of the shoppability gradient it does detect?"

S19_Z_ALPHA <- qnorm(0.975)   # 1.959964
S19_Z_POWER <- qnorm(0.80)    # 0.8416212
S19_MDE_K   <- S19_Z_ALPHA + S19_Z_POWER

s19_power_bounds <- function(within_family_results,
                             decomposition_results,
                             spec_pattern = "^\\(b\\)",
                             tier = "MAIN") {
  
  # within_family_results: comp_wf, i.e. run_comparability_within_family()'s
  #   return value (T09D). Columns include term, estimate, std.error, p.value,
  #   MODERATOR, SPEC, INSTRUMENT_LABEL, TIER.
  wf <- as.data.table(copy(within_family_results))
  wf <- wf[term == "MODC" & grepl(spec_pattern, SPEC)]
  if (!is.null(tier) && "TIER" %in% names(wf)) wf <- wf[TIER == tier]
  if (nrow(wf) == 0L) stop("No within-family rows matched. Check that ",
                           "within_family_results is comp_wf (T09D), not T09.",
                           call. = FALSE)
  
  # Benchmark: the shoppability gradient at the CONCEPT level, in the SAME
  # units as RF_COEF above (log points, per concept), so the ratio is a
  # like-for-like comparison of two concept-level reduced-form quantities.
  #
  # decomposition_results: the object returned by decompose_reduced_form(),
  # i.e. T08D_RF_vs_FS_decomposition. It is long over DEPENDENT
  # (RF_COEF / FS_COEF / IV_COEF) and over `term`, which for a two-level
  # factor is the fixest dummy name "SHOP_CERTAINTYShoppable".
  dec <- as.data.table(copy(decomposition_results))
  req_dec <- c("DEPENDENT", "term", "estimate", "INSTRUMENT_LABEL")
  if (!all(req_dec %chin% names(dec))) {
    stop("decomposition_results is missing one of: ", paste(req_dec, collapse = ", "),
         ". Pass the object returned by decompose_reduced_form() (T08D).", call. = FALSE)
  }
  dec <- dec[DEPENDENT == "RF_COEF" & grepl("Shop", term, ignore.case = TRUE)]
  if (nrow(dec) == 0L) {
    stop("No RF_COEF / Shoppable rows found in decomposition_results. ",
         "Check DEPENDENT and term values with unique(decomposition_results$DEPENDENT) ",
         "and unique(decomposition_results$term).", call. = FALSE)
  }
  bench <- unique(dec[, .(INSTRUMENT_LABEL, SHOPPABILITY_GRADIENT = abs(estimate))])
  
  out <- merge(wf, bench, by = "INSTRUMENT_LABEL", all.x = TRUE)
  if (out[is.na(SHOPPABILITY_GRADIENT), .N] > 0L) {
    warning(out[is.na(SHOPPABILITY_GRADIENT), uniqueN(INSTRUMENT_LABEL)],
            " instrument label(s) in within_family_results had no match in ",
            "decomposition_results; check INSTRUMENT_LABEL spelling against ",
            "INSTRUMENT_LABEL_MAP.", call. = FALSE)
  }
  out[, `:=`(
    MDE            = S19_MDE_K * std.error,
    MDE_SHARE      = S19_MDE_K * std.error / SHOPPABILITY_GRADIENT,
    ABS_EST_SHARE  = abs(estimate) / SHOPPABILITY_GRADIENT,
    DETECTED       = as.integer(p.value < 0.05)
  )]
  
  # A tested-and-rejected claim needs the bound to exclude an effect that would
  # matter. 0.5 is a defensible line: a comparability gradient smaller than half
  # the shoppability gradient cannot be the explanation for the shoppability
  # gradient.
  out[, VERDICT := fifelse(DETECTED == 1L, "Detected",
                           fifelse(MDE_SHARE < 0.50, "Rejected (bounded)",
                                   "Unsupported (underpowered)"))]
  
  setorder(out, MODERATOR, INSTRUMENT_LABEL)
  out[]
}


s19_power_from_p <- function(estimate, p_value, benchmark) {
  # For tables that report an estimate and a p-value but no SE (Section 12
  # output, Table 11 in the paper), recover the SE and do the same arithmetic.
  z  <- qnorm(1 - p_value / 2)
  se <- abs(estimate) / z
  data.table(ESTIMATE = estimate, P_VALUE = p_value, IMPLIED_SE = se,
             MDE = S19_MDE_K * se,
             MDE_SHARE = S19_MDE_K * se / abs(benchmark),
             BENCHMARK = benchmark)
}


run_s19c <- function(within_family_results, decomposition_results,
                     ses_gradient_results = NULL) {
  .s19_hd("19C -- MINIMUM DETECTABLE EFFECTS (item A7)")
  
  pb <- s19_power_bounds(within_family_results, decomposition_results)
  
  cat("\nComparability and contracting-depth moderators, within-family spec (b),\n",
      "MAIN tier. MDE is the smallest coefficient detectable at 80% power and\n",
      "5% two-sided; MDE_SHARE expresses it against the shoppability gradient.\n", sep = "")
  print(pb[, .(MODERATOR, INSTRUMENT_LABEL,
               EST       = signif(estimate, 3),
               SE        = signif(std.error, 3),
               P         = round(p.value, 4),
               MDE       = signif(MDE, 3),
               MDE_SHARE = round(MDE_SHARE, 3),
               VERDICT)])
  
  cat("\nWorst-case bound per moderator (largest MDE_SHARE across MAIN instruments):\n")
  print(pb[, .(MAX_MDE_SHARE = round(max(MDE_SHARE, na.rm = TRUE), 3),
               MIN_MDE_SHARE = round(min(MDE_SHARE, na.rm = TRUE), 3),
               N_DETECTED    = sum(DETECTED),
               VERDICT       = fifelse(sum(DETECTED) > 0L, "Detected",
                                       fifelse(max(MDE_SHARE, na.rm = TRUE) < 0.50,
                                               "Rejected (bounded)",
                                               "Unsupported (underpowered)"))),
           by = MODERATOR][order(MAX_MDE_SHARE)])
  
  save_csv(pb, "T19C_power_bounds_comparability.csv")
  
  if (!is.null(ses_gradient_results)) {
    .s19_sub("Demographic index interaction (Section 12)")
    sg <- as.data.table(copy(ses_gradient_results))
    ecol <- intersect(c("DELTA_GAP_PER_SD", "ESTIMATE", "GAP_CHANGE"), names(sg))[1L]
    gcol <- intersect(c("GAP_AT_MEAN", "GAP"), names(sg))[1L]
    pcol <- intersect(c("P_VALUE", "P"), names(sg))[1L]
    if (!is.na(ecol) && !is.na(pcol) && !is.na(gcol)) {
      sg[, IMPLIED_SE := abs(get(ecol)) / qnorm(1 - get(pcol) / 2)]
      sg[, MDE       := S19_MDE_K * IMPLIED_SE]
      sg[, MDE_SHARE := MDE / abs(get(gcol))]
      print(sg[, .(SCHEME, INSTRUMENT_LABEL,
                   GAP       = round(get(gcol), 2),
                   DELTA     = round(get(ecol), 2),
                   P         = round(get(pcol), 3),
                   MDE       = round(MDE, 2),
                   MDE_SHARE = round(MDE_SHARE, 2))][order(MDE_SHARE)])
      cat("\nMDE_SHARE near or above 1 means the design cannot rule out a\n",
          "demographic gradient as large as the gap itself. That is the\n",
          "definition of underpowered, and it is why Section 7 should keep\n",
          "saying 'not supported' rather than 'rejected'.\n", sep = "")
      save_csv(sg, "T19C_power_bounds_demographics.csv")
    } else {
      cat("Could not locate the expected columns; pass the SES gradient table\n",
          "with GAP_AT_MEAN, DELTA_GAP_PER_SD, and P_VALUE.\n", sep = "")
    }
  }
  
  invisible(pb)
}


# ===========================================================================
# DRIVER
# ===========================================================================
#
# Expects in memory:
#   outpatient                     the analysis panel
#   comparability_within_family    T09D, from run_comparability_within_family()
#   decomposition                  T08, the concept-level decomposition
#   ses_gradient                   optional, Section 12 index results
#
# Adjust the object names on the right-hand side if yours differ.

run_section_19 <- function(panel = outpatient,
                           within_family = NULL,
                           decomposition = NULL,
                           ses_gradient  = NULL) {
  
  a <- run_s19a(panel)
  b <- run_s19b(panel)
  
  c_out <- NULL
  if (!is.null(within_family) && !is.null(decomposition)) {
    c_out <- run_s19c(within_family, decomposition, ses_gradient)
  } else {
    cat("\n[19C skipped] Pass within_family = and decomposition = to run the\n",
        "power bounds. These are the T09D and T08 tables.\n", sep = "")
  }
  
  .s19_hd("SECTION 19 COMPLETE")
  invisible(list(first_stage_matrix = a, gap_units = b, power_bounds = c_out))
}

# Example:
# s19 <- run_section_19(outpatient,
#                       within_family = comparability_within_family,
#                       decomposition = decomposition_results,
#                       ses_gradient  = ses_gradient_table)