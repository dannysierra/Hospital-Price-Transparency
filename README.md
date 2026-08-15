# When Transparency Works: Service Shoppability, Contracting Depth, and the Price Effects of Hospital Disclosure

**Danny Sierra** · PhD Candidate, Department of Economics, Florida State University · [Ds22c@fsu.edu](mailto:Ds22c@fsu.edu)

---

## Overview

This repository holds the estimation code for my job market paper, which asks when the 2021 CMS Hospital Price Transparency Rule moves negotiated prices between hospitals and commercial insurers — and when it does not.

The paper does not test for a uniform price effect. It tests for a *conditional* one. Disclosure should discipline a negotiated rate only where the posted number is a usable basis for comparison, so the empirical strategy identifies the service and market conditions under which peer disclosure translates into pressure on price.

**The findings, in the order the paper reaches them:**

1. **The average effect is uninformative, not zero.** The pooled response is imprecise and cannot be signed. This is the first finding rather than a failed test: the average mixes services where the mechanism can operate with services where it cannot, so the average is not a quantity worth interpreting.

2. **The shoppability gradient is the result.** Shoppable services show economically meaningful price declines; non-shoppable services show a precise zero. The gradient is negative across every instrument and every classification scheme tested, and the sign does not depend on how the shoppable line is drawn.

3. **Price comparability, tested directly, is rejected as the mechanism.** Three price-dispersion measures are null once clinical family fixed effects absorb modality. The one measure that survives — contracting depth, the number of distinct payers per concept — is significant with the **opposite** sign to the ex ante prediction, which points to negotiation entrenchment rather than to comparability.

4. **Demographic heterogeneity runs against the pre-registered prediction.** The consumer shopping-capacity model predicts larger responses in more advantaged markets. The estimates do not support that, and the paper reports the rejection rather than reframing around it. See [Open items](#open-items).

The repository preserves the honest sequence: hypothesis stated, tested, rejected, alternative found. Sections that test and reject a mechanism are kept in the code rather than deleted.

---

## Paper

> **"When Transparency Works: Service Shoppability, Contracting Depth, and the Price Effects of Hospital Disclosure."**
> Danny Sierra, Florida State University, 2026. *Job market paper. Draft available on request.*

---

## Data

The negotiated-rate data are licensed and are **not** redistributed here. The code is published for methodological transparency and to document every specification decision; it is not runnable end to end without data access.

| Source | Contribution | Access |
|---|---|---|
| **Turquoise Health** | Negotiated rates parsed from hospital machine-readable files | Licensed |
| **AHA Annual Survey** | Bed counts, system membership, ownership, hospital type | Licensed |
| **CMS enforcement actions** | Warning and closure notices (exploratory instruments only) | Public |
| **American Community Survey (2022, 5-year)** | County demographics, pulled in-script via `tidycensus` | Public |
| **Census county population estimates** | Coverage denominators | Public |

**Panel.** The analysis panel is at the **hospital × clinical concept × post-month** level: roughly 1.4 million observations covering about 3,700 hospitals and 738 clinical concepts across 16 clinical families. A "concept" is a clinically coherent service defined in the project codebook, sitting between the individual billing code and the broad modality.

Note that roughly 31% of complete cases are dropped by cascading singleton removal across the county × concept and month fixed effects. This is expected, is identical across specifications, and is reported explicitly rather than left implicit — `audit_estimation_sample()` writes `QA01`, and the coverage block at the end of the script maps which counties actually contribute to identification versus merely appearing in the data.

---

## Empirical strategy

### Specification

Estimated with `fixest::feols`:

```
ln(P_ihct) = β · N_PRIOR_POSTERS_ht + γ · ln(Beds_h) + δ_mc + τ_t + ε_ihct
```

with `δ_mc` a county × concept market fixed effect (`MARKET_ID`), `τ_t` a post-month fixed effect, and standard errors two-way clustered on county and post-month. `N_PRIOR_POSTERS` counts hospitals in the same market that have already disclosed, and is instrumented by out-of-market disclosure rollout among competitor health systems.

### Four design decisions, fixed throughout

These are set once in the script header and are not re-litigated section by section.

**1. Treatment enters linearly.** Every concave transform of `N_PRIOR_POSTERS` inflates the IV coefficient by shrinking `Cov(T, Z)` in the Wald denominator, monotonically in how much it compresses the treatment. The reduced-form numerator does not move. The linear first stage runs F ≈ 43 in both interacted equations, so there is no weak-instrument problem for a transform to solve. The transform ladder is reported as robustness and is never a headline.

**2. The reduced form is the primary estimator.** IV is a ratio; the reduced form is its numerator, carrying no treatment variable, no functional-form choice, and no first-stage noise in the standard error. Because every model here is **exactly identified**, the Anderson–Rubin weak-instrument-robust test is numerically identical to the t-test on the reduced-form coefficient — so the reported RF p-values already *are* robust p-values. LIML offers nothing, since it coincides with 2SLS under exact identification.

**3. Heterogeneity is estimated by interaction, never by sample split.** The fixed effects absorb most of the instrument's effect on the treatment within any single subsample, so splitting destroys identification. A single function, `estimate_interacted()`, produces every heterogeneity result in the paper, in both categorical and continuous-moderator forms.

**4. No event study is possible.** All but two hospitals appear at exactly one posting month, so there is no within-hospital timing variation. The selection concern is addressed structurally instead: the non-shoppable group is a **within-hospital placebo** — same hospital, same posting month, same instrument value, no price response. Any selection story has to explain why selection would operate on shoppable services and not on non-shoppable ones at the same hospital in the same month.

### Instruments

Six competitor-rollout instruments are graded on how credible their exclusion restriction is and how stable their sign is across classification schemes. **The grading is by construction and by audit, never by first-stage F** — Section 7.6 shows the two strongest first stages produce the two weakest results, which is precisely why instrument selection here is a design decision rather than an outcome-based one.

| Tier | Instrument | Role |
|---|---|---|
| **MAIN** | `Competitor_only_hospitals_9m` *(primary)* | Carries every headline claim |
| **MAIN** | `Primary_strict_system_IV` | Headline |
| **MAIN** | `Competitor_outside_CBSA_hospitals_9m` | Headline |
| **CONFIRMING** | `Competitor_outside_CBSA_counties_9m` | Pooled robustness alongside MAIN |
| **CONFIRMING** | `Competitor_systems_9m` | Pooled robustness — see [Open items](#open-items) |
| **DISCREPANT** | `Competitor_outside_CBSA_systems_9m` | Reported alone, never pooled |

Only the `Z_SYS_COMPETITOR_*` family is hospital-specific and excludes the focal hospital's own system; the remaining `Z_SYS_*` measures are county-month local exposure and do not. That distinction is what the exclusion-restriction defence rests on, and it is recorded in the `EXCLUSION_NOTE` column of `T02`.

Additional instruments are carried for robustness only: a construction ladder already present in the panel (Section 13C), a 3/6/9/12-month window ladder rebuilt from source (Section 15B), and CMS enforcement measures that remain **exploratory** — enforcement plausibly moves prices directly through compliance cost and correlated conduct, so it fails the exclusion restriction on its face and no reparameterisation changes that.

### Shoppability classification

Eighteen classification schemes are built from the codebook; six are attached to the panel and carry the results. A scheme is a partition into HIGH / INTERMEDIATE / LOW, assembled from named terms, where a term is a clinical family plus an optional keyword regex on the official code description plus an optional exclusion regex.

Three rules govern the construction:

- Emergency, critical care, and everything billed as MS-DRG is forced LOW in **every** scheme. `build_schemes()` asserts this rather than assuming it.
- Residual `"_Other"` terms mean whatever is left in a family after the named sub-terms are claimed *within the same scheme*, so `"CT Other"` is scheme-dependent by design.
- Inclusion keywords are substring matches against abbreviated CPT text and over-match. `CHEST` appears inside `"Ct chest spine w/o dye"` (a thoracic spine study) and `ULTRASOUND` inside `"Colonoscopy w/ultrasound"` (a procedure adjunct). `resolve_term()` therefore honours an `exclude` regex, and a guard tests for exactly that case.

The six panel schemes:

| Column | Scheme |
|---|---|
| `SCHEME_1_CERTAINTY` | Procedural certainty (family-based; primary) |
| `SCHEME_2_THEORYV2` | Theory-based V2, MRI non-shoppable |
| `SCHEME_3_IMAGING` | Imaging vs. procedural |
| `SCHEME_4_CMS70` | CMS statutory shoppable list |
| `SCHEME_5_MDSAVE` | Upfront cash-market |
| `SCHEME_6_WITHINMOD` | High vs. low within modality |

Because several schemes induce the *same* partition of the concept universe, `deduplicate_partitions()` fingerprints each scheme–collapse variant on its estimated coefficient vector and keeps one representative per distinct partition. The paper reports the number of genuinely distinct specifications, not the number of labels.

### Inference

Shoppability is assigned to roughly 16 clinical **families**, not to 738 independent concepts. Concept-level clustered standard errors would treat 168 biopsy concepts as 168 independent observations when a single labelling decision covers all of them.

`family_permutation_test()` therefore permutes which families are called shoppable and enumerates every assignment, giving inference at the level the decision was actually made. With 16 families and 10 shoppable, C(16,10) = 8,008 assignments enumerate exactly, so the smallest attainable two-sided p-value is 1/8008. **A result at that floor is reported as being at the minimum attainable value, not as a precise small number.**

---

## Repository structure

```
.
├── .gitignore
├── README.md
└── Code/
    └── HPT_Analysis_Pipeline.R
```

This file is stage 3 of a three-stage pipeline. Stage 1 (Snowflake SQL) assembles disclosures into a cleaned price panel and builds the codebook; stage 2 (Python) constructs hospital, market, and instrument covariates. Only stage 3 is published here.

---

## Script structure

Sections 0–12 define constants and functions and load in seconds. Everything after them executes, gated by `RUN_STAGES`.

### Definitions

| Section | Contents |
|---|---|
| **0** | Paths, baseline specification, instrument registry, scheme registry, sample screens. Everything configurable lives here; nothing below hard-codes a variable name or threshold. |
| **1** | Helpers: safe coercion, coefficient extraction, first-stage diagnostics, formula builders, CSV writers, the `cache_or_run()` gate. |
| **2** | Panel loading and preparation. Every function returns data rather than mutating a global, so `outpatient` is assigned exactly once from a visible call chain. |
| **3** | Shoppability scheme construction, concept merges, scheme inheritance for constructed canonical IDs. |
| **4** | Estimators. `run_reduced_form()`, `run_first_stage()`, `run_ols()`, `run_iv()`, and `estimate_interacted()`. |
| **5** | Instrument screen (`T02`) and pooled models (`T03`). |
| **6** | Concept-level estimation. Stores RF, FS, *and* IV per concept — keeping the first stage separately is what makes the Section 8 decomposition possible. |
| **7** | Main interacted results and the transform ladder. |
| **8** | Meta-regressions, partition deduplication, family permutation inference, RF-vs-FS decomposition. |
| **9** | Price comparability as a candidate mechanism. Eight measures built leave-one-county-out, screened for degeneracy, tested row-level and within-family. |
| **10** | Comparability meta-regression: specification curve, family fixed effects, and a horse race against the shoppability label. |
| **11** | Instrument balance and placebo-outcome tests. |
| **12** | Preflight. Ten seconds; must pass before the eight-hour stage 6. |

### Run blocks and robustness

| Block | Contents |
|---|---|
| **BUILD** | Builds `schemes_long` and the panel, exports the Scheme 1 assignment, audits the estimation sample, runs preflight. |
| **Instrument audit** | Three separate questions answered separately: raw strength on a common sample, sign stability across schemes, and whether the headline result changes under each of the six instruments. Writes `QA05`–`QA08`. |
| **Stages 5–11** | The estimation stages, selected by `RUN_STAGES`. |
| **IN_SYSTEM robustness** | Re-estimates the headline with system membership controlled, closing the balance-test finding by demonstration rather than by argument. |
| **County demographics** | ACS pull via `tidycensus`, with a mandatory spot check — a wrong variable code returns a plausible column of zeros rather than an error. |
| **Section 12** | Demographic heterogeneity in the *gradient*, `Z × Shop × M`, not in the pooled response. |
| **Section 12B** | SES terciles as a categorical moderator, relaxing the linear-gradient assumption. |
| **Section 13** | System × month FE (13A), leave-one-system-out (13B), instrument construction ladder (13C), randomisation inference (13D). |
| **Section 15** | CBSA market definition (15A), instrument window ladder (15B), payer-conditional estimates (15C). |
| **Figures** | Publication figures 1–8 and robustness figures 10–20, read back from the output CSVs so they run standalone. |
| **Coverage** | County and population coverage measured on the *estimation sample*, not the full panel, plus a print-safe CONUS map. |
| **LaTeX builds** | Summary statistics, family-level scheme assignment, and shoppable-share inputs, each emitting a ready-to-`\input` `.tex`. |

---

## Output index

Results are written as CSV under `RESULT_ROOT`, prefixed by the section that produces them.

### Analysis tables

| File | Contents |
|---|---|
| `T02_instrument_first_stage_screen` | First stage and reduced form for every candidate instrument |
| `T03_pooled_OLS_RF_IV_all_outcomes` | Pooled estimates across all five price outcomes |
| `T05_concept_level_RF_FS_IV` | Concept-level RF, FS, and IV — the meta-regression input |
| `T05C_size_gradient_RF_vs_IV` | Size-gradient diagnostic |
| `T06_main_*` | Headline interacted results, MAIN tier |
| `T06D_*` | Scheme 1 procedural-exclusion sensitivity |
| `T06F` / `T06G` | All three tiers combined, estimates and heterogeneity tests |
| `T06H_*` | Pooled MAIN + CONFIRMING |
| `T06J_*` / `T06K_*` | CONFIRMING and DISCREPANT tiers |
| `T06L_tier_rerank_all_schemes` | Tier diagnostic across all six schemes, with first-stage F alongside |
| `T07` / `T07B` | Transform ladder |
| `T08*` | Meta-regressions, redundant partitions, permutation inference, decomposition, gradient by tier |
| `T09*` | Comparability measures, row-level interaction, within-family test |
| `T10` / `T10B` | Comparability meta-regression and horse race |
| `T11*` | Balance, placebo outcome, IN_SYSTEM robustness comparison |
| `T12T_*` / `T12Q_*` | Demographic interaction in the gradient; SES terciles |
| `T13A`–`T13D` | System × month FE, leave-one-system-out, construction ladder, randomisation inference |
| `T15A` / `T15C` | CBSA market definition, payer-conditional |

### Audit tables

| File | Contents |
|---|---|
| `QA01_estimation_sample_audit` | Panel rows → complete cases → estimated rows, with singleton diagnostics |
| `QA02_scheme_classification_crosswalk` | Scheme classification crosswalk |
| `QA03` / `QA04` | Scheme 1 concept assignments; procedural concepts flagged inside the shoppable group |
| `QA05`–`QA08` | Instrument audit: strength, sign stability, headline result under all six |
| `QA09_moderator_screen` | Comparability moderator screen, with the reason for every drop |
| `QA11_shoppable_shares` | Shoppable shares for the magnitude calculation |

### Figures

`fig01_headline_shoppability` · `fig02_scheme_robustness` · `fig03_transform_ladder` · `fig04_concept_distribution` · `fig05_instrument_tiers` · `fig06_decomposition` · `fig07_mechanism` · `fig08_size_gradient` · `fig10_window_ladder` · `fig11_construction` · `fig12_leave_one_out` · `fig13_sysmonth` · `fig14_randinf` · `fig15_payer` · `fig16_cbsa` · `fig17_ses_gradient` · `fig18_ses_terciles` · `fig19_demo_moderators` · `fig20_franchise`

---

## Running the code

Set `PROJECT_ROOT` in Section 0, then:

```r
# Step 1 — load Sections 0-12 (seconds)

# Step 2 — build and validate. run_preflight() must pass before anything else.
#          Run the BUILD block.

# Step 3 — main results and balance tests (~1 hour). Independent of stage 6.
RUN_STAGES <- c(7, 11)

# Step 4 — concept-level sweep (~8 hours)
RUN_STAGES <- c(6)

# Step 5 — meta-regressions and mechanism (~1 hour)
RUN_STAGES <- c(8, 9, 10)
```

Every expensive step is wrapped in `cache_or_run()`, which writes an `.rds` to `CACHE_DIR` and reloads it on later runs. Set `USE_CACHE <- FALSE` to force recomputation, or delete individual `.rds` files to invalidate selectively. **A step producing zero rows never overwrites a good cache** — it stops instead, so a silent upstream failure cannot destroy hours of completed estimation. The stage 6 loader independently verifies that cached concept-level results were estimated on the panel currently in memory.

The demographics block needs a Census API key: `census_api_key("YOUR_KEY", install = TRUE)`, run once per machine.

A few blocks in the file are labelled **interactive checkpoints**. They duplicate what the BUILD block and stage 6 loader do and exist so individual stages can be re-run without stepping through the whole file. Skip them on a cold `source()` of the entire script.

---

## Requirements

R 4.4 or later.

```r
# Core — required
install.packages(c("arrow", "data.table", "fixest", "dplyr",
                   "ggplot2", "stringr", "lubridate", "MASS"))

# Optional — checked for at the point of use
install.packages(c("R.utils",     # gzipped payer-dispersion shards
                   "tidycensus",  # ACS county demographics
                   "ggrepel", "patchwork",
                   "sf", "tigris", "ggpattern", "magick"))
```

Estimation is `fixest` throughout. `data.table` carries the panel. Figures are `ggplot2` with `patchwork`, in FSU garnet (`#782F40`), gold (`#CEB888`), and grey (`#6D6E71`).

---

## Open items

Two labelling decisions are unresolved and are recorded here rather than papered over.

**Instrument tier vocabulary.** The pipeline computes three tiers — MAIN / CONFIRMING / DISCREPANT — via `instrument_tier()`, and the tier-tagged tables (`T06F`, `T06G`, `T08F`) use them. The instrument audit tables (`QA05`–`QA08`) use a two-way MAIN / SUPPORTING split, where SUPPORTING is the union of CONFIRMING and DISCREPANT. Both are live in the code. This README uses the three-tier vocabulary.

Separately, `Competitor_systems_9m` currently sits in CONFIRMING, but its shoppability differential in `T08F` is unstable across scheme variants even though its level coefficient is sign-stable. That is the criterion that defines the DISCREPANT tier, so the assignment is under review.

**Concept count.** The panel carries 738 concepts. The concept-level sweep estimates only those clearing `MIN_SERVICE_OBS`, `MIN_SERVICE_MARKETS`, and `MIN_SERVICE_MONTHS`, so the estimable count reported in `T05` and `T06` differs from the panel count and from each other. The screens are in Section 0 and the per-table counts are in the `N_CONCEPTS` column.

**Demographic heterogeneity.** Section 12 estimates demographic variation in the shoppability *gradient*. An earlier specification interacted the instrument with demographics alone, which measures variation in the pooled response — a quantity the paper disowns, since it averages a real effect on shoppable services against a precise zero on non-shoppable ones. Any table generated under the earlier specification should not be read as evidence about the gradient.

---

## Citation

> Sierra, Danny. "When Transparency Works: Service Shoppability, Contracting Depth, and the Price Effects of Hospital Disclosure." Working Paper, 2026.