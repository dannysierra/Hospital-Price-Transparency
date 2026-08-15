# When Transparency Works: Service Shoppability, Contracting Depth, and the Price Effects of Hospital Disclosure

**Danny Sierra** · PhD Candidate, Department of Economics, Florida State University · [Ds22c@fsu.edu](mailto:Ds22c@fsu.edu)

---

## Overview

This repository holds the replication code for my job market paper, which asks when the 2021 CMS Hospital Price Transparency Rule moves negotiated prices between hospitals and commercial insurers — and when it does not.

The paper does not test for a uniform price effect. It tests for a *conditional* one. Disclosure should discipline a negotiated rate only where the posted number is a usable basis for comparison, so the empirical strategy identifies the service and market conditions under which peer disclosure translates into pressure on price.

**The findings, in the order the paper reaches them:**

1. **The average effect is uninformative, not zero.** The pooled response is imprecise and cannot be signed. This is the first finding rather than a failed test: the average mixes services where the mechanism can operate with services where it cannot, so the average is not a quantity worth interpreting.

2. **The shoppability gradient is the result.** Shoppable services show economically meaningful price declines; non-shoppable services show a precise zero. The gradient is negative across every instrument and every classification scheme tested, and the sign does not depend on how the shoppable line is drawn.

3. **Price comparability, tested directly, is rejected as the mechanism.** Three price-dispersion measures are null once clinical family fixed effects absorb modality. The one measure that survives — contracting depth, the number of distinct payers per concept — is significant with the **opposite** sign to the ex ante prediction, which points to negotiation entrenchment rather than to comparability.

4. **Demographic heterogeneity runs against the pre-registered prediction.** The consumer shopping-capacity model predicts larger responses in more advantaged markets. The estimates do not support that, and the paper reports the rejection rather than reframing around it. See [Open items](#open-items).

The repository preserves the honest sequence: hypothesis stated, tested, rejected, alternative found. Stages that test and reject a mechanism are kept in the code rather than deleted.

---

## Paper

> **"When Transparency Works: Service Shoppability, Contracting Depth, and the Price Effects of Hospital Disclosure."**
> Danny Sierra, Florida State University, 2026. *Job market paper. Draft available on request.*

---

## The pipeline

Three stages, run in order. Each writes files the next reads; nothing runs backwards.

```
   SQL (Snowflake)          Python (Jupyter)              R
   ───────────────          ────────────────              ─
   raw disclosures    →     Snowflake exports       →     analysis panels
   clinical codebook        + external controls           ↓
   cleaned prices           + treatment/instruments       estimation
   three rollup levels      = analysis-ready panels       tables & figures

   SQL/                     Notebooks/                    Code/
   01_Phase1_Codebook       HPT_Python_Pipeline.ipynb     HPT_Analysis_Pipeline.R
   02_Phase2to4_Prices
   03_Supplementary
```

| Stage | Language | What it decides |
|---|---|---|
| **1** | Snowflake SQL | Which codes exist, what they are, what clinical concept they belong to, and what each hospital charges for them |
| **2** | Python | Who had already disclosed, what the instruments are, and what the market looks like |
| **3** | R | What the effect is |

Each stage's decisions are final. Stage 2 never re-derives a concept; stage 3 never rebuilds a price. The one deliberate exception is shoppability, which is *not* assigned in SQL — the raw ingredients are exported and the eighteen classification schemes are constructed in R, so competing definitions can be tested against each other without rebuilding the warehouse.

---

## Data

The negotiated-rate data are licensed and are **not** redistributed here. The code is published for methodological transparency and to document every specification decision; it is not runnable end to end without data access.

| Source | Contribution | Access |
|---|---|---|
| **Turquoise Health** | Negotiated rates parsed from hospital machine-readable files | Licensed |
| **CMS OPPS Addendum B** | CPT/HCPCS status indicators and APCs, six quarterly vintages | Public |
| **CMS MS-DRG (IPPS Table 5)** | DRG definitions, MDC, medical/surgical type, three fiscal years | Public |
| **CMS NCCI Add-on Code edits** | Authoritative add-on/component designation | Public |
| **CMS ASC Addendum AA** | ASC-covered procedures; the objective schedulability proxy | Public |
| **AHA Annual Survey** | Bed counts, system membership, ownership, hospital type | Licensed |
| **CMS enforcement actions** | Warning and closure notices (exploratory instruments only) | Public |
| **IPUMS ACS / tidycensus** | County demographics | Public |
| **Census county population estimates** | Coverage denominators | Public |

**Study window.** July 2024 through October 2025, spanning six OPPS quarters and three MS-DRG fiscal years.

**Panel.** The analysis panel is at the **hospital × clinical concept × post-month** level: roughly 1.4 million observations covering about 3,700 hospitals and 738 clinical concepts across 16 clinical families. A "concept" is a clinically coherent service defined in the stage 1 codebook, sitting between the individual billing code and the broad modality — *MRI Brain* rather than either CPT 70552 or *MRI*.

Roughly 31% of complete cases are dropped by cascading singleton removal across the county × concept and month fixed effects. This is expected, is identical across specifications, and is reported explicitly rather than left implicit — `audit_estimation_sample()` writes `QA01`, and the coverage block maps which counties actually contribute to identification versus merely appearing in the data.

---

## Empirical strategy

### Specification

Estimated with `fixest::feols`:

```
ln(P_ihct) = β · N_PRIOR_POSTERS_ht + γ · ln(Beds_h) + δ_mc + τ_t + ε_ihct
```

with `δ_mc` a county × concept market fixed effect (`MARKET_ID`), `τ_t` a post-month fixed effect, and standard errors two-way clustered on county and post-month. `N_PRIOR_POSTERS` counts hospitals in the same market that have already disclosed, and is instrumented by out-of-market disclosure rollout among competitor health systems.

### Four design decisions, fixed throughout

These are set once in the stage 3 header and are not re-litigated section by section.

**1. Treatment enters linearly.** Every concave transform of `N_PRIOR_POSTERS` inflates the IV coefficient by shrinking `Cov(T, Z)` in the Wald denominator, monotonically in how much it compresses the treatment. The reduced-form numerator does not move. The linear first stage runs F ≈ 43 in both interacted equations, so there is no weak-instrument problem for a transform to solve. The transform ladder is reported as robustness and is never a headline.

**2. The reduced form is the primary estimator.** IV is a ratio; the reduced form is its numerator, carrying no treatment variable, no functional-form choice, and no first-stage noise in the standard error. Because every model here is **exactly identified**, the Anderson–Rubin weak-instrument-robust test is numerically identical to the t-test on the reduced-form coefficient — so the reported RF p-values already *are* robust p-values. LIML offers nothing, since it coincides with 2SLS under exact identification.

**3. Heterogeneity is estimated by interaction, never by sample split.** The fixed effects absorb most of the instrument's effect on the treatment within any single subsample, so splitting destroys identification. A single function, `estimate_interacted()`, produces every heterogeneity result in the paper.

**4. No event study is possible.** All but two hospitals appear at exactly one posting month, so there is no within-hospital timing variation. This is a structural fact about disclosure, established in stage 1: `POST_MONTH` is an observed disclosure cohort, not a panel month, and no observation is carried forward. The selection concern is addressed structurally instead — the non-shoppable group is a **within-hospital placebo**: same hospital, same posting month, same instrument value, no price response. Any selection story has to explain why selection would operate on shoppable services and not on non-shoppable ones at the same hospital in the same month.

### Instruments

Six competitor-rollout instruments, built in stage 2, are graded on how credible their exclusion restriction is and how stable their sign is across classification schemes. **The grading is by construction and by audit, never by first-stage F** — the tier rerank shows the two strongest first stages produce the two weakest results, which is precisely why instrument selection here is a design decision rather than an outcome-based one.

| Tier | Instrument | Role |
|---|---|---|
| **MAIN** | `Competitor_only_hospitals_9m` *(primary)* | Carries every headline claim |
| **MAIN** | `Primary_strict_system_IV` | Headline |
| **MAIN** | `Competitor_outside_CBSA_hospitals_9m` | Headline |
| **CONFIRMING** | `Competitor_outside_CBSA_counties_9m` | Pooled robustness alongside MAIN |
| **CONFIRMING** | `Competitor_systems_9m` | Pooled robustness — see [Open items](#open-items) |
| **DISCREPANT** | `Competitor_outside_CBSA_systems_9m` | Reported alone, never pooled |

Only the `Z_SYS_COMPETITOR_*` family is hospital-specific and excludes the focal hospital's own system; the remaining measures are county-month local exposure and do not. That distinction is what the exclusion-restriction defence rests on, and it is recorded in the `EXCLUSION_NOTE` column of `T02`. It is also why the franchise-versus-integrated evidence in `03_Supplementary_Analyses.sql` is framed as being about `Primary_strict_system_IV` specifically rather than as a general defence.

Additional instruments are carried for robustness only: a construction ladder, a 3/6/9/12-month window ladder rebuilt from source, and CMS enforcement measures that remain **exploratory** — enforcement plausibly moves prices directly through compliance cost and correlated conduct, so it fails the exclusion restriction on its face.

### Shoppability classification

Eighteen classification schemes are built in R from the stage 1 codebook; six are attached to the panel and carry the results. A scheme is a partition into HIGH / INTERMEDIATE / LOW, assembled from named terms, where a term is a clinical family plus an optional keyword regex on the official code description plus an optional exclusion regex.

Three rules govern the construction:

- Emergency, critical care, and everything billed as MS-DRG is forced LOW in **every** scheme. `build_schemes()` asserts this rather than assuming it.
- Residual `"_Other"` terms mean whatever is left in a family after the named sub-terms are claimed *within the same scheme*, so `"CT Other"` is scheme-dependent by design.
- Inclusion keywords are substring matches against abbreviated CPT text and over-match. `CHEST` appears inside `"Ct chest spine w/o dye"` (a thoracic spine study) and `ULTRASOUND` inside `"Colonoscopy w/ultrasound"` (a procedure adjunct). `resolve_term()` therefore honours an `exclude` regex, and a guard tests for exactly that case.

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
├── SQL/
│   ├── 01_Phase1_Codebook.sql          stage 1a — clinical codebook
│   ├── 02_Phase2to4_Prices.sql         stage 1b — prices and exports
│   └── 03_Supplementary_Analyses.sql   stage 1c — supporting analyses
├── Notebooks/
│   └── HPT_Python_Pipeline.ipynb       stage 2 — controls and instruments
└── Code/
    └── HPT_Analysis_Pipeline.R         stage 3 — estimation
```

---

# Stage 1 — SQL (Snowflake)

Builds the clinical codebook and every price in the study, from the raw disclosure source.

## `01_Phase1_Codebook.sql`

Four things are kept strictly separate, because collapsing them is what makes a codebook unreproducible:

| Step | What it decides | How |
|---|---|---|
| **Universe** | Which codes exist with enough support | Automatic, one raw scan |
| **Validation** | Is the code real, current, an add-on, ASC-eligible | Automatic, from official CMS files |
| **Family / concept** | What clinical category and concept a code belongs to | A small, editable rule table |
| **Shoppability** | *Deliberately not decided here* | Built in stage 3 |

No step requires per-code manual review to run. Manual judgment enters in exactly one place — `HPT_REF_MANUAL_OVERRIDES` — and only for codes explicitly listed there with a stated reason.

**The family rule table** (`HPT_REF_FAMILY_RULES`) is the core of it. Every row is a rule, not a code-level judgment call, resolved by lowest priority number:

| Priority | Rule type | Purpose |
|---|---|---|
| 0 | `CODE_PATTERN`, `RANGE` | Exclusions: quality-measure codes, dental codes, anesthesia |
| 1 | `EXACT` | Genuine exceptions where a range or keyword rule would misfire |
| 10 | `RANGE` | Clean anatomic and modality blocks |
| 20 | `KEYWORD` | Safety net, matched against the *official* description |
| 30 | `RANGE` | Coarse CPT-section fallback, so nothing is dropped unclassified |

**The concept grouping** collapses administrative variants of the same clinical service. Contrast language is stripped to a stem, so *CT chest with contrast* and *CT chest without contrast* become one concept distinguished by `CONTRAST_VARIANT`. CMS short descriptors abbreviate inconsistently — many use "dye" rather than "contrast" — so both phrasings are matched. Where the stem is empty or generic, the fallback is one concept per exact code, so nothing is silently merged on missing data.

**Key outputs:** `HPT_P1_FINAL_CODEBOOK` (full archive, so scope can be widened later without rerunning) and `HPT_P1_FINAL_CODEBOOK_SCOPED` (what feeds the price pipeline).

## `02_Phase2to4_Prices.sql`

| Phase | Contents |
|---|---|
| **2** | The single expensive raw scan; row-level cleaning flags; one canonical file per hospital × month × code type; deduplication; payer-cell construction |
| **3** | Support-adaptive price bounds; hospital × month × exact-code prices |
| **4** | County posting-cohort exposure; concept and family rollups; hospital directory; exports |

**Cleaning.** A row survives only if it is an unmodified, direct dollar amount — not a percentage, algorithm, APC-derived, or per-diem rate — that is schema-compliant, not source-flagged as an outlier, and in the correct billing class and setting tier. Each condition is a separate named flag, so the QA output shows which condition removes what.

**Trimming is scaled to support**, so a thinly-observed code is not trimmed on percentiles estimated from a handful of cells: ≥1000 payer cells trims at P0.5/P99.5, 200–999 at P1/P99, and <200 is left untrimmed. Bounds are computed once per code and applied everywhere — including in the payer-class rebuild — so no comparison can be a trimming artifact.

**Sample restrictions happen at analysis time, not here.** Every price row carries its own support metrics through every stage, so decisions of the form "use only concepts with at least N hospitals" are made in R as documented, variable screens rather than baked into SQL as an irreversible filter.

**Rollups are equal-weighted at both levels** — median of exact-code medians within a concept, then median of concept medians within a family — so a concept spanning many billing-code variants is not driven by whichever variant happens to be most widely posted, and a family covering more body parts does not mechanically outweigh one covering few. Component (add-on) codes never enter a primary rollup.

**Exports** go to four data folders (`exact_standalone`, `exact_component`, `concept`, `family`) plus a QA folder holding the codebook, hospital directory, and three phase QA summaries. Stage 2 reads expected row counts from those QA tables rather than from hardcoded constants.

## `03_Supplementary_Analyses.sql`

| Part | Contents |
|---|---|
| **A** | Payer-class derivation and payer-conditional concept prices |
| **B** | Franchise versus integrated pricing within health systems |
| **C** | County payer concentration |

Part A classifies on **plan and network, not payer name**. The major carriers all sell commercial, Medicare Advantage, and managed Medicaid under the same corporate name, so a payer-name rule would produce carrier-conditional estimates wearing a payer-conditional disguise. The key validation asks whether each carrier appears in more than one class — it should.

Part B supplies the exclusion-restriction evidence for `Primary_strict_system_IV`, the one instrument that does not exclude own-system rollout. Note that a dispersion ratio of 1.0 is *not* the benchmark for local pricing: holding the system fixed also holds fixed brand, cost structure, chargemaster vendor, and negotiating staff, so some compression is expected even under fully decentralised negotiation.

Part C constructs a concentration measure the paper otherwise states is unavailable — with the caveat, stated in the file, that contract-row share is not enrollment share.

## Running stage 1

Load the four CMS reference files into `HPT_REF_STAGE`, then run each file's sections in order. **Do not run blind in one sitting.** Each phase ends with a QA table whose status must read `PASS`, and the expensive raw scans are worth checking before continuing past them.

---

# Stage 2 — Python

`Notebooks/HPT_Python_Pipeline.ipynb`. Constructs the treatment, the instruments, and every control variable, then writes the analysis-ready panels.

Nothing here repeats work done in SQL. Raw price cleaning, payer-cell construction, exact-code trimming, and concept/family aggregation all happen upstream; this notebook starts from the finished Snowflake exports.

| Section | Contents |
|---|---|
| **0–1** | Configuration and shared helpers. Everything configurable lives in section 0 |
| **2** | Load, reconcile, and normalize the Snowflake exports; attach the hospital directory |
| **3** | Census county geography |
| **4** | Hospital disclosure events and service entry flags |
| **5** | Strict prior-poster measures at city, county, and CBSA level |
| **6** | County-level health-system instruments |
| **7** | Optional service-specific fixed-roster system instruments |
| **8–10** | ACS demographic, hospital ownership, and CMS enforcement controls |
| **11** | Merge controls, classification, and analysis flags |
| **12** | Low-memory concept-, strict-concept-, and family-level measures |
| **13** | R-ready outputs |
| **14** | County coverage map and population coverage |
| **15–16** | Final QA, inventory, data dictionaries, run summary |

**Execution is disk-backed and low-memory.** The wide exact-code, concept, and family panels are never held in memory simultaneously where that can be avoided. Section 12 caches each panel to Parquet, builds narrow lookups, releases the wide panel, then enriches one billing-code domain at a time. Holding them all at once will exhaust RAM on a 16 GB machine, and the operating system terminates the kernel without a Python traceback when it does — a failure with no error message to debug.

**Outputs** under `Data/data_final`:

| Folder | Contents |
|---|---|
| `01_R_Analysis_Panels` | Analysis-ready panels consumed by stage 3 |
| `02_Treatments_and_Instruments` | Prior-poster measures and system instruments |
| `03_Controls` | ACS, ownership, and enforcement controls |
| `04_Coverage` | County and population coverage |
| `05_QA_and_Dictionaries` | QA results and variable dictionaries |

Section 15 runs a final QA pass over every written panel. **It must report zero FAIL rows before stage 3 is run.**

Requires Python 3.11+, with `pandas`, `numpy`, and `pyarrow`. The coverage map additionally needs `geopandas` and `matplotlib`, and is skipped when `BUILD_COVERAGE_MAP` is off.

---

# Stage 3 — R

`Code/HPT_Analysis_Pipeline.R`. Every regression, diagnostic, table, and figure in the paper.

Sections 0–12 define constants and functions and load in seconds. Everything after them executes, gated by `RUN_STAGES`.

### Definitions

| Section | Contents |
|---|---|
| **0** | Paths, baseline specification, instrument registry, scheme registry, sample screens. Everything configurable lives here |
| **1** | Helpers: coefficient extraction, first-stage diagnostics, formula builders, the `cache_or_run()` gate |
| **2** | Panel loading. Every function returns data rather than mutating a global, so `outpatient` is assigned exactly once from a visible call chain |
| **3** | Shoppability scheme construction, concept merges, scheme inheritance |
| **4** | Estimators, including `estimate_interacted()` |
| **5** | Instrument screen (`T02`) and pooled models (`T03`) |
| **6** | Concept-level estimation. Stores RF, FS, *and* IV per concept |
| **7** | Main interacted results and the transform ladder |
| **8** | Meta-regressions, partition deduplication, family permutation inference, RF-vs-FS decomposition |
| **9** | Price comparability as a candidate mechanism |
| **10** | Comparability meta-regression and horse race against the shoppability label |
| **11** | Instrument balance and placebo-outcome tests |
| **12** | Preflight. Ten seconds; must pass before the eight-hour stage 6 |

Keeping the first stage separately in section 6 is what makes the section 8 decomposition possible: if the gradient appeared in IV but not in RF, it would be a denominator artefact rather than a price response.

### Run blocks and robustness

| Block | Contents |
|---|---|
| **BUILD** | Builds schemes and panel, audits the estimation sample, runs preflight |
| **Instrument audit** | Raw strength, sign stability, and whether the headline changes under each of the six instruments (`QA05`–`QA08`) |
| **Stages 5–11** | The estimation stages, selected by `RUN_STAGES` |
| **IN_SYSTEM robustness** | Re-estimates the headline with system membership controlled |
| **County demographics** | ACS pull via `tidycensus`, with a mandatory spot check |
| **Section 12 / 12B** | Demographic heterogeneity in the *gradient*; SES terciles |
| **Section 13** | System × month FE, leave-one-system-out, construction ladder, randomisation inference |
| **Section 15** | CBSA market definition, instrument window ladder, payer-conditional |
| **Figures** | Publication figures 1–8 and robustness figures 10–20, read back from output CSVs so they run standalone |
| **Coverage** | County and population coverage on the *estimation sample*, plus a print-safe map |
| **LaTeX builds** | Summary statistics, family-level scheme assignment, shoppable-share inputs |

### Running stage 3

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

Every expensive step is wrapped in `cache_or_run()`. **A step producing zero rows never overwrites a good cache** — it stops instead, so a silent upstream failure cannot destroy hours of completed estimation. The stage 6 loader independently verifies that cached concept-level results were estimated on the panel currently in memory.

A few blocks are labelled **interactive checkpoints**. They duplicate what the BUILD block and stage 6 loader do, and exist so individual stages can be re-run without stepping through the whole file. Skip them on a cold `source()`.

Requires R 4.4+. Core: `arrow`, `data.table`, `fixest`, `dplyr`, `ggplot2`, `stringr`, `lubridate`, `MASS`. Optional, checked at point of use: `R.utils`, `tidycensus`, `ggrepel`, `patchwork`, `sf`, `tigris`, `ggpattern`, `magick`.

---

## Output index

Stage 3 results are written as CSV under `RESULT_ROOT`, prefixed by the section that produces them.

### Analysis tables

| File | Contents |
|---|---|
| `T02_instrument_first_stage_screen` | First stage and reduced form for every candidate instrument |
| `T03_pooled_OLS_RF_IV_all_outcomes` | Pooled estimates across all five price outcomes |
| `T05_concept_level_RF_FS_IV` | Concept-level RF, FS, and IV — the meta-regression input |
| `T06_main_*` | Headline interacted results, MAIN tier |
| `T06D_*` | Scheme 1 procedural-exclusion sensitivity |
| `T06F` / `T06G` | All three tiers combined, estimates and heterogeneity tests |
| `T06H_*` / `T06J_*` / `T06K_*` | Pooled MAIN + CONFIRMING; CONFIRMING and DISCREPANT tiers |
| `T06L_tier_rerank_all_schemes` | Tier diagnostic across all six schemes, with first-stage F alongside |
| `T07` / `T07B` | Transform ladder |
| `T08*` | Meta-regressions, redundant partitions, permutation inference, decomposition |
| `T09*` | Comparability measures, row-level interaction, within-family test |
| `T10` / `T10B` | Comparability meta-regression and horse race |
| `T11*` | Balance, placebo outcome, IN_SYSTEM robustness |
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

## Open items

Recorded here rather than papered over.

**Instrument tier vocabulary.** The R pipeline computes three tiers — MAIN / CONFIRMING / DISCREPANT — via `instrument_tier()`, and the tier-tagged tables use them. The instrument audit tables (`QA05`–`QA08`) use a two-way MAIN / SUPPORTING split, where SUPPORTING is the union of CONFIRMING and DISCREPANT. Both are live in the code. This README uses the three-tier vocabulary.

Separately, `Competitor_systems_9m` currently sits in CONFIRMING, but its shoppability differential in `T08F` is unstable across scheme variants even though its level coefficient is sign-stable. That is the criterion defining the DISCREPANT tier, so the assignment is under review.

**Concept count.** The panel carries 738 concepts. The concept-level sweep estimates only those clearing `MIN_SERVICE_OBS`, `MIN_SERVICE_MARKETS`, and `MIN_SERVICE_MONTHS`, so the estimable count reported in `T05` and `T06` differs from the panel count and from each other. The screens are in section 0 and the per-table counts are in the `N_CONCEPTS` column.

**Demographic heterogeneity.** Section 12 estimates demographic variation in the shoppability *gradient*. An earlier specification interacted the instrument with demographics alone, which measures variation in the pooled response — a quantity the paper disowns. Any table generated under the earlier specification should not be read as evidence about the gradient.

---

## Citation

> Sierra, Danny. "When Transparency Works: Service Shoppability, Contracting Depth, and the Price Effects of Hospital Disclosure." Working Paper, 2026.