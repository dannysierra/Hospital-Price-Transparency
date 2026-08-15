/* =============================================================================

   WHEN TRANSPARENCY WORKS: SERVICE SHOPPABILITY, CONTRACTING DEPTH, AND
   THE PRICE EFFECTS OF HOSPITAL DISCLOSURE

   Replication code -- stage 1c: supplementary analyses

   Danny Sierra
   Department of Economics, Florida State University
   Ds22c@fsu.edu

   -----------------------------------------------------------------------------
   WHAT THIS FILE DOES
   -----------------------------------------------------------------------------
   Three exercises that support specific claims in the paper. None of them
   touches the raw source, so none requires the expensive scan in
   02_Phase2to4_Prices.sql. All run against tables that file already built.

     PART A  Payer-class derivation and payer-conditional concept prices.
             Supports the payer-conditional robustness check.

     PART B  Franchise versus integrated pricing within health systems.
             Supplies the exclusion-restriction evidence for the one
             instrument that does not exclude own-system rollout.

     PART C  County payer concentration.
             Constructs a concentration measure the paper otherwise states is
             unavailable. Read the caveat before using it as a test.

   -----------------------------------------------------------------------------
   REQUIRED INPUT
   -----------------------------------------------------------------------------
     O_ECON.COMMON.HPT_P2_PAYER_CELLS
     O_ECON.COMMON.HPT_P3_CODE_PRICE_BOUNDS

   Snowflake's REGEXP_LIKE full-match behaviour applies throughout this file
   exactly as described in the header of 02_Phase2to4_Prices.sql: every
   substring-intent pattern is wrapped as '.*(' || PATTERN || ').*'.
   ============================================================================= */

USE WAREHOUSE ECON_DEF_WH;
USE DATABASE O_ECON;
USE SCHEMA COMMON;


/* =============================================================================
   =============================================================================
   PART A -- PAYER-CONDITIONAL ROBUSTNESS
   =============================================================================

   THE PROBLEM. HPT_P2_PAYER_CELLS carries no payer-class column. What it has
   is CANONICAL_PAYER_NAME, CANONICAL_PLAN_NAME, and CANONICAL_NETWORK_NAME, so
   payer class has to be derived from those strings.

   THE ONE THING THAT MATTERS. Classify on PLAN and NETWORK, not on PAYER NAME.
   The major carriers all sell commercial, Medicare Advantage, and managed
   Medicaid products under the same corporate name. A payer-name rule would
   assign every product line of a carrier to a single class, and the resulting
   "payer-conditional" estimates would be carrier-conditional estimates in
   disguise -- a different and far less interesting object.

   THREE DESIGN POINTS.

   1. USE THE POOLED PRICE BOUNDS from HPT_P3_CODE_PRICE_BOUNDS, computed per
      code across all payer cells. Recomputing bounds within class would let
      the trimming rule differ by class -- managed Medicaid has far fewer cells,
      so its rule would shift from P0.5/P99.5 to no trimming -- and any
      cross-class difference could then be a trimming artifact rather than a
      pricing fact. Same bounds, applied within class, is the right comparison.

   2. DO NOT REBUILD THE TREATMENT OR THE INSTRUMENT. Both are defined at the
      hospital-posting level and are identical across payer classes. Only the
      OUTCOME varies, which makes this a concept-price rebuild plus a merge
      rather than a pipeline rerun.

   3. EXPECT MANAGED MEDICAID TO BE THIN. Check hospital counts and the first
      stage before reading anything into that class. Where the first stage
      collapses, report it as underpowered rather than as a null.
   ============================================================================= */


/* -----------------------------------------------------------------------------
   A0. INVENTORY FIRST.

   Run these two queries and read the output before trusting the rule table
   below. The rules are built on common US product-line naming conventions and
   are a starting point, not a substitute for inspecting the actual strings in
   this extract.
   ----------------------------------------------------------------------------- */

SELECT
    CANONICAL_PLAN_NAME,
    CANONICAL_NETWORK_NAME,
    COUNT(*)                          AS N_CELLS,
    COUNT(DISTINCT HOSPITAL_ID)       AS N_HOSPITALS,
    COUNT(DISTINCT CANONICAL_PAYER_NAME) AS N_CARRIERS
FROM O_ECON.COMMON.HPT_P2_PAYER_CELLS
WHERE CODE_ROLE = 'STANDALONE'
  AND MISSING_PAYER_IDENTITY_FLAG = 0
GROUP BY CANONICAL_PLAN_NAME, CANONICAL_NETWORK_NAME
ORDER BY N_CELLS DESC
LIMIT 300;

-- How much volume sits in plan and network strings that carry no product-line
-- information at all? Where PCT_UNCLASSIFIABLE is large, the payer-conditional
-- split covers the classifiable subsample only, and the paper should say so.
SELECT
    COUNT(*) AS N_CELLS,
    COUNT_IF(CANONICAL_PLAN_NAME = 'UNKNOWN_PLAN')       AS N_UNKNOWN_PLAN,
    COUNT_IF(CANONICAL_NETWORK_NAME = 'UNKNOWN_NETWORK') AS N_UNKNOWN_NETWORK,
    COUNT_IF(CANONICAL_PLAN_NAME = 'UNKNOWN_PLAN'
         AND CANONICAL_NETWORK_NAME = 'UNKNOWN_NETWORK') AS N_BOTH_UNKNOWN,
    ROUND(100 * COUNT_IF(CANONICAL_PLAN_NAME = 'UNKNOWN_PLAN'
         AND CANONICAL_NETWORK_NAME = 'UNKNOWN_NETWORK') / COUNT(*), 2) AS PCT_UNCLASSIFIABLE
FROM O_ECON.COMMON.HPT_P2_PAYER_CELLS
WHERE CODE_ROLE = 'STANDALONE' AND MISSING_PAYER_IDENTITY_FLAG = 0;


/* -----------------------------------------------------------------------------
   A1. PAYER-CLASS RULE TABLE

   Priority-ordered in the same way as HPT_REF_FAMILY_RULES: lowest PRIORITY
   wins. Order is load-bearing here.

   Medicaid is checked before Medicare because dual-eligible and
   Medicare-Medicaid products mention both programs, and because state Medicaid
   brand names are unambiguous where the word "Medicare" is not. Medicare
   supplement products are removed first of all, at priority 5, because they
   are commercially underwritten wrap policies rather than Part C and would
   otherwise be swept into Medicare Advantage by the bare MEDICARE pattern.
   ----------------------------------------------------------------------------- */

CREATE OR REPLACE TABLE O_ECON.COMMON.HPT_REF_PAYER_CLASS_RULES (
    PRIORITY        NUMBER,
    MATCH_FIELD     VARCHAR,   -- PLAN | NETWORK | PLAN_OR_NETWORK | CARRIER
    PATTERN         VARCHAR,
    PAYER_CLASS     VARCHAR,   -- MANAGED_MEDICAID | MEDICARE_ADVANTAGE | COMMERCIAL
    RULE_NOTE       VARCHAR
);

INSERT INTO O_ECON.COMMON.HPT_REF_PAYER_CLASS_RULES VALUES

    /* -- Priority 5: exclusions that would otherwise be read as Medicare -- */
    (5, 'PLAN_OR_NETWORK', 'MEDIGAP|MEDICARE SUPPLEMENT|SUPPLEMENTAL|MED SUPP',
        'COMMERCIAL',
        'Medicare supplement products are commercially underwritten wrap policies, not Part C. Must be caught before the Medicare rules below.'),

    /* -- Priority 10: managed Medicaid, including state brand names -- */
    (10, 'PLAN_OR_NETWORK',
        'MEDICAID|MCAID|(^|[^A-Z])CHIP([^A-Z]|$)|HUSKY|SOONERCARE|BADGERCARE'
        || '|MEDI-CAL|MEDICAL ASSISTANCE|TENNCARE|MASSHEALTH|APPLE HEALTH'
        || '|HEALTHCHOICE|PEACHCARE|HOOSIER|KIDMED|DENALI|AHCCCS|FAMIS'
        || '|HEALTHY BLUE|COMMUNITY PLAN|STAR KIDS|STAR PLUS|STAR\\+PLUS',
        'MANAGED_MEDICAID',
        'State Medicaid managed-care brand names plus the generic terms. Checked before Medicare because dual-eligible products name both programs.'),

    /* -- Priority 20: Medicare Advantage / Part C -- */
    (20, 'PLAN_OR_NETWORK',
        'MEDICARE ADVANTAGE|(^|[^A-Z])MA([^A-Z]|$)|PART C|MEDICARE HMO|MEDICARE PPO'
        || '|(^|[^A-Z])SNP([^A-Z]|$)|SPECIAL NEEDS|DUAL ELIGIBLE|(^|[^A-Z])DSNP([^A-Z]|$)'
        || '|MEDICARE COMPLETE|MEDICARE CHOICE|SENIOR ADVANTAGE|(^|[^A-Z])MCARE([^A-Z]|$)'
        || '|MEDICARE',
        'MEDICARE_ADVANTAGE',
        'Bare MEDICARE is last within this pattern so the more specific forms win first; supplements were already removed at priority 5.'),

    /* -- Priority 30: explicit commercial signals -- */
    (30, 'PLAN_OR_NETWORK',
        'COMMERCIAL|EXCHANGE|MARKETPLACE|(^|[^A-Z])ACA([^A-Z]|$)|EMPLOYER|GROUP'
        || '|(^|[^A-Z])PPO([^A-Z]|$)|(^|[^A-Z])HMO([^A-Z]|$)|(^|[^A-Z])EPO([^A-Z]|$)'
        || '|(^|[^A-Z])POS([^A-Z]|$)|OPEN ACCESS|CHOICE PLUS|BLUE CARD|NATIONAL',
        'COMMERCIAL',
        'Generic commercial product-line language. Reached only if no Medicaid or Medicare rule fired.')
;


/* -----------------------------------------------------------------------------
   A2. APPLY THE RULES -- one class per payer cell.
   ----------------------------------------------------------------------------- */

CREATE OR REPLACE TABLE O_ECON.COMMON.HPT_P2_PAYER_CELL_CLASS AS

WITH BASE AS (
    SELECT
        PAYER_CELL_ID,
        UPPER(COALESCE(CANONICAL_PLAN_NAME, ''))    AS PLAN_U,
        UPPER(COALESCE(CANONICAL_NETWORK_NAME, '')) AS NETWORK_U,
        UPPER(COALESCE(CANONICAL_PAYER_NAME, ''))   AS CARRIER_U
    FROM O_ECON.COMMON.HPT_P2_PAYER_CELLS
),

MATCHED AS (
    SELECT
        B.PAYER_CELL_ID,
        R.PRIORITY,
        R.PAYER_CLASS,
        R.MATCH_FIELD
    FROM BASE B
    INNER JOIN O_ECON.COMMON.HPT_REF_PAYER_CLASS_RULES R
        ON (
             (R.MATCH_FIELD IN ('PLAN','PLAN_OR_NETWORK')
                AND REGEXP_LIKE(B.PLAN_U,    '.*(' || R.PATTERN || ').*'))
          OR (R.MATCH_FIELD IN ('NETWORK','PLAN_OR_NETWORK')
                AND REGEXP_LIKE(B.NETWORK_U, '.*(' || R.PATTERN || ').*'))
          OR (R.MATCH_FIELD = 'CARRIER'
                AND REGEXP_LIKE(B.CARRIER_U, '.*(' || R.PATTERN || ').*'))
           )
    QUALIFY ROW_NUMBER() OVER (PARTITION BY B.PAYER_CELL_ID ORDER BY R.PRIORITY) = 1
)

SELECT
    B.PAYER_CELL_ID,
    COALESCE(M.PAYER_CLASS, 'UNCLASSIFIED') AS PAYER_CLASS,
    M.PRIORITY   AS MATCHED_RULE_PRIORITY,
    M.MATCH_FIELD AS MATCHED_FIELD
FROM BASE B
LEFT JOIN MATCHED M ON B.PAYER_CELL_ID = M.PAYER_CELL_ID
;


/* -- A3. AUDIT THE CLASSIFICATION before using it for anything. -- */

SELECT
    C.PAYER_CLASS,
    COUNT(*)                                  AS N_CELLS,
    ROUND(100 * COUNT(*) / SUM(COUNT(*)) OVER (), 2) AS PCT_CELLS,
    COUNT(DISTINCT P.HOSPITAL_ID)             AS N_HOSPITALS,
    COUNT(DISTINCT P.ANALYSIS_CONCEPT_ID)     AS N_CONCEPTS,
    COUNT(DISTINCT P.CANONICAL_PAYER_NAME)    AS N_CARRIERS,
    ROUND(MEDIAN(P.PAYER_CELL_RATE), 2)       AS MEDIAN_RATE
FROM O_ECON.COMMON.HPT_P2_PAYER_CELLS P
INNER JOIN O_ECON.COMMON.HPT_P2_PAYER_CELL_CLASS C USING (PAYER_CELL_ID)
WHERE P.CODE_ROLE = 'STANDALONE' AND P.MISSING_PAYER_IDENTITY_FLAG = 0
GROUP BY C.PAYER_CLASS
ORDER BY N_CELLS DESC;

-- THE KEY VALIDATION. Does each carrier appear in more than one class? It
-- should: the major carriers all sell multiple product lines. A carrier that
-- lands entirely in one class is a sign the rules are keying on carrier name
-- rather than on product line, which is the failure mode described above.
SELECT
    P.CANONICAL_PAYER_NAME,
    COUNT(DISTINCT C.PAYER_CLASS) AS N_CLASSES,
    LISTAGG(DISTINCT C.PAYER_CLASS, ', ') WITHIN GROUP (ORDER BY C.PAYER_CLASS) AS CLASSES,
    COUNT(*) AS N_CELLS
FROM O_ECON.COMMON.HPT_P2_PAYER_CELLS P
INNER JOIN O_ECON.COMMON.HPT_P2_PAYER_CELL_CLASS C USING (PAYER_CELL_ID)
WHERE P.CODE_ROLE = 'STANDALONE' AND P.MISSING_PAYER_IDENTITY_FLAG = 0
GROUP BY P.CANONICAL_PAYER_NAME
HAVING COUNT(*) >= 5000
ORDER BY N_CELLS DESC
LIMIT 40;

-- Spot check: 40 random cells per class alongside their source strings.
SELECT C.PAYER_CLASS, P.CANONICAL_PAYER_NAME, P.CANONICAL_PLAN_NAME,
       P.CANONICAL_NETWORK_NAME, C.MATCHED_RULE_PRIORITY
FROM O_ECON.COMMON.HPT_P2_PAYER_CELLS P
INNER JOIN O_ECON.COMMON.HPT_P2_PAYER_CELL_CLASS C USING (PAYER_CELL_ID)
WHERE P.CODE_ROLE = 'STANDALONE'
QUALIFY ROW_NUMBER() OVER (PARTITION BY C.PAYER_CLASS ORDER BY RANDOM()) <= 40
ORDER BY C.PAYER_CLASS;


/* -----------------------------------------------------------------------------
   A4. CLASS-CONDITIONAL CONCEPT PRICES

   Mirrors the exact-code and concept aggregations from
   02_Phase2to4_Prices.sql, with one added grouping key and the pooled bounds.
   Output is one row per hospital x month x concept x payer class.
   ----------------------------------------------------------------------------- */

CREATE OR REPLACE TABLE O_ECON.COMMON.HPT_P4_CONCEPT_PRICES_BY_PAYER_CLASS AS

WITH FLAGGED AS (
    SELECT
        P.*,
        C.PAYER_CLASS,
        IFF(P.PAYER_CELL_RATE BETWEEN B.LOWER_BOUND AND B.UPPER_BOUND, 1, 0) AS IN_BOUND_FLAG
    FROM O_ECON.COMMON.HPT_P2_PAYER_CELLS P
    INNER JOIN O_ECON.COMMON.HPT_P2_PAYER_CELL_CLASS C USING (PAYER_CELL_ID)
    INNER JOIN O_ECON.COMMON.HPT_P3_CODE_PRICE_BOUNDS B
        ON P.BILLING_CODE_TYPE = B.BILLING_CODE_TYPE AND P.BILLING_CODE = B.BILLING_CODE
    WHERE P.CODE_ROLE = 'STANDALONE'
      AND P.MISSING_PAYER_IDENTITY_FLAG = 0
      AND C.PAYER_CLASS <> 'UNCLASSIFIED'
),

EXACT_CODE AS (
    SELECT
        HOSPITAL_ID, POST_MONTH, PAYER_CLASS,
        BILLING_CODE_TYPE, BILLING_CODE,
        MIN(PROVIDER_STATE) AS PROVIDER_STATE, MIN(COUNTY_STATE) AS COUNTY_STATE,
        MIN(CBSA_CODE) AS CBSA_CODE, MIN(HEALTH_SYSTEM_ID) AS HEALTH_SYSTEM_ID,
        MEDIAN(TOTAL_BEDS) AS TOTAL_BEDS,
        MIN(ANALYSIS_SUPERFAMILY_ID) AS ANALYSIS_SUPERFAMILY_ID,
        MIN(ANALYSIS_FAMILY_ID) AS ANALYSIS_FAMILY_ID,
        MIN(ANALYSIS_CONCEPT_ID) AS ANALYSIS_CONCEPT_ID,
        COUNT_IF(IN_BOUND_FLAG = 1) AS N_PAYER_CELLS,
        COUNT(DISTINCT IFF(IN_BOUND_FLAG = 1, CANONICAL_PAYER_KEY, NULL)) AS N_DISTINCT_PAYERS,
        MEDIAN(IFF(IN_BOUND_FLAG = 1, PAYER_CELL_RATE, NULL))     AS MEDIAN_PRICE,
        MEDIAN(IFF(IN_BOUND_FLAG = 1, LN(PAYER_CELL_RATE), NULL)) AS MEDIAN_LOG_PRICE
    FROM FLAGGED
    GROUP BY HOSPITAL_ID, POST_MONTH, PAYER_CLASS, BILLING_CODE_TYPE, BILLING_CODE
    HAVING COUNT_IF(IN_BOUND_FLAG = 1) > 0
)

SELECT
    HOSPITAL_ID,
    MIN(PROVIDER_STATE) AS PROVIDER_STATE,
    MIN(COUNTY_STATE)   AS COUNTY_STATE,
    MIN(CBSA_CODE)      AS CBSA_CODE,
    MIN(HEALTH_SYSTEM_ID) AS HEALTH_SYSTEM_ID,
    MEDIAN(TOTAL_BEDS)  AS TOTAL_BEDS,
    POST_MONTH,
    PAYER_CLASS,
    BILLING_CODE_TYPE,
    MIN(ANALYSIS_SUPERFAMILY_ID) AS ANALYSIS_SUPERFAMILY_ID,
    MIN(ANALYSIS_FAMILY_ID)      AS ANALYSIS_FAMILY_ID,
    ANALYSIS_CONCEPT_ID,
    COUNT(*)                       AS N_CODES_IN_CONCEPT,
    SUM(N_PAYER_CELLS)             AS N_PAYER_CELLS_TOTAL,
    SUM(N_DISTINCT_PAYERS)         AS SUM_CODE_DISTINCT_PAYERS,
    ROUND(MEDIAN(MEDIAN_PRICE), 6) AS MEDIAN_PRICE,
    ROUND(AVG(MEDIAN_PRICE), 6)    AS MEAN_PRICE,
    ROUND(PERCENTILE_CONT(0.25) WITHIN GROUP (ORDER BY MEDIAN_PRICE), 6) AS P25_PRICE,
    ROUND(PERCENTILE_CONT(0.75) WITHIN GROUP (ORDER BY MEDIAN_PRICE), 6) AS P75_PRICE,
    MEDIAN(MEDIAN_LOG_PRICE)       AS MEDIAN_LOG_PRICE
FROM EXACT_CODE
GROUP BY HOSPITAL_ID, POST_MONTH, PAYER_CLASS, BILLING_CODE_TYPE, ANALYSIS_CONCEPT_ID
;

-- QA: coverage per class. Compare N_HOSPITALS against the pooled panel. A
-- class far below it is underpowered, which is a different finding from a null.
SELECT
    PAYER_CLASS,
    COUNT(*) AS N_ROWS,
    COUNT(DISTINCT HOSPITAL_ID) AS N_HOSPITALS,
    COUNT(DISTINCT ANALYSIS_CONCEPT_ID) AS N_CONCEPTS,
    COUNT(DISTINCT COUNTY_STATE) AS N_COUNTIES,
    ROUND(MEDIAN(MEDIAN_PRICE), 2) AS MEDIAN_PRICE
FROM O_ECON.COMMON.HPT_P4_CONCEPT_PRICES_BY_PAYER_CLASS
GROUP BY PAYER_CLASS ORDER BY N_ROWS DESC;


/* -- A5. EXPORT. In R, read this, merge the existing treatment and instrument
   columns onto it by hospital, post-month, and concept, and estimate once per
   payer class. Nothing about the instrument changes; only the outcome does. -- */

COPY INTO @O_ECON.COMMON.HPT_PY_EXPORT_STAGE/data/concept_by_payer_class/HPT_CONCEPT_PAYER_CLASS.csv.gz
FROM (
    SELECT * FROM O_ECON.COMMON.HPT_P4_CONCEPT_PRICES_BY_PAYER_CLASS
    ORDER BY PAYER_CLASS, PROVIDER_STATE, HOSPITAL_ID, POST_MONTH, ANALYSIS_CONCEPT_ID
)
FILE_FORMAT = (TYPE = CSV FIELD_OPTIONALLY_ENCLOSED_BY = '"' COMPRESSION = GZIP NULL_IF = (''))
HEADER = TRUE MAX_FILE_SIZE = 5368709120 OVERWRITE = TRUE
;

ALTER STAGE O_ECON.COMMON.HPT_PY_EXPORT_STAGE REFRESH;

SELECT RELATIVE_PATH, ROUND(SIZE/1024/1024, 2) AS SIZE_MB,
       GET_PRESIGNED_URL(@O_ECON.COMMON.HPT_PY_EXPORT_STAGE, RELATIVE_PATH, 604800) AS DOWNLOAD_URL
FROM DIRECTORY(@O_ECON.COMMON.HPT_PY_EXPORT_STAGE)
WHERE RELATIVE_PATH LIKE 'data/concept_by_payer_class/%';



/* =============================================================================
   =============================================================================
   PART B -- FRANCHISE VERSUS INTEGRATED PRICING
   =============================================================================

   WHAT THIS BEARS ON. If health systems set prices centrally across all their
   hospitals, then a hospital's own system's rollout could move its price
   through a channel other than local competitive pressure. The two competitor
   instruments exclude the focal hospital's own system by construction, so
   within-system centralisation cannot contaminate them. This evidence is about
   Primary_strict_system_IV, which does not exclude own-system rollout, and the
   paper should frame it that way rather than as a general defence.

   TWO CONSTRUCTION POINTS.

   1. CPT/HCPCS ONLY, everywhere in this part. Pooling MS-DRG inpatient
      families in makes the result look ambiguous when it is not: inpatient
      families run a substantially higher system-to-market dispersion ratio
      than outpatient ones, and they are excluded from the paper's sample
      anyway.

   2. UNIT-COUNT MATCHING for any tie-rate comparison. Tie probability falls
      mechanically as the number of compared units rises. A within-system cell
      spans however many counties the system operates in, while a within-market
      cell usually spans two or three hospitals, so an unmatched comparison
      confounds centralisation with cell size. B2a restricts both arms to
      exactly two units; B2b reports the unit-count distributions so the
      imbalance is visible rather than assumed away.

   WHAT TO REPORT. B3 is the reportable table. A standard deviation is unbiased
   in the number of observations; a tie rate is not. B2 is a supporting
   descriptive.

   HOW TO READ THE RATIO. A ratio of 1.0 is not the right benchmark for local
   pricing. Holding the system fixed also holds fixed brand, cost structure,
   chargemaster vendor, technology, and negotiating staff -- none of which are
   held fixed when comparing unrelated hospitals in a market. Some compression
   is expected even under fully decentralised negotiation.
   ============================================================================= */


/* -- B0. Base table. Outpatient CPT/HCPCS, standalone codes, clean payer
   identity, systems only. -- */

CREATE OR REPLACE TEMPORARY TABLE HPT_TMP_SYS_PRICES AS
SELECT
    HEALTH_SYSTEM_NAME, COUNTY_STATE, HOSPITAL_ID, CANONICAL_PAYER_NAME,
    ANALYSIS_FAMILY_ID, ANALYSIS_CONCEPT_ID, BILLING_CODE,
    MEDIAN(PAYER_CELL_RATE) AS MEDIAN_PRICE
FROM O_ECON.COMMON.HPT_P2_PAYER_CELLS
WHERE PAYER_CELL_RATE > 0
  AND BILLING_CODE_TYPE = 'CPT_HCPCS'
  AND CODE_ROLE = 'STANDALONE'
  AND HEALTH_SYSTEM_NAME IS NOT NULL
  AND MISSING_PAYER_IDENTITY_FLAG = 0
GROUP BY 1,2,3,4,5,6,7
;

SELECT COUNT(*) AS N_ROWS, COUNT(DISTINCT HEALTH_SYSTEM_NAME) AS N_SYSTEMS,
       COUNT(DISTINCT ANALYSIS_FAMILY_ID) AS N_FAMILIES
FROM HPT_TMP_SYS_PRICES;
-- N_FAMILIES should be the outpatient family count, all CPT. Any MDC_* family
-- appearing here means the billing-code-type filter did not apply.


/* -- B2a. TIE RATE, unit-count matched at exactly two units per cell. -- */

WITH WS AS (
    SELECT HEALTH_SYSTEM_NAME, CANONICAL_PAYER_NAME, ANALYSIS_FAMILY_ID, BILLING_CODE,
           COUNT(DISTINCT COUNTY_STATE) AS N_UNITS,
           COUNT(DISTINCT MEDIAN_PRICE) AS N_DISTINCT
    FROM HPT_TMP_SYS_PRICES GROUP BY 1,2,3,4
    HAVING COUNT(DISTINCT COUNTY_STATE) = 2
),
WM AS (
    SELECT COUNTY_STATE, CANONICAL_PAYER_NAME, ANALYSIS_FAMILY_ID, BILLING_CODE,
           COUNT(DISTINCT HOSPITAL_ID)  AS N_UNITS,
           COUNT(DISTINCT MEDIAN_PRICE) AS N_DISTINCT
    FROM HPT_TMP_SYS_PRICES GROUP BY 1,2,3,4
    HAVING COUNT(DISTINCT HOSPITAL_ID) = 2
)
SELECT 'within_system_cross_county' AS COMPARISON, ANALYSIS_FAMILY_ID,
       COUNT(*) AS N_CELLS,
       ROUND(100 * COUNT_IF(N_DISTINCT = 1) / COUNT(*), 2) AS PCT_IDENTICAL_PRICE
FROM WS GROUP BY 2
UNION ALL
SELECT 'within_market_cross_hospital', ANALYSIS_FAMILY_ID, COUNT(*),
       ROUND(100 * COUNT_IF(N_DISTINCT = 1) / COUNT(*), 2)
FROM WM GROUP BY 2
ORDER BY ANALYSIS_FAMILY_ID, COMPARISON;


/* -- B2b. How many units does each arm actually compare, unrestricted? Report
   this in the table note so the imbalance B2a corrects for is documented. -- */

WITH WS AS (
    SELECT COUNT(DISTINCT COUNTY_STATE) AS N_UNITS
    FROM HPT_TMP_SYS_PRICES
    GROUP BY HEALTH_SYSTEM_NAME, CANONICAL_PAYER_NAME, ANALYSIS_FAMILY_ID, BILLING_CODE
    HAVING COUNT(DISTINCT COUNTY_STATE) >= 2
),
WM AS (
    SELECT COUNT(DISTINCT HOSPITAL_ID) AS N_UNITS
    FROM HPT_TMP_SYS_PRICES
    GROUP BY COUNTY_STATE, CANONICAL_PAYER_NAME, ANALYSIS_FAMILY_ID, BILLING_CODE
    HAVING COUNT(DISTINCT HOSPITAL_ID) >= 2
)
SELECT 'within_system_cross_county' AS COMPARISON, COUNT(*) AS N_CELLS,
       ROUND(AVG(N_UNITS), 2) AS MEAN_UNITS, MEDIAN(N_UNITS) AS MEDIAN_UNITS,
       MAX(N_UNITS) AS MAX_UNITS FROM WS
UNION ALL
SELECT 'within_market_cross_hospital', COUNT(*), ROUND(AVG(N_UNITS), 2),
       MEDIAN(N_UNITS), MAX(N_UNITS) FROM WM;


/* -- B3. THE REPORTABLE TABLE. Dispersion ratio by family, CPT only, with a
   consistent minimum of three units in both arms. -- */

WITH SYSTEM_COUNTY AS (
    SELECT HEALTH_SYSTEM_NAME, COUNTY_STATE, CANONICAL_PAYER_NAME,
           ANALYSIS_FAMILY_ID, ANALYSIS_CONCEPT_ID,
           AVG(LN(MEDIAN_PRICE)) AS MEAN_LN_PRICE
    FROM HPT_TMP_SYS_PRICES GROUP BY 1,2,3,4,5
),
WS AS (
    SELECT ANALYSIS_FAMILY_ID, COUNT(*) AS N_CELLS,
           MEDIAN(SD_LN) AS MEDIAN_SD, AVG(N_UNITS) AS MEAN_UNITS
    FROM (
        SELECT ANALYSIS_FAMILY_ID, STDDEV(MEAN_LN_PRICE) AS SD_LN,
               COUNT(DISTINCT COUNTY_STATE) AS N_UNITS
        FROM SYSTEM_COUNTY
        GROUP BY HEALTH_SYSTEM_NAME, CANONICAL_PAYER_NAME, ANALYSIS_FAMILY_ID, ANALYSIS_CONCEPT_ID
        HAVING COUNT(DISTINCT COUNTY_STATE) >= 3
    ) WHERE SD_LN IS NOT NULL AND SD_LN < 2
    GROUP BY ANALYSIS_FAMILY_ID
),
WM AS (
    SELECT ANALYSIS_FAMILY_ID, COUNT(*) AS N_CELLS,
           MEDIAN(SD_LN) AS MEDIAN_SD, AVG(N_UNITS) AS MEAN_UNITS
    FROM (
        SELECT ANALYSIS_FAMILY_ID, STDDEV(LN(MEDIAN_PRICE)) AS SD_LN,
               COUNT(DISTINCT HOSPITAL_ID) AS N_UNITS
        FROM HPT_TMP_SYS_PRICES
        GROUP BY COUNTY_STATE, CANONICAL_PAYER_NAME, ANALYSIS_FAMILY_ID, ANALYSIS_CONCEPT_ID
        HAVING COUNT(DISTINCT HOSPITAL_ID) >= 3
    ) WHERE SD_LN IS NOT NULL AND SD_LN < 2
    GROUP BY ANALYSIS_FAMILY_ID
)
SELECT
    COALESCE(S.ANALYSIS_FAMILY_ID, M.ANALYSIS_FAMILY_ID) AS ANALYSIS_FAMILY_ID,
    S.N_CELLS AS N_CELLS_SYSTEM,  ROUND(S.MEAN_UNITS, 2) AS MEAN_COUNTIES,
    M.N_CELLS AS N_CELLS_MARKET,  ROUND(M.MEAN_UNITS, 2) AS MEAN_HOSPITALS,
    ROUND(S.MEDIAN_SD, 4) AS MEDIAN_SD_SYSTEM,
    ROUND(M.MEDIAN_SD, 4) AS MEDIAN_SD_MARKET,
    ROUND(S.MEDIAN_SD / NULLIF(M.MEDIAN_SD, 0), 3) AS RATIO_SYSTEM_TO_MARKET
FROM WS S FULL OUTER JOIN WM M USING (ANALYSIS_FAMILY_ID)
ORDER BY RATIO_SYSTEM_TO_MARKET;


/* -- B3b. Is the ratio an artifact of how many units each arm compares? Split
   the within-system arm by county count. A MEDIAN_SD flat in the county count
   means unit counts are not driving the comparison. -- */

WITH SYSTEM_COUNTY AS (
    SELECT HEALTH_SYSTEM_NAME, COUNTY_STATE, CANONICAL_PAYER_NAME,
           ANALYSIS_CONCEPT_ID, AVG(LN(MEDIAN_PRICE)) AS MEAN_LN_PRICE
    FROM HPT_TMP_SYS_PRICES GROUP BY 1,2,3,4
)
SELECT
    CASE WHEN N_UNITS = 2 THEN '2 counties'
         WHEN N_UNITS = 3 THEN '3 counties'
         WHEN N_UNITS BETWEEN 4 AND 5 THEN '4-5 counties'
         ELSE '6+ counties' END AS COUNTY_BUCKET,
    COUNT(*) AS N_CELLS, ROUND(MEDIAN(SD_LN), 4) AS MEDIAN_SD
FROM (
    SELECT STDDEV(MEAN_LN_PRICE) AS SD_LN, COUNT(DISTINCT COUNTY_STATE) AS N_UNITS
    FROM SYSTEM_COUNTY
    GROUP BY HEALTH_SYSTEM_NAME, CANONICAL_PAYER_NAME, ANALYSIS_CONCEPT_ID
    HAVING COUNT(DISTINCT COUNTY_STATE) >= 2
) WHERE SD_LN IS NOT NULL AND SD_LN < 2
GROUP BY 1 ORDER BY 1;


/* -- B4. By system, with the weight each system carries.

   Some systems do price centrally, with near-zero cross-county dispersion. The
   question is whether that matters for the design, which is what
   PCT_OF_ALL_HOSPITALS answers. A centrally priced system holding a fraction
   of a percent of hospitals is a footnote; one holding several percent belongs
   in the leave-one-system-out robustness, dropped by name, with the headline
   shown alongside. -- */

WITH SYS_SIZE AS (
    SELECT HEALTH_SYSTEM_NAME,
           COUNT(DISTINCT HOSPITAL_ID)  AS N_HOSPITALS,
           COUNT(DISTINCT COUNTY_STATE) AS N_COUNTIES
    FROM HPT_TMP_SYS_PRICES GROUP BY 1
    QUALIFY ROW_NUMBER() OVER (ORDER BY COUNT(DISTINCT HOSPITAL_ID) DESC) <= 20
),
SYSTEM_COUNTY AS (
    SELECT P.HEALTH_SYSTEM_NAME, P.COUNTY_STATE, P.CANONICAL_PAYER_NAME,
           P.ANALYSIS_CONCEPT_ID, AVG(LN(P.MEDIAN_PRICE)) AS MEAN_LN_PRICE
    FROM HPT_TMP_SYS_PRICES P INNER JOIN SYS_SIZE S USING (HEALTH_SYSTEM_NAME)
    GROUP BY 1,2,3,4
),
DISP AS (
    SELECT HEALTH_SYSTEM_NAME, COUNT(*) AS N_CELLS,
           MEDIAN(SD_LN) AS MEDIAN_SD,
           COUNT_IF(SD_LN < 0.01) / COUNT(*) AS SHARE_IDENTICAL
    FROM (
        SELECT HEALTH_SYSTEM_NAME, STDDEV(MEAN_LN_PRICE) AS SD_LN
        FROM SYSTEM_COUNTY
        GROUP BY HEALTH_SYSTEM_NAME, CANONICAL_PAYER_NAME, ANALYSIS_CONCEPT_ID
        HAVING COUNT(DISTINCT COUNTY_STATE) >= 3
    ) WHERE SD_LN IS NOT NULL
    GROUP BY HEALTH_SYSTEM_NAME
)
SELECT
    D.HEALTH_SYSTEM_NAME, S.N_HOSPITALS, S.N_COUNTIES, D.N_CELLS,
    ROUND(D.MEDIAN_SD, 4)       AS MEDIAN_SD_ACROSS_COUNTIES,
    ROUND(D.SHARE_IDENTICAL, 3) AS SHARE_ESSENTIALLY_IDENTICAL,
    ROUND(100 * S.N_HOSPITALS
          / (SELECT COUNT(DISTINCT HOSPITAL_ID) FROM HPT_TMP_SYS_PRICES), 2)
        AS PCT_OF_ALL_HOSPITALS
FROM DISP D INNER JOIN SYS_SIZE S USING (HEALTH_SYSTEM_NAME)
ORDER BY SHARE_ESSENTIALLY_IDENTICAL DESC;
-- A system above roughly 0.4 on SHARE_ESSENTIALLY_IDENTICAL is pricing
-- centrally. Read it against PCT_OF_ALL_HOSPITALS to judge whether it matters.


/* =============================================================================
   PART B2 -- FRANCHISE VERSUS INTEGRATED, STRATIFIED BY PAYER CLASS

   WHY THIS IS WORTH RUNNING beyond the pooled version above.

   The paper's transmission argument implies that compression should be UNEVEN
   across payer classes. Where rates are administratively anchored -- Medicare
   Advantage benchmarked to the Medicare fee schedule, managed Medicaid to
   state schedules -- there is little left for a system to centralise, because
   the anchor already does that work. Commercial contracts are where a system's
   own pricing policy has room to bind.

   PREDICTION: the within-system to within-market dispersion ratio should be
   LOWEST, meaning most centralised, in commercial, and HIGHEST, meaning most
   local, in the administratively anchored classes.

   If it comes back flat across classes, that is evidence the compression
   reflects shared cost structure and technology rather than deliberate central
   price-setting. That reading strengthens the exclusion argument, since shared
   inputs are absorbed by the county-by-concept fixed effects while a central
   pricing policy would not be.

   Requires Part A to have run.
   ============================================================================= */

CREATE OR REPLACE TEMPORARY TABLE HPT_TMP_SYS_PRICES_CLASS AS
SELECT
    P.HEALTH_SYSTEM_NAME,
    P.COUNTY_STATE,
    P.HOSPITAL_ID,
    P.CANONICAL_PAYER_NAME,
    C.PAYER_CLASS,
    P.ANALYSIS_FAMILY_ID,
    P.ANALYSIS_CONCEPT_ID,
    P.BILLING_CODE,
    MEDIAN(P.PAYER_CELL_RATE) AS MEDIAN_PRICE
FROM O_ECON.COMMON.HPT_P2_PAYER_CELLS P
INNER JOIN O_ECON.COMMON.HPT_P2_PAYER_CELL_CLASS C USING (PAYER_CELL_ID)
WHERE P.PAYER_CELL_RATE > 0
  AND P.BILLING_CODE_TYPE = 'CPT_HCPCS'
  AND P.CODE_ROLE = 'STANDALONE'
  AND P.HEALTH_SYSTEM_NAME IS NOT NULL
  AND P.MISSING_PAYER_IDENTITY_FLAG = 0
  AND C.PAYER_CLASS <> 'UNCLASSIFIED'
GROUP BY 1,2,3,4,5,6,7,8
;

SELECT PAYER_CLASS, COUNT(*) AS N_ROWS,
       COUNT(DISTINCT HEALTH_SYSTEM_NAME) AS N_SYSTEMS,
       COUNT(DISTINCT HOSPITAL_ID) AS N_HOSPITALS
FROM HPT_TMP_SYS_PRICES_CLASS GROUP BY 1 ORDER BY N_ROWS DESC;
-- Check the system count per class before reading anything below. A class with
-- few multi-county systems produces an unstable ratio.


WITH SYSTEM_COUNTY AS (
    SELECT HEALTH_SYSTEM_NAME, COUNTY_STATE, CANONICAL_PAYER_NAME, PAYER_CLASS,
           ANALYSIS_FAMILY_ID, ANALYSIS_CONCEPT_ID,
           AVG(LN(MEDIAN_PRICE)) AS MEAN_LN_PRICE
    FROM HPT_TMP_SYS_PRICES_CLASS
    GROUP BY 1,2,3,4,5,6
),
WS AS (
    SELECT PAYER_CLASS, COUNT(*) AS N_CELLS,
           MEDIAN(SD_LN) AS MEDIAN_SD, AVG(N_UNITS) AS MEAN_UNITS
    FROM (
        SELECT PAYER_CLASS, STDDEV(MEAN_LN_PRICE) AS SD_LN,
               COUNT(DISTINCT COUNTY_STATE) AS N_UNITS
        FROM SYSTEM_COUNTY
        GROUP BY HEALTH_SYSTEM_NAME, CANONICAL_PAYER_NAME, PAYER_CLASS,
                 ANALYSIS_FAMILY_ID, ANALYSIS_CONCEPT_ID
        HAVING COUNT(DISTINCT COUNTY_STATE) >= 3
    ) WHERE SD_LN IS NOT NULL AND SD_LN < 2
    GROUP BY PAYER_CLASS
),
WM AS (
    SELECT PAYER_CLASS, COUNT(*) AS N_CELLS,
           MEDIAN(SD_LN) AS MEDIAN_SD, AVG(N_UNITS) AS MEAN_UNITS
    FROM (
        SELECT PAYER_CLASS, STDDEV(LN(MEDIAN_PRICE)) AS SD_LN,
               COUNT(DISTINCT HOSPITAL_ID) AS N_UNITS
        FROM HPT_TMP_SYS_PRICES_CLASS
        GROUP BY COUNTY_STATE, CANONICAL_PAYER_NAME, PAYER_CLASS,
                 ANALYSIS_FAMILY_ID, ANALYSIS_CONCEPT_ID
        HAVING COUNT(DISTINCT HOSPITAL_ID) >= 3
    ) WHERE SD_LN IS NOT NULL AND SD_LN < 2
    GROUP BY PAYER_CLASS
)
SELECT
    COALESCE(S.PAYER_CLASS, M.PAYER_CLASS) AS PAYER_CLASS,
    S.N_CELLS AS N_CELLS_SYSTEM, ROUND(S.MEAN_UNITS, 2) AS MEAN_COUNTIES,
    M.N_CELLS AS N_CELLS_MARKET, ROUND(M.MEAN_UNITS, 2) AS MEAN_HOSPITALS,
    ROUND(S.MEDIAN_SD, 4) AS MEDIAN_SD_SYSTEM,
    ROUND(M.MEDIAN_SD, 4) AS MEDIAN_SD_MARKET,
    ROUND(S.MEDIAN_SD / NULLIF(M.MEDIAN_SD, 0), 3) AS RATIO_SYSTEM_TO_MARKET
FROM WS S FULL OUTER JOIN WM M USING (PAYER_CLASS)
ORDER BY RATIO_SYSTEM_TO_MARKET;
-- Compare MEAN_COUNTIES against MEAN_HOSPITALS within each class. Where they
-- diverge badly for one class, that class's ratio is not comparable to the
-- others and must be reported with its unit counts alongside.


/* Same comparison by class AND family. Use only if the class-level table looks
   stable; cells thin out quickly at this granularity, which is why both arms
   are held to a minimum of 30 cells. */

WITH SYSTEM_COUNTY AS (
    SELECT HEALTH_SYSTEM_NAME, COUNTY_STATE, CANONICAL_PAYER_NAME, PAYER_CLASS,
           ANALYSIS_FAMILY_ID, ANALYSIS_CONCEPT_ID,
           AVG(LN(MEDIAN_PRICE)) AS MEAN_LN_PRICE
    FROM HPT_TMP_SYS_PRICES_CLASS GROUP BY 1,2,3,4,5,6
),
WS AS (
    SELECT PAYER_CLASS, ANALYSIS_FAMILY_ID, COUNT(*) AS N_CELLS,
           MEDIAN(SD_LN) AS MEDIAN_SD
    FROM (
        SELECT PAYER_CLASS, ANALYSIS_FAMILY_ID, STDDEV(MEAN_LN_PRICE) AS SD_LN
        FROM SYSTEM_COUNTY
        GROUP BY HEALTH_SYSTEM_NAME, CANONICAL_PAYER_NAME, PAYER_CLASS,
                 ANALYSIS_FAMILY_ID, ANALYSIS_CONCEPT_ID
        HAVING COUNT(DISTINCT COUNTY_STATE) >= 3
    ) WHERE SD_LN IS NOT NULL AND SD_LN < 2
    GROUP BY 1,2
),
WM AS (
    SELECT PAYER_CLASS, ANALYSIS_FAMILY_ID, COUNT(*) AS N_CELLS,
           MEDIAN(SD_LN) AS MEDIAN_SD
    FROM (
        SELECT PAYER_CLASS, ANALYSIS_FAMILY_ID, STDDEV(LN(MEDIAN_PRICE)) AS SD_LN
        FROM HPT_TMP_SYS_PRICES_CLASS
        GROUP BY COUNTY_STATE, CANONICAL_PAYER_NAME, PAYER_CLASS,
                 ANALYSIS_FAMILY_ID, ANALYSIS_CONCEPT_ID
        HAVING COUNT(DISTINCT HOSPITAL_ID) >= 3
    ) WHERE SD_LN IS NOT NULL AND SD_LN < 2
    GROUP BY 1,2
)
SELECT
    COALESCE(S.PAYER_CLASS, M.PAYER_CLASS)             AS PAYER_CLASS,
    COALESCE(S.ANALYSIS_FAMILY_ID, M.ANALYSIS_FAMILY_ID) AS ANALYSIS_FAMILY_ID,
    S.N_CELLS AS N_CELLS_SYSTEM, M.N_CELLS AS N_CELLS_MARKET,
    ROUND(S.MEDIAN_SD, 4) AS MEDIAN_SD_SYSTEM,
    ROUND(M.MEDIAN_SD, 4) AS MEDIAN_SD_MARKET,
    ROUND(S.MEDIAN_SD / NULLIF(M.MEDIAN_SD, 0), 3) AS RATIO_SYSTEM_TO_MARKET
FROM WS S FULL OUTER JOIN WM M USING (PAYER_CLASS, ANALYSIS_FAMILY_ID)
WHERE S.N_CELLS >= 30 AND M.N_CELLS >= 30
ORDER BY PAYER_CLASS, RATIO_SYSTEM_TO_MARKET;


/* =============================================================================
   PART C -- COUNTY PAYER CONCENTRATION

   WHAT THIS IS FOR. The paper's theory section states a prediction it does not
   test: the disclosure effect should be amplified where insurers hold greater
   bargaining leverage, and market-level insurer concentration data are not
   available. This constructs a concentration measure from the disclosure files
   themselves.

   READ THE CAVEAT BEFORE USING IT AS THAT TEST. Contract-row share is not
   enrollment share. A payer with many negotiated line items in a county may
   simply have a broader contracted service list rather than more covered
   lives. The measure proxies negotiating breadth rather than market power, and
   it is mechanically correlated with how completely each hospital's file
   enumerates its contracts. If it is used to test the leverage prediction,
   that limitation belongs in the text and the result should be treated as
   suggestive rather than as a test of the structural parameter.

   Built from HPT_P2_PAYER_CELLS, which is already deduped and
   canonical-file-selected, so a hospital publishing overlapping files does not
   count twice.
   ============================================================================= */

CREATE OR REPLACE TABLE O_ECON.COMMON.HPT_COUNTY_PAYER_CONCENTRATION AS

WITH COUNTY_PAYER AS (
    SELECT
        COUNTY_STATE,
        CANONICAL_PAYER_NAME,
        COUNT(*)                        AS N_CONTRACT_ROWS,
        COUNT(DISTINCT HOSPITAL_ID)     AS N_HOSPITALS,
        COUNT(DISTINCT ANALYSIS_CONCEPT_ID) AS N_CONCEPTS
    FROM O_ECON.COMMON.HPT_P2_PAYER_CELLS
    WHERE PAYER_CELL_RATE > 0
      AND CODE_ROLE = 'STANDALONE'
      AND BILLING_CODE_TYPE = 'CPT_HCPCS'
      AND MISSING_PAYER_IDENTITY_FLAG = 0
      AND COUNTY_STATE IS NOT NULL
    GROUP BY 1, 2
),

SHARES AS (
    SELECT
        COUNTY_STATE,
        CANONICAL_PAYER_NAME,
        N_CONTRACT_ROWS,
        N_HOSPITALS,
        N_CONCEPTS,
        N_CONTRACT_ROWS / SUM(N_CONTRACT_ROWS) OVER (PARTITION BY COUNTY_STATE) AS SHARE,
        ROW_NUMBER() OVER (PARTITION BY COUNTY_STATE ORDER BY N_CONTRACT_ROWS DESC) AS RANK_IN_COUNTY
    FROM COUNTY_PAYER
)

SELECT
    COUNTY_STATE,
    CANONICAL_PAYER_NAME,
    N_CONTRACT_ROWS,
    N_HOSPITALS,
    N_CONCEPTS,
    ROUND(SHARE, 5) AS SHARE,
    RANK_IN_COUNTY,
    -- County-level aggregates repeated on every row, so the table merges onto
    -- the analysis panel without a second grouping step in R.
    COUNT(*)          OVER (PARTITION BY COUNTY_STATE) AS N_PAYERS_IN_COUNTY,
    ROUND(SUM(SHARE * SHARE) OVER (PARTITION BY COUNTY_STATE), 5) AS PAYER_HHI,
    ROUND(SUM(IFF(RANK_IN_COUNTY <= 3, SHARE, 0)) OVER (PARTITION BY COUNTY_STATE), 5) AS CR3
FROM SHARES
;

-- County-level summary, one row per county.
SELECT
    COUNTY_STATE,
    MAX(N_PAYERS_IN_COUNTY) AS N_PAYERS,
    MAX(PAYER_HHI)          AS PAYER_HHI,
    MAX(CR3)                AS CR3,
    MAX(IFF(RANK_IN_COUNTY = 1, CANONICAL_PAYER_NAME, NULL)) AS TOP_PAYER,
    MAX(IFF(RANK_IN_COUNTY = 1, SHARE, NULL))                AS TOP_PAYER_SHARE
FROM O_ECON.COMMON.HPT_COUNTY_PAYER_CONCENTRATION
GROUP BY COUNTY_STATE
ORDER BY PAYER_HHI DESC
LIMIT 40;

-- Distribution check. If HHI sits near 1 in most counties, the measure is
-- picking up counties with a single disclosing hospital rather than
-- concentrated insurance markets, and it should not be used.
SELECT
    CASE WHEN PAYER_HHI >= 0.50 THEN 'HHI >= 0.50 (very concentrated)'
         WHEN PAYER_HHI >= 0.25 THEN 'HHI 0.25-0.50'
         WHEN PAYER_HHI >= 0.15 THEN 'HHI 0.15-0.25'
         ELSE 'HHI < 0.15 (unconcentrated)' END AS HHI_BUCKET,
    COUNT(DISTINCT COUNTY_STATE) AS N_COUNTIES,
    ROUND(AVG(N_PAYERS_IN_COUNTY), 1) AS MEAN_PAYERS
FROM O_ECON.COMMON.HPT_COUNTY_PAYER_CONCENTRATION
GROUP BY 1 ORDER BY 1;

COPY INTO @O_ECON.COMMON.HPT_PY_EXPORT_STAGE/data/payer_concentration/HPT_COUNTY_PAYER_CONCENTRATION.csv.gz
FROM (SELECT * FROM O_ECON.COMMON.HPT_COUNTY_PAYER_CONCENTRATION
      ORDER BY COUNTY_STATE, RANK_IN_COUNTY)
FILE_FORMAT = (TYPE = CSV FIELD_OPTIONALLY_ENCLOSED_BY = '"' COMPRESSION = GZIP NULL_IF = (''))
HEADER = TRUE MAX_FILE_SIZE = 5368709120 OVERWRITE = TRUE
;

ALTER STAGE O_ECON.COMMON.HPT_PY_EXPORT_STAGE REFRESH;

SELECT RELATIVE_PATH, ROUND(SIZE/1024/1024, 2) AS SIZE_MB,
       GET_PRESIGNED_URL(@O_ECON.COMMON.HPT_PY_EXPORT_STAGE, RELATIVE_PATH, 604800) AS DOWNLOAD_URL
FROM DIRECTORY(@O_ECON.COMMON.HPT_PY_EXPORT_STAGE)
WHERE RELATIVE_PATH LIKE 'data/payer_concentration/%';
