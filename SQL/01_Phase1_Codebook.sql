/* =============================================================================

   WHEN TRANSPARENCY WORKS: SERVICE SHOPPABILITY, CONTRACTING DEPTH, AND
   THE PRICE EFFECTS OF HOSPITAL DISCLOSURE

   Replication code -- stage 1a: clinical codebook construction

   Danny Sierra
   Department of Economics, Florida State University
   Ds22c@fsu.edu

   -----------------------------------------------------------------------------
   WHAT THIS FILE DOES
   -----------------------------------------------------------------------------
   Builds the clinical codebook: which billing codes exist with enough support
   to study, whether each is a real and current code, and what clinical family
   and concept it belongs to. It is the first of three SQL/Python/R stages, and
   nothing downstream re-decides any of it.

   Four things are kept strictly separate, because collapsing them is what
   makes a codebook unreproducible:

     1. UNIVERSE       which codes exist with enough support   (automatic)
     2. VALIDATION     is the code real, current, an add-on,
                       ASC-eligible                            (automatic, from
                                                               official CMS files)
     3. FAMILY/CONCEPT what clinical category and concept a
                       code belongs to                         (a small, editable
                                                               rule table)
     4. SHOPPABILITY   deliberately not decided here          (see below)

   No step requires per-code manual review to run. Manual judgment enters in
   exactly one place, HPT_REF_MANUAL_OVERRIDES, and only for codes explicitly
   listed there.

   -----------------------------------------------------------------------------
   WHY SHOPPABILITY IS NOT ASSIGNED HERE
   -----------------------------------------------------------------------------
   The paper's central empirical exercise tests competing definitions of
   "shoppable" against each other. If this stage baked in one HIGH /
   INTERMEDIATE / LOW label, every alternative definition would require
   rebuilding the SQL pipeline, and the robustness exercise would be
   impossible in practice.

   Instead this stage exports the objective ingredients any shoppability
   definition needs -- OPPS status indicator, ASC eligibility, add-on and
   component flags, CMS-70 membership, MDC, DRG medical/surgical type, and an
   acute/emergent keyword flag -- as separate columns. Schemes are then built
   in R as functions of those columns, and eighteen of them are estimated
   without touching this file again.

   -----------------------------------------------------------------------------
   OFFICIAL REFERENCE FILES
   -----------------------------------------------------------------------------
   The study window (July 2024 through October 2025) spans six OPPS quarters
   and three MS-DRG fiscal years. All four sources are public CMS downloads.

     OPPS ADDENDUM B -- quarterly. CPT/HCPCS status indicators and APCs.
       https://www.cms.gov/medicare/medicare-fee-for-service-payment/hospitaloutpatientpps/addendum-a-and-addendum-b-updates
       One file per effective quarter: 2024-07-01, 2024-10-01, 2025-01-01,
       2025-04-01, 2025-07-01, 2025-10-01. Getting all six matters more here
       than for the other sources, because APC status and payment indicators
       change quarterly.

     MS-DRG WEIGHTS AND DEFINITIONS -- annual, effective each October 1.
       Published as Table 5 of the annual IPPS Final Rule.
       https://www.cms.gov/medicare/payment/prospective-payment-systems/acute-inpatient-pps/fy-2025-ipps-final-rule-home-page
       FY2024, FY2025, and FY2026 vintages cover Jul-Sep 2024, Oct 2024 -
       Sep 2025, and Oct 2025 respectively.

     NCCI ADD-ON CODE (AOC) EDITS -- updated quarterly, changes rarely.
       https://www.cms.gov/medicare/coding-billing/national-correct-coding-initiative-ncci-edits/medicare-ncci-add-code-edits
       Fixed-width text with an Excel version also posted. Lists every current
       add-on code and its accepted primary codes. This is the authoritative
       source for add-on status, and it replaces any description-based regex
       search for "each additional" text. One recent vintage is generally
       enough; for precision, load the vintages effective nearest 2024-07-01
       and 2025-10-01 and diff them.

     ASC COVERED PROCEDURES, ADDENDUM AA -- annual, refreshed with OPPS.
       https://www.cms.gov/medicare/medicare-fee-for-service-payment/ascpayment/asc-regulations-and-notices
       Under the CY2024 and CY2025 final rules, plus any mid-year corrections.
       Presence on this list is a CMS-determined signal that a procedure is
       plannable and schedulable, which makes it an objective proxy for
       shoppability that does not depend on the researcher's own judgment.

   Load all four into the staging tables in Section 1 with codes normalised to
   uppercase alphanumeric.

   -----------------------------------------------------------------------------
   REQUIRED INPUT
   -----------------------------------------------------------------------------
     O_ECON.COMMON.HR_SUBSET_ALL      raw disclosure source
     O_ECON.COMMON.HPT_CMS70_CODEBOOK CMS-70 required-service mapping

   -----------------------------------------------------------------------------
   OUTPUTS
   -----------------------------------------------------------------------------
   Reference tables:
     HPT_REF_OPPS_ADDENDUM_B, HPT_REF_MSDRG        loaded separately
     HPT_REF_NCCI_ADDON, HPT_REF_ASC_COVERED       loaded in Section 1
     HPT_REF_FAMILY_RULES                          the editable rule table
     HPT_REF_DRG_ACUITY_KEYWORDS
     HPT_REF_MANUAL_OVERRIDES                      empty by default

   Pipeline tables:
     HPT_P1_UNIVERSE_SUPPORT          candidate codes and their support
     HPT_P1_VALIDATED                 + official validation flags
     HPT_P1_CLASSIFIED                + family and concept assignment
     HPT_P1_FINAL_CODEBOOK            + manual overrides, full archive
     HPT_P1_FINAL_CODEBOOK_SCOPED     what actually feeds stage 1b
     HPT_P1_QA_SUMMARY

   -----------------------------------------------------------------------------
   HOW TO RUN
   -----------------------------------------------------------------------------
   Load the four reference files into HPT_REF_STAGE, then run Sections 1
   through 6 in order. Section 2 contains the single expensive scan of the raw
   source; check its QA output before continuing. Each later section ends with
   its own QA query.
   ============================================================================= */


USE WAREHOUSE ECON_DEF_WH;
USE DATABASE O_ECON;
USE SCHEMA COMMON;


/* =============================================================================
   SECTION 0 -- STUDY WINDOW CONSTANTS

   Defined once and cross-joined wherever the window is needed, so changing the
   window is a one-line edit rather than a search-and-replace.
   ============================================================================= */

CREATE OR REPLACE TABLE O_ECON.COMMON.HPT_P1_STUDY_WINDOW AS
SELECT
    '2024-07-01'::DATE AS WINDOW_START_DATE,
    '2025-10-31'::DATE AS WINDOW_END_DATE
;


/* =============================================================================
   SECTION 1 -- REFERENCE TABLES

   1A. OPPS Addendum B and MS-DRG are loaded separately into
       HPT_REF_OPPS_ADDENDUM_B and HPT_REF_MSDRG: quarterly and annual official
       vintages, codes normalised. Those loaders are unchanged by anything in
       this file.

   1B. NCCI Add-on Code (AOC) file.
   ============================================================================= */

CREATE STAGE IF NOT EXISTS O_ECON.COMMON.HPT_REF_STAGE
    DIRECTORY = (ENABLE = TRUE)
    COMMENT = 'Upload OPPS, MS-DRG, NCCI AOC, and ASC Addendum AA source files here'
;

CREATE OR REPLACE FILE FORMAT O_ECON.COMMON.HPT_REF_CSV_FORMAT
    TYPE = CSV
    SKIP_HEADER = 1
    FIELD_OPTIONALLY_ENCLOSED_BY = '"'
    TRIM_SPACE = TRUE
    EMPTY_FIELD_AS_NULL = TRUE
    NULL_IF = ('', 'NULL', 'null', 'NA', 'N/A')
    ERROR_ON_COLUMN_COUNT_MISMATCH = FALSE
    ENCODING = 'UTF8'
;

/* Upload the NCCI AOC export to the stage as ncci_addon_2024_2025.csv with
   columns in this order: EFFECTIVE_DATE, ADDON_CODE, ADDON_TYPE,
   PRIMARY_CODE_LIST (comma or pipe delimited where a Type 1 AOC lists
   several primaries). */

CREATE OR REPLACE TABLE O_ECON.COMMON.HPT_REF_NCCI_ADDON (
    EFFECTIVE_DATE      DATE,
    ADDON_CODE          VARCHAR,
    ADDON_TYPE          VARCHAR,     -- TYPE_1, TYPE_2, or TYPE_3
    PRIMARY_CODE_LIST   VARCHAR,     -- raw text; parsed below
    SOURCE_FILE         VARCHAR,
    LOADED_AT           TIMESTAMP_NTZ DEFAULT CURRENT_TIMESTAMP()
);

COPY INTO O_ECON.COMMON.HPT_REF_NCCI_ADDON (
    EFFECTIVE_DATE, ADDON_CODE, ADDON_TYPE, PRIMARY_CODE_LIST, SOURCE_FILE
)
FROM (
    SELECT
        TRY_TO_DATE($1),
        REGEXP_REPLACE(TRIM($2), '[^A-Z0-9]', ''),
        TRIM($3),
        TRIM($4),
        METADATA$FILENAME
    FROM @O_ECON.COMMON.HPT_REF_STAGE (
        FILE_FORMAT => 'O_ECON.COMMON.HPT_REF_CSV_FORMAT',
        PATTERN => '.*ncci_addon.*[.]csv([.]gz)?'
    )
)
ON_ERROR = 'CONTINUE'
;

UPDATE O_ECON.COMMON.HPT_REF_NCCI_ADDON SET
    ADDON_CODE = UPPER(ADDON_CODE),
    ADDON_TYPE = UPPER(ADDON_TYPE)
;

/* One row per add-on code, deduplicated across whatever vintages were loaded. */
CREATE OR REPLACE TABLE O_ECON.COMMON.HPT_REF_NCCI_ADDON_CURRENT AS
SELECT
    ADDON_CODE,
    MAX(EFFECTIVE_DATE)                         AS LATEST_EFFECTIVE_DATE,
    LISTAGG(DISTINCT ADDON_TYPE, ',')
        WITHIN GROUP (ORDER BY ADDON_TYPE)      AS ADDON_TYPES,
    LISTAGG(DISTINCT PRIMARY_CODE_LIST, ' | ')
        WITHIN GROUP (ORDER BY PRIMARY_CODE_LIST) AS PRIMARY_CODE_LISTS
FROM O_ECON.COMMON.HPT_REF_NCCI_ADDON
WHERE ADDON_CODE IS NOT NULL
GROUP BY ADDON_CODE
;


/* =============================================================================
   1C. ASC Covered Procedures (Addendum AA).

   Upload as asc_addendum_aa_2024_2025.csv with columns: REFERENCE_YEAR,
   EFFECTIVE_DATE, BILLING_CODE, SHORT_DESCRIPTION, PAYMENT_INDICATOR,
   SOURCE_FILE.
   ============================================================================= */

CREATE OR REPLACE TABLE O_ECON.COMMON.HPT_REF_ASC_COVERED (
    REFERENCE_YEAR      NUMBER,
    EFFECTIVE_DATE      DATE,
    BILLING_CODE        VARCHAR,
    SHORT_DESCRIPTION   VARCHAR,
    PAYMENT_INDICATOR   VARCHAR,
    SOURCE_FILE         VARCHAR,
    LOADED_AT           TIMESTAMP_NTZ DEFAULT CURRENT_TIMESTAMP()
);

COPY INTO O_ECON.COMMON.HPT_REF_ASC_COVERED (
    REFERENCE_YEAR, EFFECTIVE_DATE, BILLING_CODE, SHORT_DESCRIPTION,
    PAYMENT_INDICATOR, SOURCE_FILE
)
FROM (
    SELECT
        TRY_TO_NUMBER($1),
        TRY_TO_DATE($2),
        REGEXP_REPLACE(TRIM($3), '[^A-Z0-9]', ''),
        TRIM($4),
        TRIM($5),
        METADATA$FILENAME
    FROM @O_ECON.COMMON.HPT_REF_STAGE (
        FILE_FORMAT => 'O_ECON.COMMON.HPT_REF_CSV_FORMAT',
        PATTERN => '.*asc_addendum_aa.*[.]csv([.]gz)?'
    )
)
ON_ERROR = 'CONTINUE'
;

UPDATE O_ECON.COMMON.HPT_REF_ASC_COVERED SET
    BILLING_CODE = UPPER(BILLING_CODE),
    PAYMENT_INDICATOR = UPPER(PAYMENT_INDICATOR)
;

CREATE OR REPLACE TABLE O_ECON.COMMON.HPT_REF_ASC_COVERED_CURRENT AS
SELECT
    BILLING_CODE,
    MAX(EFFECTIVE_DATE) AS LATEST_EFFECTIVE_DATE,
    MAX_BY(SHORT_DESCRIPTION, EFFECTIVE_DATE) AS LATEST_DESCRIPTION,
    LISTAGG(DISTINCT PAYMENT_INDICATOR, ',')
        WITHIN GROUP (ORDER BY PAYMENT_INDICATOR) AS PAYMENT_INDICATORS
FROM O_ECON.COMMON.HPT_REF_ASC_COVERED
WHERE BILLING_CODE IS NOT NULL
GROUP BY BILLING_CODE
;



/* =============================================================================
   1D. FAMILY RULE TABLE

   Every row is a RULE rather than a code-level judgment call. Adding or fixing
   a family assignment is a one-row edit here, not a change to the SQL logic
   below, which is what makes the classification auditable and reproducible.

   RULE_TYPE:
     EXACT         matches one specific billing code. Used only for genuine
                   exceptions where a range or keyword rule would misfire --
                   CT-based bone-density codes that sit outside the general CT
                   numeric block, or HCPCS Level II codes that do not follow
                   CPT numbering at all.
     RANGE         matches a numeric CPT/HCPCS band. CPT is organised by
                   anatomic region within each broad section, so most imaging
                   families reduce to a manageable number of contiguous bands.
     KEYWORD       matches a regex against the OFFICIAL description, from OPPS
                   or MS-DRG rather than the free-text hospital description.
                   Needed because CPT numbering mixes modalities within one
                   range across most of the 70000-79999 radiology section, so a
                   range rule alone would misclassify.
     CODE_PATTERN  matches a regex against the BILLING CODE itself. Used for
                   whole classes of codes that are never independently priced
                   services regardless of what their description says.

   RESOLUTION: for a given code, the matching rule with the LOWEST PRIORITY
   number wins. Exclusions and EXACT overrides sit at low priority numbers,
   RANGE rules next, KEYWORD rules after that, and a coarse CPT-section
   fallback last so that nothing falls through unclassified.

   A NOTE ON PRIORITY 0 EXCLUSIONS. These match on code pattern rather than on
   description, which matters: a code whose description happens to contain a
   clinical keyword ("Anesth biopsy of nose") would otherwise be pulled into a
   clinical family it does not belong to once keyword matching is working.
   Pattern-based rules generalise to codes never individually inspected, which
   code-by-code exclusion lists do not.
   ============================================================================= */

CREATE OR REPLACE TABLE O_ECON.COMMON.HPT_REF_FAMILY_RULES (
    PRIORITY            NUMBER,
    RULE_TYPE           VARCHAR,      -- EXACT | RANGE | KEYWORD | CODE_PATTERN
    BILLING_CODE_TYPE   VARCHAR,      -- CPT_HCPCS | MS_DRG
    MATCH_CODE          VARCHAR,      -- for EXACT
    RANGE_START         VARCHAR,      -- for RANGE (numeric part only)
    RANGE_END           VARCHAR,      -- for RANGE
    KEYWORD_PATTERN     VARCHAR,      -- for KEYWORD and CODE_PATTERN
    ANALYSIS_SUPERFAMILY_ID VARCHAR,
    ANALYSIS_FAMILY_ID  VARCHAR,
    ANALYSIS_FAMILY_NAME VARCHAR,
    RULE_NOTE           VARCHAR
);

INSERT INTO O_ECON.COMMON.HPT_REF_FAMILY_RULES
    (PRIORITY, RULE_TYPE, BILLING_CODE_TYPE, MATCH_CODE, RANGE_START, RANGE_END,
     KEYWORD_PATTERN, ANALYSIS_SUPERFAMILY_ID, ANALYSIS_FAMILY_ID,
     ANALYSIS_FAMILY_NAME, RULE_NOTE)
VALUES
    /* ---- Priority 0: exclusions by code pattern. These run first so that a
       code sharing a clinical keyword with a target family never enters it. ---- */
    (0, 'CODE_PATTERN', 'CPT_HCPCS', NULL, NULL, NULL, '^[0-9]{4}F$',
        'OTHER', 'EXCLUDED_QUALITY_MEASURE', 'Excluded: Quality Measure Tracking Code',
        'CPT Category II codes exist purely for quality reporting; they are never real independently-priced services even if the raw data shows a nonzero dollar figure.'),
    (0, 'CODE_PATTERN', 'CPT_HCPCS', NULL, NULL, NULL, '^G[89][0-9]{3}$',
        'OTHER', 'EXCLUDED_QUALITY_MEASURE', 'Excluded: Quality Measure Tracking Code',
        'HCPCS G8000-G9999 codes are legacy PQRS/MIPS quality-measure documentation codes, distinct from real billable G0xxx-G7xxx service codes (e.g. G0105, G0121, G0279, which are unaffected by this rule).'),
    (0, 'CODE_PATTERN', 'CPT_HCPCS', NULL, NULL, NULL, '^D[0-9]{4}$',
        'OTHER', 'EXCLUDED_LIKELY_MISLABELED', 'Excluded: Likely Mislabeled Dental (CDT) Code',
        'CDT dental codes should never appear under a CPT/HCPCS billing_code_type; this pattern almost certainly reflects a hospital-side data-labeling error, not a real CPT/HCPCS service.'),

    /* ---- Priority 0: anesthesia is a distinct service line from the
       procedure it accompanies. Every CPT anesthesia code sits in the
       00100-01999 range, which makes this a clean range rule rather than a
       description-by-description exercise. Classified rather than dropped so
       it stays visible in the full archive, but excluded from the target
       families by the scope logic in Section 5. ---- */
    (0, 'RANGE', 'CPT_HCPCS', NULL, '00100', '01999', NULL,
        'OTHER', 'ANESTHESIA', 'Anesthesia',
        'Anesthesia service codes; a distinct billable service from the procedure they accompany, even when the description mentions the procedure by name.'),

    /* ---- Priority 1: EXACT overrides, which must beat range and keyword rules ---- */
    (1, 'EXACT', 'CPT_HCPCS', '57520', NULL, NULL, NULL,
        'OTHER', 'GYNECOLOGIC_PROCEDURE', 'Gynecologic Procedure',
        'Cervical conization; not a simple biopsy despite numeric proximity to biopsy codes.'),
    (1, 'EXACT', 'CPT_HCPCS', '76942', NULL, NULL, NULL,
        'OTHER', 'IMAGING_GUIDANCE_COMPONENT', 'Imaging Guidance (Component)',
        'Ultrasound guidance for needle placement; always a component of a parent procedure.'),
    (1, 'EXACT', 'CPT_HCPCS', '77002', NULL, NULL, NULL,
        'OTHER', 'IMAGING_GUIDANCE_COMPONENT', 'Imaging Guidance (Component)', 'Fluoroscopic guidance component.'),
    (1, 'EXACT', 'CPT_HCPCS', '77011', NULL, NULL, NULL,
        'OTHER', 'IMAGING_GUIDANCE_COMPONENT', 'Imaging Guidance (Component)', 'CT guidance component.'),
    (1, 'EXACT', 'CPT_HCPCS', '77012', NULL, NULL, NULL,
        'OTHER', 'IMAGING_GUIDANCE_COMPONENT', 'Imaging Guidance (Component)', 'CT guidance component.'),
    (1, 'EXACT', 'CPT_HCPCS', '77021', NULL, NULL, NULL,
        'OTHER', 'IMAGING_GUIDANCE_COMPONENT', 'Imaging Guidance (Component)', 'MRI guidance component.'),
    (1, 'EXACT', 'CPT_HCPCS', 'C8900', NULL, NULL, NULL,
        'IMAGING', 'MRI_MRA', 'MRI/MRA', 'HCPCS Level II facility MRA abdomen; outside CPT MRI numeric range.'),
    (1, 'EXACT', 'CPT_HCPCS', 'C8901', NULL, NULL, NULL,
        'IMAGING', 'MRI_MRA', 'MRI/MRA', 'HCPCS Level II facility MRA abdomen.'),
    (1, 'EXACT', 'CPT_HCPCS', 'C8902', NULL, NULL, NULL,
        'IMAGING', 'MRI_MRA', 'MRI/MRA', 'HCPCS Level II facility MRA abdomen.'),
    (1, 'EXACT', 'CPT_HCPCS', 'G0105', NULL, NULL, NULL,
        'GI_ENDOSCOPY', 'COLONOSCOPY_LOWER_ENDOSCOPY', 'Colonoscopy / Lower Endoscopy', 'Medicare high-risk screening colonoscopy G-code.'),
    (1, 'EXACT', 'CPT_HCPCS', 'G0121', NULL, NULL, NULL,
        'GI_ENDOSCOPY', 'COLONOSCOPY_LOWER_ENDOSCOPY', 'Colonoscopy / Lower Endoscopy', 'Medicare average-risk screening colonoscopy G-code.'),
    (1, 'EXACT', 'CPT_HCPCS', 'G0279', NULL, NULL, NULL,
        'IMAGING', 'MAMMOGRAPHY', 'Mammography', 'Screening tomosynthesis add-on G-code.'),

    /* Biopsy codes are listed as EXACT rules rather than left to the keyword
       fallback below. CMS official short descriptors abbreviate biopsy as
       "bx" inconsistently, so keyword matching alone is not reliable for a
       family this central to the analysis. */
    (1, 'EXACT', 'CPT_HCPCS', '11102', NULL, NULL, NULL, 'BIOPSY', 'BIOPSY', 'Biopsy', 'Skin biopsy, tangential.'),
    (1, 'EXACT', 'CPT_HCPCS', '11103', NULL, NULL, NULL, 'BIOPSY', 'BIOPSY', 'Biopsy', 'Skin biopsy, tangential, additional (add-on).'),
    (1, 'EXACT', 'CPT_HCPCS', '11104', NULL, NULL, NULL, 'BIOPSY', 'BIOPSY', 'Biopsy', 'Skin biopsy, punch.'),
    (1, 'EXACT', 'CPT_HCPCS', '11105', NULL, NULL, NULL, 'BIOPSY', 'BIOPSY', 'Biopsy', 'Skin biopsy, punch, additional (add-on).'),
    (1, 'EXACT', 'CPT_HCPCS', '11106', NULL, NULL, NULL, 'BIOPSY', 'BIOPSY', 'Biopsy', 'Skin biopsy, incisional.'),
    (1, 'EXACT', 'CPT_HCPCS', '11107', NULL, NULL, NULL, 'BIOPSY', 'BIOPSY', 'Biopsy', 'Skin biopsy, incisional, additional (add-on).'),
    (1, 'EXACT', 'CPT_HCPCS', '19081', NULL, NULL, NULL, 'BIOPSY', 'BIOPSY', 'Biopsy', 'Breast biopsy, stereotactic.'),
    (1, 'EXACT', 'CPT_HCPCS', '19082', NULL, NULL, NULL, 'BIOPSY', 'BIOPSY', 'Biopsy', 'Breast biopsy, stereotactic, additional (add-on).'),
    (1, 'EXACT', 'CPT_HCPCS', '19083', NULL, NULL, NULL, 'BIOPSY', 'BIOPSY', 'Biopsy', 'Breast biopsy, ultrasound-guided.'),
    (1, 'EXACT', 'CPT_HCPCS', '19084', NULL, NULL, NULL, 'BIOPSY', 'BIOPSY', 'Biopsy', 'Breast biopsy, ultrasound-guided, additional (add-on).'),
    (1, 'EXACT', 'CPT_HCPCS', '19085', NULL, NULL, NULL, 'BIOPSY', 'BIOPSY', 'Biopsy', 'Breast biopsy, MRI-guided.'),
    (1, 'EXACT', 'CPT_HCPCS', '19086', NULL, NULL, NULL, 'BIOPSY', 'BIOPSY', 'Biopsy', 'Breast biopsy, MRI-guided, additional (add-on).'),
    (1, 'EXACT', 'CPT_HCPCS', '20200', NULL, NULL, NULL, 'BIOPSY', 'BIOPSY', 'Biopsy', 'Muscle biopsy.'),
    (1, 'EXACT', 'CPT_HCPCS', '20205', NULL, NULL, NULL, 'BIOPSY', 'BIOPSY', 'Biopsy', 'Muscle biopsy, deep.'),
    (1, 'EXACT', 'CPT_HCPCS', '20220', NULL, NULL, NULL, 'BIOPSY', 'BIOPSY', 'Biopsy', 'Bone biopsy.'),
    (1, 'EXACT', 'CPT_HCPCS', '20225', NULL, NULL, NULL, 'BIOPSY', 'BIOPSY', 'Biopsy', 'Bone biopsy, deep.'),
    (1, 'EXACT', 'CPT_HCPCS', '32408', NULL, NULL, NULL, 'BIOPSY', 'BIOPSY', 'Biopsy', 'Lung/mediastinum biopsy.'),
    (1, 'EXACT', 'CPT_HCPCS', '38220', NULL, NULL, NULL, 'BIOPSY', 'BIOPSY', 'Biopsy', 'Bone marrow biopsy.'),
    (1, 'EXACT', 'CPT_HCPCS', '38221', NULL, NULL, NULL, 'BIOPSY', 'BIOPSY', 'Biopsy', 'Bone marrow biopsy.'),
    (1, 'EXACT', 'CPT_HCPCS', '38222', NULL, NULL, NULL, 'BIOPSY', 'BIOPSY', 'Biopsy', 'Bone marrow biopsy, combined.'),
    (1, 'EXACT', 'CPT_HCPCS', '38500', NULL, NULL, NULL, 'BIOPSY', 'BIOPSY', 'Biopsy', 'Lymph node biopsy.'),
    (1, 'EXACT', 'CPT_HCPCS', '38505', NULL, NULL, NULL, 'BIOPSY', 'BIOPSY', 'Biopsy', 'Lymph node biopsy, needle.'),
    (1, 'EXACT', 'CPT_HCPCS', '38510', NULL, NULL, NULL, 'BIOPSY', 'BIOPSY', 'Biopsy', 'Lymph node biopsy, open.'),
    (1, 'EXACT', 'CPT_HCPCS', '47000', NULL, NULL, NULL, 'BIOPSY', 'BIOPSY', 'Biopsy', 'Liver biopsy.'),
    (1, 'EXACT', 'CPT_HCPCS', '49180', NULL, NULL, NULL, 'BIOPSY', 'BIOPSY', 'Biopsy', 'Abdominal/retroperitoneal mass biopsy.'),
    (1, 'EXACT', 'CPT_HCPCS', '55700', NULL, NULL, NULL, 'BIOPSY', 'BIOPSY', 'Biopsy', 'Prostate biopsy.'),
    (1, 'EXACT', 'CPT_HCPCS', '57500', NULL, NULL, NULL, 'BIOPSY', 'BIOPSY', 'Biopsy', 'Cervical biopsy.'),
    (1, 'EXACT', 'CPT_HCPCS', '60100', NULL, NULL, NULL, 'BIOPSY', 'BIOPSY', 'Biopsy', 'Thyroid biopsy.'),

    /* ---- Priority 10: RANGE rules over clean anatomic and modality blocks ---- */
    (10, 'RANGE', 'CPT_HCPCS', NULL, '70450', '70470', NULL, 'IMAGING', 'CT_CTA', 'CT/CTA', 'Head/brain CT.'),
    (10, 'RANGE', 'CPT_HCPCS', NULL, '70480', '70492', NULL, 'IMAGING', 'CT_CTA', 'CT/CTA', 'Orbit/sella/maxillofacial/neck CT.'),
    (10, 'RANGE', 'CPT_HCPCS', NULL, '70496', '70498', NULL, 'IMAGING', 'CT_CTA', 'CT/CTA', 'CTA head/neck.'),
    (10, 'RANGE', 'CPT_HCPCS', NULL, '71250', '71275', NULL, 'IMAGING', 'CT_CTA', 'CT/CTA', 'Chest CT/CTA.'),
    (10, 'RANGE', 'CPT_HCPCS', NULL, '72125', '72133', NULL, 'IMAGING', 'CT_CTA', 'CT/CTA', 'Spine CT.'),
    (10, 'RANGE', 'CPT_HCPCS', NULL, '72191', '72194', NULL, 'IMAGING', 'CT_CTA', 'CT/CTA', 'Pelvis CT/CTA.'),
    (10, 'RANGE', 'CPT_HCPCS', NULL, '73200', '73206', NULL, 'IMAGING', 'CT_CTA', 'CT/CTA', 'Upper extremity CT/CTA.'),
    (10, 'RANGE', 'CPT_HCPCS', NULL, '73700', '73706', NULL, 'IMAGING', 'CT_CTA', 'CT/CTA', 'Lower extremity CT/CTA.'),
    (10, 'RANGE', 'CPT_HCPCS', NULL, '74150', '74178', NULL, 'IMAGING', 'CT_CTA', 'CT/CTA', 'Abdomen/pelvis CT/CTA.'),
    (10, 'RANGE', 'CPT_HCPCS', NULL, '74261', '74263', NULL, 'IMAGING', 'CT_CTA', 'CT/CTA', 'CT colonography.'),
    (10, 'RANGE', 'CPT_HCPCS', NULL, '75571', '75574', NULL, 'IMAGING', 'CT_CTA', 'CT/CTA', 'Cardiac CT/CTA.'),
    (10, 'RANGE', 'CPT_HCPCS', NULL, '76380', '76380', NULL, 'IMAGING', 'CT_CTA', 'CT/CTA', 'Limited/localized CT.'),

    (10, 'RANGE', 'CPT_HCPCS', NULL, '70336', '70336', NULL, 'IMAGING', 'MRI_MRA', 'MRI/MRA', 'TMJ MRI.'),
    (10, 'RANGE', 'CPT_HCPCS', NULL, '70540', '70549', NULL, 'IMAGING', 'MRI_MRA', 'MRI/MRA', 'Orbit/face/neck MRI/MRA.'),
    (10, 'RANGE', 'CPT_HCPCS', NULL, '70551', '70559', NULL, 'IMAGING', 'MRI_MRA', 'MRI/MRA', 'Brain MRI.'),
    (10, 'RANGE', 'CPT_HCPCS', NULL, '71550', '71555', NULL, 'IMAGING', 'MRI_MRA', 'MRI/MRA', 'Chest MRI/MRA.'),
    (10, 'RANGE', 'CPT_HCPCS', NULL, '72141', '72159', NULL, 'IMAGING', 'MRI_MRA', 'MRI/MRA', 'Spine MRI/MRA.'),
    (10, 'RANGE', 'CPT_HCPCS', NULL, '72195', '72198', NULL, 'IMAGING', 'MRI_MRA', 'MRI/MRA', 'Pelvis MRI/MRA.'),
    (10, 'RANGE', 'CPT_HCPCS', NULL, '73218', '73225', NULL, 'IMAGING', 'MRI_MRA', 'MRI/MRA', 'Upper extremity MRI/MRA.'),
    (10, 'RANGE', 'CPT_HCPCS', NULL, '73718', '73725', NULL, 'IMAGING', 'MRI_MRA', 'MRI/MRA', 'Lower extremity MRI/MRA.'),
    (10, 'RANGE', 'CPT_HCPCS', NULL, '74181', '74185', NULL, 'IMAGING', 'MRI_MRA', 'MRI/MRA', 'Abdomen MRI/MRA.'),
    (10, 'RANGE', 'CPT_HCPCS', NULL, '75557', '75565', NULL, 'IMAGING', 'MRI_MRA', 'MRI/MRA', 'Cardiac MRI.'),
    (10, 'RANGE', 'CPT_HCPCS', NULL, '76390', '76390', NULL, 'IMAGING', 'MRI_MRA', 'MRI/MRA', 'MR spectroscopy.'),
    (10, 'RANGE', 'CPT_HCPCS', NULL, '77046', '77049', NULL, 'IMAGING', 'MRI_MRA', 'MRI/MRA', 'Breast MRI.'),
    (10, 'RANGE', 'CPT_HCPCS', NULL, '77084', '77084', NULL, 'IMAGING', 'MRI_MRA', 'MRI/MRA', 'Whole-body/bone-marrow MRI.'),

    (10, 'RANGE', 'CPT_HCPCS', NULL, '70100', '70170', NULL, 'IMAGING', 'XRAY_FLUOROSCOPY', 'X-Ray/Fluoroscopy', 'Mandible/facial X-ray.'),
    (10, 'RANGE', 'CPT_HCPCS', NULL, '70190', '70260', NULL, 'IMAGING', 'XRAY_FLUOROSCOPY', 'X-Ray/Fluoroscopy', 'Sinus/orbit/skull X-ray.'),
    (10, 'RANGE', 'CPT_HCPCS', NULL, '70328', '70392', NULL, 'IMAGING', 'XRAY_FLUOROSCOPY', 'X-Ray/Fluoroscopy', 'TMJ/soft tissue neck/sialography.'),
    (10, 'RANGE', 'CPT_HCPCS', NULL, '71045', '71130', NULL, 'IMAGING', 'XRAY_FLUOROSCOPY', 'X-Ray/Fluoroscopy', 'Chest/rib/sternum X-ray.'),
    (10, 'RANGE', 'CPT_HCPCS', NULL, '72020', '72084', NULL, 'IMAGING', 'XRAY_FLUOROSCOPY', 'X-Ray/Fluoroscopy', 'Spine X-ray (survey/cervical/thoracic/scoliosis).'),
    (10, 'RANGE', 'CPT_HCPCS', NULL, '72100', '72120', NULL, 'IMAGING', 'XRAY_FLUOROSCOPY', 'X-Ray/Fluoroscopy', 'Lumbar spine X-ray.'),
    (10, 'RANGE', 'CPT_HCPCS', NULL, '72170', '72220', NULL, 'IMAGING', 'XRAY_FLUOROSCOPY', 'X-Ray/Fluoroscopy', 'Pelvis/sacrum X-ray.'),
    (10, 'RANGE', 'CPT_HCPCS', NULL, '73000', '73140', NULL, 'IMAGING', 'XRAY_FLUOROSCOPY', 'X-Ray/Fluoroscopy', 'Upper extremity X-ray.'),
    (10, 'RANGE', 'CPT_HCPCS', NULL, '73501', '73660', NULL, 'IMAGING', 'XRAY_FLUOROSCOPY', 'X-Ray/Fluoroscopy', 'Lower extremity/hip X-ray.'),
    (10, 'RANGE', 'CPT_HCPCS', NULL, '74018', '74022', NULL, 'IMAGING', 'XRAY_FLUOROSCOPY', 'X-Ray/Fluoroscopy', 'Abdomen X-ray.'),
    (10, 'RANGE', 'CPT_HCPCS', NULL, '74190', '74425', NULL, 'IMAGING', 'XRAY_FLUOROSCOPY', 'X-Ray/Fluoroscopy', 'GI/GU contrast fluoroscopy studies.'),

    (10, 'RANGE', 'CPT_HCPCS', NULL, '76506', '76529', NULL, 'IMAGING', 'DIAGNOSTIC_ULTRASOUND', 'Diagnostic Ultrasound', 'Ophthalmic/head ultrasound.'),
    (10, 'RANGE', 'CPT_HCPCS', NULL, '76536', '76536', NULL, 'IMAGING', 'DIAGNOSTIC_ULTRASOUND', 'Diagnostic Ultrasound', 'Thyroid/neck ultrasound.'),
    (10, 'RANGE', 'CPT_HCPCS', NULL, '76604', '76645', NULL, 'IMAGING', 'DIAGNOSTIC_ULTRASOUND', 'Diagnostic Ultrasound', 'Chest/breast ultrasound.'),
    (10, 'RANGE', 'CPT_HCPCS', NULL, '76700', '76776', NULL, 'IMAGING', 'DIAGNOSTIC_ULTRASOUND', 'Diagnostic Ultrasound', 'Abdomen/retroperitoneal ultrasound.'),
    (10, 'RANGE', 'CPT_HCPCS', NULL, '76801', '76828', NULL, 'IMAGING', 'DIAGNOSTIC_ULTRASOUND', 'Diagnostic Ultrasound', 'Obstetric ultrasound.'),
    (10, 'RANGE', 'CPT_HCPCS', NULL, '76830', '76873', NULL, 'IMAGING', 'DIAGNOSTIC_ULTRASOUND', 'Diagnostic Ultrasound', 'Pelvis/prostate/scrotum ultrasound.'),
    (10, 'RANGE', 'CPT_HCPCS', NULL, '76881', '76886', NULL, 'IMAGING', 'DIAGNOSTIC_ULTRASOUND', 'Diagnostic Ultrasound', 'Extremity/infant hip ultrasound.'),
    (10, 'RANGE', 'CPT_HCPCS', NULL, '76977', '76999', NULL, 'IMAGING', 'DIAGNOSTIC_ULTRASOUND', 'Diagnostic Ultrasound', 'Bone-density/elastography/guidance/specialized ultrasound.'),

    (10, 'RANGE', 'CPT_HCPCS', NULL, '93303', '93356', NULL, 'IMAGING', 'ECHOCARDIOGRAPHY', 'Echocardiography', 'Transthoracic/transesophageal echo and components.'),
    (10, 'RANGE', 'CPT_HCPCS', NULL, '93880', '93998', NULL, 'IMAGING', 'VASCULAR_ULTRASOUND', 'Vascular Ultrasound', 'Duplex/vascular ultrasound studies.'),
    (10, 'RANGE', 'CPT_HCPCS', NULL, '77061', '77067', NULL, 'IMAGING', 'MAMMOGRAPHY', 'Mammography', 'Diagnostic/screening mammography and tomosynthesis.'),
    (10, 'RANGE', 'CPT_HCPCS', NULL, '77078', '77086', NULL, 'IMAGING', 'BONE_DENSITY', 'Bone Density', 'CT and DXA bone-density studies (includes 77078/77079 CT-based bone density).'),

    (10, 'RANGE', 'CPT_HCPCS', NULL, '43191', '43273', NULL, 'GI_ENDOSCOPY', 'UPPER_ENDOSCOPY', 'Upper Endoscopy', 'Diagnostic and therapeutic upper GI endoscopy.'),
    (10, 'RANGE', 'CPT_HCPCS', NULL, '45300', '45398', NULL, 'GI_ENDOSCOPY', 'COLONOSCOPY_LOWER_ENDOSCOPY', 'Colonoscopy / Lower Endoscopy', 'Sigmoidoscopy/colonoscopy.'),
    (10, 'RANGE', 'CPT_HCPCS', NULL, '46600', '46615', NULL, 'GI_ENDOSCOPY', 'COLONOSCOPY_LOWER_ENDOSCOPY', 'Colonoscopy / Lower Endoscopy', 'Anoscopy.'),

    (10, 'RANGE', 'CPT_HCPCS', NULL, '99281', '99285', NULL, 'OTHER', 'EMERGENCY_DEPARTMENT', 'Emergency Department', 'ED visit levels 1-5.'),
    (10, 'RANGE', 'CPT_HCPCS', NULL, '99291', '99292', NULL, 'OTHER', 'CRITICAL_CARE', 'Critical Care', 'Initial and additional critical care time.'),
    (10, 'RANGE', 'CPT_HCPCS', NULL, '99221', '99239', NULL, 'OTHER', 'INPATIENT_OBSERVATION_EM', 'Inpatient/Observation E&M', 'Professional E/M codes; expect these to be excluded downstream as professional services.'),

    /* ---- Priority 20: KEYWORD rules on the official description. The
           radiology section mixes modalities within numeric ranges, so these
           act as a safety net for anything the range rules did not catch.
           Biopsy is included because it is too scattered numerically for a
           clean range rule. ---- */
    (20, 'KEYWORD', 'CPT_HCPCS', NULL, NULL, NULL, 'MRI|MAGNETIC RESONANCE|(^|[^A-Z])MRA([^A-Z]|$)',
        'IMAGING', 'MRI_MRA', 'MRI/MRA', 'Fallback description match for MRI/MRA.'),
    (20, 'KEYWORD', 'CPT_HCPCS', NULL, NULL, NULL, 'CT SCAN|COMPUTED TOMOGRAPH|(^|[^A-Z])CTA([^A-Z]|$)',
        'IMAGING', 'CT_CTA', 'CT/CTA', 'Fallback description match for CT/CTA.'),
    (20, 'KEYWORD', 'CPT_HCPCS', NULL, NULL, NULL, 'X-RAY|XRAY|RADIOGRAPH|FLUOROSCOP',
        'IMAGING', 'XRAY_FLUOROSCOPY', 'X-Ray/Fluoroscopy', 'Fallback description match for X-ray/fluoroscopy.'),
    (20, 'KEYWORD', 'CPT_HCPCS', NULL, NULL, NULL, 'ULTRASOUND|SONOGRAPH',
        'IMAGING', 'DIAGNOSTIC_ULTRASOUND', 'Diagnostic Ultrasound', 'Fallback description match for ultrasound.'),
    (20, 'KEYWORD', 'CPT_HCPCS', NULL, NULL, NULL, 'MAMMOGRAPH|TOMOSYNTH',
        'IMAGING', 'MAMMOGRAPHY', 'Mammography', 'Fallback description match for mammography.'),
    (20, 'KEYWORD', 'CPT_HCPCS', NULL, NULL, NULL, 'BIOPSY|(^|[^A-Z])BX([^A-Z]|$)|FINE NEEDLE ASPIRATION|ASPIRAT|(^|[^A-Z])FNA([^A-Z]|$)',
        'BIOPSY', 'BIOPSY', 'Biopsy', 'Biopsy codes are too numerically scattered across CPT sections for a range rule. Matches the CMS abbreviation "BX" as well, since official short descriptors rarely spell out "BIOPSY."'),
    (20, 'KEYWORD', 'CPT_HCPCS', NULL, NULL, NULL, 'ENDOSCOP',
        'GI_ENDOSCOPY', 'ENDOSCOPY_OTHER', 'Other Endoscopy', 'Catches endoscopy codes outside the clean upper/lower GI ranges above.'),

    /* ---- Priority 30: coarse fallback by CPT section, so nothing is dropped
           without a classification. Anything landing here with real hospital
           support is a candidate for its own explicit rule. ---- */
    (30, 'RANGE', 'CPT_HCPCS', NULL, '10000', '69999', NULL, 'OTHER', 'PROCEDURE_SURGERY', 'Procedure/Surgery', 'CPT surgery section fallback.'),
    (30, 'RANGE', 'CPT_HCPCS', NULL, '80000', '89999', NULL, 'OTHER', 'LABORATORY_PATHOLOGY', 'Laboratory/Pathology', 'CPT lab/pathology section fallback.'),
    (30, 'RANGE', 'CPT_HCPCS', NULL, '90000', '99999', NULL, 'OTHER', 'EVALUATION_MANAGEMENT', 'Evaluation & Management / Medicine', 'CPT E/M and medicine section fallback.'),
    (30, 'RANGE', 'CPT_HCPCS', NULL, '70000', '79999', NULL, 'IMAGING', 'RADIOLOGY_OTHER', 'Other Radiology', 'Radiology-section code not matched by any imaging rule above.'),

    /* ---- MS-DRG: superfamily is always inpatient. Family is refined to MDC
           granularity in Section 4 rather than one family per DRG. ---- */
    (10, 'RANGE', 'MS_DRG', NULL, '001', '999', NULL, 'INPATIENT', 'INPATIENT_DRG', 'Inpatient DRG', 'All MS-DRGs; family refined to MDC granularity in Section 4, not per-DRG.')
;


/* =============================================================================
   1E. DRG acuity keywords.

   Produces an objective flag (DRG_ACUITY_KEYWORD_TAG), not a shoppability
   label. Kept in one table so the keyword list is auditable in a single place.
   ============================================================================= */

CREATE OR REPLACE TABLE O_ECON.COMMON.HPT_REF_DRG_ACUITY_KEYWORDS (
    KEYWORD_PATTERN     VARCHAR,
    ACUITY_TAG          VARCHAR   -- ACUTE_EMERGENT | ELECTIVE_SURGICAL | MATERNITY_NEWBORN
);

INSERT INTO O_ECON.COMMON.HPT_REF_DRG_ACUITY_KEYWORDS VALUES
    ('TRAUMA|FRACTURE|POISON|TOXIC|BURN', 'ACUTE_EMERGENT'),
    ('STROKE|CEREBRAL INFARCTION|INTRACRANIAL HEMORRHAGE', 'ACUTE_EMERGENT'),
    ('ACUTE MYOCARDIAL INFARCTION|CARDIAC ARREST|SHOCK|SEPSIS|SEPTICEMIA', 'ACUTE_EMERGENT'),
    ('RESPIRATORY FAILURE|COMA|TRACHEOSTOMY|ECMO|VENTILATOR', 'ACUTE_EMERGENT'),
    ('APPENDECTOMY|APPENDIX|GASTROINTESTINAL HEMORRHAGE|OBSTRUCTION', 'ACUTE_EMERGENT'),
    ('JOINT REPLACEMENT|HIP REPLACEMENT|KNEE REPLACEMENT|SPINAL FUSION|HERNIA', 'ELECTIVE_SURGICAL'),
    ('CHOLECYSTECTOMY|MASTECTOMY|PROSTATECTOMY|HYSTERECTOMY|BARIATRIC|CATARACT', 'ELECTIVE_SURGICAL'),
    ('PACEMAKER|DEFIBRILLATOR|VALVE PROCEDURE|BYPASS', 'ELECTIVE_SURGICAL'),
    ('DELIVERY|CESAREAN|ANTEPARTUM|POSTPARTUM|NEWBORN|NEONAT', 'MATERNITY_NEWBORN')
;


/* =============================================================================
   1F. Manual override table -- the only place per-code manual judgment enters.

   Empty by default. A row belongs here only when the automatic rules get
   something wrong for a code that matters to the analysis. With the rule table
   doing the classification work, this should hold a few dozen rows at most.
   Every row carries a stated reason, so the exceptions are documented rather
   than implicit.
   ============================================================================= */

CREATE TABLE IF NOT EXISTS O_ECON.COMMON.HPT_REF_MANUAL_OVERRIDES (
    BILLING_CODE_TYPE       VARCHAR,
    BILLING_CODE            VARCHAR,
    OVERRIDE_FIELD          VARCHAR,   -- e.g. 'ANALYSIS_FAMILY_ID', 'INCLUDE_FLAG'
    OVERRIDE_VALUE          VARCHAR,
    OVERRIDE_REASON         VARCHAR,
    ADDED_AT                TIMESTAMP_NTZ DEFAULT CURRENT_TIMESTAMP()
);


/* =============================================================================
   SECTION 2 -- STAGE 1: UNIVERSE AND SUPPORT

   Every CPT/HCPCS and MS-DRG code with at least minimal support in the study
   window, from a single scan of the raw source. No curation happens here, by
   design: the universe is defined by what the data contains, not by what was
   expected to be in it.
   ============================================================================= */

CREATE OR REPLACE TRANSIENT TABLE O_ECON.COMMON.HPT_P1_UNIVERSE_SUPPORT AS

WITH RAW_NORMALIZED AS (
    SELECT
        TO_VARCHAR(PROVIDER_ID)                              AS HOSPITAL_ID,
        UPPER(TRIM(TO_VARCHAR(PROVIDER_STATE_REF)))          AS PROVIDER_STATE,
        UPPER(TRIM(TO_VARCHAR(COUNTY)))                      AS COUNTY,
        TRY_TO_TIMESTAMP_NTZ(TO_VARCHAR(INGESTED_ON))        AS INGESTED_ON,

        CASE
            WHEN UPPER(TRIM(TO_VARCHAR(BILLING_CODE_TYPE)))
                 IN ('CPT', 'HCPCS', 'CPT/HCPCS', 'CPT_HCPCS')
                THEN 'CPT_HCPCS'
            WHEN UPPER(TRIM(TO_VARCHAR(BILLING_CODE_TYPE)))
                 IN ('MS-DRG', 'MS_DRG', 'MS DRG', 'DRG')
                THEN 'MS_DRG'
            ELSE NULL
        END AS BILLING_CODE_TYPE,

        CASE
            WHEN UPPER(TRIM(TO_VARCHAR(BILLING_CODE_TYPE))) IN ('MS-DRG','MS_DRG','MS DRG','DRG')
             AND REGEXP_LIKE(REGEXP_REPLACE(UPPER(TRIM(TO_VARCHAR(BILLING_CODE))), '[^A-Z0-9]', ''), '^[0-9]{1,3}$')
                THEN LPAD(REGEXP_REPLACE(UPPER(TRIM(TO_VARCHAR(BILLING_CODE))), '[^A-Z0-9]', ''), 3, '0')
            WHEN UPPER(TRIM(TO_VARCHAR(BILLING_CODE_TYPE))) IN ('CPT','HCPCS','CPT/HCPCS','CPT_HCPCS')
             AND REGEXP_LIKE(REGEXP_REPLACE(UPPER(TRIM(TO_VARCHAR(BILLING_CODE))), '[^A-Z0-9]', ''), '^[0-9]{1,5}$')
                THEN LPAD(REGEXP_REPLACE(UPPER(TRIM(TO_VARCHAR(BILLING_CODE))), '[^A-Z0-9]', ''), 5, '0')
            ELSE REGEXP_REPLACE(UPPER(TRIM(TO_VARCHAR(BILLING_CODE))), '[^A-Z0-9]', '')
        END AS BILLING_CODE,

        NULLIF(TRIM(TO_VARCHAR(DESCRIPTION)), '')            AS RAW_DESCRIPTION,
        TRY_TO_DOUBLE(TO_VARCHAR(NEGOTIATED_DOLLAR))         AS NEGOTIATED_DOLLAR

    FROM O_ECON.COMMON.HR_SUBSET_ALL
    CROSS JOIN O_ECON.COMMON.HPT_P1_STUDY_WINDOW W

    WHERE PROVIDER_ID IS NOT NULL
      AND TRY_TO_TIMESTAMP_NTZ(TO_VARCHAR(INGESTED_ON)) BETWEEN W.WINDOW_START_DATE AND W.WINDOW_END_DATE
      AND UPPER(TRIM(TO_VARCHAR(BILLING_CODE_TYPE)))
          IN ('CPT','HCPCS','CPT/HCPCS','CPT_HCPCS','MS-DRG','MS_DRG','MS DRG','DRG')
      AND NULLIF(REGEXP_REPLACE(UPPER(TRIM(TO_VARCHAR(BILLING_CODE))), '[^A-Z0-9]', ''), '') IS NOT NULL
      AND TRY_TO_DOUBLE(TO_VARCHAR(NEGOTIATED_DOLLAR)) > 0
)

SELECT
    BILLING_CODE_TYPE,
    BILLING_CODE,
    MIN(RAW_DESCRIPTION)                          AS EXAMPLE_RAW_DESCRIPTION,
    COUNT(*)                                      AS N_POSITIVE_DOLLAR_ROWS,
    COUNT(DISTINCT HOSPITAL_ID)                   AS N_HOSPITALS,
    COUNT(DISTINCT PROVIDER_STATE)                AS N_STATES,
    COUNT(DISTINCT PROVIDER_STATE || '|' || COUNTY) AS N_COUNTIES,
    MIN(INGESTED_ON)                              AS FIRST_OBSERVED,
    MAX(INGESTED_ON)                              AS LAST_OBSERVED
FROM RAW_NORMALIZED
WHERE BILLING_CODE_TYPE IS NOT NULL
GROUP BY BILLING_CODE_TYPE, BILLING_CODE
;

-- QA. This is the one expensive scan in the file; confirm it ran sanely
-- before continuing.
SELECT
    BILLING_CODE_TYPE,
    COUNT(*)                           AS N_CANDIDATE_CODES,
    SUM(N_HOSPITALS)                   AS SUM_HOSPITAL_TOUCHES,
    COUNT_IF(N_HOSPITALS >= 20)        AS N_CODES_WITH_20PLUS_HOSPITALS,
    COUNT_IF(N_HOSPITALS < 5)          AS N_CODES_WITH_UNDER_5_HOSPITALS
FROM O_ECON.COMMON.HPT_P1_UNIVERSE_SUPPORT
GROUP BY BILLING_CODE_TYPE
;


/* =============================================================================
   SECTION 3 -- STAGE 2: OFFICIAL VALIDATION

   Joins each candidate code to the official CMS reference files. Produces
   validity flags, the official description (which every keyword rule matches
   against, rather than hospital free text), and the objective attribute
   columns that shoppability schemes are later built from.
   ============================================================================= */

CREATE OR REPLACE TABLE O_ECON.COMMON.HPT_P1_VALIDATED AS

WITH OPPS AS (
    SELECT
        UPPER(REGEXP_REPLACE(BILLING_CODE, '[^A-Z0-9]', '')) AS BILLING_CODE,
        MIN(EFFECTIVE_DATE)                                  AS FIRST_VALID_DATE,
        MAX(EFFECTIVE_DATE)                                  AS LAST_VALID_DATE,
        MAX_BY(SHORT_DESCRIPTION, EFFECTIVE_DATE)            AS OFFICIAL_DESCRIPTION,
        LISTAGG(DISTINCT UPPER(STATUS_INDICATOR), ',')
            WITHIN GROUP (ORDER BY UPPER(STATUS_INDICATOR))  AS OPPS_STATUS_INDICATORS
    FROM O_ECON.COMMON.HPT_REF_OPPS_ADDENDUM_B
    WHERE BILLING_CODE IS NOT NULL
    GROUP BY UPPER(REGEXP_REPLACE(BILLING_CODE, '[^A-Z0-9]', ''))
),

MSDRG AS (
    SELECT
        LPAD(REGEXP_REPLACE(DRG_CODE, '[^0-9]', ''), 3, '0')  AS BILLING_CODE,
        MIN(EFFECTIVE_START_DATE)                             AS FIRST_VALID_DATE,
        MAX(EFFECTIVE_END_DATE)                                AS LAST_VALID_DATE,
        MAX_BY(DRG_DESCRIPTION, EFFECTIVE_START_DATE)         AS OFFICIAL_DESCRIPTION,
        MAX_BY(MDC, EFFECTIVE_START_DATE)                     AS MDC,
        MAX_BY(MEDICAL_SURGICAL_TYPE, EFFECTIVE_START_DATE)   AS MEDICAL_SURGICAL_TYPE
    FROM O_ECON.COMMON.HPT_REF_MSDRG
    WHERE DRG_CODE IS NOT NULL
    GROUP BY LPAD(REGEXP_REPLACE(DRG_CODE, '[^0-9]', ''), 3, '0')
)

SELECT
    U.*,

    COALESCE(O.OFFICIAL_DESCRIPTION, M.OFFICIAL_DESCRIPTION) AS OFFICIAL_DESCRIPTION,

    IFF(
        (U.BILLING_CODE_TYPE = 'CPT_HCPCS' AND O.BILLING_CODE IS NOT NULL)
        OR (U.BILLING_CODE_TYPE = 'MS_DRG' AND M.BILLING_CODE IS NOT NULL),
        1, 0
    ) AS IS_OFFICIALLY_VALID_FLAG,

    IFF(
        (U.BILLING_CODE_TYPE = 'CPT_HCPCS'
             AND O.FIRST_VALID_DATE IS NOT NULL
             AND U.FIRST_OBSERVED <= O.LAST_VALID_DATE
             AND U.LAST_OBSERVED  >= O.FIRST_VALID_DATE)
        OR
        (U.BILLING_CODE_TYPE = 'MS_DRG'
             AND M.FIRST_VALID_DATE IS NOT NULL
             AND U.FIRST_OBSERVED <= M.LAST_VALID_DATE
             AND U.LAST_OBSERVED  >= M.FIRST_VALID_DATE),
        1, 0
    ) AS VALID_DURING_STUDY_WINDOW_FLAG,

    O.OPPS_STATUS_INDICATORS,
    M.MDC,
    M.MEDICAL_SURGICAL_TYPE,

    IFF(N.ADDON_CODE IS NOT NULL, 1, 0)      AS IS_NCCI_ADDON_FLAG,
    N.ADDON_TYPES,
    N.PRIMARY_CODE_LISTS                     AS NCCI_PRIMARY_CODE_LISTS,

    IFF(A.BILLING_CODE IS NOT NULL, 1, 0)    AS IS_ASC_COVERED_FLAG,
    A.PAYMENT_INDICATORS                     AS ASC_PAYMENT_INDICATORS

FROM O_ECON.COMMON.HPT_P1_UNIVERSE_SUPPORT U

LEFT JOIN OPPS O
    ON U.BILLING_CODE_TYPE = 'CPT_HCPCS' AND U.BILLING_CODE = O.BILLING_CODE

LEFT JOIN MSDRG M
    ON U.BILLING_CODE_TYPE = 'MS_DRG' AND U.BILLING_CODE = M.BILLING_CODE

LEFT JOIN O_ECON.COMMON.HPT_REF_NCCI_ADDON_CURRENT N
    ON U.BILLING_CODE_TYPE = 'CPT_HCPCS' AND U.BILLING_CODE = N.ADDON_CODE

LEFT JOIN O_ECON.COMMON.HPT_REF_ASC_COVERED_CURRENT A
    ON U.BILLING_CODE_TYPE = 'CPT_HCPCS' AND U.BILLING_CODE = A.BILLING_CODE
;

-- QA
SELECT
    BILLING_CODE_TYPE,
    COUNT(*)                                    AS N_CODES,
    COUNT_IF(IS_OFFICIALLY_VALID_FLAG = 1)      AS N_VALID,
    COUNT_IF(VALID_DURING_STUDY_WINDOW_FLAG = 1) AS N_VALID_IN_WINDOW,
    COUNT_IF(IS_NCCI_ADDON_FLAG = 1)            AS N_ADDON_CODES,
    COUNT_IF(IS_ASC_COVERED_FLAG = 1)           AS N_ASC_COVERED
FROM O_ECON.COMMON.HPT_P1_VALIDATED
GROUP BY BILLING_CODE_TYPE
;



/* =============================================================================
   SECTION 4 -- STAGE 3: FAMILY AND CONCEPT ASSIGNMENT

   Applies the rule engine, then builds the concept grouping. The concept is
   the analysis unit: it sits between the individual billing code and the broad
   modality, and collapses administrative variants of the same clinical service
   into a single object.
   ============================================================================= */

CREATE OR REPLACE TABLE O_ECON.COMMON.HPT_P1_CLASSIFIED AS

WITH NUMERIC_CODES AS (
    SELECT
        V.*,
        TRY_TO_NUMBER(V.BILLING_CODE) AS NUMERIC_CODE
    FROM O_ECON.COMMON.HPT_P1_VALIDATED V
    WHERE V.IS_OFFICIALLY_VALID_FLAG = 1
      AND V.VALID_DURING_STUDY_WINDOW_FLAG = 1
),

RULE_MATCHES AS (
    SELECT
        C.*,
        R.PRIORITY,
        R.ANALYSIS_SUPERFAMILY_ID,
        R.ANALYSIS_FAMILY_ID,
        R.ANALYSIS_FAMILY_NAME,
        R.RULE_NOTE,
        R.RULE_TYPE
    FROM NUMERIC_CODES C

    INNER JOIN O_ECON.COMMON.HPT_REF_FAMILY_RULES R
        ON C.BILLING_CODE_TYPE = R.BILLING_CODE_TYPE
       AND (
            (R.RULE_TYPE = 'EXACT'   AND C.BILLING_CODE = R.MATCH_CODE)
         OR (R.RULE_TYPE = 'RANGE'   AND C.NUMERIC_CODE IS NOT NULL
                                     AND C.NUMERIC_CODE BETWEEN TRY_TO_NUMBER(R.RANGE_START) AND TRY_TO_NUMBER(R.RANGE_END))
         OR (R.RULE_TYPE = 'KEYWORD' AND REGEXP_LIKE(UPPER(COALESCE(C.OFFICIAL_DESCRIPTION, '')), '.*(' || R.KEYWORD_PATTERN || ').*'))
         OR (R.RULE_TYPE = 'CODE_PATTERN' AND REGEXP_LIKE(C.BILLING_CODE, R.KEYWORD_PATTERN))
           )

    QUALIFY ROW_NUMBER() OVER (
        PARTITION BY C.BILLING_CODE_TYPE, C.BILLING_CODE
        ORDER BY R.PRIORITY, R.RULE_TYPE
    ) = 1
),

STEMMED AS (
    SELECT
        RM.*,

        /* Contrast and variant language is stripped so that "CT chest without
           contrast" and "CT chest with contrast" reduce to the same stem. The
           stem is what groups codes into a shared concept below, so this is
           load-bearing rather than descriptive.

           CMS official short descriptors do not consistently spell out
           "contrast" -- many MRI and CT descriptors use "dye" instead ("Mri
           brain w/dye", "Mri brain w/o dye", "Mri brain w/o & w/dye"). A
           pattern matching only "contrast" leaves every dye-phrased variant as
           its own concept instead of collapsing them. Both phrasings are
           matched here, and the CONTRAST_VARIANT column below records which
           form each code used. */
        TRIM(
            REGEXP_REPLACE(
                REGEXP_REPLACE(
                    UPPER(COALESCE(RM.OFFICIAL_DESCRIPTION, '')),
                    'WITHOUT AND WITH CONTRAST|W/O (AND|&) W/ CONTRAST|W/O[[:space:]]*(AND|&)[[:space:]]*W/[[:space:]]*DYE'
                    || '|WITHOUT CONTRAST|W/O CONTRAST|W/O[[:space:]]*DYE'
                    || '|WITH CONTRAST|W/ CONTRAST|W/[[:space:]]*DYE',
                    ''
                ),
                '[[:space:]]+', ' '
            )
        ) AS CONCEPT_STEM_TEXT,

        CASE
            WHEN REGEXP_LIKE(
                    UPPER(COALESCE(RM.OFFICIAL_DESCRIPTION,'')),
                    '.*(WITHOUT AND WITH CONTRAST|W/O (AND|&) W/ CONTRAST|W/O[[:space:]]*(AND|&)[[:space:]]*W/[[:space:]]*DYE).*'
                 )
                THEN 'WITHOUT_AND_WITH_CONTRAST'
            WHEN REGEXP_LIKE(
                    UPPER(COALESCE(RM.OFFICIAL_DESCRIPTION,'')),
                    '.*(WITHOUT CONTRAST|W/O CONTRAST|W/O[[:space:]]*DYE).*'
                 )
                THEN 'WITHOUT_CONTRAST'
            WHEN REGEXP_LIKE(
                    UPPER(COALESCE(RM.OFFICIAL_DESCRIPTION,'')),
                    '.*(WITH CONTRAST|W/ CONTRAST|W/[[:space:]]*DYE).*'
                 )
                THEN 'WITH_CONTRAST'
            ELSE NULL
        END AS CONTRAST_VARIANT,

        /* DRG family is refined from the flat INPATIENT_DRG bucket to MDC
           granularity. One family per DRG would produce dozens of
           near-singleton families, which is useless for any family-level
           analysis or for family-clustered inference. */
        CASE
            WHEN RM.BILLING_CODE_TYPE = 'MS_DRG'
                THEN 'MDC_' || COALESCE(RM.MDC, 'UNKNOWN')
            ELSE RM.ANALYSIS_FAMILY_ID
        END AS ANALYSIS_FAMILY_ID_REFINED,

        /* Objective DRG acuity flag, reported as a flag and never folded into
           a shoppability decision. */
        K.ACUITY_TAG AS DRG_ACUITY_KEYWORD_TAG

    FROM RULE_MATCHES RM

    LEFT JOIN O_ECON.COMMON.HPT_REF_DRG_ACUITY_KEYWORDS K
        ON RM.BILLING_CODE_TYPE = 'MS_DRG'
       AND REGEXP_LIKE(UPPER(COALESCE(RM.OFFICIAL_DESCRIPTION, '')), '.*(' || K.KEYWORD_PATTERN || ').*')

    QUALIFY ROW_NUMBER() OVER (
        PARTITION BY RM.BILLING_CODE_TYPE, RM.BILLING_CODE
        ORDER BY IFF(K.ACUITY_TAG IS NULL, 1, 0)
    ) = 1
)

SELECT
    S.*,

    /* THE CONCEPT GROUPING KEY. Two exact codes in the same family with the
       same stem text -- MRI brain with and without contrast, say -- collapse
       into one concept, distinguished only by CONTRAST_VARIANT.

       Where the stem is empty or generic, the fallback is one concept per
       exact code. That covers MS-DRG, where each DRG is already its own
       specific concept, and any code whose description did not come through
       OPPS or MS-DRG validation. Falling back rather than grouping means
       nothing is silently merged on missing data. */
    CASE
        WHEN S.BILLING_CODE_TYPE = 'MS_DRG'
            THEN S.ANALYSIS_FAMILY_ID_REFINED || '_DRG_' || S.BILLING_CODE
        WHEN NULLIF(S.CONCEPT_STEM_TEXT, '') IS NULL
            THEN S.ANALYSIS_FAMILY_ID_REFINED || '_' || S.BILLING_CODE
        ELSE
            S.ANALYSIS_FAMILY_ID_REFINED || '_'
            || REGEXP_REPLACE(S.CONCEPT_STEM_TEXT, '[^A-Z0-9]+', '_')
    END AS ANALYSIS_CONCEPT_ID

FROM STEMMED S
;

-- QA: how many codes fell through to the coarse fallback families? A large
-- count with real hospital support means a specific rule should be added to
-- HPT_REF_FAMILY_RULES, not that codes need individual review.
SELECT
    ANALYSIS_FAMILY_ID,
    COUNT(*)            AS N_CODES,
    SUM(N_HOSPITALS)    AS SUM_HOSPITAL_TOUCHES
FROM O_ECON.COMMON.HPT_P1_CLASSIFIED
WHERE ANALYSIS_FAMILY_ID IN ('PROCEDURE_SURGERY','RADIOLOGY_OTHER','LABORATORY_PATHOLOGY','EVALUATION_MANAGEMENT','ENDOSCOPY_OTHER')
GROUP BY ANALYSIS_FAMILY_ID
ORDER BY SUM_HOSPITAL_TOUCHES DESC
;


/* =============================================================================
   SECTION 5 -- APPLY MANUAL OVERRIDES AND FREEZE THE CODEBOOK

   Manual judgment re-enters here and nowhere else, and only for codes
   explicitly listed in HPT_REF_MANUAL_OVERRIDES.
   ============================================================================= */

CREATE OR REPLACE TABLE O_ECON.COMMON.HPT_P1_FINAL_CODEBOOK AS

SELECT
    C.BILLING_CODE_TYPE,
    C.BILLING_CODE,
    C.OFFICIAL_DESCRIPTION,

    COALESCE(
        MAX(IFF(O_FAM.OVERRIDE_FIELD = 'ANALYSIS_FAMILY_ID', O_FAM.OVERRIDE_VALUE, NULL)),
        C.ANALYSIS_FAMILY_ID_REFINED
    ) AS ANALYSIS_FAMILY_ID,

    C.ANALYSIS_SUPERFAMILY_ID,
    C.ANALYSIS_FAMILY_NAME,
    C.ANALYSIS_CONCEPT_ID,
    C.CONCEPT_STEM_TEXT,
    C.CONTRAST_VARIANT,

    C.N_HOSPITALS,
    C.N_STATES,
    C.N_COUNTIES,
    C.N_POSITIVE_DOLLAR_ROWS,
    C.FIRST_OBSERVED,
    C.LAST_OBSERVED,

    /* Objective attributes only. Shoppability schemes are built from these
       columns in R. */
    C.OPPS_STATUS_INDICATORS,
    C.IS_NCCI_ADDON_FLAG,
    C.NCCI_PRIMARY_CODE_LISTS,
    C.IS_ASC_COVERED_FLAG,
    C.ASC_PAYMENT_INDICATORS,
    C.MDC,
    C.MEDICAL_SURGICAL_TYPE,
    C.DRG_ACUITY_KEYWORD_TAG,

    IFF(CMS.BILLING_CODE IS NOT NULL, 1, 0) AS IS_CMS70_CODE,
    CMS.CMS70_SERVICE_ID,
    CMS.CMS70_SERVICE_NAME,

    MAX(IFF(O_INC.OVERRIDE_FIELD = 'INCLUDE_FLAG', O_INC.OVERRIDE_VALUE, NULL)) AS MANUAL_INCLUDE_OVERRIDE,
    LISTAGG(DISTINCT O_ANY.OVERRIDE_REASON, ' | ') AS OVERRIDE_NOTES,

    /* ---------------------------------------------------------------------
       SCOPE FLAG -- which codes feed stage 1b.

       The analysis sample is: CT, MRI, X-ray, ultrasound, mammography,
       bone density, biopsy, and endoscopy families; as many CMS-70 services
       as possible; plus inpatient MS-DRG and emergency/critical-care codes,
       which form the non-shoppable comparison group.

       THE SCOPE TEST IS ON ANALYSIS_FAMILY_ID, NOT ANALYSIS_SUPERFAMILY_ID.
       Testing the superfamily would pull in RADIOLOGY_OTHER and
       ENDOSCOPY_OTHER -- the fallback buckets holding PET scans, radiation
       therapy planning, and nuclear medicine studies -- purely because they
       share the IMAGING and GI_ENDOSCOPY superfamily labels with the target
       families. Those are not the services the design is about, so they are
       excluded by being left off the explicit list. Adding RADIOLOGY_OTHER
       here would bring them back in as a broad "other imaging" catch-all.

       SCOPE IS "TARGET FAMILY OR CMS-70", NOT FAMILY MEMBERSHIP ALONE. The
       official CMS-70 list includes lab panels and E/M office visits
       alongside imaging, and those land in the generic LABORATORY_PATHOLOGY
       and EVALUATION_MANAGEMENT fallback families under the rule table. A
       family-only filter would silently drop CMS-70 codes the design
       requires.

       INPATIENT_OBSERVATION_EM (99221-99239) is deliberately excluded despite
       sitting in the OTHER superfamily: those are professional physician E/M
       codes rather than hospital facility prices.
       --------------------------------------------------------------------- */
    CASE
        WHEN C.ANALYSIS_FAMILY_ID IN (
                'CT_CTA', 'MRI_MRA', 'XRAY_FLUOROSCOPY',
                'DIAGNOSTIC_ULTRASOUND', 'VASCULAR_ULTRASOUND', 'ECHOCARDIOGRAPHY',
                'MAMMOGRAPHY', 'BONE_DENSITY', 'BIOPSY',
                'UPPER_ENDOSCOPY', 'COLONOSCOPY_LOWER_ENDOSCOPY'
             )
            THEN 1
        WHEN C.ANALYSIS_FAMILY_ID IN ('EMERGENCY_DEPARTMENT', 'CRITICAL_CARE')
            THEN 1
        WHEN C.BILLING_CODE_TYPE = 'MS_DRG'
            THEN 1
        WHEN CMS.BILLING_CODE IS NOT NULL
            THEN 1
        ELSE 0
    END AS INCLUDE_IN_SCOPE,

    CASE
        WHEN C.ANALYSIS_FAMILY_ID IN (
                'CT_CTA', 'MRI_MRA', 'XRAY_FLUOROSCOPY',
                'DIAGNOSTIC_ULTRASOUND', 'VASCULAR_ULTRASOUND', 'ECHOCARDIOGRAPHY',
                'MAMMOGRAPHY', 'BONE_DENSITY', 'BIOPSY',
                'UPPER_ENDOSCOPY', 'COLONOSCOPY_LOWER_ENDOSCOPY'
             )
            THEN 'TARGET_IMAGING_ENDOSCOPY_BIOPSY_FAMILY'
        WHEN C.ANALYSIS_FAMILY_ID IN ('EMERGENCY_DEPARTMENT', 'CRITICAL_CARE')
            THEN 'TARGET_EMERGENCY_FAMILY'
        WHEN C.BILLING_CODE_TYPE = 'MS_DRG'
            THEN 'TARGET_INPATIENT_DRG'
        WHEN CMS.BILLING_CODE IS NOT NULL
            THEN 'CMS70_REQUIRED_OUTSIDE_TARGET_FAMILY'
        ELSE 'OUT_OF_STATED_SCOPE'
    END AS SCOPE_REASON

FROM O_ECON.COMMON.HPT_P1_CLASSIFIED C

LEFT JOIN O_ECON.COMMON.HPT_CMS70_CODEBOOK CMS
    ON C.BILLING_CODE_TYPE = CMS.BILLING_CODE_TYPE
   AND C.BILLING_CODE = CMS.BILLING_CODE
   AND CMS.MAPPING_STATUS = 'ACTIVE'

LEFT JOIN O_ECON.COMMON.HPT_REF_MANUAL_OVERRIDES O_FAM
    ON C.BILLING_CODE_TYPE = O_FAM.BILLING_CODE_TYPE AND C.BILLING_CODE = O_FAM.BILLING_CODE
LEFT JOIN O_ECON.COMMON.HPT_REF_MANUAL_OVERRIDES O_INC
    ON C.BILLING_CODE_TYPE = O_INC.BILLING_CODE_TYPE AND C.BILLING_CODE = O_INC.BILLING_CODE
LEFT JOIN O_ECON.COMMON.HPT_REF_MANUAL_OVERRIDES O_ANY
    ON C.BILLING_CODE_TYPE = O_ANY.BILLING_CODE_TYPE AND C.BILLING_CODE = O_ANY.BILLING_CODE

GROUP BY
    C.BILLING_CODE_TYPE, C.BILLING_CODE, C.OFFICIAL_DESCRIPTION,
    C.ANALYSIS_FAMILY_ID, C.ANALYSIS_FAMILY_ID_REFINED, C.ANALYSIS_SUPERFAMILY_ID, C.ANALYSIS_FAMILY_NAME,
    C.ANALYSIS_CONCEPT_ID, C.CONCEPT_STEM_TEXT, C.CONTRAST_VARIANT,
    C.N_HOSPITALS, C.N_STATES, C.N_COUNTIES, C.N_POSITIVE_DOLLAR_ROWS,
    C.FIRST_OBSERVED, C.LAST_OBSERVED, C.OPPS_STATUS_INDICATORS,
    C.IS_NCCI_ADDON_FLAG, C.NCCI_PRIMARY_CODE_LISTS, C.IS_ASC_COVERED_FLAG,
    C.ASC_PAYMENT_INDICATORS, C.MDC, C.MEDICAL_SURGICAL_TYPE,
    C.DRG_ACUITY_KEYWORD_TAG, CMS.BILLING_CODE, CMS.CMS70_SERVICE_ID, CMS.CMS70_SERVICE_NAME
;


/* =============================================================================
   SECTION 5B -- SCOPED CODEBOOK: WHAT FEEDS STAGE 1B

   HPT_P1_FINAL_CODEBOOK is retained as the full archive of every officially
   valid code, so scope can be widened later without rerunning stages 1 to 3.
   HPT_P1_FINAL_CODEBOOK_SCOPED is the narrower table that actually feeds the
   price pipeline.
   ============================================================================= */

CREATE OR REPLACE TABLE O_ECON.COMMON.HPT_P1_FINAL_CODEBOOK_SCOPED AS
SELECT *
FROM O_ECON.COMMON.HPT_P1_FINAL_CODEBOOK
WHERE INCLUDE_IN_SCOPE = 1
  AND COALESCE(MANUAL_INCLUDE_OVERRIDE, 'INCLUDE') <> 'EXCLUDE'
;

-- Scope QA: how many codes made the cut, by reason.
SELECT
    SCOPE_REASON,
    BILLING_CODE_TYPE,
    COUNT(*) AS N_CODES,
    COUNT(DISTINCT ANALYSIS_CONCEPT_ID) AS N_CONCEPTS
FROM O_ECON.COMMON.HPT_P1_FINAL_CODEBOOK_SCOPED
GROUP BY SCOPE_REASON, BILLING_CODE_TYPE
ORDER BY SCOPE_REASON, BILLING_CODE_TYPE;

-- All 70 CMS-required services should be present, even though most already
-- fall inside a target family.
SELECT COUNT(DISTINCT CMS70_SERVICE_ID) AS N_CMS70_SERVICES_IN_SCOPED_CODEBOOK
FROM O_ECON.COMMON.HPT_P1_FINAL_CODEBOOK_SCOPED
WHERE IS_CMS70_CODE = 1;


/* =============================================================================
   SECTION 5C -- CLASSIFICATION VERIFICATION

   Six queries covering different aspects of the classification. Run all six
   after any edit to HPT_REF_FAMILY_RULES: each checks something the others do
   not, so one query looking correct does not establish that the rule change
   behaved as intended.
   ============================================================================= */

-- (1) Superfamily coverage. Every target superfamily should hold codes.
SELECT
    BILLING_CODE_TYPE,
    ANALYSIS_SUPERFAMILY_ID,
    COUNT(*) AS N_CODES,
    COUNT(DISTINCT ANALYSIS_FAMILY_ID) AS N_FAMILIES
FROM O_ECON.COMMON.HPT_P1_FINAL_CODEBOOK_SCOPED
GROUP BY BILLING_CODE_TYPE, ANALYSIS_SUPERFAMILY_ID
ORDER BY BILLING_CODE_TYPE, ANALYSIS_SUPERFAMILY_ID;

-- (2) Scope breakdown by reason.
SELECT
    SCOPE_REASON,
    BILLING_CODE_TYPE,
    COUNT(*) AS N_CODES,
    COUNT(DISTINCT ANALYSIS_CONCEPT_ID) AS N_CONCEPTS
FROM O_ECON.COMMON.HPT_P1_FINAL_CODEBOOK_SCOPED
GROUP BY SCOPE_REASON, BILLING_CODE_TYPE
ORDER BY SCOPE_REASON, BILLING_CODE_TYPE;

-- (3) CMS-70 completeness. Should read exactly 70. A rule change that moves
-- this number has affected more than it was meant to.
SELECT COUNT(DISTINCT CMS70_SERVICE_ID) AS N_CMS70_SERVICES_IN_SCOPED_CODEBOOK
FROM O_ECON.COMMON.HPT_P1_FINAL_CODEBOOK_SCOPED
WHERE IS_CMS70_CODE = 1;

-- (4) Fallback bucket sizes on the FULL archive rather than the scoped table.
-- Codes claimed by a new specific rule should leave these buckets.
SELECT
    ANALYSIS_FAMILY_ID,
    COUNT(*) AS N_CODES,
    SUM(N_HOSPITALS) AS SUM_HOSPITAL_TOUCHES
FROM O_ECON.COMMON.HPT_P1_CLASSIFIED
WHERE ANALYSIS_FAMILY_ID IN ('PROCEDURE_SURGERY','RADIOLOGY_OTHER','LABORATORY_PATHOLOGY','EVALUATION_MANAGEMENT','ENDOSCOPY_OTHER')
GROUP BY ANALYSIS_FAMILY_ID
ORDER BY SUM_HOSPITAL_TOUCHES DESC;

-- (5) Concept rollup within the target families. Confirms codes sharing a stem
-- land in one concept rather than one concept per code.
SELECT
    ANALYSIS_FAMILY_ID,
    ANALYSIS_CONCEPT_ID,
    COUNT(*) AS N_CODES_IN_CONCEPT,
    LISTAGG(DISTINCT CONTRAST_VARIANT, ', ') WITHIN GROUP (ORDER BY CONTRAST_VARIANT) AS VARIANTS_PRESENT,
    LISTAGG(BILLING_CODE, ',') WITHIN GROUP (ORDER BY BILLING_CODE) AS BILLING_CODES
FROM O_ECON.COMMON.HPT_P1_FINAL_CODEBOOK_SCOPED
WHERE ANALYSIS_SUPERFAMILY_ID IN ('IMAGING','GI_ENDOSCOPY','BIOPSY')
GROUP BY ANALYSIS_FAMILY_ID, ANALYSIS_CONCEPT_ID
HAVING COUNT(*) > 1
ORDER BY N_CODES_IN_CONCEPT DESC
LIMIT 40;

-- (6) Biopsy codes in scope with their support, so the family can be checked
-- for real hospital coverage rather than a handful of rare codes.
SELECT
    BILLING_CODE,
    OFFICIAL_DESCRIPTION,
    ANALYSIS_CONCEPT_ID,
    N_HOSPITALS,
    N_STATES
FROM O_ECON.COMMON.HPT_P1_FINAL_CODEBOOK_SCOPED
WHERE ANALYSIS_SUPERFAMILY_ID = 'BIOPSY'
ORDER BY N_HOSPITALS DESC;


-- Exclusion rules: confirm the priority-0 patterns are claiming codes.
SELECT ANALYSIS_FAMILY_ID, COUNT(*) AS N_CODES
FROM O_ECON.COMMON.HPT_P1_CLASSIFIED
WHERE ANALYSIS_FAMILY_ID IN ('ANESTHESIA','EXCLUDED_QUALITY_MEASURE','EXCLUDED_LIKELY_MISLABELED')
GROUP BY ANALYSIS_FAMILY_ID;

-- DRG acuity tagging coverage.
SELECT DRG_ACUITY_KEYWORD_TAG, COUNT(*) AS N_DRGS
FROM O_ECON.COMMON.HPT_P1_FINAL_CODEBOOK_SCOPED
WHERE BILLING_CODE_TYPE = 'MS_DRG'
GROUP BY DRG_ACUITY_KEYWORD_TAG;


/* =============================================================================
   SECTION 6 -- FINAL QA
   ============================================================================= */

CREATE OR REPLACE TABLE O_ECON.COMMON.HPT_P1_QA_SUMMARY AS
SELECT
    BILLING_CODE_TYPE,
    ANALYSIS_SUPERFAMILY_ID,
    COUNT(*)                              AS N_CODES,
    COUNT(DISTINCT ANALYSIS_FAMILY_ID)    AS N_FAMILIES,
    SUM(N_HOSPITALS)                      AS SUM_HOSPITAL_TOUCHES,
    COUNT_IF(IS_CMS70_CODE = 1)           AS N_CMS70_CODES,
    COUNT_IF(IS_ASC_COVERED_FLAG = 1)     AS N_ASC_COVERED,
    COUNT_IF(IS_NCCI_ADDON_FLAG = 1)      AS N_ADDON_CODES,
    COUNT_IF(MANUAL_INCLUDE_OVERRIDE = 'EXCLUDE') AS N_MANUALLY_EXCLUDED
FROM O_ECON.COMMON.HPT_P1_FINAL_CODEBOOK
GROUP BY BILLING_CODE_TYPE, ANALYSIS_SUPERFAMILY_ID
ORDER BY BILLING_CODE_TYPE, ANALYSIS_SUPERFAMILY_ID
;

SELECT * FROM O_ECON.COMMON.HPT_P1_QA_SUMMARY;

-- CMS-70 completeness on the full archive.
SELECT COUNT(DISTINCT CMS70_SERVICE_ID) AS N_CMS70_SERVICES_COVERED
FROM O_ECON.COMMON.HPT_P1_FINAL_CODEBOOK
WHERE IS_CMS70_CODE = 1;

-- Concept rollup: confirms variants collapse into shared concepts rather than
-- defaulting to one concept per exact code.
SELECT
    ANALYSIS_FAMILY_ID,
    ANALYSIS_CONCEPT_ID,
    COUNT(*) AS N_CODES_IN_CONCEPT,
    LISTAGG(DISTINCT CONTRAST_VARIANT, ', ') WITHIN GROUP (ORDER BY CONTRAST_VARIANT) AS VARIANTS_PRESENT,
    LISTAGG(BILLING_CODE, ',') WITHIN GROUP (ORDER BY BILLING_CODE) AS BILLING_CODES
FROM O_ECON.COMMON.HPT_P1_FINAL_CODEBOOK
WHERE BILLING_CODE_TYPE = 'CPT_HCPCS'
GROUP BY ANALYSIS_FAMILY_ID, ANALYSIS_CONCEPT_ID
HAVING COUNT(*) > 1
ORDER BY N_CODES_IN_CONCEPT DESC
LIMIT 30;

-- Family size distribution, flagging any family too granular to be useful for
-- family-level analysis or family-clustered inference.
SELECT
    ANALYSIS_FAMILY_ID,
    COUNT(*) AS N_CODES
FROM O_ECON.COMMON.HPT_P1_FINAL_CODEBOOK
GROUP BY ANALYSIS_FAMILY_ID
ORDER BY N_CODES ASC
LIMIT 30;


/* =============================================================================
   NEXT STEP

   HPT_P1_FINAL_CODEBOOK_SCOPED is the input to 02_Phase2to4_Prices.sql. Its
   columns are a superset of what the price pipeline needs -- billing code
   type, code, official description, family and concept fields -- plus the
   objective attribute columns that the R stage uses to construct shoppability
   schemes:

     OPPS_STATUS_INDICATORS, IS_ASC_COVERED_FLAG, IS_NCCI_ADDON_FLAG,
     IS_CMS70_CODE, BILLING_CODE_TYPE, MEDICAL_SURGICAL_TYPE,
     DRG_ACUITY_KEYWORD_TAG, CONTRAST_VARIANT

   Before moving on, review the fallback-family QA above. Codes sitting in an
   *_OTHER bucket with substantial hospital support should get a specific rule
   in HPT_REF_FAMILY_RULES; anything still misclassified after that gets one
   row in HPT_REF_MANUAL_OVERRIDES with a stated reason.
   ============================================================================= */
