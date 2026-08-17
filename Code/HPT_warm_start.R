###############################################################################
# WARM START
#
# Loads a completed run into a fresh session without re-estimating anything.
#
# Evaluates two regions of HPT_Analysis_Pipeline.R and nothing else:
#
#   1. Everything above `# ORDER OF OPERATIONS`  -- Sections 0-12, pure
#      definitions, seconds.
#   2. Everything from `# SESSION RESTORE AND CACHE MANAGEMENT` to the end --
#      cache_status(), restore_session(), invalidate_cache().
#
# Then calls restore_session(), which reads the .rds cache and assigns every
# completed object into the global environment under the name the run blocks
# use: outpatient, schemes_long, concept_results, main, measures, and the rest.
#
# The BUILD block and every stage block are skipped. Nothing is estimated.
#
# USAGE
#   source("~/path/to/HPT_warm_start.R")
#   warm_start()                       # uses the default path below
#   warm_start("~/other/copy.R")       # or point it somewhere else
#
# The boundaries are located by marker text, not line number, so this keeps
# working when the pipeline file is edited.
###############################################################################

HPT_PIPELINE_FILE <- file.path(
  "/Users/danielsierra/Library/CloudStorage",
  "OneDrive-FloridaStateUniversity",
  "Hospital Price Transparency Paper",
  "Code", "HPT_Analysis_Pipeline.R"
)

warm_start <- function(pipeline_file = HPT_PIPELINE_FILE,
                       restore = TRUE,
                       quiet_restore = FALSE) {

  if (!file.exists(pipeline_file)) {
    stop("Pipeline file not found:\n  ", pipeline_file,
         "\nPass the correct path: warm_start(\"/path/to/HPT_Analysis_Pipeline.R\")",
         call. = FALSE)
  }

  src <- readLines(pipeline_file, warn = FALSE)

  i_run <- grep("^#\\s*ORDER OF OPERATIONS\\s*$", src)
  i_res <- grep("^#\\s*SESSION RESTORE AND CACHE MANAGEMENT\\s*$", src)

  if (length(i_run) != 1L) {
    stop("Expected exactly one '# ORDER OF OPERATIONS' marker, found ",
         length(i_run), ". The definitions/run boundary cannot be located.",
         call. = FALSE)
  }
  if (length(i_res) != 1L) {
    stop("Expected exactly one '# SESSION RESTORE AND CACHE MANAGEMENT' marker, ",
         "found ", length(i_res), ".", call. = FALSE)
  }
  if (i_res <= i_run) {
    stop("Markers are out of order; the file is not laid out as expected.",
         call. = FALSE)
  }

  defs <- src[seq_len(i_run - 1L)]
  rest <- src[seq(i_res, length(src))]

  cat("Pipeline:", pipeline_file, "\n")
  cat("Definitions: lines 1 -", i_run - 1L,
      "| restore utilities: lines", i_res, "-", length(src), "\n")

  # Section 0 opens with rm(list = ls()), which clears the global environment
  # and therefore this function's own binding. The running call is unaffected
  # because the closure is already bound in this frame, but the name would be
  # gone afterwards, so it is reinstated at the end.
  cat("\n", strrep("-", 70), "\nEvaluating definitions\n", strrep("-", 70), "\n", sep = "")
  eval(parse(text = defs), envir = globalenv())

  cat("\n", strrep("-", 70), "\nEvaluating restore utilities\n", strrep("-", 70), "\n", sep = "")
  eval(parse(text = rest), envir = globalenv())

  # Neutralised so that accidentally executing a stage block later in the
  # session is a no-op rather than an eight-hour surprise. Set it explicitly
  # when you actually want to run something.
  assign("RUN_STAGES", integer(0), envir = globalenv())

assign("warm_start", sys.function(), envir = globalenv())
assign("HPT_PIPELINE_FILE", pipeline_file, envir = globalenv())

if (restore) {
  cat("\n", strrep("-", 70), "\nRestoring cached objects\n", strrep("-", 70), "\n", sep = "")
  get("restore_session", envir = globalenv())(quiet = quiet_restore)
} else {
  cat("\nDefinitions loaded. Call restore_session() when ready.\n")
}

cat("\nRUN_STAGES set to integer(0). No stage block will execute until you",
    "\nset it deliberately. The BUILD block is unconditional -- do not run it;",
    "\nrestore_session() has already supplied `outpatient`.\n")

invisible(TRUE)
}

cat("warm_start() defined. Call warm_start() to load a completed run.\n")


###############################################################################
# load_t12_helpers()
#
# Section 12's SES/demographic-heterogeneity block (banner at "SECTION 12 --
# COUNTY DEMOGRAPHIC HETEROGENEITY...") is written as a SEPARATE self-
# contained script, not part of the Sections 0-12 definitions block that
# warm_start() loads. It lives entirely in the executed region, past
# `# ORDER OF OPERATIONS`, and its own "# 6. RUN" subsection immediately calls
# a 30-45 minute driver across schemes/moderators/instruments.
#
# This loads only the config and `.t12_*` function definitions -- everything
# from the "SECTION 12 --" banner up to (not including) "# 6. RUN" -- so
# `.t12_build_ses_index()` and friends become available without triggering
# that 30-45 minute run. Call `.t12_build_ses_index(outpatient)` yourself
# afterward; that call alone is fast (a PCA on the county cross-section).
#
# USAGE
#   load_t12_helpers()
#   ses_panel <- .t12_build_ses_index(outpatient)
###############################################################################

load_t12_helpers <- function(pipeline_file = HPT_PIPELINE_FILE) {
  if (!exists("outpatient", envir = globalenv())) {
    stop("`outpatient` is not in memory. Run warm_start() first.", call. = FALSE)
  }
  src <- readLines(pipeline_file, warn = FALSE)

  i_start <- grep("^#\\s*SECTION 12 -- COUNTY DEMOGRAPHIC HETEROGENEITY", src)
  i_end   <- grep("^#\\s*6\\. RUN\\s*$", src)

  if (length(i_start) != 1L) {
    stop("Expected exactly one 'SECTION 12 -- COUNTY DEMOGRAPHIC HETEROGENEITY' ",
         "banner, found ", length(i_start), ". Has the file been edited?", call. = FALSE)
  }
  if (length(i_end) != 1L) {
    stop("Expected exactly one '# 6. RUN' marker, found ", length(i_end),
         ". Has the file been edited?", call. = FALSE)
  }
  if (i_end <= i_start) stop("Markers out of order.", call. = FALSE)

  # Stop short of the "# ====..." rule and "# 6. RUN" header, leaving only
  # config + function definitions (last statement is .t12_driver()'s closing
  # brace).
  block <- src[seq(i_start, i_end - 4L)]

  cat("Loading .t12_ helpers: lines", i_start, "-", i_end - 4L, "\n")
  eval(parse(text = block), envir = globalenv())

  cat("Available: .t12_build_ses_index(), .t12_fit_triple(), .t12_driver(),",
      "\n  .t12_preflight(), and the T12T_* config objects.",
      "\nNothing has been run yet. Typical next step:",
      "\n  ses_panel <- .t12_build_ses_index(outpatient)\n")
  invisible(TRUE)
}

cat("load_t12_helpers() defined. Call it (after warm_start()) to get",
    ".t12_build_ses_index() without the 30-45 minute driver run.\n")
