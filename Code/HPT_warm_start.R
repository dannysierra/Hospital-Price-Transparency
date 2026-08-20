###############################################################################
# WARM START
#
# Loads a completed run into a fresh R session without re-estimating anything.
#
# ---------------------------------------------------------------------------
# THE SHORT VERSION
# ---------------------------------------------------------------------------
# You do not need this file to warm start. Two lines in a fresh session do it:
#
#     HPT_WARM_START <- TRUE
#     source("~/path/to/Code/HPT_Analysis_Pipeline.R")
#
# PARTS 1 and 2 of the pipeline evaluate normally (definitions only, seconds),
# PART 3 calls restore_session() instead of building, and every stage in PARTS
# 4 through 6 is switched off. The whole file runs in well under a minute and
# leaves outpatient, schemes_long, concept_results, main, measures, and the
# rest in the global environment under the names the run blocks use.
#
# warm_start() below is a convenience wrapper around exactly that, for when you
# would rather not remember the switch name.
#
# ---------------------------------------------------------------------------
# WHY THIS FILE CHANGED
# ---------------------------------------------------------------------------
# The previous version located the definitions/run boundary by grepping the
# pipeline for two comment strings, "# ORDER OF OPERATIONS" and "# SESSION
# RESTORE AND CACHE MANAGEMENT", then eval()'d the text between and after them.
# That worked, but it made this file depend on the pipeline's comment layout:
# reorganizing the pipeline deleted the first marker and moved the second into
# the definitions region, at which point warm_start() would have evaluated the
# entire pipeline including every estimation stage.
#
# A switch inside the pipeline cannot fail that way. The pipeline decides what
# to skip; this file just sets the flag.
#
# load_t12_helpers() still locates its block by marker, because that block is a
# self-contained script embedded inside a stage and has no switch of its own.
# Its two markers are checked below before use, and it reports clearly if the
# file has been edited past recognition rather than evaluating the wrong lines.
###############################################################################

HPT_PIPELINE_FILE <- file.path(
  Sys.getenv("HPT_ROOT", unset = path.expand(file.path(
    "~", "Library", "CloudStorage", "OneDrive-FloridaStateUniversity",
    "Hospital Price Transparency Paper"))),
  "Code", "HPT_Analysis_Pipeline.R"
)


# ===========================================================================
# warm_start()
# ===========================================================================
#
#   warm_start()                          # default path
#   warm_start("~/other/copy.R")          # or point it somewhere else
#   warm_start(restore = FALSE)           # definitions only, no cache read
#
# Returns invisibly. Prints what it loaded and what is still missing.

warm_start <- function(pipeline_file = HPT_PIPELINE_FILE, restore = TRUE) {
  
  if (!file.exists(pipeline_file)) {
    stop("Pipeline file not found:\n  ", pipeline_file,
         "\n\nEither set HPT_ROOT:\n",
         '  Sys.setenv(HPT_ROOT = "/path/to/Hospital Price Transparency Paper")\n',
         "or pass the path directly:\n",
         '  warm_start("/path/to/HPT_Analysis_Pipeline.R")',
         call. = FALSE)
  }
  
  src <- readLines(pipeline_file, warn = FALSE)
  if (!any(grepl("^HPT_WARM_START <- ", src))) {
    stop("This pipeline file has no HPT_WARM_START switch, so it is the old\n",
         "  single-block layout. Either use the reorganized pipeline, or load\n",
         "  the old one by hand.", call. = FALSE)
  }
  
  cat("Pipeline: ", pipeline_file, "\n",
      "Sourcing PARTS 1-2 (definitions) and PART 3 (restore). No estimation.\n",
      sep = "")
  
  assign("HPT_WARM_START", TRUE, envir = globalenv())
  on.exit(assign("HPT_WARM_START", FALSE, envir = globalenv()), add = TRUE)
  
  source(pipeline_file, local = FALSE, echo = FALSE)
  
  # source() re-evaluates the pipeline's own definition of these, so reinstate
  # this file's copies afterwards.
  assign("warm_start",        sys.function(),  envir = globalenv())
  assign("load_t12_helpers",  load_t12_helpers, envir = globalenv())
  assign("HPT_PIPELINE_FILE", pipeline_file,    envir = globalenv())
  
  if (!restore) {
    cat("\nrestore = FALSE, so nothing was read from the cache.",
        "\nCall restore_session() when ready.\n")
    return(invisible(TRUE))
  }
  
  cat("\nRUN_STAGES is empty and every HPT_RUN switch is FALSE, so no stage",
      "\nblock will execute. Set them deliberately to run something:",
      "\n  RUN_STAGES <- 7 ; HPT_RUN$figures <- TRUE",
      "\n\ncache_status() shows what is on disk.\n")
  invisible(TRUE)
}


# ===========================================================================
# load_t12_helpers()
# ===========================================================================
#
# Section 12's SES and demographic-heterogeneity code is written as a
# self-contained script embedded inside a stage block, so a warm start does
# not put its .t12_* functions in scope. Sourcing the whole block would kick
# off a 30-45 minute driver run across schemes, moderators, and instruments.
#
# This evaluates only the config and function definitions, stopping short of
# the block's own "# 6. RUN" subsection.
#
#   load_t12_helpers()
#   ses_panel <- .t12_build_ses_index(outpatient)      # fast: a PCA
#
# Still marker-based, because the block has no switch of its own. Both markers
# are verified before anything is evaluated.

load_t12_helpers <- function(pipeline_file = HPT_PIPELINE_FILE) {
  if (!exists("outpatient", envir = globalenv()))
    stop("`outpatient` is not in memory. Run warm_start() first.", call. = FALSE)
  if (!file.exists(pipeline_file))
    stop("Pipeline file not found:\n  ", pipeline_file, call. = FALSE)
  
  src <- readLines(pipeline_file, warn = FALSE)
  i_start <- grep("^#\\s*SECTION 12 -- COUNTY DEMOGRAPHIC HETEROGENEITY", src)
  i_end   <- grep("^#\\s*6\\. RUN\\s*$", src)
  
  if (length(i_start) != 1L || length(i_end) != 1L || i_end <= i_start) {
    stop("Could not locate the Section 12 helper block.\n",
         "  'SECTION 12 -- COUNTY DEMOGRAPHIC HETEROGENEITY' matches: ",
         length(i_start), "\n  '# 6. RUN' matches: ", length(i_end), "\n",
         "  Expected exactly one of each, in that order. The pipeline has been\n",
         "  edited past what this function can recognize; load the .t12_*\n",
         "  definitions by hand rather than risk evaluating the wrong lines.",
         call. = FALSE)
  }
  
  # Stop four lines short of "# 6. RUN" to clear its banner rule, leaving the
  # last statement as .t12_driver()'s closing brace.
  block <- src[seq(i_start, i_end - 4L)]
  
  cat("Loading .t12_ helpers from lines ", i_start, "-", i_end - 4L, "\n", sep = "")
  eval(parse(text = block), envir = globalenv())
  
  cat("Available: .t12_build_ses_index(), .t12_fit_triple(), .t12_driver(),",
      "\n  .t12_preflight(), and the T12T_* config objects.",
      "\nNothing has been run. Typical next step:",
      "\n  ses_panel <- .t12_build_ses_index(outpatient)\n")
  invisible(TRUE)
}


cat("warm_start() and load_t12_helpers() defined.\n",
    "  warm_start()        load a completed run into this session\n",
    "  cache_status()      see what is on disk before loading\n",
    "  restore_session()   read the cache directly, if definitions are already in scope\n",
    sep = "")