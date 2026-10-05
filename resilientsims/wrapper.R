#!/usr/bin/env Rscript
#' MITUS Calibration Wrapper for ResilientSims
#'
#' This script integrates MITUS TB model calibration into the ResilientSims platform.
#' It reads environment variables, processes input files, runs calibration, and
#' generates results.json according to ResilientSims specifications.
#'
#' === REQUIREMENTS SUMMARY ===
#'
#' Simulator Requirements:
#'   - Language: R
#'   - Code files: the MITUS R/ and src/ directories plus this resilientsims/ directory
#'   - Dependencies: mvtnorm, mnormt, parallel, lhs, Rcpp, MCMCpack, MASS, jsonlite
#'   - Execution command: Rscript resilientsims/wrapper.R (from the code directory)
#'   - Wrapper script: This file (handles environment setup, runs calibration)
#'
#' Configuration Requirements:
#'   Configuration must provide location-specific data files in this structure:
#'     ST/
#'       - ST_CalibDat_*.rds (calibration targets)
#'       - ST_ParamInit_*.rds (parameter initialization)
#'       - ST_StartVal_*.rds (optimization starting values)
#'       - ST_*.rds (other state-level data: immigration, mortality, population, etc.)
#'     {LOC}/  (e.g., CA/)
#'       - {LOC}_ModelInputs_*.rds
#'       - {LOC}_Param_*.rds (existing calibrated parameters, if available)
#'       - {LOC}_*.rds (location-specific data)
#'
#' File Layout:
#'   This file is the entry point: environment, runtime parameters, package
#'   install, model_load(), and the results.json/output manifest contract.
#'   Step-specific code is sourced from this directory:
#'     - calibration.R: multi-start optimization and the MAP output package
#'
#' Runtime Parameters (from JSON):
#'   - loc: Location code (e.g., "CA", "NY", "US")
#'   - n_runs: Independent optimization runs to launch (default: 15). Capped by
#'       the number of starting values in StartVal_st.
#'   - samp_i: First StartVal_st row to use; runs use rows samp_i..samp_i+n_runs-1
#'       (default: 1). In validation mode this selects the Par row to simulate.
#'   - n_parallel: Optimization runs to execute concurrently (default: 2)
#'   - n_cores: Cores passed to optim_b_st (default: 1)
#'   - TB: Include TB likelihoods (default: 1)
#'   - calib_end_year: Last year of calibration targets (default: 2021)
#'   - optimize: Run calibration (true) or validation only (false)
#'
#' Output:
#'   - results.json with per-run outcomes, the MAP selection, and summary stats
#'   - {LOC}_calibration_summary.json: the same content, as an output file that
#'     ResilientSims uploads with the other outputs (results.json itself is not)
#'   - optim_runs/: per-round optimizer results (Opt_*.rda) and per-run logs
#'   - {LOC}_Optim_all_{n_runs}_{MMDD}.rds: parameters and -log posterior for
#'     every run, in the layout the MITUS optim_data()/calib_plots_locs() tooling
#'     expects
#'   - {LOC}_Param_{YYYY-MM-DD}.rds: the MAP in model space as two identical rows,
#'     named so model_load() resolves it as the location's parameter file
#'   - {LOC}_results_1.rds: base-case projection 1950-2050 for both rows
#'     (parameter set x year x output), read by targeted testing scenarios
#'   - {LOC}_Par_optim_space.rds: the MAP in the N(0,1) optimization space
#'   - tabby2_outputs/: the Tabby2 calibration comparison files
#'
#' Expected Runtime:
#'   - Validation mode (optimize=false): ~40 seconds
#'   - Each optimization run: ~6-11 hours for a state or sub-state location
#'     (roughly 13,000 likelihood evaluations at ~2.3 seconds each)
#'   - Full calibration: roughly ceiling(n_runs / n_parallel) x per-run time
#'   - Memory: ~430 MB for the loaded model plus at most ~230 MB per concurrent run
#'
#' === ENVIRONMENT VARIABLES ===
#'   SIMULATOR_CODE_DIR       - Directory containing MITUS R source code
#'   SIMULATOR_INPUT_DIR      - Directory containing configuration data files
#'   SIMULATOR_OUTPUT_DIR     - Directory to write results
#'   SIMULATOR_RUNTIME_PARAMS_FILE - JSON file with calibration parameters

library(jsonlite)
library(parallel)

#' Directory holding this script and its companion files (Rscript passes --file=).
wrapper_dir <- function() {
  f <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(f) > 0) dirname(normalizePath(f[1])) else file.path(getwd(), "resilientsims")
}

main <- function() {
  # Read environment variables
  code_dir <- Sys.getenv("SIMULATOR_CODE_DIR", "/simulator/code")
  input_dir <- Sys.getenv("SIMULATOR_INPUT_DIR", "/simulator/input")
  output_dir <- Sys.getenv("SIMULATOR_OUTPUT_DIR", "/simulator/output")
  params_file <- Sys.getenv("SIMULATOR_RUNTIME_PARAMS_FILE", "")

  source(file.path(wrapper_dir(), "calibration.R"))

  cat("=== MITUS Calibration Wrapper ===\n")
  cat("Code directory:", code_dir, "\n")
  cat("Input directory:", input_dir, "\n")
  cat("Output directory:", output_dir, "\n")
  cat("Parameters file:", params_file, "\n\n")

  # Load runtime parameters
  if (params_file != "" && file.exists(params_file)) {
    params <- fromJSON(params_file)
    cat("Loaded parameters:\n")
    print(params)
  } else {
    stop("Runtime parameters file not found: ", params_file)
  }

  # Extract parameters with defaults
  loc <- params$loc
  if (is.null(loc) || loc == "") {
    stop("Parameter 'loc' (location code) is required")
  }

  samp_i <- if (!is.null(params$samp_i)) as.integer(params$samp_i) else 1L
  n_runs <- if (!is.null(params$n_runs)) as.integer(params$n_runs) else 15L
  n_parallel <- if (!is.null(params$n_parallel)) as.integer(params$n_parallel) else 2L
  n_cores <- if (!is.null(params$n_cores)) as.integer(params$n_cores) else 1L
  TB <- if (!is.null(params$TB)) params$TB else 1
  calib_end_year <- if (!is.null(params$calib_end_year)) as.integer(params$calib_end_year) else 2021L
  optimize <- if (!is.null(params$optimize)) params$optimize else TRUE

  if (is.na(samp_i) || samp_i < 1) stop("Parameter 'samp_i' must be a positive integer")
  if (is.na(n_runs) || n_runs < 1) stop("Parameter 'n_runs' must be a positive integer")
  if (is.na(n_parallel) || n_parallel < 1) stop("Parameter 'n_parallel' must be a positive integer")
  if (is.na(n_cores) || n_cores < 1) stop("Parameter 'n_cores' must be a positive integer")

  cat("\nCalibration settings:\n")
  cat("  Location:", loc, "\n")
  cat("  Optimization runs:", n_runs, "\n")
  cat("  First starting value:", samp_i, "\n")
  cat("  Concurrent runs:", n_parallel, "\n")
  cat("  Cores per run:", n_cores, "\n")
  cat("  TB likelihoods:", TB, "\n")
  cat("  Calibration end year:", calib_end_year, "\n")
  cat("  Run optimization:", optimize, "\n\n")

  # Install MITUS as a local package
  cat("\nInstalling MITUS package...\n")
  lib_path <- Sys.getenv("R_LIBS_USER")
  if (nchar(lib_path) == 0 || !dir.exists(lib_path)) {
    lib_path <- .libPaths()[1]
  }
  install.packages(code_dir, repos = NULL, type = "source", lib = lib_path)

  # Load MITUS package
  cat("Loading MITUS package...\n")
  library(MITUS)

  # Load MITUS data and initialize model for this location
  cat("\nLoading MITUS model for location:", loc, "\n")
  model_load(loc = loc, data_dir = input_dir)

  cat("\nMITUS model loaded successfully\n")
  cat("Global variables set:\n")
  cat("  ParamInit dimensions:", dim(ParamInit), "\n")
  cat("  StartVal dimensions:", dim(StartVal), "\n")
  if (exists("Par")) {
    cat("  Par dimensions:", dim(Par), "\n")
  }

  # Run calibration or validation
  results <- list()

  if (optimize) {
    cat("\n=== Running Multi-Start Calibration ===\n")
    cat("This may take days for a state or sub-state location...\n\n")

    calibration <- run_calibration(loc = loc, samp_i = samp_i, n_runs = n_runs,
                                   n_parallel = n_parallel, n_cores = n_cores,
                                   TB = TB, calib_end_year = calib_end_year,
                                   output_dir = output_dir)

    results$runs <- calibration$runs
    results$summary <- list(
      location = loc,
      calibration_mode = TRUE,
      n_runs = n_runs,
      first_start_value = samp_i,
      n_parallel = n_parallel,
      n_cores = n_cores,
      TB = TB,
      calib_end_year = calib_end_year,
      runs_usable = if (is.null(calibration$map)) 0 else calibration$n_usable,
      optimization_complete = !is.null(calibration$map)
    )

    if (is.null(calibration$map)) {
      cat("\nERROR: no run produced a usable parameter set;",
          "check the per-run logs in optim_runs/\n")
      results$summary$error <- "No optimization run reached a usable posterior value"
    } else {
      results$summary$map <- calibration$map[c("samp_i", "round",
                                               "neg_log_posterior", "convergence")]
      results$map_parameters_optim_space <- as.list(setNames(
        calibration$map$par, rownames(ParamInitZ)))

      map_outputs <- tryCatch({
        write_map_outputs(calibration$map, loc, output_dir)
      }, error = function(e) {
        cat("\nERROR generating MAP outputs:\n", conditionMessage(e), "\n")
        list(error = conditionMessage(e))
      })
      results$map_outputs <- map_outputs
    }

  } else {
    cat("\n=== Validation Mode (No Optimization) ===\n")
    cat("Testing model with existing parameters\n\n")

    # Run model with existing parameters
    if (!exists("Par")) {
      stop("No calibrated parameters (Par) available for location: ", loc)
    }

    result <- OutputsZint(samp_i = samp_i, ParMatrix = Par, loc = loc,
                         startyr = 1950, endyr = 2050,
                         prg_chng = def_prgchng(Par[samp_i,]),
                         ttt_list = def_ttt())

    cat("Model run completed successfully\n")

    results$summary <- list(
      location = loc,
      sample_index = samp_i,
      validation_mode = TRUE,
      model_years = c(1950, 2050)
    )
  }

  # ResilientSims keeps only results$summary from results.json and uploads the
  # files named in output_manifest, so the per-run table and MAP parameters are
  # written to their own file, before the manifest is built, to reach the output store.
  summary_file <- file.path(output_dir, paste0(loc, "_calibration_summary.json"))
  write_json(results, summary_file, pretty = TRUE, auto_unbox = TRUE, digits = NA)
  cat("\nWrote calibration summary:", basename(summary_file), "\n")

  # Create list of output files for ResilientSims
  output_files <- list.files(output_dir, recursive = TRUE, full.names = FALSE)
  results$output_manifest <- output_files[output_files != "results.json"]

  # Write results.json
  results_file <- file.path(output_dir, "results.json")
  cat("\nWriting results to:", results_file, "\n")
  write_json(results, results_file, pretty = TRUE, auto_unbox = TRUE, digits = NA)

  cat("\n=== Calibration wrapper completed successfully ===\n")
}

# Run main function
if (!interactive()) {
  main()
}
