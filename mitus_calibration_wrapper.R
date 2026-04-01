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
#'   - Code files: All MITUS R source files (119 files) + C++ src files (19 files)
#'   - Dependencies: mvtnorm, mnormt, parallel, lhs, Rcpp, MCMCpack, MASS, jsonlite
#'   - Execution command: Rscript mitus_calibration_wrapper.R
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
#' Runtime Parameters (from JSON):
#'   - loc: Location code (e.g., "CA", "NY", "US")
#'   - samp_i: Sample index (default: 1)
#'   - n_cores: Parallel cores (default: 2)
#'   - TB: Include TB likelihoods (default: 1)
#'   - optimize: Run calibration (true) or validation only (false)
#'
#' Output:
#'   - results.json with summary stats and calibration outputs
#'   - Optional: Optimization result files (Opt_*.rda)
#'   - Optional: Tabby2 calibration output files (.rds)
#'
#' Expected Runtime:
#'   - Validation mode (optimize=false): ~40 seconds
#'   - Full calibration (optimize=true): Several hours (7-9 optimization rounds)
#'
#' === ENVIRONMENT VARIABLES ===
#'   SIMULATOR_CODE_DIR       - Directory containing MITUS R source code
#'   SIMULATOR_INPUT_DIR      - Directory containing configuration data files
#'   SIMULATOR_OUTPUT_DIR     - Directory to write results
#'   SIMULATOR_RUNTIME_PARAMS_FILE - JSON file with calibration parameters

library(jsonlite)

main <- function() {
  # Read environment variables
  code_dir <- Sys.getenv("SIMULATOR_CODE_DIR", "/simulator/code")
  input_dir <- Sys.getenv("SIMULATOR_INPUT_DIR", "/simulator/input")
  output_dir <- Sys.getenv("SIMULATOR_OUTPUT_DIR", "/simulator/output")
  params_file <- Sys.getenv("SIMULATOR_RUNTIME_PARAMS_FILE", "")

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

  samp_i <- if (!is.null(params$samp_i)) params$samp_i else 1
  n_cores <- if (!is.null(params$n_cores)) params$n_cores else 2
  TB <- if (!is.null(params$TB)) params$TB else 1
  optimize <- if (!is.null(params$optimize)) params$optimize else TRUE

  cat("\nCalibration settings:\n")
  cat("  Location:", loc, "\n")
  cat("  Sample index:", samp_i, "\n")
  cat("  Cores:", n_cores, "\n")
  cat("  TB likelihoods:", TB, "\n")
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
    cat("\n=== Running Calibration Optimization ===\n")
    cat("This may take several hours...\n\n")

    # Change to output directory for optimization output files
    setwd(output_dir)

    # Run optimization
    tryCatch({
      optim_b_st(df = StartVal_st, samp_i = samp_i, n_cores = n_cores,
                 loc = loc, TB = TB)

      cat("\nOptimization completed successfully\n")

      # Load optimization results
      opt_files <- list.files(output_dir, pattern = paste0("Opt_", loc, "_r.*\\.rda$"),
                              full.names = TRUE)
      cat("Generated optimization files:\n")
      print(opt_files)

      # Load final optimization result
      if (length(opt_files) > 0) {
        # Load the last optimization round (highest number)
        opt_rounds <- sub(paste0("Opt_", loc, "_r(\\d+)_.*"), "\\1", basename(opt_files))
        final_opt_file <- opt_files[which.max(as.numeric(opt_rounds))]
        cat("\nLoading final optimization:", final_opt_file, "\n")

        load(final_opt_file)
        # optim_b_st creates objects o1, o2, ..., o8
        # Find the last one
        opt_obj_name <- sub(".*/Opt_.*_r(\\d+)_.*\\.rda", "o\\1", final_opt_file)
        final_params <- get(opt_obj_name)

        results$summary <- list(
          location = loc,
          sample_index = samp_i,
          optimization_complete = TRUE,
          final_likelihood = -final_params$value,
          convergence = final_params$convergence,
          n_iterations = final_params$counts,
          optimized_parameters = as.list(final_params$par)
        )
      }

      # Generate calibration validation plots and outputs
      cat("\nGenerating calibration validation outputs...\n")

      # Run model with optimized parameters
      ParMatrix <- matrix(final_params$par, nrow = 1)
      # Get parameter names from StartVal_st which was used in optimization
      if (!is.null(names(final_params$par))) {
        colnames(ParMatrix) <- names(final_params$par)
      } else if (exists("StartVal_st") && !is.null(colnames(StartVal_st))) {
        colnames(ParMatrix) <- colnames(StartVal_st)
      } else {
        colnames(ParMatrix) <- names(StartVal_st[samp_i,])
      }

      result <- OutputsZint(samp_i = 1, ParMatrix = ParMatrix, loc = loc,
                           startyr = 1950, endyr = 2050,
                           prg_chng = def_prgchng(ParMatrix[1,]),
                           ttt_list = def_ttt())

      # Generate calibration output files for Tabby2
      simp_date <- format(Sys.Date(), "%Y-%m-%d")
      model_calib_outputs(loc = loc, bc.array = result, samp_i = 1,
                         simp_date = simp_date)

      cat("Calibration outputs generated for Tabby2 integration\n")

    }, error = function(e) {
      cat("\nERROR during optimization:\n")
      cat(conditionMessage(e), "\n")
      results$summary <<- list(
        location = loc,
        optimization_complete = FALSE,
        error = conditionMessage(e)
      )
    })

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

  # Create list of output files for ResilientSims
  output_files <- list.files(output_dir, full.names = FALSE)
  results$output_manifest <- output_files[output_files != "results.json"]

  # Write results.json
  results_file <- file.path(output_dir, "results.json")
  cat("\nWriting results to:", results_file, "\n")
  write_json(results, results_file, pretty = TRUE, auto_unbox = TRUE)

  cat("\n=== Calibration wrapper completed successfully ===\n")
}

# Run main function
if (!interactive()) {
  main()
}
