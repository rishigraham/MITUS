#!/usr/bin/env Rscript
#' MITUS wrapper for ResilientSims
#'
#' Entry point for every MITUS simulator on the ResilientSims platform. Reads the
#' environment and runtime parameters, installs or loads the MITUS package, loads
#' the location, runs the mode selected by the entrypoint argument, and writes
#' results.json according to ResilientSims specifications.
#'
#' === MODES (first command-line argument; default "calibrate") ===
#'
#'   Rscript resilientsims/wrapper.R calibrate
#'     Multi-start calibration, then the location package and the built-in
#'     scenarios from the MAP. Hours per run; the heavy simulator.
#'   Rscript resilientsims/wrapper.R package
#'     The location package and the built-in scenarios from the configuration's
#'     existing {LOC}_Param_*.rds, without calibrating. Minutes. Refreshes the
#'     precomputed outputs after a code change and exercises the shared code.
#'   Rscript resilientsims/wrapper.R scenario
#'     One custom scenario (care cascade changes and/or targeted testing and
#'     treatment) on the configuration's parameters, exported for the hub. Seconds
#'     of model time; the interactive simulator.
#'
#' Each simulator definition fixes its mode through the entrypoint command, so the
#' runtime parameters of one mode cannot start another.
#'
#' === REQUIREMENTS SUMMARY ===
#'
#' Simulator Requirements:
#'   - Language: R
#'   - Code files: the MITUS R/ and src/ directories plus this resilientsims/ directory
#'   - Dependencies: mvtnorm, mnormt, parallel, lhs, Rcpp, MCMCpack, MASS, jsonlite
#'   - Execution command: cd $SIMULATOR_CODE_DIR && Rscript resilientsims/wrapper.R <mode>
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
#'       - {LOC}_*.rds (location-specific data)
#'       - {LOC}_Param_*.rds (calibrated parameters; required by package and scenario)
#'       - {LOC}_results_1.rds (base case; required by scenario when targeted
#'         testing is active)
#'
#' File Layout:
#'   This file is the entry point: environment, runtime parameters, package
#'   install, model_load(), mode dispatch, and the results.json/output manifest
#'   contract. Step-specific code is sourced from this directory:
#'     - calibration.R: multi-start optimization and the MAP parameters
#'     - package.R: the location package (parameters, base case, Tabby2 comparison
#'       files) and the built-in scenario export, shared by calibrate and package
#'     - scenarios.R: built-in and custom scenario runs, request validation
#'     - tabby2_export.R: the v2 JSON the ResilientHub2 explorer reads
#'
#' Runtime Parameters (from JSON):
#'   all modes
#'   - loc: Location code (e.g., "CA", "NY", "US", "SanDiego")
#'   - n_cores: Cores for one model run (default: 1)
#'   - n_parallel: Concurrent optimization runs, or concurrent built-in scenario
#'       runs (default: 2)
#'   calibrate
#'   - n_runs: Independent optimization runs to launch (default: 15). Capped by
#'       the number of starting values in StartVal_st.
#'   - samp_i: First StartVal_st row to use; runs use rows samp_i..samp_i+n_runs-1
#'       (default: 1)
#'   - TB: Include TB likelihoods (default: 1)
#'   - calib_end_year: Last year of calibration targets (default: 2021)
#'   scenario
#'   - name: Scenario name shown in the hub (required)
#'   - prg_chng: care cascade fields (any subset): start_yr, scrn_cov, IGRA_frc,
#'       ltbi_init_frc, frc_3hp, comp_3hp, frc_4r, comp_4r, frc_3hr, comp_3hr,
#'       tb_tim2tx_frc, tb_txdef_frc. Fractions as fractions; tb_tim2tx_frc as a
#'       percent of the current duration of infectiousness. Omitted fields keep
#'       the model defaults; the three regimen fractions must sum to 1.
#'   - ttt_list: targeted testing fields (any subset): NativityGrp (All|USB|NUSB),
#'       AgeGrp (All|"0 to 24"|"25 to 64"|"65+"), NRiskGrp (millions of people),
#'       FrcScrn (0-1 screened per year), StartYr, EndYr, RRprg, RRmu, RRPrev.
#'       Targeting is active when NRiskGrp and FrcScrn are both nonzero.
#'   - n_param_sets: Parameter sets to simulate (default: 2)
#'
#' Output (all modes):
#'   - results.json: summary and the output manifest
#'   - {LOC}_{mode}_summary.json: the full results, as an output file that
#'     ResilientSims uploads with the other outputs (results.json itself is not)
#' Output (calibrate):
#'   - optim_runs/: per-round optimizer results (Opt_*.rda) and per-run logs
#'   - {LOC}_Optim_all_{n_runs}_{MMDD}.rds: parameters and -log posterior for
#'     every run, in the layout the MITUS optim_data()/calib_plots_locs() tooling
#'     expects
#'   - {LOC}_Par_optim_space.rds: the MAP in the N(0,1) optimization space
#' Output (calibrate and package):
#'   - {LOC}_Param_{YYYY-MM-DD}.rds: the parameters in model space as two identical
#'     rows, named so model_load() resolves it as the location's parameter file
#'   - {LOC}_results_1.rds: base-case projection 1950-2050 for both rows
#'     (parameter set x year x output), read by targeted testing scenarios
#'   - tabby2_outputs/: the Tabby2 calibration comparison files
#'   - {LOC}_tabby2_ui_data.json: the built-in scenarios (base case, five
#'     interventions, two sensitivity analyses) as the hub's v2 series payload
#' Output (scenario):
#'   - {LOC}_scenario_results.json: the scenario's series in the same v2 schema,
#'     with the parameter file it ran against as the location basis and the
#'     request echoed under meta.custom_scenario
#'
#' Expected Runtime:
#'   - Each optimization run: ~6-11 hours for a state or sub-state location
#'     (roughly 13,000 likelihood evaluations at ~2.3 seconds each)
#'   - Full calibration: roughly ceiling(n_runs / n_parallel) x per-run time,
#'     plus a few minutes for the location package and the built-in scenarios
#'   - package: a few minutes; scenario: seconds of model time
#'   - Installing the package from source adds 1-4 minutes to every run unless
#'     the image has it preinstalled (MITUS_PREINSTALLED=1)
#'   - Memory: ~430 MB for the loaded model plus at most ~230 MB per concurrent run
#'
#' === ENVIRONMENT VARIABLES ===
#'   SIMULATOR_CODE_DIR       - Directory containing MITUS R source code
#'   SIMULATOR_INPUT_DIR      - Directory containing configuration data files
#'   SIMULATOR_OUTPUT_DIR     - Directory to write results
#'   SIMULATOR_RUNTIME_PARAMS_FILE - JSON file with runtime parameters
#'   MITUS_PREINSTALLED       - "1" when the image already has MITUS installed
#'                              from this code; skips the source install

library(jsonlite)
library(parallel)

MODES <- c("calibrate", "package", "scenario")

#' Directory holding this script and its companion files (Rscript passes --file=).
wrapper_dir <- function() {
  f <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(f) > 0) dirname(normalizePath(f[1])) else file.path(getwd(), "resilientsims")
}

#' A positive integer runtime parameter, or its default.
int_param <- function(params, name, default) {
  if (is.null(params[[name]])) return(as.integer(default))
  v <- suppressWarnings(as.integer(params[[name]]))
  if (length(v) != 1 || is.na(v) || v < 1) stop("Parameter '", name, "' must be a positive integer")
  v
}

#' The location's parameter matrix from the configuration, or stop.
configured_parameters <- function(loc, input_dir) {
  if (!exists("Par") || is.null(Par)) {
    stop("No parameter file (", loc, "_Param_*.rds) in the configuration for location: ", loc)
  }
  list(ParMatrix = par_matrix_from(Par),
       source = basename(find_loc_file(loc, "Param", data_dir = input_dir)))
}

run_calibrate <- function(params, loc, output_dir, n_cores, n_parallel) {
  samp_i <- int_param(params, "samp_i", 1)
  n_runs <- int_param(params, "n_runs", 15)
  TB <- if (!is.null(params$TB)) params$TB else 1
  calib_end_year <- int_param(params, "calib_end_year", 2021)

  cat("  Optimization runs:", n_runs, "\n")
  cat("  First starting value:", samp_i, "\n")
  cat("  TB likelihoods:", TB, "\n")
  cat("  Calibration end year:", calib_end_year, "\n")
  cat("\n=== Running Multi-Start Calibration ===\n")
  cat("This may take days for a state or sub-state location...\n\n")

  calibration <- run_calibration(loc = loc, samp_i = samp_i, n_runs = n_runs,
                                 n_parallel = n_parallel, n_cores = n_cores,
                                 TB = TB, calib_end_year = calib_end_year,
                                 output_dir = output_dir)
  results <- list(runs = calibration$runs)
  results$summary <- list(
    mode = "calibrate", location = loc, n_runs = n_runs, first_start_value = samp_i,
    n_parallel = n_parallel, n_cores = n_cores, TB = TB, calib_end_year = calib_end_year,
    runs_usable = if (is.null(calibration$map)) 0 else calibration$n_usable,
    optimization_complete = !is.null(calibration$map))

  if (is.null(calibration$map)) {
    cat("\nERROR: no run produced a usable parameter set;",
        "check the per-run logs in optim_runs/\n")
    results$summary$error <- "No optimization run reached a usable posterior value"
    return(results)
  }
  results$summary$map <- calibration$map[c("samp_i", "round", "neg_log_posterior",
                                           "convergence")]
  results$map_parameters_optim_space <- as.list(setNames(calibration$map$par,
                                                         rownames(ParamInitZ)))
  ParMatrix <- write_map_parameters(calibration$map, loc, output_dir)
  results$map_parameters <- as.list(setNames(ParMatrix[1, ], colnames(ParMatrix)))
  basis <- sprintf("MAP of %d-start calibration, %s", n_runs, format(Sys.Date(), "%Y-%m-%d"))
  c(results, build_location_package(ParMatrix, loc, output_dir, basis, n_parallel))
}

run_package <- function(params, loc, input_dir, output_dir, n_parallel) {
  cfg <- configured_parameters(loc, input_dir)
  cat("  Parameter file:", cfg$source, "\n")
  cat("\n=== Building the location package ===\n")
  results <- list(summary = list(mode = "package", location = loc,
                                 parameter_file = cfg$source, n_parallel = n_parallel))
  c(results, build_location_package(cfg$ParMatrix, loc, output_dir,
                                     basis = paste("parameter file", cfg$source), n_parallel))
}

#' Location package plus built-in scenario export; errors are recorded, not fatal.
build_location_package <- function(ParMatrix, loc, output_dir, basis, n_parallel) {
  results <- list()
  pkg <- tryCatch(write_location_package(ParMatrix, loc, output_dir), error = function(e) {
    cat("\nERROR building the location package:\n", conditionMessage(e), "\n")
    list(error = conditionMessage(e))
  })
  results$package_outputs <- pkg[setdiff(names(pkg), "base_case")]
  if (is.null(pkg$error)) {
    results$scenario_outputs <- tryCatch({
      export_builtin_scenarios(loc, ParMatrix, pkg$base_case, output_dir, basis, n_parallel)
    }, error = function(e) {
      cat("\nERROR generating the built-in scenarios:\n", conditionMessage(e), "\n")
      list(error = conditionMessage(e))
    })
  }
  results
}

run_scenario_mode <- function(params, loc, input_dir, output_dir, n_cores) {
  cfg <- configured_parameters(loc, input_dir)
  spec <- scenario_spec(params, cfg$ParMatrix)
  cat("  Parameter file:", cfg$source, "\n")
  cat("  Scenario:", spec$name, "(", spec$id, ")\n")
  cat("  Parameter sets:", min(spec$n_param_sets, nrow(cfg$ParMatrix)), "\n")
  cat("  Targeted testing active:", spec$ttt_active, "\n")
  cat("\n=== Running the scenario ===\n")

  t0 <- Sys.time()
  arr <- run_custom_scenario(loc, cfg$ParMatrix, spec, n_cores = n_cores)
  payload <- custom_scenario_payload(loc, spec, arr, cfg$source, cfg$ParMatrix)
  json_file <- paste0(loc, "_scenario_results.json")
  write_tabby2_ui_json(payload, file.path(output_dir, json_file))
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("Wrote %s: %d series in %.1f s\n", json_file, length(payload$series), elapsed))

  list(summary = list(mode = "scenario", location = loc, scenario = spec$name,
                      scenario_id = spec$id, parameter_file = cfg$source,
                      parameter_sets = dim(arr)[1], ttt_active = spec$ttt_active,
                      results_file = json_file, elapsed_seconds = round(elapsed, 1)),
       scenario = payload$meta$custom_scenario,
       scenario_description = payload$meta$scenarios[[spec$id]]$description)
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  mode <- if (length(args) >= 1) args[1] else "calibrate"
  if (!mode %in% MODES) stop("Unknown mode '", mode, "'; expected one of ", paste(MODES, collapse = ", "))

  code_dir <- Sys.getenv("SIMULATOR_CODE_DIR", "/simulator/code")
  input_dir <- Sys.getenv("SIMULATOR_INPUT_DIR", "/simulator/input")
  output_dir <- Sys.getenv("SIMULATOR_OUTPUT_DIR", "/simulator/output")
  params_file <- Sys.getenv("SIMULATOR_RUNTIME_PARAMS_FILE", "")

  for (f in c("calibration.R", "package.R", "scenarios.R", "tabby2_export.R")) {
    source(file.path(wrapper_dir(), f))
  }

  cat("=== MITUS wrapper:", mode, "===\n")
  cat("Code directory:", code_dir, "\n")
  cat("Input directory:", input_dir, "\n")
  cat("Output directory:", output_dir, "\n")
  cat("Parameters file:", params_file, "\n\n")

  if (params_file != "" && file.exists(params_file)) {
    params <- fromJSON(params_file, simplifyVector = FALSE)
    cat("Loaded parameters:\n")
    cat(toJSON(params, auto_unbox = TRUE, pretty = TRUE), "\n")
  } else {
    stop("Runtime parameters file not found: ", params_file)
  }

  loc <- params$loc
  if (is.null(loc) || loc == "") stop("Parameter 'loc' (location code) is required")
  if (!is.null(params$optimize)) {
    stop("Parameter 'optimize' is no longer used; the run mode is the entrypoint argument")
  }
  n_cores <- int_param(params, "n_cores", 1)
  n_parallel <- int_param(params, "n_parallel", 2)

  cat("\nSettings:\n")
  cat("  Mode:", mode, "\n")
  cat("  Location:", loc, "\n")
  cat("  Cores per model run:", n_cores, "\n")
  cat("  Concurrent runs:", n_parallel, "\n")

  if (identical(Sys.getenv("MITUS_PREINSTALLED"), "1")) {
    cat("\nUsing the preinstalled MITUS package\n")
  } else {
    cat("\nInstalling MITUS package...\n")
    lib_path <- Sys.getenv("R_LIBS_USER")
    if (nchar(lib_path) == 0 || !dir.exists(lib_path)) lib_path <- .libPaths()[1]
    install.packages(code_dir, repos = NULL, type = "source", lib = lib_path)
  }
  cat("Loading MITUS package...\n")
  library(MITUS)

  cat("\nLoading MITUS model for location:", loc, "\n")
  model_load(loc = loc, data_dir = input_dir)
  cat("MITUS model loaded successfully\n")
  if (exists("Par") && !is.null(Par)) cat("  Par dimensions:", dim(Par), "\n")

  results <- switch(mode,
    calibrate = run_calibrate(params, loc, output_dir, n_cores, n_parallel),
    package   = run_package(params, loc, input_dir, output_dir, n_parallel),
    scenario  = run_scenario_mode(params, loc, input_dir, output_dir, n_cores))

  # ResilientSims keeps only results$summary from results.json and uploads the
  # files named in output_manifest, so the full results are written to their own
  # file, before the manifest is built, to reach the output store.
  summary_file <- file.path(output_dir, paste0(loc, "_", mode, "_summary.json"))
  write_json(results, summary_file, pretty = TRUE, auto_unbox = TRUE, digits = NA)
  cat("\nWrote run summary:", basename(summary_file), "\n")

  output_files <- list.files(output_dir, recursive = TRUE, full.names = FALSE)
  results$output_manifest <- output_files[output_files != "results.json"]

  results_file <- file.path(output_dir, "results.json")
  cat("Writing results to:", results_file, "\n")
  write_json(results, results_file, pretty = TRUE, auto_unbox = TRUE, digits = NA)

  cat("\n=== MITUS wrapper completed successfully ===\n")
}

if (!interactive()) {
  main()
}
