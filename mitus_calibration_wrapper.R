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
#' Calibration Design:
#'   The posterior surface for state- and sub-state locations is multi-modal, so a
#'   single gradient-based descent lands in whichever local optimum its starting
#'   point drains into. This wrapper therefore runs a multi-start calibration:
#'   `n_runs` independent optimizations, each seeded from a different row of
#'   StartVal_st and each running the full BFGS/Nelder-Mead sequence of
#'   optim_b_st(). Runs are executed in parallel across `n_parallel` workers. The
#'   maximum a posteriori (MAP) parameter set -- the lowest -log posterior among
#'   the runs' final rounds -- is then used for the simulation and Tabby2 outputs.
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
#'   - optim_runs/: per-round optimizer results (Opt_*.rda) and per-run logs
#'   - {LOC}_Optim_all_{n_runs}_{MMDD}.rds: parameters and -log posterior for
#'     every run, in the layout the MITUS optim_data()/calib_plots_locs() tooling
#'     expects
#'   - {LOC}_Par_calibrated.rds, {LOC}_Par_optim_space.rds, bc.array, and the
#'     Tabby2 calibration output files (mirrored into tabby2_outputs/)
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

#' Posterior values at or above this come from the -10^12 likelihood penalty,
#' i.e. the run never reached a region where the model produced valid output.
FAIL_VALUE <- 1e11

#' Run one optimization from a single starting value, with all optimizer output
#' captured to its own log file so parallel runs do not interleave.
run_optimization <- function(samp, start_vals, loc, TB, n_cores, calib_end_year,
                            run_dir) {
  log_file <- file.path(run_dir, sprintf("optim_run_%03d.log", samp))
  cat(sprintf("[run %d] started at %s\n", samp, format(Sys.time())), file = stderr())
  t0 <- Sys.time()

  con <- file(log_file, open = "wt")
  sink(con, type = "output")
  sink(con, type = "message")
  status <- tryCatch({
    optim_b_st(df = start_vals, samp_i = samp, n_cores = n_cores, loc = loc,
               TB = TB, calib_end_year = calib_end_year)
    "completed"
  }, error = function(e) paste("error:", conditionMessage(e)))
  sink(type = "message")
  sink(type = "output")
  close(con)

  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "hours"))
  cat(sprintf("[run %d] %s after %.2f hours\n", samp, status, elapsed), file = stderr())

  list(samp_i = samp, status = status, elapsed_hours = elapsed,
       log_file = file.path(basename(run_dir), basename(log_file)))
}

#' Read every optimizer round saved by one run and return them ordered by round.
collect_rounds <- function(samp, loc, run_dir) {
  pattern <- sprintf("^Opt_%s_r([0-9]+)_%d_.*\\.rda$", loc, samp)
  files <- list.files(run_dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) return(list())

  rounds <- as.integer(sub(pattern, "\\1", basename(files)))
  files <- files[order(rounds)]
  rounds <- sort(rounds)

  collected <- lapply(seq_along(files), function(k) {
    env <- new.env()
    ok <- tryCatch({ load(files[k], envir = env); TRUE }, error = function(e) FALSE)
    obj_name <- paste0("o", rounds[k])
    if (!ok || !exists(obj_name, envir = env)) return(NULL)
    o <- get(obj_name, envir = env)
    list(round = rounds[k], value = as.numeric(o$value),
         convergence = as.integer(o$convergence), par = as.numeric(o$par))
  })
  collected[!vapply(collected, is.null, logical(1))]
}

#' Final round of a run, or NULL if it saved nothing. The optimization rounds
#' form a descending chain, so the last round a run reached is its best.
final_round <- function(rounds) {
  if (length(rounds) == 0) return(NULL)
  rounds[[length(rounds)]]
}

#' Convert a calibrated parameter vector from the N(0,1) optimization space into
#' a full model-space parameter vector, following gen_par_matrix().
optim_to_model_par <- function(par_optim) {
  names(par_optim) <- rownames(ParamInitZ)
  par_unif <- pnorm(par_optim, 0, 1)
  par_true <- par_unif
  par_true[idZ0] <- qbeta(par_unif[idZ0], shape1 = ParamInitZ[idZ0, 6],
                          shape2 = ParamInitZ[idZ0, 7])
  par_true[idZ1] <- qgamma(par_unif[idZ1], shape = ParamInitZ[idZ1, 6],
                           rate = ParamInitZ[idZ1, 7])
  par_true[idZ2] <- qnorm(par_unif[idZ2], mean = ParamInitZ[idZ2, 6],
                          sd = ParamInitZ[idZ2, 7])
  par_full <- P
  par_full[ii] <- par_true
  par_full
}

#' Store each run's calibrated parameters in the matrix layout that optim_data()
#' and calib_plots_locs() read (rows b_no_*, last column post_val).
write_optim_all <- function(finals, samp_ids, loc, n_runs, output_dir) {
  n_par <- nrow(ParamInitZ)
  n_row <- max(samp_ids)
  opt_all <- matrix(NA_real_, n_row, n_par + 1,
                    dimnames = list(paste0("b_no_", seq_len(n_row)),
                                    c(rownames(ParamInitZ), "post_val")))
  for (k in seq_along(samp_ids)) {
    if (is.null(finals[[k]])) next
    opt_all[samp_ids[k], seq_len(n_par)] <- finals[[k]]$par
    opt_all[samp_ids[k], n_par + 1] <- finals[[k]]$value
  }
  opt_file <- file.path(output_dir, sprintf("%s_Optim_all_%d_%s.rds", loc, n_runs,
                                            format(Sys.Date(), "%m%d")))
  saveRDS(opt_all, opt_file, version = 2)
  cat("Wrote optimization summary matrix:", basename(opt_file), "\n")
  opt_all
}

#' Run n_runs multi-start optimizations, collect their rounds, and pick the MAP.
run_calibration <- function(loc, samp_i, n_runs, n_parallel, n_cores, TB,
                            calib_end_year, output_dir) {
  n_avail <- nrow(StartVal_st)
  if (samp_i + n_runs - 1 > n_avail) {
    stop("n_runs=", n_runs, " starting at row ", samp_i, " needs ",
         samp_i + n_runs - 1, " starting values, but StartVal_st supplies ",
         n_avail, ". Lower n_runs or samp_i, or supply a StartVal file with ",
         "more rows (see gen_st_val_st).")
  }
  samp_ids <- samp_i + seq_len(n_runs) - 1

  run_dir <- file.path(output_dir, "optim_runs")
  dir.create(run_dir, showWarnings = FALSE, recursive = TRUE)
  # Absolute, so the per-run log paths stay valid once the workers are running
  # with run_dir as their working directory.
  run_dir <- normalizePath(run_dir)

  workers <- max(1, min(n_parallel, n_runs))
  cat("Starting", n_runs, "optimization runs from StartVal_st rows",
      paste(range(samp_ids), collapse = "-"), "on", workers, "workers\n")
  cat("Per-run optimizer output goes to", run_dir, "\n\n")
  flush.console()

  # optim_b_st() writes its Opt_*.rda files to the working directory; children
  # forked by mclapply() inherit it, and the file names carry samp_i so the
  # parallel runs cannot collide.
  old_wd <- setwd(run_dir)
  on.exit(setwd(old_wd), add = TRUE)

  t0 <- Sys.time()
  if (workers > 1 && .Platform$OS.type != "windows") {
    run_info <- mclapply(samp_ids, run_optimization, start_vals = StartVal_st,
                         loc = loc, TB = TB, n_cores = n_cores,
                         calib_end_year = calib_end_year, run_dir = run_dir,
                         mc.cores = workers, mc.preschedule = FALSE)
  } else {
    run_info <- lapply(samp_ids, run_optimization, start_vals = StartVal_st,
                       loc = loc, TB = TB, n_cores = n_cores,
                       calib_end_year = calib_end_year, run_dir = run_dir)
  }
  setwd(old_wd)
  cat(sprintf("\nAll runs finished in %.2f hours\n",
              as.numeric(difftime(Sys.time(), t0, units = "hours"))))

  # A worker that died leaves NULL or a try-error in place of its result; the
  # rounds it saved before dying are still on disk and still usable.
  rounds <- lapply(samp_ids, collect_rounds, loc = loc, run_dir = run_dir)
  finals <- lapply(rounds, final_round)

  runs <- lapply(seq_along(samp_ids), function(k) {
    info <- run_info[[k]]
    if (!is.list(info)) {
      info <- list(samp_i = samp_ids[k], status = "worker failed",
                   elapsed_hours = NA_real_,
                   log_file = file.path("optim_runs",
                                        sprintf("optim_run_%03d.log", samp_ids[k])))
    }
    entry <- list(samp_i = samp_ids[k], status = info$status,
                  elapsed_hours = round(info$elapsed_hours, 3),
                  rounds_saved = length(rounds[[k]]), log_file = info$log_file)
    if (!is.null(finals[[k]])) {
      entry$final_round <- finals[[k]]$round
      entry$neg_log_posterior <- finals[[k]]$value
      entry$convergence <- finals[[k]]$convergence
      entry$usable <- is.finite(finals[[k]]$value) && finals[[k]]$value < FAIL_VALUE
    } else {
      entry$usable <- FALSE
    }
    entry
  })

  write_optim_all(finals, samp_ids, loc, n_runs, output_dir)

  cat("\n--- Run summary (-log posterior, lower is better) ---\n")
  for (entry in runs) {
    cat(sprintf("  run %3d: %-12s rounds=%d  value=%s  convergence=%s\n",
                entry$samp_i, entry$status, entry$rounds_saved,
                if (is.null(entry$neg_log_posterior)) "none"
                else format(entry$neg_log_posterior, digits = 8),
                if (is.null(entry$convergence)) "NA" else entry$convergence))
  }

  usable <- vapply(runs, function(entry) isTRUE(entry$usable), logical(1))
  if (!any(usable)) {
    return(list(runs = runs, map = NULL, samp_ids = samp_ids))
  }
  values <- rep(Inf, length(runs))
  values[usable] <- vapply(runs[usable], function(entry) entry$neg_log_posterior,
                           numeric(1))
  map_idx <- which.min(values)
  cat(sprintf("\nMAP parameter set: run %d, round %d, -log posterior %s\n",
              runs[[map_idx]]$samp_i, runs[[map_idx]]$final_round,
              format(runs[[map_idx]]$neg_log_posterior, digits = 8)))

  list(runs = runs, samp_ids = samp_ids,
       map = list(samp_i = runs[[map_idx]]$samp_i,
                  round = runs[[map_idx]]$final_round,
                  neg_log_posterior = runs[[map_idx]]$neg_log_posterior,
                  convergence = runs[[map_idx]]$convergence,
                  par = finals[[map_idx]]$par),
       n_usable = sum(usable))
}

#' Simulate the MAP parameter set and write the Tabby2 calibration outputs.
write_map_outputs <- function(map, loc, output_dir) {
  par_model <- optim_to_model_par(map$par)

  par_optim <- matrix(map$par, nrow = 1,
                      dimnames = list(NULL, rownames(ParamInitZ)))
  saveRDS(par_optim, file.path(output_dir, paste0(loc, "_Par_optim_space.rds")),
          version = 2)

  # OutputsZint() takes parameter names from names(ParMatrix) whenever the matrix
  # has a single row, which is NULL for a matrix and strips the names param_init()
  # needs. A two-row matrix takes the colnames branch instead; samp_i=1 selects
  # the MAP row.
  ParMatrix <- matrix(par_model, nrow = 2, ncol = length(par_model), byrow = TRUE)
  colnames(ParMatrix) <- names(par_model)
  saveRDS(ParMatrix, file.path(output_dir, paste0(loc, "_Par_calibrated.rds")),
          version = 2)

  cat("\nRunning the model with the MAP parameters (1950-2050)...\n")
  bc_array <- OutputsZint(samp_i = 1, ParMatrix = ParMatrix, loc = loc,
                          startyr = 1950, endyr = 2050,
                          prg_chng = def_prgchng(ParMatrix[1, ]),
                          ttt_list = def_ttt())
  saveRDS(bc_array, file.path(output_dir, paste0(loc, "_bc_array_calibrated.rds")),
          version = 2)

  # model_calib_outputs() writes to ~/MITUS/inst/{loc}/calibration_outputs; create
  # that path and mirror the files back into the ResilientSims output directory.
  simp_date <- format(Sys.Date(), "%Y-%m-%d")
  tabby2_dir <- path.expand(file.path("~/MITUS/inst", loc, "calibration_outputs"))
  dir.create(tabby2_dir, showWarnings = FALSE, recursive = TRUE)
  model_calib_outputs(loc = loc, bc.array = bc_array, samp_i = 1,
                      simp.date = simp_date)

  mirror_dir <- file.path(output_dir, "tabby2_outputs")
  dir.create(mirror_dir, showWarnings = FALSE, recursive = TRUE)
  tabby2_files <- list.files(tabby2_dir, full.names = TRUE)
  file.copy(tabby2_files, mirror_dir, overwrite = TRUE)
  cat("Tabby2 calibration outputs generated:", length(tabby2_files), "files\n")

  list(model_years = c(1950, 2050), tabby2_files = length(tabby2_files),
       calibrated_parameters = as.list(par_model))
}

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
