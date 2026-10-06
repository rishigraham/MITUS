#' MITUS calibration for ResilientSims: multi-start optimization and MAP outputs
#'
#' Sourced by wrapper.R from this directory.
#' Expects the MITUS package loaded and model_load() already run for the location.
#'
#' Calibration Design:
#'   The posterior surface for state- and sub-state locations is multi-modal, so a
#'   single gradient-based descent lands in whichever local optimum its starting
#'   point drains into. run_calibration() therefore runs a multi-start calibration:
#'   `n_runs` independent optimizations, each seeded from a different row of
#'   StartVal_st and each running the full BFGS/Nelder-Mead sequence of
#'   optim_b_st(). Runs are executed in parallel across `n_parallel` workers. The
#'   maximum a posteriori (MAP) parameter set -- the lowest -log posterior among
#'   the runs' final rounds -- is converted to model space by
#'   write_map_parameters() and handed to the location package (package.R).

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

#' Save the MAP in optimization space and return it in model space as the
#' two-row parameter matrix the location package and scenario runs use.
write_map_parameters <- function(map, loc, output_dir) {
  par_optim <- matrix(map$par, nrow = 1,
                      dimnames = list(NULL, rownames(ParamInitZ)))
  saveRDS(par_optim, file.path(output_dir, paste0(loc, "_Par_optim_space.rds")),
          version = 2)
  par_matrix_from(optim_to_model_par(map$par))
}
