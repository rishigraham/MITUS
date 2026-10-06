#' The per-location package a calibration produces and the scenario runs consume
#'
#' Sourced by wrapper.R from this directory. Expects the MITUS package loaded and
#' model_load() already run for the location. Shared by the calibrate mode (from
#' the MAP) and the package mode (from the location's existing parameter file).

#' Parameter matrix with at least two rows. OutputsZint() takes parameter names
#' from names(ParMatrix) whenever the matrix has a single row, which is NULL for
#' a matrix and strips the names param_init() needs, and the Tabby2 scenario path
#' simulates Par[1:2,]; a single parameter set is therefore stored twice.
par_matrix_from <- function(par) {
  m <- if (is.null(dim(par))) {
    matrix(par, nrow = 1, dimnames = list(NULL, names(par)))
  } else {
    as.matrix(par)
  }
  if (nrow(m) == 1) m <- m[c(1, 1), , drop = FALSE]
  rownames(m) <- NULL
  m
}

#' Write the location's parameter file, base-case results and Tabby2 comparison
#' files. Returns the file names and the base-case array.
write_location_package <- function(ParMatrix, loc, output_dir) {
  simp_date <- format(Sys.Date(), "%Y-%m-%d")
  # {loc}_Param_{date} is the pattern model_load() resolves, so a Configuration
  # can carry the file unchanged
  param_file <- paste0(loc, "_Param_", simp_date, ".rds")
  saveRDS(ParMatrix, file.path(output_dir, param_file), version = 2)

  cat("\nRunning the base case (1950-2050)...\n")
  bc_array <- run_scenario(loc, ParMatrix, builtin = "base_case")
  # read by targeted testing scenarios in param_init()
  results_file <- paste0(loc, "_results_1.rds")
  saveRDS(bc_array, file.path(output_dir, results_file), version = 2)

  tabby2_dir <- file.path(output_dir, "tabby2_outputs")
  model_calib_outputs(loc = loc, bc.array = bc_array[1, , ], samp_i = 1,
                      simp.date = simp_date, out_dir = tabby2_dir)
  tabby2_files <- list.files(tabby2_dir)
  cat("Tabby2 calibration outputs generated:", length(tabby2_files), "files\n")

  list(param_file = param_file, results_file = results_file,
       tabby2_files = length(tabby2_files), base_case = bc_array)
}

#' Run the built-in scenarios and write the hub's {loc}_tabby2_ui_data.json.
export_builtin_scenarios <- function(loc, ParMatrix, base_case, output_dir, basis,
                                     n_parallel = 1) {
  t0 <- Sys.time()
  cat("\nRunning the", length(TABBY2_SCENARIOS), "built-in scenarios...\n")
  raw <- run_builtin_scenarios(loc, ParMatrix, base_case = base_case,
                               n_parallel = n_parallel)
  payload <- tabby2_ui_payload(loc, raw, basis)
  json_file <- paste0(loc, "_tabby2_ui_data.json")
  write_tabby2_ui_json(payload, file.path(output_dir, json_file))
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("Wrote %s: %d series in %.0f s\n", json_file, length(payload$series),
              elapsed))
  list(ui_data_file = json_file, scenarios = names(raw),
       n_series = length(payload$series), parameter_sets = dim(raw[[1]])[1],
       elapsed_seconds = round(elapsed, 1))
}
