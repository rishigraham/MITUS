#' Tabby2-style scenario runs
#'
#' Sourced by wrapper.R from this directory. Expects the MITUS package loaded and
#' model_load() already run for the location. A scenario is either one of the
#' model's built-in interventions and sensitivity analyses (an Int1..Int5 or
#' Scen1..Scen3 flag) or a custom change expressed through the prg_chng vector
#' (care cascade) and ttt_list (targeted testing and treatment) that param_init()
#' consumes. Every run returns a parameter set x year x output array.

library(parallel)

#' Built-in scenarios exported to the hub. scenario_1 ("No Transmission Within
#' the US After 2022") is left out: its implementation is commented out in
#' param_init(), so its output equals the base case.
TABBY2_SCENARIOS <- c("base_case", paste0("intervention_", 1:5), "scenario_2", "scenario_3")

#' Position of each scenario's flag in the Int1..Int5, Scen1..Scen3 argument list.
TABBY2_FLAG_INDEX <- c(base_case = 0, intervention_1 = 1, intervention_2 = 2,
                       intervention_3 = 3, intervention_4 = 4, intervention_5 = 5,
                       scenario_1 = 6, scenario_2 = 7, scenario_3 = 8)

#' Run one scenario for every parameter set and return the parameter set x year x
#' output array. `builtin` selects one of the Int1..Scen3 flags by scenario id
#' (NULL or "base_case" sets none); prg_chng and ttt_list carry custom changes.
run_scenario <- function(loc, ParMatrix, builtin = NULL,
                         prg_chng = def_prgchng(ParMatrix[1, ]), ttt_list = def_ttt(),
                         n_cores = 1, endyr = 2050) {
  flags <- rep(0, 8)
  if (!is.null(builtin) && builtin != "base_case") {
    k <- TABBY2_FLAG_INDEX[builtin]
    if (is.na(k)) stop("Unknown built-in scenario: ", builtin)
    flags[k] <- 1
  }
  out <- OutputsInt(loc = loc, ParMatrix = ParMatrix, n_cores = n_cores,
                    startyr = 1950, endyr = endyr,
                    Int1 = flags[1], Int2 = flags[2], Int3 = flags[3], Int4 = flags[4],
                    Int5 = flags[5], Scen1 = flags[6], Scen2 = flags[7], Scen3 = flags[8],
                    prg_chng = prg_chng, ttt_list = ttt_list)
  # a single parameter set comes back as a matrix; keep the array shape throughout
  if (length(dim(out)) == 2) {
    out <- array(out, dim = c(1, dim(out)), dimnames = list(NULL, NULL, colnames(out)))
  }
  out
}

#' Run the built-in scenarios, `n_parallel` at a time. A precomputed base case
#' (the {loc}_results_1 array) is reused rather than rerun.
run_builtin_scenarios <- function(loc, ParMatrix, scenarios = TABBY2_SCENARIOS,
                                  base_case = NULL, n_parallel = 1) {
  to_run <- scenarios
  if (!is.null(base_case)) to_run <- setdiff(to_run, "base_case")
  runs <- mclapply(to_run, function(s) run_scenario(loc, ParMatrix, builtin = s),
                   mc.cores = max(1, min(n_parallel, length(to_run))))
  names(runs) <- to_run
  failed <- vapply(runs, function(r) inherits(r, "try-error") || !is.array(r), logical(1))
  if (any(failed)) stop("Scenario runs failed: ", paste(to_run[failed], collapse = ", "))
  raw <- list()
  for (s in scenarios) raw[[s]] <- if (s == "base_case" && !is.null(base_case)) base_case else runs[[s]]
  raw
}
