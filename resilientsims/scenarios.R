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

PRG_CHNG_FIELDS <- c("start_yr", "scrn_cov", "IGRA_frc", "ltbi_init_frc",
                     "frc_3hp", "comp_3hp", "frc_4r", "comp_4r", "frc_3hr", "comp_3hr",
                     "tb_tim2tx_frc", "tb_txdef_frc")
TTT_FIELDS <- c("NativityGrp", "AgeGrp", "NRiskGrp", "FrcScrn", "StartYr", "EndYr",
                "RRprg", "RRmu", "RRPrev")
TTT_NATIVITY <- c("All", "USB", "NUSB")
TTT_AGE <- c("All", "0 to 24", "25 to 64", "65+")

scenario_slug <- function(name) {
  s <- gsub("^_+|_+$", "", tolower(gsub("[^A-Za-z0-9]+", "_", name)))
  if (nzchar(s)) s else "scenario"
}

#' Validate a custom scenario request and return the full prg_chng vector and
#' ttt_list with defaults filled in. Units follow param_init(): fractions as
#' fractions, tb_tim2tx_frc as a percent of the current duration of
#' infectiousness, NRiskGrp in millions of people.
scenario_spec <- function(params, ParMatrix) {
  name <- params$name
  if (is.null(name) || !nzchar(trimws(as.character(name)))) {
    stop("Parameter 'name' (scenario name) is required")
  }
  n_sets <- if (!is.null(params$n_param_sets)) as.integer(params$n_param_sets) else 2L
  if (is.na(n_sets) || n_sets < 1) stop("Parameter 'n_param_sets' must be a positive integer")

  num <- function(x, field) {
    v <- suppressWarnings(as.numeric(x))
    if (length(v) != 1 || is.na(v)) stop("Field '", field, "' must be a single number")
    v
  }
  in_range <- function(v, lo, hi, field) {
    if (v < lo || v > hi) {
      stop(sprintf("Field '%s' = %s is outside [%s, %s]", field, format(v), format(lo), format(hi)))
    }
  }

  prg <- def_prgchng(ParMatrix[1, ])
  given <- params$prg_chng
  if (length(given) > 0) {
    bad <- setdiff(names(given), PRG_CHNG_FIELDS)
    if (length(bad) > 0) stop("Unknown prg_chng field(s): ", paste(bad, collapse = ", "))
    for (f in names(given)) prg[f] <- num(given[[f]], f)
    in_range(prg["start_yr"], 2022, 2050, "start_yr")
    if (prg["scrn_cov"] <= 0) stop("Field 'scrn_cov' (multiple of current screening) must be positive")
    for (f in c("IGRA_frc", "ltbi_init_frc", "frc_3hp", "comp_3hp", "frc_4r", "comp_4r",
                "frc_3hr", "comp_3hr", "tb_txdef_frc")) in_range(prg[f], 0, 1, f)
    if (prg["tb_tim2tx_frc"] <= 0 || prg["tb_tim2tx_frc"] > 100) {
      stop("Field 'tb_tim2tx_frc' (percent of current duration) must be in (0, 100]")
    }
    regimens <- sum(prg[c("frc_3hp", "frc_4r", "frc_3hr")])
    if (abs(regimens - 1) > 1e-6) {
      stop(sprintf("Regimen fractions frc_3hp + frc_4r + frc_3hr must sum to 1 (got %.4f)", regimens))
    }
  }

  ttt <- def_ttt()
  given <- params$ttt_list
  if (length(given) > 0) {
    bad <- setdiff(names(given), TTT_FIELDS)
    if (length(bad) > 0) stop("Unknown ttt_list field(s): ", paste(bad, collapse = ", "))
    for (f in intersect(names(given), c("NativityGrp", "AgeGrp"))) ttt[[f]] <- as.character(given[[f]])
    for (f in setdiff(intersect(names(given), TTT_FIELDS), c("NativityGrp", "AgeGrp"))) {
      ttt[[f]] <- num(given[[f]], f)
    }
    if (!ttt$NativityGrp %in% TTT_NATIVITY) {
      stop("Field 'NativityGrp' must be one of ", paste(TTT_NATIVITY, collapse = ", "))
    }
    if (!ttt$AgeGrp %in% TTT_AGE) stop("Field 'AgeGrp' must be one of ", paste(TTT_AGE, collapse = ", "))
    if (ttt$NRiskGrp < 0) stop("Field 'NRiskGrp' (millions of people) must be non-negative")
    in_range(ttt$FrcScrn, 0, 1, "FrcScrn")
    in_range(ttt$StartYr, 2022, 2050, "StartYr")
    in_range(ttt$EndYr, 2022, 2050, "EndYr")
    if (ttt$StartYr > ttt$EndYr) stop("TTT StartYr must not be after EndYr")
    for (f in c("RRprg", "RRmu", "RRPrev")) {
      if (ttt[[f]] <= 0) stop("Field '", f, "' (rate ratio) must be positive")
    }
  }

  list(name = as.character(name), id = paste0("custom_", scenario_slug(name)),
       prg_chng = prg, ttt_list = ttt,
       ttt_active = ttt$NRiskGrp != 0 && ttt$FrcScrn != 0,
       n_param_sets = n_sets,
       inputs = list(prg_chng = params$prg_chng, ttt_list = params$ttt_list))
}

#' Run a validated custom scenario on the first n_param_sets parameter sets.
run_custom_scenario <- function(loc, ParMatrix, spec, n_cores = 1) {
  if (spec$ttt_active) {
    rds <- find_loc_file(loc, "results_1", data_dir = get_data_dir(), required = FALSE)
    rda <- system.file(paste0(loc, "/", loc, "_results_1.rda"), package = "MITUS")
    if (is.null(rds) && !nzchar(rda)) {
      stop("A targeted testing scenario needs the location's base-case results file (",
           loc, "_results_1.rds) in the configuration")
    }
  }
  PM <- ParMatrix[seq_len(min(spec$n_param_sets, nrow(ParMatrix))), , drop = FALSE]
  if (nrow(PM) == 1) PM <- PM[c(1, 1), , drop = FALSE]
  run_scenario(loc, PM, prg_chng = spec$prg_chng, ttt_list = spec$ttt_list, n_cores = n_cores)
}
