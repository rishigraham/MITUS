#' ResilientHub2 UI data export for Tabby2-style scenario results
#'
#' Sourced by wrapper.R from this directory. Turns scenario output arrays
#' (parameter set x year x output, see scenarios.R) into one
#' {LOC}_tabby2_ui_data.json per location in the v2 schema the hub explorer reads:
#' `meta` (scenarios, measures, groups, years, location) plus a flat `series`
#' list, each series holding one value per year 2022-2050 for one
#' scenario x measure x nativity x age group. Values are derived from the
#' 858-column model output by name and averaged over the parameter sets.

library(jsonlite)

TABBY2_YEARS_ALL <- 1950:2050
TABBY2_YEARS_OUT <- 2022:2050
TABBY2_SUMMARY_YEARS <- c(2022, 2025, 2030, 2040, 2050)

TABBY2_BANDS <- c("0_4", "5_14", "15_24", "25_34", "35_44", "45_54", "55_64",
                  "65_74", "75_84", "85_94", "95p")
TABBY2_AGE_GROUPS <- list(all = TABBY2_BANDS, `0-24` = TABBY2_BANDS[1:3],
                          `25-64` = TABBY2_BANDS[4:7], `65+` = TABBY2_BANDS[8:11])

TABBY2_SCENARIO_LABEL <- c(
  base_case      = "Base Case",
  intervention_1 = "LTBI treatment for new migrants",
  intervention_2 = "Improved LTBI treatment in the US",
  intervention_3 = "Enhanced case detection",
  intervention_4 = "Enhanced TB treatment",
  intervention_5 = "All improvements",
  scenario_2     = "No TB or LTBI in New Migrants After 2022",
  # Labelled from the Scen3 implementation in param_init() (TB and LTBI among new
  # arrivals halve every decade from 2022), which is what the Tabby2 app config
  # calls "Improving Global TB Control"; the config's positional label for the
  # third analysis ("Reduced TB Prevention Effort") does not describe this data.
  scenario_3     = "Improving Global TB Control")

TABBY2_SCENARIO_DESC <- c(
  base_case      = "Continuation of current TB policy and services.",
  intervention_1 = "Provision of LTBI testing and treatment for all new legal immigrants entering the US",
  intervention_2 = "Intensification of the current LTBI targeted testing and treatment policy for high-risk populations, doubling treatment uptake within each risk group compared to current levels, and increasing the fraction cured among individuals initiating LTBI treatment, via a 3-month Isoniazid-Rifapentine drug regimen",
  intervention_3 = "Improved detection of TB cases, such that the duration of untreated disease (time from TB incidence to the initiation of treatment) is reduced by 50%",
  intervention_4 = "Improved treatment quality for TB, such that treatment default, failure rates, and the fraction of individuals receiving an incorrect drug regimen are reduced by 50% from current levels",
  intervention_5 = "The combination of all intervention scenarios described above",
  scenario_2     = "From 2022 onwards, all individuals immigrating to the United States are free of M. TB infection, while maintaining the same total volume of immigration.",
  scenario_3     = "From 2022 onwards, the prevalence of TB disease and latent TB infection among new arrivals to the United States declines gradually, halving every 10 years, while maintaining the same total volume of immigration.")

TABBY2_SCENARIO_GROUP <- c(
  base_case = "interventions", intervention_1 = "interventions",
  intervention_2 = "interventions", intervention_3 = "interventions",
  intervention_4 = "interventions", intervention_5 = "interventions",
  scenario_2 = "analyses", scenario_3 = "analyses")

#' TB treatment completions are not exported: cSim fills them as diagnoses x the
#' monthly completion rate (TxVec[0], about 1/9) rather than the completion
#' probability, so the column is a constant multiple of initiations.
TABBY2_MEASURES <- list(
  population            = list(label = "Population", unit = "persons", group = "core"),
  tb_cases              = list(label = "TB Cases", unit = "cases/year", group = "core"),
  tb_cases_per100k      = list(label = "TB Incidence", unit = "cases per 100,000/year", group = "core"),
  tb_deaths             = list(label = "TB-Related Deaths", unit = "deaths/year", group = "core"),
  tb_deaths_per100k     = list(label = "TB Mortality", unit = "deaths per 100,000/year", group = "core"),
  ltbi_count            = list(label = "LTBI Prevalence", unit = "persons", group = "core"),
  ltbi_pct              = list(label = "LTBI Prevalence", unit = "percent", group = "core"),
  tb_infections         = list(label = "Incident M. tb Infections", unit = "infections/year", group = "core",
    description = "Incident M. TB infections each year due to transmission within the US (includes reinfection of individuals with prior infection, excludes migrants entering the US with established LTBI)"),
  tb_infections_per100k = list(label = "Incident M. tb Infections", unit = "infections per 100,000/year", group = "core"),
  ltbi_tests            = list(label = "LTBI Tests", unit = "tests/year", group = "services",
    description = "LTBI tests administered within the model's screening cascade. Scale is set by each location's calibrated LTBI-treatment-volume target, not by total real-world testing volume; per-capita levels are therefore not comparable across locations without noting their targets."),
  ltbi_tx_inits         = list(label = "LTBI Treatment Initiations", unit = "initiations/year", group = "services",
    description = "LTBI treatment initiations in the model's screening cascade; calibrated to each location's LTBI-treatment-volume target."),
  ltbi_tx_comps         = list(label = "LTBI Treatment Completions", unit = "completions/year", group = "services",
    description = "LTBI treatment initiations times the modeled completion probability (about 78% at current parameters)."),
  tb_tx_inits           = list(label = "TB Disease Treatment Initiations", unit = "initiations/year", group = "services",
    description = "Equals TB diagnoses: in the model, every diagnosed case initiates treatment."))

TABBY2_NOTES <- c(
  "scenario_1 ('No Transmission Within the US After 2022') is excluded: its implementation is commented out in current MITUS, so its data equals the base case.",
  "scenario_3 is labeled 'Improving Global TB Control' from its implementation (gradual halving of TB/LTBI among new arrivals per decade from 2022); the tabby2 app config calls its third analysis 'Reduced TB Prevention Effort', which does not describe this data.",
  "tb_tx_comps (TB treatment completions) is excluded: the model fills it as diagnoses x a monthly completion rate (~0.111), not a completion probability; the same defect is present in PPML's shipped US data.",
  "Comparison modes are computed in the front end from these absolute values; both percentage modes are ratios (base case = 100).")

#' Display label for a location code, from the state or county registry.
location_label <- function(loc) {
  # counties carry their state as a suffix; checked first because model_load()
  # registers the county in the state table, where it then looks like a state
  county_file <- system.file("extdata", "county_ID.csv", package = "MITUS")
  if (nzchar(county_file)) {
    counties <- read.csv(county_file, stringsAsFactors = FALSE)
    k <- which(counties$Code == loc)
    if (length(k) > 0) {
      stateID <- as.data.frame(get_stateID(), stringsAsFactors = FALSE)
      idx <- which(as.integer(as.character(stateID$FIPS)) ==
                     as.integer(counties$StateFIPS[k[1]]))
      suffix <- if (length(idx) > 0) paste0(", ", as.character(stateID$USPS[idx[1]])) else ""
      return(paste0(counties$Name[k[1]], suffix))
    }
  }
  info <- resolve_location(loc)
  if (is.null(info)) loc else info$name
}

#' Measure values by age group and nativity for one parameter set's
#' year x output matrix. Population columns are in millions.
measures_for <- function(m) {
  colnames(m) <- func_ResNam()
  g <- function(cols) if (length(cols) == 1) m[, cols] else rowSums(m[, cols, drop = FALSE])
  o <- list()
  for (ag in names(TABBY2_AGE_GROUPS)) {
    b <- TABBY2_AGE_GROUPS[[ag]]
    pop_all <- g(paste0("N_", b)); pop_us <- g(paste0("N_US_", b))
    cas_all <- g(paste0("NOTIF_", b)) + g(paste0("NOTIF_MORT_", b))
    cas_us  <- g(paste0("NOTIF_US_", b)) + g(paste0("NOTIF_US_MORT_", b))
    dth_us  <- g(paste0("TBMORT_US_", b)); dth_nus <- g(paste0("TBMORT_NUS_", b))
    lt_us   <- g(paste0("N_US_LTBI_", b)); lt_nus  <- g(paste0("N_FB_LTBI_", b))
    inf_us  <- g(paste0("N_newinf_USB_", b)); inf_nus <- g(paste0("N_newinf_NUSB_", b))
    svc <- function(stem) {
      usb <- g(paste0("N_", stem, "_USB_", b)); nusb <- g(paste0("N_", stem, "_NUSB_", b))
      list(all = usb + nusb, usb = usb, nusb = nusb)
    }
    tests <- svc("LtbiTests"); linits <- svc("LtbiTxInits"); lcomps <- svc("LtbiTxComps")
    tinits <- svc("TBTxInits")

    pop <- list(all = pop_all, usb = pop_us, nusb = pop_all - pop_us)
    cas <- list(all = cas_all, usb = cas_us, nusb = cas_all - cas_us)
    dth <- list(all = dth_us + dth_nus, usb = dth_us, nusb = dth_nus)
    ltb <- list(all = lt_us + lt_nus, usb = lt_us, nusb = lt_nus)
    inf <- list(all = inf_us + inf_nus, usb = inf_us, nusb = inf_nus)

    for (nat in c("all", "usb", "nusb")) {
      p <- pop[[nat]] * 1e6
      key <- function(measure) paste(ag, nat, measure, sep = "|")
      o[[key("population")]]            <- p
      o[[key("tb_cases")]]              <- cas[[nat]] * 1e6
      o[[key("tb_cases_per100k")]]      <- cas[[nat]] * 1e6 / p * 1e5
      o[[key("tb_deaths")]]             <- dth[[nat]] * 1e6
      o[[key("tb_deaths_per100k")]]     <- dth[[nat]] * 1e6 / p * 1e5
      o[[key("ltbi_count")]]            <- ltb[[nat]] * 1e6
      o[[key("ltbi_pct")]]              <- ltb[[nat]] / pop[[nat]] * 100
      o[[key("tb_infections")]]         <- inf[[nat]] * 1e6
      o[[key("tb_infections_per100k")]] <- inf[[nat]] * 1e6 / p * 1e5
      o[[key("ltbi_tests")]]            <- tests[[nat]] * 1e6
      o[[key("ltbi_tx_inits")]]         <- linits[[nat]] * 1e6
      o[[key("ltbi_tx_comps")]]         <- lcomps[[nat]] * 1e6
      o[[key("tb_tx_inits")]]           <- tinits[[nat]] * 1e6
    }
  }
  o
}

#' Series for one scenario array (parameter set x year x output), averaged over
#' the parameter sets and restricted to the export years.
scenario_series <- function(loc, loc_label, scenario_id, arr) {
  n_sets <- dim(arr)[1]
  keep <- TABBY2_YEARS_ALL %in% TABBY2_YEARS_OUT
  per_set <- lapply(seq_len(n_sets), function(i) measures_for(arr[i, , ]))
  out <- list()
  for (key in names(per_set[[1]])) {
    v <- Reduce(`+`, lapply(per_set, function(x) x[[key]])) / n_sets
    p <- strsplit(key, "|", fixed = TRUE)[[1]]
    out[[length(out) + 1]] <- list(location = loc, location_label = loc_label,
      scenario = scenario_id, measure = p[3], nativity = p[2], age_group = p[1],
      values = round(unname(v[keep]), 6))
  }
  out
}

#' Scenario metadata entries for the ids present in `raw`.
scenario_meta <- function(ids) {
  setNames(lapply(ids, function(s) list(
    id = s, label = unname(TABBY2_SCENARIO_LABEL[s]), group = unname(TABBY2_SCENARIO_GROUP[s]),
    description = unname(TABBY2_SCENARIO_DESC[s]))), ids)
}

#' Assemble the v2 payload for one location from a named list of scenario arrays.
tabby2_ui_payload <- function(loc, raw, basis, loc_label = location_label(loc)) {
  n_sets <- dim(raw[[1]])[1]
  series <- list()
  for (s in names(raw)) series <- c(series, scenario_series(loc, loc_label, s, raw[[s]]))
  list(
    meta = list(
      generated = format(Sys.Date(), "%Y-%m-%d"), model = "MITUS", schema_version = 2,
      years = TABBY2_YEARS_OUT, summary_years = TABBY2_SUMMARY_YEARS,
      locations = list(list(code = loc, label = loc_label, parameter_sets = n_sets,
                            basis = basis)),
      scenarios = scenario_meta(names(raw)),
      scenario_groups = list(interventions = "Modeled Scenarios", analyses = "Sensitivity Analyses"),
      measure_groups = list(core = "Modelled Outcomes", services = "Counts of Services"),
      nativity = list(all = "Total", usb = "US-born", nusb = "Non-US-born"),
      age_groups = names(TABBY2_AGE_GROUPS),
      measures = TABBY2_MEASURES,
      notes = TABBY2_NOTES),
    series = series)
}

write_tabby2_ui_json <- function(payload, file) {
  write_json(payload, file, auto_unbox = TRUE, digits = NA)
  invisible(file)
}
