#' Location and Year Index Utilities
#'
#' Helper functions for converting between calendar years and model indices,
#' and for discovering location-specific data files within the MITUS package.
#' The MITUS model uses 1950 as the base year (index 1).

#' Convert a calendar year to a model index
#'@name year_to_idx
#'@param year Calendar year (numeric)
#'@param base_year Base year of the model (default: 1950)
#'@return Integer index (1-based, where base_year = 1)
#'@export
year_to_idx<-function(year,base_year=1950){
  year-base_year+1
}

#' Convert a calendar year and month to a monthly model index
#'@name year_month_to_idx
#'@param year Calendar year (numeric)
#'@param month Month (1-12)
#'@param base_year Base year of the model (default: 1950)
#'@return Integer monthly index (1-based, where base_year January = 1)
#'@export
year_month_to_idx<-function(year,month,base_year=1950){
  (year-base_year)*12+month
}

#' Find a location-specific data file
#'
#' Searches for RDS files matching {loc}_{filetype}[_.]*.rds. When data_dir
#' is provided, searches there first (e.g., an external input directory).
#' Falls back to the installed package inst/ directory via system.file().
#' Supports dated filenames (ST_CalibDat_04-20-22.rds) and plain filenames
#' (SanDiego_CalibDat.rds). If multiple matches, returns most recently modified.
#'@name find_loc_file
#'@param loc Location code (e.g., "CA", "SanDiego", "US", "ST")
#'@param filetype File type pattern (e.g., "ModelInputs", "CalibDat")
#'@param fallback_loc Optional fallback location to try if no file found
#'@param data_dir Optional external directory to search before package inst/
#'@param required If TRUE (default), stop with error when no file found
#'@return Full file path to the matched RDS file
#'@export
find_loc_file<-function(loc,filetype,fallback_loc=NULL,data_dir=NULL,required=TRUE){
  pattern<-paste0("^",loc,"_?",filetype,"([_.].*)?\\.rds$")
  fb_pattern<-if(!is.null(fallback_loc)) paste0("^",fallback_loc,"_?",filetype,"([_.].*)?\\.rds$")

  # helper: search a directory for matching files
  search_dir<-function(base_dir,prefix,pat){
    d<-file.path(base_dir,prefix)
    if(!dir.exists(d)) return(NULL)
    files<-list.files(d,pattern=pat,full.names=TRUE,ignore.case=TRUE)
    if(length(files)==0) return(NULL)
    if(length(files)>1){
      warning(paste0("Multiple files matched for ",prefix,"/",filetype,
                     "; using most recently modified: ",basename(files[which.max(file.mtime(files))])))
    }
    files[which.max(file.mtime(files))]
  }

  # search external data_dir first (e.g., SIMULATOR_INPUT_DIR)
  if(!is.null(data_dir)){
    f<-search_dir(data_dir,loc,pattern)
    if(!is.null(f)) return(f)
    if(!is.null(fallback_loc)){
      f<-search_dir(data_dir,fallback_loc,fb_pattern)
      if(!is.null(f)) return(f)
    }
  }
  # fall back to installed package inst/ directory
  pkg_dir<-system.file(loc,package="MITUS")
  if(nchar(pkg_dir)>0){
    files<-list.files(pkg_dir,pattern=pattern,full.names=TRUE,ignore.case=TRUE)
    if(length(files)>0){
      if(length(files)>1){
        warning(paste0("Multiple files matched for ",loc,"/",filetype,
                       "; using most recently modified: ",basename(files[which.max(file.mtime(files))])))
      }
      return(files[which.max(file.mtime(files))])
    }
  }
  if(!is.null(fallback_loc)){
    fb_pkg_dir<-system.file(fallback_loc,package="MITUS")
    if(nchar(fb_pkg_dir)>0){
      files<-list.files(fb_pkg_dir,pattern=fb_pattern,full.names=TRUE,ignore.case=TRUE)
      if(length(files)>0){
        if(length(files)>1){
          warning(paste0("Multiple files matched for ",fallback_loc,"/",filetype,
                         "; using most recently modified: ",basename(files[which.max(file.mtime(files))])))
        }
        return(files[which.max(file.mtime(files))])
      }
    }
  }
  if(required){
    searched<-loc
    if(!is.null(fallback_loc)) searched<-paste0(loc," or ",fallback_loc)
    stop(paste0("Required data file not found: ",filetype," for location '",
                searched,"'"))
  }
  return(NULL)
}

#' Resolve a location code to its lookup entry
#'
#' Checks stateID first, then falls back to county_ID.csv for sub-state
#' locations. Pure lookup — does not modify global state. Use
#' \code{register_location} to add county entries to stateID at runtime.
#'@name resolve_location
#'@param loc Location code (USPS abbreviation, county code, or FIPS)
#'@return List with name, fips, code, loc_type, st. NULL if not found.
#'@export
resolve_location<-function(loc){
  data("stateID",package="MITUS")
  StateID<-as.data.frame(stateID)
  if(loc=="US"){
    return(list(name="United States",fips="0",code="US",
                loc_type="national",st=0))
  }
  # check stateID by USPS code (includes any previously registered counties)
  idx<-which(StateID$USPS==loc)
  if(length(idx)>0){
    return(list(name=StateID$Name[idx],fips=StateID$FIPS[idx],
                code=loc,loc_type="state",st=idx))
  }
  # check county_ID.csv
  county_file<-system.file("extdata","county_ID.csv",package="MITUS")
  if(nchar(county_file)>0 && file.exists(county_file)){
    countyID<-read.csv(county_file,stringsAsFactors=FALSE)
    cidx<-which(countyID$Code==loc | as.character(countyID$FIPS)==loc)
    if(length(cidx)>0){
      county<-countyID[cidx[1],]
      return(list(name=county$Name,fips=county$FIPS,
                  code=county$Code,loc_type="county",st=NA,
                  state_fips=county$StateFIPS))
    }
  }
  return(NULL)
}

#' Register a location in the global stateID table
#'
#' Adds a location to stateID so downstream \code{which(StateID$USPS==loc)}
#' lookups work. Call once during model_load, not during calibration iterations.
#' Idempotent — does nothing if already registered.
#'@name register_location
#'@param loc_info List as returned by \code{resolve_location}
#'@return Row index (st) of the location in stateID (invisible)
#'@export
register_location<-function(loc_info){
  data("stateID",package="MITUS")
  StateID<-as.data.frame(stateID)
  idx<-which(StateID$USPS==loc_info$code)
  if(length(idx)>0) return(invisible(idx))
  new_row<-data.frame(Name=loc_info$name,FIPS=as.character(loc_info$fips),
                      USPS=loc_info$code,stringsAsFactors=FALSE)
  StateID<-rbind(StateID,new_row)
  stateID<<-as.matrix(StateID)
  return(invisible(nrow(StateID)))
}

#' Merge a flat county-level CalibDat overlay into a country-wide CalibDat base
#'
#' MITUS's calibration code expects \code{CalibDat} to be the country-wide
#' state-indexed structure where many fields are length-51 lists indexed by
#' position in \code{stateID}. A sub-state location (e.g. a county) provides
#' its data as a flat single-location overlay; this function inserts that
#' overlay into the country-wide structure at row \code{st} so the existing
#' \code{[[st]]} dereference in likelihood code resolves correctly.
#'
#' Slots in \code{state_indexed_fields} are populated from the overlay when
#' present and left as the base value (typically NULL beyond length 51) when
#' absent — i.e. no silent fallback to another state's data. Non-state-indexed
#' fields are replaced with the overlay's value when the overlay provides them.
#'
#' Issues a warning if the country-wide CalibDat contains length-51 list
#' fields not in \code{state_indexed_fields}, indicating possible drift if
#' MITUS adds new state-indexed fields after this generalization. The hardcoded
#' field list should be updated in that case.
#'
#'@name merge_county_calibdat
#'@param base Country-wide CalibDat (e.g., contents of \code{ST_CalibDat_*.rds})
#'@param overlay Flat single-location CalibDat with canonical field names
#'@param st Integer row index assigned to the location in \code{stateID}
#'@return Merged CalibDat suitable for use as both \code{CalibDat} and \code{CalibDatState}
#'@export
merge_county_calibdat<-function(base,overlay,st){
  state_indexed_fields<-c(
    "cases_yr_st","cases_yr_ag_nat_st","cases_yr_ag_nat_st_5yr",
    "hr_cases","homeless_pop","TLTBI_volume_state",
    "tbdeaths","pop_00_17","pop_00_19","pop_50_10"
  )
  out<-base
  # 1. Slot county data into known state-indexed fields
  for(fld in state_indexed_fields){
    if(!is.null(overlay[[fld]])&&!is.null(out[[fld]])){
      out[[fld]][[st]]<-overlay[[fld]]
    }
  }
  # 2. Replace non-state-indexed fields with overlay values
  for(fld in names(overlay)){
    if(!(fld %in% state_indexed_fields)){
      out[[fld]]<-overlay[[fld]]
    }
  }
  # 3. Drift check: warn about length-51 list fields not in our merge list
  is_state_list<-function(x) is.list(x)&&!is.data.frame(x)&&length(x)==51
  base_state_lists<-vapply(base,is_state_list,logical(1))
  unmerged<-setdiff(names(base)[base_state_lists],state_indexed_fields)
  if(length(unmerged)>0){
    warning(paste0(
      "merge_county_calibdat: country-wide CalibDat contains length-51 list ",
      "fields not in state_indexed_fields: ",paste(unmerged,collapse=", "),
      ". County data is missing for these fields. Update state_indexed_fields ",
      "in R/location_utils.R."
    ))
  }
  out
}
