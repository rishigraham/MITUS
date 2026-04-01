#'loads data for the geography of interest
#'@name model_load
#'@param loc two letter postal abbreviation for states; US for national
#'@return void
model_load<-function(loc="US"){
#' add loc as a global variable
# loc<<-loc
  library(mnormt)
  library(parallel)
  library(lhs)
  library(Rcpp)
  library(MCMCpack)
  library(MASS)

  loc_info <- resolve_location(loc)
  if (is.null(loc_info) && loc != "US") {
    stop(paste0("Location '", loc, "' not found in stateID or county_ID.csv"))
  }
  # Register county locations in stateID so downstream st-based lookups work
  if (!is.null(loc_info) && loc_info$loc_type == "county") {
    register_location(loc_info)
  }

#'load necessary datasets
#'Model Input
if (loc=="US"){
  CalibDat<<-readRDS(find_loc_file("US", "CalibDat"))
  CalibDatCases<<-CalibDat
  ## ParamInit and StartVal were replaced on 8/1 to clean up, but 07-07 are last calibrated version
  ParamInit<<-as.data.frame(readRDS(find_loc_file("US", "ParamInit")))
  StartVal<<-readRDS(find_loc_file("US", "StartVal"))
  Inputs<<-readRDS(find_loc_file("US", "Inputs"))
  Opt <<- readRDS(find_loc_file("US", "Optim", required = FALSE))
  Par <<- readRDS(find_loc_file("US", "Param", required = FALSE))
} else {
  CalibDat<<-CalibDatState<<-readRDS(find_loc_file(loc, "CalibDat", fallback_loc="ST"))
  CalibDatCases<<-CalibDat
  ParamInit_st<<-ParamInit<<-readRDS(find_loc_file("ST", "ParamInit", fallback_loc=loc))
  StartVal_st<<-StartVal<<-readRDS(find_loc_file("ST", "StartVal", fallback_loc=loc))
  Inputs<<-readRDS(find_loc_file(loc, "ModelInputs"))
  # Load DeathByAge data once (used by calibration likelihood functions)
  dba_file <- find_loc_file(loc, "deathbyAge", fallback_loc = "ST", required = FALSE)
  if (!is.null(dba_file)) {
    DeathByAge <<- readRDS(dba_file)
  }
  par_file <- find_loc_file(loc, "Param", required = FALSE)
  if (!is.null(par_file)) {
    Par <<- readRDS(par_file)
  }
  #last input change was to update the RR active TB by age in immigrants
}
if (loc=="US"){
  LgtCurveY2 <- function(StYr,Endyr,EndVal) { z <- log(1/0.005-1)
  zz  <- seq(-z*(1+2*(StYr-1950)/(Endyr-StYr)),z*(1+2*(2019-Endyr)/(Endyr-StYr)),by=(2*z)/(Endyr-StYr))
  zz  <- as.numeric(EndVal)/(1+exp(-zz));    zz  }
  ImptWeights <- LgtCurveY2(2000,2019,0.95)+0.05
  names(ImptWeights) <- 1950:2019

  wts <<- ImptWeights
  P  <<- ParamInit[,1]
  names(P) <<- rownames(ParamInit)

  ii <<-  ParamInit[,5]==1
  ParamInitZ <<- ParamInit[ParamInit$Calib==1,]
  idZ0 <<- ParamInitZ[,4]==0
  idZ1 <<- ParamInitZ[,4]==1
  idZ2 <<- ParamInitZ[,4]==2
} else {
  LgtCurveY2 <- function(StYr,Endyr,EndVal) { z <- log(1/0.005-1)
  zz  <- seq(-z*(1+2*(StYr-1950)/(Endyr-StYr)),z*(1+2*(2019-Endyr)/(Endyr-StYr)),by=(2*z)/(Endyr-StYr))
  zz  <- as.numeric(EndVal)/(1+exp(-zz));    zz  }

  ImptWeights <- LgtCurveY2(2000,2019,0.95)+0.05
  names(ImptWeights) <- 1950:2019
  wts <<- ImptWeights
  W <- wts[year_to_idx(1993):year_to_idx(2018)];  W["2016"] <- 4
  wtZ <<-W

  #creation of background parameters
  #elements of P will be replaced from either the StartVals in the case
  #of optimization or the user inputted dataset

  P  <<- ParamInit_st[,1]
  names(P) <<- rownames(ParamInit_st)
  ii <<-  ParamInit_st[,5]==1
  ParamInitZ <<- ParamInit_st[ParamInit_st$Calib==1,]
  idZ0 <<- ParamInitZ[,4]==0
  idZ1 <<- ParamInitZ[,4]==1
  idZ2 <<- ParamInitZ[,4]==2
  ParamInit<<-ParamInit_st
}
  prgchng<<-def_prgchng(P)

return(invisible(NULL))
}

#'loads data for the geography of interest
#'@name model_load_demo
#'@param loc two letter postal abbreviation for states; US for national
#'@return void
model_load_demo<-function(loc="US"){
  library(mnormt)
  library(parallel)
  library(lhs)
  library(Rcpp)
  library(MCMCpack)
  library(MASS)

  loc_info <- resolve_location(loc)
  if (is.null(loc_info) && loc != "US") {
    stop(paste0("Location '", loc, "' not found in stateID or county_ID.csv"))
  }

  #'lazy load necessary datasets
  #'Model Input
  if (loc=="US"){
    CalibDat<<-readRDS(find_loc_file("US", "CalibDat"))
    ParamInit_demo<<-readRDS(find_loc_file("US", "ParamInitdemo"))
    StartVal_demo<<-readRDS(find_loc_file("US", "StartValdemo"))
    Inputs<<-readRDS(find_loc_file("US", "Inputs"))

  } else {
    CalibDat<<-CalibDatState<<-readRDS(find_loc_file(loc, "CalibDat", fallback_loc="ST"))
    ParamInit_st<<-ParamInit<<-readRDS(find_loc_file("ST", "ParamInit"))
    StartVal_st<<-StartVal<<-readRDS(find_loc_file("ST", "StartVal"))
    Inputs<<-readRDS(find_loc_file(loc, "ModelInputs"))
  }

  if (loc=="US")
  {
    wts <<- CalibDat[["ImptWeights"]]

  } else {
    wts <<- CalibDatState[["ImptWeights"]]
  }
  #'creation of background parameters
  #'elements of P will be replaced from either the StartVals in the case
  #'of optimization or the user inputted dataset
  #'
  P  <<- ParamInit_demo[,1]
  names(P) <<- rownames(ParamInit_demo)

  ii <<-  ParamInit_demo[,5]==1
  ParamInitZ <<- ParamInit_demo[ParamInit_demo$Calib==1,]
  idZ0 <<- ParamInitZ[,4]==0
  idZ1 <<- ParamInitZ[,4]==1
  idZ2 <<- ParamInitZ[,4]==2

  #'format the calibration data
  # targets[["ParamInitZ"]]<-ParamInitZ
  # targets[["idZ0"]]<-idZ0
  # targets[["idZ1"]]<-idZ1
  # targets[["idZ2"]]<-idZ2
  # targets[["P"]]<-P
  ## ALSO LOAD IN THE BASELINE PROGRAM CHANGE VECTOR
  prgchng<<-def_prgchng(P)
  return(invisible(NULL))
}
