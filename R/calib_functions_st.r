#'Calibration Functions for use with the tb_model.cpp model for state level
#'This script creates several individual log likelihood functions
#'for the calibration of the State Level TB model in tb_model.cpp
#'These llikelihood functions are called in IMIS_functions.R
#'takes in the outputs and calibration data and creates likelihood functions

# Calibration year index constants (year_to_idx: yr-1950+1)
idx_5yr_end<-c(49,54,59,64,70)   # 1998,2003,2008,2013,2019
idx_5yr_start<-c(45,50,55,60,65) # 1994,1999,2004,2009,2014
idx_decades<-c(11,21,31,41,51,61,70) # 1960,1970,1980,1990,2000,2010,2019

# Accessors for the calibration data, keyed on the names each table carries.

#' Columns of cases_nat_st_5yr holding the 5-year case counts
#' (e.g. "X1996.2000".."X2016.2020"), oldest bucket first.
#' @noRd
nat_5yr_bucket_cols <- function(d) {
  grep("^X[0-9]{4}\\.[0-9]{4}$", colnames(d), value = TRUE)
}

#' Age-group columns of cases_yr_ag_nat_st_5yr (e.g. "X0.4".."X85p"), youngest
#' age group first.
#' @noRd
ag_5yr_cols <- function(d) {
  grep("^X[0-9]", colnames(d), value = TRUE)
}

#' Rows of DeathByAge[[st]] for the given calendar years, NA where a year is
#' absent. DeathByAge is an unnamed matrix with the year in column 1.
#' @noRd
dba_year_rows <- function(m, years) {
  match(years, m[, 1])
}

#' Column of DeathByAge[[st]] holding deaths totalled across age groups.
#' @noRd
dba_total_col <- function(m) ncol(m)

#'Total Diagnosed Cases 1953-2016
#'Motivation: Normal, mean centered with CI = +/- 5% of the mean
#'@name notif_tot_lLik_st
#'@param V vector of total notifi
#'cations 1953-2014
#'@return likelihood
notif_tot_lLik_st <- function(V,st) { # V = simulator NOTIF_ALL+NOTIF_MORT_ALL, all years
  raw       <- CalibDatCases[["cases_yr_st"]][[st]]
  years_all <- as.integer(raw[, "year"])
  keep      <- years_all %in% as.integer(names(wts))
  years     <- years_all[keep]
  notif_tot <- as.numeric(raw[keep, "cases"])
  wt_slice  <- wts[as.character(years)]
  ### the most recent year of data carries double weight
  wt_slice[length(wt_slice)] <- wt_slice[length(wt_slice)]*2
  V_slice   <- V[year_to_idx(years)]
  ok        <- !is.na(notif_tot) & notif_tot != 0
  adj_1 <- sum((dnorm(notif_tot, notif_tot, notif_tot*0.1/1.96, log=T) * wt_slice)[ok])
  #notif tot is in real scale must scale outputs up
  sum((dnorm(notif_tot, V_slice*1e6, notif_tot*0.1/1.96, log=T) * wt_slice)[ok]) - adj_1
}

### ### ### TOTAL DIAGNOSED CASES 1953-1993  ### ### ### ### ### ### D
notif_decline_lLik_st <- function(V, st=st) {
  notif_tot     <- CalibDatCases[["cases_yr_st"]][[st]][-28,2];
  # V = vector of total notifications 1953-1993
  notif_decline  <- CalibDatCases[["cases_prop_change_53_94"]]
  notif_tot2     <- cumprod(notif_decline)/prod(notif_decline)*notif_tot[1]
  adj_1b         <- sum(dnorm(notif_tot2,notif_tot2,notif_tot2*0.2/1.96,log=T)*wts[year_to_idx(1953):year_to_idx(1993)])
  #notif tot is in real scale must scale outputs up
  sum(dnorm(notif_tot2,V*1e6,notif_tot2*0.2/1.96,log=T)*wts[year_to_idx(1953):year_to_idx(1993)]) - adj_1b
}

### ### ### CASES FB DISTRIBUTION 1993-2018  ### ### ### ### ### ###  D
# Motivatioµn: dirichlet-multinomial, multinomial data with additional non-sampling biases
notif_fb_lLik_st <- function(V,st,rho=0.005) { # V = table of notifications by fb 1993-2016 (row=24 years, col=fb then us)
  notif_age_fb0     <- CalibDatCases[["cases_yr_ag_nat_st"]][[st]][,,"nusb"]
  notif_age_us0     <- CalibDatCases[["cases_yr_ag_nat_st"]][[st]][,,"usb"]
  notif_fb      <- cbind(notif_age_fb0[,"total"],notif_age_us0[,"total"])
  adj_3         <- sum(dDirMult(M=notif_fb+0.01,n=notif_fb,Rho=rho)*wts[year_to_idx(1993):year_to_idx(2018)])
  #scale does not matter for dirichlet llikelihood
  (sum(dDirMult(M=V,n=notif_fb,Rho=rho)*wts[year_to_idx(1993):year_to_idx(2018)]) - adj_3)*2
}
###############################################################################################

notif_fb_5yr_lLik_st <- function(V,st,rho=0.0005) { # V = table of notifications by fb 1993-2016 (row=24 years, col=fb then us)
  ### Read in the nativity stratified data
  cn              <- CalibDatCases$cases_nat_st_5yr
  bkt             <- nat_5yr_bucket_cols(cn)
  notif_age_fb0   <- cn[cn$State.Code==st & cn$usb==0,bkt]
  notif_age_us0   <- cn[cn$State.Code==st & cn$usb==1,bkt]
  notif_fb      <- cbind(t(notif_age_fb0),t(notif_age_us0))
  adj_3         <- sum(dDirMult(M=notif_fb,n=notif_fb,Rho=rho)*wts[idx_5yr_end])
  V2<-matrix(0,5,2)
  V2[1,]<-colSums(V[1:5,]);V2[2,]<-colSums(V[6:10,]);V2[3,]<-colSums(V[11:15,]); V2[4,]<-colSums(V[16:20,]); V2[5,]<-colSums(V[21:25,])
  #scale does not matter for dirichlet llikelihood
  (sum(dDirMult(M=V2,n=notif_fb,Rho=rho)*wts[idx_5yr_end]) - adj_3)*2
}

#'FB Diagnosed Cases 1953-2016
#'Motivation: Normal, mean centered with CI = +/- 5% of the mean
#'@param V vector of total notifications 1953-2014
#'@return likelihood
notif_fb_5yr_lik <- function(V,st=st) {
  cn         <- CalibDatCases$cases_nat_st_5yr
  notif_fb   <- as.numeric(unlist(cn[cn$State.Code==st & cn$usb==0,nat_5yr_bucket_cols(cn)]))
  adj_1         <- sum(dnorm(notif_fb,notif_fb,notif_fb*0.05/1.96,log=T)*wts[idx_5yr_end])
  V2<-rep(0,5)
  V2[1]<-sum(V[1:5]);V2[2]<-sum(V[6:10]);V2[3]<-sum(V[11:15]); V2[4]<-sum(V[16:20]); V2[5]<-sum(V[21:25])
  (sum(dnorm(notif_fb,V2*1e6,notif_fb*0.05/1.96,log=T)*wts[idx_5yr_end]) - adj_1)
}

#'US Diagnosed Cases 1953-2016
#'Motivation: Normal, mean centered with CI = +/- 5% of the mean
#'@param V vector of total notifications 1953-2014
#'@return likelihood

notif_us_5yr_lik <- function(V,st=st) {
  cn         <- CalibDatCases$cases_nat_st_5yr
  notif_us   <- as.numeric(unlist(cn[cn$State.Code==st & cn$usb==1,nat_5yr_bucket_cols(cn)]))
  adj_1         <- sum(dnorm(notif_us,notif_us,notif_us*0.05/1.96,log=T)*wts[idx_5yr_end])
  V2<-rep(0,5)
  V2[1]<-sum(V[1:5]);V2[2]<-sum(V[6:10]);V2[3]<-sum(V[11:15]); V2[4]<-sum(V[16:20]); V2[5]<-sum(V[21:25])
  (sum(dnorm(notif_us,V2*1e6,notif_us*0.05/1.96,log=T)*wts[idx_5yr_end]) - adj_1)
}


### ### ### CASES FB DISTRIBUTION SLOPES OVER PAST 5 year  ### ### ### ### ### ###
notif_fbus_slp_lLik_st <- function(V,st) {
  # notif_fbus_slp5     <- as.numeric(CalibDatCases[["case_change_5"]][st,2:3])
  cn              <- CalibDatCases$cases_nat_st_5yr
  ### slope is taken across the two most recent 5-year buckets
  last2           <- tail(nat_5yr_bucket_cols(cn), 2)
  notif_age_fb0   <- cn[cn$State.Code==st & cn$usb==0,last2]
  notif_age_us0   <- cn[cn$State.Code==st & cn$usb==1,last2]
  tot_case_nat<-cbind(t(notif_age_fb0),t(notif_age_us0))
  #calculate the slopes
  notif_fbus_slp5<-apply(log(tot_case_nat),2,function(x) lm(x~I(1:2))$coef[2])
  # if (st ==35){
  #   adj_3a              <- dnorm(notif_fbus_slp5[1],notif_fbus_slp5[1],0.005,log=T)# V = table of notifications by fb 2011-2016 (row=6 years, col=fb then us)
  # } else{
    adj_3a              <- sum(dnorm(notif_fbus_slp5,notif_fbus_slp5,0.005,log=T))# V = table of notifications by fb 2011-2016 (row=6 years, col=fb then us)
  # }
  V1a <- matrix(c(sum(V[1:5,1]),sum(V[1:5,2]), sum(V[6:10,1]), sum(V[6:10,2])),2,2, byrow = TRUE)
  V2 <- apply(log(V1a),2,function(x) lm(x~I(1:2))$coef[2])
  if (st==35){
    # return(dnorm(notif_fbus_slp5[1],V2[1],0.005,log=T) - adj_3a)
    notif_fbus_slp5[2] <- -0.1
  }
    return(sum(dnorm(notif_fbus_slp5,V2,0.005,log=T)) - adj_3a)
}

### ### ### CASES HR DISTRIBUTION 1993-2014  ### ### ### ### ### ### D
# Motivation: dirichlet-multinomial, multinomial data with additional non-sampling biases

# notif_hr_lLik_st <- function(V,st,rho=0.005) { # V = table of notifications by tx history (row=97:16, col=n then e)
#   notif_hr0     <- CalibDatCases[["hr_cases"]][[st]]
#   notif_hr      <- cbind(notif_hr0[,1],1-notif_hr0[,1])#*notif_us_hr0[,2]
#   adj_5b           <- sum(dDirMult(M=notif_hr+0.01,n=notif_hr,Rho=rho)*wts[c(year_to_idx(1994), year_to_idx(1999), year_to_idx(2004), year_to_idx(2009), year_to_idx(2014))])
#   V2 <- rbind(colSums(V[1:5,]),colSums(V[6:10,]),colSums(V[11:15,]),colSums(V[16:20,]), colSums(V[21:25,]))
#   #scale does not matter for dirichlet llikelihood
#   sum(dDirMult(M=V2,n=notif_hr,Rho=rho)*wts[c(year_to_idx(1994), year_to_idx(1999), year_to_idx(2004), year_to_idx(2009), year_to_idx(2014))]) - adj_5b
# }

###############################################################################################
### ### ### US CASES AGE DISTRIBUTION 5yrs 1993-2016  ### ### ### ### ### ### D
# Motivation: dirichlet-multinomial, multinomial data with additional non-sampling biases
notif_age_us_5yr_lLik_st <- function(V,st,rho=0.1) { # V = table of us notifications by age 1993-2016 (row=24 years, col=11 ages)
  ### Read in the age and nativity stratified data
  ag5                  <- CalibDatCases$cases_yr_ag_nat_st_5yr[[st]]
  notif_age_us_5yr     <- as.data.frame(ag5[ag5[,"usb"]==1,ag_5yr_cols(ag5)])
  ### Format the model estimates to match the calibration data
  V2 <- V[,-11]; V2[,10] <- V2[,10]+V[,11]
  V3<-matrix(0,5,10)
  V3[1,]<-colSums(V2[1:5,]);V3[2,]<-colSums(V2[6:10,]);V3[3,]<-colSums(V2[11:15,]); V3[4,]<-colSums(V2[16:20,]); V3[5,]<-colSums(V2[21:25,])
  ### Check for missing calibration data
  if (sum(is.na(notif_age_us_5yr > 0))) {
    ### Initialize the likelihood value to zero
    tot_lik <- 0; adj_2a <-0
    ### If data is missing we need to create a single new bucket for these data and estimates
    ### Read in the total US cases to use as a total
    cn              <- CalibDatCases$cases_nat_st_5yr
    notif_age_us0   <- cn[cn$State.Code==st & cn$usb==1,nat_5yr_bucket_cols(cn)]
    for (i in 1:nrow(notif_age_us_5yr)){
      print(i)
      ### Check that this row has an NA
      if (sum(is.na(notif_age_us_5yr[i,])) > 0){
        ### Calibration data
        ### Which of the age groups are NA
        index <- which(is.na(notif_age_us_5yr[i,]))
        ### Remove these age groups from the dataframe
        known_cases <- notif_age_us_5yr[i,-index]
        ### Set the total of missing cases equal to the difference of total cases and known cases
        missing_cases <- notif_age_us0[i] - sum(known_cases)
        ### Add missing cases to the dataframe
        tot_cases <- c(unlist(known_cases), unlist(missing_cases))
        ### Model estimates (match to estimates)
        ### Create a sum of the estimates that correspond to missing age groups
        known_est   <- V3[i,-index]
        missing_est <- sum(V3) - sum(known_est)
        tot_est <- c(known_est, missing_est)
      } else {
        tot_est <- V3[i,]
        tot_cases <- notif_age_us_5yr[i,]
      }
      adj_2a  <- adj_2a +  dDirMult(M=tot_cases+0.01,n=tot_cases,Rho=rho)*wts[idx_5yr_end][i]
      print(adj_2a)
      tot_lik <- tot_lik + dDirMult(M=tot_est*1e6,n=tot_cases,Rho=rho)*wts[idx_5yr_end][i]
      }
  } else {
    adj_2a            <- sum(dDirMult(M=notif_age_us_5yr+0.01,n=notif_age_us_5yr,Rho=rho)*wts[idx_5yr_end])
    #scale does not matter for dirichlet llikelihood
    tot_lik <- sum(dDirMult(M=V3*1e6,n=notif_age_us_5yr,Rho=rho)*wts[idx_5yr_end]) - adj_2a
  }
  return(tot_lik)
}

### ### ### NUS CASES AGE DISTRIBUTION 5yrs 1993-2016  ### ### ### ### ### ### D
# Motivation: dirichlet-multinomial, multinomial data with additional non-sampling biases
notif_age_nus_5yr_lLik_st <- function(V,st,rho=0.1) { # V = table of us notifications by age 1993-2016 (row=24 years, col=11 ages)
  ### Read in the age and nativity stratified data
  ag5                   <- CalibDatCases$cases_yr_ag_nat_st_5yr[[st]]
  notif_age_nus_5yr     <- as.data.frame(ag5[ag5[,"usb"]==0,ag_5yr_cols(ag5)])
  ### Format the model estimates to match the calibration data
  V2 <- V[,-11]; V2[,10] <- V2[,10]+V[,11]
  V3<-matrix(0,5,10)
  V3[1,]<-colSums(V2[1:5,]);V3[2,]<-colSums(V2[6:10,]);V3[3,]<-colSums(V2[11:15,]); V3[4,]<-colSums(V2[16:20,]); V3[5,]<-colSums(V2[21:25,])
  ### Check for missing calibration data
  if (sum(is.na(notif_age_nus_5yr)) > 0) {
    ### Initialize the likelihood value to zero
    tot_lik <- 0; adj_2a <-0
    ### If data is missing we need to create a single new bucket for these data and estimates
    ### Read in the total US cases to use as a total
    cn               <- CalibDatCases$cases_nat_st_5yr
    notif_age_nus0   <- cn[cn$State.Code==st & cn$usb==0,nat_5yr_bucket_cols(cn)]
    for (i in 1:nrow(notif_age_nus_5yr)){
      ### Check that this row has an NA
      if (sum(is.na(notif_age_nus_5yr[i,])) > 0){
        ### Calibration data
        ### Which of the age groups are NA
        index <- which(is.na(notif_age_nus_5yr[i,]))
        ### Remove these age groups from the dataframe
        known_cases <- notif_age_nus_5yr[i,-index]
        ### Set the total of missing cases equal to the difference of total cases and known cases
        missing_cases <- notif_age_nus0[i] - sum(known_cases)
        ### Add missing cases to the dataframe
        tot_cases <- c(unlist(known_cases), unlist(missing_cases))
        ### Model estimates (match to estimates)
        ### Create a sum of the estimates that correspond to missing age groups
        known_est   <- V3[i,-index]
        missing_est <- sum(V3) - sum(known_est)
        tot_est <- c(known_est, missing_est)
      } else {
        tot_est <- V3[i,]
        tot_cases <- notif_age_nus_5yr[i,]
      }
      adj_2a  <- adj_2a +  dDirMult(M=tot_cases+0.01,n=tot_cases,Rho=rho)*wts[idx_5yr_end][i]
      tot_lik <- tot_lik + dDirMult(M=tot_est*1e6,n=tot_cases,Rho=rho)*wts[idx_5yr_end][i]
    }
  } else {
    adj_2a            <- sum(dDirMult(M=notif_age_nus_5yr+0.01,n=notif_age_nus_5yr,Rho=rho)*wts[idx_5yr_end])
    #scale does not matter for dirichlet llikelihood
    tot_lik <- sum(dDirMult(M=V3*1e6,n=notif_age_nus_5yr,Rho=rho)*wts[idx_5yr_end]) - adj_2a
  }
  return(tot_lik)
}

##smoothed estimates
notif_hr_lLik_st <- function(V,st) { # V = table of notifications by tx history (row=97:16, col=n then e)
  sid <- get_stateID()
  hr  <- CalibDatCases[["hr_cases_sm"]]
  rws <- which(hr[,"st"]==sid[st,1])
  ### smoothed high-risk fraction x total cases, per 5-year bucket
  notif_hr_5yr <- hr[rws,"fract_sm"] * hr[rws,"ss"]
  # notif_5yr          <-    c(sum(CalibDatCases[["cases_yr_st"]][[st]][2:6,2]),
  #                             sum(CalibDatCases[["cases_yr_st"]][[st]][7:11,2]),
  #                             sum(CalibDatCases[["cases_yr_st"]][[st]][12:16,2]),
  #                             sum(CalibDatCases[["cases_yr_st"]][[st]][17:21,2]),
  #                             sum(CalibDatCases[["cases_yr_st"]][[st]][22:26,2]))
  # notif_hr_5yr      <- notif_hr*notif_5yr
  adj_6           <- sum(dnorm(notif_hr_5yr,notif_hr_5yr,notif_hr_5yr*0.1/1.96,log=T)*wts[idx_5yr_start])
  sum(dnorm(notif_hr_5yr,V*1e6,notif_hr_5yr*0.1/1.96,log=T)*wts[idx_5yr_start]) - adj_6
}

#' DISTRIBUTION OF CASES RECENT TRANSMISSION VS NO RECENT TRANSMISSION
#'@param V distribution of recent and not recent transmission
#'@return likelihood
recent_trans_dist_lLik_st  <- function(V,st) {
  rct_trans_dist        <- CalibDat[["rct_cases_sm"]][st,"fract_sm"]
  adj_13          <- dbeta(rct_trans_dist,rct_trans_dist*1000,(1-rct_trans_dist)*1000,log=T)
  dbeta(rct_trans_dist,V[1]*1000,V[2]*1000,log=T)  - adj_13
}

### ### ### CASES FB RECENT ENTRY DISTRIBUTION 1993-2013  ### ### ### ### ### ### D
# Motivation: should be a normal distribution because it is based on a model result
notif_fb_rec_lLik_st<-function(V,st){
  sid <- get_stateID()
  rt  <- CalibDatState[["rt_fb_cases_sm"]]
  notif_rec<-rt[which(rt[,"st"]==sid[st,1]),"p_recent"]
  cases_ag_nat <- CalibDatState[["cases_yr_ag_nat_st"]][[st]]
  notif_fb          <-  rbind(sum(cases_ag_nat[2:6,"total","nusb"]),
                              sum(cases_ag_nat[7:11,"total","nusb"]),
                              sum(cases_ag_nat[12:16,"total","nusb"]),
                              sum(cases_ag_nat[17:21,"total","nusb"]),
                              sum(cases_ag_nat[22:26,"total","nusb"]))
  notif_fb_rec      <- notif_rec*notif_fb
  adj_6           <- sum(dnorm(notif_fb_rec,notif_fb_rec,notif_fb_rec*0.1/1.96,log=T)*wts[idx_5yr_start])
  sum(dnorm(notif_fb_rec,V*1e6,notif_fb_rec*0.1/1.96,log=T)*wts[idx_5yr_start]) - adj_6
}

# notif_fb_rec_lLik_st <- function(V,st,rho=0.005) { # V = table of notifications by rec 1993-2014 (row=22 years, col=pos then neg)
#   notif_rec<-CalibDatState[["rt_fb_cases_sm"]][which(CalibDatState[["rt_fb_cases_sm"]][,1]==stateID[st,1]),9]
#   notif_fb          <-  rbind(sum(CalibDatState[["cases_yr_ag_nat_st"]][[st]][2:6,12,"nusb"]),
#                         sum(CalibDatState[["cases_yr_ag_nat_st"]][[st]][7:11,12,"nusb"]),
#                         sum(CalibDatState[["cases_yr_ag_nat_st"]][[st]][12:16,12,"nusb"]),
#                         sum(CalibDatState[["cases_yr_ag_nat_st"]][[st]][17:21,12,"nusb"]),
#                         sum(CalibDatState[["cases_yr_ag_nat_st"]][[st]][22:26,12,"nusb"]))
#   notif_fb_rec      <- cbind(notif_rec*notif_fb, (1-notif_rec)*notif_fb)
#   adj_6             <- sum(dDirMult(M=notif_fb_rec,n=notif_fb_rec,Rho=rho)*wts[c(year_to_idx(1994), year_to_idx(1999), year_to_idx(2004), year_to_idx(2009), year_to_idx(2014))])
#   #scale does not matter for dirichlet llikelihood
#   (sum(dDirMult(M=V,n=notif_fb_rec,Rho=rho)*wts[c(year_to_idx(1994), year_to_idx(1999), year_to_idx(2004), year_to_idx(2009), year_to_idx(2014))]) - adj_6)*5
#   }

### ### ### TREATMENT OUTCOMES 1993-2012  ### ### ### ### ### ### D
# Motivation: dirichlet-multinomial, multinomial data with additional non-sampling biases

tx_outcomes_lLik_st <- function(V,rho=0.01) { # V = simulator tx-outcome output for cols 132:134, all years
  raw           <- CalibDatState[["tx_outcomes"]]
  years_all     <- as.integer(raw[, "year"])
  keep          <- years_all %in% as.integer(names(wts))
  years         <- years_all[keep]
  raw_k         <- raw[keep, ]
  tx_outcomes   <- cbind(1 - rowSums(raw_k[, 2:3]), raw_k[, 2], raw_k[, 3]) * raw_k[, 4]
  wt_slice      <- wts[as.character(years)]
  V_slice       <- V[year_to_idx(years), , drop = FALSE]
  adj_11 <- sum(dDirMult(M = tx_outcomes + 0.01, n = tx_outcomes, Rho = 0.01) * wt_slice)
  #scale does not matter for dirichlet llikelihood
  sum(dDirMult(M = V_slice, n = tx_outcomes, Rho = rho) * wt_slice) - adj_11
  }

### ### ### TOTAL LTBI TREATMENT INITS 2002  ### ### ### ### ### ### D

tltbi_tot_lLik_st   <- function(V,st) { # V = total TLTBI inits in 2002 (scalar)
  # Motivation: norm, mean centered with CI = +/- 10% of mean
  tltbi_vol        <- CalibDatState[["TLTBI_volume_state"]][[st]]
  adj_12           <- dnorm(tltbi_vol[1],tltbi_vol[1],diff(tltbi_vol[2:3])/1.96,log=T)
  dnorm(tltbi_vol[1],V*1e6,diff(tltbi_vol[2:3])/1.96,log=T) - adj_12
  }

### ### ### DISTRIBUTION OF LTBI TREATMENT INITS 2002  ### ### ### ### ### ###  D
tltbi_dist_lLik_st  <- function(V) {
  ### fractions of LTBI treatment inits that are foreign-born / high-risk
  TLTBI_dist       <- CalibDatState[["TLTBI_dist"]][c("FB","HR")]
  adj_13           <- sum( dbeta(TLTBI_dist,TLTBI_dist*100,(1-TLTBI_dist)*100,log=T) )# V = dist TLTBI inits in 2002 (vector fraction FB, HR, HV in 2002)
  sum( dbeta(TLTBI_dist,V*100,(1-V)*100,log=T) ) - adj_13  }

### ### ### LTBI PREVALENCE BY AGE 2011, US  ### ### ### ### ### ### D
# Motivation: additional prior on LTBI, using beta densities parameterized to Miramontes/Hill results

ltbi_us_11_lLik_st <- function(V) { # V = LTBI in US pop 2011 (row=11 ages, col= ltbi, non-ltbi)
  ltbi_us_11      <- CalibDatState[["LTBI_prev_US_11_IGRA"]]
  adj_15          <- sum( dbeta(ltbi_us_11[,2]/rowSums(ltbi_us_11[,2:3]),ltbi_us_11[,2],ltbi_us_11[,3],log=T) )
  V[9,] <- colSums(V[9:11,])
  (sum( dbeta(V[2:9,1]/rowSums(V[2:9,]),ltbi_us_11[,2],ltbi_us_11[,3],log=T) ) - adj_15)*2  }

ltbi_us_11_dp_lLik_st <- function(V) { # V = LTBI in US pop 2011 (row=11 ages, col= ltbi, non-ltbi)
  ltbi_us_11_dp      <- CalibDatState[["LTBI_prev_US_11_DoubPos"]]
  adj_15dp           <- sum( dbeta(ltbi_us_11_dp[,2]/rowSums(ltbi_us_11_dp[,2:3]),ltbi_us_11_dp[,2],ltbi_us_11_dp[,3],log=T) )

  V[9,] <- colSums(V[9:11,])
  (sum( dbeta(V[2:9,1]/rowSums(V[2:9,]),ltbi_us_11_dp[,2],ltbi_us_11_dp[,3],log=T) ) - adj_15dp)*2  }

### ### ### LTBI PREVALENCE BY AGE 2011, FB  ### ### ### ### ### ### D
# Motivation: multinomial adjusted to match effective sample size due to survey weighting
ltbi_fb_11_lLik_st <- function(V) { # V = LTBI in FB pop 2011 (row=11 ages, col= ltbi, non-ltbi)

  ltbi_fb_11      <- CalibDatState[["LTBI_prev_FB_11_IGRA"]]
  adj_16          <- sum( dbeta(ltbi_fb_11[,2]/rowSums(ltbi_fb_11[,2:3]),ltbi_fb_11[,2],ltbi_fb_11[,3],log=T) )
  V[9,] <- colSums(V[9:11,])
  (sum( dbeta(V[2:9,1]/rowSums(V[2:9,]),ltbi_fb_11[,2],ltbi_fb_11[,3],log=T) ) - adj_16)*2  }


ltbi_fb_11_dp_lLik_st <- function(V) { # V = LTBI in FB pop 2011 (row=11 ages, col= ltbi, non-ltbi)
  ltbi_fb_11_dp      <- CalibDatState[["LTBI_prev_FB_11_DoubPos"]]
  adj_16dp           <- sum( dbeta(ltbi_fb_11_dp[,2]/rowSums(ltbi_fb_11_dp[,2:3]),ltbi_fb_11_dp[,2],ltbi_fb_11_dp[,3],log=T) )
  V[9,] <- colSums(V[9:11,])
  (sum( dbeta(V[2:9,1]/rowSums(V[2:9,]),ltbi_fb_11_dp[,2],ltbi_fb_11_dp[,3],log=T) ) - adj_16dp)*2  }
### ### ### Total TB DEATHS 1999-2016 ### ### ### ### ### ### D

tbdeaths_lLik_st <- function(V,st) { # V = simulator TB-death output for cols 227:237, all years
  raw       <- CalibDatState[["tbdeaths"]][[st]]
  years_all <- as.integer(raw[, "Year"])
  keep      <- years_all %in% as.integer(names(wts))
  years     <- years_all[keep]
  tb_deaths <- as.numeric(raw[keep, "Deaths"])
  wt_slice  <- wts[as.character(years)]
  V_slice   <- V[year_to_idx(years), , drop = FALSE]
  V2        <- rowSums(V_slice) * 1e6
  ok        <- !is.na(tb_deaths)
  adj_19 <- sum((dnorm(tb_deaths, tb_deaths, tb_deaths*0.1/1.96, log=T) * wt_slice)[ok])
  sum((dnorm(tb_deaths, V2, tb_deaths*0.1/1.96, log=T) * wt_slice)[ok]) - adj_19
}
### ### ### ANN DECLINE IN TB DEATHS 1968-2015  ### ### ### ### ### ### D

tbdeaths_decline_lLik_st <- function(V) { # V = vector of tb deaths 1968-2015
  tbdeaths_decline      <- CalibDatState[["deaths_ann_decline_68_15"]]
  adj_19a            <- sum(dnorm(tbdeaths_decline,tbdeaths_decline,0.015/1.96,log=T))
  V2 <- (1-(V[48]/V[1])^(1/47))
  sum(dnorm(tbdeaths_decline,V2,0.015/1.96,log=T)) - adj_19a
}
### ### ### TB DEATHS AGE DISTRIBUTION 1999-2016  ### ### ### ### ### ### D
# Motivation: dirichlet-multinomial, multinomial data with additional non-sampling biases

tb_dth_age_lLik_st <- function(V,rho=0.005) { # V = simulator TB-death output for cols 227:237, all years
  raw           <- CalibDatState[["tbdeaths_age_yr"]]
  years_all     <- as.integer(raw[, "year"])
  keep          <- years_all %in% as.integer(names(wts))
  years         <- years_all[keep]
  tb_deaths_age <- raw[keep, -1]
  wt_slice      <- wts[as.character(years)]
  V_slice       <- V[year_to_idx(years), , drop = FALSE]
  adj_19b <- sum(dDirMult(M = tb_deaths_age + 0.005, n = tb_deaths_age, Rho = rho) * wt_slice)
  V2 <- V_slice[, -11]; V2[, 10] <- V2[, 10] + V_slice[, 11]
  sum(dDirMult(M = V2, n = tb_deaths_age, Rho = rho) * wt_slice) - adj_19b
}
### ### ### TOTAL POP EACH DECADE, FOR FB  ### ### ### ### ### ###  D
# Motivation: norm, mean centered with CI = +/- 2 million wts[1+0:6*10]

tot_pop_yr_fb_lLik_st <- function(V,st) { # V = total pop (rows=year, cols=us, fb)
  tot_pop_yr      <- CalibDatState[["pop_50_10"]][[st]]
  tot_pop_yr_fb   <- tot_pop_yr[tot_pop_yr[,2]==0,]
  #get 2019 population
  pop_ag_11_190  <- CalibDatState[["pop_00_19"]][[st]][,c("age_group","usb","2019")]
  #get 2019 fb population
  pop_ag_11_19nus <-sum(pop_ag_11_190[pop_ag_11_190[,2]==0,3][-11])
  #append the foreign born population
  tot_pop_yr_fb   <- c(tot_pop_yr_fb[,-c(1:2)], pop_ag_11_19nus)
  # if (loc != "HI" & loc != "AK"){
  adj_17          <- sum(dnorm(tot_pop_yr_fb[-1],tot_pop_yr_fb[-1],tot_pop_yr_fb[7]*0.05/1.96,log=T)*wts[idx_decades])
  #total population is in real numbers so we need to scale up output
  sum(dnorm(tot_pop_yr_fb[-1],V[idx_decades]*1e6,tot_pop_yr_fb[7]*0.05/1.96,log=T)*wts[idx_decades]) - adj_17}
  # else{
  #   adj_17          <- sum(dnorm(tot_pop_yr_fb,tot_pop_yr_fb,tot_pop_yr_fb[7]*0.1/1.96,log=T)*wts[c(1+1:6*10,68)])
  #   #total population is in real numbers so we need to scale up output
  #   sum(dnorm(tot_pop_yr_fb[-1],V[c(11,21,31,41,51,61,68)]*1e6,tot_pop_yr_fb[7]*0.1/1.96,log=T)*wts[c(1+1:6*10,68)]) - adj_17
  # }
# }
### D
tot_pop_yr_us_lLik_st_00_10 <- function(V,st) {
  tot_pop_yr      <- CalibDatState[["pop_50_10"]][[st]]
  # V = total pop (rows=year, cols=us, fb)
  tot_pop_yr_us  <- tot_pop_yr[tot_pop_yr[,2]==1,3]
  # tot_pop_yr_us<-as.matrix(tot_pop_yr_us)
  # tot_pop_yr_us<-colSums(as.matrix(tot_pop_yr_us)[,-c(1:2)])
  adj_17b        <- sum(dnorm(tot_pop_yr_us[6:7],tot_pop_yr_us[6:7],tot_pop_yr_us[7]*0.05/1.96,log=T)*wts[c(year_to_idx(2000), year_to_idx(2010))])
  sum(dnorm(tot_pop_yr_us[6:7],V[c(year_to_idx(2000), year_to_idx(2010))]*1e6,tot_pop_yr_us[7]*0.05/1.96,log=T)*wts[c(year_to_idx(2000), year_to_idx(2010))]) - adj_17b  } # CI = +/- 2mil

tot_pop_yr_us_lLik_st <- function(V,st) {

  tot_pop_yr      <- CalibDatState[["pop_50_10"]][[st]]
  tot_pop_yr_us   <- tot_pop_yr[tot_pop_yr[,2]==1,]
  #get 2019 population
  pop_ag_11_190  <- CalibDatState[["pop_00_19"]][[st]][,c("age_group","usb","2019")]
  #get 2019 fb population
  pop_ag_11_19us <-sum(pop_ag_11_190[pop_ag_11_190[,2]==1,3][-11])
  #append the foreign born population
  tot_pop_yr_us   <- c(tot_pop_yr_us[,-c(1:2)], pop_ag_11_19us)
  adj_17          <- sum(dnorm(tot_pop_yr_us[-1],tot_pop_yr_us[-1],tot_pop_yr_us[7]*0.05/1.96,log=T)*wts[idx_decades])
  #total population is in real numbers so we need to scale up output
  sum(dnorm(tot_pop_yr_us[-1],V[idx_decades]*1e6,tot_pop_yr_us[7]*0.05/1.96,log=T)*wts[idx_decades]) - adj_17}

### ### ### TOTAL POP AGE DISTRIBUTION 2017  ### ### ### ### ### ### D
# Motivation: reported estimates represent pseudo-data for a multinomial likelihood, with ESS = 500
tot_pop19_ag_fb_lLik_st <- function(V,st,ESS=500) { # V =  US pop in 2014 (row=11 ages, col= us, fb)
  pop_ag_11_190  <- CalibDatState[["pop_00_19"]][[st]][,c("age_group","usb","2019")]
  pop_ag_11_19us <-pop_ag_11_190[pop_ag_11_190[,2]==1,3][-11]
  pop_ag_11_19nus <-pop_ag_11_190[pop_ag_11_190[,2]==0,3][-11]

  pop_ag_11_19   <- cbind(pop_ag_11_19us/sum(pop_ag_11_19us)+.001, pop_ag_11_19nus/sum(pop_ag_11_19nus)+.001)
  adj_18         <- (sum(log(pop_ag_11_19[,1])*pop_ag_11_19[,1])+sum(log(pop_ag_11_19[,2])*pop_ag_11_19[,2]))*ESS
  V1 <- rbind(V[1:9,],V[10,]+V[11,])
  V2<-cbind(V1[,1]/sum(V1[,1]),V1[,2]/sum(V1[,2]))
  (sum(log(V2[,1])*pop_ag_11_19[,1])+sum(log(V2[,2])*pop_ag_11_19[,2]))*ESS - adj_18
}

### ### ### TOTAL POP AGE DISTRIBUTION 2017  ### ### ### ### ### ### D
# Motivation: reported estimates represent pseudo-data for a multinomial likelihood, with ESS = 500
tot_pop1719_ag_fb_lLik_st <- function(V,st,ESS=500) { # V =  US pop in 2014 (row=11 ages, col= us, fb)
  pop_ag_11_190  <- CalibDatState[["pop_00_19"]][[st]][,c("age_group","usb","2017","2018","2019")]
  pop_ag_11_19us <- rowSums(pop_ag_11_190[pop_ag_11_190[,2]==1,3:5])[-11]
  #top code at 75p
  # pop_ag_11_19us[9] <- pop_ag_11_19us[9] + pop_ag_11_19us[10]
  # pop_ag_11_19us <- pop_ag_11_19us[-10]
  pop_ag_11_19nus <-rowSums(pop_ag_11_190[pop_ag_11_190[,2]==0,3:5])[-11]
  #top code at 75p
  pop_ag_11_19   <- cbind(pop_ag_11_19us/sum(pop_ag_11_19us), pop_ag_11_19nus/sum(pop_ag_11_19nus))
  adj_18         <- (sum(log(pop_ag_11_19[,1])*pop_ag_11_19[,1])+sum(log(pop_ag_11_19[,2])*pop_ag_11_19[,2]))*ESS
  V1 <- rbind(V[1:9,],V[10,]+V[11,])
  V2<-cbind(V1[,1]/sum(V1[,1]),V1[,2]/sum(V1[,2]))
  (sum(log(V2[,1])*pop_ag_11_19[,1])+sum(log(V2[,2])*pop_ag_11_19[,2]))*ESS - adj_18
}

#' TOTAL ALL-CAUSE DEATHS, MOST RECENT AVAILABLE YEAR
#' Motivation: norm, mean centered with CI = +/- 10% of mean
#'@name dth_tot_lLik_st
#'@param V simulator all-cause deaths by age (rows = years, 11 age columns)
#'@return likelihood
dth_tot_lLik_st <- function(V,st) {
  dba           <- DeathByAge[[st]]
  yrs           <- as.integer(dba[, 1])
  year          <- max(yrs[yrs %in% as.integer(names(wts))])
  ST_deaths_tot <- dba[dba_year_rows(dba, year), dba_total_col(dba)]
  wt_slice      <- wts[as.character(year)]
  V_slice       <- rowSums(V[year_to_idx(year), , drop = FALSE])
  adj_20a       <- sum(dnorm(ST_deaths_tot, ST_deaths_tot, ST_deaths_tot*0.1/1.96, log=T) * wt_slice)
  sum(dnorm(ST_deaths_tot, V_slice*1e6, ST_deaths_tot*0.1/1.96, log=T) * wt_slice) - adj_20a
}

#' ALL-CAUSE DEATHS AGE DISTRIBUTION, MOST RECENT TWO AVAILABLE YEARS
#' Motivation: dirichlet-multinomial, multinomial data with additional non-sampling biases
#'@param V simulator all-cause deaths by age (rows = years, 11 age columns)
#'@param rho correlation parameter
#'@return likelihood
tot_dth_age_lLik_st <- function(V,st,rho=0.01) {
  dba      <- DeathByAge[[st]]
  yrs      <- as.integer(dba[, 1])
  yrs      <- tail(sort(yrs[yrs %in% as.integer(names(wts))]), 2)
  tda      <- dba[dba_year_rows(dba, yrs), -c(1, dba_total_col(dba)), drop = FALSE]
  wt_slice <- wts[as.character(yrs)]
  Vy       <- V[year_to_idx(yrs), , drop = FALSE]
  V2       <- Vy[, -11, drop = FALSE]; V2[, 10] <- V2[, 10] + Vy[, 11]
  adj_20b  <- sum(dDirMult(M=tda+0.1, n=tda, Rho=rho) * wt_slice)
  sum(dDirMult(M=V2, n=tda, Rho=rho) * wt_slice) - adj_20b
}

#' Mortality Risk Group Distribution 1999-2014
#' Motivation: dirichlet-multinomial, multinomial data with additional non-sampling biases
#'@name mort_dist_lLik_st
#'@param V table of mort_dist 1999-2014 (row=16 years, col=11 ages)
#'@param rho correlation parameter
#'@return likelihood
mort_dist_lLik_st <- function(V,rho=0.1) {
  md     <- rowSums(dist_gen)
  mort_dist     <-matrix(md,length(year_to_idx(2015):year_to_idx(2016)),4, byrow = TRUE)
  adj_21        <- sum(dDirMult(M=mort_dist+0.01,n=mort_dist,Rho=0.1)*wts[year_to_idx(2015):year_to_idx(2016)])
  tot_lik<-0
  for(ag in 1:11){
    V1<-V[,(1:4)+4*(ag-1)]
    x<-sum(dDirMult(M=(V1*1e6),n=mort_dist,Rho=rho)*wts[year_to_idx(2015):year_to_idx(2016)]) - adj_21
    tot_lik<-tot_lik+x
    # print(x)
  }
  tot_lik<-tot_lik
  return(tot_lik)
}

#' Mortality Risk Group Distribution 1999-2014
#' Motivation: dnorm
#'@name mort_dist_lLik_norm_st
#'@param V table of mort_dist 1999-2014 (row=16 years, col=11 ages)
#'@param rho correlation parameter
#'@return likelihood
mort_dist_lLik_norm_st <- function(V) {
  md     <- rowSums(dist_gen)
  mort_dist     <-matrix(md,length(year_to_idx(2000):year_to_idx(2016)),4, byrow = TRUE)
  adj_21b        <- sum(dnorm(mort_dist,mort_dist,mort_dist*0.1/1.96, log=T)*wts[year_to_idx(2000):year_to_idx(2016)])
  tot_lik<-0
  for(ag in 1:11){
    V1<-V[,(1:4)+4*(ag-1)]
    x<-sum(dnorm(mort_dist, V1,mort_dist*0.1/1.96, log=T)*wts[year_to_idx(2000):year_to_idx(2016)]) - adj_21b
    tot_lik<-tot_lik+x
    # print(x)
  }
  tot_lik<-tot_lik
  return(tot_lik)
}
### ### ### HOMELESS POP 2010  ### ### ### ### ### ### names(CalibDatState)
# Motivation: norm, mean centered with CI = +/- 25% of mean

homeless_10_lLik_st <- function(V,st) { # V = homeless pop in 2010 (scalar)
  homeless_pop      <- CalibDatState[["homeless_pop"]][[st]]
  adj_23b          <- dnorm(homeless_pop[1],homeless_pop[1],diff(homeless_pop[2:3])/2/1.96,log=T)
  dnorm(homeless_pop[1],V,diff(homeless_pop[2:3])/2/1.96,log=T) - adj_23b   }

### ### ### LIKELIHOOD FOR BORGDORFF ESTIMATES  ### ### ### ### ### ###
borgdorff_lLik_st <- function(Par_list,N_red=1) {  # Par_list = list(Mpfast[,c(1,3,2,4)], Mrslow[,c(1,3,2,4)], rfast, rRecov)
  ss_borgdorff   <- 854.5463/4
  datB           <- CalibDatState[["borgdorff_data"]]
  adj_24         <- sum(diff(-datB[,2])*ss_borgdorff*log(diff(-datB[,2])) + (1-diff(-datB[,2]))*ss_borgdorff*log(1-diff(-datB[,2])))
  zz <- tryCatch({
    Mpfast     <- Par_list[[1]];
    # ORpfastRF <- Par[2];
    Mrslow     <- Par_list[[2]]*12;
    # RRrSlowRF <- Par[4];
    rfast     <- Par_list[[3]]*12;
    rRecov    <- Par_list[[4]]*12;
    pfast_v <- rslow_v <- rep(NA,4)
    pfast_v <- Mpfast[3,]
    rslow_v <- Mrslow[3,]
    p0 <- matrix(NA,4,10)  # this added
    for (i in 1:4){
      p0[i,] <- pfast_v[i] *(1-(1-exp(-(rfast+rRecov)*datB[,1]))*(rfast/(rfast+rRecov))) +
        (1-pfast_v[i])*(1-(1-exp(-(rslow_v[i]+rRecov)*datB[,1]))*(rslow_v[i]/(rslow_v[i]+rRecov)))
    }
    # p1 <- as.numeric(t(p0)%*%v21a[17,41:44])

    p1 <- as.numeric(t(p0)%*%colSums(dist_gen))
    p <- 1-(1-p1)/(1-p1)[nrow(datB)]
    sum(diff(-datB[,2])*ss_borgdorff*log(diff(-p)) + (1-diff(-datB[,2]))*ss_borgdorff*log(1-diff(-p)))/N_red - adj_24/N_red
  },error=function(e) -Inf )
  if(is.nan(zz)) { zz = -10^4
  } else {
    if(zz== -Inf) zz = -10^4 }
  zz }

### ### ### LIKELIHOOD FOR FEREBEE ESTIMATES  ### ### ### ### ### ###
ferebee_lLik_st <- function(Par_list,N_red=4) {
  datF         <- CalibDatState[["ferebee_data"]]
  adj_25       <- sum(datF[,3]*log(datF[,3]/datF[,2]) + (datF[,2]-datF[,3])*log(1-datF[,3]/datF[,2]))
  n_yr_F       <- nrow(datF)
  zz <- tryCatch({

    Mpfast     <- Par_list[[1]];
    # ORpfastRF <- Par[2];
    Mrslow     <- Par_list[[2]]*12;
    # RRrSlowRF <- Par[4];
    rfast     <- Par_list[[3]]*12;
    rRecov    <- Par_list[[4]]*12;
    pfast_v <- rslow_v <- rep(NA,4)
    pfast_v <- Mpfast[3,]
    rslow_v <- Mrslow[3,]
    p0 <- matrix(NA,4,11)  # this added
    for (i in 1:4){
      p0[i,] <- pfast_v[i] *(1-(1-exp(-(rfast+rRecov)*(0:n_yr_F)))*(rfast/(rfast+rRecov))) +
        (1-pfast_v[i])*(1-(1-exp(-(rslow_v[i]+rRecov)*(0:n_yr_F)))*(rslow_v[i]/(rslow_v[i]+rRecov)))
    }
    p1 <- as.numeric(t(p0)%*%colSums(dist_gen))
    r2 <- -log(1-diff(-p1)/p1[-(n_yr_F+1)])/1
    sum((datF[,3]*log(r2) + (datF[,2]-datF[,3])*log(1-r2)))/N_red - adj_25/N_red
  },error=function(e) -Inf )
  if(is.nan(zz)) {
    zz = -10^4
  } else {
    if(zz== -Inf) zz = -10^4 }
  zz  }


### ### ### LIKELIHOOD FOR SUTHERLAND ESTIMATES  ### ### ### ### ### ###
sutherland_lLik_st <- function(Par_list,N_red=4) {
  datS            <- CalibDatState[["sutherland_data"]]
  datSz           <- datS; datSz[datS[,3]==0,3] <- 0.01
  adj_26          <- sum(datSz[,3]*log(datSz[,3]/datSz[,2]) + (datSz[,2]-datSz[,3])*log(1-datSz[,3]/datSz[,2]))
  n_yr_S          <- nrow(datS)
  zz <- tryCatch({
    Mpfast     <- Par_list[[1]];
    # ORpfastRF <- Par[2];
    Mrslow     <- Par_list[[2]]*12;
    # RRrSlowRF <- Par[4];
    rfast     <- Par_list[[3]]*12;
    rRecov    <- Par_list[[4]]*12;
    s<-pfast_v <- rslow_v <- rep(NA,4)
    pfast_v <- Mpfast[3,]
    rslow_v <- Mrslow[3,]
    p0 <- matrix(NA,4,16)  # this added
    for (i in 1:4){
      p0[i,] <- pfast_v[i] *(1-(1-exp(-(rfast+rRecov)*(0:n_yr_S)))*(rfast/(rfast+rRecov))) +
        (1-pfast_v[i])*(1-(1-exp(-(rslow_v[i]+rRecov)*(0:n_yr_S)))*(rslow_v[i]/(rslow_v[i]+rRecov)))
    }
    p1 <- as.numeric(t(p0)%*%colSums(dist_gen))
    r2 <- -log(1-diff(-p1)/p1[-(n_yr_S+1)])/1
    sum((datS[,3]*log(r2) + (datS[,2]-datS[,3])*log(1-r2)))/N_red - adj_26/N_red
  },error=function(e) -Inf )
  if(is.nan(zz)) {
    zz = -10^4
  } else {
    if(zz== -Inf) zz = -10^4 }
  zz  }


### ### ### LIKELIHOOD FOR TIEMERSMA ESTS ### ### ### ### ### ###
tiemersma_lLik_st <- function(Par) { # Par= c(rSlfCur,muIp)
  adj_27         <- dnorm(3.0,3.0,0.5/1.96,log=T)+dbeta(0.9,0.9*50,(1-0.9)*50,log=T)
  rSlfCur <- Par[1]*12;
  muIp    <- Par[2]*12 #muIp is the value for all TB active disease (not just smear pos)

  dur <- (1/(rSlfCur+muIp))
  cf_sp <- muIp/(rSlfCur+muIp);

  l1 <- dnorm(3.0,dur,0.5/1.96,log=T)
  l2 <- dbeta(0.9,cf_sp*50,(1-cf_sp)*50,log=T)
  l1+l2 - adj_27
}
