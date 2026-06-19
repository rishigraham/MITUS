make_base <- function() {
  list(
    cases_yr_st = lapply(seq_len(51), function(i) matrix(i, nrow = 2, ncol = 2)),
    cases_nat_st_5yr = data.frame(
      State = rep(c("Alabama","California"), each = 3),
      State.Code = rep(c(1, 5), each = 3),
      usb = rep(c(0, 1, 3), 2),
      X1996.2000 = 1:6,
      stringsAsFactors = FALSE
    ),
    hr_cases_sm = data.frame(
      st = factor(rep(c("Alabama","California"), each = 5),
                  levels = c("Alabama","California")),
      yr = rep(1:5, 2),
      hr = (1:10) / 10
    ),
    rt_fb_cases_sm = data.frame(
      st = factor(rep(c("Alabama","California"), each = 5),
                  levels = c("Alabama","California")),
      yr = rep(1:5, 2),
      p_recent = (1:10) / 100
    ),
    deaths_ann_decline_68_15 = 0.05
  )
}

test_that("state-indexed list fields receive overlay at slot st", {
  base <- make_base()
  overlay <- list(cases_yr_st = matrix(99, nrow = 2, ncol = 2))
  out <- merge_county_calibdat(base, overlay, st = 52)
  expect_equal(out$cases_yr_st[[52]], overlay$cases_yr_st)
  expect_equal(out$cases_yr_st[[5]], base$cases_yr_st[[5]])
})

test_that("by-row state-keyed tables append overlay rows", {
  base <- make_base()
  sd_name <- "San Diego County"
  overlay <- list(
    cases_nat_st_5yr = data.frame(
      State = sd_name, State.Code = 52, usb = c(0, 1, 3), X1996.2000 = 7:9,
      stringsAsFactors = FALSE
    ),
    hr_cases_sm = data.frame(
      st = factor(rep(sd_name, 5),
                  levels = c(levels(base$hr_cases_sm$st), sd_name)),
      yr = 1:5, hr = (11:15) / 10
    ),
    rt_fb_cases_sm = data.frame(
      st = factor(rep(sd_name, 5),
                  levels = c(levels(base$rt_fb_cases_sm$st), sd_name)),
      yr = 1:5, p_recent = (11:15) / 100
    )
  )
  out <- merge_county_calibdat(base, overlay, st = 52)

  expect_equal(nrow(out$cases_nat_st_5yr), nrow(base$cases_nat_st_5yr) + 3)
  expect_equal(nrow(out$hr_cases_sm), nrow(base$hr_cases_sm) + 5)
  expect_equal(nrow(out$rt_fb_cases_sm), nrow(base$rt_fb_cases_sm) + 5)

  expect_equal(
    nrow(out$cases_nat_st_5yr[out$cases_nat_st_5yr$State.Code == 52, ]), 3
  )
  expect_equal(nrow(out$hr_cases_sm[out$hr_cases_sm$st == sd_name, ]), 5)
  expect_equal(nrow(out$rt_fb_cases_sm[out$rt_fb_cases_sm$st == sd_name, ]), 5)

  expect_equal(
    nrow(out$cases_nat_st_5yr[out$cases_nat_st_5yr$State.Code == 5, ]), 3
  )
})

test_that("other non-indexed fields are replaced by overlay value", {
  base <- make_base()
  overlay <- list(deaths_ann_decline_68_15 = 0.10)
  out <- merge_county_calibdat(base, overlay, st = 52)
  expect_equal(out$deaths_ann_decline_68_15, 0.10)
})

test_that("missing overlay fields leave base values untouched", {
  base <- make_base()
  overlay <- list()
  out <- merge_county_calibdat(base, overlay, st = 52)
  expect_identical(out, base)
})
