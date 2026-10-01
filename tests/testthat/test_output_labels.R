test_that("NM x I output column labels match the written values", {
  # The 16-column NM x I blocks are written with NM as the fastest index. The same run
  # carries subtotals by M (N_US_M*, N_NUS_M*, N_ag_k_M*) and by I (N_US_I*, N_NUS_I*)
  # from separate loops, so summing a block's columns by label must reproduce them.
  loc <- "CA"
  model_load(loc)
  M <- OutputsZint(samp_i = 1, ParMatrix = P, loc = loc, startyr = 1950, endyr = 2020,
                   prg_chng = def_prgchng(P), ttt_list = def_ttt())
  M <- as.matrix(M)
  col_sum <- function(pattern) rowSums(M[, grep(pattern, colnames(M)), drop = FALSE])
  tol <- 1e-8

  for (j in 1:4) {
    expect_equal(col_sum(sprintf("^N_NM%d_I[1-4]$", j)),
                 unname(M[, sprintf("N_US_M%d", j)] + M[, sprintf("N_NUS_M%d", j)]), tolerance = tol)
    expect_equal(col_sum(sprintf("^N_NM[1-4]_I%d$", j)),
                 unname(M[, sprintf("N_US_I%d", j)] + M[, sprintf("N_NUS_I%d", j)]), tolerance = tol)
  }

  ages <- c("0-4", "5-14", "15-24", "25-34", "35-44", "45-54", "55-64", "65-74", "75-84", "85-94", "95p")
  for (a in seq_along(ages)) {
    for (j in 1:4) {
      pattern <- sprintf("^(%%_)?%s_NM%d_I[1-4]$", ages[a], j)   # two 0-4 labels carry a "%_" prefix
      expect_length(grep(pattern, colnames(M)), 4)
      expect_equal(col_sum(pattern), unname(M[, sprintf("N_ag_%d_M%d", a, j)]), tolerance = tol)
    }
  }
})
