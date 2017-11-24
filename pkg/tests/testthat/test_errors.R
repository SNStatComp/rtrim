context("formula interface errors")

test_that("invalid formulas", {

  df <- data.frame(count=1, site=2, year=3, month=4, habitat=5, acid=5)

  # expect_error(trim(count + onzin, df),       "object .* not found")
  # expect_error(trim(count ~ site + year),     "no data given", fixed=TRUE)
  # expect_error(trim(count ~ site + year, pi), "data should be a data frame", fixed=TRUE)
  # expect_error(trim(count ~ site, df),        "model should have form 'count ~ site + year ...", fixed=TRUE)
  # expect_error(trim(count + site ~ year, df), "model should have form 'count ~ ...'", fixed=TRUE)
  # expect_error(trim(count ~ site - year, df), "Model contains unallowed operator: -", fixed=TRUE)
  # expect_error(trim(count ~ site + year + habitat : acid, df), "covariates should be included using the '+' operator (':' found)", fixed=TRUE)
})

