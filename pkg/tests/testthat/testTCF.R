
context("Test tcf parsing")

f <- tempfile()
writeLines("
FILE F:\\TRIM\\Skylark.dat
TITLE Skylark.dat
NTIMES 8
NCOVARS 2
LABELS
 HABITAT
 Cov2
End
Comment Hello Bird
MISSING -1
WEIGHT Present
WEIGHTING off
SERIALCOR on
OVERDISP on
BASETIME 1
MODEL 3
Covariates 2
outputfiles F
RUN", con=f)


test_that("read_tcf parses tcf files", {
  x <- read_tcf(f)
  expect_equal(x@origin,f)
  expect_equal(x@file,"F:/TRIM/Skylark.dat")
  expect_equal(x@title,"Skylark.dat")
  expect_equal(x@ntimes,8L)
  expect_equal(x@ncovars,2L)
  expect_equal(x@labels,c("HABITAT","Cov2"))
  expect_equal(x@weight,TRUE)
  expect_equal(x@weighting,FALSE)
  expect_equal(x@comment,"Hello Bird")
  expect_equal(x@serialcor,TRUE)
  expect_equal(x@overdisp,TRUE)
  expect_equal(x@basetime,1L)
  expect_equal(x@model,3L)
  expect_equal(x@covariates,2L)
  expect_equal(x@changepoints,integer(0))
  expect_equal(x@stepwise,FALSE)
  expect_equal(x@outputfiles,"F")
  expect_equal(x@run,TRUE)
})
tryCatch(unlink(f),error=function(e)cat(sprintf("Could not unlinke temporary file %s",f)))