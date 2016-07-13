
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
END
COMMENT Hello Bird
MISSING -1
WEIGHT Present
WEIGHTING off
SERIALCOR on
OVERDISP on
BASETIME 1
MODEL 3
COVARIATES 2
OUTPUTFILES F
RUN", con=f)


test_that("read_tcf parses tcf files", {
  x <- read_tcf(f)[[1]]
  expect_equal(x$file,"F:/TRIM/Skylark.dat")
  expect_equal(x$title,"Skylark.dat")
  expect_equal(x$ntimes,8L)
  expect_equal(x$ncovars,2L)
  expect_equal(x$labels,c("HABITAT","Cov2"))
  expect_equal(x$missing, -1L)
  expect_equal(x$weight,"Present")
  expect_equal(x$weighting,"off")
  expect_equal(x$comment,"Hello Bird")
  expect_equal(x$serialcor,"on")
  expect_equal(x$overdisp,"on")
  expect_equal(x$basetime,1L)
  expect_equal(x$model,3L)
  expect_equal(x$covariates,2L)
  expect_equal(x$changepoints,integer(0))
  expect_equal(x$stepwise,"off")
  expect_equal(x$outputfiles,"F")
})

tryCatch(unlink(f),error=function(e)cat(sprintf("Could not unlinke temporary file %s",f)))

test_that("parsing multi-model files",{})




