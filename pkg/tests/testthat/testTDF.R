testthat::context("Reading TRIM input files")

f <- tempfile()
# happy flow
writeLines("1 1992 186 1.0000 1
1 1993 60 1.0000 1
1 1994 39 1.0000 1
1 1995 3 1.0000 1
1 1996 3 1.0000 1
1 1997 0 1.0000 1
1 1998 0 1.0000 1
1 1999 0 1.0000 1
1 2000 0 1.0000 1
1 2001 -1 1.0000 1",con=f)

test_that("reading happy flow",{
  dat <- read_tdf(file=f,missing_code=-1L,nsnif=5)
  expect_equal(dat[1,],
      data.frame(
        site=1
        ,time=1992
        ,count=186
        ,weight=1
        ,covar01=1
      ))
  expect_equal(dat[10,3],NA_integer_)  
})




