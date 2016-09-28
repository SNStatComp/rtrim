

context("data checkers")

test_that("basic assertions",{
  expect_error(assert_positive(-1,"foo"),regexp = "foo")
  expect_error(assert_positive(0,"foo"),regexp = "foo")
  expect_true(assert_positive(1,"foo"))  

  expect_error(assert_increasing(c(3,1,4)) )
  expect_true(assert_increasing(1:2))
  
  expect_error(assert_sufficient_counts(count=0:2,index = 1:3))
  expect_true(assert_sufficient_counts(count=1:3,index=1:3))
})




test_that("Sufficient data for piecewise linear trend model (Model 2)",{
  
  ## no covariates
  d <- data.frame(count=rep(0:10,2), site=rep(1:2,11), time=rep(1:11,each=2))
  expect_true(assert_plt_model(d$count, d$time, changepoints=c(1,4,10),covars=list()))
  expect_true(assert_plt_model(d$count, d$time, changepoints=integer(0),covars=list()))
  d$count[21:22] <- 0
  expect_error(assert_plt_model(d$count, d$time, changepoints=c(1,4,10),covars=list()),regexp="10")
  
  ## with covariates
  d <- data.frame(
    count = c(0,1,7,3)
    , time = rep(1:2,each=2)
    , cov=c("a","b","a","a")
  )
  expect_true(assert_plt_model(d$count,d$time,changepoints=1, covars=list(cov=d$cov)))
  d <- data.frame(
    count = c(2,1,7,3,0,2)
    , time = rep(1:3,each=2)
    , cov=c("a","b","a","b","a","b")
  )
  expect_error(
    assert_plt_model(d$count,d$time,changepoints=c(1,2), covars=list(cov=d$cov))
  )
})


test_that("Sufficient data for model 3 with covariates",{
  
  d <- data.frame(count = rep(0:1,2), time = rep(1:2,each=2), cov = rep(1:2,2))
  expect_error(assert_covariate_counts(d$count, d$time, d["cov"]),regexp="cov")
  d <- data.frame(count = rep(7:8,2),time = rep(1:2,each=2), cov = rep(1:2,2))
  expect_true(assert_covariate_counts(d$count, d$time, d["cov"]))
  
  # with errors in two covariates
  
  d <- data.frame(count = rep(0:1,2), time = rep(1:2,each=2)
                  , covA = rep(1:2,2)
                  , covB = rep(letters[1:2],2))
  expect_error(assert_covariate_counts(d$count, d$time, d[3:4]),regexp="covB")
  expect_true(assert_covariate_counts(d$count, d$time, data.frame()))
  
  
})


context("Autodelete")

test_that("autodelete w/o covariates",{
  d <- data.frame(count = rep(1,10),time=1:10)
  # case nothing to delete.
  expect_equal(
    autodelete(count=d$count, time=d$time, changepoints=c(4,7),covars=list())
    , c(4,7)    
  )
  # case something to delete
  d[5:7,1] <- 0
  expect_equal(
    autodelete(count = d$count, time = d$time
               , changepoints = c(4,7), covars=list())
    , 7
    )
})

test_that("autodelete with covariates",{
  d <- data.frame(
    time = 1:10
    , covar = rep(letters[1:2], times=5)
    , count = rep(1,10)
  )
  # case all fine
  expect_equal(
    autodelete(count = d$count, time = d$time, changepoints=c(4,7),covars=list(cov=d$covar))
  , c(4,7)  
  )  
  # case delete 7
  d$count[9] <- 0  
  expect_equal(
    autodelete(count = d$count, time = d$time, changepoints=c(4,7),covars=list(cov=d$covar))
    , 4  
  )  
  # changepoints with explicit 1 at the beginning.
  expect_equal(
    autodelete(count = d$count, time = d$time, changepoints=c(1,4,7),covars=list(cov=d$covar))
    ,c(1, 4)  
  )  
  
})






