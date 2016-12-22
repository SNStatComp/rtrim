

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


test_that("check_observations",{
  # model 3, no covariates
  d <- data.frame(time = 1:3, count=0:2)
  out <- check_observations(d,model=3)
  expect_false(out$sufficient)
  expect_equal(out$errors,list(time=1))

  # model 3 with covariates
  d <- data.frame(count = rep(0:1,2), time = rep(1:2,each=2), cov = rep(1:2,2))
  out <- check_observations(d, model=3,covars="cov")
  expect_false(out$sufficient)
  expect_equal(out$errors,list(cov=data.frame(time=factor(1:2),cov=factor(c(1,1),levels=1:2))) )

  # model 2, no covariates
  d <- data.frame(time = 1:10, count = c(rep(1,7),rep(0,3)))
  out <- check_observations(d,model=2,changepoints = c(4,7))
  expect_false(out$sufficient)
  expect_equal(out$errors$changepoint, 7)
  
  # model 2, with covariates
  
  d <- data.frame(
    time=1:4
    , X = c("A","A","A","B")
    , count = c(1,1,1,0)
  )
  
  out <- check_observations(d,model=2, covars = "X",changepoints=1)
  expect_false(out$sufficient)
  expect_equal(out$errors$X, data.frame(changepoint=factor(1),X=factor("B",levels=c("A","B")) ))
})


test_that("Sufficient data for piecewise linear trend model (Model 2)",{

  ## no covariates
  d <- data.frame(count=rep(0:10,2), site=rep(1:2,11), time=rep(1:11,each=2))
  expect_true(assert_plt_model(d$count, d$time, changepoints=c(1,4,10),covars=list()))
  expect_true(assert_plt_model(d$count, d$time, changepoints=integer(0),covars=list()))
  #d$count[21:22] <- 0 # BUG: these are AFTER changepoint 10, so don't trigger an error
  d$count[19:22] <- 0 # corrected
  expect_error(assert_plt_model(d$count, d$time, changepoints=c(1,4,10),covars=list()),regexp="10")

  ## with covariates
  d <- data.frame(
    count = c(0,1,7,3)
    , time = rep(1:2,each=2)
    , cov=c("a","b","a","a")
  )
  #expect_true(assert_plt_model(d$count,d$time,changepoints=1, covars=list(cov=d$cov)))
  expect_error(assert_plt_model(d$count,d$time,changepoints=1, covars=list(cov=d$cov)))

  d <- data.frame(
    count = c(2,1,7,3,0,2) # BUG: passes, because the count '7' provided enough data for changepoint '2' / cat 'a'
    # count = c(2,1,0,3,0,2) # fixed
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
    , 4
    )
})

test_that("autodelete with covariates",{
  d <- data.frame(
    time = rep(1:10, times=2)
    , covar = rep(letters[1:2], each=10)
    , count = rep(1, 20)
  )

  # case all fine
  # BUG: there is no real time 1 for class b
  expect_equal(
    autodelete(count = d$count, time = d$time, changepoints=c(4,7),covars=list(cov=d$covar))
  , c(4,7)
  )

   #case delete 7 
   d$count[8:10] <- 0
   expect_equal(
     autodelete(count = d$count, time = d$time, changepoints=c(4,7),covars=list(cov=d$covar))
     , 4
   )
  #
  
   expect_equal(
     autodelete(count = d$count, time = d$time, changepoints=c(1,4,7),covars=list(cov=d$covar))
     ,c(1, 4)
   )

})






