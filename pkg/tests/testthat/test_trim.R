
# compare TRIM model m against trim output string to.
trimtest <- function(m, to){
  
  # data basics
  expect_equal(m$nsite, get_n_site(to))
  expect_equal(m$ntime, get_n_time(to))
    
  # time index check
  tgt <- get_time_indices(to)
  out <- index(m,"both")
  for ( i in seq_len(ncol(tgt)) ){
    expect_true( max(abs(out[,i]-tgt[,i]), na.rm=TRUE) < 1e-4
      , info=sprintf("Time index column %d",i)
    )
  }
  
  # time totals check
  tgt <- get_time_totals(to)
  out <- totals(m,"both")
  for ( i in seq_len(ncol(tgt)) ){
    expect_true( max(abs(out[,i]-tgt[,i]), na.rm=TRUE) < 1e-4
      , info=sprintf("Time totals column %d",i)
    )
  }
  
  # overall slope
  tgt <- get_overal_imputed_slope(to)
  out <- overall(m,"imputed")
  # std errors are not tested here...
  expect_true(abs(out$coef[2,1] - tgt[1]) < 1e-4)
  expect_true(abs(out$coef[2,3] - tgt[3]) < 1e-4)

  
  # but here, with a somewhat higher tolerance:
  expect_true( max( abs(out$coef[2,c(2,4)] - tgt[c(2,4)]) ) < 1e-3)

  # goodness-of-fit
  tgt <- get_gof(to)
  out <- gof(m)
  expect_equal(out$chi2$chi2, tgt$chi2$chi2, tol = 1e-3, info="chi2 value")
  expect_equal(out$chi2$df, tgt$chi2$df, info="chi2 df")
  expect_equal(out$chi2$p, tgt$chi2$p, tol= 1e-4, info="chi2 p-value") 
  expect_equal(out$LR$LR, tgt$LR$LR, tol=1e-3, info="Likelihood ratio")
  expect_equal(out$LR$df, tgt$LR$df,info="Likelihood ratio df")  
  expect_equal(abs(out$AIC), abs(tgt$AIC), tol=1e-4, info="AIC value") 
  
  # wald test
  tgt <- get_wald(to)
  out <- wald(m)
  expect_equal(out$W,tgt$W,tol=1e-1, info="Wald test value")
  expect_true(all(out$df==tgt$df), info="Wald test df")
  expect_equal(out$model,tgt$model, info="model type")
  expect_equal(out$p, tgt$p, tol=1e-4, info="Wald test p-value")
  
  # coefficients
  tgt <- get_coef(to)
  out <- coefficients(m,which="both")
  
  expect_equal(out$model, tgt$model)
  if ( nrow(tgt$coef) == 1){
    for ( i in 3:6 ){
      expect_equal(out$coef[1,i], tgt$coef[,i-2], tol=1e-4
         , info=sprintf("Coefficients column %d",i)
     )
    }
  } else {
    for ( i in 1:6 ){
      expect_equal(out$coef[,i],tgt$coef[,i],tol=1e-4
         , info=sprintf("Coefficients column %d",i))
    }
  }
}


context("TRIM: no covariates, no changepoints")

test_that("model 2",{
  tc <- read_tcf("outfiles/skylark-1d.tcf")
  m <- trim(tc)
  to <- read_tof("outfiles/skylark-1d.out")
  trimtest(m, to)
})


context("TRIM: no covariates, with changepoints")

test_that("model 2",{
  tc <- read_tcf("outfiles/skylark-1e.tcf")
  m <- trim(tc)
  to <- read_tof("outfiles/skylark-1e.out")
  trimtest(m,to)
})






context("Output printers")
test_that("S3 output printers", {
  data(skylark)
  m2 <- trim(count ~ time + site, data=skylark, model=2, overdisp=TRUE, serialcor = TRUE)
  expect_output(print(m2))
  expect_output(print(coef(m2)))
  expect_output(print(wald(m2)))
  expect_output(print(overall(m2)))
  expect_output(print(totals(m2)))
  expect_output(print(index(m2)))
})


