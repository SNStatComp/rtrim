context("TRIM models without covariates")

test_that("skylark-1d model",{
  tc <- read_tcf("outfiles/skylark-1d.tcf")
  m <- trim(tc)
  to <- read_tof("outfiles/skylark-1d.out")
  
  # time index check
  tgt <- get_time_indices(to)
  out <- index(m,"both")
  for ( i in seq_len(ncol(tgt)) ){
    expect_true( max(abs(out[,i]-tgt[,i]), na.rm=TRUE) < 1e-4
      , info=sprintf("Time index column %d",i)
    )
  }
  
  tgt <- get_time_totals(to)
  out <- totals(m,"both")
  for ( i in seq_len(ncol(tgt)) ){
    expect_true( max(abs(out[,i]-tgt[,i]), na.rm=TRUE) < 1e-4
      , info=sprintf("Time totals column %d",i)
    )
  }
  
  tgt <- get_overal_imputed_slope(to)
  out <- overall(m,"imputed")
  # std errors are not tested here...
  expect_true( max( abs(out$coef[2,c(1,3)] - tgt[c(1,3)]) ) < 1e-4)
  # but here, with a somewhat higher tolerance:
  expect_true( max( abs(out$coef[2,c(2,4)] - tgt[c(2,4)]) ) < 1e-3)

  tgt <- get_gof(to)
  out <- gof(m)
  expect_true(abs(out$chi2$chi2 - tgt[[1]][1]) < 1e-3)
  expect_true(out$chi2$df == tgt[[1]][2])
  expect_true(abs(out$chi2$p - tgt[[1]][3]) < 1e-5)    
  
  
  
})


