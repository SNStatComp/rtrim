# ##################################################### Stepwise refinement ####

# =============================================== Main refinement prodecure ====

#' TRIM stepwise refinement
#'
#' @param count a numerical vector of count data.
#' @param time.id a numerical vector time points for each count data point.
#' @param site.id a numerical vector time points for each count data point.
#' @param covars an optional list of covariates
#' @param model a model type selector
#' @param serialcor a flag indication of autocorrelation has to be taken into account.
#' @param overdisp a flag indicating of overdispersion has to be taken into account.
#' @param changepoints a numerical vector change points (only for Model 2)
#'
#' @return a list of class \code{trim}, that contains all output, statistiscs, etc.
#'   Usually this information is retrieved by a set of postprocessing functions
#'
#' @keywords internal
trim_refine <- function(count, time.id, site.id, covars=list(),
                        model=2, serialcor=FALSE, overdisp=FALSE,
                        changepoints=integer(0),weights=numeric(0))
{
  org_cp = changepoints
  ncp <- length(org_cp)
  active <- rep(TRUE, ncp) # Keeps track which original changepoints are active or not

  # Always start with an estimation using all proposed changepoints
  cur_cp <- org_cp
  z <- trim_workhorse(count, time.id, site.id, covars, model, serialcor, overdisp, cur_cp, weights)

  # # Hack: remove all except the first changepoints
  # n <- length(org_cp)
  # active[2:n] <- FALSE
  # cur_cp <- org_cp[active]
  # z <- trim_workhorse(count, time.id, site.id, covars, model, serialcor, overdisp, cur_cp)

  for (iter in 1:100) {
    # Phase 1: can one of the changepoints be removed?
    # (Only applies for 2 or more changepoints)
    if (sum(active)>=2) {
      W <- wald(z)
      max_p = max(W$dslope$p)
      if (max_p > 0.2) {
        i = which.max(W$dslope$p)
        del_cp <- cur_cp[i]
        del_p  <- max_p
        rprintf("\n>>> Deleted changepoint %d (p = %.4f) <<<\n\n", del_cp, del_p)
        # remove from original changepoints
        i = which(org_cp == del_cp)
        active[i] = FALSE
        removed <- TRUE
      } else removed <- FALSE
    } else removed <- FALSE

    # If a changepoint has been removed, we'll need to re-estimate the model
    if (removed) {
      cur_cp = org_cp[active]
      z <- trim_workhorse(count, time.id, site.id, covars, model, serialcor, overdisp, cur_cp, weights)
    }

    # Phase 2: try to re-insert previously removed changepoints
    alpha <- z$alpha
    beta  <- z$beta
    nacp <- sum(active) # Number of active changpoints
    p <- numeric(ncp)
    for (i in 1:ncp) if (active[i]==FALSE) {
      # deleted changepoints
      num.active.before = ifelse(i==1, 0, sum(active[1:(i-1)]))
      num.active.after  = ifelse(i==ncp, 0, sum(active[(i+1):ncp]))
      # cast beta in a matrix with beta0 in first column, covars in other columns
      beta_t <- matrix(as.vector(beta), nacp)
      if (num.active.before==0) {
        stop("Should never happen")
      } else if (num.active.after==0) {
        beta_last = beta_t[num.active.before, ,drop=FALSE]
        beta_t = rbind(beta_t, beta_last)
      } else {
        beta1 = beta[1:num.active.before, ,drop=FALSE]
        beta_last = beta_t[num.active.before, ,drop=FALSE]
        beta2 = beta[(nacp-num.active.after+1):nacp, ,drop=FALSE]
        beta_t = rbind(beta1, beta_last, beta2)
      }
      beta_t <- matrix(beta_t) # Reshape into column vector
      active_t <- active # Create list of current (test) changepoints
      active_t[i] <- TRUE
      cp_t <- org_cp[active_t]
      idx <- c(rep(FALSE,num.active.before), TRUE, rep(FALSE, num.active.after))
      p[i] <- Score(z, alpha, beta_t, cp_t, idx)
    } else p[i] = 1.0 # active changepoints

    # A changepoint is re-inserted if the minimum signficance is lower than a
    # specified threshold
    min_p <- min(p)
    if (min_p < 0.15) {
      i = which.min(p)
      active[i] <- TRUE
      ins_cp <- org_cp[i]
      rprintf("\n>>> Re-inserted changepoint %d (p=%.4f) <<<\n\n", ins_cp, min_p)
      added <- TRUE
    } else added <- FALSE

    # If a changepoint has been re-inserted, we'll need to re-estimate the model
    if (added) {
      cur_cp = org_cp[active]
      z <- trim_workhorse(count, time.id, site.id, covars, model, serialcor, overdisp, cur_cp, weights)
    }

    # Finished refinement?
    if (removed==FALSE & added==FALSE) {
      rprintf("Stepwise refinement ready.\n")
      break
    }
  }

  # Return last estimated model
  z
}

# ======================================================= Score computation ====

Score <- function(z, alpha, beta, changepoints, index) {
  # Unpack TRIM variables
  f  <- z$f
  rho <- z$rho
  sig2 <- z$sig2
  covars <- z$covars
  cvmat <- z$cvmat
  nsite <- z$nsite
  ntime <- z$ntime
  if (is.null(z$wt)) {
    wt <- matrix(1.0, nsite, ntime)
  } else {
    wt <- z$wt
  }

  ncp    <- length(changepoints)
  ncovar <- length(covars)

  # Define covar aux info
  if (ncovar>0) {
    nclass <- numeric(ncovar)
    for (i in 1:ncovar) {
      if (is.factor(covars[[i]])) {
        nclass[i] <- nlevels(covars[[i]] )
      } else {
        nclass[i] <- max(covars[[i]])
      }
    }
  }

  # Define B0, the global B for model 2
  B0 <- matrix(0, ntime, ncp)
  for (i in 1:ncp) {
    cp1  <-  changepoints[i]
    cp2  <-  ifelse(i<ncp, changepoints[i+1], ntime)
    if (cp1>1) B0[1:(cp1-1), i] <- 0
    B0[cp1:cp2, i] <- 0:(cp2-cp1)
    if (cp2<ntime) B0[(cp2+1):ntime,i] <- B0[cp2,i]
  }

  # Prepare for covariate contributions to $B$
  if (ncovar>0) {
    cvmask <- list()
    for (cv in 1:ncovar) {
      cvmask[[cv]] = list()
      for (cls in 2:nclass[cv]) {
        cvmask[[cv]][[cls]] <- list()
        for (i in 1:nsite) {
          cvmask[[cv]][[cls]][[i]] <- which(cvmat[[cv]][i, ]!=cls)
        }
      }
    }
  }

  # Global R
  if (rho==0.0) {
    Rg <- diag(1, ntime) # default (no autocorrelation) value
  } else {
    Rg <- rho ^ abs(row(diag(ntime)) - col(diag(ntime)))
  }

  # Compute score matrix
  U_b <- 0
  i_b <- 0
  for (i in 1:nsite) {
    # Select observations
    f_i <- f[i, ]
    observed <- is.finite(f_i)
    f_i <- f_i[observed]

    # Define B
    B <- B0
    if (ncovar>0) { # add a copy of $B$ for each covar class
      for (cv in ncovar) {
        for (cls in 2:nclass[cv]) {
          Btmp <- B0
          mask <- cvmask[[cv]][[cls]][[i]]
          if (length(mask)>0) Btmp[mask, ] = 0
          B <- cbind(B, Btmp)
        }
      }
    }
    # Compute mu
    mu = exp(alpha[i] + B %*% beta) / wt[i, ]
    mu_i = mu[observed]

    d_mu_i <- diag(mu_i, length(mu_i)) # Length argument guarantees diag creation

    # Compute $V$ and $\Omega$
    if (rho==0.0 && sig2==1.0) { # ML
      V_i <- sig2 * d_mu_i
    } else { # GEE
      idx <- which(observed)
      R_i <- Rg[idx,idx]
      V_i <- sig2 * sqrt(d_mu_i) %*% R_i %*% sqrt(d_mu_i)
    }
    V_inv <- solve(V_i)

    Omega <- d_mu_i %*% V_inv %*% d_mu_i

    B_i = B[observed, ,drop=FALSE] # recyle index for covariates
    U_b <- U_b + t(B_i) %*% d_mu_i %*% V_inv %*% (f_i - mu_i)

    nobs_i = length(f_i)
    ones <- matrix(1, nobs_i, 1)
    d_i <- as.numeric(t(ones) %*% Omega %*% ones)
    i_b <- i_b - t(B_i) %*% (Omega - (Omega %*% ones %*% t(ones) %*% Omega) / d_i) %*% B_i
  }
  V <- solve(-i_b)
  S <- t(U_b) %*% V %*% U_b

  df <- length(beta) / length(changepoints) # Number of beta-blocks (baseline+covariates)
  p <- 1 - pchisq(S, df=df)
  p
}
