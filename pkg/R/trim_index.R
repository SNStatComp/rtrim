# ################################################################# Indices ####


# ============================================= Internal workhorse function ====

.index <- function(tt, var_tt, b) {

  # Time index $\tau_j$ is defined as time totals $\Mu$, normalized by the time total for the base
  # year, i.e.\,
  # $$ \tau_j = \Mu_j / \Mu_b $$
  # where $b\in[1\ldots J]$ indicates the base year.
  tau <- tt / tt[b]

  # Uncertainty is again quantified as a standard error $\sqrt{\var{\cdot}}$,
  # approximated using the delta method, now extended for the multivariate case:
  # \begin{equation}
  #   \var{\tau_j} = \var{f(\Mu_b,\Mu_j)} = d^T V(\Mu_b,\Mu_j) d \label{var_tau}
  # \end{equation}
  # where $d$ is a vector containing the partial derivatives of $f(\Mu_b,\Mu_j)$
  # \begin{equation}
  #   d = \begin{pmatrix} -\Mu_j \Mu_b^{-2} \\ \Mu_b^{-1} \end{pmatrix}
  # \end{equation}
  # and $V$ the covariance matrix of $\Mu_b$ and $\Mu_j$:
  # \begin{equation}
  #   V(\Mu_b,\Mu_j) = \begin{pmatrix*}[l]
  #     \var{\Mu_b} & \cov{\Mu_b, \Mu_j} \\
  #     \cov{\Mu_b, \Mu_j} & \var{\Mu_j}
  #   \end{pmatrix*}
  # \end{equation}
  # Note that for the base year $b$, where $\tau_b\equiv1$, Eqn~\eqref{var_tau} results in
  # $\var{\tau_b}=0$, which is also expected conceptually because $\tau_b$ is not an estimate but an exact and fixed result.
  J <- length(tt)
  var_tau <- numeric(J)
  for (j in 1:J) {
    d <- matrix(c(-tt[j] / tt[b]^2, 1/tt[b]))
    V <- var_tt[c(b,j), c(b,j)]
    var_tau[j] <- t(d) %*% V %*% d
  }

  out = list(tau=tau, var_tau=var_tau)
}



# ========================================================== User interface ====

#' Extract time-indices from TRIM output
#'
#' @param x TRIM output structure (output of a call to \code{\link{trim}})
#' @param base Base time point, for which the index is 1
#' @param which Selector to distinguish between time indices based on the imputed data (default),
#' the modelled data, or both.
#'
#' @return a data frame containing indices and their uncertainty (expressed as standard error)
#' @export
#'
#' @examples
#' z <- trim(tcf,dat);
#' index(z) # prints the indices for the imputed data
#' print(index(z,which="imputed")) # idem
#' print(index(z,4)) # using the 4th time point as reference
#' index(z, which="both") # mimics classic TRIM
#' SE <- index(z)$se_imp # Extract standard error for the imputed data
index <- function(trm, base=1, which=c("imputed","model","both")) {
  stopifnot(class(trm)=="trim")

  # Computation and output is user-configurable
  which <- match.arg(which)
  if (which=="model") {
    # Call workhorse function to do the actual computation
    mod <- .index(trm$tt_mod, trm$var_tt_mod, base)
    # Store results in a data frame
    out = data.frame(time  = 1:trm$ntime,
                     model = mod$tau,
                     se_mod = sqrt(mod$var_tau))
  } else if (which=="imputed") {
    # Idem, using the imputed time totals instead
    imp <- .index(trm$tt_imp, trm$var_tt_imp, base)
    out = data.frame(time    = 1:trm$ntime,
                     imputed = imp$tau,
                     se_imp  = sqrt(imp$var_tau))
  } else if (which=="both") {
    # Idem, using both modelled and imputed time totals.
    mod <- .index(trm$tt_mod, trm$var_tt_mod, base)
    imp <- .index(trm$tt_imp, trm$var_tt_imp, base)
    out = data.frame(time    = 1:trm$ntime,
                     model   = mod$tau,
                     se_mod  = sqrt(mod$var_tau),
                     imputed = imp$tau,
                     se_imp  = sqrt(imp$var_tau))
  } else stop("Can't happen") # because other cases are catched by match.arg()

  out
}


