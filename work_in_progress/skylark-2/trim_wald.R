# \title{TRIM code documentation}
# \author{Patrick Bogaart}
# \date{\today}

# ############################################################## Wald tests ####

# ================================================================== Theory ====

# TRIM provides a number of tests for the significance of groups of parameters.
# These so called Wald-tests are based on the estimated covariance matrix of the parameters,
# and since this covariance matrix takes the overdispersion and serial
# correlation into account (if specified), these tests are valid not only if
# the counts are assumed to be independent Poisson observations but also if
# $\sigma$ and/or $\rho$ are estimated.

# The form of the Wald-statistic for testing simultaneously whether several parameters are different from zero is
# \begin{equation}
#   W = {\hat\theta}^T \left[\var{\hat\theta}\right]^{-1} \hat\theta
#   \label{Wald-vector}
# \end{equation}
# with $\hat\theta$ a vector containing the parameter estimates to be tested,
# and $\var{\hat\theta}$ the covariance matrix of $\hat\theta$.
# For the univariate case (i.e., $\theta$ is a scalar),
# Eqn~\eqref{Wald-vector} reduces to the univariate case
# \begin{equation}
#   W = {\hat\theta}^2 / \var{\hat\theta}
#   \label{Wald-scalar}
# \end{equation}

# The following Wald-tests are performed by TRIM
# \begin{itemize}
# \item Test for the significance of the slope parameter (model 2).
# \item Tests for the significance of changes in slope (model 2).
# \item Test for the significance of the deviations from a linear trend (model 3).
# \item Tests for the significance of the effect of each covariate (models 2 and 3).
# \end{itemize}
# The Wald-tests are asymptotically $\chi^2_{df}$ distributed,
# with the number of degrees of freedom equal to the rank of the covariance matrix $\var{\hat\theta}$.
# The hypothesis that the tested parameters are zero is rejected for large values of the test-statistic and small values of the associated significance probabilities (denoted by $p$),
# so parameters are significantly different from zero if $p$ is smaller than some chosen significance level (customary choices are 0.01, 0.05 and 0.10).


# ============================================================= Computation ====

wald <- function(x) UseMethod("wald")

#' Test significance of TRIM coefficients by using a Wald test
#'
#' @param x TRIM output structure (i.e., output of a call to \code{trim})
#'
#' @return a model-dependent list of Wald statistics
#' @export
#'
#' @examples
#' z2 <- trim(..., model=2)
#' print(wald(z2))  # print info on significance of slope parameters
#' z3 <- trim(..., model=3)
#' print(wald(z3))  # print info on significance of deviations from linear trend
wald.trim <- function(z)
{
  # Collect TRIM output variables that we need here.
  model  <- z$model
  ntime  <- z$ntime
  nbeta  <- length(z$beta)
  nbeta0 <- z$nbeta0
  covars <- z$covars
  ncovar <- z$ncovar
  nclass <- z$nclass
  beta <- z$beta
  var_beta <- z$var_beta

  if (model==3 && ncovar==0) {
    gstar <- z$gstar
    var_gstar <- z$var_gstar
  }

  wald <- list() # Create empty output object

  # Test for the significance of the slope parameter ---------------------------

  # This test applies to the case where model 2 is used without covariates or changepoints.
  # There thus is a single $\beta$ representing the trend for all sites and throughout the whole period,
  # and the univariate approach Eqn~\eqref{Wald-scalar} applies.
  if (model==2 && nbeta==1 && ncovar==0) {
    theta <- as.numeric(beta)
    var_theta <- as.numeric(var_beta)
    W  <- theta^2 / var_theta # Compute the Wald statistic by \eqref{Wald-scalar}
    df <- 1 # Degrees of freedom
    p  <- 1 - pchisq(W, df=df) # $p$-value, based on $W$ being $\chi^2$ distributed.
    wald$slope <- list(W=W, df=df, p=p)
  }

  # Tests for the significance of changes in slope -----------------------------

  # When model 2 is used with changepoints, $\beta$ is now a vector slope parameters,
  # and the Wald test is used to test if these slopes significantly change after a changepoint.
  # Thus, the Wald test is not applied to individual slope magnitudes $\beta_i$,
  # but on the \emph{change in} slope $\beta'_i$, where
  # $$ \beta'_i = \beta_i - \beta_{i-1} $$ and $$ \beta'_1 = \beta_1 $$
  # The vector $\beta'$ can be obtained from the linear tranformation
  # $ \beta' = A \beta$
  # where transformation matrix $A$ is a simple banded matrix structured as

  ## $$ A = \begin{pmatrix*}[r]
  ##    1 &  0 & 0 & \cdots & 0 & 0\\
  ##   -1 &  1 & 0 & \cdots & 0 & 0\\
  ##    0 & -1 & 1 & \cdots & 0 & 0\\
  ##    \vdots & \vdots & \vdots & \ddots & \vdots & \vdots\\
  ##   0 & 0 & 0 & \cdots & 1 & 0 \\
  ##   0 & 0 & 0 & \cdots & -1 & 1
  ##   \end{pmatrix*} $$

  # $$ A = \begin{pmatrix*}[r]
  #    1 &  0 & 0 & \cdots \\
  #   -1 &  1 & 0 & \cdots \\
  #    0 & -1 & 1 & \cdots \\
  #    \vdots & \vdots & \vdots & \ddots
  #   \end{pmatrix*} $$
  # that is,
  # \begin{equation}
  #   A_{i,j} = \begin{cases}
  #      1 &\quad\text{for $i=j$},\\
  #     -1 &\quad\text{for $i=j+1$},\\
  #      0 &\quad\text{otherwise}.
  #   \end{cases}
  # \end{equation}
  else if (model==2 && nbeta>1 && ncovar==0) {
    A <- diag(nbeta)           # Start with a diagonal matrix
    idx <- row(A)==(col(A)+1)  # The band just below the diagonal
    A[idx] <- -1
    dbeta <- A %*% beta
    # The covariance matrix of $\beta'$, $V^{\beta'}$, can be obtained from $V^\beta$ by
    # applying the Taylor (???) delta (???) method
    # $$ V^{\beta'} = A V^\beta A^T $$
    # and the resulting diagonal elements can be taken as the variance of the corresponding $\beta'$:
    # $$ \var{\beta'_i} = V^{\beta'}_{i,i} $$
    var_dbeta <- A %*% var_beta %*% t(A)
    # Note that the Wald test is applied to each change point individually.
    theta <- as.numeric(dbeta)
    var_theta <- diag(var_dbeta)
    W <- theta^2 / var_theta # Eqn~\eqref{Wald-scalar}
    df <- 1 # degrees of freedom
    p  <- 1 - pchisq(W, df=df) # $p$-value, based on $W$ being $\chi^2$ distributed.
    wald$dslope <- list(W=W, df=df, p=p)
  }

  # A variant is the same test for  multiple covariates
  else if (model==2 && nbeta>1 && ncovar>0) {
    # Again compute the transformation matrix
    A <-diag(nbeta0)
    idx <- row(A)==(col(A)+1) # band just below the diagonal
    A[idx] <- -1
    # Parameter vector $\beta$ and it's covariance matrix $V^\beta$ consists of `blocks'
    # representing either the baseline slopes $\beta_0$ or the impact of covariates $\beta_k$.
    # If we have $n$ changepoints, then these blocks are $n\times 1$ and $n\times n$ for
    # $\beta$ and $V^\beta$, respectively. First create an index for the first block.
    nblock = sum(nclass-1)+1
    stopifnot(nblock*nbeta0==nbeta)
    idx0 = 1:nbeta0
    # Again, $\beta'$ is easily computed from $\beta$, except that the transformation
    # $\beta' = A\beta$ is applied to each block
    dbeta = matrix(0, nbeta, 1)
    for (i in 1:nblock) {
      idx = idx0 + nbeta0*(i-1)
      dbeta[idx, ] <- A %*% beta[idx, ]
    }
    # idem for the covariance matrix
    var_dbeta <- matrix(0, nbeta, nbeta)
    for (i in 1:nblock) {
      ridx = seq(to=i*nbeta0, len=nbeta0) # ((i-1)*nbeta0+1) : (i*nbeta0)
      for (j in 1:nblock) {
        cidx <- seq(to=j*nbeta0, len=nbeta0)
        var_dbeta[ridx,cidx] <- A %*% var_beta[ridx,cidx] %*% t(A)
      }
    }

    # Compute a Wald statistic for each changepoint
    W = numeric(nbeta0)
    for (b in 1:nbeta0) {
      idx <- seq(from=b, by=nbeta0, len=nblock)
      theta = dbeta[idx]
      var_theta = var_dbeta[idx,idx]
      W[b] = t(theta) %*% solve(var_theta) %*% theta
    }
    df <- nblock
    p <- 1 - pchisq(W, df)
    wald$dslope <- list(W=W, df=df, p=p)
  }


  # Test for the significance of the deviations from a linear trend ------------

  # For Model 3, we use the Wald test to test if the residuals around the overall
  # trend (i.e., the $\gamma_j^\ast$) significantly differ from 0.
  # For this case, the vectorized Wald equation~\eqref{Wald-vector} is used.
  # Note that this test is only performed when there are no covariates.
  else if (model==3 && ncovar==0) {
    J <- ntime
    theta <- matrix(gstar) # Column vector of all $J$ $\gamma^\ast$.
    var_theta <- var_gstar[-1,-1] # Covariance matrix; drop the $\beta^\ast$ terms.

    # We now have $J$ equations, but due to the double contraints 2 of them are linear
    # dependent on the others. Let's confirm this:
    eig <- eigen(var_theta)$values
    stopifnot(sum(eig<1e-7)==2)

    # Shrink $\theta$ and it's covariance matrix to remove the dependent equations.
    theta <- theta[3:J]
    var_theta <- var_theta[3:J, 3:J]

    W <- t(theta) %*% solve(var_theta) %*% theta # Eqn~\eqref{Wald-vector}
    W <- as.numeric(W) # Convert from $1\times1$ matrix to proper atomic
    df <- J-2 # degrees of freedom
    p  <- 1 - pchisq(W, df=df) # $p$-value, based on $W$ being $\chi^2$ distributed.

    wald$deviations <- list(W=W, df=df, p=p)
  }
  else if (model==3 && ncovar>0) {
    # pass
  } else stop("Can't happen")

  # Tests for the significance of the effect of each covariate -----------------

  # As explained in ????, for both models 2 and 3, the covariate effects are modelled as additions $\beta_k$ to the
  # baseline slope parameters $\beta_0$ which represent the first class of all covariates.
  if (ncovar>0) {
    wald$covar <- data.frame(Covariate=names(covars), W=0, df=0, p=0)
    size <- (nclass-1)*nbeta0    # size of covariate block witin total $\beta$ vector
    last <- cumsum(size)+nbeta0  # last element of covariate block
    first <- last - size + 1     # first element
    for (i in 1:length(covars)) {
      idx <- first[i] : last[i]
      theta <- matrix(beta[idx])
      var_theta <- var_beta[idx, idx]
      W = t(theta) %*% solve(var_theta) %*% theta
      wald$covar$W[i]     <- W
      wald$covar$df[i]    <- nbeta0
      wald$covar$p[i]     <- 1 - pchisq(W, df=nbeta0)
    }
  }

  # Output results in a list with type 'trim.wald' to enable specialized further processing.
  class(wald) <- "trim.wald"
  wald
}


# ================================================================ Printing ====

# A simple printing function is provided that mimics the output of TRIM for Windows.

print.trim.wald <- function(w) {
  stopifnot(class(w)=="trim.wald")

  if (!is.null(w$covar)) {
    printf("Wald test for significance of covariates\n")
    print(w$covar, row.names=FALSE)
    printf("\n")
  }

  if (!is.null(w$slope)) {
    printf("Wald test for significance of slope parameter\n")
    printf("  Wald = %.2f, df=%d, p=%f\n", w$slope$W, w$slope$df, w$slope$p)
  } else if (!is.null(w$dslope)) {
    printf("Wald test for significance of changes in slope\n")
    df = data.frame(Changepoint = 1:length(w$dslope$W),
                    Wald_test = w$dslope$W, df = w$dslope$df, p = w$dslope$p)
    print(df, row.names=FALSE)
  } else if (!is.null(w$deviations)) {
    printf("Wald test for significance of deviations from linear trend\n")
    printf("  Wald = %.2f, df=%d, p=%f\n", w$deviations$W, w$deviations$df, w$deviations$p)
  }
}
