# #################################################### Parameter estimation ####

# This Section describes the core TRIM function, which estimates the TRIM parameters.

alpha_method <- 1
graph_debug <- FALSE
compatible <- FALSE

# ##################################################### Estimation function ####

#' TRIM estimation function
#'
#' @param count a numerical vector of count data.
#' @param site.id a numerical vector time points for each count data point.
#' @param time.id a numerical vector time points for each count data point.
#' @param covars an optional data frame with covariates
#' @param model a model type selector (1, 2 or 3)
#' @param serialcor a flag indication of autocorrelation has to be taken into account.
#' @param overdisp a flag indicating of overdispersion has to be taken into account.
#' @param changepoints a numerical vector change points (only for Model 2)
#' @param stepwise a flag indicating stepwise refinement of changepoints is to be used.
#' @param autodelete a flag indicating auto-deletion of changepoints with too little observations.
#' @param weights a numerical vector of weights.
#'
#' @return a list of class \code{trim}, that contains all output, statistiscs, etc.
#'   Usually this information is retrieved by a set of postprocessing functions
#'
#' @keywords internal
trim_estimate <- function(count, site.id, time.id, covars=data.frame()
                         , model=2, serialcor=FALSE, overdisp=FALSE
                         , changepoints=integer(0)
                         , autodelete=TRUE, weights=numeric(0)
                         , stepwise=FALSE, covin=list()
                         , ...)
{
  call <- sys.call(1)

  # kick out missing/zero sites
  useful <- count>0
  ok = rep(TRUE, length(count))
  sites = unique(site.id)
  del.sites = character(0)
  nkickout = 0
  for (site in sites) {
    idx = site.id==site
    if (!any(count[idx]>0, na.rm=TRUE)) {
      ok[idx] = FALSE
      nkickout = nkickout+1
      del.sites <- c(del.sites, site)
    }
  }
  if (nkickout>0) {
    count = count[ok]
    site.id = site.id[ok]
    time.id = time.id[ok]
    if (length(weights)>0) weights = weights[ok]
    rprintf("Removed %d %s without observations: (%s)\n", nkickout,
            ifelse(nkickout==1, "site","sites"), paste0(del.sites, collapse=", "))
  }

  # Handle "auto" changepoints
  if (is.character(changepoints)) {
    if (changepoints %in% c("all","auto")) {
      if (changepoints == "auto") stepwise=TRUE
      J <- length(unique(time.id))
      changepoints <- 1 : (J-1)
    }
  }

  if (isTRUE(stepwise) && model != 2){
    stop(sprintf("Stepwise removal only works for model 2"), call.=FALSE)
  }

  t1 <- Sys.time()
  if (isTRUE(stepwise)) {
    m <- trim_refine(count, site.id, time.id, covars, model, serialcor
          , overdisp, changepoints, weights)
  } else {
    # data input checks: throw error if not enough counts available.
    if (model == 2 && length(changepoints)>0 && autodelete){
      changepoints <- autodelete(count=count, time=time.id
        , changepoints = changepoints, covars=covars)
    } else if (model == 2){
      assert_plt_model(count = count, time = time.id
              , changepoints = changepoints, covars = covars)

    } else if (model == 3){
      assert_sufficient_counts(count = count, index = time.id)
      assert_covariate_counts(count = count, time = time.id, covars=covars)
    }

    # compute actual model
    m <- trim_workhorse(count, site.id, time.id, covars, model
          , serialcor, overdisp, changepoints, weights, covin
          , ...)
  }
  t2 <- Sys.time()
  m$dt <- difftime(t2,t1)
  rprintf("Running trim took %8.4f %s\n",m$dt,attr(m$dt,"units"))
  m$call <- call
  m
}


# ###################################################### Workhorse function ####

#' TRIM workhorse function
#'
#' @param count a numerical vector of count data.
#' @param site.id a numerical vector time points for each count data point.
#' @param time.id an numerical vector time points for each count data point.
#' @param covars an optional data frame with covariates
#' @param model a model type selector
#' @param serialcor a flag indication of autocorrelation has to be taken into account.
#' @param overdisp a flag indicating of overdispersion has to be taken into account.
#' @param changepoints a numerical vector change points (only for Model 2)
#' @param weights a numerical vector of weights.
#' @param covin ...
#' @param conv_crit convergence criterion.
#' @param max_iter maximum number of iterations allowed.
#' @param soft: specifies if trim on error returns an error code (TRUE) or just stops with a message (FALSE)
#'
#' @return a list of class \code{trim}, that contains all output, statistiscs, etc.
#'   Usually this information is retrieved by a set of postprocessing functions
#'
#'
#' @keywords internal
trim_workhorse <- function(count, site.id, time.id, covars=data.frame(),
                         model=2, serialcor=FALSE, overdisp=FALSE,
                         changepoints=integer(0), weights=numeric(0),
                         covin = list(),
                         conv_crit=1e-5, max_iter=200, max_sub_step=7,
                         soft=FALSE, debug=FALSE)
{

  # =========================================================== Preparation ====

  # Check the arguments. \verb!count! should be a vector of numerics.
  stopifnot(class(count) %in% c("integer","numeric"))
  n = length(count)

  # \verb!site.id! should be a vector of numbers, strings or factors
  stopifnot(class(site.id) %in% c("integer","character","factor"))
  stopifnot(length(site.id)==n)

  # \verb!time.id! should be an ordered factor, or a vector of consecutive years or numbers
  # Note the use of "any" because of multiple classes for ordered factors
  stopifnot(any(class(time.id) %in% c("integer","numeric","factor")))
  if (any(class(time.id) %in% c("integer","numeric"))) {
    check = unique(diff(sort(unique(time.id))))
    stopifnot(check==1 && length(check)==1)
  }
  stopifnot(length(time.id)==n)
  # Convert the time points to a factor

  # \verb!covars! should be a list where each element (if any) is a vector
  stopifnot(class(covars)=="data.frame")
  ncovar = length(covars)
  use.covars <- ncovar>0
  if (use.covars) {
    # convert to numerical values
    icovars <- vector("list", ncovar)
    for (i in 1:ncovar) if (any(is.na(covars[[i]])))
      stop(sprintf('NA values not allowed for covariate "%s".', names(covars)[i]), call.=FALSE)
    for (i in 1:ncovar) icovars[[i]] = as.integer((covars[[i]]))

    # Also, each covariate $i$ should be a number (ID) ranging $1\ldots nclass_i$
    nclass <- numeric(ncovar)
    for (i in 1:ncovar) {
      cv <- icovars[[i]] # The vector of covariate class ID's
      stopifnot(min(cv)==1) # Assert lower end of range
      nclass[i] = max(cv)  # Upper end of range
      #stopifnot(nclass[i]>1) # Assert upper end
      stopifnot(length(unique(cv))==nclass[i]) # Assert the range is contiguous
    }
  } else {
    nclass <- 0
  }

  # remove covariates that have only a single class
  while (any(nclass==1)) {
    idx <-  which(nclass==1)[1]
    warning(sprintf("Removing covariate \"%s\" which has only one class.", names(covars)[idx])
            , call. = FALSE)
    covars <-  covars[-idx]
    icovars <- icovars[-idx]
    nclass <- nclass[-idx]
    ncovar <- ncovar-1
    if (ncovar==0) {
      warning("No covariates left", call. = FALSE)
      use.covars <- FALSE
    }
  }

  # \verb!model! should be in the range 1 to 3
  stopifnot(model %in% 1:3)

  # Weights should be either absent, or aligned with the counts
  if (length(weights)>0) {
    use.weights <- TRUE
    stopifnot(length(weights)==length(count))
  } else use.weights <- FALSE

  # User-specified covariance
  use.covin <- length(covin)>0

  # Convert site and time to factors, if they're not yet
  if (any(class(time.id) %in% c("integer","numeric"))) time.id <- ordered(time.id)
  ntime = length(levels(time.id))

  #org.site.id <- site.id # Remember the original values for output purposes.
  if (class(site.id) %in% c("integer","numeric")) site.id <- factor(site.id)
  nsite = length(levels(site.id))

  # Create observation matrix $f$.
  # Convert the data from a vector representation to a matrix representation.
  # It's OK to have missing site/time combinations; these will automatically
  # translate to NA values.
  f <- matrix(0, nsite, ntime) # ??? check if we should not use NA instead of 0!!!
  rows <- as.integer(site.id) # `site.id' is a factor, thus this results in $1\ldots I$.
  cols <- as.integer(time.id) # idem, $1 \ldots J$.
  idx <- (cols-1)*nsite+rows   # Create column-major linear index from row/column subscripts.
  f[idx] <- count    # ... such that we can paste all data into the right positions
  ff = f # allow browser() inspection

  # TRIM is not intended for extrapolation. Therefore, issue a warning if the first or last
  # time points do not contain positive observations.
  totals <- colSums(f, na.rm=TRUE)
  if (sum(totals)==0) stop("No positive observations in the data.")
  if (totals[1]==0) {
    n = which(totals>0)[1] - 1
    warning(sprintf("Data starts with %d years without positive observations.", n))
  }
  totals <- rev(totals)
  if (totals[1]==0) {
    n = which(totals>0)[1] - 1
    warning(sprintf("Data ends with %d years without positive observations.", n))
  }

  # Create similar matrices for all covariates
  if (use.covars) {
    cvmat <- list()
    for (i in 1:ncovar) {
      cv = icovars[[i]]
      m <- matrix(NA, nsite, ntime)
      m[idx] <- cv
      if (any(is.na(m))) stop(sprintf('(implicit) NA values in covariate "%s".', names(covars)[i]), call.=FALSE)
      cvmat[[i]] <- m
    }
  } else {
    cvmat <- NULL
  }

  # idem for the weights
  wt <- matrix(1.0, nsite, ntime)
  if (use.weights) wt[idx] <- weights

  # We often need some specific subset of the data, e.g.\ all observations for site 3.
  # These are conveniently found by combining the following indices:
  observed <- is.finite(f)  # Flags observed (TRUE) / missing (FALSE) data
  positive <- is.finite(f) & f > 0.0 # Flags useful ($f_{i,j}>0$) observations
  site <- as.vector(row(f)) # Internal site identifiers are the row numbers of the original matrix.
  time <- as.vector(col(f)) # Idem for time points.
  nobs <- rowSums(observed) # Number of actual observations per site
  npos <- rowSums(positive) # Number of useful ($f_{i,j}>0$) observations per site.

  # Check if the covin matrices all have the right size.
  if (use.covin) {
    for (i in 1:nsite) {
      stopifnot(nrow(covin[[i]])==nobs[i])
      stopifnot(ncol(covin[[i]])==nobs[i])
    }
  }

  # For model 2, we do not allow for changepoints $<1$ or $\geq J$.
  # At the same time, a changepoint $1$ must be present
  if (model==2) {
    if (length(changepoints)==0) {
      use.changepoints <- FALSE # Pretend we're not using changepoints at all
      changepoints <- 1L        # but internally use them nevertheless
    } else {
      use.changepoints <- TRUE
      years <- as.integer(levels(time.id))
      if (all(changepoints %in% years)) {
        # Convert changepoints in years  to 1..J
        changepoints <- match(changepoints, years)
      }
      stopifnot(all(changepoints>=1L))
      stopifnot(all(changepoints<ntime))
      stopifnot(all(diff(changepoints)>0))
    }
  }

  if (model!=2 && length(changepoints) > 0) stop(sprintf("Changepoints cannot be specified for model %d", model))

  # We make use of the generic model structure
  # $$ \log\mu = A\alpha + B\beta $$
  # where design matrices $A$ and $B$ both have $IJ$ rows.
  # For efficiency reasons the model estimation algorithm works on a per-site basis.
  # There is thus no need to store these full matrices. Instead, $B$ is constructed as a
  # smaller matrix that is valid for any site, and $A$ is not used at all.

  # Create matrix $B$, which is model-dependent.
  if (model==1) {
    # Model 1 has not really any beta's, but the code runs easier if we have a fake beta
    # with a fioxed value of 1. Therefore, create a corresponding B
    J <- ntime
    B <- matrix(0, J, 1)
  } else if (model==2) {
    ncp  <-  length(changepoints)
    J <- ntime
    B <- matrix(0, J, ncp)
    for (i in 1:ncp) {
      cp1  <-  changepoints[i]
      cp2  <-  ifelse(i<ncp, changepoints[i+1], J)
      if (cp1>1) B[1:(cp1-1), i]  <-  0
      B[cp1:cp2,i]  <-  0:(cp2-cp1)
      if (cp2<J) B[(cp2+1):J,i]  <-  B[cp2,i]
    }
  } else if (model==3) {
    # Model 3 in it's canonical form uses a single time parameter $\gamma$ per time step,
    # so design matrix $B$ is essentially a $J\times$J identity matrix.
    # Note, hoewever, that by definition $\gamma_1=0$, so effectively there are $J-1$ $\gamma$-values to consider.
    # As a consequence, the first column is deleted.
    B <- diag(ntime) # Construct $J\times$J identity matrix
    B <- B[ ,-1]     # Remove first column
  }

  # For some purposes (e.g. vcov() ), we do need a dummy first column in B
  Bfull <- cbind(0, B)

  # optionally add covariates. Each covar class (except class 1) adds an extra copy of $B$,
  # where rows are cleared if that site/time combi does not participate in
  if (use.covars) {
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
    # The amount of extra parameter sets is the total amount of covariate classes
    # minus the number of covariates (because class 1 does not add extra params)
    num.extra.beta.sets <- sum(nclass-1)
  }

  # When we use covariates, $B$ is site-specific. We thus define a function to make
  # the proper $B$ for each site $i$

  B0 <- B # The "standard" B
  rm(B)
  make.B <- function(i, debug=FALSE) {
    if (debug) { printf("make.B(%i): B0:", i); str(B0) }
    if (use.covars) {
      # Model 2 with covariates. Add a copy of B for each covar class
      Bfinal <- B0
      for (cv in 1:ncovar) {
        if (debug) printf("adding covar %d\n", cv)
        for (cls in 2:nclass[cv]) {
          if (debug) printf("adding class %d\n", cls)
          Btmp <- B0
          mask <- cvmask[[cv]][[cls]][[i]]
          if (length(mask)>0) Btmp[mask, ] = 0
          Bfinal <- cbind(Bfinal, Btmp)
        }
      }
    } else {
      Bfinal = B0
    }
    if (debug) { printf("make.B(%i): Bfinal:", i); str(Bfinal)}
    Bfinal
  }

  # ================================== Setup parameters and state variables ====

  # Parameter $\alpha$ has a unique value for each site.
  # alpha <- matrix(0, nsite,1) # Store as column vector
  alpha <- matrix(log(rowSums(f, na.rm=TRUE)/ntime))
  # alpha <- matrix(log(rowMeans(f*wt, na.rm=TRUE)));

  # Parameter $\beta$ is model dependent.
  if (model==1) {
    nbeta <- 1
  } else if (model==2) {
    # For model 2 we have one $\beta$ per change points
    nbeta <- length(changepoints)
  } else if (model==3) {
    # For model 3, we have one $\beta$ per time $j>1$
    nbeta = ntime-1
  }
  # If we have covariates, $\beta$'s are repeated for each covariate class $>1$.
  nbeta0 <- nbeta # Number of `baseline' (i.e., without covariates) $\beta$'s.
  if (use.covars) {
    nbeta <- nbeta0 * (sum(nclass-1)+1)
  }

  # Now that we know much parameters are requested, check if we do have enough observations.
  # For this, we ignore the `0' observations
  nalpha <- length(alpha)
  if (sum(npos) < (nalpha+nbeta)) {
    msg <- sprintf("Not enough positive observations (%d) to specify %d parameters (%d alpha + %d beta)", sum(npos), nalpha+nbeta, nalpha, nbeta)
    if (soft) return(list(error=msg)) else stop(msg, call.=FALSE)
  }

  # All $\beta_j$ are initialized at 0, to reflect no trend (model 1 or 2) or no time effects (model 3)
  beta <- matrix(0, nbeta,1) # Store as column vector

  # Variable $\mu$ holds the estimated counts.
  mu <- matrix(0, nsite, ntime)

  # Setup error handling
  err.out <- NULL

  # ====================================================== Model estimation ====

  # TRIM estimates the model parameters $\alpha$ and $\beta$ in an iterative fashion,
  # so separate functions are defined for the updates of these and other variables needed.

  #3 Site-parameters $\alpha$
  # Update $\alpha_i$ using:
  # $$ \alpha_i^t = \log z_i' f_i - \log z_i' \exp(B_i \beta^{t-1}) $$
  # where vector $z$ contains just ones if autocorrelation and overdispersion are
  # ignored (i.e., Maximum Likelihood, ML),
  # or weights, when these are taken into account (i.e., Generalized Estimating
  # Equations, GEE).
  # In this case,
  # $$ z = \mu V^{-1} $$
  # with $V$ a covariance matrix (see Section~\ref{covariance}).

  update_alpha <- function(method=c("ML","GEE")) {
    for (i in 1:nsite) {
      obs = observed[i, ]
      if (alpha_method==1) {
        B = make.B(i)
        f_i <- matrix(f[i, obs]) # cast as column vector
        wt_i <- matrix(wt[i, obs]) # caset as column vector
        B_i <- B[obs, , drop=FALSE]
        if (method=="ML") { # no covariance; $V_i = \diag{mu}$
          z_t <- matrix(1, 1, nobs[i])
        } else if (method=="GEE") { # Use covariance
          mu_i = mu[i, obs]
          z_t <- mu_i %*% V_inv[[i]] # define correlation weights
        } else stop("Can't happen")
        if (z_t %*% f_i > 0) { # Application of method 1 is possible
          alpha[i] <<- log(z_t %*% f_i) - log(z_t %*% exp(B_i %*% beta - log(wt_i)))
        } else { # Fall back to method 2
          sumf <- sum(f_i)
          sumu <- sum(mu_i)
          dalpha <- if (sumf/sumu > 1e-7) log(sumf/sumu) else 0.0
          alpha[i] <<- alpha[i] + dalpha;
        }
        #if (z_t %*% f_i < 0) z_t = matrix(1, 1, nobs[i]) # alternative hack
      } else { #method 2: classic TRIM
        f_i <- f[i, obs]
        mu_i <- mu[i, obs]
        sumf <- sum(f_i)
        sumu <- sum(mu_i)
        dalpha <- if (sumf/sumu > 1e-7) log(sumf/sumu) else 0.0
        alpha[i] <<- alpha[i] + dalpha;
      }
    }
    #printf("\n\n** %f ** \n\n", max(abs(alpha1-alpha2)))
    #alpha <<- alpha1 # works better in some cases, i.e. testset 10104_0.tcf
    #if (!all(is.finite(alpha))) alpha <- alpha2
    if (any(!is.finite(alpha))) stop("non-finite alpha problem")
  }

  # ----------------------------------------------- Time parameters $\beta$ ----

  # Estimates for parameters $\beta$ are improved by computing a change in $\beta$ and
  # adding that to the previous values:
  # $$ \beta^t = \beta^{t-1} - (i_b)^{-1} U_b^\ast \label{beta}$$
  # where $i_b$ is a derivative matrix (see Section~\ref{Hessian})
  # and $U_b^\ast$ is a Fisher Scoring matrix (see Section~\ref{Scoring}).
  # Note that the `improvement' as defined by \eqref{beta} can actually results in a decrease in model fit.
  # These cases are identified by measuring the model Likelihood Ratio (Eqn~\eqref{LR}).
  # If this measure increases, then smaller adjustment steps are applied.
  # This process is repeated until an actually improvement is found.

  max_dbeta <- 0.0

  update_beta <- function(method=c("ML","GEE"))
  {
    # Compute the proposed change in $\beta$.
    dbeta  <-  -solve(i_b) %*% U_b
    max_dbeta <<- max(abs(dbeta))
    if (!all(is.finite(dbeta))) stop("non-finite beta problem")
    # This is the maximum update; if it results in an \emph{increased} likelihood ratio,
    # then we have to take smaller steps. First record the original state and likelihood.
    beta0  <- beta
    lik0  <- likelihood()
    stepsize <- 1.0
    for (subiter in 1:max_sub_step) {
      beta  <<- beta0 + stepsize*dbeta
      if (any(!is.finite(beta))) stop("non-finite beta problem")
      if (any(beta >  7)) stop("Model non-estimable due to excessive high beta values", call.=FALSE)
      if (any(beta < -7)) stop("Model non-estimable due to excessive small beta values", call.=FALSE)

      update_mu(fill=FALSE)
      update_alpha(method)
      if (!is.null(err.out)) return()
      update_mu(fill=FALSE)

      lik <- likelihood()
      likc <- (1+conv_crit) * lik0 # threshold value
      if (lik>0 && lik < likc) break else stepsize <- stepsize / 2 # Stop or try again
      # condition lik>0 added to prevent crash due to extreme lik(mu(dbeta(i_b))) as in testcase 412_21
    }
    subiter
  }

  # --------------------- Covariance and autocorrelation \label{covariance} ----

  # Covariance matrix $V_i$ is defined by
  # \begin{equation}
  #   V_i = \sigma^2 \sqrt{\diag{\mu}} R \sqrt{\diag{\mu}} \label{V1}
  # \end{equation}
  # where $\sigma^2$ is a dispersion parameter (Section~\ref{sig2})
  # and $R$ is an (auto)correlation matrix.
  # Both of these two elements are optional.
  # If the counts are perfectly Possion distributed, $\sigma^2=1$,
  # and if autocorrelation is disabled (i.e.\ counts are independent),
  # Eqn~\eqref{V1} reduces to
  # \begin{equation}
  #   V_i = \sigma^2 \diag{\mu} \label{V2}
  # \end{equation}
  V_inv  <- vector("list", nsite) # Create storage space for $V_i^{-1}$.
  Omega <- vector("list", nsite)

  update_V <- function(method=c("ML","GEE")) {
    for (i in 1:nsite) {
      obs <- observed[i, ]
      mu_i <- mu[i, obs]
      f_i  <- f[i, obs]
      d_mu_i <- diag(mu_i, length(mu_i)) # Length argument guarantees diag creation
      if (method=="ML") {
        V_i <- sig2 * d_mu_i
      } else if (method=="GEE") {
        idx <- which(obs)
        R_i <- Rg[idx,idx]
        V_i <- sig2 * sqrt(d_mu_i) %*% R_i %*% sqrt(d_mu_i)
      } else stop("Can't happen")
      if (any(abs(diag(V_i))<1e-12)) browser()
      V_inv[[i]] <<- solve(V_i) # Store $V^{-1}# for later use
      Omega[[i]] <<- d_mu_i %*% V_inv[[i]] %*% d_mu_i # idem for $\Omega_i$
    }
  }

  # The (optional) autocorrelation structure for any site $i$ is stored in $n_i\times n_i$ matrix $R_i$.
  # In case there are no missing values, $n_i=J$, and the `full' or `generic' autocorrelation matrix $R$ is expressed
  # as
  # \begin{equation}
  # R = \begin{pmatrix}
  #   1          & \rho       & \rho^2     & \cdots & \rho^{J-1} \\
  #   \rho       & 1          & \rho       & \cdots & \rho^{J-2} \\
  #   \vdots     & \vdots     & \vdots     & \ddots & \vdots     \\
  #   \rho^{J-1} & \rho^{J-2} & \rho^{J-3} & \cdots & 1
  #   \end{pmatrix}
  # \end{equation}
  # where $\rho$ is the lag-1 autocorrelation.
  Rg <- diag(1, ntime) # default (no autocorrelation) value
  update_R <- function() {
    Rg <<- rho ^ abs(row(diag(ntime)) - col(diag(ntime)))
  }

  # Lag-1 autocorrelation parameter $\rho$ is estimated as
  # \begin{equation}
  #   \hat{\rho} = \frac{1}{n_{i,j,j+1}\hat{\sigma}^2} \left(\Sum_i^I\Sum_j^{J-1} r_{i,j}r_{i,j+1}) \right)
  # \end{equation}
  # where the summation is over observed pairs $i,j$--$i,j+1$, and $n_{i,j,j+1}$ is the number of pairs involved.
  # Again, both $\rho$ and $R$ are computes $\rho$ in a stepwise per-site fashion.
  # Also, site-specific autocorrelation matrices $R_i$ are formed by removing the rows and columns from $R$
  # corresponding with missing observations.
  rho <- 0.0 # default value (ML)
  update_rho <- function() {
    # First estimate $\rho$
    rho   <-  0.0
    count <-  0
    for (i in 1:nsite) {
      for (j in 1:(ntime-1)) {
        if (observed[i,j] && observed[i,j+1]) { # short-circuit AND intended
          rho <- rho + r[i,j] * r[i,j+1]
          count <- count+1
        }
      }
    }
    rho <<- rho / (count * sig2) # compute and store in outer environment
    # if (rho < 0) rho <<- 0.0 # Don't allow negative autocorrelation
  }

  # ------------------------------------------ Overdispersion. \label{sig2} ----

  # Dispersion parameter $\sigma^2$ is estimated as
  # \begin{equation}
  #   \hat{\sigma}^2 = \frac{1}{n_f - n_\alpha - n_\beta} \sum_{i,j} r_{ij}^2
  # \end{equation}
  # where the $n$ terms are the number of observations, $\alpha$'s and $\beta$'s, respectively.
  # Summation is over the observed $i,j$ only.
  # and $r_{ij}$ are Pearson residuals (Section~\ref{r})
  sig2 <- 1.0 # default value (Maximum Likelihood case)
  update_sig2 <- function() {
    df <- sum(nobs) - length(alpha) - length(beta) # degrees of freedom
    sig2 <<- if (df>0) sum(r^2, na.rm=TRUE) / df else 1.0
    if (!is.finite(sig2)) stop("Overdispersion problem")
  }

  # -------------------------------------------- Pearson residuals\label{r} ----

  # Deviations between measured and estimated counts are quantified by the
  # Pearson residuals $r_{ij}$, given by
  # \begin{equation}
  #   r_{ij} = (f_{ij} - \mu_{ij}) / \sqrt{\mu_{ij}}
  # \end{equation}
  r <- matrix(0, nsite, ntime)
  update_r <- function() {
    r[observed] <<- (f[observed]-mu[observed]) / sqrt(mu[observed])
  }

  # -------------------------------------------- Derivatives and GEE scores ----

  # \label{Hessian}\label{Scoring}
  # Derivative matrix $i_b$ is defined as
  # \begin{equation}
  #   -i_b = \sum_i B_i' \left(\Omega_i - \frac{1}{d_i}\Omega_i z_i z_i' \Omega_i\right) B_i \label{i_b}
  # \end{equation}
  # where
  # \begin{equation}
  #   \Omega_i = \diag{\mu_i} V_i^{-1} \diag{\mu_i} \label{Omega_i}
  # \end{equation}
  # with $V_i$ the covariance matrix for site $i$, and
  # \begin{equation}
  #   d_i = z_i' \Omega_i z_i \label{d_i}
  # \end{equation}
  i_b <- 0
  U_b <- 0



  update_U_i <- function() {
    i_b <<- 0 # Also store in outer environment for later retrieval
    U_b <<- 0
    for (i in 1:nsite) {
      obs <- observed[i, ]
      B = make.B(i)
      mu_i <- mu[i, obs]
      f_i  <- f[i, obs]
      d_mu_i <- diag(mu_i, length(mu_i)) # Length argument guarantees diag creation
      ones <- matrix(1, nobs[i], 1)
      d_i <- as.numeric(t(ones) %*% Omega[[i]] %*% ones) # Could use sum(Omega) as well...
      B_i <- B[obs, ,drop=FALSE] # recyle index for e.g. covariates in $B$
      i_b <<- i_b - t(B_i) %*% (Omega[[i]] - (Omega[[i]] %*% ones %*% t(ones) %*% Omega[[i]]) / d_i) %*% B_i
      U_b <<- U_b + t(B_i) %*% d_mu_i %*% V_inv[[i]] %*% (f_i - mu_i)
    }
    if (any(abs(colSums(i_b))< 1e-12)) stop("Data does not contain enough information to estimate model.", call.=FALSE)
    # invertable <- class(try(solve(i_b), silent=T))=="matrix"
    # if (!invertable) {
    #   browser()
    # }
  }

  # ------------------------------------------------------- Count estimates ----

  # Let's not forget to provide a function to update the modelled counts $\mu_{ij}$:
  # $$ \mu^t = \exp(A\alpha^t + B\beta^{t-1} - \log w) $$
  # where it is noted that we do not use matrix $A$. Instead, the site-specific
  # parameters $\alpha_i$ are used directly:
  # $$ \mu_i^t = \exp(\alpha_i^t + B\beta^{t-1} - \log w) $$
  update_mu <- function(fill) {
    for (i in 1:nsite) {
      B = make.B(i)
      if (use.weights) {
        mu[i, ] <<- (exp(alpha[i] + B %*% beta) / wt[i, ])
      } else {
        mu[i, ] <<- exp(alpha[i] + B %*% beta)
      }
    }
    # clear estimates for non-observed cases, if required.
    if (!fill) mu[!observed] <<- 0.0
  }

  # ------------------------------------------------------------ Likelihood ----

  likelihood <- function() {
    # lik <- 2*sum(f*log(f/mu), na.rm=TRUE)
    lik <- 2*sum(f[positive]*log(f[positive]/mu[positive]))
    # ok = f>1e-6 & mu>1e-6
    # lik <- 2*sum(f[ok]*log(f[ok]/mu[ok]))
    lik
  }

  # ----------------------------------------------------------- Convergence ----

  # The parameter estimation algorithm iterates until convergence is reached.
  # `convergence' here is defined in a multivariate way: we demand convergence in
  # model paramaters $\alpha$ and $\beta$, model estimates $\mu$ and likelihood measure $L$.
  new_par <- new_cnt <- new_lik <- new_rho <- new_sig <- NULL
  old_par <- old_cnt <- old_lik <- old_rho <- old_sig <- NULL

  check_convergence <- function(iter) {

    # Collect new data for convergence test
    # (Store in outer environment to make them persistent)
    new_par <<- c(as.vector(alpha), as.vector(beta))
    new_cnt <<- as.vector(mu)
    new_lik <<- likelihood()
    new_rho <<- rho
    new_sig <<- sig2

    if (iter>1) {
      max_par_change <- max(abs(new_par - old_par))
      max_cnt_change <- max(abs(new_cnt - old_cnt))
      max_lik_change <- max(abs(new_lik - old_lik))
      rho_change <- abs(new_rho - old_rho)
      sig_change <- abs(new_sig - old_sig)
      beta_change <- max_dbeta
      conv_par <- max_par_change < conv_crit
      conv_cnt <- max_cnt_change < conv_crit
      conv_lik <- max_lik_change < conv_crit
      conv_rho <- rho_change < conv_crit
      conv_sig <- sig_change < conv_crit
      conv_beta <- beta_change < conv_crit
      # convergence <- conv_par && conv_cnt && conv_lik
      convergence <- conv_rho && conv_sig && conv_beta
      # rprintf(" Max change: %10e %10e %10e ", max_par_change, max_cnt_change, max_lik_change)
      rprintf(" Max change: %10e %10e %10e ", rho_change, sig_change, beta_change)
    } else {
      convergence = FALSE
    }

    # Today's new stats are tomorrow's old stats
    old_par <<- new_par
    old_cnt <<- new_cnt
    old_lik <<- new_lik
    old_rho <<- new_rho
    old_sig <<- new_sig

    convergence
  }

  # --------------------------------------------- Main estimation procedure ----

  # Now we have all the building blocks ready to start the iteration procedure.
  # We start `smooth', with a couple of Maximum Likelihood iterations
  # (i.e., not considering $\sigma^2\neq1$ or $\rho>0$), after which we move to on
  # GEE iterations if requested.
  method    <- "ML" # start with Maximum Likelihood
  final_method <- ifelse(serialcor || overdisp, "GEE", "ML") # optionally move on to GEE
  if (use.covin) final_method <- "ML" # fall back during upscaling

  update_mu(fill=FALSE)
  # browser()

  for (iter in 1:max_iter) {
    if (compatible && iter==4) method <- final_method
    rprintf("Iteration %d (%s)", iter, method)
    update_alpha(method)
    if (!is.null(err.out)) return(err.out)
    update_mu(fill=FALSE)
    if (method=="GEE") {
      update_r()
      if (serialcor || overdisp)  update_sig2() # hack
      if (serialcor) update_rho()
      update_R()
      if (!overdisp) sig2 <- 1.0 # hack
    }
    update_V(method)
    update_U_i() # update Score $U_b$ and Fisher Information $i_b$
    if (model>1) {
      subiters <- update_beta(method)
      if (!is.null(err.out)) return(err.out)
      rprintf(", %d subiters", subiters)
    }
    rprintf(", lik=%.3f", likelihood())
    if (overdisp)  rprintf(", sig^2=%.3f", sig2)
    if (serialcor) rprintf(", rho=%.3f;", rho)

    convergence <- check_convergence(iter)

    if (graph_debug) {
      nn = colSums(observed)
      fp = colSums(f, na.rm=TRUE) / nn
      mp = colSums(mu, na.rm=TRUE)/ nn
      plot(1:ntime, fp, type='b', ylim=range(c(fp,mp), na.rm=TRUE))
      points(1:ntime, mp, col="red")
      Sys.sleep(0.1)
    }
    if (convergence && method==final_method) {
      rprintf("\n")
      break
    } else if (convergence) {
      rprintf("\nChanging ML --> GEE\n")
      method = "GEE"
    } else {
      rprintf("\n")
    }
  }

  # If we reach the preset maximum number of iterations, we clearly have not reached
  # convergence.
  converged <- iter < max_iter
  convergence_msg <- sprintf("%s reached after %d iterations",
                             ifelse(convergence,"Convergence","No convergence"),
                             iter)
  rprintf("%s\n", convergence_msg)

  # Run the final model
  update_mu(fill=TRUE)

  # Covariance matrix
  if (model>1) V <- -solve(i_b)

  if (use.covin) {
    UUT <- matrix(0,nbeta,nbeta)
    for (i in 1:nsite) {
      B <- make.B(i)
      # mu_i <- mu[site==i & observed]
      # n_i <- length(mu_i)
      # d_mu_i <- diag(mu_i, n_i) # Length argument guarantees diag creation
      # OM <- Omega[[i]]
      # d_i <- sum(OM) # equivalent with z' Omega z, as in the TRIM manual
      Bi <- B[observed[site==i], ,drop=FALSE]
      var_fi <- covin[[i]]
      vi <- rowSums(var_fi)
      vipp <- sum(var_fi)
      mui <- mu[site==i & observed]
      muip <- sum(mui)
      term1 <- var_fi
      term2 <- vi %*% t(mui) / muip
      term3 <- mui %*% t(vi) / muip
      term4 <- vipp * mui %*% t(mui) / muip^2
      UUT <- UUT + t(Bi) %*% (term1-term2-term3+term4) %*% Bi
    }
    V <- V %*% UUT %*% V
  }

  # Variance of the $\beta$'s is given by the covariance matrix
  if (model>1) var_beta = V

  # ============================================================ Imputation ====

  # The imputation process itself is trivial: just replace all missing observations
  # $f_{i,j}$ by the model-based estimates $\mu_{i,j}$.
  imputed <- ifelse(observed, f, mu)


  # ============================================= Output and postprocessing ====

  # Measured, modelled and imputed count data are stored in a TRIM output object,
  # together with parameter values and other usefull information.

  # Convert time point back to their original (numerical) values
  time.id <- as.numeric(levels(time.id))
  site.id <- factor(levels(site.id))

  z <- list(title=title, f=f, nsite=nsite, ntime=ntime, time.id=time.id,
            site.id = site.id,
            nbeta0=nbeta0, covars=covars, ncovar=ncovar, cvmat=cvmat,
            model=model, changepoints=changepoints, converged=converged,
            mu=mu, imputed=imputed, alpha=alpha)
  if (model>1) {
    z$beta <- beta
    z$var_beta <- var_beta
  }

  z$method <- ifelse(use.covin, "Pseudo ML", final_method)
  z$convergence <- convergence_msg

  # mark the original request for changepoints to allow wald(z) to perform a dslope
  # test in case of a single changepoint. (otherwise we cannot distinguish it from an
  # automatic changepoint 1)
  if (model==2) z$ucp <- use.changepoints

  if (use.covars) {
    z$ncovar <- ncovar # todo: eliminate?
    z$nclass <- nclass
  }
  if (use.weights) {
    z$wt = wt
  } else {
    z$wt = NULL
  }

  class(z) <- "trim"

  # Several kinds of statistics can now be computed, and added to this output object.

  # ------------------------------------- Overdispersion and Autocorrelation ---

  # remove if not required: makes summary printing and extracting more elegant.
  z$sig2 <- if (overdisp) sig2 else NULL
  z$rho  <- if (serialcor) rho else  NULL

  #-------------------------------------------- Coefficients and uncertainty ---

  if (model==1) {
    z$coefficients <- NULL
  }

  if (model==2) {
    se_beta  <- sqrt(diag(var_beta))

    ncp <- length(changepoints)
    from_cp <- changepoints
    upto_cp <- if (ncp==1) ntime else c(changepoints[2:ncp], ntime)
    coefs = data.frame(
      from   = time.id[from_cp],
      upto   = time.id[upto_cp],
      add    = beta,
      se_add = se_beta,
      mul    = exp(beta),
      se_mul = exp(beta) * se_beta
    )

    if (use.covars) {
      # Add some prefix columns with covariate and factor ID
      # Note that we have to specify all covariate levels here, to prevent them
      # from being to converted to NA later
      prefix <- data.frame(covar = factor("baseline",levels=c("baseline", names(covars))),
                           cat   = 0)
      coefs <- cbind(prefix, coefs)
      idx = 1:nbeta0
      for (i in 1:ncovar) {
        for (j in 2:nclass[i]) {
          idx <- idx + nbeta0
          coefs$covar[idx]  <- names(covars)[i]
          coefs$cat[idx]    <- j
        }
      }
    }

    z$coefficients <- coefs
  } # if model==2

  if (model==3) {
    # Model coefficients are output in two types; as additive parameters:
    # $$ \log\mu_{ij} = \alpha_i + \gamma_j $$
    # and as multiplicative parameters:
    # $$ \mu_{ij} = a_i g_j $$
    # where $a_i=e^{\alpha_i}$ and $g_j = e^{\gamma_j}$.

    # For the first time point, $\gamma_1=0$ by definition.
    # So we have to add values of 0 for the baseline case and each covariate category $>1$, if any.
    #gamma <- matrix(beta, nrow=nbeta0) # Each covariate category in a column
    #gamma <- rbind(0, gamma) # Add row of 0's for first time point
    #gamma <- matrix(gamma, ncol=1) # Cast back into a column vector
    gamma <- beta
    g     <- exp(gamma)

    # Parameter uncertainty is expressed as standard errors.
    # For the additive parameters $\gamma$, the variance is estimated as
    # $$ \var{\gamma} = (-i_b)^{-1} $$
    var_gamma <-  -solve(i_b) # BUG: should be var_beta (inc covin effect)
    # Finally, we compute the standard error as $\se{\gamma} = \sqrt{\diag{\var{\gamma}}}$
    se_gamma  <-  sqrt(diag(var_gamma))

    # The standard error of the multiplicative parameters $g_j$ is opproximated by
    # using the delta method, which is based on a Taylor expansion:
    # \begin{equation}
    #   \var{f(\theta)} = \left(f'(\theta)\right)^2 \var{\theta}
    # \end{equation}
    # which for $f(\theta)=e^\theta$ translates to
    # $$ \var{g} = \var{e^{\gamma}} = e^{2\gamma} \var{\gamma} $$
    # leading to
    # $$ \se{g} = e^{\gamma} \se{\gamma} = g \se{\gamma} $$
    se_g <- g * se_gamma

    # Baseline coefficients.
    # Note that, because $\gamma_1\equiv0$, it was not estimated,
    # and as a results $j=1$ was not incuded in $i_b$, nor in $\var{gamma}$ as computed above.
    # We correct this by adding the `missing' 0 (or 1 for multiplicative parameters) during output
    idx = 1:nbeta0
    coefs <- data.frame(
      time   = 1:ntime,
      add    = c(0, gamma[idx]),
      se_add = c(0, se_gamma[idx]),
      mul    = c(1, g[idx]),
      se_mul = c(0, g[idx] * se_gamma[idx])
    )

    # Covariate categories ($>1$)
    if (use.covars) {
      prefix = data.frame(covar="baseline", cat=0)
      coefs <- cbind(prefix, coefs)
      for (i in 1:ncovar) {
        for (j in 2:nclass[i]) {
          idx <- idx + nbeta0
          df <- data.frame(
            covar  = names(covars)[i],
            cat    = j,
            time   = 1:ntime,
            add    = c(0, gamma[idx]),
            se_add = c(0, se_gamma[idx]),
            mul    = c(1, g[idx]),
            se_mul = c(0, g[idx] * se_gamma[idx])
          )
          coefs <- rbind(coefs, df)
        }
      }
    }

    z$coefficients <- coefs
  } # if model==3


  # ----------------------------------------------------------- Time totals ----

  # Recompute Score matrix $i_b$ with final $\mu$'s
  saved.ib = i_b
  ib <- 0
  for (i in 1:nsite) {
    B <- make.B(i)
    mu_i <- mu[site==i & observed]
    n_i <- length(mu_i)
    d_mu_i <- diag(mu_i, n_i) # Length argument guarantees diag creation
    OM <- Omega[[i]]
    d_i <- sum(OM) # equivalent with z' Omega z, as in the TRIM manual
    B_i <- B[observed[site==i], ,drop=FALSE]
    om <- colSums(OM)
    OMzzOM <- om %*% t(om) # equivalent with OM z z' OM, as in the TRIM manual
    term <- t(B_i) %*% (OM - (OMzzOM) / d_i) %*% B_i
    ib <- ib - term
  }

  # Matrices E and F take missings into account
  E <- -ib # replace by covin hessian

  rep.rows <- function(row, n) matrix(rep(row, each=n), nrow=n)
  rep.cols <- function(col, n) matrix(rep(col, times=n), ncol=n)

  var_model_tt <- function(observed_only=FALSE, mask=NULL) {

    if (!use.covin) { # use computed covariance

      F <- matrix(0, nsite, nbeta)
      d <- numeric(nsite)
      for (i in 1:nsite) {
        B <- make.B(i)
        d[i] <- sum(Omega[[i]])
        w_i <- colSums(Omega[[i]])
        B_i <- B[observed[site==i], ,drop=FALSE]
        F_i <- (t(w_i) %*% B_i) / d[i]
        F[i, ] <- F_i
      }

      # Matrices G and H are for all (weighted) mu's
      wmu <- if (use.weights) wt * mu else mu
      # Zero-out unobserved $i,j$ to facilitate the computation of the variance of the imputed counts
      if (observed_only) wmu[!observed] = 0.0
      # Zero-out unselected sites to facilitate the computation of the variance of covariate categories
      if (!is.null(mask)) wmu[!mask]     = 0.0

      GddG <- matrix(0, ntime,ntime)
      for (i in 1:nsite) {
        for (j in 1:ntime) for (k in 1:ntime) {
          GddG[j,k] <- GddG[j,k] + wmu[i,j]*wmu[i,k]/d[i]
        }
      }

      GF <- matrix(0, ntime, nbeta)
      for (i in 1:nsite) {
        for (j in 1:ntime) for (k in 1:nbeta)  {
          GF[j,k] <- GF[j,k] + wmu[i,j] * F[i,k]
        }
      }

      H <- matrix(0, ntime, nbeta)
      for (i in 1:nsite) {
        B <- make.B(i)
        for (k in 1:nbeta) for (j in 1:ntime) {
          H[j,k]  <- H[j,k] + B[j,k] * wmu[i,j]
        }
      }

      GFminH <- GF - H

      # All building blocks are ready. Use them to compute the variance
      if (model==1) V <- GddG
      else          V <- GddG + GFminH %*% solve(E) %*% t(GFminH)

    } else { # Use input covariance
      if (!is.null(mask)) stop("Alas, covariates+upscaling not implemented yet.");
      # First loop op sites to compute $GF - H$.
      GFminH = matrix(0, ntime, nbeta)
      for (i in 1:nsite) { # First loop
        u <-  mu[i, ]
        wu <- wt[i, ] * mu[i, ]
        obs <- observed[i, ] # observed j's
        if (observed_only) wu[!obs] <- 0.0
        w <- mu[i, obs]
        d <- sum(w)
        Bi = make.B(i)
        Bo = Bi[obs, ]
        w_mat = rep.cols(w, nbeta)
        F = colSums(Bo * w_mat / d)
        wu_mat <- rep.cols(wu, nbeta)
        F_mat <- rep.rows(F, ntime)
        GFminH <- GFminH + wu_mat * (F_mat - Bi)
      }
      M2 <- GFminH %*% solve(E)

      # second loop
      V <- matrix(0, ntime, ntime)
      for (i in 1:nsite) {
        u <-  mu[i, ]
        wu <- wt[i, ] * mu[i, ]
        obs <- observed[i, ] # observed j's
        if (observed_only) wu[!obs] <- 0.0
        w <- mu[i, obs]
        d <- sum(w)
        M1 <- rep.cols(wu/d, nobs[i])
        Bi = make.B(i)
        Bo = Bi[obs, ]
        w_mat = rep.cols(w, nbeta)
        F = colSums(Bo * w_mat / d)
        F_mat <- rep.rows(F, nobs[i])
        M3 <- t(F_mat - Bo)
        M4 = M1 + M2 %*% M3
        Vi <- M4 %*% covin[[i]] %*% t(M4)
        V <- V + Vi
      }
    }
    V
  }

  var_tt_mod <- var_model_tt(observed_only = FALSE)

  # To compute the variance of the time totals of the imputed data, we first
  # substract the contribution due to te observations, as computed by above scheme,
  # and replace it by the contribution due to the observations, as resulting from the
  # covariance matrix.

  var_observed_tt <- function(mask=NULL) {
    # Variance due to observations
    V = matrix(0, ntime, ntime)
    if (!use.covin) {
      wwmu <- if (use.weights) wt * wt * mu else mu
      wwmu[!observed] <- 0 # # erase estimated $\mu$'s
      if (!is.null(mask)) wmu[!mask] <- 0.0 # Erase 'other' covariate categories also.
      for (i in 1:nsite) {
        if (serialcor) {
          srdu = sqrt(diag(wwmu[i, ]))
          Vi = sig2 * srdu %*% Rg %*% srdu
        } else {
          Vi = sig2 * diag(wwmu[i, ])
        }
        V = V + Vi
      }
    } else {
      for (i in 1:nsite) {
        Vi = covin[[i]]
        obs <- observed[i, ]
        wt1 = rep.rows(wt[i, obs], nobs[i])
        wt2 = rep.cols(wt[i, obs], nobs[i])
        V[obs,obs] <- V[obs,obs] + wt1 * wt2 * Vi
      }
    }
    V
  }

  # Combine
  var_tt_obs_old <- var_model_tt(observed_only=TRUE)
  var_tt_obs_new <- var_observed_tt()
  var_tt_imp = var_tt_mod - var_tt_obs_old + var_tt_obs_new

  # Time totals of the model, and it's standard error
  wmu <- if (use.weights) wt * mu else mu
  tt_mod    <- colSums(wmu)
  se_tt_mod <- round(sqrt(diag(var_tt_mod)))

  wimp <- if (use.weights) wt * imputed  else imputed #kan eleganter

  tt_imp     <- colSums(wimp)
  se_tt_imp <- round(sqrt(diag(var_tt_imp)))

  # Store in TRIM output
  z$tt_mod <- tt_mod
  z$tt_imp <- tt_imp
  z$var_tt_mod <- var_tt_mod
  z$var_tt_imp <- var_tt_imp

  z$time.totals <- data.frame(
    time    = time.id,
    fitted  = round(tt_mod),
    se_fit  = se_tt_mod,
    imputed = round(tt_imp),
    se_imp  = se_tt_imp
  )

  # Indices for covariates. First baseline
  if (use.covars) {
    z$covar_tt <- list()
    for (i in 1:ncovar) {
      cvname <- names(covars)[i]
      cvtt = vector("list", nclass[i])
      for (j in 1:nclass[i]) {
        classname <- ifelse(is.factor(covars[[i]]), levels(covars[[i]])[j], j)
        mask <- cvmat[[i]]==j
        # Time totals (fitted)
        mux <- wmu
        mux[!mask] <- 0.0
        tt_mod <- colSums(mux)
        # Corresponding variance
        var_tt_mod <- var_model_tt(mask=mask)
        # Time-total (imputed)
        impx <- wimp
        impx[!mask] <- 0.0
        tt_imp <- colSums(impx)
        # Variance
        var_tt_obs_old <- var_model_tt(observed_only=TRUE, mask=mask)
        var_tt_obs_new <- var_observed_tt(mask)
        var_tt_imp = var_tt_mod - var_tt_obs_old + var_tt_obs_new
        # Store
        df <- list(covariate=cvname, class=j, mod=tt_mod, var_mod=var_tt_mod, imp=tt_imp, var_imp=var_tt_imp)
        cvtt[[j]] <- df
      }
      z$covar_tt[[cvname]] <- cvtt
    }
  } else z$covar_tt <- NULL

  # ----------------------------------------- Reparameterisation of Model 3 ----

  # Here we consider the reparameterization of the time-effects model in terms of
  # a model with a linear trend and deviations from this linear trend for each time point.
  # The time-effects model is given by
  # \begin{equation}
  #   \log\mu_{ij}=\alpha_i+\gamma_j,
  # \end{equation}
  # with $\gamma_j$ the effect for time $j$ on the log-expected counts and $\gamma_1=0$. This reparameterization can be expressed as
  # \begin{equation}
  #   \log\mu_{ij}=\alpha^*_i+\beta^*d_j+\gamma^*_j,
  # \end{equation}
  # with $d_j=j-\bar{j}$ and $\bar{j}$ the mean of the integers $j$ representing
  # the time points.
  # The parameter $\alpha^*_i$ is the intercept and the parameter $\beta^*$ is
  # the slope of the least squares regression line through the $J$ log-expected
  # time counts in site $i$ and  $\gamma^*_j$ can be seen as the residuals of this
  # linear fit.
  # From regression theory we have that the `residuals'"'  $\gamma^*_j$ sum to zero
  # and are orthogonal to the explanatory variable, i.e.
  # \begin{equation}
  #   \sum_j\gamma^*_j = 0 \quad \text{and} \quad \sum_jd_j\gamma^*_j = 0. \label{constraints}
  # \end{equation}
  # Using these constraints we obtain the equations:
  # \begin{gather}
  #   \log\mu_{ij}           = \alpha^*_i+\beta^*d_j+\gamma^*_j=\alpha_i+\gamma_j  \label{repar1}\\
  #   \sum_j \log\mu_{ij}    = J\alpha^*_j = J\alpha_i+\sum_j\gamma_j \label{repar2}\\
  #   \sum_j d_j\log\mu_{ij} = \beta^*\sum_jd^2_j = \sum_jd_j\gamma_j \label{repar3},
  # \end{gather}
  # where \eqref{repar1} is the re-parameterization equation itself and \eqref{repar2}
  # and \eqref{repar3} are obtained by using the constraints~\eqref{constraints}

  # From \eqref{repar2} we have that $\alpha^*_i=\alpha_i+\frac{1}{J}\sum_j\gamma_j$.
  # Now, by using the equations \eqref{repar1} thru \eqref{repar3} and defining
  # $D=\sum_jd^2_j$, we can express the parameters $\beta^*$ and $\gamma^*$ as
  # functions of the parameters $\gamma$ as follows:
  # \begin{align}
  #   \label{betaster}
  #   \beta^* &=\frac{1}{D}\sum_jd_j\gamma_j,\\ \nonumber
  #   \label{gammaster}
  #   \gamma^*_j &= \alpha_i+\gamma_j-\alpha^*_i-\beta^*d_j  \quad (\text{using (5)})\\ \nonumber
  #   &=\alpha_i-\left( \alpha_i+\frac{1}{J}\sum_j\gamma_j\right) +\gamma_j-d_j\frac{1}{D}\sum_jd_j\gamma_j \\
  #   &=\gamma_j-\frac{1}{J}\sum_j\gamma_j-d_j\frac{1}{D}\sum_jd_j\gamma_j.
  # \end{align}
  # Since $\beta^*$ and $\gamma^*_j$ are linear functions of the parameters $\gamma_j$
  # they can be expressed in matrix notation by
  # \begin{equation}
  #   \left ( \begin{array} {c}
  #          \beta^* \\
  #          \boldsymbol{\gamma}^*
  #   \end{array} \right ) = \mathbf{T}\boldsymbol{\gamma},
  # \end{equation}
  # with $\boldsymbol{\gamma}^*=(\gamma^*_1,\ldots,\gamma^*_J)^T$,
  # $\boldsymbol{\gamma}=(\gamma_1,\ldots,\gamma_J)^T$ and $\mathbf{T}$
  # the $(J+1) \times J$ transformation matrix that transforms
  # $\boldsymbol{\gamma}$ to  $\left (\beta^*,(\boldsymbol{\gamma}^*)^T\right)^T$.
  # From \eqref{betaster} and \eqref{gammaster} it follows that the elements of
  # $\mathbf{T}$ are given by:
  # \begin{align}
  #   \label{matrixT} \nonumber
  #   &\mathbf{T}_{(1,j)}=\frac{d_j}{D} &\quad (i=1,j=1,\ldots,J)\\ \nonumber
  #   &\mathbf{T}_{(i,j)}=1-\frac{1}{J}-\frac{1}{D}d_{i-1}d_j &\quad(i=2,\ldots,J+1,j=1,\ldots,J,i-1=j)\\ \nonumber
  #   &\mathbf{T}_{(i,j)}=-\frac{1}{J}-\frac{1}{D}d_{i-1}d_j &\quad(i=2,\ldots,J+1,j=1,\ldots,J,i-1 \neq j)
  # \end{align}

  if (model==3 && !use.covars) {

    TT <- matrix(0, ntime+1, ntime)
    J <- ntime
    j <- 1:J; d <- j - mean(j) # i.e, $ d_j = j-\frac{1}{J}\sum_j j$
    D <- sum(d^2)             # i.e., $ D = \sum_j d_j^2$
    TT[1, ] <- d / D
    for (i in 2:(J+1)) for (j in 1:J) {
      if (i-1 == j) {
        TT[i,j] <- 1 - (1/J) - d[i-1]*d[j]/D
      } else {
        TT[i,j] <-   - (1/J) - d[i-1]*d[j]/D
      }
    }

    gstar <- TT %*% c(0, gamma) # Add the implicit $\gamma_1=0$
    bstar <- gstar[1]
    gstar <- gstar[2:(J+1)]

    # The covariance matrix of the transformed parameter vector can now be obtained
    # from the covariance matrix $\mathbf{T}\boldsymbol{\gamma}$ of $\boldsymbol{\gamma}$ as
    # \begin{equation}
    #   V\left( \begin{array} {c} \beta^* \\\boldsymbol{\gamma}^* \end{array} \right)
    #   = \mathbf{T}V(\boldsymbol{\gamma})\mathbf{T}^T
    # \end{equation}

    var_gstar <- TT %*% rbind(0,cbind(0,var_gamma)) %*% t(TT) # Again, $\gamma_1=0$
    se_bstar  <- sqrt(diag(var_gstar))[1]
    se_gstar  <- sqrt(diag(var_gstar))[2:(ntime+1)]

    z$gstar <- gstar
    z$var_gstar <- var_gstar

    z$linear.trend <- data.frame(
      Additive       = bstar,
      std.err        = se_bstar,
      Multiplicative = exp(bstar),
      std.err.       = exp(bstar) * se_bstar,
      row.names      = "Slope",
      check.names    = FALSE)

    # Deviations from the linear trend
    z$deviations <- data.frame(
      Time       = 1:ntime,
      Additive   = gstar,
      std.err.   = se_gstar,
      Multiplicative = exp(gstar),
      std.err.   = exp(gstar) * se_gstar,
      check.names = FALSE
    )

  }


  # ======================================================== Return results ====

  # The TRIM result is returned to the user\ldots
  rprintf("(Exiting workhorse function)\n")
  z
}
# \ldots which ends the main TRIM function.
