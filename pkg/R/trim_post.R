# ########################################### TRIM postprocessing functions ####

# ================================================================= Summary ====

# ----------------------------------------------------------------- extract ----

#' Extract summary information for a TRIM job
#'
#' @param object TRIM output structure (i.e., output of a call to \code{trim})
#' @param ... Currently unused
#'
#' @return a list of type \code{trim.summary} containing elements for estimation method
#'   (\code{est.method}), overdispersion (\code{sig2}) and autocorrelation (\code{rho})
#' @export
#'
#' @family analyses
#' @seealso \code{\link{trim}}
#' @examples
#'
#' data(skylark)
#' z <- trim(skylark, count ~ time + site,model=2,overdisp=TRUE)
#' summary(z)
#' # extract autocorrelation strength
#' rho <- summary(z)$rho
summary.trim <- function(object,...) {
  x <- object
  if (is.finite(x$sig2) || is.finite(x$rho)) {
    out = list(est.method="Generalised Estimating Equations")
  } else {
    out = list(est.method="Maximum Likelihood")
  }
  out$sig2 <- x$sig2
  out$rho  <- x$rho
  class(out) <- "trim.summary"
  out
}

# ------------------------------------------------------------------- print ----

#' Print method for \code{trim.summary}
#'
#' @export
#' @param x An object of class \code{trim.summary}
#' @keywords internal
print.trim.summary <- function(x,...) {
  printf("\nEstimation method = %s\n", x$est.method)
  if (is.finite(x$sig2)) printf("  Estimated Overdispersion     = %f\n", x$sig2)
  if (is.finite(x$rho))  printf("  Estimated Serial Correlation = %f\n", x$rho)
}


# ============================================================ Coefficients ====

# ----------------------------------------------------------------- Extract ----

#' Extract TRIM model coefficients.
#'
#' @param object TRIM output structure (i.e., output of a call to \code{trim})
#' @param which What coefficients to return.
#' @param ... currently unused
#'
#' @return a list containing the model type (element \code{model}, 1,2 or 3), and
#' element \code{coef}, which is a data frame containing the actual coefficients.
#' @export
#'
#' @family analyses
#' @examples
#' data(skylark)
#' z <- trim(skylark, count ~ time + site,model=2,overdisp=TRUE)
#' summary(z)
#' # extract autocorrelation strength
#' rho <- summary(z)$rho
coef.trim <- function(object, which=c("additive","multiplicative","both"),...) {

  # Craft a custom output
  which <- match.arg(which)

  # Last 4 columns contain the additive and multiplicative parameters.
  # Select the appropriate subset from these, and all columns before these 4.
  n = ncol(object$coefficients)
  stopifnot(n>=4)
  if (which=="additive") {
    cols <- 1:(n-2)
  } else if (which=="multiplicative") {
    cols <- c(1:(n-4), (n-1):n)
  } else if (which=="both") {
    cols = 1:n
  }
  out <- object$coefficients[cols]
  out
}


# ============================================================= Time totals ====

# ----------------------------------------------------------------- Extract ----

#' Extract time-totals from TRIM output
#'
#' @param x TRIM output structure (i.e., output of a call to \code{trim})
#' @param which Selector to distinguish between time totals based on the imputed data (default),
#' the modelled data, or both.
#'
#' @return a structure of class \code{trim.totals}, which is a list with as single
#' element the data frame \code{totals}.
#' @export
#'
#' @family analyses
#' @examples
#' data(skylark)
#' z <- trim(count ~ time + site, data=skylark, model=2);
#' totals(z)
#' #print(totals(z,"imputed")) # idem
#' #totals(z, "both") # mimics classic TRIM
#' #SE <- totals(z)$totals$std.err
totals <- function(x, which=c("imputed","model","both")) {
  stopifnot(class(x)=="trim")

  # Select output columns from the pre-computed time totals
  which <- match.arg(which)
  totals <- switch(which
    , model   = x$time.totals[c(1,2,3)]
    , imputed = x$time.totals[c(1,4,5)]
    , both    = x$time.totals
  )

  # wrap the time.index field in a list and make it an S3 class
  # (so that it becomes recognizable as a TRIM time-indices)
  class(totals) <- c("trim.totals","data.frame")
  totals
}

#------------------------------------------------------------------ Export ----

export <- function(x, species, stratum) UseMethod("export")

export.trim.totals <- function(x, species, stratum) {
  stopifnot(class(x)=="trim.totals")

  # Create extra columns to be put before the actual time totals
  df1 = data.frame(species=species, stratum=stratum)
  df2 = x$totals
  df = cbind(df1, df2)
  print(df, row.names=FALSE)
}

#------------------------------------------------------------------- Print ----

#' Print method for trim time totals
#'
#' @param tt \code{trim.totals} object
#' @param ... currently unused
#'
#' @export
#' @keywords internal
print.trim.totals <- function(tt,...) {
  printf("Time totals\n")
  print.data.frame(tt, row.names=FALSE)
}


# ============================================== Variance-Covariance matrix ====

# ----------------------------------------------------------------- extract ----

#' Extract variance-covariance matrix from TRIM output
#'
#' @param x TRIM output structure (i.e., output of a call to \code{trim})
#' @param which Selector to distinguish between variance-covariance based on the
#' imputed data (default), or the modelled data.
#'
#' @return a JxJ matrix, where J is the number or time points.
#' @export
#'
#' @family analyses
#' @examples
#' data(skylark)
#' z <- trim(count ~ time + site, data=skylark, model=2);
#' totals(z)
#' vcv1 <- varcovar(z)       # Use imputed data
#' vcv2 <- varcovar(z,"mod") # Use modelled data
varcovar <- function(x, which=c("imputed","model")) {
  stopifnot(inherits(x,"trim"))

  which <- match.arg(which)
  vcv <- switch(which
    , model   = x$var_tt_mod
    , imputed = x$var_tt_imp
  )
  vcv
}



# =========================================== Reparameterisation of Model 3 ====

# ----------------------------------------------------------------- extract ----

#' Extract coefficients of the reparameterisation of TRIM model 3
#'
#' @param x TRIM output structure (i.e., output of a call to \code{trim})
#'
#' @return a list with elements
#'   \code{trend}, containing additive and multiplicative parameters of the linear trend,
#'   and \code{dev}, containing the deviations from that trend.
#' @export
#'
#' @family analyses
#' @examples
#' data(skylark)
#' z <- trim(count ~ time + site, data=skylark, model=3)
#' # print coefficients of the linear trend
#' print(linear(z))
#' # get linear trend magnitude
#' slope <- linear(z)$trend$Multiplicative
#' # ... and the deviations from that trend
#' devs  <- linear(z)$dev$Multiplicative
linear <- function(x) {
  stopifnot(inherits(x,"trim"))
  if (x$model !=3 ){
    message("Cannot extract linear coefficients from TRIM model %d",x$model)
    return(NULL)
  }

  structure(list(trend=x$linear.trend, dev=x$deviations), class="trim.linear")
}

# ------------------------------------------------------------------- print ----


#' Print an object of class trim.linear
#'
#' @param x An object of class \code{trim.linear}
#'
#' @export
#' @keywords internal
print.trim.linear <- function(x,...) {

  printf("Linear Trend + Deviations for Each Time\n")
  print(x$trend, row.names=TRUE)
  printf("\n")
  print(x$dev, row.names=FALSE)
}




