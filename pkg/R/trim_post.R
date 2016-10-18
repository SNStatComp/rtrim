# ########################################### TRIM postprocessing functions ####

# ================================================================= Summary ====

# ----------------------------------------------------------------- extract ----


#' Print summary information for a TRIM job
#'
#' Print a summary of a \code{\link{trim}} object.
#'
#' @param object an object of class \code{\link{trim}}.
#' @param ... Currently unused
#'
#' @return \code{NULL}, invisibly.
#' @export
#'
#' @family analyses
#' @seealso \code{\link{trim}}
#' @examples
#'
#' data(skylark)
#' z <- trim(skylark, count ~ time + site,model=2,overdisp=TRUE)
#' summary(z)
summary.trim <- function(object,...) {

  cl <- paste(capture.output(print(object$call)),collapse="\n")
  printf("Call:\n%s\n",cl)

  printf("\nCoefficients:\n")
  print(coef.trim(object,"both"))
  printf("\n")

  printf(" Overdispersion    : %8.4f\n",object$rho)
  printf(" Serial Correlation: %8.4f\n",object$sig2)
  printf("\n")

  print(gof(object))
  invisible(NULL)
}


#' Extract serial correlation from TRIM object
#'
#' @param x An object of class \code{\link{trim}}
#'
#' @return The serial correlation coefficient if computed, otherwise \code{NULL}.
#'
#' @export
#' @family analyses
serial_correlation <- function(x){
  stopifnot(inherits(x,"trim"))
  x$rho
}


#' Extract overdispersion from trim object
#'
#' @param x An object of class \code{\link{trim}}
#'
#' @return The overdispersion value if computed, otherwise \code{NULL}.
#'
#' @export
#' @family analyses
#'
overdispersion <- function(x){
  stopifnot(inherits(x,"trim"))
  x$sig2
}


# ============================================================ Coefficients ====

# ----------------------------------------------------------------- Extract ----

#' Extract TRIM model coefficients.
#'
#' @section Details:
#' 
#' Extract the site, growth or time effect parameters computed with
#' \code{\link{trim}}.
#' 
#' @section Additive versus multiplicative representation:
#' 
#' In the simplest cases (no covariates, no change points), the trim
#' Model 2 and Model 3 can be summarized as follows:
#' 
#' \itemize{
#' \item{Model 2: \eqn{\ln\mu_{ij}=\alpha_i + \beta\times(j-1)} }
#' \item{Model 3: \eqn{\ln\mu_{ij}=\alpha_i + \gamma_j}.}
#' }
#' 
#' Here, \eqn{\mu_{ij}} is the estimated number of counts at site \eqn{i}, time 
#' \eqn{j}. The parameters \eqn{\alpha_i}, \eqn{\beta} and \eqn{\gamma_j} are 
#' refererred to as coefficients in the additive representation. Bby
#' exponentiating both sides of the above equations, alternative representations
#' can be written down. Explicitly, one can show that
#'
#' \itemize{
#' \item{Model 2: \eqn{\mu_{ij}= a_ib^{(j-1)} = b\mu_{ij-1}}, where \eqn{a_i=e^{\alpha_i}} and \eqn{b=e^\beta}.}
#' \item{Model 3: \eqn{\mu_{ij}=a_ic_j}, where \eqn{a_i=e^{\alpha_i}}, \eqn{c_1=1} and \eqn{c_j=e^{\gamma_j}} for \eqn{j>1}.}
#' }
#'
#' The parameters \eqn{a_i}, \eqn{b} and \eqn{c_j} are referred to as
#' coefficients in the \emph{multiplicative representation}.
#'
#' @param object TRIM output structure (i.e., output of a call to \code{trim})
#' @param which What coefficients to return.
#' @param ... currently unused
#'
#' @return A \code{data.frame} containing coefficients and their standard errors.
#' Depending on the requested type of coefficients column names are \code{add}
#' and \code{se_add} for additive coefficients and/or \code{mul} and \code{se_mul}
#' for multiplicative coefficients. For model 2, the output has columns
#' \code{from} and \code{upto}, indicating the time slices for which the coefficients
#' are valid. For model 3, a column \code{time} is present, indicating to which
#' time point each (set of) coefficient(s) pertain.
#'
#'
#' @export
#'
#' @family analyses
#' @examples
#' data(skylark)
#' z <- trim(skylark, count ~ time + site,model=2,overdisp=TRUE)
#' coefficients(z)
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
#' @param which select what totals to compute (see \code{Details} section).
#'
#' @return a \code{data.frame} with subclass \code{trim.totals} 
#'  (for pretty-printing). The columns are \code{time}, \code{model}
#'  and \code{se_mod} (for standard error), and/or \code{imputed}
#'  and \code{se_imp}, depending on the selection.
#' 
#' @section Details:
#' 
#' The idea of \code{TRIM} is to impute those site-time combinations where
#' no counts are available. Time-totals (i.e. summed over sites) can be obtained 
#' for two cases:
#' 
#' \itemize{
#' \item{\code{"imputed"}: Time totals are computed after replacing missing values with values predicted by the model}.
#' \item{\code{"model"}: Time totals are computed after replacing both missing values and observed values with
#' values predicted by the model.}
#' }
#' 
#' 
#' @return A \code{data.frame} with time totals and their standard errors.
#' 
#' @export
#'
#' @family analyses
#' @examples
#' data(skylark)
#' z <- trim(count ~ time + site, data=skylark, model=2,cp=c(3,5))
#' totals(z)
#' 
#' totals(z, "both") # mimics classic TRIM
#' 
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
print.trim.totals <- function(x,...) {
  printf("Time totals\n")
  print.data.frame(x, row.names=FALSE)
}


# ============================================== Variance-Covariance matrix ====

# ----------------------------------------------------------------- extract ----

#' Extract variance-covariance matrix from TRIM output
#'
#' @param x TRIM output structure (i.e., output of a call to \code{trim})
#' @param which Selector to distinguish between variance-covariance based on the
#' imputed data (default), or the modelled data.
#'
#' @return a JxJ matrix, where J is the number of time points.
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

#' Extract coefficients of the reparameterisation of trim model 3
#'
#'
#'
#' @param x an object of class \code{\link{trim}}
#'
#' @return a list of class \code{trim.linear} with elements \code{trend}, 
#'   containing additive and multiplicative parameters of the overall linear
#'   trend, and \code{dev}, containing the deviations from that trend.
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




