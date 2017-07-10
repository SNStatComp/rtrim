# ########################################### TRIM postprocessing functions ####

# ================================================================= Print ======


#' print a 'trim' object
#'
#' @param x a \code{\link{trim}} object
#' @param ... currently unused
#'
#' @export
#' @keywords internal
print.trim <- function(x,...){
  cat("Call:\n")
  print(x$call)
  cat("\n",x$convergence,"\n")
  cat("\nCoefficients:\n")
  print(coef.trim(x))
}

# ================================================================= Summary ====

# ----------------------------------------------------------------- extract ----


#' Summary information for a TRIM job
#'
#' Print a summary of a \code{\link{trim}} object.
#'
#' @param object an object of class \code{\link{trim}}.
#' @param ... Currently unused
#'
#' @return A \code{list} of class \code{trim.summary} containing the call that
#'   created the object, the model code, the coefficients (in additive and
#'   multiplicative form) , the goodness of fit parameters,the overdispersion
#'   and the serial correlation parameters (if computed).
#'
#' @export
#'
#' @family analyses
#' @seealso \code{\link{trim}}
#' @examples
#'
#' data(skylark)
#' z <- trim(count ~ site + time, data=skylark, model=2, overdisp=TRUE)
#'
#' summary(z)
summary.trim <- function(object,...) {

  structure(list(
    call = object$call
    , coefficients = coef.trim(object)
    , gof = gof(object)
    , overdispersion = overdispersion(object)
    , serialcorrelation = serial_correlation(object)
    , model = object$model
    , method = object$method
    , convergence = object$convergence
  ),class="trim.summary")
}


#' @export
#' @keywords internal
print.trim.summary <- function(x,...){

  cl <- paste(capture.output(print(x$call)),collapse="\n")
  printf("Call:\n%s\n",cl)
  printf("\n")

  printf("Model  : %d\n", x$model)
  printf("Method : %s (%s)\n", x$method, x$convergence)

  if (x$model>1) {
    printf("\nCoefficients:\n")
    print(x$coefficients)
    printf("\n")
  }

  printf(" Overdispersion     : %.4f\n",x$overdispersion)
  printf(" Serial Correlation : %.4f\n",x$serialcorrelation)
  printf("\n")

  print(x$gof)
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
#' refererred to as coefficients in the additive representation. By
#' exponentiating both sides of the above equations, alternative representations
#' can be written down. Explicitly, one can show that
#'
#' \itemize{
#' \item{Model 2: \eqn{\mu_{ij}= a_ib^{(j-1)} = b\mu_{ij-1}}, where \eqn{a_i=e^{\alpha_i}} and \eqn{b=e^\beta}.}
#' \item{Model 3: \eqn{\mu_{ij}=a_ic_j}, where \eqn{a_i=e^{\alpha_i}}, \eqn{c_1=1} and \eqn{c_j=e^{\gamma_j}} for \eqn{j>1}.}
#' }
#'
#' The parameters \eqn{a_i}, \eqn{b} and \eqn{c_j} are referred to as
#' coefficients in the \emph{multiplicative form}.
#'
#' @section Trend and deviation (Model 3 only):
#'
#' The equation for Model 3
#'
#' \eqn{\ln\mu_{ij}  = \alpha_i + \gamma_j},
#'
#' can also be written as an overall slope resulting from a linear regression of
#' the \eqn{\mu_{ij}} over time,  plus site- and time effects that
#' record deviations from this overall slope.  In such a reparametrisation
#' the previous equation can be written as
#'
#' \eqn{\ln\mu_{ij} = \alpha_i^* + \beta^*d_j + \gamma_j^*,}
#'
#' where \eqn{d_j} equals \eqn{j} minus the mean over all \eqn{j} (i.e. if \eqn{j=1,2,\ldots,J}
#' then \eqn{d_j = j-(J+1)/2}). It is not hard to show that
#' \itemize{
#' \item{The \eqn{\alpha_i^*} are the mean \eqn{\ln\mu_{ij}} per site}
#' \item{The \eqn{\gamma_j^*} must sum to zero.}
#' }
#' The coefficients \eqn{\alpha_i^*} and \eqn{\gamma_j^*} are obtained by
#' setting \code{representation="deviations"}. If \code{representation="trend"},
#' the overall trend parameters \eqn{\beta^*} and \eqn{\alpha^*} from the overall
#' slope defined by \eqn{\alpha^* + \beta^*d_j} is returned.
#'
#' Finally, note that both the overall slope and the deviations can be written
#' in multiplicative form as well.
#'
#'
#' @param object TRIM output structure (i.e., output of a call to \code{trim})
#' @param representation \code{[character]} Choose the coefficient
#'   representation. Options \code{"trend"} and \code{"deviations"} are for model 3 only.
#' @param ... currently unused
#'
#' @return A \code{data.frame} containing coefficients and their standard errors,
#' both in additive and multiplicative form.
#'
#' @export
#'
#' @family analyses
#' @examples
#' data(skylark)
#' z <- trim(count ~ site + time, data=skylark, model=2, overdisp=TRUE)
#' coefficients(z)
coef.trim <- function(object,
    representation=c("standard","trend","deviations"),...) {


  representation <- match.arg(representation)

  if (representation %in% c("deviations","trend") && object$model != 3){
    stop(
      sprintf("Cannot extract  %s from TRIM model %d\n",representation,object$model)
      , call.=TRUE)
  }

  switch(representation
    , "standard" = object$coefficients
    , "deviations" = setNames(object$deviations,c("time","add","se_add","mul","se_mul"))
    , "trend" = setNames(object$linear.trend,c("add","se_add","mul","se_mul"))
  )

}


# ============================================================= Time totals ====

# ----------------------------------------------------------------- Extract ----

#' Extract time-totals from TRIM output
#'
#' @param x TRIM output structure (i.e., output of a call to \code{trim})
#' @param which select what totals to compute (see \code{Details} section).
#'
#' @return A \code{data.frame} with subclass \code{trim.totals}
#'  (for pretty-printing). The columns are \code{time}, \code{fitted}
#'  and \code{se_fit} (for standard error), and/or \code{imputed}
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
#' \item{\code{"fitted"}: Time totals are computed after replacing both missing values and observed values with
#' values predicted by the model.}
#' }
#'
#' @export
#'
#' @family analyses
#' @examples
#' data(skylark)
#' z <- trim(count ~ site + time, data=skylark, model=2, changepoints=c(3,5))
#' totals(z)
#'
#' totals(z, "both") # mimics classic TRIM
#'
totals <- function(x, which=c("imputed","fitted","both")) {
  stopifnot(class(x)=="trim")

  # Select output columns from the pre-computed time totals
  which <- match.arg(which)
  totals <- switch(which
    , fitted  = x$time.totals[c(1,2,3)]
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

#------------------------------------------------------------------ Plot -----

plot.trim.totals <- function(t1, ..., leg.pos="topleft") {

  # Create custom palette based on Color Brewer Set 1
  brewer_set1 <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#999999")
  opaque <- brewer_set1
  aqua   <- brewer_set1
  for (i in 1:9) aqua[i] <- adjustcolor(aqua[i], 0.2)


  # Build a list of time-totals with optional titles
  tt = list(t1)
  optional = list(...)

  # cat("tt pre:\n")
  # str(tt)

  # cat("optional:\n")
  # str(optional)

  nopt = length(optional)
  for (i in seq_len(nopt)) {
    x = optional[[i]]
    if ("character" %in% class(x)) {
      attr(tt[[length(tt)]], "tag") <- x
    } else if ("trim.totals" %in% class(x)) {
      tt[[length(tt)+1]] <- x
    } else {
      stop(sprintf("Invalid data type for optional argument %d: %s", i, class(x)))
    }
  }

  # cat("tt post:\n")
  # str(tt)

  # cat("leg.pos:\n")
  # str(leg.pos)

  # First pass to compute total range
  n = length(tt)
  for (i in 1:n) {
    x = tt[[i]][[1]] # Time point or years
    y = tt[[i]][[2]] # imputed or fitted
    s = tt[[i]][[3]] # Standard error
    ylo = y-s
    yhi = y+s
    if (i==1) {
      xrange <- range(x)
      yrange <- range(ylo, yhi)
    } else {
      xrange <- range(xrange, range(x))
      yrange <- range(yrange, range(ylo, range(yhi)))
    }
  }

  # Ensure y-axis starts at 0.0
  yrange <- range(0.0, yrange)

  # empty plot for correct axes
  plot(xrange, yrange, type='n', xlab="Time point", ylab="Time totals")

  # Second pass: plot them
  for (i in 1:n) {
    x = tt[[i]][[1]] # Time point or years
    y = tt[[i]][[2]] # imputed or fitted
    s = tt[[i]][[3]] # Standard error
    ylo = y-s
    yhi = y+s

    xx = c(x, rev(x))
    ci = c(ylo, rev(yhi))

    polygon(xx,ci, col=aqua[i], border=NA)
    lines(x,y, col=opaque[i])
  }

  # third pass: legend
  nnamed  = 0
  nnoname = 0
  for (i in 1:n) {
    s <- attr(tt[[i]],"tag")
    if (is.null(s)) {
      nnoname <- nnoname + 1
      s <- sprintf("<unnamed> %d", nnoname)
    } else {
      nnamed = nnamed + 1
    }
    if (i==1) {
      leg.colors <- opaque[i]
      leg.names  <- s
    } else {
      leg.colors <- c(leg.colors, opaque[i])
      leg.names <- c(leg.names, s)
    }
  }
  if (n>1 | nnamed>0) {
    legend(leg.pos, legend=leg.names, col=leg.colors, lty=1, lwd=2, bty='n', inset=0.02, y.intersp=1.5);
  }
}

# ============================================== Variance-Covariance matrix ====

# ----------------------------------------------------------------- extract ----

#' Extract variance-covariance matrix from TRIM output
#'
#' @param object TRIM output structure (i.e., output of a call to \code{trim})
#' @param which \code{[character]} Selector to distinguish between variance-covariance based on the
#' imputed data (default), or the modelled data.
#' @param ... Arguments to pass to or from other methods (currently unused)
#'
#' @return a JxJ matrix, where J is the number of time points.
#' @export
#'
#' @family analyses
#' @examples
#' data(skylark)
#' z <- trim(count ~ site + time, data=skylark, model=2);
#' totals(z)
#' vcv1 <- vcov(z)       # Use imputed data
#' vcv2 <- vcov(z,"model") # Use modelled data
vcov.trim <- function(object, which = c("imputed","model"), ... ) {
  stopifnot(inherits(object,"trim"))
  which <- match.arg(which)

  vcv <- switch(which
    , model   = object$var_tt_mod
    , imputed = object$var_tt_imp
  )
  vcv
}




# ================================================================= Results ====

# Function \verb!result()! collects and combines the observed, modelled, and imputed
# counts. These results are presented as a data frame, which is readily exported to
# a file by the user.

#' collect observed, modelled, and imputed counts from TRIM output
#'
#' @param z TRIM output structure (i.e., output of a call to \code{trim})
#'
#' @return A \code{data.frame}, one row per site-time combination, with columns for
#' site, time, observed counts, modelled counts and imputed counts.
#' Missing observations are marked as \code{NA}.
#'
#' @export
#'
#' @family analyses
#' @examples
#' data(skylark)
#' z <- trim(count ~ site + time, data=skylark, model=2);
#' out <- results(z)
results <- function(z) {
  stopifnot(inherits(z,"trim"))

  out <- data.frame(
    site = rep(z$site.id, each=z$ntime),
    time = rep(z$time.id, times=z$nsite),
    observed = as.vector(t(z$f)),
    fitted   = as.vector(t(z$mu)),
    imputed  = as.vector(t(z$imputed))
  )
  class(out) <- c("trim.results","data.frame")
  out
}

plot.trim.results <- function(z, ...) {
  sites = levels(z$site)
  nsite = nlevels(z$site)
  hues = seq(0, 360, length.out = nsite+1)[1:nsite]
  colors = hcl(hues, 100, 65) # C and L Similar to ggplot
  # hues = seq(0, 1, length.out = nsite+1)[1:nsite]
  # colors = hsv(hues, 0.5, 1)
  xrange = range(z$time)
  yrange = range(z$observed, z$modelled, na.rm=TRUE)
  plot(xrange,yrange, type='n', xlab="Time", ylab="Counts")
  for (i in 1: nsite) {
    df = z[z$site == sites[i], ]
    points(df$time, df$observed, pch=16, col=colors[i])
    lines(df$time, df$modelled, col=colors[i])
  }
}

# ================================================================== Advice ====

#' Give advice on further refinement of TRIM models
#'
#' @param z an object of class \code{\link{trim}}.
#'
#' @export
#'
#' @family analyses
#' @seealso \code{\link{trim}}
#' @examples
#'
#' data(skylark)
#' z <- trim(count ~ site + time, data=skylark, model=2, overdisp=TRUE)
#' now_what(z)
now_what <- function(z) {
  stopifnot(inherits(z,"trim"))
  Wald <- wald(z)
  advice_given <- FALSE

  if (!is.null(Wald$dslope)) {
    p = Wald$dslope$p
    if (any(p > 0.2)) {
      ntot = length(p)
      ndel = sum(p > 0.2)
      worst = which.max(p)
      printf("%d out of %d changepoints appear to be insignificant;", ndel, ntot)
      printf("changepoint #%d would be the first candidate to remove.\n", worst)
      advice_given <- TRUE
    }
  }

  if (z$model==3) {
    LR_p <- gof(z)$LR$p
    if (LR_p < 0.05) {
      printf("Model 3 has a bad fit (%g < 0.05); Try a different model.\n", LR_p)
      advice_given <- TRUE
    }
  }

  if (!advice_given) rprintf("Model appears to be adequate; no suggestions for further improvement.\n")
}

