# TRIM postprocessing functions

#===============================================================================
#                                                                Goodness-of-fit
#===============================================================================

#-------------------------------------------------------------------------------
#                                                                        extract


#' Extract TRIM goodness-of-fit information.
#'
#' TRIM computes three goodness-of-fit measures:
#' \itemize{
#'   \item Chi-squared
#'   \item Likelihood ration
#'   \item Akaike Information content
#' }
#'
#' @param x TRIM output structure (i.e., output of a call to \code{trim})
#'
#' @return a list of type "trim.gof", containing elements \code{chi2}, \code{LR}
#' and \code{AIC}, for Chi-squared, Likelihoof Ratio and Akaike information content,
#' respectively.
#' @export
#'
#' @family analyses
#' 
#' @examples
#' data(skylark)
#' z <- trim(skylark, count ~ time + site, model=2)
#' # prettyprint GOF information
#' gof(z)
#' 
#' # get individual elements, e.g. p-value
#' L <- gof(z)
#' LR_p <- L$LR$p # get p-value for likelihood ratio
#' 
gof <- function(x) UseMethod("gof")

#' @export
#' @rdname gof
gof.trim <- function(x) {
  structure(list(chi2 = x$chi2, LR=x$LR, AIC=x$AIC), class="trim.gof")
}

#-------------------------------------------------------------------------------
#                                                                          print

#' Print method for \code{trim.gof}
#'
#' @export
#' @param x a \code{trim.gof} object
#' @keywords internal
print.trim.gof <- function(x,...) {
   
  # print welcome message
  cat(sprintf("GOODNESS OF FIT\n"))

  # print $\Chi^2$ results
  with(x$chi2,
       printf("%24s =%.2f, df=%d, p=%.4f\n", "Chi-square", chi2, df, p))

  # print Likelihood ratio
  with(x$LR,
       printf("%24s =%.2f, df=%d, p=%.4f\n", "Likelihood Ratio", LR, df, p))

  # Akaike
  with(x,
       printf("%24s =%.2f\n", "AIC (up to a constant)", AIC))
}


#===============================================================================
#                                                                        Summary
#===============================================================================

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


#===============================================================================
#                                                                   Coefficients
#===============================================================================

#-------------------------------------------------------------------------------
#                                                                        Extract

#' Extract TRIM model coefficients.
#'
#' @param object TRIM output structure (i.e., output of a call to \code{trim})
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
coef.trim <- function(object, which=c("additive","multiplicative","both")) {
  stopifnot(class(object)=="trim")

  # Craft a custom output
  which <- match.arg(which)
  if (which=="additive") {
    if (object$use.changepoints) out <- cbind(object$coef$int, object$coef$add)
    else                    out <- object$coef$add
  } else if (which=="multiplicative") {
    if (object$use.changepoints) out <- cbind(object$coef$int, object$coef$mul)
    else                    out <- object$coef$mul
  } else if (which=="both") {
    if (object$use.changepoints) out <- cbind(object$coef$int, object$coef$add, object$coef$mul)
    else                    out <- cbind(object$coef$add, object$coef$mul)
  } else stop(sprintf("Invalid options which=%s", which))

  structure(list(model=object$model, coef=out), class="trim.coef")
}

#-------------------------------------------------------------------------------
#                                                                          Print

#' Print method for \code{trim.coef}
#'
#' @export
#' @param x An object of class \code{trim.coef}
#' @keywords internal
print.trim.coef <- function(x,...) {
  stopifnot(class(x)=="trim.coef")

  printf("Parameter estimates\n")
  if (x$model==2) {
    print(x$coef)
  } else if (x$model==3) {
    printf("\nParameters for each time point\n")
    print(x$coef, row.names=FALSE)
  } else stop("Can't happen")

}


#===============================================================================
#                                                                    Time totals
#===============================================================================

#-------------------------------------------------------------------------------
#                                                                        Extract

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
  which <- match.arg(which)
  
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

#-------------------------------------------------------------------------------
#                                                                         Export

export <- function(x, species, stratum) UseMethod("export")

export.trim.totals <- function(x, species, stratum) {
  stopifnot(class(x)=="trim.totals")

  # Create extra columns to be put before the actual time totals
  df1 = data.frame(species=species, stratum=stratum)
  df2 = x$totals
  df = cbind(df1, df2)
  print(df, row.names=FALSE)
}

#-------------------------------------------------------------------------------
#                                                                         Print

print.trim.totals <- function(x,...) {
  printf("Time totals\n")
  print(x$totals, row.names=FALSE)
}


#===============================================================================
#                                                                   Time indices
#===============================================================================

#-------------------------------------------------------------------------------
#                                                                        Extract

#' Extract time-indices from TRIM output
#'
#' @param x TRIM output structure (i.e., output of a call to \code{trim})
#' @param which Selector to distinguish between time indices based on the imputed data (default),
#' the modelled data, or both.
#'
#' @return a structure of class \code{trim.index}, which is a \code{data.frame}
#' with an extra class label.
#' 
#' @export
#'
#' @family analyses
#' @examples
#' data(skylark)
#' z <- trim(count ~ time + site, data=skylark, model=2)
#' index(z) 
#' # mimic classic TRIM:
#' index(z, "both") 
#' # Extract standard error for the imputed data
#' SE <- index(z)$std.err 
index <- function(x, which=c("imputed","model","both")) {
  stopifnot(inherits(x,"trim"))

  # Select output columns from the pre-computed time totals
  which <- match.arg(which)
  idx <- switch(which
    , model   = x$time.index[c(1,2,3)]
    , imputed = x$time.index[c(1,4,5)]
    , both    = x$time.index
  )

  # wrap the time.index field in a list and make it an S3 class
  # (so that it becomes recognizable as a TRIM time-indices)
  class(idx) <- c("trim.index","data.frame")
  idx
}

#-------------------------------------------------------------------------------
#                                                                         Export

export <- function(x, species, stratum) UseMethod("export")

export.trim.index <- function(x, species, stratum) {

  # Create extra columns to be put before the actual time indices
  df1 = data.frame(species=species, stratum=stratum)
  df2 = x$idx
  df = cbind(df1, df2)
  print(df, row.names=FALSE)
}

#-------------------------------------------------------------------------------
#                                                                         Print

#' Print an object of class trim.index
#' 
#' @param x An object of class \code{trim.index}
#' 
#' @export
#' @keywords internal
print.trim.index <- function(x,...) {
  printf("Time indices\n")
  print.data.frame(x, row.names=FALSE)
}

#===============================================================================
#                                                  Reparameterisation of Model 3
#===============================================================================

#-------------------------------------------------------------------------------
#                                                                        extract

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

#-------------------------------------------------------------------------------
#                                                                          print

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

#===============================================================================
#                                                                      Wald test
#===============================================================================

#-------------------------------------------------------------------------------
#                                                                        Extract

#' Extract information of Wald tests as peformed by TRIM
#'
#' @param x TRIM output structure (i.e., output of a call to \code{trim})
#'
#' @return a model-dependent list of Wald statistics
#' @export
#'
#' @family analyses
#' @examples
#' data(skylark)
#' z2 <- trim(count ~ time + site, data=skylark, model=2)
#' # print info on significance of slope parameters
#' wald(z2)  
#' z3 <- trim(count ~ time + site, data=skylark, model=3)
#' # print info on significance of deviations from linear trend
#' wald(z3)
wald <- function(x) {
  structure(x$wald, class="trim.wald")
}

#-------------------------------------------------------------------------------
#                                                                          Print

#' Print an object of class trim.wald
#' 
#' @param x An object of class \code{trim.wald}
#' 
#' @export
#' @keywords internal
print.trim.wald <- function(x,...) {
  if (x$model==2) {
    if (length(x$W)==1) {
      printf("Wald test for significance of slope parameter\n")
      printf("  Wald = %.2f, df=%d, p=%f\n", x$W, x$df, x$p)
    } else {
      printf("Wald test for significance of changes in slope\n")
      df = data.frame(Changepoint = 1:length(x$W),
                      Wald_test = x$W,
                      df = x$df,
                      p = x$p)
      print(df, row.names=FALSE)
    }
  } else if (x$model==3) {
    printf("Wald test for significance of deviations from linear trend\n")
    printf("  Wald = %.2f, df=%d, p=%f\n", x$W, x$df, x$p)
  } else stop("Can't happen")

}

#===============================================================================
#                                                                  Overall slope
#===============================================================================


#-------------------------------------------------------------------------------
#                                                                        Extract

#' Extract information of TRIM overall slope
#'
#' @param x TRIM output object
#' @param which Choose between "imputed" or "model" counts
#'
#' @return a list containing information overall slope coefficients (\code{coef}),
#'   the p-value of the overall slope (\code{p}),
#'   the size-effect (\code{effect}),
#' @export
#'
#' @family analyses
#' @examples
#' data(skylark)
#' z <- trim(count ~ time + site, data=skylark, model=2)
#' overall(z)
overall <- function(x, which=c("imputed","model")) {
  stopifnot(class(x)=="trim")
  which <- match.arg(which)
  if (which=="imputed") {
     structure(x$overall$imp, class="trim.overall", src="imputed")
  } else {
     structure(x$overall$mod, class="trim.overall", src="model")
  }
}

#-------------------------------------------------------------------------------
#                                                                          Print


#' Print an object of class trim.overall
#' 
#' @param x An object of class \code{trim.overall}
#' 
#' @export
#' @keywords internal
print.trim.overall <- function(x,...) {
  stopifnot(class(x)=="trim.overall")

  # Compute 95% confidence interval of multiplicative slope
  bhat <- x$coef[[3]][2] # multiplicative trend (i.e., not log-transformed)
  berr <- x$coef[[4]][2] # corresponding standard error
  alpha <- 0.05
  df <- x$J-2
  tval <- qt((1-alpha/2), df)
  blo <- bhat - tval * berr
  bhi <- bhat + tval * berr
  # Compute effect size
  change <- bhat ^ (x$J-1) - 1
  # Build an informative string
  src = attr(x, "src")
  info <- sprintf("p=%f, conf.int (mul)=[%.f, %f], change=%.2f%%", x$p, blo, bhi, 100*change)
  printf("Overall slope (%s): %s\n", src, info)
  print(x$coef, row.names=TRUE)}

#-------------------------------------------------------------------------------
#                                                                           Plot

#' Plot overall slope
#'
#' @param x An object of class \code{trim.overall} (returned by \code{\link{overall}})
#' @param imputed Toggle to show imputed counts 
#' @param ... Further options passed to \code{\link[graphics]{plot}}
#'
#' @family analyses
#' @export
plot.trim.overall <- function(x, imputed=TRUE, ...) {
  X <- x
  title <- attr(X, "title")

  J <- X$J

  # Collect all data for plotting: time-totals
  j <- 1:J
  ydata <- X$tt

  # error bars
  y0 = ydata - X$err
  y1 = ydata + X$err

  # Trend line
  a <- X$coef[[1]][1] # intercept
  b <- X$coef[[1]][2] # slope
  x <- seq(1, J, length.out=100)
  ytrend <- exp(a + b*x)

  # Confidence band
  xconf <- c(x, rev(x))
  alpha <- 0.05
  df <- J - 2
  t <- qt((1-alpha/2), df)
  dx2 <- (x-mean(j))^2
  sumdj2 <- sum((j-mean(j))^2)
  dy <- t * sqrt((X$SSR/(J-2))*(1/J + dx2/sumdj2))
  ylo <- exp(a + b*x - dy)
  yhi <- exp(a + b*x + dy)
  yconf <- c(ylo, rev(yhi))

  # Compute the total range of all plot elements
  xrange = range(x)
  yrange = range(range(yconf), range(y0), range(y1))

  # Now plot layer-by-layer
  cbred <- rgb(228,26,28, maxColorValue = 255)
  cbblue <- rgb(55,126,184, maxColorValue = 255)
  plot(xrange, yrange, type='n', xlab="Time point", ylab="Count"
       , main=title,las=1,...)
  polygon(xconf, yconf, col=gray(0.9), lty=0)
  lines(x, ytrend, col=cbred, lwd=3)
  segments(j,y0, j,y1, lwd=3, col=gray(0.5))
  points(j, ydata, col=cbblue, type='b', pch=16, lwd=3)

}

#===============================================================================
#                                                                           Plot
#===============================================================================


# Plotting of data is easier if we convert it to a data table.
# Here is a function that does the conversion.
mat2df <- function(m, src=NA) {
  nsite <- nrow(m)
  ntime <- ncol(m)
  df <- data.frame(
    Site  = factor(rep(1:nsite, times=ntime)),
    Time  = rep(1:ntime, each=nsite),
    Count = as.vector(m)
  )
  df <- df[order(df$Site, df$Time), ]   # Sort by site, then by time
  if (!is.na(src)) df$Source=src        # set optional data source
  return(df)
}

