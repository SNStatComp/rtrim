# ########################################################### Overall slope ####

#' Compute overall slope
#'
#' The overal slope represents the total growth over the piecewise linear model.
#'
#'
#' @param x an object of class \code{\link{trim}}.
#' @param which \code{[character]} Choose between \code{"imputed"} or
#'   \code{"model"} counts.
#' @param cp \code{[numeric]} Change points for which to compute the overall slope.
#'
#' @section Details:
#' 
#' The overall slope represents the mean growth or decline over a period of time.
#' This can be determined over the whole time period for which the modelis fitted (this is the default)
#' or may be computed over time slices that can be defined with the \code{cp} parameter.
#' The values for \code{cp} do not depend on \code{changepoints} that were used when 
#' specifying the \code{trim} model (See also the example below).
#'
#'
#' @return a list of class \code{trim.overall} containing overall slope
#'   coefficients (\code{coef}), the p-value of the overall slope (\code{p}),
#'   and the size-effect (\code{effect}).
#' @export
#'
#' @family analyses
#' @examples
#' 
#' # obtain the overall slope accross all change points.
#' data(skylark)
#' z <- trim(count ~ time + site, data=skylark, model=2)
#' overall(z)
#' plot(overall(z))
#' 
#' # Overall is a list, you can get information out if it using the $ syntax,
#' # for example
#' L <- overall(z)
#' L$coef
#' 
#' # Obtain the slope from changepoint to changepoint
#' z <- trim(count ~ time + site, data=skylark, model=2,changepoints=c(1,4,6))
#' # slope from time point 1 to 5
#' overall(z,cp=c(1,5,7))
overall <- function(x, which=c("imputed","model"), cp=numeric(0)) {
  stopifnot(class(x)=="trim")
  which = match.arg(which)

  # extract vars from TRIM output
  tt_mod <- x$tt_mod
  tt_imp <- x$tt_imp
  var_tt_mod <- x$var_tt_mod
  var_tt_imp <- x$var_tt_imp
  ntime <- x$ntime

  # The overall slope is computed for both the modeled and the imputed $\Mu$'s.
  # So we define a function to do the actual work

  .compute.overall.slope <- function(tt, var_tt, src) {
    J <- length(tt)
    stopifnot(nrow(var_tt)==J && ncol(var_tt)==J)

    # Use Ordinary Least Squares (OLS) to estimate slope parameter $\beta$
    X <- cbind(1, seq_len(J)) # design matrix
    y <- matrix(log(tt))
    bhat <- solve(t(X) %*% X) %*% t(X) %*% y # OLS estimate of $b = (\alpha,\beta)^T$
    yhat <- X %*% bhat

    # Apply the sandwich method to take heteroskedasticity into account
    dvtt <- 1/tt # derivative of $\log{\Mu}$
    Om <- diag(dvtt) %*% var_tt %*% diag(dvtt) # $\var{log{\Mu}}$
    var_beta <- solve(t(X) %*% X) %*% t(X) %*% Om %*% X %*% solve(t(X) %*% X)
    b_err <- sqrt(diag(var_beta))

    # Compute the $p$-value, using the $t$-distribution
    df <- J - 2
    t_val <- bhat[2] / b_err[2]
    if (df>0) p <- 2 * pt(abs(t_val), df, lower.tail=FALSE)
    else      p <- NA

    # Also compute effect size as relative change during the monitoring period.
    effect <- abs(yhat[J] - yhat[1]) / yhat[1]

    # Reverse-engineer the SSR (sum of squared residuals) from the standard error
    j <- 1:J
    D <- sum((j-mean(j))^2)
    SSR <- b_err[2]^2 * D * (J-2)

    # Export the results
    df <- data.frame(
      add       = bhat,
      se_add    = b_err,
      mul       = exp(bhat),
      se_mul    = exp(bhat) * b_err,
      row.names = c("intercept","slope")
    )
    list(src=src, coef=df,p=p, effect=effect, J=J, tt=tt, err=sqrt(diag(var_tt)), SSR=SSR)
  }

  if (which=="imputed") {
    tt     <- tt_imp
    var_tt <- var_tt_imp
    src = "imputed"
  } else if (which=="model") {
    tt     <- tt_mod
    var_tt <- var_tt_mod
    src = "imputed"
  }

  if (length(cp)==0) {
    # Normal overall slope
    out <- .compute.overall.slope(tt, var_tt, src)
    out$type <- "normal" # mark output as 'normal' overall slope
  } else {
    # overall slope per changepoint
    J <- length(tt)
    ncp <- length(cp)
    cpx <- c(cp, J) # Extend list of overall changepoints with final year
    coef.collector = data.frame()
    for (i in 1:ncp) {
      idx <- cpx[i] : cpx[i+1]
      tmp <- .compute.overall.slope(tt[idx], var_tt[idx,idx], src)
      slope <- tmp$coef[2,] # slope is on second row of output dataframe
      prefix <- data.frame(from=x$time.id[cpx[i]], upto=x$time.id[cpx[i+1]])
      coef.collector <- rbind(coef.collector, cbind(prefix, slope))
    }
    out <- list(src=src,coef=coef.collector, type="changept") #store 'n' mark
  }

  structure(out, class="trim.overall")
}


# ------------------------------------------------------------------- Print ----

#' Print an object of class trim.overall
#'
#' @param x An object of class \code{trim.overall}
#'
#' @export
#' @keywords internal
print.trim.overall <- function(x,...) {
  if (x$type=="normal") {
    # Compute 95\% confidence interval of multiplicative slope
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
    info <- sprintf("p=%f, conf.int (mul)=[%f, %f], change=%.2f%%", x$p, blo, bhi, 100*change)
    printf("Overall slope (%s): %s\n", x$src, info)
    print(x$coef, row.names=TRUE)
  } else if (x$type=="changept") {
    printf("Overall slope (%s) for changepoints\n", x$src)
    print(x$coef, row.names=FALSE)
  } else stop("Can't happen")
}

#--------------------------------------------------------------------- Plot ----

#' Plot overall slope
#'
#' Creates a plot of the overall slope, its 95\% confidence band, the
#' total population per time and their 95\% confidence intervals.
#' 
#' @param x An object of class \code{trim.overall} (returned by \code{\link{overall}})
#' @param imputed Toggle to show imputed counts
#' @param ... Further options passed to \code{\link[graphics]{plot}}
#'
#' @family analyses
#' 
#' @examples 
#' data(skylark)
#' m <- trim(count ~ time + site, data=skylark, model=2)
#' plot(overall(m))
#' 
#' @export
plot.trim.overall <- function(x, imputed=TRUE, ...) {
  X <- x
  title <- if (is.null(list(...)$main)){
    attr(X, "title")
  } else {
    list(...)$main
  }

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
  plot(xrange, yrange, type='n', xlab="Time point", ylab="Count", las=1, main=title,...)
  polygon(xconf, yconf, col=gray(0.9), lty=0)
  lines(x, ytrend, col=cbred, lwd=3)
  segments(j,y0, j,y1, lwd=3, col=gray(0.5))
  points(j, ydata, col=cbblue, type='b', pch=16, lwd=3)

}
