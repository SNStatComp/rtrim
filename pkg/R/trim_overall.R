# ########################################################### Overall slope ####

#' Compute overall slope
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
  which = match.arg(which)

  # extract vars from TRIM output
  tt_mod <- x$tt_mod
  tt_imp <- x$tt_imp
  var_tt_mod <- x$var_tt_mod
  var_tt_imp <- x$var_tt_imp
  ntime <- x$ntime

  # The overall slope is computed for both the modeled and the imputed $\Mu$'s.
  # So we define a function to do the actual work

  .compute.overall.slope <- function(tt, var_tt) {
    stopifnot(length(tt)==ntime)
    J <- ntime

    # Use Ordinary Least Squares (OLS) to estimate slope parameter $\beta$
    X <- cbind(1, seq_len(ntime)) # design matrix
    y <- matrix(log(tt))
    bhat <- solve(t(X) %*% X) %*% t(X) %*% y # OLS estimate of $b = (\alpha,\beta)^T$
    yhat <- X %*% bhat

    # Apply the sandwich method to take heteroskedasticity into account
    dvtt <- 1/tt_mod # derivative of $\log{\Mu}$
    Om <- diag(dvtt) %*% var_tt %*% diag(dvtt) # $\var{log{\Mu}}$
    var_beta <- solve(t(X) %*% X) %*% t(X) %*% Om %*% X %*% solve(t(X) %*% X)
    b_err <- sqrt(diag(var_beta))

    # Compute the $p$-value, using the $t$-distribution
    df <- ntime - 2
    t_val <- bhat[2] / b_err[2]
    p <- 2 * pt(abs(t_val), df, lower.tail=FALSE)

    # Also compute effect size as relative change during the monitoring period.
    effect <- abs(yhat[J] - yhat[1]) / yhat[1]

    # Reverse-engineer the SSR (sum of squared residuals) from the standard error
    j <- 1:J
    D <- sum((j-mean(j))^2)
    SSR <- b_err[2]^2 * D * (J-2)

    # Export the results
    df <- data.frame(
      Additive       = bhat,
      std.err.       = b_err,
      Multiplicative = exp(bhat),
      std.err.       = exp(bhat) * b_err,
      row.names      = c("Intercept","Slope"),
      check.names    = FALSE
    )
    list(coef=df,p=p, effect=effect, J=J, tt=tt, err=x$time.totals[[3]], SSR=SSR)
  }

  if (which=="imputed") {
    out = .compute.overall.slope(tt_imp, var_tt_imp)
    out$src = "imputed"
  } else if (which=="model") {
    out = .compute.overall.slope(tt_mod, var_tt_mod)
    out$src = "model"
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
  info <- sprintf("p=%f, conf.int (mul)=[%.f, %f], change=%.2f%%", x$p, blo, bhi, 100*change)
  printf("Overall slope (%s): %s\n", x$src, info)
  print(x$coef, row.names=TRUE)
}

#--------------------------------------------------------------------- Plot ----

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
  plot(xrange, yrange, type='n', xlab="Time point", ylab="Count", main=title)
  polygon(xconf, yconf, col=gray(0.9), lty=0)
  lines(x, ytrend, col=cbred, lwd=3)
  segments(j,y0, j,y1, lwd=3, col=gray(0.5))
  points(j, ydata, col=cbblue, type='b', pch=16, lwd=3)

}
