# ########################################################### Overall slope ####

#' Compute overall slope
#'
#' @param x TRIM output object
#' @param which Choose between "imputed" or "model" counts
#' @param cp \code{[numeric]} Change points for which to compute the overall slope.
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
overall <- function(x, which=c("imputed","model"), cp=numeric(0)) {
  stopifnot(class(x)=="trim")
  which = match.arg(which)
  #browser()

  # extract vars from TRIM output
  tt_mod <- x$tt_mod
  tt_imp <- x$tt_imp
  var_tt_mod <- x$var_tt_mod
  var_tt_imp <- x$var_tt_imp
  ntime <- x$ntime

  # if (base>0) { # use index instead
  #     browser()
  #     tt     = tt_mod
  #     var_tt = var_tt_mod
  #     b = base
  #
  #     tau <- tt / tt[b]
  #     J <- length(tt)
  #     var_tau <- numeric(J)
  #     for (j in 1:J) {
  #       d <- matrix(c(-tt[j] / tt[b]^2, 1/tt[b]))
  #       V <- var_tt[c(b,j), c(b,j)]
  #       var_tau[j] <- t(d) %*% V %*% d
  #     }
  #
  #     tt_mod <- tt_imp <- tau
  #     var_tt_mod <- var_tt_imp <- diag(var_tau)
  # }


  # The overall slope is computed for both the modeled and the imputed $\Mu$'s.
  # So we define a function to do the actual work

  .compute.overall.slope <- function(tt, var_tt, src) {
    # tpt = time points, either 1..J or year1..yearn
    J <- length(tt)
    stopifnot(nrow(var_tt)==J && ncol(var_tt)==J)

    # handle zero time totals (might happen in imputed TT)
    problem = tt<1e-6
    log_tt = log(tt)
    log_tt[problem] <- 0.0
    alt_tt <- exp(log_tt)

    # Use Ordinary Least Squares (OLS) to estimate slope parameter $\beta$
    X <- cbind(1, 1:J) # design matrix
    y <- matrix(log_tt)
    #y[tt<1e-6] = 0.0 # Handle zero (or very low) counts
    bhat <- solve(t(X) %*% X) %*% t(X) %*% y # OLS estimate of $b = (\alpha,\beta)^T$
    yhat <- X %*% bhat

    # Apply the sandwich method to take heteroskedasticity into account
    dvtt <- 1/alt_tt # derivative of $\log{\Mu}$
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
    src = "model"
  }

  if (length(cp)==0) {
    # Normal overall slope
    out <- .compute.overall.slope(tt, var_tt, src)
    out$type <- "normal" # mark output as 'normal' overall slope
    out$timept <- x$time.id # export tiem points for proper plotting
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
  tpt = X$timept

  # Collect all data for plotting: time-totals
  j <- 1:J
  ydata <- X$tt

  # error bars
  y0 = ydata - X$err
  y1 = ydata + X$err

  # Trend line
  a <- X$coef[[1]][1] # intercept
  b <- X$coef[[1]][2] # slope
  x <- seq(1, J, length.out=100) # continue timepoint 1..J
  ytrend <- exp(a + b*x)

  # Confidence band
  xcont <- seq(min(tpt), max(tpt), len=100) # continue year1..yearn
  xconf <- c(xcont, rev(xcont))
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
  xrange = range(xcont)
  yrange = range(range(yconf), range(y0), range(y1))

  # Now plot layer-by-layer (using ColorBrewer colors)
  cbred <- rgb(228,26,28, maxColorValue = 255)
  cbblue <- rgb(55,126,184, maxColorValue = 255)
  plot(xrange, yrange, type='n', xlab="Time point", ylab="Count", main=title)
  polygon(xconf, yconf, col=gray(0.9), lty=0)
  lines(xcont, ytrend, col=cbred, lwd=3)
  segments(tpt,y0, tpt,y1, lwd=3, col=gray(0.5))
  points(tpt, ydata, col=cbblue, type='b', pch=16, lwd=3)
}
