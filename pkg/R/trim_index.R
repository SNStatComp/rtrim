# ################################################################# Indices ####


# ============================================= Internal workhorse function ====

.old_index <- function(tt, var_tt, b) {

  # Time index $\tau_j$ is defined as time totals $\Mu$, normalized by the time total for the base
  # year, i.e.\,
  # $$ \tau_j = \Mu_j / \Mu_b $$
  # where $b\in[1\ldots J]$ indicates the base year.
  tau <- tt / tt[b]

  # Don't attempt to compute the index variance if time total variance is not given
  if (is.null(var_tt)) return(list(tau=tau, var_tau=NULL))

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
    if (j==b) {
      # SE in the base year is always 0. In principle this is also the result
      # from the algorihtm below applied to $j=b$, but in practice this may result
      # in very small negative values, causing NaN's when send to sqrt().
      var_tau[j] <- 0.0
    } else {
      d <- matrix(c(-tt[j] / tt[b]^2, 1/tt[b]))
      V <- var_tt[c(b,j), c(b,j)]
      var_tau[j] <- t(d) %*% V %*% d
    }
  }

  out <- list(tau=tau, var_tau=var_tau)
  out
}


.index <- function(tt, var_tt, base, level=NULL, sig2=NULL, formal=TRUE) {

  nbase <- length(base)

  if (nbase==1) {
    scale <- 1 / tt[base]
  } else {
    scale <- 1 / mean(tt[base])
  }
  tau <- tt * scale

  J <- length(tt)
  var_tau <- numeric(J)
  d <- matrix(0, nbase+1, 1)
  for (j in 1:J) {

    if (nbase==1 && !formal) {
      # compute scaled standard error
      #var_tau[j] <- var_tt[j,j] / tt[base]^2
      SE <- sqrt(var_tt[j,j])
      SE_scaled  <- SE * scale
      var_tau[j] <- SE_scaled^2
      next
    }

    if (nbase==1 && j==base) {
      var_tau[j] <- 0.0
    } else {
      d[1:nbase] <- -nbase * tt[j] / sum(tt[base])^2
      d[nbase+1] <- nbase / sum(tt[base])
      idx <- c(base, j)
      V <- var_tt[idx,idx]
      var_tau[j] <- t(d) %*% V %*% d
    }
  }

  out <- list(tau=tau, var_tau=var_tau)

  # Optionally add confidence interval
  if (!is.null(level)) {
    mul <- ci_multipliers(lambda=tt, sig2=sig2, level=level)
    out$lo_tau <- tau - mul$lo * sqrt(var_tau)
    out$hi_tau <- tau + mul$hi * sqrt(var_tau)
  }

  out
}

# ========================================================== User interface ====

#' Extract time-indices from TRIM output.
#'
#' Indices are obtained by dividing the modelled or imputed time totals by a reference value.
#' Most commonly, the time totals for the starting year are used as reference.
#' As a result, the index value for this year will be 1.0, with a standard error of 0.0 by definition.\cr
#' Alternatively, a range of years can be used as reference. In this case, the mean time totals for this range will be used as
#' reference, and the standard errors will be larger than 0.0.\cr
#' Starting with \code{rtrim} 2.2, an additional method can be selected,
#' which uses a simpler scaling approach to standard errors of the indices
#'
#' @param x an object of class \code{\link{trim}}
#' @param which (character) Selector to distinguish between time indices based on the imputed data (default),
#' the fitted model, or both.
#' @param covars (logical) Switch to compute indices for covariate categories as well.
#' @param base (integer or numeric) One or more years, used as as reference for the index.
#' If just a single year is given, the time total of the corresponding year will be uses as a reference.
#' If a range of years is given, the average of the corresponding time totals will be used as reference.\cr
#' Alternatively, the reference year(s) can be identified using their rank number,
#' i.e. \code{base=1} always refers to the starting year, \code{base=2} to the second year, etc.
#' @param level (numeric) the confidence interval required.
#' Must be in the range 0 to 1. A value of 0.95 results in 95\% confidence intervals.
#' The default value of NULL results in no confidence interval to be computed.
#' @param method (character) Method selector.
#' Options are \code{"formal"} (default) to use a formal computation of standard errors,
#' resulting in \eqn{\text{SE}=0} for the reference year,
#' and \code{"scaled"} to use a simpler approach, based on linear scaling of the time-totals SE.
#' @param long (logical) Switch to return 'long' output
#' (default is 'wide', as in rtrim versions < 2.2)
#'
#' @return A data frame containing indices and their uncertainty expressed as
#'   standard error. Depending on the chosen output, columns \code{fitted}
#'   and \code{se_fit}, and/or \code{imputed} and \code{se_imp} are present.
#'   If \code{covars} is \code{TRUE}, additional indices are computed for the
#'   individual covariate categories. In this case additional columns
#'   \code{covariate} and \code{category} are present. The overall indices are
#'   marked as covariate `Overall' and category 0.\cr
#'   In case \code{long=TRUE} a long table is returned, and a different naming convention is used.
#'   e.g., imputed/fitted info is in column \code{series},
#'   and standard error are always in column \code{SE}.

#' @export
#'
#' @family analyses
#'
#' @examples
#'
#' data(skylark)
#' z <- trim(count ~ site + time, data=skylark, model=2)
#' index(z)
#' # mimic classic TRIM:
#' index(z, "both")
#' # Extract standard errors for the imputed data
#' SE <- index(z,"imputed")$se_mod
#' # Include covariates
#' skylark$Habitat <- factor(skylark$Habitat) # hack
#' z <- trim(count ~ site + time + Habitat, data=skylark, model=2)
#' ind <- index(z, covars=TRUE)
#' plot(ind)
#' # Use alternative base year
#' index(z, base=3)
#' # Use average of first 5 years as reference for indexing
#' index(z, base=1:5)
#' # Prevent SE=0 for the reference year
#' index(z, method="scaled")
index <- function(x, which=c("imputed","fitted","both"), covars=FALSE, base=1, level=NULL, method=c("formal","scaled"), long=FALSE) {
  stopifnot(inherits(x,"trim"))

  # Match base to actual time points
  if (base[1] %in% x$time.id) {
    # First base time point is an actual year. Check that all the others are
    stopifnot(all(base %in% x$time.id))
    # Then convert year to time point
    for (i in seq_along(base)) base[i] <- which(base[i] == x$time.id)
  } else if (base<1 || base>x$nyear) {
    msg <- sprintf("Invalid base year %d. Must be either %d--%d or %d--%d",
                   base, 1, x$nyear, x$time.id[1], x$time.id[x$nyear])
    stop(msg)
  }

  # Start with overall indices (i.e. ignoring covariate categories, if applicable)
  # Computation and output is user-configurable
  which  <- match.arg(which)
  method <- match.arg(method)
  formal <- method=="formal"

  if (which=="fitted") {                # -- FITTED --
    # Call workhorse function to do the actual computation
    mod <- .index(x$tt_mod, x$var_tt_mod, base, level, x$sig2, formal)
    # Store results in a data frame
    if (long==TRUE) {                   # long format
      out <- data.frame(
        variable = "index",
        series   = "fitted",
        year     = x$time.id,
        value    = mod$tau,
        SE       = sqrt(mod$var_tau)
      )
    } else {                            # wide format
      out <- data.frame(
        time  = x$time.id,
        fitted = mod$tau,
        se_fit = sqrt(mod$var_tau))
    }
    # Add confidence interval?
    if (!is.null(level)) {
      out$lo <- mod$lo
      out$hi <- mod$hi
    }
  } else if (which=="imputed") {        # -- IMPUTED --
    # Idem, using the imputed time totals instead
    imp <- .index(x$tt_imp, x$var_tt_imp, base, level, x$sig2, formal)
    if (long==TRUE) {
      out <- data.frame(
        variable = "index",
        series   = "imputed",
        year     = x$time.id,
        value    = imp$tau,
        SE       = sqrt(imp$var_tau)
      )
    } else {
      out = data.frame(time    = x$time.id,
                       imputed = imp$tau,
                       se_imp  = sqrt(imp$var_tau))
    }
    if (!is.null(level)) {
      out$lo <- imp$lo
      out$hi <- imp$hi
    }
  }
  else if (which=="both") {           # -- BOTH --
    # Idem, using both modelled and imputed time totals.
    mod <- .index(x$tt_mod, x$var_tt_mod, base, level, x$sig2, formal)
    imp <- .index(x$tt_imp, x$var_tt_imp, base, level, x$sig2, formal)
    if (long==TRUE) {
      out1 <- data.frame(
        variable = "index",
        series   = "fitted",
        year     = x$time.id,
        value    = mod$tau,
        SE       = sqrt(mod$var_tau)
      )
      out2 <- data.frame(
        variable = "index",
        series   = "imputed",
        year     = x$time.id,
        value    = imp$tau,
        SE       = sqrt(imp$var_tau)
      )
      if (!is.null(level)) {
        out1$lo <- mod$lo
        out1$hi <- mod$hi
        out2$lo <- imp$lo
        out2$hi <- imp$hi
      }
      out <- rbind(out1, out2)
    } else { # Wide format
      if (!is.null(level)) stop("Confidence intervals can only be output for either imputed or fitted indices, but not for both")
      out = data.frame(time    = x$time.id,
                       fitted   = mod$tau,
                       se_fit  = sqrt(mod$var_tau),
                       imputed = imp$tau,
                       se_imp  = sqrt(imp$var_tau))
    }
  } else stop("Can't happen") # because other cases are catched by match.arg()

  # Add indices for covariate categories, if applicable
  if (covars) {
    # create first, "overall" section by pasting a "left" covar part with the "right" index data block
    out <- cbind(data.frame(covariate="Overall", category="(none)", stringsAsFactors=TRUE), out)

    tt <- x$covar_tt
    index <- list()
    ncovar <- length(tt)
    for (i in seq_len(ncovar)) {
      tti = tt[[i]]
      nclass <- length(tti)
      name <- names(tt)[i] # covariate name
      index[[name]] <- vector("list", nclass)
      for (j in seq_len(nclass)) {
        ttij <- tti[[j]]
        catname <- as.character(ttij$class) # old code; todo: fix properly by setting a cat name earlier
        catname <- levels(x$covars[[name]])[j]
        # construct left block (ensure proper # of rows)
        n <- length(x$time.id)
        df <- data.frame(covariate = rep(ttij$covariate, n),
                         category  = rep(catname, n),
                         stringsAsFactors=TRUE)
        if (!long) df$time <- x$time.id
        # Compute model index+variance
        if (which %in% c("fitted","both")) {
          idx <- .index(ttij$mod, ttij$var_mod, base, level, x$sig2, formal)
          if (long) {
            df2 <- data.frame(variable="index",series="fitted", year=x$time.id, value=idx$tau, SE=sqrt(idx$var_tau))
          } else{
            df2 <- data.frame(fitted=idx$tau, se_fit=sqrt(idx$var_tau))
          }
          if (!is.null(level)) {
            df2$lo <- idx$lo
            df2$hi <- idx$hi
          }
          if (long) {
            out <- rbind(out, cbind(df, df2))
          } else {
            df <- cbind(df, df2)
          }
        }
        # Idem for imputed index + variance
        if (which %in% c("imputed","both")) {
          idx <- .index(ttij$imp, ttij$var_imp, base, level, x$sig2, formal)
          if (long) {
            df2 <- data.frame(variable="index",series="imputed", year=x$time.id, value=idx$tau, SE=sqrt(idx$var_tau))
          } else {
            df2 <- data.frame(imputed=idx$tau, se_imp=sqrt(idx$var_tau))
          }
          if (!is.null(level)) {
            df2$lo <- idx$lo
            df2$hi <- idx$hi
          }
          if (long) {
            out <- rbind(out, cbind(df, df2))
          } else {
            df <- cbind(df, df2)
          }
        }
        if (!long) out <- rbind(out, df)
      }
    }
  }
  class(out) <- c("trim.index", class(out)) # prepend new class type
  out
}


#' Plot time-indices from trim output.
#'
#' Uncertainty ranges exressed as standard errors are always plotted.
#' Confidence intervals are plotted when they are present in the \code{trim.index} object, i.e. when requested for in the call to \code{\link{index}}.
#'
#' @param x an object of class \code{trim.index}, as resulting from e.g. a call to \code{\link{index}}.
#' @param ... additional \code{trim.index} objects, or parameters that will be passed on to \code{\link[graphics]{plot}}.
#' @param names   optional character vector with names for the various series.
#' @param covar \code{[character]} the name of a covariate to include in the plot.
#'   If set to \code{"auto"} (the default), the first (or only) covariate is used.
#'   If set to \code{"none"} plotting of covariates is suppressed and only the overall index is shown.
#' @param xlab a title for the x-axis. The default value is "auto" will be changed to "Time Point" if the time ID's start at 1, and to "Year" otherwise.
#' @param ylab a title for the y-axis. The default value is "Index".
#' @param pct  Switch to show the index values as percent instead as fraction (i.e., for the base year it will be 100 instead of 1)
#' @param band Defines if the uncertainty band will be plotted using standard errors ("se") or confidence intervals ("ci").
#'
#' @export
#'
#' @family analyses
#' @family graphical post-processing
#'
#' @examples
#'
#' # Simple example
#' data(skylark2)
#' z <- trim(count ~ site + year, data=skylark2, model=3)
#' idx <- index(z)
#' plot(idx)
#'
#' # Example with user-modified title, and different y-axis scaling
#' plot(idx, main="Skylark", pct=TRUE)
#'
#' # Using covariates:
#' z <- trim(count ~ site + year + habitat, data=skylark2, model=3)
#' idx <- index(z, covars=TRUE)
#' plot(idx)
#'
#' # Suppressing the plotting of covariate indices:
#' plot(idx, covar="none")
#'
plot.trim.index <- function(x, ..., names=NULL, covar="auto", xlab="auto", ylab="Index", pct=FALSE, band="se") {

  # This function operates in 3 modes:
  # 1 : single index
  # 2 : multiple indices
  # 3 : single index, with covariate categories

  # Assume a simple job
  mode <- 1
  z <- x
  zz <- list(z)
  nidx <- 1

  # First determine the mode. Check for mode 2
  # ellipsis <- list(...)
  ellipsis <- as.list(substitute(list(...)))[-1L]
  n <- length(ellipsis) # number of ellipsis arguments
  if (n>0) {
    keep <- rep(TRUE, n) # records which (named!) arguments to keep for passing to plot()
    if (is.null(names(ellipsis))) { # none has names
      named <- rep(FALSE, n)
    } else {                        # some have names
      named <- nchar(names(ellipsis)) > 0
    }
    for (i in seq_along(ellipsis)) {
      if (named[i]) next # skip over named arguments that are captured in the ...
      item <- ellipsis[[i]]
      if (is.symbol(item)) item <- eval(item) # needed to convert symbol -> data.frame
      if (inherits(item, "trim.index")) {
        mode <- 2
        nidx <- nidx + 1
        zz[[nidx]] <- item
        keep[i] <- FALSE # additional index data sets are removed from the ellipsis argument
      } else if (inherits(item, "character")) {
        # todo: check if this arguments immediately follows an index argument
        attr(zz[[nidx]], "tag") <- item
        keep[i] <- FALSE
      } else stop("Unknown type for unnamed argument")
    }
    ellipsis <- ellipsis[keep]
  }

  # Now we know how much independent indices are involved.
  nidx <- length(zz)


  # Check for mode 3 (not allowed in cojunction with mode 2)
  if (covar=="auto" && !"covariate" %in% names(z)) covar <- "none"
  use.covars <- covar!="none"
  if (use.covars) {
    if (mode==2) stop("Plotting of indices with covariates is not allowed when multiple indices.")
    if (!"covariate" %in% names(z))
      stop("No covariate info in index data")
    if (covar=="auto") covar <- levels(z$covariate)[2] # Skip "Overall"
    if (!covar %in% levels(z$covariate))
      stop(sprintf('Covariate "%s" not present in index data', covar))
    mode <- 3
  }

  # Create custom palette based on Color Brewer Set 1
  brewer_set1 <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#999999")
  opaque <- brewer_set1
  aqua   <- brewer_set1
  if (mode==1 || mode==2) {
    # one or multiple indices: different colors, but same transparency
    for (i in seq_along(aqua)) aqua[i] <- adjustcolor(aqua[i], 0.5)
  } else if (mode==3) {
    # covariates: main index is darker than the covariates
    aqua[1] <- adjustcolor(aqua[1], 0.5)
    for (i in 2:9) aqua[i] <- adjustcolor(aqua[i], 0.25) # covariates plot lighter
  }

  # get index/stderr columns. Prefer the imputed ones
  idx_col <- integer(nidx)
  for (i in seq_len(nidx)) {
    zi <- zz[[i]]
    if ("imputed" %in% names(zi)) {
      idx_col[i] <- which(names(zi)=="imputed")
    } else {
      idx_col[i] <- which(names(zi)=="fitted")
    }
  }
  err_col <- idx_col +1

  # Handle the "pct" switch
  yscale <- ifelse(pct, 100L, 1L)

  # set up index series
  series <- list()
  nseries <- 0

  if (mode==1) {
    z <- zz[[1]] # get from list (it might have been named...)
    x <- z$time
    y <- z[[idx_col]] * yscale
    err <- z[[err_col]] * yscale
    yslo = y - err # 's' means std err
    yshi = y + err
    nseries <- 1
    yclo <- z$lo * yscale # 'c' means confidence interval, might be NULL which is OK
    ychi <- z$hi * yscale
    if (band=="ci") {
      if (is.null(yclo) || is.null(ychi)) stop("No confidence interval present")
      yslo <- yclo
      yshi <- ychi
      yclo <- ychi <- NULL
    }
    name <- attr(z, "tag") # might be NULL
    series[[1]] <- list(x=x, y=y,
                        yslo=yslo, yshi=yshi, yclo=yclo, ychi=ychi,
                        fill=aqua[1], stroke=opaque[1], lty="solid", name=name)
  } else if (mode==2) {
    # Create or handle names
    if (is.null(names)) names <- sprintf("<no name> #%d", 1:nidx) # default names
    if (length(names)!=nidx) stop("Number of names is not equal to number of series")
    for (i in seq.int(nidx)) {
      zi <- zz[[i]]
      name <- attr(zi, "tag")
      if (!is.null(name)) names[i] <- name
    }

    # Now create the series
    for (i in seq.int(nidx)) {
      zi <- zz[[i]]
      x <- zi$time
      y <- zi[[idx_col[i]]] * yscale
      err <- zi[[err_col[i]]] * yscale
      yslo = y - err
      yshi = y + err
      yclo <- zi$lo * yscale # 'c' means confidence interval, might be NULL which is OK
      ychi <- zi$hi * yscale
      if (band=="ci") {
        if (is.null(yclo) || is.null(ychi)) stop("No confidence interval present")
        yslo <- yclo
        yshi <- ychi
        yclo <- ychi <- NULL
      }
      nseries <- nseries + 1
      series[[i]] <- list(x=x, y=y,
                          yslo=yslo,yshi=yshi, yclo=yclo, ychi=ychi,
                          fill=aqua[i], stroke=opaque[i], lty="solid", name=names[i])
    }
  } else if (mode==3) {
    # Setup series for covariate categories

    # First split in overall and covar-cats
    z <- zz[[1]]
    rows    <- z$covariate=="Overall"
    overall <- z[rows, ]
    rows    <- z$covariate==covar
    other   <- z[rows, ]

    # First setup overall index...
    x <- overall$time
    y <- overall[[idx_col]] * yscale
    err <- overall[[err_col]] * yscale
    yslo = y - err
    yshi = y + err
    yclo <- overall$lo * yscale # 'c' means confidence interval, might be NULL which is OK
    ychi <- overall$hi * yscale
    if (band=="ci") {
      if (is.null(yclo) || is.null(ychi)) stop("No confidence interval present")
      yslo <- yclo
      yshi <- ychi
      yclo <- ychi <- NULL
    }
    nseries <- 1
    name <- attr(overall, "tag")
    if (is.null(name)) name <- "Overall"
    series[[nseries]] <- list(x=x,y=y,
                              yslo=yslo, yshi=yshi, yclo=yclo, ychi=ychi,
                              fill=aqua[1], stroke=opaque[1], lty="solid", name=name)

    # ... then the covariate categories
    cats <- levels(factor(other$category)) # Determine the *actual* categories
    for (cat in cats) {
      nseries <- nseries+1
      rows <- other$category==cat
      x <- overall$time
      y <- other[rows, idx_col] * yscale
      err <- other[rows, err_col] * yscale
      yslo = y - err
      yshi = y + err
      yclo <- other$lo[rows] * yscale # 'c' means confidence interval, might be NULL which is OK
      ychi <- other$hi[rows] * yscale
      if (band=="ci") {
        if (is.null(yclo) || is.null(ychi)) stop("No confidence interval present")
        yslo <- yclo
        yshi <- ychi
        yclo <- ychi <- NULL
      }
      name <- sprintf("%s: %s", covar, cat)
      series[[nseries]] <- list(x=x, y=y,
                                yslo=yslo, yshi=yshi, yclo=yclo, ychi=ychi,
                                fill=aqua[nseries], stroke=opaque[nseries], lty="dashed", name=name)
    }
  } else stop("Can't happen.")

  # Determine axis labels iff automatic (just use the last 'x' defined)
  if (xlab=="auto") {
    xlab <- ifelse(x[1]==1, "Time Point", "Year")
  }

  ### analysis ###

  xrange <- range(x)

  # Compute the y-range of all series: std.err ranges
  yrange <- range(series[[1]]$yslo, series[[1]]$yshi)
  if (nseries>1) for (i in 2:nseries) {
    yrange <- range(series[[i]]$yslo, series[[i]]$yshi, yrange)
  }
  # ... CI ranges
  for (i in 1:nseries) {
    yrange <- range(series[[i]]$yclo, series[[i]]$ychi, yrange)
  }
  # ... some extras
  yrange <- range(yrange, 0.0) # honest scaling
  yrange <- range(yrange, 1.1*yscale) # include index=1 (or 100%)


  ### plotting ###

  # Plot 'empty' overall index.
  # We do need some special tricks to pass the plot-specific elements of ...
  par(las=1)
  args <- ellipsis
  args$x <- NULL
  args$y <- NULL
  args$type='n'
  args <- c(list(x=NULL,y=NULL, type='n', xlim=xrange, ylim=yrange, xlab=xlab, ylab=ylab), ellipsis)
  # plot(NULL, NULL, type='n', xlim=xrange, ylim=yrange, xlab=xlab, ylab=ylab, ...) # won't work
  do.call(plot, args) # does work
  abline(h=yscale, lty="dashed")

  # set up legend
  leg.names <- "Overall"
  leg.colors <- opaque[1]

  # plot all series, in reverse order to highlight the first series.
  # First the error bars...
  for (i in rev(1:nseries)) {
    ser <- series[[i]]
    xx <- c(ser$x, rev(ser$x))
    yy <- c(ser$yslo, rev(ser$yshi))
    polygon(xx, yy, col=ser$fill, border=NA)
    segments(ser$x, ser$yslo, ser$x, ser$yshi, col="white", lwd=2)
  }
  # Optionally confidence intervals
  for (i in rev(1:nseries)) {
    ser <- series[[i]]
    if (length(ser$yclo)==0) next # No CI: skip series
    lines(ser$x, ser$yclo, col=ser$stroke, lwd=1, lty="dashed")
    lines(ser$x, ser$ychi, col=ser$stroke, lwd=1, lty="dashed")
  }
  # ... and then the indices themselves
  for (i in rev(1:nseries)) {
    ser <- series[[i]]
    lines(ser$x, ser$y, col=ser$stroke, lwd=2, lty=ser$lty)
    points(ser$x, ser$y, col=ser$stroke, pch=16)
  }

  # Plot legend if appropriate
  if (nseries>1) {
    leg.names <- leg.colors <- character(nseries)
    leg.lty <- character(nseries)
    for (i in 1:nseries) {
      leg.names[i]  <- series[[i]]$name
      leg.colors[i] <- series[[i]]$stroke
      leg.lty[i]    <- series[[i]]$lty
    }

    # compute optimal legend position
    xx <- numeric(0)
    yy <- numeric(0)
    for (i in 1:nseries) {
      ser <- series[[i]]
      xx <- c(xx, ser$x)
      yy <- c(yy, ser$y)
    }
    mx <- mean(x)
    my <- mean(y)
    if (mx < mean(xrange)) {
      if (my < mean(yrange)) leg.pos <- "topright"
      else                   leg.pos <- "bottomright"
    } else {
      if (my < mean(yrange)) leg.pos <- "topleft"
      else                   leg.pos <- "bottomleft"
    }

    legend(leg.pos, legend=leg.names, col=leg.colors, lty=leg.lty, lwd=2, bty='n', inset=0.02, y.intersp=1.5);
  }
}


