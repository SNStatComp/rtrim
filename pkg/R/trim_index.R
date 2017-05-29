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


.index <- function(tt, var_tt, base) {

  nbase <- length(base)

  if (nbase==1) {
    tau <- tt / tt[base]
  } else {
    tau <- tt / mean(tt[base])
  }

  J <- length(tt)
  var_tau <- numeric(J)
  d <- matrix(0, nbase+1, 1)
  for (j in 1:J) {

    if (nbase==1 && j==base) {
      var_tau[j] <- 0.0
    } else {
      d[1:nbase] <- -nbase * tt[j] / sum(tt[base])^2
      d[nbase+1] <- nbase / sum(tt[base])
    }

    idx <- c(base, j)
    V <- var_tt[idx,idx]
    var_tau[j] <- t(d) %*% V %*% d
  }

  out <- list(tau=tau, var_tau=var_tau)
}

# ========================================================== User interface ====

#' Extract time-indices from trim output. Infices are obtained by dividing the modelled or imputed time totals by a reference value.
#' Most commonly, the time totals for the first time point are used as reference.
#' As a result, the index value for this first time point will be 1.0, with a standard error of 0.0 by definition.
#'
#' @param x an object of class \code{\link{trim}}
#' @param which \code{[character]} Selector to distinguish between time indices based on the imputed data (default),
#' the fitted model, or both.
#' @param covars \code{[logical]} Switch to compute indices for covariate categories as well.
#' @param base \code{[integer|numeric]} One or more base time point, used as as reference for the index.
#' If just a single number is given, the time total of the correspondong time point will be uses as  reference.
#' If a range of numbers is given, the average of the corresponding time totals will be used as reference.
#' The base time points can be given in the interval 1...J, or,
#' if the time points are proper years, say year1...yearn, the base year can be given.
#' So, if the data range 2000...2016, \code{base=2} and \code{base=2001} are equivalent.
#'
#' @return A data frame containing indices and their uncertainty expressed as
#'   standard error. Depending on the chosen output, columns \code{fitted}
#'   and \code{se_fit}, and/or \code{imputed} and \code{se_imp} are present.
#'   If \code{covars} is \code{TRUE}, additional indices are computed for the
#'   individual covariate categories. In this case additional columns
#'   \code{covariate} and \code{category} are present. The overall indices are
#'   marked as covariate `Overall' and category 0.
#'
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
#' z <- trim(count ~ site + time + Habitat, data=skylark, model=2)
#' ind <- index(z, covars=TRUE)
#' plot(ind)
#' # Use alternative base year
#' index(z, base=3)
#' # Use average of first 5 years as reference for indexing
#' index(z, base=1:5)
index <- function(x, which=c("imputed","fitted","both"), covars=FALSE, base=1) {
  stopifnot(inherits(x,"trim"))

  # Match base to actual time points
  if (base[1] %in% x$time.id) {
    # First base time point is an actual year. Check that all the others are
    stopifnot(all(base %in% x$time.id))
    # Then convert year to time point
    for (i in seq_along(base)) base[i] <- which(base[i] == x$time.id)
  }

  # Start with overall indices (i.e. ignoring covariate categories, if applicable)
  # Computation and output is user-configurable
  which <- match.arg(which)
  if (which=="fitted") {
    # Call workhorse function to do the actual computation
    mod <- .index(x$tt_mod, x$var_tt_mod, base)
    # Store results in a data frame
    out <- data.frame(time  = x$time.id,
                      fitted = mod$tau,
                      se_fit = sqrt(mod$var_tau))
  } else if (which=="imputed") {
    # Idem, using the imputed time totals instead
    imp <- .index(x$tt_imp, x$var_tt_imp, base)
    out = data.frame(time    = x$time.id,
                     imputed = imp$tau,
                     se_imp  = sqrt(imp$var_tau))
  } else if (which=="both") {
    # Idem, using both modelled and imputed time totals.
    mod <- .index(x$tt_mod, x$var_tt_mod, base)
    imp <- .index(x$tt_imp, x$var_tt_imp, base)
    out = data.frame(time    = x$time.id,
                     fitted   = mod$tau,
                     se_fit  = sqrt(mod$var_tau),
                     imputed = imp$tau,
                     se_imp  = sqrt(imp$var_tau))
  } else stop("Can't happen") # because other cases are catched by match.arg()

  # Add indices for covariate categories, if applicable
  if (covars) {
    out <- cbind(data.frame(covariate="Overall", category=0), out)

    tt <- x$covar_tt
    index <- list()
    ncovar <- length(tt)
    for (i in seq_len(ncovar)) {
      tti = tt[[i]]
      nclass <- length(tti)
      name <- names(tt)[i]
      index[[name]] <- vector("list", nclass)
      for (j in seq_len(nclass)) {
        ttij <- tti[[j]]
        df = data.frame(covariate=ttij$covariate, category=ttij$class, time=x$time.id)
        # Compute model index+variance
        if (which %in% c("fitted","both")) {
          idx <- .index(ttij$mod, ttij$var_mod, base)
          df2 <- data.frame(fitted=idx$tau, se_fit=sqrt(idx$var_tau))
          df <- cbind(df, df2)
        }
        # Idem for imputed index + variance
        if (which %in% c("imputed","both")) {
          idx <- .index(ttij$imp, ttij$var_imp, base)
          df2 <- data.frame(imputed=idx$tau, se_imp=sqrt(idx$var_tau))
          df <- cbind(df, df2)
        }
        out <- rbind(out, df)
      }
    }
  }
  class(out) <- c("trim.index","data.frame")
  out
}


#' Plot time-indices from trim output
#'
#' @param x an object of class \code{trim.index}, as resulting from e.g. a call to \code{\link{index}}.
#' @param covar \code{[character]} the name of a covariate to include in the plot.
#'   If set to \code{"auto"} (the default), the first (or only) covariate is used.
#'   If set to \code{"none"} plotting of covariates is suppressed and only the overall index is shown.
#' @param ... Further options passed to \code{\link[graphics]{plot}}
#'
#' @export
#'
#' @family analyses
#'
#' @examples
#'
#' data(skylark)
#' z <- trim(count ~ site + time + Habitat, data=skylark, model=2)
#' idx <- index(z, covars=TRUE)
#' plot(idx, covar="Habitat", main="Skylark")
#'
plot.trim.index <- function(x, covar="auto", ...) {
  z <- x # hack
  # Create custom palette based on Color Brewer Set 1
  brewer_set1 <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#999999")
  opaque <- brewer_set1
  aqua   <- brewer_set1
  aqua[1] <- adjustcolor(aqua[1], 0.3)
  for (i in 2:9) aqua[i] <- adjustcolor(aqua[i], 0.1)

  # Use covariates in *this* plot?
  if (covar=="auto" && !"covariate" %in% names(z)) covar <- "none"
  use.covars <- covar!="none"
  if (use.covars) {
    if (!"covariate" %in% names(z))
      stop("No covariate info in index data")
    if (covar=="auto") covar <- levels(z$covariate)[2] # Skip "Overall"
    if (!covar %in% levels(z$covariate))
      stop(sprintf('Covariate "%s" not present in index data', covar))
  }

  # get index/stderr columns. Prefer the imputed ones
  if ("imputed" %in% names(z)) {
    idx_col <- which(names(z)=="imputed")
  } else {
    idx_col <- which(names(z)=="fitted")
  }
  err_col <- idx_col +1

  if ("covariate" %in% names(z)) {
    # split in overall and covar-cats
    rows <- z$covariate=="Overall"
    overall <- z[rows, ]
    rows <- z$covariate==covar
    other   <- z[rows, ]
  } else {
    overall <- z
  }

  x = overall$time
  y1 = overall[[idx_col]]
  ylo1 = y1 - overall[[err_col]]
  yhi1 = y1 + overall[[err_col]]
  yrange1 = range(y1, ylo1, yhi1)


  if (use.covars) {
    other$category <- factor(other$category)
    yrange2 <- range(other[[idx_col]] - other[[err_col]])
    yrange3 <- range(other[[idx_col]] + other[[err_col]])
    yrange <- range(yrange1, yrange2, yrange3)
  } else {
    yrange <- yrange1
  }

  # Plot 'empty' overall index
  par(las=1)
  plot(x, y1, type='n', ylim=yrange, xlab="Time point", ylab="Index", ...)
  abline(h=1.0, lty="dashed")
  xx = c(x, rev(x))

  # set up legend
  leg.names <- "Overall"
  leg.colors <- opaque[1]

  if (use.covars) {
    # Plot covar cat indices
    cidx <- 2 # color index to use (first element is for overall index)
    for (cat in levels(other$category)) {
      rows <- other$category==cat
      y = other[rows, idx_col]
      ylo = y - other[rows, err_col]
      yhi = y + other[rows, err_col]
      yy = c(ylo, rev(yhi))
      polygon(xx,yy,col=aqua[cidx], border=NA)
      lines(x, y, col=opaque[cidx], lwd=2)
      segments(x,ylo, x,yhi, col="white", lwd=2)
      points(x,y, col=opaque[cidx], pch=16)
      # Append legend info
      leg.names <- c(leg.names, sprintf("%s cat. %s", covar, cat))
      leg.colors  <- c(leg.colors, opaque[cidx])
      # Move on
      cidx <- cidx+1
    }
  }

  # Finally plot overall index
  yy = c(ylo1, rev(yhi1))
  polygon(xx,yy, col=aqua[1], border=NA)
  lines(x, y1, col=opaque[1], lwd=3)
  segments(x,ylo1, x,yhi1, col="white", lwd=1)
  points(x,y1, col=opaque[1], pch=16, cex=1)

  if (use.covars)
    legend("topleft", legend=leg.names, col=leg.colors, lty=1, lwd=2, bty='n', inset=0.02, y.intersp=1.5);
}

plot.midx <- function(idx1, ..., main="", leg.pos="topleft")
{
  # Create custom palette based on Color Brewer Set 1
  brewer_set1 <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#999999")
  opaque <- brewer_set1
  aqua   <- brewer_set1
  for (i in 1:9) aqua[i] <- adjustcolor(aqua[i], 0.3)

  # Build a list of indices with optional titles
  idx = list(idx1)
  optional = list(...)

  nopt = length(optional)
  for (i in seq_len(nopt)) {
    x = optional[[i]]
    if ("character" %in% class(x)) {
      attr(idx[[length(idx)]], "tag") <- x
    } else if ("trim.index" %in% class(x)) {
      idx[[length(idx)+1]] <- x
    } else {
      stop(sprintf("Invalid data type for optional argument %d: %s", i, class(x)))
    }
  }

  # First pass to compute total range
  n = length(idx)
  for (i in 1:n) {
    x = idx[[i]][[1]] # Time point or years
    y = idx[[i]][[2]] # imputed or fitted index
    s = idx[[i]][[3]] # Standard error
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

  # empty plot for correct axes
  plot(xrange, yrange, type='n', xlab="Time point", ylab="Index", main=main)

  # Second pass: SE panels
  for (i in 1:n) {
    x = idx[[i]][[1]] # Time point or years
    y = idx[[i]][[2]] # imputed or fitted index
    s = idx[[i]][[3]] # Standard error
    ylo = y-s
    yhi = y+s

    xx = c(x, rev(x))
    ci = c(ylo, rev(yhi))

    polygon(xx,ci, col=aqua[i], border=NA)
  }

  # Third pass: SE lines
  for (i in 1:n) {
    x = idx[[i]][[1]] # Time point or years
    y = idx[[i]][[2]] # imputed or fitted index
    s = idx[[i]][[3]] # Standard error
    ylo = y-s
    yhi = y+s

    segments(x,ylo, x,yhi, col="white", lwd=1)
  }

  # Fourth pass: lines+points
  for (i in 1:n) {
    x = idx[[i]][[1]] # Time point or years
    y = idx[[i]][[2]] # imputed or fitted index

    lines(x,y, col=opaque[i], lwd=2)
    points(x,y, col=opaque[i], pch=16)
  }

  # Fifth pass: legend
  nnamed  = 0
  nnoname = 0
  for (i in 1:n) {
    s <- attr(idx[[i]],"tag")
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