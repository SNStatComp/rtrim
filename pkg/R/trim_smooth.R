run_kalman <- function(x,y,Q11, h=1, smooth=FALSE) {
  # setup
  a00 <- matrix(0, 2,1)
  P00 <- 10000*max(y) * diag(2)
  Q <- matrix(c(Q11,0,0,0), 2,2)
  T <- matrix(c(2,1,-1,0), 2,2) # Transition matrix
  z <- matrix(c(1,0), 2,1) # observation coefficients


  N <- length(y)
  nu <- numeric(N)
  f <- numeric(N)

  if (length(h)==1) h <- rep(h, N)

  ap <- a <- array(0, c(2,1,N)) # State vector (history)
  Pp <- P <- array(0, c(2,2,N)) # Variance-covariance matrix (history)

  for (i in 1:N) {
    # 1. Predict
    aprev <- if (i==1) a00 else a[,,i-1]
    Pprev <- if (i==1) P00 else P[,,i-1]
    ap[,,i] <- T %*% aprev                                                  # D.7
    Pp[,,i] <- T %*% Pprev %*% t(T) + Q                                     # D.8
    # 2. Innovations
    nu[i] <- y[i] - t(z) %*% ap[,,i]
    f[i] <- t(z) %*% Pp[,,i] %*% z + h[i]
    # 3. Update
    a[,,i] <- ap[,,i] + Pp[,,i] %*% z %*% nu[i] / f[i]                       # D.12
    P[,,i] <- Pp[,,i] - Pp[,,i] %*% z %*% t(z) %*% Pp[,,i] / f[i] # D.13
  }

  # Compute likelihood
  Ns <- ceiling(0.15*N) # Skip first time steps
  idx <- (Ns+1) : N
  sigma2 <- mean(nu[idx]^2/f[idx])                                         # D.24
  loglik <- sum(log(sigma2 * f[idx]))                                     # D.25

  out <- list(
    rawdata   = list(x=x, y=y, se=sqrt(h)),
    yy_filter = a[1,1,],
    sigma2=sigma2, loglik=loglik,
    internal  = list(Q=Q, a=a, P=P, Pp=Pp))

  # Optional smoothing
  if (smooth) {
    as    <- array(0, c(2,1,N)) # Smoothed state vector history
    Ps    <- array(0, c(2,2,N)) # Smoothed variance matrix history
    Pstar <- array(0, c(2,2,N))

    as[,,N] <- a[,,N]
    Ps[,,N] <- P[,,N]
    for (i in (N-1) : 1) {
      Pstar[,,i] <- P[,,i] %*% t(T) %*% solve(Pp[,,i+1])                           # D.16
      as[,,i] <- a[,,i] + (Pstar[,,i] %*% (as[,,i+1] - T %*% a[,,i]))              # D.14
      Ps[,,i] <- P[,,i] + Pstar[,,i] %*% (Ps[,,i+1] - Pp[,,i+1]) %*% t(Pstar[,,i]) # D.15
      # Note: the last Pp above is not required (could be just P)
    }
    # out$as <- as
    out$smoothed <- list(x=x, y=as[1,1,], se=sqrt(Ps[1,1, ]*sigma2))
    out$internal$Ps    <- Ps    # required for lagdiff
    out$internal$Pstar <- Pstar # idem
  }
  out
}

fit_kalman <- function(x,y,h, optimize=TRUE) {
  f <- function(Q, x,y,h) {
    out <- run_kalman(x,y, Q,h, smooth=FALSE)
    out$loglik
  }
  if (optimize) {
    Qhi  <- 1
    for (i in 1:100) {
      print(i)
      Qopt <- optimize(f, c(0,Qhi), x=x,y=y,h=h)$minimum
      if (Qopt < 0.8*Qhi) {
        cat(sprintf("*** Qopt = %.3f %.3f ***\n", Qopt, Qhi))
        break
      }
      Qhi <- Qhi * 2
    }
  } else Qopt <- 1
  out <- run_kalman(x,y, Qopt, h, smooth=TRUE)
  out
}

smooth <- function(x, which=c("imputed","fitted")) {
  stopifnot(class(x)=="trim")
  which = match.arg(which)

  # extract vars from TRIM output
  years <- x$time.id
  tt_mod <- x$tt_mod
  tt_imp <- x$tt_imp
  var_tt_mod <- x$var_tt_mod
  var_tt_imp <- x$var_tt_imp
  J <- ntime <- x$ntime

  if (which=="imputed") {
    tt     <- tt_imp
    var_tt <- var_tt_imp
    src = "imputed"
  } else if (which=="fitted") {
    tt     <- tt_mod
    var_tt <- var_tt_mod
    src = "fitted"
  }

  J = length(tt)
  out <- fit_kalman(years, tt, diag(var_tt), optimize=TRUE)
  #out <- fit_kalman(years, tt, 1, optimize=FALSE)
  out$type="normal"
  out$tt = tt
  out$err = sqrt(diag(var_tt))
  out$internal <- NULL # remove trendspotter internal info
  structure(out, class="trim.smooth")
}

#--------------------------------------------------------------------- Plot ----

#' Plot overall slope
#'
#' Creates a plot of the overall slope, its 95\% confidence band, the
#' total population per time and their 95\% confidence intervals.
#'
#' @param x An object of class \code{trim.overall} (returned by \code{\link{overall}})
#' @param imputed \code{[logical]} Toggle to show imputed counts
#' @param ... Further options passed to \code{\link[graphics]{plot}}
#'
#' @family analyses
#'
#' @examples
#' data(skylark)
#' m <- trim(count ~ site + time, data=skylark, model=2)
#' plot(overall(m))
#'
#' @export
plot.trim.smooth <- function(x, imputed=TRUE, ...) {
  trend <- x
  title <- if (is.null(list(...)$main)){
    attr(trend, "title")
  } else {
    list(...)$main
  }

  make_cband <- function(x) {
    ylo <- x$y - x$se * 1.96
    yhi <- x$y + x$se * 1.96
    xconf <- c(x$x, rev(x$x))
    yconf <- c(ylo, rev(yhi))
    x$cb <- cbind(xconf, yconf)
    x
  }

  # Convert SE to confidence band
  trend$rawdata  <- make_cband(trend$rawdata)
  trend$smoothed <- make_cband(trend$smoothed)

  xrange  <- range(trend$rawdata$x)
  yrange1 <- range(trend$rawdata$cb[,2])
  yrange2 <- range(trend$smoothed$cb[,2])
  yrange <- range(0.0, yrange1, yrange2)  # include y=0 !

  # Now plot layer-by-layer (using ColorBrewer colors)
  cbred <- rgb(228,26,28, maxColorValue = 255)
  cbblue <- rgb(55,126,184, maxColorValue = 255)
  plot(xrange, yrange, type='n', xlab="Year", ylab="Count", las=1, main=title,...)
  # raw data in red; trend in blue
  polygon(trend$rawdata$cb, col=adjustcolor(cbred, 0.2), lty=0)
  polygon(trend$smoothed$cb, col=adjustcolor(cbblue, 0.2), lty=0)

  points(trend$rawdata$x, trend$rawdata$y,  col=cbred, pch=16)
  lines(trend$smoothed$x, trend$smoothed$y, col=cbblue, lwd=2)
  trend
}
