.colorize <- function(m, pal1, pal2=NULL, idx2=NULL, log=FALSE, na.col) {
  # m : matrix
  # pal1 : palette
  # pal2 : optional second palette
  # idx2 : logical matrix of grid cells where pal2 should be used
  # log  : flag to set log-transformation of values in m
  # na.col : color to use for NA cells

  # Internal functions

  .normalize <- function(x) {
    minx <- min(x, na.rm=TRUE)
    maxx <- max(x, na.rm=TRUE)
    (x-minx) / (maxx-minx)
  }

  # (next lines are obsolete: 0 will be painted white)
  # # test if all data points are >0 iff log=TRUE
  # if (log) {
  #   nzero <- sum(m==0, na.rm=TRUE)
  #   if (nzero>0) {
  #     msg <- sprintf("Can't make a heat map because %d zero elements; consider setting log=FALSE.", nzero)
  #     warning(msg)
  #     # stop(msg, call.=FALSE)
  #   }
  # }

  rgb.pal.1 <- col2rgb(pal1) / 255.0
  variant <- ifelse(is.null(pal2), 1, 2) # 1: only pal1; 2: pal1+pal2
  if (variant==2) {
    rgb.pal.2 <- col2rgb(pal2) / 255.0
    stopifnot(!is.null(idx2))
  }
  if (class(na.col)=="character") na.col <- col2rgb(na.col)/255.0
  ncolor <- length(pal1)
  nr <- nrow(m)
  nc <- ncol(m)
  img <- array(0, c(nr, nc, 3L))

  # determine special values
  if (log) {
    na_idx   <- is.na(m)
    zero_idx <- is.finite(m) & m==0
    ok_idx  <- is.finite(m)  & m>0
  } else {
    na_idx   <- is.na(m)
    ok_idx   <- is.finite(m)
  }

  # transform
  if (log) {
    m[zero_idx] <- NA
    m <- log10(m)
  }

  idx1 <- is.finite(m)

  # try to uniformly distribute colors. Probably there's only one largest value.
  eps <- 1e-7
  im <- 1L + as.integer((ncolor-eps) * .normalize(m))
  for (i in 1:3) {
    layer <- matrix(na.col[i], nr, nc)
    if (variant==1) {
      layer[idx1] <- rgb.pal.1[i, im[ok_idx]]
    } else {
      layer[idx1] <- rgb.pal.1[i, im[ok_idx]]
      layer[idx2] <- rgb.pal.2[i, im[idx2]]
    }
    # mark zeros as WHITE
    if (log) {
      layer[zero_idx] <- 1.0 # white
    }
    img[ , ,i] <- layer
  }
  img
}


# heatmap <- function(x) UseMethod("heatmap")


#' Plot a heatmap representation of observed and/or imputed counts.
#'
#' This function organizes the observed and/or imputed counts into a matrix where
#' rows represent sites and columns represent time points.
#' A bitmap image is constructed in which each pixel corresponds to an element of this matrix.
#' Each pixel is colored according the correspondong count status, and the type of heatmap plot requested ('data', 'imputed' or 'fitted').
#'
#' The 'imputed' heatmap uses the most elaborate color scheme:
#' Site/time combinations that are observed are colored red, the higher the count, the darker the red.
#' Site/time combinations that are imputed are colored blue, the higher the estimate, the darker the blue.
#'
#' For the 'data' heatmap, missing site/time combinations are colored gray.
#'
#' For the 'fitted' heatmap, all site/time combinations are colored blue.
#'
#' By default, all counts are log-transformed prior to colorization, and observed counts of 0 are indicates as white pixels.
#'
#' @param z output of a call to \code{\link{trim}}.
#' @param what the type of heatmap to be plotted: 'data' (default), 'imputed' or 'fitted'.
#' @param log flag to indicate whether the count should be log-transformed first.
#' @param xlab x-axis label. The default value "auto" will evaluate to either "Year" or "Time point"
#' @param ylab y-axis label
#' @param ... other parameters to be passed to \code{\link[graphics]{plot}}
#'
#' @export
#' @family graphical post-processing
#'
#' @examples
#' data(skylark2)
#' z <- trim(count ~ site + year, data=skylark2, model=3)
#' heatmap(z,"imputed")
#'
heatmap <- function(z, what=c("data","imputed","fitted"), log=TRUE, xlab="auto", ylab="Site #", ...) {
  # Create an RGB heatmap image

  # Define the color palettes to use: RColorBrewer "Reds" (9) and "Blues" (9)
  reds  <- c("#FFF5F0", "#FEE0D2", "#FCBBA1", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D", "#A50F15", "#67000D")
  blues <- c("#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#08519C", "#08306B")

  na.col=gray(0.85)

  # In case of monthly data, convert the 3D arrays to 2D matrices
  I <- z$nsite
  J <- z$nyear
  M <- z$nmonth
  if (M > 1) {
    obs <- matrix( aperm(z$f,       c(1,3,2)), I, J*M)
    imp <- matrix( aperm(z$imputed, c(1,3,2)), I, J*M)
    fit <- matrix( aperm(z$mu,      c(1,3,2)), I, J*M)
  } else {
    obs <- z$f
    imp <- z$imputed
    fit <- z$mu
  }

  what <- match.arg(what)
  if (what=="data") {
    img <- .colorize(obs, reds, log=log, na.col=na.col)
  } else if (what=="imputed") {
    img <- .colorize(imp, reds, blues, is.na(obs), log=log, na.col=na.col)
  } else if (what=="fitted") {
    img <- .colorize(fit, blues, log=log, na.col=na.col)
  } else stop("Can't happen.")

  # and plot it
  nc <- dim(img)[2] # number of columns
  xx = c(0, nc)+0.5
  xx = c(min(z$time.id)-0.5, max(z$time.id)+0.5)
  # if (nc != J)                 stop("Problem with # of time points")
  # if (nc != length(z$time.id)) stop("Problem with # of time points")
  # xx <- range(z$time.id) + 0.5
  if (xlab=="auto") xlab <- ifelse(z$time.id[1]==1, "Time point", "Year")

  nr <- dim(img)[1] # number of rows ...
  yy = c(0, nr)+0.5

  plot(xx,yy, type='n', ylim=rev(yy), xlab=xlab, ylab=ylab, ...)
  rasterImage(img, xx[1], yy[1], xx[2], yy[2], interpolate=FALSE)

  # Draw grid for monthly only
  if (M>1) {
    xgrid <- (xx[1]+1) : (xx[2]-1)
    ygrid1 <- rep(yy[1], length(xgrid))
    ygrid2 <- rep(yy[2], length(xgrid))
    segments(xgrid,ygrid1, xgrid,ygrid2, col="black", lwd=1)
  }
  rect(xx[1], yy[1], xx[2], yy[2], border="black")
}
