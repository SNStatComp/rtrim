# TRIM postprocessing functions

#===============================================================================
#                                                                   Data summary
#===============================================================================

#-------------------------------------------------------------------------------
#                                                                        extract

summary.TRIMdata <- function(x)
{
  stopifnot(class(x)=="TRIMdata")

  # Collect covariate data
  #   covars <- list()
  #   covar_cols = ifelse(weight, 5,4) : ncol(x$df)
  #   for (col in covar_cols) {
  #     name <- names(x$df)[col]
  #     nlevels = nlevels(x$df[[col]])
  #   }
  covar_cols <- ifelse(x$weight, 5,4) : ncol(x$df)
  covars <- data.frame(col = covar_cols,
                       name = names(x$df)[covar_cols],
                       levels = sapply(x$df[covar_cols], nlevels)
  )

  # Create summary output
  out <- list(ncols=ncol(x$df), file=x$file, nsite=x$nsite, ntime=x$ntime,
              missing=x$missing, covars=covars,
              nzero=x$nzero, npos=x$npos, nobs=x$nobs, nmis=x$nmis,
              ncount=x$ncount, totcount=x$totcount)
  class(out) <- "summary.TRIMdata"
  out
}

dominant_sites <- function(x, threshold=10) {
  stopifnot(class(x)=="TRIMdata")

  # Compute site totals
  ST <- ddply(x$df, .(site), summarize, total=sum(count, na.rm=TRUE))
  ST$percent <- 100 * ST$total / x$totcount

  # Dominate sites: more than 10% of total observations
  DOM <- subset(ST, percent>threshold)

  # Create summary output
  out <- list(sites=DOM, threshold=threshold)
  class(out) <- "trim.dom"
  out
}

average <- function(x)
{
  stopifnot(class(x)=="TRIMdata")

  # Number of observations and mean count for each  time point can be computed
  # directly using the plyr tools.
  out <- ddply(x$df, .(time), summarise,
               observations=sum(is.finite(count)),
               average=mean(count, na.rm=TRUE))
  out$index <- out$average / out$average[1]
  out
}

#-------------------------------------------------------------------------------
#                                                                          print

print.summary.TRIMdata <- function(x)
{
  stopifnot(class(x)=="summary.TRIMdata")

  printf("\nThe following %d variables have been read from file: %s\n", x$ncols, x$file)
  printf("1. %-20s number of values: %d\n", "Site", x$nsite)
  printf("2. %-20s number of values: %d\n", "Time", x$ntime)
  printf("3. %-20s missing = %d\n", "Count", x$missing)
  # TODO: weight column
  for (i in 1:nrow(x$covars)) {
    printf("%d. %-20s number of values: %d\n", x$covars$col[i], x$covars$name[i], x$covars$levels[i])
  }

  printf("\n")
  printf("Number of observed zero counts      %8d\n", x$nzero)
  printf("Number of observed positive counts  %8d\n", x$npos)
  printf("Total number of observed counts     %8d\n", x$nobs)
  printf("Number of missing counts            %8d\n", x$nmis)
  printf("Total number of counts              %8d\n", x$ncount)
  printf("\n")
  printf("Total count                         %8d\n", x$totcount)

}

print.trim.dom <- function(dom) {
  stopifnot(class(dom)=="trim.dom")

  printf("\nSites containing more than %d%% of the total count:\n", dom$threshold)
  print(dom$sites, row.names=FALSE)
}

#===============================================================================
#                                                                Goodness-of-fit
#===============================================================================

#-------------------------------------------------------------------------------
#                                                                        extract


gof <- function(x) UseMethod("gof")

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
#' and \code{AIC}, for Chi-squared, Likelihoof Ratio and Akaike informatiuon content,
#' respectively.
#' @export
#'
#' @examples
#' z <- trim(tcf)
#' gof(z) # prints information on the goodness of fit.
#' LR_p <- gof(z)$LR$p # get p-value for likelihood ratio
gof.trim <- function(x) {
  stopifnot(class(x)=="trim")
  structure(list(chi2 = x$chi2, LR=x$LR, AIC=x$AIC), class="trim.gof")
}

#-------------------------------------------------------------------------------
#                                                                          print

print.trim.gof <- function(x) {
  stopifnot(class(x)=="trim.gof")
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
#' @param x TRIM output structure (i.e., output of a call to \code{trim})
#'
#' @return a list of type 'trim.summary' containing elements for estimation method
#'   (\code{est.method}), overdispersion (\code{sig2}) and autocorrelation (\code{rho})
#' @export
#'
#' @examples
#' z <- trim(...)
#' summary(z) # prints some summary info
#' rho <- summary(z)$rho # extract autocorrelation strength
summary.trim <- function(x) {
  stopifnot(class(x)=="trim")

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

print.trim.summary <- function(x) {
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
#' @param x TRIM output structure (i.e., output of a call to \code{trim})
#'
#' @return a list containing the model type (element \code{model}, 1,2 or 3), and
#' element \code{coef}, which is a data frame containing the actual coefficients.
#' @export
#'
#' @examples
#' z <- trim(...)
#' print(coef(z)) # print all coefficients
#' mul <- coef(z)$coef[3] # Extract multiplicative coefficients.
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

print.trim.coef <- function(x) {
  stopifnot(class(x)=="trim.coef")

  printf("Parameter estimates\n")
  if (x$model==2) {
    if (names(x$coef)[1]=="from") printf("Slope for Time Intervals\n")
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
#' @examples
#' z <- trim(tcf,dat);
#' totals(z) # prints the time-totals for the imputed data
#' print(totals(z,"imputed")) # idem
#' totals(z, "both") # mimics classic TRIM
#' SE <- totals(z)$totals$std.err
totals <- function(x, which=c("imputed","model","both")) {
  stopifnot(class(x)=="trim")

  # Select output columns from the pre-computed time totals
  which <- match.arg(which)
  if (which=="model") {
    totals = x$time.totals[c(1,2,3)]
  } else if (which=="imputed") {
    totals = x$time.totals[c(1,4,5)]
  } else if (which=="both") {
    totals <- x$time.totals
  } else stop(sprintf("Invalid options which=%s", which))

  # wrap the time.index field in a list and make it an S3 class
  # (so that it becomes recognizable as a TRIM time-indices)
  structure(list(totals=totals), class="trim.totals")
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

print.trim.totals <- function(x) {
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
#' @return a structure of class \code{trim.index}, which is a list with as single
#' element the data frame \code{idx}.
#' @export
#'
#' @examples
#' z <- trim(tcf,dat);
#' index(z) # prints the time-totals for the imputed data
#' print(index(z,"imputed")) # idem
#' index(z, "both") # mimics classic TRIM
#' SE <- index(z)$idx$std.err # Extract standard error for the imputed data
#' @examples
index <- function(x, which=c("imputed","model","both")) {
  stopifnot(class(x)=="trim")

  # Select output columns from the pre-computed time totals
  which <- match.arg(which)
  if (which=="model") {
    idx = x$time.index[c(1,2,3)]
  } else if (which=="imputed") {
    idx = x$time.index[c(1,4,5)]
  } else if (which=="both") {
    idx <- x$time.index
  } else stop(sprintf("Invalid options which=%s", which))

  # wrap the time.index field in a list and make it an S3 class
  # (so that it becomes recognizable as a TRIM time-indices)
  structure(list(idx=idx), class="trim.index")
}

#-------------------------------------------------------------------------------
#                                                                         Export

export <- function(x, species, stratum) UseMethod("export")

export.trim.index <- function(x, species, stratum) {
  stopifnot(class(x)=="trim.index")

  # Create extra columns to be put before the actual time indices
  df1 = data.frame(species=species, stratum=stratum)
  df2 = x$idx
  df = cbind(df1, df2)
  print(df, row.names=FALSE)
}

#-------------------------------------------------------------------------------
#                                                                         Print

print.trim.index <- function(x) {
  stopifnot(class(x)=="trim.index")
  printf("Time indices\n")
  print(x$idx, row.names=FALSE)
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
#' @examples
#' z <- trim(...)
#' print(linear(z)) # print coefficients of the linear trend
#' slope <- linear(z)$trend$Multiplicative # get linear trend magnitude
#' devs  <- linear(z)$dev$Multiplicative   # ... and the deviations from that trend
linear <- function(x) {
  stopifnot(class(x)=="trim")
  stopifnot(x$model==3)

  structure(list(trend=x$linear.trend, dev=x$deviations), class="trim.linear")
}

#-------------------------------------------------------------------------------
#                                                                          print

print.trim.linear <- function(x) {
  stopifnot(class(x)=="trim.linear")

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
#' @examples
#' z2 <- trim(..., model=2)
#' print(wald(z2))  # print info on significance of slope parameters
#' z3 <- trim(..., model=3)
#' print(wald(z3))  # print info on significance of deviations from linear trend
wald <- function(x) {
  structure(x$wald, class="trim.wald")
}

#-------------------------------------------------------------------------------
#                                                                          Print

print.trim.wald <- function(x) {
  stopifnot(class(x)=="trim.wald")

  if (x$model==2) {
    if (length(W)==1) {
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
#' @param which Selector to choose between "imputed" or "model" counts
#'
#' @return a list containing information overall slope coefficients (\code{coef}),
#'   the p-value of the overall slope (\code{p}),
#'   the size-effect (\code{effect}),
#' @export
#'
#' @examples
overall <- function(x, which=c("imputed","model")) {
  stopifnot(class(x)=="trim")
  which = match.arg(which)
  if (which=="imputed") {
     out=structure(x$overall$imp, class="trim.overall", src="imputed")
  } else if (which=="model") {
     out=structure(x$overall$mod, class="trim.overall", src="model")
  }
  out
}

#-------------------------------------------------------------------------------
#                                                                          Print

print.trim.overall <- function(x) {
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

plot.trim.overall <- function(X, imputed=TRUE, ...) {
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

plot.trim <- function(x) {
  #Prepare for plotting. First convert the estimations $\mu$ to a data frame.
  #Because these estimations are for specific time points, we call it the `discrete' model.
  Discrete <- rbind(
    mat2df(x$data,  "Data"),
    mat2df(x$mu,    "Model")
  )
  Discrete <- subset(Discrete, is.finite(Count))

  # Similarly create a data frame for estimations applied to continuous time (the 'continuous' model)
  if (x$model!=3) {
    ctime <- seq(1, x$ntime, length.out=100)
    CModel <- data.frame()
    for (i in 1:x$nsite) {
      if (x$model==1) {
        tmp <- data.frame(Site=i, Time=ctime, Count=exp(out$alpha[i]))
      } else if (out$model==2) {
        tmp <- data.frame(Site=i, Time=ctime, Count=exp(out$alpha[i] + out$beta*(ctime-1.0)))
      }
      CModel <- rbind(CModel, tmp)
    }
    CModel$Site <- factor(CModel$Site)
  }

  # Plot data and model (both discrete and continuous)
  g <- ggplot(Discrete, aes(x=Time, y=Count, colour=Site)) + theme_bw()
  g <- g + geom_point(aes(shape=Source), size=4)
  g <- g + scale_shape_manual(values=c(1,20))
  if (model!=3) g <- g + geom_path(data=CModel)
  g <- g + scale_x_continuous(breaks=1:x$ntime)
  g <- g + labs(shape="")
  if (nchar(title)>0) g <- g + labs(title=title)

  # Add the overall trend (based on the imputed)
  intercept <- x$overall.slope$imp$coef[[1]][1]
  slope     <- x$overall.slope$imp$coef[[1]][2]
  t <- seq(1, x$ntime, length.out=100)
  y <- exp(intercept + slope*t)
  trend <- data.frame(Time=t, Count=y, Site=NA)
  g <- g + geom_path(data=trend)

  print(g)

}


#-------------------------------------------------------------------------------
#                                                                       Plotting


plot_data_df <- function(df, title="") {
  ntime <- max(df$Time)
  # remove NA rows
  df <- subset(df, is.finite(Count))
  g <- ggplot(df, aes(x=Time,y=Count,colour=Site)) + theme_bw()
  g <- g + geom_path(linetype="dashed")
  g <- g + geom_point(size=4, shape=20)
  g <- g + scale_x_continuous(breaks=1:ntime)
  if (nchar(title)) g <- g + labs(title=title)
  print(g)
}


plot_data_mat <- function(m, ...) {
  df <- mat2df(m)
  plot_data_df(df, ...)
}


plot_output <- function(out, title="") {
}

#plot_output(Model)

