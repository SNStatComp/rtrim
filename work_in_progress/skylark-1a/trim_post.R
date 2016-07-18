# TRIM postprocessing functions

#===============================================================================
#                                                                Goodness-of-fit
#===============================================================================

#-------------------------------------------------------------------------------
#                                                                        extract

# Create generic gooedness-of-fit function
gof <- function(x) UseMethod("gof")

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
#                                                                        Summary
#===============================================================================

summary.trim <- function(z) {
  stopifnot(class(z)=="trim")

  if (is.finite(z$sig2) || is.finite(z$rho)) {
    out = list(est.method="Generalised Estimating Equations")
  } else {
    out = list(est.method="Maximum Likelihood")
  }
  out$sig2 <- z$sig2
  out$rho  <- z$rho
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

coef.trim <- function(x) {
  stopifnot(class(x)=="trim")

  structure(list(model=x$model, coef=x$coefficients), class="trim.coef")
}

#-------------------------------------------------------------------------------
#                                                                          Print

print.trim.coef <- function(x) {
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
#                                                                    Time Totals
#===============================================================================

#-------------------------------------------------------------------------------
#                                                                        Extract

totals <- function(x) {
  stopifnot(class(x)=="trim")

  # wrap the time.index field in a list and make it an S3 class
  # (so that it becomes recognizable as a TRIM time-indices)
  structure(list(totals=x$time.totals), class="trim.totals")
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

index <- function(x) {
  stopifnot(class(x)=="trim")

  # wrap the time.index field in a list and make it an S3 class
  # (so that it becomes recognizable as a TRIM time-indices)
  structure(list(idx=x$time.index), class="trim.index")
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

linear <- function(x) {
  structure(list(trend=x$linear.trend, deviations=x$deviations), class="trim.linear")
}

#-------------------------------------------------------------------------------
#                                                                          print

print.trim.linear <- function(x) {
  printf("Linear Trend + Deviations for Each Time\n")
  print(x$trend,      row.names=TRUE)
  printf("\n")
  print(x$deviations, row.names=FALSE)
}

#===============================================================================
#                                                                      Wald test
#===============================================================================

#-------------------------------------------------------------------------------
#                                                                        Extract

wald <- function(x) {
  structure(x$wald, class="trim.wald")
}

#-------------------------------------------------------------------------------
#                                                                          Print

print.trim.wald <- function(x) {
  if (x$model==2) {
    printf("Wald test for significance of slope parameter\n")
  } else if (x$model==3) {
    printf("Wald test for significance of deviations from linear trend\n")
  } else stop("Can't happen")
  printf("  Wald = %.2f, df=%d, p=%f\n", x$W, x$df, x$p)
}

#===============================================================================
#                                                                  Overall slope
#===============================================================================


#-------------------------------------------------------------------------------
#                                                                        Extract

overall <- function(x, imputed=TRUE) {
  if (imputed) out=structure(x$overall$imp, class="trim.overall", src="imputed", title=x$title)
  else         out=structure(x$overall$mod, class="trim.overall", src="model",   title=x$title)
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

