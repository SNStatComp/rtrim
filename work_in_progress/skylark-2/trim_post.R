# ########################################### TRIM postprocessing functions ####

# ============================================================ Data summary ====

# ----------------------------------------------------------------- extract ----

summary.TRIMdata <- function(x)
{
  stopifnot(class(x)=="TRIMdata")

  # Collect covariate data
  ##   covars <- list()
  ##   covar_cols = ifelse(weight, 5,4) : ncol(x$df)
  ##   for (col in covar_cols) {
  ##     name <- names(x$df)[col]
  ##     nlevels = nlevels(x$df[[col]])
  ##   }
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

# ------------------------------------------------------------------- print ----

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


# ================================================================= Summary ====


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

# ============================================================ Coefficients ====

# ----------------------------------------------------------------- Extract ----

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
coef.trim <- function(x, which=c("additive","multiplicative","both")) {
  stopifnot(class(x)=="trim")
  which <- match.arg(which)
  # Last 4 columns contain the additive and multiplicative parameters.
  # Select the appropriate subset from these, and all columns before these 4.
  n = ncol(x$coefficients)
  stopifnot(n>=4)
  if (which=="additive") {
    cols <- 1:(n-2)
  } else if (which=="multiplicative") {
    cols <- c(1:(n-4), (n-1):n)
  } else if (which=="both") {
    cols = 1:n
  } else stop(sprintf("Invalid options which=%s", which))
  out <- x$coefficients[cols] # ??? does not return correct function output
  out
}

# ============================================================= Time totals ====

# ----------------------------------------------------------------- Extract ----

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

  totals
}



# =========================================== Reparameterisation of Model 3 ====

# ----------------------------------------------------------------- extract ----

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

# ------------------------------------------------------------------- print ----

print.trim.linear <- function(x) {
  stopifnot(class(x)=="trim.linear")

  printf("Linear Trend + Deviations for Each Time\n")
  print(x$trend, row.names=TRUE)
  printf("\n")
  print(x$dev, row.names=FALSE)
}



# ================================================================ Plotting ====


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
  # Prepare for plotting. First convert the estimations $\mu$ to a data frame.
  # Because these estimations are for specific time points, we call it the `discrete' model.
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


# ---------------------------------------------------------------- Plotting ----


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


