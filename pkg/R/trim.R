#' Estimate TRIM model parameters
#'
#'
#' Compute TRIM model parameters.
#'
#'
#' @section Models:
#'
#' The purpose of \code{trim} is to estimate population totals over time,
#' based on a set of counts \eqn{f_{ij}} at sites \eqn{i=1,2,\ldots,I}
#' and times \eqn{j=1,2,\ldots,J}. If no count data is available at
#' site and time \eqn{(i,j)}, a value \eqn{\mu_{ij}} will be imputed.
#'
#' In \bold{Model 1}, the imputed values are modeled as
#'
#' \eqn{\ln\mu_{ij} = \alpha_i,}
#'
#' where \eqn{\alpha_i} is the site effect. This model implies that the counts
#' vary accross sites, not over time. The model-based \link[=totals]{time totals} are equal to
#' each time point and the model-based \link[=index]{indices} are all equal to one.
#'
#' In \bold{Model 2}, the imputed values are modeled as
#'
#' \eqn{\ln\mu_{ij} = \alpha_i + \beta\times(j-1).}
#'
#' Here, \eqn{\alpha_i} is the log-count of site \eqn{i} averaged over time and
#' \eqn{\beta} is the mean growth factor that is shared by all sites over all of
#' time. The assumption of a constant growth rate may be relaxed by passing
#' a number of \code{changepoints} that indicate at what times the growth
#' rate is allowed to change. Using a \code{\link[=wald.trim]{wald}} test
#' one can investigate whether the changes in slope at the changepoints are
#' significant. Setting \code{stepwise=TRUE} makes \code{trim} automatically
#' remove changepoints where the slope does not change significantly.
#'
#' In \bold{Model 3}, the imputed values are modeled as
#'
#' \eqn{\ln\mu_{ij}=\alpha_i + \gamma_j},
#'
#' where \eqn{\gamma_j} is the deviatiation of log-counts at time \eqn{j},
#' averaged over all sites. To make this model identifiable, the value of
#' \eqn{\gamma_1=0} by definition. Model 3 can be shown to be equivalent to
#' Model 2 with a changepoint at every time point. Using a
#' \code{\link[=wald.trim]{wald}} test, one can estimate whether the collection
#' of deviations \eqn{\gamma_i} make the model differ significantly from an
#' overall linear trend (Model 2 without changepoints).
#'
#' The parameters \eqn{\alpha_i}, \eqn{\beta} and \eqn{\gamma_j} are referred to
#' as the \emph{additive representation} of the coefficients. Once computed,
#' they can be represented and extracted in several representations, using the
#' \code{\link[=coef.trim]{coefficients}} function. (See also the examples
#' below).
#'
#' Other model parameters can be extracted using functions such as
#' \code{\link{gof}} (for goodness of fit), \code{\link[=summary.trim]{summary}}
#' or \code{\link{totals}}. Refer to the `See also' section for an overview.
#'
#'
#' @section {Using covariates}:
#'
#' In the basic case of Models 2 and 3, the growth parameter \eqn{\beta} does
#' not vary accross sites. If auxiliary information is available (for instance
#' a classification of the type of soil or vegetation), the effect of these
#' variables on the per-site growth rate can be taken into account.
#'
#' For \bold{Model 2 with covariates} the growth factor \eqn{\beta} is
#' replaced with a factor
#'
#' \eqn{\beta_0 + \sum_{k=1}^K z_{ijk}\beta_k}.
#'
#' Here, \eqn{\beta_0} is referred to as the \emph{baseline} and \eqn{z_{ijk}} is a
#' dummy variable that combines dummy variables for all covariates. Since a
#' covariate with \eqn{L} classes is modeled by \eqn{L-1} dummy variables, the
#' value of \eqn{K} is equal to the sum of the numbers of categories for all
#' covariates minus the number of covariates. Observe that this model allows for
#' a covariate to change over time at a certain sites. It is therefore possible
#' to include situations for example where a site turns from farmland to rural
#' area. The \code{\link[=coef.trim]{coefficients}} function will report every
#' individual value of \eqn{\beta}. With a \code{\link[=wald.trim]{wald}} test,
#' the significance of contributions of covariates can be tested.
#'
#' For \bold{Model 3 with covariates} the parameter \eqn{\gamma_j} is replaced by
#'
#' \eqn{\gamma_{j0} + \sum_{k=1}^Kz_{ijk}\gamma_{jk}.}
#'
#' Again, the \eqn{\gamma_{j0}} are referred to as baseline parameters and the
#' \eqn{\gamma_{jk}} record mean differences in log-counts within a set of sites
#' with equal values for the covariates. All coefficients can be extracted with
#' \code{\link[=coef.trim]{coefficients}} and the significance of covariates can
#' be investigated with the \code{\link[=wald.trim]{wald}} test.
#'
#'
#' @section Estimation options:
#'
#' In the simplest case, the counts at different times and sites are considered
#' independently Poisson distributed. The (often too strict) assumption that
#' counts are independent over time may be dropped, so correlation between time
#' points at a certain site can be taken into account. The assumption of being
#' Poisson distributed can be relaxed as well. In general, the
#' variance-covariance structure of counts \eqn{f_{ij}} at site \eqn{i} for time
#' \eqn{j} is modeled as
#' \itemize{
#' \item{\eqn{\textrm{var}(f_{ij}) = \sigma^2\mu_{ij} }}
#' \item{\eqn{\textrm{cor}(f_{ij},f_{i,j+1}) = \rho} },
#' }
#' where \eqn{\sigma} is called the \emph{overdispersion}, \eqn{\mu_{ij}} is
#' the estimated count for site \eqn{i}, time \eqn{j} and \eqn{\rho} is called
#' the \emph{serial correlation}.
#'
#' If \eqn{\sigma=1}, a pure Poisson distribution is assumed to model the
#' counts. Setting \code{overdispersion = TRUE} makes \code{trim} relax this
#' condition. Setting  \code{serialcor=TRUE} allows \code{trim} to assume a
#' non-zero correlation between adjacent time points, thus relaxing the
#' assumption of independence over time.
#'
#' @section Demands on data:
#' The data set must contain sufficient counts to be able to estimate the model.
#' In particular
#' \itemize{
#' \item{For model 2 without covariates there must be at least one observation
#'   for each time segment defined by the change points.}
#' \item{For model 2 with covariates there must be at least one observation for
#'   every value of each covariate, at each time segment defined by the change
#'   points.}
#' \item{For model 3 without covariates there must be at least one observation
#'   for each time point.}
#' \item{For model 3 with covariates there must be at least one observation for
#'   every value of each covariate, at each time point.}
#' }
#'
#' The function \code{\link{check_observations}} identifies cases where too few
#' observations are present to compute a model. Setting the option
#' \code{autodelete=TRUE} (Model 2 only) makes \code{trim} remove changepoints
#' such that at each time piece sufficient counts are available to estimate the
#' model.
#'
#'
#'
#' @param x a \code{\link{trimcommand}}, a \code{data.frame}, or a \code{formula}
#'  If \code{x} is a \code{formula}, the dependent variable (left-hand-side)
#'  is treated as the 'counts' variable. The first and second independent variable
#'  are treated as the 'site' and 'time' variable, \bold{in that specific order}. All
#'  other variables are treated as covariates.
#' @param ... Currently unused
#'
#'
#' @export
#'
#' @family analyses
#' @family modelspec
#' @seealso \href{../doc/rtrim_for_TRIM_users.html}{rtrim for TRIM users}, \href{../doc/Skylark_example.html}{TRIM by example}
#'
#' @examples
#' data(skylark)
#' skylark$Habitat <- factor(skylark$Habitat) # hack
#' m <- trim(count ~ site + time, data=skylark, model=2)
#' summary(m)
#' coefficients(m)
#'
#' # An example using weights
#' # set up some random weights (one for each site)
#' w <- runif(55, 0.1, 0.9)
#' # match weights to sites
#' weights <- w[skylark$site]
#' # run model
#' m <- trim(count ~ site + time, data=skylark, model=3, weights=weights)
#'
#' # An example using change points, a covariate, and overdispersion
#' # 1 is added as cp automatically
#'
#' cp <- c(2,6)
#' m <- trim(count ~ site + time + Habitat, data=skylark, model=2, changepoints=cp, overdisp=TRUE)
#' coefficients(m)
#' # check significance of changes in slope
#' wald(m)
#' plot(overall(m))
#'
trim <- function(x,...){
  UseMethod('trim')
}

#' @rdname trim
#' @export
trim.trimcommand <- function(x,...){
  call <- sys.call()

  dat <- read_tdf(x)
  covars <- x$labels[x$covariates]

  if (isTRUE(x$weighting)) { wgt <- dat$weight }
  else                     { wgt <- numeric(0) }

  if (isTRUE(x$covin)) covin <- read_icv(x)
  else                 covin <- list()

  # Create 'automatic' changepoint #1
  if (x$model==2 && length(x$changepoints)==0) x$changepoints=1L

  trim_estimate(count=dat$count
      , site.id = dat$site
      , time.id = dat$time
      , covars = dat[covars]
      , model = x$model
      , serialcor = x$serialcor
      , overdisp = x$overdisp
      , changepoints = x$changepoints
      , stepwise = x$stepwise
      , weights = wgt
      , covin = covin
      , ...)
}



#' @param count.id \code{[character]} name of the column holding species counts
#' @param site.id \code{[character]} name of the column holding the site id
#' @param time.id \code{[character]} name of the column holding the time of counting
#' @param covars \code{[character]} name(s) of column(s) holding covariates
#' @param model \code{[numeric]} TRIM model type 1, 2, or 3.
#' @param weights \code{[numeric]} Optional vector of site weights. The length of
#' \code{weights} must be equal to the number of rows in the data.
#' @param serialcor \code{[logical]} Take serial correlation into account (See `Estimation details')
#' @param overdisp \code{[logical]} Take overdispersion into account (See `Estimation options').
#' @param changepoints \code{[numeric]} Indices for changepoints (`Models').
#' @param stepwise \code{[logical]} Perform stepwise refinement of changepoints.
#' @param autodelete \code{[logical]} Auto-delete changepoints when number of observations is too small. (See
#'  `Demands on data').
#'
#' @rdname trim
#' @export
trim.data.frame <- function(x, count.id = "count", site.id="site", time.id="time", season.id=NULL
                            , covars=character(0)
                            , model = 2
                            , weights=numeric(0)
                            , serialcor=FALSE
                            , overdisp=FALSE
                            , changepoints=ifelse(model==2, 1L, integer(0))
                            , stepwise=FALSE
                            , autodelete=FALSE, ...) {

  if (nrow(x)==0) stop("Empty data frame")

  if (length(changepoints)>0 && is.na(changepoints)) changepoints <- integer(0) # fix

  # # handle 'bare' weights specifier (column name, not as string)
  # if (class(substitute(weights))=="name") {
  #   weights <- deparse(substitute(weights))
  # }

  # Handle character-type weights specifier
  if (is.character(weights)) {
    if (!weights %in% names(x)) {
      msg <- sprintf("Weights column \"%s\" not present in data frame.", weights)
      stop(msg, call.=FALSE)
    }
    weights <- x[[weights]]
  }

  # Check covariates
  if (length(covars)>0) {
    for (covar in covars) {

      # Does the covariate exist?
      if (!covar %in% names(x)) {
        msg <- sprintf("Covariate column \"%s\" not present in data frame.", covar)
        stop(msg, call.=FALSE)
      }

      # # is is a factor? If not, raise an error
      # cls <- class(x[[covar]])
      # if (!"factor" %in% cls) {
      #   msg <- sprintf("Covariate \"%s\" is not a factor but of class \"%s\".", covar, paste(cls,collapse=","))
      #   stop(msg, call.=FALSE)
      # }

      # Is it not a factor? Convert it.
      cls <- class(x[[covar]])
      if (!"factor" %in% cls) {
        rprintf("Converting covariate \"%s\" from class \"%s\" to a factor.\n", covar, paste(cls,collapse=","))
        x[[covar]] <- factor(x[[covar]])
      }

    }
  }

  stopifnot(is.numeric(model), model %in% 1:4)
  stopifnot(isTRUE(serialcor)||!isTRUE(serialcor))
  stopifnot(isTRUE(overdisp)||!isTRUE(overdisp))
  stopifnot(isTRUE(stepwise)||!isTRUE(stepwise))
  stopifnot(all(weights>0), length(weights) == 0  || length(weights) == nrow(x))

  # estimate the model and return
  month <- if (is.null(season.id)) NULL else x[[season.id]]
  trim_estimate(
    count   = x[[count.id]],
    site.id = x[[site.id]],
    time.id = x[[time.id]],
    month = month,
    covars  = x[covars],
    model   = model,
    serialcor =serialcor, overdisp=overdisp, changepoints = changepoints,
    stepwise  = stepwise, autodelete = autodelete,
    weights   = weights,
    ...
  )
}

#' @rdname trim
#' @param data \code{[data.frame]} Data containing at least counts, sites, and times
#' @export
trim.formula <- function(x, data=NULL, model=2, weights=numeric(0)
          , serialcor=FALSE
          , overdisp=FALSE
          , changepoints=ifelse(model==2, 1L, integer(0))
          , stepwise=FALSE
          , autodelete=FALSE, ...) {
  # Parameter 'data' MUST be a data.frame
  if (is.null(data)) stop("No data given")
  stopifnot(inherits(data,"data.frame"))
  # Parse formula 'x' and make an interpretation based on the data type
  # L <- parse_formula(x, names(data))
  # Formula 'x' is organized as operator / operand1 / operand2.
  # The left-hand side of the formula is therefore the second list element
  lhs <- all.vars(x[[2]])
  if ( length(lhs) != 1)
    stop(sprintf("Expected precisely one dependent variable, got %s",pr(lhs)))
  # Idem for the right hand side, which is the 3th list element
  rhs <- all.vars(x[[3]])
  if ( length(rhs) < 2 )
    stop(sprintf("Expected at least two dependent variables, got %s",pr(rhs)))
  # Combine and check if all formula vars exist as columns in the data.frame
  all_vars   <- c(lhs,rhs)
  valid_vars <- all_vars %in% names(data)
  if (!all(valid_vars)){
    stop(sprintf("Variables %s not found in data", pr(all_vars[!valid_vars])))
  }
  # Now, unpack the formula. First 3 terms are mandatory: count ~ site + year
  count.id <- lhs
  site.id  <- rhs[1]
  time.id  <- rhs[2]
  # season/month MAY be present as a fourth term. Recognize it by its data type
  if (length(rhs)>2) {
    term <- rhs[3]
    if ("integer" %in% class(data[[term]]) || "numeric" %in% class(data[[term]])) {
      season.id <- term
      covar.id  <- rhs[-(1:3)]
    } else {
      season.id <- NULL
      covar.id  <- rhs[-(1:2)]
    }
  } else {
    season.id <- NULL
    covar.id  <- character(0)
  }

  # Fix changepoints=NA --> int(0)
  # TODO: try to remove
  if (length(changepoints)>0 && is.na(changepoints)) changepoints <- integer(0) # fix

  # pass on to trim.data.frame
  trim.data.frame(x=data
      , count.id=count.id, site.id=site.id, time.id=time.id, season.id=season.id, covars = covar.id
      , model=model, weights=weights
      , serialcor=serialcor, overdisp=overdisp, changepoints=changepoints
      , stepwise=stepwise, autodelete=autodelete, ...)
}


