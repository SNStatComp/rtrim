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
#' \eqn{\ln\mu_{ij}=\alpha_i + \beta_j},
#'
#' where \eqn{\beta_j} is the deviatiation of log-counts at time \eqn{j},
#' averaged over all sites. To make this model identifiable, the value of
#' \eqn{\beta_1=0} by definition. Model 3 can be shown to be equivalent to
#' Model 2 with a changepoint at every time point. Using a
#' \code{\link[=wald.trim]{wald}} test, one can estimate whether the collection
#' of deviations \eqn{\beta_i} make the model differ significantly from an
#' overall linear trend (Model 2 without changepoints).
#'
#' The parameters \eqn{\alpha_i} and \eqn{\gamma_j} are referred to
#' as the \emph{additive representation} of the coefficients. Once computed,
#' they can be represented and extracted in several representations, using the
#' \code{\link[=coef.trim]{coefficients}} function. (See also the examples
#' below).
#'
#' Other model parameters can be extracted using functions such as
#' \code{\link{gof}} (for goodness of fit), \code{\link[=summary.trim]{summary}}
#' or \code{\link{totals}}. Refer to the `See also' section for an overview.
#'
#' @section {Using yearly and monthly counts}
#'
#' In many data sets will use use only yearly count data, in which case the
#' time \eqn{j} will reflect the year number.
#' An extension of \code{trim} is to use monthly (or any other sub-yearly) count data,
#' in combination with index computations on the yearly time scale.
#'
#' In this case, counts are given as \eqn{f_{i,j,m}} with \eqn{m=1,2,\ldots,M} the month number.
#' As before, \eqn{\mu_{i,j,m}} will be imputed in case of missing counts.
#'
#' The contibution of month factors to the model is always similar to the way year factors are used in Model 3,
#' that is,
#'
#' \eqn{\ln\mu_{i,j,m} = \alpha_i + \beta\times(j-1) + \gamma_m}
#' for Model 2, and
#'  \eqn{\ln\mu_{i,j,m} = \alpha_i + \beta_j + \gamma_m}
#' for Model 3.
#'
#' For the same reason why \eqn{\beta_1=0} for Model 3, \eqn{\gamma_1=0} in case of montly parameters.
#'
#' #' @section {Using covariates}:
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
#'
#' For \bold{Model 3 with covariates} the parameter \eqn{\beta_j} is replaced by
#'
#' \eqn{\beta_{j0} + \sum_{k=1}^Kz_{ijk}\beta_{jk}.}
#'
#' Again, the \eqn{\beta_{j0}} are referred to as baseline parameters and the
#' \eqn{\beta_{jk}} record mean differences in log-counts within a set of sites
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
#' \item{For montly data, there must be at least one observation for every month.}
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
#' m <- trim(count ~ site + time, data=skylark, model=2)
#' summary(m)
#' coefficients(m)
#'
#' # An example using weights
#' # set up some random weights (one for each site)
#' w <- runif(55, 0.1, 0.9)
#' # match weights to sites
#' skylark$weights <- w[skylark$site]
#' # run model
#' m <- trim(count ~ site + time, data=skylark, weights="weights", model=3)
#'
#' # An example using change points, a covariate, and overdispersion
#' # 1 is added as cp automatically
#' cp <- c(2,6)
#' m <- trim(count ~ site + time + Habitat, data=skylark, model=2, changepoints=cp, overdisp=TRUE)
#' coefficients(m)
#' # check significance of changes in slope
#' wald(m)
#' plot(overall(m))
#'
trim <- function(x, ...){
  UseMethod('trim')
}

#' @rdname trim
#' @export
trim.trimcommand <- function(tcf, ...) {
  call <- sys.call()

  dat <- read_tdf(tcf)
  covars <- tcf$labels[tcf$covariates]

  if (isTRUE(tcf$weighting)) { wgt <- dat$weight }
  else                       { wgt <- NULL }

  if (isTRUE(tcf$covin)) covin <- read_icv(x)
  else                   covin <- list()

  # Create 'automatic' changepoint #1
  if (tcf$model==2 && length(tcf$changepoints)==0) tcf$changepoints=1L

  trim_estimate(count=dat$count
      , site=dat$site
      , year=dat$time
      , month=NULL
      , weights=wgt
      , covars=dat[covars]
      , model=tcf$model
      , changepoints=tcf$changepoints
      , overdisp=tcf$overdisp
      , serialcor=tcf$serialcor
      , stepwise=tcf$stepwise
      , autodelete=tcf$autodelete
      , ...)
}



#' @param count_col \code{[character]} name of the column holding species counts
#' @param site_col \code{[character]} name of the column holding the site id
#' @param year_col \code{[character]} name of the column holding the time of counting
#' @param month_col \\code{[character]} optional name of the column holding the season of counting
#' @param model \code{[numeric]} TRIM model type 1, 2, or 3.
#' @param weights_col \code{[numeric]} Optional vector of site weights. The length of
#' @param covar_cols \code{[character]} name(s) of column(s) holding covariates
#' @param serialcor \code{[logical]} Take serial correlation into account (See `Estimation details')
#' @param overdisp \code{[logical]} Take overdispersion into account (See `Estimation options').
#' @param changepoints \code{[numeric]} Indices for changepoints (`Models').
#' @param stepwise \code{[logical]} Perform stepwise refinement of changepoints.
#' @param autodelete \code{[logical]} Auto-delete changepoints when number of observations is too small. (See
#'  `Demands on data').
#'
#' @rdname trim
#' @export
trim.data.frame <- function(data, count_col="count", site_col="site", year_col="year", month_col=NULL
                            , weights_col=NULL, covar_cols=NULL
                            , model=2, changepoints=ifelse(model==2, 1L, integer(0))
                            , overdisp=FALSE, serialcor=FALSE
                            , stepwise=FALSE, autodelete=TRUE, ...) {

  # check data source
  stopifnot(class(data)=="data.frame")
  stopifnot(nrow(data)>0)

  # check data columns
  stopifnot(count_col %in% names(data))
  count <- data[[count_col]]

  stopifnot(site_col %in% names(data))
  site <- data[[site_col]]

  stopifnot(year_col %in% names(data))
  year <- data[[year_col]]

  if (is.null(month_col)) {
    month <- NULL
  } else {
    stopifnot(month_col %in% names(data))
    month <- data[[month_col]]
  }

  if (is.null(weights_col)) {
    weights <- NULL
  } else {
    stopifnot(weights_col %in% names(data))
    weights <- data[[weights_col]]
  }

  if (is.null(covar_cols)) {
    covars <- data.frame()
  } else if (length(covar_cols)==0) {
    covars <- data.frame()
  } else {
    for (covar_col in covar_cols) {
      stopifnot(covar_col %in% names(data))
    }
    covars <- data[covar_cols]
  }

  # Check model specification and parameters

  stopifnot(is.numeric(model), model %in% 1:4)

  if (length(changepoints)>0 && is.na(changepoints)) changepoints <- integer(0) # fix

  stopifnot(is.logical(overdisp))
  stopifnot(is.logical(serialcor))
  stopifnot(is.logical(autodelete))
  stopifnot(is.logical(stepwise))

  # proceed with the (internal) workhorse function
  # browser()
  trim_estimate(count=count, site=site, year=year, month=month, weights=weights,
                covars=covars, model=model, changepoints=changepoints,
                overdisp=overdisp, serialcor=serialcor,
                stepwise=stepwise, autodelete=autodelete,
                ...)
}


#' @rdname trim
#' @param data \code{[data.frame]} Data containing at least counts, sites, and times
#' @export
trim.formula <- function(f, data=NULL, weights=NULL, monthly=FALSE, ...)
{
  # internal function
  unpack_formula <- function(f) {
    if (length(f)==1) {
      out <- as.character(f)
    } else if (length(f)==3) {
      opr <- unpack_formula(f[[1]]) # operator
      op1 <- unpack_formula(f[[2]]) # operand 1
      op2 <- unpack_formula(f[[3]]) # operand 2
      out <- c(op1, opr, op2)
    } else stop(sprintf("Unexpected: formula (element) of length %s", length(f)))
    out
  }

  # argument 'data' MUST be data frame
  if (is.null(data)) stop("no data given")
  if (class(data)!="data.frame") stop("data should be a data frame")

  # first check if all elements of the formula are known.
  vars <- all.vars(f)
  known <- vars %in% names(data)
  if (!all(known)) {
    unknown <- vars[!known]
    msg <- sprintf("Model contains unkown variable: %s", paste(unknown, collapse=", "))
    stop(msg)
  }

  terms <- all.names(f)
  operators <- terms[!(terms %in% vars)]
  allowed <- operators %in% c("~","+",":")
  if (!all(allowed)) {
    forbidden <- operators[!allowed]
    msg <- sprintf("Model contains unallowed operator: %s", paste(forbidden, collapse=", "))
    stop(msg)
  }

  terms <- unpack_formula(f)
  # first part should be LHS ~ site + year ...
  if (length(terms)<5)   stop("model should have form 'count ~ site + year ...")
  if (terms[[2]] != "~") stop("model should have form 'count ~ ...'")
  if (terms[[4]] != "+") stop("model should have form 'count ~ site + year ...")

  count_col <- terms[1]
  site_col  <- terms[3]
  year_col  <- terms[5]
  # printf("found count ID: %s\n", count)
  # printf("found site ID: %s\n", site)
  # printf("found year ID: %s\n", year)
  terms <- terms[-(1:5)]

  # optionally: month specifier
  if (length(terms)>0 && terms[1]==":") {
    month_col <- terms[2]
    # printf("found month ID: %s\n", month)
    terms <- terms[-(1:2)]
  }
  else month_col <- NULL

  # optionally: covariates
  covar_cols <- character(0)
  while (length(terms)>0) {
    if (terms[1]!="+") {
      msg <- sprintf("covariates should be included using the '+' operator ('%s' found)", terms[1])
      stop(msg)
    }
    covar <- terms[2]
    covar_cols <- c(covar_cols, covar)
    terms <- terms[-(1:2)]
  }
  if (length(covar_cols)==0) covar_cols <- NULL
  # if (length(covars)>0) {
  #   printf("found covars: %s\n", paste(covars, collapse=", "))
  # }


  ########################
  # # Parameter 'data' MUST be a data.frame
  # if (is.null(data)) stop("No data given")
  # stopifnot(inherits(data,"data.frame"))
  #
  # # Parse formula 'f'.
  # # As a first check, all terms of the formula must be in the data frame.
  # all_vars <- all.vars(f)
  # for (v in all_vars) {
  #   if (!v %in% names(data)) {
  #     msg <- sprintf("model term '%s' not found in data.", v)
  #     hint <- sprintf("(data columns are: '%s').", paste(names(data), collapse="','"))
  #     stop(paste(msg, hint))
  #   }
  # }
  #
  # # formula's are stored in reverse polish notation, i.e.,
  # # z ~ x + y is stored as ~(z,(x,y)) or general: (~, lhs, rhs)
  # lhs <- all.vars(f[[2]])
  # if (length(lhs)!=1) stop("Model should have 1 LHS term")
  # count_col <- lhs
  #
  # rhs <- all.vars(f[[3]])
  # nrhs <- length(rhs)
  # if (nrhs<2) stop("Model should have at least 2 RHS term")
  # site_col <- rhs[1]
  # year_col <- rhs[2]
  # nused <- 2L
  #
  # # if TRIM is used in montly mode, a third RHS value identifies the month column
  # if (monthly) {
  #   if (nrhs<3) stop(" Montly models require at least 3 RHS terms")
  #   month_col <- rhs[3]
  #   nused <- nused + 1L
  # } else {
  #   month_col <- NULL
  # }
  #
  # # Any additional columns specify covariates.
  # if (nrhs > nused) {
  #   covar_cols <- rhs[(nused+1) : nrhs]
  # } else {
  #   covar_cols <- NULL
  # }

  weights_col <- weights # rename for consistency in call below

  # pass on to trim.data.frame
  trim.data.frame(data=data, count_col=count_col, site_col=site_col, year_col=year_col,
      month_col=month_col, weights_col=weights_col, covar_cols=covar_cols,
      ...)
}
