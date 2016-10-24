#' Estimate TRIM model parameters
#'
#'
#' Compute TRIM model parameters as specified in the
#' \href{https://www.cbs.nl/NR/rdonlyres/2E9912EB-534B-4A32-AD22-17A73402C083/0/trim3man.pdf}{TRIM3 manual}.
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
#'
#' @param x a \code{\link{trimcommand}}, a \code{data.frame}, or a \code{formula}
#' @param ... Currently unused
#'
#'
#' @export
#'
#' @family modelspec
#' @seealso \href{../doc/rtrim_for_TRIM_users.html}{rtrim for TRIM users}, \code{\link{summary.trim}}
#'
#' @examples
#' data(skylark)
#' m <- trim(count ~ time + site, data=skylark, model=2)
#' gof(m)
#'
#' # An example using weights
#' # set up some random weights (one for each site)
#' w <- runif(55, 0.1, 0.9)
#' # match weights to sites
#' weights <- w[skylark$site]
#' # run model
#' m <- trim(count ~ time + site, data=skylark, model=3)
#'
#' # An example using change points, a covariate, and overdispersion
#' # 1 is added as cp automatically
#' cp <- c(2,6)
#' m <- trim(count ~ time + site + Habitat, data=skylark, model=2, changepoints=cp, overdisp=TRUE)
#' plot(overall(m))
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
  else             { wgt <- numeric(0) }

  if (isTRUE(x$covin)) covin <- read_icv(x)
  else                 covin <- list()

  out <- trim_estimate(count=dat$count
                      , time.id = dat$time
                      , site.id = dat$site
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

#' @param formula \code{[formula]} The dependent variable (left-hand-side)
#'  is treated as the 'counts' variable. The first and second independent variable
#'  are treated as the 'time' and 'site' variable, in that specific order. All
#'  other variables are treated as covariates.
#' @param model \code{[numeric]} TRIM model type 1, 2, or 3.
#' @param weights \code{[numeric]} Optional vector of site weigts. The length of
#' \code{weights} must be equal to the number of rows in the data.
#' @param serialcor \code{[logical]} Take serial correlation into account.
#' @param overdisp \code{[logical]} Take overdispersion into account.
#' @param changepoints \code{[numeric]} Indices for changepoints.
#' @param stepwise \code{[logical]} Perform stepwise refinement of changepoints.
#' @param autodelete \code{[logical]} Auto-delete changepoints when number of observations is too small.
#'
#' @rdname trim
#' @export
trim.data.frame <- function(x, formula, model = 2, weights=numeric(0)
  , serialcor=FALSE, overdisp=FALSE, changepoints=integer(0), stepwise=FALSE
  , autodelete=FALSE, ...){

  # argument parsing
  L <- parse_formula(formula,vars=names(x))

  stopifnot(is.numeric(model),model %in% 2:3)
  stopifnot(isTRUE(serialcor)||!isTRUE(serialcor))
  stopifnot(isTRUE(overdisp)||!isTRUE(overdisp))
  stopifnot(isTRUE(stepwise)||!isTRUE(stepwise))
  stopifnot(all(weights>0), length(weights) != nrow(x))

  # estimate the model and return
  m <- trim_estimate(
    count = x[[L$count]]
    , time.id = x[[L$time]]
    , site.id = x[[L$site]]
    , covars = x[L$cov]
    , model = model
    , serialcor=serialcor
    , overdisp=overdisp
    , changepoints = changepoints
    , stepwise = stepwise
    , autodelete = autodelete
    , weights = weights
  )
}

#' @rdname trim
#' @param data \code{[data.frame]} Data containing at least counts, times, and sites.
#' @export
trim.formula <- function(x, data, model=c(1,2,3), weights=numeric(0)
          , serialcor=FALSE, overdisp=FALSE, changepoints=integer(0), stepwise=FALSE
          , autodelete=FALSE, ...){
  stopifnot(inherits(data,"data.frame"))
  trim.data.frame(x=data, formula=x, model=model, weights=weights
      , serialcor=serialcor, overdisp=overdisp, changepoints=changepoints
      , stepwise=stepwise, autodelete=autodelete)
}


parse_formula <- function(x, vars){
  lhs <- all.vars(x[[2]])
  if ( length(lhs) != 1)
    stop(sprintf("Expected precisely one dependent variable, got %s",pr(lhs)))
  rhs <- all.vars(x[[3]])
  if ( length(rhs) < 2 )
    stop(sprintf("Expected at least two dependent variables, got %s",pr(rhs)))
  all_vars <- c(lhs,rhs)
  valid_vars <- all_vars %in% vars
  if (!all(valid_vars)){
    stop(sprintf("Variables %s not found in data"),pr(all_vars[!valid_vars]))
  }
  list(count = lhs, time = rhs[1], site=rhs[2], cov=rhs[-(1:2)])
}




