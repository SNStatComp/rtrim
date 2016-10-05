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
trim <- function(x,...){
  UseMethod('trim')
}

#' @rdname trim
#' @export
trim.trimcommand <- function(x,...){
  dat <- read_tdf(x)
  covars <- x$labels[x$covariates]

  if (isTRUE(x$weighting)) { wgt <- dat$weight }
  else             { wgt <- numeric(0) }

  out = trim_estimate(count=dat$count
                      , time.id = dat$time
                      , site.id = dat$site
                      , covars = dat[covars]
                      , model = x$model
                      , serialcor = x$serialcor
                      , overdisp = x$overdisp
                      , changepoints = x$changepoints
                      , stepwise = x$stepwise
                      , weights = wgt)
  out
}

#' @param formula \code{[formula]} The dependent variable (left-hand-side)
#'  is treated as the 'counts' variable. The first and second independent variable
#'  are treated as the 'time' and 'site' variable, in that specific order. All
#'  other variables are treated as covariates.
#' @param model TRIM model type.
#' @param weights \code{[numeric]} Optional vector of weights.
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
  #if (missing(weights)) weights <- rep(1,nrow(x))
  stopifnot(is.numeric(model),model %in% 2:3)
  stopifnot(isTRUE(serialcor)||!isTRUE(serialcor))
  stopifnot(isTRUE(overdisp)||!isTRUE(overdisp))
  stopifnot(isTRUE(stepwise)||!isTRUE(stepwise))

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




