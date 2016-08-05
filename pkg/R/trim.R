#' Estimate TRIM model
#'
#' @param x a \code{\link{trimcommand}} object or a \code{data.frame}
#'
#' @export
trim <- function(x,...){
  UseMethod('trim')
}

#' @rdname trim
#' @export
trim.trimcommand <- function(x,...){
  dat <- read_tdf(x)
  trim_estimate(count=dat$count
      , time.id = dat$time
      , site.id = dat$site
      , model = x$model
      , serialcor = x$serialcor
      , overdisp = x$overdisp
      , changepoints = x$changepoints)
}




