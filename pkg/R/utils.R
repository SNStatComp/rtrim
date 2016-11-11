
#' Set verbosity of trim model functions
#'
#' @section Details:
#' 
#' Control how much output \code{\link{trim}} writes to the screen while
#' fitting the model. By default, \code{trim} only returns the output
#' and does not write any progress to the screen. After calling
#' \code{set_trim_verbose(TRUE)}, \code{trim} will write information
#' about running iterations and convergence to the screen during optmization.
#' 
#'
#'
#' @param verbose \code{[logical]} toggle verbosity. \code{TRUE} means: be
#' verbose, \code{FALSE} means be quiet (this is the default).
#' @family modelspec
#' @export
set_trim_verbose <- function(verbose=FALSE){
  stopifnot(isTRUE(verbose)|!isTRUE(verbose))
  options(trim_verbose=verbose)
}

# Convenience function for console output during runs
rprintf <- function(fmt,...) { if(getOption("trim_verbose")) cat(sprintf(fmt,...)) }

# Similar, but for object/summary printing
printf <- function(fmt,...) {cat(sprintf(fmt,...))}
  
