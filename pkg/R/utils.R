
 #' Set verbosity of trim model functions
 #'
 #' @section Details:
 #' 
 #' Control how much output \code{\link{trim}} writes to the screen while
 #' fitting the model. 
 #'
 #'
 #' @param verbose \code{[logical]} toggle verbosity. \code{TRUE} means: be
 #' verbose, \code{FALSE} means be quiet (this is the default).
 #'
 #' @export
 set_trim_verbose <- function(verbose=FALSE){
   stopifnot(isTRUE(verbose)|!isTRUE(verbose))
   options(trim_verbose=verbose)
 }
 set_trim_verbose(TRUE)
 
 # Convenience function for console output during runs
 rprintf <- function(fmt,...) { if(getOption("trim_verbose")) cat(sprintf(fmt,...)) }
 
 # Similar, but for object/summary printing
 printf <- function(fmt,...) {cat(sprintf(fmt,...))}
  
