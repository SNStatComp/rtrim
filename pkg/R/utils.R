
 #' Set verbosity of trim model functions
 #'
 #' @param verbose \code{[logical]} toggle verbosity.
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
  
