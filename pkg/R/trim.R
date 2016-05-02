#' @include tcf.R


#' Compute the TRIM model
#' 
#' @name trim
#' @title Compute TRIM model
#' 
#' 
#' @param x \code{[character|TRIMcommand|data.frame]} An R-object (see Details)
#' @param y [optional] R-object
#' @param ... Extra parameters, to be passed to underlying methods.
#' 
#' 
#' @export
setGeneric("trim",function(x,y,...) standardGeneric("trim"))


# methods for the TRIM-like workflow

#' @rdname trim
setMethod("trim",signature = c("TRIMcommand","ANY"), function(x,y=NULL,...){
  x <- tc_merge(x,...)
  dat <- read_tdf(data_file(x))
  do_trim(dat=dat, cmd=cmd)
  ## possibly do some output object building.
})

#' @rdname trim
setMethod("trim",signature=c("character","ANY"),function(x,y=NULL,...){
  cmd <- read_tcf(x,...)
  getMethod("trim",c("TRIMcommand","ANY"))(cmd,y,...)
})

#' @rdname trim
setMethod("trim",signature=c("data.frame","TRIMcommand"), function(x,y,...){
  y <- tc_merge(y,...)
  do_trim(dat=x, cmd=y)
  ## possibly some output object building
})



tc_merge <- function(cmd,...){
  # merge optins in ... with TRIMcommand object cmd
  cmd
}




