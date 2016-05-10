#' @include tcf.R


#' Compute the TRIM model
#' 
#' 
#' @description Compute the TRIM model
#' 
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
  dat <- read_tdf(trim_file(x))
  do_trim(dat=dat, cmd=x)
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



# merge optins in ... with TRIMcommand object cmd
# an error is generated when arguments in ... are not slots
# in TRIMcommand. Assigned values will be coerced to the correct
# type when possible.
tc_merge <- function(cmd,...){
  slot_types <- getSlots("TRIMcommand")
  slot_names <- names(slot_types)
  L <- list(...)
  
  in_names <- names(L)
  
  ii <- in_names %in% slot_names
  if ( !all(ii) ){
    w <- sprintf("Arguments %s are not valid TRIMcommand slots and will be ignored",
            paste(in_names[!ii],collapse=", "))
    warning(w)
    in_names <- in_names[ii]
  }
  for ( n in in_names ){
    slot(cmd, n) <- as(L[n],slot_types[n])
  }
  cmd
}




