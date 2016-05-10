
#' Get or set values for TRIMcommand objects
#'
#' @param x [TRIMcommand] object, see \\code{\\link{TRIMcommand}}
#' @param value (optional) Value to assign to the slot.
#'
#' @keywords internal
#' @rdname trim-options
trim_origin <- function(x){
  slot(x,"origin")
}

#' @rdname trim-options
`trim_origin<-` <- function(x,value){
  slot(x,"origin") <- value
}

#' @rdname trim-options
trim_file <- function(x){
  slot(x,"file")
}

#' @rdname trim-options
`trim_file<-` <- function(x,value){
  slot(x,"file") <- value
}

#' @rdname trim-options
trim_title <- function(x){
  slot(x,"title")
}

#' @rdname trim-options
`trim_title<-` <- function(x,value){
  slot(x,"title") <- value
}

#' @rdname trim-options
trim_ntimes <- function(x){
  slot(x,"ntimes")
}

#' @rdname trim-options
`trim_ntimes<-` <- function(x,value){
  slot(x,"ntimes") <- value
}

#' @rdname trim-options
trim_ncovars <- function(x){
  slot(x,"ncovars")
}

#' @rdname trim-options
`trim_ncovars<-` <- function(x,value){
  slot(x,"ncovars") <- value
}

#' @rdname trim-options
trim_labels <- function(x){
  slot(x,"labels")
}

#' @rdname trim-options
`trim_labels<-` <- function(x,value){
  slot(x,"labels") <- value
}

#' @rdname trim-options
trim_weight <- function(x){
  slot(x,"weight")
}

#' @rdname trim-options
`trim_weight<-` <- function(x,value){
  slot(x,"weight") <- value
}

#' @rdname trim-options
trim_comment <- function(x){
  slot(x,"comment")
}

#' @rdname trim-options
`trim_comment<-` <- function(x,value){
  slot(x,"comment") <- value
}

#' @rdname trim-options
trim_weighting <- function(x){
  slot(x,"weighting")
}

#' @rdname trim-options
`trim_weighting<-` <- function(x,value){
  slot(x,"weighting") <- value
}

#' @rdname trim-options
trim_serialcor <- function(x){
  slot(x,"serialcor")
}

#' @rdname trim-options
`trim_serialcor<-` <- function(x,value){
  slot(x,"serialcor") <- value
}

#' @rdname trim-options
trim_overdisp <- function(x){
  slot(x,"overdisp")
}

#' @rdname trim-options
`trim_overdisp<-` <- function(x,value){
  slot(x,"overdisp") <- value
}

#' @rdname trim-options
trim_basetime <- function(x){
  slot(x,"basetime")
}

#' @rdname trim-options
`trim_basetime<-` <- function(x,value){
  slot(x,"basetime") <- value
}

#' @rdname trim-options
trim_model <- function(x){
  slot(x,"model")
}

#' @rdname trim-options
`trim_model<-` <- function(x,value){
  slot(x,"model") <- value
}

#' @rdname trim-options
trim_covariates <- function(x){
  slot(x,"covariates")
}

#' @rdname trim-options
`trim_covariates<-` <- function(x,value){
  slot(x,"covariates") <- value
}

#' @rdname trim-options
trim_changepoints <- function(x){
  slot(x,"changepoints")
}

#' @rdname trim-options
`trim_changepoints<-` <- function(x,value){
  slot(x,"changepoints") <- value
}

#' @rdname trim-options
trim_stepwise <- function(x){
  slot(x,"stepwise")
}

#' @rdname trim-options
`trim_stepwise<-` <- function(x,value){
  slot(x,"stepwise") <- value
}

#' @rdname trim-options
trim_outputfiles <- function(x){
  slot(x,"outputfiles")
}

#' @rdname trim-options
`trim_outputfiles<-` <- function(x,value){
  slot(x,"outputfiles") <- value
}

#' @rdname trim-options
trim_run <- function(x){
  slot(x,"run")
}

#' @rdname trim-options
`trim_run<-` <- function(x,value){
  slot(x,"run") <- value
}
