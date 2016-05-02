
#' Get or set values for TRIMcommand objects
#'
#' @param x [TRIMcommand] object, see \\code{\\link{TRIMcommand}}
#' @param value (optional) Value to assign to the slot.
#'
#' @keywords internal
#' @rdname trim-options
origin <- function(x){
  slot(x,"origin")
}

#' @rdname trim-options
`origin<-` <- function(x,value){
  slot(x,"origin") <- value
}

#' @rdname trim-options
data_file <- function(x){
  slot(x,"data_file")
}

#' @rdname trim-options
`data_file<-` <- function(x,value){
  slot(x,"data_file") <- value
}

#' @rdname trim-options
title <- function(x){
  slot(x,"title")
}

#' @rdname trim-options
`title<-` <- function(x,value){
  slot(x,"title") <- value
}

#' @rdname trim-options
ntimes <- function(x){
  slot(x,"ntimes")
}

#' @rdname trim-options
`ntimes<-` <- function(x,value){
  slot(x,"ntimes") <- value
}

#' @rdname trim-options
ncovars <- function(x){
  slot(x,"ncovars")
}

#' @rdname trim-options
`ncovars<-` <- function(x,value){
  slot(x,"ncovars") <- value
}

#' @rdname trim-options
labels <- function(x){
  slot(x,"labels")
}

#' @rdname trim-options
`labels<-` <- function(x,value){
  slot(x,"labels") <- value
}

#' @rdname trim-options
weight <- function(x){
  slot(x,"weight")
}

#' @rdname trim-options
`weight<-` <- function(x,value){
  slot(x,"weight") <- value
}

#' @rdname trim-options
comment <- function(x){
  slot(x,"comment")
}

#' @rdname trim-options
`comment<-` <- function(x,value){
  slot(x,"comment") <- value
}

#' @rdname trim-options
weighting <- function(x){
  slot(x,"weighting")
}

#' @rdname trim-options
`weighting<-` <- function(x,value){
  slot(x,"weighting") <- value
}

#' @rdname trim-options
serialcor <- function(x){
  slot(x,"serialcor")
}

#' @rdname trim-options
`serialcor<-` <- function(x,value){
  slot(x,"serialcor") <- value
}

#' @rdname trim-options
overdisp <- function(x){
  slot(x,"overdisp")
}

#' @rdname trim-options
`overdisp<-` <- function(x,value){
  slot(x,"overdisp") <- value
}

#' @rdname trim-options
basetime <- function(x){
  slot(x,"basetime")
}

#' @rdname trim-options
`basetime<-` <- function(x,value){
  slot(x,"basetime") <- value
}

#' @rdname trim-options
model <- function(x){
  slot(x,"model")
}

#' @rdname trim-options
`model<-` <- function(x,value){
  slot(x,"model") <- value
}

#' @rdname trim-options
covariates <- function(x){
  slot(x,"covariates")
}

#' @rdname trim-options
`covariates<-` <- function(x,value){
  slot(x,"covariates") <- value
}

#' @rdname trim-options
changepoints <- function(x){
  slot(x,"changepoints")
}

#' @rdname trim-options
`changepoints<-` <- function(x,value){
  slot(x,"changepoints") <- value
}

#' @rdname trim-options
stepwise <- function(x){
  slot(x,"stepwise")
}

#' @rdname trim-options
`stepwise<-` <- function(x,value){
  slot(x,"stepwise") <- value
}

#' @rdname trim-options
outputdata_files <- function(x){
  slot(x,"outputdata_files")
}

#' @rdname trim-options
`outputdata_files<-` <- function(x,value){
  slot(x,"outputdata_files") <- value
}

#' @rdname trim-options
run <- function(x){
  slot(x,"run")
}

#' @rdname trim-options
`run<-` <- function(x,value){
  slot(x,"run") <- value
}
