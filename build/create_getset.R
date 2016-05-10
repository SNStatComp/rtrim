
devtools::load_all('pkg')
library(whisker)

get_template <- "
#' @rdname trim-options
trim_{{slot}} <- function(x){
  slot(x,\"{{slot}}\")
}
"

set_template <- "
#' @rdname trim-options
`trim_{{slot}}<-` <- function(x,value){
  slot(x,\"{{slot}}\") <- value
}
"

start_block <- "
#' Get or set values for TRIMcommand objects
#'
#' @param x [TRIMcommand] object, see \\\\code{\\\\link{TRIMcommand}}
#' @param value (optional) Value to assign to the slot.
#'
#' @keywords internal"
cat(start_block, file="build/getset.R")

slot_names <- slotNames("TRIMcommand")
for ( nm in slot_names ){
  cat(whisker.render(get_template,list(slot=nm)), file="build/getset.R",append=TRUE)
  cat(whisker.render(set_template,list(slot=nm)),file="build/getset.R",append=TRUE)
}

file.copy("build/getset.R","pkg/R",overwrite=TRUE)


