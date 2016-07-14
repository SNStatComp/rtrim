
#' Read TRIM data files
#' 
#' Read data files intended for the original TRIM programme.
#' 
#' @section The TRIM data file format:
#' 
#' TRIM input data is stored in an \code{ASCII} encoded file where columns
#' are separated by one or more spaces. In the original format, the following
#' columns and restrictions applied:
#'\tabular{lll}{ 
#' \bold{Variable}        \tab \bold{Values} \tab \bold{Required/optional} \cr
#' Site identifier        \tab \code{integer < 1E8}. \tab required\cr
#' Time-point identifier  \tab \code{integer < 1e4}  \tab required\cr
#' Count                  \tab \code{integer <2e9} or missing code \tab required\cr
#' Weight                 \tab \code{real > 0.001}\tab optional\cr
#' Category of 1st covariate  \tab \code{integer} in \code{1,2,...,90}\tab optional\cr
#' \eqn{\vdots}           \tab\tab\cr
#' Category of last covariate \tab\tab\cr
#' } 
#' 
#' The restrictions on the values applied because of the internal representation
#' of data in TRIM. In the current implementation, all columns are read into R's
#' \code{integer} or \code{numeric} format. This means that the largest possible
#' integer is now given by \code{.Machine$integer.max} (see
#' \code{\link[base]{.Machine}}).
#' 
#' 
#' @param x a filename or a \code{\link{trimbatch}} object
#' @param missing \code{[integer]} missing value indicator
#' @param weight \code{[logical]} indicate presence of a weight column
#' @param ncovars \code{[logical]} The number of covariates in the file
#' @param labels \code{[character]} (optional) specify labels for the covariates. 
#'     Defaults to \code{cov<i>} (\code{i=1,2,...,ncovars}) if none are specified.
#' @param ... (unused)
#'
#' @return A \code{data.frame} with the following columns
#' \tabular{lll}{
#' \bold{Variable}    \tab\bold{status}   \tab \bold{R class}\cr
#' \code{site}        \tab requiered \tab \code{integer}\cr
#' \code{time}        \tab required  \tab \code{integer}\cr
#' \code{count}       \tab required  \tab \code{numeric}\cr
#' \code{weight}      \tab optional  \tab \code{numeric}\cr
#' \code{<covariate1>}\tab optional\tab \code{integer}\cr
#' \code{...}\tab\tab\cr
#' \code{<covariateN>}\tab optional\tab \code{integer}\cr
#' }
#' Missing values are translated to \code{\link{NA}}. 
#'
#'
#' @export
read_tdf <- function(x,...){
  UseMethod("read_tdf")
}

#' @rdname read_tdf
read_tdf.character <- function(x, missing = -1, weight = FALSE, ncovars=0, labels=character(0)){
  tdfread(file=x, missing=missing, weight=weight,ncovars=ncovars, labels=labels) 
}


#' @rdname read_tdf
read_tdf.trimbatch <- function(x,...){
  tc <- x[[1]]
  tdfread(tc$file, missing = tc$missing, weight = tc$weight, ncovars = tc$ncovars, labels=tc$labels)
}


# workhorse function for the S3 interfaces
tdfread <- function(file, missing, weight, ncovars, labels){
 
  if ( ncovars > 0 && length(labels) == 0 ){
    labels <- paste0("cov",seq_len(ncovars))
  } else if ( ncovars != length(labels)) {
    stop(sprintf("Length of 'labels' (%d) unequal to 'ncovars' (%d)",length(labels),ncovars))
  }
   
  colclasses <- c(site = "integer", time = "integer", count="numeric")
  if (weight) colclasses['weight'] <- "numeric"
  # add labels and names for covariates
  colclasses <- c(colclasses, setNames(rep("integer",ncovars), labels))
  
  
  # by default, one or more blanks (space, tab) are used as separators
  tab <- tryCatch(
    read.table(file, header=FALSE, colClasses=colclasses, col.names = names(colclasses))
    , error=function(e) snifreport(file, colclasses))
  if ( nrow(tab) > 0 ) tab[tab == missing] <- NA
  tab
}


snifreport <- function(file, colclasses){
  if (!file.exists(file)) stop(sprintf("Could not find file %s",file))
  ncl <- length(colclasses)
  lns <- readLines(file,n=5)
  cls <- paste(paste0(names(colclasses),"<",colclasses,">"),collapse=" ")
  msg <- sprintf("\n\rExpected %s columns: %s\nStart of file looks like this:\n",ncl,cls)
  msg <- paste0(msg,paste(sprintf("\r%s\n",lns), collapse=""))
  stop(msg, call.=FALSE)
}


