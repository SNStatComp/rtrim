
#' Read TRIM data files
#'
#' Read data files intended for the original TRIM programme.
#'
#' @section The TRIM data file format:
#'
#' TRIM input data is stored in a \code{ASCII} encoded file where headerless columns
#' are separated by one or more spaces. Below are the columns as \code{read_tdf} expects
#' them.
#'
#' \tabular{lll}{
#' \bold{Variable}    \tab\bold{status}   \tab \bold{R type}\cr
#' \code{site}        \tab requiered \tab \code{integer}\cr
#' \code{time}        \tab required  \tab \code{integer}\cr
#' \code{count}       \tab required  \tab \code{numeric}\cr
#' \code{weight}      \tab optional  \tab \code{numeric}\cr
#' \code{<covariate1>}\tab optional\tab \code{integer}\cr
#' \code{...}\tab\tab\cr
#' \code{<covariateN>}\tab optional\tab \code{integer}\cr
#' }
#'
#'
#' @param x a filename or a \code{\link{trimcommand}} object
#' @param missing \code{[integer]} missing value indicator.
#'      Missing values are translated to \code{\link[base]{NA}}.
#' @param weight \code{[logical]} indicate presence of a weight column
#' @param ncovars \code{[logical]} The number of covariates in the file
#' @param labels \code{[character]} (optional) specify labels for the covariates.
#'     Defaults to \code{cov<i>} (\code{i=1,2,...,ncovars}) if none are specified.
#' @param ... (unused)
#'
#' @return A \code{data.frame}.
#'
#' @family modelspec
#' @export
read_tdf <- function(x,...){
  UseMethod("read_tdf")
}

#' @rdname read_tdf
#' @export
read_tdf.character <- function(x, missing = -1, weight = FALSE, ncovars=0, labels=character(0),...){
  tdfread(file=x, missing=missing, weight=weight,ncovars=ncovars, labels=labels)
}


#' @rdname read_tdf
#' @export
read_tdf.trimcommand <- function(x,...){
  tdfread(x$file, missing = x$missing, weight = x$weight, ncovars = x$ncovars, labels=x$labels)
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
  if (nrow(tab)==0) stop(sprintf("file \"%s\" appears to be empty", file))
  if (nrow(tab) > 0) tab[tab == missing] <- NA
  tab
}


snifreport <- function(file, colclasses){
  if (!file.exists(file)) stop(sprintf("Could not find file %s",file))
  ncl <- length(colclasses)
  lns <- readLines(file,n=5)
  if (length(lns)==0) stop(sprintf("file \"%s\" appears to be empty", file))
  cls <- paste(paste0(names(colclasses),"<",colclasses,">"),collapse=" ")
  msg <- sprintf("\n\rExpected %s columns: %s\nStart of file looks like this:\n",ncl,cls)
  msg <- paste0(msg,paste(sprintf("\r%s\n",lns), collapse=""))
  stop(msg, call.=FALSE)
}


