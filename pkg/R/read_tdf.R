
#' Read TRIM data file
#'
#' @param file \code{[character]} Input data file. See details section for spec.
#' @param missing_code \code{[integer]} Code for missing counts (see Details).
#' @param nsnif \code{[integer]} Number of lines read to determine input format.
#'
#'
#' @section Details:
#' 
#' TRIM input data is specified in the original TRIM manual as follows.
#' 
#' The data file is an ASCII file containing one line (a record) for each 
#' combination of Site and Time. So, for \eqn{I} sites and \eqn{J} time- points, the 
#' number of records is \eqn{I\times J}. Each record contains the following variables
#' (the order is important!), separated by one or more spaces.
#' 
#' 
#' \tabular{lll}{
#' Variable               \tab Values \tab Required/optional \cr
#' Site identifier        \tab \code{integer < 1E8}. \tab required\cr
#' Time-point identifier  \tab \code{integer < 1e4}  \tab required\cr
#' Count                  \tab \code{integer <2e9} or missing code \tab required\cr
#' Weight                 \tab \code{real > 0.001}\tab optional\cr
#' Category of 1st covariate          \tab \code{integer} in \code{1,2,...,90}\tab optional\cr
#' \eqn{\vdots}           \tab\tab\cr
#' Category of last covariate         \tab\tab\cr
#' }
#' The missing code (see section 3.2) must be a integer in the range \code{(-32767...32767)}
#' and should be chosen outside the range of observed counts. Zero will usually not be
#' outside the range, but a negative number such as -1 will always be outside the range
#' of observed counts.
#'
#'
#' @export
read_tdf <- function(file, missing_code=-1L, nsnif=10L){
  # snif the file structure
  lines <- readLines(con=file, n=nsnif, warn=FALSE)
  L <- strsplit(x=lines,split=" +")
  len <- sapply(L,length)
  ncol <- unique(len)
  if (length(ncol) !=1 ) {
    stop(sprintf("Detected different numbers of columns in first %d rows",length(L)))
  }
  if (ncol==0){
    warning("This file contains no records")
    return(NULL)
  } else if( ncol < 4){
    stop(sprintf("A TRIM data input file must contain at least 4 columns, found %s",len))
  }
  
  col_classes <- c("integer","integer","integer","numeric")
  ncat <- len-4
  col_classes <- c(col_classes, rep("integer",ncol-4))
  col_names <- c("site","time","count","weight")
  col_names <- c(col_names,sprintf("covar%02d",seq_len(ncol-4))) 
  dat <- read.table(file=file
    , header=FALSE
    , sep=""
    , colClasses=col_classes
    , col.names = col_names
  )
  
  within(dat, count[count==missing_code] <- NA)
  
}
