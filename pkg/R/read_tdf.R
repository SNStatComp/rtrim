
#' Read TRIM data file
#'
#' @param file \code{[character]} Input data file. See details section for spec.
#' @param missing_code \code{[integer]} Code for missing counts (see Details).
#' @param snif \code{[integer]} Number of lines read to determine input format.
#' @param weight \code{[logical]} Is there a weight column present?
#' @param strict \code{[logical]} Check data against TRIM requirements? (see Details).
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
#' If \code{strict==TRUE}, the data are checked against the demands in the table.
#' Otherwise they are just treated as R-native types with less restrictions
#' on the size of integers.
#'
#' @return A \code{data.frame} with columns named \code{site}, \code{time},
#' \code{count}, and optionally \code{weight} and \code{covar01...covarNN}
#' for covariate labels.
#'
#' @export
#read_tdf <- function(file, missing_code=-1L, snif=10L, weight=FALSE, strict=FALSE) {
read_tdf <- function(tcf=NULL,
                     file="nofile.dat", missing_code=-1L, weight=FALSE, covars=character(0),
                     snif=10L, strict=FALSE, dbg=TRUE)
{
  # Extended argument parsing
  if (!is.null((tcf))) {

    if (missing(file)) {
      if (dbg) cat(sprintf("Using TCF file = \"%s\"\n", tcf@file))
      file <- tcf@file
    }
    
    if (missing(missing_code)) {
      if (dbg) cat(sprintf("Using TCF missing = %d\n", tcf@missing))
      missing_code <- tcf@missing
    }
    
    if (missing(weight)) {
      if (dbg) cat(sprintf("Using TCF weight = %s\n", tcf@weight))
      weight <- tcf@weight
    }
    
    if (missing(covars)) {
      covars <- tcf@labels
      if (dbg) cat(sprintf("Using TCF covars = [%s]\n", paste(covars,collapse=", ")))
    }
    
  } else if (dbg) {
    
    if (missing(file)) sprintf("Using default: file=\"%s\"\n", file)
    else               sprintf("Using specified: file=\"%s\"\n", file)
    
    if (missing(file)) sprintf("Using default: missing=%d\n", missing_code)
    else               sprintf("Using specified: missing=%d\n", missing_code)

    if (missing(file)) sprintf("Using default: weight=%s\n", weight)
    else               sprintf("Using specified: weight=%s\n", weight)
    
    if (missing(covars)) sprintf("Using default: covars=[%s]\n", paste(covars,collapse=", "))
    else               sprintf("Using specified: covars=[%s]\n", paste(covars,collapse=", "))
  } else stop("Can't happen")
  
  # handle file origin
  if (!is.null((tcf))) file <- file.path(dirname(tcf@origin), file) 
  
  # snif the file structure
  lines <- readLines(con=file, n=snif, warn=FALSE)
  lines <- trimws(lines) # remove leading/trailing whitespace which fcks up splitting
  L <- strsplit(x=lines,split=" +")
  len <- sapply(L,length)
  ncol <- unique(len)
  if (length(ncol)==0 || ncol==0){
    warning("This file contains no records")
    return(NULL)
  }
  
  # Check correct number of columns
  if (length(ncol) !=1 ) stop(sprintf("Detected different numbers of columns in first %d rows",length(L)))
  mincol <- 3 # site, time, count
  if (weight) mincol <- mincol+1
  ncovars <- length(covars)
  mincol <- mincol + ncovars
  if (ncol < mincol) stop(sprintf("TRIM data input file must contain at least %d columns, found %d", mincol, ncol))

  # required columns
  col_classes <- c("integer", "integer", "numeric")
  col_names   <- c("site",    "time",    "count")
  
  # optional column: weight
  if (weight) {
    col_classes <- c(col_classes, "numeric")
    col_names   <- c(col_names,   "weight")
  }
  
  # Specified covariate columns
  col_classes <- c(col_classes, rep("integer",ncovars))
  col_names   <- c(col_names, covars)
  
  # Unspecified covariate columns
  if (ncol > mincol) {
    nextra = ncol - mincol
    col_classes <- c(col_classes, rep("integer", nextra))
    col_names   <- c(col_names,   sprintf("covar%02d",seq_len(nextra)+ncovars))
    ncovars <- ncovars + nextra
  }
  
  # last check
  stopifnot(length(col_classes)==ncol)
  stopifnot(length(col_names)==ncol)
  

  dat <- read.table(file=file
    , header=FALSE
    , sep=""
    , colClasses=col_classes
    , col.names = col_names
  )
  
  dat <- within(dat, count[count==missing_code] <- NA)
  if (strict) check_tdf(dat)
  
  # Convert site, time and covariates to factors
  dat$site = factor(dat$site)
  dat$time = ordered(dat$time)
  covar_cols = ifelse(weight, 5,4) : ncol # beware: ifelse(weight, 5:ncol, 4:ncol) does not work
  for (col in covar_cols) {
    dat[ ,col] = factor(dat[ ,col])
  }

  if (dbg) {
    print(head(dat, n=3))
    cat(sprintf("...\n"))
    print(tail(dat, n=3))
  }
  
  dat
}

# check against the TRIM requirements
check_tdf <- function(x, weight){
  stopifnot(all(x$site)<1e8)
  stopifnot(all(x$time<1e4))
  stopifnot(all(x$count<2e9), all(x$count>=0))
  if (weight) stopifnot(all(x$weight > 0.001))
  nvar <- ifelse(weight,4,3)
  ncov <- ncol(x) - nvar
  icov <- seq_len(ncov) + nvar
  stopifnot(all(x[I] <=90), all(x[I]>0) )
  x
}


