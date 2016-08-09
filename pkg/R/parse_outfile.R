
# library(rtrim)
# data(skylark)
# set_trim_verbose(FALSE)
# m <- trim(count ~ time + site, data=skylark, model=2)
# 
# r <- readLines("pkg/tests/testthat/outfiles/skylark-1d.out")
# r <- paste(r,collapse="\n")
# cat(r)

#' Read a TRIM3 output file
#'
#' @param file \code{[character]} filename
#' 
#' @return A character string of class \code{tof}
#' 
#' @family parse_output
#' 
#' @keywords internal
read_tof <- function(file){
  r <- readLines(file)
  r <- paste(r,collapse="\n")
  class(r) <- c("tof","character")
  r
}

print.tof <- function(x,...){
  cat(x)
}

#' Extract time indices from \code{tof} object
#'
#' @param x An object of class \code{tof}
#' 
#' @return A \code{data.frame}
#' @family parse_output
#' @keywords internal
get_time_indices <- function(x){
  stopifnot(inherits(x,"tof"))
  tab <- get_indextab(x,"TIME INDICES")
  # in the first row, absence of std error is misread.
  tab[1,4] <- tab[1,3]
  tab[1,3] <- NA  
  tab
}

#' Extract time totals from \code{tof} object
#'
#' @inheritParams get_time_indices
#' 
#' @return A \code{data.frame}
#' @family parse_output
#' @keywords internal
get_time_totals <- function(x){
  stopifnot(inherits(x,"tof"))
  get_indextab(x,"TIME TOTALS")
}

get_indextab <- function(x, label){
  re <- paste0(label,".*?\n[[:blank:]]*\n")
  # find and extracte table title, headings and table.
  m <- regexpr(re,x,ignore.case=TRUE)
  s <- regmatches(x,m)
  # remove double line endings
  s <- gsub("\n\n","\n",s)
  # remove title line
  re2 <- paste0(label,".*?\n+")
  s <- gsub(re2,"",s,ignore.case=TRUE)
  # now we have a fixed-width table, readable with std. functions
  read.table(text=s, header=TRUE, strip.white=TRUE
        , fill=TRUE, colClasses=c("integer",rep("numeric",4)))
}

#' Extract estimation method from \code{tof} object
#'
#' @inheritParams get_time_indices
#'
#' @return \code{character} string
#' @family parse_output
#' @keywords internal
get_estimation_method <- function(x){
  stopifnot(inherits(x,"tof"))
  m <- regexpr("ESTIMATION METHOD.*?\n",x)
  s <- regmatches(x,m)
  s <- gsub("\n","",s)
  L <- strsplit(s,"[[:blank:]]+=[[:blank:]]+")
  L[[1]][length(L[[1]])]
}


#' Extract overall model slope \code{tof} object
#'
#' @inheritParams get_time_indices
#'
#' @return \code{character} string
#' @family parse_output
#' @keywords internal
get_overal_model_slope <- function(x){
  stopifnot(inherits(x,"tof"))
  get_slope(x,"MODEL")
}


#' Extract overall imputed slope \code{tof} object
#'
#' @inheritParams get_time_indices
#'
#' @return \code{character} string
#' @family parse_output
#' @keywords internal
get_overal_imputed_slope <- function(x){
  stopifnot(inherits(x,"tof"))
  get_slope(x,"IMPUTED")
}

get_slope <- function(x,label){
  re <- paste0("OVERALL SLOPE ",label,".*?\n[[:blank:]]*\n")
  mm <- regexpr(re,x)
  s <- regmatches(x,mm)
  s <- trimws(gsub("\n\n","\n",s)[[1]])
  L <- trimws(strsplit(s,"\n")[[1]])
  labels <- strsplit(L[2],"[[:blank:]]+")[[1]]
  values <- as.numeric(strsplit(L[3],"[[:blank:]]+")[[1]])
  setNames(values,labels)
}


get_oneliner <- function(x,label){
  re <- paste0("ESTIMATED ",label,".*?\\n")
  mm <- regexpr(re,x,ignore.case=TRUE)
  s <-  regmatches(x,mm)
  s <- trimws(gsub("\n","",s))
  L <- trimws(strsplit(s,"=")[[1]])
  as.numeric(L[2])
}

get_overdispersion <- function(x){
  stopifnot(inherits(x,"tof"))
  get_oneliner(x,"OVERDISPERSION")
}

get_serial_correlation <- function(x){
  stopifnot(inherits(x,"tof"))
  get_oneliner(x,"SERIAL CORRELATION")
}


to_vector <- function(x){
  sapply(x, function(u){
    mm <- regexpr("[^0-9]+",u,ignore.case=TRUE)
    lab <- trimws(regmatches(u,mm))
    mm <- regexpr("[0-9]+\\.?[0-9]*",u)
    num <- as.numeric(regmatches(u,mm))
    setNames(num,lab)
  },USE.NAMES=FALSE)
}

#' Extract goodness of fit from \code{tof} object
#'
#' @inheritParams get_time_indices
#'
#' @return \code{character} string
#' @family parse_output
#' @keywords internal
get_gof <- function(x){
  re <- "GOODNESS OF FIT.*?\\n[[:blank:]]*\\n"
  mm <- regexpr(re,x,ignore.case=TRUE)
  s <- regmatches(x,mm)
  s <- gsub("\n\n","\n",s)
  s <- trimws(strsplit(s,"\\n")[[1]][-1])
  L <- strsplit(s,",")
  lapply(L,to_vector)
}



