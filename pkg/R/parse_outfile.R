
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


#' Extract nr of times from \code{tof} object
#'
#' @param x An object of class \code{tof}
#' 
#' @return \code{numeric}
#' @family parse_output
#' @keywords internal
get_n_site <- function(x){
  mm <- regexpr("Site[[:blank:]]{5,}.*?\n",x)
  s <- regmatches(x,mm)
  get_num(s)
}


#' Extract nr of sites from \code{tof} object
#'
#' @param x An object of class \code{tof}
#' 
#' @return \code{numeric}
#' @family parse_output
#' @keywords internal
get_n_time <- function(x){
  mm <- regexpr("Time[[:blank:]]{5,}.*?\n",x)
  s <- regmatches(x,mm)
  get_num(s)
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

# extract numbers from text.
get_num <- function(x){
  mm <- regexpr("[0-9]+\\.?[0-9]*",x)
  as.numeric(regmatches(x,mm))
}


#' Extract goodness of fit from \code{tof} object
#'
#' @inheritParams get_time_indices
#'
#' @return List of type \code{trim.gof}
#' @family parse_output
#' @keywords internal
get_gof <- function(x){
  re <- "GOODNESS OF FIT.*?\\n[[:blank:]]*\\n"
  mm <- regexpr(re,x,ignore.case=TRUE)
  s <- regmatches(x,mm)
  s <- gsub("\n\n","\n",s)
  s <- trimws(strsplit(s,"\\n")[[1]][-1])
  L <- strsplit(s,",")
  K <- lapply(L,function(x) sapply(x, get_num, USE.NAMES=FALSE))
  names(K) <- c("chi2","LR","AIC")
  K$chi2 <- setNames(as.list(K$chi2),c("chi2","df","p"))
  K$LR <- setNames(as.list(K$LR),c("LR","df","p"))
  K$AIC <- unname(K$AIC)
  structure(K,class="trim.gof")
}

#' Extract Wald test parameters from \code{tof} object
#'
#' @inheritParams get_time_indices
#'
#' @return list of type \code{trim.wald}
#' @family parse_output
#' @keywords internal
get_wald <- function(x){
  stopifnot(inherits(x,"tof"))
  
  model <- get_model(x)
  
  # Wald test parameters.
  # TODO: implement for model=1, model=3
  L <- switch(as.character(model)
    , "2" = {
      re <- "Wald-test[[:blank:]]{5,}.*?\\n[[:blank:]]*\\n"
      mm <- regexpr(re,x,ignore.case=TRUE)   
      s <- regmatches(x,mm)
      num <- sapply(strsplit(s,",")[[1]],get_num,USE.NAMES=FALSE)
      list(model=model, W = num[1],df=num[2], p=num[3])
  })
  structure(L,class="trim.wald")
}


#' Extract model type from \code{tof} object
#'
#' @inheritParams get_time_indices
#'
#' @return Model number, either 1, 2 or 3.
#' @family parse_output
#' @keywords internal
get_model <- function(x){
  stopifnot(inherits(x,"tof"))
  if (grep("WALD-TEST FOR SIGNIFICANCE OF SLOPE",x)){ 
    2L
  } else if(grep("WALD-TEST FOR SIGNIFICANCE OF CHANGES IN SLOPE",x)) {
    3L
  } else {
    1L
  }
}

#' Extract model coefficients from \code{tof} object
#'
#' @inheritParams get_time_indices
#'
#' @return list of class \code{trim.coef}
#' @family parse_output
#' @keywords internal
get_coef <- function(x){
  stopifnot(inherits(x,"tof"))
  model <- get_model(x)
  
  # find parameter estimates block
  re <- "PARAMETER ESTIMATES[[:blank:]]*\\n[[:blank:]]*\\n.*?\\n[[:blank:]]*\\n"
  mm <- regexpr(re,x)
  s <- regmatches(x,mm)
  # remove first two lines
  s <- sub(".*?\\n\\n","",s)
  coef <- read.table(text=s, header=TRUE)
  names(coef)[4] <- "std.err"
  structure(list(model=model, coef=coef),class="trim.coef")
}



