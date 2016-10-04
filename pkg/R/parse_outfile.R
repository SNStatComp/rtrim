
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

#' Read a TRIM3 variance-covariance output file
#'
#' @param file \\code{[character]} filename
#'
#' @return A matrix of class \code{[numeric]}
#'
#' @family parse_output
#'
#' @keywords internal
read_vcv <- function(file) {
  df <- read.table(file) # Read as data.frame
  m  <- as.matrix(df)    # Convert to matrix...
  out <- unname(m)       # ... and remove dimname attributes to prevent test_equal problem
}

#' Extract TRIM version used for output
#'
#' @return \code{character}
#' @family parse_output
#' @keywords internal
get_version <- function(x){
  re <- "TRIM (\\d\\.\\d+) :  TRend analysis.*"
  m <- regexec(re, x)
  a <- m[[1]][2]
  b <- attr(m[[1]],"match.length")[2]
  version <- substr(x, a, a+b-1)
  version
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
  tab[1,c(3,5)] <- 0
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
  v <- get_version(x)
  # str(x)
  # str(v)
  if (v=="3.61") re <- paste0("OVERALL SLOPE ",label,".*?\n[[:blank:]]*(\n|$)")
  else           re <- paste0("OVERALL SLOPE ",label,".*intercept.*?\n[[:blank:]]*(\n|$)")
  mm <- regexpr(re,x)
  s <- regmatches(x,mm)
  s <- trimws(gsub("\n\n","\n",s)[[1]])
  L <- trimws(strsplit(s,"\n")[[1]])
  labels <- strsplit(L[2],"[[:blank:]]+")[[1]]
  values <- as.numeric(strsplit(L[3],"[[:blank:]]+")[[1]])
  # str(mm)
  # str(s)
  # str(labels)
  # str(values)
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
  # Read to second block if necessary
  m <- regexpr("STEPWISE SELECTION OF CHANGEPOINTS.*",x)
  if (m[1] != -1) x <- regmatches(x, m)
  # OK, proceed with target
  get_oneliner(x,"OVERDISPERSION")
}

get_serial_correlation <- function(x){
  stopifnot(inherits(x,"tof"))
  # Read to second block if necessary
  m <- regexpr("STEPWISE SELECTION OF CHANGEPOINTS.*",x)
  if (m[1] != -1) x <- regmatches(x, m)
  # OK, proceed with target
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
  # Read to second block if necessary
  m <- regexpr("STEPWISE SELECTION OF CHANGEPOINTS.*",x)
  if (m[1] != -1) x <- regmatches(x, m)
  # OK, proceed with target
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
  # Read to second block if necessary
  m <- regexpr("STEPWISE SELECTION OF CHANGEPOINTS.*",x)
  if (m[1] != -1) x <- regmatches(x, m)

  re <- "WALD-TEST FOR SIGNIFICANCE OF DEVIATIONS FROM LINEAR TREND.*?\\n\\n"
  s <- get_str(re,x)
  deviations <- if(length(s)>0){
    s <- strip_line(s)
    u <- as.list(get_num(strsplit(s,",")[[1]]))
    setNames(u,c("W","df","p"))
  } else {NULL}

  re <- "WALD-TEST FOR SIGNIFICANCE OF SLOPE PARAMETER.*?\\n\\n"
  s <- get_str(re,x)
  slope <- if (length(s)>0){
    s <- strip_line(s)
    u <- as.list(get_num(strsplit(s,",")[[1]]))
    setNames(u,c("W","df","p"))
  } else {NULL}

  re <- "WALD-TEST FOR SIGNIFICANCE OF CHANGES IN SLOPE.*?\\n\\n"
  s <- get_str(re,x)
  dslope <- if (length(s)>0){
    s <- strip_line(s)
    u <- read.table(text=s,header=TRUE)
    list(W=u[,2],df=u[1,3],p=u[,4])
  } else {NULL}

  re <- "WALD-TEST FOR SIGNIFICANCE OF COVARIATES.*?\\n\\n"
  s <- get_str(re,x)
  covar <- if (length(s)>0){
    s <- strip_line(s)
    u <- read.table(text=s,header=TRUE)
    setNames(u,c("Covariate","W","df","p"))
  } else {NULL}

  structure(
    list(deviations=deviations, slope=slope, dslope=dslope, covar=covar )
    ,class="trim.wald")

}


get_str<-function(re,x,...){
  mm <- regexpr(re,x,...)
  regmatches(x,mm)
}

strip_line<-function(x,n=1){
  x <- sub(".*?\\n","",x)
  if (n>1) strip_line(x,n=n-1) else x
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
  if (grepl("RESULTS FOR MODEL: Linear Trend",x,ignore.case=TRUE)){
    2L
  } else if(grepl("RESULTS FOR MODEL: Effects for each time point",x,ignore.case=TRUE)) {
    3L
  } else {
    1L
  }
}

#' Extract model coefficients from \code{tof} object
#'
#' @inheritParams get_time_indices
#' @param covars \code{[character]} names of covariates.
#' @return list of class \code{trim.coef}
#' @family parse_output
#' @keywords internal
get_coef <- function(x,covars){
  stopifnot(inherits(x,"tof"))
  model <- get_model(x)
  hascov <- has_covariates(x)
  hascp <- has_changepoints(x)

  if (model==2 & !hascp & !hascov ){ # e.g. skylark-1d
    re <- "PARAMETER ESTIMATES[[:blank:]]*\\n[[:blank:]]*\\n.*?\\n[[:blank:]]*\\n"
    s <- get_str(re,x)
    s <- strip_line(s,2)
    s <- replace_with_space("Slope",s)
    out <- read.table(text=s,header=TRUE)
    names(out) <- c("add","se_add","mul","se_mul")
    fromto <- range(get_time_indices(x)$Time)
    out$from <- fromto[1]
    out$upto <- fromto[2]
    out <- out[c("from","upto","add","se_add","mul","se_mul")]
  }

  if (model==2 & hascp & !hascov){ # e.g. skylark-1e, 1f
    re <- "PARAMETER ESTIMATES[[:blank:]]*\\n[[:blank:]]*\\n.*?\\n[[:blank:]]*\\n"
    s <- get_str(re,x)
    s <- strip_line(s,3)
    out <- read.table(text=s,header=TRUE)
    names(out)[3:6] <- c("add","se_add","mul","se_mul")
  }

  if (model==2 & hascp & hascov){ # e.g. skylark-2a
    re <- "PARAMETER ESTIMATES.*?\\n\\n\\n"
    s <- get_str(re,x)
    s <- strip_line(s,3)
    s <- sub("\\n\\n\\n","",s)
    u <- strsplit(s,"\\n\\n")[[1]]
    out <- do.call(rbind,lapply(u,coef_m2_covars))
    out$from <- shd(out$from)
    out$upto <- shd(out$upto)
    out <- out[with(out,order(covar,cat)),,drop=FALSE]
    out$covar <- sub("Constant","baseline",out$covar)
    i <- grepl("Covariate",out$covar)
    covnum <- get_num(out$covar[i])
    out$covar[i] <- covars[covnum]
    out
  }



  if (model == 3 & !hascov){ # e.g. skylark-1a,b,c
    re <- "Parameters for each time point[[:blank:]]*\\n[[:blank:]]*\\n.*?\\n[[:blank:]]*\\n"
    s <-get_str(re,x,ignore.case=TRUE)
    # remove first two lines
    s <- sub(".*?\\n\\n","",s)
    out <- read.table(text=s, header=TRUE, strip.white=TRUE, fill=TRUE)
    out[1,c(3,5)] <- 0
    out[1,4] <- 1
    names(out) <- c("time","add","se_add","mul","se_mul")
  }

  if (model==3 & hascov){ # e.g. skylark-2b
    re <- "PARAMETER ESTIMATES.*?\\n\\n\\n"
    s <- get_str(re,x)
    s <- strip_line(s,3)
    s <- sub("\\n\\n\\n","",s)
    u <- strsplit(s,"\\n\\n")[[1]]
    out <- do.call(rbind,lapply(u,coef_m3_covars))
    out$covar <- sub("Constant","baseline",out$covar)
    i <- grepl("Covariate",out$covar)
    covnum <- get_num(out$covar[i])
    out$covar[i] <- covars[covnum]
    out
  }
  out

}


# get_covariates <- function(x){
#   lines <- strsplit(x,"\n")[[1]]
#   lines <- lines[grepl("number of values",lines)]
#   lines <- lines[!grepl("(Site)|(Time)|(Count)",lines)]
#   lines <- trimws(sub("number of values.*","",lines))
#   trimws(sub("[0-9]+\\.","",lines))
# }

has_covariates <- function(x){
  grepl("Covariate",x)
}


has_changepoints <- function(x){
  grepl("from +upto",x)
}

replace_with_space <- function(pattern,x){
  str <- get_str(pattern,x)
  replace <- paste0(rep(" ", nchar(str)),collapse="")
  sub(str,replace,x)
}

## Helper functions for parsing coefficients/model 3
# sequential hotdeck imputor
shd <- function(x){
  if (length(x)==1) return(x)
  val <- x[1]
  for ( i in 2:length(x)){
    if (is.na(x[i])){
      x[i] <- val
    } else {
      val <- x[i]
    }
  }
  x
}

# model 2, with covariates
coef_m2_covars <- function(x){
  if (grepl("from upto",x)){
    x <- strip_line(x)
    y <- get_str("[0-9].*?\\n",x)
    fromto <- sapply(trimws(strsplit(y," +")[[1]]), get_num, USE.NAMES=FALSE)
    x <- strip_line(x,2)
    beta_i <- "Constant"
    x <- sub("Constant","        ",x)
    cat <- 0
  } else { # Covariate
    fromto <- c(NA,NA)
    beta_i <- trimws(get_str(".*?\\n",x))
    x <- strip_line(x,2)
    st <- get_str("Category.*?[0-9]+?",x)
    cat <- get_num(st)
    replace <- paste(rep(" ",nchar(st)),collapse="")
    x <- sub(st,replace,x)
  }
   out <- read.table(text=x)
   names(out) <- c("add","se_add","mul","se_mul")
   out$from <- fromto[1]
   out$upto <- fromto[2]
   out$covar <- beta_i
   out$cat <- cat
   out[c("covar","cat","from","upto","add","se_add","mul","se_mul")]
}

# model 3 with covariates
coef_m3_covars <- function(x){
  beta_i <- trimws(get_str(".*?\\n",x))
  x <- strip_line(x,2)

  if (grepl("Category",x)){
    cat <- get_num(get_str("Category.*?\n",x))
    x <- strip_line(x,1)
  } else {
    cat <- 0
  }

  x <- gsub("Time","    ",x)
  out <- read.table(text=x,header=FALSE)
  names(out) <- c("time","add","se_add","mul","se_mul")
  out$cat <- cat
  out$covar <- beta_i
  out[c("covar","cat","time","add","se_add","mul","se_mul")]
}



