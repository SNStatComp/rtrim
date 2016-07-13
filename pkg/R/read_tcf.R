#' Create a TRIMCommand object
#'
#' Define a TRIM calculation through a TRIMCommand object. 
#' 
#' @section Description:
#' 
#' If no parameters are passed, a default TRIMCommand is returned. All parameters listed here
#' are optional.
#' 
#' @param ... Options in the form of \code{key=value}. See below for all options.
#' 
#' @section Options:
#' 
#' \itemize{
#' \item{ \code{file}    \code{[character]} name of file containing training data.}
#' \item{ \code{title}   \code{[character]} A string to be printed in the output file.}
#' \item{ \code{ntimes}  \code{[character]} Number of time points.}
#' \item{ \code{ncovars} \code{[character]} Number of covariates.}
#' \item{ \code{labels}  \code{[character]} Covariate label.}
#' \item{ \code{missing} \code{[integer]} Missing value indicator.}
#' \item{ \code{weight}  \code{[logical]} Whether a weight column is present in the \code{file}.}
#' \item{ \code{comment} \code{[character]} A string to be printed in the output file.}
#' \item{ \code{weighting} \code{[logical]} Whether weights are to be used in the model.}
#' \item{ \code{serialcor} \code{[logical]} Whether serial correlation is assumed in the model.}
#' \item{ \code{overdist} \code{[logical]} Whether overdispersion is taken into account by the model.}
#' \item{ \code{basetime} \code{[integer]} Position of the base time point (must be positive).}
#' \item{ \code{model}    \code{[integer]} What model to use (1, 2 or 3).}
#' \item{ \code{covariates} \code{[integer]} Number of covariates to include.}
#' \item{ \code{changepoints} \code{[integer]} Positions of the change points to include.}
#' \item{ \code{stepwise} \code{[logical]} Whether stepwise selection of the changepoints is to be used.}
#' \item{ \code{outputfiles} \code{[character]} Type of outputfile to generate ('F' and/or 'S')}
#'}
#'
#' @export
new_TRIMCommand <- function(...){
  # decide on default values
  tc <- list(
    file           = character(0)
    , title        = character(0)
    , ntimes       = integer(0)
    , ncovars      = integer(0)
    , labels       = character(0)
    , missing      = integer(0)
    , weight       = character(0)
    , comment      = character(0)
    , weighting    = character(0)
    , serialcor    = character(0)
    , overdisp     = character(0)
    , basetime     = integer(0)
    , model        = integer(0)
    , covariates   = integer(0)
    , changepoints = integer(0)
    , stepwise     = "off"
    , outputfiles  = character(0)
  )
  class(tc) <- c("TRIMCommand","list")
  L <- list(...)
  for ( nm in names(L) ){
    if (! nm %in% names(tc) ) stop(sprintf("'%s' is not a valid TRIM keyword",nm))
    if (nm == "file") L[[nm]] <- convert_path(L[[nm]])
    # convert and set
    if (length(L[[nm]])>0) tc[[nm]] <- as(L[[nm]], class(tc[[nm]]))
  }
  tc
}

#' Read a TRIM command file
#'
#' Read TRIM Command Files, compatible with the Windows TRIM programme.
#'
#' @section TRIM Command file format:
#' 
#' TRIM command files are text files that specify a TRIM job, where a job
#' consists of one or more models to be computed on a single data input file.
#' TRIM command files are commonly stored with the extension \code{.tcf}, but
#' this is not a strict requirement.
#'
#' A TRIM command file is built up out of a sequence of commands and optionally 
#' comments. A command can be single-line or multi-line. In the case of a 
#' single-line command, the line starts with a keyword, followed by one or more 
#' spaces, followed by one or more option values. If there are multiple option 
#' values, these are separated one or more spaces. A multi-line command starts 
#' with a keyword, followed by one option value on each consecutive line. The 
#' end of a multiline command is indicated with the keyword \code{END} on a new
#' line. Currently, the only multi-line command is \code{LABELS}.
#' 
#' The keyword \code{RUN} (at the beginning of a single line) ends the
#' specification of a single model. After this a new model can be specified.
#' Parameters not specified in the next model will be copied from the previous
#' one.
#' 
#' @section Encoding issues:
#'
#' To read files containing non-ASCII characters encoded in a format that is not
#' native to your system, specifiy the \code{encoding} option. This causes R to 
#' re-encode to native encoding upon reading. Input encodings supported for your
#' system can be listed by calling \code{\link[base]{iconvlist}()}. For more 
#' information on Encoding in R, see \code{\link[base]{Encoding}}.
#' 
#' @section Note on filenames:
#' 
#' If the \code{file} is specified using backslashes to separate directories
#' (Windows style), this will be converted to a filename using forward slashes
#' (POSIX style, as used by R).
#' 
#' 
#'
#' @param file Location of tcf file.
#' @param encoding The encoding in which the file is stored.
#'
#' @return An object of class \code{TRIMCommand}
#' 
#' @seealso \code{\link{new_TRIMCommand}}
#' @export
read_tcf <- function(file, encoding=getOption("encoding")){
  con <- file(description = file, encoding=encoding)
  tcf <- paste(readLines(con), collapse="\n") 
  close(con)
  
  tcflist <- trimws(strsplit(tcf,"(\\n|^)RUN")[[1]])
  L <- vector(mode="list",length=length(tcflist))
  L[[1]] <- tc_from_char(tcflist[[1]])
  for ( i in 1+seq_along(L[-1]) ){
    L[[i]] <- tc_from_char(tcflist[[i]], default = L[[i-1]])
  }
  class(L) <- "TRIMCommandList"
  L
}


print.TRIMCommandList <- function(x,...){
  y <- x[[1]]
  cat(sprintf("TRIMCommandList: %s\n",y$title))
  cat(sprintf("file: %s (%d means missing)\n"
              , y$file, y$missing))
  cat(sprintf("Weights %s, %d covariates labeled %s\n",y$weight, y$ncovars
              ,paste0("",paste(y$labels,collapse=", "))))

  cat("\nModel parameter overview:\n") 
  models <- data.frame(
    Nr = seq_along(x)
    , comment = sapply(x,`[[`,"comment")
    , weighting = sapply(x,`[[`,"weighting")
    , serialcor = sapply(x, `[[`,"serialcor")
    , overdisp  = sapply(x, `[[`,"overdisp")
    , basetime  = sapply(x, `[[`,"basetime")
    , model     = sapply(x, `[[`, "model")
    , covariates = sapply(x, function(m) paste(m$covariates,collapse=", ") )
    , changepoints = sapply(x, function(m) paste(m$changepoints, collapse=", "))
    , stepwise = sapply(x, `[[`, "stepwise")
    , outputfiles = sapply(x, `[[`,"outputfiles")
    , stringsAsFactors=FALSE
  )
  print(models)
}


convert_path <- function(x){
  if (isTRUE(grepl("\\\\",x)) ){
    y <- gsub("\\\\","/",x)
    y
  } else {
    x
  }
}


shortfilename <- function(x){
  if ( identical(x, character(0)) || nchar(x) <= 20 ) return(x)
  st <- substr(x,1,3)
  n <- nchar(x)
  en <- substr(x,n-13,n)
  paste0(st,"...",en)
}

print.TRIMCommand <- function(x,pretty=FALSE,...){
  if (!pretty){
    cat("Object of class TRIMcommand:\n")
    for ( nm in names(x) ){
      cat(sprintf("%12s: %s\n",nm,paste0("",paste(x[[nm]]),collapse=", ")) )
    }
  } else {
    cat(sprintf("TRIMcommand: %s (%s)\n",x$title,x$comment))
    cat(sprintf("  %d time points from %s (%d represents missing)\n",x$ntimes, shortfilename(x$file),x$missing))
    cat(sprintf("  Model %d with serial correlation %s and overdispersion %s\n",x$model,x$serialcor, x$overdisp))
    cat(sprintf("  Weights are %s and turned %s for the model\n", tolower(x$weight), x$weighting))
    cat(sprintf("  Data contains %d covariates labeled %s\n", x$ncovars, paste0("",paste(x$labels,collapse=", ")) ))
    cat(sprintf("  Covatiates used in the model include %s\n",paste(x$covariates,collapse=", ")))
    if (length(x$changepoints)==0){ 
      cat(sprintf("  No changepoints defined."))
    } else {
      cat(sprintf("  Changepoints defined at %s.",paste0(x$changepoints,collapse=", ")))
    }
    cat(" Stepwise turned %s\n",paste0("",x$stepwise))
  }
}

summary.TRIMcommand <- function(x,...){
  print(x,pretty=TRUE)
}

key_regex <- function(trimkey){
  re <- paste0("(\\n|^)", trimkey,".+?")
  re <- if (trimkey == "LABELS"){
    paste0(re,"END")
  } else {
    paste0(re,"(\\n|$)")
  }
  re
}

extract_keyval <- function(trimkey, x){
  trimkey <- toupper(trimkey)
  re <- key_regex(trimkey)
  m <- regexpr(re,x,ignore.case=TRUE)      # Fetch key location
  s <- regmatches(x,m)                     # extract substring
  re <- paste0("(",trimkey,")|(END)")      # remove key and whitespace
  s <- trimws(gsub(re,"",s))               # remove key and return
  # split values when relevant
  if (trimkey != "COMMENT" && length(s)>0 && nchar(s)>0) {
    s <- unlist(strsplit(s, split="([[:blank:]]|\n)+"))
  }
  s
}


tc_from_char <- function(x, default = new_TRIMCommand()){
  L <- lapply(names(default), extract_keyval, x)
  L <- setNames(L, names(default))
  for ( i in seq_along(L))
    if (length(L[[i]])==0) L[[i]] <- default[[i]]
  do.call(new_TRIMCommand, L)
}


setNames <- function (object = nm, nm) {
  names(object) <- nm
  object
}

