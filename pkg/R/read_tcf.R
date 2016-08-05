
#' Create a trimcommand object
#'
#'
#' @section Description:
#' 
#' A \code{trimcommand} object stores a single TRIM model, including the
#' specification of the data file. Normally, such an object is defined by
#' reading a legacy TRIM command file.
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
#' @seealso \code{\link{read_tcf}}
#' 
#' @export
trimcommand <- function(...){
  # decide on default values
  tc <- list(
    file           = character(0)
    , title        = character(0)
    , ntimes       = integer(0)
    , ncovars      = integer(0)
    , labels       = character(0)
    , missing      = integer(0)
    , weight       = logical(0)
    , comment      = character(0)
    , weighting    = logical(0)
    , serialcor    = logical(0)
    , overdisp     = logical(0)
    , basetime     = integer(0)
    , model        = integer(0)
    , covariates   = integer(0)
    , changepoints = 1L
    , stepwise     = logical(0)
    , outputfiles  = character(0)
  )
  class(tc) <- c("trimcommand","list")
  L <- list(...)
  for ( nm in names(L) ){
    if (! nm %in% names(tc) ) stop(sprintf("'%s' is not a valid TRIM keyword",nm))
    if (nm == "file") L[[nm]] <- convert_path(L[[nm]])
    # convert and set
    if (length(L[[nm]])>0) tc[[nm]] <- as_rtrim(L[[nm]], tc[[nm]])
  }
  tc
}


trimbatch <- function(...){
  tc <- list(trimcommand(...))
  class(tc) <- c("trimbatch", "list")
  tc
}


# Add a model to a trimbatch object
#
# Set up multiple models in a trimbatch object.
# 
#
# @param x a \code{trimbatch} object.
# @param ... model parameters (see \code{\link{trimbatch}}). Unspecified parameters
#    are copied from the last model in the list.
# @export
#add_model <- function(x,...){
#  UseMethod("add_model")
#}

# @export 
# @rdname add_model


# convert from character representation of tcf to rtrim internal representation.
as_rtrim <- function(value, template){
  if ( inherits(template, "logical") ){
    if ( tolower(value) %in% c("present","on") ) TRUE else FALSE
  } else {
    as(value,class(template))
  }
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
#' A TRIM command file consists of two parts. The first part describes the
#' data file to be read, the second part describes the model(s) to be run. A
#' TRIM command file can only contain a single data specification part, but multiple
#' models may be specified.
#' 
#' Each command starts on a new line with a keyword, followed by at least
#' one space and at least one option value, where multiple option values are
#' separated by spaces. All commands must be written on a single line, except
#' the \code{LABELS} command (to set labels for covariates). The latter command
#' starts with \code{LABELS} on a single line, followed by a newline, followed
#' by a new label on each following line. The keyword \code{END} (at the beginning of a line) 
#' signals the end of the labels command.
#' 
#' The keyword \code{RUN} (at the beginning of a single line) ends the
#' specification of a single model. After this a new model can be specified.
#' Parameters not specified in the current model will be copied from the previous
#' one.
#' 
#' @section TRIM commands:
#' 
#' The commands are identical to those in the original TRIM software. Commands
#' that represent a simple toggle (on/off, present/absent) are translated to a
#' \code{logical} upon reading. Below we give commands in upper case, but the
#' commands are parsed case insensitively.
#' 
#' \tabular{ll}{
#' \bold{Data}\tab\cr
#' \code{FILE}   \tab data filename and path.\cr
#' \code{TITLE}  \tab A title (appears in output when exported).\cr
#' \code{NTIMES} \tab [positive integer] Number of time points in data file.\cr
#' \code{NCOVARS}\tab [nonnegative integer] Number of covariates in data file.\cr
#' \code{LABELS} \tab Covariate labels (multiline command). \cr
#' \code{END}    \tab Signals end of \code{LABELS} command.\cr
#' \code{MISSING}\tab missing value indicator.\cr
#' \code{WEIGHT} \tab [\code{present}, \code{absent}] Indicates whether weights are present in the data file [translated to \code{logical}].\cr
#' \bold{Model}  \tab\cr
#' \code{COMMENT}  \tab A comment for the current model.\cr
#' \code{WEIGHTING}\tab [\code{on},\code{off}] Switch use of weights for current model [translated to \code{logical}].\cr
#' \code{SERIALCOR}\tab [\code{on},\code{off}] Switch use of serial correlation for current model [translated to \code{logical}].\cr
#' \code{OVERDISP}\tab [\code{on},\code{off}] Switch use of overdispersion for current model [translated to \code{logical}].\cr
#' \code{BASETIME}\tab [integer] Index of base time-point.\cr
#' \code{MODEL}\tab [\code{1}, \code{2}, \code{3}] Choose the current model\cr
#' \code{COVARIATES}\tab [integers] indices of covariates to use (1st covariate has index 1)\cr
#' \code{CHANGEPOINTS} \tab [integers] indices of changepoints\cr
#' \code{RUN}\tab Signals end of current model specification.
#' }
#' 
#' 
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
#' @param file Location of TRIM command file.
#' @param encoding The encoding in which the file is stored.
#' @param simplify Return a single \code{trimcommand} object if only one 
#'   model is specified in the TRIM command file.
#' 
#' 
#' @return A trimcommand object, or in the case of multiple models in a single
#' TRIM command file, a \code{list} of \code{trimcommand} objects. In the
#' latter case, a useful summary can be printed with \code{\link{summary.trimbatch}}.
#' 
#' @seealso \code{\link{read_tcf}}, \code{\link{trimcommand}}
#' @export
read_tcf <- function(file, encoding=getOption("encoding"),simplify=TRUE){
  con <- file(description = file, encoding=encoding)
  tcf <- paste(readLines(con), collapse="\n") 
  close(con)
  
  tcflist <- trimws(strsplit(tcf,"(\\n|^)RUN")[[1]])
  L <- vector(mode="list",length=length(tcflist))
  L[[1]] <- tc_from_char(tcflist[[1]])
  for ( i in 1+seq_along(L[-1]) ){
    L[[i]] <- tc_from_char(tcflist[[i]], default = L[[i-1]])
  }
  class(L) <- c("trimbatch","list")
  if (length(L) == 1 && simplify ) L[[1]] else L 
}

#' summarize a trimbatch object
#'
#' @export
#' @keywords internal
#' @param x an object
#' @param ... options (ignored)
summary.trimbatch <- function(object,...){
  y <- object[[1]]
  cat(sprintf("trimbatch: %s\n",pr(y$title, len=Inf)))
  cat(sprintf("file: %s (%s means missing)\n"
              , pr(y$file, len=50), pr(y$missing)))
  cat(sprintf("Weights %s, %s covariates labeled %s\n",pr(y$weight), pr(y$ncovars)
              ,paste0("",paste(y$labels,collapse=", "))))

  cat("\nModel parameter overview:\n") 
  oneliner(object)
}


convert_path <- function(x){
  if (isTRUE(grepl("\\\\",x)) ){
    y <- gsub("\\\\","/",x)
    y
  } else {
    x
  }
}




#' print a trimcommand object
#'
#' @export
#' @keywords internal
#' @param x an R object
#' @param ... optional parameters (ignored)
print.trimcommand <- function(x,...){
  cat("Object of class trimcommand:\n")
  for ( nm in names(x) ){
    cat(sprintf("%12s: %s\n",nm,paste0("",paste(x[[nm]]),collapse=", ")) )
  }
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
  m <- regexpr(re,x,ignore.case=TRUE)           # Fetch key location
  s <- regmatches(x,m)                          # extract substring
  re <- paste0("(",trimkey,")|(END)")           # remove key and whitespace
  s <- trimws(gsub(re,"",s,ignore.case = TRUE)) # remove key and return
  # split values when relevant
  if (trimkey != "COMMENT" && length(s)>0 && nchar(s)>0) {
    s <- unlist(strsplit(s, split="([[:blank:]]|\n)+"))
  }
  s
}


tc_from_char <- function(x, default = trimcommand()){
  L <- lapply(names(default), extract_keyval, x)
  L <- setNames(L, names(default))
  for ( i in seq_along(L))
    if (length(L[[i]])==0) L[[i]] <- default[[i]]
  do.call(trimcommand, L)
}


setNames <- function (object = nm, nm) {
  names(object) <- nm
  object
}


# Utilities for pretty printing

shortstr<- function(x,len=12){
  if ( identical(x, character(0)) || nchar(x) <= len ) return(x)
  st <- substr(x,1,4)
  n <- nchar(x)
  en <- substr(x,n-5,n)
  paste0(st,"..",en)
}

pr <- function(s, ...){
  a <- if ( length(s) == 0 ) "<none>" else paste(as.character(s),collapse=", ")
  shortstr(a, ...)
}

oneliner <- function(x){
  cat(sprintf("%10s %9s %9s %8s %8s %6s %6s %12s %8s %8s\n"
              ,"comment","weighting","serialcor","overdisp"
              ,"basetime","model","covars","changepoints","stepwise","outfiles"))
  for ( i in seq_along(x)){
    tc <- x[[i]]
    cat(sprintf("%10s %9s %9s %8s %8s %6s %6s %12s %8s %8s\n"
      , pr(tc$comment)
      , pr(tc$weighting)
      , pr(tc$serialcor)
      , pr(tc$overdisp)
      , pr(tc$basetime)
      , pr(tc$model)
      , pr(tc$covariates)
      , pr(tc$changepoints)
      , pr(tc$stepwise)
      , pr(tc$outputfiles)
      ))
  }
  
}


