#' Stores TRIM command options
#'
#'
#' @slot origin Where did the data for this object come from (either a filname or 'commandline').
#' @slot file    \code{[character]} name of file containing training data.
#' @slot title   \code{[character]} A string to be printed in the output file.
#' @slot ntimes  \code{[character]} Number of time points.
#' @slot ncovars \code{[character]} Number of covariates.
#' @slot labels  \code{[character]} Covariate label.
#' @slot missing \code{[integer]} Missing value indicator.
#' @slot weight  \code{[logical]} Whether a weight column is present in the \code{file}.
#' @slot comment \code{[character]} A string to be printed in the output file.
#' @slot weighting \code{[logical]} Whether weights are to be used in the model.
#' @slot serialcor \code{[logical]} Whether serial correlation is assumed in the model.
#' @slot overdist \code{[logical]} Whether overdispersion is taken into account by the model.
#' @slot basetime \code{[integer]} Position of the base time point (must be positive).
#' @slot model    \code{[integer]} What model to use (1, 2 or 3).
#' @slot covariates \code{[integer]} Number of covariates to include.
#' @slot changepoints \code{[integer]} Positions of the change points to include.
#' @slot stepwise \code{[logical]} Whether stepwise selection of the changepoints is to be used.
#' @slot outputfiles \code{[character]} Type of outputfile to generate ('F' and/or 'S')
#' @slot run \code{[logical]}, IGNORED (run the file)
#'
#' @rdname TRIMcommand
#' 
#' @seealso \code{\link{read_tcf}}
#'
TRIMcommand <- setClass(Class="TRIMcommand"
  , slots=c(
     origin       = "character"
    , file         = "character"
    , title        = "character"
    , ntimes       = "integer"
    , ncovars      = "integer"
    , labels       = "character"
    , missing      = "integer"
    , weight       = "logical"
    , comment      = "character"
    , weighting    = "logical"
    , serialcor    = "logical"
    , overdisp     = "logical"
    , basetime     = "integer"
    , model        = "integer"
    , covariates   = "integer"
    , changepoints = "integer"
    , stepwise     = "logical"
    , outputfiles  = "character"
    , run          = "logical"
    )
)


extract_tcf_key <- function(x,key,endkey=NULL,type=c('character','integer','logical')) {
  type <- match.arg(type)
  re <- ifelse(is.null(endkey)
               , paste0(key,".+?(\\n|$)") # PWB: fixed erroneous neglectance of last item
               , paste0(key,".+?END"))
  
  
  m <- regexpr(re,x,ignore.case=TRUE)         # Fetch key-value
  s <- trimws(regmatches(x,m))                # ... remove surrounding whitespace
  s = trimws(gsub(key,"",s,ignore.case=TRUE)) # Remove key and more  whitespace
  
  # remove endkey if appropriate
  if (!is.null(endkey)) s <- trimws(gsub(endkey,"",s,ignore.case=TRUE))
  
  # try splitting into tokens
  if (length(s)>0 && nchar(s)>0) {
    L <- strsplit(s, split="([[:blank:]]|\n)+")
    s <- unlist(L)
  }
  
  if (type == 'character') {
    out <- ifelse(length(s)==0, NA_character_, s)
  } else if (type == 'integer') {
    out <- ifelse(length(s)==0, NA_integer_, as.integer(s))
  } else if (type == "logical") {
    out <- ifelse(length(s)==0, NA, ifelse(length(s) == 1 && tolower(s) %in% c("on","present"),TRUE, FALSE))
  } else {
    stop("Can't happen")
  }
  
  out
}


#' Read a TRIM command file
#'
#' @section Details:
#' 
#' Read Trim Command Files, compatible with the Windows TRIM programme.
#' 
#'
#' @param file Location of tcf file.
#' 
#' @return An object of class \code{\link{TRIMcommand}}
#' @export
read_tcf <- function(file){

  # TODO: interpret changepoints properly
  #tcf <- readChar(file,file.info(file)$size)
  tcf <- paste(readLines(file), collapse="\n") # PWB: fixed CRLF issue

  TRIMcommand(
    origin = file
    , file         = convert_path(extract_tcf_key(tcf,"FILE"))
    , title        = extract_tcf_key(tcf, "TITLE")
    , ntimes       = extract_tcf_key(tcf, "NTIMES" , type="integer")
    , ncovars      = extract_tcf_key(tcf, "NCOVARS", type="integer")
    , labels       = extract_tcf_key(tcf, "LABELS", endkey="END")
    , missing      = extract_tcf_key(tcf, "MISSING", type="integer")
    , weight       = extract_tcf_key(tcf, "WEIGHT", type="logical")
    , comment      = extract_tcf_key(tcf, "COMMENT")
    , weighting    = extract_tcf_key(tcf, "WEIGHTING", type="logical")
    , serialcor    = extract_tcf_key(tcf, "SERIALCOR", type="logical")
    , overdisp     = extract_tcf_key(tcf, "OVERDISP", type="logical")
    , basetime     = extract_tcf_key(tcf, "BASETIME", type="integer")
    , model        = extract_tcf_key(tcf, "MODEL", type="integer")
    , covariates   = extract_tcf_key(tcf, "COVARIATES", type="integer")
    , changepoints = extract_tcf_key(tcf, "CHANGEPOINTS", type="integer")
    , stepwise     = extract_tcf_key(tcf,"STEPWISE", type="logical")
    , outputfiles  = extract_tcf_key(tcf,"OUTPUTFILES")
    , run          = grepl("RUN",tcf)
    )
}


setMethod("show",signature = "TRIMcommand",function(object){
  cat(sprintf("Object of class '%s' with the following parameters:\n",class(object)))
  for (n in slotNames(object)){
    cat(sprintf("%-12s: %s\n",n,as.character(paste(slot(object,n),collapse=", ") ) ))
  }
})

convert_path <- function(x){
  if (grepl("\\\\",x)){
    y <- gsub("\\\\","/",x)
    message(sprintf("Converting DOS style path '%s'\nto POSIX compliant path '%s'\n",x,y))
    y
  } else {
    x
  }
}

### some tests
#x <- read_tcf("test.tcf")



