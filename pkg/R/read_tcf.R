
new_TRIMcommand <- function(...){
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
  )
  class(tc) <- c("TRIMCommand","list")
  L <- list(...)
  for ( nm in names(L) ){
    if (! nm %in% names(tc) ) stop(sprintf("'%s' is not a valid TRIM keyword",nm))
    if (nm == "file") L[[nm]] <- convert_path(L[[nm]])
    # convert and set
    tc[[nm]] <- as(L[[nm]], class(tc[[nm]]))
  }
  tc
}


convert_path <- function(x){
  if (grepl("\\\\",x)){
    y <- gsub("\\\\","/",x)
#    message(sprintf("Converting DOS style path '%s'\nto POSIX compliant path '%s'\n",x,y))
    y
  } else {
    x
  }
}


shortfilename <- function(x){
  if (nchar(x)<=20) return(x)
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


tc_from_char <- function(x){
  tc <- new_TRIMcommand()
  L <- lapply(names(tc), extract_keyval, x)
  L <- setNames(L, names(tc))
  do.call(new_TRIMcommand, L)
}

read_tcf <- function(file, encoding=getOption("encoding")){
  con <- file(description = file, encoding=encoding)
  tcf <- paste(readLines(con), collapse="\n") 
  close(con)
  
  tcflist <- trimws(strsplit(tcf,"(\\n|^)RUN"))
  L <- lapply(tcflist, tc_from_char)
  class(L) <- "TRIMCommandList"
  L
}

print.TRIMCommandList <- function(x,pretty=TRUE,...){
  cat(sprintf("TRIMCommandList with %d models:\n",length(x)) )
  i <- 0
  for ( tc in x ){
    i <- i+1
    cat(sprintf("%02d: ",i))
    print(tc, pretty=pretty)
  }
}

summary.TRIMCommandList <- function(x,...){
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
    , stringsAsFactors=FALSE
  )
  print(models)
}


