read_icv <- function(x,...){
  UseMethod("read_icv")
}

read_icv.character <- function(x, J=0, ...){
  icv_read(filenm=x, J=J)
}


read_icv.trimcommand <- function(x,...){
  basename <- strsplit(x$file,"\\.")[[1]][1]
  filenm <- paste0(basename,".ICV")
  icv_read(filenm, x$ntimes)
}


# Workhorse reader
icv_read <- function(filenm, J=0)
{
  m = unname(as.matrix(read.table(filenm)))
  if (J>0) { # perform a check
    stopifnot(ncol(m)==J) # COVIN should have J columns
    stopifnot(nrow(m)%%J == 0) #... and I*J rows
  } else J <- ncol(m) # autodetect
  nsite <- nrow(m) %/% J
  covin <- vector("list", nsite)
  idx = 1:J
  for (i in 1:nsite) {
    covin[[i]] <- m[idx, ]
    idx = idx + J
  }
  covin
}