# Some basic assertions to test whether models can be run.

# capture how an object is printed in a string.
print_and_capture <- function(x){
  paste(capture.output(print(x)),collapse="\n")
}

# all x positive or an error
assert_positive <- function(x, varname){
  if (any(x <= 0)){
    i <- which(x<=0)
    stop(sprintf("Found zero or less counts for %s %s",varname, paste(names(x[i]),collapse=", ")),call.=FALSE)
  }
  invisible(TRUE)
}

# x strictly increasting or an error
assert_increasing <- function(x, varname){
  if( !all(diff(x)>0) ){
    stop(sprintf(
      "%s not ordered or containing duplicates",varname
    ), call.=FALSE)
  }
  invisible(TRUE)
}


# sufficient data per index (index=time for model 3, pieces for model 2)
assert_sufficient_counts <- function(count, index){
  time_totals <- tapply(X = count, INDEX = index, FUN = sum, na.rm=TRUE)
  assert_positive(time_totals, names(index))
}

# Get an indicator for the pieces in 'piecewise linear model'
# that are encoded in changepoints.
pieces_from_changepoints <- function(time, changepoints) {
  # convert time from (possibly non-contiguous) years to time points 1..ntime
  tpt <- as.integer(ordered(time))

  pieces <- integer(length(tpt)) + 1
  C <- changepoints
  if (C[length(C)] != max(tpt)) C <- append(C,max(tpt))

  for ( i in seq_along(C[-1])){
    # j <- seq(C[i] + 1, C[i+1])
    j <- seq(C[i], C[i+1] - 1)
    pieces[tpt %in% j] <- C[i]
  }
  pieces
}

# sufficient data for piecewise linear trend model
assert_plt_model <- function(count, time, changepoints, covars){
  assert_increasing(changepoints, "change points")
  if (length(changepoints)==0) changepoints <- 1
  # label the pieces in piecewise linear regression
  pieces <- pieces_from_changepoints(time, changepoints)
  if (length(covars)==0){
    assert_sufficient_counts(count, list(changepoint=pieces))
  } else {
    assert_covariate_counts(count, pieces, covars, timename="changepoint")
  }
}


# get a list of errors: for which time (pieces) and covariate values
# are there zero counts? Result is an empty list or a named list
# of matrices with columns 'timename', value (of the covariate)
get_cov_count_errlist <- function(count, time, covars, timename="time"){
  ERR <- list()

  for ( i in seq_along(covars) ){
    covname <- names(covars)[i]
    cov <- covars[[i]]
    index <- list(time=time,value = cov)
    names(index)[1] <- timename
    tab <- tapply(count, INDEX=index, FUN=sum, na.rm=TRUE)
    err <- which(tab <= 0,arr.ind=TRUE)
    dimnames(err) <- setNames(list(NULL,colnames(err)),c("",covname))
    if (length(err) > 0){
      err[,2] <- colnames(tab)[err[,2]]
      ERR[[covname]] <- err
    }
  }

  ERR
}

# count: vector of counts
# time: vector of time points
# covar: list of covariate vectors
assert_covariate_counts <- function(count, time, covars, timename="time"){
  err <- get_cov_count_errlist(count, time, covars, timename=timename)
  if ( length(err)>0 )
    stop("Zero observations for the following cases:\n"
         , gsub("\\$.*?\n","",print_and_capture(err))
         , call.=FALSE)
  invisible(TRUE)
}


# Return the first changepoint to delete (if any).
# returns the value of the CP, or -1 when nothing
# needs to be deleted.
get_deletion <- function(count, time, changepoints, covars){
  if ( changepoints[1] != 1) changepoints <- c(1,changepoints)
  pieces <- pieces_from_changepoints(time=time, changepoints=changepoints)
  out <- -1
  if ( length(covars)> 0){
    err <- get_cov_count_errlist(count, pieces, covars)
    if ( length(err)>0){
      e <- err[[1]]
      out <- changepoints[as.numeric(e[1,1])]
    }
  } else {
    tab <- tapply(count, list(pieces=pieces), sum,na.rm=TRUE)
    j <- tab <= 0
    if (any(j)){
      out <- as.numeric(names(tab)[which(j)[1]])
    }
  }
  out
}

autodelete <- function(count, time, changepoints, covars){

  out <- get_deletion(count, time, changepoints, covars)
  while ( out > 0){
    rprintf("Auto-deleting change point %d\n",as.integer(out))
    # delete changepoint
    changepoints <- changepoints[changepoints != out]
    out <- get_deletion(count, time, changepoints, covars)
  }
  changepoints
}

