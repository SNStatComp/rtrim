# Some basic assertions to test whether models can be run.

print_and_capture <- function(x){
  paste(capture.output(print(x)),collapse="\n")
}


assert_positive <- function(x, varname){
  if (any(x <= 0)){
    i <- which(x<=0)
    stop(sprintf("Found zero or less counts for %s %s",varname, paste(names(x[i]),collapse=", ")),call.=FALSE)
  }
  invisible(TRUE)
}

assert_increasing <- function(x, varname){
  if( !all(diff(x)>0) ){
    stop(sprintf(
      "%s not ordered or containing duplicates",varname
    ), call.=FALSE)
  }
  invisible(TRUE)
}


# sufficient data per index (index=time for model 3)
assert_sufficient_counts <- function(count, index){
  time_totals <- tapply(X = count, INDEX = index, FUN = sum, na.rm=TRUE)  
  assert_positive(time_totals, names(index))
}

# sufficient data for piecewise linear trend model
assert_plt_model <- function(count, time, changepoints){
  assert_increasing(changepoints, "change points")
  if (length(changepoints)==0) changepoints <- 1
  # label the pieces in piecewise linear regression
  pieces <- integer(length(time)) + 1
  C <- changepoints
  if (C[length(C)] != max(time)) C <- append(C,max(time))
  
  for ( i in seq_along(C[-1])){
    j <- seq(C[i] + 1, C[i+1])
    pieces[time %in% j] <- C[i]
  }
  assert_sufficient_counts(count, list(changepoint=pieces))
}

# count: vector of counts
# time: vector of time points
# covar: list of covariate vectors
assert_covariate_counts <- function(count, time, covars){
  ERR <- list()
  for ( i in seq_along(covars) ){
    covname <- names(covars)[i]
    cov <- covars[[i]]
    index <- list(time=time,value = cov)
    tab <- tapply(count, INDEX=index, FUN=sum, na.rm=TRUE)
    err <- which(tab <= 0,arr.ind=TRUE)
    dimnames(err) <- setNames(list(NULL,colnames(err)),c("",covname))
    if (length(err) > 0){
      ERR[[covname]] <- err
    }
  }
  if (length(ERR)>0) 
    stop("Zero observations for the following cases:\n",gsub("\\$.*?\n","",print_and_capture(ERR)),call.=FALSE)
  invisible(TRUE)
}



