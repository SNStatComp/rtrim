# Some basic assertions to test whether models can be run.

#' Check whether there are sufficient observations to run a model
#'
#'
#'
#'
#' @param x A \code{\link{trimcommand}} object, a \code{data.frame}, or the location of a TRIM command file.
#' @param ... Parameters passed to other methods.
#' 
#' @family modelspec
#' 
#' @export
check_observations <- function(x,...){ 
  UseMethod("check_observations")
}

#' @param model \code{[numeric]} Model 1, 2 or 3?
#' @param covar \code{[character|numeric]} column index of covariates in \code{x}
#' @param time.id \code{[character|numeric]} column index of time points in \code{x}
#' @param count.id \code{[character|numeric]} column index of the counts in \code{x}
#' @param changepoints \code{[numeric]} Changepoints (model 2 only)
#' @param eps \code{[numeric]} Numbers whose absolute magnitude are lesser than \code{eps} are considered zero.
#' 
#' @return A \code{list} with two components. The component \code{sufficient} takes the value
#' \code{TRUE} or \code{FALSE} depending on whether sufficient counts have been found.
#' The component \code{errors} is a \code{list}, of which the structure depends on the chosen model,
#' that indicates under what conditions insufficient data is present to estimate the model.
#' 
#' \itemize{
#' \item{For model 3 without covariates, \code{$errors} is a list whose single element is a vector of time
#' points with insufficient counts}.
#' \item{For model 3 with covariates, \code{$errors} is a named list with an element for each covariate
#' for which insufficients counts are encountered. Each element is a two-column \code{data.frame}. The
#' first column indicates the time point, the second column indicates for what covariate value insufficient
#' counts are found.}
#' \item{Model 2: TODO}
#' }
#' 
#' 
#' 
#' @export
#' @rdname check_observations
check_observations.data.frame <- function(x, model, covar = list()
  , changepoints = numeric(0), time.id="time",count.id="count", eps=1e-8, ...){
 
  stopifnot(model %in% 1:3)
  
  out <- list()
  
  if (model==3 && length(covar) == 0 ){
    time_totals <- tapply(X = x[,count.id], INDEX = x[,time.id], FUN = sum, na.rm=TRUE)
    ii <- time_totals <= eps
    out$sufficient <- !any(ii)
    out$errors <- setNames(list(x[ii,time.id]),time.id)
  } else if (model == 3 && length(covar>0)) {
    out$errors <- get_cov_count_errlist(x[,count.id],x[,time.id],covars=x[covar],timename=time.id)
    out$sufficient <- length(out$errors) == 0
  } else {
    warning("Data checking is currently implemented for mode 3 only. Returning empty list")
  }
  
  out
}


#' @export
#' @rdname check_observations
check_observations.trimcommand <- function(x,...){
  dat <- read_tdf(x$file)
  check_observations.data.frame(x=dat,model=x$model, covar=x$labels[x$covariates]
                                , changepoints = x$changepoints)
}

#' @export
#' @rdname check_observations
check_observarions.character <- function(x,...){
  tc <- read_tcf(x)
  check_observations.trimcommand(tc,...)
}


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
    df <- as.data.frame(as.table(tab))
    err <- df[df$Freq==0, 1:2]
    if (nrow(err) > 0){
      names(err)[2] <- covname
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

