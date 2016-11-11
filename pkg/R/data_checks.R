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
#' @param covars \code{[character|numeric]} column index of covariates in \code{x}
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
#' first column indicates the time point, the second column indicates for which covariate value insufficient
#' counts are found.}
#' \item{For Model 2, without covariates \code{$errors} is a list with a single
#' element \code{changepoints}. It points out what changepoints lead to a time
#' slice with zero observations.}
#' \item{For Model 2, with covariates \code{$errors} is a named list with an
#' element for each covariate for which inssufficients counts are encountered.
#' Each element is a two-column \code{data.frame}, The first colum indicates the
#' changepoint, the second column indicates for which covariate value
#' insufficient counts are found.}
#' }
#'
#'
#'
#' @export
#' @rdname check_observations
check_observations.data.frame <- function(x, model, covars = list()
  , changepoints = numeric(0), time.id="time",count.id="count", eps=1e-8, ...){

  stopifnot(model %in% 1:3)

  out <- list()
  if (model==3 && length(covars) == 0 ){
    time_totals <- tapply(X = x[,count.id], INDEX = x[,time.id], FUN = sum, na.rm=TRUE)
    ii <- time_totals <= eps
    out$sufficient <- !any(ii)
    out$errors <- setNames(list(x[ii,time.id]),time.id)
  } else if (model == 3 && length(covars>0)) {
    out$errors <- get_cov_count_errlist(x[,count.id],x[,time.id],covars=x[covars],timename=time.id)
    out$sufficient <- length(out$errors) == 0
  } else if ( model == 2 ) {
    pieces <- pieces_from_changepoints(x[,time.id],changepoints)
    ok <- pieces > 0 # allow zero counts for changepoint 0
    if ( length(covars) == 0){
      time_totals <- tapply(X=x[ok,count.id],INDEX=pieces[ok], FUN = sum, na.rm=TRUE)
      ii <- time_totals <= eps
      out$sufficient <- !any(ii)
      out$errors <- list(changepoint = as.numeric(names(time_totals))[ii])
    } else {
      out$errors <- get_cov_count_errlist(x[ok,count.id], pieces[ok], x[ok,covars,drop=FALSE], timename="changepoint")
      out$sufficient <- length(out$errors) == 0
    }
  }

  out
}


#' @export
#' @rdname check_observations
check_observations.trimcommand <- function(x,...){
  dat <- read_tdf(x$file)
  check_observations.data.frame(x=dat,model=x$model, covars=x$labels[x$covariates]
                                , changepoints = x$changepoints)
}

#' @export
#' @rdname check_observations
check_observations.character <- function(x,...){
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
assert_sufficient_counts <- function(count, index) {
  time_totals <- tapply(X=count, INDEX=index, FUN=sum, na.rm=TRUE)
  assert_positive(time_totals, names(index))
}

# Get an indicator for the pieces in 'piecewise linear model'
# that are encoded in changepoints.
pieces_from_changepoints <- function(time, changepoints) {
  if (length(changepoints)==0) return(rep(0,length(time))) # nothing to do

  # convert time from (possibly non-contiguous) years to time points 1..ntime
  tpt <- as.integer(ordered(time))

  pieces <- integer(length(tpt))
  C <- changepoints
  if (C[length(C)] != max(tpt)) C <- append(C,max(tpt))

  for ( i in seq_along(C[-1])){
    j <- seq(C[i] + 1, C[i+1])
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
  ok = pieces>0 # Allow zero observations for changepoint 0
  if (length(covars)==0){
    assert_sufficient_counts(count[ok], list(changepoint=pieces[ok]))
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

    # df$Freq[is.na(df$Freq)] <- 0 # replace NA -> 0
    #
    # # Allow no-data at cp 0
    # CP0 = levels(df$time)[1]
    # err <- df[df$Freq==0 & df$time!=CP0, 1:2]
    err <- df[df$Freq==0, 1:2]
    row.names(err) <- NULL
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
get_deletion <- function(count, time, changepoints, covars) {
  # if ( changepoints[1] != 1) changepoints <- c(1,changepoints)
  out <- -1
  if (length(changepoints)==1) return(out) # Never propose to delete a lonely changepoint
  pieces <- pieces_from_changepoints(time=time, changepoints=changepoints)

  if ( length(covars)> 0){
    err <- get_cov_count_errlist(count, pieces, covars,timename="piece")
    if ( length(err)>0){
      e <- err[[1]]
      out <- as.numeric(as.character(e[1,1]))
    }
  } else {
    tab <- tapply(count, list(pieces=pieces), sum,na.rm=TRUE)
    j <- tab <= 0
    if (any(j)){
      wj = which(j)[1]
      out <- as.numeric(names(tab)[min(wj+1, length(tab))])
    }
  }
  out
}

autodelete <- function(count, time, changepoints, covars) {
  out <- get_deletion(count, time, changepoints, covars)
  while (out > 0) {
    cp = as.integer(out)
    yr = time[cp]
    if (cp==yr) rprintf("Auto-deleting change point %d\n", cp)
    else        rprintf("Auto-deleting change point %d (%d)\n", cp, yr)
    # delete changepoint
    changepoints <- changepoints[changepoints != out]
    out <- get_deletion(count, time, changepoints, covars)
  }
  changepoints
}

