
.empty_sites <- function(count, site) {
  # convert site to factor; remember the original value
  # (allowed: integer; numeric; character; factor)
  if (inherits(site, "integer")) {
    site <- factor(site)
  }

}



# Some basic assertions to test whether models can be run.

#' Check whether there are sufficient observations to run a model
#'
#' @param x A \code{\link{trimcommand}} object, a \code{data.frame}, or the location of a TRIM command file.
#' @param ... Parameters passed to other methods.
#'
#' @family modelspec
#'
#' @export
check_observations <- function(x, ...){
  UseMethod("check_observations")
}

#' @param model \code{[numeric]} Model 1, 2 or 3?
#' @param count_col \code{[character|numeric]} column index of the counts in \code{x}
#' @param year_col \code{[character|numeric]} column index of years or time points in \code{x}
#' @param month_col \code{[character|numeric]} optional column index of months in \code{x}
#' @param covars \code{[character|numeric]} column index of covariates in \code{x}
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
check_observations.data.frame <- function(x, model, count_col="count", year_col="year", month_col=NULL,
                                          covars = character(0), changepoints=numeric(0), eps=1e-8, ...) {

  if (!isTRUE(model %in% 1:3)) stop("model must be 1, 2, or 3")

  if (!(count_col %in% names(x))) stop(sprintf("Column %s not found in data.frame", count_col))
  if (!(year_col %in% names(x))) stop(sprintf("Column %s not found in data.frame", year_col))
  if (!is.null(month_col)) {
    if (!(month_col %in% names(x))) stop(sprintf("Column %s not found in data.frame", month_col))
  }
  for (cv in covars) {
    if (!(cv %in% names(x))) stop(sprintf("Column %s not found in data.frame", cv))
  }

  out <- list()
  if (model==3 && length(covars) == 0) {
    # model 3, simple mode: annual or annual + monthly
    if (is.null(month_col)) { # annual only
      yt <- tapply(X=x[,count_col], INDEX = x[,year_col], FUN=sum, na.rm=TRUE) # total counts per year
      ii <- yt < eps
      out$sufficient <- !any(ii)
      if (!out$sufficient) {
        out$errors <- list()
        out$errors[[year_col]] <- names(ii)[ii]
      }
    } else { # annual + monthly
      yt <- tapply(X=x[,count_col], INDEX = x[,year_col], FUN=sum, na.rm=TRUE) # total counts per year
      mt <- tapply(X=x[,count_col], INDEX = x[,month_col], FUN=sum, na.rm=TRUE) # per month
      ii <- yt < eps
      jj <- mt < eps
      out$sufficient <- !any(ii) & !any(jj)
      if (!out$sufficient) {
        out$errors <- list()
        if (any(ii)) out$errors[[year_col]] <- names(ii)[ii]
        if (any(jj)) out$errors[[month_col]] <- names(jj)[jj]
      }
    }
  } else if (model == 3 && length(covars>0)) {
    out$errors <- get_cov_count_errlist(x[,count_col],x[,year_col],covars=x[covars],timename=year_col)
    out$sufficient <- length(out$errors) == 0
  } else if ( model == 2 ) {
    pieces <- pieces_from_changepoints(x[,year_col],changepoints)
    ok <- pieces > 0 # allow zero counts for changepoint 0
    if ( length(covars) == 0){
      time_totals <- tapply(X=x[ok,count_col],INDEX=pieces[ok], FUN = sum, na.rm=TRUE)
      ii <- time_totals <= eps
      out$sufficient <- !any(ii)
      out$errors <- list(changepoint = as.numeric(names(time_totals))[ii])
    } else {
      out$errors <- get_cov_count_errlist(x[ok,count_col], pieces[ok], x[ok,covars,drop=FALSE], timename="changepoint")
      out$sufficient <- length(out$errors) == 0
    }
  }

  out
}


#' @export
#' @rdname check_observations
check_observations.trimcommand <- function(x, ...){
  dat <- read_tdf(x$file)
  check_observations.data.frame(x=dat,model=x$model, covars=x$labels[x$covariates]
                                , changepoints = x$changepoints, ...)
}

#' @export
#' @rdname check_observations
check_observations.character <- function(x, ...){
  tc <- read_tcf(x)
  check_observations.trimcommand(tc, ...)
}


# capture how an object is printed in a string.
print_and_capture <- function(x){
  paste(capture.output(print(x)),collapse="\n")
}

# all x positive or an error
assert_positive <- function(x, varname) {
  if (any(x <= 0)){
    i <- which(x<=0)
    msg <- if (is.null(varname)) sprintf("Found zero or less counts for %s", paste(names(x[i]),collapse=", "))
           else                  sprintf("Found zero or less counts for %s %s",varname, paste(names(x[i]),collapse=", "))
    stop(msg,call.=FALSE)
  }
  invisible(TRUE)
}


# sufficient data per index (index=time for model 3, pieces for model 2)
assert_sufficient_counts <- function(count, index) {
  time_totals <- tapply(X=count, INDEX=index, FUN=sum, na.rm=TRUE)
  assert_positive(time_totals, names(index))
}

#' Link the time vector (assuming years) to the pieces (as in 'piecewise linear model) that are encoded in the changepoints.
#'
#' @param year integer vector of years (or time points 1..J)
#' @param changepoints integer vector of suggested changepoints in the range 1..J-1
#'
#' @returns a vector of the changepoints (as in 1..J-1) for each year.
pieces_from_changepoints <- function(year, changepoints, dbg=F) {
  if (dbg) {
    cat("pieces_from_changepoints()\n")
  }
  # convert actual time from (possibly non-contiguous) years to time points 1..J
  jj <- as.integer(ordered(year))
  J <- max(jj)

  # Changepoints must be converted to 1..J-1 if not already so
  if (length(changepoints)==0) {
    # case 0: no changepoints
    cpts <- 0L
  } else if (min(changepoints)>=1 && max(changepoints)<J) {
    # case 1: already in 1..J-1
    cpts <- changepoints
    if (cpts[1] > 0) cpts <- c(0L, cpts) # add prefix 0
  } else if (all(changepoints %in% year)) {
    # case 2: actual years (used); convert to 0, 1..J-1
    str(year)
    unique_years <- sort(unique(year))
    cpts <- match(changepoints, unique_years)
    cpts <- c(0L, cpts) # add prefix 0
  } else {
    stop("Invalid changepoints specified")
  }

  # Assign each time point to the corresponding change point
  pieces <- integer(length(year)) # Allocate memory;
  # N.B.: we actually use the default value of 0L for cases
  # of no changepoints, or before the first changepoint
  for (cpt in cpts) {
    idx <- which(jj >= cpt)
    pieces[idx] <- cpt
  }

  # Ready
  if (dbg) rprintf("  returning pieces: %s\n", vfmt(pieces, 20))
  return(pieces)
}

if (F) {
  pfc <- pieces_from_changepoints

  rprintf("--- Test 1\n")
  year <- c(1,2,3)
  cpts <- c(1,2)
  pieces_from_changepoints(year, cpts, dbg=T)
  rprintf("  Should be:        1, 2, 2\n")

  rprintf("--- Test 2\n")
  year <- c(1,2,3,4)
  cpts <- c(1,3)
  pieces_from_changepoints(year, cpts, dbg=T)
  rprintf("  Should be:        1, 1, 3, 3\n")

  rprintf("--- Test 3\n")
  year <- c(1,2,3,4)
  cpts <- c(1,3)
  pieces_from_changepoints(year, cpts, dbg=T)
  rprintf("  Should be:        1, 1, 3, 3\n")

  stop("intended")
}


old_pieces_from_changepoints <- function(year, changepoints, dbg=F) {
  if (dbg) cat("pieces_from_changepoints()\n")
  # convert actual time from (possibly non-contiguous) years to time points 1..J
  jj <- as.integer(ordered(year))
  J <- max(jj)

  # Changepoints must be converted to 1..J-1 if not already so
  if (length(changepoints)==0) {
    # case 0: no changepoints
    cpts <- integer(0)
  } else if (min(changepoints)>=0 && max(changepoints)<J) {
    # case 1: already in 0..J-1 or 1..J-1
    cpts <- changepoints
  } else if (all(changepoints %in% year)) {
    # case 2: actual years (used); convert to 1..J-1
    str(year)
    unique_years <- sort(unique(year))
    cpts <- match(changepoints, unique_years)
  } else {
    stop("Invalid changepoints specified")
  }

  # Assign each time point to the corresponding change point
  pieces <- integer(length(year)) # Allocate memory;
  # N.B.: we actually use the default value of 0L for cases
  # of no changepoints, or before the first changepoint
  for (cpt in cpts) {
    idx <- which(jj > cpt)
    pieces[idx] <- cpt
  }

  # Ready
  if (dbg) rprintf("  returning pieces: %s\n", vfmt(pieces, 20))
  return(pieces)
}


## Check model 2

# sufficient data for piecewise linear trend model
assert_plt_model <- function(count, time, changepoints, covars){

  # First check if the changepoints, are strictly increasing
  if (!all(diff(changepoints)>0)) {
    msg <- "changepoints not ordered, or containing duplicates"
    stop(msg, call.=FALSE)
  }

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
  for (i in seq_along(covars)) { # For all covariates
    covname <- names(covars)[i]
    cov <- covars[[i]]
    index <- list(time=time, value=cov)
    names(index)[1] <- timename
    tab <- tapply(count, INDEX=index, FUN=sum, na.rm=TRUE)
    df <- as.data.frame(as.table(tab))
    # df[,1] = as.integer(df[,1]) # time or piece chareacter->integer
    # allow no-pos-data on time pt 0
    if (timename=="changepoint") {
      idx <- df$changepoint != "0"
      df <- df[idx, ]
    }

    df$Freq[is.na(df$Freq)] <- 0 # replace NA -> 0
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
# time: vector of time point or piece ID
# covar: list of covariate vectors
assert_covariate_counts <- function(count, time, covars, timename="time"){
  err <- get_cov_count_errlist(count, time, covars, timename=timename)
  if ( length(err)>0 )
    stop("Zero observations for the following cases:\n"
         , gsub("\\$.*?\n","",print_and_capture(err))
         , call.=FALSE)
  invisible(TRUE)
}


#' Return the index of the first changepoint to delete
#'
#' @param count integer vector of counts (NA, 0 or positive)
#' @param year integer vector of years (or time points 1..J)
#' @param changepoints integer vector of suggested changepoints in the range 1..J-1
#' @param covars optional covariance info.
#' @param dbg debug flag
#'
#' @returns the index of the first changepoint to delete; or 0 if there are none to delete
get_deletion <- function(count, year, changepoints, covars=NULL, dbg=dbg) {
  if (dbg) rprintf("\nget_deletion()\n")

  .count_pos <- function(x) {
    # count number of positive elements in vector x.
    ispos <- x > 0
    sum(ispos, na.rm=TRUE)
  }

  # if ( changepoints[1] != 1) changepoints <- c(1,changepoints)
  out <- 0L
  if (length(changepoints)==1) return(out) # Never propose to delete a lonely changepoint
  #todo: this must always be chanepoint '0'
  # link the time data to the corresponding change points
  pieces <- pieces_from_changepoints(year, changepoints, dbg=dbg)

  if (length(covars) > 0){
    err <- get_cov_count_errlist(count, pieces, covars,timename="piece")
    if (length(err) > 0){
      # extract for the first covariant ([[1]]),
      # the first column, representing th piece (second [[1]]).
      # These are chanepoints as factor, so with as.integer() we get their position.
      # however, a '0' changepoints was added earlier, so we have to extract it to find the correct index
      out <- as.integer(err[[1]][[1]][1])-1L
      # e <- err[[1]]
      # cat("e:"); str(e); str(e[1,1]); str(as.integer(e[1,1]))
      # out <- as.numeric(as.character(e[1,1]))
    }
  } else { # no covars to deal with
    # count the number of positive counts per piece
    tab <- tapply(count, list(pieces=pieces), .count_pos) # was: using sum, na.rm=TRUE
    if (dbg) {
      rprintf("Tabulating pieces:\n")
      print(tab)
      rprintf("---\n")
    }
    # look for this first piece without positive counts
    j <- tab <= 0
    if (any(j)){
      # return the index of that piece/changepoint
      idx = which(j)[1]
      out <- unname(idx)
    }
  }
  if (dbg) rprintf("  Returning: %d\n", out)
  return(out)
}

if (F) {
  rprintf("--- Test 1\n")
  count <- c(1,2,1,2)
  year  <- c(1,2,3,4)
  cpts  <- c(1,2,3)
  out <- get_deletion(count, year, cpts, dbg=T)
  #expect_equal(out, 0)
  rprintf("  Should be: 0\n")

  rprintf("--- Test 2\n")
  count <- c(1,0,2,1)
  year  <- c(1,2,3,4)
  cpts  <- c(1,2,3)
  get_deletion(count, year, cpts, dbg=T)
  rprintf("  Should be: 2\n")

  rprintf("--- Test 3\n")
  count <- c(1,2,1,0)
  year  <- c(1,2,3,4)
  cpts  <- c(1,2,3)
  get_deletion(count, year, cpts, dbg=T)
  rprintf("  Should be: 2\n")

  # rprintf("--- Test 2\n")
  # year <- c(1,2,3,4)
  # cpts <- c(1,3)
  # pieces_from_changepoints(year, cpts, dbg=T)
  # rprintf("  Should be:        1, 1, 3, 3\n")
  #
  # rprintf("--- Test 3\n")
  # year <- c(1,2,3,4)
  # cpts <- c(1,3)
  # pieces_from_changepoints(year, cpts, dbg=T)
  # rprintf("  Should be:        1, 1, 3, 3\n")

  stop("intended")
}

old_get_deletion <- function(count, time, changepoints, covars, dbg=dbg) {
  if (dbg) rprintf("get_deletion()\n")

  .count_pos <- function(x) {
    # count number of positive elements in vector x.
    ispos <- x > 0
    sum(ispos, na.rm=TRUE)
  }

  # if ( changepoints[1] != 1) changepoints <- c(1,changepoints)
  out <- 0L
  if (length(changepoints)==1) return(out) # Never propose to delete a lonely changepoint
  # link the time data to the corresponding change points
  pieces <- pieces_from_changepoints(time, changepoints, dbg=dbg)

  if (length(covars) > 0){
    err <- get_cov_count_errlist(count, pieces, covars,timename="piece")
    if (length(err) > 0){
      # extract for the first covariant ([[1]]),
      # the first column, representing th piece (second [[1]]).
      # These are chanepoints as factor, so with as.integer() we get their position.
      # however, a '0' changepoints was added earlier, so we have to extract it to find the correct index
      out <- as.integer(err[[1]][[1]][1])-1L
      # e <- err[[1]]
      # cat("e:"); str(e); str(e[1,1]); str(as.integer(e[1,1]))
      # out <- as.numeric(as.character(e[1,1]))
    }
  } else { # no covars to deal with
    # count the number of positive counts per piece
    tab <- tapply(count, list(pieces=pieces), .count_pos) # was: using sum, na.rm=TRUE
    if (dbg) {
      rprintf("Tabulating pieces:\n")
      print(tab)
      rprintf("---\n")
    }
    # look for this first piece without positive counts
    j <- tab <= 0
    if (any(j)){
      # return the index of that piece/changepoint
      idx = which(j)[1]
      out <- unname(idx)
    }
  }
  if (dbg) rprintf("  Returning: %d\n", out)
  return(out)
}

#' Autodelete function
#'
#' @param count integer vector of counts (NA, 0 or positive)
#' @param time integer vector of years (or time points 1..J)
#' @param changepoints integer vector of suggested changepoints in the range 1..J-1
#' @param covars optional covariance info.
#' @param dbg debug flag
#'
#' @returns a vector with new changepoints
#' @export
#'
#' @examples
org_autodelete <- function(count, year, changepoints, covars=NULL, dbg=T) {
  if (dbg) {
    cat("autodelete()\n")
    cat("  count:"); str(count)
    cat("  year: ");  str(year)
    cat("  cpts: ");  str(changepoints)
  }

  # # convert actual time from (possibly non-contiguous) years to time points 1..J
  # timept <- as.integer(ordered(year))
  # J <- max(timept)
  #
  # # Changepoints must be converted to 0, 1..J-1 if not already so
  # if (length(changepoints)==0) {
  #   # case 0: no changepoints
  #   cpts <- 0L
  # } else if (min(changepoints)>=0 && max(changepoints)<J) {
  #   # case 1: already in 0..J-1 or 1..J-1
  #   cpts <- changepoints
  #   if (cpts[1]>0) cpts <- c(0L, cpts) # add prefix 0
  # } else if (all(changepoints %in% year)) {
  #   # case 2: actual years (used); convert to 1..J-1
  #   unique_years <- sort(unique(year))
  #   cpts <- match(changepoints, unique_years)
  #   cpts <- c(0L, cpts) # add prefix 0
  # } else {
  #   stop("Invalid changepoints specified:", vfmt(changepoints, 20))
  # }


  # get the index of the first changepoint to delete
  idx <- get_deletion(count, year, changepoints, covars, dbg=dbg)
  stopifnot(is.integer(idx))
  niter <- 1L
  while (idx > 0L) {
    # !! add 1 to cpmpensate for the implict changepoint '0'
    idx <- idx + 1L

    # Delete this changepoint.

    # first a check
    if (idx > length(changepoints)) {
      # 260109: Deze aanpak crasht bij trailing no-pos-data cases (zie usercase 26010108 Adriaan)
      # omdat het laatste piece een index van ncpt+1 heeft.
      msg <- sprintf("Can't happen: index %d too large in autodelete()", idx)
      stop(msg, call.=FALSE)
    }

    # move on to actual deletion
    if (dbg) rprintf("Auto-deleting change point %d : %d\n", idx, changepoints[idx])
    changepoints <- changepoints[-idx] # was: changepoints[changepoints != out]
    if (dbg) rprintf("  Changepoints remaining: %s\n", vfmt(changepoints, 20))

    # get the index of the next changepoint to delete
    idx <- get_deletion(count, year, changepoints, covars, dbg=dbg)

    # prevent infinite loops.
    niter <- niter + 1L
    if (niter > 100L) stop("Infinite loop in autodelete()")
  }
  if (dbg) rprintf("  Returning: %s\n", vfmt(changepoints))
  return(changepoints)
}

old_autodelete <- function(count, year, changepoints, covars=NULL, dbg=T) {
  if (dbg) {
    cat("autodelete()\n")
    cat("  count:"); str(count)
    cat("  year: ");  str(year)
    cat("  cpts: ");  str(changepoints)
  }
  # get the index of the first changepoint to delete
  idx <- get_deletion(count, year, changepoints, covars, dbg=dbg)
  stopifnot(is.integer(idx))
  niter <- 1L
  while (idx > 0L) {
    # Delete this changepoint.

    # first a check
    if (idx > length(changepoints)) {
      # 260109: Deze aanpak crasht bij trailing no-pos-data cases (zie usercase 26010108 Adriaan)
      # omdat het laatste piece een index van ncpt+1 heeft.
      msg <- sprintf("Can't happen: index %d too large in autodelete()", idx)
      stop(msg, call.=FALSE)
    }

    # move on to actual deletion
    if (dbg) rprintf("Auto-deleting change point %d : %d\n", idx, changepoints[idx])
    changepoints <- changepoints[-idx] # was: changepoints[changepoints != out]

    # get the index of the next changepoint to delete
    idx <- get_deletion(count, year, changepoints, covars, dbg=dbg)

    # prevent infinite loops.
    niter <- niter + 1L
    if (niter > 100L) stop("Infinite loop in autodelete()")
  }
  if (dbg) rprintf("  Returning: %s\n", vfmt(changepoints))
  return(changepoints)
}

prev_autodelete <- function(count, year, changepoints=NULL, covars=NULL, min_obs=3L, dbg=T) {
  if (dbg) rprintf("\nAutodelete()\n")
  # complete fresh take on autodelete.
  # we now analyse if for each piece there will be at least two time points with
  # positive counts

  min_piece <- 3L

  # if there are no changepoints, then we're OK
  if (is.null(changepoints)) return(0L)

  # Step 1: convert years to time points 1..J
  # convert actual time from (possibly non-contiguous) years to time points 1..J
  unique_years <- sort(unique(year))
  J <- length(unique_years)
  timept <- match(year, unique_years)
  if (dbg) rprintf("  count:  %s\n", vfmt(count,   20))
  if (dbg) rprintf("  years:  %s\n", vfmt(year,   20))
  if (dbg) rprintf("  timept: %s\n", vfmt(timept, 20))
  stopifnot(min(timept)==1L)
  stopifnot(max(timept)==J)

  # Step 2: convert changepoints from years to 1..(J-1) if nescessary
  if (min(changepoints)>=1L & max(changepoints)<J) {
    # not needed
    cpts <- changepoints
  } else if (min(changepoints)>=min(unique_years) & max(changepoints)<=max(unique_years)) {
    cpts <- match(changepoints, unique_years)
  } else {
    stop("Invalid changepoints:", vfmt(changepoints, 20))
  }
  if (dbg) rprintf("  changepoints:  %s\n", vfmt(changepoints,   20))
  if (dbg) rprintf("  cpts:          %s\n", vfmt(cpts, 20))

  # Step 3: count for each time points how many positive counts there are
  # count_ok <- ifelse(is.na(v), FALSE, v>0) # usable or not
  # for (j in 1:J) {
  #   idx <- timept==j
  # }
  .count_pos <- function(x) {
    # count number of positive elements in vector x.
    ispos <- x > 0
    sum(ispos, na.rm=TRUE)
  }
  pos_counts <- unname(tapply(count, list(j=timept), .count_pos)) # > 0
  pos_pos    <-  pos_counts > 0

  if (dbg) rprintf("  pos. counts:   %s\n", vfmt(pos_counts, 20))

  max_iter <- 2 * length(changepoints)
  for (iter in 1:max_iter) {
    if (dbg) rprintf("  Iteration %d\n", iter)

    # Step 4: for each piece we need at least two years with usable data
    P <- length(cpts) # number of pieces (last piece will be from the )
    if (P==0) stop("No changepoints left!")

    pfil_yr  <- integer(P) # records for each piece the # of years with data
    pfil_obs <- integer(P) # records for each piece the # of positive observations
    pok      <- logical(P)

    # note:
    # - pfil_yr  must be at least 2 (two years with data)
    # - pfil_obs must be at least 3 (three data points)

    for (p in 1:P) {
      cp1 <- cpts[p] # starting time point of the piece
      if (p<P) { # intermediate piece
        cp2 <- cpts[p+1] #ending time point of the piece
      } else if (cp1<J) { # change point before end of time series
        cp2 <- J
      } else if (cp1==J) { # change point at end of time series
        cp2 <- J
      } else {
        stop("Can't happen")
      }
      y1 <- unique_years[cp1]
      y2 <- unique_years[cp2]
      pfil_yr[p]  <- sum(pos_pos[cp1:cp2])
      pfil_obs[p] <- sum(pos_counts[cp1:cp2])
      ok <- pfil_yr[p] >= 2 & pfil_obs[p] >= min_obs
      pok[p] <- ok
      if (dbg)  rprintf("    piece %d (%d:%d - %d:%d) #yr=%d #obs=%d : %s\n",
                        p, cp1, y1, cp2, y2, pfil_yr[p], pfil_obs[p],
                        ifelse(ok,"ok", "NOT ok"))
    }
    if (dbg) rprintf("    pfil_yr:  %s\n", vfmt(pfil_yr,  20))
    if (dbg) rprintf("    pfil_obs: %s\n", vfmt(pfil_obs, 20))
    if (dbg) rprintf("    pok: %s\n", vfmt(pok, 20))

    if (all(pok)) { # All changepoints are OK
      if (dbg) rprintf("    All changepoints OK: %s\n", vfmt(changepoints, 20))
      return(changepoints)
    }

    # Not all OK; delete one of the changepoints

    # find the last piece that has insufficient data
    idx <- tail(which(!pok), 1)
    if (idx < length(changepoints)) idx <- idx + 1L # remove *tailing* changepoint of the piece
    # remove it,
    if (dbg) rprintf("    Removing changepoint #%d: %d\n", idx, changepoints[idx])
    changepoints <- changepoints[-idx]
    cpts <- cpts[-idx]
    if (dbg) rprintf("    Changepoints remaining: %s\n", vfmt(changepoints, 20))

    # next iteration
  }

  stop("Taking too long")
}


check_model3 <- function(count, site, year, month=NULL) {
  ndata <- length(count)
  use_months <- !is.null(month)

  # sites must be factors at this stage!
  stopifnot(inherits(site,"factor"))
  I <-  nsite  <- nlevels(site)
  site_nr <- as.integer(site) # 1, 2, ..., I

  # same for years
  stopifnot(inherits(year,"factor"))
  J <- nyear <- nlevels(year)
  year_nr <- as.integer(year) # 1, 2, ..., J

  if (use_months) {
    stopifnot(inherits(month,"factor"))
    M <- nmonth <- nlevels(month)
    month_nr <- as.integer(month) # 1, 2, ... M
  }

  # Create observation matrix $f$.
  # Convert the data from a vector representation to a matrix representation.
  # It's OK to have missing site/time combinations; these will automatically
  # translate to NA values.
  if (use_months) {
    f <- array(NA, dim=c(nsite, nyear, nmonth))
    for (m in 1:M) {
      fm <- matrix(NA, nsite, nyear)
      midx <- month_nr == m # month factor -> 1,2,3,etc
      rows <- site_nr[midx]
      cols <- year_nr[midx]
      idx <- (cols-1)*nsite+rows
      fm[idx] <- count[midx]
      f[ , ,m] <- fm
    }
    dimnames(f)[[1]] <- levels(site)
    dimnames(f)[[2]] <- levels(year)
    dimnames(f)[[3]] <- levels(month)
  } else {
    f <- matrix(NA, nsite, nyear)
    rows <- site_nr # works because site_nr = 1...I
    cols <- year_nr # idem
    idx <- (cols-1)*nsite+rows   # Create column-major linear index from row/column subscripts.
    f[idx] <- count    # ... such that we can paste all data into the right positions
    dimnames(f)[[1]] <- levels(site)
    dimnames(f)[[2]] <- levels(year)
  }

  # Availability matrix (i.e/, availability for beta parameters)
  if (use_months) {
    avail <- array(FALSE, dim=c(nsite, nyear, nmonth))
    avail[f > 0] <- TRUE
  } else {
    avail <- matrix(FALSE, nsite, nyear)
    avail[f > 0] <- TRUE
  }

  # without months: simple check/adjustment of availability.
  # first sites, then years.
  if (!use_months) {
    # first check/adjust for sites
    for (i in 1:nsite) {
      n <- sum(avail[i, ])
      print(n)
      if (n==0) {
        msg <- sprintf("No data available for site #%d (%s)", i, levels(site)[i])
        stop(msg)
      }
      if (n==1) {
        msg <- sprintf("Single pos.obs for site #%d (%s); removing availability for another purposes\n",
                i, levels(site)[i])
        warning(msg)
        avail[i, ] <- FALSE
      }
      # ok; more than 1 pos.obs for this site.
    }
    # then check for years
    for (j in 1:nyear) {
      n <- sum(avail[i, j])
      if (n==0) {
        msg <- sprintf("No data available for year #%d (%s)", j, levels(year)[j])
        stop(msg)
      }
    }
  }

  # with months:
  # - 1:
  # - 2: first check for dual responsibilities (site month)
  # - 3: then remove avaiablity for one-obs sites / monthts
  # - 4:finally check years

  if (use_months) {

    # step 1: check for empty sites, years, and months
    for (i in 1:nsite) {
      n <- sum(avail[i,,])
      if (n==0) {
        msg <- sprintf("Empty site #%d (%s)", i, levels(site)[i])
        stop(msg)
      }
    }
    for (j in 1:nyear) {
      n <- sum(avail[,j,])
      if (n==0) {
        msg <- sprintf("Empty year #%d (%s)", j, levels(year)[j])
        stop(msg)
      }
    }
    for (m in 1:nmonth) {
      n <- sum(avail[,,m])
      if (n==0) {
        msg <- sprintf("Empty month #%d (%s)", m, levels(month)[m])
        stop(msg)
      }
    }

    # step 2: dual responsibilities
    for (i in 1:nsite) {
      ni <- sum(avail[i,,])
      if (ni==1) {
        for (m in 1:nmonth) {
          nim <- sum(avail[i,,m])
          if (nim==1) {
            msg <- sprintf("Single pos.obs for site #%d (%s) *and* month #%d (%s)",
                           i, levels(site)[i], m, levels(month)[m])
            stop(msg)
          }
        }
      }
    }

    # Step 3: remove availability for single-pos.obs sites and monthts
    for (i in 1:nsite) {
      n <- sum(avail[i,,])
      if (n==1) {
        msg <- sprintf("Removing availability for site #%d (%s)\n", i, levels(site)[i])
        rprintf(msg)
        avail[i,,] <- FALSE
      }
    }
    for (m in 1:nmonth) {
      n <- sum(avail[,,m])
      if (n==1) {
        msg <- sprintf("Removing availability for month #%d (%s)\n", m, levels(month)[m])
        rprintf(msg)
        avail[,,m] <- FALSE
      }
    }

    # Step 4: re-test years
    for (j in 1:nyear) {
      n <- sum(avail[,j,])
      if (n==0) {
        msg <- sprintf("Not sufficient data available for year #%d (%s)", j, levels(year)[j])
        stop(msg)
      }
    }
  }
  rprintf("OK\n")
  return(0)
}

rprintf <- function(fmt, ...) cat(sprintf(fmt, ...))

# #should work
# count <- 1:5
# site <-  rep("A", 5)
# year <- 1:5
# check_model3(count, factor(site), factor(year))
#
# # should crash at site 2
# count <- c(1,1,0,0)
# site <- c(1,1,2,2)
# year <- c(1,2,1,2)
# check_model3(count, factor(site), factor(year))

# # should crash at year 2
# count <- c(1,0,1,0)
# site <- c(1,1,2,2)
# year <- c(1,2,1,2)
# check_model3(count, factor(site), factor(year))

# # should crash at site1 / month 1 ((dual resp))
# s1m1 <- data.frame(site=1, month=1, count=c(1,0,0), year=1:3)
# s1m2 <- data.frame(site=1, month=2, count=c(0,0,0), year=1:3)
# s2m1 <- data.frame(site=2, month=1, count=c(0,0,0), year=1:3)
# s2m2 <- data.frame(site=2, month=2, count=c(1,1,1), year=1:3)
# df <- rbind(s1m1, s1m2,s2m1,s2m2)
# print(df)
# check_model3(df$count, factor(df$site), factor(df$year), month=factor(df$month))

# stop("Intended")

check_model2 <- function(count, site, year, month=NULL, changepoints=0, autodelete=FALSE, min_obs, dbg=T) {
  ndata <- length(count)
  use_months <- !is.null(month)

  # prepare output
  out <- list(
    ok <- TRUE,
    msg <- 0,
    deletion <- 0L
  )

  # sites must be factors at this stage!
  stopifnot(inherits(site,"factor"))
  I <-  nsite  <- nlevels(site)
  site_nr <- as.integer(site) # 1, 2, ..., I

  # same for years
  stopifnot(inherits(year,"factor"))
  J <- nyear <- nlevels(year)
  year_nr <- as.integer(year) # 1, 2, ..., J

  if (use_months) {
    stopifnot(inherits(month,"factor"))
    M <- nmonth <- nlevels(month)
    month_nr <- as.integer(month) # 1, 2, ... M
  }

  # Step 2: convert changepoints from years to 1..(J-1) if necessary
  uyear <- as.integer(levels(year)) # unique years
  if (min(changepoints)>=1L & max(changepoints)<J) {
    # not needed; already in OK format
    cpt <- changepoints
    change_years <- FALSE # to assist in reporting
  } else if (min(changepoints)>=min(uyear) & max(changepoints) < max(uyear)) {
    cpt <- match(changepoints, uyear)
    change_years <- TRUE
  } else {
    stop("Invalid changepoints:", vfmt(changepoints, 20))
  }
  # if (dbg) rprintf("  changepoints:  %s\n", vfmt(changepoints, 20))
  # if (dbg) rprintf("  cpts:          %s\n", vfmt(cpt, 20))
  ncpt <- length(cpt)

  # year numbers j to piece number p
  # year2piece <- integer(J)
  # for (p in 1:ncpt) {
  #   idx <- year_nr > cpt[p]
  #   year2piece[idx] <- p
  # }
  # if (dbg) rprintf("  year2piece:    %s\n", vfmt(year2piece, 20))
  # explicit start-stop indices (except for [piece '0'])
  pieces <- data.frame(changepoint=changepoints, cpt=cpt, start=cpt+1L)
  #pieces$stop <- ifelse(ncpt==1, J, c(cpt[2:ncpt], J))
  if (ncpt==1) {
    pieces$stop <- J
  } else {
    pieces$stop <- c(cpt[2:ncpt], J)
  }
  pieces$nobs <- 0L
  pieces$ok <- FALSE

  # Create observation matrix $f$.
  # Convert the data from a vector representation to a matrix representation.
  # It's OK to have missing site/time combinations; these will automatically
  # translate to NA values.
  if (use_months) {
    f <- array(NA, dim=c(nsite, nyear, nmonth))
    for (m in 1:M) {
      fm <- matrix(NA, nsite, nyear)
      midx <- month_nr == m # month factor -> 1,2,3,etc
      rows <- site_nr[midx]
      cols <- year_nr[midx]
      idx <- (cols-1)*nsite+rows
      fm[idx] <- count[midx]
      f[ , ,m] <- fm
    }
    dimnames(f)[[1]] <- levels(site)
    dimnames(f)[[2]] <- levels(year)
    dimnames(f)[[3]] <- levels(month)
  } else {
    f <- matrix(NA, nsite, nyear)
    rows <- site_nr # works because site_nr = 1...I
    cols <- year_nr # idem
    idx <- (cols-1)*nsite+rows   # Create column-major linear index from row/column subscripts.
    f[idx] <- count    # ... such that we can paste all data into the right positions
    dimnames(f)[[1]] <- levels(site)
    dimnames(f)[[2]] <- levels(year)
  }

  # Availability matrix (i.e/, availability for beta parameters)
  if (use_months) {
    avail <- array(FALSE, dim=c(nsite, nyear, nmonth))
    avail[f > 0] <- TRUE
  } else {
    avail <- matrix(FALSE, nsite, nyear)
    avail[f > 0] <- TRUE
  }

  # Now proceed to the actual checking.

  # without months:
  # - Check if sites are not empty.
  # - Remove availability for single-obs sites
  # - Check pieces; if any not ok; issue error, or propose a deletion

  if (!use_months) {
    # first check/adjust for sites
    for (i in 1:nsite) {
      n <- sum(avail[i, ])
      print(n)
      if (n==0) {
        msg <- sprintf("No data available for site #%d (%s)", i, levels(site)[i])
        stop(msg)
      }
      if (n==1) {
        j <- which(avail[i, ])
        msg <- sprintf("Removing availability of single pos.obs for site #%d (%s) / year #%d (%s)\n",
                       i, levels(site)[i], j, levels(year)[j])
        if (dbg) rprintf(msg)
        avail[i, ] <- FALSE
      }
      # ok; more than 1 pos.obs for this site.
    }

    # Check all pieces; no further response yet.
    for (p in 1:ncpt) {
      j1 <- pieces$start[p]
      j2 <- pieces$stop[p]
      n <- sum(avail[,j1:j2])
      pieces$nobs[p] <- n
      pieces$ok[p] <- n >= min_obs
    }
  }

  # with months:
  # - 1: check empty sites and months
  # - 2: first check for dual responsibilities (site month)
  # - 3: then remove availability for one-obs sites / months
  # - 4: finally check years

  if (use_months) {

    # step 1: check for empty sites, and months
    for (i in 1:nsite) {
      n <- sum(avail[i,,])
      if (n==0) {
        msg <- sprintf("Empty site #%d (%s)", i, levels(site)[i])
        stop(msg)
      }
    }
    for (m in 1:nmonth) {
      n <- sum(avail[,,m])
      if (n==0) {
        msg <- sprintf("Empty month #%d (%s)", m, levels(month)[m])
        stop(msg)
      }
    }

    # step 2: dual responsibilities
    for (i in 1:nsite) {
      ni <- sum(avail[i,,])
      if (ni==1) {
        for (m in 1:nmonth) {
          nim <- sum(avail[i,,m])
          if (nim==1) {
            msg <- sprintf("Single pos.obs for site #%d (%s) *and* month #%d (%s)",
                           i, levels(site)[i], m, levels(month)[m])
            stop(msg)
          }
        }
      }
    }

    savail <- mavail <- pavail <- avail

    # Step 3: remove availability for single-pos.obs sites and months
    for (i in 1:nsite) {
      n <- sum(avail[i,,])
      if (n==1) {
        msg <- sprintf("Removing availability of site #%d (%s)\n", i, levels(site)[i])
        if (dbg) rprintf(msg)
        mavail[i,,] <- FALSE
        pavail[i,,] <- FALSE
      }
    }
    for (m in 1:nmonth) {
      n <- sum(avail[,,m])
      if (n==1) {
        msg <- sprintf("Removing availability of month #%d (%s)\n", m, levels(month)[m])
        if (dbg) rprintf(msg)
        savail[,,m] <- FALSE
        pavail[,,m] <- FALSE
      }
    }

    # Step 4: re-test individual months and sites
    for (i in 1:nsite) {
      n <- sum(savail[i,,])
      if (n==0) {
        msg <- sprintf("Empty site #%d (%s)", i, levels(site)[i])
        stop(msg)
      }
    }
    for (m in 1:nmonth) {
      n <- sum(mavail[,,m])
      if (n==0) {
        msg <- sprintf("Empty month #%d (%s)", m, levels(month)[m])
        stop(msg)
      }
    }


    # step 5: test changepoints
    # Check all pieces; no further response yet.
    for (p in 1:ncpt) {
      j1 <- pieces$start[p]
      j2 <- pieces$stop[p]
      n <- sum(pavail[ ,j1:j2, ])
      pieces$nobs[p] <- n
      pieces$ok[p] <- n >= min_obs
    }
  }

  # Treatment of pieces that are not OK is similar for with/out months

  # everything OK?
  if (all(pieces$ok)) return(list(ok=TRUE))

  # Not OK; Set up error messages.
  idx <- which(!pieces$ok)
  if (length(idx)==1) {
    msg <- sprintf("Not enough data for changepoint %s.", pieces$changepoint[idx])
  } else {
    msg <- sprintf("Not enough data for changepoints %s.", vfmt(pieces$changepoint[idx], 20))
  }

  # Suggest to to remove the LAST problematic changepoint
  suggestion <- tail(idx, 1)
  if (suggestion==1) {
    if (nrow(pieces)==1) stop("Autodelete problem: should not happen")
    suggestion <- suggestion + 1L
  }

  # Either issue an error or return to caller
  if (autodelete) {
    return(list(ok=FALSE, msg=msg, deletion=suggestion))
  } else {
    stop(msg)
  }
}

check_model2_simple <- function(count, site, year, month=NULL, changepoints=0, autodelete=FALSE, min_obs, dbg=T) {
  ndata <- length(count)
  use_months <- !is.null(month)

  # prepare output
  out <- list(
    ok <- TRUE,
    msg <- 0,
    deletion <- 0L
  )

  # sites must be factors at this stage!
  stopifnot(inherits(site,"factor"))
  I <-  nsite  <- nlevels(site)
  site_nr <- as.integer(site) # 1, 2, ..., I

  # same for years
  stopifnot(inherits(year,"factor"))
  J <- nyear <- nlevels(year)
  year_nr <- as.integer(year) # 1, 2, ..., J

  if (use_months) {
    stopifnot(inherits(month,"factor"))
    M <- nmonth <- nlevels(month)
    month_nr <- as.integer(month) # 1, 2, ... M
  }

  # Step 2: convert changepoints from years to 1..(J-1) if necessary
  uyear <- as.integer(levels(year)) # unique years
  if (min(changepoints)>=1L & max(changepoints)<J) {
    # not needed; already in OK format
    cpt <- changepoints
    change_years <- FALSE # to assist in reporting
  } else if (min(changepoints)>=min(uyear) & max(changepoints) < max(uyear)) {
    cpt <- match(changepoints, uyear)
    change_years <- TRUE
  } else {
    stop("Invalid changepoints:", vfmt(changepoints, 20))
  }
  if (dbg) rprintf("  changepoints:  %s\n", vfmt(changepoints, 20))
  if (dbg) rprintf("  cpts:          %s\n", vfmt(cpt, 20))
  ncpt <- length(cpt)

  # year numbers j to piece number p
  # year2piece <- integer(J)
  # for (p in 1:ncpt) {
  #   idx <- year_nr > cpt[p]
  #   year2piece[idx] <- p
  # }
  # if (dbg) rprintf("  year2piece:    %s\n", vfmt(year2piece, 20))
  # explicit start-stop indices (except for [piece '0'])
  pieces <- data.frame(changepoint=changepoints, cpt=cpt, start=cpt+1L)
  if (ncpt==1) {
    pieces$stop <- J
  } else {
    pieces$stop <- c(cpt[2:ncpt], J)
  }
  #pieces$stop <- ifelse(ncpt==1, J, c(cpt[2:ncpt], J))
  pieces$nobs <- 0L
  pieces$ok <- FALSE
  print(pieces)
  stopifnot(all(pieces$stop >= pieces$start))

  # Create observation matrix $f$.
  # Convert the data from a vector representation to a matrix representation.
  # It's OK to have missing site/time combinations; these will automatically
  # translate to NA values.
  if (use_months) {
    f <- array(NA, dim=c(nsite, nyear, nmonth))
    for (m in 1:M) {
      fm <- matrix(NA, nsite, nyear)
      midx <- month_nr == m # month factor -> 1,2,3,etc
      rows <- site_nr[midx]
      cols <- year_nr[midx]
      idx <- (cols-1)*nsite+rows
      fm[idx] <- count[midx]
      f[ , ,m] <- fm
    }
    dimnames(f)[[1]] <- levels(site)
    dimnames(f)[[2]] <- levels(year)
    dimnames(f)[[3]] <- levels(month)
  } else {
    f <- matrix(NA, nsite, nyear)
    rows <- site_nr # works because site_nr = 1...I
    cols <- year_nr # idem
    idx <- (cols-1)*nsite+rows   # Create column-major linear index from row/column subscripts.
    f[idx] <- count    # ... such that we can paste all data into the right positions
    dimnames(f)[[1]] <- levels(site)
    dimnames(f)[[2]] <- levels(year)
  }

  # Availability matrix (i.e/, availability for beta parameters)
  if (use_months) {
    avail <- array(FALSE, dim=c(nsite, nyear, nmonth))
    avail[f > 0] <- TRUE
  } else {
    avail <- matrix(FALSE, nsite, nyear)
    avail[f > 0] <- TRUE
  }

  # Now proceed to the actual checking.

  # without months:
  # - Check if sites are not empty.
  # - Remove availability for single-obs sites
  # - Check pieces; if any not ok; issue error, or propose a deletion

  if (!use_months) {
    # first check/adjust for sites
    for (i in 1:nsite) {
      n <- sum(avail[i, ])
      print(n)
      if (n==0) {
        msg <- sprintf("No data available for site #%d (%s)", i, levels(site)[i])
        stop(msg)
      }
      if (n==1) {
        j <- which(avail[i, ])
        msg <- sprintf("Removing availability of single pos.obs for site #%d (%s) / year #%d (%s)\n",
                       i, levels(site)[i], j, levels(year)[j])
        if (dbg) rprintf(msg)
        avail[i, ] <- FALSE
      }
      # ok; more than 1 pos.obs for this site.
    }

    # Check all pieces; no further response yet.
    for (p in 1:ncpt) {
      j1 <- pieces$start[p]
      j2 <- pieces$stop[p]
      n <- sum(avail[,j1:j2])
      pieces$nobs[p] <- n
      pieces$ok[p] <- n >= min_obs
    }
  }

  # with months:
  # - 1: check empty sites and months
  # - 2: first check for dual responsibilities (site month)
  # - 3: then remove availability for one-obs sites / months
  # - 4: finally check years

  if (use_months) {

    # step 1: check for empty sites, and months
    for (i in 1:nsite) {
      n <- sum(avail[i,,])
      if (n==0) {
        msg <- sprintf("Empty site #%d (%s)", i, levels(site)[i])
        stop(msg)
      }
    }
    for (m in 1:nmonth) {
      n <- sum(avail[,,m])
      if (n==0) {
        msg <- sprintf("Empty month #%d (%s)", m, levels(month)[m])
        stop(msg)
      }
    }

    # step 2: dual responsibilities
    for (i in 1:nsite) {
      ni <- sum(avail[i,,])
      if (ni==1) {
        for (m in 1:nmonth) {
          nim <- sum(avail[i,,m])
          if (nim==1) {
            msg <- sprintf("Single pos.obs for site #%d (%s) *and* month #%d (%s)",
                           i, levels(site)[i], m, levels(month)[m])
            stop(msg)
          }
        }
      }
    }


    # Step 3: remove availability for single-pos.obs sites and months
    for (i in 1:nsite) {
      n <- sum(avail[i,,])
      if (n==1) {
        msg <- sprintf("Removing availability of site #%d (%s)\n", i, levels(site)[i])
        if (dbg) rprintf(msg)
        avail[i,,] <- FALSE
      }
    }
    for (m in 1:nmonth) {
      n <- sum(avail[,,m])
      if (n==1) {
        msg <- sprintf("Removing availability of month #%d (%s)\n", m, levels(month)[m])
        if (dbg) rprintf(msg)
        avail[,,m] <- FALSE
      }
    }

    # # Step 4: re-test individual months and sites
    # for (i in 1:nsite) {
    #   n <- sum(avail[i,,])
    #   if (n==0) {
    #     msg <- sprintf("Empty site #%d (%s)", i, levels(site)[i])
    #     stop(msg)
    #   }
    # }
    # for (m in 1:nmonth) {
    #   n <- sum(avail[,,m])
    #   if (n==0) {
    #     msg <- sprintf("Empty month #%d (%s)", m, levels(month)[m])
    #     stop(msg)
    #   }
    # }


    # step 5: test changepoints
    # Check all pieces; no further response yet.
    for (p in 1:ncpt) {
      j1 <- pieces$start[p]
      j2 <- pieces$stop[p]
      n <- sum(avail[ ,j1:j2, ])
      pieces$nobs[p] <- n
      pieces$ok[p] <- n >= min_obs
    }
  }
  print(pieces)

  # Treatment of pieces that are not OK is similar for with/out months

  # everything OK?
  if (all(pieces$ok)) return(list(ok=TRUE))

  # Not OK; Set up error messages.
  idx <- which(!pieces$ok)
  if (length(idx)==1) {
    msg <- sprintf("Not enough data for changepoint %s.", pieces$changepoint[idx])
  } else {
    msg <- sprintf("Not enough data for changepoints %s.", vfmt(pieces$changepoint[idx], 20))
  }

  # Suggest to to remove the LAST problematic changepoint
  suggestion <- tail(idx, 1)
  if (suggestion==1) {
    if (nrow(pieces)==1) stop("Autodelete problem: should not happen")
    suggestion <- suggestion + 1L
  }

  # Either issue an error or return to caller
  if (autodelete) {
    return(list(ok=FALSE, msg=msg, deletion=suggestion))
  } else {
    stop(msg)
  }
}

autodelete <- function (count, site, year, month, changepoints, min_obs=1L, dbg=TRUE) {

  # prepare sites
  if (inherits(site, "numeric")) site <- as.integer(site)
  if (inherits(site, "integer")) {
    site <- factor(site)                # Convert to a factor
    site_id <- as.integer(levels(site)) # Original values
  } else if (inherits(site,"character")) {
    site <- factor(site)
    site_id <- levels(site)
  } else if (inherits(site,"factor")) {
    site <- factor(site) # Refactor to get rid of unused levels
    site_id <- levels(site)
  } else {
    msg <- sprintf("Invalid site class: %s", paste0(class(site), collapse=","))
    stop(msg, call.=FALSE)
  }

  # prepare years
  # Years must be integers or numerics with a constant step size.
  # Convert numerics that are integers in disguise
  if (inherits(year, "numeric")) year <- as.integer(year)
  if (inherits(year,"integer")) {
    # Integers are allowed iff they have a constant step size
    delta <- unique(diff(sort(unique(year))))
    if (length(delta)!=1L) stop("Years don't have a constant interval")
    year <- ordered(year)
    year_id <- as.integer(levels(year))
  } else  if (inherits(year, "numeric")) {
    # Idem for numerics
    delta <- unique(diff(sort(unique(year))))
    if (length(delta)!=1L) stop("Years don't have a constant interval")
    year <- ordered(year)
    year_id <- as.numeric(levels(year))
  } else {
    msg <- sprintf("Invalid year class: %s", paste0(class(year), collapse=","))
    stop(msg, call.=FALSE)
  }

  # prepare months
  if (is.null(month)) {
    use_months <- FALSE
  } else {
    use_months <- TRUE

    # Convert numerics that are integers in disguise
    if (inherits(month, "numeric")) month <- as.integer(month)
    if (inherits(month,"integer")) {
      month <- ordered(month)
      month_id  <- as.integer(levels(month)) # the original month identifiers
    } else if (inherits(month, "character")) {
      use.months <- TRUE
      month <- ordered(month, levels=unique(month)) # use order of appearance
      month_id <- levels(month)
    } else if (inherits(month,"ordered")) {
      use.months <- TRUE
      month <- ordered(month) # Get rid of unused levels
      month_id <- ordered(levels(month), levels(month))
    } else if (inherits(month,"factor")) {
      use.months <- TRUE
      month <- factor(month) # Get rid of unused levels
      month_id <- ordered(levels(month), levels(month))
    } else {
      msg <- sprintf("Invalid year class: %s", paste0(class(year), collapse=","))
      stop(msg, call.=FALSE)
    }
  }

  for (iter in 1:1000) {
    if (dbg) rprintf("Autodelete() iteration %d\n", iter)

    out <- check_model2(count, site, year, month, changepoints, autodelete=TRUE, min_obs=min_obs, dbg=dbg)
    if (out$ok) break
    # not ok; issue message and remove a changepoint
    if (dbg) rprintf("%s\n", out$msg)
    idx <- out$deletion
    if (dbg) rprintf("Removing changepoint #%d (%s)\n", idx, changepoints[idx])
    changepoints <- changepoints[-idx]
  }
  if (dbg) rprintf("Autodeletion complete\n")
  return(changepoints)
}

Main <- function(run=F) {
  if (!run) return()

  # meest minimale dataset
  count <- c(0, 0, 1, 1)
  site  <- c(1, 1, 1, 1)
  year  <- c(1, 2, 3, 4) + 2000
  cpts  <- c(1, 2, 3) + 2000
  autodelete(count, site, year, NULL, cpts)

  # count <- c(1, 0, 1, 2, 0, NA)
  # year  <- c(1, 2, 4, 1, 2, 4) + 1999
  # cpts  <- c(1,2)
  # new_autodelete(count, year, cpts)
}

rprintf <- function(fmt, ...) cat(sprintf(fmt,...))
vfmt <- function(v, maxlen=9) {
  n <- length(v)
  if (n>maxlen) {
    s <- sprintf("%s ... %s", toString(v[1:(maxlen-4)]), toString(v[(n-2):n]))
  } else {
    s <- toString(v)
  }
  paste0("[", s, "]")
}

Main(F)

#
# load("../tests/testthat/testdata/131183.RData")
#
# count <- df$count
# year <- df$year
# cpts <- sort(unique(df$year))
# J <- max(as.integer(ordered(year)))
# cpts <- cpts[1:J-1]
#
# out <- autodelete(count, year, cpts)
# print(out)

#trim(count ~ site + year, data=df, model=2, overdisp=TRUE, serialcor=TRUE, changepoints="all", autodelete=TRUE)

# # test code 260109-1 - Triggert probleem ana het eind.
# cat("\n\n--- testing with time 1... --- should be 4\n")
# time  <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
# count <- c(1, 1, 1, 0, 1, 1, 1, 1, 1,  2)
# cpts  <- c(1, 2, 3, 4, 5, 6, 7, 8, 9    )
# #cpts <- c(0, 1,4,7)
# out <- autodelete(count, time, cpts)
# print(out)

# new test code 2025
# test code
# cat("\n\n--- testing with time 1... --- should be 4\n")
# time  <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
# count <- c(1, 1, 1, 1, 0, 0, 0, 1, 9,  9)
# cpts  <- c(1,       4,       7          ) # cpt 7 is removed to provide data to cpt 4
# out <- autodelete(count, time, cpts)
# print(out)
#
# cat("\n\n--- testing with time 1... --- should be 4\n")
# time  <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
# count <- c(1, 1, 1, 1,NA,NA,NA, 1, 1,  1)
# cpts  <- c(1,       4,       7          ) # cpt 7 is removed to provide data to cpt 4
# out <- autodelete(count, time, cpts)
# print(out)

# # test code
# cat("\n\n--- testing with time 1... --- should be 4\n")
# time  <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
# count <- c(1, 1, 1, 1, 0, 0, 0, 1, 1,  1)
# cpts  <- c(         4,       7          ) # cpt 7 is removed to provide data to cpt 4
# out <- autodelete(count, time, cpts)
# print(out)
#
# cat("\n\n--- testing with time 10... --- should be 4\n")
# time  <- 10:19
# count <- c(1, 1, 1, 1, 0, 0, 0, 1, 1,  1)
# cpts  <- c(         4,       7          ) # cpt 7 is removed to provide data to cpt 4
# out <- autodelete(count, time, cpts)
# print(out)
#
# cat("\n\n--- testing with time 10... --- should be 1,4\n")
# time  <- 10:19
# count <- c(1, 1, 1, 1, 0, 0, 0, 1, 1,  1)
# cpts  <- c(1,       4,       7          ) # cpt 7 is removed to provide data to cpt 4
# out <- autodelete(count, time, cpts)
# print(out)
#
# cat("\n\n--- testing with cpts in years --- should be 13\n")
# cpts  <- c(13, 16) # same, but in years
# out <- autodelete(count, time, cpts)
# print(out)
#
# cat("\n\n--- with covars / all fine --- should be 4,7\n")
# time  <- rep(1:10, times=2)
# covar <- rep(letters[1:2], each=10)
# count <- rep(1, 20)
# out <- autodelete(count, time, c(4,7), covars=list(cov=covar))
# print(out)
#
# cat("\n\n--- with covars / delete 7 ---should be 4\n")
# count[8:10] <- 0
# out <- autodelete(count, time, c(4,7), covars=list(cov=covar))
# print(out)
#
# cat("\n\n--- with covars / delete 7 --- should be 1,4 \n")
# count[8:10] <- 0
# out <- autodelete(count, time, c(1,4,7), covars=list(cov=covar))
# print(out)
