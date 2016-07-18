#rm(list=ls())
library(stringr)

#printf <- function(...) cat(sprintf(...))

read_out = function(filename, debug=FALSE) {
  lines = trimws(readLines(filename))

  DBG=debug

  eat.white = function() {
    # remove empty lines
    while (lines[1]=="") lines <<- lines[-1]
  }

  eat <- function(pat, debug=DBG) {
    if (debug) cat("\n\n")

    pat = gsub("%d", "(\\\\d+)",          pat) # create RE for integer
    pat = gsub("%f", "(-?\\\\d+\\\\.\\\\d+)", pat) # idem floationg point
    pat = gsub("%s", "([[:print:]]+)",  pat) # idem string
    pat = gsub("\\s+","\\\\s+",         pat) # make all whitespace flexible
    pat <- paste0("^", pat, "$")

    eat.white() # remove empty lines

    # pop the current line
    curline = lines[1]
    lines <<- lines[-1]

    if (debug) {
      printf("Current line: %s\n", curline)
      printf("Pattern     : %s\n", pat)
    }

    # first assert that we do have the pattern
    stopifnot(grepl(pat, curline))

    out=""
    # extract requested groups, only for non-empty patterns
    if (pat!="^$") {
      m = str_match(curline, pat)
      if (debug) print(m)
      stopifnot(is.na(m[1,1])==FALSE)
      if (ncol(m)>1) out = m[1,2:ncol(m)]
    }
    if (debug) printf("output: %s\n", out)
    out
  }

  optional <- function(pat, debug=DBG) {
    if (debug) cat("\n\n")

    pat = gsub("%d", "(\\\\d+)",          pat) # create RE for integer
    pat = gsub("%f", "(-?\\\\d+\\\\.\\\\d+)", pat) # idem floationg point
    pat = gsub("%s", "([[:print:]]+)",  pat) # idem string
    pat = gsub("\\s+","\\\\s+",         pat) # make all whitespace flexible
    pat <- paste0("^", pat, "$")

    eat.white() # remove empty lines

    if (debug) {
      printf("Current line: %s\n", lines[1])
      printf("Pattern     : %s\n", pat)
    }

    if (grepl(pat, lines[1])) {
      if (debug) printf("Optional pattern found.\n")
      # pop the current line
      curline = lines[1]
      lines <<- lines[-1]
    } else {
      printf('Optional pattern not found.\n')
      return("")
    }

    # first assert that we do have the pattern
    stopifnot(grepl(pat, curline))

    out=""
    # extract requested groups, only for non-empty patterns
    if (pat!="^$") {
      m = str_match(curline, pat)
      if (debug) print(m)
      stopifnot(is.na(m[1,1])==FALSE)
      if (ncol(m)>1) out = m[1,2:ncol(m)]
    }
    if (debug) printf("output: %s\n", out)
    out
  }

  eat.table <- function(my.names=NULL, header=TRUE, fix=FALSE) {
    while (lines[1]=="") lines <<- lines[-1] # remove empty lines
    n = grep("^$", lines)[1] # detect table size
    if (is.null(my.names)) {
      df = read.table(textConnection(lines[1:n]), header=header, fill=TRUE)
    } else {
      df = read.table(textConnection(lines[1:n]), header=FALSE, fill=TRUE, skip=1)
      names(df) <- my.names
    }
    # fix line 1 (time indices, mainly)
    if (fix) {
      print(df)
      print(class(df))
      if (all.equal(is.finite(as.vector(df[1, ])), c(TRUE,TRUE,TRUE,FALSE,FALSE))) {
        df[1,4] = df[1,3]
        df[1,3] = 0
        df[1,5] = 0
      }
    }

    lines <<- lines[-(1:n)]
    df
  }

  z <- list()

  # Header
  eat("TRIM 3.61.*")
  printf("done\n")
  eat("STATISTICS NETHERLANDS")
  eat("Date.*")
  eat("Title : %s") -> title
  optional("Comment: %s") -> comment

  # Datafile and variables
  eat("The following %d variables have been read from file:") -> z$ncols
  eat("%s") -> data_path
  eat("\\d+\\. Site  number of values: %d") -> z$nsite
  eat("\\d+\\. Time  number of values: %d") -> z$ntime
  eat("\\d+\\. Count  missing = %d") -> z$missing
  eat(".*")
  eat(".*")

  # data stats
  eat("Number of sites without positive counts \\(removed\\)  %d") -> nerrsite
  eat("Number of observed zero counts  %d") -> z$nzero
  eat("Number of observed positive counts  %d") -> z$npos
  eat("Total number of observed counts  %d") -> z$nobs
  eat("Number of missing counts  %d") -> z$nmis
  eat("Total number of counts  %d") -> z$ncount
  eat("Total count  %f") -> z$totcount

  # Dominant sites
  z$dom = list()
  eat("Sites containing more than %d% of the total count") -> z$dom$threshold
  eat.table(c("site.nr", "obs.tot","percent")) -> z$dom$sites

  eat("Time Point Averages")
  eat.table() -> z$avg

  eat("RESULTS FOR MODEL\\: (.+)") -> model_description
  if (grepl("Linear Trend", model_description)) {
    z$model <- 2
  } else if (grepl("Effects for Each Time Point", model_description)) {
    z$model <- 3
  } else
    stop("Unknown model")
  eat("-+")
  eat("ESTIMATION METHOD = (.+)") -> method
  eat("Total time used: .* seconds")

  optional("Estimated Overdispersion = %f") -> sig2
  optional("Estimated Serial Correlation = %f") -> rho

  eat("GOODNESS OF FIT")
  z$gof = list()
  eat("Chi-square  %f, df %d, p %f") -> z$gof$chi2
  eat("Likelihood Ratio  %f, df  %d, p %f") -> z$gof$LR
  eat("AIC \\(up to a constant\\)  %f") -> z$gof$AIC

  eat("WALD.*")
  eat("Wald-Test  %f, df %d, p  %f") -> z$wald

  eat("PARAMETER ESTIMATES")
  if (z$model==2) {
    eat.table() -> z$coef
  } else if (z$model==3) {
    eat("Parameters for Each Time Point")
    coef = eat.table()
    coef[1,4] = coef[1,3]; coef[1,c(3,5)] <- 0 # fix empty spaces
    z$coef = coef

    eat("Linear Trend \\+ Deviations for Each Time")
    z$linear = list()
    eat.table() -> z$linear$trend
    eat("Time Deviations")
    eat.table(header=FALSE) -> z$linear$deviations
  }
  eat("Time INDICES")
  eat.table() -> z$time.idx
  r = which(is.na(z$time.idx[[5]])) # problem rows
  z$time.idx[r,2:5] <- c(1, 0, 1, 0)

  eat("TIME TOTALS")
  eat.table() -> z$time.tot

  z$overall = list()
  z$overall$mod = list()
  eat("OVERALL SLOPE MODEL: %s \\(p<%f\\) %s") -> z$overall$mod$txt
  eat.table() ->  z$overall$mod$par

  z$overall$imp = list()
  eat("OVERALL SLOPE IMPUTED\\+\\(recommended\\): %s \\(p<%f\\) %s") -> z$overall$imp$txt
  eat.table() ->  z$overall$imp$par

  z
}