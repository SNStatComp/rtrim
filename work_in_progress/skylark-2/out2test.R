#rm(list=ls())
library(stringr)

#printf <- function(...) cat(sprintf(...))

read_out = function(filename, tcf, debug=FALSE) {
  lines = trimws(readLines(filename))

  DBG=debug

  eat.white = function() {
    # remove empty lines
    while (lines[1]=="") lines <<- lines[-1]
  }

  eat <- function(pat, debug=DBG) {
    if (debug) printf("\n\nEating: %s\n", pat)

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
    trimws(out)
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
      return(NA)
    }

    # first assert that we do have the pattern
    stopifnot(grepl(pat, curline))

    out=curline
    # extract requested groups, only for non-empty patterns
    if (pat!="^$") {
      m = str_match(curline, pat)
      if (debug) print(m)
      stopifnot(is.na(m[1,1])==FALSE)
      if (ncol(m)>1) out = m[1,2:ncol(m)]
    }
    if (debug) printf("output: %s\n", out)
    trimws(out)
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

  eat.list <- function() {
    # Eat a single-line list of integer numbers, eg "1 1 2 3 5 8"
    while (lines[1]=="") lines <<- lines[-1] # remove empty lines
    lst <- scan(textConnection(lines[1]), what=integer())
    lines <<- lines[-1]
    lst
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
  if (tcf$weight) eat("\\d+\\. weight")
  z$dat.covar.name <- character(tcf$ncovars) # Covariate names as in the data
  z$dat.covar.size <- integer(tcf$ncovars)
  for (i in 1:tcf$ncovars) {
    eat("\\d+\\. %s  number of values: %d") -> tmp
    z$dat.covar.name[i] = tmp[1]
    z$dat.covar.size[i] = as.integer(tmp[2])
  }

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

  if (tcf$model==1) {
    stop("Model 1 not supported")
  } else if (tcf$model==2) {
    eat("RESULTS FOR MODEL\\: Linear Trend")
    z$model <- 2
  } else if (tcf$model==3) {
    eat("RESULTS FOR MODEL\\: Effects for Each Time Point")
    z$model <- 3
  } else stop("Can't happen")
  eat("-+")

  # Actual covariates
  ncovar <- length(tcf$covariates)
  if (ncovar>0) {
    eat("Effects of covariate\\(s\\)")
    idx = integer(ncovar)
    for (i in 1:ncovar) {
      eat("%s") -> name
      stopifnot(name %in% z$dat.covar.name)
      idx[i] = which(name == z$dat.covar.name)
    }
    z$mod.covar.name = z$dat.covar.name[idx]
    z$mod.covar.size = z$dat.covar.size[idx]
  }

  if (length(tcf$changepoints)>0) {
    eat("Changes in Slope at Timepoints")
    eat.list() -> changepoints
    use.changepoints <- TRUE
    num.changepoints <- length(changepoints)
  } else {
    use.changepoints <- FALSE
  }

  eat("ESTIMATION METHOD = (.+)") -> method
  eat("Total time used: .* seconds")

  optional("Estimated Overdispersion = %f") -> sig2
  optional("Estimated Serial Correlation = %f") -> rho

  eat("GOODNESS OF FIT")
  z$gof = list()
  eat("Chi-square  %f, df %d, p %f") -> z$gof$chi2
  eat("Likelihood Ratio  %f, df  %d, p %f") -> z$gof$LR
  eat("AIC \\(up to a constant\\)  %f") -> z$gof$AIC

  z$wald = list()
  if (ncovar>0) {
    eat("WALD-TEST FOR SIGNIFICANCE OF COVARIATES")
    eat.table() -> z$wald$covar
  }

  if (z$model == 2) {
    if (use.changepoints) {
      eat("WALD-TEST FOR SIGNIFICANCE OF CHANGES IN SLOPE")
      eat.table() -> tmp
      tmp <- tmp[-1] # remove "Changepoints" column
      names(tmp) <- c("W","df","p")
      z$wald$dslope <- tmp
    } else {
      eat("WALD-TEST FOR SIGNIFICANCE OF SLOPE PARAMETER")
      eat("Wald-Test  %f, df %d, p  %f") -> tmp
      z$wald$slope <- as.list(as.numeric(tmp))
      names(z$wald$slope) <- c("W","df","p")
    }
  }

  if (z$model==3 && ncovar==0) {
    eat("WALD-TEST FOR SIGNIFICANCE OF DEVIATIONS FROM LINEAR TREND")
    eat("Wald-Test  %f, df %d, p  %f") -> tmp
    z$wald$deviations <- as.list(as.numeric(tmp))
    names(z$wald$deviations) <- c("W","df","p")
  }

  # # Use only actual covariates
  # covar.names = c(z$covar1[1], z$covar2[1])
  # covar.nclass = as.numeric(c(z$covar1[2], z$covar2[2]))
  #


  eat("PARAMETER ESTIMATES")
  if (z$model==2) {
    if (ncovar>0) {
      eat("Effects of Covariates.*")
      ncp  = num.changepoints
      ncat = sum(z$mod.covar.size - 1)
      npar = ncp * (ncat+1)
      # print(z$dat.covar.name)
      # print(z$dat.covar.size)
      # print(z$mod.covar.name)
      # print(z$mod.covar.size)
      # printf("ncp=%d ncat=%d npar=%d\n", ncp, ncat, npar)
      # stop("intended")
      coefs = matrix(0, npar, 8)
      r=1
      for (j in 1:ncp) {
        eat("from upto")
        eat("%d %d") -> fromto
        eat("Additive      std.err.   Multiplicative  std.err.")
        eat("Constant %f %f %f %f") -> pars
        coefs[r,1:2] <- c(0,0)
        coefs[r,3:4] <- as.numeric(fromto)
        coefs[r,5:8] <- as.numeric(pars)
        r <- r+1
        for (cv in 1:ncovar) {
          eat("Covariate %d") -> tmp
          coefs[r,1] <- as.numeric(tmp)
          eat("--*")
          for (k in 2:z$dat.covar.size[cv]) {
            eat("Category %d %f %f %f %f") -> tmp
            coefs[r,2] <- as.numeric(tmp[1]) # category
            coefs[r,3:4] <- as.numeric(fromto)
            coefs[r,5:8] <- as.numeric(tmp[2:5]) # params
            r <- r+1
          }
        }
      }
      # sort to align with R-TRIM output
      idx = order(coefs[,1], coefs[,2], coefs[,3])
      coefs <- coefs[idx, ]
    } else {
      if (use.changepoints) eat("Slope for Time Intervals")
      eat.table() -> coefs
      if (use.changepoints==FALSE) {
        prefix = data.frame(from=1, upto=as.numeric(z$ntime)) # prevent ntime->factor
        coefs = cbind(prefix, coefs)
      }
    }
  } else if (z$model==3) {
    if (ncovar==0) {
      eat("Parameters for Each Time Point")
      coefs = eat.table()
      coefs[1,4] = coefs[1,3]; coefs[1,c(3,5)] <- 0 # fix empty spaces

      eat("Linear Trend \\+ Deviations for Each Time")
      z$linear = list()
      eat.table() -> z$linear$trend
      eat("Time Deviations")
      eat.table(header=FALSE) -> z$linear$dev
    } else {
      eat("Additive std\\.err\\. Multiplicative std\\.err\\.")
      col_names <- c("time","add","se_add","mul","se_mul")
      eat("Constant")
      eat("--+")
      n = tcf$ntime
      # First row has a prefix "Time" which has to be removed
      df1 <- read.table(textConnection(lines[1]),   header=FALSE)[,-1] # skip prefix
      df2 <- read.table(textConnection(lines[2:n]), header=FALSE)
      lines <- lines[-(1:n)]
      names(df1) <- col_names
      names(df2) <- col_names
      coefs = cbind(data.frame(covar="baseline",cat=0), rbind(df1,df2))
      # Covariate parameters
      for (cv in 1:ncovar) {
        cv.name = z$mod.covar.name[cv]
        eat("Covariate %d") -> tmp
        stopifnot(tmp==cv) # reported convariates should be 1, 2, ...
        eat("--+")
        for (i in 2:z$mod.covar.size[cv]) {
          eat("Category %d") -> tmp
          stopifnot(tmp==i) # reported catagories should be 1, 2, ...
          cat.name <- sprintf("cat_%s", tmp)
          # Again read coefficients in two steps
          df1 <- read.table(textConnection(lines[1]), header=FALSE)[,-1]
          df2 <- read.table(textConnection(lines[2:n]), header=FALSE)
          lines <- lines[-(1:n)]
          names(df1) <- col_names
          names(df2) <- col_names
          prefix = data.frame(covar=cv.name, cat=i)
          coefs = rbind(coefs, cbind(prefix, rbind(df1,df2)))
        }
      }
    }
  }
  z$coefficients = coefs

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