# ######################################################### Goodness of fit ####

# The goodness-of-fit of the model is assessed using three statistics:
# Chi-squared, Likelihood Ratio and Aikaike Information Content.

# ============================================================= Computation ====

# Here we define `gof' as a S3 generic function
#' Extract TRIM goodness-of-fit information.
#'
#' \code{\link{trim}} computes three goodness-of-fit measures:
#' \itemize{
#'   \item Chi-squared
#'   \item Likelihood ratio
#'   \item Akaike Information content
#' }
#'
#' @param x an object of class \code{\link{trim}}. 
#'
#' @return a list of type "trim.gof", containing elements \code{chi2}, \code{LR}
#' and \code{AIC}, for Chi-squared, Likelihoof Ratio and Akaike informatiuon content,
#' respectively.
#' @export
#' 
#' @family analyses
#' @examples 
#' data(skylark)
#' z <- trim(count ~ time + site, data=skylark, model=2)
#' # prettyprint GOF information
#' gof(z)
#' 
#' # get individual elements, e.g. p-value
#' L <- gof(z)
#' LR_p <- L$LR$p # get p-value for likelihood ratio
#' 
gof <- function(x) UseMethod("gof")

# Here is a simple wrapper function for TRIM output lists.
#' @export
#' @rdname gof
gof.trim <- function(x) {
  stopifnot(class(x)=="trim")
  gof.numeric(x$f, x$mu, x$alpha, x$beta)
}

# Here is the workhorse function

gof.numeric <- function(f, mu, alpha, beta) {
  observed <- is.finite(f)

  # The $\chi^2$ (Chi-square) statistic is given by
  # \begin{equation}
  #   \chi^2 = \sum_{ij}\frac{f_{i,j}-\mu_{i,j}}{\mu_{i,j}}
  # \end{equation}
  # where the summation is over the observed $i,j$'s only.
  # Significance is assessed by comparing against a $\chi^2$ distribution with
  # $df$ degrees of freedom, equal to the number of observations
  # minus the total number of parameters involved, i.e.\
  # $df = n_f - n_\alpha - n_\beta$.
  chi2 <- sum((f-mu)^2/mu, na.rm=TRUE)
  df   <- sum(observed) - length(alpha) - length(beta)
  p    <-  1 - pchisq(chi2, df=df)
  chi2 <- list(chi2=chi2, df=df, p=p) # store in a list

  # Similarly, the \emph{Likelihood ratio} (LR) is computed as
  # \begin{equation}
  #   \operatorname{LR} = 2\sum_{ij}f_{ij} \log\frac{f_{i,j}}{\mu_{i,j}} \label{LR}
  # \end{equation}
  # and again compared against a $\chi^2$ distribution.
  LR <- 2 * sum(f * log(f / mu), na.rm=TRUE)
  df <- sum(observed) - length(alpha) - length(beta)
  p  <- 1 - pchisq(LR, df=df)
  LR <- list(LR=LR, df=df, p=p)

  # The Akaike Information Content (AIC) is related to the LR as:
  AIC <- LR$LR - 2*LR$df

  # Output all goodness-of-fit measures in a single list
  structure(list(chi2=chi2, LR=LR, AIC=AIC), class="trim.gof")
}

# ================================================================ Printing ====

# A simple printing function is provided that mimics TRIM for Windows output.

#' Print method for \code{trim.gof}
#'
#' @export
#' @param x a \code{trim.gof} object
#' @keywords internal
print.trim.gof <- function(x,...) {
  # print welcome message
  cat(sprintf("Goodness of fit:\n"))

  # print $\chi^2$ results
  with(x$chi2,
       printf("%24s = %.2f, df=%d, p=%.4f\n", "Chi-square", chi2, df, p))

  # idem, Likelihood ratio
  with(x$LR,
       printf("%24s = %.2f, df=%d, p=%.4f\n", "Likelihood Ratio", LR, df, p))

  # idem, Akaike Information Content
  with(x,
       printf("%24s = %.2f\n", "AIC (up to a constant)", AIC))
}
