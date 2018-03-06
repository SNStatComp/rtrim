#' Trend and Indices for Monitoring Data
#'
#' The TRIM model is used to estimate species populations based on frequent
#' (annual) counts at a varying collection of sites. The model is able to take
#' account of missing data by imputing prior to estimation of population totals.
#' The current package is a complete re-implementation of the \code{Delphi}
#' based
#' \href{https://www.cbs.nl/en-gb/society/nature-and-environment/indices-and-trends--trim--}{TRIM}
#' software developed at Statistics Netherlands by Jeroen Pannekoek.
#'
#' @section Getting started:
#'
#' Users of the original TRIM software can get started by following the
#' \href{../doc/rtrim_for_TRIM_users.html}{rtrim for trim users} or
#' following the \href{../doc/Skylark_example.html}{rtrim by example} manuals.
#'
#' Users of rtrim version 1.*.* may be interested in the many extensions, mainly described in the vignette
#' \href{../doc/rtrim_2_extensions.html}{rtrim 2 extensions}
#'
#' Users who would like to have more insight what is going on under the hood, may take a look in
#' \href{../doc/TRIM_methods_v2.pdf}{Models and statistical mehods in rtrim},
#' \href{../doc/TRIM_confidence_intervals.html}{rtrim confidence intervals}
#' and
#' \href{../doc/taming_overdispersion.html}{Taming overdispersion}.
#'
#' Enjoy!
#'
#' @name rtrim-package
#' @docType package
#' @import methods
#' @importFrom utils read.table head tail str capture.output
#' @importFrom grDevices  gray rgb hcl adjustcolor col2rgb
#' @importFrom graphics lines plot points polygon segments title abline legend par rasterImage rect
#' @importFrom stats pchisq qchisq pt qt qnorm time setNames quantile qgamma
#'
{}

.onLoad <- function(libname, pkgname){
  options(trim_verbose=FALSE)
}

