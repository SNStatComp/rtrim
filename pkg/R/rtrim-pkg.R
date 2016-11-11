#' Trend and Indices for Monitoring Data
#'
#' The TRIIM model is used to estimate species populations based on frequent 
#' (annual) counts at a varying collection of sites. The model is able to take 
#' account of missing data by imputing prior to estimation of population totals.
#' The current package is a complete re-implementation of the \code{Delphi}
#' based 
#' \href{https://www.cbs.nl/en-gb/society/nature-and-environment/indices-and-trends--trim--}{TRIM}
#' software by J. Pannekoek.
#'
#'
#' @section Getting started:
#'
#' Users of the original TRIM software can get started by following the \href{../doc/rtrim_for_TRIM_users.html}{rtrim for trim users} manual.
#'
#'
#' @name rtrim-package
#' @docType package
#' @import methods
#' @importFrom utils read.table head tail str capture.output
#' @importFrom grDevices  gray rgb hcl adjustcolor
#' @importFrom graphics lines plot points polygon segments title abline legend par
#' @importFrom stats pchisq pt qt time setNames
#'
{}

.onLoad <- function(libname, pkgname){
  options(trim_verbose=FALSE)
}

