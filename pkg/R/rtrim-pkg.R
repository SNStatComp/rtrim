#' Trend and Indices for Monitoring Data
#' 
#' Population totals are estimated based on ab
#'
#'
#'
#' @name rtrim-package
#' @docType package
#' @import methods
#' @importFrom utils read.table head tail str capture.output
#' @importFrom grDevices  gray rgb hcl
#' @importFrom graphics lines plot points polygon segments title
#' @importFrom stats pchisq pt qt time setNames
#'
{}

.onLoad <- function(libname, pkgname){
  options(trim_verbose=FALSE)
}

