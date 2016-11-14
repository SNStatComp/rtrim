#' @name skylark
#' @aliases skylark2
#' @title Skylark population data
#' 
#' @usage data(skylark); data(skylark2)
#' 
#' @description 
#' 
#' The Skylark dataset that was included with the original TRIM software.
#' 
#' The dataset can be loaded in two forms. The \bold{\code{skylark}} dataset is
#' exactly equal to the data set in the original TRIM software:
#' 
#' \tabular{lll}{
#'  \bold{Column}      \tab \bold{Type}    \tab \bold{Description}\cr
#'  \code{site}        \tab \code{integer} \tab Site number\cr
#'  \code{time}        \tab \code{integer} \tab Time point coded as integer sequence\cr
#'  \code{count}       \tab \code{numeric} \tab Counted skylarks\cr
#'  \code{Habitat}     \tab \code{integer} \tab Habitat type (1, 2)\cr
#'  \code{Deposition}  \tab \code{integer} \tab Deposition type (1, 2, 3, 4)
#' } 
#' 
#' The current implementation is more flexible and allows time points to be coded as years
#' and covariates as factors. The \bold{\code{skylark2}} data set looks as follows.
#' 
#' \tabular{lll}{
#'  \bold{Column}      \tab \bold{Type}    \tab \bold{Description}\cr
#'  \code{site}        \tab \code{factor}  \tab Site number\cr
#'  \code{year}        \tab \code{integer} \tab Time point coded as year\cr
#'  \code{count}       \tab \code{integer} \tab Counted skylarks\cr
#'  \code{Habitat}     \tab \code{factor}  \tab Habithat type (\code{dunes}, \code{heath})\cr
#'  \code{Deposition}  \tab \code{integer} \tab Deposition type (1, 2, 3, 4)\cr
#'  \code{Weight}      \tab \code{numeric} \tab Site weight
#' } 
#' 
#' 
#'
#' @docType data
#' @format \code{.RData}
#'
NULL
