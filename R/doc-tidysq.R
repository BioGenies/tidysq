#' tidysq: tidy analysis of biological sequences
#'
#' @description The \code{tidysq} package is a toolbox for the analysis of 
#' biological sequences in a tidy way.
#' @author Michal Burdukiewicz, Dominik Rafacz, Mateusz Bąkała, Leon Eyrich Jessen
#' @docType package
#' @importFrom stats rnorm setNames
#' @importFrom tibble as_tibble tibble
#' @import checkmate
#' @import vctrs
#' @name tidysq-package
#' @aliases tidysq
NULL

#' @useDynLib tidysq, .registration=TRUE
#' @importFrom Rcpp sourceCpp
NULL