#' tidysq: tidy analysis of biological sequences
#'
#' @description The \code{tidysq} package is a toolbox for the analysis of 
#' biological sequences in a tidy way.
#' @author Michal Burdukiewicz, Dominik Rafacz, Leon Eyrich Jessen
#' @docType package
#' @importFrom stats na.omit rnorm setNames
#' @importFrom utils download.file installed.packages
#' @importFrom tibble as_tibble tibble
#' @import vctrs
#' @name tidysq-package
#' @aliases tidysq
NULL

#' @useDynLib tidysq
#' @importFrom Rcpp sourceCpp
#' @importFrom RcppParallel RcppParallelLibs
NULL