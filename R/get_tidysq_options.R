#' Obtain current state of tidysq options
#' 
#' @return a \code{\link[=namedList-class]{named list}} containing values of 
#' currently set options associated with the \pkg{tidysq} package. 
#' 
#' @details To see full list of options and their meanings, 
#' see \code{\link{tidysq-options}} man page.
#' 
#' @seealso \code{\link{tidysq-options}}
get_tidysq_options <- function() {
  options()[c("tidysq_NA_letter",
              "tidysq_on_warning",
              "tidysq_pillar_max_width",
              "tidysq_print_max_sequences",
              "tidysq_print_use_color",
              "tidysq_safe_mode")]
}