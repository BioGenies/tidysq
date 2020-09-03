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
  options()[c("tidysq_a_bite_na",
              "tidysq_a_cln_sub_letters",
              "tidysq_a_no_given_enc",
              "tidysq_a_typify_small_cap_let",
              "tidysq_g_fast_mode",
              "tidysq_p_max_pillar_width",
              "tidysq_p_max_sequences",
              "tidysq_p_na_char",
              "tidysq_p_use_color")]
}