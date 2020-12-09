#' Apply function to each sequence
#' 
#' Applies given function to each sequence. Sequences are passed to function as character vectors
#' (or numeric, if type of \code{sq} is \strong{enc}) or single character strings, depending on 
#' parameter.
#' 
#' @template x
#' @param fun [\code{function(1)}]\cr
#'  A function to apply to each sequence in \code{sq} object; it should
#'  take a character vector, numeric vector or single character string as an input.
#' @template three-dots
#' @param paste_char [\code{logical(1)}]\cr
#'  A value indicating in which form sequences should be
#'  passed to the function \code{fun}; if \code{FALSE} (default), they will be treated as character
#'  vectors, if \code{TRUE}, they will be pasted into a single string.
#' @param use_na_letter [\code{logical(1)}]\cr
#'  A value indicating whether to use a printing character
#'  to represent \code{\link{NA}} values; if \code{TRUE}, letter from option "tidysq_p_na_letter"
#'  will be used instead of \code{NA} values (default value for this option is "!", for details
#'  see \code{\link{tidysq-options}}), otherwise just \code{NA} values will be used; default value
#'  for this parameter is equal to \code{paste_char} value; \code{use_na_letter} cannot be
#'  \code{FALSE} if \code{paste_char} is \code{TRUE}.
#' 
#' @return A list of values returned by function for each sequence in corresponding order.
#' 
#' @seealso \code{\link[=sq-class]{sq}} \code{\link[base]{lapply}}
#' @export 
#' @export
sqapply <- function(x, fun, ...,
                    single_string = FALSE, 
                    NA_letter = getOption("tidysq_NA_letter")) {
  assert_class(x, "sq")
  assert_function(fun)
  assert_flag(single_string)
  assert_string(NA_letter, min.chars = 1)
  
  CPP_apply_R_function(x, 
                       function(sequence) fun(sequence, ...), 
                       single_string, NA_letter)
}