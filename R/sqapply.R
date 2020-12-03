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
sqapply <- function(x, fun, ..., paste_char = FALSE,
                    use_na_letter = paste_char) {
  assert_class(x, "sq")
  assert_flag(paste_char)
  assert_flag(use_na_letter)
  assert_false(paste_char && use_na_letter)
  
  na_letter <- getOption("tidysq_NA_letter")
  type <- sq_type(x)
  .apply_sq(x, if (paste_char) "string" else "char", "none", function(s) {
    if (!use_na_letter) s[s == na_letter] <- NA
    if (type == "enc") s <- as.numeric(s)
    fun(s, ...)
  })
}