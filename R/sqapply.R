#' Apply function to each sequence
#' 
#' Applies given function to each sequence. Sequences are passed to function as character vectors
#' (or numeric, if type of \code{sq} is \strong{enc}) or single character strings, depending on 
#' parameter.
#' 
#' @inheritParams reverse
#' @param fun a \code{\link{function}} to apply to each sequence in \code{sq} object; it should
#' take a character vector, numeric vector or single character string as an input.
#' @param ... another arguments passed to \code{fun}.
#' @param paste_char a \code{\link{logical}} value indicating in which form sequences should be
#' passed to the function \code{fun}; if \code{FALSE} (default), they will be treated as character
#' vectors, if \code{TRUE}, they will be pasted into a single string.
#' @param use_na_char a \code{\link{logical}} value indicating whether to use a printing character
#' to represent \code{\link{NA}} values; if \code{TRUE}, letter from option "tidysq_p_na_char"
#' will be used instead of \code{NA} values (default value for this option is "!", for details
#' see \code{\link{tidysq-options}}), otherwise just \code{NA} values will be used; default value
#' for this parameter is equal to \code{paste_char} value; \code{use_na_char} cannot be 
#' \code{FALSE} if \code{paste_char} is \code{TRUE}.
#' 
#' @return A list of values returned by function for each sequence in corresponding order.
#' 
#' @examples 
#' sq_ami <- construct_sq(c("YCYWIFTSRIK", "GGWGDVKCG", "KTKHIEQKL"))
#' 
#' # count how many times "K" appears in each sequence:
#' sqapply(sq_ami, function(sequence) sum(sequence == "K"))
#' 
#' # duplicate each element of sequence, then paste it and construct new sq
#' duplicated <- sqapply(sq_ami, function(sequence) paste(rep(sequence, each = 2), collapse = ""))
#' construct_sq(unlist(duplicated))
#' 
#' @seealso \code{\link{sq}} \code{\link[base]{lapply}}
#' @export 
sqapply <- function(sq, fun, ..., paste_char = FALSE, 
                    use_na_char = paste_char) {
  .validate_sq(sq)
  .check_logical(paste_char, "'paste_char'", single_elem = TRUE)
  .check_logical(use_na_char, "'use_na_char'", single_elem = TRUE)
  .check_paste_or_na(paste_char, use_na_char)
  
  na_char <- na_character(alphabet(sq))
  type <- .get_sq_type(sq)
  .apply_sq(sq, if (paste_char) "string" else "char", "none", function(s) {
    if (!use_na_char) s[s == na_char] <- NA
    if (type == "enc") s <- as.numeric(s)
    fun(s, ...)
  })
}