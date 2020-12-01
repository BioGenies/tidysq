#' Generate random sequences
#'
#' @templateVar alph_null_ok FALSE
#' 
#' Generates a \code{\link[=sq-class]{sq}} object with specified number of sequences of given length
#' and given type.
#' 
#' @param n a positive \code{\link{integer}} value - number of sequences to generate.
#' @param len a positive \code{\link{integer}} value - length of each sequence if \code{sd} not 
#' specified and mean length of sequences if \code{sd} specified
#' @template alphabet
#' @param sd a positive \code{\link{numeric}} value; if specified, gives standard deviation of
#' length of generated sequences.
#' @param use_gap - a \code{\link{logical}} value; if \code{TRUE}, sequences will be generated
#' with random gaps inside.
#' @return An object of class \code{sq} with type as specified.
#' 
#' Sequences are generated using \code{\link{sample}} function. There is no possibility of 
#' generating a sequence of length 0, even if \code{sd} is given. Letter '*' is not used 
#' in generating \strong{ami} sequences.
#' 
#' @seealso \code{\link{construct_sq}} \code{\link[=sq-class]{sq}}
#' @importFrom stringi stri_rand_strings stri_paste
#' @export
random_sq <- function(n, len, alphabet, sd = NULL, use_gap = FALSE) {
  assert_count(n)
  assert_count(len)
  assert_character(alphabet, any.missing = FALSE, min.len = 0, unique = TRUE, null.ok = TRUE)
  assert_number(sd, null.ok = TRUE)
  assert_flag(use_gap)
  
  if (length(alphabet) == 1) {
    type <- interpret_type(alphabet)
    if (type == "unt")
      stop("method 'random_sq' cannot take 'unt' as alphabet type", call. = FALSE)
    else
      alphabet <- get_standard_alphabet(type)
  } else {
    alphabet <- sq_alphabet(alphabet, "atp")
  }
  
  if (!is.null(sd)) {
    # TODO: consider using other distribution than normal maybe?
    len <- round(rnorm(n, len, sd))
    len <- ifelse(len <= 0, 1, len)
  }
  
  CPP_random_sq(n, len, alphabet, use_gap)
}
