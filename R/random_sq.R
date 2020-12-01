#' Generate random sequences
#'
#' @templateVar alph_null_ok FALSE
#' 
#' @description Generates an \code{\link[=sq-class]{sq}} object with specified
#' number of sequences of given length and alphabet.
#' 
#' @param n [\code{integer(1)}]\cr
#'  A number of sequences to generate - must be non-negative.
#' @param len [\code{integer(1)}]\cr
#'  Length of each sequence if \code{sd} not specified and mean length of
#'  sequences if \code{sd} specified - must be non-negative.
#' @template alphabet
#' @param sd [\code{integer(1)}]\cr
#'  If specified, gives standard deviation of length of generated sequences -
#'  must be non-negative.
#' @param use_gap [\code{logical(1)}]\cr
#'  If \code{TRUE}, sequences will be generated with random gaps inside
#'  (commonly denoted as "\code{-}").
#'
#' @return An object of class \code{sq} with type as specified.
#'
#' @details
#' Letter '*' is not used in generating \strong{ami} sequences. If parameter
#' \code{sd} is passed, then all generated negative values are replaced with 0s.
#'
#' @family io_functions
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
    len <- ifelse(len < 0, 0, len)
  }
  
  CPP_random_sq(n, len, alphabet, use_gap)
}
