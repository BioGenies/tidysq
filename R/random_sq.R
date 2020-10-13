#' Generate random sequences
#' 
#' Generates a \code{\link{sq}} object with specified number of sequences of given length 
#' and given type.
#' 
#' @param n a positive \code{\link{integer}} value - number of sequences to generate.
#' @param len a positive \code{\link{integer}} value - length of each sequence if \code{sd} not 
#' specified and mean length of sequences if \code{sd} specified
#' @param type a type of generated sq object; possible values are "ami", "dna" and "rna"
#' (see section \emph{sq types} in \code{\link{sq}} documentation for details).
#' @param is_clean a \code{\link{logical}} value - if \code{TRUE}, letters will be drawn from
#' clean alphabet, if \code{FALSE} - ambiguous letters might be also generated.
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
#' @examples 
#' # setting seed for reproducibility
#' set.seed(16)
#' 
#' # generating random sequences
#' random_sq(10, 10, "ami", TRUE)
#' random_sq(25, 18, "rna", TRUE, sd = 6)
#' random_sq(50, 8, "dna", FALSE, sd = 3)
#' random_sq(6, 100, "ami", TRUE, use_gap = TRUE)
#' @seealso \code{\link{construct_sq}} \code{\link{sq}}
#' @importFrom stringi stri_rand_strings stri_paste
#' @export
random_sq <- function(n, len, type, sd = NULL, use_gap = FALSE) {
  assert_count(n)
  assert_count(len)
  assert_sq_type(type)
  assert_number(sd, null.ok = TRUE)
  assert_flag(use_gap)
  
  alph <- get_standard_alph(type)
  if (!use_gap) alph <- .skip_characters(alph, "-")
  if (type == "ami") alph <- .skip_characters(alph, "*")
  
  alph_regex <- stri_paste("[", stri_paste(alph, collapse = ""), "]")
  if (!is.null(sd)) {
    # TODO: consider using other distribution than normal maybe?
    len <- round(rnorm(n, len, sd))
    len <- ifelse(len <= 0, 1, len)
  }
  
  sq <- stri_rand_strings(n, len, alph_regex)
  sq(sq, type)
}
