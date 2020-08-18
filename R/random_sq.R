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
#' @export
random_sq <- function(n, len, type, is_clean, sd = NULL, use_gap = FALSE) {
  .check_integer(n, "'n'", single_elem = TRUE)
  .check_integer(len, "'len'", single_elem = TRUE)
  .check_type(type)
  .check_logical(is_clean, "'is_clean'", single_elem = TRUE)
  .check_numeric(sd, "'sd'", allow_null = TRUE, single_elem = TRUE)
  .check_logical(use_gap, "'use_gap'", single_elem = TRUE)
  
  alph <- .get_standard_alph(type, is_clean)
  if (!use_gap) alph <- .skip_characters(alph, "-")
  if (type == "ami") alph <- .skip_characters(alph, "*")
  
  if (is.null(sd))
    sq <- sapply(1:n, function(i) paste0(sample(alph, len, replace = TRUE), collapse = ""))
  else {
    len <- round(rnorm(n, len, sd))
    len <- ifelse(len <= 0, 1, len)
    sq <- sapply(1:n, function(i) paste0(sample(alph, len[i], replace = TRUE), collapse = ""))
  }
  construct_sq(sq, type, is_clean)
}
