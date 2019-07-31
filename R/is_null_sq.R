#' Test if a sequence is empty
#' 
#' Test a sq object for presence of empty sequences
#' 
#' @param sq \code{\link{sq}} object to be tested
#'  
#' @return a logical vector of the same length as input sq, 
#' indicating which elements are \code{NULL sq}, i.e., an empty sequence.
#' 
#' @details This function allows identification of empty sequences 
#' represented by the \code{NULL sq} values in the sq object. It returns
#' a logical for every element of the sq object - \code{TRUE} if
#' its value is \code{NULL sq} and \code{FALSE} otherwise. 
#' \code{NULL sq} values may be introduced as a result of 
#' \code{\link{clean}} function in place of sequences containing
#' ambiguous elements. 
#'
#' @examples 
#' # Creating an object to work on:
#' sq <- construct_sq(c("ACGATTAGACG", "", "GACGANTCCAGNTAC"), type = "nuc")
#' 
#' # Testing for presence of empty sequences:
#' is_null_sq(sq)
#' 
#' # Testing for presence of empty sequences after cleaning - sequence 
#' # containing ambiguous elements is replaced by NULL sq:
#' cln_sq <- clean(sq)
#' is_null_sq(cln_sq)
#' 
#' @seealso sq clean clnsq
#' @export

is_null_sq <- function(sq) {
  validate_sq(sq)
  sapply(sq, function(s) identical(s, as.raw(0)))
}