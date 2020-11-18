#' Test if a sequence is empty
#' 
#' Test a sq object for presence of empty sequences
#' 
#' @param sq a\code{\link{sq}} object to be tested
#'  
#' @return A logical vector of the same length as input sq, 
#' indicating which elements are \code{NULL sq}, i.e., an empty sequence
#' of length 0.
#' 
#' @details This function allows identification of empty sequences (that
#' have length 0) represented by the \code{NULL sq} values in the sq object. 
#' It returns a logical for every element of the sq object - \code{TRUE} if
#' its value is \code{NULL sq} and \code{FALSE} otherwise. \code{NULL sq} 
#' values may be introduced as a result of \code{\link{clean}} and 
#' \code{\link{remove_na}} functions. The first one replaces sequences 
#' containing ambiguous elements with \code{NULL sq} values, whereas the 
#' latter replaces \code{NA} values (which may be introduced by 
#' \code{\link{bite}}) with \code{NULL sq}.
#'
#' @seealso \code{\link{sq}} \code{\link{clean}}
#' @export
is_null_sq <- function(x) {
  assert_class(x, "sq")
  get_sq_lengths(x) == 0
}