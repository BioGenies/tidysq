#' Reverse sequence
#' 
#' @description Reverse given list of sequences.
#' 
#' @param x a \code{\link{sq}} object.
#' 
#' @details The \code{reverse} function reverses each sequence in supplied 
#' \code{\link{sq}} object (e.q. transforms "MIAANYTWIL" to "LIWTYNAAIM"). 
#' Empty sequences are left with no effect. This operation does not change 
#' the type of the input object nor its alphabet.
#' 
#' Since the function \code{reverse} returns a \code{\link{sq}} object, the 
#' \code{\link[=sq-print]{print}} function is implicitly called.
#' 
#' @return A \code{\link{sq}} object of the same type as input object but 
#' each sequence is reversed.
#' 
#' @seealso \code{\link{sq}} \code{\link{clean}} \code{\link{sq-print}}
#' 
#' @export
reverse <- function(x) {
  assert_class(x, "sq")
  ret <- .apply_sq(x, "int", "int", rev)
  vec_restore(ret, x)
}