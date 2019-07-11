#' Reverse sequence
#' 
#' Reverse given list of sequences. 
#' 
#' @param sq sq object
#' 
#' @return reversed sq object (with preserved type)
#' 
#' @examples 
#' reverse(construct_sq("ACTAGAGTGATAGA", type = "nuc"))
#' reverse(construct_sq(c("fafasfasfFSA", "ygagayagfa", "adsDaf"), type = "ami"))
#'
#' @export
reverse <- function(sq) {
  validate_sq(sq)
  ret <- lapply(sq, rev)
  .set_class_alph(ret, sq)
}