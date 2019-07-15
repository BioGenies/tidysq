#' Reverse sequence
#' 
#' Reverse given list of sequences.
#' 
#' @param sq \code{\link{sq}} object
#' 
#' @return \code{\link{sq}} object of the same type as input sq, where each element is reversed.
#' 
#' @examples 
#' # Reverse just one sequence
#' reverse(construct_sq("ACTAGAGTGATAGA", type = "nuc"))
#' 
#' # Reverse list of sequences
#' reverse(construct_sq(c("fafasfasfFSA", "ygagayagfa", "adsDaf"), type = "ami"))
#' 
#' # Reverse list of uncleaned sequences
#' reverse(construct_sq(c("PASJIFEHF", "hvfisxxx", "xxxer")))
#' 
#' # Reverse list of cleaned sequences - only_elements = FALSE
#' reverse(clean(construct_sq(c("fafasfasfFSA", "ygagayagfa", "adxxaf", "xxx"), type = "ami")))
#'
#' # Reverse list of cleaned sequences - only_elements = TRUE
#' reverse(clean(construct_sq(c("fafasfasfFSA", "ygagayagfa", "adxxaf", "xxx"), type = "ami"), only_elements = TRUE))
#' 
#' @export
reverse <- function(sq) {
  validate_sq(sq)
  ret <- lapply(sq, rev)
  .set_class_alph(ret, sq)
}