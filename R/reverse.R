#' Reverse sequence
#' 
#' @description Reverse given list of sequences.
#' 
#' @param sq a \code{\link{sq}} object.
#' 
#' @details The \code{reverse} function reverses each sequence in supplied 
#' \code{\link{sq}} object (e.q. transforms "MIAANYTWIL" to "LIWTYNAAIM"). 
#' Empty sequences are left with no effect. This operation does not change 
#' the type of the input object nor its alphabet.
#' 
#' Since the function \code{reverse} returns a \code{\link{sq}} object, the 
#' \code{\link{print.sq}} function is implicitly called.
#' 
#' @return A \code{\link{sq}} object of the same type as input object but 
#' each sequence is reversed.
#' 
#' @examples 
#' # Creating sq objects using construct_sq:
#' sq_ami <- construct_sq(c("MIAANYTWIL","TIAALGNIIYRAIE", 
#'                          "NYERTGHLI", "MAYXXXIALN"), type = "ami")
#' sq_nuc <- construct_sq(c("ATGCAGGA", "GACCGAACGAN", 
#'                          ""), type = "nuc")
#' sq_unt <- construct_sq(c("ATGCAGGA!", "TGACGAGCTTA", "", "TIAALGNIIYRAIE"))
#' 
#' # Reverse sequences:
#' reverse(sq_ami)
#' reverse(sq_nuc)
#' reverse(sq_unt)
#' 
#' # Reverse cleaned sequences:
#' reverse(clean(sq_ami))
#' reverse(clean(sq_nuc))
#' 
#' @seealso \code{\link{sq}} \code{\link{clean}} \code{\link{print.sq}}
#' 
#' @export
reverse <- function(sq) {
  validate_sq(sq)
  alph_size <- .get_alph_size(.get_alph(sq))
  ret <- .apply_sq(sq, "int", "int", rev)
  .set_class_alph(ret, sq)
}