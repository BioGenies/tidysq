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
#' \code{\link[=sq-print]{print}} function is implicitly called.
#' 
#' @return A \code{\link{sq}} object of the same type as input object but 
#' each sequence is reversed.
#' 
#' @examples 
#' # Creating sq objects using construct_sq:
#' sq_ami <- construct_sq(c("MIAANYTWIL","TIAALGNIIYRAIE", 
#'                          "NYERTGHLI", "MAYXXXIALN"), type = "ami")
#' sq_dna <- construct_sq(c("ATGCAGGA", "GACCGAACGAN", ""), type = "dna")
#' sq_unt <- construct_sq(c("ATGCAGGA!", "TGACGAGCTTA", "", "TIAALGNIIYRAIE"))
#' 
#' # Reverse sequences:
#' reverse(sq_ami)
#' reverse(sq_dna)
#' reverse(sq_unt)
#' 
#' # Reverse cleaned sequences:
#' reverse(clean(sq_ami))
#' reverse(clean(sq_dna))
#' 
#' @seealso \code{\link{sq}} \code{\link{clean}} \code{\link{print.sq}}
#' 
#' @export
reverse <- function(sq) {
  .validate_sq(sq)
  ret <- .apply_sq(sq, "int", "int", rev)
  vec_restore(ret, sq)
}