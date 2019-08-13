#' Transform encoded sequence to a list
#' 
#' @description Transform encoded by \code{\link{encode}} sequence, stored in
#' \code{\link{sq}} object, to a list.
#' 
#' @param indices \code{encsq} a \code{\link{sq}} object which was encoded 
#' using \code{\link{encode}} function.
#' 
#' In a new object you can check what value is assigned to each letter.
#' 
#' @return A named lists with sequences and theirs encoding.
#' 
#' @details Function is used to transform an \code{\link{sq}} object with 
#' \code{\link{encsq}} class to a named list.
#' 
#' Each nucleic or amino acid sequence and assigned encodings.
#' 
#' @examples 
#' 
#' # Create sq object with sequences containing letters from 
#' stadard alphabet and extended alphabet to work on:
#' 
#' sq_nuc <- construct_sq(c("TATGAATTAGCTGTCTTTGCTGCTTTGGTTATCTATGA", 
#'                          "CTTTGGTTATCTAGCTGTATGA", "TATCTAGCTGTATG", 
#'                          "CTGCTG", "CTTAGA", "CCCT", "CTGAATGT"), 
#'                        type = "nuc")
#'                        
#' sq_ami <- construct_sq(c("NYMITGGREEYERTVIYRAIALNAANYTWIL", 
#'                          "TIAALGNIIYRAIE", "NYERTGHLI", "MAYNNNIALN", 
#'                          "MN", "NAAAT"), type = "ami")
#'                          
#' sq_nuc_ex <- construct_sq(c("TATCTAGCTGTATG", "CUGCUG", "CUUAGA", 
#'                             "CCCT", "CUGAAUGU"))
#'                             
#' sq_ami_ex <- construct_sq(c("MAYUOUONNNIALN", "UUMXBZONO", 
#'                             "NAAGAT"))
#' 
#' 
#' # Create encoding for standard, extended alphabet and iport from 
#' other sources:
#' 
#' enc_nuc <- c(A = 1, C = 2, G = 2, T = 2)
#' 
#' enc_ami <- c(A = 5, C = 5, D = 6, E = 6, F = 4, 
#'              G = 1, H = 5, I = 3, K = 2, L = 3, 
#'              M = 5, N = 6, P = 2, Q = 6, R = 2, 
#'              S = 6, T = 6, V = 3, W = 4, Y = 4)
#'              
#' enc_nuc_ex <- c(A = 1, C = 2, G = 2, T = 2, U = 3)
#' 
#' enc_ami_ex <- c(A = 5, C = 5, D = 6, E = 6, F = 4, 
#'                 G = 1, H = 5, I = 3, K = 2, L = 3, 
#'                 M = 5, N = 6, P = 2, Q = 6, R = 2, 
#'                 S = 6, T = 6, V = 3, W = 4, Y = 4, 
#'                 U = -0.1, O = 0.753, X = -53.95, 
#'                 B = 7.77, Z = 0)
#'                 
#' data("aaprop")
#' enc_aa <- aaprop[20,]
#' 
#' 
#' # Encode sequences and assign it to a variable:
#' 
#' e1 <- encode(sq_nuc, enc_nuc)
#' e2 <- encode(sq_ami, enc_ami)
#' e3 <- encode(sq_nuc, enc_nuc_ex)
#' e4 <- encode(sq_ami, enc_ami_ex)
#' e5 <- encode(sq_nuc_ex, enc_nuc)
#' e6 <- encode(sq_ami_ex, enc_ami)
#' e7 <- encode(sq_nuc_ex, enc_nuc_ex)
#' e8 <- encode(sq_ami_ex, enc_ami_ex)
#' e9 <- encode(sq_nuc, c(A = 1, G = 0.02))
#' e10 <- encode(sq_ami, c(A = 5, H = 5, I = 0.3, K = -2, 
#'                         L = -3.1, M = 5, N = 6))
#' e_aa <- encode(sq_ami, enc_aa)
#' 
#' 
#' # Transform \code{\link{encode}} result to a list:
#' 
#' ## Sequence with standard alphabet, encoding with standard alphabet
#' encsq_to_list(e1)
#' encsq_to_list(e2)
#' 
#' ## Sequence with standard alphabet, encoding with extended alphabet
#' encsq_to_list(e3)
#' encsq_to_list(e4)
#' 
#' ## Sequence with extended alphabet, encoding with standard alphabet
#' encsq_to_list(e5)
#' encsq_to_list(e6)
#' 
#' ## Sequence with extended alphabet, encoding with extended alphabet
#' encsq_to_list(e7)
#' encsq_to_list(e8)
#' 
#' ## Sequence with standard alphabet, encoding without all letters assigned
#' encsq_to_list(e9)
#' encsq_to_list(e10)
#' 
#' ## Import encoding from other sources
#' encsq_to_list(e_aa)
#' 
#' 
#' # What encoding is assigned to a letter
#' 
#' ee1 <- encsq_to_list(e1)
#' ee1[[1]]["T"]
#' 
#' @seealso sq encsq encode
#' 
#' @export
encsq_to_list <- function(encsq) {
  validate_sq(encsq, type = "enc")
  
  alph <- .get_alph(encsq)
  .apply_sq(encsq, "int", "none", function(s) alph[s])
}