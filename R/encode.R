#' Encode sequences as numeric values
#' 
#' @description Encode sequences, stored in \code{\link{sq}} object, with 
#' numeric values defined by user in named vector.
#' 
#' @inheritParams reverse
#' @param encoding \code{encoding} a named \code{\link{numeric}} vector, that consists of 
#' values assigned to nucleic or amino acid letters. Letters without assigned encoding will 
#' be shown as \code{\link[=sq]{NA}}. Names of vector should be unique and should be elements
#' of alphabet of \code{sq} object.
#' 
#' @return \code{\link{sq}} object of type\strong{enc}.
#' 
#' @details Each position in a sequence is replaced by a numeric value 
#' assigned to that letter. 
#' 
#' Sometimes for research purposes one wants to replace letters by values, 
#' described by various properties of nucleic or amino acids. It can be a 
#' sparse encoding, hydrophobicity scale, heat capacities, entropies, chemical 
#' shift index, probability matrix (BLOSUM, PAM), 
#' sequence profile or reduced alphabet.
#' 
#' The newly constructed sequence will have a type \strong{enc} (see details on types 
#' in \code{\link{sq}}), which represents encoded sequences.
#' 
#' The named vector (ex. \code{c(G = 1, K = 2, P = 2)}) should have all letters 
#' assigned, otherwise unassigned letters will be shown as \code{NA}. If any letter that
#' appears in an alphabet appears in at least one of sequences, user will be informed about it.
#' Default action is a warning printed in the console, but it can be changed via setting
#' "tidysq_a_no_given_enc" (see details at \code{\link{tidysq-options}}).
#' 
#' In fact the only thing that is replaced is an alphabet - letters are substituted
#' with values assigned to them. The internal structure of the object remains unchanged.
#' 
#' If one wants to access numeric values of encoded sequences, they may use 
#' \code{\link{as.matrix}} or \code{\link{encsq_to_list}}.
#' 
#' @examples 
#' 
#' # Create sq object with sequences containing letters from 
#' # standard alphabet to work on:
#' 
#' sq_dna <- construct_sq(c("TATGAATTAGCTGTCTTTGCTGCTTTGGTTATCTATGA", 
#'                          "CTTTGGTTATCTAGCTGTATGA", "TATCTAGCTGTATG", 
#'                          "CTGCTG", "CTTAGA", "CCCT", "CTGAATGT"), 
#'                        type = "dna")
#' 
#' sq_ami <- construct_sq(c("NYMITGGREEYERTVIYRAIALNAANYTWIL", 
#'                         "TIAALGNIIYRAIE", "NYERTGHLI", 
#'                         "MAYNNNIALN", "MN", "NAAAT"), 
#'                         type = "ami")
#' 
#' 
#' # Create an object, called sq, with sequences containing letters from 
#' # standard and extended alphabet to work on:
#' 
#' sq_rna_ex <- construct_sq(c("VAHCHAGDUGBBVG", "CUGCVB", "DUUBGA", "CCCU", 
#'                             "CUHAABBU"), type = "rna")
#' 
#' sq_ami_ex <- construct_sq(c("MAYUOUONNNIALN", "UUMXBZONO", "NAAGAT"), type = "ami")
#' 
#' 
#' # Create encoding for a standard alphabet 
#' 
#' enc_dna <- c(A = 1, C = 2, G = 2, T = 2)
#' enc_rna <- c(A = 7, C = 5, G = 3, U = 2)
#' enc_ami <- c(A = 5, C = 5, D = 6, E = 6, F = 4, 
#'              G = 1, H = 5, I = 3, K = 2, L = 3, 
#'              M = 5, N = 6, P = 2, Q = 6, R = 2, 
#'              S = 6, T = 6, V = 3, W = 4, Y = 4)
#'                         
#'              
#' # Create encoding for an extended alphabet
#' 
#' enc_dna_ex <- c(A = 1, C = 1, G = 2, T = 3,
#'                 B = 5, D = 8, H = 13, V = 21)
#' enc_rna_ex <- c(A = 1, C = 1, G = 2, U = 3,
#'                 B = 5, D = 8, H = 13, V = 21)
#' enc_ami_ex <- c(A = 5, C = 5, D = 6, E = 6, F = 4, 
#'                 G = 1, H = 5, I = 3, K = 2, L = 3, 
#'                 M = 5, N = 6, P = 2, Q = 6, R = 2, 
#'                 S = 6, T = 6, V = 3, W = 4, Y = 4, 
#'                 U = -0.1, O = 0.753, X = -53.95, B = 7.77, Z = 0)
#'
#' 
#' # Encode with simplified standard alphabet sequences without 
#' # non-standard letters
#' 
#' encode(sq_dna, enc_dna)
#' encode(sq_ami, enc_ami)
#' 
#' 
#' # Encode with simplified extended alphabet sequences without 
#' # non-standard letters
#' 
#' encode(sq_dna, enc_dna_ex)
#' encode(sq_ami, enc_ami_ex)
#' 
#' 
#' # Encode with simplified standard alphabet sequences with 
#' # non-standard letters
#' 
#' encode(sq_rna_ex, enc_rna)
#' encode(sq_ami_ex, enc_ami)
#' 
#' 
#' # Encode with simplified extended alphabet sequences with 
#' # non-standard letters
#' 
#' encode(sq_rna_ex, enc_rna_ex)
#' encode(sq_ami_ex, enc_ami_ex)
#' 
#' 
#' # Encode without assigning all letters
#' 
#' encode(sq_dna, c(A = 1, G = 0.02))
#' encode(sq_ami, c(A = 5, H = 5, I = 0.3, K = -2, L = -3.1, M = 5, N = 6))
#' 
#' 
#' # Import encoding from \code{AAindex database}
#' 
#' data("AAindex_norm")
#' enc_aa <- AAindex_norm[20, ]
#' 
#' encode(sq_ami, enc_aa)
#' 
#' 
#' @seealso \code{\link{sq}} \code{\link{as.matrix}} or \code{\link{encsq_to_list}}
#' 
#' @export
encode <- function(sq, encoding) {
  .validate_sq(sq)
  .check_isnt_missing(encoding, "'encoding'")
  .check_is_named(encoding, "'encoding'")
  .check_numeric(encoding, "'encoding'", allow_zero = TRUE, allow_negative = TRUE, 
                 allow_na = TRUE, allow_nan = TRUE, allow_inf = TRUE)
  .check_is_unique(names(encoding), "'encoding'")
  
  type <- .get_sq_type(sq)
  if (type %in% c("ami", "dna", "rna"))
    names(encoding) <- toupper(names(encoding))
  
  alph <- alphabet(sq)
  alph_size <- .get_alph_size(alph)
  is_given <- alph %in% names(encoding)
  if (!all(is_given)) {
    ind <- (1:length(alph))[!is_given]
    for (s in sq) {
      if (any(C_unpack_ints(s, alph_size) %in% ind)) {
        .handle_opt_txt("tidysq_a_no_given_enc",
                        "there are letters in the alphabet of 'sq' that appear in sequences, but were not given in 'encoding' - assuming NA")
        break
      }
    }
    non_given <- rep(NA_real_, length(ind))
    names(non_given) <- alph[ind]
    encoding <- c(encoding, non_given)
  }
  
  new_list_of(sq,
              ptype = raw(),
              alphabet = encoding[alph],
              class = c("encsq", "sq"))
}
