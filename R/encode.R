#' Encode sequences
#' 
#' @description Function encodes sequences using numeric values defined by user in named vector.
#' 
#' @param sq \code{\link{sq}} object.
#' @param indices \code{encoding} a named vector, that consists of nucleic or amino acid letters with their numeric encoding.
#' Letters without assigned encoding will be shown as \code{NA}.
#' 
#' 
#' @return \code{\link{encsq}} object of the same type as input sq with encoded alphabet.
#' 
#' @details Each position in a sequence is replaced by a numeric value assigned to that letter. 
#' 
#' Sometimes for research purposes one wants to replace letters by various values, described by physio-chemical properties of nucleic or amino acids. 
#' It can be a hydrophobicity scale, heat capacities, entropies, chemical shift index or probability matrix (BLOSUM, PAM).
#' 
#' The newly constructed sequence will have a new class \code{\link{encsq})}, representing sequence encoded with custom alphabet.
#' 
#' The named vector (ex. \code{c(G = 1, K = 2, P = 2)}) should have all letters assigned, otherwise unasigned letters will be shown as \code{NA}.
#' 
#' All replaced letters will have the numeric type.
#' 
#' 
#' @examples 
#' 
#' # Create object, called sq, with sequences containing letters from stadard alphabet to work on:
#' 
#' sq_nuc <- construct_sq(c("TATGAATTAGCTGTCTTTGCTGCTTTGGTTATCTATGA", "CTTTGGTTATCTAGCTGTATGA", 
#'                         "TATCTAGCTGTATG", "CTGCTG", "CTTAGA", "CCCT", "CTGAATGT"), type = "nuc")
#' 
#' sq_ami <- construct_sq(c("NYMITGGREEYERTVIYRAIALNAANYTWIL", "TIAALGNIIYRAIE", "NYERTGHLI", 
#'                         "MAYNNNIALN", "MN", "NAAAT"), type = "ami")
#' 
#'                
#' # Create object, called sq, with sequences containing letters from stadard and extended alphabet to work on:
#'    
#' sq_nuc_ex <- construct_sq(c("TATCTAGCTGTATG", "CUGCUG", "CUUAGA", "CCCT", "CUGAAUGU"))
#' 
#' sq_ami_ex <- construct_sq(c("MAYUOUONNNIALN", "UUMXBZONO", "NAAGAT"))       
#' 
#'        
#' # Create encoding for standard alphabet 
#' 
#' enc_nuc <- c(A = 1, C = 2, G = 2, T = 2)
#' enc_ami <- c(A = 5, C = 5, D = 6, E = 6, F = 4, 
#'              G = 1, H = 5, I = 3, K = 2, L = 3, 
#'              M = 5, N = 6, P = 2, Q = 6, R = 2, 
#'              S = 6, T = 6, V = 3, W = 4, Y = 4)                   
#'                         
#'              
#' # Create encoding for extended alphabet  
#' 
#' enc_nuc_ex <- c(A = 1, C = 2, G = 2, T = 2, U = 3)
#' enc_ami_ex <- c(A = 5, C = 5, D = 6, E = 6, F = 4, 
#'              G = 1, H = 5, I = 3, K = 2, L = 3, 
#'              M = 5, N = 6, P = 2, Q = 6, R = 2, 
#'              S = 6, T = 6, V = 3, W = 4, Y = 4, 
#'              U = -0.1, O = 0.753, X = -53.95, B = 7.77, Z = 0)
#'                      
#'                                               
#' # Encode with simplified standard alphabet sequences without non-standard letters
#' 
#' encode(sq_nuc, enc_nuc)
#' encode(sq_ami, enc_ami)
#' 
#' 
#' # Encode with simplified extended alphabet sequences without non-standard letters
#' 
#' encode(sq_nuc, enc_nuc_ex)
#' encode(sq_ami, enc_ami_ex)
#' 
#' 
#' # Encode with simplified standard alphabet sequences with non-standard letters
#' 
#' encode(sq_nuc_ex, enc_nuc)
#' encode(sq_ami_ex, enc_ami)
#' 
#' 
#' # Encode with simplified extended alphabet sequences with non-standard letters
#' 
#' encode(sq_nuc_ex, enc_nuc_ex)
#' encode(sq_ami_ex, enc_ami_ex)
#' 
#' 
#' # Encode without assigning all letters
#' 
#' encode(sq_nuc, c(A = 1, G = 2))
#' encode(sq_ami, c(A = 5, H = 5, I = 3, K = 2, L = 3, M = 5, N = 6))
#' 
#' 
#' # Use created encoding from \code{AAindex database}
#' 
#' data("aaprop")
#' enc_aa <- aaprop[20,]
#' 
#' encode(sq_ami, enc_aa)
#' 
#' 
#' @seealso sq encsq
#' 
#' @exportClass encsq
#' @export
encode <- function(sq, encoding) {
  validate_sq(sq)
  type <- .get_sq_type(sq)
  
  if (!is.numeric(encoding)) {
    stop("'encoding' should be a numeric vector")
  }
  
  if (type %in% c("ami", "nuc")) {
    names(encoding) <- toupper(names(encoding))
    has_gap <- "-" %in% names(encoding)
    has_end <- "*" %in% names(encoding)
    
    if (!(has_gap && has_end)) {
      .handle_opt_txt("tidysq_encode_nogap_action",
                      "'encoding' don't contain values for '-' or '*' or both, assuming NA")
    }
    if (!has_gap) encoding <- c(encoding, `-` = NA)
    if (!has_end) encoding <- c(encoding, `*` = NA)
  }
  
  alph <- .get_alph(sq)
  if (!all(alph %in% names(encoding))) {
    if (length(encoding) == length(alph)) {
      names(encoding) <- alph
      .handle_opt_txt("tidysq_encode_noname_action",
                      "'encoding' is unnamed, assuming values corresponds to letters in order")
    } else stop("'encoding' should be a named vector which names are superset of alphabet of 'sq'")
  }
  
  sq <- .set_alph(sq, encoding[alph])
  .set_class(sq, "enc", FALSE)
}






sq_nuc <- construct_sq(c("TATGAATTAGCTGTCTTTGCTGCTTTGGTTATCTATGA", "CTTTGGTTATCTAGCTGTATGA", "TATCTAGCTGTATG", "CTGCTG", "CTTAGA", "CCCT", "CTGAATGT"), type = "nuc")

sq_ami <- construct_sq(c("NYMITGGREEYERTVIYRAIALNAANYTWIL", "TIAALGNIIYRAIE", "NYERTGHLI", "MAYNNNIALN", "MN", "NAAAT"), type = "ami")

sq_nuc_ex <- construct_sq(c("TATCTAGCTGTATG", "CUGCUG", "CUUAGA", "CCCT", "CUGAAUGU"))

sq_ami_ex <- construct_sq(c("MAYUOUONNNIALN", "UUMXBZONO", "NAAGAT"))  

 enc_nuc <- c(A = 1, C = 2, G = 2, T = 2)
 enc_ami <- c(A = 5, C = 5, D = 6, E = 6, F = 4, 
              G = 1, H = 5, I = 3, K = 2, L = 3, 
              M = 5, N = 6, P = 2, Q = 6, R = 2, 
             S = 6, T = 6, V = 3, W = 4, Y = 4)                   
                         
              
 # Create encoding for extended alphabet  
 
 enc_nuc_ex <- c(A = 1, C = 2, G = 2, T = 2, U = 3)
 enc_ami_ex <- c(A = 5, C = 5, D = 6, E = 6, F = 4, 
              G = 1, H = 5, I = 3, K = 2, L = 3, 
             M = 5, N = 6, P = 2, Q = 6, R = 2, 
             S = 6, T = 6, V = 3, W = 4, Y = 4, 
              U = -0.1, O = 0.753, X = -53.95, B = 7.77, Z = 0)

