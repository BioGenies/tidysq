#' Encode sequences using simplified or custom alphabets into numeric values
#' 
#' @description Function allows to encode whole sequences using simplified or user-defined alphabets.
#' Each letter must be must be assigned to it's numeric value. 
#' 
#' @param sq \code{\link{sq}} object.
#' @param indices \code{encoding} named or numeric vector. 
#' 
#' Named vector must consist of nucleic or amino acid letters with their numeric encoding.
#' All standard nucleic/amino acid letters must be assigned, otherwise an error will be introduced.
#' 
#' Non-standard letters can only be assigned in named vector.
#' 
#' Values in numeric vector must correspond to letters in alphabetic order. This will work only for standard alphabet.
#' For nucleic acids - A, C, G, T. For amino acids - A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y.
#' 
#' 
#' @return \code{\link{encsq}} object of the same type as input sq with encoded alphabet.
#' 
#' @details \code{encode}
#' 
#'   allows to replace ambigous/extraordinary 
#' letters in nucleic or amino acid sequence with user-defined or IUPAC symbols. 
#' Letters can also be replaced with \code{NA} values, so that they can be later 
#' removed, from the sequence, by \code{clean} function.
#' 
#' \code{substitute_letters} can be used to replace default amino acid letters 
#' with encodings. They can be user-defined or be derived from various simplified alphabets.
#' 
#' One letter of the alphabet may be replaced by a multiple characters. 
#' 
#' The alphabet characters to be replaced need to be written in capital letters and must originate from default alphabets, otherwise error will be introduced.
#' Multiple string of letters to be substituted (ex. \item{c(AHG = "replacement")}) will also produce an error.
#' 
#' Replacing multiple letters with the same symbol (ex. \item{c(A = "replacement1", H  = "replacement1", G = "replacement1")}) is allowed.
#' 
#' Created sequence will be deprived of \code{\link{cln})} subtype, if the original sequence possessed it.
#' This will also occur when the letter to be replaced will not be found in the sequence. It remain unchanged but will lose subclass.
#' 
#' The newly constructed will have a new class \code{\link{cln})}, representing atypical alphabet.
#' 
#' 
#' @examples 
#' 
#' 
#' @seealso sq atpsq
#' enc (encoded, created during encoding
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






library(tidysq)
sq_ami <- construct_sq(c("NYMITGGREEYERTVIYRAIALNAANYTWIL", "TIAALGNIIYRAIE", "NYERTGHLI", "MAYNNNIALN", "MN", "NAAAT"), type = "ami")
sq_ami <- construct_sq(c("MN", "NAAAT"), type = "ami")

sq_nuc <- construct_sq(c("ATGC", "AAAUU"), type = "nuc")


library(AmyloGram)

AG_enc_raw <- unlist(AmyloGram_model[["enc"]])

enc_AG <- as.numeric(substr(names(AG_enc_raw), 1, 1))
# enc_AG <- as.numeric(names(AG_enc_raw))
names(enc_AG) <- toupper(AG_enc_raw)
enc_AG 


encode(sq_ami, enc_AG) 
zz <- encode(sq_ami, enc_AG)
encode(sq_ami, c(A = 1))

encode(sq_ami, c(seq(1, 20)))

substitute_letters(sq_ami, enc_AG)

sq_ami <- construct_sq(c("MN", "NAAAT"), type = "ami")
sq_ami1 <- construct_sq(c("AACDUO", "KLM"), type = "ami")

encode(sq_ami, c(A = 1, C= 1, D=1, E=1, F=1, G=1, H=1, I=1, K=1, L=1, M=1, N=1, P=1, Q=1, R=1, S=1, T=1, V=1, W=1, Y=1, U=2, O=3))
encode(sq_ami1, c(A = 1, C= 1, D=1, E=1, F=1, G=1, H=1, I=1, K=1, L=1, M=1, N=1, P=1, Q=1, R=1, S=1, T=1, V=1, W=1, Y=1, U=2, O=3))
encode(sq_ami, seq(1:20))
encode(sq_ami1, seq(1:20))


sq_nuc <- construct_sq(c("ATGC", "AAAUU"), type = "nuc")

encode(sq_nuc, c(1, 2, 3, 4))
encode(sq_nuc, c(1,2,3,4,5))
encode(sq_nuc, c(A=1, T=4, G=3, C=2, U=5))

