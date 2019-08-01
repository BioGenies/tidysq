#' Substitute default amino acid or nucleic acid alphabet with a custom encoding
#' 
#' @description Replace default amino acid or nucleic acid letters in a sequence, 
#' stored in \code{\link{sq}} object, with a custom encoding. 
#' Selected letters in the amino/nucleic acid sequence are replaced 
#' by a user-defined symbols.
#' 
#' 
#' @param sq \code{\link{sq}} object.
#' @param indices \code{encoding} vector of letters to be replaced together with their replacements.
#' One letter can be replaced with multiple symbols. 
#' To perform substitution create a named vector \item{c(A = Ala, H = His, amino_or_nucleic_acid_symbol = replacement)}.
#' 
#' @return \code{\link{atpsq}} object of the same type as input sq with replaced alphabet, defined by user.
#' 
#' @details \code{substitute_letters} allows to replace desired letters in the amino acid or nucleic acid sequence.
#' One letter of the alphabet may be replaced by a multiple characters. The function allows to replace single, 
#' multiple letters or even a whole alphabet.
#' 
#' Sometimes one needs to introduce artificial amino acids or nucleotides into the sequence, replacing
#' others or ambiguous ones. 
#' 
#' The alphabet characters to be replaced need to be written in capital letters and must originate from default alphabets, otherwise error will be introduced.
#' This will occur even when the letter to be replaced won't occur in the sequence.
#' 
#' Created sequence will be deprived of \code{\link{cln})} subtype, if the original sequence possessed it.
#' 
#' 
#' @examples 
#' # Creating object, called sq to work on:
#'
#' sq_nuc <- construct_sq(c("TATGAATTAGCTGTCTTTGCTGCTTTGGTTATCTATGA", "CTTTGGTTATCTAGCTGTATGA", 
#'                         "TATCTAGCTGTATG", "CTGCTG", "CTTAGA", "CCCT", "CTGAATGT"), type = "nuc")
#' 
#' sq_ami <- construct_sq(c("NYMITGGREEYERTVIYRAIALNAANYTWIL", "TIAALGNIIYRAIE", "NYERTGHLI", 
#'                         "MAYNNNIALN", "MN", "NAAAT"), type = "ami")
#'                      
#' 
#' # Replace single letter of alphabet with single character encoding:
#' 
#' substitute_letters(sq_nuc, c(T = "t", A = "a", C = "H", G = "Z"))
#' substitute_letters(sq_nuc, c(T = 1, A = 2, C = 3, G = 4))
#' 
#' substitute_letters(sq_ami, c(M = "m", Q = "g", R = "#", D = "$"))
#' substitute_letters(sq_ami, c(M = "2", Q = "5", R = "9", D = "7"))
#' 
#' 
#' # Replace single letter of alphabet with multiple character encoding:
#' 
#' substitute_letters(sq_nuc, c(T = "th", A = "ad", C = "cy", G = "gu"))
#' substitute_letters(sq_nuc, c(T = 111, A = 222, C = 333, G = 444))
#' 
#' substitute_letters(sq_ami, c(M = "Met", Q = "Gln", R = "Arg", D = "Asp"))
#' substitute_letters(sq_ami, c(M = "222", Q = "555", R = "999", D = "777"))
#' 
#' 
#' # Use created encoding
#' 
#' enc_nuc <- c(T = "t", A = "a", C = "c", G = "g")
#' enc_ami <- c(M = "Met", Q = "Gln", R = "Arg", D = "Asp", H = "His", K = "Lys", A = "Ala")
#' 
#' substitute_letters(sq_nuc, enc_nuc)
#' substitute_letters(sq_ami, enc_ami)
#' 
#' 
#' # Use created encoding from other package (ex. \code{\link[AmyloGram]{myloGram_model}})
#' 
#' AG_enc_raw <- unlist(AmyloGram_model[["enc"]])
#' 
#' enc_AG <- substr(names(AG_enc_raw), 1, 1)
#' names(enc_AG) <- toupper(AG_enc_raw)
#' enc_AG 
#' 
#' sq_ami <- construct_sq(c("NYMITGGREEYERTVIYRAIALNAANYTWIL", "TIAALGNIIYRAIE", "NYERTGHLI", 
#'                         "MAYNNNIALN", "MN", "NAAAT"), type = "ami")
#' substitute_letters(sq_ami, enc_AG)
#' 
#' @seealso sq atpsq
#' @export

substitute_letters <- function(sq, encoding) {
  validate_sq(sq)
  
  alph <- .get_alph(sq)
  
  if (!all(names(encoding) %in% alph)) {
    stop("all names of 'encoding' has to be letters from alphabet (elements of 'alphabet' attribute of 'sq')")
  }
  
  names(alph) <- 1:length(alph)
  alph_inds <- !(alph %in% names(encoding))
  names(encoding) <- match(names(encoding), alph)
  
  transl_table <- c(alph[alph_inds], encoding)
  new_alph <- na.omit(unique(transl_table))
  inds_func <- match(transl_table, new_alph)
  names(inds_func) <- names(transl_table)
  
  ret <- .recode_sq(sq, alph, new_alph, inds_func)
  if (.is_cleaned(sq)) {
    .handle_opt_txt("tidysq_subsitute_letters_cln",
                    "column passed to muatting had 'cln' subtype, output column doesn't have it")
  }
  class(ret) <- c("atpsq", "sq")
  attr(ret, "alphabet") <- new_alph[!is.na(new_alph)]
  ret
}