#' Substitute letters in a sequence
#' 
#' @description 1) Replace ambigous/extraordinary letters in nucleic or 
#' amino acid sequence, stored in \code{\link{sq}} object, with the ones 
#' that are compliant with the IUPAC standard, ones that are user-defined 
#' or with \code{NA} values.
#' 
#' 2) Replace default amino acid letters in a sequence with a custom encoding 
#' to create simplified alphabets.
#' 
#' The function is only used to replace letters in the alphabet. 
#' It cannot be used to merge surrounding characters.
#' 
#' 
#' @param sq \code{\link{sq}} object.
#' @param encoding a vector of letters to be replaced together with their replacements.
#' One letter can be replaced with multiple symbols. 
#' To perform substitution create a named vector ex. 
#' \code{c(A = Ala, H = His, amino_or_nucleic_acid_symbol = replacement)}.
#' 
#' @return a \code{\link{sq}} object with \strong{atp}  type with replaced alphabet, 
#' defined by user.
#' 
#' @details \code{substitute_letters} allows to replace ambigous/extraordinary 
#' letters in nucleic or amino acid sequence with user-defined or IUPAC 
#' symbols. Letters can also be replaced with \code{NA} values, so that they 
#' can be later removed, from the sequence, by \code{clean} function.
#' 
#' \code{substitute_letters} can be used to replace default amino acid letters 
#' with encodings. They can be user-defined or be derived from various 
#' simplified alphabets.
#' 
#' One letter of the alphabet may be replaced by a multiple characters. 
#' 
#' The alphabet characters to be replaced need to be written in capital letters
#' and must originate from default alphabets, otherwise error will be 
#' introduced.
#' 
#' Multiple string of letters to be substituted 
#' (ex. \code{c(AHG = "replacement")}) will also produce an error.
#' 
#' Replacing multiple letters with the same symbol 
#' (ex. \code{c(A = "rep1", H  = "rep1", G = "rep1")}) is allowed.
#' 
#' Created sequence will be deprived of \strong{cln} subtype, 
#' if the original sequence possessed it. This will also occur when
#' the letter to be replaced will not be found in the sequence. 
#' It remain unchanged but will lose subclass.
#' 
#' The newly constructed will have a new type \strong{atp}, 
#' representing sequences with atypical alphabet.
#' 
#' All replaced letters will have the character type. 
#' Multiple symbol replacement will be recognized as one letter and 
#' cannot be separated in future operations into single letters. 
#' 
#' @examples 
#' # Creating sq object to work on:
#'
#' sq_nuc <- construct_sq(c("TATGAATTAGCTGTCTTTGCTGCTTTGGTTATCTATGA", 
#'                          "CTTTGGTTATCTAGCTGTATGA", "TATCTAGCTGTATG", 
#'                          "CTGCTG", "CTTAGA", "CCCT", "CTGAATGT"), 
#'                        type = "nuc")
#' 
#' sq_ami <- construct_sq(c("NYMITGGREEYERTVIYRAIALNAANYTWIL", 
#'                         "TIAALGNIIYRAIE", "NYERTGHLI", 
#'                         "MAYNNNIALN", "MN", "NAAAT"), 
#'                         type = "ami")
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
#' sub_nuc <- c(T = "t", A = "a", C = "c", G = "g")
#' sub_ami <- c(M = "Met", Q = "Gln", R = "Arg", D = "Asp", 
#'              H = "His", K = "Lys", A = "Ala")
#' 
#' substitute_letters(sq_nuc, sub_nuc)
#' substitute_letters(sq_ami, sub_ami)
#' 
#' 
#' # Use created encoding from other package 
#' # (ex. \code{\link[AmyloGram]{AmyloGram_model}})
#' 
#' library(AmyloGram)
#' 
#' AG_sub_raw <- unlist(AmyloGram_model[["enc"]])
#' 
#' sub_AG <- substr(names(AG_sub_raw), 1, 1)
#' names(sub_AG) <- toupper(AG_sub_raw)
#' sub_AG 
#' 
#' substitute_letters(sq_ami, sub_AG)
#' 
#' @seealso sq atpsq
#' @export

substitute_letters <- function(sq, encoding) {
  validate_sq(sq)
  alph <- .get_alph(sq)
  .check_isnt_missing(encoding, "'encoding'")
  .check_isnt_null(encoding, "'encoding'")
  .check_is_named(encoding, "'encoding'")
  .check_enc_names_in_alph(encoding, alph)
  .check_is_unique(names(encoding), "names of 'encoding'")
  if (is.numeric(encoding)) {
    .check_integer(encoding, "if is numeric, 'encoding'", allow_na = TRUE)
    name <- names(encoding)
    encoding <- as.character(encoding)
    names(encoding) <- name
  } else .check_character(encoding, "if is character, 'encoding'", allow_na = TRUE)
  
  inds_fun <- alph
  inds_fun[match(names(encoding), alph)] <- encoding
  new_alph <- na.omit(unique(inds_fun))
  names(inds_fun) <- as.character(1:length(alph))
  inds_fun <- match(inds_fun, new_alph)
  inds_fun[is.na(inds_fun)] <- .get_na_val(new_alph)
  
  ret <- .apply_sq(sq, "int", "int", function(s) {
    inds_fun[s]
  }, new_alph)
  if (.is_cleaned(sq)) {
    .handle_opt_txt("tidysq_subsitute_letters_cln",
                    "'sq' object passed to substitute_letters had 'cln' subtype, output doesn't have it")
  }

  ret <- .set_alph(ret, new_alph)
  .set_class(ret, "atp")
}