#' Substitute default amino acid or nucleic acid alphabet with a custom encoding
#' 
#' @description Replace default amino acid or nucleic acid letters in a sequence, 
#' stored in \code{\link{sq}} object, with a custom encoding. 
#' Selected letters in the amino acid or nucleic acid sequence are replaced 
#' by a user-defined symbols.
#' 
#' 
#' @param sq \code{\link{sq}} object.
#' @param indices \code{encoding} vector of letters to be replaced together with their replacements.
#' 
#' @return \code{\link{atpsq}} object of the same type as input sq with replaced alphabet.
#' 
#' @details Function allows to replace default alphabet encoding with letters desired by the user.
#' One letter of the alphabet may be replaced by a string of characters.
#' 
#' Sometimes one needs to replace default amino/nucleic acid alphabet with custom one. 
#' Such an example could be the use of simplified amino acid alphabet, which take into account 
#' different physicochemical properties of amino acids to cluster them into fewer groups.
#' The simplification preserves the informative character of the alphabet while reducing the 
#' number of required operations when using it in futher steps of pipeline, such as machine learning.
#' 
#' 
#' #' @details Amino acids and nucleic acid sequences are represented as \code{\link{sq}} 
#' object in the \code{\link{tidysq}} package. Often one needs to get only a 
#' single letter or the sequence of a defined range from the original sequences. 
#' A subsequence is a sequence that can be derived from the original sequence 
#' by trimming some elements (letters) without changing the order of the 
#' remaining elements. To obtain a subsequence from each sequence contained in 
#' the \code{\link{sq}} object with the same indices. This is for example 
#' useful to extract an user-defined region from a sequence. 
#' 
#' The usage of \code{bite} follows the normal R conventions. For details 
#' refer to the R documentation (see 
#' \url{https://cran.r-project.org/doc/manuals/r-release/R-intro.html#Index-vectors}). 
#' Extracting indices not present in the sequence results in introducing 
#' \code{\link{NA}} (‘Not Available’ / Missing Values) values. 
#' Information about it is printed on console depending on value of option 
#' 'tidysq_bite_na_action' - it can be either a warning (default), error, 
#' message or no information (you can check details in \code{\link{sq-options})}. 
#' \code{NA} values can be removed by using \code{\link{remove_na}} function.
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