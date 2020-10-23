#' Test sq object for presence of given motifs
#' 
#' @description Test if elements of a \code{\link{sq}} object contain given motifs
#' 
#' @param x a \code{\link{sq}} object to be tested.
#' @param y a \code{character} vector of motifs to be searched for.
#' 
#' @return A \code{\link{logical}} vector of the same length as input \code{sq}, indicating 
#' which elements contain all given motifs
#' 
#' @details This function allows testing if elements of a \code{sq} object contain 
#' a given motif or motifs, which includes start sequences and stop codons. It returns 
#' a \code{logical} for every element of the \code{sq} object - \code{TRUE} if it 
#' contains the motif and \code{FALSE} otherwise. When multiple motifs are searched, 
#' \code{TRUE} will be returned only for sequences that contain all the given motifs. 
#' 
#' This function only indicates if a motif is present within a sequence, to 
#' find all motifs, and their positions within sequences use 
#' \code{\link{find_motifs}}.
#' 
#' @section Allowed and forbidden letters and characters details:
#' Note: if a sq object contains characters: ^$?=()\.|+*{}[] in its alphabet, 
#' search for motifs cannot be performed and an error will be displayed (with 
#' exception of sq objects of type ami - in their alphabet there is '*' letter 
#' and it can be contained in sought motif). To search for motifs with those 
#' characters, you have to replace them first using 
#' \code{\link{substitute_letters}}. 
#' 
#' In case of sq objects of type \strong{ami}, \strong{dna} and \strong{rna},
#' motifs have to consist of upper case letters from amino acid, DNA and RNA alphabets
#' respectively. Use of lower case letters will return an error. Two additional
#' characters are allowed: '^' and '$' indicating the beginning and the end of a sequence 
#' respectively. Moreover, notice that '*' character may be used in amino acid 
#' motifs, as it is a part of the amino acid alphabet. If a motif contains 
#' ambiguous letters, all possible matches will be searched for, e.g., amino 
#' acid motif "MAJ" (where "J" is an ambiguous letter indicating L, or I) will 
#' find motifs: "MAJ", "MAL" and "MAI". 
#' 
#' Detailed list of all letters corresponding to each ambiguous letter may be found at
#' \code{\link{aminoacids_df}} and \code{\link{nucleotides_df}}.
#' 
#' @examples 
#' # Creating objects to work on:
#' sq_dna <- construct_sq(c("CTGAATGCAGTACCGTAAT", "ATGCCGTAAATGCCAT", 
#'                          "CAGACCANNNATAG"), type = 'dna')
#' sq_ami <- construct_sq(c("AGNTYIKFGGAYTI", "MATEGILIAADGYTWIL", 
#'                          "MIPADHICAANGIENAGIK"), type = 'ami')
#'                          
#' # Test if DNA sequences contain start codon 'ATG':
#' sq_dna %has% 'ATG'
#' 
#' # Test if DNA sequences begin with start codon 'ATG'
#' sq_dna %has% '^ATG'
#' 
#' # Test if DNA sequences end with 'TAG' (one of the stop codons):
#' sq_dna %has% 'TAG$'
#' 
#' # Test if amino acid sequences contain motif of two alanines followed by
#' # aspartic acid or asparagine ('AAB' motif will match 'AAB', 'AAD' and 'AAN'):
#' sq_ami %has% 'AAB'
#' 
#' # Test if amino acid sequences contain two motifs:
#' sq_ami %has% c("AAXG", "mat")
#' 
#' @seealso \code{\link{sq}} \code{\link{substitute_letters}} \code{\link{find_motifs}}
#' @export
`%has%` <- function(x, y) {
  assert_character(y, any.missing = FALSE, min.len = 1)
  
  UseMethod("%has%")
}

#' @export
`%has%.default` <- function(x, y)
  stop("operator '%has%' is not overloaded for this type of objects")

#' @export
`%has%.sq` <- function(x, y) {
  assert_alph_regex_friendly(alphabet(x))
  
  x <- as.character(x)
  
  ret <- do.call(cbind, lapply(y, function(s) grepl(s, x)))
  apply(ret, 1, all)
}

#' @export
`%has%.amisq` <- function(x, y) {
  y <- toupper(y)
  assert_motifs_for_type(y, "ami")
  
  y <- .replace_ami_motif(y)
  x <- as.character(x)
  
  ret <- do.call(cbind, lapply(y, function(s) grepl(s, x)))
  apply(ret, 1, all)
}

#' @export
`%has%.dnasq` <- function(x, y) {
  y <- toupper(y)
  assert_motifs_for_type(y, "dna")
  
  y <- .replace_dna_motif(y)
  x <- as.character(x)
  
  ret <- do.call(cbind, lapply(y, function(s) grepl(s, x)))
  apply(ret, 1, all)
}

#' @export
`%has%.rnasq` <- function(x, y) {
  y <- toupper(y)
  assert_motifs_for_type(y, "rna")
  
  y <- .replace_rna_motif(y)
  x <- as.character(x)
  
  ret <- do.call(cbind, lapply(y, function(s) grepl(s, x)))
  apply(ret, 1, all)
}

#' @importFrom stringi stri_replace_all_regex
.replace_ami_motif <- function(motif) {
  stri_replace_all_regex(
    motif,
    c("\\*", "B", "J", "Z", "X"),
    c("\\\\*", "[BDN]", "[JIL]", "[ZEQ]", "[A-Z]"),
    vectorize_all = FALSE
  )
}

#' @importFrom stringi stri_replace_all_regex
.replace_dna_motif <- function(motif) {
  stri_replace_all_regex(
    motif,
    c("W", "S", "M", "K", "R", "Y", "B", "D", "H", "V", "N"),
    c("[WAT]", "[SCG]", "[MAC]", "[KGT]", "[RAG]", "[YCT]",
      "[BSKYCTG]", "[DWKRATG]", "[HWMYACT]", "[VSMRACG]",
      "[NBDHVWSMKRYACTG]"),
    vectorize_all = FALSE
  )
}

#' @importFrom stringi stri_replace_all_regex
.replace_rna_motif <- function(motif) {
  stri_replace_all_regex(
    motif,
    c("W", "S", "M", "K", "R", "Y", "B", "D", "H", "V", "N"),
    c("[WAU]", "[SCG]", "[MAC]", "[KGU]", "[RAG]", "[YCU]",
      "[BSKYCUG]", "[DWKRAUG]", "[HWMYACU]", "[VSMRACG]",
      "[NBDHVWSMKRYACUG]"),
    vectorize_all = FALSE
  )
}
