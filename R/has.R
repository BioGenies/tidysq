#' Test sq object for presence of given motifs
#' 
#' @description Test if elements of a sq object contain given motifs
#' 
#' @param x \code{\link{sq}} object to be tested
#' @param y \code{character} vector of motifs to be searched for
#' 
#' @return a logical vector of the same length as input sq, indicating 
#' which elements contain all of given motifs
#' 
#' @details This function allows testing if elements of a sq object contain 
#' a given motif or motifs. It returns a logical for every element of the sq 
#' object - \code{TRUE} if it contains the motif and \code{FALSE} otherwise. 
#' In case of search for multiple motifs, \code{TRUE} will be returned only 
#' for sequences that contain all of the given motifs. 
#' 
#' Note if a sq object contains characters: ^$?=()\.|+*{}[] in its alphabet, 
#' search for motifs cannot be performed and an error will be displayed (with 
#' exception of sq objects of type ami - in their alphabet there is "*" letter 
#' and it can be contained in sought motif"). To search for motifs with those 
#' characters, you have to replace them first using 
#' \code{\link{substitute_letters}}. 
#' 
#' In case of sq objects of type 'ami' and 'nuc', motifs have to consist 
#' of letters from amino acid and nucleotide alphabet respectively. 
#' However, two additional characters are allowed: '^' and '$' indicating
#' the beginning and the end of a sequence respectively. Moreover, notice
#' that '*' character may be used in amino acid motifs, as it is a part of 
#' the amino acid alphabet. If a motif contains ambiguous letters, all 
#' possible matches will be tested, e.g., amino acid motif "MAJ" (where "J" 
#' is an ambiguous letter indicating L or I) will return \code{TRUE} for 
#' sequences containing "MAJ", "MAL" and "MAI". If a motif contains lower 
#' case letters, they will be converted to upper case.  
#' 
#' This function only indicates if a motif is present within a sequence, to 
#' find all motifs and their positions within sequences use 
#' \code{\link{find_motifs}}.
#' 
#' @examples 
#' # Creating objects to work on:
#' sq_nuc <- construct_sq(c("CTGAATGCAGTACCGTAAT", "ATGCCGTAAATGCCAT", 
#'                          "CAGACCANNNATAG"), type = 'nuc')
#' sq_ami <- construct_sq(c("AGNTYIKFGGAYTI", "MATEGILIAADGYTWIL", 
#'                          "MIPADHICAANGIENAGIK"), type = 'ami')
#'                          
#' # Test if nucleotide sequences contain start codon 'ATG':
#' sq_nuc %has% 'ATG'
#' 
#' # Test if nucleotide sequences begin with start codon 'ATG'
#' sq_nuc %has% '^ATG'
#' 
#' # Test if nucleotide sequences end with 'TAG' (one of the stop codons):
#' sq_nuc %has% 'TAG$'
#' 
#' # Test if amino acid sequences contain motif of two alanines followed by
#' # aspartic acid or asparagine ('AAB' motif will match 'AAB', 'AAD' and 'AAN'):
#' sq_ami %has% 'AAB'
#' 
#' # Test if amino acid sequences contain two motifs:
#' sq_ami %has% c("AAXG", "mat")
#' 
#' @seealso sq substitute_letters find_motifs
#' @export
`%has%` <- function(x, y) {
  UseMethod("%has%")
}

#' @exportMethod `%has%`
#' @export
`%has%.default` <- function(x, y) {
  stop("operator '%has%' is not overloaded for this type of objects")
}

#' @exportMethod `%has%` sq
#' @export
`%has%.sq` <- function(x, y) {
  if (!is.character(y)) {
    stop("object which you're looking for in 'sq' object needs to be a character vector")
  }
  alph <- .get_alph(x)
  if (!all(unlist(strsplit(y, "")) %in% c(alph, "^", "$"))) {
    stop("motifs that you're searching for in 'sq' object needs to consist of letters from alphabet of 'sq'")
  }
  if (any(alph %in% c("^", "$", "?", "(", "=", ")", "\\", ".", "|", "+", "*", "{", "}", "[", "]"))) {
    stop("you cannot search for motifs if any of those characters: ^$?=()\\.|+*{}[] are elements of 'sq' alphabet; if you want to use them, please substitute those letters with some other using 'substitute_letters'")
  }
  
  
  x <- as.character(x)
  
  ret <- sapply(y, function(s) grepl(s, x))
  ret <- apply(ret, 1, all)
  ret
}

#' @exportMethod `%has%` amisq
#' @export
`%has%.amisq` <- function(x, y) {
  if (!is.character(y)) {
    stop("object which you're looking for in 'sq' object needs to be a character vector")
  }
  y <- strsplit(toupper(y), "")
  if (!all(unlist(y) %in% c(.get_standard_alph("ami", FALSE), "^", "$"))) {
    stop("motifs that you're searching for in 'sq' object needs to consist of letters from aminoacids alphabet and optionally '^' or '$' characters")
  }
  y <- lapply(y, function(s) replace(s, s == "B", "[BDN]"))
  y <- lapply(y, function(s) replace(s, s == "J", "[JIL]"))
  y <- lapply(y, function(s) replace(s, s == "X", "[A-Z]"))
  y <- lapply(y, function(s) replace(s, s == "Z", "[ZEQ]"))
  y <- sapply(y, function(s) paste(s, collapse = ""))
  
  alph <- .get_alph(x)
  x <- as.character(x)
  
  ret <- sapply(y, function(s) grepl(s, x))
  ret <- apply(ret, 1, all)
  ret
}

#' @exportMethod `%has%` nucsq
#' @export
`%has%.nucsq` <- function(x, y) {
  if (!is.character(y)) {
    stop("object which you're looking for in 'sq' object needs to be a character vector")
  }
  y <- strsplit(toupper(y), "")
  if (!all(unlist(y) %in% c(.get_standard_alph("nuc", FALSE), "^", "$"))) {
    stop("motifs that you're searching for in 'sq' object needs to consist of letters from nucleotides alphabet and optionally '^' or '$' characters")
  }
  y <- lapply(y, function(s) replace(s, s == "W", "[WATU]"))
  y <- lapply(y, function(s) replace(s, s == "S", "[SCG]"))
  
  y <- lapply(y, function(s) replace(s, s == "M", "[MAC]"))
  y <- lapply(y, function(s) replace(s, s == "K", "[KGTU]"))
  y <- lapply(y, function(s) replace(s, s == "R", "[RAG]"))
  y <- lapply(y, function(s) replace(s, s == "Y", "[YCTU]"))
  
  y <- lapply(y, function(s) replace(s, s == "B", "[BCTGU]"))
  y <- lapply(y, function(s) replace(s, s == "D", "[DATGU]"))
  y <- lapply(y, function(s) replace(s, s == "H", "[HACTU]"))
  y <- lapply(y, function(s) replace(s, s == "V", "[VACG]"))
  
  y <- lapply(y, function(s) replace(s, s == "N", "[ACTGUWSMKRYBDHVN]"))
  
  y <- sapply(y, function(s) paste(s, collapse = ""))
  
  alph <- .get_alph(x)
  x <- as.character(x)
  
  ret <- sapply(y, function(s) grepl(s, x))
  ret <- apply(ret, 1, all)
  ret
}