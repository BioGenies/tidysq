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
  y <- lapply(y, function(s) replace(s, s == "*", "\\*"))
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