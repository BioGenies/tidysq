#' Import sq objects from other objects
#' 
#' Creates \code{\link[sq]{sq object}} from \code{object} of class from another package.
#' Currently supported packages are \code{ape} with its formats (\code{AAbin} and \code{DNAbin}),
#' \code{Bioconductor} (\code{AAStringSet}, \code{DNAStringSet}) and
#' \code{seqinr} (\code{SeqFastaAA}, \code{SeqFastadna}).
#' 
#' @param object - an object of one of classes: \code{AAbin}, \code{DNAbin}, \code{AAStringSet}, 
#' \code{DNAStringSet}, \code{SeqFastaAA}, \code{SeqFastadna}
#' @return a \code{tibble} with \code{sq} column of \code{\link{sq}} type representing the same 
#' sequences as given object; the object has a type corresponding to the input type; if given
#' sequences had names, output tibble has also another column \code{name} with those names
#' 
#' @details 
#' Providing object of class other than specified will result in error.
#' 
#' @examples 
#' ## ape example
#' library(ape)
#' ape_nuc <- as.DNAbin(list(one = c("C", "T", "C", "A"), two = c("T", "G", "A", "G", "G")))
#' import_sq(ape_nuc)
#' 
#' ## Biostrings example
#' library(Biostrings)
#' Biostrings_nuc <- DNAStringSet(c(one = "CTCA", two = "TGAGG"))
#' import_sq(Biostrings_nuc)
#' 
#' ## seqinr example
#' seqinr_nuc <- as.SeqFastadna(list(one = c("C", "T", "C", "A"), 
#'                                   two = c("T", "G", "A", "G", "G")))
#' import_sq(seqinr_nuc)
#' 
#' @seealso \code{\link{export_sq}} \code{\link{sq}}
#' @export
import_sq <- function(object) {
  if (class(object) %in% c("AAbin", "DNAbin")) {
    sq <- as.character(object)
    if (is.matrix(sq)) {
      name <- rownames(sq)
      sq <- construct_sq(apply(sq, 1, function(s) paste(s, collapse = "")))
    } else if (is.list(sq)) {
      name <- names(sq)
      sq <- construct_sq(sapply(sq, function(s) paste(s, collapse = "")))
    }
    
    if (is.null(name))
      tibble(sq = sq)
    else 
      tibble(name = name, sq = sq)
  } else if (class(object) %in% c("AAStringSet", "DNAStringSet")) {
    sq <- construct_sq(as.character(object))
    name <- names(object)
    
    if (is.null(name))
      tibble(sq = sq)
    else 
      tibble(name = name, sq = sq)
  } else if (class(object) %in% c("SeqFastaAA", "SeqFastadna")) {
    sq <- sapply(object, function(s) paste(s, collapse = ""))
    name <- names(sq)
    if (is.null(sq))
      name <- sapply(object, function(s) attr(s, "name"))
    
    sq <- construct_sq(sq)
    if (is.null(name))
      tibble(sq = sq)
    else 
      tibble(name = name, sq = sq)
  } else stop("this function cannot handle objects with class as given", call. = FALSE)
}