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
  } else if (class(object) %in% c("SeqFastaAA", "SeqFastaDNA")) {
    sq <- sapply(object, function(s) paste(s, collapse = ""))
    name <- names(sq)
    if (is.null(sq))
      name <- sapply(object, function(s) attr(s, "name"))
    
    if (is.null(name))
      tibble(sq = sq)
    else 
      tibble(name = name, sq = sq)
  } else stop("this function cannot handle objects with class as given")
}