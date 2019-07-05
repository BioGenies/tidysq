#' @export
simplify <- function(sqtbl, encoding) {
  validate_sqtibble(sqtbl)
  sqtypes <- extract_sq_types(sqtbl)
  if (nrow(sqtbl) == 0) {
    stop("you cannot encode sqtibble with less than one row")
  }
  if (length(unique(sqtypes)) != 1 ||
      any(sqtypes == "sim") ||
      any(sqtypes == "unt")) {
    stop("'sqtbl' should contain only sequences of the same type - either 'aa' or 'nuc'")
  } 
  
  if (!is.character(encoding) ) {
    stop("'encoding' should be a named vector, where names are upper latin letters or '-' and elements are lower latin letters and '-'")
  }
  if ((sqtypes[1] == "aa" && !setequal(names(encoding), aminoacids_df[,"one"])) ||
      (sqtypes[1] == "nuc" && !setequal(names(encoding), nucleotides_df[,"one"]))) {
    stop("names of encoding should be equal to alphabet that is going to be encoded")
  }
  
  sqtbl[["sq"]] <- lapply(sqtbl[["sq"]], function(sq) {
    levels(sq) <- encoding[levels(sq)] 
    class(sq) <- c("simsq", "sq", "factor")
    sq
  })
  set_sqcol(sqtbl)  
}

