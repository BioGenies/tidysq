#' @exportClass simsq
#' @export
simplify <- function(sqtbl, encoding) {
  validate_sqtibble(sqtbl)
  sqtypes <- extract_sq_types(sqtbl)
  if (nrow(sqtbl) == 0) {
    stop("you cannot simplify sqtibble with less than one row")
  }
  if (length(unique(sqtypes)) != 1 ||
      any(sqtypes == "sim") ||
      any(sqtypes == "unt")) {
    stop("'sqtbl' should contain only sequences of the same type - either 'aa' or 'nuc'")
  } 
  ind_cln <- extract_is_clean(sqtbl)
  if (any(ind_cln) && !all(ind_cln)) {
    warning("'sqtbl' contain at least one cleaned sequence and at least one uncleaned")
  }
  
  if (!is.character(encoding) ) {
    stop("'encoding' should be a named vector, where names are upper latin letters or '-' and elements are lower latin letters and '-'")
  }
  if (all(ind_cln)) {
    if ((sqtypes[1] == "aa" && !setequal(names(encoding), aminoacids_df[!aminoacids_df[["amb"]],"one"])) ||
        (sqtypes[1] == "nuc" && !setequal(names(encoding), nucleotides_df[!nucleotides_df[["amb"]],"one"]))) {
      stop("names of encoding should be equal to alphabet that is going to be encoded (cleaned alphabet)")
    }
  } else {
    if ((sqtypes[1] == "aa" && !setequal(names(encoding), aminoacids_df[,"one"])) ||
        (sqtypes[1] == "nuc" && !setequal(names(encoding), nucleotides_df[,"one"]))) {
      stop("names of encoding should be equal to alphabet that is going to be encoded")
    }
  }
  
  sqtbl[["sq"]] <- lapply(sqtbl[["sq"]], function(sq) {
    levels(sq) <- encoding[levels(sq)] 
    class(sq) <- c("simsq", "sq", "factor")
    sq
  })
  set_sqcol(sqtbl)  
}

