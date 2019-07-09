#' @export
bite <- function(sq, indices) {
  validate_sq(sq)

  if (!(is.numeric(indices) && 
        floor(indices) == indices)) {
    stop("'indicies' has to be an integer vector")
  }
  
  had_na <- any(sapply(sq, function(s) any(is.na(s))))
  ret <- lapply(sq, function(s) s[indices])
  has_na <- any(sapply(ret, function(s) any(is.na(s))))
  if (has_na & !had_na) {
    handle_opt_txt("tidysq_bite_na_action",
                   "some sequences are subsetted with index bigger than length - NA's introduced")
  }
  .set_class_alph(ret, sq)
}
