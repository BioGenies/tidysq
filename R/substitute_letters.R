#' @export
substitute_letters <- function(sq, encoding) {
  validate_sq(sq)
  
  alph <- .get_alph(sq)
  
  if (!all(names(encoding) %in% alph)) {
    stop("all names of 'encoding' has to be letters from alphabet (elements of 'alphabet' attribute of 'sq')")
  }
  
  inds_fun <- alph
  inds_fun[match(names(encoding), alph)] <- encoding
  new_alph <- na.omit(unique(inds_fun))
  names(inds_fun) <- as.character(1:length(alph))
  inds_fun <- match(inds_fun, new_alph)
  inds_fun[is.na(inds_fun)] <- .get_na_val(new_alph)
  
  ret <- .apply_sq(sq, "int", "int", function(s) {
    inds_fun[s]
  })
  if (.is_cleaned(sq)) {
    .handle_opt_txt("tidysq_subsitute_letters_cln",
                    "column passed to muatting had 'cln' subtype, output column doesn't have it")
  }
  class(ret) <- c("atpsq", "sq")
  attr(ret, "alphabet") <- new_alph
  ret
}