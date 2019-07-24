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