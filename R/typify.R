#' @export
typify <- function(sq, dest_type) {
  validate_sq(sq)
  if (missing(dest_type) || 
      !(dest_type %in% c("ami", "nuc"))) {
    stop("'dest_type' should be either 'ami' or 'nuc'")
  }
  
  type <- .get_sq_type(sq)
  if (type %in% c('ami', 'nuc')) {
    return(sq)
  }
  
  alph <- .get_alph(sq)
  up_alph <- unique(toupper(alph))
  dest_alph <- if (dest_type == "ami") aminoacids_df[, "one"] else nucleotides_df[, "one"]
  
  if (!all(up_alph %in% dest_alph)) {
    stop("some sequences have levels that are invalid for given 'dest_type'; you can check them with 'get_invalid_letters' function and fix them with 'substitute_letters'")
  }
  
  if (!(length(alph) == length(up_alph))) {
    .handle_opt_txt("tidysq_typify_small_cap_let",
                    "in 'alphabet' attribute of 'sq' some letters show up as both lower and capital")
  }
  
  inds_func <- match(toupper(alph), dest_alph)
  names(inds_func) <- as.character(1:length(alph))
  
  ret <- lapply(sq, function(s) {
    unname(inds_func[as.character(s)])
  })
  
  class(ret) <- c(paste0(dest_type, "sq"), "sq")
  attr(ret, "alphabet") <- dest_alph
  ret
}