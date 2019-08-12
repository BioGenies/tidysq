#' @export
random_sq <- function(n, len, type, is_clean, sd = NULL, use_gap = FALSE) {
  .check_integer(n, "'n'", single_elem = TRUE)
  .check_integer(len, "'len'", single_elem = TRUE)
  .check_type(type)
  .check_logical(is_clean, "'is_clean'", single_elem = TRUE)
  .check_numeric(sd, "'sd'", allow_null = TRUE, single_elem = TRUE)
  .check_logical(use_gap, "'use_gap'", single_elem = TRUE)
  
  alph <- .get_standard_alph(type, is_clean)
  if (!use_gap) alph <- setdiff(alph, "-")
  if (type == "ami") alph <- setdiff(alph, "*")
  
  if (is.null(sd))
    sq <- sapply(1:n, function(i) paste0(sample(alph, len, replace = TRUE), collapse = ""))
  else {
    len <- round(rnorm(n, len, sd))
    len <- ifelse(len <= 0, 1, len)
    sq <- sapply(1:n, function(i) paste0(sample(alph, len[i], replace = TRUE), collapse = ""))
  }
  construct_sq(sq, type, is_clean)
}