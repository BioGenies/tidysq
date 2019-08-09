#' @export
random_sq <- function(n, len, type, is_clean, sd = NULL, use_gap = FALSE) {
  .check_n_is_int(n)
  .check_len_is_int(len)
  .check_type_in_ami_nuc(type)
  .check_is_clean_in_TRUE_FALSE(is_clean)
  .check_sd_is_numeric_or_NULL(sd)
  .check_use_gap_in_TRUE_FALSE(use_gap)
  
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