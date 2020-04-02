#' @importFrom stringi stri_split_regex
#' @importFrom tibble tibble
#' @export
digest <- function(ami_sq, name, digest_pattern) {
  validate_sq(ami_sq)
  .check_character(name, "name")
  .check_eq_lens(ami_sq, name, "ami_sq", "name")
  .check_character(digest_pattern, "digest_pattern")
  .check_sq_has_type(ami_sq, "ami_sq", "ami")
  .check_is_clean(ami_sq, "ami_sq")
  
  
  .check_digest_pattern(digest_pattern)
  digest_pattern = .regexify_pattern(digest_pattern)
  detected <- stri_split_regex(as.character(ami_sq), digest_pattern)
  tibble(name = unlist(lapply(1:length(name), function(i) rep(name[i], length(detected[[i]])))),
         peptides = construct_sq_ami(unlist(detected), TRUE))
}