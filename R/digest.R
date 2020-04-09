#' Digest proteins at specified cleavage sites
#' 
#' TEMPORARY DOCS
#' 
#' @param ami_sq a \code{\link{sq}} object with \strong{ami} type and \strong{cln} subtype.
#' @param name a non-\code{NULL} \code{character} vector without \code{\link{NA}} values, 
#' containing names of the sequences in the sq. It has to be of the same length 
#' as the \code{sq}. 
#' @param digest_pattern a \code{character} protease cleavage site pattern.
#' 
#' @return A \code{\link[tibble]{tibble}} with following columns:
#'  \item{name}{name of the sequence}
#'  \item{peptide}{fragment of the sequence, result of digesting}
#'  
#' @details TEMPORARY DOCS
#' 
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
         peptides = construct_sq(unlist(detected), "ami", TRUE))
}