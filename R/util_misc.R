.get_alph_size <- function(alph) {
  ceiling(log2(length(alph) + 2))
}

.get_na_val <- function(alph) {
  2 ^ .get_alph_size(alph) - 1
}

.get_alph <- function(sq) {
  attr(sq, "alphabet")
}

.get_real_alph <- function(sq) {
  unique(unlist(strsplit(sq, "")))
}

.get_sq_subclass <- function(sq) {
  intersect(class(sq), c("amisq", "nucsq", "dnasq", "rnasq", "untsq", "atpsq", "encsq"))
}

.get_sq_type <- function(sq) {
  sqclasses <- intersect(class(sq), c("amisq", "dnasq", "rnasq", "untsq", "atpsq", "encsq"))
  dict <- c(amisq = "ami", dnasq = "dna", rnasq = "rna", untsq = "unt", atpsq = "atp", encsq = "enc")
  dict[sqclasses]
}

.get_standard_alph <- function(type, is_clean) {
       if (type == "ami" &&  is_clean) aminoacids_df[!aminoacids_df[["amb"]], "one"]
  else if (type == "ami" && !is_clean) aminoacids_df[, "one"]
  else if (type == "dna" &&  is_clean) nucleotides_df[nucleotides_df[["dna"]], "one"]
  else if (type == "dna" && !is_clean) nucleotides_df[nucleotides_df[["dna"]] | nucleotides_df[["amb"]], "one"]
  else if (type == "rna" &&  is_clean) nucleotides_df[nucleotides_df[["rna"]], "one"]
  else if (type == "rna" && !is_clean) nucleotides_df[nucleotides_df[["rna"]] | nucleotides_df[["amb"]], "one"]
}

.is_cleaned <- function(sq) {
  "clnsq" %in% class(sq)
}

.set_class <- function(sq, type, is_clean = FALSE) {
  class(sq) <- c(if (is_clean) "clnsq" else NULL, paste0(type, "sq"), "sq", "list")
  sq
}

.set_alph <- function(sq, alph) {
  attr(sq, "alphabet") <- alph
  sq
}

.set_original_length <- function(sq, orig_lengths) {
  for (index in 1:length(sq)) {
    attr(sq[[index]], "original_length") <- orig_lengths[index]
  }
  sq
}

.set_class_alph <- function(new_sq, sq) {
  class(new_sq) <- class(sq)
  attr(new_sq, "alphabet") <- .get_alph(sq)
  new_sq
}

.construct_sq_s <- function(sq, alph, classes) {
  sq <- .pack_to_sq(sq, alph)
  attr(sq, "alphabet") <- alph
  class(sq) <- classes
  sq
}

.guess_ami_is_clean <- function(real_alph) {
  if (all(real_alph %in% .get_standard_alph("ami", TRUE)))
    TRUE
  else if (all(real_alph %in% .get_standard_alph("ami", FALSE)))
    FALSE
  else stop("there are letters that aren't in IUPAC standard! (see: aminoacids_df)")
}

.guess_dna_is_clean <- function(real_alph) {
  if (all(real_alph %in% .get_standard_alph("dna", TRUE)))
    TRUE
  else if (all(real_alph %in% .get_standard_alph("dna", FALSE)))
    FALSE
  else stop("there are letters that aren't in IUPAC standard! (see: nucleotides_df)")
}

.guess_rna_is_clean <- function(real_alph) {
  if (all(real_alph %in% .get_standard_alph("rna", TRUE)))
    TRUE
  else if (all(real_alph %in% .get_standard_alph("rna", FALSE)))
    FALSE
  else stop("there are letters that aren't in IUPAC standard! (see: nucleotides_df)")
}

.guess_sq_type_subtype <- function(sq) {
  real_alph <- toupper(.get_real_alph(sq))
  .guess_type_subtype_by_alph(real_alph)
}

.guess_type_subtype_by_alph <- function(alph) {
  # TODO: any better idea for it?
  possib_rets <- list(list(type = "dna", is_clean = TRUE),
                      list(type = "rna", is_clean = TRUE),
                      list(type = "ami", is_clean = TRUE),
                      list(type = "dna", is_clean = FALSE),
                      list(type = "rna", is_clean = FALSE),
                      list(type = "ami", is_clean = FALSE))
  for (ret in possib_rets)
    if (all(alph %in% .get_standard_alph(ret[["type"]], ret[["is_clean"]]))) return(ret)
  list(type = "unt", is_clean = NULL)  
}

# TODO: check what it does
.merge_ind <- function(res_ind, begs) {
  n <- length(res_ind)
  m <- length(begs)
  ret <- logical(n + m)
  act_res <- 1
  act_beg <- 1
  act_out <- 1
  while (act_res <= n &&
         act_beg <= m) {
    if (res_ind[act_res] < begs[act_beg]) {
      ret[act_out] <- TRUE
      act_res = act_res + 1
    } else {
      ret[act_out] <- FALSE
      act_beg = act_beg + 1
    }
    act_out = act_out + 1
  }
  ret[act_out:(n+m)] <- (act_res <= n)
  ret
}

.get_readable_file <- function(file) {
  if (!file.exists(file)) {
    tmp <- tempfile()
    download.file(file, tmp)
    tmp
  } else normalizePath(file)
}

#' @importFrom stringi stri_replace_all_regex stri_match_all_regex stri_split_regex
.regexify_pattern <- function(digest_pattern) {
  ami_alph <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", 
                "P", "Q", "R", "S", "T", "V", "W", "Y")
  negated <- stri_match_all_regex(digest_pattern, "(?<=\\<)[ACDEFGHIKLMNPQRSTVWY]+(?=\\>)")[[1]]
  pattern <- stri_split_regex(digest_pattern, "\\<[ACDEFGHIKLMNPQRSTVWY]+\\>")[[1]]
  if (length(negated) == 1 && is.na(negated[1])) {
    negated <- NULL
  } else {
    negated <- paste0("[^", sapply(strsplit(as.character(negated), ""), 
                                   function(neg_set) paste(base::setdiff(ami_alph, neg_set), collapse = "")), "]")
  }
  pattern <- paste0(pattern, c(negated, ""), collapse = "")
  
  sides <- stri_split_regex(pattern, "\\.")[[1]]
  sides[[1]] <- ifelse(nchar(sides[[1]]) > 1, paste0("(?<=", sides[[1]], ")"), "")
  sides[[2]] <- ifelse(nchar(sides[[2]]) > 1, paste0("(?=", sides[[2]], ")"), "")

  paste0(sides[[1]], sides[[2]])
}