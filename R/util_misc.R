.get_sq_subclass <- function(sq) {
  intersect(class(sq), c("amisq", "nucsq", "dnasq", "rnasq", "untsq", "atpsq", "encsq"))
}

.get_sq_type <- function(sq) {
  sqclasses <- intersect(class(sq), c("amisq", "dnasq", "rnasq", "untsq", "atpsq", "encsq"))
  dict <- c(amisq = "ami", dnasq = "dna", rnasq = "rna", untsq = "unt", atpsq = "atp", encsq = "enc")
  dict[sqclasses]
}

.is_cleaned <- function(sq) {
  "clnsq" %in% class(sq)
}

.set_original_length <- function(sq, orig_lengths) {
  if (length(sq) == 0) return(sq)
  for (index in 1:length(sq)) {
    attr(sq[[index]], "original_length") <- orig_lengths[index]
  }
  sq
}

# TODO: verify if all calls in code don't pass "list" in classes vector
.construct_sq_s <- function(sq, alph, classes) {
  sq <- .pack_to_sq(sq, alph)
  new_list_of(sq,
              ptype = raw(),
              alphabet = alph,
              class = classes)
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
  if (test_file_exists(file)) {
    normalizePath(file)
  } else {
    tmp <- tempfile()
    download.file(file, tmp)
    tmp
  }
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