.is_cleaned <- function(x) {
  "clnsq" %in% class(x)
}

.set_original_length <- function(x, orig_lengths) {
  if (length(x) == 0) return(x)
  for (index in 1:length(x)) {
    attr(x[[index]], "original_length") <- orig_lengths[index]
  }
  x
}

# TODO: verify if all calls in code don't pass "list" in classes vector
.construct_sq_s <- function(x, alph, classes) {
  x <- .pack_to_sq(x, alph)
  new_list_of(x,
              ptype = raw(),
              alphabet = alph,
              class = classes)
}

# TODO: outsource
.get_readable_file <- function(file) {
  if (test_file_exists(file)) {
    normalizePath(file)
  } else {
    tmp <- tempfile()
    download.file(file, tmp)
    tmp
  }
}
