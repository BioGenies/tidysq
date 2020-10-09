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
