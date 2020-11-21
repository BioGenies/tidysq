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
