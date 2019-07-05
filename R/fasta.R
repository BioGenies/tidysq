#' @export
read_fasta <- function(file, type = "unt") {
  if (!is.character(file) ||
      !(length(file) == 1)) {
    stop("'file' has to be a string giving file to read from")
  }
  if (!file.exists(file)) {
    stop("'file' doesn't exists")
  }
  if (!(type %in% c("aa", "nuc", "unt"))) {
    stop("'type' needs to be one of 'aa', 'nuc', 'unt' (default)")
  }
  
  #used from biogram
  all_lines <- readLines(file)
  prot_id <- cumsum(grepl("^>", all_lines))
  all_prots <- split(all_lines, prot_id)
  
  sq_list <- lapply(all_prots, function(prot)
    unlist(strsplit(prot[-1], split = "")))
  
  names_vec <- sub(">", "", sapply(all_prots, function(prot) prot[1]), fixed = TRUE)
  
  sq_list <- lapply(sq_list, function(sq) construct_sq(sq, type))
  construct_sqtibble(names_vec, unname(sq_list))
}
