#' @export
read_fasta <- function(file, type = "unt") {
  if (!is.character(file) ||
      !(length(file) == 1)) {
    stop("'file' has to be a string giving file to read from")
  }
  if (!file.exists(file)) {
    stop("'file' doesn't exists")
  }
  if (!(type %in% c("ami", "nuc", "unt"))) {
    stop("'type' needs to be one of 'ami', 'nuc', 'unt' ('unt' is default)")
  }
  
  #used from biogram
  all_lines <- readLines(file)
  prot_id <- cumsum(grepl("^>", all_lines))
  all_prots <- split(all_lines, prot_id)
  
  sq_list <- sapply(all_prots, function(prot) prot[-1])
  
  names_vec <- sub(">", "", sapply(all_prots, function(prot) prot[1]), fixed = TRUE)
  
  construct_sqtibble(sq_list, names_vec, type)
}
