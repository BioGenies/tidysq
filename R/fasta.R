#' @export
read_fasta <- function(file, type = "amb") {
  all_lines <- readLines(file)
  prot_id <- cumsum(grepl("^>", all_lines))
  all_prots <- split(all_lines, prot_id)
  
  sq_list <- lapply(all_prots, function(ith_seq)
    unlist(strsplit(ith_seq[-1], split = "")))
  
  names_vec <- sub(">", "", sapply(all_prots, function(ith_seq) ith_seq[1]), fixed = TRUE)
  
  # Return output sqtibble
  sq_list <- lapply(sq_list, function(sq) construct_sq(sq, type))
  construct_sqtibble(names_vec, sq_list)
}
