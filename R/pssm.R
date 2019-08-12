#' Create empty Position Specific Scoring Matrix
#'
#' @description Create an empty Position Specific Scoring Matrix 
#' for a specified alphabet.
#'
#' @param num_pos the number of positions in a sequence.
#' @param type \code{character} type of a sequence
#' @return A numeric matrix with rows equal to \code{num_pos}
#' and number of columns equal to the number of characters in the specified 
#' alphabet.
#' @examples
#' # empty nucleotide PSSM
#' empty_pssm(num_pos = 9, type = "nuc")
#' 
#' # empty aminoacid PSSM
#' empty_pssm(num_pos = 9, type = "ami")
#' @export
empty_pssm <- function(num_pos, type) {
  .check_type(type)
  .check_integer(num_pos, "'num_pos'", single_elem = TRUE)
  .create_empty_pssm(num_pos, type)
}


#' @export
compute_pssm <- function(sq, mode = "freqs", background_dist = "All") {
  validate_sq(sq)
  .check_type(.get_sq_type(sq), "type of 'sq'")
  .check_is_clean(sq)
  .check_all_elem(length(unique(.get_lens(sq))) != 1, "'sq'", "have to have equal length")
  # add mode checking
  
  switch(mode,
         counts = .compute_pssm_counts(sq),
         freqs = .compute_pssm_freqs(sq),
         shannon = .compute_pssm_shannon(sq),
         information = .compute_pssm_kl(sq),
         `kullback-leibler` = .compute_pssm_kl(sq, background_dist))
}

#' Score a sequence against a PSSM using sum of positional scores.
#'
#' @inheritParams bite
#' @inheritParams flatten_pssm
#' @return A vector of scores, one for each input sequence.
#' @examples
#' pssm1 <- empty_pssm(num_pos = 9, type = "nuc")
#' pssm1[1:9, 1:6] <- rnorm(54)
#' score_pssm(as.sq(c("GGGCTGGCG","CTTCAGACT")), pssm1)
#' @export
score_pssm <- function(sq, pssm) {
  validate_sq(sq)
  #.check_is_pssm(pssm)
  #.check_sq_type_in_ami_nuc(.get_sq_type(sq))
  #.check_types_match(sq, pssm)
  #.check_sq_lens_eq(sq)
  len <- .get_lens(sq[1])
  n <- length(sq)
  #.check_lens_eq_num_pos(len, nrow(pssm))
  sq <- as.matrix(sq)
  #.check_matrix_no_na(sq)
  #.check_matrix_no_star(sq)
  sq_long  <- as.character(t(sq))
  s_vec <- pssm[cbind(rep(1:len, n), sq_long)]
  s_mat <- matrix(s_vec, ncol = len, byrow = TRUE)
  rowSums(s_mat)
}

#' @export
scale_pssm <- function(pssm) {
  #.check_is_pssm_or_numeric(pssm)
  (pssm - min(pssm, na.rm = TRUE)) / (max(pssm, na.rm = TRUE) - min(pssm, na.rm = TRUE))
}

#' 'Flatten' PSSM to a single row
#'
#' 'Flatten' PSSM to a single row, such that e.g. m x n = 9 x 20 becomes 1 x 180
#'
#' @param pssm a position specific scoring matrix
#' @return A numeric matrix with 1 row and m * n columns
#' @examples
#' flatten_pssm(empty_pssm(num_pos = 9, type = "nuc"))
#' @export
flatten_pssm <- function(pssm) {
  #.check_is_pssm(pssm)
  out_col_names <- sapply(colnames(pssm), 
                          function(name) paste(rownames(pssm), name, sep = "_")) 
  ret <- as.numeric(pssm)
  names(ret) <- out_col_names
  ret
  
}