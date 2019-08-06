#' @export
empty_pssm <- function(num_pos, type) {
  .check_type_in_ami_nuc(type)
  .check_num_pos_proper_int(num_pos)
  .create_empty_pssm(num_pos, type)
}


#' @export
compute_pssm <- function(sq, mode, background_dist = "All") {
  validate_sq(sq)
  .check_sq_type_in_ami_nuc(.get_sq_type(sq))
  .check_sq_is_clean(sq)
  .check_sq_lens_eq(sq)
  .check_mode_is_proper(mode)
  
  switch(mode,
         counts = .compute_pssm_counts(sq),
         freqs = .compute_pssm_freqs(sq),
         shannon = .compute_pssm_shannon(sq),
         information = .compute_pssm_kl(sq),
         `kullback-leibler` = .compute_pssm_kl(sq, background_dist))
}

#' @export
score_pssm <- function(sq, pssm) {
  validate_sq(sq)
  .check_is_pssm(pssm)
  .check_sq_type_in_ami_nuc(.get_sq_type(sq))
  .check_types_match(sq, pssm)
  .check_sq_lens_eq(sq)
  len <- .get_lens(sq[1])
  n <- length(sq)
  .check_lens_eq_num_pos(len, nrow(pssm))
  sq <- as.matrix(sq)
  .check_matrix_no_na(sq)
  .check_matrix_no_star(sq)
  sq_long  <- as.character(sq)
  pssm_long <- pssm[, sq_long]
  s_vec <- colSums(pssm_long * matrix(rep(diag(len), n), nrow = len))
  s_mat <- matrix(s_vec, ncol = len, byrow = TRUE)
  rowSums(s_mat)
}

#' @export
scale_pssm <- function(pssm) {
  .check_is_pssm_or_numeric(pssm)
  (pssm - min(pssm, na.rm = TRUE)) / (max(pssm, na.rm = TRUE) - min(pssm, na.rm = TRUE))
}

#' @export
flatten_pssm <- function(pssm) {
  .check_is_pssm(pssm)
  out_col_names <- sapply(colnames(pssm), 
                          function(name) paste(rownames(pssm), name, sep = "_")) 
  ret <- as.numeric(pssm)
  names(ret) <- out_col_names
  ret
  
}