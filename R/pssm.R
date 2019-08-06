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

score_pssm <- function(sq, pssm) {}

scale_pssm <- function(pssm) {}

flatten_pssm <- function(pssm) {}