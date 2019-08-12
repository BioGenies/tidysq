.create_empty_pssm <- function(num_pos, type) {
  alph <- setdiff(.get_standard_alph(type, TRUE), "*")
  matrix(0, num_pos, length(alph), dimnames = list(1:num_pos, alph))
}

.compute_pssm_counts <- function(sq) {
  num_pos <- max(.get_lens(sq))
  alph <- setdiff(.get_alph(sq), "*")
  ret <- .create_empty_pssm(num_pos, .get_sq_type(sq))
  
  sq <- as.matrix(sq)
  .check_matrix_no_na(sq)
  .check_matrix_no_star(sq)
  
  for (let in alph) {
    ret[, let] <- colSums(sq == let)
  }
  ret
}

.compute_pssm_freqs <- function(sq) {
  pssm <- .compute_pssm_counts(sq)
  pssm / rowSums(pssm)
}

.compute_pssm_shannon <- function(sq) {
  pssm <- .compute_pssm_freqs(sq)
  ifelse(pssm == 0, 0, -pssm * log2(pssm))
}

.compute_pssm_kl <- function(sq, background_dist = NULL) {
  pssm <- .compute_pssm_freqs(sq)
  
  if (is.null(background_dist)) 
    q <- 1 / if (.get_sq_type(sq) == "ami") 20 else 4
  else {
    .check_back_dist_is_proper(background_dist)
    background_dist <- c(as.numeric(BGFREQS[background_dist, ]), 0)
    q <- matrix(rep(background_dist, each = nrow(pssm)), nrow(pssm), ncol(pssm))
  }
  ifelse(pssm == 0, 0, pssm * log2(pssm / q))
}