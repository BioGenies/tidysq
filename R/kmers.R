#' @export
count_kmers <- function(sq, dsts, position = FALSE) {
  validate_sq(sq)
  
  if (!is.list(dsts) ||
      !all(sapply(dsts, function(dst) is.numeric(dst))) ||
      any(sapply(dsts, function(dst) any(floor(dst) != dst)))) {
    stop("'dsts' has to be a list of integer vector with non-negative values")
  }
  
  if (!is.logical(position) ||
      is.na(position)) {
    stop("'position' should be either TRUE or FALSE")
  }
  
  alph <- .get_alph(sq)
  
  sqmatrix <- as.matrix(sq)
  
  n_patterns <- length(dsts)
  do.call(cbind, lapply(1:n_patterns, function(i) {
    count_pattern_kmers(sqmatrix, dsts[[i]], position, alph)
  }))
}

#' @export
extract_kmers <- function(sq, len, dsts, position = FALSE) {
  validate_sq(sq)
  .check_sq_lens_eq(sq)
  .check_len_is_int(len)
  #.check_dsts_is_numeric(dsts)
  #.check_dsts_aprop_length(dsts)
  
  if (!is.logical(position) ||
      is.na(position)) {
    stop("'position' should be either TRUE or FALSE")
  }
  
  sqmatrix <- as.matrix(sq)
  # length of sequence
  len_seq <- ncol(sqmatrix)
  # number of sequences
  n_seqs <- nrow(sqmatrix)
  
  # look for n-gram indices for d
  ngram_ind <- get_ngrams_ind(len_seq, len, dsts)
  
  max_grams <- calc_max_grams(len_seq, len, ngram_ind)
  
  # extract n-grams from sequene
  res <- t(vapply(1L:n_seqs, function(i) {
    grams <- seq2ngrams_helper(sqmatrix[i, ], ind = ngram_ind, max_grams)
    paste(grams, paste0(attr(ngram_ind, "d"), collapse = "."), 
          sep = "_")
  }, rep("a", max_grams)))
  if (max_grams == 1)
    res <- t(res)
  
  # add position information if requested
  if (position)
    res <- do.call(cbind, lapply(1L:ncol(res), function(pos_id)
      paste0(pos_id, "_", res[, pos_id])))
  
  res
} 

#' @import slam
count_pattern_kmers <- function(sqmatrix, dst, position, alph) {
  len_sq <- ncol(sqmatrix)
  num_sq <- nrow(sqmatrix)
  len <- sum(dst) + length(dst) + 1
  possib_ngrams <- create_ngrams(len, alph)
  ngram_ind <- get_ngrams_ind(len_sq, len, dst)
  max_grams <- calc_max_ngrams(len_sq, len, ngram_ind)
  grams <- vapply(1L:num_sq, 
                  function(i) seq2ngrams_helper(sqmatrix[i, ], 
                                                ind = ngram_ind, 
                                                max_grams), 
                  rep("a", max_grams))
  if (!is.matrix(grams)) {
    grams <- matrix(grams, ncol = num_sq)
  }
  if (position) {
    pos_possib_ngrams <- create_ngrams(len, alph, max_grams)
    res <- do.call(cbind, 
                   c(lapply(possib_ngrams,
                            function(ngram) as.simple_triplet_matrix(t(
                              vapply(1L:num_sq,
                                     function(sq) grams[, sq] == ngram, 
                                     rep(0, max_grams)))))))
    colnames(res) <- pos_possib_ngrams
  }
  else {
    res <- do.call(cbind, 
                   lapply(possib_ngrams, 
                          function(ngram) as.simple_triplet_matrix(
                            vapply(1L:num_sq,
                                   function(sq) sum(grams[, sq] == ngram, na.rm = TRUE), 
                                   0))))
    colnames(res) <- possib_ngrams
  }
  colnames(res) <- paste(colnames(res), paste0(attr(ngram_ind, 
                                                    "dst"), collapse = "."), sep = "_")
  res
}
create_ngrams <- function(len, alph, possible_grams = NULL) {
  grid_list <- lapply(1:len, function(i) alph)
  res <- apply(expand.grid(grid_list), 1, function(x) paste(x,
      collapse = "."
    ))
  if (!is.null(possible_grams)) {
    res <- as.vector(sapply(res, function(i) paste(1:possible_grams,
        i,
        sep = "_"
      )))
  }
  res
}
get_ngrams_ind <- function(len_seq, len, dst) {
  ind <- lapply(1:len, function(i) (1 + i - 1):(len_seq - len +
      i))
  if (len > 1) {
    if (length(dst) == 1 && len > 2) {
      dst <- rep(dst, len - 1)
    }
    if (sum(dst) > 0) {
      ind[-1] <- lapply(1L:length(dst), function(i) ind[[i +
          1]] + sum(dst[1L:i]))
      not_taken <- ind[[1]][(length(ind[[1]]) - sum(dst) +
        1):length(ind[[1]])]
      ind <- lapply(ind, function(i) i[-not_taken])
    }
    attr(ind, "dst") <- dst
  }
  else {
    attr(ind, "dst") <- 0
  }
  ind
}
calc_max_ngrams <- function(len_sq, len, ngram_ind) {
  max_grams <- len_sq - len - sum(attr(ngram_ind, "dst")) + 1
  #this check should be moved to validation of input data 
   if (max_grams < 1) {
    stop("n-gram too long.")
  }
  max_grams
}

calc_max_grams <- function(len_seq, len, ngram_ind){
  # use attr(ngram_ind, "d") instead of d because of distance recycling
  max_grams <- len_seq - len - sum(attr(ngram_ind, "d")) + 1
  if (max_grams < 1)
    stop("n-gram too long.")
  max_grams
}
seq2ngrams_helper <- function(sq, ind, max_grams) {
  if (length(ind) > 1) {
    element_matrix <- vapply(ind, function(i) sq[i], rep(
      sq[1],
      max_grams
    ))
    if (max_grams == 1) {
      grams <- paste(element_matrix, collapse = ".")
    }
    else {
      grams <- vapply(1L:nrow(element_matrix), function(ith_row) paste(element_matrix[ith_row, ], collapse = "."), "a")
    }
  }
  else {
    grams <- as.character(sq)
  }
  grams
}