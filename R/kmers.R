#' @export
count_kmers <- function(sq, dists, alph = NULL, position = FALSE) {
  validate_sq(sq)
  .check_isnt_missing(dists, "'dists'")
  .check_isnt_null(dists, "'dists'")
  if (is.list(dists)) 
    .check_list_dists(dists)
  else  {
    .check_integer(dists, "'dists'", allow_zero_len = TRUE, allow_zero = TRUE)
    dists <- list(dists)
  }
  .check_dists_prop_len(sq, dists)
  .check_isnt_missing(alph, "'alph'")
  if (is.null(alph)) alph <- .get_alph(sq)
  else .check_alph_is_subset(sq, alph)
  .check_logical(position, "'position'", single_elem = TRUE)
  sqmatrix <- as.matrix(sq)
  n_patterns <- length(dists)
  do.call(cbind, lapply(1:n_patterns, function(i) {
    count_pattern_kmers(sqmatrix, dists[[i]], position, alph)
  }))
}

#' @export
extract_kmers <- function(sq, dists) {
  validate_sq(sq)
  .check_integer(dists, "'dists'", allow_zero_len = TRUE, allow_zero = TRUE)
  .check_dists_prop_len(sq, dists)
  
  lens <- .get_lens(sq)
  n <- length(sq)
  kmers_len <- sum(dists) + length(dists) + 1
  
  sqmatrix <- as.matrix(sq)
  
  kmers <- do.call(c, lapply(1:n, function(i) {
    do.call(c, lapply(1:(lens[i] - kmers_len + 1), function(ind) 
      paste(sqmatrix[i, ind:(ind + kmers_len - 1)], collapse = "")))
  }))
  
  ret_sq <- sq[do.call(c, lapply(1:n, function(i) rep(i, lens[i] - kmers_len + 1)))]
  
  ngram_ind <- lapply(lens, get_ngrams_ind, len = kmers_len, dst = dists)
  max_grams <- sapply(1:n, function(i) calc_max_grams(lens[i], kmers_len, ngram_ind[i]))
  
  ret_ids <- do.call(c, lapply(1:n, function(i) {
    grams <- na.omit(seq2ngrams_helper(sqmatrix[i, ], ind = ngram_ind[[i]], max_grams[i]))
    paste(grams, paste0(attr(ngram_ind[[i]], "d"), collapse = "."),
          sep = "_")
  }))
  
  tibble(sq = ret_sq, kmer_id = ret_ids, kmer = construct_sq(kmers))
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