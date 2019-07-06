#' @export
count_kmers <- function(sqtbl, lens, dsts, position = FALSE) {
  validate_sqtibble(sqtbl)
  if (nrow(sqtbl) == 0) {
    stop("you cannot count kmers in sqtibble with less than one row")
  }
  sqtypes <- extract_sq_types(sqtbl)
  if (length(unique(sqtypes)) != 1 ||
      any(sqtypes == "unt")) {
    stop("'sqtbl' should contain only sequences of the same type - either 'aa', 'nuc' or 'sim'")
  } 
  
  if (length(unique(lengths(sqtbl[["sq"]]))) != 1 ||
      any(sqtypes == "unt")) {
    stop("'sqtbl' should have at least one sequence")
  } 
  if (!sqtypes[1] == "sim" &&
      !all(extract_is_clean(sqtbl))) {
    stop("'aa' or 'nuc' sequences should be cleaned")
  }
  if (sqtypes[1] == "sim" &&
      !(length(unique(lapply(sqtbl[["sq"]], function(sq) levels(sq))))) == 1) {
    stop("some sequences have other alphabets than others")
  }
  
  
  if (!is.numeric(lens) ||
      floor(lens) != lens ||
      any(lens <= 0)) {
    stop("'lens' has to be an integer vector with positive values")
  }
  if (!is.list(dsts) ||
      !all(sapply(dsts, function(dst) is.numeric(dst))) ||
      any(sapply(dsts, function(dst) any(floor(dst) != dst)))) {
    stop("'dsts' has to be a list of integer vector with non-negative values")
  }
  if (length(lens) != length(dsts)) {
    stop("'dsts' has to have equal length to 'lens'")
  }
  if (!all(lens > lengths(dsts) | 
          (lens == 1 & sapply(dsts, function(dst) length(dst) == 1 && dst ==0)))) {
    stop("each of 'dsts' element should have length less than value of corresponding element of 'lens' or be equal to 1")
  }
  
  if (!is.logical(position) ||
      is.na(position)) {
    stop("'position' should be either TRUE or FALSE")
  }
  
  alphabet <- levels(sqtbl[["sq"]][[1]])
  
  sqmatrix <- as.matrix(do.call(rbind, lapply(sqtbl[["sq"]], function(sq) as.character(sq))))
  
  n_patterns <- length(lens)
  ret <- do.call(cbind, lapply(1:n_patterns, function(i) {
    count_pattern_kmers(sqmatrix, lens[i], dsts[[i]], position, alphabet)
  }))
  # alternatively:
  # d_count_all(sqmatrix, lens, dsts, position, alphabet)
  ret
}

#' @import slam
count_pattern_kmers <- function(sqmatrix, len, dst, position, alphabet) {
  len_sq <- ncol(sqmatrix)
  num_sq <- nrow(sqmatrix)
  possib_ngrams <- create_ngrams(len, alphabet)
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
    pos_possib_ngrams <- create_ngrams(len, alphabet, max_grams)
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
create_ngrams <- function(len, alphabet, possible_grams = NULL) {
  grid_list <- lapply(1:len, function(i) alphabet)
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

#following functions are not yet ready to be used !

d_count_all <- function(sqmatrix, lens, dsts, position, alphabet) {
  #we assume here, that position of dsts and lens is ordered
  nlens <- lens + sapply(dsts, sum)
  dsts_by_ngrams_lens <- split(dsts, nlens)
  
  ngrams_by_nlens <- d_create_bgrams_structure(lens, dsts_by_nlens, alphabet)
  #...
}

d_create_ngrams_structure <- function(dsts, alphabet) {
  pat_strct <- d_create_patterns_structure(dsts)
  
  max_ng_len <- length(pat_strct[["pats"]])
  
  ngrams_strct <- vector("list", max_ng_len)
  nonvirt_ind <- pat_strct[["nvi"]]
  
  ngrams_strct[[1]] <- alphabet
  
  
  for (i in 2:max_ng_len) {
    pats_i <- pat_strct[["pats"]][[i]]
    for (j in 1:length(ngrams_strct[[i-1]])) {
      for (k in 1:length(pats_i)) {
        if (substr(pats_i, i, i) == "_") {
          ngrams_strct[[i]] <- c(ngrams_strct[[i]], paste0(ngrams_strct[[i-1]][j], "_"))
        } else {
          #...
        }
      }
    }
  }
}

d_create_patterns_structure <- function(dsts) {
  #we assume, that lens and dsts are correct
  ng_lens <- sapply(dsts, sum) + lengths(dsts) + 1
  max_ng_len <- max(ng_lens)
  
  pats_by_lens <- lapply(1:max_ng_len,
         function(i) {
           lapply(dsts[ng_lens == i], d_create_single_pattern)
         })
  nonvirt_ind <- lengths(pats_by_lens)
  #kids_refs <- vector("list", max_ng_len-1)
  parents_refs <- vector("list", max_ng_len-1)
  
  for (i in max_ng_len:2) {
    parent_pats <- lapply(pats_by_lens[[i]], function(pat) substr(pat, 1, i-1))
    parent_inds <- match(parent_pats, pats_by_lens[[i-1]])
    pats_by_lens[[i-1]] <- c(pats_by_lens[[i-1]], unique(parent_pats[is.na(parent_inds)]))
    parent_inds <- match(parent_pats, pats_by_lens[[i-1]])
    n <- length(parent_inds)
    
    parents_refs[[i]] <- parent_inds
    
    #unique(parent_inds)[unique(parent_inds)] <- lapply(unique(parent_inds), function(ind) (1:n)[parent_inds == ind])
  }
  
  list(pats = pats_by_lens, refs = parents_refs, nvi = nonvirt_ind)
}

d_create_single_pattern <- function(dsts) {
  if (length(dsts) == 0) {
    "#"
  } else {
    n <- length(dsts)
    ret <- rep("_", sum(dsts) + n + 1)
    ret[c(1, cumsum(dsts) + 2:(n + 1))] <- "#"
    paste(ret, collapse = "")
  }
}

