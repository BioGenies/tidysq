#' Count k-mers in sequences
#'
#' Counts all k-mers or position-specific k-mers present in the input sequence(s).
#'
#' @param sq a vector or matrix describing sequence(s). 
#' @param lens \code{integer} vector of lengths of k-mers.
#' @param dists \code{integer} vector of distances between elements of k-mer (0 means 
#' consecutive elements). See Details.
#' @param alph \code{integer}, \code{numeric} or \code{character} vector of all
#' possible elements of the alphabet.
#' @param position \code{logical}, if \code{TRUE} position-specific k-mers are counted.
#' 
#' @details \code{ns} vector and \code{ds} vector must have equal length. Elements of 
#' \code{ds} vector are used as equivalents of \code{d} parameter for respective values 
#' of \code{ns}. For example, if \code{ns} is \code{c(4, 4, 4)}, the \code{ds} must be a list of 
#' length 3. Each element of the \code{ds} list must have length 3 or 1, as appropriate
#' for a \code{d} parameter in \code{count_ngrams} function.
#' 
#' A \code{dists} vector should be always \code{n} - 1 in length.
#' For example when \code{n} = 3, \code{d} = c(1,2) means A_A__A. For \code{n} = 4, 
#' \code{d} = c(2,0,1) means A__AA_A. If vector \code{d} has length 1, it is recycled to
#' length \code{n} - 1.
#' 
#' k-mer names follow a specific convention and have three parts for position-specific
#' k-mers and two parts otherwise. The parts are separated by \code{_}. The \code{.} symbol
#' is used to separate elements within a part. The general naming scheme is 
#' \code{POSITION_KMER_DISTANCE}. The optional \code{POSITION} part of the name indicates
#' the actual position of the k-mer in the sequence(s) and will be present 
#' only if \code{pos} = \code{TRUE}. This part is always a single integer. The \code{KMER}
#' part of the name is a sequence of elements in the k-mer. For example, \code{4.2.2}
#' indicates the k-mer 422 (e.g. TCC). The \code{DISTANCE} part of the name is a vector of
#' distance(s). For example, \code{0.0} indicates zero distances (continuous k-mers), while
#' \code{1.2} represents distances for the k-mer A_A__A.
#' 
#' Examples of k-mer names:
#' \itemize{
#' \item{46_4.4.4_0.1 : trigram 44_4 on position 46}
#' \item{12_2.1_2     : bigram 2__1 on position 12}
#' \item{8_1.1.1_0.0  : continuous trigram 111 on position 8}
#' \item{1.1.1_0.0    : continuous trigram 111 without position information}
#' }
#' 
#' @return a \code{\link[slam]{simple_triplet_matrix}} where columns represent
#' k-mers and rows sequences. See \code{Details} for specifics of the naming convention.
#' 
#' @note By default, the counted k-mer data is stored in a memory-saving format.
#' To convert an object to a 'classical' matrix use the \code{\link[base]{as.matrix}}
#' function. See examples for further information.
#' @export
#' @seealso 
#' Extract k-mers from sequence(s): \code{\link{extract_kmers}}.
#' @examples 
#' set.seed(1410)
#' random_sqs <- random_sq(5, 20, type = "nuc", is_clean = TRUE)
#' 
#' # count 1-mers without position information for nucleotides
#' count_kmers(random_sqs, 1, 0, alph = c("A", "C", "G", "T"), pos = FALSE)
#' 
#' # count position-specific 3-mers from multiple nucleotide sequences
#' kmers <- count_kmers(random_sqs, 3, c(0, 0, 0), 
#'                      alph = c("A", "C", "G", "T"), pos = TRUE)
#' 
#' # output results of the k-mer counting to screen
#' as.matrix(kmers)
#' 
#' # count multiple k-mers of different lengths and distances
#' count_kmers(random_sqs, c(1, 2, 2), list(0, 0, 1), 
#'             alph = c("A", "C", "G", "T"), pos = FALSE)
#' @export
count_kmers <- function(sq, lens = 0, dists = list(0), alph = NULL, position = FALSE) {
  validate_sq(sq)
  .check_isnt_null(dists, "'dists'")
  if (is.list(dists)) 
    .check_list_dists(dists)
  else  {
    .check_integer(dists, "'dists'", allow_zero_len = TRUE, allow_zero = TRUE)
    dists <- list(dists)
  }
  .check_dists_prop_len(sq, dists)
  
  if(length(lens) != length(dists)) {
    stop("lens vector and distances vector must have equal length.", call. = FALSE)
  }
  
  .check_isnt_missing(alph, "'alph'")
  if (is.null(alph)) alph <- .get_alph(sq)
  else .check_alph_is_subset(sq, alph)
  .check_logical(position, "'position'", single_elem = TRUE)
  sqmatrix <- as.matrix(sq)
  n_patterns <- length(dists)

  do.call(cbind, lapply(1:n_patterns, function(i) {
    count_pattern_kmers(sqmatrix, lens[i], dists[[i]], position, alph)
  }))
}

#' Extract k-mers from sequences
#'
#' Extracts vector of k-mers present in sequence(s).
#'
#' @inheritParams count_kmers
#' @details A format of \code{d} vector is discussed in Details of 
#' \code{\link{count_kmers}}.
#' @return A \code{character} matrix of k-mers, where every row corresponds to a
#' different sequence.
#' @examples 
#' set.seed(1410)
#' random_sqs <- random_sq(3, 5, type = "nuc", is_clean = TRUE)
#' 
#' # 1-mers
#' extract_kmers(random_sqs, lens = 1, dists = list(0))
#' 
#' # 2-mers with and without gap
#' extract_kmers(random_sqs, lens = c(2, 2), dists = list(0, 1))
#' @export
extract_kmers <- function(sq, lens = 0, dists = list(0)) {
  validate_sq(sq)
  #.check_integer(dists, "'dists'", allow_zero_len = TRUE, allow_zero = TRUE)
  #.check_dists_prop_len(sq, dists)
  sq_lens <- .get_lens(sq)
  sqmatrix <- as.matrix(sq)
  
  do.call(rbind, lapply(1L:nrow(sqmatrix), function(ith_seq) {
    ngram_ind <- lapply(1L:length(lens), function(i) 
      get_ngrams_ind(len_seq = sq_lens[ith_seq], len = lens[i], dst = dists[[i]]))
    
    max_grams <- lapply(1L:length(lens), function(i) 
      calc_max_grams(len_seq = sq_lens[ith_seq], len = lens[i], ngram_ind = ngram_ind[[i]]))
    
    
    kmers_vector <- do.call(c, lapply(1L:length(lens), function(i) {
      grams <- na.omit(seq2ngrams_helper(sqmatrix[ith_seq, ], ind = ngram_ind[[i]], max_grams[i]))
      paste(grams, paste0(attr(ngram_ind[[i]], "d"), collapse = "."),
            sep = "_")
    }))
    # TODO: add position and the kmer in the sq format
    tibble(origin_sq = sq[ith_seq], kmer = kmers_vector)
  }))
} 

#' @import slam
count_pattern_kmers <- function(sqmatrix, len, dst, position, alph) {
  len_sq <- ncol(sqmatrix)
  num_sq <- nrow(sqmatrix)

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
    stop("k-mer too long.")
  }
  max_grams
}

calc_max_grams <- function(len_seq, len, ngram_ind){
  # use attr(ngram_ind, "d") instead of d because of distance recycling
  max_grams <- len_seq - len - sum(attr(ngram_ind, "d")) + 1
  if (max_grams < 1)
    stop("k-mer too long.")
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