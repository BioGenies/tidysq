#' @importFrom stringi stri_sub stri_locate_all_regex stri_count_regex
#' @export
find_motifs <- function(sq, name, motifs) {
  validate_sq(sq)
  if (missing(name) ||
      !is.character(name) ||
      any(is.null(name) | is.na(name)) ||
      (length(name) != length(sq))) {
    stop("'name' has to be a non-NULL character vector without NA's, of lenght equal to length of sq")
  }
  
  if (!is.character(motifs) ||
      any(is.na(motifs)) ||
      length(motifs) < 1) {
    stop("'motifs' needs to be a non-NULL character vector without NA's")
  }
  
  type <- .get_sq_type(sq)
  
  sq_c <- sq
  motifs_c <- motifs
  motifs_l <- nchar(motifs) - stri_count_regex(motifs, "[\\\\$]")
  motifs <- ifelse(motifs_l == 1,
                   motifs, 
                   paste0(stri_sub(motifs, 1, 1), "(?=", stri_sub(motifs, 2), ")"))
  motifs <- strsplit(motifs, "")
  
  alph <- .get_alph(sq)
  
  if (type == "ami") {
    motifs <- lapply(motifs, toupper)
    if (!all(unlist(strsplit(motifs_c, "")) %in% c(.get_standard_alph("ami", FALSE), "^", "$"))) {
      stop("motifs that you're searching for in 'sq' object needs to consist of letters from aminoacids alphabet and optionally '^' or '$' characters")
    }
    
    motifs <- lapply(motifs, function(motif) replace(motif, motif == "*", "\\*"))
    motifs <- lapply(motifs, function(motif) replace(motif, motif == "B", "[BDN]"))
    motifs <- lapply(motifs, function(motif) replace(motif, motif == "J", "[JIL]"))
    motifs <- lapply(motifs, function(motif) replace(motif, motif == "X", "[A-Z]"))
    motifs <- lapply(motifs, function(motif) replace(motif, motif == "Z", "[ZEQ]"))
    motifs <- sapply(motifs, function(motif) paste(motif, collapse = ""))
  } else if (type == "nuc") {
    motifs <- lapply(motifs, toupper)
    if (!all(unlist(strsplit(motifs_c, "")) %in% c(.get_standard_alph("nuc", FALSE), "^", "$"))) {
      stop("motifs that you're searching for in 'sq' object needs to consist of letters from nucleotides alphabet and optionally '^' or '$' characters")
    }
    motifs <- lapply(motifs, function(motif) replace(motif, motif == "W", "[WATU]"))
    motifs <- lapply(motifs, function(motif) replace(motif, motif == "S", "[SCG]"))
    
    motifs <- lapply(motifs, function(motif) replace(motif, motif == "M", "[MAC]"))
    motifs <- lapply(motifs, function(motif) replace(motif, motif == "K", "[KGTU]"))
    motifs <- lapply(motifs, function(motif) replace(motif, motif == "R", "[RAG]"))
    motifs <- lapply(motifs, function(motif) replace(motif, motif == "Y", "[YCTU]"))
    
    motifs <- lapply(motifs, function(motif) replace(motif, motif == "B", "[BCTGU]"))
    motifs <- lapply(motifs, function(motif) replace(motif, motif == "D", "[DATGU]"))
    motifs <- lapply(motifs, function(motif) replace(motif, motif == "H", "[HACTU]"))
    motifs <- lapply(motifs, function(motif) replace(motif, motif == "V", "[VACG]"))
    
    motifs <- lapply(motifs, function(motif) replace(motif, motif == "N", "[ACTGUWSMKRYBDHVN]"))
    
    motifs <- sapply(motifs, function(motif) paste(motif, collapse = ""))
    
  } else if (any(alph %in% c("^", "$", "?", "(", "=", ")", "\\", ".", "|", "+", "*", "{", "}", "[", "]"))) {
    stop("you cannot search for motifs if any of those characters: ^$?=()\\.|+*{}[] are elements of 'sq' alphabet; if you want to use them, please substitute those letters with some other using 'substitute_letters'")
  }
  
  
  sq <- as.character(sq)
  
  n <- length(sq)
  m <- length(motifs)
  
  match_inds <- lapply(1:m, function(i) {
    ret <- stri_locate_all_regex(sq, motifs[i])
    ret <- cbind(do.call(rbind, ret), s_ind = unlist(lapply(1:n, function(i) rep(i, nrow(ret[[i]])))))
    ret <- ret[!is.na(ret[,"start"]), , drop = FALSE]
    ret[, "end"] <- ret[, "end"] + rep(motifs_l[i], nrow(ret)) - 1
    ret
  })
  
  matched_inds <- cbind(do.call(rbind, match_inds)) 
  sought <-  unlist(lapply(1:m, function(j) rep(motifs_c[j], nrow(match_inds[[j]]))))
  
  sq_col <- sq_c[matched_inds[, "s_ind"]]
  nm_col <- name[matched_inds[, "s_ind"]]
  found <- stri_sub(sq[matched_inds[, "s_ind"]], from = matched_inds[, "start"], to = matched_inds[, "end"])
  
  sq_col <- .set_class_alph(sq_col, sq_c)
  
  tibble(name = nm_col, 
         sq = sq_col, 
         sought = sought, 
         found = .construct_sq_s(found, alph, class(sq_c)), 
         start = matched_inds[, "start"], 
         end = matched_inds[, "end"])
}