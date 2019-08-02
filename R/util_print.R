#' @importFrom pillar type_sum
#' @exportMethod type_sum amisq
#' @export
type_sum.amisq <- function(x) {
  "ami"
}

#' @importFrom pillar type_sum
#' @exportMethod type_sum nucsq
#' @export
type_sum.nucsq <- function(x) {
  "nuc"
}

#' @importFrom pillar type_sum
#' @exportMethod type_sum untsq
#' @export
type_sum.untsq <- function(x) {
  "unt"
}

#' @importFrom pillar type_sum
#' @exportMethod type_sum atpsq
#' @export
type_sum.atpsq <- function(x) {
  "atp"
}

#' @importFrom pillar type_sum
#' @exportMethod type_sum encsq
#' @export
type_sum.encsq <- function(x) {
  "enc"
}

#' @importFrom pillar type_sum
#' @exportMethod type_sum clnsq
#' @export
type_sum.clnsq <- function(x) {
  paste0("(c)", NextMethod())
}

#' @importFrom crayon blue
#' @importFrom crayon silver
#' @importFrom crayon green
#' @importFrom crayon col_nchar
#' @exportMethod print sq
#' @export
print.sq <- function(x,  
                     max_sequences = getOption("tidysq_max_print_sequences"),
                     use_color = getOption("tidysq_colorful_sq_print"), 
                     letters_sep = NULL) {
  
  alph <- .get_alph(x)
  
  #if parameter is NULL and all letters are lenght one, no space
  if (is.null(letters_sep)) {
    letters_sep <- if (all(nchar(alph) == 1)) "" else " "
  }
  p_width <- getOption("width")
  
  #select at most max_sequences to print
  num_lines <- min(max_sequences, length(x))
  sq <- x[1:num_lines]
  
  #cut sq object so that we don't need to debitify long sequences
  # 6 is minimum lenght of p_lens and p_inds, 8 is byte lenght
  sq_cut <- .cut_sq(sq, ceiling((p_width - 6) / (8 * (nchar(letters_sep) + 1))))
  sq_cut <- .debitify_sq(sq_cut, "char")
  
  #color NA's
  na_char <- .get_na_char()
  if (use_color) sq_cut <- lapply(sq_cut, function(s) {
    s[s == na_char] <- silver(na_char)
    s
  })
  
  #max index number width
  inds_width <- nchar(num_lines) + 2
  
  #indices to print
  p_inds <- format(paste0("[", 1:num_lines, "]"), 
                   width = inds_width, justify = "right")
  
  #lengths of sequences
  lens <- .get_lens(sq)
  
  #max lenght number width
  lens_width <- max(nchar(lens)) + 2
  
  #lengths to print
  p_lens <- paste0("<", lens, ">")
  if (use_color) p_lens <- blue(p_lens)
  
  needs_dots <- rep(FALSE, num_lines)
  
  for (i in 1:num_lines) {
    s <- sq_cut[[i]]
    # we count how much characters can we print by counting cumulative extent
    cum_lens <- cumsum(col_nchar(s)) + (0:(length(s) - 1)) * nchar(letters_sep)
    
    #max lenght of this line is p_width minus lenghts of lens and inds
    res_lens <- p_width - col_nchar(p_lens[i]) - nchar(p_inds[i]) - 2
    
    #we remove characters we cannot print
    s <- s[cum_lens < res_lens]
    n <- length(s)
    
    #if printed sequence is shorter than original, we need also space for dots
    if (n < lens[i]) {
      s <- s[cum_lens[1:n] < res_lens - 3]
      needs_dots[i] <- TRUE
    }
    sq_cut[[i]] <- s
  }
  
  #paste sequene
  p_body <- sapply(sq_cut, function(s) paste(s, collapse = letters_sep))
  if (use_color) p_body <- green(p_body)
  
  #dots
  p_dots <- ifelse(needs_dots, "...", "")
  if (use_color) p_dots <- silver(p_dots)
  
  #spaces between sequence and lens
  p_spaces <- sapply(p_width - col_nchar(p_body) - nchar(p_inds) - 
                       col_nchar(p_lens) - col_nchar(p_dots) - 2, 
                     function(l) paste0(rep(" ", l), collapse = ""))
  
  #paste and cat everything
  p_body <- paste0(p_inds, " ", p_body, p_dots, p_spaces, " ", p_lens, collapse = "\n")
  header <- .get_print_header(sq)
  footer <- .get_print_footer(x, num_lines)
  cat(header, p_body, footer, sep = "\n")
}

#' @exportMethod print encsq
#' @export
print.encsq <- function(x, ...) {
  sqclass <- "enc (encoded values) sequences vector:\n"
  na_char <- .get_na_char()

  alph <- .get_alph(x)
  decoded <- .apply_sq(x, "int", "none", function(s) alph[s])
  decoded <- sapply(decoded, function(s) ifelse(length(s) == 0, 
                                                "<NULL sq>", 
                                                paste(ifelse(is.na(s), na_char, s), collapse = " ")))
  max_width <- max(nchar(1:length(x)))
  inds <- paste0("[", 1:length(x), "] ")
  cat(sqclass, paste0(format(inds, width = max_width + 3, justify = "right"), 
                      decoded, 
                      collapse = "\n"), 
      "\n", sep = "")
}

#' @importFrom crayon col_nchar
#' @importFrom crayon green
#' @importFrom crayon silver
#' @importFrom pillar pillar_shaft
#' @importFrom pillar new_pillar_shaft
#' @importFrom pillar get_max_extent
#' @importFrom pillar get_extent
#' @exportMethod pillar_shaft sq
#' @export
pillar_shaft.sq <- function(x, ...) {
  x <- .debitify_sq(x, "string")
  if (.get_color_opt()) {
    na_char <- .get_na_char()
    x <- gsub(na_char, silver(na_char), x)
  }
  if (.get_color_opt()) x <- green(x)
    
  longest_str <- get_max_extent(x)
  min_str_width <- if (longest_str >= 6) 6 else longest_str
  
  opt <- .get_print_length()
  
  new_pillar_shaft(x,
                   width = min(longest_str + col_nchar(longest_str) + 3, 
                               opt + col_nchar(longest_str) + 6),
                   min_width = col_nchar(longest_str) + min_str_width + 3,
                   class = "pillar_shaft_sq",
                   align = "left")
}

#' @importFrom crayon col_nchar
#' @importFrom crayon blue
#' @importFrom crayon silver
#' @importFrom crayon col_substring
#' @importFrom pillar new_ornament
#' @exportMethod format pillar_shaft_sq
#' @export
format.pillar_shaft_sq <- function(x, width, ...) {
  if (width < attr(x, "min_width")) {
    stop("need at least width ", attr(x, "min_width"), ", requested ", 
         width, ".", call. = FALSE)
  } 
  
  s_opt <- if (.get_color_opt()) silver else identity
  b_opt <- if (.get_color_opt()) blue else identity
  
  seq_lengths <- col_nchar(x)
  nums_lengths <- nchar(seq_lengths)
  available_seq_space <- width - nums_lengths - 3
  inds_need_dots <- available_seq_space < seq_lengths
  available_seq_space[inds_need_dots] <- available_seq_space[inds_need_dots] - 3
  row <- col_substring(x, 1, available_seq_space)
  row[inds_need_dots] <- paste0(row[inds_need_dots], s_opt("..."))
  whitespace_lenghts <- available_seq_space - seq_lengths
  whitespace_lenghts[whitespace_lenghts < 0] <- 0
  whitespace_vec <- sapply(whitespace_lenghts, 
                           function(n) paste(rep(" ", n), collapse = ""))
  row <- paste0(row, whitespace_vec, b_opt(" <"), b_opt(seq_lengths), b_opt(">"))
  
  new_ornament(row, width = width, align = "left")
}

#' @importFrom pillar pillar_shaft
#' @importFrom pillar new_pillar_shaft
#' @importFrom pillar get_max_extent
#' @exportMethod pillar_shaft encsq
#' @export
pillar_shaft.encsq <- function(x, ...) {
  alph <- .get_alph(x)
  x <- .apply_sq(x, "int", "none", function(s) alph[s])
  
  x_min <- sapply(x, function(x) paste(format(x, digits = 1, nsmall = 1, scientific = FALSE), collapse = ""))
  
  longest_str <- get_max_extent(x_min)
  
  opt <- .get_print_length()
  
  new_pillar_shaft(x,
                   width = min(longest_str, 
                               opt + 3),
                   min_width = 7,
                   class = "pillar_shaft_encsq",
                   align = "left")
}

#' @importFrom crayon cyan
#' @importFrom crayon silver
#' @importFrom pillar new_ornament
#' @exportMethod format pillar_shaft_encsq
#' @export
format.pillar_shaft_encsq <- function(x, width, ...) {
  if (width < attr(x, "min_width")) {
    stop("need at least width ", attr(x, "min_width"), ", requested ", 
         width, ".", call. = FALSE)
  } 
  
  s_opt <- if (.get_color_opt()) cyan else identity
  
  x_t <- sapply(x, function(xth) paste(format(xth, digits = 1, nsmall = 1, scientific = FALSE), collapse = " "))
  if (max(nchar(x_t)) > width) {
    x_t <- paste0(substr(x_t, 1, width - 3), silver("..."))
  }
  
  row <- s_opt(x_t)
  
  new_ornament(row, width = width, align = "left")
}

.cut_sq <- function(sq, num_oct) {
  alph_size <- tidysq:::.get_alph_size(tidysq:::.get_alph(sq))
  ret <- lapply(sq, function(s) if (length(s) < num_oct * alph_size) s 
         else s[1:(num_oct * alph_size)])
  .set_class_alph(ret, sq)
}

.get_print_header <- function(sq) {
  type <- .get_sq_type(sq)
  if (length(type) != 1) {
    "sq (improper subtype!):"
  } else {
    type_msg <- switch(type,
                       ami = "ami (amino acids)",
                       nuc = "nuc (nucleotides)",
                       unt = "unt (unspecified type)",
                       atp = "atp (atypical alphabet)")
    clean_msg <- if (.is_cleaned(sq)) ", cln (cleaned)" else ""
    paste0(type_msg, clean_msg, " sequences list:")
  }
}

.get_print_footer <- function(sq, num_lines) {
  if (length(sq) > num_lines) 
    paste0("printed ", num_lines, " out of ", length(sq), "")
  else ""
}
