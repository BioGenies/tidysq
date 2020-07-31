.print_sq <- function(x, alph, max_sequences, use_color, letters_sep, body_color) {
  # select at most max_sequences to print
  num_lines <- min(max_sequences, length(x))
  if (num_lines == 0) {
    cat(.get_print_empty_sq(x))
    return()
  }
  sq <- x[1:num_lines]
  
  p_width <- getOption("width")
  
  # cut sq object so that we don't need to debitify long sequences
  # 6 is minimum length of p_lens and p_inds, 8 is byte length
  sq_cut <- .cut_sq(sq, ceiling((p_width - 6) / (8 * (nchar(letters_sep) + 1))))
  sq_cut <- .unpack_from_sq(sq_cut, "int")
  sq_cut <- lapply(sq_cut, function(s) alph[s])
  
  # lengths of sequences
  lens <- lengths(sq)
  
  # max width of index number
  inds_width <- nchar(num_lines) + 2
  
  # indices to print
  p_inds <- format(paste0("[", 1:num_lines, "]"), 
                   width = inds_width, justify = "right")
  
  p_seqs <- .get_p_seqs(sq_cut, lens, letters_sep, body_color, p_width - inds_width - 1, use_color)
  
  cat(
    .get_print_nonempty_header(sq),  # print header
    paste0(p_inds, " ", p_seqs, collapse = "\n"),  # print body
    if (length(x) > num_lines) .get_print_footer(x, num_lines),  # print footer (if needed)
    sep = "\n"
  )
}

.cut_sq <- function(sq, num_oct) {
  # note: num_oct should be greater than 0, not that it'd make sense to use 0 or less
  alph_size <- .get_alph_size(.get_alph(sq))
  ret <- lapply(sq, function(s) {
    if (length(s) <= num_oct * alph_size) s 
    else structure(s[1:(num_oct * alph_size)], original_length = attr(s, "original_length"))
  })
  .set_class_alph(ret, sq)
}

.get_print_header <- function(sq) {
  type <- .get_sq_type(sq)
  if (length(type) != 1) {
    "sq (improper subtype!)"
  } else {
    paste0(vec_ptype_full(sq), " sequences list")
  }
}

.get_print_nonempty_header <- function(sq) {
  paste0(.get_print_header(sq), ":")
}

.get_print_empty_sq <- function(sq) {
  paste0(.get_print_header(sq), " of length 0")
}

.get_print_footer <- function(sq, num_lines) {
  paste0("printed ", num_lines, " out of ", length(sq), "")
}

#' @importFrom crayon col_nchar
#' @importFrom crayon blue
#' @importFrom crayon silver
#' @importFrom utils tail
.get_p_seqs <- function(x, lens, letters_sep, body_color, width, use_color = .get_color_opt()) {
  # max width of length number
  lens_width <- max(nchar(lens)) + 2
  
  x <- mapply(function(sq, len) {
    if (len == 0) {
      structure("<NULL>", dots = "")
    } else {
      # we count how much characters can we print by counting cumulative extent
      cum_lens <- cumsum(col_nchar(sq)) + (0:(length(sq) - 1)) * nchar(letters_sep)
      # max length of this line is its width minus the lens_width
      res_lens <- width - lens_width - 1
      
      # if total length is greater than reserved space, we have to cut it down
      if (tail(cum_lens, n = 1) > res_lens) {
        # if printed sequence is shorter than original, we also need space for dots
        # find last index that allows sq to fit into reserved space
        n <- Position(identity, cum_lens <= res_lens - 3, right = TRUE)
        sq <- sq[1:n]
        attr(sq, "dots") <- "..."
      } else {
        attr(sq, "dots") <- ""
      }
      sq
    }
  }, x, lens, SIMPLIFY = FALSE)
  
  # strings formatted (without colors) to be printed
  # sequences
  p_body <- sapply(x, paste, collapse = letters_sep)
  # dots
  p_dots <- sapply(x, attr, "dots")
  # lengths
  p_lens <- paste0("<", lens, ">")
  # spaces between sequences and lens
  p_spaces <- sapply(width - col_nchar(p_body) - col_nchar(p_lens) - col_nchar(p_dots),
                     function(l) paste0(rep(" ", l), collapse = ""))
  
  # add colors to text if necessary
  if (use_color) {
    p_lens <- blue(p_lens)
    p_body <- mapply(function(content, len)
      if (len == 0) silver(content) else body_color(content),
      p_body, lens)
    p_dots <- silver(p_dots)
  }
  
  # return pasted everything as a vector of strings, each for one (printed) sequence
  paste0(p_body, p_dots, p_spaces, p_lens)
}

#' @importFrom crayon col_nchar
#' @importFrom pillar new_pillar_shaft
.pillar_shaft_sq <- function(x, alph, letters_sep, body_color) {
  p_width <- getOption("width")
  
  # cut sq object so that we don't need to debitify long sequences
  # 6 is minimum length of p_lens and p_inds, 8 is byte length
  sq_cut <- .cut_sq(x, ceiling((p_width - 6) / (8 * (nchar(letters_sep) + 1))))
  sq_cut <- .unpack_from_sq(sq_cut, "int")
  sq_cut <- lapply(sq_cut, function(s) alph[s])
  
  # maximum length of length numbers
  lens <- lengths(x)
  max_len_width <- max(nchar(lens))
  
  max_str_width <- max(sapply(sq_cut, function(s) sum(col_nchar(s))))
  min_str_width <- if (max_str_width >= 6) 6 else max_str_width
  
  new_pillar_shaft(sq_cut,
                   width = max_len_width + 3 + min(max_str_width, .get_print_length() + 3),
                   min_width = max_len_width + 3 + min_str_width,
                   class = "pillar_shaft_sq",
                   align = "left",
                   lens = lens,
                   letters_sep = letters_sep,
                   body_color = body_color)
}
