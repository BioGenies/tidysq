#' @importFrom cli col_silver col_green
#' @importFrom crayon col_nchar
#' @export
format.sq <- function(x, ...,
                      max_sequences = getOption("tidysq_p_max_sequences"),
                      use_color = getOption("tidysq_p_use_color"),
                      letters_sep = NULL) {
  assert_count(max_sequences)
  assert_flag(use_color)
  assert_string(letters_sep, null.ok = TRUE)
  
  # color NA's
  alph <- alphabet(x)
  attr(alph, "na_letter") <- if (use_color)
    col_silver(getOption("tidysq_NA_letter")) else
      getOption("tidysq_NA_letter")
  
  # if parameter is NULL and all letters are length one, no space
  if (is.null(letters_sep))
    letters_sep <- if (all(col_nchar(alph[!is.na(alph)]) == 1)) "" else "\u00a0"
  
  .format_sq(x, max_sequences, use_color, letters_sep, col_green)
}

#' @importFrom cli col_cyan
#' @export
format.encsq <- function(x, ...,
                         max_sequences = getOption("tidysq_p_max_sequences"),
                         use_color = getOption("tidysq_p_use_color"),
                         letters_sep = NULL,
                         digits = 2) {
  assert_count(max_sequences)
  assert_flag(use_color)
  assert_string(letters_sep, null.ok = TRUE)
  assert_count(digits)
  
  alphabet(x) <- format(alphabet(x), digits = digits, scientific = FALSE)
  
  # if parameter is NULL default sep is space
  if (is.null(letters_sep))
    letters_sep <- "\u00a0"
  
  .format_sq(x, max_sequences, use_color, letters_sep, col_cyan)
}

#' @importFrom pillar new_ornament
#' @export
format.pillar_shaft_sq <- function(x, width, ...) {
  if (width < attr(x, "min_width"))
    stop("need at least width ", attr(x, "min_width"), ", requested ", width, ".", call. = FALSE)
  
  # TODO: maybe get rid of those intermediate variables?
  lens <- attr(x, "lens")
  letters_sep <- attr(x, "letters_sep")
  body_color <- attr(x, "body_color")
  align <- attr(x, "align")
  
  p_seqs <- .get_p_seqs(x, lens, letters_sep, body_color, width)
  
  # cat everything
  new_ornament(p_seqs, width = width, align = align)
}

.format_sq <- function(x, max_sequences, use_color, letters_sep, body_color) {
  # select at most max_sequences to print
  num_lines <- min(max_sequences, length(x))
  x <- x[1:num_lines]
  
  p_width <- getOption("width")
  
  # cut sq object so that we don't need to debitify long sequences
  # 6 is minimum length of p_lens and p_inds, 8 is byte length
  sq_cut <- .cut_sq(x, ceiling((p_width - 6) / (8 * (nchar(letters_sep) + 1))))
  sq_cut <- unpack(sq_cut, "INTS")
  sq_cut <- lapply(sq_cut, function(s) alphabet(x)[s])
  
  # lengths of sequences
  lens <- get_sq_lengths(x)
  
  # max width of index number
  inds_width <- nchar(num_lines) + 2
  
  # indices to print
  p_inds <- format(paste0("[", 1:num_lines, "]"), 
                   width = inds_width, justify = "right")
  
  p_seqs <- .get_p_seqs(sq_cut, lens, letters_sep, body_color, p_width - inds_width - 1, use_color)
  
  paste0(p_inds, "\u00a0", p_seqs, collapse = "\n")
}

.cut_sq <- function(x, num_oct) {
  # note: num_oct should be greater than 0, not that it'd make sense to use 0 or less
  alph_size <- .get_alph_size(alphabet(x))
  ret <- lapply(x, function(s) {
    if (length(s) <= num_oct * alph_size) s
    else structure(s[1:(num_oct * alph_size)], original_length = attr(s, "original_length"))
  })
  vec_restore(ret, x)
}

#' @importFrom cli col_blue col_silver
#' @importFrom crayon col_nchar
#' @importFrom utils tail
.get_p_seqs <- function(x, lens, letters_sep, body_color, width, use_color = .get_color_opt()) {
  # max width of length number
  lens_width <- max(nchar(lens)) + 2
  
  x <- mapply(function(sequence, len) {
    if (len == 0) {
      structure("<NULL>", dots = "")
    } else {
      # we count how much characters can we print by counting cumulative extent
      cum_lens <- cumsum(col_nchar(sequence)) + (0:(length(sequence) - 1)) * nchar(letters_sep)
      # max length of this line is its width minus the lens_width
      res_lens <- width - lens_width - 1
      
      # if total length is greater than reserved space, we have to cut it down
      if (tail(cum_lens, n = 1) > res_lens) {
        # if printed sequence is shorter than original, we also need space for dots
        # find last index that allows sq to fit into reserved space
        n <- Position(identity, cum_lens <= res_lens - 3, right = TRUE)
        sequence <- sequence[1:n]
        attr(sequence, "dots") <- "..."
      } else {
        attr(sequence, "dots") <- ""
      }
      sequence
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
                     function(l) paste0(rep("\u00a0", l), collapse = ""))
  
  # add colors to text if necessary
  if (use_color) {
    p_lens <- col_blue(p_lens)
    p_body <- mapply(function(content, len)
      if (len == 0) col_silver(content) else body_color(content),
      p_body, lens)
    p_dots <- col_silver(p_dots)
  }
  
  # return pasted everything as a vector of strings, each for one (printed) sequence
  paste0(p_body, p_dots, p_spaces, p_lens)
}

#' @importFrom crayon col_nchar
#' @importFrom pillar new_pillar_shaft
.pillar_shaft_sq <- function(x, letters_sep, body_color) {
  p_width <- getOption("width")
  
  # cut sq object so that we don't need to debitify long sequences
  # 6 is minimum length of p_lens and p_inds, 8 is byte length
  sq_cut <- .cut_sq(x, ceiling((p_width - 6) / (8 * (nchar(letters_sep) + 1))))
  sq_cut <- unpack(sq_cut, "INTS")
  sq_cut <- lapply(sq_cut, function(s) alphabet(x)[s])
  
  # maximum length of length numbers
  lens <- get_sq_lengths(x)
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
