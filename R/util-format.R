# TODO: issue #59

#' @importFrom cli col_green
#' @export
format.sq <- function(x, ...,
                      max_sequences = getOption("tidysq_print_max_sequences"),
                      use_color = getOption("tidysq_print_use_color"),
                      NA_letter = getOption("tidysq_NA_letter"),
                      letters_sep = NULL) {
  assert_count(max_sequences)
  assert_flag(use_color)
  assert_string(NA_letter)
  assert_string(letters_sep, null.ok = TRUE)
  
  internal_format_sq(x, max_sequences, use_color, NA_letter, letters_sep, col_green)
}

#' @importFrom cli col_cyan
#' @export
format.encsq <- function(x, ...,
                         max_sequences = getOption("tidysq_print_max_sequences"),
                         use_color = getOption("tidysq_print_use_color"),
                         NA_letter = getOption("tidysq_NA_letter"),
                         letters_sep = NULL,
                         digits = 2) {
  assert_count(max_sequences)
  assert_flag(use_color)
  assert_string(letters_sep, null.ok = TRUE)
  assert_count(digits)
  assert_string(NA_letter)
  
  alphabet(x) <- format(alphabet(x), digits = digits, scientific = FALSE)
  
  internal_format_sq(x, max_sequences, use_color, NA_letter, letters_sep, col_cyan)
}

#' @importFrom pillar new_ornament
#' @export
format.pillar_shaft_sq <- function(x, width, ...) {
  if (width < attr(x, "min_width"))
    stop("need at least width ", attr(x, "min_width"), ", requested ", width, ".", call. = FALSE)

  lens <- attr(x, "lens")
  letters_sep <- attr(x, "letters_sep")
  body_color <- attr(x, "body_color")
  align <- attr(x, "align")
  
  p_seqs <- format_sequences_and_lengths(x, lens, letters_sep, width, TRUE)
  
  # cat everything
  new_ornament(p_seqs, width = width, align = align)
}

internal_format_sq <- function(x, max_sequences, use_color, NA_letter, letters_sep, body_color = identity) {
  # select at most max_sequences to print
  x <- reduce_num_sequences(x, max_sequences)

  max_inds_width <- nchar(length(x)) + 2 # max width of index number is the longest number + space for "[]"
  
  console_width <- getOption("width")
  NA_letter <- if (use_color) col_silver(NA_letter) else NA_letter
  letters_sep <- choose_letters_sep(letters_sep, alphabet(x))
  lens <- get_sq_lengths(x) # lengths of sequences
  
  # cut sq object so that we don't need to unpack long sequences
  x <- precut_body_as_unpackeds(x, console_width, letters_sep, NA_letter, body_color)
  
  paste0(format_indices(length(x), max_inds_width),
         "\u00a0",  #unbreakable space
         format_sequences_and_lengths(x, lens, letters_sep,
                                      console_width - (max_inds_width + 1),  #subtracting chars for inds and space
                                      use_color),
         collapse = "\n")
}

#' @importFrom crayon col_nchar
#' @importFrom pillar new_pillar_shaft
pillar_shaft_sq <- function(x, letters_sep, NA_letter, body_color = identity, max_pillar_width) {
  # maximum length of length numbers
  lens <- get_sq_lengths(x)
  max_len_width <- max(nchar(lens))

  NA_letter <- col_silver(NA_letter)
  
  # cut sq object so that we don't need to unpack long sequences
  x <- precut_body_as_unpackeds(x, max_pillar_width, letters_sep, NA_letter, body_color)

  
  max_str_width <- max(sapply(x, function(s) sum(col_nchar(s))))
  min_str_width <- if (max_str_width >= 6) 6 else max_str_width
  
  new_pillar_shaft(x,
                   width = max_len_width + 3 + min(max_str_width, max_pillar_width + 3),
                   min_width = max_len_width + 3 + min_str_width,
                   class = "pillar_shaft_sq",
                   align = "left",
                   lens = lens,
                   letters_sep = letters_sep,
                   body_color = body_color)
}

choose_letters_sep <- function(letters_sep, alph) {
  # \u00a0 is hard space
  if (!is.null(letters_sep)) letters_sep
  else if (is.numeric(alph)) "\u00a0"
  else if (all(nchar(alph) == 1)) ""
  else "\u00a0"
}

format_indices <- function(num_lines, max_inds_width)
  format(paste0("[", 1:num_lines, "]"), width = max_inds_width, justify = "right")


cut_sq <- function(x, print_width, letters_sep_nchar) {
  # note: num_bytes should be greater than 0, not that it'd make sense to use 0 or less
  # (max) number of bytes that fit = (max) number of bits that fit / 8  rounded up
  # (max) number of bits that fit = (max) number of letters that fit * bits per letter (a.k.a alphabet size)
  # (max) number of letters that fit = available space (a.k.a print_width - 6) / (min) chars per letter (a.k.a 1 (min letter width) + letters_sep_nchar))
  alph_size <- size(alphabet(x))
  num_bytes <- ceiling((print_width - 6) * alph_size / (8 * (letters_sep_nchar + 1)))
  ret <- lapply(x, function(s) {
    # if sequence is shorter than number of bytes that fit - leave it
    # otherwise cut it and calculate how many letters fit
    if (length(s) <= num_bytes) s
    else structure(s[1:(num_bytes)],
                   original_length = min(attr(s, "original_length"), floor(num_bytes * 8 / alph_size)))
  })
  vec_restore(ret, x)
}

precut_body_as_unpackeds <- function(x, print_width, letters_sep, NA_letter, body_color) {
  sq_cut <- cut_sq(x, print_width, nchar(letters_sep))
  alphabet(sq_cut) <- vec_restore(body_color(alphabet(x)), alphabet(x))
  unpack(sq_cut, "STRINGS", NA_letter)
}

#' @importFrom cli col_blue col_silver
#' @importFrom crayon col_nchar
#' @importFrom utils tail
format_sequences_and_lengths <- function(x, lens, letters_sep, width, use_color) {
  lens_width <- max(nchar(lens)) + 2 # max width of length number

  col_dots <- if (use_color) col_silver else identity
  col_lens <- if (use_color) col_blue else identity

  x <- mapply(function(sequence, len) {
    if (len == 0) {
      structure(col_dots("<NULL>"), dots = "")
    } else {
      # we count how much characters can we print by counting cumulative extent
      cumulative_lens <- cumsum(col_nchar(sequence)) + (0:(length(sequence) - 1)) * nchar(letters_sep)
      # max length of this line is its width minus the lens_width
      res_lens <- width - lens_width - 1

      # if total length is greater than reserved space, we have to cut it down
      if (tail(cumulative_lens, n = 1) > res_lens) {
        # if printed sequence is shorter than original, we also need space for dots
        # find last index that allows sq to fit into reserved space
        n <- Position(identity, cumulative_lens <= res_lens - 3, right = TRUE)
        sequence <- sequence[1:n]
        attr(sequence, "dots") <- "..."
      } else {
        attr(sequence, "dots") <- ""
      }
      sequence
    }
  }, x, lens, SIMPLIFY = FALSE)

  # strings formatted to be printed
  elements <- sapply(x, paste, collapse = letters_sep)
  dots <- col_dots(sapply(x, attr, "dots"))
  lengths <- col_lens(paste0("<", lens, ">"))
  spaces <- sapply(width - col_nchar(elements) - col_nchar(lengths) - col_nchar(dots),
                   function(l) paste0(rep("\u00a0", l), collapse = ""))

  # return pasted everything as a vector of strings, each for one (printed) sequence
  paste0(elements, dots, spaces, lengths)
}

reduce_num_sequences <- function(x, max_sequences) {
  x <- x[1:min(max_sequences, length(x))]
}