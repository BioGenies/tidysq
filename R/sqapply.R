#' @export 
sqapply <- function(sq, fun, ..., paste_char = FALSE, 
                    use_na_char = if (paste_char) TRUE else FALSE) {
  validate_sq(sq)
  .check_logical(paste_char, "'paste_char'", single_elem = TRUE)
  .check_logical(use_na_char, "'use_na_char'", single_elem = TRUE)
  .check_paste_or_na(paste_char, use_na_char)
  
  na_char <- .get_na_char()
  type <- .get_sq_type(sq)
  .apply_sq(sq, if (paste_char) "string" else "char", "none", function(s) {
    if (!use_na_char) s[s == na_char] <- NA
    if (type == "enc") s <- as.numeric(s)
    fun(s, ...)
  })
}