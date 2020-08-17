# alphabet assignment ----

alphabet <- function(sq)
  attr(sq, "alphabet")

`alphabet<-` <- function(sq, value) {
  attr(sq, "alphabet") <- value
  sq
}

# alphabet creation ----

sq_alphabet <- function(alph, na_char = .get_na_char()) {
  new_vctr(
    alph,
    na_character = na_char,
    class = c("sq_alphabet", "character")
  )
}

sq_alphabet_ptype <- function()
  sq_alphabet(character())

.skip_characters <- function(alph, chars)
  vec_restore(setdiff(alph, chars), alph)

# alphabet reading ----

`[.sq_alphabet` <- function(x, i) {
  ret <- vec_data(x)[i]
  ret[i == .get_na_val(x)] <- na_character(x)
  ret
}

na_character <- function(alph)
  attr(alph, "na_character")

`na_character<-` <- function(alph, value) {
  attr(alph, "na_character") <- value
  alph
}

# various internal methods put together (to check!) ----

.get_alph_size <- function(alph) {
  ceiling(log2(length(alph) + 2))
}

.get_na_val <- function(alph) {
  2 ^ .get_alph_size(alph) - 1
}

.get_real_alph <- function(str_sq) {
  new_vctr(
    C_get_real_alph(str_sq),
    na_character = .get_na_char(),
    class = c("sq_alphabet", "character")
  )
}

.get_standard_alph <- function(type, is_clean) {
  new_vctr(
    if (type == "ami" &&  is_clean)
      aminoacids_df[!aminoacids_df[["amb"]], "one"]
    else if (type == "ami" && !is_clean)
      aminoacids_df[, "one"]
    else if (type == "dna" &&  is_clean)
      nucleotides_df[nucleotides_df[["dna"]], "one"]
    else if (type == "dna" && !is_clean)
      nucleotides_df[nucleotides_df[["dna"]] | nucleotides_df[["amb"]], "one"]
    else if (type == "rna" &&  is_clean)
      nucleotides_df[nucleotides_df[["rna"]], "one"]
    else if (type == "rna" && !is_clean)
      nucleotides_df[nucleotides_df[["rna"]] | nucleotides_df[["amb"]], "one"],
    na_character = .get_na_char(),
    class = c("sq_alphabet", "character")
  )
}

.guess_ami_is_clean <- function(real_alph) {
  if (all(real_alph %in% .get_standard_alph("ami", TRUE)))
    TRUE
  else if (all(real_alph %in% .get_standard_alph("ami", FALSE)))
    FALSE
  else stop("there are letters that aren't in IUPAC standard! (see: aminoacids_df)")
}

.guess_dna_is_clean <- function(real_alph) {
  if (all(real_alph %in% .get_standard_alph("dna", TRUE)))
    TRUE
  else if (all(real_alph %in% .get_standard_alph("dna", FALSE)))
    FALSE
  else stop("there are letters that aren't in IUPAC standard! (see: nucleotides_df)")
}

.guess_rna_is_clean <- function(real_alph) {
  if (all(real_alph %in% .get_standard_alph("rna", TRUE)))
    TRUE
  else if (all(real_alph %in% .get_standard_alph("rna", FALSE)))
    FALSE
  else stop("there are letters that aren't in IUPAC standard! (see: nucleotides_df)")
}
