alphabet <- function(sq)
  attr(sq, "alphabet")

`alphabet<-` <- function(sq, value) {
  attr(sq, "alphabet") <- value
  sq
}

`[.sq_alphabet` <- function(x, i)
  if (i == .get_na_val(x)) na_character(x) else vec_data(x)[i]

na_character <- function(alph)
  attr(alph, "na_character")

`na_character<-` <- function(alph, value) {
  attr(alph, "na_character") <- value
  alph
}
