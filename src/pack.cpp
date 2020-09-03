#include <Rcpp.h>
#include <tidysq_generic.h>

// [[Rcpp::interfaces(cpp, r)]]

unsigned short C_get_alph_size(Rcpp::CharacterVector alph);
Rcpp::RawVector C_match(Rcpp::CharacterVector letters,
                        Rcpp::CharacterVector alph);
Rcpp::RawVector C_match(Rcpp::RawVector letters,
                        Rcpp::CharacterVector alph);

// [[Rcpp::export]]
Rcpp::RawVector C_pack_raws(Rcpp::RawVector unpacked,
                          const unsigned short alph_size) {
  return tidysq::pack_raws_internal(unpacked, alph_size);
}
  
// [[Rcpp::export]]
Rcpp::RawVector  C_pack_ints(Rcpp::IntegerVector unpacked,
                          const unsigned short alph_size) {
  Rcpp::RawVector ret(unpacked);
  return C_pack_raws(ret, alph_size);
}

//[[Rcpp::export]]
Rcpp::RawVector C_pack_chars(Rcpp::CharacterVector unpacked,
                           Rcpp::CharacterVector alph) {
  unsigned short alph_size = C_get_alph_size(alph);
  Rcpp::RawVector ret = C_match(unpacked, alph);
  return C_pack_raws(ret, alph_size);
}

//[[Rcpp::export]]
Rcpp::RawVector C_pack_string(Rcpp::RawVector unpacked,
                            Rcpp::CharacterVector alph) {
  unsigned short alph_size = C_get_alph_size(alph);
  Rcpp::RawVector ret = C_match(unpacked, alph);
  return C_pack_raws(ret, alph_size);
}