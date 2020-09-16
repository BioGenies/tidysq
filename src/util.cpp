#include <Rcpp.h>
#include <cmath>
#include <iostream>
#include <tidysq_generic.h>

// [[Rcpp::interfaces(cpp, r)]]
// [[Rcpp::export]]
unsigned short C_get_alph_size(Rcpp::CharacterVector alph) {
  return tidysq::get_alph_size_internal(alph);
}

// [[Rcpp::interfaces(cpp, r)]]
// [[Rcpp::export]]
std::vector<unsigned char> unpack_raws_to_std_vector(Rcpp::RawVector v, unsigned short alphabet_size) {
    return tidysq::unpack_raws_internal<Rcpp::RawVector, std::vector<unsigned char>>(v, alphabet_size);
}

// [[Rcpp::export]]
unsigned int C_get_out_len(Rcpp::RawVector packed,
                           const unsigned short alph_size) {
    return tidysq::get_out_len_internal<Rcpp::RawVector>(packed, alph_size);
}

// [[Rcpp::export]]
unsigned short C_get_na_val(const unsigned short alph_size) {
  return (1 << alph_size) - 1;
}

//' Match index of letter in alphabet for each letter
Rcpp::RawVector C_match(Rcpp::CharacterVector letters,
                        Rcpp::CharacterVector alph) {
  Rcpp::RawVector ret(letters.size());
  for (int i = 0; i < letters.size(); i++) {
    ret[i] = (1 << C_get_alph_size(alph)) - 1;
    for (int j = 0; j < alph.size(); j ++) {
      if (alph[j] == letters[i]) {
        ret[i] = j + 1;
        break;
      }
    }
  }

  return ret;
}

//' Match index of letter in alphabet for each letter
Rcpp::RawVector C_match(Rcpp::RawVector letters,
                        Rcpp::CharacterVector alph) {
  Rcpp::RawVector flat_alph(alph.size());
  for (int i = 0; i < alph.size(); i++) {
    flat_alph[i] = alph[i][0];
  }
  
  Rcpp::RawVector ret(letters.size());
  for (int i = 0; i < letters.size(); i++) {
    ret[i] = (1 << C_get_alph_size(alph)) - 1;
    for (int j = 0; j < alph.size(); j ++) {
      if (flat_alph[j] == letters[i]) {
        ret[i] = j + 1;
        break;
      }
    }
  }

  return ret;
}

//' Match index of letter in alphabet for each letter
Rcpp::CharacterVector C_match(Rcpp::RawVector letters,
                              Rcpp::CharacterVector alph,
                              Rcpp::CharacterVector na_letter) {
  Rcpp::CharacterVector ret(letters.size());
  unsigned char na_val = C_get_na_val(C_get_alph_size(alph));
  for (int i = 0; i < letters.size(); i++) {
    if (letters[i] == na_val) {
      ret[i] = na_letter[0];
    } else {
      ret[i] = alph[letters[i] - 1];
    }
  }

  return ret;
}