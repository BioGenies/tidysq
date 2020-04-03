#include <Rcpp.h>
#include <cmath>
#include <iostream>

// [[Rcpp::interfaces(cpp, r)]]

// [[Rcpp::export]]
unsigned short get_alph_size(Rcpp::CharacterVector alph) {
  return ceil(log2(alph.size() + 2));
}

unsigned short get_na_val(const unsigned short alph_size) {
  return (1 << alph_size) - 1;
}

// [[Rcpp::export]]
Rcpp::RawVector match_chars(Rcpp::CharacterVector letters,
                            Rcpp::CharacterVector alph) {
  Rcpp::RawVector ret(letters.size());
  for (int i = 0; i < letters.size(); i++) {
    ret[i] = (1 << get_alph_size(alph)) - 1;
    for (int j = 0; j < alph.size(); j ++) {
      if (alph[j] == letters[i]) {
        ret[i] = j + 1;
        break;
      }
    }
  }

  return ret;
}

// [[Rcpp::export]]
Rcpp::RawVector match_char(Rcpp::RawVector letters,
                           Rcpp::CharacterVector alph) {
  Rcpp::RawVector flat_alph(alph.size());
  for (int i = 0; i < alph.size(); i++) {
    flat_alph[i] = alph[i][0];
  }
  
  Rcpp::RawVector ret(letters.size());
  for (int i = 0; i < letters.size(); i++) {
    ret[i] = (1 << get_alph_size(alph)) - 1;
    for (int j = 0; j < alph.size(); j ++) {
      if (flat_alph[j] == letters[i]) {
        ret[i] = j + 1;
        break;
      }
    }
  }

  return ret;
}

// [[Rcpp::export]]
Rcpp::CharacterVector match_raws(Rcpp::RawVector letters,
                                 Rcpp::CharacterVector alph,
                                 Rcpp::CharacterVector na_char) {
  Rcpp::CharacterVector ret(letters.size());
  unsigned char na_val = get_na_val(get_alph_size(alph));
  for (int i = 0; i < letters.size(); i++) {
    if (letters[i] == na_val) {
      ret[i] = na_char[0];
    } else {
      ret[i] = alph[letters[i] - 1];
    }
  }

  return ret;
}