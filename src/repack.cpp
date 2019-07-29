#include <Rcpp.h>

Rcpp::RawVector unpack_unt(Rcpp::RawVector PACKED, const char na_char,
                           Rcpp::RawVector alph, const unsigned short ALPH_SIZE);

Rcpp::RawVector pack_unt(Rcpp::RawVector UNPACKED, Rcpp::RawVector alph, 
                         const unsigned short ALPH_SIZE);


inline char match(char letter, Rcpp::RawVector alph, 
                  char na_char, Rcpp::RawVector letters_fun) {
  for (int i = 0; i < alph.size(); i++) {
    if (letter == alph[i])
      return letters_fun[i];
  }
  return na_char;
}

// [[Rcpp::export]]
Rcpp::RawVector repack(Rcpp::RawVector PACKED, 
                       char na_char,
                       Rcpp::RawVector alph, 
                       Rcpp::RawVector new_alph,
                       unsigned short int alph_size,
                       unsigned short int new_alph_size,
                       Rcpp::RawVector letters_fun) {
  Rcpp::RawVector unpacked = unpack_unt(PACKED, na_char, alph, alph_size); 
  for (int i = 0; i < unpacked.size(); i++) {
    std::cout << "before match: " << unpacked[i] << "  ";
    unpacked[i] = match(unpacked[i], alph, na_char, letters_fun);
    std::cout << "after: " << unpacked[i];
  }
  return pack_unt(unpacked, new_alph, new_alph_size);
}