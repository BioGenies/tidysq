#include <Rcpp.h>
#include <string>
#include <iostream>

// [[Rcpp::interfaces(cpp, r)]]

// [[Rcpp::export]]
Rcpp::CharacterVector C_get_real_alph(Rcpp::CharacterVector sq) {
  std::vector<std::string> sq_conv = Rcpp::as<std::vector<std::string>>(sq);
  unsigned int input_size = sq_conv.size();
  
  std::vector<char> chars;
  
  // reads passed string in order (is parallelization possible?)
  for (unsigned int i = 0; i < input_size; i++) {
    for (char& c : sq_conv[i]) {
      std::vector<char>::iterator it = std::find(chars.begin(), chars.end(), c);
      if (it == chars.end()) {
        chars.push_back(c);
      }
    }
  }
  return Rcpp::wrap(chars);
}
