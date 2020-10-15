#include <Rcpp.h>
#include <algorithm>
#include <set>
#include <iterator>

#include "tidysq/exports.h"
#include "tidysq/ops/internal/util.h"

using namespace tidysq;

// [[Rcpp::export]]
Rcpp::List CPP_find_invalid_letters(Rcpp::List x,
                                    Rcpp::StringVector dest_type) {
  const Sq<RCPP> sq = importFromR(x, "!");
  Alphabet alph = sq.alphabet();
  std::vector<std::string> dest_alph =
          util::standard_letters_for_type(util::sq_type_for_abbr(dest_type));
  
  std::vector<LetterValue> invalid_indices;
  for (LetterValue i = 0; i < alph.length(); ++i) {
    if (std::none_of(dest_alph.begin(), dest_alph.end(),
                     [alph, i](std::string letter){ return alph[i] == letter; })) {
      invalid_indices.push_back(i);
    }
  }
  
  Rcpp::List ret = Rcpp::List::create();
  for (LetterValue i = 0; i < sq.length(); ++i) {
    ret.push_back(Rcpp::StringVector());
  }
  
  if (invalid_indices.empty()) {
    return ret;
  }
  
  return Rcpp::List::create(invalid_indices);
}
