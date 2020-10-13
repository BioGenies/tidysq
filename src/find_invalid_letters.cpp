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
    util::getStandardLetters(util::get_sq_type_from_abbr(dest_type));
  
  std::vector<LetValue> invalid_indices;
  for (LetValue i = 0; i < alph.length(); ++i) {
    if (std::none_of(dest_alph.begin(), dest_alph.end(),
                     [alph, i](std::string letter){ return alph[i] == letter; })) {
      invalid_indices.push_back(i);
    }
  }
  return Rcpp::List::create(invalid_indices);
}
