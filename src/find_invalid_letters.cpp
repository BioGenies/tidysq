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
  const Alphabet alph = sq.alphabet();
  const std::vector<std::string> dest_alph =
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
  
  // Process only whenever we know there are invalid letters
  // Might remove that, as iterating over invalid_indices does basically the same thing
  if (!invalid_indices.empty()) {
    for (LenSq i = 0; i < sq.length(); ++i) {
      const Sequence<RCPP> sequence = sq[i];
      std::vector<LetterValue> invalid_found = {};
      for (LetterValue index : invalid_indices) {
        // We need input iterator to make it work
        if (std::any_of(sequence.begin(), sequence.end(), [=](ElementPacked elem){ return elem == index; })) {
          ret[i].push_back(alph[index]);
        }
      }
    }
  }
  
  return ret;
}
