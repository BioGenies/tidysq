#include <Rcpp.h>

#include "tidysq/tidysq-includes.h"
#include "tidysq/substitute_letters.h"

using namespace tidysq;

// [[Rcpp::export]]
Rcpp::List CPP_substitute_letters(const Rcpp::List &x,
                                  const Rcpp::StringVector &encoding,
                                  const Rcpp::StringVector &NA_letter) {
    std::unordered_map<Letter, Letter> encoding_map;
    std::vector<std::string> std_encoding = util::convert_string_vector(encoding);
    std::vector<std::string> names = util::convert_string_vector(Rcpp::StringVector(encoding.names()));
    for (int i = 0; i < encoding.length(); ++i) {
        encoding_map.insert_or_assign(names[i], std_encoding[i]);
    }
    return export_to_R(substitute_letters<RCPP_IT>(import_from_R(x, NA_letter), encoding_map));
}