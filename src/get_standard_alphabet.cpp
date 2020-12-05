#include <Rcpp.h>

#include "tidysq.h"

using namespace tidysq;

// [[Rcpp::export]]
Rcpp::StringVector CPP_get_standard_alphabet(const std::string& dest_type) {
    return export_to_R(Alphabet(util::sq_type_for_sq_type_abbr(dest_type)));
}
