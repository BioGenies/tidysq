#include <Rcpp.h>

#include "tidysq.h"

using namespace tidysq;

// [[Rcpp::export]]
Rcpp::StringVector CPP_get_standard_alphabet(const Rcpp::StringVector& dest_type) {
    return export_to_R(Alphabet(util::sq_type_for_sq_type_abbr(util::convert_to_scalar(dest_type))));
}
