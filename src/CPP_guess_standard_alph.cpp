#include "tidysq/exports.h"

using namespace tidysq;

//[[Rcpp::export]]
Rcpp::StringVector CPP_guess_standard_alph(const Rcpp::StringVector &alph,
                                           const Rcpp::StringVector &NA_letter) {
    return export_to_R(
            Alphabet(util::convert_string_vector(alph), util::get_scalar_string_value(NA_letter)));
}