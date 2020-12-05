#include "tidysq.h"

using namespace tidysq;

//[[Rcpp::export]]
Rcpp::StringVector CPP_guess_standard_alph(const std::vector<std::string> &alph,
                                           const tidysq::Letter &NA_letter) {
    return export_to_R(
            Alphabet(alph, NA_letter));
}