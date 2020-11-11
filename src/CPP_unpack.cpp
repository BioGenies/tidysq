#include <Rcpp.h>

#include "tidysq/tidysq-includes.h"

using namespace tidysq;

//[[Rcpp::export]]
Rcpp::List CPP_unpack_RAWS(const Rcpp::List& sq, const Rcpp::StringVector& NA_letter) {
    return export_to_R(
            import_from_R(sq, NA_letter)
                    .unpack<RCPP, RAWS>());
}

//[[Rcpp::export]]
Rcpp::List CPP_unpack_INTS(const Rcpp::List& sq, const Rcpp::StringVector& NA_letter) {
    return export_to_R(
            import_from_R(sq, NA_letter)
                    .unpack<RCPP, INTS>());
}

//[[Rcpp::export]]
Rcpp::List CPP_unpack_STRINGS(const Rcpp::List& sq, const Rcpp::StringVector& NA_letter) {
    return export_to_R(
            import_from_R(sq, NA_letter)
                    .unpack<RCPP, STRINGS>());
}

//[[Rcpp::export]]
Rcpp::StringVector CPP_unpack_STRING(const Rcpp::List& sq, const Rcpp::StringVector& NA_letter) {
    return export_to_R(
            import_from_R(sq, NA_letter)
                    .unpack<RCPP, STRING>());
}