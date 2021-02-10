#include <Rcpp.h>

#include "tidysq.h"

using namespace tidysq;

//[[Rcpp::export]]
Rcpp::List CPP_unpack_RAWS(const Rcpp::List& sq, const tidysq::Letter& NA_letter) {
    return export_to_R(import_sq_from_R(sq, NA_letter).unpack<RCPP_IT, RAWS_PT>());
}

//[[Rcpp::export]]
Rcpp::List CPP_unpack_INTS(const Rcpp::List& sq, const tidysq::Letter& NA_letter) {
    return export_to_R(import_sq_from_R(sq, NA_letter).unpack<RCPP_IT, INTS_PT>());
}

//[[Rcpp::export]]
Rcpp::List CPP_unpack_STRINGS(const Rcpp::List& sq, const tidysq::Letter& NA_letter) {
    return export_to_R(import_sq_from_R(sq, NA_letter).unpack<RCPP_IT, STRINGS_PT>());
}

//[[Rcpp::export]]
Rcpp::StringVector CPP_unpack_STRING(const Rcpp::List& sq, const tidysq::Letter& NA_letter) {
    return export_to_R(import_sq_from_R(sq, NA_letter).unpack<RCPP_IT, STRING_PT>());
}