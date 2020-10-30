#include <Rcpp.h>

#include "tidysq/exports.h"

using namespace tidysq;

//[[Rcpp::export]]
Rcpp::List CPP_unpack_RAWS(Rcpp::List sq, Rcpp::StringVector NA_letter) {
  return importFromR(sq, NA_letter)
  .unpack<RCPP, RAWS>()
  .exportToR();
}

//[[Rcpp::export]]
Rcpp::List CPP_unpack_INTS(Rcpp::List sq, Rcpp::StringVector NA_letter) {
  return importFromR(sq, NA_letter)
  .unpack<RCPP, INTS>()
  .exportToR();
}

//[[Rcpp::export]]
Rcpp::List CPP_unpack_STRINGS(Rcpp::List sq, Rcpp::StringVector NA_letter) {
  return importFromR(sq, NA_letter)
  .unpack<RCPP, STRINGS>()
  .exportToR();
}

//[[Rcpp::export]]
Rcpp::StringVector CPP_unpack_STRING(Rcpp::List sq, Rcpp::StringVector NA_letter) {
  return importFromR(sq, NA_letter)
  .unpack<RCPP, STRING>()
  .exportToR();
}