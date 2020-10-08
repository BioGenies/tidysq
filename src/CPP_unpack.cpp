#include <Rcpp.h>

#include "tidysq/exports.h"

using namespace tidysq;

//[[Rcpp::export]]
Rcpp::List CPP_unpack_RAWS(Rcpp::List sq) {
  return importFromR(sq)
  .unpack<RCPP, RAWS>()
  .exportToR();
}

//[[Rcpp::export]]
Rcpp::List CPP_unpack_INTS(Rcpp::List sq) {
  return importFromR(sq)
  .unpack<RCPP, INTS>()
  .exportToR();
}

//[[Rcpp::export]]
Rcpp::List CPP_unpack_STRINGS(Rcpp::List sq) {
  return importFromR(sq)
  .unpack<RCPP, STRINGS>()
  .exportToR();
}

////[[Rcpp::export]]
//Rcpp::StringVector CPP_unpack_STRING(Rcpp::List sq) {
//  return Sq<RCPP>(sq)
//  .unpack<RCPP, STRING>()
//  .exportToR();
//}