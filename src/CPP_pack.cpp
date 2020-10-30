#include <Rcpp.h>

#include "tidysq/exports.h"

using namespace tidysq;

//[[Rcpp::export]]
Rcpp::List CPP_pack_RAWS(Rcpp::List proto, Rcpp::StringVector alphabet, Rcpp::StringVector NA_letter) {
  return importProtoFromR<RAWS>(proto, alphabet, NA_letter)
  .pack<RCPP>()
  .exportToR();
}

//[[Rcpp::export]]
Rcpp::List CPP_pack_INTS(Rcpp::List proto, Rcpp::StringVector alphabet, Rcpp::StringVector NA_letter) {
  return importProtoFromR<INTS>(proto, alphabet, NA_letter)
  .pack<RCPP>()
  .exportToR();
}

//[[Rcpp::export]]
Rcpp::List CPP_pack_STRINGS(Rcpp::List proto, Rcpp::StringVector alphabet, Rcpp::StringVector NA_letter) {
  return importProtoFromR<STRINGS>(proto, alphabet, NA_letter)
  .pack<RCPP>()
  .exportToR();
}

//[[Rcpp::export]]
Rcpp::List CPP_pack_STRING(Rcpp::StringVector proto, Rcpp::StringVector alphabet, Rcpp::StringVector NA_letter) {
  return importProtoFromR<STRING>(proto, alphabet, NA_letter)
  .pack<RCPP>()
  .exportToR();
}