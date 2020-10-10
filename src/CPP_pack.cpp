#include <Rcpp.h>

#include "tidysq/exports.h"

using namespace tidysq;

//[[Rcpp::export]]
Rcpp::List CPP_pack_RAWS(Rcpp::List proto, Rcpp::StringVector alphabet) {
  return importProtoFromR<RAWS>(proto, alphabet)
  .pack<RCPP>()
  .exportToR();
}

//[[Rcpp::export]]
Rcpp::List CPP_pack_INTS(Rcpp::List proto, Rcpp::StringVector alphabet) {
  return importProtoFromR<INTS>(proto, alphabet)
  .pack<RCPP>()
  .exportToR();
}

//[[Rcpp::export]]
Rcpp::List CPP_pack_STRINGS(Rcpp::List proto, Rcpp::StringVector alphabet) {
  return importProtoFromR<STRINGS>(proto, alphabet)
  .pack<RCPP>()
  .exportToR();
}

//[[Rcpp::export]]
Rcpp::List CPP_pack_STRING(Rcpp::StringVector proto, Rcpp::StringVector alphabet) {
  return importProtoFromR<STRING>(proto, alphabet)
  .pack<RCPP>()
  .exportToR();
}