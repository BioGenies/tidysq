#include <Rcpp.h>

#include "tidysq/tidysq.h"

using namespace tidysq;

//[[Rcpp::export]]
Rcpp::List CPP_pack_RAWS(Rcpp::List proto, Rcpp::StringVector alphabet) {
  return SqProto<RCPP, RAWS>(proto, Alphabet(alphabet))
  .pack<RCPP>()
  .exportToR();
}

//[[Rcpp::export]]
Rcpp::List CPP_pack_INTS(Rcpp::List proto, Rcpp::StringVector alphabet) {
  return SqProto<RCPP, INTS>(proto, Alphabet(alphabet))
  .pack<RCPP>()
  .exportToR();
}

//[[Rcpp::export]]
Rcpp::List CPP_pack_STRINGS(Rcpp::List proto, Rcpp::StringVector alphabet) {
  return SqProto<RCPP, STRINGS>(proto, Alphabet(alphabet))
  .pack<RCPP>()
  .exportToR();
}

//[[Rcpp::export]]
Rcpp::List CPP_pack_STRING(Rcpp::StringVector proto, Rcpp::StringVector alphabet) {
  return SqProto<RCPP, STRING>(proto, Alphabet(alphabet))
  .pack<RCPP>()
  .exportToR();
}