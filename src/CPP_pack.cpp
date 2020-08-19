#include <Rcpp.h>
#include "tidysq/types/SqRCPP.h"
#include "tidysq/types/SqProtoRCPP.h"
#include "tidysq/types/AlphabetRCPP.h"

using namespace tidysq;

//' @export
//[[Rcpp::export]]
Rcpp::List CPP_pack_RAWS(Rcpp::List proto, Rcpp::StringVector alphabet) {
  return SqProto<RCPP, RAWS>(proto, Alphabet<RCPP>(alphabet))
  .pack<RCPP>()
  .exportToR();
}

//' @export
//[[Rcpp::export]]
Rcpp::List CPP_pack_INTS(Rcpp::List proto, Rcpp::StringVector alphabet) {
  return SqProto<RCPP, INTS>(proto, Alphabet<RCPP>(alphabet))
  .pack<RCPP>()
  .exportToR();
}

//' @export
//[[Rcpp::export]]
Rcpp::List CPP_pack_STRINGS(Rcpp::List proto, Rcpp::StringVector alphabet) {
  return SqProto<RCPP, STRINGS>(proto, Alphabet<RCPP>(alphabet))
  .pack<RCPP>()
  .exportToR();
}

//' @export
//[[Rcpp::export]]
Rcpp::List CPP_pack_STRING(Rcpp::StringVector proto, Rcpp::StringVector alphabet) {
  return SqProto<RCPP, STRING>(proto, Alphabet<RCPP>(alphabet))
  .pack<RCPP>()
  .exportToR();
}