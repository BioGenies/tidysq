#include <Rcpp.h>
#include "tidysq/types/SqRCPP.h"
#include "tidysq/types/SqProtoRCPP.h"
#include "tidysq/types/AlphabetRCPP.h"

using namespace tidysq;

//' @export
//[[Rcpp::export]]
void CPP_unpack_RAWS(Rcpp::List sq) {
  Sq<RCPP> a(sq);
  //return //Sq<RCPP>(sq)
  //.unpack<RCPP, RAWS>()
  //.exportToR();
}

//' @export
//[[Rcpp::export]]
Rcpp::List CPP_unpack_INTS(Rcpp::List sq) {
  return Sq<RCPP>(sq)
  .unpack<RCPP, INTS>()
  .exportToR();
}

//' @export
//[[Rcpp::export]]
Rcpp::List CPP_unpack_STRINGS(Rcpp::List sq) {
  return Sq<RCPP>(sq)
  .unpack<RCPP, STRINGS>()
  .exportToR();
}

//' @export
//[[Rcpp::export]]
Rcpp::StringVector CPP_unpack_STRING(Rcpp::List sq) {
  return Sq<RCPP>(sq)
  .unpack<RCPP, STRING>()
  .exportToR();
}