#include <Rcpp.h>
#include "tidysq/types/SqRCPP.h"
#include "tidysq/types/SqProtoRCPP.h"
#include "tidysq/types/AlphabetRCPP.h"

using namespace Rcpp;
using namespace tidysq;

//' @export
//[[Rcpp::export]]
List tmpPack(List raws, StringVector alphabet) {
  return SqProto<RCPP>(raws, Alphabet<RCPP>(alphabet))
    .pack<RCPP>()
    .exportToR();
}

//' @export
//[[Rcpp::export]]
List tmpPack2(StringVector raws, StringVector alphabet) {
  return SqProto<RCPP, STRING>(raws, Alphabet<RCPP>(alphabet))
  .pack<RCPP>()
  .exportToR();
}

//' @export
//[[Rcpp::export]]
List tmpUnpack(List raws) {
  return Sq<RCPP>(raws)
    .unpack<RCPP, RAWS>()
    .exportToR();
}

//' @export
//[[Rcpp::export]]
List tmpUnpack2(List raws) {
  return Sq<RCPP>(raws)
  .unpack<RCPP, INTS>()
  .exportToR();
}


//' @export
//[[Rcpp::export]]
List tmpUnpack3(List raws) {
  return Sq<RCPP>(raws)
  .unpack<RCPP, STRINGS>()
  .exportToR();
}